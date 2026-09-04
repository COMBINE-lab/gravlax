//! Static-table rANS for heavy-tailed unsigned value streams.
//!
//! On noisy chunk streams (class backreferences, weights, and representative/multimapper offsets),
//! zstd encoded about 5–25% above the order-0 entropy of the decoded values: match modeling finds
//! little repetition, and varint byte boundaries hide the value distribution. A memoryless coder
//! using the observed value distribution closes that gap. Streams where zstd beats order-0
//! entropy (anchors, layouts, shapes, and patterns) remain varint+zstd.
//!
//! Symbolization: values < 128 are their own symbol; larger values escape to a bit-length class
//! (8..=64) with the low `len-1` bits carried verbatim in a separate bitstream. Frequencies are
//! normalized to 12 bits; the coder is single-state byte-wise rANS (L = 2^23).

use anyhow::{bail, Result};

pub const NSYM: usize = 128 + 57; // direct 0..127, escapes for bit lengths 8..=64
const SCALE_BITS: u32 = 12;
const SCALE: u32 = 1 << SCALE_BITS;
const LOW: u64 = 1 << 23;
// Chunks are 4 Mb genomic bins and are expected to stay far below this. A corrupt `n` must not
// turn a tiny payload into an unbounded allocation or decode loop.
pub const MAX_VALUES_PER_STREAM: usize = 1 << 25;

fn sym_of(v: u64) -> (usize, u32, u64) {
    // (symbol, extra bit count, extra bits)
    if v < 128 {
        (v as usize, 0, 0)
    } else {
        let len = 64 - v.leading_zeros(); // >= 8
        (128 + (len - 8) as usize, len - 1, v - (1u64 << (len - 1)))
    }
}

fn val_of(sym: usize, extra: u64) -> u64 {
    if sym < 128 {
        sym as u64
    } else {
        let len = sym as u32 - 128 + 8;
        (1u64 << (len - 1)) + extra
    }
}

pub struct Table {
    pub freq: [u32; NSYM],
    cum: [u32; NSYM + 1],
    lookup: Vec<u16>, // SCALE entries -> symbol
}

impl Table {
    /// Normalize observed counts to SCALE, guaranteeing every present symbol a nonzero slot.
    /// An empty stream receives a canonical all-zero-symbol table; the table is serialized with
    /// the archive but is never consulted while decoding its zero values.
    pub fn from_counts(counts: &[u64; NSYM]) -> Result<Table> {
        let total: u64 = counts.iter().sum();
        if total == 0 {
            let mut freq = [0u32; NSYM];
            freq[0] = SCALE;
            return Ok(Table::from_freqs(freq));
        }
        let mut freq = [0u32; NSYM];
        let mut assigned = 0u32;
        for i in 0..NSYM {
            if counts[i] > 0 {
                freq[i] = (((counts[i] as u128) * SCALE as u128 / total as u128) as u32).max(1);
                assigned += freq[i];
            }
        }
        // Repair to exactly SCALE by charging/crediting the largest bucket.
        let imax = (0..NSYM).max_by_key(|&i| freq[i]).unwrap();
        let diff = SCALE as i64 - assigned as i64;
        let nf = freq[imax] as i64 + diff;
        if nf < 1 {
            bail!("cannot normalize rans table");
        }
        freq[imax] = nf as u32;
        Ok(Table::from_freqs(freq))
    }

    pub fn from_freqs(freq: [u32; NSYM]) -> Table {
        let mut cum = [0u32; NSYM + 1];
        for i in 0..NSYM {
            cum[i + 1] = cum[i] + freq[i];
        }
        let mut lookup = vec![0u16; SCALE as usize];
        for s in 0..NSYM {
            for slot in cum[s]..cum[s + 1] {
                lookup[slot as usize] = s as u16;
            }
        }
        Table { freq, cum, lookup }
    }

    pub fn serialize(&self, out: &mut Vec<u8>) {
        for f in self.freq {
            crate::archive::put_varint(out, f as u64);
        }
    }

    pub fn deserialize(c: &mut crate::format::Cursor) -> Result<Table> {
        let mut freq = [0u32; NSYM];
        let mut sum = 0u64;
        for f in freq.iter_mut() {
            let value = c.varint()?;
            if value > SCALE as u64 {
                bail!("rans frequency {value} exceeds scale {SCALE}");
            }
            *f = value as u32;
            sum += value;
        }
        if sum != SCALE as u64 {
            bail!("rans table sums to {sum}, expected {SCALE}");
        }
        Ok(Table::from_freqs(freq))
    }
}

/// Count symbols for table building.
pub fn count(values: &[u64], counts: &mut [u64; NSYM]) {
    for &v in values {
        counts[sym_of(v).0] += 1;
    }
}

/// Encode: payload = [n varint][sym bytes len varint][sym bytes][extra bitstream].
pub fn encode(values: &[u64], t: &Table, out: &mut Vec<u8>) {
    crate::archive::put_varint(out, values.len() as u64);
    // Symbols encode in reverse so decode streams forward.
    let mut body: Vec<u8> = Vec::with_capacity(values.len());
    let mut x: u64 = LOW;
    let mut extras_bits: Vec<u8> = Vec::new(); // bit-packed, MSB first
    let (mut bitbuf, mut nbits) = (0u64, 0u32);
    for &v in values {
        let (_, nb, eb) = sym_of(v);
        if nb > 0 {
            bitbuf = (bitbuf << nb) | eb;
            nbits += nb;
            while nbits >= 8 {
                extras_bits.push((bitbuf >> (nbits - 8)) as u8);
                nbits -= 8;
            }
        }
    }
    if nbits > 0 {
        extras_bits.push((bitbuf << (8 - nbits)) as u8);
    }
    for &v in values.iter().rev() {
        let (s, _, _) = sym_of(v);
        let f = t.freq[s] as u64;
        let x_max = ((LOW >> SCALE_BITS) << 8) * f;
        while x >= x_max {
            body.push(x as u8);
            x >>= 8;
        }
        x = ((x / f) << SCALE_BITS) + (x % f) + t.cum[s] as u64;
    }
    crate::archive::put_varint(out, (body.len() + 8) as u64);
    out.extend_from_slice(&x.to_le_bytes());
    // body was emitted low-byte-first during reverse encode; decoder pops from the end.
    out.extend_from_slice(&body);
    out.extend_from_slice(&extras_bits);
}

pub fn decode(payload: &[u8], t: &Table) -> Result<Vec<u64>> {
    decode_limited(payload, t, MAX_VALUES_PER_STREAM)
}

/// Decode while enforcing a caller-supplied upper bound from surrounding format metadata.
pub fn decode_limited(payload: &[u8], t: &Table, max_values: usize) -> Result<Vec<u64>> {
    let mut c = crate::format::Cursor::new(payload);
    let n = usize::try_from(c.varint()?).map_err(|_| anyhow::anyhow!("rans value count is too large"))?;
    if n > max_values {
        bail!("rans value count {n} exceeds caller limit {max_values}");
    }
    let sym_len = usize::try_from(c.varint()?).map_err(|_| anyhow::anyhow!("rans symbol stream is too large"))?;
    let start = c.position();
    let Some(sym_end) = start.checked_add(sym_len) else { bail!("rans payload length overflow") };
    if sym_end > payload.len() || sym_len < 8 {
        bail!("rans payload truncated");
    }
    let mut x = u64::from_le_bytes(payload[start..start + 8].try_into().unwrap());
    let mut body = &payload[start + 8..sym_end];
    let extras = &payload[sym_end..];
    let (mut bitpos, mut out) = (0usize, Vec::with_capacity(n));
    let mut syms = Vec::with_capacity(n);
    for _ in 0..n {
        let slot = (x & (SCALE as u64 - 1)) as usize;
        let s = t.lookup[slot] as usize;
        let f = t.freq[s] as u64;
        x = f * (x >> SCALE_BITS) + (slot as u64) - t.cum[s] as u64;
        while x < LOW {
            let Some((&b, rest)) = body.split_last() else { bail!("rans body underrun") };
            body = rest;
            x = (x << 8) | b as u64;
        }
        syms.push(s);
    }
    for s in syms {
        if s < 128 {
            out.push(s as u64);
        } else {
            // Byte-at-a-time MSB-first bit gather (bit-identical to the old per-bit loop, ~5x
            // fewer iterations — escapes carry 7..=63 bits and dominate the position streams).
            let mut nb = (s as u32 - 128 + 8 - 1) as usize;
            let mut eb = 0u64;
            while nb > 0 {
                let Some(&byte) = extras.get(bitpos >> 3) else {
                    bail!("rans extra-bit stream underrun");
                };
                let avail = 8 - (bitpos & 7);
                let take = avail.min(nb);
                let chunk = (byte >> (avail - take)) & (((1u16 << take) - 1) as u8);
                eb = (eb << take) | chunk as u64;
                bitpos += take;
                nb -= take;
            }
            out.push(val_of(s, eb));
        }
    }
    let expected_extra_bytes = bitpos.div_ceil(8);
    if extras.len() != expected_extra_bytes {
        bail!(
            "rans extra-bit stream has {} bytes, expected {expected_extra_bytes}",
            extras.len()
        );
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_mixed_values() {
        let mut vals: Vec<u64> = Vec::new();
        for i in 0..50_000u64 {
            vals.push(i % 7);
            if i % 3 == 0 {
                vals.push(1000 + (i * 37) % 90_000);
            }
            if i % 1000 == 0 {
                vals.push(u32::MAX as u64 + i);
            }
        }
        let mut counts = [0u64; NSYM];
        count(&vals, &mut counts);
        let t = Table::from_counts(&counts).unwrap();
        let mut payload = Vec::new();
        encode(&vals, &t, &mut payload);
        let back = decode(&payload, &t).unwrap();
        assert_eq!(back, vals);
    }

    #[test]
    fn table_roundtrip() {
        let mut counts = [0u64; NSYM];
        counts[0] = 1000;
        counts[5] = 30;
        counts[130] = 7;
        let t = Table::from_counts(&counts).unwrap();
        let mut buf = Vec::new();
        t.serialize(&mut buf);
        let t2 = Table::deserialize(&mut crate::format::Cursor::new(&buf)).unwrap();
        assert_eq!(t.freq, t2.freq);
    }

    #[test]
    fn empty_stream_encodes() {
        let counts = [0u64; NSYM];
        let t = Table::from_counts(&counts).unwrap();
        assert_eq!(t.freq[0], SCALE);
        assert!(t.freq[1..].iter().all(|frequency| *frequency == 0));
        let mut payload = Vec::new();
        encode(&[], &t, &mut payload);
        assert_eq!(decode(&payload, &t).unwrap(), Vec::<u64>::new());
    }

    #[test]
    fn truncated_extra_bits_are_rejected() {
        let vals = [1000u64, 90_000];
        let mut counts = [0u64; NSYM];
        count(&vals, &mut counts);
        let t = Table::from_counts(&counts).unwrap();
        let mut payload = Vec::new();
        encode(&vals, &t, &mut payload);
        payload.pop();
        assert!(decode(&payload, &t).unwrap_err().to_string().contains("extra-bit"));
    }

    #[test]
    fn trailing_extra_bytes_are_rejected() {
        let vals = [1000u64];
        let mut counts = [0u64; NSYM];
        count(&vals, &mut counts);
        let t = Table::from_counts(&counts).unwrap();
        let mut payload = Vec::new();
        encode(&vals, &t, &mut payload);
        payload.push(0);
        assert!(decode(&payload, &t).is_err());
    }
}
