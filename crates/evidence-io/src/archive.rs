//! Column-oriented archive encoder.
//!
//! Fields are written to separate streams and compressed independently. Each stream is more
//! homogeneous than an interleaved record, and the resulting per-stream sizes show which fields
//! account for the encoded bytes.
//!
//! Two structural choices account for much of the compression:
//!
//! * **Shape interning.** A molecule's geometry (block lengths, junction offsets) is stored
//!   relative to its start and interned. 10x 3' reads are overwhelmingly a single ungapped block
//!   of one length, so the shape stream should collapse to near nothing.
//! * **Positional delta coding.** Molecules are sorted by (chrom, start) so positions become small
//!   gaps rather than 32-bit absolutes.
//!
//! Whether the UMI stream is emitted is a parameter rather than a constant. The UMI is over half
//! the record, and it is needed only when a future annotation must redo per-gene UMI collapse.

use crate::{Molecule, Strand};
use anyhow::Result;
use rustc_hash::FxHashMap;
use std::io::Write;

/// Geometry of a molecule's placement, relative to its start position.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Shape {
    /// `(offset_from_start, len)` for each aligned block.
    pub blocks: Vec<(u32, u32)>,
}

impl Shape {
    pub fn of(p: &crate::Placement) -> Shape {
        let s = p.start();
        Shape {
            blocks: p.blocks.iter().map(|b| (b.start - s, b.len())).collect(),
        }
    }
}

/// Append a LEB128 varint.
pub fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let byte = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(byte);
            return;
        }
        out.push(byte | 0x80);
    }
}

/// Zigzag so small negative deltas stay small.
pub fn put_svarint(out: &mut Vec<u8>, v: i64) {
    put_varint(out, ((v << 1) ^ (v >> 63)) as u64);
}

#[derive(Default)]
struct Streams {
    chrom: Vec<u8>,
    pos_delta: Vec<u8>,
    shape_id: Vec<u8>,
    flags: Vec<u8>,
    cb_delta: Vec<u8>,
    umi: Vec<u8>,
    n_reads: Vec<u8>,
    shape_dict: Vec<u8>,
}

/// Compressed size of each stream, in bytes.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct SizeReport {
    pub molecules: u64,
    pub chrom: u64,
    pub pos_delta: u64,
    pub shape_id: u64,
    pub flags: u64,
    pub cb_delta: u64,
    pub umi: u64,
    pub n_reads: u64,
    pub shape_dict: u64,
    pub distinct_shapes: u64,
    pub total: u64,
}

impl SizeReport {
    pub fn bits_per_molecule(&self) -> f64 {
        if self.molecules == 0 {
            0.0
        } else {
            8.0 * self.total as f64 / self.molecules as f64
        }
    }
}

fn zstd_len(buf: &[u8], level: i32) -> Result<u64> {
    let mut enc = zstd::Encoder::new(Vec::new(), level)?;
    enc.write_all(buf)?;
    Ok(enc.finish()?.len() as u64)
}

/// Encode molecules and report the compressed size of each stream.
///
/// `molecules` is sorted in place by (chrom, start, cb) so the delta streams are meaningful;
/// callers that need the original order should clone first.
pub fn encode(molecules: &mut [Molecule], with_umi: bool, level: i32) -> Result<SizeReport> {
    molecules.sort_by_key(|m| {
        let p = m.placements.first();
        (
            p.map(|p| p.chrom).unwrap_or(u32::MAX),
            p.map(|p| p.start()).unwrap_or(0),
            m.cb,
        )
    });

    let mut s = Streams::default();
    let mut shapes: FxHashMap<Shape, u32> = FxHashMap::default();
    let mut shape_order: Vec<Shape> = Vec::new();

    let (mut last_chrom, mut last_pos, mut last_cb) = (u32::MAX, 0u32, 0u32);

    for m in molecules.iter() {
        let Some(p) = m.placements.first() else { continue };

        // Position deltas reset at each chromosome boundary, so a new chromosome does not emit one
        // enormous negative gap.
        if p.chrom != last_chrom {
            put_varint(&mut s.chrom, p.chrom as u64);
            last_chrom = p.chrom;
            last_pos = 0;
        } else {
            put_varint(&mut s.chrom, 0);
        }
        put_svarint(&mut s.pos_delta, p.start() as i64 - last_pos as i64);
        last_pos = p.start();

        let shape = Shape::of(p);
        let id = match shapes.get(&shape) {
            Some(&id) => id,
            None => {
                let id = shape_order.len() as u32;
                shapes.insert(shape.clone(), id);
                shape_order.push(shape);
                id
            }
        };
        put_varint(&mut s.shape_id, id as u64);

        // Strand in bit 0; bit 1 records that this molecule kept more than one placement, which
        // replay must know about because ambiguity is preserved rather than resolved at ingest.
        let f = (matches!(p.strand, Strand::Reverse) as u8) | ((m.placements.len() > 1) as u8) << 1;
        s.flags.push(f);

        put_svarint(&mut s.cb_delta, m.cb as i64 - last_cb as i64);
        last_cb = m.cb;

        if with_umi {
            s.umi.extend_from_slice(&m.umi.to_le_bytes());
        }
        put_varint(&mut s.n_reads, m.n_reads as u64);
    }

    for shape in &shape_order {
        put_varint(&mut s.shape_dict, shape.blocks.len() as u64);
        for (off, len) in &shape.blocks {
            put_varint(&mut s.shape_dict, *off as u64);
            put_varint(&mut s.shape_dict, *len as u64);
        }
    }

    let mut r = SizeReport {
        molecules: molecules.len() as u64,
        chrom: zstd_len(&s.chrom, level)?,
        pos_delta: zstd_len(&s.pos_delta, level)?,
        shape_id: zstd_len(&s.shape_id, level)?,
        flags: zstd_len(&s.flags, level)?,
        cb_delta: zstd_len(&s.cb_delta, level)?,
        umi: if with_umi { zstd_len(&s.umi, level)? } else { 0 },
        n_reads: zstd_len(&s.n_reads, level)?,
        shape_dict: zstd_len(&s.shape_dict, level)?,
        distinct_shapes: shape_order.len() as u64,
        total: 0,
    };
    r.total = r.chrom + r.pos_delta + r.shape_id + r.flags + r.cb_delta + r.umi + r.n_reads + r.shape_dict;
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Block, Placement};

    fn mol(cb: u32, umi: u32, chrom: u32, start: u32) -> Molecule {
        Molecule {
            cb,
            umi,
            placements: vec![Placement {
                chrom,
                strand: Strand::Forward,
                blocks: vec![Block { start, end: start + 91 }],
                junctions: vec![],
                nm: 0,
                score: 0,
                nh: 1,
                clip: (0, 0),
            }],
            n_reads: 2,
            residual: None,
        }
    }

    #[test]
    fn varint_roundtrip_boundaries() {
        for v in [0u64, 1, 127, 128, 300, u32::MAX as u64] {
            let mut b = Vec::new();
            put_varint(&mut b, v);
            // Decode inline to confirm the encoding is self-terminating and exact.
            let (mut got, mut shift) = (0u64, 0);
            for byte in &b {
                got |= ((byte & 0x7f) as u64) << shift;
                shift += 7;
            }
            assert_eq!(got, v, "varint roundtrip failed for {v}");
        }
    }

    #[test]
    fn svarint_keeps_small_negatives_small() {
        let mut b = Vec::new();
        put_svarint(&mut b, -1);
        assert_eq!(b.len(), 1, "-1 must encode in one byte, not ten");
    }

    #[test]
    fn identical_geometry_interns_to_one_shape() {
        // The core structural bet: 10x 3' molecules share one geometry, so the shape stream should
        // collapse regardless of how many molecules there are.
        let mut ms: Vec<_> = (0..1000).map(|i| mol(i, i * 7, 0, 1000 + i * 100)).collect();
        let r = encode(&mut ms, true, 3).unwrap();
        assert_eq!(r.distinct_shapes, 1);
        assert_eq!(r.molecules, 1000);
    }

    #[test]
    fn dropping_umi_removes_the_umi_stream_and_shrinks_total() {
        let mut a: Vec<_> = (0..2000).map(|i| mol(i, i.wrapping_mul(2654435761), 0, 1000 + i * 50)).collect();
        let mut b = a.clone();
        let with = encode(&mut a, true, 3).unwrap();
        let without = encode(&mut b, false, 3).unwrap();
        assert!(with.umi > 0);
        assert_eq!(without.umi, 0);
        assert!(
            without.total < with.total,
            "dropping the UMI must shrink the archive: {} vs {}",
            without.total,
            with.total
        );
    }

    #[test]
    fn sorted_positions_beat_scattered_ones() {
        // Delta coding is only a win if input is position-sorted; encode() sorts, so a shuffled
        // input must compress identically to a pre-sorted one.
        let mut sorted: Vec<_> = (0..1000).map(|i| mol(i, i, 0, 1000 + i * 100)).collect();
        let mut shuffled = sorted.clone();
        shuffled.reverse();
        let a = encode(&mut sorted, true, 3).unwrap();
        let b = encode(&mut shuffled, true, 3).unwrap();
        assert_eq!(a.pos_delta, b.pos_delta);
    }

    #[test]
    fn chromosome_change_resets_position_delta() {
        let mut ms = vec![mol(0, 1, 0, 200_000_000), mol(1, 2, 1, 1000)];
        let r = encode(&mut ms, true, 3).unwrap();
        // Two molecules on different chromosomes must not cost a huge negative delta.
        assert!(r.pos_delta < 64, "pos_delta stream unexpectedly large: {}", r.pos_delta);
    }
}
