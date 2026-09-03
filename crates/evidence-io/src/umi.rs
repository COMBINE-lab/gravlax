//! 2-bit UMI packing. Current 10x chemistries use fixed 10 or 12 bp UMIs, so a `u32` holds either
//! exactly (and supports any ACGT-only sequence up to 16 bp).
//!
//! `N` bases cannot be represented; sequences containing them are rejected rather than silently
//! coerced, because a coerced UMI would collapse into the wrong molecule.

/// Pack an ACGT-only sequence of at most 16 bases into a `u32`. Returns `None` on any other base.
pub fn pack(seq: &[u8]) -> Option<u32> {
    if seq.len() > 16 {
        return None;
    }
    let mut v: u32 = 0;
    for &b in seq {
        let code = match b {
            b'A' | b'a' => 0u32,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        v = (v << 2) | code;
    }
    Some(v)
}

/// Unpack `len` bases previously packed by [`pack`].
pub fn unpack(mut v: u32, len: usize) -> Vec<u8> {
    let mut out = vec![0u8; len];
    for i in (0..len).rev() {
        out[i] = b"ACGT"[(v & 0b11) as usize];
        v >>= 2;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_12bp() {
        let s = b"ACGTACGTACGT";
        let p = pack(s).unwrap();
        assert_eq!(unpack(p, 12), s);
    }

    #[test]
    fn roundtrip_10bp_for_5prime_v2() {
        let s = b"ACGTACGTAC";
        let p = pack(s).unwrap();
        assert_eq!(unpack(p, 10), s);
    }

    #[test]
    fn distinct_umis_pack_distinctly() {
        assert_ne!(pack(b"AAAAAAAAAAAA"), pack(b"AAAAAAAAAAAC"));
    }

    #[test]
    fn leading_a_is_not_lost() {
        // A packs to 0, so a naive implementation loses leading As; length is carried separately.
        assert_eq!(unpack(pack(b"AAAC").unwrap(), 4), b"AAAC");
    }

    #[test]
    fn rejects_n_rather_than_coercing() {
        assert_eq!(pack(b"ACGTNACGTACG"), None);
    }

    #[test]
    fn rejects_overlong() {
        assert_eq!(pack(b"ACGTACGTACGTACGTA"), None);
    }
}
