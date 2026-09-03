//! CIGAR → aligned blocks and splice junctions.
//!
//! Block semantics follow STAR's, which matters because replay must reproduce its gene assignment.
//! In STAR an alignment is a list of blocks whose separators are tagged by `canonSJ`: values `>= 0`
//! are real splice junctions, `-3` is a mate jump, and `-1`/`-2` are indels. The classification
//! code in `alignToTranscript` advances the transcript exon cursor **only** at `canonSJ >= 0`, and
//! `alignToTranscriptMinOverlap` explicitly walks across indel separators to "expand the block".
//!
//! So the assignment-relevant view is: **split blocks at `N` only.** A deletion extends the current
//! block (it consumes reference), an insertion is skipped (it does not), and soft/hard clips are
//! recorded separately as residual-sequence bookkeeping for E3.

use evidence_io::{Block, Junction};

/// A CIGAR operation, in the subset that appears in RNA-seq alignments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// `M`/`=`/`X`: consumes query and reference.
    Match(u32),
    /// `I`: consumes query only.
    Ins(u32),
    /// `D`: consumes reference only; absorbed into the surrounding block.
    Del(u32),
    /// `N`: consumes reference only; splits blocks and yields a junction.
    Skip(u32),
    /// `S`: consumes query only; residual sequence.
    SoftClip(u32),
    /// `H`: consumes neither.
    HardClip(u32),
    /// `P`: consumes neither.
    Pad(u32),
}

/// Decompose an alignment starting at 0-based `pos` into blocks, junctions, and (5', 3') clip
/// lengths.
pub fn blocks_and_junctions(pos: u32, ops: &[Op]) -> (Vec<Block>, Vec<Junction>, (u16, u16)) {
    let mut blocks = Vec::new();
    let mut junctions = Vec::new();
    let mut cur: Option<Block> = None;
    let mut refpos = pos;

    for op in ops {
        match *op {
            Op::Match(n) => {
                cur = Some(match cur {
                    Some(b) => Block { start: b.start, end: refpos + n },
                    None => Block { start: refpos, end: refpos + n },
                });
                refpos += n;
            }
            Op::Del(n) => {
                // Reference-consuming, so it extends the block rather than breaking it.
                if let Some(b) = cur {
                    cur = Some(Block { start: b.start, end: b.end + n });
                }
                refpos += n;
            }
            Op::Skip(n) => {
                if let Some(b) = cur.take() {
                    blocks.push(b);
                }
                junctions.push(Junction { donor: refpos, acceptor: refpos + n });
                refpos += n;
            }
            // Query-only or no-op: no effect on genomic block structure.
            Op::Ins(_) | Op::SoftClip(_) | Op::HardClip(_) | Op::Pad(_) => {}
        }
    }
    if let Some(b) = cur {
        blocks.push(b);
    }

    let lead = leading_clip(ops.iter());
    let trail = leading_clip(ops.iter().rev());
    (blocks, junctions, (lead, trail))
}

/// Clip length at whichever end the iterator starts from. Hard clips sit outside soft clips, so
/// both are summed until the first aligned operation.
fn leading_clip<'a, I: Iterator<Item = &'a Op>>(ops: I) -> u16 {
    let mut n = 0u32;
    for op in ops {
        match *op {
            Op::SoftClip(k) | Op::HardClip(k) => n += k,
            Op::Pad(_) => {}
            _ => break,
        }
    }
    n.min(u16::MAX as u32) as u16
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ungapped_match_is_one_block() {
        let (b, j, c) = blocks_and_junctions(100, &[Op::Match(91)]);
        assert_eq!(b, vec![Block { start: 100, end: 191 }]);
        assert!(j.is_empty());
        assert_eq!(c, (0, 0));
    }

    #[test]
    fn skip_splits_blocks_and_yields_junction() {
        // 50M 150N 41M — the canonical spliced 10x read.
        let (b, j, _) = blocks_and_junctions(100, &[Op::Match(50), Op::Skip(150), Op::Match(41)]);
        assert_eq!(
            b,
            vec![Block { start: 100, end: 150 }, Block { start: 300, end: 341 }]
        );
        assert_eq!(j, vec![Junction { donor: 150, acceptor: 300 }]);
    }

    #[test]
    fn deletion_extends_block_rather_than_splitting_it() {
        // STAR walks across indel separators when matching exons, so a D must not create a
        // second block; otherwise replay would spuriously demand an extra exon boundary.
        let (b, j, _) = blocks_and_junctions(100, &[Op::Match(20), Op::Del(3), Op::Match(20)]);
        assert_eq!(b, vec![Block { start: 100, end: 143 }]);
        assert!(j.is_empty());
    }

    #[test]
    fn insertion_does_not_consume_reference() {
        let (b, _, _) = blocks_and_junctions(100, &[Op::Match(20), Op::Ins(5), Op::Match(20)]);
        assert_eq!(b, vec![Block { start: 100, end: 140 }]);
    }

    #[test]
    fn soft_clips_are_recorded_and_do_not_shift_coordinates() {
        let (b, _, c) =
            blocks_and_junctions(100, &[Op::SoftClip(7), Op::Match(80), Op::SoftClip(4)]);
        assert_eq!(b, vec![Block { start: 100, end: 180 }]);
        assert_eq!(c, (7, 4));
    }

    #[test]
    fn hard_and_soft_clips_sum_at_each_end() {
        let (_, _, c) = blocks_and_junctions(
            100,
            &[Op::HardClip(5), Op::SoftClip(3), Op::Match(80), Op::SoftClip(2), Op::HardClip(1)],
        );
        assert_eq!(c, (8, 3));
    }

    #[test]
    fn multiple_junctions_are_ordered() {
        let (b, j, _) = blocks_and_junctions(
            10,
            &[Op::Match(10), Op::Skip(90), Op::Match(10), Op::Skip(90), Op::Match(10)],
        );
        assert_eq!(b.len(), 3);
        assert_eq!(
            j,
            vec![
                Junction { donor: 20, acceptor: 110 },
                Junction { donor: 120, acceptor: 210 },
            ]
        );
    }

    #[test]
    fn deletion_adjacent_to_junction_keeps_junction_coordinates_exact() {
        // A D immediately before an N must not shift the donor site.
        let (b, j, _) = blocks_and_junctions(
            100,
            &[Op::Match(20), Op::Del(2), Op::Skip(100), Op::Match(30)],
        );
        assert_eq!(b[0], Block { start: 100, end: 122 });
        assert_eq!(j, vec![Junction { donor: 122, acceptor: 222 }]);
    }
}
