//! Extraction of annotation-independent evidence from BAM alignments.

pub mod cigar;

use evidence_io::{Placement, Strand};

/// Build a [`Placement`] from the fields a BAM record exposes.
///
/// Tag values are passed in rather than read here so this stays testable without constructing BAM
/// records.
pub fn placement_from_alignment(
    chrom: u32,
    pos: u32,
    reverse: bool,
    ops: &[cigar::Op],
    nm: u16,
    score: i32,
    nh: u16,
) -> Placement {
    let (blocks, junctions, clip) = cigar::blocks_and_junctions(pos, ops);
    Placement {
        chrom,
        strand: if reverse { Strand::Reverse } else { Strand::Forward },
        blocks,
        junctions,
        nm,
        score,
        nh,
        clip,
    }
}

/// Whether two placements would lead STAR to the same gene assignment.
///
/// Compares only the fields `alignToTranscript` and `geneCountsAddAlign` actually read. This
/// determines whether annotation-free and annotation-aware alignments carry equivalent evidence
/// for gene assignment.
pub fn same_assignment_evidence(a: &Placement, b: &Placement) -> bool {
    a.chrom == b.chrom
        && a.strand == b.strand
        && a.nh == b.nh
        && a.blocks == b.blocks
        && a.junctions == b.junctions
}

/// How two placements differ. Every disagreement receives a concrete category rather than being
/// reported only as part of an aggregate count.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Divergence {
    Identical,
    /// Different chromosome or strand: a wholly different locus.
    Locus,
    /// Same locus and junctions, but block boundaries shifted (e.g. trimmed or extended ends).
    BlockBoundary,
    /// Splice structure differs: the case sjdb insertion is expected to cause.
    Junction,
    /// Multimapping status changed, which flips assignability regardless of coordinates.
    Multiplicity,
    /// Aligned in one configuration and absent in the other.
    Presence,
}

pub fn classify(a: Option<&Placement>, b: Option<&Placement>) -> Divergence {
    let (a, b) = match (a, b) {
        (Some(a), Some(b)) => (a, b),
        (None, None) => return Divergence::Identical,
        _ => return Divergence::Presence,
    };
    if same_assignment_evidence(a, b) {
        return Divergence::Identical;
    }
    if a.chrom != b.chrom || a.strand != b.strand {
        return Divergence::Locus;
    }
    // Junctions are checked before multiplicity: a changed splice structure is the mechanism we
    // care about, and it frequently changes NH as a side effect.
    if a.junctions != b.junctions {
        return Divergence::Junction;
    }
    if a.nh != b.nh {
        return Divergence::Multiplicity;
    }
    Divergence::BlockBoundary
}

#[cfg(test)]
mod tests {
    use super::*;
    use evidence_io::{Block, Junction};

    fn p(blocks: Vec<Block>, junctions: Vec<Junction>, nh: u16) -> Placement {
        Placement {
            chrom: 0,
            strand: Strand::Forward,
            blocks,
            junctions,
            nm: 0,
            score: 0,
            nh,
            clip: (0, 0),
        }
    }

    #[test]
    fn identical_evidence_ignores_score_and_edits() {
        let mut a = p(vec![Block { start: 10, end: 60 }], vec![], 1);
        let mut b = a.clone();
        a.nm = 5;
        a.score = 90;
        b.nm = 0;
        b.score = 12;
        // STAR reads neither field when assigning genes, so these must not count as divergence.
        assert!(same_assignment_evidence(&a, &b));
        assert_eq!(classify(Some(&a), Some(&b)), Divergence::Identical);
    }

    #[test]
    fn junction_change_is_detected() {
        let a = p(
            vec![Block { start: 10, end: 60 }, Block { start: 200, end: 240 }],
            vec![Junction { donor: 60, acceptor: 200 }],
            1,
        );
        let b = p(
            vec![Block { start: 10, end: 60 }, Block { start: 300, end: 340 }],
            vec![Junction { donor: 60, acceptor: 300 }],
            1,
        );
        assert_eq!(classify(Some(&a), Some(&b)), Divergence::Junction);
    }

    #[test]
    fn multiplicity_change_is_detected() {
        let a = p(vec![Block { start: 10, end: 60 }], vec![], 1);
        let b = p(vec![Block { start: 10, end: 60 }], vec![], 3);
        assert_eq!(classify(Some(&a), Some(&b)), Divergence::Multiplicity);
    }

    #[test]
    fn block_boundary_shift_is_detected() {
        let a = p(vec![Block { start: 10, end: 60 }], vec![], 1);
        let b = p(vec![Block { start: 12, end: 60 }], vec![], 1);
        assert_eq!(classify(Some(&a), Some(&b)), Divergence::BlockBoundary);
    }

    #[test]
    fn missing_on_one_side_is_presence() {
        let a = p(vec![Block { start: 10, end: 60 }], vec![], 1);
        assert_eq!(classify(Some(&a), None), Divergence::Presence);
        assert_eq!(classify(None, None), Divergence::Identical);
    }

    #[test]
    fn different_chromosome_is_locus() {
        let a = p(vec![Block { start: 10, end: 60 }], vec![], 1);
        let mut b = a.clone();
        b.chrom = 7;
        assert_eq!(classify(Some(&a), Some(&b)), Divergence::Locus);
    }
}
