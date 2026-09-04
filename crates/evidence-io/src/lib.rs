//! Core annotation-independent evidence types.
//!
//! The design constraint that shapes everything here: STAR assigns reads to genes using only
//! aligned blocks, splice junctions, strand and multimapping count (see
//! `Transcriptome_classifyAlign.cpp::alignToTranscript`). Nothing in the assignment path reads
//! sequence, quality or alignment score. Conditional on fixed genome placements and the specified
//! assignment consumer, the fields below at level [`Ev::E1`] retain the geometry that consumer
//! reads. The deployed archive applies an additional, empirically evaluated molecule-level
//! reduction. Richer levels exist for consumers whose decisions *do* depend on sequence and to
//! bound what a future annotation might need.

pub mod alignment_provenance;
pub mod archive;
pub mod format;
pub mod genome;
pub mod rans;
pub mod terminal_tail;
pub mod umi;

use serde::{Deserialize, Serialize};

/// Evidence richness levels. Higher levels are supersets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Ev {
    /// cell, UMI, genomic start/end, strand.
    E0,
    /// + aligned blocks and splice junctions.
    E1,
    /// + mismatch/indel summaries and variant witnesses.
    E2,
    /// + residual sequence for anything the genome does not explain.
    E3,
    /// near-lossless: everything needed to reproduce fresh mapping decisions.
    E4,
}

/// A half-open genomic block `[start, end)`, 0-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Block {
    pub start: u32,
    pub end: u32,
}

impl Block {
    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start
    }
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

/// An intron inferred from a CIGAR `N` operation: `[donor, acceptor)` in genomic coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Junction {
    pub donor: u32,
    pub acceptor: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

/// One genomic placement of a read or molecule. This is the unit that
/// [`replay`](../replay/index.html) tests against an annotation.
///
/// `blocks` are maximal aligned runs separated by junctions; CIGAR `D`/`I` do **not** split a
/// block, matching STAR's own notion (`canonSJ >= 0` marks a real junction, negative values mark
/// indels and mate jumps).
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Placement {
    pub chrom: u32,
    pub strand: Strand,
    pub blocks: Vec<Block>,
    pub junctions: Vec<Junction>,
    /// Edit distance (BAM `NM`), retained from E2 up.
    pub nm: u16,
    /// Alignment score (BAM `AS`), retained from E2 up.
    pub score: i32,
    /// Number of reported alignments for this read (BAM `NH`). STAR refuses to assign a gene to
    /// any read with `NH > 1` in the default `Gene` path, so this is assignment-relevant.
    pub nh: u16,
    /// Soft-clipped lengths at the 5' and 3' ends of the alignment.
    pub clip: (u16, u16),
}

impl Placement {
    #[inline]
    pub fn start(&self) -> u32 {
        self.blocks.first().map(|b| b.start).unwrap_or(0)
    }
    #[inline]
    pub fn end(&self) -> u32 {
        self.blocks.last().map(|b| b.end).unwrap_or(0)
    }
    #[inline]
    pub fn aligned_len(&self) -> u32 {
        self.blocks.iter().map(|b| b.len()).sum()
    }

    /// The tuple used to compare assignment-relevant evidence across mapping configurations.
    /// It deliberately excludes `nm`, `score`, and `clip` because gene assignment reads only
    /// the returned fields.
    pub fn assignment_key(&self) -> (u32, Strand, &[Block], &[Junction], u16) {
        (self.chrom, self.strand, &self.blocks, &self.junctions, self.nh)
    }
}

/// A molecule: a cell barcode, a UMI, and the evidence for one or more genomic placements.
///
/// Multiple placements are retained rather than resolved. Reads sharing `(cb, umi)` but with
/// incompatible placements stay as separate hypotheses; forcing an early merge would bake in a
/// decision that should remain available to the annotation supplied at replay time.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    /// Interned cell-barcode id.
    pub cb: u32,
    /// 2-bit packed UMI (12 bp for 10x 3' v3).
    pub umi: u32,
    pub placements: Vec<Placement>,
    /// Raw reads collapsed into this molecule. The storage thesis lives here: this is the factor
    /// by which a molecule archive beats a read archive.
    pub n_reads: u32,
    /// Residual sequence not explained by the genome (E3+): soft-clips, unmapped mates.
    pub residual: Option<Vec<u8>>,
}

impl Molecule {
    /// Project this record down to a lower evidence level, for the E0..E4 ablation. Returns a copy
    /// with the fields above `level` cleared, so archive-size and replay-accuracy can be measured
    /// on exactly the same molecules.
    pub fn project(&self, level: Ev) -> Molecule {
        let mut m = self.clone();
        if level < Ev::E3 {
            m.residual = None;
        }
        for p in &mut m.placements {
            if level < Ev::E2 {
                p.nm = 0;
                p.score = 0;
                p.clip = (0, 0);
            }
            if level < Ev::E1 {
                // E0 keeps only the span: collapse blocks and drop junctions. This is expected to
                // fail junction concordance in replay, which is the point of measuring it.
                let (s, e) = (p.start(), p.end());
                p.blocks = vec![Block { start: s, end: e }];
                p.junctions.clear();
            }
        }
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn placement() -> Placement {
        Placement {
            chrom: 0,
            strand: Strand::Forward,
            blocks: vec![
                Block { start: 100, end: 150 },
                Block { start: 300, end: 340 },
            ],
            junctions: vec![Junction { donor: 150, acceptor: 300 }],
            nm: 2,
            score: 88,
            nh: 1,
            clip: (3, 0),
        }
    }

    #[test]
    fn placement_geometry() {
        let p = placement();
        assert_eq!(p.start(), 100);
        assert_eq!(p.end(), 340);
        assert_eq!(p.aligned_len(), 90);
    }

    #[test]
    fn e0_collapses_blocks_and_drops_junctions() {
        let m = Molecule {
            cb: 1,
            umi: 2,
            placements: vec![placement()],
            n_reads: 4,
            residual: Some(vec![b'A']),
        };
        let e0 = m.project(Ev::E0);
        assert_eq!(e0.placements[0].blocks, vec![Block { start: 100, end: 340 }]);
        assert!(e0.placements[0].junctions.is_empty());
        assert!(e0.residual.is_none());

        // E1 must preserve the splice structure that E0 destroys.
        let e1 = m.project(Ev::E1);
        assert_eq!(e1.placements[0].blocks.len(), 2);
        assert_eq!(e1.placements[0].junctions.len(), 1);
        assert_eq!(e1.placements[0].nm, 0, "E1 must not retain edit info");
    }

    #[test]
    fn e4_projection_is_identity() {
        let m = Molecule {
            cb: 1,
            umi: 2,
            placements: vec![placement()],
            n_reads: 4,
            residual: Some(vec![b'A']),
        };
        assert_eq!(m.project(Ev::E4), m);
    }
}
