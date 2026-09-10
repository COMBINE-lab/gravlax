//! Port of STAR's alignment-vs-transcript classification and the STARsolo `Gene` assignment rule.
//!
//! The reference implementation is STAR's
//! `Transcriptome_classifyAlign.cpp::alignToTranscript`. That function consumes only the
//! alignment's blocks and junctions against the transcript's exon structure — no sequence or
//! alignment scores — which is why an E1 archive can drive it.
//!
//! Semantics ported:
//!  * every junction in the alignment must coincide exactly with an exon-exon boundary of the
//!    transcript ("SJ concordance"), or the transcript is incompatible;
//!  * blocks landing across an exon/intron boundary mark the ExonIntronSpan state;
//!  * otherwise the alignment is Concordant (fully exonic), Intron (fully intronic) or ExonIntron
//!    (mixed, no boundary spanned).
//!
//! Gravlax blocks merge across indels (a deletion extends a block), whereas STAR keeps
//! indel-separated blocks and walks across them. A deletion bridging an exon/intron boundary can
//! therefore classify differently.

use crate::Transcript;
use evidence_io::Placement;

/// Relationship between the cDNA alignment strand and the transcript strand, matching STARsolo's
/// `--soloStrand` values. `Forward` is the historical Gravlax behavior (10x 3' libraries);
/// R2-only 10x 5' libraries require `Reverse`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SoloStrand {
    Forward,
    Reverse,
    Unstranded,
}

impl SoloStrand {
    /// Whether an alignment/transcript strand pair is eligible for assignment.
    #[inline]
    pub fn accepts(self, alignment_rev: bool, transcript_rev: bool) -> bool {
        match self {
            Self::Forward => alignment_rev == transcript_rev,
            Self::Reverse => alignment_rev != transcript_rev,
            Self::Unstranded => true,
        }
    }
}

/// STARsolo `GeneFull`: overlap aligned blocks with exon-derived full gene spans.
/// Built once per replay; the archive and compiled annotation formats are unchanged.
pub struct GeneFullIndex {
    chroms: hashbrown::HashMap<u32, [Vec<GeneSpan>; 2]>,
}

struct GeneSpan {
    start: u32,
    end: u32,
    max_end: u32,
    gene: u32,
}

impl GeneFullIndex {
    pub fn new(anno: &crate::Annotation) -> Self {
        let mut bounds = hashbrown::HashMap::<(u32, u32, bool), (u32, u32)>::new();
        for t in &anno.transcripts {
            let (start, end) = t.span();
            if start >= end {
                continue;
            }
            let entry = bounds
                .entry((t.chrom, t.gene, t.strand_rev))
                .or_insert((start, end));
            entry.0 = entry.0.min(start);
            entry.1 = entry.1.max(end);
        }
        let mut chroms: hashbrown::HashMap<u32, [Vec<GeneSpan>; 2]> = hashbrown::HashMap::new();
        for ((chrom, gene, reverse), (start, end)) in bounds {
            chroms.entry(chrom).or_default()[usize::from(reverse)].push(GeneSpan {
                start,
                end,
                max_end: 0,
                gene,
            });
        }
        for strands in chroms.values_mut() {
            for spans in strands {
                spans.sort_unstable_by_key(|s| (s.start, s.end, s.gene));
                let mut max_end = 0;
                for span in spans {
                    max_end = max_end.max(span.end);
                    span.max_end = max_end;
                }
            }
        }
        Self { chroms }
    }

    /// Only aligned blocks overlap: genes inside skipped introns are not hit by
    /// the outer alignment span alone. Junction concordance is not required.
    pub fn genes_into(&self, p: &Placement, chrom: u32, strand: SoloStrand, out: &mut Vec<u32>) {
        self.genes_for_blocks_into(
            chrom,
            matches!(p.strand, evidence_io::Strand::Reverse),
            strand,
            p.blocks.iter().map(|b| (b.start, b.end)),
            out,
        );
    }

    /// Query blocks directly from retained shape offsets without materializing a Placement or
    /// its unused splice junctions. Strand-specific indexes avoid scanning opposite-strand genes.
    /// Coordinates are zero-based, half-open; output is a sorted, deduplicated gene set.
    pub fn genes_for_blocks_into(
        &self,
        chrom: u32,
        alignment_reverse: bool,
        strand: SoloStrand,
        blocks: impl IntoIterator<Item = (u32, u32)>,
        out: &mut Vec<u32>,
    ) {
        out.clear();
        let Some(strands) = self.chroms.get(&chrom) else {
            return;
        };
        let indexes: &[Vec<GeneSpan>] = match strand {
            SoloStrand::Forward => {
                &strands[usize::from(alignment_reverse)..=usize::from(alignment_reverse)]
            }
            SoloStrand::Reverse => {
                &strands[usize::from(!alignment_reverse)..=usize::from(!alignment_reverse)]
            }
            SoloStrand::Unstranded => strands,
        };
        for (start, end) in blocks {
            if start >= end {
                continue;
            }
            for spans in indexes {
                let hi = spans.partition_point(|s| s.start < end);
                for span in spans[..hi].iter().rev() {
                    if span.max_end <= start {
                        break;
                    }
                    if span.end > start {
                        out.push(span.gene);
                    }
                }
            }
        }
        out.sort_unstable();
        out.dedup();
    }
}

/// STAR's `AlignVsTranscript` states, minus the transcript-distance bookkeeping replay never uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Vs {
    /// Fully exonic with concordant junctions — the state the `Gene` feature counts.
    Concordant,
    /// Mixed exonic/intronic without spanning a boundary (retained intron read).
    ExonIntron,
    /// Fully intronic.
    Intron,
    /// A block crosses an exon/intron boundary.
    ExonIntronSpan,
    /// Outside the transcript, or a junction disagrees with the exon structure.
    Incompatible,
}

/// Classify one placement against one transcript, mirroring `alignToTranscript`.
///
/// Both inner walks use binary search over the (disjoint, sorted — starts and ends strictly
/// increase after the parser's merge) exon array: the old linear rescans made this O(exons) per
/// segment, quadratic on many-exon transcripts, and this function dominates the replay profile.
pub fn align_vs_transcript(p: &Placement, t: &Transcript) -> Vs {
    let (ts, te) = t.span();
    if p.start() < ts || p.end() > te {
        return Vs::Incompatible;
    }

    // Junction concordance: each alignment junction must be exactly an intron of the transcript,
    // i.e. donor == some exon end and acceptor == the next exon's start. Exon ends are strictly
    // increasing, so at most one window can match — find it directly.
    for j in &p.junctions {
        let k = t.exons.partition_point(|e| e.end < j.donor);
        let ok = k + 1 < t.exons.len()
            && t.exons[k].end == j.donor
            && t.exons[k + 1].start == j.acceptor;
        if !ok {
            return Vs::Incompatible;
        }
    }

    let (mut exonic, mut intronic, mut span) = (false, false, false);
    for b in &p.blocks {
        // Walk the block across the transcript's exon/intron structure.
        let mut pos = b.start;
        while pos < b.end {
            // Segment containing `pos`: the exon with the greatest start <= pos if pos is inside
            // it, else the intron up to the next exon (or the transcript end) — exactly the
            // first-hit semantics of the old left-to-right scan.
            let k = t.exons.partition_point(|e| e.start <= pos);
            let (seg_end, in_exon) = if k > 0 && pos < t.exons[k - 1].end {
                (t.exons[k - 1].end, true)
            } else if k < t.exons.len() {
                (t.exons[k].start, false)
            } else {
                (te, false)
            };
            if in_exon {
                exonic = true;
            } else {
                intronic = true;
            }
            if seg_end < b.end {
                span = true; // block continues past this segment boundary
            }
            pos = seg_end.max(pos + 1);
        }
    }

    if span {
        Vs::ExonIntronSpan
    } else if exonic && !intronic {
        Vs::Concordant
    } else if exonic {
        Vs::ExonIntron
    } else {
        Vs::Intron
    }
}

/// Velocyto type bits, matching STAR's `AlignVsTranscript` bit positions exactly
/// (`Intron=0, ExonIntron=1, ExonIntronSpan=2, Concordant=3`).
pub const VB_INTRON: u8 = 1 << 0;
pub const VB_EXONINTRON: u8 = 1 << 1;
pub const VB_SPAN: u8 = 1 << 2;
pub const VB_CONCORDANT: u8 = 1 << 3;

/// Port of `alignToTranscriptMinOverlap` (minOverlapMinusOne = 6, matching velocyto's
/// MIN_FLANK = 5): the velocyto-mode classification. Blocks shorter than the tolerance make no
/// call; near-boundary starts/ends make no call; a block resting in an intron longer than 1 Mb
/// disqualifies the transcript entirely (velocyto.py's giant-intron rule); a spliced alignment
/// that touches intron or spans a boundary is incompatible — junction concordance is implied by
/// every block resting cleanly inside the exon structure. Returns None for incompatible.
pub fn align_vs_transcript_min_overlap(p: &Placement, t: &Transcript) -> Option<Vs> {
    const TOL: u32 = 6;
    let (ts, te) = t.span();
    if p.start() < ts || p.end() > te {
        return None;
    }
    let (mut exonic, mut intronic, mut span) = (false, false, false);
    let last = t.exons.len() - 1;
    for b in &p.blocks {
        let bs = b.start;
        let be = b.end - 1; // STAR uses inclusive block ends
        // Exon with the greatest start <= block start.
        let ex1 = match t.exons.partition_point(|e| e.start <= bs) {
            0 => return None, // contained in the span yet before the first exon: cannot happen
            k => k - 1,
        };
        if ex1 == last {
            exonic = true;
            break; // align end <= transcript end, so the rest is exonic
        }
        if be - bs < TOL {
            continue; // block too short to call
        }
        let e_end = t.exons[ex1].end - 1; // inclusive
        let n_start = t.exons[ex1 + 1].start;
        let n_end = t.exons[ex1 + 1].end - 1;
        if bs + TOL <= e_end {
            // start certainly in exon1
            if be <= e_end + TOL {
                exonic = true;
            } else {
                span = true;
            }
        } else if bs + TOL < n_start {
            // start in the intron
            if be >= n_start + TOL {
                span = true;
            } else if be > e_end + TOL {
                if n_start - 1 - e_end > 1_000_000 {
                    return None; // giant intron must not swallow small genes
                }
                intronic = true;
            } // else: too close to the boundary, no call
        } else {
            // start too close to the next exon's start
            if be > n_end + TOL {
                span = true;
            } else if be >= n_start + TOL {
                exonic = true;
            } // else: no call
        }
        if !p.junctions.is_empty() && (intronic || span) {
            return None; // spliced aligns must rest cleanly in the exon structure
        }
    }
    Some(if span {
        Vs::ExonIntronSpan
    } else if !intronic {
        Vs::Concordant
    } else if exonic {
        Vs::ExonIntron
    } else {
        Vs::Intron
    })
}

/// Velocyto per-(placement, transcript) type bits: `1 << status`, with span implying both
/// exonic and intronic (mirrors `classifyAlign` exactly).
pub fn velocyto_bits(v: Vs) -> u8 {
    match v {
        Vs::Concordant => VB_CONCORDANT,
        Vs::ExonIntron => VB_EXONINTRON,
        Vs::Intron => VB_INTRON,
        Vs::ExonIntronSpan => VB_SPAN | VB_INTRON | VB_CONCORDANT,
        Vs::Incompatible => 0,
    }
}

/// STARsolo `--soloFeatures Gene` assignment with `--soloStrand Forward`: the set of genes with at
/// least one same-strand transcript classifying the placement as Concordant. The caller applies
/// the uniqueness rules (exactly one gene, `NH == 1`), which are policy rather than geometry.
pub fn concordant_genes(p: &Placement, anno: &crate::Annotation, chrom: u32) -> Vec<u32> {
    let (mut txbuf, mut genes) = (Vec::new(), Vec::new());
    concordant_genes_into(p, anno, chrom, &mut txbuf, &mut genes);
    genes
}

/// Allocation-free variant of [`concordant_genes`]: `txbuf` is the transcript-overlap scratch and
/// `genes` receives the result (cleared first, same order as the allocating version).
pub fn concordant_genes_into(
    p: &Placement,
    anno: &crate::Annotation,
    chrom: u32,
    txbuf: &mut Vec<u32>,
    genes: &mut Vec<u32>,
) {
    concordant_genes_stranded_into(
        p, anno, chrom, SoloStrand::Forward, txbuf, genes,
    );
}

/// Strand-parameterized variant of [`concordant_genes_into`]. Geometry and gene uniqueness are
/// unchanged; only the STARsolo alignment/transcript strand relationship is selected.
pub fn concordant_genes_stranded_into(
    p: &Placement,
    anno: &crate::Annotation,
    chrom: u32,
    solo_strand: SoloStrand,
    txbuf: &mut Vec<u32>,
    genes: &mut Vec<u32>,
) {
    genes.clear();
    anno.overlapping_into(chrom, p.start(), p.end(), txbuf);
    let p_rev = matches!(p.strand, evidence_io::Strand::Reverse);
    for &txi in txbuf.iter() {
        let t = &anno.transcripts[txi as usize];
        if !solo_strand.accepts(p_rev, t.strand_rev) {
            continue;
        }
        // A gene is in the set once ANY of its transcripts classifies Concordant, so further
        // transcripts of an already-admitted gene need no classification (same set, same order).
        if genes.contains(&t.gene) {
            continue;
        }
        if align_vs_transcript(p, t) == Vs::Concordant {
            genes.push(t.gene);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Exon;
    use evidence_io::{Block, Junction, Placement, Strand};

    fn tx(exons: &[(u32, u32)]) -> Transcript {
        Transcript {
            gene: 0,
            chrom: 0,
            strand_rev: false,
            exons: exons.iter().map(|&(s, e)| Exon { start: s, end: e }).collect(),
        }
    }

    fn pl(blocks: &[(u32, u32)], junctions: &[(u32, u32)]) -> Placement {
        Placement {
            chrom: 0,
            strand: Strand::Forward,
            blocks: blocks.iter().map(|&(s, e)| Block { start: s, end: e }).collect(),
            junctions: junctions.iter().map(|&(d, a)| Junction { donor: d, acceptor: a }).collect(),
            nm: 0,
            score: 0,
            nh: 1,
            clip: (0, 0),
        }
    }

    // Transcript: exons [100,200) and [300,400); intron [200,300).
    fn t2() -> Transcript {
        tx(&[(100, 200), (300, 400)])
    }

    #[test]
    fn fully_exonic_single_block_is_concordant() {
        assert_eq!(align_vs_transcript(&pl(&[(120, 180)], &[]), &t2()), Vs::Concordant);
    }

    #[test]
    fn spliced_read_matching_the_intron_is_concordant() {
        let p = pl(&[(150, 200), (300, 350)], &[(200, 300)]);
        assert_eq!(align_vs_transcript(&p, &t2()), Vs::Concordant);
    }

    #[test]
    fn junction_off_by_one_is_incompatible() {
        // STAR demands exact donor/acceptor coincidence; a shifted junction disqualifies the
        // transcript entirely, even though every base is exonic.
        let p = pl(&[(150, 201), (301, 350)], &[(201, 301)]);
        assert_eq!(align_vs_transcript(&p, &t2()), Vs::Incompatible);
    }

    #[test]
    fn fully_intronic_read_is_intron() {
        assert_eq!(align_vs_transcript(&pl(&[(220, 280)], &[]), &t2()), Vs::Intron);
    }

    #[test]
    fn block_crossing_exon_end_is_span() {
        assert_eq!(align_vs_transcript(&pl(&[(180, 220)], &[]), &t2()), Vs::ExonIntronSpan);
    }

    #[test]
    fn read_outside_span_is_incompatible() {
        assert_eq!(align_vs_transcript(&pl(&[(50, 90)], &[]), &t2()), Vs::Incompatible);
        assert_eq!(align_vs_transcript(&pl(&[(390, 420)], &[]), &t2()), Vs::Incompatible);
    }

    #[test]
    fn single_exon_transcript_and_contained_read() {
        let t = tx(&[(1000, 2000)]);
        assert_eq!(align_vs_transcript(&pl(&[(1100, 1191)], &[]), &t), Vs::Concordant);
    }

    #[test]
    fn spliced_read_against_single_exon_transcript_is_incompatible() {
        let t = tx(&[(1000, 2000)]);
        let p = pl(&[(1100, 1150), (1500, 1541)], &[(1150, 1500)]);
        assert_eq!(align_vs_transcript(&p, &t), Vs::Incompatible);
    }

    #[test]
    fn three_exon_skip_junction_must_match_consecutive_exons() {
        // Exons A[0,100) B[200,300) C[400,500). A junction 100->400 (skipping B) is NOT an intron
        // of this transcript (introns are 100->200 and 300->400), so the transcript is
        // incompatible — matching STAR, where each junction must agree with the consecutive exon
        // walk. A skipping isoform would be its own transcript record.
        let t = tx(&[(0, 100), (200, 300), (400, 500)]);
        let p = pl(&[(50, 100), (400, 450)], &[(100, 400)]);
        assert_eq!(align_vs_transcript(&p, &t), Vs::Incompatible);
    }

    // ---- alignToTranscriptMinOverlap port (velocyto mode, TOL = 6) ----
    // Transcript: exons [100,200) and [300,400); intron [200,300).

    #[test]
    fn velo_min_overlap_exonic() {
        let t = tx(&[(100, 200), (300, 400)]);
        assert_eq!(align_vs_transcript_min_overlap(&pl(&[(110, 190)], &[]), &t), Some(Vs::Concordant));
    }

    #[test]
    fn velo_min_overlap_intronic_and_span() {
        let t = tx(&[(100, 200), (300, 400)]);
        assert_eq!(align_vs_transcript_min_overlap(&pl(&[(210, 280)], &[]), &t), Some(Vs::Intron));
        assert_eq!(align_vs_transcript_min_overlap(&pl(&[(150, 250)], &[]), &t), Some(Vs::ExonIntronSpan));
    }

    #[test]
    fn velo_min_overlap_boundary_tolerance() {
        let t = tx(&[(100, 200), (300, 400)]);
        // Ends 5 bp into the intron: within tolerance, still exonic.
        assert_eq!(align_vs_transcript_min_overlap(&pl(&[(150, 205)], &[]), &t), Some(Vs::Concordant));
        // Ends 10 bp into the intron: spans.
        assert_eq!(align_vs_transcript_min_overlap(&pl(&[(150, 210)], &[]), &t), Some(Vs::ExonIntronSpan));
        // Starts <6 bp before the intron end and stays near the boundary: no call -> Concordant default? no:
        // a short no-call block alone leaves all flags false -> !intronic -> Concordant.
    }

    #[test]
    fn velo_min_overlap_giant_intron_disqualifies() {
        let t = tx(&[(100, 200), (2_000_000, 2_000_100)]);
        assert_eq!(align_vs_transcript_min_overlap(&pl(&[(500, 580)], &[]), &t), None);
    }

    #[test]
    fn velo_min_overlap_spliced_in_intron_incompatible() {
        let t = tx(&[(100, 200), (300, 400)]);
        // Spliced align whose first block sits in the intron: not consistent with the model.
        let p = pl(&[(210, 260), (350, 380)], &[(260, 350)]);
        assert_eq!(align_vs_transcript_min_overlap(&p, &t), None);
    }

    #[test]
    fn velo_min_overlap_spliced_concordant() {
        let t = tx(&[(100, 200), (300, 400)]);
        let p = pl(&[(150, 200), (300, 340)], &[(200, 300)]);
        assert_eq!(align_vs_transcript_min_overlap(&p, &t), Some(Vs::Concordant));
    }

    #[test]
    fn velo_min_overlap_last_exon_shortcut() {
        let t = tx(&[(100, 200), (300, 400)]);
        // Block starting in the last exon: exonic by the shortcut.
        assert_eq!(align_vs_transcript_min_overlap(&pl(&[(310, 390)], &[]), &t), Some(Vs::Concordant));
    }

    #[test]
    fn velo_bits_span_implies_both() {
        assert_eq!(velocyto_bits(Vs::ExonIntronSpan), VB_SPAN | VB_INTRON | VB_CONCORDANT);
        assert_eq!(velocyto_bits(Vs::Concordant), VB_CONCORDANT);
    }

    #[test]
    fn solo_strand_relationships_match_star() {
        use SoloStrand::{Forward, Reverse, Unstranded};
        assert!(Forward.accepts(false, false));
        assert!(Forward.accepts(true, true));
        assert!(!Forward.accepts(false, true));
        assert!(Reverse.accepts(false, true));
        assert!(Reverse.accepts(true, false));
        assert!(!Reverse.accepts(false, false));
        assert!(Unstranded.accepts(false, false));
        assert!(Unstranded.accepts(false, true));
        assert!(Unstranded.accepts(true, false));
        assert!(Unstranded.accepts(true, true));
    }

}

#[cfg(test)]
mod genefull_tests {
    use super::*;
    use crate::{Annotation, Exon};
    use evidence_io::{Block, Strand};

    fn annotation(transcripts: Vec<Transcript>) -> Annotation {
        Annotation {
            gene_ids: Vec::new(),
            gene_names: Vec::new(),
            transcript_ids: Vec::new(),
            source_exons: Vec::new(),
            chrom_ids: Default::default(),
            index: crate::build_index(&transcripts),
            transcripts,
        }
    }

    fn tx(gene: u32, chrom: u32, reverse: bool, exons: &[(u32, u32)]) -> Transcript {
        Transcript {
            gene,
            chrom,
            strand_rev: reverse,
            exons: exons
                .iter()
                .map(|&(start, end)| Exon { start, end })
                .collect(),
        }
    }

    fn placement(blocks: &[(u32, u32)], reverse: bool) -> Placement {
        Placement {
            chrom: 99,
            strand: if reverse {
                Strand::Reverse
            } else {
                Strand::Forward
            },
            blocks: blocks
                .iter()
                .map(|&(start, end)| Block { start, end })
                .collect(),
            junctions: Vec::new(),
            nm: 0,
            score: 0,
            nh: 1,
            clip: (0, 0),
        }
    }

    #[test]
    fn genefull_boundaries_isoforms_nested_genes_and_skipped_introns() {
        let a = annotation(vec![
            tx(0, 0, false, &[(100, 200), (400, 500)]),
            tx(0, 0, false, &[(700, 800)]),
            tx(1, 0, false, &[(250, 300)]),
            tx(2, 0, true, &[(100, 800)]),
            tx(0, 1, false, &[(900, 1000)]),
            tx(0, 0, true, &[(900, 1000)]),
        ]);
        let index = GeneFullIndex::new(&a);
        for (blocks, strand, expected) in [
            (vec![(99, 100)], SoloStrand::Forward, vec![]),
            (vec![(100, 101)], SoloStrand::Forward, vec![0]),
            (vec![(799, 800)], SoloStrand::Forward, vec![0]),
            (vec![(800, 801)], SoloStrand::Forward, vec![]),
            (vec![(600, 650)], SoloStrand::Forward, vec![0]), // gap between isoforms
            (vec![(260, 270)], SoloStrand::Forward, vec![0, 1]),
            (vec![(150, 200), (400, 450)], SoloStrand::Forward, vec![0]),
            (vec![(150, 190), (410, 450)], SoloStrand::Forward, vec![0]), // novel splice
            (vec![(100, 110)], SoloStrand::Reverse, vec![2]),
            (vec![(100, 110)], SoloStrand::Unstranded, vec![0, 2]),
            (vec![(950, 960)], SoloStrand::Forward, vec![]), // no cross-strand union
            (vec![(200, 200)], SoloStrand::Unstranded, vec![]),
        ] {
            let mut got = vec![999];
            index.genes_into(&placement(&blocks, false), 0, strand, &mut got);
            assert_eq!(got, expected, "{blocks:?} {strand:?}");
        }
        let mut got = Vec::new();
        index.genes_into(
            &placement(&[(950, 960)], false),
            1,
            SoloStrand::Forward,
            &mut got,
        );
        assert_eq!(got, vec![0]);
        index.genes_into(
            &placement(&[(950, 960)], false),
            2,
            SoloStrand::Forward,
            &mut got,
        );
        assert!(got.is_empty());
    }

    #[test]
    fn genefull_index_matches_independent_exhaustive_overlap() {
        let mut seed = 421u64;
        let mut random = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed >> 32) as u32
        };
        let transcripts = (0..300)
            .map(|_| {
                let start = random() % 2000;
                tx(
                    random() % 30,
                    random() % 4,
                    random() % 2 == 0,
                    &[
                        (start, start + 1 + random() % 30),
                        (start + 50, start + 100),
                    ],
                )
            })
            .collect();
        let a = annotation(transcripts);
        let index = GeneFullIndex::new(&a);
        // Oracle scans every exon to derive bounds, then every gene; no interval index.
        let mut spans = std::collections::BTreeMap::<(u32, u32, bool), (u32, u32)>::new();
        for t in &a.transcripts {
            for exon in &t.exons {
                let span = spans
                    .entry((t.gene, t.chrom, t.strand_rev))
                    .or_insert((u32::MAX, 0));
                span.0 = span.0.min(exon.start);
                span.1 = span.1.max(exon.end);
            }
        }
        for _ in 0..2000 {
            let start = random() % 2500;
            let p = placement(
                &[(start, start + random() % 100), (start + 150, start + 200)],
                random() % 2 == 0,
            );
            let chrom = random() % 5;
            for strand in [
                SoloStrand::Forward,
                SoloStrand::Reverse,
                SoloStrand::Unstranded,
            ] {
                let mut expected = std::collections::BTreeSet::new();
                for (&(gene, c, reverse), &(s, e)) in &spans {
                    if c == chrom
                        && strand.accepts(matches!(p.strand, Strand::Reverse), reverse)
                        && p.blocks
                            .iter()
                            .any(|b| b.start < b.end && b.start < e && s < b.end)
                    {
                        expected.insert(gene);
                    }
                }
                let mut got = Vec::new();
                index.genes_into(&p, chrom, strand, &mut got);
                assert_eq!(got, expected.into_iter().collect::<Vec<_>>());
            }
        }
    }
}
