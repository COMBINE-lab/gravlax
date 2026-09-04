//! Paired annotation replay and consequence explanation.
//!
//! The archive is decoded once.  Every evidence row is evaluated against both annotations under
//! the same strand and multimapper policy, but the two sides retain independent class aggregation
//! and final 1MM UMI collapse.  The resulting count delta therefore remains exact even when a
//! changed class also changes which neighbouring classes collapse.  The class/state ledger is
//! complete; molecule records and changed rows are bounded deterministic witnesses. Causes name
//! observed state transitions (including an explicit annotation-order tie-break artifact) and
//! deliberately do not claim unique minimal counterfactual attribution.

use super::{decode_chunk, StreamingReplayArchive};
use crate::rows::{Extracted, MolRec, PatAlt, Row, SAME_SHAPE};
use anno::assign::SoloStrand;
use anyhow::{bail, Context, Result};
use evidence_io::{Block, Junction, Placement, Strand};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
#[cfg(test)]
use std::fs::{self, OpenOptions};
#[cfg(test)]
use std::io::{BufWriter, Write};
use std::path::Path;
#[cfg(test)]
use std::sync::atomic::{AtomicU64, Ordering};

pub const REPORT_SCHEMA: &str = "gravlax.annotation-comparison.v1";

/// How stable annotation identifiers are joined across the two sides.  No gene-symbol join is
/// performed.  `Unversioned` removes exactly one terminal `.<digits>` suffix and rejects any
/// within-annotation collision created by that normalization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GeneKeyPolicy {
    Unversioned,
    Exact,
}

/// Replay policy shared by the two annotations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompareOptions {
    pub solo_strand: SoloStrand,
    pub gene_key_policy: GeneKeyPolicy,
    /// Maximum number of molecule witnesses in the whole report.  The bounded scan retains the
    /// earliest directly changed record for each of the first distinct changed classes observed
    /// in canonical archive order.
    pub max_molecule_witnesses: usize,
    /// Maximum changed evidence rows retained inside any one molecule witness.
    pub max_row_transitions_per_molecule: usize,
}

impl Default for CompareOptions {
    fn default() -> Self {
        Self {
            solo_strand: SoloStrand::Forward,
            gene_key_policy: GeneKeyPolicy::Unversioned,
            max_molecule_witnesses: 10_000,
            max_row_transitions_per_molecule: 32,
        }
    }
}

/// Machine-readable statement of what is exact and where interpretation stops.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExactnessBoundary {
    /// Counts and signed deltas come from two complete, independent final collapses.
    pub final_count_deltas_are_exact: bool,
    /// Every changed class and its two independently collapsed states is present in the ledger.
    pub class_state_ledger_is_complete: bool,
    /// Molecule records and their rows are explicitly bounded witnesses rather than a full ledger.
    pub molecule_witnesses_are_bounded: bool,
    /// Causes name observed algorithmic state transitions, not a unique minimal counterfactual
    /// cause. The annotation-order tie-break cause is explicitly non-biological.
    pub causes_are_observed_state_transitions_not_counterfactual: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReplaySummary {
    pub assigned_rows: u64,
    pub selected_classes: u64,
    pub final_gene_umis: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ArchiveIdentity {
    pub archive_version: u32,
    /// Root of the authenticated v2 directory, observed from the same open reader used to decode
    /// evidence.  Legacy v1 archives intentionally report `None`.
    pub rooted_content_commitment_hex: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CountDelta {
    pub cell: u32,
    pub comparison_gene_id: String,
    pub gene_id_before: Option<String>,
    pub gene_id_after: Option<String>,
    pub before: u32,
    pub after: u32,
    pub delta: i64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GeneSupport {
    pub comparison_gene_id: String,
    pub gene_id: String,
    pub weight: u32,
}

/// The independently resolved state of one UMI class on one annotation side.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ClassState {
    pub gene_support: Vec<GeneSupport>,
    pub selected_comparison_gene_id: Option<String>,
    pub selected_gene_id: Option<String>,
    pub selected_weight: u32,
    /// True exactly when this class contributes one final gene/UMI count.
    pub counted: bool,
    /// The class that survives greedy 1MM collapse, or `None` when unassigned.
    pub canonical_class: Option<u32>,
    /// All adjacent classes assigned to the same gene on this side, sorted by class id.
    pub same_gene_neighbors: Vec<u32>,
}

impl ClassState {
    fn final_contribution(&self) -> Option<&str> {
        self.counted
            .then_some(self.selected_comparison_gene_id.as_deref())
            .flatten()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
#[expect(
    clippy::enum_variant_names,
    reason = "the full variant names are part of the serialized annotation-comparison schema"
)]
pub enum TransitionCause {
    CandidateSetChanged,
    RowAssignmentChanged,
    ClassSupportChanged,
    /// Equal comparison-key support tied for the maximum on both sides, but annotation-local gene
    /// order selected different winners. This is an explicit replay-method artifact, not a
    /// biological structural change.
    AnnotationOrderTieBreakChanged,
    ClassWinnerChanged,
    CollapseNeighborhoodChanged,
    CollapseOutcomeChanged,
    FinalContributionChanged,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ClassTransitionKind {
    GainedFinalCount,
    LostFinalCount,
    ReassignedFinalCount,
    ChangedWithoutFinalCountDelta,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceKind {
    UniqueRepresentative,
    MultimapperSignature,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GeneIdentity {
    /// Cross-annotation key selected by [`GeneKeyPolicy`].
    pub comparison_gene_id: String,
    /// Identifier exactly as stored in this annotation.
    pub gene_id: String,
}

/// A row whose annotation-dependent candidate set or singleton assignment changed.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RowTransition {
    pub row_index: u32,
    pub kind: EvidenceKind,
    pub chrom: String,
    pub pos: u32,
    pub strand_reverse: bool,
    pub shape: u32,
    pub pattern: Option<u32>,
    pub weight: u32,
    pub before_candidates: Vec<GeneIdentity>,
    pub after_candidates: Vec<GeneIdentity>,
    pub before_singleton: Option<GeneIdentity>,
    pub after_singleton: Option<GeneIdentity>,
    pub causes: Vec<TransitionCause>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ClassEvidenceSummary {
    pub molecule_records: u64,
    pub rows: u64,
    pub changed_rows: u64,
    pub molecule_witnesses: u64,
    pub omitted_molecule_witnesses: u64,
    pub changed_row_witnesses: u64,
    pub omitted_changed_row_witnesses: u64,
}

/// One class-level transition after both sides have completed their own global collapse.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ClassTransition {
    pub cell: u32,
    pub umi_class: u32,
    pub kind: ClassTransitionKind,
    pub evidence: ClassEvidenceSummary,
    pub before: ClassState,
    pub after: ClassState,
    pub causes: Vec<TransitionCause>,
}

/// One archive molecule record belonging to a changed class.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MoleculeWitness {
    /// Ordinal in canonical archive chunk/molecule order.
    pub ordinal: u64,
    pub cell: u32,
    pub umi_class: u32,
    pub chrom: String,
    pub anchor: u32,
    pub rows: u32,
    /// Only rows whose candidate set or singleton assignment changed are repeated here.
    pub changed_rows: Vec<RowTransition>,
    pub changed_rows_total: u64,
    pub changed_rows_omitted: u64,
    pub before_class: ClassState,
    pub after_class: ClassState,
    pub causes: Vec<TransitionCause>,
}

/// Deterministic paired-replay report.  Every vector is sorted by its documented stable key.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnnotationComparison {
    pub schema: String,
    pub archive_identity: ArchiveIdentity,
    /// Human-readable 16-base cell barcodes in archive dictionary order. Numeric `cell` fields
    /// throughout the report are indexes into this vector.
    pub cell_barcodes: Vec<String>,
    pub gene_key_policy: GeneKeyPolicy,
    pub exactness: ExactnessBoundary,
    pub archive_passes: u32,
    pub decoded_chunks: u64,
    pub decoded_molecules: u64,
    pub evidence_rows: u64,
    /// Exact number of molecule records belonging to a changed class.
    pub changed_molecule_records: u64,
    pub unchanged_molecule_records: u64,
    pub molecule_witnesses_omitted: u64,
    pub before: ReplaySummary,
    pub after: ReplaySummary,
    /// Nonzero deltas sorted by `(cell, comparison_gene_id)`.
    pub count_deltas: Vec<CountDelta>,
    /// Changed classes sorted by `(cell, umi_class)`.
    pub class_transitions: Vec<ClassTransition>,
    /// Bounded, deterministic witnesses sorted by archive ordinal.
    pub molecule_witnesses: Vec<MoleculeWitness>,
}

struct GeneKeyIndex {
    by_gene: Vec<GeneIdentity>,
    original_by_key: BTreeMap<String, String>,
}

impl GeneKeyIndex {
    fn new(annotation: &anno::Annotation, policy: GeneKeyPolicy, side: &str) -> Result<Self> {
        let mut by_gene = Vec::with_capacity(annotation.gene_ids.len());
        let mut original_by_key = BTreeMap::new();
        for gene_id in &annotation.gene_ids {
            let key = comparison_gene_key(gene_id, policy).to_string();
            if key.is_empty() {
                bail!(
                    "{side} annotation gene id {gene_id:?} has an empty comparison key under \
                     {policy:?} matching"
                );
            }
            if let Some(previous) = original_by_key.insert(key.clone(), gene_id.clone()) {
                bail!(
                    "{side} annotation gene-id normalization collision: {previous:?} and \
                     {gene_id:?} both map to {key:?}"
                );
            }
            by_gene.push(GeneIdentity {
                comparison_gene_id: key,
                gene_id: gene_id.clone(),
            });
        }
        Ok(Self {
            by_gene,
            original_by_key,
        })
    }

    fn identity(&self, gene: u32) -> Result<&GeneIdentity> {
        self.by_gene
            .get(gene as usize)
            .with_context(|| format!("annotation gene index {gene} has no stable id"))
    }
}

struct ComparisonGeneKeys {
    before: GeneKeyIndex,
    after: GeneKeyIndex,
}

impl ComparisonGeneKeys {
    fn new(
        before: &anno::Annotation,
        after: &anno::Annotation,
        policy: GeneKeyPolicy,
    ) -> Result<Self> {
        Ok(Self {
            before: GeneKeyIndex::new(before, policy, "before")?,
            after: GeneKeyIndex::new(after, policy, "after")?,
        })
    }
}

fn comparison_gene_key(gene_id: &str, policy: GeneKeyPolicy) -> &str {
    if policy == GeneKeyPolicy::Exact {
        return gene_id;
    }
    let Some((prefix, suffix)) = gene_id.rsplit_once('.') else {
        return gene_id;
    };
    if !suffix.is_empty() && suffix.bytes().all(|byte| byte.is_ascii_digit()) {
        prefix
    } else {
        gene_id
    }
}

/// Compare two source- or compiled-annotation models while decoding the `.aie` archive once.
pub fn compare_archive(
    archive: &Path,
    before: &anno::Annotation,
    after: &anno::Annotation,
    options: CompareOptions,
) -> Result<AnnotationComparison> {
    // Validate the join before opening/decoding the archive: ambiguous identifiers are a policy
    // error, not something that should yield a partial scientific comparison.
    let gene_keys = ComparisonGeneKeys::new(before, after, options.gene_key_policy)?;
    let stream = StreamingReplayArchive::open(archive)?;
    compare_stream(&stream, before, after, &gene_keys, options)
}

/// Compute the complete report first, then publish JSON atomically without replacing a path that
/// already exists.  A replay, validation, serialization, or publication error leaves no partial
/// destination file.
#[cfg(test)]
pub fn write_comparison_json_noclobber(
    archive: &Path,
    before: &anno::Annotation,
    after: &anno::Annotation,
    options: CompareOptions,
    output: &Path,
) -> Result<AnnotationComparison> {
    let report = compare_archive(archive, before, after, options)?;
    publish_json_noclobber(output, &report)?;
    Ok(report)
}

#[derive(Debug, Clone)]
struct AssignmentTuple {
    cell: u32,
    class: u32,
    gene: u32,
    weight: u32,
}

#[derive(Debug, Default, Clone)]
struct DirectClassEvidence {
    molecule_records: u64,
    rows: u64,
    changed_rows: u64,
    causes: BTreeSet<TransitionCause>,
}

#[derive(Debug)]
struct MoleculeDraft {
    ordinal: u64,
    cell: u32,
    class: u32,
    chrom: String,
    anchor: u32,
    rows: u32,
    changed_rows: Vec<RowTransition>,
    changed_rows_total: u64,
    direct_causes: BTreeSet<TransitionCause>,
}

#[derive(Default)]
struct BatchPart {
    before: Vec<AssignmentTuple>,
    after: Vec<AssignmentTuple>,
    witness_candidates: BTreeMap<u32, MoleculeDraft>,
    class_evidence: BTreeMap<u32, DirectClassEvidence>,
}

struct PairClassifier<'a> {
    x: &'a Extracted,
    before: &'a anno::Annotation,
    after: &'a anno::Annotation,
    before_chrom: Vec<Option<u32>>,
    after_chrom: Vec<Option<u32>>,
    strand: SoloStrand,
    place: Placement,
    before_tx: Vec<u32>,
    after_tx: Vec<u32>,
    before_genes: Vec<u32>,
    after_genes: Vec<u32>,
    before_alt: Vec<u32>,
    after_alt: Vec<u32>,
}

impl<'a> PairClassifier<'a> {
    fn new(
        x: &'a Extracted,
        before: &'a anno::Annotation,
        after: &'a anno::Annotation,
        strand: SoloStrand,
    ) -> Self {
        let before_chrom = x
            .chrom_names
            .iter()
            .map(|name| before.chrom_ids.get(name).copied())
            .collect();
        let after_chrom = x
            .chrom_names
            .iter()
            .map(|name| after.chrom_ids.get(name).copied())
            .collect();
        Self {
            x,
            before,
            after,
            before_chrom,
            after_chrom,
            strand,
            place: Placement {
                chrom: 0,
                strand: Strand::Forward,
                blocks: Vec::new(),
                junctions: Vec::new(),
                nm: 0,
                score: 0,
                nh: 1,
                clip: (0, 0),
            },
            before_tx: Vec::new(),
            after_tx: Vec::new(),
            before_genes: Vec::new(),
            after_genes: Vec::new(),
            before_alt: Vec::new(),
            after_alt: Vec::new(),
        }
    }

    fn classify(&mut self, row: &Row) -> Result<()> {
        self.before_genes.clear();
        self.after_genes.clear();
        if row.pattern == u32::MAX {
            let bc = self.before_chrom.get(row.chrom as usize).copied().flatten();
            let ac = self.after_chrom.get(row.chrom as usize).copied().flatten();
            if bc.is_none() && ac.is_none() {
                return Ok(());
            }
            let shape = self.x.shapes.get(row.shape as usize).with_context(|| {
                format!("row shape {} is outside the shape dictionary", row.shape)
            })?;
            placement_into_checked(
                &mut self.place,
                row.chrom,
                row.pos,
                row.strand_rev,
                shape,
                1,
            )?;
            if let Some(chrom) = bc {
                anno::assign::concordant_genes_stranded_into(
                    &self.place,
                    self.before,
                    chrom,
                    self.strand,
                    &mut self.before_tx,
                    &mut self.before_genes,
                );
            }
            if let Some(chrom) = ac {
                anno::assign::concordant_genes_stranded_into(
                    &self.place,
                    self.after,
                    chrom,
                    self.strand,
                    &mut self.after_tx,
                    &mut self.after_genes,
                );
            }
        } else {
            let pattern = self.x.patterns.get(row.pattern as usize).with_context(|| {
                format!(
                    "row pattern {} is outside the pattern dictionary",
                    row.pattern
                )
            })?;
            for alt in pattern {
                self.classify_alternative(row, alt)?;
            }
        }
        Ok(())
    }

    fn classify_alternative(&mut self, row: &Row, alt: &PatAlt) -> Result<()> {
        let bc = self.before_chrom.get(alt.chrom as usize).copied().flatten();
        let ac = self.after_chrom.get(alt.chrom as usize).copied().flatten();
        // Gene replay's SkipAlt policy does not inspect geometry on a chromosome absent from an
        // annotation.  Preserve that policy independently on each side.
        if bc.is_none() && ac.is_none() {
            return Ok(());
        }
        let apos = i64::from(row.pos)
            .checked_add(alt.offset)
            .and_then(|p| u32::try_from(p).ok())
            .with_context(|| {
                format!(
                    "multimapper alternative coordinate {} + {} is outside u32",
                    row.pos, alt.offset
                )
            })?;
        let shape_id = if alt.shape == SAME_SHAPE {
            row.shape
        } else {
            alt.shape
        };
        let shape = self.x.shapes.get(shape_id as usize).with_context(|| {
            format!("alternative shape {shape_id} is outside the shape dictionary")
        })?;
        placement_into_checked(
            &mut self.place,
            alt.chrom,
            apos,
            row.strand_rev != alt.strand_flip,
            shape,
            2,
        )?;
        if let Some(chrom) = bc {
            anno::assign::concordant_genes_stranded_into(
                &self.place,
                self.before,
                chrom,
                self.strand,
                &mut self.before_tx,
                &mut self.before_alt,
            );
            extend_unique(&mut self.before_genes, &self.before_alt);
        }
        if let Some(chrom) = ac {
            anno::assign::concordant_genes_stranded_into(
                &self.place,
                self.after,
                chrom,
                self.strand,
                &mut self.after_tx,
                &mut self.after_alt,
            );
            extend_unique(&mut self.after_genes, &self.after_alt);
        }
        Ok(())
    }
}

fn extend_unique(dst: &mut Vec<u32>, src: &[u32]) {
    for &value in src {
        if !dst.contains(&value) {
            dst.push(value);
        }
    }
}

fn placement_into_checked(
    place: &mut Placement,
    chrom: u32,
    pos: u32,
    strand_reverse: bool,
    shape: &evidence_io::archive::Shape,
    nh: u16,
) -> Result<()> {
    if shape.blocks.is_empty() {
        bail!("archive shape has no blocks");
    }
    place.chrom = chrom;
    place.strand = if strand_reverse {
        Strand::Reverse
    } else {
        Strand::Forward
    };
    place.nh = nh;
    place.blocks.clear();
    place.junctions.clear();
    for &(offset, length) in &shape.blocks {
        if length == 0 {
            bail!("archive shape contains a zero-length block");
        }
        let start = pos
            .checked_add(offset)
            .context("placement block start overflow")?;
        let end = start
            .checked_add(length)
            .context("placement block end overflow")?;
        place.blocks.push(Block { start, end });
    }
    for blocks in place.blocks.windows(2) {
        place.junctions.push(Junction {
            donor: blocks[0].end,
            acceptor: blocks[1].start,
        });
    }
    Ok(())
}

fn compare_stream(
    stream: &StreamingReplayArchive,
    before: &anno::Annotation,
    after: &anno::Annotation,
    gene_keys: &ComparisonGeneKeys,
    options: CompareOptions,
) -> Result<AnnotationComparison> {
    let archive_identity = ArchiveIdentity {
        archive_version: stream.reader.archive_version(),
        rooted_content_commitment_hex: stream.reader.content_commitment().map(|root| root.to_hex()),
    };
    let mut before_tuples = Vec::with_capacity(stream.n_mols);
    let mut after_tuples = Vec::with_capacity(stream.n_mols);
    let mut witness_candidates: BTreeMap<u32, MoleculeDraft> = BTreeMap::new();
    let mut direct = vec![DirectClassEvidence::default(); stream.x.n_classes as usize];
    let mut decoded_molecules = 0usize;

    let mut chunk_ordinals = Vec::with_capacity(stream.chunks.len());
    let mut ordinal = 0u64;
    for chunk in &stream.chunks {
        chunk_ordinals.push(ordinal);
        ordinal = ordinal
            .checked_add(u64::from(chunk.n_mols))
            .context("archive molecule ordinal overflow")?;
    }

    let batch_size = rayon::current_num_threads().max(1) * 2;
    for (batch_no, batch) in stream.chunks.chunks(batch_size).enumerate() {
        let first = batch_no * batch_size;
        let decoded: Vec<Vec<MolRec>> = batch
            .par_iter()
            .enumerate()
            .map(|(j, info)| {
                let i = first + j;
                let (compressed, raw_len) = stream.reader.read_compressed_at(&format!("c{i}"))?;
                let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                decode_chunk(&raw, info, Some(&stream.cell_of_class), &stream.rans_tables)
            })
            .collect::<Result<_>>()?;
        decoded_molecules = decoded_molecules
            .checked_add(decoded.iter().map(Vec::len).sum::<usize>())
            .context("decoded molecule count overflow")?;

        let parts: Vec<BatchPart> = decoded
            .par_iter()
            .enumerate()
            .map(|(j, mols)| {
                classify_chunk(
                    &stream.x,
                    before,
                    after,
                    gene_keys,
                    mols,
                    chunk_ordinals[first + j],
                    options,
                )
            })
            .collect::<Result<_>>()?;
        for mut part in parts {
            before_tuples.append(&mut part.before);
            after_tuples.append(&mut part.after);
            // A BTreeMap gives deterministic class lookup, but class order is not archive order.
            // Merge by ordinal so a partially filled global budget retains the genuinely earliest
            // distinct changed classes rather than whichever class IDs happen to sort first.
            let mut part_candidates: Vec<_> = part.witness_candidates.into_values().collect();
            part_candidates.sort_by_key(|candidate| (candidate.ordinal, candidate.class));
            for candidate in part_candidates {
                let class = candidate.class;
                if let Some(old) = witness_candidates.get_mut(&class) {
                    if witness_is_better(&candidate, old) {
                        *old = candidate;
                    }
                } else if witness_candidates.len() < options.max_molecule_witnesses {
                    witness_candidates.insert(class, candidate);
                }
            }
            for (class, evidence) in part.class_evidence {
                let dst = direct.get_mut(class as usize).with_context(|| {
                    format!("molecule references class {class} beyond class table")
                })?;
                dst.molecule_records = dst
                    .molecule_records
                    .checked_add(evidence.molecule_records)
                    .context("class molecule-record count overflow")?;
                dst.rows = dst
                    .rows
                    .checked_add(evidence.rows)
                    .context("class row count overflow")?;
                dst.changed_rows = dst
                    .changed_rows
                    .checked_add(evidence.changed_rows)
                    .context("class changed-row count overflow")?;
                dst.causes.extend(evidence.causes);
            }
        }
    }
    if decoded_molecules != stream.n_mols {
        bail!(
            "molecule count mismatch: {decoded_molecules} decoded vs {} in meta",
            stream.n_mols
        );
    }
    build_report(
        &stream.x,
        &stream.cell_of_class,
        gene_keys,
        before_tuples,
        after_tuples,
        witness_candidates,
        direct,
        stream.chunks.len(),
        decoded_molecules,
        archive_identity,
        options,
    )
}

fn classify_chunk(
    x: &Extracted,
    before: &anno::Annotation,
    after: &anno::Annotation,
    gene_keys: &ComparisonGeneKeys,
    mols: &[MolRec],
    ordinal_base: u64,
    options: CompareOptions,
) -> Result<BatchPart> {
    let mut classifier = PairClassifier::new(x, before, after, options.solo_strand);
    let mut part = BatchPart::default();
    for (molecule_offset, mol) in mols.iter().enumerate() {
        if mol.umi_class >= x.n_classes {
            bail!(
                "molecule references class {} beyond {} classes",
                mol.umi_class,
                x.n_classes
            );
        }
        let class_evidence = part.class_evidence.entry(mol.umi_class).or_default();
        class_evidence.molecule_records = class_evidence
            .molecule_records
            .checked_add(1)
            .context("class molecule-record count overflow")?;
        let mut changed_rows = Vec::new();
        let mut changed_rows_total = 0u64;
        let mut direct_causes = BTreeSet::new();
        let mut row_index = 0u32;
        let mut visit = |row: Row, kind: EvidenceKind| -> Result<()> {
            class_evidence.rows = class_evidence
                .rows
                .checked_add(1)
                .context("class row count overflow")?;
            classifier.classify(&row)?;
            let before_genes = classifier.before_genes.as_slice();
            let after_genes = classifier.after_genes.as_slice();
            if let [gene] = before_genes {
                part.before.push(AssignmentTuple {
                    cell: row.cell,
                    class: row.umi_class,
                    gene: *gene,
                    weight: row.weight,
                });
            }
            if let [gene] = after_genes {
                part.after.push(AssignmentTuple {
                    cell: row.cell,
                    class: row.umi_class,
                    gene: *gene,
                    weight: row.weight,
                });
            }

            let mut causes = BTreeSet::new();
            if !candidate_sets_equal(
                &gene_keys.before,
                before_genes,
                &gene_keys.after,
                after_genes,
            )? {
                causes.insert(TransitionCause::CandidateSetChanged);
            }
            if singleton_key(&gene_keys.before, before_genes)?
                != singleton_key(&gene_keys.after, after_genes)?
            {
                causes.insert(TransitionCause::RowAssignmentChanged);
            }
            if !causes.is_empty() {
                let before_candidates = stable_candidates(&gene_keys.before, before_genes)?;
                let after_candidates = stable_candidates(&gene_keys.after, after_genes)?;
                let before_singleton = singleton_id(&before_candidates);
                let after_singleton = singleton_id(&after_candidates);
                changed_rows_total = changed_rows_total
                    .checked_add(1)
                    .context("molecule changed-row count overflow")?;
                class_evidence.changed_rows = class_evidence
                    .changed_rows
                    .checked_add(1)
                    .context("class changed-row count overflow")?;
                class_evidence.causes.extend(causes.iter().copied());
                direct_causes.extend(causes.iter().copied());
                let chrom = x
                    .chrom_names
                    .get(row.chrom as usize)
                    .cloned()
                    .unwrap_or_else(|| format!("#{}", row.chrom));
                if changed_rows.len() < options.max_row_transitions_per_molecule {
                    changed_rows.push(RowTransition {
                        row_index,
                        kind,
                        chrom,
                        pos: row.pos,
                        strand_reverse: row.strand_rev,
                        shape: row.shape,
                        pattern: (row.pattern != u32::MAX).then_some(row.pattern),
                        weight: row.weight,
                        before_candidates,
                        after_candidates,
                        before_singleton,
                        after_singleton,
                        causes: causes.into_iter().collect(),
                    });
                }
            }
            row_index = row_index
                .checked_add(1)
                .context("molecule row count exceeds u32")?;
            Ok(())
        };
        for chain in &mol.chains {
            for &(pos, shape) in &chain.reps {
                visit(
                    Row {
                        cell: mol.cell,
                        umi_class: mol.umi_class,
                        chrom: mol.chrom,
                        pos,
                        strand_rev: mol.strand_rev,
                        shape,
                        weight: chain.weight,
                        pattern: u32::MAX,
                    },
                    EvidenceKind::UniqueRepresentative,
                )?;
            }
        }
        for &(pos, shape, pattern, weight) in &mol.mms {
            visit(
                Row {
                    cell: mol.cell,
                    umi_class: mol.umi_class,
                    chrom: mol.chrom,
                    pos,
                    strand_rev: mol.strand_rev,
                    shape,
                    weight,
                    pattern,
                },
                EvidenceKind::MultimapperSignature,
            )?;
        }
        let chrom = x
            .chrom_names
            .get(mol.chrom as usize)
            .cloned()
            .unwrap_or_else(|| format!("#{}", mol.chrom));
        if options.max_molecule_witnesses > 0 && changed_rows_total > 0 {
            let candidate = MoleculeDraft {
                ordinal: ordinal_base
                    .checked_add(molecule_offset as u64)
                    .context("archive molecule ordinal overflow")?,
                cell: mol.cell,
                class: mol.umi_class,
                chrom,
                anchor: mol.anchor(),
                rows: row_index,
                changed_rows,
                changed_rows_total,
                direct_causes,
            };
            if let Some(old) = part.witness_candidates.get_mut(&mol.umi_class) {
                if witness_is_better(&candidate, old) {
                    *old = candidate;
                }
            } else if part.witness_candidates.len() < options.max_molecule_witnesses {
                part.witness_candidates.insert(mol.umi_class, candidate);
            }
        }
    }
    Ok(part)
}

fn witness_is_better(candidate: &MoleculeDraft, current: &MoleculeDraft) -> bool {
    let candidate_direct = candidate.changed_rows_total > 0;
    let current_direct = current.changed_rows_total > 0;
    (candidate_direct && !current_direct)
        || (candidate_direct == current_direct && candidate.ordinal < current.ordinal)
}

fn stable_candidates(index: &GeneKeyIndex, genes: &[u32]) -> Result<Vec<GeneIdentity>> {
    let mut ids = Vec::with_capacity(genes.len());
    for &gene in genes {
        ids.push(index.identity(gene)?.clone());
    }
    ids.sort_by(|left, right| left.comparison_gene_id.cmp(&right.comparison_gene_id));
    ids.dedup_by(|left, right| left.comparison_gene_id == right.comparison_gene_id);
    Ok(ids)
}

fn singleton_id(candidates: &[GeneIdentity]) -> Option<GeneIdentity> {
    (candidates.len() == 1).then(|| candidates[0].clone())
}

fn candidate_sets_equal(
    left_index: &GeneKeyIndex,
    left: &[u32],
    right_index: &GeneKeyIndex,
    right: &[u32],
) -> Result<bool> {
    if left.len() != right.len() {
        return Ok(false);
    }
    for &left_gene in left {
        let key = left_index.identity(left_gene)?.comparison_gene_id.as_str();
        let mut found = false;
        for &right_gene in right {
            if right_index.identity(right_gene)?.comparison_gene_id == key {
                found = true;
                break;
            }
        }
        if !found {
            return Ok(false);
        }
    }
    Ok(true)
}

fn singleton_key<'a>(index: &'a GeneKeyIndex, genes: &[u32]) -> Result<Option<&'a str>> {
    match genes {
        [gene] => Ok(Some(index.identity(*gene)?.comparison_gene_id.as_str())),
        _ => Ok(None),
    }
}

#[cfg(test)]
fn candidate_keys(candidates: &[GeneIdentity]) -> Vec<&str> {
    candidates
        .iter()
        .map(|identity| identity.comparison_gene_id.as_str())
        .collect()
}

fn support_keys(support: &[GeneSupport]) -> Vec<(&str, u32)> {
    support
        .iter()
        .map(|entry| (entry.comparison_gene_id.as_str(), entry.weight))
        .collect()
}

fn has_tied_selected_maximum(state: &ClassState) -> bool {
    let Some(maximum) = state.gene_support.iter().map(|entry| entry.weight).max() else {
        return false;
    };
    state.selected_weight == maximum
        && state
            .gene_support
            .iter()
            .filter(|entry| entry.weight == maximum)
            .count()
            > 1
}

fn annotation_order_tie_break_changed(before: &ClassState, after: &ClassState) -> bool {
    before.selected_comparison_gene_id != after.selected_comparison_gene_id
        && support_keys(&before.gene_support) == support_keys(&after.gene_support)
        && has_tied_selected_maximum(before)
        && has_tied_selected_maximum(after)
}

#[derive(Debug, Clone, Default)]
struct InternalOutcome {
    support: Vec<(u32, u32)>,
    selected_gene: Option<u32>,
    selected_weight: u32,
    canonical_class: Option<u32>,
    same_gene_neighbors: Vec<u32>,
}

struct SideCollapse {
    counts: BTreeMap<(u32, String), u32>,
    outcomes: Vec<InternalOutcome>,
    assigned_rows: u64,
    selected_classes: u64,
    final_gene_umis: u64,
}

fn collapse_side(
    x: &Extracted,
    class_cells: &[u32],
    gene_keys: &GeneKeyIndex,
    mut tuples: Vec<AssignmentTuple>,
) -> Result<SideCollapse> {
    if class_cells.len() != x.n_classes as usize {
        bail!(
            "cell-of-class table has {} values for {} classes",
            class_cells.len(),
            x.n_classes
        );
    }
    let assigned_rows = u64::try_from(tuples.len()).context("assigned-row count exceeds u64")?;
    tuples.sort_unstable_by_key(|t| (t.cell, t.class, t.gene, t.weight));
    let mut outcomes = vec![InternalOutcome::default(); x.n_classes as usize];
    let mut kept: Vec<(u32, u32, u32, u32)> = Vec::new(); // cell, gene, class, weight
    let mut i = 0usize;
    while i < tuples.len() {
        let cell = tuples[i].cell;
        let class = tuples[i].class;
        let expected_cell = *class_cells
            .get(class as usize)
            .with_context(|| format!("assignment references class {class} beyond class table"))?;
        if cell != expected_cell {
            bail!("class {class} belongs to cell {expected_cell}, not assignment cell {cell}");
        }
        let mut support = Vec::new();
        while i < tuples.len() && tuples[i].cell == cell && tuples[i].class == class {
            let gene = tuples[i].gene;
            gene_keys.identity(gene)?;
            let mut weight = 0u32;
            while i < tuples.len()
                && tuples[i].cell == cell
                && tuples[i].class == class
                && tuples[i].gene == gene
            {
                weight = weight
                    .checked_add(tuples[i].weight)
                    .context("per-class gene evidence weight exceeds u32")?;
                i += 1;
            }
            support.push((gene, weight));
        }
        let &(best_gene, best_weight) = support
            .iter()
            .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.cmp(&a.0)))
            .context("assigned class has no gene support")?;
        outcomes[class as usize].support = support;
        outcomes[class as usize].selected_gene = Some(best_gene);
        outcomes[class as usize].selected_weight = best_weight;
        kept.push((cell, best_gene, class, best_weight));
    }

    let mut adjacency = vec![Vec::<u32>::new(); x.n_classes as usize];
    for &(a, b) in &x.edges {
        if a >= x.n_classes || b >= x.n_classes {
            bail!(
                "UMI edge ({a}, {b}) references a class beyond {}",
                x.n_classes
            );
        }
        if class_cells[a as usize] != class_cells[b as usize] {
            bail!("UMI edge ({a}, {b}) crosses cells");
        }
        adjacency[a as usize].push(b);
        adjacency[b as usize].push(a);
    }
    for neighbors in &mut adjacency {
        neighbors.sort_unstable();
        neighbors.dedup();
    }

    kept.sort_unstable();
    let mut local_counts: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    let mut k0 = 0usize;
    while k0 < kept.len() {
        let cell = kept[k0].0;
        let gene = kept[k0].1;
        let mut k1 = k0 + 1;
        while k1 < kept.len() && kept[k1].0 == cell && kept[k1].1 == gene {
            k1 += 1;
        }
        let mut order: Vec<(u32, u32)> = kept[k0..k1]
            .iter()
            .map(|entry| (entry.2, entry.3))
            .collect();
        order.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        let mut by_class: Vec<(u32, usize)> = order
            .iter()
            .enumerate()
            .map(|(rank, (class, _))| (*class, rank))
            .collect();
        by_class.sort_unstable();
        let mut canonical_rank = vec![0usize; order.len()];
        let mut roots = 0u32;
        for (rank, &(class, _)) in order.iter().enumerate() {
            let mut eligible = Vec::new();
            let mut target_rank: Option<usize> = None;
            for &neighbor in &adjacency[class as usize] {
                if let Ok(position) = by_class.binary_search_by_key(&neighbor, |entry| entry.0) {
                    let neighbor_rank = by_class[position].1;
                    eligible.push(neighbor);
                    if neighbor_rank < rank {
                        target_rank =
                            Some(target_rank.map_or(neighbor_rank, |old| old.min(neighbor_rank)));
                    }
                }
            }
            eligible.sort_unstable();
            eligible.dedup();
            outcomes[class as usize].same_gene_neighbors = eligible;
            if let Some(target) = target_rank {
                canonical_rank[rank] = canonical_rank[target];
            } else {
                canonical_rank[rank] = rank;
                roots = roots
                    .checked_add(1)
                    .context("collapsed root count exceeds u32")?;
            }
            outcomes[class as usize].canonical_class = Some(order[canonical_rank[rank]].0);
        }
        local_counts.insert((cell, gene), roots);
        k0 = k1;
    }

    let mut counts: BTreeMap<(u32, String), u32> = BTreeMap::new();
    for ((cell, gene), count) in local_counts {
        let comparison_gene_id = gene_keys.identity(gene)?.comparison_gene_id.clone();
        let entry = counts.entry((cell, comparison_gene_id)).or_default();
        *entry = entry
            .checked_add(count)
            .context("stable-gene count exceeds u32")?;
    }
    let selected_classes = outcomes
        .iter()
        .filter(|outcome| outcome.selected_gene.is_some())
        .count() as u64;
    let final_gene_umis = counts.values().try_fold(0u64, |sum, &count| {
        sum.checked_add(u64::from(count))
            .context("final gene/UMI total exceeds u64")
    })?;
    Ok(SideCollapse {
        counts,
        outcomes,
        assigned_rows,
        selected_classes,
        final_gene_umis,
    })
}

fn public_state(outcome: &InternalOutcome, gene_keys: &GeneKeyIndex) -> Result<ClassState> {
    let mut support_by_id: BTreeMap<String, (String, u32)> = BTreeMap::new();
    for &(gene, weight) in &outcome.support {
        let identity = gene_keys.identity(gene)?;
        let entry = support_by_id
            .entry(identity.comparison_gene_id.clone())
            .or_insert_with(|| (identity.gene_id.clone(), 0));
        entry.1 = entry
            .1
            .checked_add(weight)
            .context("stable-gene support exceeds u32")?;
    }
    let gene_support = support_by_id
        .into_iter()
        .map(|(comparison_gene_id, (gene_id, weight))| GeneSupport {
            comparison_gene_id,
            gene_id,
            weight,
        })
        .collect();
    let selected_identity = outcome
        .selected_gene
        .map(|gene| gene_keys.identity(gene).cloned())
        .transpose()?;
    Ok(ClassState {
        gene_support,
        selected_comparison_gene_id: selected_identity
            .as_ref()
            .map(|identity| identity.comparison_gene_id.clone()),
        selected_gene_id: selected_identity.map(|identity| identity.gene_id),
        selected_weight: outcome.selected_weight,
        counted: false,
        canonical_class: outcome.canonical_class,
        same_gene_neighbors: outcome.same_gene_neighbors.clone(),
    })
}

fn state_for_class(
    class: u32,
    outcome: &InternalOutcome,
    gene_keys: &GeneKeyIndex,
) -> Result<ClassState> {
    let mut state = public_state(outcome, gene_keys)?;
    state.counted = state.canonical_class == Some(class);
    Ok(state)
}

#[allow(clippy::too_many_arguments)]
fn build_report(
    x: &Extracted,
    class_cells: &[u32],
    gene_keys: &ComparisonGeneKeys,
    before_tuples: Vec<AssignmentTuple>,
    after_tuples: Vec<AssignmentTuple>,
    mut witness_candidates: BTreeMap<u32, MoleculeDraft>,
    direct: Vec<DirectClassEvidence>,
    decoded_chunks: usize,
    decoded_molecules: usize,
    archive_identity: ArchiveIdentity,
    options: CompareOptions,
) -> Result<AnnotationComparison> {
    let before = collapse_side(x, class_cells, &gene_keys.before, before_tuples)?;
    let after = collapse_side(x, class_cells, &gene_keys.after, after_tuples)?;

    let mut all_count_keys = BTreeSet::new();
    all_count_keys.extend(before.counts.keys().cloned());
    all_count_keys.extend(after.counts.keys().cloned());
    let mut count_deltas = Vec::new();
    for (cell, comparison_gene_id) in all_count_keys {
        let key = (cell, comparison_gene_id.clone());
        let before_count = before.counts.get(&key).copied().unwrap_or(0);
        let after_count = after.counts.get(&key).copied().unwrap_or(0);
        if before_count != after_count {
            count_deltas.push(CountDelta {
                cell,
                gene_id_before: gene_keys
                    .before
                    .original_by_key
                    .get(&comparison_gene_id)
                    .cloned(),
                gene_id_after: gene_keys
                    .after
                    .original_by_key
                    .get(&comparison_gene_id)
                    .cloned(),
                comparison_gene_id,
                before: before_count,
                after: after_count,
                delta: i64::from(after_count) - i64::from(before_count),
            });
        }
    }

    let mut transitions_by_class: Vec<Option<(ClassState, ClassState, Vec<TransitionCause>)>> =
        vec![None; x.n_classes as usize];
    let mut class_transitions = Vec::new();
    for class in 0..x.n_classes {
        let before_state =
            state_for_class(class, &before.outcomes[class as usize], &gene_keys.before)?;
        let after_state =
            state_for_class(class, &after.outcomes[class as usize], &gene_keys.after)?;
        let mut causes = direct[class as usize].causes.clone();
        if support_keys(&before_state.gene_support) != support_keys(&after_state.gene_support) {
            causes.insert(TransitionCause::ClassSupportChanged);
        }
        if annotation_order_tie_break_changed(&before_state, &after_state) {
            causes.insert(TransitionCause::AnnotationOrderTieBreakChanged);
        }
        if before_state.selected_comparison_gene_id != after_state.selected_comparison_gene_id {
            causes.insert(TransitionCause::ClassWinnerChanged);
        }
        if before_state.same_gene_neighbors != after_state.same_gene_neighbors {
            causes.insert(TransitionCause::CollapseNeighborhoodChanged);
        }
        if before_state.counted != after_state.counted
            || before_state.canonical_class != after_state.canonical_class
        {
            causes.insert(TransitionCause::CollapseOutcomeChanged);
        }
        if before_state.final_contribution() != after_state.final_contribution() {
            causes.insert(TransitionCause::FinalContributionChanged);
        }
        if causes.is_empty() {
            continue;
        }
        let kind = match (
            before_state.final_contribution(),
            after_state.final_contribution(),
        ) {
            (None, Some(_)) => ClassTransitionKind::GainedFinalCount,
            (Some(_), None) => ClassTransitionKind::LostFinalCount,
            (Some(left), Some(right)) if left != right => ClassTransitionKind::ReassignedFinalCount,
            _ => ClassTransitionKind::ChangedWithoutFinalCountDelta,
        };
        let causes: Vec<_> = causes.into_iter().collect();
        transitions_by_class[class as usize] =
            Some((before_state.clone(), after_state.clone(), causes.clone()));
        class_transitions.push(ClassTransition {
            cell: class_cells[class as usize],
            umi_class: class,
            kind,
            evidence: ClassEvidenceSummary {
                molecule_records: direct[class as usize].molecule_records,
                rows: direct[class as usize].rows,
                changed_rows: direct[class as usize].changed_rows,
                molecule_witnesses: 0,
                omitted_molecule_witnesses: direct[class as usize].molecule_records,
                changed_row_witnesses: 0,
                omitted_changed_row_witnesses: direct[class as usize].changed_rows,
            },
            before: before_state,
            after: after_state,
            causes,
        });
    }
    class_transitions.sort_by_key(|transition| (transition.cell, transition.umi_class));

    let changed_molecule_records = class_transitions.iter().try_fold(0u64, |sum, transition| {
        sum.checked_add(transition.evidence.molecule_records)
            .context("changed-class molecule-record count exceeds u64")
    })?;
    let mut molecule_witnesses = Vec::new();
    for transition in &mut class_transitions {
        if molecule_witnesses.len() >= options.max_molecule_witnesses {
            break;
        }
        let Some(draft) = witness_candidates.remove(&transition.umi_class) else {
            continue;
        };
        let Some((before_state, after_state, class_causes)) =
            transitions_by_class[draft.class as usize].as_ref()
        else {
            continue;
        };
        let mut causes: BTreeSet<_> = class_causes
            .iter()
            .copied()
            .filter(|cause| {
                !matches!(
                    cause,
                    TransitionCause::CandidateSetChanged | TransitionCause::RowAssignmentChanged
                )
            })
            .collect();
        causes.extend(draft.direct_causes);
        let changed_row_witnesses = u64::try_from(draft.changed_rows.len())
            .context("changed-row witness count exceeds u64")?;
        let changed_rows_omitted = draft
            .changed_rows_total
            .checked_sub(changed_row_witnesses)
            .context("row witnesses exceed changed rows")?;
        transition.evidence.molecule_witnesses = 1;
        transition.evidence.omitted_molecule_witnesses = transition
            .evidence
            .molecule_records
            .checked_sub(1)
            .context("class witness count exceeds molecule records")?;
        transition.evidence.changed_row_witnesses = changed_row_witnesses;
        transition.evidence.omitted_changed_row_witnesses = transition
            .evidence
            .changed_rows
            .checked_sub(changed_row_witnesses)
            .context("class row witnesses exceed changed rows")?;
        molecule_witnesses.push(MoleculeWitness {
            ordinal: draft.ordinal,
            cell: draft.cell,
            umi_class: draft.class,
            chrom: draft.chrom,
            anchor: draft.anchor,
            rows: draft.rows,
            changed_rows: draft.changed_rows,
            changed_rows_total: draft.changed_rows_total,
            changed_rows_omitted,
            before_class: before_state.clone(),
            after_class: after_state.clone(),
            causes: causes.into_iter().collect(),
        });
    }
    molecule_witnesses.sort_by_key(|molecule| molecule.ordinal);

    let evidence_rows = direct.iter().try_fold(0u64, |sum, evidence| {
        sum.checked_add(evidence.rows)
            .context("evidence-row total exceeds u64")
    })?;
    let decoded_molecules_u64 =
        u64::try_from(decoded_molecules).context("decoded molecule count exceeds u64")?;
    let unchanged_molecule_records = decoded_molecules_u64
        .checked_sub(changed_molecule_records)
        .context("changed molecule count exceeds decoded count")?;
    let witness_count =
        u64::try_from(molecule_witnesses.len()).context("molecule witness count exceeds u64")?;
    let molecule_witnesses_omitted = changed_molecule_records
        .checked_sub(witness_count)
        .context("molecule witnesses exceed changed records")?;

    Ok(AnnotationComparison {
        schema: REPORT_SCHEMA.to_string(),
        archive_identity,
        cell_barcodes: x
            .cells
            .iter()
            .map(|packed| {
                String::from_utf8(evidence_io::umi::unpack(*packed, 16))
                    .context("archive cell barcode was not valid ASCII")
            })
            .collect::<Result<Vec<_>>>()?,
        gene_key_policy: options.gene_key_policy,
        exactness: ExactnessBoundary {
            final_count_deltas_are_exact: true,
            class_state_ledger_is_complete: true,
            molecule_witnesses_are_bounded: true,
            causes_are_observed_state_transitions_not_counterfactual: true,
        },
        archive_passes: 1,
        decoded_chunks: u64::try_from(decoded_chunks).context("chunk count exceeds u64")?,
        decoded_molecules: decoded_molecules_u64,
        evidence_rows,
        changed_molecule_records,
        unchanged_molecule_records,
        molecule_witnesses_omitted,
        before: ReplaySummary {
            assigned_rows: before.assigned_rows,
            selected_classes: before.selected_classes,
            final_gene_umis: before.final_gene_umis,
        },
        after: ReplaySummary {
            assigned_rows: after.assigned_rows,
            selected_classes: after.selected_classes,
            final_gene_umis: after.final_gene_umis,
        },
        count_deltas,
        class_transitions,
        molecule_witnesses,
    })
}

#[cfg(test)]
fn publish_json_noclobber(output: &Path, report: &AnnotationComparison) -> Result<()> {
    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);
    let parent = output
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let file_name = output
        .file_name()
        .with_context(|| format!("output {} has no file name", output.display()))?
        .to_string_lossy();
    let mut temporary = None;
    for _ in 0..100 {
        let id = NEXT_TEMP.fetch_add(1, Ordering::Relaxed);
        let candidate = parent.join(format!(".{file_name}.{}.{}.tmp", std::process::id(), id));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(file) => {
                temporary = Some((candidate, file));
                break;
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error).with_context(|| format!("creating {}", candidate.display()))
            }
        }
    }
    let (temporary_path, file) = temporary.context("could not allocate a temporary output path")?;
    let result = (|| -> Result<()> {
        let mut writer = BufWriter::new(file);
        serde_json::to_writer_pretty(&mut writer, report)?;
        writer.write_all(b"\n")?;
        writer.flush()?;
        writer.get_ref().sync_all()?;
        drop(writer);
        fs::hard_link(&temporary_path, output).with_context(|| {
            format!(
                "publishing {} without replacing an existing path",
                output.display()
            )
        })?;
        Ok(())
    })();
    match result {
        Ok(()) => {
            // Publishing has already installed a complete hard link. An orphaned private
            // staging name is preferable to returning an error that implies no output exists.
            let _ = fs::remove_file(&temporary_path);
            Ok(())
        }
        Err(error) => {
            let _ = fs::remove_file(&temporary_path);
            Err(error)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archivecmd::write_archive;
    use crate::rows::{replay_rows_stranded, MolChain};
    use evidence_io::archive::Shape;
    use smallvec::smallvec;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEST: AtomicU64 = AtomicU64::new(0);

    struct TestDir(PathBuf);

    impl TestDir {
        fn new(label: &str) -> Self {
            let id = NEXT_TEST.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "gravlax-annotation-compare-{label}-{}-{id}",
                std::process::id()
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }

        fn path(&self, name: &str) -> PathBuf {
            self.0.join(name)
        }
    }

    impl Drop for TestDir {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).ok();
        }
    }

    fn annotation(dir: &TestDir, name: &str, body: &str) -> anno::Annotation {
        let path = dir.path(name);
        fs::write(&path, body).unwrap();
        anno::Annotation::from_gtf(&path).unwrap()
    }

    fn exon(chrom: &str, start: u32, end: u32, gene: &str, transcript: &str) -> String {
        format!(
            "{chrom}\ttest\texon\t{start}\t{end}\t.\t+\t.\tgene_id \"{gene}\"; \
             transcript_id \"{transcript}\"; gene_name \"{gene}\";\n"
        )
    }

    fn molecule(class: u32, pos: u32, weight: u32) -> MolRec {
        MolRec {
            cell: 0,
            umi_class: class,
            chrom: 0,
            strand_rev: false,
            chains: smallvec![MolChain {
                weight,
                reps: smallvec![(pos, 0)],
            }],
            mms: smallvec![],
        }
    }

    fn two_rep_molecule(class: u32, first: u32, second: u32, weight: u32) -> MolRec {
        MolRec {
            cell: 0,
            umi_class: class,
            chrom: 0,
            strand_rev: false,
            chains: smallvec![MolChain {
                weight,
                reps: smallvec![(first, 0), (second, 0)],
            }],
            mms: smallvec![],
        }
    }

    fn extracted(mols: Vec<MolRec>, edges: Vec<(u32, u32)>, n_classes: u32) -> Extracted {
        Extracted {
            mols,
            edges,
            cells: vec![0x4141_4141],
            shapes: vec![Shape {
                blocks: vec![(0, 20)],
            }],
            patterns: Vec::new(),
            n_classes,
            chrom_names: vec!["chr1".to_string()],
        }
    }

    fn write_test_archive(dir: &TestDir, name: &str, x: &Extracted, chunk_bp: u32) -> PathBuf {
        let path = dir.path(name);
        // The archive codec learns one global table per stream and rejects an entirely empty
        // training stream.  A scientifically inert record on an annotation-absent chromosome
        // exercises all streams without changing either replay's assignments or counts.
        let padding_chrom = x.chrom_names.len() as u32;
        let padding_pattern = x.patterns.len() as u32;
        let mut mols = x.mols.clone();
        mols.push(MolRec {
            cell: 0,
            umi_class: x.n_classes,
            chrom: padding_chrom,
            strand_rev: false,
            chains: smallvec![MolChain {
                weight: 1,
                reps: smallvec![(0, 0), (1, 0)],
            }],
            mms: smallvec![(0, 0, padding_pattern, 1)],
        });
        let mut patterns = x.patterns.clone();
        patterns.push(vec![PatAlt {
            chrom: padding_chrom,
            offset: 0,
            strand_flip: false,
            shape: SAME_SHAPE,
        }]);
        let mut chrom_names = x.chrom_names.clone();
        chrom_names.push("__codec_padding__".to_string());
        let padded = Extracted {
            mols,
            edges: x.edges.clone(),
            cells: x.cells.clone(),
            shapes: x.shapes.clone(),
            patterns,
            n_classes: x.n_classes + 1,
            chrom_names,
        };
        write_archive(&padded, &path, 1, chunk_bp, None).unwrap();
        path
    }

    fn reference_deltas(
        x: &Extracted,
        before: &anno::Annotation,
        after: &anno::Annotation,
        policy: GeneKeyPolicy,
    ) -> (Vec<CountDelta>, u64, u64, u64, u64) {
        let keys = ComparisonGeneKeys::new(before, after, policy).unwrap();
        let (before_counts, before_rows, before_total) =
            replay_rows_stranded(x, before, SoloStrand::Forward);
        let (after_counts, after_rows, after_total) =
            replay_rows_stranded(x, after, SoloStrand::Forward);
        let mut b = BTreeMap::<(u32, String), u32>::new();
        let mut a = BTreeMap::<(u32, String), u32>::new();
        for ((cell, gene), count) in before_counts {
            let key = keys
                .before
                .identity(gene)
                .unwrap()
                .comparison_gene_id
                .clone();
            *b.entry((cell, key)).or_default() += count;
        }
        for ((cell, gene), count) in after_counts {
            let key = keys
                .after
                .identity(gene)
                .unwrap()
                .comparison_gene_id
                .clone();
            *a.entry((cell, key)).or_default() += count;
        }
        let mut all = BTreeSet::new();
        all.extend(b.keys().cloned());
        all.extend(a.keys().cloned());
        let deltas = all
            .into_iter()
            .filter_map(|(cell, comparison_gene_id)| {
                let key = (cell, comparison_gene_id.clone());
                let before_count = b.get(&key).copied().unwrap_or(0);
                let after_count = a.get(&key).copied().unwrap_or(0);
                (before_count != after_count).then(|| CountDelta {
                    cell,
                    gene_id_before: keys
                        .before
                        .original_by_key
                        .get(&comparison_gene_id)
                        .cloned(),
                    gene_id_after: keys.after.original_by_key.get(&comparison_gene_id).cloned(),
                    comparison_gene_id,
                    before: before_count,
                    after: after_count,
                    delta: i64::from(after_count) - i64::from(before_count),
                })
            })
            .collect();
        (deltas, before_rows, before_total, after_rows, after_total)
    }

    fn reverse_kind(kind: ClassTransitionKind) -> ClassTransitionKind {
        match kind {
            ClassTransitionKind::GainedFinalCount => ClassTransitionKind::LostFinalCount,
            ClassTransitionKind::LostFinalCount => ClassTransitionKind::GainedFinalCount,
            ClassTransitionKind::ReassignedFinalCount => ClassTransitionKind::ReassignedFinalCount,
            ClassTransitionKind::ChangedWithoutFinalCountDelta => {
                ClassTransitionKind::ChangedWithoutFinalCountDelta
            }
        }
    }

    #[test]
    fn identical_annotations_are_zero_and_output_is_deterministic() {
        let dir = TestDir::new("identity");
        let gtf = format!(
            "{}{}",
            exon("chr1", 101, 180, "GA.1", "TA"),
            exon("chr1", 401, 480, "GX.1", "TX")
        );
        let a = annotation(&dir, "a.gtf", &gtf);
        let x = extracted(vec![molecule(0, 100, 2), molecule(1, 400, 1)], vec![], 2);
        let archive = write_test_archive(&dir, "x.aie", &x, 200);
        let first = compare_archive(&archive, &a, &a, CompareOptions::default()).unwrap();
        let second = compare_archive(&archive, &a, &a, CompareOptions::default()).unwrap();

        assert_eq!(first, second);
        assert_eq!(
            serde_json::to_vec(&first).unwrap(),
            serde_json::to_vec(&second).unwrap()
        );
        assert!(first.count_deltas.is_empty());
        assert!(first.class_transitions.is_empty());
        assert!(first.molecule_witnesses.is_empty());
        assert_eq!(
            first.cell_barcodes,
            vec![String::from_utf8(evidence_io::umi::unpack(x.cells[0], 16)).unwrap()]
        );
        assert_eq!(first.changed_molecule_records, 0);
        assert_eq!(first.unchanged_molecule_records, 3);
        assert_eq!(first.archive_passes, 1);
        assert_eq!(
            first.archive_identity.archive_version,
            evidence_io::format::VERSION
        );
        assert_eq!(
            first
                .archive_identity
                .rooted_content_commitment_hex
                .as_ref()
                .unwrap()
                .len(),
            64
        );
    }

    #[test]
    fn signed_deltas_match_independent_replays_and_swap_is_antisymmetric() {
        let dir = TestDir::new("antisymmetry");
        let a = annotation(&dir, "a.gtf", &exon("chr1", 101, 180, "GA.1", "TA"));
        let b = annotation(&dir, "b.gtf", &exon("chr1", 101, 180, "GB.7", "TB"));
        let x = extracted(vec![molecule(0, 100, 1)], vec![], 1);
        let archive = write_test_archive(&dir, "x.aie", &x, 500);
        let ab = compare_archive(&archive, &a, &b, CompareOptions::default()).unwrap();
        let ba = compare_archive(&archive, &b, &a, CompareOptions::default()).unwrap();
        let (expected, br, bt, ar, at) = reference_deltas(&x, &a, &b, GeneKeyPolicy::Unversioned);

        assert_eq!(ab.count_deltas, expected);
        assert_eq!(
            (ab.before.assigned_rows, ab.before.final_gene_umis),
            (br, bt)
        );
        assert_eq!((ab.after.assigned_rows, ab.after.final_gene_umis), (ar, at));
        assert_eq!(ab.count_deltas.len(), ba.count_deltas.len());
        for (left, right) in ab.count_deltas.iter().zip(&ba.count_deltas) {
            assert_eq!(left.comparison_gene_id, right.comparison_gene_id);
            assert_eq!(left.before, right.after);
            assert_eq!(left.after, right.before);
            assert_eq!(left.delta, -right.delta);
            assert_eq!(left.gene_id_before, right.gene_id_after);
            assert_eq!(left.gene_id_after, right.gene_id_before);
        }
        for (left, right) in ab.class_transitions.iter().zip(&ba.class_transitions) {
            assert_eq!((left.cell, left.umi_class), (right.cell, right.umi_class));
            assert_eq!(left.before, right.after);
            assert_eq!(left.after, right.before);
            assert_eq!(left.causes, right.causes);
            assert_eq!(left.kind, reverse_kind(right.kind));
        }
        for (left, right) in ab.molecule_witnesses.iter().zip(&ba.molecule_witnesses) {
            assert_eq!(left.ordinal, right.ordinal);
            assert_eq!(left.before_class, right.after_class);
            assert_eq!(left.after_class, right.before_class);
            assert_eq!(left.causes, right.causes);
            for (lr, rr) in left.changed_rows.iter().zip(&right.changed_rows) {
                assert_eq!(lr.before_candidates, rr.after_candidates);
                assert_eq!(lr.after_candidates, rr.before_candidates);
                assert_eq!(lr.causes, rr.causes);
            }
        }
    }

    #[test]
    fn reordered_identical_gtf_exposes_annotation_order_tie_break_artifact() {
        let dir = TestDir::new("annotation-order-tie");
        let gene_a = exon("chr1", 101, 180, "GA", "TA");
        let gene_b = exon("chr1", 301, 380, "GB", "TB");
        let before = annotation(&dir, "before.gtf", &format!("{gene_a}{gene_b}"));
        let after = annotation(&dir, "after.gtf", &format!("{gene_b}{gene_a}"));
        // Both retained rows have identical singleton assignments on the two sides. The only
        // difference is that equal support is resolved by each annotation's local gene order.
        let x = extracted(vec![molecule(0, 100, 1), molecule(0, 300, 1)], vec![], 1);
        let archive = write_test_archive(&dir, "x.aie", &x, 500);
        let report = compare_archive(&archive, &before, &after, CompareOptions::default()).unwrap();
        let (expected, before_rows, before_total, after_rows, after_total) =
            reference_deltas(&x, &before, &after, GeneKeyPolicy::Unversioned);

        assert_eq!(report.count_deltas, expected);
        assert_eq!(report.before.assigned_rows, before_rows);
        assert_eq!(report.after.assigned_rows, after_rows);
        assert_eq!(report.before.final_gene_umis, before_total);
        assert_eq!(report.after.final_gene_umis, after_total);
        assert_eq!(report.count_deltas.len(), 2);
        assert_eq!(report.count_deltas[0].comparison_gene_id, "GA");
        assert_eq!(report.count_deltas[0].delta, -1);
        assert_eq!(report.count_deltas[1].comparison_gene_id, "GB");
        assert_eq!(report.count_deltas[1].delta, 1);

        assert_eq!(report.class_transitions.len(), 1);
        let transition = &report.class_transitions[0];
        assert_eq!(transition.kind, ClassTransitionKind::ReassignedFinalCount);
        assert_eq!(
            transition.before.gene_support,
            transition.after.gene_support
        );
        assert_eq!(transition.before.selected_weight, 1);
        assert_eq!(transition.after.selected_weight, 1);
        assert_eq!(transition.evidence.changed_rows, 0);
        assert!(transition
            .causes
            .contains(&TransitionCause::AnnotationOrderTieBreakChanged));
        assert!(transition
            .causes
            .contains(&TransitionCause::ClassWinnerChanged));
        assert!(transition
            .causes
            .contains(&TransitionCause::FinalContributionChanged));
        assert!(!transition
            .causes
            .contains(&TransitionCause::ClassSupportChanged));
        assert!(report.molecule_witnesses.is_empty());
    }

    #[test]
    fn independent_final_collapse_exposes_nonlinear_compensation() {
        let dir = TestDir::new("collapse");
        let before_gtf = format!(
            "{}{}",
            exon("chr1", 101, 140, "GA", "TA"),
            exon("chr1", 301, 340, "GA", "TA")
        );
        let after_gtf = exon("chr1", 301, 340, "GA", "TA");
        let before = annotation(&dir, "before.gtf", &before_gtf);
        let after = annotation(&dir, "after.gtf", &after_gtf);
        let x = extracted(
            vec![molecule(0, 100, 2), molecule(1, 300, 1)],
            vec![(0, 1)],
            2,
        );
        let archive = write_test_archive(&dir, "x.aie", &x, 500);
        let report = compare_archive(&archive, &before, &after, CompareOptions::default()).unwrap();

        // Losing class 0's assignment makes class 1 a root, so the exact final GA count remains
        // one.  Subtracting pre-collapse class assignments would miss this compensation.
        assert!(report.count_deltas.is_empty());
        assert_eq!(report.before.final_gene_umis, 1);
        assert_eq!(report.after.final_gene_umis, 1);
        assert_eq!(report.class_transitions.len(), 2);
        let class0 = report
            .class_transitions
            .iter()
            .find(|t| t.umi_class == 0)
            .unwrap();
        let class1 = report
            .class_transitions
            .iter()
            .find(|t| t.umi_class == 1)
            .unwrap();
        assert_eq!(class0.kind, ClassTransitionKind::LostFinalCount);
        assert_eq!(class1.kind, ClassTransitionKind::GainedFinalCount);
        assert!(class0
            .causes
            .contains(&TransitionCause::CollapseNeighborhoodChanged));
        assert!(class1
            .causes
            .contains(&TransitionCause::CollapseOutcomeChanged));
        assert_eq!(class1.before.canonical_class, Some(0));
        assert_eq!(class1.after.canonical_class, Some(1));
    }

    #[test]
    fn multimapper_candidate_union_is_explained_not_subtracted() {
        let dir = TestDir::new("multimapper");
        let before_gtf = format!(
            "{}{}",
            exon("chr1", 101, 140, "GA", "TA"),
            exon("chr1", 301, 340, "GA", "TA")
        );
        let after_gtf = format!(
            "{}{}",
            exon("chr1", 101, 140, "GA", "TA"),
            exon("chr1", 301, 340, "GB", "TB")
        );
        let before = annotation(&dir, "before.gtf", &before_gtf);
        let after = annotation(&dir, "after.gtf", &after_gtf);
        let mut x = extracted(Vec::new(), vec![], 1);
        x.patterns.push(vec![
            PatAlt {
                chrom: 0,
                offset: 0,
                strand_flip: false,
                shape: SAME_SHAPE,
            },
            PatAlt {
                chrom: 0,
                offset: 200,
                strand_flip: false,
                shape: SAME_SHAPE,
            },
        ]);
        x.mols.push(MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: false,
            chains: smallvec![],
            mms: smallvec![(100, 0, 0, 1)],
        });
        let archive = write_test_archive(&dir, "x.aie", &x, 500);
        let report = compare_archive(&archive, &before, &after, CompareOptions::default()).unwrap();

        assert_eq!(report.count_deltas.len(), 1);
        assert_eq!(report.count_deltas[0].comparison_gene_id, "GA");
        assert_eq!(report.count_deltas[0].delta, -1);
        let row = &report.molecule_witnesses[0].changed_rows[0];
        assert_eq!(row.kind, EvidenceKind::MultimapperSignature);
        assert_eq!(candidate_keys(&row.before_candidates), vec!["GA"]);
        assert_eq!(candidate_keys(&row.after_candidates), vec!["GA", "GB"]);
        assert!(row.causes.contains(&TransitionCause::CandidateSetChanged));
        assert!(row.causes.contains(&TransitionCause::RowAssignmentChanged));
    }

    #[test]
    fn cross_chunk_two_representatives_and_witness_bounds_are_stable() {
        let dir = TestDir::new("cross-chunk");
        let before_gtf = format!(
            "{}{}{}{}",
            exon("chr1", 101, 140, "GA", "TA"),
            exon("chr1", 121, 160, "GA", "TA"),
            exon("chr1", 601, 640, "GC", "TC"),
            exon("chr1", 1101, 1140, "GB", "TB")
        );
        let after_gtf = format!(
            "{}{}",
            exon("chr1", 601, 640, "GC", "TC"),
            exon("chr1", 1101, 1140, "GB", "TB")
        );
        let before = annotation(&dir, "before.gtf", &before_gtf);
        let after = annotation(&dir, "after.gtf", &after_gtf);
        // Class 0 is introduced in chunk 0, appears with two representatives there, and is used
        // again by back-reference in chunk 2.  Class aggregation must remain global.
        let x = extracted(
            vec![
                two_rep_molecule(0, 100, 120, 3),
                molecule(1, 600, 1),
                molecule(0, 1100, 4),
            ],
            vec![],
            2,
        );
        let archive = write_test_archive(&dir, "x.aie", &x, 500);
        let options = CompareOptions {
            max_molecule_witnesses: 1,
            max_row_transitions_per_molecule: 1,
            ..CompareOptions::default()
        };
        let report = compare_archive(&archive, &before, &after, options).unwrap();
        let repeated = compare_archive(&archive, &before, &after, options).unwrap();
        let (expected, br, bt, ar, at) =
            reference_deltas(&x, &before, &after, GeneKeyPolicy::Unversioned);

        assert_eq!(report, repeated);
        assert_eq!(report.decoded_chunks, 4);
        assert_eq!(report.archive_passes, 1);
        assert_eq!(report.count_deltas, expected);
        assert_eq!(
            (report.before.assigned_rows, report.before.final_gene_umis),
            (br, bt)
        );
        assert_eq!(
            (report.after.assigned_rows, report.after.final_gene_umis),
            (ar, at)
        );
        let class0 = report
            .class_transitions
            .iter()
            .find(|t| t.umi_class == 0)
            .unwrap();
        assert_eq!(
            class0.before.selected_comparison_gene_id.as_deref(),
            Some("GA")
        );
        assert_eq!(class0.before.selected_weight, 6);
        assert_eq!(
            class0.after.selected_comparison_gene_id.as_deref(),
            Some("GB")
        );
        assert_eq!(class0.after.selected_weight, 4);
        assert_eq!(report.changed_molecule_records, 2);
        assert_eq!(report.molecule_witnesses.len(), 1);
        assert_eq!(report.molecule_witnesses_omitted, 1);
        assert_eq!(report.molecule_witnesses[0].ordinal, 0);
        assert_eq!(report.molecule_witnesses[0].changed_rows_total, 2);
        assert_eq!(report.molecule_witnesses[0].changed_rows.len(), 1);
        assert_eq!(report.molecule_witnesses[0].changed_rows_omitted, 1);
        assert_eq!(class0.evidence.changed_rows, 2);
        assert_eq!(class0.evidence.changed_row_witnesses, 1);
        assert_eq!(class0.evidence.omitted_changed_row_witnesses, 1);

        let zero = compare_archive(
            &archive,
            &before,
            &after,
            CompareOptions {
                max_molecule_witnesses: 0,
                max_row_transitions_per_molecule: 0,
                ..CompareOptions::default()
            },
        )
        .unwrap();
        assert!(zero.molecule_witnesses.is_empty());
        assert_eq!(
            zero.molecule_witnesses_omitted,
            zero.changed_molecule_records
        );
        assert!(zero.class_transitions.iter().all(|transition| {
            transition.evidence.molecule_witnesses == 0
                && transition.evidence.changed_row_witnesses == 0
        }));
    }

    #[test]
    fn global_witness_cap_prefers_earlier_ordinal_over_lower_class_id() {
        let dir = TestDir::new("witness-ordinal");
        let before = annotation(&dir, "before.gtf", &exon("chr1", 101, 800, "GA", "TA"));
        let after = annotation(&dir, "after.gtf", &exon("chr1", 5001, 5100, "GB", "TB"));
        // Class 0 is introduced by an unchanged row, then class 1 consumes one global witness
        // slot in the first chunk. In the next chunk, class 2 changes before a changed back-reference
        // to lower-ID class 0. With one slot left, archive order must win over the BTreeMap's
        // class-key order.
        let x = extracted(
            vec![
                molecule(0, 0, 1),
                molecule(1, 100, 1),
                molecule(2, 600, 1),
                molecule(0, 700, 1),
            ],
            vec![],
            3,
        );
        let archive = write_test_archive(&dir, "x.aie", &x, 500);
        let options = CompareOptions {
            max_molecule_witnesses: 2,
            max_row_transitions_per_molecule: 1,
            ..CompareOptions::default()
        };

        let report = compare_archive(&archive, &before, &after, options).unwrap();
        let repeated = compare_archive(&archive, &before, &after, options).unwrap();
        let (expected, ..) = reference_deltas(&x, &before, &after, GeneKeyPolicy::Unversioned);

        assert_eq!(report, repeated);
        assert_eq!(report.count_deltas, expected);
        assert_eq!(report.class_transitions.len(), 3);
        assert_eq!(report.changed_molecule_records, 4);
        assert_eq!(report.molecule_witnesses_omitted, 2);
        assert_eq!(
            report
                .molecule_witnesses
                .iter()
                .map(|witness| (witness.ordinal, witness.umi_class))
                .collect::<Vec<_>>(),
            vec![(1, 1), (2, 2)]
        );
        let later_low_id = report
            .class_transitions
            .iter()
            .find(|transition| transition.umi_class == 0)
            .unwrap();
        assert_eq!(later_low_id.evidence.molecule_witnesses, 0);
        assert_eq!(later_low_id.evidence.omitted_molecule_witnesses, 2);
    }

    #[test]
    fn unversioned_gene_keys_prevent_release_suffix_artifacts_and_reject_collisions() {
        let dir = TestDir::new("gene-keys");
        let before = annotation(&dir, "before.gtf", &exon("chr1", 101, 180, "ENSG1.4", "TA"));
        let after = annotation(&dir, "after.gtf", &exon("chr1", 101, 180, "ENSG1.5", "TB"));
        let x = extracted(vec![molecule(0, 100, 1)], vec![], 1);
        let archive = write_test_archive(&dir, "x.aie", &x, 500);

        let normalized =
            compare_archive(&archive, &before, &after, CompareOptions::default()).unwrap();
        assert!(normalized.count_deltas.is_empty());
        assert!(normalized.class_transitions.is_empty());

        let before_two = annotation(
            &dir,
            "before-two.gtf",
            &format!(
                "{}{}",
                exon("chr1", 101, 180, "ENSG1.4", "TA2"),
                exon("chr1", 301, 380, "ENSG1.4", "TA2")
            ),
        );
        let after_one = annotation(
            &dir,
            "after-one.gtf",
            &exon("chr1", 101, 180, "ENSG1.5", "TB2"),
        );
        let x_two = extracted(vec![molecule(0, 100, 1), molecule(1, 300, 1)], vec![], 2);
        let archive_two = write_test_archive(&dir, "two.aie", &x_two, 500);
        let normalized_delta = compare_archive(
            &archive_two,
            &before_two,
            &after_one,
            CompareOptions::default(),
        )
        .unwrap();
        assert_eq!(normalized_delta.count_deltas.len(), 1);
        assert_eq!(normalized_delta.count_deltas[0].comparison_gene_id, "ENSG1");
        assert_eq!(
            normalized_delta.count_deltas[0].gene_id_before.as_deref(),
            Some("ENSG1.4")
        );
        assert_eq!(
            normalized_delta.count_deltas[0].gene_id_after.as_deref(),
            Some("ENSG1.5")
        );
        assert_eq!(normalized_delta.count_deltas[0].delta, -1);

        let exact = compare_archive(
            &archive,
            &before,
            &after,
            CompareOptions {
                gene_key_policy: GeneKeyPolicy::Exact,
                ..CompareOptions::default()
            },
        )
        .unwrap();
        assert_eq!(exact.count_deltas.len(), 2);
        assert_eq!(exact.count_deltas[0].comparison_gene_id, "ENSG1.4");
        assert_eq!(
            exact.count_deltas[0].gene_id_before.as_deref(),
            Some("ENSG1.4")
        );
        assert_eq!(exact.count_deltas[0].gene_id_after, None);
        assert_eq!(exact.count_deltas[1].comparison_gene_id, "ENSG1.5");

        let collision_gtf = format!(
            "{}{}",
            exon("chr1", 101, 180, "ENSG2.1", "TC1"),
            exon("chr1", 301, 380, "ENSG2.2", "TC2")
        );
        let collision = annotation(&dir, "collision.gtf", &collision_gtf);
        let missing_archive = dir.path("does-not-exist.aie");
        let error = compare_archive(
            &missing_archive,
            &collision,
            &after,
            CompareOptions::default(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("normalization collision"));

        let empty_key = annotation(
            &dir,
            "empty-key.gtf",
            &exon("chr1", 101, 180, ".1", "TEMPTY"),
        );
        let error = compare_archive(
            &missing_archive,
            &empty_key,
            &after,
            CompareOptions::default(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("empty comparison key"));
    }

    #[test]
    fn output_publication_is_atomic_and_never_overwrites() {
        let dir = TestDir::new("publication");
        let a = annotation(&dir, "a.gtf", &exon("chr1", 101, 180, "GA", "TA"));
        let x = extracted(vec![molecule(0, 100, 1)], vec![], 1);
        let archive = write_test_archive(&dir, "x.aie", &x, 500);
        let output = dir.path("comparison.json");

        fs::write(&output, b"sentinel\n").unwrap();
        assert!(write_comparison_json_noclobber(
            &archive,
            &a,
            &a,
            CompareOptions::default(),
            &output,
        )
        .is_err());
        assert_eq!(fs::read(&output).unwrap(), b"sentinel\n");

        let invalid_archive = dir.path("invalid.aie");
        fs::write(&invalid_archive, b"not an archive").unwrap();
        let absent = dir.path("absent.json");
        assert!(write_comparison_json_noclobber(
            &invalid_archive,
            &a,
            &a,
            CompareOptions::default(),
            &absent,
        )
        .is_err());
        assert!(!absent.exists());
    }
}
