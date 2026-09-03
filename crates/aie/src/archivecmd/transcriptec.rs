//! Annotation-conditional transcript equivalence classes from the retained `.aie` quotient.
//!
//! This module deliberately stops at a deterministic compatibility relation.  It does not run
//! an abundance model and a singleton compatible-transcript set is not evidence that an entire
//! isoform was observed or phased.  For one stored record, compatible transcripts are unioned
//! over alternative placements. A record with no compatible transcript makes the class
//! explicitly no-compatible; otherwise record-level sets are intersected for the global `(cell,
//! UMI-class)` identity, including occurrences in different archive chunks.

use super::{decode_chunk, StreamingReplayArchive};
use crate::rows::{placement_from_parts_into, Extracted, MolRec, SAME_SHAPE};
use anyhow::{bail, Context, Result};
use evidence_io::{archive::Shape, Placement, Strand};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub const TRANSCRIPT_EQUIVALENCE_SCHEMA: &str = "gravlax.transcript-equivalence.v1";

/// Explicit interpretation limits attached to every result.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TranscriptEquivalenceSemantics {
    pub compatibility: &'static str,
    pub alternative_placements: &'static str,
    pub record_reduction: &'static str,
    pub exact_scope: &'static str,
    pub abundance_inferred: bool,
    pub full_isoform_phasing_claimed: bool,
}

impl Default for TranscriptEquivalenceSemantics {
    fn default() -> Self {
        Self {
            compatibility: "exonic block compatibility with exact annotated junctions under the selected strand policy",
            alternative_placements: "union compatible transcripts within each retained multimapper record",
            record_reduction: "fail closed when any retained record has no compatible transcript; otherwise intersect record-level sets across the global UMI class",
            exact_scope: "the retained archive quotient and supplied annotation",
            abundance_inferred: false,
            full_isoform_phasing_claimed: false,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TranscriptStrandPolicy {
    Forward,
    Reverse,
    Unstranded,
}

/// Coordinate window used only to decide whether a retained placement belongs to a query.
/// Compatibility remains an independent transcript-level decision.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub struct TranscriptSelectorLocus {
    pub contig: String,
    pub start: u32,
    pub end: u32,
}

impl From<anno::assign::SoloStrand> for TranscriptStrandPolicy {
    fn from(value: anno::assign::SoloStrand) -> Self {
        match value {
            anno::assign::SoloStrand::Forward => Self::Forward,
            anno::assign::SoloStrand::Reverse => Self::Reverse,
            anno::assign::SoloStrand::Unstranded => Self::Unstranded,
        }
    }
}

/// Required scientific identity and optional scope for transcript-EC derivation.
#[derive(Clone, Debug)]
pub struct TranscriptEquivalenceOptions {
    pub solo_strand: anno::assign::SoloStrand,
    /// BLAKE3 of the exact annotation resource content used for this derivation.
    pub annotation_content_digest: [u8; 32],
    /// Optional annotation transcript-index universe. Compatibility is evaluated only against
    /// these transcripts; indices are validated before the archive is opened.
    pub selected_transcript_indices: Option<BTreeSet<u32>>,
    /// Optional retained-placement windows. A class is selector-relevant when at least one
    /// retained placement block overlaps one of these windows, even if that placement is
    /// incompatible with every selected transcript.
    pub selector_loci: Option<Vec<TranscriptSelectorLocus>>,
    /// Scoped feature/locus callers can suppress classes that never have a compatible record in
    /// the selected universe. Classes with any selected-compatible record remain fail-closed if
    /// another retained record has no compatible transcript.
    pub include_classes_without_compatible_evidence: bool,
}

impl TranscriptEquivalenceOptions {
    pub fn new(solo_strand: anno::assign::SoloStrand, annotation_content_digest: [u8; 32]) -> Self {
        Self {
            solo_strand,
            annotation_content_digest,
            selected_transcript_indices: None,
            selector_loci: None,
            include_classes_without_compatible_evidence: true,
        }
    }

    pub fn with_transcript_universe(
        mut self,
        transcript_indices: impl IntoIterator<Item = u32>,
    ) -> Self {
        self.selected_transcript_indices = Some(transcript_indices.into_iter().collect());
        self
    }

    pub fn with_selector_loci(
        mut self,
        loci: impl IntoIterator<Item = TranscriptSelectorLocus>,
    ) -> Self {
        let mut loci: Vec<_> = loci.into_iter().collect();
        loci.sort();
        loci.dedup();
        self.selector_loci = Some(loci);
        self
    }

    pub fn only_classes_with_compatible_evidence(mut self) -> Self {
        self.include_classes_without_compatible_evidence = false;
        self
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TranscriptEquivalenceScope {
    pub annotation_content_blake3: String,
    /// `None` denotes the full annotation transcript universe.
    pub selected_transcript_ids: Option<Vec<String>>,
    /// `None` denotes an unselected full-archive derivation. When present, these canonical
    /// windows are the exact retained-placement predicate behind `selector_relevant`.
    pub selector_loci: Option<Vec<TranscriptSelectorLocus>>,
    pub include_classes_without_compatible_evidence: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TranscriptArchiveIdentity {
    pub archive_version: u32,
    pub rooted_commitment: Option<[u8; 32]>,
}

/// One interned, deterministic set of compatible transcript identifiers.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TranscriptSetSummary {
    /// Content-derived identity, independent of result ordering and filtering.
    pub ec_id: String,
    pub transcript_ids: Vec<String>,
    pub gene_ids: Vec<String>,
    pub class_count: u64,
    pub cell_count: u64,
    pub ambiguous: bool,
    pub complete_class_count: u64,
}

/// The transcript-set assignment and evidence-quality flags for one archive UMI class.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct UmiClassTranscriptSet {
    pub umi_class: u32,
    pub cell_id: u32,
    /// At least one retained placement block overlaps the explicit query selector. Always true
    /// for derivations without selector windows.
    pub selector_relevant: bool,
    /// Absent for fail-closed no-compatible outcomes and conflicts with an empty intersection.
    pub ec_id: Option<String>,
    pub retained_record_count: u64,
    pub represented_alignment_count: u64,
    pub compatible_record_count: u64,
    pub unmatched_record_count: u64,
    pub ambiguous: bool,
    /// At least one retained record has zero compatible transcripts (or the class has no retained
    /// records). This fail-closed state is distinct from a conflict among nonempty record sets.
    pub no_compatible_transcript: bool,
    pub conflict: bool,
    /// True only when every retained record has a compatible transcript, the nonempty sets have
    /// a nonempty intersection, and no middle placements were omitted by span-extreme reduction.
    /// This is completeness within the archive quotient, not full-transcript phasing.
    pub complete_within_archive_quotient: bool,
    /// False when a chain with more than two source alignments is represented by two span
    /// extremes. One representative at any weight is complete because all placements coincide.
    pub retained_representatives_complete: bool,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TranscriptCellSummary {
    pub cell_id: u32,
    pub barcode: String,
    pub class_count: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TranscriptEquivalenceTotals {
    pub archive_classes_scanned: u64,
    pub classes: u64,
    pub classes_filtered_out: u64,
    pub cells: u64,
    pub transcript_sets: u64,
    pub assigned_classes: u64,
    pub unassigned_classes: u64,
    pub ambiguous_classes: u64,
    pub no_compatible_transcript_classes: u64,
    pub conflicting_classes: u64,
    pub complete_classes: u64,
    pub incomplete_classes: u64,
    pub retained_records: u64,
    pub represented_alignments: u64,
}

/// Versioned, serialization-stable transcript equivalence result.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TranscriptEquivalenceReport {
    pub schema: &'static str,
    pub archive_version: u32,
    pub archive_root_blake3: Option<String>,
    pub strand_policy: TranscriptStrandPolicy,
    pub scope: TranscriptEquivalenceScope,
    pub semantics: TranscriptEquivalenceSemantics,
    pub transcript_sets: Vec<TranscriptSetSummary>,
    pub classes: Vec<UmiClassTranscriptSet>,
    pub cells: Vec<TranscriptCellSummary>,
    pub totals: TranscriptEquivalenceTotals,
}

impl TranscriptEquivalenceReport {
    /// Deterministic JSON: the result contains only sorted vectors and no serialized maps.
    pub fn to_pretty_json(&self) -> Result<String> {
        serde_json::to_string_pretty(self).context("serializing transcript equivalence report")
    }
}

/// Load either an AIC or GTF and derive transcript equivalence classes from an `.aie` archive.
pub fn derive_transcript_equivalence_classes(
    archive: &Path,
    annotation: &Path,
    annotation_identity: anno::intent::AnnotationIdentity,
    options: &TranscriptEquivalenceOptions,
) -> Result<TranscriptEquivalenceReport> {
    let bound = anno::intent::BoundAnnotation::from_path(annotation, annotation_identity)
        .map_err(anyhow::Error::new)?;
    let observed_digest = bound
        .identity()
        .digest
        .as_deref()
        .context("bound annotation did not report its observed content digest")?;
    let expected_digest = format!(
        "blake3:{}",
        blake3::Hash::from_bytes(options.annotation_content_digest).to_hex()
    );
    if observed_digest != expected_digest {
        bail!(
            "annotation content digest mismatch for {}: expected {}, observed {}",
            annotation.display(),
            expected_digest,
            observed_digest
        );
    }
    derive_transcript_equivalence_classes_with_annotation(archive, bound.annotation(), options)
}

/// Typed entry point for callers that already hold the compiled annotation model.
pub fn derive_transcript_equivalence_classes_with_annotation(
    archive: &Path,
    annotation: &anno::Annotation,
    options: &TranscriptEquivalenceOptions,
) -> Result<TranscriptEquivalenceReport> {
    validate_transcript_identifiers(annotation)?;
    validate_options(annotation, options)?;
    let archive_path = archive.to_path_buf();
    let replay = StreamingReplayArchive::open(&archive_path).with_context(|| {
        format!(
            "opening transcript-equivalence archive {}",
            archive.display()
        )
    })?;
    let archive_identity = TranscriptArchiveIdentity {
        archive_version: replay.reader.archive_version(),
        rooted_commitment: replay.reader.content_commitment().map(|root| root.digest),
    };
    let mut accumulator = ClassAccumulator::new(replay.x.n_classes)?;
    let batch_size = rayon::current_num_threads().max(1) * 2;
    let mut decoded_molecules = 0usize;

    for (batch_number, batch) in replay.chunks.chunks(batch_size).enumerate() {
        let first = batch_number * batch_size;
        let decoded: Vec<Vec<MolRec>> = batch
            .par_iter()
            .enumerate()
            .map(|(offset, info)| {
                let chunk_number = first + offset;
                let (compressed, raw_len) = replay
                    .reader
                    .read_compressed_at(&format!("c{chunk_number}"))?;
                let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                decode_chunk(&raw, info, Some(&replay.cell_of_class), &replay.rans_tables)
            })
            .collect::<Result<_>>()?;
        decoded_molecules = decoded_molecules
            .checked_add(decoded.iter().map(Vec::len).sum::<usize>())
            .context("decoded molecule count overflow")?;
        accumulator.add_chunks(&decoded, &replay.x, annotation, options)?;
    }

    if decoded_molecules != replay.n_mols {
        bail!(
            "molecule count mismatch: {decoded_molecules} decoded vs {} in archive metadata",
            replay.n_mols
        );
    }
    accumulator.finish(&replay.x.cells, annotation, options, &archive_identity)
}

/// In-memory counterpart used by scientific callers and focused regression tests.
pub fn derive_transcript_equivalence_classes_from_extracted(
    extracted: &Extracted,
    annotation: &anno::Annotation,
    options: &TranscriptEquivalenceOptions,
    archive_identity: &TranscriptArchiveIdentity,
) -> Result<TranscriptEquivalenceReport> {
    validate_transcript_identifiers(annotation)?;
    validate_options(annotation, options)?;
    let mut accumulator = ClassAccumulator::new(extracted.n_classes)?;
    accumulator.add_chunks(
        std::slice::from_ref(&extracted.mols),
        extracted,
        annotation,
        options,
    )?;
    accumulator.finish(&extracted.cells, annotation, options, archive_identity)
}

fn validate_transcript_identifiers(annotation: &anno::Annotation) -> Result<()> {
    if annotation.transcript_ids.len() != annotation.transcripts.len() {
        bail!(
            "transcript identifier metadata has {} entries for {} transcripts",
            annotation.transcript_ids.len(),
            annotation.transcripts.len()
        );
    }
    if !annotation.transcripts.is_empty() && annotation.transcript_ids.iter().all(Option::is_none) {
        bail!(
            "transcript equivalence classes require stable transcript IDs, but this annotation has none; AIC v1 omitted transcript IDs, so recompile the source GTF to AIC v2 or pass the GTF directly"
        );
    }
    let mut seen = BTreeSet::new();
    for (index, id) in annotation.transcript_ids.iter().enumerate() {
        let id = id.as_deref().with_context(|| {
            format!(
                "transcript equivalence classes require a stable ID for annotation transcript {index}"
            )
        })?;
        if id.is_empty() {
            bail!("annotation transcript {index} has an empty transcript ID");
        }
        if !seen.insert(id) {
            bail!("annotation contains duplicate transcript ID {id:?}");
        }
    }
    Ok(())
}

fn validate_options(
    annotation: &anno::Annotation,
    options: &TranscriptEquivalenceOptions,
) -> Result<()> {
    if let Some(selected) = &options.selected_transcript_indices {
        for &index in selected {
            if index as usize >= annotation.transcripts.len() {
                bail!(
                    "selected transcript index {index} is beyond the annotation's {} transcripts",
                    annotation.transcripts.len()
                );
            }
        }
    }
    if let Some(loci) = &options.selector_loci {
        for (index, locus) in loci.iter().enumerate() {
            if locus.contig.is_empty() {
                bail!("selector locus {index} has an empty contig");
            }
            if locus.start >= locus.end {
                bail!(
                    "selector locus {index} on {} has start {} greater than or equal to end {}",
                    locus.contig,
                    locus.start,
                    locus.end
                );
            }
        }
    }
    Ok(())
}

struct CompatibilityContext<'a> {
    extracted: &'a Extracted,
    annotation: &'a anno::Annotation,
    archive_to_annotation_chrom: Vec<Option<u32>>,
    solo_strand: anno::assign::SoloStrand,
    selected_transcript_indices: Option<&'a BTreeSet<u32>>,
    selector_windows_by_archive_chrom: Option<Vec<Vec<(u32, u32)>>>,
}

impl<'a> CompatibilityContext<'a> {
    fn new(
        extracted: &'a Extracted,
        annotation: &'a anno::Annotation,
        options: &'a TranscriptEquivalenceOptions,
    ) -> Result<Self> {
        let archive_to_annotation_chrom = extracted
            .chrom_names
            .iter()
            .map(|name| annotation.chrom_ids.get(name).copied())
            .collect();
        let selector_windows_by_archive_chrom = if let Some(loci) = &options.selector_loci {
            let archive_chroms: BTreeMap<&str, usize> = extracted
                .chrom_names
                .iter()
                .enumerate()
                .map(|(index, name)| (name.as_str(), index))
                .collect();
            let mut windows = vec![Vec::new(); extracted.chrom_names.len()];
            for locus in loci {
                let &chrom = archive_chroms.get(locus.contig.as_str()).with_context(|| {
                    format!(
                        "selector locus contig {:?} is absent from the archive chromosome dictionary",
                        locus.contig
                    )
                })?;
                windows[chrom].push((locus.start, locus.end));
            }
            for chrom_windows in &mut windows {
                chrom_windows.sort_unstable();
                chrom_windows.dedup();
            }
            Some(windows)
        } else {
            None
        };
        Ok(Self {
            extracted,
            annotation,
            archive_to_annotation_chrom,
            solo_strand: options.solo_strand,
            selected_transcript_indices: options.selected_transcript_indices.as_ref(),
            selector_windows_by_archive_chrom,
        })
    }

    fn placement_is_selector_relevant(
        &self,
        archive_chrom: u32,
        placement: &Placement,
    ) -> Result<bool> {
        let Some(windows_by_chrom) = &self.selector_windows_by_archive_chrom else {
            return Ok(true);
        };
        let windows = windows_by_chrom
            .get(archive_chrom as usize)
            .with_context(|| {
                format!(
                    "archive placement references chromosome {archive_chrom} beyond its chromosome dictionary"
                )
            })?;
        Ok(placement.blocks.iter().any(|block| {
            windows
                .iter()
                .any(|&(start, end)| block.start < end && block.end > start)
        }))
    }
}

struct CompatibilityScratch {
    placement: Placement,
    overlapping_transcripts: Vec<u32>,
    current_candidates: Vec<u32>,
    alternative_union: Vec<u32>,
}

impl Default for CompatibilityScratch {
    fn default() -> Self {
        Self {
            placement: Placement {
                chrom: 0,
                strand: Strand::Forward,
                blocks: Vec::new(),
                junctions: Vec::new(),
                nm: 0,
                score: 0,
                nh: 1,
                clip: (0, 0),
            },
            overlapping_transcripts: Vec::new(),
            current_candidates: Vec::new(),
            alternative_union: Vec::new(),
        }
    }
}

fn validate_shape_at(shape: &Shape, position: u32) -> Result<()> {
    if shape.blocks.is_empty() {
        bail!("archive placement shape has no blocks");
    }
    let mut previous_end = 0u32;
    for (index, &(offset, length)) in shape.blocks.iter().enumerate() {
        if length == 0 {
            bail!("archive placement shape block {index} has zero length");
        }
        if index > 0 && offset < previous_end {
            bail!("archive placement shape blocks overlap or are out of order");
        }
        let relative_end = offset
            .checked_add(length)
            .context("archive placement shape coordinate overflow")?;
        position
            .checked_add(relative_end)
            .context("archive placement genomic coordinate overflow")?;
        previous_end = relative_end;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn placement_candidates(
    archive_chrom: u32,
    position: u32,
    strand_rev: bool,
    shape_id: u32,
    nh: u16,
    context: &CompatibilityContext<'_>,
    scratch: &mut CompatibilityScratch,
) -> Result<bool> {
    scratch.current_candidates.clear();
    let shape = context
        .extracted
        .shapes
        .get(shape_id as usize)
        .with_context(|| format!("archive placement references missing shape {shape_id}"))?;
    validate_shape_at(shape, position)?;
    placement_from_parts_into(
        &mut scratch.placement,
        archive_chrom,
        position,
        strand_rev,
        shape,
        nh,
    );
    let selector_relevant =
        context.placement_is_selector_relevant(archive_chrom, &scratch.placement)?;
    let annotation_chrom = context
        .archive_to_annotation_chrom
        .get(archive_chrom as usize)
        .with_context(|| {
            format!(
                "archive placement references chromosome {archive_chrom} beyond its chromosome dictionary"
            )
        })?;
    let Some(annotation_chrom) = *annotation_chrom else {
        return Ok(selector_relevant);
    };
    context.annotation.overlapping_into(
        annotation_chrom,
        scratch.placement.start(),
        scratch.placement.end(),
        &mut scratch.overlapping_transcripts,
    );
    for &transcript_index in &scratch.overlapping_transcripts {
        if context
            .selected_transcript_indices
            .is_some_and(|selected| !selected.contains(&transcript_index))
        {
            continue;
        }
        let transcript = context
            .annotation
            .transcripts
            .get(transcript_index as usize)
            .with_context(|| {
                format!("annotation overlap index references transcript {transcript_index}")
            })?;
        if context
            .solo_strand
            .accepts(strand_rev, transcript.strand_rev)
            && anno::assign::align_vs_transcript(&scratch.placement, transcript)
                == anno::assign::Vs::Concordant
        {
            scratch.current_candidates.push(transcript_index);
        }
    }
    scratch.current_candidates.sort_unstable();
    scratch.current_candidates.dedup();
    Ok(selector_relevant)
}

fn checked_alternative_position(anchor: u32, offset: i64) -> Result<u32> {
    let position = i64::from(anchor)
        .checked_add(offset)
        .context("multimapper alternative coordinate overflow")?;
    u32::try_from(position).context("multimapper alternative coordinate is negative or exceeds u32")
}

#[derive(Debug)]
struct MoleculeObservation {
    umi_class: u32,
    cell_id: u32,
    candidates: Option<Vec<u32>>,
    retained_record_count: u64,
    represented_alignment_count: u64,
    compatible_record_count: u64,
    unmatched_record_count: u64,
    conflict: bool,
    selector_relevant: bool,
    retained_representatives_complete: bool,
}

struct CandidateReduction {
    candidates: Option<Vec<u32>>,
    retained_record_count: u64,
    compatible_record_count: u64,
    unmatched_record_count: u64,
    conflict: bool,
}

impl CandidateReduction {
    fn new() -> Self {
        Self {
            candidates: None,
            retained_record_count: 0,
            compatible_record_count: 0,
            unmatched_record_count: 0,
            conflict: false,
        }
    }

    fn observe(&mut self, candidates: &[u32]) -> Result<()> {
        self.retained_record_count = self
            .retained_record_count
            .checked_add(1)
            .context("retained record count overflow")?;
        if candidates.is_empty() {
            self.unmatched_record_count = self
                .unmatched_record_count
                .checked_add(1)
                .context("unmatched record count overflow")?;
            return Ok(());
        }
        self.compatible_record_count = self
            .compatible_record_count
            .checked_add(1)
            .context("compatible record count overflow")?;
        self.candidates = Some(match self.candidates.take() {
            None => candidates.to_vec(),
            Some(previous) => {
                let intersection = intersect_sorted(&previous, candidates);
                if previous.len() > 0 && intersection.is_empty() {
                    self.conflict = true;
                }
                intersection
            }
        });
        Ok(())
    }
}

fn classify_molecule(
    molecule: &MolRec,
    context: &CompatibilityContext<'_>,
    scratch: &mut CompatibilityScratch,
) -> Result<MoleculeObservation> {
    let mut reduction = CandidateReduction::new();
    let mut represented_alignment_count = 0u64;
    let mut retained_representatives_complete = true;
    let mut selector_relevant = context.selector_windows_by_archive_chrom.is_none();

    for chain in &molecule.chains {
        if chain.weight == 0 {
            bail!("UMI class {} has a zero-weight chain", molecule.umi_class);
        }
        represented_alignment_count = represented_alignment_count
            .checked_add(u64::from(chain.weight))
            .context("represented alignment count overflow")?;
        if chain.reps.is_empty() {
            retained_representatives_complete = false;
            continue;
        }
        if chain.reps.len() > 2 || chain.reps.len() > chain.weight as usize {
            bail!(
                "UMI class {} has {} representatives for a weight-{} chain",
                molecule.umi_class,
                chain.reps.len(),
                chain.weight
            );
        }
        if chain.reps.len() == 2 && chain.weight > 2 {
            retained_representatives_complete = false;
        }
        for &(position, shape_id) in &chain.reps {
            selector_relevant |= placement_candidates(
                molecule.chrom,
                position,
                molecule.strand_rev,
                shape_id,
                1,
                context,
                scratch,
            )?;
            reduction.observe(&scratch.current_candidates)?;
        }
    }

    for &(anchor, anchor_shape, pattern_id, weight) in &molecule.mms {
        if weight == 0 {
            bail!(
                "UMI class {} has a zero-weight multimapper record",
                molecule.umi_class
            );
        }
        represented_alignment_count = represented_alignment_count
            .checked_add(u64::from(weight))
            .context("represented alignment count overflow")?;
        let alternatives = context
            .extracted
            .patterns
            .get(pattern_id as usize)
            .with_context(|| format!("archive record references missing pattern {pattern_id}"))?;
        scratch.alternative_union.clear();
        for alternative in alternatives {
            let position = checked_alternative_position(anchor, alternative.offset)?;
            let shape_id = if alternative.shape == SAME_SHAPE {
                anchor_shape
            } else {
                alternative.shape
            };
            selector_relevant |= placement_candidates(
                alternative.chrom,
                position,
                molecule.strand_rev != alternative.strand_flip,
                shape_id,
                2,
                context,
                scratch,
            )?;
            scratch
                .alternative_union
                .extend_from_slice(&scratch.current_candidates);
        }
        scratch.alternative_union.sort_unstable();
        scratch.alternative_union.dedup();
        reduction.observe(&scratch.alternative_union)?;
    }

    if reduction.retained_record_count == 0 {
        retained_representatives_complete = false;
    }
    Ok(MoleculeObservation {
        umi_class: molecule.umi_class,
        cell_id: molecule.cell,
        candidates: reduction.candidates,
        retained_record_count: reduction.retained_record_count,
        represented_alignment_count,
        compatible_record_count: reduction.compatible_record_count,
        unmatched_record_count: reduction.unmatched_record_count,
        conflict: reduction.conflict,
        selector_relevant,
        retained_representatives_complete,
    })
}

fn intersect_sorted(left: &[u32], right: &[u32]) -> Vec<u32> {
    let mut output = Vec::with_capacity(left.len().min(right.len()));
    let (mut i, mut j) = (0usize, 0usize);
    while i < left.len() && j < right.len() {
        match left[i].cmp(&right[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                output.push(left[i]);
                i += 1;
                j += 1;
            }
        }
    }
    output
}

struct ClassState {
    cell_id: u32,
    candidates: Option<Vec<u32>>,
    retained_record_count: u64,
    represented_alignment_count: u64,
    compatible_record_count: u64,
    unmatched_record_count: u64,
    conflict: bool,
    selector_relevant: bool,
    retained_representatives_complete: bool,
}

impl ClassState {
    fn from_observation(observation: MoleculeObservation) -> Self {
        Self {
            cell_id: observation.cell_id,
            candidates: observation.candidates,
            retained_record_count: observation.retained_record_count,
            represented_alignment_count: observation.represented_alignment_count,
            compatible_record_count: observation.compatible_record_count,
            unmatched_record_count: observation.unmatched_record_count,
            conflict: observation.conflict,
            selector_relevant: observation.selector_relevant,
            retained_representatives_complete: observation.retained_representatives_complete,
        }
    }

    fn merge(&mut self, observation: MoleculeObservation, umi_class: u32) -> Result<()> {
        if self.cell_id != observation.cell_id {
            bail!(
                "UMI class {umi_class} appears in cells {} and {}",
                self.cell_id,
                observation.cell_id
            );
        }
        self.retained_record_count = self
            .retained_record_count
            .checked_add(observation.retained_record_count)
            .context("retained record count overflow")?;
        self.represented_alignment_count = self
            .represented_alignment_count
            .checked_add(observation.represented_alignment_count)
            .context("represented alignment count overflow")?;
        self.compatible_record_count = self
            .compatible_record_count
            .checked_add(observation.compatible_record_count)
            .context("compatible record count overflow")?;
        self.unmatched_record_count = self
            .unmatched_record_count
            .checked_add(observation.unmatched_record_count)
            .context("unmatched record count overflow")?;
        self.conflict |= observation.conflict;
        self.selector_relevant |= observation.selector_relevant;
        self.retained_representatives_complete &= observation.retained_representatives_complete;
        if let Some(incoming) = observation.candidates {
            self.candidates = Some(match self.candidates.take() {
                None => incoming,
                Some(previous) => {
                    let intersection = intersect_sorted(&previous, &incoming);
                    if !previous.is_empty() && !incoming.is_empty() && intersection.is_empty() {
                        self.conflict = true;
                    }
                    intersection
                }
            });
        }
        Ok(())
    }
}

struct ClassAccumulator {
    classes: Vec<Option<ClassState>>,
}

impl ClassAccumulator {
    fn new(n_classes: u32) -> Result<Self> {
        let n = usize::try_from(n_classes).context("class count exceeds usize")?;
        Ok(Self {
            classes: std::iter::repeat_with(|| None).take(n).collect(),
        })
    }

    fn add_chunks(
        &mut self,
        chunks: &[Vec<MolRec>],
        extracted: &Extracted,
        annotation: &anno::Annotation,
        options: &TranscriptEquivalenceOptions,
    ) -> Result<()> {
        let context = CompatibilityContext::new(extracted, annotation, options)?;
        let observations: Vec<MoleculeObservation> = chunks
            .par_iter()
            .flat_map_iter(|chunk| chunk.iter())
            .map_init(CompatibilityScratch::default, |scratch, molecule| {
                classify_molecule(molecule, &context, scratch)
            })
            .collect::<Result<_>>()?;
        for observation in observations {
            let class_index = observation.umi_class as usize;
            let declared_classes = self.classes.len();
            let slot = self.classes.get_mut(class_index).with_context(|| {
                format!(
                    "archive molecule references class {} beyond declared class count {}",
                    observation.umi_class, declared_classes
                )
            })?;
            match slot {
                Some(state) => state.merge(observation, class_index as u32)?,
                None => *slot = Some(ClassState::from_observation(observation)),
            }
        }
        Ok(())
    }

    fn finish(
        self,
        cells: &[u32],
        annotation: &anno::Annotation,
        options: &TranscriptEquivalenceOptions,
        archive_identity: &TranscriptArchiveIdentity,
    ) -> Result<TranscriptEquivalenceReport> {
        let archive_classes_scanned = self.classes.len();
        let mut resolved_classes = Vec::with_capacity(self.classes.len());
        for (umi_class, state) in self.classes.into_iter().enumerate() {
            let state = state.with_context(|| {
                format!("archive declares UMI class {umi_class}, but no record carries that class")
            })?;
            if !options.include_classes_without_compatible_evidence
                && state.compatible_record_count == 0
            {
                continue;
            }
            let no_compatible_transcript = state.retained_record_count == 0
                || state.compatible_record_count == 0
                || state.unmatched_record_count > 0;
            let mut transcript_ids = Vec::new();
            let mut gene_ids = BTreeSet::new();
            let reported_candidates = if no_compatible_transcript || state.conflict {
                &[][..]
            } else {
                state.candidates.as_deref().unwrap_or_default()
            };
            for transcript_index in reported_candidates {
                let transcript = annotation
                    .transcripts
                    .get(*transcript_index as usize)
                    .with_context(|| {
                        format!("candidate transcript index {transcript_index} is out of range")
                    })?;
                let transcript_id = annotation.transcript_ids[*transcript_index as usize]
                    .as_ref()
                    .context("validated transcript ID unexpectedly missing")?;
                transcript_ids.push(transcript_id.clone());
                let gene_id = annotation
                    .gene_ids
                    .get(transcript.gene as usize)
                    .with_context(|| {
                        format!(
                            "annotation transcript {transcript_id:?} references missing gene {}",
                            transcript.gene
                        )
                    })?;
                gene_ids.insert(gene_id.clone());
            }
            transcript_ids.sort();
            transcript_ids.dedup();
            let ambiguous = transcript_ids.len() > 1;
            let unassigned = transcript_ids.is_empty();
            let complete_within_archive_quotient = !unassigned
                && !no_compatible_transcript
                && !state.conflict
                && state.retained_representatives_complete;
            resolved_classes.push(ResolvedClass {
                umi_class: u32::try_from(umi_class).context("UMI class exceeds u32")?,
                cell_id: state.cell_id,
                selector_relevant: state.selector_relevant,
                transcript_ids,
                gene_ids: gene_ids.into_iter().collect(),
                retained_record_count: state.retained_record_count,
                represented_alignment_count: state.represented_alignment_count,
                compatible_record_count: state.compatible_record_count,
                unmatched_record_count: state.unmatched_record_count,
                ambiguous,
                no_compatible_transcript,
                conflict: state.conflict,
                complete_within_archive_quotient,
                retained_representatives_complete: state.retained_representatives_complete,
            });
        }

        assemble_report(
            resolved_classes,
            cells,
            annotation,
            options,
            archive_identity,
            archive_classes_scanned,
        )
    }
}

struct ResolvedClass {
    umi_class: u32,
    cell_id: u32,
    selector_relevant: bool,
    transcript_ids: Vec<String>,
    gene_ids: Vec<String>,
    retained_record_count: u64,
    represented_alignment_count: u64,
    compatible_record_count: u64,
    unmatched_record_count: u64,
    ambiguous: bool,
    no_compatible_transcript: bool,
    conflict: bool,
    complete_within_archive_quotient: bool,
    retained_representatives_complete: bool,
}

#[derive(Default)]
struct TranscriptSetAccumulator {
    gene_ids: Vec<String>,
    cells: BTreeSet<u32>,
    class_count: u64,
    complete_class_count: u64,
}

fn update_length_framed(hasher: &mut blake3::Hasher, bytes: &[u8]) -> Result<()> {
    let length = u64::try_from(bytes.len()).context("content-address frame exceeds u64")?;
    hasher.update(&length.to_le_bytes());
    hasher.update(bytes);
    Ok(())
}

fn transcript_ec_id(
    annotation_content_digest: &[u8; 32],
    sorted_transcript_ids: &[String],
) -> Result<String> {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"gravlax.transcript-equivalence-class.v1\0");
    update_length_framed(&mut hasher, annotation_content_digest)?;
    let count = u64::try_from(sorted_transcript_ids.len())
        .context("transcript equivalence class exceeds u64 members")?;
    hasher.update(&count.to_le_bytes());
    for transcript_id in sorted_transcript_ids {
        update_length_framed(&mut hasher, transcript_id.as_bytes())?;
    }
    Ok(hasher.finalize().to_hex().to_string())
}

fn assemble_report(
    resolved_classes: Vec<ResolvedClass>,
    packed_cells: &[u32],
    annotation: &anno::Annotation,
    options: &TranscriptEquivalenceOptions,
    archive_identity: &TranscriptArchiveIdentity,
    archive_classes_scanned: usize,
) -> Result<TranscriptEquivalenceReport> {
    let mut sets: BTreeMap<Vec<String>, TranscriptSetAccumulator> = BTreeMap::new();
    let mut cell_class_counts: BTreeMap<u32, u64> = BTreeMap::new();
    for class in &resolved_classes {
        if class.cell_id as usize >= packed_cells.len() {
            bail!(
                "UMI class {} references cell {} beyond the {}-cell dictionary",
                class.umi_class,
                class.cell_id,
                packed_cells.len()
            );
        }
        if !class.transcript_ids.is_empty() {
            let set = sets.entry(class.transcript_ids.clone()).or_default();
            if set.gene_ids.is_empty() {
                set.gene_ids = class.gene_ids.clone();
            } else if set.gene_ids != class.gene_ids {
                bail!("identical transcript sets resolved to inconsistent gene sets");
            }
            set.class_count = set
                .class_count
                .checked_add(1)
                .context("class count overflow")?;
            set.cells.insert(class.cell_id);
            set.complete_class_count = set
                .complete_class_count
                .checked_add(u64::from(class.complete_within_archive_quotient))
                .context("complete class count overflow")?;
        }
        let cell_count = cell_class_counts.entry(class.cell_id).or_default();
        *cell_count = cell_count
            .checked_add(1)
            .context("cell class count overflow")?;
    }

    let mut set_ids = BTreeMap::new();
    let mut id_contents: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut transcript_sets = Vec::with_capacity(sets.len());
    for (transcript_ids, aggregate) in sets {
        let ec_id = transcript_ec_id(&options.annotation_content_digest, &transcript_ids)?;
        if let Some(previous) = id_contents.insert(ec_id.clone(), transcript_ids.clone()) {
            if previous != transcript_ids {
                bail!("BLAKE3 transcript-equivalence identity collision");
            }
        }
        set_ids.insert(transcript_ids.clone(), ec_id.clone());
        transcript_sets.push(TranscriptSetSummary {
            ec_id,
            ambiguous: transcript_ids.len() > 1,
            transcript_ids,
            gene_ids: aggregate.gene_ids,
            class_count: aggregate.class_count,
            cell_count: u64::try_from(aggregate.cells.len()).context("cell count exceeds u64")?,
            complete_class_count: aggregate.complete_class_count,
        });
    }
    transcript_sets.sort_by(|left, right| left.ec_id.cmp(&right.ec_id));

    let mut totals = TranscriptEquivalenceTotals {
        archive_classes_scanned: u64::try_from(archive_classes_scanned)
            .context("archive class count exceeds u64")?,
        classes: u64::try_from(resolved_classes.len()).context("class count exceeds u64")?,
        classes_filtered_out: u64::try_from(
            archive_classes_scanned
                .checked_sub(resolved_classes.len())
                .context("filtered class count underflow")?,
        )
        .context("filtered class count exceeds u64")?,
        cells: u64::try_from(cell_class_counts.len()).context("cell count exceeds u64")?,
        transcript_sets: u64::try_from(transcript_sets.len())
            .context("transcript set count exceeds u64")?,
        assigned_classes: 0,
        unassigned_classes: 0,
        ambiguous_classes: 0,
        no_compatible_transcript_classes: 0,
        conflicting_classes: 0,
        complete_classes: 0,
        incomplete_classes: 0,
        retained_records: 0,
        represented_alignments: 0,
    };
    let mut classes = Vec::with_capacity(resolved_classes.len());
    for class in resolved_classes {
        let unassigned = class.transcript_ids.is_empty();
        let ec_id = if unassigned {
            None
        } else {
            Some(
                set_ids
                    .get(&class.transcript_ids)
                    .context("transcript set was not interned")?
                    .clone(),
            )
        };
        totals.assigned_classes = totals
            .assigned_classes
            .checked_add(u64::from(!unassigned))
            .context("assigned class count overflow")?;
        totals.unassigned_classes = totals
            .unassigned_classes
            .checked_add(u64::from(unassigned))
            .context("unassigned class count overflow")?;
        totals.ambiguous_classes = totals
            .ambiguous_classes
            .checked_add(u64::from(class.ambiguous))
            .context("ambiguous class count overflow")?;
        totals.no_compatible_transcript_classes = totals
            .no_compatible_transcript_classes
            .checked_add(u64::from(class.no_compatible_transcript))
            .context("no-compatible-transcript class count overflow")?;
        totals.conflicting_classes = totals
            .conflicting_classes
            .checked_add(u64::from(class.conflict))
            .context("conflicting class count overflow")?;
        totals.complete_classes = totals
            .complete_classes
            .checked_add(u64::from(class.complete_within_archive_quotient))
            .context("complete class count overflow")?;
        totals.incomplete_classes = totals
            .incomplete_classes
            .checked_add(u64::from(!class.complete_within_archive_quotient))
            .context("incomplete class count overflow")?;
        totals.retained_records = totals
            .retained_records
            .checked_add(class.retained_record_count)
            .context("retained record total overflow")?;
        totals.represented_alignments = totals
            .represented_alignments
            .checked_add(class.represented_alignment_count)
            .context("represented alignment total overflow")?;
        classes.push(UmiClassTranscriptSet {
            umi_class: class.umi_class,
            cell_id: class.cell_id,
            selector_relevant: class.selector_relevant,
            ec_id,
            retained_record_count: class.retained_record_count,
            represented_alignment_count: class.represented_alignment_count,
            compatible_record_count: class.compatible_record_count,
            unmatched_record_count: class.unmatched_record_count,
            ambiguous: class.ambiguous,
            no_compatible_transcript: class.no_compatible_transcript,
            conflict: class.conflict,
            complete_within_archive_quotient: class.complete_within_archive_quotient,
            retained_representatives_complete: class.retained_representatives_complete,
        });
    }

    let cells = cell_class_counts
        .into_iter()
        .map(|(cell_id, class_count)| {
            let barcode =
                String::from_utf8(evidence_io::umi::unpack(packed_cells[cell_id as usize], 16))
                    .expect("packed cell barcode always decodes to ASCII");
            TranscriptCellSummary {
                cell_id,
                barcode,
                class_count,
            }
        })
        .collect();

    let set_class_total = transcript_sets.iter().try_fold(0u64, |sum, set| {
        sum.checked_add(set.class_count)
            .context("transcript-set class total overflow")
    })?;
    if set_class_total != totals.assigned_classes
        || totals.classes + totals.classes_filtered_out != totals.archive_classes_scanned
        || totals.assigned_classes + totals.unassigned_classes != totals.classes
        || totals.complete_classes + totals.incomplete_classes != totals.classes
    {
        bail!("transcript equivalence conservation invariant failed");
    }

    let selected_transcript_ids = options
        .selected_transcript_indices
        .as_ref()
        .map(|selected| {
            let mut ids = selected
                .iter()
                .map(|&index| {
                    annotation.transcript_ids[index as usize]
                        .clone()
                        .context("validated selected transcript ID unexpectedly missing")
                })
                .collect::<Result<Vec<_>>>()?;
            ids.sort();
            Ok::<Vec<String>, anyhow::Error>(ids)
        })
        .transpose()?;

    Ok(TranscriptEquivalenceReport {
        schema: TRANSCRIPT_EQUIVALENCE_SCHEMA,
        archive_version: archive_identity.archive_version,
        archive_root_blake3: archive_identity
            .rooted_commitment
            .map(|digest| blake3::Hash::from_bytes(digest).to_hex().to_string()),
        strand_policy: options.solo_strand.into(),
        scope: TranscriptEquivalenceScope {
            annotation_content_blake3: format!(
                "blake3:{}",
                blake3::Hash::from_bytes(options.annotation_content_digest).to_hex()
            ),
            selected_transcript_ids,
            selector_loci: options.selector_loci.clone(),
            include_classes_without_compatible_evidence: options
                .include_classes_without_compatible_evidence,
        },
        semantics: TranscriptEquivalenceSemantics::default(),
        transcript_sets,
        classes,
        cells,
        totals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rows::{MolChain, PatAlt};
    use evidence_io::archive::Shape;
    use smallvec::{smallvec, SmallVec};
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_FILE: AtomicU64 = AtomicU64::new(0);

    struct TestFile(PathBuf);

    impl TestFile {
        fn new(extension: &str) -> Self {
            let serial = NEXT_FILE.fetch_add(1, Ordering::Relaxed);
            Self(std::env::temp_dir().join(format!(
                "gravlax-transcriptec-{}-{serial}.{extension}",
                std::process::id()
            )))
        }

        fn path(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for TestFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    fn exon(chrom: &str, start: u32, end: u32, strand: char, gene: &str, tx: &str) -> String {
        format!(
            "{chrom}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gene}\"; gene_name \"{gene}\"; transcript_id \"{tx}\";\n"
        )
    }

    fn annotation(records: &[String]) -> anno::Annotation {
        let gtf = TestFile::new("gtf");
        std::fs::write(gtf.path(), records.concat()).unwrap();
        anno::Annotation::from_gtf(gtf.path()).unwrap()
    }

    fn packed_cells(n: usize) -> Vec<u32> {
        [
            b"AAAAAAAAAAAAAAAA".as_slice(),
            b"CCCCCCCCCCCCCCCC".as_slice(),
            b"GGGGGGGGGGGGGGGG".as_slice(),
        ]
        .into_iter()
        .take(n)
        .map(|barcode| evidence_io::umi::pack(barcode).unwrap())
        .collect()
    }

    fn one_block(length: u32) -> Shape {
        Shape {
            blocks: vec![(0, length)],
        }
    }

    fn unique_molecule(
        cell: u32,
        class: u32,
        chrom: u32,
        strand_rev: bool,
        position: u32,
        shape: u32,
    ) -> MolRec {
        MolRec {
            cell,
            umi_class: class,
            chrom,
            strand_rev,
            chains: smallvec![MolChain {
                weight: 1,
                reps: smallvec![(position, shape)]
            }],
            mms: SmallVec::new(),
        }
    }

    fn extracted(
        mols: Vec<MolRec>,
        shapes: Vec<Shape>,
        patterns: Vec<Vec<PatAlt>>,
        n_classes: u32,
        n_cells: usize,
        chrom_names: &[&str],
    ) -> Extracted {
        Extracted {
            mols,
            edges: Vec::new(),
            cells: packed_cells(n_cells),
            shapes,
            patterns,
            n_classes,
            chrom_names: chrom_names.iter().map(|name| (*name).to_owned()).collect(),
        }
    }

    fn fixture_options(solo_strand: anno::assign::SoloStrand) -> TranscriptEquivalenceOptions {
        TranscriptEquivalenceOptions::new(
            solo_strand,
            *blake3::hash(b"transcriptec-test-annotation").as_bytes(),
        )
    }

    fn fixture_archive_identity() -> TranscriptArchiveIdentity {
        TranscriptArchiveIdentity {
            archive_version: evidence_io::format::VERSION,
            rooted_commitment: Some(*blake3::hash(b"transcriptec-test-archive").as_bytes()),
        }
    }

    fn derive_fixture(
        extracted: &Extracted,
        annotation: &anno::Annotation,
        solo_strand: anno::assign::SoloStrand,
    ) -> TranscriptEquivalenceReport {
        derive_transcript_equivalence_classes_from_extracted(
            extracted,
            annotation,
            &fixture_options(solo_strand),
            &fixture_archive_identity(),
        )
        .unwrap()
    }

    fn transcript_ids(report: &TranscriptEquivalenceReport, class: u32) -> Vec<String> {
        let class = report
            .classes
            .iter()
            .find(|entry| entry.umi_class == class)
            .unwrap();
        let Some(ec_id) = &class.ec_id else {
            return Vec::new();
        };
        report
            .transcript_sets
            .iter()
            .find(|set| &set.ec_id == ec_id)
            .unwrap()
            .transcript_ids
            .clone()
    }

    #[test]
    fn multimapper_candidates_are_unioned_across_alternative_placements() {
        let annotation = annotation(&[
            exon("chr1", 101, 200, '+', "G2", "T2"),
            exon("chr2", 301, 400, '-', "G1", "T1"),
        ]);
        let archive = extracted(
            vec![MolRec {
                cell: 0,
                umi_class: 0,
                chrom: 0,
                strand_rev: false,
                chains: SmallVec::new(),
                mms: smallvec![(110, 0, 0, 3)],
            }],
            vec![one_block(20)],
            vec![vec![
                PatAlt {
                    chrom: 0,
                    offset: 0,
                    strand_flip: false,
                    shape: SAME_SHAPE,
                },
                PatAlt {
                    chrom: 1,
                    offset: 200,
                    strand_flip: true,
                    shape: SAME_SHAPE,
                },
            ]],
            1,
            1,
            &["chr1", "chr2"],
        );

        let report = derive_fixture(&archive, &annotation, anno::assign::SoloStrand::Forward);
        assert_eq!(transcript_ids(&report, 0), ["T1", "T2"]);
        assert_eq!(report.transcript_sets[0].gene_ids, ["G1", "G2"]);
        assert!(report.classes[0].ambiguous);
        assert!(!report.classes[0].conflict);
        assert!(report.classes[0].complete_within_archive_quotient);
        assert_eq!(report.classes[0].retained_record_count, 1);
        assert_eq!(report.classes[0].represented_alignment_count, 3);
    }

    #[test]
    fn both_span_extreme_representatives_refine_the_set_and_flag_incompleteness() {
        let annotation = annotation(&[
            exon("chr1", 101, 180, '+', "G", "T1"),
            exon("chr1", 101, 130, '+', "G", "T2"),
        ]);
        let chain = |weight| MolChain {
            weight,
            reps: smallvec![(110, 0), (150, 0)],
        };
        let archive = extracted(
            vec![
                MolRec {
                    cell: 0,
                    umi_class: 0,
                    chrom: 0,
                    strand_rev: false,
                    chains: smallvec![chain(2)],
                    mms: SmallVec::new(),
                },
                MolRec {
                    cell: 0,
                    umi_class: 1,
                    chrom: 0,
                    strand_rev: false,
                    chains: smallvec![chain(3)],
                    mms: SmallVec::new(),
                },
            ],
            vec![one_block(10)],
            Vec::new(),
            2,
            1,
            &["chr1"],
        );

        let report = derive_fixture(&archive, &annotation, anno::assign::SoloStrand::Forward);
        assert_eq!(transcript_ids(&report, 0), ["T1"]);
        assert_eq!(transcript_ids(&report, 1), ["T1"]);
        assert!(report.classes[0].retained_representatives_complete);
        assert!(report.classes[0].complete_within_archive_quotient);
        assert!(!report.classes[1].retained_representatives_complete);
        assert!(!report.classes[1].complete_within_archive_quotient);
        assert_eq!(report.totals.retained_records, 4);
        assert_eq!(report.totals.represented_alignments, 5);
    }

    #[test]
    fn strand_policy_is_applied_to_transcript_candidates() {
        let annotation = annotation(&[
            exon("chr1", 101, 200, '+', "GP", "TP"),
            exon("chr1", 101, 200, '-', "GM", "TM"),
        ]);
        let archive = extracted(
            vec![unique_molecule(0, 0, 0, false, 110, 0)],
            vec![one_block(20)],
            Vec::new(),
            1,
            1,
            &["chr1"],
        );

        let forward = derive_fixture(&archive, &annotation, anno::assign::SoloStrand::Forward);
        let reverse = derive_fixture(&archive, &annotation, anno::assign::SoloStrand::Reverse);
        let unstranded =
            derive_fixture(&archive, &annotation, anno::assign::SoloStrand::Unstranded);

        assert_eq!(transcript_ids(&forward, 0), ["TP"]);
        assert_eq!(transcript_ids(&reverse, 0), ["TM"]);
        assert_eq!(transcript_ids(&unstranded, 0), ["TM", "TP"]);
        assert!(unstranded.classes[0].ambiguous);
    }

    #[test]
    fn ec_identity_is_content_addressed_and_selected_universe_filters_unrelated_classes() {
        let annotation = annotation(&[
            exon("chr1", 101, 200, '+', "G1", "T1"),
            exon("chr1", 101, 200, '+', "G2", "T2"),
        ]);
        let archive = extracted(
            vec![
                unique_molecule(0, 0, 0, false, 110, 0),
                unique_molecule(0, 1, 0, false, 400, 0),
            ],
            vec![one_block(20)],
            Vec::new(),
            2,
            1,
            &["chr1"],
        );
        let digest = *blake3::hash(b"annotation-content-a").as_bytes();
        let selected = TranscriptEquivalenceOptions::new(anno::assign::SoloStrand::Forward, digest)
            .with_transcript_universe([0]);
        let selected_report = derive_transcript_equivalence_classes_from_extracted(
            &archive,
            &annotation,
            &selected,
            &fixture_archive_identity(),
        )
        .unwrap();
        let scoped_report = derive_transcript_equivalence_classes_from_extracted(
            &archive,
            &annotation,
            &selected.clone().only_classes_with_compatible_evidence(),
            &fixture_archive_identity(),
        )
        .unwrap();
        assert_eq!(transcript_ids(&selected_report, 0), ["T1"]);
        assert_eq!(transcript_ids(&scoped_report, 0), ["T1"]);
        assert_eq!(
            selected_report.classes[0].ec_id,
            scoped_report.classes[0].ec_id
        );
        assert_eq!(scoped_report.classes.len(), 1);
        assert_eq!(scoped_report.totals.archive_classes_scanned, 2);
        assert_eq!(scoped_report.totals.classes_filtered_out, 1);
        assert_eq!(
            scoped_report.scope.selected_transcript_ids,
            Some(vec!["T1".into()])
        );
        assert_eq!(scoped_report.classes[0].ec_id.as_ref().unwrap().len(), 64);

        let different_identity = TranscriptEquivalenceOptions::new(
            anno::assign::SoloStrand::Forward,
            *blake3::hash(b"annotation-content-b").as_bytes(),
        )
        .with_transcript_universe([0])
        .only_classes_with_compatible_evidence();
        let rebound = derive_transcript_equivalence_classes_from_extracted(
            &archive,
            &annotation,
            &different_identity,
            &fixture_archive_identity(),
        )
        .unwrap();
        assert_ne!(scoped_report.classes[0].ec_id, rebound.classes[0].ec_id);
    }

    #[test]
    fn selector_overlap_is_independent_of_transcript_compatibility() {
        let annotation = annotation(&[exon("chr1", 101, 130, '+', "G1", "T1")]);
        let archive = extracted(
            vec![
                unique_molecule(0, 0, 0, false, 105, 0),
                unique_molecule(0, 1, 0, false, 120, 0),
                unique_molecule(0, 2, 0, false, 400, 0),
            ],
            vec![one_block(20)],
            Vec::new(),
            3,
            1,
            &["chr1"],
        );
        let options = fixture_options(anno::assign::SoloStrand::Forward)
            .with_transcript_universe([0])
            .with_selector_loci([TranscriptSelectorLocus {
                contig: "chr1".into(),
                start: 100,
                end: 130,
            }]);
        let report = derive_transcript_equivalence_classes_from_extracted(
            &archive,
            &annotation,
            &options,
            &fixture_archive_identity(),
        )
        .unwrap();

        assert!(report.classes[0].selector_relevant);
        assert_eq!(transcript_ids(&report, 0), ["T1"]);
        assert!(report.classes[1].selector_relevant);
        assert!(report.classes[1].no_compatible_transcript);
        assert!(report.classes[1].ec_id.is_none());
        assert!(!report.classes[2].selector_relevant);
        assert_eq!(
            report.scope.selector_loci,
            Some(vec![TranscriptSelectorLocus {
                contig: "chr1".into(),
                start: 100,
                end: 130,
            }])
        );
    }

    #[test]
    fn selector_contig_missing_from_archive_fails_closed() {
        let annotation = annotation(&[exon("chr2", 101, 130, '+', "G1", "T1")]);
        let archive = extracted(
            vec![unique_molecule(0, 0, 0, false, 105, 0)],
            vec![one_block(20)],
            Vec::new(),
            1,
            1,
            &["chr1"],
        );
        let options = fixture_options(anno::assign::SoloStrand::Forward)
            .with_transcript_universe([0])
            .with_selector_loci([TranscriptSelectorLocus {
                contig: "chr2".into(),
                start: 100,
                end: 130,
            }]);
        let error = derive_transcript_equivalence_classes_from_extracted(
            &archive,
            &annotation,
            &options,
            &fixture_archive_identity(),
        )
        .unwrap_err()
        .to_string();

        assert!(error.contains("chr2"), "{error}");
        assert!(
            error.contains("absent from the archive chromosome dictionary"),
            "{error}"
        );
    }

    #[test]
    fn empty_records_conflicts_conservation_and_output_are_deterministic() {
        let annotation = annotation(&[
            exon("chr1", 101, 130, '+', "G1", "T1"),
            exon("chr1", 201, 230, '+', "G2", "T2"),
        ]);
        let mols = vec![
            unique_molecule(0, 0, 0, false, 105, 0),
            unique_molecule(0, 0, 0, false, 400, 0),
            unique_molecule(0, 1, 0, false, 105, 0),
            unique_molecule(0, 1, 0, false, 205, 0),
            unique_molecule(1, 2, 0, false, 400, 0),
        ];
        let archive = extracted(
            mols.clone(),
            vec![one_block(10)],
            Vec::new(),
            3,
            2,
            &["chr1"],
        );
        let mut reversed = mols;
        reversed.reverse();
        let reversed_archive =
            extracted(reversed, vec![one_block(10)], Vec::new(), 3, 2, &["chr1"]);

        let report = derive_fixture(&archive, &annotation, anno::assign::SoloStrand::Forward);
        let reordered = derive_fixture(
            &reversed_archive,
            &annotation,
            anno::assign::SoloStrand::Forward,
        );
        assert_eq!(report, reordered);
        assert_eq!(
            report.to_pretty_json().unwrap(),
            reordered.to_pretty_json().unwrap()
        );

        assert!(transcript_ids(&report, 0).is_empty());
        assert_eq!(report.classes[0].unmatched_record_count, 1);
        assert!(report.classes[0].no_compatible_transcript);
        assert!(!report.classes[0].conflict);
        assert!(!report.classes[0].complete_within_archive_quotient);
        assert!(transcript_ids(&report, 1).is_empty());
        assert!(report.classes[1].conflict);
        assert!(!report.classes[1].no_compatible_transcript);
        assert!(transcript_ids(&report, 2).is_empty());
        assert!(report.classes[2].no_compatible_transcript);
        assert!(!report.classes[2].conflict);

        assert_eq!(report.totals.classes, 3);
        assert_eq!(report.totals.cells, 2);
        assert_eq!(report.totals.assigned_classes, 0);
        assert_eq!(report.totals.unassigned_classes, 3);
        assert_eq!(report.totals.no_compatible_transcript_classes, 2);
        assert_eq!(report.totals.conflicting_classes, 1);
        assert_eq!(
            report
                .transcript_sets
                .iter()
                .map(|set| set.class_count)
                .sum::<u64>(),
            report.totals.assigned_classes
        );
        assert_eq!(report.cells[0].class_count, 2);
        assert_eq!(report.cells[1].class_count, 1);
        assert!(report.transcript_sets.is_empty());
        assert!(report.classes.iter().all(|class| class.ec_id.is_none()));
        assert!(!report.semantics.abundance_inferred);
        assert!(!report.semantics.full_isoform_phasing_claimed);
    }

    fn read_u32(bytes: &[u8], at: &mut usize) -> u32 {
        let end = *at + 4;
        let value = u32::from_le_bytes(bytes[*at..end].try_into().unwrap());
        *at = end;
        value
    }

    fn skip_string(bytes: &[u8], at: &mut usize) {
        let len = read_u32(bytes, at) as usize;
        *at += len;
        assert!(*at <= bytes.len());
    }

    /// Convert a compiler-produced AIC v2 to the exact v1 layout by removing the trailing
    /// identifier metadata. This keeps the test independent of annotation crate internals.
    fn write_aic_v1(v2: &Path, v1: &Path) {
        const HEADER: usize = 8 + 4 + 8 + 32;
        let bytes = std::fs::read(v2).unwrap();
        assert_eq!(&bytes[..8], b"GRVLXAIC");
        let payload = &bytes[HEADER..];
        let mut at = 0usize;
        let genes = read_u32(payload, &mut at) as usize;
        for _ in 0..genes * 2 {
            skip_string(payload, &mut at);
        }
        let chroms = read_u32(payload, &mut at) as usize;
        for _ in 0..chroms {
            skip_string(payload, &mut at);
        }
        let transcripts = read_u32(payload, &mut at) as usize;
        for _ in 0..transcripts {
            at += 4 + 4 + 1;
            let exons = read_u32(payload, &mut at) as usize;
            at += exons * 8;
        }
        let indexed_chroms = read_u32(payload, &mut at) as usize;
        for _ in 0..indexed_chroms {
            at += 4;
            let entries = read_u32(payload, &mut at) as usize;
            at += entries * 16;
        }
        assert!(at < payload.len());
        let legacy_payload = &payload[..at];
        let mut legacy = Vec::new();
        legacy.extend_from_slice(b"GRVLXAIC");
        legacy.extend_from_slice(&1u32.to_le_bytes());
        legacy.extend_from_slice(&(legacy_payload.len() as u64).to_le_bytes());
        legacy.extend_from_slice(blake3::hash(legacy_payload).as_bytes());
        legacy.extend_from_slice(legacy_payload);
        std::fs::write(v1, legacy).unwrap();
    }

    #[test]
    fn archive_scan_intersects_globally_across_chunks_and_aic_v1_fails_clearly() {
        let records = vec![
            exon("chr1", 1, 100, '+', "G1", "T1"),
            exon("chr1", 1, 100, '+', "G2", "T2"),
            exon("chr1", 201, 300, '+', "G2", "T2"),
        ];
        let gtf = TestFile::new("gtf");
        std::fs::write(gtf.path(), records.concat()).unwrap();
        let annotation = anno::Annotation::from_gtf(gtf.path()).unwrap();
        let aic = TestFile::new("aic");
        annotation.write_compiled(aic.path()).unwrap();
        let legacy_aic = TestFile::new("v1.aic");
        write_aic_v1(aic.path(), legacy_aic.path());

        let archive_data = extracted(
            vec![
                MolRec {
                    cell: 0,
                    umi_class: 0,
                    chrom: 0,
                    strand_rev: false,
                    chains: smallvec![MolChain {
                        weight: 2,
                        reps: smallvec![(10, 0), (20, 0)],
                    }],
                    // Populate every archive value stream while keeping this record compatible
                    // with both first-exon transcripts.
                    mms: smallvec![(10, 0, 0, 1)],
                },
                MolRec {
                    cell: 0,
                    umi_class: 0,
                    chrom: 0,
                    strand_rev: false,
                    chains: smallvec![MolChain {
                        weight: 2,
                        reps: smallvec![(210, 0), (220, 0)],
                    }],
                    mms: SmallVec::new(),
                },
            ],
            vec![one_block(20)],
            vec![vec![PatAlt {
                chrom: 0,
                offset: 0,
                strand_flip: false,
                shape: SAME_SHAPE,
            }]],
            1,
            1,
            &["chr1"],
        );
        let archive = TestFile::new("aie");
        super::super::write_archive(&archive_data, &archive.0, 1, 100, None).unwrap();

        let options = TranscriptEquivalenceOptions::new(
            anno::assign::SoloStrand::Forward,
            *blake3::hash(&std::fs::read(aic.path()).unwrap()).as_bytes(),
        );
        let report = derive_transcript_equivalence_classes(
            archive.path(),
            aic.path(),
            anno::intent::AnnotationIdentity::new("test-assembly", "test-annotation").unwrap(),
            &options,
        )
        .unwrap();
        assert_eq!(transcript_ids(&report, 0), ["T2"]);
        assert_eq!(report.classes[0].retained_record_count, 5);
        assert_eq!(report.totals.classes, 1);
        assert_eq!(report.archive_version, evidence_io::format::VERSION);
        assert!(report.archive_root_blake3.is_some());

        let legacy_options = TranscriptEquivalenceOptions::new(
            anno::assign::SoloStrand::Forward,
            *blake3::hash(&std::fs::read(legacy_aic.path()).unwrap()).as_bytes(),
        );
        let error = derive_transcript_equivalence_classes(
            archive.path(),
            legacy_aic.path(),
            anno::intent::AnnotationIdentity::new("test-assembly", "legacy-annotation").unwrap(),
            &legacy_options,
        )
        .unwrap_err()
        .to_string();
        assert!(error.contains("AIC v1 omitted transcript IDs"), "{error}");
    }
}
