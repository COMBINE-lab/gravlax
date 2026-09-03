//! Public CLI adapter for paired annotation comparison.
//!
//! The scientific replay lives in `archivecmd::annotationcompare`.  This module owns only input
//! binding, result schemas, rendering, and the no-clobber output boundary.

use crate::archivecmd::annotationcompare::{
    compare_archive, AnnotationComparison, ClassState, ClassTransitionKind, CompareOptions,
    GeneKeyPolicy, TransitionCause,
};
use anno::assign::SoloStrand;
use anno::intent::{AnnotationIdentity, BoundAnnotation};
use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use gravlax_output::{
    publish_file_no_clobber, write_table, AnnotationProvenance, DataType, Durability, Field,
    OutputFormat, Producer, Provenance, ResultContext, ResultEnvelope, ScalarValue, TableSchema,
    TypedTable,
};
use serde::Serialize;
use serde_json::json;
use std::collections::BTreeMap;
use std::fmt::Write as FmtWrite;
use std::io::Write as IoWrite;
use std::path::{Path, PathBuf};

const RESULT_SCHEMA: &str = "gravlax.annotation.compare.v1";
const COUNT_DELTAS_SCHEMA: &str = "gravlax.annotation.compare.count-deltas.v1";
const CLASS_TRANSITIONS_SCHEMA: &str = "gravlax.annotation.compare.class-transitions.v1";
const CONTRIBUTING_CAUSES_SCHEMA: &str = "gravlax.annotation.compare.contributing-causes.v1";
const WITNESSES_SCHEMA: &str = "gravlax.annotation.compare.witnesses.v1";

#[derive(Debug, Parser)]
pub struct Args {
    /// Evidence archive replayed once against both annotations.
    #[arg(value_name = "ARCHIVE")]
    pub archive: PathBuf,

    /// Annotation used as the before side of the comparison (GTF or compiled AIC).
    #[arg(long, value_name = "A")]
    pub annotation_a: PathBuf,

    /// Annotation used as the after side of the comparison (GTF or compiled AIC).
    #[arg(long, value_name = "B")]
    pub annotation_b: PathBuf,

    /// Reference assembly shared by both annotations.
    #[arg(long, value_name = "ASSEMBLY")]
    pub assembly: String,

    /// Human-readable provenance label for annotation A.
    #[arg(long, value_name = "LABEL")]
    pub annotation_a_label: String,

    /// Human-readable provenance label for annotation B.
    #[arg(long, value_name = "LABEL")]
    pub annotation_b_label: String,

    /// Expected BLAKE3 content digest for annotation A (`blake3:<hex>`).
    #[arg(long, value_name = "DIGEST")]
    pub annotation_a_digest: Option<String>,

    /// Expected BLAKE3 content digest for annotation B (`blake3:<hex>`).
    #[arg(long, value_name = "DIGEST")]
    pub annotation_b_digest: Option<String>,

    /// Cross-annotation gene identifier matching policy.
    #[arg(long, value_enum, default_value_t = GeneKeyArg::Unversioned)]
    pub gene_key: GeneKeyArg,

    /// Strand convention used for both annotation assignments.
    #[arg(long, value_enum, default_value_t = SoloStrandArg::Forward)]
    pub solo_strand: SoloStrandArg,

    /// Maximum molecule witnesses retained in the whole report.
    #[arg(long, default_value_t = 10_000)]
    pub max_molecule_witnesses: usize,

    /// Maximum changed evidence rows retained in one molecule witness.
    #[arg(long, default_value_t = 32)]
    pub max_row_transitions_per_molecule: usize,

    /// Permit annotations whose observed content digests are identical.
    #[arg(long)]
    pub allow_identical: bool,

    /// Output representation. TSV requires one table selected with `--table`.
    #[arg(long, value_enum, default_value_t = Format::Text)]
    pub format: Format,

    /// Table emitted by TSV output; invalid with text or JSON output.
    #[arg(long, value_enum)]
    pub table: Option<Table>,

    /// Write atomically to a new file instead of standard output.
    #[arg(short = 'o', long, value_name = "FILE")]
    pub output: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum GeneKeyArg {
    Unversioned,
    Exact,
}

impl GeneKeyArg {
    fn policy(self) -> GeneKeyPolicy {
        match self {
            Self::Unversioned => GeneKeyPolicy::Unversioned,
            Self::Exact => GeneKeyPolicy::Exact,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Unversioned => "unversioned",
            Self::Exact => "exact",
        }
    }
}

impl std::fmt::Display for GeneKeyArg {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(self.name())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum SoloStrandArg {
    Forward,
    Reverse,
    Unstranded,
}

impl SoloStrandArg {
    fn policy(self) -> SoloStrand {
        match self {
            Self::Forward => SoloStrand::Forward,
            Self::Reverse => SoloStrand::Reverse,
            Self::Unstranded => SoloStrand::Unstranded,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
            Self::Unstranded => "unstranded",
        }
    }
}

impl std::fmt::Display for SoloStrandArg {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(self.name())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Format {
    Text,
    Json,
    Tsv,
}

impl std::fmt::Display for Format {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self {
            Self::Text => "text",
            Self::Json => "json",
            Self::Tsv => "tsv",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Table {
    CountDeltas,
    ClassTransitions,
    ContributingCauses,
    Witnesses,
}

#[derive(Debug, Serialize)]
struct CompareResultData {
    summary: CompareSummary,
    semantics: CompareSemantics,
    count_deltas: TypedTable,
    class_transitions: TypedTable,
    contributing_causes: TypedTable,
    witnesses: TypedTable,
}

#[derive(Debug, Serialize)]
struct CompareSummary {
    annotation_a_label: String,
    annotation_b_label: String,
    gene_key_policy: String,
    solo_strand: String,
    archive_version: u32,
    archive_root_blake3: Option<String>,
    archive_passes: u32,
    archive_cells: u64,
    decoded_chunks: u64,
    decoded_molecules: u64,
    evidence_rows: u64,
    annotation_a_final_gene_umis: u64,
    annotation_b_final_gene_umis: u64,
    count_delta_rows: u64,
    class_transition_rows: u64,
    changed_molecule_records: u64,
    unchanged_molecule_records: u64,
    molecule_witness_rows: u64,
    molecule_witnesses_omitted: u64,
}

#[derive(Debug, Serialize)]
struct CompareSemantics {
    final_count_deltas: &'static str,
    class_transition_ledger: &'static str,
    contributing_causes: &'static str,
    annotation_order_tie_break: &'static str,
    molecule_witnesses: &'static str,
    final_count_deltas_are_exact: bool,
    class_transition_ledger_is_complete: bool,
    contributing_causes_are_nonexclusive: bool,
    contributing_causes_are_additive_attributions: bool,
    annotation_order_tie_break_is_biological_change: bool,
    molecule_witnesses_are_bounded: bool,
}

pub fn run(args: Args) -> Result<()> {
    validate_output_contract(&args)?;
    preflight_output(args.output.as_deref())?;

    let before = bind_annotation(
        &args.annotation_a,
        &args.assembly,
        &args.annotation_a_label,
        args.annotation_a_digest.as_deref(),
    )
    .context("bind annotation A")?;
    let after = bind_annotation(
        &args.annotation_b,
        &args.assembly,
        &args.annotation_b_label,
        args.annotation_b_digest.as_deref(),
    )
    .context("bind annotation B")?;

    let before_digest = before
        .identity()
        .digest
        .clone()
        .context("bound annotation A has no observed digest")?;
    let after_digest = after
        .identity()
        .digest
        .clone()
        .context("bound annotation B has no observed digest")?;
    let identical = before_digest == after_digest;
    if identical && !args.allow_identical {
        bail!(
            "annotation A and annotation B have identical observed content digest {}; pass \
             --allow-identical to run the explicit A/A control",
            before_digest
        );
    }

    let options = CompareOptions {
        solo_strand: args.solo_strand.policy(),
        gene_key_policy: args.gene_key.policy(),
        max_molecule_witnesses: args.max_molecule_witnesses,
        max_row_transitions_per_molecule: args.max_row_transitions_per_molecule,
    };
    let comparison = compare_archive(
        &args.archive,
        before.annotation(),
        after.annotation(),
        options,
    )
    .context("compare annotations")?;

    let context = result_context(&args, &comparison, &before_digest, &after_digest, identical)?;
    let data = result_data(&args, &comparison)?;
    let bytes = render(&args, &comparison, &context, data)?;
    emit(args.output.as_deref(), &bytes)
}

fn validate_output_contract(args: &Args) -> Result<()> {
    match (args.format, args.table) {
        (Format::Tsv, None) => bail!("--format tsv requires --table"),
        (Format::Text | Format::Json, Some(_)) => {
            bail!("--table is only valid with --format tsv")
        }
        _ => Ok(()),
    }
}

fn bind_annotation(
    path: &Path,
    assembly: &str,
    label: &str,
    digest: Option<&str>,
) -> Result<BoundAnnotation> {
    let identity = AnnotationIdentity::new(assembly, label)?;
    let identity = match digest {
        Some(digest) => identity.with_digest(digest)?,
        None => identity,
    };
    Ok(BoundAnnotation::from_path(path, identity)?)
}

fn result_context(
    args: &Args,
    comparison: &AnnotationComparison,
    before_digest: &str,
    after_digest: &str,
    identical: bool,
) -> Result<ResultContext> {
    let archive = match &comparison.archive_identity.rooted_content_commitment_hex {
        Some(root) => format!("aie-directory-root-v2:{root}"),
        None => format!(
            "aie-v{}-unrooted:{}",
            comparison.archive_identity.archive_version,
            args.archive.display()
        ),
    };
    let mut parameters = BTreeMap::new();
    parameters.insert("archive_path".into(), json!(args.archive));
    parameters.insert("annotation_a_path".into(), json!(args.annotation_a));
    parameters.insert("annotation_b_path".into(), json!(args.annotation_b));
    parameters.insert("gene_key".into(), json!(args.gene_key.name()));
    parameters.insert("solo_strand".into(), json!(args.solo_strand.name()));
    parameters.insert(
        "max_molecule_witnesses".into(),
        json!(args.max_molecule_witnesses),
    );
    parameters.insert(
        "max_row_transitions_per_molecule".into(),
        json!(args.max_row_transitions_per_molecule),
    );

    let mut warnings = Vec::new();
    if comparison
        .archive_identity
        .rooted_content_commitment_hex
        .is_none()
    {
        warnings.push(
            "legacy v1 archive has no rooted content commitment; seal or rewrite it as v2 for \
             portable archive provenance"
                .into(),
        );
    }
    if identical {
        warnings.push(
            "annotations have identical observed content digests; --allow-identical enabled this \
             A/A control"
                .into(),
        );
    }

    let context = ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            archives: vec![archive],
            assembly: Some(args.assembly.clone()),
            annotation: None,
            annotation_digest: None,
            annotations: vec![
                AnnotationProvenance {
                    role: "before".into(),
                    assembly: args.assembly.clone(),
                    annotation: args.annotation_a_label.clone(),
                    digest: before_digest.into(),
                },
                AnnotationProvenance {
                    role: "after".into(),
                    assembly: args.assembly.clone(),
                    annotation: args.annotation_b_label.clone(),
                    digest: after_digest.into(),
                },
            ],
            parameters,
        },
        warnings,
    };
    context.validate()?;
    Ok(context)
}

fn result_data(args: &Args, comparison: &AnnotationComparison) -> Result<CompareResultData> {
    Ok(CompareResultData {
        summary: CompareSummary {
            annotation_a_label: args.annotation_a_label.clone(),
            annotation_b_label: args.annotation_b_label.clone(),
            gene_key_policy: args.gene_key.name().into(),
            solo_strand: args.solo_strand.name().into(),
            archive_version: comparison.archive_identity.archive_version,
            archive_root_blake3: comparison
                .archive_identity
                .rooted_content_commitment_hex
                .clone(),
            archive_passes: comparison.archive_passes,
            archive_cells: comparison.cell_barcodes.len() as u64,
            decoded_chunks: comparison.decoded_chunks,
            decoded_molecules: comparison.decoded_molecules,
            evidence_rows: comparison.evidence_rows,
            annotation_a_final_gene_umis: comparison.before.final_gene_umis,
            annotation_b_final_gene_umis: comparison.after.final_gene_umis,
            count_delta_rows: comparison.count_deltas.len() as u64,
            class_transition_rows: comparison.class_transitions.len() as u64,
            changed_molecule_records: comparison.changed_molecule_records,
            unchanged_molecule_records: comparison.unchanged_molecule_records,
            molecule_witness_rows: comparison.molecule_witnesses.len() as u64,
            molecule_witnesses_omitted: comparison.molecule_witnesses_omitted,
        },
        semantics: CompareSemantics {
            final_count_deltas: "Exact signed B-minus-A counts obtained after two independent final gene/UMI collapses.",
            class_transition_ledger: "Complete changed-class state ledger; class rows are not additive count-delta contributions.",
            contributing_causes: "Non-exclusive observed state changes, not unique or additive counterfactual attributions.",
            annotation_order_tie_break: "An annotation_order_tie_break_changed cause identifies equal maximum comparison-key support resolved differently only because annotation-local gene order changed; it preserves exact replay parity and is not a biological structural change.",
            molecule_witnesses: "Deterministically selected bounded examples; total and omitted counts remain explicit.",
            final_count_deltas_are_exact: true,
            class_transition_ledger_is_complete: true,
            contributing_causes_are_nonexclusive: true,
            contributing_causes_are_additive_attributions: false,
            annotation_order_tie_break_is_biological_change: false,
            molecule_witnesses_are_bounded: true,
        },
        count_deltas: count_delta_table(comparison)?,
        class_transitions: class_transition_table(comparison)?,
        contributing_causes: cause_table(comparison)?,
        witnesses: witness_table(comparison)?,
    })
}

fn count_delta_table(comparison: &AnnotationComparison) -> Result<TypedTable> {
    let schema = TableSchema::new(
        COUNT_DELTAS_SCHEMA,
        vec![
            Field::new("cell", DataType::UInt64),
            Field::new("cell_barcode", DataType::String),
            Field::new("comparison_gene_id", DataType::String),
            Field::new("annotation_a_gene_id", DataType::String).nullable(),
            Field::new("annotation_b_gene_id", DataType::String).nullable(),
            Field::new("annotation_a_count", DataType::UInt64),
            Field::new("annotation_b_count", DataType::UInt64),
            Field::new("signed_delta_b_minus_a", DataType::Int64),
        ],
    )?;
    let mut rows = Vec::with_capacity(comparison.count_deltas.len());
    for delta in &comparison.count_deltas {
        rows.push(vec![
            ScalarValue::UInt64(delta.cell.into()),
            cell_barcode(comparison, delta.cell)?,
            ScalarValue::String(delta.comparison_gene_id.clone()),
            optional_string(&delta.gene_id_before),
            optional_string(&delta.gene_id_after),
            ScalarValue::UInt64(delta.before.into()),
            ScalarValue::UInt64(delta.after.into()),
            ScalarValue::Int64(delta.delta),
        ]);
    }
    Ok(TypedTable::new(schema, rows)?)
}

fn class_transition_table(comparison: &AnnotationComparison) -> Result<TypedTable> {
    let schema = TableSchema::new(
        CLASS_TRANSITIONS_SCHEMA,
        vec![
            Field::new("cell", DataType::UInt64),
            Field::new("cell_barcode", DataType::String),
            Field::new("umi_class", DataType::UInt64),
            Field::new("transition_kind", DataType::String),
            Field::new("molecule_records", DataType::UInt64),
            Field::new("evidence_rows", DataType::UInt64),
            Field::new("changed_evidence_rows", DataType::UInt64),
            Field::new("annotation_a_selected_comparison_gene_id", DataType::String).nullable(),
            Field::new("annotation_a_selected_gene_id", DataType::String).nullable(),
            Field::new("annotation_a_selected_weight", DataType::UInt64),
            Field::new("annotation_a_counted", DataType::Boolean),
            Field::new("annotation_a_canonical_class", DataType::UInt64).nullable(),
            Field::new("annotation_a_gene_support", DataType::Json),
            Field::new("annotation_a_same_gene_neighbors", DataType::Json),
            Field::new("annotation_b_selected_comparison_gene_id", DataType::String).nullable(),
            Field::new("annotation_b_selected_gene_id", DataType::String).nullable(),
            Field::new("annotation_b_selected_weight", DataType::UInt64),
            Field::new("annotation_b_counted", DataType::Boolean),
            Field::new("annotation_b_canonical_class", DataType::UInt64).nullable(),
            Field::new("annotation_b_gene_support", DataType::Json),
            Field::new("annotation_b_same_gene_neighbors", DataType::Json),
            Field::new("contributing_cause_count", DataType::UInt64),
            Field::new("molecule_witnesses", DataType::UInt64),
            Field::new("omitted_molecule_witnesses", DataType::UInt64),
            Field::new("changed_row_witnesses", DataType::UInt64),
            Field::new("omitted_changed_row_witnesses", DataType::UInt64),
        ],
    )?;
    let mut rows = Vec::with_capacity(comparison.class_transitions.len());
    for transition in &comparison.class_transitions {
        let mut row = vec![
            ScalarValue::UInt64(transition.cell.into()),
            cell_barcode(comparison, transition.cell)?,
            ScalarValue::UInt64(transition.umi_class.into()),
            ScalarValue::String(transition_kind(transition.kind).into()),
            ScalarValue::UInt64(transition.evidence.molecule_records),
            ScalarValue::UInt64(transition.evidence.rows),
            ScalarValue::UInt64(transition.evidence.changed_rows),
        ];
        row.extend(class_state_cells(&transition.before)?);
        row.extend(class_state_cells(&transition.after)?);
        row.extend([
            ScalarValue::UInt64(transition.causes.len() as u64),
            ScalarValue::UInt64(transition.evidence.molecule_witnesses),
            ScalarValue::UInt64(transition.evidence.omitted_molecule_witnesses),
            ScalarValue::UInt64(transition.evidence.changed_row_witnesses),
            ScalarValue::UInt64(transition.evidence.omitted_changed_row_witnesses),
        ]);
        rows.push(row);
    }
    Ok(TypedTable::new(schema, rows)?)
}

fn class_state_cells(state: &ClassState) -> Result<Vec<ScalarValue>> {
    Ok(vec![
        optional_string(&state.selected_comparison_gene_id),
        optional_string(&state.selected_gene_id),
        ScalarValue::UInt64(state.selected_weight.into()),
        ScalarValue::Boolean(state.counted),
        optional_u32(state.canonical_class),
        json_cell(&state.gene_support)?,
        json_cell(&state.same_gene_neighbors)?,
    ])
}

fn cause_table(comparison: &AnnotationComparison) -> Result<TypedTable> {
    let schema = TableSchema::new(
        CONTRIBUTING_CAUSES_SCHEMA,
        vec![
            Field::new("cell", DataType::UInt64),
            Field::new("cell_barcode", DataType::String),
            Field::new("umi_class", DataType::UInt64),
            Field::new("transition_kind", DataType::String),
            Field::new("contributing_cause", DataType::String),
            Field::new("nonexclusive", DataType::Boolean),
            Field::new("additive_count_attribution", DataType::Boolean),
        ],
    )?;
    let mut rows = Vec::new();
    for transition in &comparison.class_transitions {
        for cause in &transition.causes {
            rows.push(vec![
                ScalarValue::UInt64(transition.cell.into()),
                cell_barcode(comparison, transition.cell)?,
                ScalarValue::UInt64(transition.umi_class.into()),
                ScalarValue::String(transition_kind(transition.kind).into()),
                ScalarValue::String(transition_cause(*cause).into()),
                ScalarValue::Boolean(true),
                ScalarValue::Boolean(false),
            ]);
        }
    }
    Ok(TypedTable::new(schema, rows)?)
}

fn witness_table(comparison: &AnnotationComparison) -> Result<TypedTable> {
    let schema = TableSchema::new(
        WITNESSES_SCHEMA,
        vec![
            Field::new("archive_ordinal", DataType::UInt64),
            Field::new("cell", DataType::UInt64),
            Field::new("cell_barcode", DataType::String),
            Field::new("umi_class", DataType::UInt64),
            Field::new("chrom", DataType::String),
            Field::new("anchor", DataType::UInt64),
            Field::new("evidence_rows", DataType::UInt64),
            Field::new("changed_rows_total", DataType::UInt64),
            Field::new("changed_rows_omitted", DataType::UInt64),
            Field::new("annotation_a_selected_comparison_gene_id", DataType::String).nullable(),
            Field::new("annotation_a_selected_gene_id", DataType::String).nullable(),
            Field::new("annotation_a_counted", DataType::Boolean),
            Field::new("annotation_a_canonical_class", DataType::UInt64).nullable(),
            Field::new("annotation_b_selected_comparison_gene_id", DataType::String).nullable(),
            Field::new("annotation_b_selected_gene_id", DataType::String).nullable(),
            Field::new("annotation_b_counted", DataType::Boolean),
            Field::new("annotation_b_canonical_class", DataType::UInt64).nullable(),
            Field::new("contributing_causes", DataType::Json),
            Field::new("changed_row_witnesses", DataType::Json),
        ],
    )?;
    let mut rows = Vec::with_capacity(comparison.molecule_witnesses.len());
    for witness in &comparison.molecule_witnesses {
        let causes = witness
            .causes
            .iter()
            .map(|cause| transition_cause(*cause))
            .collect::<Vec<_>>();
        rows.push(vec![
            ScalarValue::UInt64(witness.ordinal),
            ScalarValue::UInt64(witness.cell.into()),
            cell_barcode(comparison, witness.cell)?,
            ScalarValue::UInt64(witness.umi_class.into()),
            ScalarValue::String(witness.chrom.clone()),
            ScalarValue::UInt64(witness.anchor.into()),
            ScalarValue::UInt64(witness.rows.into()),
            ScalarValue::UInt64(witness.changed_rows_total),
            ScalarValue::UInt64(witness.changed_rows_omitted),
            optional_string(&witness.before_class.selected_comparison_gene_id),
            optional_string(&witness.before_class.selected_gene_id),
            ScalarValue::Boolean(witness.before_class.counted),
            optional_u32(witness.before_class.canonical_class),
            optional_string(&witness.after_class.selected_comparison_gene_id),
            optional_string(&witness.after_class.selected_gene_id),
            ScalarValue::Boolean(witness.after_class.counted),
            optional_u32(witness.after_class.canonical_class),
            json_cell(&causes)?,
            json_cell(&witness.changed_rows)?,
        ]);
    }
    Ok(TypedTable::new(schema, rows)?)
}

fn optional_string(value: &Option<String>) -> ScalarValue {
    value.as_ref().map_or(ScalarValue::Null, |value| {
        ScalarValue::String(value.clone())
    })
}

fn cell_barcode(comparison: &AnnotationComparison, cell: u32) -> Result<ScalarValue> {
    let barcode = comparison
        .cell_barcodes
        .get(cell as usize)
        .with_context(|| format!("comparison references absent cell id {cell}"))?;
    Ok(ScalarValue::String(barcode.clone()))
}

fn optional_u32(value: Option<u32>) -> ScalarValue {
    value.map_or(ScalarValue::Null, |value| ScalarValue::UInt64(value.into()))
}

fn json_cell<T: Serialize>(value: &T) -> Result<ScalarValue> {
    let value = serde_json::to_value(value)?;
    debug_assert!(!value.is_null());
    Ok(ScalarValue::Json(value))
}

fn transition_kind(kind: ClassTransitionKind) -> &'static str {
    match kind {
        ClassTransitionKind::GainedFinalCount => "gained_final_count",
        ClassTransitionKind::LostFinalCount => "lost_final_count",
        ClassTransitionKind::ReassignedFinalCount => "reassigned_final_count",
        ClassTransitionKind::ChangedWithoutFinalCountDelta => "changed_without_final_count_delta",
    }
}

fn transition_cause(cause: TransitionCause) -> &'static str {
    match cause {
        TransitionCause::CandidateSetChanged => "candidate_set_changed",
        TransitionCause::RowAssignmentChanged => "row_assignment_changed",
        TransitionCause::ClassSupportChanged => "class_support_changed",
        TransitionCause::AnnotationOrderTieBreakChanged => "annotation_order_tie_break_changed",
        TransitionCause::ClassWinnerChanged => "class_winner_changed",
        TransitionCause::CollapseNeighborhoodChanged => "collapse_neighborhood_changed",
        TransitionCause::CollapseOutcomeChanged => "collapse_outcome_changed",
        TransitionCause::FinalContributionChanged => "final_contribution_changed",
    }
}

fn render(
    args: &Args,
    comparison: &AnnotationComparison,
    context: &ResultContext,
    data: CompareResultData,
) -> Result<Vec<u8>> {
    match args.format {
        Format::Text => Ok(render_text(args, comparison, context).into_bytes()),
        Format::Json => {
            let envelope = ResultEnvelope::new(RESULT_SCHEMA, context.clone(), data)?;
            let mut bytes = serde_json::to_vec(&envelope)?;
            bytes.push(b'\n');
            Ok(bytes)
        }
        Format::Tsv => {
            let table = selected_table(&data, args.table.expect("validated TSV table"));
            let mut bytes = Vec::new();
            write_table(
                &mut bytes,
                &table.schema,
                table.rows.clone(),
                OutputFormat::Tsv,
                context,
            )?;
            Ok(bytes)
        }
    }
}

fn selected_table(data: &CompareResultData, table: Table) -> &TypedTable {
    match table {
        Table::CountDeltas => &data.count_deltas,
        Table::ClassTransitions => &data.class_transitions,
        Table::ContributingCauses => &data.contributing_causes,
        Table::Witnesses => &data.witnesses,
    }
}

fn render_text(args: &Args, comparison: &AnnotationComparison, context: &ResultContext) -> String {
    let positive: i64 = comparison
        .count_deltas
        .iter()
        .map(|delta| delta.delta.max(0))
        .sum();
    let negative: i64 = comparison
        .count_deltas
        .iter()
        .map(|delta| (-delta.delta).max(0))
        .sum();
    let mut text = String::new();
    let _ = writeln!(
        text,
        "Annotation comparison: {} (A/before) -> {} (B/after)",
        args.annotation_a_label, args.annotation_b_label
    );
    match &comparison.archive_identity.rooted_content_commitment_hex {
        Some(root) => {
            let _ = writeln!(
                text,
                "Archive: v{} root {}",
                comparison.archive_identity.archive_version, root
            );
        }
        None => {
            let _ = writeln!(
                text,
                "Archive: legacy v{} without a rooted content commitment",
                comparison.archive_identity.archive_version
            );
        }
    }
    let _ = writeln!(
        text,
        "Final gene/UMIs: A={} B={} (+{} / -{} across {} changed cell-gene rows)",
        comparison.before.final_gene_umis,
        comparison.after.final_gene_umis,
        positive,
        negative,
        comparison.count_deltas.len()
    );
    let _ = writeln!(
        text,
        "Changed UMI classes: {}; molecule witnesses: {} retained, {} omitted",
        comparison.class_transitions.len(),
        comparison.molecule_witnesses.len(),
        comparison.molecule_witnesses_omitted
    );
    let _ = writeln!(
        text,
        "Semantics: final count deltas are exact after two independent final collapses."
    );
    let _ = writeln!(
        text,
        "Causes are non-exclusive observed state changes, not additive count-delta attributions; witnesses are bounded."
    );
    let _ = writeln!(
        text,
        "Annotation-order tie-break causes are replay-method artifacts, not biological structural changes."
    );
    let _ = writeln!(
        text,
        "Policy: gene-key={} solo-strand={} archive-passes={}",
        args.gene_key, args.solo_strand, comparison.archive_passes
    );
    if comparison.count_deltas.is_empty() {
        let _ = writeln!(text, "No final gene/UMI count differences.");
    } else {
        let _ = writeln!(text, "Largest signed B-minus-A changes:");
        let mut deltas = comparison.count_deltas.iter().collect::<Vec<_>>();
        deltas.sort_by_key(|delta| {
            (
                std::cmp::Reverse(delta.delta.unsigned_abs()),
                delta.cell,
                &delta.comparison_gene_id,
            )
        });
        for delta in deltas.iter().take(10) {
            let barcode = comparison
                .cell_barcodes
                .get(delta.cell as usize)
                .map(String::as_str)
                .unwrap_or("<invalid-cell>");
            let _ = writeln!(
                text,
                "  cell {} ({}) {}: {} -> {} ({:+})",
                barcode,
                delta.cell,
                delta.comparison_gene_id,
                delta.before,
                delta.after,
                delta.delta
            );
        }
        if deltas.len() > 10 {
            let _ = writeln!(
                text,
                "  ... {} more rows (use --format json or tsv)",
                deltas.len() - 10
            );
        }
    }
    for warning in &context.warnings {
        let _ = writeln!(text, "Warning: {warning}");
    }
    text
}

fn preflight_output(output: Option<&Path>) -> Result<()> {
    let Some(output) = output else {
        return Ok(());
    };
    if output.exists() {
        bail!("refusing to overwrite existing output {}", output.display());
    }
    let parent = output_parent(output);
    if !parent.is_dir() {
        bail!("output parent {} is not a directory", parent.display());
    }
    Ok(())
}

fn emit(output: Option<&Path>, bytes: &[u8]) -> Result<()> {
    match output {
        Some(path) => {
            let outcome = publish_file_no_clobber(path, Durability::File, |writer| {
                writer.write_all(bytes)?;
                Ok(())
            })
            .with_context(|| {
                format!(
                    "publish comparison without overwriting existing output {}",
                    path.display()
                )
            })?;
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
            Ok(())
        }
        None => {
            let stdout = std::io::stdout();
            let mut lock = stdout.lock();
            lock.write_all(bytes)
                .context("write comparison to stdout")?;
            lock.flush().context("flush comparison stdout")
        }
    }
}

fn output_parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotation_order_tie_break_has_a_stable_public_cause_name() {
        assert_eq!(
            transition_cause(TransitionCause::AnnotationOrderTieBreakChanged),
            "annotation_order_tie_break_changed"
        );
    }
}
