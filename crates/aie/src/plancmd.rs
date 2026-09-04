//! Versioned, declarative analysis plans.
//!
//! Plans resolve named project resources into absolute inputs and project-contained outputs. The
//! resolved command list is persisted before execution, making the exact invocation inspectable
//! without introducing a second implementation of any analysis operation.

use crate::projectcmd::{
    ensure_project_directory, find_project, find_project_from, resolve_project_output,
    validate_resource_id, ProjectAnnotationIdentity, ProjectContext, ResourceKind,
};
use anno::intent::{
    AnnotationIdentity, FeatureKind, IdentifierQuery, IntentResolver, MatchBasis, Strand,
};
use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use gravlax_output::{install_open_file_no_clobber, Durability, LOGICAL_OUTPUT_MAP_ENV};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command as ProcessCommand, Stdio};

pub const PLAN_SCHEMA_VERSION: u32 = 1;
pub const RESOLVED_PLAN_SCHEMA_VERSION: u32 = 6;
const SUPPORTED_KINDS: &[&str] = &[
    "inspect-archive",
    "query-region",
    "query-junction",
    "query-junctions",
    "query-jset",
    "query-events",
    "query-apa",
    "compare-annotations",
    "query-transcript-ecs",
    "federate-junction",
    "collection-region",
    "collection-junction",
    "collection-jset",
    "cohort-events",
    "cohort-splice-graph",
    "compile-annotation",
    "ingest-archive",
    "replay-rows",
    "extend-annotation",
];

#[derive(Parser)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Validate and resolve a YAML or JSON plan without running it.
    Check(CheckArgs),
    /// Resolve, persist, and execute a plan through existing `aie` commands.
    Run(RunArgs),
}

#[derive(clap::Args)]
struct CheckArgs {
    /// Versioned YAML or JSON plan.
    plan: PathBuf,
    /// Project directory or manifest; otherwise search upward from the plan.
    #[arg(long)]
    project: Option<PathBuf>,
    /// Explain every resource, output, and command resolution.
    #[arg(long)]
    explain: bool,
    /// Emit the fully resolved plan as JSON.
    #[arg(long)]
    json: bool,
}

#[derive(clap::Args)]
struct RunArgs {
    /// Versioned YAML or JSON plan.
    plan: PathBuf,
    /// Project directory or manifest; otherwise search upward from the plan.
    #[arg(long)]
    project: Option<PathBuf>,
    /// Explain every resource, output, and command resolution before execution.
    #[arg(long)]
    explain: bool,
    /// Resolve and display the plan without writing or executing anything.
    #[arg(long)]
    dry_run: bool,
    /// Skip a completed step only when its exact plan digest and every output identity match the
    /// persisted completion record. Existing unverified outputs are an error.
    #[arg(long)]
    resume: bool,
}

#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum OutputFormat {
    #[default]
    Human,
    Tsv,
    Json,
}

/// Explicit opt-in to the uniform result/report contract.
///
/// This is deliberately separate from the older per-command `format` and `output` fields.  A
/// source-plan v1 document that omits this object therefore compiles to the same legacy command
/// line as before, while a new plan can request the normalized streaming contract without
/// overloading an established field.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum UniformFormat {
    Text,
    Tsv,
    Json,
}

impl UniformFormat {
    fn as_arg(self) -> &'static str {
        match self {
            Self::Text => "text",
            Self::Tsv => "tsv",
            Self::Json => "json",
        }
    }
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct UniformOutput {
    pub format: UniformFormat,
    /// Project-contained destination.  Omit it to stream the schema-bearing result to stdout.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub output: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum SoloStrand {
    #[default]
    Forward,
    Reverse,
    Unstranded,
}

#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ComparisonGeneKey {
    #[default]
    Unversioned,
    Exact,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum AnnotationComparisonTable {
    CountDeltas,
    ClassTransitions,
    ContributingCauses,
    Witnesses,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum TranscriptEcTable {
    Catalog,
    Counts,
    Membership,
}

#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum QueryAggregation {
    #[default]
    Auto,
    Cell,
    Group,
    Bulk,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct QueryScope {
    /// Named `cells` project resource.
    #[serde(default)]
    pub cells: Option<String>,
    /// Named `groups` project resource.
    #[serde(default)]
    pub groups: Option<String>,
    #[serde(default)]
    pub aggregation: QueryAggregation,
}

/// A biological identifier may use a compact scalar form (`feature: TP53`) or carry explicit
/// expectations that are checked against immutable project metadata. The qualified form prevents
/// an old plan from silently being interpreted against a different assembly or annotation label.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(untagged)]
pub enum FeatureRequest {
    Identifier(String),
    Qualified(QualifiedFeatureRequest),
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct QualifiedFeatureRequest {
    pub identifier: String,
    #[serde(default)]
    pub assembly: Option<String>,
    #[serde(default)]
    pub annotation: Option<String>,
}

impl FeatureRequest {
    fn identifier(&self) -> &str {
        match self {
            Self::Identifier(identifier) => identifier,
            Self::Qualified(request) => &request.identifier,
        }
    }

    fn expected_assembly(&self) -> Option<&str> {
        match self {
            Self::Identifier(_) => None,
            Self::Qualified(request) => request.assembly.as_deref(),
        }
    }

    fn expected_annotation(&self) -> Option<&str> {
        match self {
            Self::Identifier(_) => None,
            Self::Qualified(request) => request.annotation.as_deref(),
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum EventType {
    AltAcceptor,
    AltDonor,
    Cassette,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ApaStrand {
    Forward,
    Reverse,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AnalysisPlan {
    pub schema_version: u32,
    #[serde(default)]
    pub name: Option<String>,
    pub steps: Vec<PlanStep>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
pub enum PlanStep {
    InspectArchive {
        id: String,
        archive: String,
        #[serde(default)]
        verify_content: bool,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
    },
    QueryRegion {
        id: String,
        archive: String,
        #[serde(default)]
        locus: Option<String>,
        #[serde(default)]
        feature: Option<FeatureRequest>,
        #[serde(default = "default_query_top")]
        top: usize,
        #[serde(default)]
        annotation: Option<String>,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
        #[serde(default)]
        scope: QueryScope,
    },
    QueryJunction {
        id: String,
        archive: String,
        locus: String,
        #[serde(default = "default_query_top")]
        top: usize,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
        #[serde(default)]
        scope: QueryScope,
    },
    QueryJunctions {
        id: String,
        archive: String,
        #[serde(default)]
        locus: Option<String>,
        #[serde(default)]
        feature: Option<FeatureRequest>,
        #[serde(default)]
        either: bool,
        #[serde(default = "default_min_support")]
        min_support: u64,
        #[serde(default)]
        with_cells: bool,
        #[serde(default)]
        min_cells: usize,
        #[serde(default)]
        annotation: Option<String>,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
        #[serde(default)]
        scope: QueryScope,
    },
    QueryJset {
        id: String,
        archive: String,
        include: Vec<String>,
        exclude: Vec<String>,
        #[serde(default = "default_query_top")]
        top: usize,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
        #[serde(default)]
        scope: QueryScope,
    },
    QueryEvents {
        id: String,
        archive: String,
        #[serde(default)]
        locus: Option<String>,
        #[serde(default)]
        feature: Option<FeatureRequest>,
        #[serde(default)]
        event_types: Vec<EventType>,
        #[serde(default = "default_events_min_support")]
        min_support: u64,
        #[serde(default = "default_min_informative")]
        min_informative: usize,
        #[serde(default = "default_max_events")]
        max_events: usize,
        #[serde(default = "default_query_top")]
        top: usize,
        #[serde(default)]
        annotation: Option<String>,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
        #[serde(default)]
        scope: QueryScope,
    },
    QueryApa {
        id: String,
        archive: String,
        #[serde(default)]
        locus: Option<String>,
        #[serde(default)]
        feature: Option<FeatureRequest>,
        #[serde(default)]
        annotation: Option<String>,
        #[serde(default = "default_site_gap")]
        site_gap: u32,
        #[serde(default)]
        strand: Option<ApaStrand>,
        #[serde(default)]
        tsv: bool,
        #[serde(default)]
        groups: Option<String>,
        #[serde(default)]
        genome: Option<String>,
        #[serde(default)]
        drop_ip: bool,
        #[serde(default)]
        permute: usize,
        #[serde(default = "default_seed")]
        seed: u64,
        #[serde(default)]
        plot: Option<PathBuf>,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
    },
    CompareAnnotations {
        id: String,
        archive: String,
        annotation_a: String,
        annotation_b: String,
        #[serde(default)]
        gene_key: ComparisonGeneKey,
        #[serde(default)]
        solo_strand: SoloStrand,
        #[serde(default = "default_max_molecule_witnesses")]
        max_molecule_witnesses: usize,
        #[serde(default = "default_max_row_transitions_per_molecule")]
        max_row_transitions_per_molecule: usize,
        #[serde(default)]
        allow_identical: bool,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        table: Option<AnnotationComparisonTable>,
        #[serde(default)]
        output: Option<PathBuf>,
    },
    QueryTranscriptEcs {
        id: String,
        archive: String,
        #[serde(default)]
        locus: Option<String>,
        #[serde(default)]
        feature: Option<FeatureRequest>,
        annotation: String,
        #[serde(default)]
        solo_strand: SoloStrand,
        #[serde(default)]
        scope: QueryScope,
        #[serde(default)]
        emit_membership: bool,
        #[serde(default = "default_max_ecs")]
        max_ecs: usize,
        #[serde(default = "default_max_memberships")]
        max_memberships: usize,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        table: Option<TranscriptEcTable>,
        #[serde(default)]
        output: Option<PathBuf>,
    },
    FederateJunction {
        id: String,
        archives: Vec<String>,
        locus: String,
        #[serde(default = "default_federate_top")]
        top: usize,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
    },
    CollectionRegion {
        id: String,
        collection: String,
        #[serde(default)]
        locus: Option<String>,
        #[serde(default)]
        feature: Option<FeatureRequest>,
        #[serde(default)]
        annotation: Option<String>,
        #[serde(default = "default_federate_top")]
        top: usize,
        #[serde(default)]
        explain_routing: bool,
        #[serde(default)]
        verify_content: bool,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
    },
    CollectionJunction {
        id: String,
        collection: String,
        locus: String,
        #[serde(default)]
        min_support: u64,
        #[serde(default = "default_federate_top")]
        top: usize,
        #[serde(default)]
        explain_routing: bool,
        #[serde(default)]
        verify_content: bool,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
    },
    CollectionJset {
        id: String,
        collection: String,
        include: Vec<String>,
        exclude: Vec<String>,
        #[serde(default)]
        min_support: u64,
        #[serde(default = "default_federate_top")]
        top: usize,
        #[serde(default)]
        explain_routing: bool,
        #[serde(default)]
        verify_content: bool,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
    },
    CohortEvents {
        id: String,
        /// Cohort sample ID -> named archive resource.
        samples: BTreeMap<String, String>,
        /// Cohort sample ID -> named groups resource.
        #[serde(default)]
        groups: BTreeMap<String, String>,
        #[serde(default)]
        locus: Option<String>,
        #[serde(default)]
        feature: Option<FeatureRequest>,
        #[serde(default)]
        event_types: Vec<EventType>,
        #[serde(default = "default_events_min_support")]
        min_support: u64,
        #[serde(default = "default_min_samples")]
        min_samples: usize,
        #[serde(default = "default_min_informative")]
        min_informative: usize,
        #[serde(default)]
        min_row_informative: usize,
        #[serde(default = "default_max_events")]
        max_events: usize,
        #[serde(default)]
        annotation: Option<String>,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
    },
    CohortSpliceGraph {
        id: String,
        #[serde(default)]
        locus: Option<String>,
        #[serde(default)]
        feature: Option<FeatureRequest>,
        /// Annotation used only to resolve `feature`; the delegated splice-graph command consumes
        /// the resulting explicit locus.
        #[serde(default)]
        annotation: Option<String>,
        /// Named `design` project resource.
        design: String,
        #[serde(default)]
        contrast: Option<String>,
        #[serde(default)]
        counts_only: bool,
        #[serde(default = "default_min_support")]
        min_support: u64,
        #[serde(default = "default_min_samples")]
        min_edge_samples: usize,
        #[serde(default = "default_min_sample_umis")]
        min_sample_umis: usize,
        #[serde(default = "default_min_replicates")]
        min_replicates: usize,
        #[serde(default = "default_min_umis")]
        min_path_umis: usize,
        #[serde(default = "default_min_path_samples")]
        min_path_samples: usize,
        #[serde(default = "default_max_events")]
        max_paths: usize,
        #[serde(default)]
        format: OutputFormat,
        #[serde(default)]
        output: Option<PathBuf>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_output: Option<UniformOutput>,
    },
    CompileAnnotation {
        id: String,
        annotation: String,
        output: PathBuf,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_report: Option<UniformOutput>,
    },
    IngestArchive {
        id: String,
        bam: String,
        whitelist: String,
        output: PathBuf,
        #[serde(default)]
        genome: Option<String>,
        #[serde(default = "default_locus_gap")]
        locus_gap: u32,
        #[serde(default = "default_zstd_level")]
        zstd_level: i32,
        #[serde(default = "default_chunk_mb")]
        chunk_mb: u32,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_report: Option<UniformOutput>,
    },
    ReplayRows {
        id: String,
        archive: String,
        annotation: String,
        barcodes: String,
        out_dir: PathBuf,
        #[serde(default = "default_locus_gap")]
        locus_gap: u32,
        #[serde(default)]
        velocity: bool,
        #[serde(default)]
        solo_strand: SoloStrand,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_report: Option<UniformOutput>,
    },
    ExtendAnnotation {
        id: String,
        archive: String,
        annotation: String,
        out_gtf: PathBuf,
        #[serde(default)]
        report: Option<PathBuf>,
        #[serde(default)]
        genome: Option<String>,
        #[serde(default = "default_max_extend")]
        max_extend: u32,
        #[serde(default = "default_evidence_gap")]
        evidence_gap: u32,
        #[serde(default = "default_min_umis")]
        min_umis: usize,
        #[serde(default = "default_min_cells")]
        min_cells: usize,
        #[serde(default = "default_site_gap")]
        site_gap: u32,
        #[serde(default = "default_min_extension")]
        min_extension: u32,
        #[serde(default)]
        clip_any_strand: bool,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        uniform_report: Option<UniformOutput>,
    },
}

fn default_query_top() -> usize {
    20
}
fn default_max_molecule_witnesses() -> usize {
    10_000
}
fn default_max_row_transitions_per_molecule() -> usize {
    32
}
fn default_max_ecs() -> usize {
    100_000
}
fn default_max_memberships() -> usize {
    1_000_000
}
fn default_federate_top() -> usize {
    5
}
fn default_min_support() -> u64 {
    1
}
fn default_events_min_support() -> u64 {
    2
}
fn default_min_informative() -> usize {
    1
}
fn default_min_samples() -> usize {
    2
}
fn default_max_events() -> usize {
    100_000
}
fn default_locus_gap() -> u32 {
    2_000
}
fn default_zstd_level() -> i32 {
    19
}
fn default_chunk_mb() -> u32 {
    4
}
fn default_max_extend() -> u32 {
    10_000
}
fn default_evidence_gap() -> u32 {
    2_000
}
fn default_min_umis() -> usize {
    5
}
fn default_min_cells() -> usize {
    3
}
fn default_site_gap() -> u32 {
    24
}
fn default_min_extension() -> u32 {
    50
}
fn default_seed() -> u64 {
    1
}
fn default_min_sample_umis() -> usize {
    10
}
fn default_min_replicates() -> usize {
    2
}
fn default_min_path_samples() -> usize {
    2
}

impl PlanStep {
    pub fn id(&self) -> &str {
        match self {
            Self::InspectArchive { id, .. }
            | Self::QueryRegion { id, .. }
            | Self::QueryJunction { id, .. }
            | Self::QueryJunctions { id, .. }
            | Self::QueryJset { id, .. }
            | Self::QueryEvents { id, .. }
            | Self::QueryApa { id, .. }
            | Self::CompareAnnotations { id, .. }
            | Self::QueryTranscriptEcs { id, .. }
            | Self::FederateJunction { id, .. }
            | Self::CollectionRegion { id, .. }
            | Self::CollectionJunction { id, .. }
            | Self::CollectionJset { id, .. }
            | Self::CohortEvents { id, .. }
            | Self::CohortSpliceGraph { id, .. }
            | Self::CompileAnnotation { id, .. }
            | Self::IngestArchive { id, .. }
            | Self::ReplayRows { id, .. }
            | Self::ExtendAnnotation { id, .. } => id,
        }
    }

    pub fn kind(&self) -> &'static str {
        match self {
            Self::InspectArchive { .. } => "inspect-archive",
            Self::QueryRegion { .. } => "query-region",
            Self::QueryJunction { .. } => "query-junction",
            Self::QueryJunctions { .. } => "query-junctions",
            Self::QueryJset { .. } => "query-jset",
            Self::QueryEvents { .. } => "query-events",
            Self::QueryApa { .. } => "query-apa",
            Self::CompareAnnotations { .. } => "compare-annotations",
            Self::QueryTranscriptEcs { .. } => "query-transcript-ecs",
            Self::FederateJunction { .. } => "federate-junction",
            Self::CollectionRegion { .. } => "collection-region",
            Self::CollectionJunction { .. } => "collection-junction",
            Self::CollectionJset { .. } => "collection-jset",
            Self::CohortEvents { .. } => "cohort-events",
            Self::CohortSpliceGraph { .. } => "cohort-splice-graph",
            Self::CompileAnnotation { .. } => "compile-annotation",
            Self::IngestArchive { .. } => "ingest-archive",
            Self::ReplayRows { .. } => "replay-rows",
            Self::ExtendAnnotation { .. } => "extend-annotation",
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ResolvedPlan {
    pub schema_version: u32,
    pub plan_schema_version: u32,
    pub producer: ResolvedProducer,
    pub name: String,
    pub source_path: PathBuf,
    pub source_digest: String,
    pub project_name: String,
    pub project_root: PathBuf,
    pub project_manifest: PathBuf,
    pub project_manifest_digest: String,
    pub resources: BTreeMap<String, ResolvedResource>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub embedded_resources: BTreeMap<String, ResolvedEmbeddedResource>,
    pub steps: Vec<ResolvedStep>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedProducer {
    pub name: String,
    pub version: String,
    /// Versioned identity of the resolver/execution contract, independent of file schemas.
    pub plan_engine: String,
    /// Exact bytes of the executable that compiled and runs the resolved command list.
    pub executable_identity: ResolvedIdentity,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ResolvedResource {
    pub kind: ResourceKind,
    pub path: PathBuf,
    pub external: bool,
    pub bytes: u64,
    pub identity: ResolvedIdentity,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_identity: Option<ProjectAnnotationIdentity>,
}

/// A path embedded inside a registered, human-editable resource such as a cohort design.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedEmbeddedResource {
    pub owner_resource: String,
    pub sample: String,
    pub role: String,
    pub declared_path: String,
    pub kind: ResourceKind,
    pub path: PathBuf,
    pub external: bool,
    pub bytes: u64,
    pub identity: ResolvedIdentity,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct ResolvedIdentity {
    pub scheme: String,
    pub digest: String,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ResolvedOutputKind {
    File,
    Directory,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ResolvedOutput {
    /// Stable name used by downstream `step:<id>:<name>` references.
    pub name: String,
    pub path: PathBuf,
    pub kind: ResolvedOutputKind,
    /// Semantic input role this artifact may satisfy in a later step.
    pub resource_kind: ResourceKind,
    /// Deterministic project-private target passed to the child command before no-clobber install.
    pub staging_path: PathBuf,
    /// Present for annotation outputs whose biological semantics and identity metadata are
    /// inherited from a registered source annotation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_semantics: Option<ResolvedAnnotationSemantics>,
    /// Declared coordinate assembly inherited by a derived archive when it can be established
    /// from a registered genome input.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedAnnotationSemantics {
    pub assembly: String,
    pub annotation: String,
    /// Existing annotation snapshot used to resolve identifiers before the derived annotation is
    /// produced. The identity binds this exact source content.
    pub source_path: PathBuf,
    pub source_identity: ResolvedIdentity,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedStepInput {
    pub reference: String,
    pub producer_step: String,
    pub output_name: String,
    pub path: PathBuf,
    pub kind: ResolvedOutputKind,
    pub resource_kind: ResourceKind,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedPreparedInput {
    pub path: PathBuf,
    pub bytes: u64,
    pub identity: ResolvedIdentity,
    pub content: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ResolvedStep {
    pub id: String,
    pub kind: String,
    /// Arguments after the `aie` executable. No shell is involved in execution.
    pub args: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stdout: Option<PathBuf>,
    /// Present only when the source plan explicitly selected the uniform result/report contract.
    /// `stdout: null` alone means inherited child stdout; this record makes the machine-readable
    /// format and publication boundary unambiguous.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uniform_io: Option<ResolvedUniformIo>,
    pub outputs: Vec<ResolvedOutput>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub input_resources: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub step_inputs: Vec<ResolvedStepInput>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub embedded_resources: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub prepared_inputs: Vec<ResolvedPreparedInput>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub biological_intent: Option<ResolvedBiologicalIntent>,
    /// Annotation inputs with explicit scientific roles. For a prior-step compiled annotation,
    /// `expected_command_identity` is unavailable until that producer runs; its source identity
    /// and the producer completion identity still bind the dataflow without fabricating a digest.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub annotation_inputs: Vec<ResolvedAnnotationInput>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_comparison: Option<ResolvedAnnotationComparisonIntent>,
    pub output_schema_ids: Vec<String>,
    pub io_estimate: ResolvedIoEstimate,
    pub explanation: Vec<String>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ResolvedUniformIoKind {
    Result,
    Report,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ResolvedUniformPublication {
    Stdout,
    AtomicNoClobberFile,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedUniformIo {
    pub kind: ResolvedUniformIoKind,
    pub format: UniformFormat,
    pub publication: ResolvedUniformPublication,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ResolvedAnnotationRole {
    A,
    B,
    Query,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedAnnotationInput {
    pub role: ResolvedAnnotationRole,
    pub resource: String,
    pub annotation_path: PathBuf,
    pub source_path: PathBuf,
    pub assembly: String,
    pub annotation: String,
    pub source_identity: ResolvedIdentity,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected_command_identity: Option<ResolvedIdentity>,
    pub compatibility: Vec<ResolvedAssemblyCompatibility>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedAnnotationComparisonIntent {
    pub annotation_a_resource: String,
    pub annotation_b_resource: String,
    pub assembly: String,
    pub gene_key: ComparisonGeneKey,
    pub solo_strand: SoloStrand,
    pub final_count_delta_semantics: String,
    pub transition_evidence_semantics: String,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedBiologicalIntent {
    pub requested: String,
    pub resolved_kind: FeatureKind,
    pub stable_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_name: Option<String>,
    pub matched_by: MatchBasis,
    pub gene_ids: Vec<String>,
    pub transcript_ids: Vec<String>,
    pub annotation_resource: String,
    pub annotation_path: PathBuf,
    pub assembly: String,
    pub annotation: String,
    /// Digest of the exact annotation snapshot used for identifier resolution.
    pub annotation_digest: String,
    pub contig: String,
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    /// Explicit 0-based, half-open locus delegated to the scientific command.
    pub locus: String,
    pub compatibility: Vec<ResolvedAssemblyCompatibility>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum AssemblyCompatibilityStatus {
    Verified,
    Unverified,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedAssemblyCompatibility {
    pub resource: String,
    pub kind: ResourceKind,
    pub status: AssemblyCompatibilityStatus,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub declared_assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chromosome_digest: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genome_digest: Option<String>,
    pub note: String,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum IoEstimateBound {
    WholeSelectedFiles,
    KnownInputsOnly,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResolvedIoEstimate {
    /// Exact size of known selected input artifacts at check time, deduplicated by canonical path.
    pub known_selected_input_bytes: u64,
    pub known_selected_input_files: usize,
    pub unknown_prior_step_outputs: usize,
    /// A deliberately conservative bound. Zero is not a prediction; it expresses that routing,
    /// caching, and operating-system reads are not modeled by the plan checker.
    pub read_bytes_lower_bound: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_bytes_upper_bound: Option<u64>,
    pub bound: IoEstimateBound,
    pub note: String,
}

pub const COMPLETION_SCHEMA_VERSION: u32 = 2;

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StepCompletion {
    pub schema_version: u32,
    pub resolved_plan_digest: String,
    pub step_id: String,
    pub step_digest: String,
    pub inputs: Vec<CompletedStepInput>,
    pub outputs: Vec<CompletedOutput>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct CompletedStepInput {
    pub producer_step: String,
    pub output_name: String,
    pub path: PathBuf,
    pub kind: ResolvedOutputKind,
    pub bytes: u64,
    pub identity: ResolvedIdentity,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct CompletedOutput {
    pub path: PathBuf,
    pub kind: ResolvedOutputKind,
    pub bytes: u64,
    pub identity: ResolvedIdentity,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ResumeDecision {
    Run,
    SkipVerified,
}

struct Resolver<'a> {
    project: &'a ProjectContext,
    all_step_ids: BTreeSet<String>,
    prior_outputs: BTreeMap<String, Vec<ResolvedOutput>>,
    current_step: String,
    resources: BTreeMap<String, ResolvedResource>,
    embedded_resources: BTreeMap<String, ResolvedEmbeddedResource>,
    identity_cache: BTreeMap<(ResourceKind, PathBuf), (u64, ResolvedIdentity)>,
    step_resources: BTreeSet<String>,
    step_inputs: BTreeMap<String, ResolvedStepInput>,
    step_embedded_resources: BTreeSet<String>,
}

#[derive(Clone)]
struct AnnotationResolutionInput {
    command_path: PathBuf,
    source_path: PathBuf,
    metadata: ProjectAnnotationIdentity,
    source_identity: ResolvedIdentity,
    command_identity: Option<ResolvedIdentity>,
}

impl<'a> Resolver<'a> {
    fn new(project: &'a ProjectContext, all_step_ids: BTreeSet<String>) -> Self {
        Self {
            project,
            all_step_ids,
            prior_outputs: BTreeMap::new(),
            current_step: String::new(),
            resources: BTreeMap::new(),
            embedded_resources: BTreeMap::new(),
            identity_cache: BTreeMap::new(),
            step_resources: BTreeSet::new(),
            step_inputs: BTreeMap::new(),
            step_embedded_resources: BTreeSet::new(),
        }
    }

    fn begin_step(&mut self, step_id: &str) {
        self.current_step = step_id.to_owned();
        self.step_resources.clear();
        self.step_inputs.clear();
        self.step_embedded_resources.clear();
    }

    fn finish_step(&mut self) -> (Vec<String>, Vec<ResolvedStepInput>, Vec<String>) {
        (
            std::mem::take(&mut self.step_resources)
                .into_iter()
                .collect(),
            std::mem::take(&mut self.step_inputs)
                .into_values()
                .collect(),
            std::mem::take(&mut self.step_embedded_resources)
                .into_iter()
                .collect(),
        )
    }

    fn register_outputs(&mut self, step: &ResolvedStep) {
        self.prior_outputs
            .insert(step.id.clone(), step.outputs.clone());
    }

    fn resource(&mut self, name: &str, accepted: &[ResourceKind]) -> Result<PathBuf> {
        if name.starts_with("step:") {
            return self.step_output_resource(name, accepted);
        }
        validate_resource_id(name)?;
        self.step_resources.insert(name.to_owned());
        let (kind, path) = self.project.resolve_resource(name, accepted)?;
        require_utf8_path(&path, &format!("path for resource `{name}`"))?;
        let external = self.project.manifest.resources[name].external;
        let annotation_identity = self.project.manifest.resources[name]
            .annotation_identity
            .clone();
        let assembly = self.project.manifest.resources[name].assembly.clone();
        let identity_key = (kind, path.clone());
        let (bytes, identity) = if let Some(identity) = self.identity_cache.get(&identity_key) {
            identity.clone()
        } else {
            let identity = resolve_resource_identity(kind, &path)?;
            self.identity_cache.insert(identity_key, identity.clone());
            identity
        };
        self.resources
            .entry(name.to_owned())
            .or_insert_with(|| ResolvedResource {
                kind,
                path: path.clone(),
                external,
                bytes,
                identity,
                assembly,
                annotation_identity,
            });
        Ok(path)
    }

    /// Resolve an annotation artifact and the stable source snapshot that can be inspected during
    /// plan checking. A compiled annotation emitted by an earlier step inherits its semantics from
    /// that step's source annotation, so compile-then-query plans remain resolvable before output
    /// files exist.
    fn annotation_resolution(&mut self, name: &str) -> Result<AnnotationResolutionInput> {
        let command_path = self.resource(name, &[ResourceKind::Annotation])?;
        if name.starts_with("step:") {
            let input = self.step_inputs.get(name).with_context(|| {
                format!("internal error resolving annotation reference `{name}`")
            })?;
            let output = self
                .prior_outputs
                .get(&input.producer_step)
                .and_then(|outputs| {
                    outputs
                        .iter()
                        .find(|output| output.name == input.output_name)
                })
                .with_context(|| format!("annotation producer for `{name}` is unavailable"))?;
            let semantics = output.annotation_semantics.clone().with_context(|| {
                format!(
                    "step output `{name}` has no inherited annotation identity; register the source annotation with --assembly and --annotation-label"
                )
            })?;
            return Ok(AnnotationResolutionInput {
                command_path,
                source_path: semantics.source_path,
                metadata: ProjectAnnotationIdentity {
                    assembly: semantics.assembly,
                    annotation: semantics.annotation,
                },
                source_identity: semantics.source_identity,
                command_identity: None,
            });
        }

        let resolved = self
            .resources
            .get(name)
            .with_context(|| format!("internal error resolving annotation resource `{name}`"))?;
        let metadata = resolved.annotation_identity.clone().with_context(|| {
            format!(
                "annotation resource `{name}` has no scientific identity; re-register it with --assembly and --annotation-label"
            )
        })?;
        Ok(AnnotationResolutionInput {
            command_path: command_path.clone(),
            source_path: command_path,
            metadata,
            source_identity: resolved.identity.clone(),
            command_identity: Some(resolved.identity.clone()),
        })
    }

    fn inherited_annotation_semantics(&self, name: &str) -> Option<ResolvedAnnotationSemantics> {
        if name.starts_with("step:") {
            let input = self.step_inputs.get(name)?;
            return self
                .prior_outputs
                .get(&input.producer_step)?
                .iter()
                .find(|output| output.name == input.output_name)?
                .annotation_semantics
                .clone();
        }
        let resolved = self.resources.get(name)?;
        let metadata = resolved.annotation_identity.as_ref()?;
        Some(ResolvedAnnotationSemantics {
            assembly: metadata.assembly.clone(),
            annotation: metadata.annotation.clone(),
            source_path: resolved.path.clone(),
            source_identity: resolved.identity.clone(),
        })
    }

    fn resolved_annotation_input(
        &mut self,
        role: ResolvedAnnotationRole,
        name: &str,
    ) -> Result<ResolvedAnnotationInput> {
        let input = self.annotation_resolution(name)?;
        if input.source_identity.scheme != "full-file-blake3-v1" {
            bail!(
                "annotation `{name}` requires a full-file BLAKE3 source identity, found {}",
                input.source_identity.scheme
            );
        }
        Ok(ResolvedAnnotationInput {
            role,
            resource: name.to_owned(),
            annotation_path: input.command_path,
            source_path: input.source_path,
            assembly: input.metadata.assembly,
            annotation: input.metadata.annotation,
            source_identity: input.source_identity,
            expected_command_identity: input.command_identity,
            compatibility: Vec::new(),
        })
    }

    fn verify_annotation_coordinate_resource(
        &self,
        resource_name: &str,
        annotation: &mut ResolvedAnnotationInput,
        explanation: &mut Vec<String>,
    ) -> Result<()> {
        let (kind, path, declared_assembly) = if resource_name.starts_with("step:") {
            let input = self.step_inputs.get(resource_name).with_context(|| {
                format!("internal error resolving coordinate input `{resource_name}`")
            })?;
            let output = self
                .prior_outputs
                .get(&input.producer_step)
                .and_then(|outputs| {
                    outputs
                        .iter()
                        .find(|output| output.name == input.output_name)
                })
                .with_context(|| {
                    format!("coordinate producer for `{resource_name}` is unavailable")
                })?;
            (
                output.resource_kind,
                output.path.clone(),
                output.assembly.clone(),
            )
        } else {
            let resource = self.resources.get(resource_name).with_context(|| {
                format!("internal error resolving coordinate resource `{resource_name}`")
            })?;
            (
                resource.kind,
                resource.path.clone(),
                resource.assembly.clone(),
            )
        };
        let compatibility = inspect_declared_assembly_compatibility(
            resource_name,
            kind,
            &path,
            declared_assembly,
            &annotation.assembly,
        )?;
        explanation.push(format!(
            "assembly compatibility `{resource_name}` -> {}: {}",
            match compatibility.status {
                AssemblyCompatibilityStatus::Verified => "verified",
                AssemblyCompatibilityStatus::Unverified => "unverified",
            },
            compatibility.note
        ));
        annotation.compatibility.push(compatibility);
        Ok(())
    }

    fn verify_coordinate_resource(
        &self,
        name: &str,
        intent: &mut ResolvedBiologicalIntent,
        explanation: &mut Vec<String>,
    ) -> Result<()> {
        let (kind, path, assembly) = if name.starts_with("step:") {
            let input = self
                .step_inputs
                .get(name)
                .with_context(|| format!("internal error resolving coordinate input `{name}`"))?;
            let output = self
                .prior_outputs
                .get(&input.producer_step)
                .and_then(|outputs| {
                    outputs
                        .iter()
                        .find(|output| output.name == input.output_name)
                })
                .with_context(|| format!("coordinate producer for `{name}` is unavailable"))?;
            (
                output.resource_kind,
                output.path.clone(),
                output.assembly.clone(),
            )
        } else {
            let resource = self.resources.get(name).with_context(|| {
                format!("internal error resolving coordinate resource `{name}`")
            })?;
            (
                resource.kind,
                resource.path.clone(),
                resource.assembly.clone(),
            )
        };
        let compatibility = inspect_coordinate_compatibility(name, kind, &path, assembly, intent)?;
        explanation.push(format!(
            "assembly compatibility `{name}` -> {}: {}",
            match compatibility.status {
                AssemblyCompatibilityStatus::Verified => "verified",
                AssemblyCompatibilityStatus::Unverified => "unverified",
            },
            compatibility.note
        ));
        intent.compatibility.push(compatibility);
        Ok(())
    }

    fn verify_embedded_archives(
        &self,
        intent: &mut ResolvedBiologicalIntent,
        explanation: &mut Vec<String>,
    ) -> Result<()> {
        for key in &self.step_embedded_resources {
            let resource = &self.embedded_resources[key];
            if resource.kind != ResourceKind::Archive {
                continue;
            }
            let compatibility =
                inspect_coordinate_compatibility(key, resource.kind, &resource.path, None, intent)?;
            explanation.push(format!(
                "assembly compatibility `{key}` -> {}: {}",
                match compatibility.status {
                    AssemblyCompatibilityStatus::Verified => "verified",
                    AssemblyCompatibilityStatus::Unverified => "unverified",
                },
                compatibility.note
            ));
            intent.compatibility.push(compatibility);
        }
        Ok(())
    }

    fn step_output_resource(
        &mut self,
        reference: &str,
        accepted: &[ResourceKind],
    ) -> Result<PathBuf> {
        let mut parts = reference.split(':');
        if parts.next() != Some("step") {
            unreachable!("step output references are routed by their prefix");
        }
        let producer_step = parts.next().unwrap_or_default();
        let requested_output = parts.next();
        if producer_step.is_empty() || parts.next().is_some() {
            bail!(
                "invalid step output reference `{reference}`; use step:<id> or step:<id>:<output-name>"
            );
        }
        let Some(outputs) = self.prior_outputs.get(producer_step) else {
            if self.all_step_ids.contains(producer_step) {
                bail!(
                    "step `{}` references `{reference}` before its producer; step outputs may only be consumed by later declarations",
                    self.current_step
                );
            }
            bail!("step output reference `{reference}` names unknown step `{producer_step}`");
        };
        let output = if let Some(output_name) = requested_output {
            if output_name.is_empty() {
                bail!("invalid step output reference `{reference}`; output name must not be empty");
            }
            outputs
                .iter()
                .find(|output| output.name == output_name)
                .with_context(|| {
                    let names = outputs
                        .iter()
                        .map(|output| output.name.as_str())
                        .collect::<Vec<_>>()
                        .join(", ");
                    format!(
                        "step `{producer_step}` has no output `{output_name}`; available outputs: {names}"
                    )
                })?
        } else {
            if outputs.len() != 1 {
                let names = outputs
                    .iter()
                    .map(|output| output.name.as_str())
                    .collect::<Vec<_>>()
                    .join(", ");
                bail!("step output reference `{reference}` is ambiguous; choose one of: {names}");
            }
            &outputs[0]
        };
        if !accepted.contains(&output.resource_kind) {
            let expected = accepted
                .iter()
                .map(|kind| kind.as_str())
                .collect::<Vec<_>>()
                .join(" or ");
            bail!(
                "step output `{reference}` is {}, but this input requires {expected}",
                output.resource_kind.as_str()
            );
        }
        let resolved = ResolvedStepInput {
            reference: reference.to_owned(),
            producer_step: producer_step.to_owned(),
            output_name: output.name.clone(),
            path: output.path.clone(),
            kind: output.kind,
            resource_kind: output.resource_kind,
        };
        self.step_inputs.insert(reference.to_owned(), resolved);
        Ok(output.path.clone())
    }

    #[expect(
        clippy::too_many_arguments,
        reason = "the arguments preserve the distinct source-resource provenance fields"
    )]
    fn embedded_resource(
        &mut self,
        key: String,
        owner_resource: &str,
        sample: &str,
        role: &str,
        declared_path: &str,
        kind: ResourceKind,
        path: PathBuf,
    ) -> Result<()> {
        require_utf8_path(
            &path,
            &format!("embedded {role} path for sample `{sample}`"),
        )?;
        let identity_key = (kind, path.clone());
        let (bytes, identity) = if let Some(identity) = self.identity_cache.get(&identity_key) {
            identity.clone()
        } else {
            let identity = resolve_resource_identity(kind, &path)?;
            self.identity_cache.insert(identity_key, identity.clone());
            identity
        };
        let embedded = ResolvedEmbeddedResource {
            owner_resource: owner_resource.to_owned(),
            sample: sample.to_owned(),
            role: role.to_owned(),
            declared_path: declared_path.to_owned(),
            kind,
            external: !path.starts_with(&self.project.root),
            path,
            bytes,
            identity,
        };
        if self
            .embedded_resources
            .insert(key.clone(), embedded)
            .is_some()
        {
            bail!("duplicate embedded resource key `{key}`");
        }
        self.step_embedded_resources.insert(key);
        Ok(())
    }

    fn output(&self, path: &Path, kind: ResolvedOutputKind) -> Result<ResolvedOutput> {
        self.typed_output(path, kind, "output", ResourceKind::File)
    }

    fn typed_output(
        &self,
        path: &Path,
        kind: ResolvedOutputKind,
        name: &str,
        resource_kind: ResourceKind,
    ) -> Result<ResolvedOutput> {
        let path = self.project.resolve_output(path)?;
        require_utf8_path(&path, "resolved output path")?;
        Ok(ResolvedOutput {
            name: name.to_owned(),
            path,
            kind,
            resource_kind,
            staging_path: PathBuf::new(),
            annotation_semantics: None,
            assembly: None,
        })
    }
}

fn inspect_coordinate_compatibility(
    resource_name: &str,
    kind: ResourceKind,
    path: &Path,
    declared_assembly: Option<String>,
    intent: &ResolvedBiologicalIntent,
) -> Result<ResolvedAssemblyCompatibility> {
    if let Some(assembly) = &declared_assembly {
        if assembly != &intent.assembly {
            bail!(
                "assembly mismatch: feature uses `{}` but coordinate resource `{resource_name}` is registered as `{assembly}`",
                intent.assembly
            );
        }
    }
    let mut chromosome_digest = None;
    let mut genome_digest = None;
    let mut locus_checked = false;
    if kind == ResourceKind::Archive && path.is_file() {
        if let Ok(archive) = crate::archivecmd::LazyArchive::open(path) {
            if !archive
                .chrom_names
                .iter()
                .any(|contig| contig == &intent.contig)
            {
                bail!(
                    "coordinate resource `{resource_name}` has no contig `{}` required by feature `{}`",
                    intent.contig,
                    intent.requested
                );
            }
            locus_checked = true;
            chromosome_digest = Some(format!("blake3:{}", archive.chrom_digest));
            genome_digest = archive.genome_sig.map(|signature| signature.digest);
        }
    }
    let (status, note) = if declared_assembly.is_some() {
        let suffix = if locus_checked {
            "; archive chromosome dictionary also contains the resolved contig"
        } else {
            "; no parseable chromosome/genome signature was available for an additional check"
        };
        (
            AssemblyCompatibilityStatus::Verified,
            format!("declared assembly matches `{}`{suffix}", intent.assembly),
        )
    } else {
        let evidence = if locus_checked {
            "the archive chromosome dictionary contains the resolved contig, but that does not identify an assembly"
        } else {
            "no declared assembly or parseable chromosome/genome signature establishes an assembly identity"
        };
        (
            AssemblyCompatibilityStatus::Unverified,
            format!("{evidence}; compatibility remains unverified"),
        )
    };
    Ok(ResolvedAssemblyCompatibility {
        resource: resource_name.to_owned(),
        kind,
        status,
        declared_assembly,
        chromosome_digest,
        genome_digest,
        note,
    })
}

fn inspect_declared_assembly_compatibility(
    resource_name: &str,
    kind: ResourceKind,
    path: &Path,
    declared_assembly: Option<String>,
    annotation_assembly: &str,
) -> Result<ResolvedAssemblyCompatibility> {
    if let Some(assembly) = &declared_assembly {
        if assembly != annotation_assembly {
            bail!(
                "assembly mismatch: annotation uses `{annotation_assembly}` but coordinate resource `{resource_name}` is registered as `{assembly}`"
            );
        }
    }
    let mut chromosome_digest = None;
    let mut genome_digest = None;
    if kind == ResourceKind::Archive && path.is_file() {
        if let Ok(archive) = crate::archivecmd::LazyArchive::open(path) {
            chromosome_digest = Some(format!("blake3:{}", archive.chrom_digest));
            genome_digest = archive.genome_sig.map(|signature| signature.digest);
        }
    }
    let (status, note) = if declared_assembly.is_some() {
        (
            AssemblyCompatibilityStatus::Verified,
            format!("declared assembly matches `{annotation_assembly}`"),
        )
    } else {
        (
            AssemblyCompatibilityStatus::Unverified,
            "no declared coordinate-resource assembly establishes compatibility; the result must retain this unverified status".to_owned(),
        )
    };
    Ok(ResolvedAssemblyCompatibility {
        resource: resource_name.to_owned(),
        kind,
        status,
        declared_assembly,
        chromosome_digest,
        genome_digest,
        note,
    })
}

fn resolve_resource_identity(kind: ResourceKind, path: &Path) -> Result<(u64, ResolvedIdentity)> {
    if kind == ResourceKind::Archive {
        let before = fs::metadata(path)
            .with_context(|| format!("reading archive metadata from {}", path.display()))?;
        let reader = evidence_io::format::SectionReader::open(path)
            .with_context(|| format!("reading archive identity from {}", path.display()))?;
        if let Some(commitment) = reader.content_commitment() {
            let after = fs::metadata(path)
                .with_context(|| format!("re-reading archive metadata from {}", path.display()))?;
            ensure_file_unchanged(&before, &after, path)?;
            return Ok((
                after.len(),
                ResolvedIdentity {
                    scheme: "aie-directory-root-v2".to_owned(),
                    digest: commitment.to_hex(),
                },
            ));
        }
        // Legacy archives do not authenticate their directory. Their only exact persistent
        // identity is the complete encoded file.
    } else if kind == ResourceKind::Collection {
        let before = fs::metadata(path)
            .with_context(|| format!("reading collection metadata from {}", path.display()))?;
        let digest = crate::collectioncmd::native_collection_identity(path)
            .with_context(|| format!("reading collection identity from {}", path.display()))?;
        let after = fs::metadata(path)
            .with_context(|| format!("re-reading collection metadata from {}", path.display()))?;
        ensure_file_unchanged(&before, &after, path)?;
        return Ok((
            after.len(),
            ResolvedIdentity {
                scheme: "aicollection-directory-root-v1".to_owned(),
                digest,
            },
        ));
    }
    full_file_identity(path)
}

fn full_file_identity(path: &Path) -> Result<(u64, ResolvedIdentity)> {
    let mut file =
        fs::File::open(path).with_context(|| format!("opening {} for identity", path.display()))?;
    let before = file
        .metadata()
        .with_context(|| format!("reading metadata for {}", path.display()))?;
    if !before.is_file() {
        bail!("identity input is not a regular file: {}", path.display());
    }
    let mut hasher = blake3::Hasher::new();
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        let read = file
            .read(&mut buffer)
            .with_context(|| format!("hashing {}", path.display()))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    let after = file
        .metadata()
        .with_context(|| format!("re-reading metadata for {}", path.display()))?;
    ensure_file_unchanged(&before, &after, path)?;
    Ok((
        after.len(),
        ResolvedIdentity {
            scheme: "full-file-blake3-v1".to_owned(),
            digest: hasher.finalize().to_hex().to_string(),
        },
    ))
}

fn ensure_file_unchanged(before: &fs::Metadata, after: &fs::Metadata, path: &Path) -> Result<()> {
    let changed = before.len() != after.len()
        || before.modified().ok() != after.modified().ok()
        || !same_file(before, after);
    if changed {
        bail!(
            "file changed while its identity was being resolved: {}",
            path.display()
        );
    }
    Ok(())
}

#[cfg(unix)]
fn same_file(before: &fs::Metadata, after: &fs::Metadata) -> bool {
    use std::os::unix::fs::MetadataExt;
    before.dev() == after.dev() && before.ino() == after.ino()
}

#[cfg(not(unix))]
fn same_file(_before: &fs::Metadata, _after: &fs::Metadata) -> bool {
    true
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Check(args) => run_check(args),
        Command::Run(args) => run_plan(args),
    }
}

fn run_check(args: CheckArgs) -> Result<()> {
    let project = project_for_plan(&args.plan, args.project.as_deref())?;
    let resolved = resolve_plan(&args.plan, &project)?;
    if args.json {
        println!("{}", serde_json::to_string_pretty(&resolved)?);
        if args.explain {
            eprint!("{}", explain_plan(&resolved));
        }
    } else if args.explain {
        print!("{}", explain_plan(&resolved));
        println!("valid: {} supported step(s)", resolved.steps.len());
    } else {
        println!(
            "valid plan `{}`: {} supported step(s) in project `{}`",
            resolved.name,
            resolved.steps.len(),
            resolved.project_name
        );
    }
    Ok(())
}

fn run_plan(args: RunArgs) -> Result<()> {
    let project = project_for_plan(&args.plan, args.project.as_deref())?;
    let resolved = resolve_plan(&args.plan, &project)?;
    let plan_digest = resolved_plan_digest(&resolved)?;
    let decisions = if args.resume {
        resume_decisions(&resolved, &plan_digest)?
    } else {
        ensure_outputs_absent(&resolved)?;
        vec![ResumeDecision::Run; resolved.steps.len()]
    };
    if args.explain || args.dry_run {
        print!("{}", explain_plan(&resolved));
    }
    if args.resume {
        explain_resume(&resolved, &decisions);
    }
    if args.dry_run {
        println!("dry run: no snapshot written and no steps executed");
        return Ok(());
    }
    let snapshot = persist_resolved_plan(&resolved)?;
    println!("resolved plan: {}", snapshot.display());
    execute_plan(&resolved, &plan_digest, &decisions)?;
    println!("completed plan `{}`", resolved.name);
    Ok(())
}

fn project_for_plan(plan: &Path, explicit: Option<&Path>) -> Result<ProjectContext> {
    if explicit.is_some() {
        find_project(explicit)
    } else {
        find_project_from(plan)
    }
}

/// Parse a YAML or JSON plan. The extension is part of the contract so a file is never silently
/// interpreted using a more permissive format than its author intended.
pub fn load_plan(path: &Path) -> Result<(AnalysisPlan, Vec<u8>)> {
    let bytes = fs::read(path).with_context(|| format!("reading plan {}", path.display()))?;
    let extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    validate_declared_step_kinds(&bytes, &extension, path)?;
    let plan = match extension.as_str() {
        "yaml" | "yml" => serde_yaml::from_slice(&bytes)
            .with_context(|| format!("parsing YAML plan {}", path.display()))?,
        "json" => serde_json::from_slice(&bytes)
            .with_context(|| format!("parsing JSON plan {}", path.display()))?,
        _ => bail!(
            "plan must use a .yaml, .yml, or .json extension: {}",
            path.display()
        ),
    };
    Ok((plan, bytes))
}

#[derive(Deserialize)]
struct KindEnvelope {
    steps: Vec<KindOnly>,
}

#[derive(Deserialize)]
struct KindOnly {
    kind: String,
}

fn validate_declared_step_kinds(bytes: &[u8], extension: &str, path: &Path) -> Result<()> {
    let envelope: KindEnvelope = match extension {
        "yaml" | "yml" => serde_yaml::from_slice(bytes)
            .with_context(|| format!("reading step kinds from {}", path.display()))?,
        "json" => serde_json::from_slice(bytes)
            .with_context(|| format!("reading step kinds from {}", path.display()))?,
        _ => return Ok(()),
    };
    for (index, step) in envelope.steps.iter().enumerate() {
        if !SUPPORTED_KINDS.contains(&step.kind.as_str()) {
            bail!(
                "unsupported plan step kind `{}` at steps[{}]; supported kinds: {}",
                step.kind,
                index,
                SUPPORTED_KINDS.join(", ")
            );
        }
    }
    Ok(())
}

/// Validate and compile a plan into the existing CLI command surfaces.
pub fn resolve_plan(path: &Path, project: &ProjectContext) -> Result<ResolvedPlan> {
    let source_path =
        fs::canonicalize(path).with_context(|| format!("resolving plan {}", path.display()))?;
    require_utf8_path(&project.root, "project root")?;
    require_utf8_path(&source_path, "plan path")?;
    if !source_path.starts_with(&project.root) {
        bail!(
            "plan {} is outside project {}; keep plans inside the workspace",
            source_path.display(),
            project.root.display()
        );
    }
    let (plan, source_bytes) = load_plan(&source_path)?;
    resolve_plan_loaded(source_path, plan, source_bytes, project)
}

/// Compile an already-parsed plan without creating a source file. The synthetic source name is
/// used only in provenance/display and must be one plain UTF-8 filename. The source digest binds
/// the deterministic JSON serialization of the supplied model.
pub fn resolve_plan_model(
    plan: AnalysisPlan,
    synthetic_source_name: &str,
    project: &ProjectContext,
) -> Result<ResolvedPlan> {
    if synthetic_source_name.is_empty()
        || synthetic_source_name.chars().any(char::is_control)
        || Path::new(synthetic_source_name)
            .file_name()
            .and_then(|name| name.to_str())
            != Some(synthetic_source_name)
    {
        bail!("synthetic plan source must be one nonempty printable UTF-8 filename");
    }
    require_utf8_path(&project.root, "project root")?;
    let source_path = project.root.join("plans").join(synthetic_source_name);
    let source_bytes = serde_json::to_vec(&plan).context("serializing in-memory analysis plan")?;
    resolve_plan_loaded(source_path, plan, source_bytes, project)
}

fn resolve_plan_loaded(
    source_path: PathBuf,
    plan: AnalysisPlan,
    source_bytes: Vec<u8>,
    project: &ProjectContext,
) -> Result<ResolvedPlan> {
    if plan.schema_version != PLAN_SCHEMA_VERSION {
        bail!(
            "unsupported plan schema version {}; this aie supports version {}",
            plan.schema_version,
            PLAN_SCHEMA_VERSION
        );
    }
    if plan.steps.is_empty() {
        bail!("plan must contain at least one step");
    }
    let name = plan.name.unwrap_or_else(|| {
        source_path
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or("analysis")
            .to_owned()
    });
    validate_plan_name(&name)?;

    let mut ids = BTreeSet::new();
    for step in &plan.steps {
        validate_step_id(step.id())?;
        if !ids.insert(step.id().to_owned()) {
            bail!("duplicate plan step id `{}`", step.id());
        }
    }
    let mut resolver = Resolver::new(project, ids);
    let mut steps = Vec::with_capacity(plan.steps.len());
    for step in &plan.steps {
        let resolved = compile_step(step, &mut resolver)
            .with_context(|| format!("resolving step `{}` ({})", step.id(), step.kind()))?;
        resolver.register_outputs(&resolved);
        steps.push(resolved);
    }
    validate_output_collisions(&steps)?;

    let executable = std::env::current_exe().context("locating the aie executable")?;
    let (_, executable_identity) = full_file_identity(&executable)
        .context("hashing the aie executable for resolved-plan provenance")?;
    Ok(ResolvedPlan {
        schema_version: RESOLVED_PLAN_SCHEMA_VERSION,
        plan_schema_version: PLAN_SCHEMA_VERSION,
        producer: ResolvedProducer {
            // The crates.io package is named `gravlax`, but the executable and
            // its result/provenance wire identity remain `aie`.
            name: "aie".to_owned(),
            version: env!("CARGO_PKG_VERSION").to_owned(),
            plan_engine: "aie-declarative-plan-v2".to_owned(),
            executable_identity,
        },
        name,
        source_path,
        source_digest: blake3::hash(&source_bytes).to_hex().to_string(),
        project_name: project.manifest.name.clone(),
        project_root: project.root.clone(),
        project_manifest: project.manifest_path.clone(),
        project_manifest_digest: project.manifest_digest.clone(),
        resources: resolver.resources,
        embedded_resources: resolver.embedded_resources,
        steps,
    })
}

fn compile_step(step: &PlanStep, resolver: &mut Resolver<'_>) -> Result<ResolvedStep> {
    resolver.begin_step(step.id());
    let mut args = Vec::<String>::new();
    let mut explanation = Vec::<String>::new();
    let mut outputs = Vec::<ResolvedOutput>::new();
    let mut prepared_inputs = Vec::<ResolvedPreparedInput>::new();
    let mut biological_intent = None;
    let mut annotation_inputs = Vec::<ResolvedAnnotationInput>::new();
    let mut annotation_comparison = None;
    let mut stdout = None;
    let mut uniform_io = None;

    match step {
        PlanStep::InspectArchive {
            archive,
            verify_content,
            format,
            output,
            uniform_output,
            ..
        } => {
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "inspect-archive",
            )?;
            if uniform_output.is_none() && *format == OutputFormat::Tsv {
                bail!("inspect-archive supports human or json format, not tsv");
            }
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            args.extend(["inspect-archive".to_owned(), path_arg(&archive_path)]);
            if *verify_content {
                args.push("--verify-content".to_owned());
            }
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                if *format == OutputFormat::Json {
                    args.push("--json".to_owned());
                }
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::QueryRegion {
            archive,
            locus,
            feature,
            top,
            annotation,
            format,
            output,
            uniform_output,
            scope,
            ..
        } => {
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "query-region",
            )?;
            let (locus, intent) = resolve_locus_or_feature(
                "query-region",
                locus,
                feature,
                annotation,
                resolver,
                &mut explanation,
            )?;
            biological_intent = intent;
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            if let Some(intent) = biological_intent.as_mut() {
                resolver.verify_coordinate_resource(archive, intent, &mut explanation)?;
            }
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            args.extend([
                "query".to_owned(),
                path_arg(&archive_path),
                "region".to_owned(),
                locus,
                "--top".to_owned(),
                top.to_string(),
            ]);
            if let Some(annotation) = annotation {
                let annotation_path = resolver.resource(annotation, &[ResourceKind::Annotation])?;
                explain_resource(&mut explanation, "annotation", annotation, &annotation_path);
                args.extend(["--gtf".to_owned(), path_arg(&annotation_path)]);
            }
            resolve_query_scope(scope, resolver, &mut args, &mut explanation)?;
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                push_format(&mut args, *format);
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::QueryJunction {
            archive,
            locus,
            top,
            format,
            output,
            uniform_output,
            scope,
            ..
        } => {
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "query-junction",
            )?;
            validate_locus(locus)?;
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            args.extend([
                "query".to_owned(),
                path_arg(&archive_path),
                "junction".to_owned(),
                locus.clone(),
                "--top".to_owned(),
                top.to_string(),
            ]);
            resolve_query_scope(scope, resolver, &mut args, &mut explanation)?;
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                push_format(&mut args, *format);
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::QueryJunctions {
            archive,
            locus,
            feature,
            either,
            min_support,
            with_cells,
            min_cells,
            annotation,
            format,
            output,
            uniform_output,
            scope,
            ..
        } => {
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "query-junctions",
            )?;
            let (locus, intent) = resolve_locus_or_feature(
                "query-junctions",
                locus,
                feature,
                annotation,
                resolver,
                &mut explanation,
            )?;
            biological_intent = intent;
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            if let Some(intent) = biological_intent.as_mut() {
                resolver.verify_coordinate_resource(archive, intent, &mut explanation)?;
            }
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            args.extend([
                "query".to_owned(),
                path_arg(&archive_path),
                "junctions".to_owned(),
                locus,
                "--min-support".to_owned(),
                min_support.to_string(),
            ]);
            if *either {
                args.push("--either".to_owned());
            }
            if *with_cells {
                args.push("--with-cells".to_owned());
            }
            if *min_cells != 0 {
                args.extend(["--min-cells".to_owned(), min_cells.to_string()]);
            }
            if let Some(annotation) = annotation {
                let annotation_path = resolver.resource(annotation, &[ResourceKind::Annotation])?;
                explain_resource(&mut explanation, "annotation", annotation, &annotation_path);
                args.extend(["--gtf".to_owned(), path_arg(&annotation_path)]);
            }
            resolve_query_scope(scope, resolver, &mut args, &mut explanation)?;
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                push_format(&mut args, *format);
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::QueryJset {
            archive,
            include,
            exclude,
            top,
            format,
            output,
            uniform_output,
            scope,
            ..
        } => {
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "query-jset",
            )?;
            validate_junction_sets(include, exclude)?;
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            args.extend([
                "query".to_owned(),
                path_arg(&archive_path),
                "jset".to_owned(),
            ]);
            push_repeated(&mut args, "--include", include);
            push_repeated(&mut args, "--exclude", exclude);
            args.extend(["--top".to_owned(), top.to_string()]);
            resolve_query_scope(scope, resolver, &mut args, &mut explanation)?;
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                push_format(&mut args, *format);
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::QueryEvents {
            archive,
            locus,
            feature,
            event_types,
            min_support,
            min_informative,
            max_events,
            top,
            annotation,
            format,
            output,
            uniform_output,
            scope,
            ..
        } => {
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "query-events",
            )?;
            let (locus, intent) = resolve_locus_or_feature(
                "query-events",
                locus,
                feature,
                annotation,
                resolver,
                &mut explanation,
            )?;
            biological_intent = intent;
            if *max_events == 0 {
                bail!("max_events must be greater than zero");
            }
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            if let Some(intent) = biological_intent.as_mut() {
                resolver.verify_coordinate_resource(archive, intent, &mut explanation)?;
            }
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            args.extend([
                "query".to_owned(),
                path_arg(&archive_path),
                "events".to_owned(),
                locus,
                "--min-support".to_owned(),
                min_support.to_string(),
                "--min-informative".to_owned(),
                min_informative.to_string(),
                "--max-events".to_owned(),
                max_events.to_string(),
                "--top".to_owned(),
                top.to_string(),
            ]);
            push_event_types(&mut args, event_types);
            if let Some(annotation) = annotation {
                let annotation_path = resolver.resource(annotation, &[ResourceKind::Annotation])?;
                explain_resource(&mut explanation, "annotation", annotation, &annotation_path);
                args.extend(["--gtf".to_owned(), path_arg(&annotation_path)]);
            }
            resolve_query_scope(scope, resolver, &mut args, &mut explanation)?;
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                push_format(&mut args, *format);
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::QueryApa {
            archive,
            locus,
            feature,
            annotation,
            site_gap,
            strand,
            tsv,
            groups,
            genome,
            drop_ip,
            permute,
            seed,
            plot,
            output,
            uniform_output,
            ..
        } => {
            reject_legacy_result_options(uniform_output, *tsv, output.is_some(), "query-apa")?;
            if uniform_output.is_some() && plot.is_some() {
                bail!("query-apa uniform_output cannot be combined with plot");
            }
            let (locus, intent) = resolve_locus_or_feature(
                "query-apa",
                locus,
                feature,
                annotation,
                resolver,
                &mut explanation,
            )?;
            biological_intent = intent;
            if *site_gap == 0 {
                bail!("site_gap must be greater than zero");
            }
            if *drop_ip && genome.is_none() {
                bail!("query-apa drop_ip requires a named genome resource");
            }
            if *permute > 0 && groups.is_none() {
                bail!("query-apa permute requires a named groups resource");
            }
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            if let Some(intent) = biological_intent.as_mut() {
                resolver.verify_coordinate_resource(archive, intent, &mut explanation)?;
            }
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            args.extend([
                "query".to_owned(),
                path_arg(&archive_path),
                "apa".to_owned(),
                locus,
                "--site-gap".to_owned(),
                site_gap.to_string(),
                "--seed".to_owned(),
                seed.to_string(),
            ]);
            if let Some(strand) = strand {
                args.extend([
                    "--strand".to_owned(),
                    match strand {
                        ApaStrand::Forward => "+",
                        ApaStrand::Reverse => "-",
                    }
                    .to_owned(),
                ]);
            }
            if *tsv {
                args.push("--tsv".to_owned());
            }
            if let Some(groups) = groups {
                let groups_path =
                    resolver.resource(groups, &[ResourceKind::Groups, ResourceKind::Metadata])?;
                explain_resource(&mut explanation, "groups", groups, &groups_path);
                args.extend(["--groups".to_owned(), path_arg(&groups_path)]);
            }
            if let Some(genome) = genome {
                let genome_path = resolver.resource(genome, &[ResourceKind::Genome])?;
                if let Some(intent) = biological_intent.as_mut() {
                    resolver.verify_coordinate_resource(genome, intent, &mut explanation)?;
                }
                explain_resource(&mut explanation, "genome", genome, &genome_path);
                args.extend(["--genome".to_owned(), path_arg(&genome_path)]);
            }
            if *drop_ip {
                args.push("--drop-ip".to_owned());
            }
            if *permute > 0 {
                args.extend(["--permute".to_owned(), permute.to_string()]);
            }
            if let Some(plot) = plot {
                let resolved = resolver.typed_output(
                    plot,
                    ResolvedOutputKind::File,
                    "plot",
                    ResourceKind::File,
                )?;
                require_extension(&resolved.path, "svg", "query-apa plot output")?;
                explanation.push(format!("plot -> {}", resolved.path.display()));
                args.extend(["--plot".to_owned(), path_arg(&resolved.path)]);
                outputs.push(resolved);
            }
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::CompareAnnotations {
            archive,
            annotation_a,
            annotation_b,
            gene_key,
            solo_strand,
            max_molecule_witnesses,
            max_row_transitions_per_molecule,
            allow_identical,
            format,
            table,
            output,
            ..
        } => {
            validate_bundle_format(
                *format,
                table.map(annotation_comparison_table_arg),
                "compare-annotations",
            )?;
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            let mut annotation_a_input =
                resolver.resolved_annotation_input(ResolvedAnnotationRole::A, annotation_a)?;
            let mut annotation_b_input =
                resolver.resolved_annotation_input(ResolvedAnnotationRole::B, annotation_b)?;
            if annotation_a_input.assembly != annotation_b_input.assembly {
                bail!(
                    "annotation assembly mismatch: `{annotation_a}` is `{}` but `{annotation_b}` is `{}`",
                    annotation_a_input.assembly,
                    annotation_b_input.assembly
                );
            }
            if annotation_a_input.source_identity == annotation_b_input.source_identity
                && !allow_identical
            {
                bail!(
                    "annotation resources `{annotation_a}` and `{annotation_b}` resolve to identical source content; set allow_identical: true only for an intentional zero-delta control"
                );
            }
            validate_comparison_annotation_source(
                &annotation_a_input.source_path,
                *gene_key,
                "annotation A",
            )?;
            validate_comparison_annotation_source(
                &annotation_b_input.source_path,
                *gene_key,
                "annotation B",
            )?;
            resolver.verify_annotation_coordinate_resource(
                archive,
                &mut annotation_a_input,
                &mut explanation,
            )?;
            resolver.verify_annotation_coordinate_resource(
                archive,
                &mut annotation_b_input,
                &mut explanation,
            )?;
            explain_resource(
                &mut explanation,
                "annotation A",
                annotation_a,
                &annotation_a_input.annotation_path,
            );
            explain_resource(
                &mut explanation,
                "annotation B",
                annotation_b,
                &annotation_b_input.annotation_path,
            );
            args.extend([
                "compare-annotations".to_owned(),
                path_arg(&archive_path),
                "--annotation-a".to_owned(),
                path_arg(&annotation_a_input.annotation_path),
                "--annotation-b".to_owned(),
                path_arg(&annotation_b_input.annotation_path),
                "--assembly".to_owned(),
                annotation_a_input.assembly.clone(),
                "--annotation-a-label".to_owned(),
                annotation_a_input.annotation.clone(),
                "--annotation-b-label".to_owned(),
                annotation_b_input.annotation.clone(),
            ]);
            push_expected_annotation_digest(
                &mut args,
                "--annotation-a-digest",
                annotation_a_input.expected_command_identity.as_ref(),
            )?;
            push_expected_annotation_digest(
                &mut args,
                "--annotation-b-digest",
                annotation_b_input.expected_command_identity.as_ref(),
            )?;
            args.extend([
                "--gene-key".to_owned(),
                comparison_gene_key_arg(*gene_key).to_owned(),
                "--solo-strand".to_owned(),
                solo_strand_arg(*solo_strand).to_owned(),
                "--max-molecule-witnesses".to_owned(),
                max_molecule_witnesses.to_string(),
                "--max-row-transitions-per-molecule".to_owned(),
                max_row_transitions_per_molecule.to_string(),
            ]);
            if *allow_identical {
                args.push("--allow-identical".to_owned());
            }
            push_named_format(&mut args, *format);
            if let Some(table) = table {
                args.extend([
                    "--table".to_owned(),
                    annotation_comparison_table_arg(*table).to_owned(),
                ]);
            }
            resolve_capability_output(output, resolver, &mut args, &mut outputs, &mut explanation)?;
            annotation_comparison = Some(ResolvedAnnotationComparisonIntent {
                annotation_a_resource: annotation_a.clone(),
                annotation_b_resource: annotation_b.clone(),
                assembly: annotation_a_input.assembly.clone(),
                gene_key: *gene_key,
                solo_strand: *solo_strand,
                final_count_delta_semantics: "exact signed B-minus-A difference between two independent complete annotation-specific assignment and UMI-collapse reductions".to_owned(),
                transition_evidence_semantics: "archive UMI-class transitions and non-exclusive contributing causes; these are not additive or uniquely attributable contributions to final count deltas".to_owned(),
            });
            annotation_inputs.extend([annotation_a_input, annotation_b_input]);
            explanation.push(
                "execution -> one full archive scan feeding two independent reductions; no locus or output filter is used as an I/O bound".to_owned(),
            );
            explanation.push(
                "interpretation -> final count deltas are exact; class transitions, causes, and bounded witnesses are explanatory evidence, not a decomposition of the nonlinear delta".to_owned(),
            );
        }
        PlanStep::QueryTranscriptEcs {
            archive,
            locus,
            feature,
            annotation,
            solo_strand,
            scope,
            emit_membership,
            max_ecs,
            max_memberships,
            format,
            table,
            output,
            ..
        } => {
            validate_bundle_format(
                *format,
                table.map(transcript_ec_table_arg),
                "query-transcript-ecs",
            )?;
            if *max_ecs == 0 {
                bail!("max_ecs must be greater than zero");
            }
            if *max_memberships == 0 {
                bail!("max_memberships must be greater than zero");
            }
            if *table == Some(TranscriptEcTable::Membership) && !emit_membership {
                bail!("query-transcript-ecs table membership requires emit_membership: true");
            }
            let annotation_name = Some(annotation.clone());
            let (resolved_locus, intent) = resolve_locus_or_feature(
                "query-transcript-ecs",
                locus,
                feature,
                &annotation_name,
                resolver,
                &mut explanation,
            )?;
            biological_intent = intent;
            if let Some(intent) = &biological_intent {
                if intent.resolved_kind != FeatureKind::Gene {
                    bail!(
                        "query-transcript-ecs feature `{}` resolved as {:?}; transcript equivalence queries require a gene feature or an explicit locus",
                        intent.requested,
                        intent.resolved_kind
                    );
                }
            }
            let mut annotation_input =
                resolver.resolved_annotation_input(ResolvedAnnotationRole::Query, annotation)?;
            validate_transcript_ec_annotation_source(&annotation_input.source_path)?;
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            if let Some(intent) = biological_intent.as_mut() {
                resolver.verify_coordinate_resource(archive, intent, &mut explanation)?;
                annotation_input.compatibility = intent.compatibility.clone();
            } else {
                resolver.verify_annotation_coordinate_resource(
                    archive,
                    &mut annotation_input,
                    &mut explanation,
                )?;
            }
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            explain_resource(
                &mut explanation,
                "annotation",
                annotation,
                &annotation_input.annotation_path,
            );
            args.extend([
                "query".to_owned(),
                path_arg(&archive_path),
                "transcript-ecs".to_owned(),
                "--annotation-file".to_owned(),
                path_arg(&annotation_input.annotation_path),
                "--assembly".to_owned(),
                annotation_input.assembly.clone(),
                "--annotation-label".to_owned(),
                annotation_input.annotation.clone(),
            ]);
            push_expected_annotation_digest(
                &mut args,
                "--annotation-digest",
                annotation_input.expected_command_identity.as_ref(),
            )?;
            if let Some(intent) = &biological_intent {
                args.extend(["--feature".to_owned(), format!("gene:{}", intent.stable_id)]);
            } else {
                args.extend(["--locus".to_owned(), resolved_locus]);
            }
            args.extend([
                "--solo-strand".to_owned(),
                solo_strand_arg(*solo_strand).to_owned(),
                "--max-ecs".to_owned(),
                max_ecs.to_string(),
                "--max-memberships".to_owned(),
                max_memberships.to_string(),
            ]);
            resolve_query_scope(scope, resolver, &mut args, &mut explanation)?;
            if *emit_membership {
                args.push("--emit-membership".to_owned());
            }
            push_named_format(&mut args, *format);
            if let Some(table) = table {
                args.extend([
                    "--table".to_owned(),
                    transcript_ec_table_arg(*table).to_owned(),
                ]);
            }
            resolve_capability_output(output, resolver, &mut args, &mut outputs, &mut explanation)?;
            annotation_inputs.push(annotation_input);
            explanation.push(
                "execution -> full archive scan so stored alternative placements and every occurrence of a global UMI class are considered; no execution read upper bound is claimed".to_owned(),
            );
            explanation.push(
                "interpretation -> transcript equivalence classes are deterministic compatibility sets from retained archive geometry, not abundance estimates, isoform calls, or molecule-to-molecule phasing".to_owned(),
            );
        }
        PlanStep::FederateJunction {
            archives,
            locus,
            top,
            output,
            uniform_output,
            ..
        } => {
            reject_legacy_result_options(
                uniform_output,
                false,
                output.is_some(),
                "federate-junction",
            )?;
            validate_locus(locus)?;
            if archives.len() < 2 {
                bail!("federate-junction requires at least two archives");
            }
            args.push("federate".to_owned());
            let mut seen_archives = BTreeMap::<PathBuf, String>::new();
            for archive in archives {
                let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
                if let Some(previous) = seen_archives.insert(archive_path.clone(), archive.clone())
                {
                    bail!(
                        "federate-junction resources `{previous}` and `{archive}` resolve to the same archive {}",
                        archive_path.display()
                    );
                }
                explain_resource(&mut explanation, "archive", archive, &archive_path);
                args.push(path_arg(&archive_path));
            }
            args.extend([locus.clone(), "--top".to_owned(), top.to_string()]);
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::CollectionRegion {
            collection,
            locus,
            feature,
            annotation,
            top,
            explain_routing,
            verify_content,
            format,
            output,
            uniform_output,
            ..
        } => {
            let (locus, intent) = resolve_locus_or_feature(
                "collection-region",
                locus,
                feature,
                annotation,
                resolver,
                &mut explanation,
            )?;
            biological_intent = intent;
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "collection-region",
            )?;
            if uniform_output.is_none() {
                require_human_or_json(*format, "collection-region")?;
            }
            let collection_path = resolver.resource(collection, &[ResourceKind::Collection])?;
            if let Some(intent) = biological_intent.as_mut() {
                resolver.verify_coordinate_resource(collection, intent, &mut explanation)?;
            }
            explain_resource(&mut explanation, "collection", collection, &collection_path);
            args.extend([
                "collection".to_owned(),
                "region".to_owned(),
                path_arg(&collection_path),
                locus,
                "--top".to_owned(),
                top.to_string(),
            ]);
            push_collection_flags(
                &mut args,
                *explain_routing,
                *verify_content,
                if uniform_output.is_some() {
                    OutputFormat::Human
                } else {
                    *format
                },
            );
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::CollectionJunction {
            collection,
            locus,
            min_support,
            top,
            explain_routing,
            verify_content,
            format,
            output,
            uniform_output,
            ..
        } => {
            validate_locus(locus)?;
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "collection-junction",
            )?;
            if uniform_output.is_none() {
                require_human_or_json(*format, "collection-junction")?;
            }
            let collection_path = resolver.resource(collection, &[ResourceKind::Collection])?;
            explain_resource(&mut explanation, "collection", collection, &collection_path);
            args.extend([
                "collection".to_owned(),
                "junction".to_owned(),
                path_arg(&collection_path),
                locus.clone(),
                "--min-support".to_owned(),
                min_support.to_string(),
                "--top".to_owned(),
                top.to_string(),
            ]);
            push_collection_flags(
                &mut args,
                *explain_routing,
                *verify_content,
                if uniform_output.is_some() {
                    OutputFormat::Human
                } else {
                    *format
                },
            );
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::CollectionJset {
            collection,
            include,
            exclude,
            min_support,
            top,
            explain_routing,
            verify_content,
            format,
            output,
            uniform_output,
            ..
        } => {
            validate_junction_sets(include, exclude)?;
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "collection-jset",
            )?;
            if uniform_output.is_none() {
                require_human_or_json(*format, "collection-jset")?;
            }
            let collection_path = resolver.resource(collection, &[ResourceKind::Collection])?;
            explain_resource(&mut explanation, "collection", collection, &collection_path);
            args.extend([
                "collection".to_owned(),
                "jset".to_owned(),
                path_arg(&collection_path),
            ]);
            push_repeated(&mut args, "--include", include);
            push_repeated(&mut args, "--exclude", exclude);
            args.extend([
                "--min-support".to_owned(),
                min_support.to_string(),
                "--top".to_owned(),
                top.to_string(),
            ]);
            push_collection_flags(
                &mut args,
                *explain_routing,
                *verify_content,
                if uniform_output.is_some() {
                    OutputFormat::Human
                } else {
                    *format
                },
            );
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::CohortEvents {
            samples,
            groups,
            locus,
            feature,
            event_types,
            min_support,
            min_samples,
            min_informative,
            min_row_informative,
            max_events,
            annotation,
            format,
            output,
            uniform_output,
            ..
        } => {
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "cohort-events",
            )?;
            let (locus, intent) = resolve_locus_or_feature(
                "cohort-events",
                locus,
                feature,
                annotation,
                resolver,
                &mut explanation,
            )?;
            biological_intent = intent;
            if samples.len() < 2 {
                bail!("cohort-events requires at least two named samples");
            }
            if *min_samples == 0 || *min_samples > samples.len() {
                bail!(
                    "min_samples must be between 1 and the {} configured samples",
                    samples.len()
                );
            }
            if *max_events == 0 {
                bail!("max_events must be greater than zero");
            }
            args.extend(["cohort".to_owned(), "events".to_owned(), locus]);
            let mut seen_archives = BTreeMap::<PathBuf, String>::new();
            for (sample_id, resource_name) in samples {
                validate_cohort_id(sample_id, "sample")?;
                let archive_path = resolver.resource(resource_name, &[ResourceKind::Archive])?;
                if let Some(intent) = biological_intent.as_mut() {
                    resolver.verify_coordinate_resource(resource_name, intent, &mut explanation)?;
                }
                if let Some(previous_sample) =
                    seen_archives.insert(archive_path.clone(), sample_id.clone())
                {
                    bail!(
                        "cohort samples `{previous_sample}` and `{sample_id}` resolve to the same archive {}",
                        archive_path.display()
                    );
                }
                explanation.push(format!(
                    "sample `{sample_id}` uses archive `{resource_name}` -> {}",
                    archive_path.display()
                ));
                args.extend([
                    "--sample".to_owned(),
                    format!("{sample_id}={}", archive_path.display()),
                ]);
            }
            for (sample_id, resource_name) in groups {
                validate_cohort_id(sample_id, "groups sample")?;
                if !samples.contains_key(sample_id) {
                    bail!("groups entry `{sample_id}` has no matching cohort sample");
                }
                let groups_path = resolver.resource(
                    resource_name,
                    &[ResourceKind::Groups, ResourceKind::Metadata],
                )?;
                explanation.push(format!(
                    "groups for `{sample_id}` use `{resource_name}` -> {}",
                    groups_path.display()
                ));
                args.extend([
                    "--groups".to_owned(),
                    format!("{sample_id}={}", groups_path.display()),
                ]);
            }
            push_event_types(&mut args, event_types);
            args.extend([
                "--min-support".to_owned(),
                min_support.to_string(),
                "--min-samples".to_owned(),
                min_samples.to_string(),
                "--min-informative".to_owned(),
                min_informative.to_string(),
                "--min-row-informative".to_owned(),
                min_row_informative.to_string(),
                "--max-events".to_owned(),
                max_events.to_string(),
            ]);
            if let Some(annotation) = annotation {
                let annotation_path = resolver.resource(annotation, &[ResourceKind::Annotation])?;
                explain_resource(&mut explanation, "annotation", annotation, &annotation_path);
                args.extend(["--gtf".to_owned(), path_arg(&annotation_path)]);
            }
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                push_format(&mut args, *format);
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::CohortSpliceGraph {
            locus,
            feature,
            annotation,
            design,
            contrast,
            counts_only,
            min_support,
            min_edge_samples,
            min_sample_umis,
            min_replicates,
            min_path_umis,
            min_path_samples,
            max_paths,
            format,
            output,
            uniform_output,
            ..
        } => {
            let (locus, intent) = resolve_locus_or_feature(
                "cohort-splice-graph",
                locus,
                feature,
                annotation,
                resolver,
                &mut explanation,
            )?;
            biological_intent = intent;
            reject_legacy_result_options(
                uniform_output,
                *format != OutputFormat::Human,
                output.is_some(),
                "cohort-splice-graph",
            )?;
            if uniform_output.is_none() {
                require_human_or_json(*format, "cohort-splice-graph")?;
            }
            if *counts_only == contrast.is_some() {
                bail!("cohort-splice-graph requires exactly one of contrast or counts_only: true");
            }
            if *min_support == 0 {
                bail!("min_support must be greater than zero");
            }
            for (name, value) in [
                ("min_edge_samples", *min_edge_samples),
                ("min_sample_umis", *min_sample_umis),
                ("min_replicates", *min_replicates),
                ("min_path_umis", *min_path_umis),
                ("min_path_samples", *min_path_samples),
                ("max_paths", *max_paths),
            ] {
                if value == 0 {
                    bail!("{name} must be greater than zero");
                }
            }
            let design_path =
                resolver.resource(design, &[ResourceKind::Design, ResourceKind::Metadata])?;
            let external_design = resolver.project.manifest.resources[design].external;
            let canonical_design =
                resolve_design_paths(resolver, design, &design_path, external_design)?;
            if let Some(intent) = biological_intent.as_mut() {
                resolver.verify_embedded_archives(intent, &mut explanation)?;
            }
            if *min_edge_samples > canonical_design.sample_count {
                bail!(
                    "min_edge_samples cannot exceed the {} design samples",
                    canonical_design.sample_count
                );
            }
            if let Some(contrast) = contrast {
                let (left, right) = validate_contrast(contrast)?;
                let expected = BTreeSet::from([left.to_owned(), right.to_owned()]);
                let observed = canonical_design
                    .condition_counts
                    .keys()
                    .cloned()
                    .collect::<BTreeSet<_>>();
                if observed != expected {
                    bail!(
                        "design conditions {:?} do not exactly match contrast {left}:{right}",
                        observed
                    );
                }
                for condition in [left, right] {
                    let replicates = canonical_design.condition_counts[condition];
                    if replicates < *min_replicates {
                        bail!(
                            "condition `{condition}` has {replicates} design sample(s), fewer than min_replicates {}",
                            min_replicates
                        );
                    }
                }
            }
            let prepared_design =
                resolved_design_input(resolver, design, canonical_design.content)?;
            explain_resource(&mut explanation, "design", design, &design_path);
            explanation.push(format!(
                "canonical design -> {}",
                prepared_design.path.display()
            ));
            args.extend([
                "cohort".to_owned(),
                "splice-graph".to_owned(),
                locus,
                "--design".to_owned(),
                path_arg(&prepared_design.path),
                "--min-support".to_owned(),
                min_support.to_string(),
                "--min-edge-samples".to_owned(),
                min_edge_samples.to_string(),
                "--min-sample-umis".to_owned(),
                min_sample_umis.to_string(),
                "--min-replicates".to_owned(),
                min_replicates.to_string(),
                "--min-path-umis".to_owned(),
                min_path_umis.to_string(),
                "--min-path-samples".to_owned(),
                min_path_samples.to_string(),
                "--max-paths".to_owned(),
                max_paths.to_string(),
            ]);
            prepared_inputs.push(prepared_design);
            if let Some(contrast) = contrast {
                args.extend(["--contrast".to_owned(), contrast.clone()]);
            } else {
                args.push("--counts-only".to_owned());
            }
            if uniform_output.is_some() {
                uniform_io = resolve_uniform_io(
                    uniform_output,
                    ResolvedUniformIoKind::Result,
                    resolver,
                    &mut args,
                    &mut outputs,
                    &mut explanation,
                )?;
            } else {
                if *format == OutputFormat::Json {
                    args.push("--json".to_owned());
                }
                resolve_stdout(
                    output,
                    resolver,
                    &mut stdout,
                    &mut outputs,
                    &mut explanation,
                )?;
            }
        }
        PlanStep::CompileAnnotation {
            annotation,
            output,
            uniform_report,
            ..
        } => {
            let annotation_path = resolver.resource(annotation, &[ResourceKind::Annotation])?;
            require_extension(&annotation_path, "gtf", "compile-annotation input")?;
            explain_resource(&mut explanation, "annotation", annotation, &annotation_path);
            let mut resolved = resolver.typed_output(
                output,
                ResolvedOutputKind::File,
                "annotation",
                ResourceKind::Annotation,
            )?;
            resolved.annotation_semantics = resolver.inherited_annotation_semantics(annotation);
            require_extension(&resolved.path, "aic", "compile-annotation output")?;
            explanation.push(format!("output -> {}", resolved.path.display()));
            args.extend([
                "compile-annotation".to_owned(),
                path_arg(&annotation_path),
                "--out".to_owned(),
                path_arg(&resolved.path),
            ]);
            outputs.push(resolved);
            uniform_io = resolve_uniform_io(
                uniform_report,
                ResolvedUniformIoKind::Report,
                resolver,
                &mut args,
                &mut outputs,
                &mut explanation,
            )?;
        }
        PlanStep::IngestArchive {
            bam,
            whitelist,
            output,
            genome,
            locus_gap,
            zstd_level,
            chunk_mb,
            uniform_report,
            ..
        } => {
            if !(-7..=22).contains(zstd_level) {
                bail!("zstd_level must be between -7 and 22");
            }
            if *chunk_mb == 0 {
                bail!("chunk_mb must be greater than zero");
            }
            let bam_path = resolver.resource(bam, &[ResourceKind::Bam])?;
            let whitelist_path = resolver.resource(whitelist, &[ResourceKind::Whitelist])?;
            explain_resource(&mut explanation, "bam", bam, &bam_path);
            explain_resource(&mut explanation, "whitelist", whitelist, &whitelist_path);
            let mut resolved = resolver.typed_output(
                output,
                ResolvedOutputKind::File,
                "archive",
                ResourceKind::Archive,
            )?;
            require_extension(&resolved.path, "aie", "ingest-archive output")?;
            explanation.push(format!("output -> {}", resolved.path.display()));
            args.extend([
                "ingest-archive".to_owned(),
                path_arg(&bam_path),
                "--whitelist".to_owned(),
                path_arg(&whitelist_path),
                "--out".to_owned(),
                path_arg(&resolved.path),
                "--locus-gap".to_owned(),
                locus_gap.to_string(),
                "--zstd-level".to_owned(),
                zstd_level.to_string(),
                "--chunk-mb".to_owned(),
                chunk_mb.to_string(),
            ]);
            if let Some(genome) = genome {
                let genome_path = resolver.resource(genome, &[ResourceKind::Genome])?;
                resolved.assembly = resolver
                    .resources
                    .get(genome)
                    .and_then(|resource| resource.assembly.clone());
                explain_resource(&mut explanation, "genome", genome, &genome_path);
                args.extend(["--genome".to_owned(), path_arg(&genome_path)]);
            }
            outputs.push(resolved);
            uniform_io = resolve_uniform_io(
                uniform_report,
                ResolvedUniformIoKind::Report,
                resolver,
                &mut args,
                &mut outputs,
                &mut explanation,
            )?;
        }
        PlanStep::ReplayRows {
            archive,
            annotation,
            barcodes,
            out_dir,
            locus_gap,
            velocity,
            solo_strand,
            uniform_report,
            ..
        } => {
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            let annotation_path = resolver.resource(annotation, &[ResourceKind::Annotation])?;
            let barcodes_path = resolver.resource(barcodes, &[ResourceKind::Barcodes])?;
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            explain_resource(&mut explanation, "annotation", annotation, &annotation_path);
            explain_resource(&mut explanation, "barcodes", barcodes, &barcodes_path);
            let resolved = resolver.typed_output(
                out_dir,
                ResolvedOutputKind::Directory,
                "directory",
                ResourceKind::File,
            )?;
            explanation.push(format!("output directory -> {}", resolved.path.display()));
            args.extend([
                "replay-rows".to_owned(),
                path_arg(&archive_path),
                "--gtf".to_owned(),
                path_arg(&annotation_path),
                "--barcodes".to_owned(),
                path_arg(&barcodes_path),
                "--out-dir".to_owned(),
                path_arg(&resolved.path),
                "--locus-gap".to_owned(),
                locus_gap.to_string(),
                "--solo-strand".to_owned(),
                solo_strand_arg(*solo_strand).to_owned(),
            ]);
            if *velocity {
                args.push("--velocity".to_owned());
            }
            outputs.push(resolved);
            uniform_io = resolve_uniform_io(
                uniform_report,
                ResolvedUniformIoKind::Report,
                resolver,
                &mut args,
                &mut outputs,
                &mut explanation,
            )?;
        }
        PlanStep::ExtendAnnotation {
            archive,
            annotation,
            out_gtf,
            report,
            genome,
            max_extend,
            evidence_gap,
            min_umis,
            min_cells,
            site_gap,
            min_extension,
            clip_any_strand,
            uniform_report,
            ..
        } => {
            let archive_path = resolver.resource(archive, &[ResourceKind::Archive])?;
            let annotation_path = resolver.resource(annotation, &[ResourceKind::Annotation])?;
            require_extension(&annotation_path, "gtf", "extend-annotation input")?;
            explain_resource(&mut explanation, "archive", archive, &archive_path);
            explain_resource(&mut explanation, "annotation", annotation, &annotation_path);
            let out_gtf = resolver.typed_output(
                out_gtf,
                ResolvedOutputKind::File,
                "annotation",
                ResourceKind::Annotation,
            )?;
            require_extension(&out_gtf.path, "gtf", "extend-annotation output")?;
            explanation.push(format!("extended annotation -> {}", out_gtf.path.display()));
            args.extend([
                "extend".to_owned(),
                path_arg(&archive_path),
                "--gtf".to_owned(),
                path_arg(&annotation_path),
                "--out-gtf".to_owned(),
                path_arg(&out_gtf.path),
                "--max-extend".to_owned(),
                max_extend.to_string(),
                "--evidence-gap".to_owned(),
                evidence_gap.to_string(),
                "--min-umis".to_owned(),
                min_umis.to_string(),
                "--min-cells".to_owned(),
                min_cells.to_string(),
                "--site-gap".to_owned(),
                site_gap.to_string(),
                "--min-extension".to_owned(),
                min_extension.to_string(),
            ]);
            outputs.push(out_gtf);
            if let Some(report) = report {
                let report = resolver.typed_output(
                    report,
                    ResolvedOutputKind::File,
                    "report",
                    ResourceKind::File,
                )?;
                explanation.push(format!("report -> {}", report.path.display()));
                args.extend(["--report".to_owned(), path_arg(&report.path)]);
                outputs.push(report);
            }
            if let Some(genome) = genome {
                let genome_path = resolver.resource(genome, &[ResourceKind::Genome])?;
                explain_resource(&mut explanation, "genome", genome, &genome_path);
                args.extend(["--genome".to_owned(), path_arg(&genome_path)]);
            }
            if *clip_any_strand {
                args.push("--clip-any-strand".to_owned());
            }
            uniform_io = resolve_uniform_io(
                uniform_report,
                ResolvedUniformIoKind::Report,
                resolver,
                &mut args,
                &mut outputs,
                &mut explanation,
            )?;
        }
    }

    assign_staging_paths(step.id(), &mut args, stdout.as_deref(), &mut outputs)?;
    validate_compiled_command(&args)?;
    let output_schema_ids = output_schema_ids(step);
    for schema in &output_schema_ids {
        explanation.push(format!("output schema -> {schema}"));
    }
    explanation.push(format!("command -> {}", display_command(&args)));
    let (input_resources, step_inputs, embedded_resources) = resolver.finish_step();
    let io_estimate = estimate_step_io(
        step,
        resolver,
        &input_resources,
        &step_inputs,
        &embedded_resources,
        &prepared_inputs,
    )?;
    explanation.push(explain_io_estimate(&io_estimate));
    Ok(ResolvedStep {
        id: step.id().to_owned(),
        kind: step.kind().to_owned(),
        args,
        stdout,
        uniform_io,
        outputs,
        input_resources,
        step_inputs,
        embedded_resources,
        prepared_inputs,
        biological_intent,
        annotation_inputs,
        annotation_comparison,
        output_schema_ids,
        io_estimate,
        explanation,
    })
}

fn resolve_locus_or_feature(
    step_kind: &str,
    locus: &Option<String>,
    feature: &Option<FeatureRequest>,
    annotation: &Option<String>,
    resolver: &mut Resolver<'_>,
    explanation: &mut Vec<String>,
) -> Result<(String, Option<ResolvedBiologicalIntent>)> {
    match (locus, feature) {
        (Some(locus), None) => {
            validate_locus(locus)?;
            explanation.push(format!(
                "locus -> {locus} (caller-supplied, 0-based half-open)"
            ));
            Ok((locus.clone(), None))
        }
        (None, Some(feature)) => {
            let annotation_name = annotation.as_ref().with_context(|| {
                format!("{step_kind} feature resolution requires a named annotation resource")
            })?;
            let identifier = feature.identifier();
            let query = IdentifierQuery::parse(identifier)
                .with_context(|| format!("invalid biological feature `{identifier}`"))?;
            let input = resolver.annotation_resolution(annotation_name)?;
            if let Some(expected) = feature.expected_assembly() {
                validate_feature_expectation(expected, "feature assembly")?;
                if expected != input.metadata.assembly {
                    bail!(
                        "feature expects assembly `{expected}`, but annotation resource `{annotation_name}` is registered as `{}`",
                        input.metadata.assembly
                    );
                }
            }
            if let Some(expected) = feature.expected_annotation() {
                validate_feature_expectation(expected, "feature annotation")?;
                if expected != input.metadata.annotation {
                    bail!(
                        "feature expects annotation `{expected}`, but annotation resource `{annotation_name}` is registered as `{}`",
                        input.metadata.annotation
                    );
                }
            }
            if input.source_identity.scheme != "full-file-blake3-v1" {
                bail!(
                    "annotation resolution requires a full-file BLAKE3 identity, found {}",
                    input.source_identity.scheme
                );
            }
            let expected_digest = format!("blake3:{}", input.source_identity.digest);
            let identity = AnnotationIdentity::new(
                input.metadata.assembly.clone(),
                input.metadata.annotation.clone(),
            )?
            .with_digest(expected_digest)?;
            let intent_resolver = IntentResolver::from_path(&input.source_path, identity)
                .with_context(|| {
                    format!(
                        "loading annotation `{annotation_name}` for biological feature resolution"
                    )
                })?;
            let resolved = intent_resolver.resolve(&query).with_context(|| {
                format!(
                    "resolving biological feature `{identifier}` against annotation `{annotation_name}`"
                )
            })?;
            if resolved.loci.len() != 1 {
                bail!(
                    "biological feature `{identifier}` resolves to {} distinct loci; use a stable identifier with one locus or supply an explicit locus",
                    resolved.loci.len()
                );
            }
            let genomic = &resolved.loci[0];
            let explicit_locus = format!("{}:{}-{}", genomic.contig, genomic.start, genomic.end);
            validate_locus(&explicit_locus)?;
            let annotation_digest = resolved
                .identity
                .digest
                .clone()
                .context("annotation resolver did not return its bound content digest")?;
            explanation.push(format!(
                "feature `{identifier}` -> {} {}:{}-{} ({:?}); {} / {}; {}",
                resolved.stable_id,
                genomic.contig,
                genomic.start,
                genomic.end,
                genomic.strand,
                resolved.identity.assembly,
                resolved.identity.annotation,
                annotation_digest
            ));
            Ok((
                explicit_locus.clone(),
                Some(ResolvedBiologicalIntent {
                    requested: identifier.to_owned(),
                    resolved_kind: resolved.kind,
                    stable_id: resolved.stable_id,
                    display_name: resolved.display_name,
                    matched_by: resolved.matched_by,
                    gene_ids: resolved.gene_ids,
                    transcript_ids: resolved.transcript_ids,
                    annotation_resource: annotation_name.clone(),
                    annotation_path: input.command_path,
                    assembly: resolved.identity.assembly,
                    annotation: resolved.identity.annotation,
                    annotation_digest,
                    contig: genomic.contig.clone(),
                    start: genomic.start,
                    end: genomic.end,
                    strand: genomic.strand,
                    locus: explicit_locus,
                    compatibility: Vec::new(),
                }),
            ))
        }
        (Some(_), Some(_)) => {
            bail!("{step_kind} must specify exactly one of `locus` or `feature`, not both")
        }
        (None, None) => bail!("{step_kind} must specify exactly one of `locus` or `feature`"),
    }
}

/// Resolve one project-registered biological feature without constructing or writing a plan.
/// Explorer and other read-only interfaces use the same fail-closed path as `plan check`.
pub fn resolve_project_feature(
    project: &ProjectContext,
    annotation: &str,
    feature: &FeatureRequest,
) -> Result<ResolvedBiologicalIntent> {
    let mut resolver = Resolver::new(project, BTreeSet::new());
    resolver.begin_step("feature-lookup");
    let (_, intent) = resolve_locus_or_feature(
        "feature lookup",
        &None,
        &Some(feature.clone()),
        &Some(annotation.to_owned()),
        &mut resolver,
        &mut Vec::new(),
    )?;
    intent.context("feature lookup did not produce a biological intent")
}

fn validate_feature_expectation(value: &str, label: &str) -> Result<()> {
    if value.trim() != value
        || value.is_empty()
        || value.len() > 200
        || value.chars().any(char::is_control)
    {
        bail!("{label} must be 1-200 printable characters without surrounding whitespace");
    }
    Ok(())
}

fn output_schema_ids(step: &PlanStep) -> Vec<String> {
    let one = |schema: &str| vec![schema.to_owned()];
    let many = |schemas: &[&str]| schemas.iter().map(|schema| (*schema).to_owned()).collect();
    match step {
        PlanStep::InspectArchive { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.archive.inspect-report.v1",
                    "gravlax.archive.section-accounting.v1",
                ])
            } else {
                one("gravlax.archive.identity.v1")
            }
        }
        PlanStep::QueryRegion { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.query.region.result.v1",
                    "gravlax.query.region.counts.v1",
                ])
            } else {
                one("gravlax.query.region.v1")
            }
        }
        PlanStep::QueryJunction {
            scope,
            uniform_output,
            ..
        } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.query.junction.result.v1",
                    "gravlax.query.junction.counts.v1",
                ])
            } else {
                one(
                    if scope.cells.is_some()
                        || scope.groups.is_some()
                        || scope.aggregation != QueryAggregation::Auto
                    {
                        "gravlax.query.junction.v2"
                    } else {
                        "gravlax.query.junction.v1"
                    },
                )
            }
        }
        PlanStep::QueryJunctions {
            scope,
            with_cells,
            min_cells,
            uniform_output,
            ..
        } => {
            if uniform_output.is_some() {
                let mut schemas = many(&[
                    "gravlax.query.junctions.result.v1",
                    "gravlax.query.junctions.junctions.v1",
                ]);
                if *with_cells
                    || *min_cells != 0
                    || scope.cells.is_some()
                    || scope.groups.is_some()
                    || scope.aggregation != QueryAggregation::Auto
                {
                    schemas.push("gravlax.query.junctions.counts.v1".to_owned());
                }
                schemas
            } else {
                one(
                    if scope.cells.is_some()
                        || scope.groups.is_some()
                        || scope.aggregation != QueryAggregation::Auto
                    {
                        "gravlax.query.junctions.v2"
                    } else {
                        "gravlax.query.junctions.v1"
                    },
                )
            }
        }
        PlanStep::QueryJset { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.query.jset.result.v1",
                    "gravlax.query.jset.junctions.v1",
                    "gravlax.query.jset.counts.v1",
                ])
            } else {
                one("gravlax.query.jset.v1")
            }
        }
        PlanStep::QueryEvents { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.query.events.result.v1",
                    "gravlax.query.events.events.v1",
                    "gravlax.query.events.components.v1",
                    "gravlax.query.events.counts.v1",
                ])
            } else {
                one("gravlax.query.events.v1")
            }
        }
        PlanStep::QueryApa {
            tsv,
            groups,
            uniform_output,
            ..
        } => {
            if uniform_output.is_some() {
                let mut schemas =
                    many(&["gravlax.query.apa.result.v1", "gravlax.query.apa.sites.v1"]);
                if groups.is_some() {
                    schemas.extend(many(&[
                        "gravlax.query.apa.group-counts.v1",
                        "gravlax.query.apa.group-test.v1",
                    ]));
                }
                schemas
            } else {
                one(if *tsv {
                    "gravlax.query.apa.tsv.v1"
                } else {
                    "gravlax.query.apa.summary.v1"
                })
            }
        }
        PlanStep::CompareAnnotations { format, table, .. } => match (format, table) {
            (OutputFormat::Tsv, Some(table)) => one(annotation_comparison_table_schema(*table)),
            (OutputFormat::Json, None) => vec![
                "gravlax.annotation.compare.v1".to_owned(),
                "gravlax.annotation.compare.count-deltas.v1".to_owned(),
                "gravlax.annotation.compare.class-transitions.v1".to_owned(),
                "gravlax.annotation.compare.contributing-causes.v1".to_owned(),
                "gravlax.annotation.compare.witnesses.v1".to_owned(),
            ],
            _ => one("gravlax.annotation.compare.v1"),
        },
        PlanStep::QueryTranscriptEcs {
            format,
            table,
            emit_membership,
            ..
        } => match (format, table) {
            (OutputFormat::Tsv, Some(table)) => one(transcript_ec_table_schema(*table)),
            (OutputFormat::Json, None) => {
                let mut schemas = vec![
                    "gravlax.query.transcript-ecs.v1".to_owned(),
                    "gravlax.query.transcript-ecs.catalog.v1".to_owned(),
                    "gravlax.query.transcript-ecs.counts.v1".to_owned(),
                ];
                if *emit_membership {
                    schemas.push("gravlax.query.transcript-ecs.membership.v1".to_owned());
                }
                schemas
            }
            _ => one("gravlax.query.transcript-ecs.v1"),
        },
        PlanStep::FederateJunction { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.federate.junction.result.v1",
                    "gravlax.federate.junction.archives.v1",
                    "gravlax.federate.junction.counts.v1",
                ])
            } else {
                one("gravlax.federate.junction.v1")
            }
        }
        PlanStep::CollectionRegion { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.collection.region.result.v1",
                    "gravlax.collection.region.samples.v1",
                    "gravlax.collection.region.cells.v1",
                ])
            } else {
                one("gravlax.collection.region.v2")
            }
        }
        PlanStep::CollectionJunction { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.collection.junction.result.v1",
                    "gravlax.collection.junction.samples.v1",
                    "gravlax.collection.junction.cells.v1",
                ])
            } else {
                one("gravlax.collection.junction.v2")
            }
        }
        PlanStep::CollectionJset { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.collection.jset.result.v1",
                    "gravlax.collection.jset.requests.v1",
                    "gravlax.collection.jset.samples.v1",
                    "gravlax.collection.jset.cells.v1",
                ])
            } else {
                one("gravlax.collection.jset.v2")
            }
        }
        PlanStep::CohortEvents { uniform_output, .. } => {
            if uniform_output.is_some() {
                many(&[
                    "gravlax.cohort.events.result.v1",
                    "gravlax.cohort.events.samples.v1",
                    "gravlax.cohort.events.events.v1",
                    "gravlax.cohort.events.components.v1",
                    "gravlax.cohort.events.counts.v1",
                ])
            } else {
                one("gravlax.cohort.events.v1")
            }
        }
        PlanStep::CohortSpliceGraph {
            contrast,
            uniform_output,
            ..
        } => {
            if uniform_output.is_some() {
                let mut schemas = many(&[
                    "gravlax.cohort.splice-graph.result.v1",
                    "gravlax.cohort.splice-graph.samples.v1",
                    "gravlax.cohort.splice-graph.nodes.v1",
                    "gravlax.cohort.splice-graph.edges.v1",
                    "gravlax.cohort.splice-graph.paths.v1",
                    "gravlax.cohort.splice-graph.path-counts.v1",
                    "gravlax.cohort.splice-graph.edge-counts.v1",
                ]);
                if contrast.is_some() {
                    schemas.extend(many(&[
                        "gravlax.cohort.splice-graph.tests.v1",
                        "gravlax.cohort.splice-graph.skipped-tests.v1",
                    ]));
                }
                schemas
            } else {
                one("gravlax.cohort.splice-graph.v1")
            }
        }
        PlanStep::CompileAnnotation { uniform_report, .. } => {
            let mut schemas = one("gravlax.annotation.compiled.v2");
            if uniform_report.is_some() {
                schemas.extend(many(&[
                    "gravlax.annotation.compile.result.v1",
                    "gravlax.annotation.compile.artifacts.v1",
                ]));
            }
            schemas
        }
        PlanStep::IngestArchive { uniform_report, .. } => {
            let mut schemas = one("gravlax.archive.v2");
            if uniform_report.is_some() {
                schemas.extend(many(&[
                    "gravlax.archive.ingest-report.v1",
                    "gravlax.archive.section-accounting.v1",
                ]));
            }
            schemas
        }
        PlanStep::ReplayRows {
            velocity,
            uniform_report,
            ..
        } => {
            let mut schemas = many(&[
                "gravlax.replay.gene-matrix.v1",
                "gravlax.replay.mex-artifact.v1",
            ]);
            if *velocity {
                schemas.push("gravlax.replay.velocity-matrices.v1".to_owned());
            }
            if uniform_report.is_some() {
                schemas.extend(many(&[
                    "gravlax.archive.replay-report.v1",
                    "gravlax.archive.replay-artifact-files.v1",
                ]));
            }
            schemas
        }
        PlanStep::ExtendAnnotation {
            report,
            uniform_report,
            ..
        } => {
            let mut schemas = one("gravlax.annotation.gtf.v1");
            if report.is_some() {
                schemas.push("gravlax.annotation.extension-report.v1".to_owned());
            }
            if uniform_report.is_some() {
                schemas.extend(many(&[
                    "gravlax.extend.result.v1",
                    "gravlax.extend.artifacts.v1",
                    "gravlax.extend.genes.v1",
                ]));
            }
            schemas
        }
    }
}

fn annotation_comparison_table_schema(table: AnnotationComparisonTable) -> &'static str {
    match table {
        AnnotationComparisonTable::CountDeltas => "gravlax.annotation.compare.count-deltas.v1",
        AnnotationComparisonTable::ClassTransitions => {
            "gravlax.annotation.compare.class-transitions.v1"
        }
        AnnotationComparisonTable::ContributingCauses => {
            "gravlax.annotation.compare.contributing-causes.v1"
        }
        AnnotationComparisonTable::Witnesses => "gravlax.annotation.compare.witnesses.v1",
    }
}

fn transcript_ec_table_schema(table: TranscriptEcTable) -> &'static str {
    match table {
        TranscriptEcTable::Catalog => "gravlax.query.transcript-ecs.catalog.v1",
        TranscriptEcTable::Counts => "gravlax.query.transcript-ecs.counts.v1",
        TranscriptEcTable::Membership => "gravlax.query.transcript-ecs.membership.v1",
    }
}

fn estimate_step_io(
    step: &PlanStep,
    resolver: &Resolver<'_>,
    input_resources: &[String],
    step_inputs: &[ResolvedStepInput],
    embedded_resources: &[String],
    prepared_inputs: &[ResolvedPreparedInput],
) -> Result<ResolvedIoEstimate> {
    let mut known = BTreeMap::<PathBuf, u64>::new();
    for name in input_resources {
        let resource = resolver
            .resources
            .get(name)
            .with_context(|| format!("internal error estimating resource `{name}`"))?;
        known.insert(resource.path.clone(), resource.bytes);
    }
    for name in embedded_resources {
        let resource = resolver
            .embedded_resources
            .get(name)
            .with_context(|| format!("internal error estimating embedded resource `{name}`"))?;
        known.insert(resource.path.clone(), resource.bytes);
    }
    for input in prepared_inputs {
        known.insert(input.path.clone(), input.bytes);
    }
    let known_selected_input_bytes = known.values().try_fold(0u64, |total, bytes| {
        total
            .checked_add(*bytes)
            .context("selected input byte estimate overflow")
    })?;
    let route_dependent = matches!(
        step,
        PlanStep::CollectionRegion { .. }
            | PlanStep::CollectionJunction { .. }
            | PlanStep::CollectionJset { .. }
    );
    let unavailable_reason = if route_dependent {
        "collection member archives are selected dynamically"
    } else if !step_inputs.is_empty() {
        "one or more inputs are produced by earlier steps and have no size at check time"
    } else {
        "command read multiplicity, caching, decompression, and operating-system I/O are not modeled"
    };
    let read_bytes_upper_bound = None;
    let bound = IoEstimateBound::KnownInputsOnly;
    let note = format!(
        "exact known selected file sizes only; no execution read upper bound is claimed because {unavailable_reason}"
    );
    Ok(ResolvedIoEstimate {
        known_selected_input_bytes,
        known_selected_input_files: known.len(),
        unknown_prior_step_outputs: step_inputs.len(),
        read_bytes_lower_bound: 0,
        read_bytes_upper_bound,
        bound,
        note,
    })
}

fn explain_io_estimate(estimate: &ResolvedIoEstimate) -> String {
    let upper = estimate
        .read_bytes_upper_bound
        .map_or_else(|| "unavailable".to_owned(), |bytes| bytes.to_string());
    format!(
        "selected-input estimate -> {} known byte(s) across {} file(s), {} prior-step input(s); read bound 0..{upper} byte(s): {}",
        estimate.known_selected_input_bytes,
        estimate.known_selected_input_files,
        estimate.unknown_prior_step_outputs,
        estimate.note
    )
}

fn assign_staging_paths(
    step_id: &str,
    args: &mut [String],
    stdout: Option<&Path>,
    outputs: &mut [ResolvedOutput],
) -> Result<()> {
    for output in outputs {
        output.staging_path = staging_output_path(&output.path, output.kind, step_id)?;
        if stdout == Some(output.path.as_path()) {
            continue;
        }
        let final_arg = path_arg(&output.path);
        let staged_arg = path_arg(&output.staging_path);
        let mut replacements = 0usize;
        for argument in args.iter_mut() {
            if *argument == final_arg {
                *argument = staged_arg.clone();
                replacements += 1;
            }
        }
        if replacements != 1 {
            bail!(
                "internal plan compiler error: expected one output argument for {}, found {replacements}",
                output.path.display()
            );
        }
    }
    Ok(())
}

fn staging_output_path(output: &Path, kind: ResolvedOutputKind, step_id: &str) -> Result<PathBuf> {
    let file_name = output
        .file_name()
        .and_then(|name| name.to_str())
        .context("resolved output file name is not valid UTF-8")?;
    let staging_name = if kind == ResolvedOutputKind::File {
        match (
            output.file_stem().and_then(|name| name.to_str()),
            output.extension().and_then(|name| name.to_str()),
        ) {
            (Some(stem), Some(extension)) => {
                format!(".{stem}.aie-stage-{step_id}.{extension}")
            }
            _ => format!(".{file_name}.aie-stage-{step_id}"),
        }
    } else {
        format!(".{file_name}.aie-stage-{step_id}")
    };
    Ok(output.with_file_name(staging_name))
}

fn validate_compiled_command(args: &[String]) -> Result<()> {
    let argv = std::iter::once("aie".to_owned()).chain(args.iter().cloned());
    crate::Cli::try_parse_from(argv)
        .map(|_| ())
        .map_err(|error| anyhow::anyhow!("compiled command does not match this aie CLI: {error}"))
}

fn resolve_stdout(
    output: &Option<PathBuf>,
    resolver: &Resolver<'_>,
    stdout: &mut Option<PathBuf>,
    outputs: &mut Vec<ResolvedOutput>,
    explanation: &mut Vec<String>,
) -> Result<()> {
    if let Some(output) = output {
        let resolved = resolver.output(output, ResolvedOutputKind::File)?;
        explanation.push(format!("stdout -> {}", resolved.path.display()));
        *stdout = Some(resolved.path.clone());
        outputs.push(resolved);
    } else {
        explanation.push("stdout -> terminal".to_owned());
    }
    Ok(())
}

fn resolve_uniform_io(
    requested: &Option<UniformOutput>,
    kind: ResolvedUniformIoKind,
    resolver: &Resolver<'_>,
    args: &mut Vec<String>,
    outputs: &mut Vec<ResolvedOutput>,
    explanation: &mut Vec<String>,
) -> Result<Option<ResolvedUniformIo>> {
    let Some(requested) = requested else {
        return Ok(None);
    };
    let (format_flag, output_flag, output_name, label) = match kind {
        ResolvedUniformIoKind::Result => ("--format", "--output", "result", "uniform result"),
        ResolvedUniformIoKind::Report => (
            "--report-format",
            "--report-output",
            "uniform-report",
            "uniform report",
        ),
    };
    args.extend([format_flag.to_owned(), requested.format.as_arg().to_owned()]);
    let (publication, output) = if let Some(output) = &requested.output {
        let resolved = resolver.typed_output(
            output,
            ResolvedOutputKind::File,
            output_name,
            ResourceKind::File,
        )?;
        explanation.push(format!("{label} -> {}", resolved.path.display()));
        args.extend([output_flag.to_owned(), path_arg(&resolved.path)]);
        let path = resolved.path.clone();
        outputs.push(resolved);
        (ResolvedUniformPublication::AtomicNoClobberFile, Some(path))
    } else {
        explanation.push(format!("{label} -> stdout"));
        (ResolvedUniformPublication::Stdout, None)
    };
    Ok(Some(ResolvedUniformIo {
        kind,
        format: requested.format,
        publication,
        output,
    }))
}

fn reject_legacy_result_options(
    requested: &Option<UniformOutput>,
    legacy_format: bool,
    legacy_output: bool,
    kind: &str,
) -> Result<()> {
    if requested.is_some() && (legacy_format || legacy_output) {
        bail!("{kind} uniform_output cannot be combined with legacy format or output fields");
    }
    Ok(())
}

fn resolve_capability_output(
    output: &Option<PathBuf>,
    resolver: &Resolver<'_>,
    args: &mut Vec<String>,
    outputs: &mut Vec<ResolvedOutput>,
    explanation: &mut Vec<String>,
) -> Result<()> {
    if let Some(output) = output {
        let resolved = resolver.typed_output(
            output,
            ResolvedOutputKind::File,
            "result",
            ResourceKind::File,
        )?;
        explanation.push(format!("result -> {}", resolved.path.display()));
        args.extend(["--output".to_owned(), path_arg(&resolved.path)]);
        outputs.push(resolved);
    } else {
        explanation.push("result -> terminal".to_owned());
    }
    Ok(())
}

fn resolve_query_scope(
    scope: &QueryScope,
    resolver: &mut Resolver<'_>,
    args: &mut Vec<String>,
    explanation: &mut Vec<String>,
) -> Result<()> {
    if scope.cells.is_some() && scope.groups.is_some() {
        bail!("query scope cannot specify both cells and groups");
    }
    if scope.aggregation == QueryAggregation::Group && scope.groups.is_none() {
        bail!("group aggregation requires a named groups resource");
    }
    if let Some(cells) = &scope.cells {
        let path = resolver.resource(cells, &[ResourceKind::Cells, ResourceKind::Metadata])?;
        explain_resource(explanation, "cells", cells, &path);
        args.extend(["--cells".to_owned(), path_arg(&path)]);
    }
    if let Some(groups) = &scope.groups {
        let path = resolver.resource(groups, &[ResourceKind::Groups, ResourceKind::Metadata])?;
        explain_resource(explanation, "groups", groups, &path);
        args.extend(["--groups".to_owned(), path_arg(&path)]);
    }
    let aggregation = match scope.aggregation {
        QueryAggregation::Auto => None,
        QueryAggregation::Cell => Some("cell"),
        QueryAggregation::Group => Some("group"),
        QueryAggregation::Bulk => Some("bulk"),
    };
    if let Some(aggregation) = aggregation {
        args.extend(["--agg".to_owned(), aggregation.to_owned()]);
        explanation.push(format!("aggregation -> {aggregation}"));
    }
    Ok(())
}

fn validate_junction_sets(include: &[String], exclude: &[String]) -> Result<()> {
    if include.is_empty() || exclude.is_empty() {
        bail!("junction-set steps require nonempty include and exclude lists");
    }
    let mut included = BTreeSet::new();
    for locus in include {
        validate_locus(locus)?;
        if !included.insert(locus) {
            bail!("duplicate inclusion junction `{locus}`");
        }
    }
    let mut excluded = BTreeSet::new();
    for locus in exclude {
        validate_locus(locus)?;
        if included.contains(locus) {
            bail!("junction `{locus}` appears in both include and exclude lists");
        }
        if !excluded.insert(locus) {
            bail!("duplicate exclusion junction `{locus}`");
        }
    }
    Ok(())
}

fn push_repeated(args: &mut Vec<String>, flag: &str, values: &[String]) {
    for value in values {
        args.extend([flag.to_owned(), value.clone()]);
    }
}

fn push_event_types(args: &mut Vec<String>, event_types: &[EventType]) {
    for event_type in event_types {
        let value = match event_type {
            EventType::AltAcceptor => "alt-acceptor",
            EventType::AltDonor => "alt-donor",
            EventType::Cassette => "cassette",
        };
        args.extend(["--event-type".to_owned(), value.to_owned()]);
    }
}

fn require_human_or_json(format: OutputFormat, kind: &str) -> Result<()> {
    if format == OutputFormat::Tsv {
        bail!("{kind} supports human or json format, not tsv");
    }
    Ok(())
}

fn push_collection_flags(
    args: &mut Vec<String>,
    explain_routing: bool,
    verify_content: bool,
    format: OutputFormat,
) {
    if explain_routing {
        args.push("--explain".to_owned());
    }
    if verify_content {
        args.push("--verify-content".to_owned());
    }
    if format == OutputFormat::Json {
        args.push("--json".to_owned());
    }
}

fn validate_cohort_id(value: &str, label: &str) -> Result<()> {
    let valid = !value.is_empty()
        && value.len() <= 64
        && value.chars().all(|character| {
            character.is_ascii_alphanumeric() || matches!(character, '.' | '_' | '-')
        });
    if !valid {
        bail!("invalid {label} ID `{value}`; use 1-64 ASCII letters, digits, '.', '_' or '-'");
    }
    Ok(())
}

fn validate_contrast(contrast: &str) -> Result<(&str, &str)> {
    let (left, right) = contrast
        .split_once(':')
        .with_context(|| format!("contrast must be CONDITION_A:CONDITION_B, got `{contrast}`"))?;
    validate_cohort_id(left, "contrast condition")?;
    validate_cohort_id(right, "contrast condition")?;
    if left == right {
        bail!("contrast conditions must be distinct");
    }
    Ok((left, right))
}

/// A named design remains human-editable, but plan checking still verifies every embedded path.
/// Internal designs are portable and therefore project-contained; an explicitly external design
/// may refer to external read-only archives and cell lists.
struct CanonicalDesign {
    content: String,
    sample_count: usize,
    condition_counts: BTreeMap<String, usize>,
}

fn resolve_design_paths(
    resolver: &mut Resolver<'_>,
    design_name: &str,
    design: &Path,
    allow_external: bool,
) -> Result<CanonicalDesign> {
    let project_root = resolver.project.root.clone();
    let text = fs::read_to_string(design)
        .with_context(|| format!("reading cohort design {}", design.display()))?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .map(|line| line.trim_end_matches('\r'))
        .context("cohort design is empty")?;
    if header != "sample\tcondition\tarchive\tcells" {
        bail!("cohort design header must be exactly: sample<TAB>condition<TAB>archive<TAB>cells");
    }
    let mut canonical_design = String::from("sample\tcondition\tarchive\tcells\n");
    let base = design.parent().context("cohort design has no parent")?;
    let mut samples = BTreeSet::new();
    let mut archive_paths = BTreeSet::new();
    let mut condition_counts = BTreeMap::<String, usize>::new();
    for (index, raw_line) in lines.enumerate() {
        let line_no = index + 2;
        let line = raw_line.trim_end_matches('\r');
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 4 || fields.iter().any(|field| field.is_empty()) {
            bail!("cohort design line {line_no} must have four nonempty tab-separated fields");
        }
        validate_cohort_id(fields[0], "design sample")?;
        validate_cohort_id(fields[1], "design condition")?;
        if !samples.insert(fields[0]) {
            bail!("cohort design contains duplicate sample `{}`", fields[0]);
        }
        *condition_counts.entry(fields[1].to_owned()).or_default() += 1;
        let mut canonical_paths = Vec::with_capacity(2);
        for (role, value) in [("archive", fields[2]), ("cells", fields[3])] {
            if role == "cells" && value == "." {
                canonical_paths.push(".".to_owned());
                continue;
            }
            let path = Path::new(value);
            let path = if path.is_absolute() {
                path.to_owned()
            } else {
                base.join(path)
            };
            let resolved = fs::canonicalize(&path).with_context(|| {
                format!(
                    "resolving cohort design {role} on line {line_no}: {}",
                    path.display()
                )
            })?;
            if !allow_external && !resolved.starts_with(&project_root) {
                bail!(
                    "cohort design {role} on line {line_no} escapes project {}; register the design itself with --external to allow external read-only paths: {}",
                    project_root.display(),
                    resolved.display()
                );
            }
            if !resolved.is_file() {
                bail!(
                    "cohort design {role} on line {line_no} is not a file: {}",
                    resolved.display()
                );
            }
            if role == "archive" && !archive_paths.insert(resolved.clone()) {
                bail!(
                    "cohort design reuses resolved archive {} as more than one biological sample",
                    resolved.display()
                );
            }
            require_utf8_path(
                &resolved,
                &format!("cohort design {role} on line {line_no}"),
            )?;
            let canonical_path = path_arg(&resolved);
            require_tsv_field(
                &canonical_path,
                &format!("resolved cohort design {role} path on line {line_no}"),
            )?;
            let kind = if role == "archive" {
                ResourceKind::Archive
            } else {
                ResourceKind::Cells
            };
            let key = format!("{design_name}/{}/{role}", fields[0]);
            resolver.embedded_resource(key, design_name, fields[0], role, value, kind, resolved)?;
            canonical_paths.push(canonical_path);
        }
        canonical_design.push_str(fields[0]);
        canonical_design.push('\t');
        canonical_design.push_str(fields[1]);
        canonical_design.push('\t');
        canonical_design.push_str(&canonical_paths[0]);
        canonical_design.push('\t');
        canonical_design.push_str(&canonical_paths[1]);
        canonical_design.push('\n');
    }
    if samples.len() < 2 {
        bail!("cohort design requires at least two samples");
    }
    Ok(CanonicalDesign {
        content: canonical_design,
        sample_count: samples.len(),
        condition_counts,
    })
}

fn resolved_design_input(
    resolver: &Resolver<'_>,
    design_name: &str,
    content: String,
) -> Result<ResolvedPreparedInput> {
    let digest = blake3::hash(content.as_bytes()).to_hex().to_string();
    let path = resolver
        .project
        .root
        .join(".aie/resolved-inputs")
        .join(format!("{design_name}-{}.tsv", &digest[..12]));
    require_utf8_path(&path, "canonical cohort design path")?;
    Ok(ResolvedPreparedInput {
        path,
        bytes: content.len() as u64,
        identity: ResolvedIdentity {
            scheme: "full-file-blake3-v1".to_owned(),
            digest,
        },
        content,
    })
}

fn explain_resource(lines: &mut Vec<String>, role: &str, name: &str, path: &Path) {
    lines.push(format!("{role} `{name}` -> {}", path.display()));
}

fn push_format(args: &mut Vec<String>, format: OutputFormat) {
    match format {
        OutputFormat::Human => {}
        OutputFormat::Tsv => args.push("--tsv".to_owned()),
        OutputFormat::Json => args.push("--json".to_owned()),
    }
}

fn push_named_format(args: &mut Vec<String>, format: OutputFormat) {
    let value = match format {
        OutputFormat::Human => "text",
        OutputFormat::Tsv => "tsv",
        OutputFormat::Json => "json",
    };
    args.extend(["--format".to_owned(), value.to_owned()]);
}

fn validate_bundle_format(format: OutputFormat, table: Option<&str>, kind: &str) -> Result<()> {
    match (format, table) {
        (OutputFormat::Tsv, None) => {
            bail!("{kind} format tsv requires an explicit table selection")
        }
        (OutputFormat::Human | OutputFormat::Json, Some(_)) => {
            bail!("{kind} table selection is only valid with format tsv")
        }
        _ => Ok(()),
    }
}

fn comparison_gene_key_arg(value: ComparisonGeneKey) -> &'static str {
    match value {
        ComparisonGeneKey::Unversioned => "unversioned",
        ComparisonGeneKey::Exact => "exact",
    }
}

fn annotation_comparison_table_arg(value: AnnotationComparisonTable) -> &'static str {
    match value {
        AnnotationComparisonTable::CountDeltas => "count-deltas",
        AnnotationComparisonTable::ClassTransitions => "class-transitions",
        AnnotationComparisonTable::ContributingCauses => "contributing-causes",
        AnnotationComparisonTable::Witnesses => "witnesses",
    }
}

fn transcript_ec_table_arg(value: TranscriptEcTable) -> &'static str {
    match value {
        TranscriptEcTable::Catalog => "catalog",
        TranscriptEcTable::Counts => "counts",
        TranscriptEcTable::Membership => "membership",
    }
}

fn push_expected_annotation_digest(
    args: &mut Vec<String>,
    flag: &str,
    identity: Option<&ResolvedIdentity>,
) -> Result<()> {
    let Some(identity) = identity else {
        return Ok(());
    };
    if identity.scheme != "full-file-blake3-v1" {
        bail!(
            "annotation command identity must use full-file-blake3-v1, found {}",
            identity.scheme
        );
    }
    args.extend([flag.to_owned(), format!("blake3:{}", identity.digest)]);
    Ok(())
}

fn validate_comparison_annotation_source(
    path: &Path,
    gene_key: ComparisonGeneKey,
    label: &str,
) -> Result<()> {
    let annotation = anno::Annotation::from_path(path)
        .with_context(|| format!("loading {label} from {}", path.display()))?;
    if gene_key == ComparisonGeneKey::Unversioned {
        let mut normalized = BTreeMap::<String, &str>::new();
        for gene_id in &annotation.gene_ids {
            let key = unversioned_gene_id(gene_id).to_owned();
            if let Some(previous) = normalized.insert(key.clone(), gene_id) {
                bail!(
                    "{label} gene-id normalization collision: `{previous}` and `{gene_id}` both map to `{key}`"
                );
            }
        }
    }
    Ok(())
}

fn unversioned_gene_id(gene_id: &str) -> &str {
    let Some((prefix, suffix)) = gene_id.rsplit_once('.') else {
        return gene_id;
    };
    if !suffix.is_empty() && suffix.bytes().all(|byte| byte.is_ascii_digit()) {
        prefix
    } else {
        gene_id
    }
}

fn validate_transcript_ec_annotation_source(path: &Path) -> Result<()> {
    let annotation = anno::Annotation::from_path(path).with_context(|| {
        format!(
            "loading transcript-equivalence annotation from {}",
            path.display()
        )
    })?;
    if annotation.transcript_ids.len() != annotation.transcripts.len() {
        bail!(
            "transcript identifier metadata has {} entries for {} transcripts",
            annotation.transcript_ids.len(),
            annotation.transcripts.len()
        );
    }
    if !annotation.transcripts.is_empty() && annotation.transcript_ids.iter().all(Option::is_none) {
        bail!(
            "transcript equivalence classes require stable transcript IDs, but this annotation has none; recompile the source GTF to AIC v2 or use the GTF directly"
        );
    }
    let mut seen = BTreeSet::new();
    for (index, identifier) in annotation.transcript_ids.iter().enumerate() {
        let identifier = identifier.as_deref().with_context(|| {
            format!(
                "transcript equivalence classes require a stable ID for annotation transcript {index}"
            )
        })?;
        if identifier.is_empty() {
            bail!("annotation transcript {index} has an empty transcript ID");
        }
        if !seen.insert(identifier) {
            bail!("annotation contains duplicate transcript ID `{identifier}`");
        }
    }
    Ok(())
}

fn solo_strand_arg(strand: SoloStrand) -> &'static str {
    match strand {
        SoloStrand::Forward => "forward",
        SoloStrand::Reverse => "reverse",
        SoloStrand::Unstranded => "unstranded",
    }
}

fn require_utf8_path(path: &Path, label: &str) -> Result<()> {
    if path.to_str().is_none() {
        bail!("{label} must use a valid UTF-8 path so it can be preserved in a resolved plan");
    }
    Ok(())
}

fn path_arg(path: &Path) -> String {
    path.to_str()
        .expect("resolved plan paths are validated as UTF-8")
        .to_owned()
}

fn require_tsv_field(value: &str, label: &str) -> Result<()> {
    if value.chars().any(char::is_control) {
        bail!("{label} contains a control character and cannot be represented safely in TSV");
    }
    Ok(())
}

fn validate_output_collisions(steps: &[ResolvedStep]) -> Result<()> {
    let mut outputs = Vec::<(&str, &'static str, &Path, ResolvedOutputKind)>::new();
    for step in steps {
        for output in &step.outputs {
            for (label, path) in [
                ("output", output.path.as_path()),
                ("staging output", output.staging_path.as_path()),
            ] {
                for (other_step, other_label, other_path, other_kind) in &outputs {
                    let same = path == *other_path;
                    let nested_in_other = *other_kind == ResolvedOutputKind::Directory
                        && path.starts_with(other_path);
                    let contains_other = output.kind == ResolvedOutputKind::Directory
                        && other_path.starts_with(path);
                    if same || nested_in_other || contains_other {
                        bail!(
                            "{label} collision between step `{}` and its {other_label} in step `{}` at {}",
                            step.id,
                            other_step,
                            path.display()
                        );
                    }
                }
                outputs.push((step.id.as_str(), label, path, output.kind));
            }
        }
    }
    Ok(())
}

fn validate_plan_name(name: &str) -> Result<()> {
    if name.trim().is_empty() || name.len() > 200 || name.chars().any(char::is_control) {
        bail!("plan name must be 1-200 printable characters");
    }
    Ok(())
}

fn validate_step_id(id: &str) -> Result<()> {
    validate_resource_id(id).map_err(|_| {
        anyhow::anyhow!(
            "invalid step id `{id}`; use 1-64 ASCII letters, digits, '.', '_' or '-', starting with a letter or digit"
        )
    })
}

fn validate_locus(locus: &str) -> Result<()> {
    if locus.chars().any(char::is_whitespace) {
        bail!("locus must not contain whitespace: `{locus}`");
    }
    let (chromosome, range) = locus
        .rsplit_once(':')
        .with_context(|| format!("invalid locus `{locus}`; expected chromosome:start-end"))?;
    let (start, end) = range
        .split_once('-')
        .with_context(|| format!("invalid locus `{locus}`; expected chromosome:start-end"))?;
    let start = start
        .parse::<u32>()
        .with_context(|| format!("invalid locus start in `{locus}`"))?;
    let end = end
        .parse::<u32>()
        .with_context(|| format!("invalid locus end in `{locus}`"))?;
    if chromosome.is_empty() || start >= end {
        bail!("invalid locus `{locus}`; chromosome must be nonempty and start < end");
    }
    Ok(())
}

fn require_extension(path: &Path, extension: &str, label: &str) -> Result<()> {
    if path.extension().and_then(|value| value.to_str()) != Some(extension) {
        bail!("{label} must end in .{extension}: {}", path.display());
    }
    Ok(())
}

fn explain_plan(plan: &ResolvedPlan) -> String {
    let mut explanation = String::new();
    explanation.push_str(&format!("plan: {}\n", plan.name));
    explanation.push_str(&format!("source: {}\n", plan.source_path.display()));
    explanation.push_str(&format!(
        "project: {} ({})\n",
        plan.project_name,
        plan.project_root.display()
    ));
    for (index, step) in plan.steps.iter().enumerate() {
        explanation.push_str(&format!(
            "step {}: {} ({})\n",
            index + 1,
            step.id,
            step.kind
        ));
        for line in &step.explanation {
            explanation.push_str(&format!("  {line}\n"));
        }
    }
    explanation
}

fn display_command(args: &[String]) -> String {
    let mut rendered = String::from("aie");
    for argument in args {
        rendered.push(' ');
        rendered.push_str(&display_argument(argument));
    }
    rendered
}

fn display_argument(argument: &str) -> String {
    if !argument.is_empty()
        && argument
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '/' | '.' | '_' | '-' | ':' | '='))
    {
        argument.to_owned()
    } else {
        format!("{argument:?}")
    }
}

fn resolved_plan_bytes(plan: &ResolvedPlan) -> Result<Vec<u8>> {
    let mut bytes = serde_json::to_vec_pretty(plan)?;
    bytes.push(b'\n');
    Ok(bytes)
}

pub fn resolved_plan_digest(plan: &ResolvedPlan) -> Result<String> {
    Ok(blake3::hash(&resolved_plan_bytes(plan)?)
        .to_hex()
        .to_string())
}

pub fn resolved_step_digest(step: &ResolvedStep) -> Result<String> {
    Ok(blake3::hash(&serde_json::to_vec(step)?)
        .to_hex()
        .to_string())
}

fn ensure_outputs_absent(plan: &ResolvedPlan) -> Result<()> {
    for step in &plan.steps {
        for output in &step.outputs {
            for (label, path) in [
                ("output", &output.path),
                ("stale staging output", &output.staging_path),
            ] {
                match fs::symlink_metadata(path) {
                    Ok(_) => bail!(
                        "refusing to overwrite {label} from step `{}`: {}; use --resume only for outputs with an exact completion record",
                        step.id,
                        path.display()
                    ),
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                    Err(error) => {
                        return Err(error)
                            .with_context(|| format!("inspecting {label} {}", path.display()));
                    }
                }
            }
        }
    }
    Ok(())
}

fn resume_decisions(plan: &ResolvedPlan, plan_digest: &str) -> Result<Vec<ResumeDecision>> {
    let mut decisions = Vec::with_capacity(plan.steps.len());
    for step in &plan.steps {
        for output in &step.outputs {
            match fs::symlink_metadata(&output.staging_path) {
                Ok(_) => bail!(
                    "cannot resume step `{}` while stale staging output exists: {}",
                    step.id,
                    output.staging_path.display()
                ),
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                Err(error) => {
                    return Err(error).with_context(|| {
                        format!(
                            "inspecting staging output {}",
                            output.staging_path.display()
                        )
                    });
                }
            }
        }
        if step.outputs.is_empty() {
            decisions.push(ResumeDecision::Run);
            continue;
        }
        let mut present = Vec::with_capacity(step.outputs.len());
        for output in &step.outputs {
            present.push(match fs::symlink_metadata(&output.path) {
                Ok(_) => true,
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => false,
                Err(error) => {
                    return Err(error)
                        .with_context(|| format!("inspecting output {}", output.path.display()));
                }
            });
        }
        if present.iter().all(|value| !value) {
            decisions.push(ResumeDecision::Run);
            continue;
        }
        if !present.iter().all(|value| *value) {
            bail!(
                "cannot resume step `{}`: only some declared outputs exist; remove or restore the step outputs as a set",
                step.id
            );
        }
        verify_step_completion(plan, plan_digest, step)?;
        decisions.push(ResumeDecision::SkipVerified);
    }
    Ok(decisions)
}

fn verify_step_completion(
    plan: &ResolvedPlan,
    plan_digest: &str,
    step: &ResolvedStep,
) -> Result<()> {
    let record_path = completion_record_path(plan, plan_digest, &step.id);
    let record = load_completion_record(&record_path).with_context(|| {
        format!(
            "cannot trust existing outputs for step `{}` without an exact completion record at {}",
            step.id,
            record_path.display()
        )
    })?;
    let expected_step_digest = resolved_step_digest(step)?;
    if record.schema_version != COMPLETION_SCHEMA_VERSION
        || record.resolved_plan_digest != plan_digest
        || record.step_id != step.id
        || record.step_digest != expected_step_digest
    {
        bail!(
            "completion record does not match the resolved plan and step for `{}`",
            step.id
        );
    }
    let observed = completed_outputs(&step.outputs)?;
    if record.outputs != observed {
        bail!(
            "cannot resume step `{}`: one or more output identities differ from its completion record",
            step.id
        );
    }
    let observed_inputs = completed_step_inputs(plan, step)?;
    if record.inputs != observed_inputs {
        bail!(
            "cannot resume step `{}`: one or more upstream step-output identities differ from its completion record",
            step.id
        );
    }
    Ok(())
}

fn explain_resume(plan: &ResolvedPlan, decisions: &[ResumeDecision]) {
    for (step, decision) in plan.steps.iter().zip(decisions) {
        match decision {
            ResumeDecision::Run if step.outputs.is_empty() => println!(
                "resume: step `{}` has no declared artifacts and will run",
                step.id
            ),
            ResumeDecision::Run => {
                println!("resume: step `{}` has no outputs yet and will run", step.id)
            }
            ResumeDecision::SkipVerified => println!(
                "resume: step `{}` has exact verified outputs and will be skipped",
                step.id
            ),
        }
    }
}

pub fn completion_record_path(plan: &ResolvedPlan, plan_digest: &str, step_id: &str) -> PathBuf {
    plan.project_root
        .join(".aie/completions")
        .join(plan_digest)
        .join(format!("{step_id}.json"))
}

fn load_completion_record(path: &Path) -> Result<StepCompletion> {
    let metadata = fs::symlink_metadata(path)
        .with_context(|| format!("reading completion record metadata {}", path.display()))?;
    if !metadata.file_type().is_file() {
        bail!(
            "completion record is not a regular file: {}",
            path.display()
        );
    }
    let bytes =
        fs::read(path).with_context(|| format!("reading completion record {}", path.display()))?;
    let record: StepCompletion = serde_json::from_slice(&bytes)
        .with_context(|| format!("parsing completion record {}", path.display()))?;
    Ok(record)
}

fn completed_outputs(outputs: &[ResolvedOutput]) -> Result<Vec<CompletedOutput>> {
    outputs.iter().map(completed_output).collect()
}

fn resolved_step_input_output<'a>(
    plan: &'a ResolvedPlan,
    input: &ResolvedStepInput,
) -> Result<&'a ResolvedOutput> {
    let producer = plan
        .steps
        .iter()
        .find(|candidate| candidate.id == input.producer_step)
        .with_context(|| {
            format!(
                "step input `{}` names missing producer `{}`",
                input.reference, input.producer_step
            )
        })?;
    let output = producer
        .outputs
        .iter()
        .find(|candidate| candidate.name == input.output_name)
        .with_context(|| {
            format!(
                "step input `{}` names missing output `{}` from producer `{}`",
                input.reference, input.output_name, input.producer_step
            )
        })?;
    if output.path != input.path
        || output.kind != input.kind
        || output.resource_kind != input.resource_kind
    {
        bail!(
            "step input `{}` does not match its producer's resolved output",
            input.reference
        );
    }
    Ok(output)
}

fn completed_step_inputs(
    plan: &ResolvedPlan,
    step: &ResolvedStep,
) -> Result<Vec<CompletedStepInput>> {
    step.step_inputs
        .iter()
        .map(|input| {
            let output = resolved_step_input_output(plan, input)?;
            let completed = completed_output(output)?;
            Ok(CompletedStepInput {
                producer_step: input.producer_step.clone(),
                output_name: input.output_name.clone(),
                path: completed.path,
                kind: completed.kind,
                bytes: completed.bytes,
                identity: completed.identity,
            })
        })
        .collect()
}

fn verify_step_inputs(plan: &ResolvedPlan, plan_digest: &str, step: &ResolvedStep) -> Result<()> {
    for name in &step.input_resources {
        let resource = plan.resources.get(name).with_context(|| {
            format!(
                "resolved step `{}` names missing resource `{name}`",
                step.id
            )
        })?;
        verify_resource_identity(
            &format!("project resource `{name}`"),
            resource.kind,
            &resource.path,
            resource.bytes,
            &resource.identity,
        )?;
    }
    for key in &step.embedded_resources {
        let resource = plan.embedded_resources.get(key).with_context(|| {
            format!(
                "resolved step `{}` names missing embedded resource `{key}`",
                step.id
            )
        })?;
        let owner = plan
            .resources
            .get(&resource.owner_resource)
            .with_context(|| {
                format!(
                    "embedded resource `{key}` names missing owner `{}`",
                    resource.owner_resource
                )
            })?;
        let declared = Path::new(&resource.declared_path);
        let candidate = if declared.is_absolute() {
            declared.to_owned()
        } else {
            owner
                .path
                .parent()
                .context("embedded-resource owner has no parent")?
                .join(declared)
        };
        let resolved = fs::canonicalize(&candidate).with_context(|| {
            format!(
                "re-resolving embedded {} for sample `{}` from {}",
                resource.role,
                resource.sample,
                candidate.display()
            )
        })?;
        if resolved != resource.path {
            bail!(
                "embedded {} for sample `{}` now resolves to {}, not the snapshotted {}",
                resource.role,
                resource.sample,
                resolved.display(),
                resource.path.display()
            );
        }
        if !owner.external && !resolved.starts_with(&plan.project_root) {
            bail!(
                "embedded {} for sample `{}` escaped the project after resolution",
                resource.role,
                resource.sample
            );
        }
        verify_resource_identity(
            &format!(
                "embedded {} for sample `{}`",
                resource.role, resource.sample
            ),
            resource.kind,
            &resource.path,
            resource.bytes,
            &resource.identity,
        )?;
    }
    let mut verified_producers = BTreeSet::new();
    for input in &step.step_inputs {
        resolved_step_input_output(plan, input)?;
        if verified_producers.insert(input.producer_step.as_str()) {
            let producer = plan
                .steps
                .iter()
                .find(|candidate| candidate.id == input.producer_step)
                .expect("resolved_step_input_output found this producer");
            verify_step_completion(plan, plan_digest, producer).with_context(|| {
                format!(
                    "verifying upstream step `{}` before consumer `{}`",
                    input.producer_step, step.id
                )
            })?;
        }
    }
    Ok(())
}

fn verify_resource_identity(
    label: &str,
    kind: ResourceKind,
    path: &Path,
    expected_bytes: u64,
    expected_identity: &ResolvedIdentity,
) -> Result<()> {
    let (observed_bytes, observed_identity) =
        resolve_resource_identity(kind, path).with_context(|| format!("re-verifying {label}"))?;
    if observed_bytes != expected_bytes || &observed_identity != expected_identity {
        bail!("{label} changed after the resolved plan was created");
    }
    Ok(())
}

fn materialize_prepared_inputs(plan: &ResolvedPlan, step: &ResolvedStep) -> Result<()> {
    if step.prepared_inputs.is_empty() {
        return Ok(());
    }
    let directory =
        ensure_project_directory(&plan.project_root, Path::new(".aie/resolved-inputs"))?;
    for prepared in &step.prepared_inputs {
        if prepared.path.parent() != Some(directory.as_path()) {
            bail!(
                "prepared input for step `{}` is outside the resolved-input directory: {}",
                step.id,
                prepared.path.display()
            );
        }
        let content = prepared.content.as_bytes();
        let observed = ResolvedIdentity {
            scheme: "full-file-blake3-v1".to_owned(),
            digest: blake3::hash(content).to_hex().to_string(),
        };
        if prepared.bytes != content.len() as u64 || prepared.identity != observed {
            bail!(
                "prepared input metadata is inconsistent for step `{}`",
                step.id
            );
        }
        install_immutable_file(&prepared.path, content, "prepared input")?;
    }
    Ok(())
}

fn install_immutable_file(path: &Path, bytes: &[u8], label: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            if !metadata.file_type().is_file() {
                bail!("{label} is not a regular file: {}", path.display());
            }
            if fs::read(path)? != bytes {
                bail!("existing {label} has different content: {}", path.display());
            }
            return Ok(());
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error).with_context(|| format!("inspecting {label} {}", path.display()));
        }
    }
    let temporary = temporary_output_path(path)?;
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .with_context(|| format!("creating temporary {label} {}", temporary.display()))?;
    if let Err(error) = file.write_all(bytes).and_then(|_| file.sync_all()) {
        let _ = fs::remove_file(&temporary);
        return Err(error).with_context(|| format!("writing {label} {}", temporary.display()));
    }
    match install_open_file_no_clobber(&file, &temporary, path, Durability::File) {
        Ok(outcome) => {
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
            Ok(())
        }
        Err(install_error) => {
            let metadata = fs::symlink_metadata(path)
                .with_context(|| format!("inspecting concurrent {label} {}", path.display()))?;
            if !metadata.file_type().is_file() {
                bail!(
                    "concurrent {label} is not a regular file: {}",
                    path.display()
                );
            }
            let existing = fs::read(path)
                .with_context(|| format!("reading concurrent {label} {}", path.display()))?;
            if existing == bytes {
                Ok(())
            } else {
                Err(anyhow::Error::new(install_error)).with_context(|| {
                    format!("concurrent {label} has different content: {}", path.display())
                })
            }
        }
    }
}

fn completed_output(output: &ResolvedOutput) -> Result<CompletedOutput> {
    let metadata = fs::symlink_metadata(&output.path)
        .with_context(|| format!("reading completed output {}", output.path.display()))?;
    if metadata.file_type().is_symlink() {
        bail!(
            "completed output must not be a symlink: {}",
            output.path.display()
        );
    }
    match output.kind {
        ResolvedOutputKind::File => {
            if !metadata.is_file() {
                bail!("expected file output at {}", output.path.display());
            }
            // Resume is an exact physical-output check, not merely a logical archive identity.
            // A rooted archive commits to payload digests without reading payload bytes; using
            // that cheap identity here would trust an archive whose encoded payload was corrupted
            // after completion. The explicit `--resume` scan therefore hashes the whole artifact.
            let (bytes, identity) = full_file_identity(&output.path)?;
            Ok(CompletedOutput {
                path: output.path.clone(),
                kind: output.kind,
                bytes,
                identity,
            })
        }
        ResolvedOutputKind::Directory => {
            if !metadata.is_dir() {
                bail!("expected directory output at {}", output.path.display());
            }
            let (bytes, digest) = directory_tree_identity(&output.path)?;
            Ok(CompletedOutput {
                path: output.path.clone(),
                kind: output.kind,
                bytes,
                identity: ResolvedIdentity {
                    scheme: "directory-tree-blake3-v1".to_owned(),
                    digest,
                },
            })
        }
    }
}

fn directory_tree_identity(root: &Path) -> Result<(u64, String)> {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"gravlax-directory-tree-v1\0");
    let mut bytes = 0u64;
    hash_directory_entries(root, root, &mut hasher, &mut bytes)?;
    Ok((bytes, hasher.finalize().to_hex().to_string()))
}

fn hash_directory_entries(
    root: &Path,
    directory: &Path,
    hasher: &mut blake3::Hasher,
    bytes: &mut u64,
) -> Result<()> {
    let mut entries = fs::read_dir(directory)
        .with_context(|| format!("reading output directory {}", directory.display()))?
        .collect::<std::io::Result<Vec<_>>>()?;
    entries.sort_by_key(|entry| entry.file_name());
    for entry in entries {
        let path = entry.path();
        let relative = path
            .strip_prefix(root)
            .expect("directory entry is below root");
        let name = relative.as_os_str().as_encoded_bytes();
        let metadata = fs::symlink_metadata(&path)?;
        if metadata.file_type().is_symlink() {
            bail!("output directory contains a symlink: {}", path.display());
        } else if metadata.is_dir() {
            hasher.update(b"d");
            hasher.update(&(name.len() as u64).to_le_bytes());
            hasher.update(name);
            hash_directory_entries(root, &path, hasher, bytes)?;
        } else if metadata.is_file() {
            hasher.update(b"f");
            hasher.update(&(name.len() as u64).to_le_bytes());
            hasher.update(name);
            hasher.update(&metadata.len().to_le_bytes());
            *bytes = bytes
                .checked_add(metadata.len())
                .context("directory output byte count overflow")?;
            let mut file = fs::File::open(&path)?;
            let mut buffer = vec![0u8; 1024 * 1024];
            loop {
                let read = file.read(&mut buffer)?;
                if read == 0 {
                    break;
                }
                hasher.update(&buffer[..read]);
            }
        } else {
            bail!(
                "output directory contains a special file: {}",
                path.display()
            );
        }
    }
    Ok(())
}

/// Store a deterministic, content-addressed JSON snapshot beneath `.aie/resolved-plans`.
pub fn persist_resolved_plan(plan: &ResolvedPlan) -> Result<PathBuf> {
    let bytes = resolved_plan_bytes(plan)?;
    let digest = blake3::hash(&bytes).to_hex().to_string();
    let slug = snapshot_slug(&plan.name);
    let directory = ensure_project_directory(&plan.project_root, Path::new(".aie/resolved-plans"))?;
    let path = directory.join(format!("{slug}-{}.resolved.json", &digest[..12]));
    if path.exists() {
        ensure_existing_snapshot(&path, &bytes)?;
        return Ok(path);
    }
    let temporary = temporary_output_path(&path)?;
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .with_context(|| format!("creating temporary snapshot {}", temporary.display()))?;
    if let Err(error) = file.write_all(&bytes).and_then(|_| file.sync_all()) {
        let _ = fs::remove_file(&temporary);
        return Err(error).context("writing resolved-plan snapshot");
    }
    match install_open_file_no_clobber(&file, &temporary, &path, Durability::File) {
        Ok(outcome) => {
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
        }
        Err(install_error) => {
            if let Err(check_error) = ensure_existing_snapshot(&path, &bytes) {
                return Err(anyhow::Error::new(install_error)).with_context(|| {
                    format!("installing snapshot {} failed; {check_error:#}", path.display())
                });
            }
        }
    }
    Ok(path)
}

fn ensure_existing_snapshot(path: &Path, expected: &[u8]) -> Result<()> {
    let metadata = fs::symlink_metadata(path)
        .with_context(|| format!("inspecting existing snapshot {}", path.display()))?;
    if !metadata.file_type().is_file() {
        bail!(
            "resolved-plan snapshot is not a regular file: {}",
            path.display()
        );
    }
    let existing = fs::read(path)?;
    if existing != expected {
        bail!("resolved-plan digest collision at {}", path.display());
    }
    Ok(())
}

fn snapshot_slug(name: &str) -> String {
    let mut slug = String::new();
    for character in name.chars() {
        if slug.len() >= 48 {
            break;
        }
        if character.is_ascii_alphanumeric() || matches!(character, '-' | '_') {
            slug.push(character.to_ascii_lowercase());
        } else if !slug.ends_with('-') {
            slug.push('-');
        }
    }
    let slug = slug.trim_matches('-');
    if slug.is_empty() {
        "analysis".to_owned()
    } else {
        slug.to_owned()
    }
}

fn persist_step_completion(
    plan: &ResolvedPlan,
    plan_digest: &str,
    step: &ResolvedStep,
) -> Result<PathBuf> {
    let record = StepCompletion {
        schema_version: COMPLETION_SCHEMA_VERSION,
        resolved_plan_digest: plan_digest.to_owned(),
        step_id: step.id.clone(),
        step_digest: resolved_step_digest(step)?,
        inputs: completed_step_inputs(plan, step)?,
        outputs: completed_outputs(&step.outputs)?,
    };
    let mut bytes = serde_json::to_vec_pretty(&record)?;
    bytes.push(b'\n');
    if plan_digest.len() != 64
        || !plan_digest
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        bail!("resolved plan digest is not canonical lowercase BLAKE3 hex");
    }
    let relative_directory = Path::new(".aie/completions").join(plan_digest);
    let directory = ensure_project_directory(&plan.project_root, &relative_directory)?;
    let path = directory.join(format!("{}.json", step.id));
    let temporary = temporary_output_path(&path)?;
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .with_context(|| format!("creating temporary completion {}", temporary.display()))?;
    if let Err(error) = file.write_all(&bytes).and_then(|_| file.sync_all()) {
        let _ = fs::remove_file(&temporary);
        return Err(error).context("writing step completion record");
    }
    match install_open_file_no_clobber(&file, &temporary, &path, Durability::File) {
        Ok(outcome) => {
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
        }
        Err(install_error) => {
            let metadata = fs::symlink_metadata(&path).with_context(|| {
                format!("inspecting concurrent completion {}", path.display())
            })?;
            if !metadata.file_type().is_file() || fs::read(&path)? != bytes {
                return Err(anyhow::Error::new(install_error)).with_context(|| {
                    format!(
                        "concurrent completion record differs from the exact step completion: {}",
                        path.display()
                    )
                });
            }
        }
    }
    Ok(path)
}

fn execute_plan(
    plan: &ResolvedPlan,
    plan_digest: &str,
    decisions: &[ResumeDecision],
) -> Result<()> {
    let executable = std::env::current_exe().context("locating the aie executable")?;
    for (index, (step, decision)) in plan.steps.iter().zip(decisions).enumerate() {
        verify_step_inputs(plan, plan_digest, step)?;
        if *decision == ResumeDecision::SkipVerified {
            // Earlier steps may be long-running. Re-read the completion and every output at the
            // actual skip point instead of trusting the preflight observation indefinitely.
            verify_step_completion(plan, plan_digest, step)?;
            println!(
                "skipping step {}/{} `{}` (exact completion verified)",
                index + 1,
                plan.steps.len(),
                step.id
            );
            continue;
        }
        materialize_prepared_inputs(plan, step)?;
        println!(
            "running step {}/{} `{}` ({})",
            index + 1,
            plan.steps.len(),
            step.id,
            step.kind
        );
        prepare_output_paths(plan, &step.outputs)?;
        let mut command = ProcessCommand::new(&executable);
        command
            .args(&step.args)
            .current_dir(&plan.project_root)
            // Artifact-producing commands receive staging paths so the plan runner can publish
            // each output without clobbering. Reports, however, are part of the public result
            // contract and must describe the final logical destinations. Always override this
            // internal variable, including with `{}`, so a parent environment cannot inject path
            // aliases into a planned child.
            .env(
                LOGICAL_OUTPUT_MAP_ENV,
                logical_output_map_json(&step.outputs)?,
            );
        if let Some(output) = &step.stdout {
            let resolved = step
                .outputs
                .iter()
                .find(|candidate| &candidate.path == output)
                .context("stdout output is missing from the resolved output list")?;
            let file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&resolved.staging_path)
                .with_context(|| {
                    format!("creating staged stdout {}", resolved.staging_path.display())
                })?;
            command.stdout(Stdio::from(file));
        }
        if let Err(error) = verify_producer_executable(plan, &executable) {
            cleanup_staging_outputs(plan, &step.outputs);
            return Err(error).with_context(|| {
                format!(
                    "re-verifying the aie executable before plan step `{}`",
                    step.id
                )
            });
        }
        let status = match command.status() {
            Ok(status) => status,
            Err(error) => {
                cleanup_staging_outputs(plan, &step.outputs);
                return Err(error).with_context(|| format!("starting plan step `{}`", step.id));
            }
        };
        if !status.success() {
            cleanup_staging_outputs(plan, &step.outputs);
            bail!(
                "plan step `{}` failed with {} (resolved plan remains at .aie/resolved-plans)",
                step.id,
                status
            );
        }
        install_staged_outputs(plan, &step.outputs)?;
        if !step.outputs.is_empty() {
            let completion = persist_step_completion(plan, plan_digest, step)?;
            println!("completion: {}", completion.display());
        }
    }
    Ok(())
}

fn logical_output_map_json(outputs: &[ResolvedOutput]) -> Result<String> {
    let mut paths = BTreeMap::new();
    for output in outputs {
        let staging = output
            .staging_path
            .to_str()
            .context("resolved staging output path is not valid UTF-8")?;
        let logical = output
            .path
            .to_str()
            .context("resolved logical output path is not valid UTF-8")?;
        if paths.insert(staging, logical).is_some() {
            bail!(
                "resolved plan contains duplicate staging output path {}",
                output.staging_path.display()
            );
        }
    }
    serde_json::to_string(&paths).context("serializing logical output path map")
}

fn verify_producer_executable(plan: &ResolvedPlan, executable: &Path) -> Result<()> {
    let (_, observed) = full_file_identity(executable)?;
    if observed != plan.producer.executable_identity {
        bail!(
            "the aie executable changed after this plan was resolved; check the plan again before running it"
        );
    }
    Ok(())
}

fn prepare_output_paths(plan: &ResolvedPlan, outputs: &[ResolvedOutput]) -> Result<()> {
    for output in outputs {
        for (label, path) in [
            ("plan output", &output.path),
            ("staging output", &output.staging_path),
        ] {
            ensure_output_parent(&plan.project_root, path)?;
            revalidate_output_path(&plan.project_root, path)?;
            match fs::symlink_metadata(path) {
                Ok(_) => bail!("refusing to overwrite {label} {}", path.display()),
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                Err(error) => {
                    return Err(error)
                        .with_context(|| format!("inspecting {label} {}", path.display()));
                }
            }
        }
    }
    Ok(())
}

fn ensure_output_parent(project_root: &Path, output: &Path) -> Result<()> {
    let parent = output.parent().context("plan output has no parent")?;
    let relative = parent.strip_prefix(project_root).with_context(|| {
        format!(
            "output parent is outside project {}: {}",
            project_root.display(),
            parent.display()
        )
    })?;
    if !relative.as_os_str().is_empty() {
        ensure_project_directory(project_root, relative)?;
    }
    Ok(())
}

fn revalidate_output_path(project_root: &Path, expected: &Path) -> Result<()> {
    let relative = expected.strip_prefix(project_root).with_context(|| {
        format!(
            "resolved output escaped project {}: {}",
            project_root.display(),
            expected.display()
        )
    })?;
    let observed = resolve_project_output(project_root, relative)?;
    if observed != expected {
        bail!(
            "output path resolution changed after plan validation: {} now resolves to {}",
            expected.display(),
            observed.display()
        );
    }
    Ok(())
}

fn install_staged_outputs(plan: &ResolvedPlan, outputs: &[ResolvedOutput]) -> Result<()> {
    install_staged_outputs_with(&plan.project_root, outputs, |_, _| Ok(()))
}

fn install_staged_outputs_with<F>(
    project_root: &Path,
    outputs: &[ResolvedOutput],
    mut before_install: F,
) -> Result<()>
where
    F: FnMut(usize, &ResolvedOutput) -> Result<()>,
{
    // Validate the complete set before exposing any final path. Cross-path all-or-nothing
    // publication is not generally available from the filesystem: if a later destination loses
    // a race, preserve every complete output already installed. Removing it by pathname could
    // instead delete an object another process installed concurrently.
    let prepared = outputs
        .iter()
        .map(|output| prepare_staged_output(project_root, output))
        .collect::<Result<Vec<_>>>()?;
    let mut installed = Vec::with_capacity(outputs.len());
    for (index, prepared_output) in prepared.iter().enumerate() {
        let output = prepared_output.output;
        before_install(index, output)?;
        if let Err(error) = install_staged_output(project_root, prepared_output) {
            let installed_paths = installed
                .iter()
                .map(|installed: &&ResolvedOutput| installed.path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ");
            let staging_paths = outputs
                .iter()
                .filter(|candidate| candidate.staging_path.exists())
                .map(|candidate| candidate.staging_path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ");
            bail!(
                "{error:#}; a concurrent destination change prevented complete publication; \
                 preserved {} already-installed output(s) [{}] and staging artifacts [{}] for \
                 inspection rather than risking deletion of another process's data",
                installed.len(),
                installed_paths,
                staging_paths,
            );
        }
        installed.push(output);
    }
    #[cfg(not(target_os = "linux"))]
    for output in outputs {
        if output.kind != ResolvedOutputKind::Directory {
            continue;
        }
        if let Err(error) = fs::remove_dir_all(&output.staging_path) {
            let installed_paths = installed
                .iter()
                .map(|installed| installed.path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ");
            bail!(
                "removing staged output {} failed: {error}; all final outputs were installed \
                 completely and preserved [{}] rather than risking deletion of concurrently \
                 replaced paths",
                output.staging_path.display(),
                installed_paths,
            );
        }
    }
    Ok(())
}

struct PreparedStagedOutput<'a> {
    output: &'a ResolvedOutput,
    held: Option<fs::File>,
    metadata: Option<fs::Metadata>,
}

fn open_staged_path(path: &Path) -> Result<fs::File> {
    let mut options = OpenOptions::new();
    options.read(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.custom_flags(libc::O_NOFOLLOW);
    }
    options
        .open(path)
        .with_context(|| format!("opening staged output {}", path.display()))
}

fn prepare_staged_output<'a>(
    project_root: &Path,
    output: &'a ResolvedOutput,
) -> Result<PreparedStagedOutput<'a>> {
    revalidate_output_path(project_root, &output.path)?;
    match fs::symlink_metadata(&output.path) {
        Ok(_) => bail!(
            "refusing to replace output created while the step ran: {}",
            output.path.display()
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error)
                .with_context(|| format!("inspecting output {}", output.path.display()));
        }
    }
    revalidate_output_path(project_root, &output.staging_path)?;
    let metadata = fs::symlink_metadata(&output.staging_path).with_context(|| {
        format!(
            "step did not create declared staging output {}",
            output.staging_path.display()
        )
    })?;
    if metadata.file_type().is_symlink() {
        bail!(
            "step created a symlink instead of an output: {}",
            output.staging_path.display()
        );
    }
    match output.kind {
        ResolvedOutputKind::File if metadata.is_file() => {
            let held = open_staged_path(&output.staging_path)?;
            let held_metadata = held.metadata()?;
            let path_metadata = fs::symlink_metadata(&output.staging_path)?;
            if !same_file(&held_metadata, &path_metadata)
                || held_metadata.len() != path_metadata.len()
            {
                bail!(
                    "staging output changed while it was opened: {}",
                    output.staging_path.display()
                );
            }
            Ok(PreparedStagedOutput {
                output,
                held: Some(held),
                metadata: Some(held_metadata),
            })
        }
        ResolvedOutputKind::Directory if metadata.is_dir() => {
            validate_staged_directory(&output.staging_path)?;
            #[cfg(target_os = "linux")]
            {
                let held = open_staged_path(&output.staging_path)?;
                let held_metadata = held.metadata()?;
                let path_metadata = fs::symlink_metadata(&output.staging_path)?;
                if !held_metadata.is_dir()
                    || !path_metadata.is_dir()
                    || !same_file(&held_metadata, &path_metadata)
                {
                    bail!(
                        "staging directory changed while it was opened: {}",
                        output.staging_path.display()
                    );
                }
                Ok(PreparedStagedOutput {
                    output,
                    held: Some(held),
                    metadata: Some(held_metadata),
                })
            }
            #[cfg(not(target_os = "linux"))]
            Ok(PreparedStagedOutput {
                output,
                held: None,
                metadata: None,
            })
        }
        ResolvedOutputKind::File => bail!(
            "step created a non-file for declared output {}",
            output.path.display()
        ),
        ResolvedOutputKind::Directory => bail!(
            "step created a non-directory for declared output {}",
            output.path.display()
        ),
    }
}

fn install_staged_output(
    project_root: &Path,
    prepared: &PreparedStagedOutput<'_>,
) -> Result<()> {
    let output = prepared.output;
    revalidate_output_path(project_root, &output.path)?;
    match output.kind {
        ResolvedOutputKind::File => {
            let outcome = install_open_file_no_clobber(
                prepared
                    .held
                    .as_ref()
                    .expect("prepared file output retains its descriptor"),
                &output.staging_path,
                &output.path,
                Durability::Flush,
            )?;
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
            Ok(())
        }
        ResolvedOutputKind::Directory => {
            install_directory_no_clobber(&output.staging_path, &output.path)?;
            #[cfg(target_os = "linux")]
            {
                let expected = prepared
                    .metadata
                    .as_ref()
                    .expect("prepared directory output retains metadata");
                let observed = fs::symlink_metadata(&output.path)?;
                if !observed.is_dir() || !same_file(expected, &observed) {
                    bail!(
                        "installed directory is not the validated staging directory; preserving \
                         {} because a concurrent process may have replaced a path",
                        output.path.display()
                    );
                }
            }
            Ok(())
        }
    }
}

#[cfg(target_os = "linux")]
fn install_directory_no_clobber(staging: &Path, output: &Path) -> Result<()> {
    rename_no_replace(staging, output).with_context(|| {
        format!(
            "atomically installing output directory without overwrite {}",
            output.display()
        )
    })
}

#[cfg(target_os = "linux")]
fn rename_no_replace(source: &Path, destination: &Path) -> Result<()> {
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;

    let source = CString::new(source.as_os_str().as_bytes()).context("source path contains NUL")?;
    let destination = CString::new(destination.as_os_str().as_bytes())
        .context("destination path contains NUL")?;
    // Call the kernel entry point directly. glibc exports a `renameat2` wrapper, but musl does
    // not, so referring to `libc::renameat2` leaves a fully static musl build with an unresolved
    // symbol even though the syscall is available on every supported Linux kernel.
    let result = unsafe {
        libc::syscall(
            libc::SYS_renameat2,
            libc::AT_FDCWD,
            source.as_ptr(),
            libc::AT_FDCWD,
            destination.as_ptr(),
            libc::RENAME_NOREPLACE,
        )
    };
    if result == 0 {
        Ok(())
    } else {
        Err(std::io::Error::last_os_error()).context("renameat2(RENAME_NOREPLACE)")
    }
}

#[cfg(not(target_os = "linux"))]
fn install_directory_no_clobber(staging: &Path, output: &Path) -> Result<()> {
    fs::create_dir(output).with_context(|| {
        format!(
            "creating output directory without overwrite {}",
            output.display()
        )
    })?;
    if let Err(error) = link_directory_contents(staging, output) {
        let cleanup = fs::remove_dir_all(output);
        if let Err(cleanup) = cleanup {
            bail!("{error:#}; removing partial output also failed: {cleanup}");
        }
        return Err(error);
    }
    Ok(())
}

fn validate_staged_directory(directory: &Path) -> Result<()> {
    let entries = fs::read_dir(directory)
        .with_context(|| format!("reading staged directory {}", directory.display()))?;
    for entry in entries {
        let path = entry?.path();
        let metadata = fs::symlink_metadata(&path)?;
        if metadata.file_type().is_symlink() {
            bail!(
                "staged output directory contains a symlink: {}",
                path.display()
            );
        } else if metadata.is_dir() {
            validate_staged_directory(&path)?;
        } else if !metadata.is_file() {
            bail!(
                "staged output directory contains a special file: {}",
                path.display()
            );
        }
    }
    Ok(())
}

#[cfg(not(target_os = "linux"))]
fn link_directory_contents(source: &Path, destination: &Path) -> Result<()> {
    let mut entries = fs::read_dir(source)
        .with_context(|| format!("reading staged directory {}", source.display()))?
        .collect::<std::io::Result<Vec<_>>>()?;
    entries.sort_by_key(|entry| entry.file_name());
    for entry in entries {
        let source_path = entry.path();
        let destination_path = destination.join(entry.file_name());
        let metadata = fs::symlink_metadata(&source_path)?;
        if metadata.is_dir() {
            fs::create_dir(&destination_path)?;
            link_directory_contents(&source_path, &destination_path)?;
        } else {
            fs::hard_link(&source_path, &destination_path)?;
        }
    }
    Ok(())
}

fn cleanup_staging_outputs(plan: &ResolvedPlan, outputs: &[ResolvedOutput]) {
    for output in outputs {
        if revalidate_output_path(&plan.project_root, &output.staging_path).is_err() {
            continue;
        }
        let Ok(metadata) = fs::symlink_metadata(&output.staging_path) else {
            continue;
        };
        if metadata.is_dir() && !metadata.file_type().is_symlink() {
            let _ = fs::remove_dir_all(&output.staging_path);
        } else {
            let _ = fs::remove_file(&output.staging_path);
        }
    }
}

fn temporary_output_path(output: &Path) -> Result<PathBuf> {
    let file_name = output
        .file_name()
        .context("plan output has no file name")?
        .to_string_lossy();
    for attempt in 0u32..1_000 {
        let temporary =
            output.with_file_name(format!(".{file_name}.tmp.{}-{attempt}", std::process::id()));
        if !temporary.exists() {
            return Ok(temporary);
        }
    }
    bail!(
        "could not allocate temporary output beside {}",
        output.display()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locus_validation_is_strict() {
        assert!(validate_locus("chr1:10-20").is_ok());
        assert!(validate_locus("chr1:20-10").is_err());
        assert!(validate_locus("chr1:ten-20").is_err());
        assert!(validate_locus("chr1 10-20").is_err());
        assert!(validate_locus("chr1:4294967296-4294967297").is_err());
    }

    #[test]
    fn canonical_design_paths_must_be_safe_tsv_fields() {
        assert!(require_tsv_field("/data/sample.aie", "path").is_ok());
        assert!(require_tsv_field("/data/with\ttab.aie", "path").is_err());
        assert!(require_tsv_field("/data/with\nnewline.aie", "path").is_err());
    }

    #[test]
    fn child_logical_output_map_is_exact_and_versioned_by_its_environment_name() {
        let outputs = vec![ResolvedOutput {
            name: "annotation".to_owned(),
            path: PathBuf::from("/project/results/genes.aic"),
            kind: ResolvedOutputKind::File,
            resource_kind: ResourceKind::Annotation,
            staging_path: PathBuf::from("/project/results/.genes.aie-stage-compile.aic"),
            annotation_semantics: None,
            assembly: None,
        }];
        assert_eq!(LOGICAL_OUTPUT_MAP_ENV, "GRAVLAX_LOGICAL_OUTPUT_MAP_V1");
        assert_eq!(logical_output_map_json(&[]).unwrap(), "{}");
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&logical_output_map_json(&outputs).unwrap())
                .unwrap(),
            serde_json::json!({
                "/project/results/.genes.aie-stage-compile.aic":
                    "/project/results/genes.aic"
            })
        );
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn directory_install_is_complete_and_never_replaces_an_existing_target() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gravlax-directory-install-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        let staging = root.join("staging");
        let output = root.join("output");
        fs::create_dir(&staging).unwrap();
        fs::write(staging.join("artifact"), b"complete").unwrap();
        install_directory_no_clobber(&staging, &output).unwrap();
        assert!(!staging.exists());
        assert_eq!(fs::read(output.join("artifact")).unwrap(), b"complete");

        let blocked_staging = root.join("blocked-staging");
        let blocked_output = root.join("blocked-output");
        fs::create_dir(&blocked_staging).unwrap();
        fs::write(blocked_staging.join("new"), b"new").unwrap();
        fs::create_dir(&blocked_output).unwrap();
        fs::write(blocked_output.join("existing"), b"existing").unwrap();
        assert!(install_directory_no_clobber(&blocked_staging, &blocked_output).is_err());
        assert!(blocked_staging.join("new").is_file());
        assert_eq!(
            fs::read(blocked_output.join("existing")).unwrap(),
            b"existing"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn later_install_race_preserves_concurrent_replacement_and_complete_partial_set() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gravlax-output-set-race-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        let first_path = root.join("first.txt");
        let second_path = root.join("second.txt");
        let first_staging = root.join(".first.stage");
        let second_staging = root.join(".second.stage");
        fs::write(&first_staging, b"first from runner").unwrap();
        fs::write(&second_staging, b"second from runner").unwrap();
        let outputs = vec![
            ResolvedOutput {
                name: "first".to_owned(),
                path: first_path.clone(),
                kind: ResolvedOutputKind::File,
                resource_kind: ResourceKind::File,
                staging_path: first_staging.clone(),
                annotation_semantics: None,
                assembly: None,
            },
            ResolvedOutput {
                name: "second".to_owned(),
                path: second_path.clone(),
                kind: ResolvedOutputKind::File,
                resource_kind: ResourceKind::File,
                staging_path: second_staging.clone(),
                annotation_semantics: None,
                assembly: None,
            },
        ];

        let error = install_staged_outputs_with(&root, &outputs, |index, _| {
            if index == 1 {
                fs::remove_file(&first_path)?;
                fs::write(&first_path, b"concurrent replacement")?;
                fs::write(&second_path, b"concurrent destination")?;
            }
            Ok(())
        })
        .unwrap_err()
        .to_string();

        assert!(error.contains("preserved 1 already-installed output"));
        assert_eq!(fs::read(&first_path).unwrap(), b"concurrent replacement");
        assert_eq!(fs::read(&second_path).unwrap(), b"concurrent destination");
        assert!(!first_staging.exists());
        assert!(!second_staging.exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn plan_file_install_uses_held_descriptor_when_staging_path_is_replaced() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gravlax-plan-staging-race-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        let output_path = root.join("result.txt");
        let staging_path = root.join(".result.stage");
        let displaced_path = root.join("held-produced-file");
        fs::write(&staging_path, b"produced by plan step").unwrap();
        let outputs = vec![ResolvedOutput {
            name: "result".to_owned(),
            path: output_path.clone(),
            kind: ResolvedOutputKind::File,
            resource_kind: ResourceKind::File,
            staging_path: staging_path.clone(),
            annotation_semantics: None,
            assembly: None,
        }];

        install_staged_outputs_with(&root, &outputs, |_, _| {
            fs::rename(&staging_path, &displaced_path)?;
            fs::write(&staging_path, b"concurrent replacement")?;
            Ok(())
        })
        .unwrap();

        assert_eq!(fs::read(&output_path).unwrap(), b"produced by plan step");
        assert_eq!(fs::read(&displaced_path).unwrap(), b"produced by plan step");
        assert_eq!(fs::read(&staging_path).unwrap(), b"concurrent replacement");
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn snapshot_names_are_portable() {
        assert_eq!(snapshot_slug("My plan: v1"), "my-plan-v1");
        assert_eq!(snapshot_slug("***"), "analysis");
    }

    #[test]
    fn unsupported_kind_has_a_targeted_error() {
        let source = br#"schema_version: 1
steps:
  - id: nope
    kind: magical-query
"#;
        let error = validate_declared_step_kinds(source, "yaml", Path::new("plan.yaml"))
            .unwrap_err()
            .to_string();
        assert!(error.contains("unsupported plan step kind `magical-query`"));
    }

    #[test]
    fn completion_identity_detects_archive_payload_corruption() {
        use std::io::{Read as _, Seek as _, SeekFrom, Write as _};

        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-completion-identity-{}-{nonce}.aie",
            std::process::id()
        ));
        let mut writer = evidence_io::format::SectionWriter::create(&path, 1).unwrap();
        writer
            .section(
                "payload",
                b"payload bytes that compress into a nonempty frame",
            )
            .unwrap();
        writer.finish().unwrap();

        let reader = evidence_io::format::SectionReader::open(&path).unwrap();
        let (_, offset, _, compressed_len) = &reader.entries()[0];
        assert!(*compressed_len > 0);
        let payload_offset = *offset + 1 + "payload".len() as u64 + 16;
        let (_, logical_before) = resolve_resource_identity(ResourceKind::Archive, &path).unwrap();
        let output = ResolvedOutput {
            name: "archive".to_owned(),
            path: path.clone(),
            kind: ResolvedOutputKind::File,
            resource_kind: ResourceKind::Archive,
            staging_path: path.with_extension("stage.aie"),
            annotation_semantics: None,
            assembly: None,
        };
        let physical_before = completed_output(&output).unwrap();

        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .open(&path)
            .unwrap();
        file.seek(SeekFrom::Start(payload_offset)).unwrap();
        let mut byte = [0u8; 1];
        file.read_exact(&mut byte).unwrap();
        byte[0] ^= 0x01;
        file.seek(SeekFrom::Start(payload_offset)).unwrap();
        file.write_all(&byte).unwrap();
        file.sync_all().unwrap();

        let (_, logical_after) = resolve_resource_identity(ResourceKind::Archive, &path).unwrap();
        let physical_after = completed_output(&output).unwrap();
        assert_eq!(logical_before, logical_after);
        assert_ne!(physical_before.identity, physical_after.identity);
        fs::remove_file(path).unwrap();
    }
}
