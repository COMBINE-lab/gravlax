//! Loopback-only, read-only project Explorer and scientific plan builder.
//!
//! The browser is deliberately a thin renderer. Feature resolution, typed plan construction,
//! command compilation, route explanations, and resource identities all come from the Rust
//! backend. The HTTP layer can inspect and export those values, but cannot execute a plan or
//! modify a project.

use crate::projectcmd::{load_project, ProjectContext, ResourceKind};
use anyhow::{bail, Context, Result};
use clap::Parser;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{Read, Write};
use std::net::{Ipv4Addr, SocketAddr, TcpListener, TcpStream};
use std::path::{Component, Path, PathBuf};
use std::sync::Arc;
use std::time::{Duration, UNIX_EPOCH};

const PROJECT_MANIFEST: &str = "aie-project.yaml";
const MAX_REQUEST_BYTES: usize = 16 * 1024;
const MAX_QUERY_BYTES: usize = 12 * 1024;
const MAX_QUERY_FIELDS: usize = 256;
const MAX_QUERY_KEY_BYTES: usize = 64;
const MAX_QUERY_VALUE_BYTES: usize = 4 * 1024;
const MAX_ARTIFACTS: usize = 10_000;
const MAX_STRUCTURED_PREVIEW_BYTES: u64 = 16 * 1024 * 1024;
const EXPLORER_HTML: &[u8] = include_bytes!("../assets/explorer.html");

#[derive(Parser, Debug)]
pub struct Args {
    /// Project directory or `aie-project.yaml`; defaults to discovery from the current directory.
    #[arg(long, value_name = "PATH")]
    pub project: Option<PathBuf>,

    /// Loopback TCP port. Use 0 to ask the operating system for an available port.
    #[arg(long, default_value_t = 8787)]
    pub port: u16,
}

#[derive(Clone, Debug)]
pub struct ArtifactSummary {
    pub id: String,
    pub kind: &'static str,
    pub bytes: u64,
    pub modified_unix_seconds: Option<u64>,
    pub media_type: &'static str,
    pub status: &'static str,
    pub detail: Option<String>,
    pub data: Value,
}

impl ArtifactSummary {
    fn value(&self) -> Value {
        json!({
            "id": self.id,
            "kind": self.kind,
            "bytes": self.bytes,
            "modified_unix_seconds": self.modified_unix_seconds,
            "media_type": self.media_type,
            "status": self.status,
            "detail": self.detail,
            "data": self.data,
            "download_url": format!("/api/v1/artifact?path={}", percent_encode(&self.id)),
        })
    }
}

pub struct OpenArtifact {
    pub summary: ArtifactSummary,
    pub file: File,
}

#[derive(Clone, Debug)]
struct FeatureLookup {
    annotation: String,
    feature: String,
    assembly: Option<String>,
    annotation_label: Option<String>,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
enum PlannerView {
    AnnotationComparison,
    Region,
    SpliceEvents,
    Junction,
    JunctionSet,
    TerminalBoundary,
}

impl PlannerView {
    fn parse(value: &str) -> Result<Self> {
        match value {
            "annotation-comparison" => Ok(Self::AnnotationComparison),
            "region" => Ok(Self::Region),
            "splice-events" => Ok(Self::SpliceEvents),
            "junction" => Ok(Self::Junction),
            "junction-set" => Ok(Self::JunctionSet),
            "terminal-boundary" => Ok(Self::TerminalBoundary),
            _ => bail!("view must be annotation-comparison, region, splice-events, junction, junction-set, or terminal-boundary"),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::AnnotationComparison => "annotation-comparison",
            Self::Region => "region",
            Self::SpliceEvents => "splice-events",
            Self::Junction => "junction",
            Self::JunctionSet => "junction-set",
            Self::TerminalBoundary => "terminal-boundary",
        }
    }
}

#[derive(Clone, Debug)]
struct PlannerRequest {
    view: PlannerView,
    archive: String,
    annotation_a: Option<String>,
    annotation_b: Option<String>,
    comparison_gene_key: crate::plancmd::ComparisonGeneKey,
    comparison_solo_strand: crate::plancmd::SoloStrand,
    max_molecule_witnesses: usize,
    max_row_transitions_per_molecule: usize,
    allow_identical: bool,
    annotation: Option<String>,
    feature: Option<String>,
    assembly: Option<String>,
    annotation_label: Option<String>,
    junction: Option<String>,
    include: Vec<String>,
    exclude: Vec<String>,
    cells: Option<String>,
    groups: Option<String>,
    aggregation: crate::plancmd::QueryAggregation,
    site_gap: u32,
    strand: Option<String>,
}

/// Read-only boundary between the Explorer transport/UI and project or plan implementations.
///
/// Artifact ids are backend-owned opaque strings.  `open_artifact` must repeat its own access
/// checks: callers must never infer a filesystem path directly from an HTTP request.
trait ExplorerBackend: Send + Sync {
    fn overview(&self) -> Result<Value>;
    fn artifacts(&self) -> Result<Vec<ArtifactSummary>>;
    fn open_artifact(&self, id: &str) -> Result<OpenArtifact>;
    fn resolve_feature(&self, request: &FeatureLookup) -> Result<Value>;
    fn preview_plan(&self, request: &PlannerRequest) -> Result<Value>;
}

#[derive(Debug)]
pub struct FilesystemBackend {
    root: PathBuf,
    project: Option<ProjectContext>,
}

impl FilesystemBackend {
    pub fn discover(selected: &Path) -> Result<Self> {
        if !selected.exists() {
            bail!(
                "Explorer project path {} does not exist",
                selected.display()
            );
        }
        let selected = selected
            .canonicalize()
            .with_context(|| format!("resolving Explorer project path {}", selected.display()))?;
        let (root, manifest) = if selected.is_file() {
            if selected
                .file_name()
                .is_none_or(|name| name != PROJECT_MANIFEST)
            {
                bail!(
                    "Explorer expects a project directory or {PROJECT_MANIFEST}, not {}",
                    selected.display()
                );
            }
            let root = selected
                .parent()
                .context("project manifest has no parent directory")?
                .to_path_buf();
            (root, Some(selected))
        } else {
            let mut cursor = Some(selected.as_path());
            let mut found = None;
            while let Some(directory) = cursor {
                let candidate = directory.join(PROJECT_MANIFEST);
                if candidate.is_file() {
                    found = Some(candidate);
                    break;
                }
                cursor = directory.parent();
            }
            if let Some(manifest) = found {
                (manifest.parent().unwrap().to_path_buf(), Some(manifest))
            } else {
                (selected, None)
            }
        };
        let project = manifest
            .as_deref()
            .map(load_project)
            .transpose()
            .context("loading the Explorer project")?;
        Ok(Self { root, project })
    }

    pub fn root(&self) -> &Path {
        &self.root
    }

    fn project(&self) -> Result<&ProjectContext> {
        self.project
            .as_ref()
            .context("Explorer scientific planning requires an aie-project.yaml")
    }

    fn qualified_feature(
        &self,
        annotation: &str,
        feature: &str,
        expected_assembly: Option<&str>,
        expected_annotation: Option<&str>,
    ) -> Result<crate::plancmd::FeatureRequest> {
        let project = self.project()?;
        let resource = project
            .manifest
            .resources
            .get(annotation)
            .with_context(|| format!("project resource `{annotation}` is not registered"))?;
        if resource.kind != ResourceKind::Annotation {
            bail!(
                "project resource `{annotation}` is {}, not annotation",
                resource.kind.as_str()
            );
        }
        let identity = resource.annotation_identity.as_ref().with_context(|| {
            format!(
                "annotation resource `{annotation}` has no assembly/release identity; register it with project add --assembly and --annotation-label"
            )
        })?;
        if expected_assembly.is_some_and(|expected| expected != identity.assembly) {
            bail!(
                "annotation resource `{annotation}` is registered for {}, not {}",
                identity.assembly,
                expected_assembly.unwrap()
            );
        }
        if expected_annotation.is_some_and(|expected| expected != identity.annotation) {
            bail!(
                "annotation resource `{annotation}` is registered as {}, not {}",
                identity.annotation,
                expected_annotation.unwrap()
            );
        }
        if feature.trim().is_empty()
            || feature.len() > 1_024
            || feature.chars().any(char::is_control)
        {
            bail!("feature must contain 1-1024 printable characters");
        }
        Ok(crate::plancmd::FeatureRequest::Qualified(
            crate::plancmd::QualifiedFeatureRequest {
                identifier: feature.to_owned(),
                assembly: Some(identity.assembly.clone()),
                annotation: Some(identity.annotation.clone()),
            },
        ))
    }

    fn relative_id(&self, path: &Path) -> Result<String> {
        let relative = path
            .strip_prefix(&self.root)
            .with_context(|| format!("{} is outside the Explorer project", path.display()))?;
        let mut parts = Vec::new();
        for component in relative.components() {
            match component {
                Component::Normal(part) => parts.push(part.to_string_lossy().into_owned()),
                _ => bail!("artifact path contains a non-normal component"),
            }
        }
        Ok(parts.join("/"))
    }

    fn classify(id: &str) -> Option<&'static str> {
        if id == PROJECT_MANIFEST {
            Some("project")
        } else if id == "plans" || id.starts_with("plans/") {
            Some("plan")
        } else if id == ".aie/resolved-plans" || id.starts_with(".aie/resolved-plans/") {
            Some("resolved-plan")
        } else if id == "results"
            || id.starts_with("results/")
            || id == ".aie/results"
            || id.starts_with(".aie/results/")
        {
            Some("result")
        } else {
            None
        }
    }

    fn allowed_extension(kind: &str, path: &Path) -> bool {
        let extension = path
            .extension()
            .and_then(|value| value.to_str())
            .unwrap_or("");
        match kind {
            "project" => path
                .file_name()
                .is_some_and(|name| name == PROJECT_MANIFEST),
            "plan" => matches!(extension, "yaml" | "yml" | "json"),
            "resolved-plan" => extension == "json",
            // Everything beneath the project-owned results directory is an explicit artifact;
            // archives and matrices can be streamed without teaching the UI each file format.
            "result" => true,
            _ => false,
        }
    }

    fn collect_directory(
        &self,
        directory: &Path,
        depth: usize,
        out: &mut Vec<ArtifactSummary>,
    ) -> Result<()> {
        if depth > 8 || !directory.is_dir() {
            return Ok(());
        }
        let mut entries: Vec<_> = std::fs::read_dir(directory)
            .with_context(|| {
                format!(
                    "reading Explorer artifact directory {}",
                    directory.display()
                )
            })?
            .collect::<std::io::Result<_>>()?;
        entries.sort_by_key(|entry| entry.file_name());
        for entry in entries {
            if out.len() >= MAX_ARTIFACTS {
                bail!("Explorer artifact catalogue exceeds the {MAX_ARTIFACTS}-file safety limit");
            }
            let file_type = entry.file_type()?;
            if file_type.is_symlink() {
                continue;
            }
            let path = entry.path();
            if file_type.is_dir() {
                self.collect_directory(&path, depth + 1, out)?;
            } else if file_type.is_file() {
                if let Some(summary) = self.summarize(&path)? {
                    out.push(summary);
                }
            }
        }
        Ok(())
    }

    fn summarize(&self, path: &Path) -> Result<Option<ArtifactSummary>> {
        let id = self.relative_id(path)?;
        let Some(kind) = Self::classify(&id) else {
            return Ok(None);
        };
        if !Self::allowed_extension(kind, path) {
            return Ok(None);
        }
        let metadata = path.metadata()?;
        let modified_unix_seconds = metadata
            .modified()
            .ok()
            .and_then(|modified| modified.duration_since(UNIX_EPOCH).ok())
            .map(|duration| duration.as_secs());
        let (status, detail, data) = match kind {
            "project" => (
                "valid",
                None,
                self.project.as_ref().map_or(Value::Null, |project| {
                    json!({
                        "name": project.manifest.name,
                        "schema_version": project.manifest.schema_version,
                        "resources": project.manifest.resources.len(),
                    })
                }),
            ),
            "plan" => match &self.project {
                Some(_) if metadata.len() > MAX_STRUCTURED_PREVIEW_BYTES => (
                    "unresolved",
                    Some(format!(
                        "plan exceeds the {} MiB validation-preview limit",
                        MAX_STRUCTURED_PREVIEW_BYTES / 1024 / 1024
                    )),
                    Value::Null,
                ),
                Some(project) => match crate::plancmd::resolve_plan(path, project) {
                    Ok(resolved) => (
                        "valid",
                        None,
                        json!({
                            "name": resolved.name,
                            "source_digest": resolved.source_digest,
                            "resources": resolved.resources.len(),
                            "steps": resolved.steps.len(),
                            "step_kinds": resolved.steps.iter().map(|step| &step.kind).collect::<Vec<_>>(),
                        }),
                    ),
                    Err(error) => ("invalid", Some(format!("{error:#}")), Value::Null),
                },
                None => (
                    "unresolved",
                    Some(format!("no {PROJECT_MANIFEST} is available")),
                    Value::Null,
                ),
            },
            "resolved-plan" if metadata.len() > MAX_STRUCTURED_PREVIEW_BYTES => (
                "snapshot",
                Some(format!(
                    "snapshot exceeds the {} MiB structured-preview limit",
                    MAX_STRUCTURED_PREVIEW_BYTES / 1024 / 1024
                )),
                Value::Null,
            ),
            "resolved-plan" => {
                match std::fs::read(path)
                    .ok()
                    .and_then(|bytes| serde_json::from_slice::<Value>(&bytes).ok())
                {
                    Some(value) => (
                        "snapshot",
                        None,
                        json!({
                            "schema_version": value.get("schema_version"),
                            "name": value.get("name"),
                            "source_digest": value.get("source_digest"),
                            "steps": value.get("steps").and_then(Value::as_array).map(Vec::len),
                        }),
                    ),
                    None => (
                        "invalid",
                        Some("resolved plan is not valid JSON".to_owned()),
                        Value::Null,
                    ),
                }
            }
            "result" => ("artifact", None, Value::Null),
            _ => unreachable!(),
        };
        Ok(Some(ArtifactSummary {
            id,
            kind,
            bytes: metadata.len(),
            modified_unix_seconds,
            media_type: media_type(path),
            status,
            detail,
            data,
        }))
    }

    fn resolved_artifact_path(&self, id: &str) -> Result<PathBuf> {
        let requested = Path::new(id);
        if requested.is_absolute()
            || requested
                .components()
                .any(|component| !matches!(component, Component::Normal(_)))
        {
            bail!("invalid Explorer artifact id");
        }
        let kind = Self::classify(id).context("artifact is outside the Explorer catalogue")?;
        let joined = self.root.join(requested);
        if !Self::allowed_extension(kind, &joined) {
            bail!("artifact type is not available in Explorer");
        }
        let canonical = joined
            .canonicalize()
            .with_context(|| format!("resolving Explorer artifact {id}"))?;
        if !canonical.starts_with(&self.root) || !canonical.is_file() {
            bail!("artifact resolves outside the Explorer project");
        }
        Ok(canonical)
    }
}

fn planner_scope(request: &PlannerRequest) -> crate::plancmd::QueryScope {
    crate::plancmd::QueryScope {
        cells: request.cells.clone(),
        groups: request.groups.clone(),
        aggregation: request.aggregation,
    }
}

fn planner_feature(
    backend: &FilesystemBackend,
    request: &PlannerRequest,
) -> Result<Option<crate::plancmd::FeatureRequest>> {
    match (&request.annotation, &request.feature) {
        (Some(annotation), Some(feature)) => backend
            .qualified_feature(
                annotation,
                feature,
                request.assembly.as_deref(),
                request.annotation_label.as_deref(),
            )
            .map(Some),
        (None, None) => Ok(None),
        _ => bail!("annotation and feature must be supplied together"),
    }
}

fn build_planner_model(
    backend: &FilesystemBackend,
    request: &PlannerRequest,
) -> Result<crate::plancmd::AnalysisPlan> {
    let feature = planner_feature(backend, request)?;
    let scope = planner_scope(request);
    let annotation = request.annotation.clone();
    let result_path = PathBuf::from(format!("results/explorer-{}.json", request.view.as_str()));
    let uniform_json = || {
        Some(crate::plancmd::UniformOutput {
            format: crate::plancmd::UniformFormat::Json,
            output: Some(result_path.clone()),
        })
    };
    let step = match request.view {
        PlannerView::AnnotationComparison => crate::plancmd::PlanStep::CompareAnnotations {
            id: "explore".to_owned(),
            archive: request.archive.clone(),
            annotation_a: request
                .annotation_a
                .clone()
                .context("annotation-comparison view requires annotation A")?,
            annotation_b: request
                .annotation_b
                .clone()
                .context("annotation-comparison view requires annotation B")?,
            gene_key: request.comparison_gene_key,
            solo_strand: request.comparison_solo_strand,
            max_molecule_witnesses: request.max_molecule_witnesses,
            max_row_transitions_per_molecule: request.max_row_transitions_per_molecule,
            allow_identical: request.allow_identical,
            format: crate::plancmd::OutputFormat::Json,
            table: None,
            output: Some(result_path.clone()),
        },
        PlannerView::Region => crate::plancmd::PlanStep::QueryRegion {
            id: "explore".to_owned(),
            archive: request.archive.clone(),
            locus: None,
            feature,
            top: 20,
            annotation,
            format: crate::plancmd::OutputFormat::Human,
            output: None,
            uniform_output: uniform_json(),
            scope,
        },
        PlannerView::SpliceEvents => crate::plancmd::PlanStep::QueryEvents {
            id: "explore".to_owned(),
            archive: request.archive.clone(),
            locus: None,
            feature,
            event_types: Vec::new(),
            min_support: 2,
            min_informative: 1,
            max_events: 100_000,
            top: 20,
            annotation,
            format: crate::plancmd::OutputFormat::Human,
            output: None,
            uniform_output: uniform_json(),
            scope,
        },
        PlannerView::Junction => crate::plancmd::PlanStep::QueryJunction {
            id: "explore".to_owned(),
            archive: request.archive.clone(),
            locus: request
                .junction
                .clone()
                .context("junction view requires a junction")?,
            top: 20,
            format: crate::plancmd::OutputFormat::Human,
            output: None,
            uniform_output: uniform_json(),
            scope,
        },
        PlannerView::JunctionSet => crate::plancmd::PlanStep::QueryJset {
            id: "explore".to_owned(),
            archive: request.archive.clone(),
            include: request.include.clone(),
            exclude: request.exclude.clone(),
            top: 20,
            format: crate::plancmd::OutputFormat::Human,
            output: None,
            uniform_output: uniform_json(),
            scope,
        },
        PlannerView::TerminalBoundary => {
            if request.cells.is_some() {
                bail!(
                    "terminal-boundary view supports all archive cells or a groups resource; the underlying APA query does not accept a cells-list scope"
                );
            }
            if request.aggregation != crate::plancmd::QueryAggregation::Auto {
                bail!("terminal-boundary view uses the APA query's all-cell or group aggregation and does not accept an aggregation override");
            }
            let annotation_name = request
                .annotation
                .as_deref()
                .context("terminal-boundary view requires annotation")?;
            let feature_request = feature
                .as_ref()
                .context("terminal-boundary view requires feature")?;
            let resolved = crate::plancmd::resolve_project_feature(
                backend.project()?,
                annotation_name,
                feature_request,
            )?;
            let strand = match request.strand.as_deref() {
                Some("+") => Some(crate::plancmd::ApaStrand::Forward),
                Some("-") => Some(crate::plancmd::ApaStrand::Reverse),
                None => Some(match resolved.strand {
                    anno::intent::Strand::Forward => crate::plancmd::ApaStrand::Forward,
                    anno::intent::Strand::Reverse => crate::plancmd::ApaStrand::Reverse,
                }),
                Some(_) => unreachable!("planner query parsing validates strand"),
            };
            crate::plancmd::PlanStep::QueryApa {
                id: "explore".to_owned(),
                archive: request.archive.clone(),
                locus: None,
                feature,
                annotation,
                site_gap: request.site_gap,
                strand,
                tsv: false,
                groups: request.groups.clone(),
                genome: None,
                drop_ip: false,
                permute: 0,
                seed: 1,
                plot: None,
                output: None,
                uniform_output: uniform_json(),
            }
        }
    };
    Ok(crate::plancmd::AnalysisPlan {
        schema_version: crate::plancmd::PLAN_SCHEMA_VERSION,
        name: Some(format!("explorer-{}", request.view.as_str())),
        steps: vec![step],
    })
}

fn shell_argument(argument: &str) -> String {
    if !argument.is_empty()
        && argument
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"/._-:=".contains(&byte))
    {
        argument.to_owned()
    } else {
        format!("'{}'", argument.replace('\'', "'\"'\"'"))
    }
}

fn logical_export_arguments(step: &crate::plancmd::ResolvedStep) -> Result<Vec<String>> {
    let mut arguments = step.args.clone();
    for output in &step.outputs {
        let staged = output.staging_path.to_str().with_context(|| {
            format!(
                "resolved staging path is not valid UTF-8: {}",
                output.staging_path.display()
            )
        })?;
        let destination = output.path.to_str().with_context(|| {
            format!(
                "resolved output path is not valid UTF-8: {}",
                output.path.display()
            )
        })?;
        for argument in &mut arguments {
            if argument == staged {
                *argument = destination.to_owned();
            }
        }
    }
    Ok(arguments)
}

fn cli_export(arguments: &[String]) -> String {
    let mut out = String::from("aie");
    for argument in arguments {
        out.push(' ');
        out.push_str(&shell_argument(argument));
    }
    out.push('\n');
    out
}

fn python_export(arguments: &[String]) -> Result<String> {
    let arguments = serde_json::to_string_pretty(arguments)?;
    Ok(format!(
        "from gravlax import Client\n\naie = Client()\narguments = {arguments}\nresult = aie.run(arguments)\nprint(result.stdout, end=\"\")\n"
    ))
}

fn planner_route(
    request: &PlannerRequest,
    step: &crate::plancmd::ResolvedStep,
) -> Result<(String, String, String)> {
    let mut scope = if let Some(groups) = &request.groups {
        format!(
            "Only barcodes listed in groups resource `{groups}` contribute; output is aggregated by its group labels."
        )
    } else if let Some(cells) = &request.cells {
        format!(
            "Only barcodes listed in cells resource `{cells}` contribute; other archive cells are excluded."
        )
    } else {
        "All cells represented in the selected archive contribute.".to_owned()
    };
    match request.aggregation {
        crate::plancmd::QueryAggregation::Auto => {}
        crate::plancmd::QueryAggregation::Cell => {
            scope.push_str(" Output remains at per-cell resolution.");
        }
        crate::plancmd::QueryAggregation::Group => {
            scope.push_str(" Output is reduced to the registered group labels.");
        }
        crate::plancmd::QueryAggregation::Bulk => {
            scope.push_str(" Selected evidence is reduced to one bulk total.");
        }
    }
    Ok(match request.view {
        PlannerView::AnnotationComparison => {
            let comparison = step.annotation_comparison.as_ref().context(
                "resolved annotation-comparison step has no comparison intent disclosure",
            )?;
            let execution = step
                .explanation
                .iter()
                .find_map(|line| line.strip_prefix("execution -> "))
                .context("resolved annotation-comparison step has no execution disclosure")?;
            (
                execution.to_owned(),
                "All archived molecule classes contribute to both annotation reductions; this preview has no locus, cell, group, or output filter.".to_owned(),
                format!(
                    "Final count deltas: {}. Transition evidence: {}.",
                    comparison.final_count_delta_semantics,
                    comparison.transition_evidence_semantics
                ),
            )
        }
        PlannerView::Region => (
            "The range index selects overlapping archive chunks, then exact molecule anchors in the resolved 0-based half-open interval are reduced by cell or group.".to_owned(),
            scope,
            "This is an anchor-region query: evidence whose reported alternatives never anchor in the interval is outside the result. The annotation resolves the feature and labels the view; it does not change archived alignments.".to_owned(),
        ),
        PlannerView::SpliceEvents => (
            "The junction catalogue proposes coordinate-defined alternative-donor, alternative-acceptor, and cassette events; the union of their postings is decoded once for exact class-level usage.".to_owned(),
            scope,
            "Candidates below catalogue support 2 or exact informative support 1 are excluded. Molecules supporting neither side are not in the informative denominator; a class supporting both sides is reported explicitly.".to_owned(),
        ),
        PlannerView::Junction => (
            "The exact donor/acceptor pair is looked up in the junction catalogue; its postings select only supporting chunks before class-deduplicated counts are reduced.".to_owned(),
            scope,
            "Coordinates are exact 0-based splice boundaries. Molecules without the requested junction and cells outside the selected scope are excluded.".to_owned(),
        ),
        PlannerView::JunctionSet => (
            "All requested junction postings are unioned, decoded once, and converted to per-class inclusion/exclusion masks before cell or group reduction.".to_owned(),
            scope,
            "Classes supporting neither side are excluded from informative usage; inclusion-only, exclusion-only, and both-side evidence remain distinct. Duplicate or cross-side junctions are rejected.".to_owned(),
        ),
        PlannerView::TerminalBoundary => (
            "The range index selects the resolved feature window; transcript-oriented 3'-most aligned molecule boundaries are clustered into strand-aware sites with the requested gap.".to_owned(),
            scope,
            "Aligned fragment boundaries are evidence locations, not asserted cleavage sites. This plan does not consult genome sequence, so it reports rather than filters possible internal priming; evidence outside the feature interval or selected strand is excluded.".to_owned(),
        ),
    })
}

impl ExplorerBackend for FilesystemBackend {
    fn overview(&self) -> Result<Value> {
        let artifact_counts = self.artifacts()?.into_iter().fold(
            std::collections::BTreeMap::<&'static str, usize>::new(),
            |mut counts, artifact| {
                *counts.entry(artifact.kind).or_default() += 1;
                counts
            },
        );
        let resources = self.project.as_ref().map_or_else(Vec::new, |project| {
            project
                .manifest
                .resources
                .iter()
                .map(
                    |(name, resource)| match project.resolve_resource(name, &[resource.kind]) {
                        Ok((_, path)) => json!({
                            "name": name,
                            "kind": resource.kind.as_str(),
                            "declared_path": resource.path,
                            "resolved_path": path,
                            "external": resource.external,
                            "bytes": path.metadata().map(|metadata| metadata.len()).ok(),
                            "assembly": resource.assembly,
                            "annotation_identity": resource.annotation_identity,
                            "status": "ok",
                        }),
                        Err(error) => json!({
                            "name": name,
                            "kind": resource.kind.as_str(),
                            "declared_path": resource.path,
                            "external": resource.external,
                            "assembly": resource.assembly,
                            "annotation_identity": resource.annotation_identity,
                            "status": "invalid",
                            "error": format!("{error:#}"),
                        }),
                    },
                )
                .collect()
        });
        Ok(json!({
            "schema": "gravlax.explorer.overview.v1",
            "read_only": true,
            "project": {
                "root": self.root,
                "manifest": self.project.as_ref().map(|project| &project.manifest_path),
                "configured": self.project.is_some(),
                "display_name": self.project.as_ref().map(|project| project.manifest.name.as_str())
                    .or_else(|| self.root.file_name().and_then(|name| name.to_str())),
                "resources": resources,
            },
            "artifact_counts": artifact_counts,
        }))
    }

    fn artifacts(&self) -> Result<Vec<ArtifactSummary>> {
        let mut out = Vec::new();
        if let Some(project) = &self.project {
            if let Some(summary) = self.summarize(&project.manifest_path)? {
                out.push(summary);
            }
        }
        for directory in [
            self.root.join("plans"),
            self.root.join(".aie/resolved-plans"),
            self.root.join("results"),
            self.root.join(".aie/results"),
        ] {
            self.collect_directory(&directory, 0, &mut out)?;
        }
        out.sort_by(|left, right| left.id.cmp(&right.id));
        Ok(out)
    }

    fn open_artifact(&self, id: &str) -> Result<OpenArtifact> {
        let path = self.resolved_artifact_path(id)?;
        let summary = self
            .summarize(&path)?
            .context("artifact is not in the Explorer catalogue")?;
        let file = File::open(&path).with_context(|| format!("opening Explorer artifact {id}"))?;
        let confirmed = path
            .canonicalize()
            .with_context(|| format!("re-resolving Explorer artifact {id} after open"))?;
        if confirmed != path || !confirmed.starts_with(&self.root) {
            bail!("artifact changed path while Explorer was opening it");
        }
        let opened = file.metadata()?;
        let named = confirmed.metadata()?;
        if opened.len() != summary.bytes || !same_artifact_file(&opened, &named) {
            bail!("artifact changed while Explorer was opening it");
        }
        Ok(OpenArtifact { summary, file })
    }

    fn resolve_feature(&self, request: &FeatureLookup) -> Result<Value> {
        let feature = self.qualified_feature(
            &request.annotation,
            &request.feature,
            request.assembly.as_deref(),
            request.annotation_label.as_deref(),
        )?;
        let resolved =
            crate::plancmd::resolve_project_feature(self.project()?, &request.annotation, &feature)
                .with_context(|| {
                    format!(
                        "resolving `{}` against annotation resource `{}`",
                        request.feature, request.annotation
                    )
                })?;
        Ok(json!({
            "schema": "gravlax.explorer.feature-resolution.v1",
            "read_only": true,
            "resolved": resolved,
        }))
    }

    fn preview_plan(&self, request: &PlannerRequest) -> Result<Value> {
        let plan = build_planner_model(self, request)?;
        let mut yaml = serde_yaml::to_string(&plan).context("serializing Explorer plan as YAML")?;
        if !yaml.ends_with('\n') {
            yaml.push('\n');
        }
        let mut plan_json =
            serde_json::to_string_pretty(&plan).context("serializing Explorer plan as JSON")?;
        plan_json.push('\n');
        let resolved = crate::plancmd::resolve_plan_model(
            plan,
            &format!("explorer-{}.yaml", request.view.as_str()),
            self.project()?,
        )?;
        let step = resolved
            .steps
            .first()
            .context("Explorer plan resolver returned no step")?;
        if resolved.steps.len() != 1 {
            bail!("Explorer plan resolver returned more than one step");
        }
        let export_arguments = logical_export_arguments(step)?;
        let command = cli_export(&export_arguments);
        let python = python_export(&export_arguments)?;
        let (route, scope, exclusions) = planner_route(request, step)?;
        let comparison_labels = step
            .annotation_comparison
            .as_ref()
            .map(|comparison| {
                let annotation_a = step
                    .annotation_inputs
                    .iter()
                    .find(|input| input.role == crate::plancmd::ResolvedAnnotationRole::A)
                    .context("resolved annotation comparison has no annotation A input")?;
                let annotation_b = step
                    .annotation_inputs
                    .iter()
                    .find(|input| input.role == crate::plancmd::ResolvedAnnotationRole::B)
                    .context("resolved annotation comparison has no annotation B input")?;
                Ok::<_, anyhow::Error>((
                    annotation_a.annotation.clone(),
                    annotation_b.annotation.clone(),
                    comparison.assembly.clone(),
                ))
            })
            .transpose()?;
        let coordinate_summary =
            if let Some((annotation_a, annotation_b, assembly)) = &comparison_labels {
                format!("A: {annotation_a} -> B: {annotation_b} · {assembly}")
            } else if let Some(intent) = &step.biological_intent {
                intent.locus.clone()
            } else if let Some(junction) = &request.junction {
                junction.clone()
            } else {
                format!(
                    "include [{}]; exclude [{}]",
                    request.include.join(", "),
                    request.exclude.join(", ")
                )
            };
        let title = if let Some((annotation_a, annotation_b, _)) = &comparison_labels {
            format!("{annotation_a} -> {annotation_b} · annotation comparison")
        } else {
            let subject = step
                .biological_intent
                .as_ref()
                .map(|intent| {
                    intent
                        .display_name
                        .clone()
                        .unwrap_or_else(|| intent.stable_id.clone())
                })
                .unwrap_or_else(|| coordinate_summary.clone());
            format!("{} · {}", subject, request.view.as_str())
        };
        Ok(json!({
            "schema": "gravlax.explorer.plan-preview.v1",
            "read_only": true,
            "view": request.view,
            "title": title,
            "coordinate_summary": coordinate_summary,
            "resolved_intent": step.biological_intent,
            "route": {
                "path": route,
                "scope": scope,
                "exclusions": exclusions,
            },
            "exports": {
                "yaml": yaml,
                "json": plan_json,
                "command": command,
                "python": python,
            },
            "provenance": {
                "resolved_plan_schema_version": resolved.schema_version,
                "plan_schema_version": resolved.plan_schema_version,
                "plan_source_digest": resolved.source_digest,
                "project_manifest_digest": resolved.project_manifest_digest,
                "producer": resolved.producer,
                "resources": resolved.resources,
            },
            "resolved_step": step,
        }))
    }
}

#[cfg(unix)]
fn same_artifact_file(left: &std::fs::Metadata, right: &std::fs::Metadata) -> bool {
    use std::os::unix::fs::MetadataExt;

    left.dev() == right.dev()
        && left.ino() == right.ino()
        && left.len() == right.len()
        && left.modified().ok() == right.modified().ok()
}

#[cfg(not(unix))]
fn same_artifact_file(left: &std::fs::Metadata, right: &std::fs::Metadata) -> bool {
    left.len() == right.len() && left.modified().ok() == right.modified().ok()
}

fn media_type(path: &Path) -> &'static str {
    match path.extension().and_then(|value| value.to_str()) {
        Some("json" | "jsonl") => "application/json; charset=utf-8",
        Some("yaml" | "yml") => "application/yaml; charset=utf-8",
        Some("tsv" | "csv" | "txt" | "md" | "mtx") => "text/plain; charset=utf-8",
        Some("gz") => "application/gzip",
        _ => "application/octet-stream",
    }
}

fn percent_encode(value: &str) -> String {
    let mut out = String::new();
    for byte in value.bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b'~') {
            out.push(byte as char);
        } else {
            out.push_str(&format!("%{byte:02X}"));
        }
    }
    out
}

fn percent_decode(value: &str) -> Result<String> {
    let bytes = value.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut index = 0;
    while index < bytes.len() {
        match bytes[index] {
            b'%' => {
                if index + 2 >= bytes.len() {
                    bail!("truncated percent escape");
                }
                let high = (bytes[index + 1] as char)
                    .to_digit(16)
                    .context("invalid percent escape")?;
                let low = (bytes[index + 2] as char)
                    .to_digit(16)
                    .context("invalid percent escape")?;
                out.push((high * 16 + low) as u8);
                index += 3;
            }
            b'+' => {
                out.push(b' ');
                index += 1;
            }
            byte => {
                out.push(byte);
                index += 1;
            }
        }
    }
    if out.contains(&0) {
        bail!("NUL is not allowed in a query value");
    }
    String::from_utf8(out).context("query value is not UTF-8")
}

#[derive(Debug)]
struct Request {
    method: String,
    target: String,
    host: String,
}

fn read_request(stream: &mut TcpStream) -> Result<Request> {
    stream.set_read_timeout(Some(Duration::from_secs(5)))?;
    let mut bytes = Vec::new();
    let mut chunk = [0u8; 1024];
    while !bytes.windows(4).any(|window| window == b"\r\n\r\n") {
        let read = stream.read(&mut chunk)?;
        if read == 0 {
            bail!("connection closed before the HTTP headers completed");
        }
        bytes.extend_from_slice(&chunk[..read]);
        if bytes.len() > MAX_REQUEST_BYTES {
            bail!("HTTP headers exceed {MAX_REQUEST_BYTES} bytes");
        }
    }
    let text = std::str::from_utf8(&bytes).context("HTTP headers are not UTF-8")?;
    let mut lines = text.split("\r\n");
    let mut request_line = lines
        .next()
        .context("missing HTTP request line")?
        .split_whitespace();
    let method = request_line
        .next()
        .context("missing HTTP method")?
        .to_owned();
    let target = request_line
        .next()
        .context("missing HTTP target")?
        .to_owned();
    let version = request_line.next().context("missing HTTP version")?;
    if request_line.next().is_some() || !matches!(version, "HTTP/1.0" | "HTTP/1.1") {
        bail!("unsupported HTTP request line");
    }
    if !target.starts_with('/') || target.starts_with("//") {
        bail!("Explorer accepts only origin-form HTTP targets");
    }
    let mut host = None;
    for line in lines {
        if line.is_empty() {
            break;
        }
        let (name, value) = line.split_once(':').context("malformed HTTP header")?;
        if name.eq_ignore_ascii_case("host") {
            if host.is_some() {
                bail!("duplicate Host header");
            }
            host = Some(value.trim().to_owned());
        }
        if name.eq_ignore_ascii_case("content-length") && value.trim() != "0" {
            bail!("Explorer does not accept request bodies");
        }
        if name.eq_ignore_ascii_case("transfer-encoding") {
            bail!("Explorer does not accept transfer-encoded requests");
        }
    }
    Ok(Request {
        method,
        target,
        host: host.unwrap_or_default(),
    })
}

fn valid_loopback_host(host: &str) -> bool {
    let host = host.trim().to_ascii_lowercase();
    host == "localhost"
        || host.starts_with("localhost:")
        || host == "127.0.0.1"
        || host.starts_with("127.0.0.1:")
        || host == "[::1]"
        || host.starts_with("[::1]:")
}

enum ResponseBody {
    Bytes(Vec<u8>),
    File(File),
}

struct Response {
    status: u16,
    reason: &'static str,
    media_type: &'static str,
    headers: Vec<(String, String)>,
    length: u64,
    body: ResponseBody,
}

impl Response {
    fn bytes(status: u16, reason: &'static str, media_type: &'static str, body: Vec<u8>) -> Self {
        Self {
            status,
            reason,
            media_type,
            headers: Vec::new(),
            length: body.len() as u64,
            body: ResponseBody::Bytes(body),
        }
    }

    fn json(value: Value) -> Result<Self> {
        Ok(Self::bytes(
            200,
            "OK",
            "application/json; charset=utf-8",
            serde_json::to_vec_pretty(&value)?,
        ))
    }

    fn error(status: u16, reason: &'static str, message: impl Into<String>) -> Self {
        let body = serde_json::to_vec(&json!({
            "schema": "gravlax.explorer.error.v1",
            "error": {"status": status, "message": message.into()},
        }))
        .expect("JSON error serialization cannot fail");
        Self::bytes(status, reason, "application/json; charset=utf-8", body)
    }
}

fn explorer_error_response(error: &anyhow::Error) -> Response {
    use anno::intent::ResolutionError;

    let Some(resolution) = error.downcast_ref::<ResolutionError>() else {
        return Response::error(400, "Bad Request", format!("{error:#}"));
    };
    let (status, reason, code, data) = match resolution {
        ResolutionError::Ambiguous { query, candidates } => (
            409,
            "Conflict",
            "ambiguous-identifier",
            json!({"query": query, "candidates": candidates}),
        ),
        ResolutionError::NotFound { query, identity } => (
            404,
            "Not Found",
            "identifier-not-found",
            json!({"query": query, "identity": identity}),
        ),
        ResolutionError::IdentifierMetadataUnavailable {
            query,
            identity,
            missing,
        } => (
            422,
            "Unprocessable Content",
            "identifier-metadata-unavailable",
            json!({"query": query, "identity": identity, "missing": missing}),
        ),
        ResolutionError::EmptyIdentifier => (400, "Bad Request", "empty-identifier", Value::Null),
    };
    let body = serde_json::to_vec(&json!({
        "schema": "gravlax.explorer.error.v1",
        "error": {
            "status": status,
            "code": code,
            "message": resolution.to_string(),
            "data": data,
        },
    }))
    .expect("JSON resolution-error serialization cannot fail");
    Response::bytes(status, reason, "application/json; charset=utf-8", body)
}

fn target_parts(target: &str) -> (&str, Option<&str>) {
    target
        .split_once('?')
        .map_or((target, None), |(path, query)| (path, Some(query)))
}

type QueryParameters = BTreeMap<String, Vec<String>>;

fn query_parameters(query: Option<&str>) -> Result<QueryParameters> {
    let query = query.unwrap_or_default();
    if query.len() > MAX_QUERY_BYTES {
        bail!("Explorer query exceeds {MAX_QUERY_BYTES} bytes");
    }
    let mut out = QueryParameters::new();
    let mut fields = 0usize;
    for field in query.split('&').filter(|field| !field.is_empty()) {
        fields += 1;
        if fields > MAX_QUERY_FIELDS {
            bail!("Explorer query exceeds {MAX_QUERY_FIELDS} fields");
        }
        let (key, value) = field.split_once('=').unwrap_or((field, ""));
        let key = percent_decode(key)?;
        let value = percent_decode(value)?;
        if key.is_empty() || key.len() > MAX_QUERY_KEY_BYTES {
            bail!("query keys must contain 1-{MAX_QUERY_KEY_BYTES} UTF-8 bytes");
        }
        if value.len() > MAX_QUERY_VALUE_BYTES {
            bail!("query value `{key}` exceeds {MAX_QUERY_VALUE_BYTES} UTF-8 bytes");
        }
        out.entry(key).or_default().push(value);
    }
    Ok(out)
}

fn allow_query_keys(parameters: &QueryParameters, allowed: &[&str]) -> Result<()> {
    let allowed = allowed.iter().copied().collect::<BTreeSet<_>>();
    if let Some(key) = parameters
        .keys()
        .find(|key| !allowed.contains(key.as_str()))
    {
        bail!("unknown Explorer query parameter `{key}`");
    }
    Ok(())
}

fn one_query_parameter(parameters: &QueryParameters, name: &str) -> Result<Option<String>> {
    let Some(values) = parameters.get(name) else {
        return Ok(None);
    };
    if values.len() != 1 {
        bail!("duplicate {name} query parameter");
    }
    Ok(Some(values[0].clone()))
}

fn required_query_parameter(parameters: &QueryParameters, name: &str) -> Result<String> {
    let value = one_query_parameter(parameters, name)?
        .with_context(|| format!("missing {name} query parameter"))?;
    if value.trim().is_empty() {
        bail!("{name} query parameter must not be empty");
    }
    Ok(value)
}

fn repeated_query_parameter(parameters: &QueryParameters, name: &str) -> Vec<String> {
    parameters.get(name).cloned().unwrap_or_default()
}

fn feature_lookup_from_query(parameters: &QueryParameters) -> Result<FeatureLookup> {
    allow_query_keys(
        parameters,
        &["annotation", "feature", "assembly", "annotation-label"],
    )?;
    Ok(FeatureLookup {
        annotation: required_query_parameter(parameters, "annotation")?,
        feature: required_query_parameter(parameters, "feature")?,
        assembly: one_query_parameter(parameters, "assembly")?,
        annotation_label: one_query_parameter(parameters, "annotation-label")?,
    })
}

fn planner_request_from_query(parameters: &QueryParameters) -> Result<PlannerRequest> {
    allow_query_keys(
        parameters,
        &[
            "view",
            "archive",
            "annotation-a",
            "annotation-b",
            "gene-key",
            "solo-strand",
            "max-molecule-witnesses",
            "max-row-transitions-per-molecule",
            "allow-identical",
            "annotation",
            "feature",
            "assembly",
            "annotation-label",
            "junction",
            "include",
            "exclude",
            "cells",
            "groups",
            "aggregation",
            "site-gap",
            "strand",
        ],
    )?;
    let view = PlannerView::parse(&required_query_parameter(parameters, "view")?)?;
    let cells = one_query_parameter(parameters, "cells")?;
    let groups = one_query_parameter(parameters, "groups")?;
    if cells.is_some() && groups.is_some() {
        bail!("cells and groups are mutually exclusive query scopes");
    }
    let aggregation = match one_query_parameter(parameters, "aggregation")?.as_deref() {
        None | Some("auto") => crate::plancmd::QueryAggregation::Auto,
        Some("cell") => crate::plancmd::QueryAggregation::Cell,
        Some("group") => crate::plancmd::QueryAggregation::Group,
        Some("bulk") => crate::plancmd::QueryAggregation::Bulk,
        Some(_) => bail!("aggregation must be auto, cell, group, or bulk"),
    };
    let annotation = one_query_parameter(parameters, "annotation")?;
    let annotation_a = one_query_parameter(parameters, "annotation-a")?;
    let annotation_b = one_query_parameter(parameters, "annotation-b")?;
    let comparison_gene_key = match one_query_parameter(parameters, "gene-key")?.as_deref() {
        None | Some("unversioned") => crate::plancmd::ComparisonGeneKey::Unversioned,
        Some("exact") => crate::plancmd::ComparisonGeneKey::Exact,
        Some(_) => bail!("gene-key must be unversioned or exact"),
    };
    let comparison_solo_strand = match one_query_parameter(parameters, "solo-strand")?.as_deref() {
        None | Some("forward") => crate::plancmd::SoloStrand::Forward,
        Some("reverse") => crate::plancmd::SoloStrand::Reverse,
        Some("unstranded") => crate::plancmd::SoloStrand::Unstranded,
        Some(_) => bail!("solo-strand must be forward, reverse, or unstranded"),
    };
    let max_molecule_witnesses = one_query_parameter(parameters, "max-molecule-witnesses")?
        .map(|value| {
            value
                .parse::<usize>()
                .context("max-molecule-witnesses must be a nonnegative integer")
        })
        .transpose()?
        .unwrap_or(10_000);
    let max_row_transitions_per_molecule =
        one_query_parameter(parameters, "max-row-transitions-per-molecule")?
            .map(|value| {
                value
                    .parse::<usize>()
                    .context("max-row-transitions-per-molecule must be a nonnegative integer")
            })
            .transpose()?
            .unwrap_or(32);
    let allow_identical = match one_query_parameter(parameters, "allow-identical")?.as_deref() {
        None | Some("false") => false,
        Some("true") => true,
        Some(_) => bail!("allow-identical must be true or false"),
    };
    let feature = one_query_parameter(parameters, "feature")?;
    let assembly = one_query_parameter(parameters, "assembly")?;
    let annotation_label = one_query_parameter(parameters, "annotation-label")?;
    let junction = one_query_parameter(parameters, "junction")?;
    let include = repeated_query_parameter(parameters, "include");
    let exclude = repeated_query_parameter(parameters, "exclude");
    let site_gap = one_query_parameter(parameters, "site-gap")?
        .map(|value| value.parse::<u32>().context("site-gap must be an integer"))
        .transpose()?
        .unwrap_or(24);
    let strand = one_query_parameter(parameters, "strand")?;
    if !(1..=10_000).contains(&site_gap) {
        bail!("site-gap must be between 1 and 10000 bp");
    }
    if strand
        .as_deref()
        .is_some_and(|value| !matches!(value, "+" | "-"))
    {
        bail!("strand must be + or -");
    }
    let feature_view = matches!(
        view,
        PlannerView::Region | PlannerView::SpliceEvents | PlannerView::TerminalBoundary
    );
    let comparison_view = view == PlannerView::AnnotationComparison;
    if comparison_view {
        if annotation_a
            .as_ref()
            .is_none_or(|value| value.trim().is_empty())
            || annotation_b
                .as_ref()
                .is_none_or(|value| value.trim().is_empty())
        {
            bail!("annotation-comparison view requires annotation-a and annotation-b");
        }
        if cells.is_some() || groups.is_some() || parameters.contains_key("aggregation") {
            bail!("annotation-comparison view always compares the full archive and does not accept cells, groups, or aggregation");
        }
    } else if [
        "annotation-a",
        "annotation-b",
        "gene-key",
        "solo-strand",
        "max-molecule-witnesses",
        "max-row-transitions-per-molecule",
        "allow-identical",
    ]
    .iter()
    .any(|name| parameters.contains_key(*name))
    {
        bail!(
            "{} view does not accept annotation-comparison fields",
            view.as_str()
        );
    }
    if feature_view && (annotation.is_none() || feature.is_none()) {
        bail!("{} view requires annotation and feature", view.as_str());
    }
    if !feature_view
        && (annotation.is_some()
            || feature.is_some()
            || assembly.is_some()
            || annotation_label.is_some())
    {
        bail!(
            "{} view does not accept feature-resolution fields",
            view.as_str()
        );
    }
    match view {
        PlannerView::AnnotationComparison => {}
        PlannerView::Junction => {
            if junction
                .as_ref()
                .is_none_or(|value| value.trim().is_empty())
            {
                bail!("junction view requires one exact junction");
            }
            if !include.is_empty() || !exclude.is_empty() {
                bail!("junction view does not accept junction-set fields");
            }
        }
        PlannerView::JunctionSet => {
            if include.is_empty() || exclude.is_empty() {
                bail!("junction-set view requires at least one include and one exclude junction");
            }
            if junction.is_some() {
                bail!("junction-set view does not accept the junction field");
            }
        }
        _ => {
            if junction.is_some() || !include.is_empty() || !exclude.is_empty() {
                bail!("{} view does not accept junction fields", view.as_str());
            }
        }
    }
    if view != PlannerView::TerminalBoundary
        && (parameters.contains_key("site-gap") || strand.is_some())
    {
        bail!("site-gap and strand are only valid for terminal-boundary view");
    }
    Ok(PlannerRequest {
        view,
        archive: required_query_parameter(parameters, "archive")?,
        annotation_a,
        annotation_b,
        comparison_gene_key,
        comparison_solo_strand,
        max_molecule_witnesses,
        max_row_transitions_per_molecule,
        allow_identical,
        annotation,
        feature,
        assembly,
        annotation_label,
        junction,
        include,
        exclude,
        cells,
        groups,
        aggregation,
        site_gap,
        strand,
    })
}

fn route(request: &Request, backend: &dyn ExplorerBackend) -> Result<Response> {
    if !valid_loopback_host(&request.host) {
        return Ok(Response::error(
            421,
            "Misdirected Request",
            "Explorer accepts only loopback Host headers",
        ));
    }
    if !matches!(request.method.as_str(), "GET" | "HEAD") {
        let mut response = Response::error(405, "Method Not Allowed", "Explorer is read-only");
        response.headers.push(("Allow".into(), "GET, HEAD".into()));
        return Ok(response);
    }
    let (path, query) = target_parts(&request.target);
    let parameters = query_parameters(query)?;
    match path {
        "/" | "/index.html" => {
            allow_query_keys(&parameters, &[])?;
            Ok(Response::bytes(
                200,
                "OK",
                "text/html; charset=utf-8",
                EXPLORER_HTML.to_vec(),
            ))
        }
        "/favicon.ico" => {
            allow_query_keys(&parameters, &[])?;
            Ok(Response::bytes(
                204,
                "No Content",
                "image/x-icon",
                Vec::new(),
            ))
        }
        "/api/v1/health" => {
            allow_query_keys(&parameters, &[])?;
            Response::json(json!({
                "schema": "gravlax.explorer.health.v1",
                "ok": true,
                "read_only": true,
            }))
        }
        "/api/v1/overview" => {
            allow_query_keys(&parameters, &[])?;
            Response::json(backend.overview()?)
        }
        "/api/v1/artifacts" => {
            allow_query_keys(&parameters, &[])?;
            Response::json(json!({
                "schema": "gravlax.explorer.artifacts.v1",
                "artifacts": backend.artifacts()?.iter().map(ArtifactSummary::value).collect::<Vec<_>>(),
            }))
        }
        "/api/v1/artifact" => {
            allow_query_keys(&parameters, &["path"])?;
            let id = required_query_parameter(&parameters, "path")?;
            let artifact = backend.open_artifact(&id)?;
            let mut response = Response {
                status: 200,
                reason: "OK",
                media_type: artifact.summary.media_type,
                headers: vec![(
                    "Content-Disposition".into(),
                    format!(
                        "inline; filename*=UTF-8''{}",
                        percent_encode(
                            Path::new(&artifact.summary.id)
                                .file_name()
                                .and_then(|name| name.to_str())
                                .unwrap_or("artifact")
                        )
                    ),
                )],
                length: artifact.summary.bytes,
                body: ResponseBody::File(artifact.file),
            };
            response.headers.push((
                "X-Gravlax-Artifact".into(),
                percent_encode(&artifact.summary.id),
            ));
            Ok(response)
        }
        "/api/v1/resolve" => {
            let lookup = feature_lookup_from_query(&parameters)?;
            Response::json(backend.resolve_feature(&lookup)?)
        }
        "/api/v1/plan-preview" => {
            let request = planner_request_from_query(&parameters)?;
            Response::json(backend.preview_plan(&request)?)
        }
        _ => Ok(Response::error(
            404,
            "Not Found",
            "Explorer route not found",
        )),
    }
}

fn write_response(stream: &mut TcpStream, mut response: Response, head: bool) -> Result<()> {
    stream.set_write_timeout(Some(Duration::from_secs(30)))?;
    let mut headers = format!(
        "HTTP/1.1 {} {}\r\nContent-Type: {}\r\nContent-Length: {}\r\nConnection: close\r\nX-Content-Type-Options: nosniff\r\nReferrer-Policy: no-referrer\r\nCache-Control: no-store\r\n",
        response.status, response.reason, response.media_type, response.length
    );
    if response.media_type.starts_with("text/html") {
        headers.push_str(
            "Content-Security-Policy: default-src 'self'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'\r\n",
        );
    }
    for (name, value) in response.headers {
        if !name.contains(['\r', '\n']) && !value.contains(['\r', '\n']) {
            headers.push_str(&format!("{name}: {value}\r\n"));
        }
    }
    headers.push_str("\r\n");
    stream.write_all(headers.as_bytes())?;
    if !head {
        match &mut response.body {
            ResponseBody::Bytes(bytes) => stream.write_all(bytes)?,
            ResponseBody::File(file) => {
                std::io::copy(file, stream)?;
            }
        }
    }
    stream.flush()?;
    Ok(())
}

fn handle_connection(mut stream: TcpStream, backend: &dyn ExplorerBackend) -> Result<()> {
    let request = match read_request(&mut stream) {
        Ok(request) => request,
        Err(error) => {
            write_response(
                &mut stream,
                Response::error(400, "Bad Request", error.to_string()),
                false,
            )?;
            return Ok(());
        }
    };
    let head = request.method == "HEAD";
    let response = match route(&request, backend) {
        Ok(response) => response,
        Err(error) => explorer_error_response(&error),
    };
    write_response(&mut stream, response, head)
}

fn serve(
    listener: TcpListener,
    backend: Arc<dyn ExplorerBackend>,
    request_limit: Option<usize>,
) -> Result<()> {
    let mut handled = 0usize;
    for connection in listener.incoming() {
        let stream = connection.context("accepting Explorer connection")?;
        // Binding to 127.0.0.1 is the primary boundary; repeat the peer check so that a future
        // listener refactor cannot silently widen the server's reach.
        if !stream.peer_addr()?.ip().is_loopback() {
            continue;
        }
        if let Err(error) = handle_connection(stream, backend.as_ref()) {
            eprintln!("Explorer request failed: {error:#}");
        }
        handled += 1;
        if request_limit.is_some_and(|limit| handled >= limit) {
            break;
        }
    }
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    let selected = args
        .project
        .unwrap_or(std::env::current_dir().context("finding the current directory")?);
    let backend = Arc::new(FilesystemBackend::discover(&selected)?);
    // No user-selectable host exists: Explorer is loopback-only by construction.
    let listener = TcpListener::bind((Ipv4Addr::LOCALHOST, args.port))
        .with_context(|| format!("binding Explorer to 127.0.0.1:{}", args.port))?;
    let address: SocketAddr = listener.local_addr()?;
    println!("Gravlax Explorer (read-only)");
    println!("  project: {}", backend.root().display());
    println!("  open:    http://{address}/");
    println!("Press Ctrl-C to stop.");
    serve(listener, backend, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Scratch(PathBuf);

    impl Scratch {
        fn new() -> Self {
            let nonce = std::time::SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let path = std::env::temp_dir().join(format!(
                "gravlax-explorer-test-{}-{nonce}",
                std::process::id()
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            std::fs::remove_dir_all(&self.0).ok();
        }
    }

    fn fixture() -> (Scratch, FilesystemBackend) {
        let scratch = Scratch::new();
        std::fs::write(
            scratch.0.join(PROJECT_MANIFEST),
            "schema_version: 1\nname: test\nresources:\n  sample:\n    kind: archive\n    path: data/sample.aie\n    assembly: GRCh38.p14\n  genes:\n    kind: annotation\n    path: data/genes.gtf\n    annotation_identity:\n      assembly: GRCh38.p14\n      annotation: TEST 1\n  genes_next:\n    kind: annotation\n    path: data/genes-next.gtf\n    annotation_identity:\n      assembly: GRCh38.p14\n      annotation: TEST 2\n  unidentified:\n    kind: annotation\n    path: data/genes.gtf\n  selected:\n    kind: cells\n    path: data/cells.txt\n  clusters:\n    kind: groups\n    path: data/groups.tsv\n",
        )
        .unwrap();
        std::fs::create_dir(scratch.0.join("data")).unwrap();
        std::fs::write(
            scratch.0.join("data/genes-next.gtf"),
            concat!(
                "chr1\tX\texon\t101\t220\t.\t+\t.\tgene_id \"ENSG1.5\"; transcript_id \"ENST1.3\"; exon_id \"ENSE1.8\"; gene_name \"ALPHA\";\n",
                "chr1\tX\texon\t301\t400\t.\t+\t.\tgene_id \"ENSG1.5\"; transcript_id \"ENST1.3\"; exon_id \"ENSE2.2\"; gene_name \"ALPHA\";\n",
            ),
        )
        .unwrap();
        evidence_io::format::SectionWriter::create(&scratch.0.join("data/sample.aie"), 1)
            .unwrap()
            .finish()
            .unwrap();
        std::fs::write(
            scratch.0.join("data/genes.gtf"),
            concat!(
                "chr1\tX\texon\t101\t200\t.\t+\t.\tgene_id \"ENSG1.4\"; transcript_id \"ENST1.2\"; exon_id \"ENSE1.7\"; gene_name \"ALPHA\";\n",
                "chr1\tX\texon\t301\t400\t.\t+\t.\tgene_id \"ENSG1.4\"; transcript_id \"ENST1.2\"; exon_id \"ENSE2.1\"; gene_name \"ALPHA\";\n",
                "chr2\tX\texon\t501\t600\t.\t-\t.\tgene_id \"ENSG2.1\"; transcript_id \"ENST2.8\"; exon_id \"ENSE3.1\"; gene_name \"DUP\";\n",
                "chr3\tX\texon\t701\t800\t.\t+\t.\tgene_id \"ENSG3.1\"; transcript_id \"ENST3.1\"; exon_id \"ENSE4.1\"; gene_name \"DUP\";\n",
            ),
        )
        .unwrap();
        std::fs::write(scratch.0.join("data/cells.txt"), "AAAC\n").unwrap();
        std::fs::write(scratch.0.join("data/groups.tsv"), "AAAC\tT\n").unwrap();
        std::fs::create_dir(scratch.0.join("plans")).unwrap();
        std::fs::write(
            scratch.0.join("plans/replay.yaml"),
            "schema_version: 1\nname: inspect-sample\nsteps:\n  - id: inspect\n    kind: inspect-archive\n    archive: sample\n    format: json\n",
        )
        .unwrap();
        std::fs::create_dir_all(scratch.0.join(".aie/resolved-plans")).unwrap();
        std::fs::write(
            scratch.0.join(".aie/resolved-plans/replay.json"),
            b"{\"schema\":\"gravlax.resolved-plan.v1\"}\n",
        )
        .unwrap();
        std::fs::create_dir(scratch.0.join("results")).unwrap();
        std::fs::write(scratch.0.join("results/report.tsv"), b"feature\tcount\n").unwrap();
        let backend = FilesystemBackend::discover(&scratch.0).unwrap();
        (scratch, backend)
    }

    fn request(method: &str, target: &str, host: &str) -> Request {
        Request {
            method: method.into(),
            target: target.into(),
            host: host.into(),
        }
    }

    fn planner(view: PlannerView) -> PlannerRequest {
        let feature_view = matches!(
            view,
            PlannerView::Region | PlannerView::SpliceEvents | PlannerView::TerminalBoundary
        );
        PlannerRequest {
            view,
            archive: "sample".to_owned(),
            annotation_a: (view == PlannerView::AnnotationComparison).then(|| "genes".to_owned()),
            annotation_b: (view == PlannerView::AnnotationComparison)
                .then(|| "genes_next".to_owned()),
            comparison_gene_key: crate::plancmd::ComparisonGeneKey::Unversioned,
            comparison_solo_strand: crate::plancmd::SoloStrand::Forward,
            max_molecule_witnesses: 10_000,
            max_row_transitions_per_molecule: 32,
            allow_identical: false,
            annotation: feature_view.then(|| "genes".to_owned()),
            feature: feature_view.then(|| "ALPHA".to_owned()),
            assembly: feature_view.then(|| "GRCh38.p14".to_owned()),
            annotation_label: feature_view.then(|| "TEST 1".to_owned()),
            junction: (view == PlannerView::Junction).then(|| "chr1:200-300".to_owned()),
            include: (view == PlannerView::JunctionSet)
                .then(|| vec!["chr1:200-300".to_owned()])
                .unwrap_or_default(),
            exclude: (view == PlannerView::JunctionSet)
                .then(|| vec!["chr1:200-350".to_owned()])
                .unwrap_or_default(),
            cells: None,
            groups: None,
            aggregation: crate::plancmd::QueryAggregation::Auto,
            site_gap: 24,
            strand: None,
        }
    }

    fn response_json(response: Response) -> Value {
        match response.body {
            ResponseBody::Bytes(bytes) => serde_json::from_slice(&bytes).unwrap(),
            ResponseBody::File(_) => panic!("expected in-memory JSON response"),
        }
    }

    fn file_snapshot(root: &Path) -> Vec<(PathBuf, Vec<u8>)> {
        fn visit(root: &Path, directory: &Path, out: &mut Vec<(PathBuf, Vec<u8>)>) {
            let mut entries = std::fs::read_dir(directory)
                .unwrap()
                .collect::<std::io::Result<Vec<_>>>()
                .unwrap();
            entries.sort_by_key(|entry| entry.file_name());
            for entry in entries {
                let kind = entry.file_type().unwrap();
                if kind.is_dir() {
                    visit(root, &entry.path(), out);
                } else if kind.is_file() {
                    out.push((
                        entry.path().strip_prefix(root).unwrap().to_path_buf(),
                        std::fs::read(entry.path()).unwrap(),
                    ));
                }
            }
        }
        let mut out = Vec::new();
        visit(root, root, &mut out);
        out
    }

    #[test]
    fn catalogue_contains_exact_project_and_plan_artifacts() {
        let (_scratch, backend) = fixture();
        let artifacts = backend.artifacts().unwrap();
        let ids: Vec<_> = artifacts
            .iter()
            .map(|artifact| artifact.id.as_str())
            .collect();
        assert_eq!(
            ids,
            [
                ".aie/resolved-plans/replay.json",
                "aie-project.yaml",
                "plans/replay.yaml",
                "results/report.tsv"
            ]
        );
        let plan = artifacts
            .iter()
            .find(|artifact| artifact.id == "plans/replay.yaml")
            .unwrap();
        assert_eq!(plan.status, "valid", "{:?}", plan.detail);
        assert_eq!(plan.data["steps"], 1);
        let mut opened = backend.open_artifact("plans/replay.yaml").unwrap();
        let mut bytes = Vec::new();
        opened.file.read_to_end(&mut bytes).unwrap();
        assert!(bytes.ends_with(b"format: json\n"));
    }

    #[test]
    fn overview_exposes_named_resource_roles_status_and_identity() {
        let (_scratch, backend) = fixture();
        let overview = backend.overview().unwrap();
        let resources = overview["project"]["resources"].as_array().unwrap();
        let annotation = resources
            .iter()
            .find(|resource| resource["name"] == "genes")
            .unwrap();
        assert_eq!(annotation["kind"], "annotation");
        assert_eq!(annotation["status"], "ok");
        assert_eq!(annotation["annotation_identity"]["assembly"], "GRCh38.p14");
        let archive = resources
            .iter()
            .find(|resource| resource["name"] == "sample")
            .unwrap();
        assert_eq!(archive["assembly"], "GRCh38.p14");
        assert!(archive["bytes"].as_u64().unwrap() > 0);
    }

    #[test]
    fn feature_lookup_is_exact_content_bound_and_fails_closed() {
        let (_scratch, backend) = fixture();
        let value = backend
            .resolve_feature(&FeatureLookup {
                annotation: "genes".to_owned(),
                feature: "transcript:ENST1".to_owned(),
                assembly: Some("GRCh38.p14".to_owned()),
                annotation_label: Some("TEST 1".to_owned()),
            })
            .unwrap();
        assert_eq!(value["resolved"]["stable_id"], "ENST1.2");
        assert_eq!(value["resolved"]["locus"], "chr1:100-400");
        assert_eq!(value["resolved"]["strand"], "forward");
        let digest = value["resolved"]["annotation_digest"].as_str().unwrap();
        assert!(digest.starts_with("blake3:"));
        assert_eq!(digest.len(), "blake3:".len() + 64);

        let ambiguous_error = backend
            .resolve_feature(&FeatureLookup {
                annotation: "genes".to_owned(),
                feature: "DUP".to_owned(),
                assembly: None,
                annotation_label: None,
            })
            .unwrap_err();
        let ambiguous = format!("{ambiguous_error:#}");
        assert!(ambiguous.contains("ambiguous"), "{ambiguous}");
        let ambiguity_response = explorer_error_response(&ambiguous_error);
        assert_eq!(ambiguity_response.status, 409);
        let ambiguity = response_json(ambiguity_response);
        assert_eq!(ambiguity["error"]["code"], "ambiguous-identifier");
        assert_eq!(
            ambiguity["error"]["data"]["candidates"]
                .as_array()
                .unwrap()
                .len(),
            2
        );

        let unidentified = backend
            .resolve_feature(&FeatureLookup {
                annotation: "unidentified".to_owned(),
                feature: "ALPHA".to_owned(),
                assembly: None,
                annotation_label: None,
            })
            .unwrap_err()
            .to_string();
        assert!(unidentified.contains("no assembly/release identity"));

        let mismatch = backend
            .resolve_feature(&FeatureLookup {
                annotation: "genes".to_owned(),
                feature: "ALPHA".to_owned(),
                assembly: Some("GRCh37".to_owned()),
                annotation_label: Some("TEST 1".to_owned()),
            })
            .unwrap_err()
            .to_string();
        assert!(mismatch.contains("GRCh38.p14"));
        assert!(mismatch.contains("GRCh37"));
    }

    #[test]
    fn typed_plan_exports_share_one_resolved_command_and_do_not_write() {
        let (scratch, backend) = fixture();
        let before = file_snapshot(&scratch.0);
        let value = backend.preview_plan(&planner(PlannerView::Region)).unwrap();
        assert_eq!(value["schema"], "gravlax.explorer.plan-preview.v1");
        assert_eq!(value["resolved_intent"]["stable_id"], "ENSG1.4");
        assert_eq!(value["resolved_intent"]["locus"], "chr1:100-400");
        assert_eq!(value["resolved_step"]["uniform_io"]["kind"], "result");
        assert_eq!(value["resolved_step"]["uniform_io"]["format"], "json");
        assert_eq!(
            value["resolved_step"]["uniform_io"]["publication"],
            "atomic-no-clobber-file"
        );
        assert!(value["exports"]["yaml"]
            .as_str()
            .unwrap()
            .contains("uniform_output:"));
        assert!(value["exports"]["yaml"]
            .as_str()
            .unwrap()
            .contains("results/explorer-region.json"));
        assert_eq!(
            value["resolved_intent"]["compatibility"][0]["status"],
            "verified"
        );
        assert!(
            !value["provenance"]["resources"]["sample"]["identity"]["scheme"]
                .as_str()
                .unwrap()
                .is_empty()
        );
        assert!(
            !value["provenance"]["resources"]["sample"]["identity"]["digest"]
                .as_str()
                .unwrap()
                .is_empty()
        );
        let json_plan: crate::plancmd::AnalysisPlan =
            serde_json::from_str(value["exports"]["json"].as_str().unwrap()).unwrap();
        let yaml_plan: crate::plancmd::AnalysisPlan =
            serde_yaml::from_str(value["exports"]["yaml"].as_str().unwrap()).unwrap();
        assert_eq!(
            serde_json::to_value(&json_plan).unwrap(),
            serde_json::to_value(yaml_plan).unwrap()
        );
        let recompiled = crate::plancmd::resolve_plan_model(
            json_plan,
            "recompiled-export.yaml",
            backend.project().unwrap(),
        )
        .unwrap();
        assert_eq!(
            serde_json::to_value(&recompiled.steps[0].args).unwrap(),
            value["resolved_step"]["args"]
        );
        let logical_arguments = logical_export_arguments(&recompiled.steps[0]).unwrap();
        let arguments = serde_json::to_string_pretty(&logical_arguments).unwrap();
        assert!(value["exports"]["python"]
            .as_str()
            .unwrap()
            .contains(&arguments));
        let expected_command = format!("aie {}\n", logical_arguments.join(" "));
        assert_eq!(value["exports"]["command"], expected_command);
        assert!(!value["exports"]["command"]
            .as_str()
            .unwrap()
            .contains(".aie-stage-"));
        assert!(value["exports"]["command"]
            .as_str()
            .unwrap()
            .contains("results/explorer-region.json"));
        assert_eq!(
            backend.preview_plan(&planner(PlannerView::Region)).unwrap(),
            value,
            "the same project and request must yield byte-stable exports and provenance"
        );
        assert_eq!(file_snapshot(&scratch.0), before);
    }

    #[test]
    fn annotation_comparison_preview_is_canonical_content_bound_and_read_only() {
        let (scratch, backend) = fixture();
        let before = file_snapshot(&scratch.0);
        let mut request = planner(PlannerView::AnnotationComparison);
        request.comparison_gene_key = crate::plancmd::ComparisonGeneKey::Exact;
        request.comparison_solo_strand = crate::plancmd::SoloStrand::Unstranded;
        request.max_molecule_witnesses = 0;
        request.max_row_transitions_per_molecule = 7;

        let value = backend.preview_plan(&request).unwrap();
        assert_eq!(value["view"], "annotation-comparison");
        assert_eq!(value["resolved_step"]["kind"], "compare-annotations");
        assert_eq!(
            value["resolved_step"]["output_schema_ids"],
            json!([
                "gravlax.annotation.compare.v1",
                "gravlax.annotation.compare.count-deltas.v1",
                "gravlax.annotation.compare.class-transitions.v1",
                "gravlax.annotation.compare.contributing-causes.v1",
                "gravlax.annotation.compare.witnesses.v1",
            ])
        );
        assert_eq!(
            value["resolved_step"]["annotation_comparison"]["gene_key"],
            "exact"
        );
        assert_eq!(
            value["resolved_step"]["annotation_comparison"]["solo_strand"],
            "unstranded"
        );
        assert!(
            value["resolved_step"]["annotation_comparison"]["final_count_delta_semantics"]
                .as_str()
                .unwrap()
                .contains("B-minus-A")
        );
        assert!(
            value["resolved_step"]["annotation_comparison"]["transition_evidence_semantics"]
                .as_str()
                .unwrap()
                .contains("not additive")
        );
        let inputs = value["resolved_step"]["annotation_inputs"]
            .as_array()
            .unwrap();
        assert_eq!(inputs.len(), 2);
        assert_eq!(inputs[0]["role"], "a");
        assert_eq!(inputs[0]["annotation"], "TEST 1");
        assert_eq!(inputs[1]["role"], "b");
        assert_eq!(inputs[1]["annotation"], "TEST 2");
        assert_ne!(
            inputs[0]["source_identity"]["digest"],
            inputs[1]["source_identity"]["digest"]
        );
        for input in inputs {
            assert_eq!(input["compatibility"][0]["status"], "verified");
        }
        assert!(value["route"]["path"]
            .as_str()
            .unwrap()
            .contains("full archive scan"));
        assert!(value["route"]["exclusions"]
            .as_str()
            .unwrap()
            .contains("not additive"));
        assert!(value["resolved_step"]["io_estimate"]["read_bytes_upper_bound"].is_null());
        assert!(value["exports"]["yaml"]
            .as_str()
            .unwrap()
            .contains("kind: compare-annotations"));
        assert!(value["exports"]["command"]
            .as_str()
            .unwrap()
            .contains("--max-molecule-witnesses 0"));
        assert!(!value["exports"]["command"]
            .as_str()
            .unwrap()
            .contains(".aie-stage-"));
        assert!(value["exports"]["python"]
            .as_str()
            .unwrap()
            .contains("compare-annotations"));
        assert!(value["exports"]["yaml"]
            .as_str()
            .unwrap()
            .contains("results/explorer-annotation-comparison.json"));

        let plan: crate::plancmd::AnalysisPlan =
            serde_json::from_str(value["exports"]["json"].as_str().unwrap()).unwrap();
        let recompiled = crate::plancmd::resolve_plan_model(
            plan,
            "recompiled-comparison.yaml",
            backend.project().unwrap(),
        )
        .unwrap();
        assert_eq!(
            serde_json::to_value(&recompiled.steps[0].args).unwrap(),
            value["resolved_step"]["args"]
        );
        assert_eq!(file_snapshot(&scratch.0), before);
    }

    #[test]
    fn every_planner_view_has_an_explicit_route_and_exclusion_contract() {
        let (_scratch, backend) = fixture();
        for view in [
            PlannerView::AnnotationComparison,
            PlannerView::Region,
            PlannerView::SpliceEvents,
            PlannerView::Junction,
            PlannerView::JunctionSet,
            PlannerView::TerminalBoundary,
        ] {
            let value = backend.preview_plan(&planner(view)).unwrap();
            assert!(!value["route"]["path"].as_str().unwrap().is_empty());
            assert!(!value["route"]["scope"].as_str().unwrap().is_empty());
            assert!(!value["route"]["exclusions"].as_str().unwrap().is_empty());
            assert_eq!(value["resolved_step"]["id"], "explore");
            assert!(value["exports"]["json"].is_string());
            if view != PlannerView::AnnotationComparison {
                let expected_schemas = match view {
                    PlannerView::Region => serde_json::json!([
                        "gravlax.query.region.result.v1",
                        "gravlax.query.region.counts.v1"
                    ]),
                    PlannerView::SpliceEvents => serde_json::json!([
                        "gravlax.query.events.result.v1",
                        "gravlax.query.events.events.v1",
                        "gravlax.query.events.components.v1",
                        "gravlax.query.events.counts.v1"
                    ]),
                    PlannerView::Junction => serde_json::json!([
                        "gravlax.query.junction.result.v1",
                        "gravlax.query.junction.counts.v1"
                    ]),
                    PlannerView::JunctionSet => serde_json::json!([
                        "gravlax.query.jset.result.v1",
                        "gravlax.query.jset.junctions.v1",
                        "gravlax.query.jset.counts.v1"
                    ]),
                    PlannerView::TerminalBoundary => serde_json::json!([
                        "gravlax.query.apa.result.v1",
                        "gravlax.query.apa.sites.v1"
                    ]),
                    PlannerView::AnnotationComparison => unreachable!(),
                };
                assert_eq!(
                    value["resolved_step"]["output_schema_ids"],
                    expected_schemas
                );
                assert_eq!(value["resolved_step"]["uniform_io"]["kind"], "result");
                assert_eq!(value["resolved_step"]["uniform_io"]["format"], "json");
                let destination = format!("results/explorer-{}.json", view.as_str());
                assert!(value["resolved_step"]["uniform_io"]["output"]
                    .as_str()
                    .unwrap()
                    .ends_with(&destination));
                assert!(value["exports"]["yaml"]
                    .as_str()
                    .unwrap()
                    .contains(&destination));
            }
        }
    }

    #[test]
    fn path_traversal_never_reaches_the_filesystem() {
        let (_scratch, backend) = fixture();
        assert!(backend.open_artifact("../aie-project.yaml").is_err());
        assert!(backend.open_artifact("plans/../../etc/passwd").is_err());
        assert!(backend.open_artifact("not-catalogued.txt").is_err());
    }

    #[cfg(unix)]
    #[test]
    fn symlinked_artifacts_are_never_catalogued_or_served() {
        use std::os::unix::fs::symlink;

        let (scratch, backend) = fixture();
        let outside = scratch.0.parent().unwrap().join(format!(
            "{}-outside",
            scratch.0.file_name().unwrap().to_string_lossy()
        ));
        std::fs::write(&outside, "secret\n").unwrap();
        symlink(&outside, scratch.0.join("results/escape.txt")).unwrap();
        let ids = backend
            .artifacts()
            .unwrap()
            .into_iter()
            .map(|artifact| artifact.id)
            .collect::<Vec<_>>();
        assert!(!ids.iter().any(|id| id == "results/escape.txt"));
        assert!(backend.open_artifact("results/escape.txt").is_err());
        std::fs::remove_file(outside).unwrap();
    }

    #[test]
    fn mutating_methods_and_non_loopback_hosts_are_rejected() {
        let (_scratch, backend) = fixture();
        let post = route(
            &request("POST", "/api/v1/artifacts", "localhost:8787"),
            &backend,
        )
        .unwrap();
        assert_eq!(post.status, 405);
        let foreign = route(&request("GET", "/", "example.org"), &backend).unwrap();
        assert_eq!(foreign.status, 421);
    }

    #[test]
    fn scientific_endpoints_are_get_only_bounded_and_deny_unknown_fields() {
        let (_scratch, backend) = fixture();
        let resolved = route(
            &request(
                "GET",
                "/api/v1/resolve?annotation=genes&feature=ALPHA&assembly=GRCh38.p14&annotation-label=TEST+1",
                "localhost:8787",
            ),
            &backend,
        )
        .unwrap();
        assert_eq!(resolved.status, 200);
        assert_eq!(response_json(resolved)["resolved"]["stable_id"], "ENSG1.4");
        let preview = route(
            &request(
                "HEAD",
                "/api/v1/plan-preview?view=junction&archive=sample&junction=chr1%3A200-300",
                "127.0.0.1:8787",
            ),
            &backend,
        )
        .unwrap();
        assert_eq!(preview.status, 200);
        assert_eq!(
            response_json(preview)["resolved_step"]["kind"],
            "query-junction"
        );

        assert!(route(
            &request(
                "GET",
                "/api/v1/resolve?annotation=genes&feature=ALPHA&feature=ENSG1.4",
                "localhost:8787",
            ),
            &backend,
        )
        .is_err());
        assert!(route(
            &request(
                "GET",
                "/api/v1/resolve?annotation=genes&feature=ALPHA&execute=true",
                "localhost:8787",
            ),
            &backend,
        )
        .is_err());
        let oversized = format!(
            "/api/v1/resolve?annotation=genes&feature={}",
            "A".repeat(MAX_QUERY_VALUE_BYTES + 1)
        );
        assert!(route(&request("GET", &oversized, "localhost:8787"), &backend).is_err());
    }

    #[test]
    fn annotation_comparison_query_contract_fails_closed() {
        let (_scratch, backend) = fixture();
        let response = route(
            &request(
                "GET",
                "/api/v1/plan-preview?view=annotation-comparison&archive=sample&annotation-a=genes&annotation-b=genes_next&gene-key=unversioned&solo-strand=forward&max-molecule-witnesses=0&max-row-transitions-per-molecule=0",
                "localhost:8787",
            ),
            &backend,
        )
        .unwrap();
        assert_eq!(response.status, 200);
        let preview = response_json(response);
        assert_eq!(preview["resolved_step"]["kind"], "compare-annotations");
        assert_eq!(
            preview["resolved_step"]["annotation_inputs"][0]["compatibility"][0]["status"],
            "verified"
        );

        for target in [
            "/api/v1/plan-preview?view=annotation-comparison&archive=sample&annotation-a=genes",
            "/api/v1/plan-preview?view=annotation-comparison&archive=sample&annotation-a=genes&annotation-b=genes_next&cells=selected",
            "/api/v1/plan-preview?view=annotation-comparison&archive=sample&annotation-a=genes&annotation-b=genes_next&max-molecule-witnesses=-1",
            "/api/v1/plan-preview?view=region&archive=sample&annotation=genes&feature=ALPHA&annotation-a=genes",
        ] {
            assert!(
                route(&request("GET", target, "localhost:8787"), &backend).is_err(),
                "comparison contract unexpectedly accepted {target}"
            );
        }

        let identical = "/api/v1/plan-preview?view=annotation-comparison&archive=sample&annotation-a=genes&annotation-b=genes";
        assert!(route(&request("GET", identical, "localhost:8787"), &backend).is_err());
        let intentional_control = format!("{identical}&allow-identical=true");
        assert_eq!(
            route(
                &request("GET", &intentional_control, "localhost:8787"),
                &backend,
            )
            .unwrap()
            .status,
            200
        );
    }

    #[test]
    fn embedded_ui_uses_safe_dom_text_paths_only() {
        let html = std::str::from_utf8(EXPLORER_HTML).unwrap();
        assert!(html.contains("/api/v1/plan-preview"));
        assert!(html.contains("/api/v1/resolve"));
        assert!(html.contains("textContent"));
        for unsafe_api in [
            "innerHTML",
            "outerHTML",
            "insertAdjacentHTML",
            "document.write(",
        ] {
            assert!(!html.contains(unsafe_api), "UI contains {unsafe_api}");
        }
    }

    #[test]
    fn artifact_route_preserves_bytes_and_type() {
        let (scratch, backend) = fixture();
        let response = route(
            &request(
                "GET",
                "/api/v1/artifact?path=plans%2Freplay.yaml",
                "127.0.0.1:8787",
            ),
            &backend,
        )
        .unwrap();
        assert_eq!(response.status, 200);
        assert_eq!(response.media_type, "application/yaml; charset=utf-8");
        assert_eq!(
            response.length,
            std::fs::metadata(scratch.0.join("plans/replay.yaml"))
                .unwrap()
                .len()
        );
    }

    #[test]
    fn server_listener_is_loopback_and_serves_one_request() {
        let (_scratch, backend) = fixture();
        let listener = TcpListener::bind((Ipv4Addr::LOCALHOST, 0)).unwrap();
        let address = listener.local_addr().unwrap();
        assert!(address.ip().is_loopback());
        let thread =
            std::thread::spawn(move || serve(listener, Arc::new(backend), Some(1)).unwrap());
        let mut client = TcpStream::connect(address).unwrap();
        client
            .write_all(
                format!(
                    "GET /api/v1/health HTTP/1.1\r\nHost: 127.0.0.1:{}\r\n\r\n",
                    address.port()
                )
                .as_bytes(),
            )
            .unwrap();
        let mut response = String::new();
        client.read_to_string(&mut response).unwrap();
        thread.join().unwrap();
        assert!(response.starts_with("HTTP/1.1 200 OK"));
        assert!(response.contains("gravlax.explorer.health.v1"));
    }

    #[test]
    fn percent_codec_roundtrips_artifact_ids() {
        let id = ".aie/resolved-plans/a plan.json";
        assert_eq!(percent_decode(&percent_encode(id)).unwrap(), id);
        assert!(percent_decode("%ZZ").is_err());
    }
}
