//! Project workspaces for repeatable `aie` analyses.
//!
//! A project is deliberately small: one human-editable, versioned manifest plus conventional
//! `plans/`, `results/`, and `.aie/resolved-plans/` directories. Internal resource paths are
//! relative to the project root so the workspace remains movable; explicitly external inputs are
//! stored as canonical absolute, read-only paths. Plan outputs always stay inside the project.

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Component, Path, PathBuf};

pub const PROJECT_FILE: &str = "aie-project.yaml";
pub const PROJECT_SCHEMA_VERSION: u32 = 1;

#[derive(Parser)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Create a portable project workspace.
    Init(InitArgs),
    /// Register a named input file in the project manifest.
    Add(AddArgs),
    /// Show the project and the resolved status of its resources.
    Show(ShowArgs),
}

#[derive(clap::Args)]
struct InitArgs {
    /// Directory to initialize.
    #[arg(default_value = ".")]
    directory: PathBuf,
    /// Human-readable project name (defaults to the directory name).
    #[arg(long)]
    name: Option<String>,
}

#[derive(clap::Args)]
struct AddArgs {
    /// Stable name used by plans, for example `brain-v1`.
    name: String,
    /// Existing input file. Outside files require `--external`.
    path: PathBuf,
    /// Resource type. `auto` infers it from the file name.
    #[arg(long, value_enum, default_value_t = ResourceKindArg::Auto)]
    kind: ResourceKindArg,
    /// Register a canonical absolute path outside the project as a read-only input. Internal
    /// resources remain relative and keep the project movable.
    #[arg(long)]
    external: bool,
    /// Project directory or manifest; otherwise search upward from the current directory.
    #[arg(long)]
    project: Option<PathBuf>,
    /// Replace an existing resource with this name.
    #[arg(long)]
    replace: bool,
    /// Reference assembly (for example `GRCh38.p14`). Annotation resources also require
    /// `--annotation-label`; archives, collections, and genomes may carry the assembly alone.
    #[arg(long)]
    assembly: Option<String>,
    /// Immutable annotation release/label (for example `GENCODE 49`). Must be supplied together
    /// with `--assembly` and may only be used for annotations.
    #[arg(long)]
    annotation_label: Option<String>,
}

#[derive(clap::Args)]
struct ShowArgs {
    /// Project directory or manifest; otherwise search upward from the current directory.
    #[arg(long)]
    project: Option<PathBuf>,
    /// Emit a versioned JSON document.
    #[arg(long)]
    json: bool,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ResourceKind {
    Archive,
    Collection,
    Annotation,
    Genome,
    Bam,
    Barcodes,
    Whitelist,
    Groups,
    Cells,
    Design,
    Metadata,
    File,
}

impl ResourceKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Archive => "archive",
            Self::Collection => "collection",
            Self::Annotation => "annotation",
            Self::Genome => "genome",
            Self::Bam => "bam",
            Self::Barcodes => "barcodes",
            Self::Whitelist => "whitelist",
            Self::Groups => "groups",
            Self::Cells => "cells",
            Self::Design => "design",
            Self::Metadata => "metadata",
            Self::File => "file",
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum ResourceKindArg {
    Auto,
    Archive,
    Collection,
    Annotation,
    Genome,
    Bam,
    Barcodes,
    Whitelist,
    Groups,
    Cells,
    Design,
    Metadata,
    File,
}

impl ResourceKindArg {
    fn resolve(self, path: &Path) -> ResourceKind {
        match self {
            Self::Auto => infer_resource_kind(path),
            Self::Archive => ResourceKind::Archive,
            Self::Collection => ResourceKind::Collection,
            Self::Annotation => ResourceKind::Annotation,
            Self::Genome => ResourceKind::Genome,
            Self::Bam => ResourceKind::Bam,
            Self::Barcodes => ResourceKind::Barcodes,
            Self::Whitelist => ResourceKind::Whitelist,
            Self::Groups => ResourceKind::Groups,
            Self::Cells => ResourceKind::Cells,
            Self::Design => ResourceKind::Design,
            Self::Metadata => ResourceKind::Metadata,
            Self::File => ResourceKind::File,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ProjectResource {
    pub kind: ResourceKind,
    pub path: PathBuf,
    #[serde(default, skip_serializing_if = "is_false")]
    pub external: bool,
    /// Optional declared assembly for coordinate-bearing non-annotation resources. This is a
    /// scientific compatibility label, not a substitute for the computed content identity.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
    /// Scientific identity required when a plan resolves a biological identifier. The exact
    /// content digest is intentionally computed from the resource bytes during plan resolution
    /// and is not duplicated as editable manifest metadata.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_identity: Option<ProjectAnnotationIdentity>,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ProjectAnnotationIdentity {
    pub assembly: String,
    pub annotation: String,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ProjectManifest {
    pub schema_version: u32,
    pub name: String,
    #[serde(default)]
    pub resources: BTreeMap<String, ProjectResource>,
}

#[derive(Clone, Debug)]
pub struct ProjectContext {
    /// Canonical project root.
    pub root: PathBuf,
    /// Canonical path to `aie-project.yaml`.
    pub manifest_path: PathBuf,
    pub manifest: ProjectManifest,
    /// BLAKE3 of the exact manifest bytes parsed into `manifest`.
    pub manifest_digest: String,
}

impl ProjectContext {
    /// Resolve a named resource and enforce its declared role. `file` is an intentional escape
    /// hatch for inputs whose more specific type is not known yet.
    pub fn resolve_resource(
        &self,
        name: &str,
        accepted: &[ResourceKind],
    ) -> Result<(ResourceKind, PathBuf)> {
        let resource = self
            .manifest
            .resources
            .get(name)
            .with_context(|| format!("project resource `{name}` is not registered"))?;
        if resource.kind != ResourceKind::File && !accepted.contains(&resource.kind) {
            let expected = accepted
                .iter()
                .map(|kind| kind.as_str())
                .collect::<Vec<_>>()
                .join(" or ");
            bail!(
                "project resource `{name}` is {}, but this field requires {expected}",
                resource.kind.as_str()
            );
        }
        let path = resolve_resource_path(&self.root, resource)
            .with_context(|| format!("resolving project resource `{name}`"))?;
        if !path.is_file() {
            bail!(
                "project resource `{name}` is not a file: {}",
                path.display()
            );
        }
        Ok((resource.kind, path))
    }

    /// Resolve a prospective output below the project root. Existing symlink ancestors are
    /// canonicalized before containment is checked.
    pub fn resolve_output(&self, path: &Path) -> Result<PathBuf> {
        resolve_project_output(&self.root, path)
    }
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Init(args) => run_init(args),
        Command::Add(args) => run_add(args),
        Command::Show(args) => run_show(args),
    }
}

fn run_init(args: InitArgs) -> Result<()> {
    fs::create_dir_all(&args.directory)
        .with_context(|| format!("creating project directory {}", args.directory.display()))?;
    let root = fs::canonicalize(&args.directory)
        .with_context(|| format!("resolving project directory {}", args.directory.display()))?;
    let manifest_path = root.join(PROJECT_FILE);
    if manifest_path.exists() {
        bail!("project already exists at {}", manifest_path.display());
    }
    let default_name = root
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("aie-project")
        .to_owned();
    let name = args.name.unwrap_or(default_name);
    validate_project_name(&name)?;
    let manifest = ProjectManifest {
        schema_version: PROJECT_SCHEMA_VERSION,
        name,
        resources: BTreeMap::new(),
    };
    for directory in ["plans", "results", ".aie/resolved-plans"] {
        ensure_project_directory(&root, Path::new(directory))?;
    }
    // Install the manifest last: its presence is the marker of a successfully initialized
    // workspace, while harmless empty directories can be retried.
    write_new_manifest(&manifest_path, &manifest)?;
    println!(
        "initialized project `{}` at {}",
        manifest.name,
        root.display()
    );
    println!("manifest: {}", manifest_path.display());
    println!("next: aie project add <name> <path>");
    Ok(())
}

fn run_add(args: AddArgs) -> Result<()> {
    validate_resource_id(&args.name)?;
    let discovered = find_project(args.project.as_deref())?;
    let _manifest_lock = lock_project_manifest(&discovered.root)?;
    // Reload only after acquiring the stable sidecar lock. Concurrent add processes therefore
    // merge against the latest committed manifest instead of dropping one update.
    let mut project = load_project(&discovered.manifest_path)?;
    let source = fs::canonicalize(&args.path)
        .with_context(|| format!("resolving resource {}", args.path.display()))?;
    if !source.is_file() {
        bail!("resource is not a file: {}", source.display());
    }
    let stored_path = if args.external {
        source.clone()
    } else {
        let relative = source.strip_prefix(&project.root).with_context(|| {
            format!(
                "resource {} is outside project {}; pass --external to register this canonical read-only input",
                source.display(),
                project.root.display()
            )
        })?;
        let relative = clean_relative_path(relative, "resource path")?;
        if is_reserved_project_path(&relative) {
            bail!("project metadata cannot be registered as an analysis resource");
        }
        relative
    };
    if project.manifest.resources.contains_key(&args.name) && !args.replace {
        bail!(
            "resource `{}` already exists; pass --replace to update it",
            args.name
        );
    }
    let kind = args.kind.resolve(&source);
    let (assembly, annotation_identity) = scientific_metadata_from_args(
        kind,
        args.assembly.as_deref(),
        args.annotation_label.as_deref(),
    )?;
    project.manifest.resources.insert(
        args.name.clone(),
        ProjectResource {
            kind,
            path: stored_path,
            external: args.external,
            assembly,
            annotation_identity,
        },
    );
    save_manifest(
        &project.manifest_path,
        &project.manifest,
        &project.manifest_digest,
    )?;
    println!(
        "registered `{}` as {} -> {}",
        args.name,
        kind.as_str(),
        source.display()
    );
    Ok(())
}

#[derive(Serialize)]
struct ShownProject<'a> {
    schema_version: u32,
    name: &'a str,
    root: &'a Path,
    manifest: &'a Path,
    resources: Vec<ShownResource<'a>>,
}

#[derive(Serialize)]
struct ShownResource<'a> {
    name: &'a str,
    kind: ResourceKind,
    path: &'a Path,
    external: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    assembly: Option<&'a str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    annotation_identity: Option<&'a ProjectAnnotationIdentity>,
    resolved_path: Option<PathBuf>,
    status: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
}

fn run_show(args: ShowArgs) -> Result<()> {
    let project = find_project(args.project.as_deref())?;
    let resources = project
        .manifest
        .resources
        .iter()
        .map(
            |(name, resource)| match resolve_resource_path(&project.root, resource) {
                Ok(path) if path.is_file() => ShownResource {
                    name,
                    kind: resource.kind,
                    path: &resource.path,
                    external: resource.external,
                    assembly: resource.assembly.as_deref(),
                    annotation_identity: resource.annotation_identity.as_ref(),
                    resolved_path: Some(path),
                    status: "ok",
                    error: None,
                },
                Ok(path) => ShownResource {
                    name,
                    kind: resource.kind,
                    path: &resource.path,
                    external: resource.external,
                    assembly: resource.assembly.as_deref(),
                    annotation_identity: resource.annotation_identity.as_ref(),
                    resolved_path: Some(path),
                    status: "invalid",
                    error: Some("not a file".to_owned()),
                },
                Err(error) => ShownResource {
                    name,
                    kind: resource.kind,
                    path: &resource.path,
                    external: resource.external,
                    assembly: resource.assembly.as_deref(),
                    annotation_identity: resource.annotation_identity.as_ref(),
                    resolved_path: None,
                    status: "invalid",
                    error: Some(format!("{error:#}")),
                },
            },
        )
        .collect::<Vec<_>>();
    let shown = ShownProject {
        schema_version: PROJECT_SCHEMA_VERSION,
        name: &project.manifest.name,
        root: &project.root,
        manifest: &project.manifest_path,
        resources,
    };
    if args.json {
        println!("{}", serde_json::to_string_pretty(&shown)?);
    } else {
        println!("project: {}", shown.name);
        println!("root: {}", shown.root.display());
        println!("manifest: {}", shown.manifest.display());
        println!("resources: {}", shown.resources.len());
        for resource in shown.resources {
            let detail = resource
                .resolved_path
                .as_deref()
                .unwrap_or(resource.path)
                .display();
            if let Some(error) = resource.error {
                println!(
                    "  {} ({}{}) -> {} [{}: {}]",
                    resource.name,
                    resource.kind.as_str(),
                    if resource.external { ", external" } else { "" },
                    detail,
                    resource.status,
                    error
                );
            } else {
                println!(
                    "  {} ({}{}) -> {} [{}]",
                    resource.name,
                    resource.kind.as_str(),
                    if resource.external { ", external" } else { "" },
                    detail,
                    resource.status
                );
                if let Some(identity) = resource.annotation_identity {
                    println!(
                        "    annotation identity: {} / {}",
                        identity.assembly, identity.annotation
                    );
                }
                if let Some(assembly) = resource.assembly {
                    println!("    assembly: {assembly}");
                }
            }
        }
    }
    Ok(())
}

/// Find and load a project. An explicit path may name either the manifest or its directory. With
/// no explicit path, the search walks upward from the current directory.
pub fn find_project(explicit: Option<&Path>) -> Result<ProjectContext> {
    let manifest_path = if let Some(path) = explicit {
        if path.is_dir() {
            path.join(PROJECT_FILE)
        } else {
            path.to_owned()
        }
    } else {
        return find_project_from(&std::env::current_dir().context("reading current directory")?);
    };
    load_project(&manifest_path)
}

/// Search `start` and its ancestors for a project manifest.
pub fn find_project_from(start: &Path) -> Result<ProjectContext> {
    let mut directory = if start.is_dir() {
        fs::canonicalize(start)
            .with_context(|| format!("resolving project search path {}", start.display()))?
    } else {
        let parent = start
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        fs::canonicalize(parent)
            .with_context(|| format!("resolving project search path {}", parent.display()))?
    };
    loop {
        let candidate = directory.join(PROJECT_FILE);
        if candidate.is_file() {
            return load_project(&candidate);
        }
        if !directory.pop() {
            bail!("no {PROJECT_FILE} found; run `aie project init` or pass --project");
        }
    }
}

/// Load a project from a manifest path. The project root and manifest path are canonicalized;
/// individual resources are resolved lazily so `project show` can report stale entries together.
pub fn load_project(manifest_path: &Path) -> Result<ProjectContext> {
    if !manifest_path.is_file() {
        bail!("project manifest not found: {}", manifest_path.display());
    }
    let manifest_path = fs::canonicalize(manifest_path)
        .with_context(|| format!("resolving project manifest {}", manifest_path.display()))?;
    if manifest_path.file_name().and_then(|name| name.to_str()) != Some(PROJECT_FILE) {
        bail!("project manifest must be named {PROJECT_FILE}");
    }
    let root = manifest_path
        .parent()
        .context("project manifest has no parent directory")?
        .to_owned();
    let bytes = fs::read(&manifest_path)
        .with_context(|| format!("reading project manifest {}", manifest_path.display()))?;
    let manifest_digest = blake3::hash(&bytes).to_hex().to_string();
    let manifest: ProjectManifest = serde_yaml::from_slice(&bytes)
        .with_context(|| format!("parsing project manifest {}", manifest_path.display()))?;
    if manifest.schema_version != PROJECT_SCHEMA_VERSION {
        bail!(
            "unsupported project schema version {}; this aie supports version {}",
            manifest.schema_version,
            PROJECT_SCHEMA_VERSION
        );
    }
    validate_project_name(&manifest.name)?;
    for (name, resource) in &manifest.resources {
        validate_resource_id(name)?;
        if resource.external {
            if !resource.path.is_absolute() {
                bail!("external resource `{name}` must use an absolute path");
            }
        } else {
            clean_relative_path(&resource.path, &format!("path for resource `{name}`"))?;
        }
        validate_scientific_metadata(
            resource.kind,
            resource.assembly.as_deref(),
            resource.annotation_identity.as_ref(),
        )
        .with_context(|| format!("validating project resource `{name}`"))?;
    }
    Ok(ProjectContext {
        root,
        manifest_path,
        manifest,
        manifest_digest,
    })
}

pub fn validate_resource_id(name: &str) -> Result<()> {
    let valid = !name.is_empty()
        && name.len() <= 64
        && name != "."
        && name != ".."
        && name
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphanumeric())
        && name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.'));
    if !valid {
        bail!(
            "invalid resource name `{name}`; use 1-64 ASCII letters, digits, '.', '_' or '-', starting with a letter or digit"
        );
    }
    Ok(())
}

pub fn clean_relative_path(path: &Path, label: &str) -> Result<PathBuf> {
    if path.as_os_str().is_empty() {
        bail!("{label} must not be empty");
    }
    let mut clean = PathBuf::new();
    for component in path.components() {
        match component {
            Component::Normal(part) => clean.push(part),
            Component::CurDir if clean.as_os_str().is_empty() => {}
            Component::CurDir => {}
            Component::ParentDir => bail!("{label} must not contain '..': {}", path.display()),
            Component::RootDir | Component::Prefix(_) => {
                bail!(
                    "{label} must be relative to the project: {}",
                    path.display()
                )
            }
        }
    }
    if clean.as_os_str().is_empty() {
        bail!("{label} must name a path below the project root");
    }
    Ok(clean)
}

/// Create or validate a project-contained directory one component at a time. Existing symlinks
/// are rejected even when they point back inside the project: metadata and output roots should
/// have one unambiguous physical location, and creation must not follow a link outside first.
pub fn ensure_project_directory(root: &Path, relative: &Path) -> Result<PathBuf> {
    let relative = clean_relative_path(relative, "project directory")?;
    let mut cursor = root.to_owned();
    for component in relative.components() {
        let Component::Normal(component) = component else {
            unreachable!("clean_relative_path returned only normal components");
        };
        cursor.push(component);
        match fs::symlink_metadata(&cursor) {
            Ok(metadata) => validate_directory_component(&cursor, &metadata)?,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                match fs::create_dir(&cursor) {
                    Ok(()) => {}
                    Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
                        let metadata = fs::symlink_metadata(&cursor).with_context(|| {
                            format!(
                                "inspecting concurrently created directory {}",
                                cursor.display()
                            )
                        })?;
                        validate_directory_component(&cursor, &metadata)?;
                    }
                    Err(error) => {
                        return Err(error)
                            .with_context(|| format!("creating directory {}", cursor.display()));
                    }
                }
            }
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("inspecting directory {}", cursor.display()));
            }
        }
    }
    let resolved = fs::canonicalize(&cursor)
        .with_context(|| format!("resolving directory {}", cursor.display()))?;
    if !resolved.starts_with(root) {
        bail!(
            "project directory escapes through a symlink: {}",
            cursor.display()
        );
    }
    Ok(cursor)
}

fn validate_directory_component(path: &Path, metadata: &fs::Metadata) -> Result<()> {
    if metadata.file_type().is_symlink() {
        bail!(
            "project directory must not be a symlink: {}",
            path.display()
        );
    }
    if !metadata.is_dir() {
        bail!(
            "project directory component is not a directory: {}",
            path.display()
        );
    }
    Ok(())
}

fn resolve_existing_project_path(root: &Path, relative: &Path) -> Result<PathBuf> {
    let relative = clean_relative_path(relative, "project path")?;
    let joined = root.join(relative);
    let resolved = fs::canonicalize(&joined)
        .with_context(|| format!("path does not exist: {}", joined.display()))?;
    if !resolved.starts_with(root) {
        bail!(
            "path escapes the project through a symlink: {}",
            joined.display()
        );
    }
    Ok(resolved)
}

fn resolve_resource_path(root: &Path, resource: &ProjectResource) -> Result<PathBuf> {
    if resource.external {
        if !resource.path.is_absolute() {
            bail!(
                "external resource path must be absolute: {}",
                resource.path.display()
            );
        }
        fs::canonicalize(&resource.path)
            .with_context(|| format!("external path does not exist: {}", resource.path.display()))
    } else {
        resolve_existing_project_path(root, &resource.path)
    }
}

pub(crate) fn resolve_project_output(root: &Path, relative: &Path) -> Result<PathBuf> {
    let relative = clean_relative_path(relative, "output path")?;
    if is_reserved_project_path(&relative) {
        bail!(
            "output path is reserved for project metadata: {}",
            relative.display()
        );
    }
    let output = root.join(relative);
    // `Path::exists` follows symlinks and therefore misses a broken symlink. Walk every existing
    // component with `symlink_metadata` first so neither a broken final link nor a link in a
    // not-yet-created subtree can redirect a writer outside the workspace.
    let components = output
        .strip_prefix(root)
        .expect("output was constructed below root")
        .iter()
        .collect::<Vec<_>>();
    let mut cursor = root.to_owned();
    for (index, component) in components.iter().enumerate() {
        cursor.push(component);
        match fs::symlink_metadata(&cursor) {
            Ok(metadata) if metadata.file_type().is_symlink() => {
                let resolved = fs::canonicalize(&cursor).with_context(|| {
                    format!("output contains a broken symlink: {}", cursor.display())
                })?;
                if !resolved.starts_with(root) {
                    bail!(
                        "output escapes the project through a symlink: {}",
                        output.display()
                    );
                }
                if index + 1 < components.len() && !resolved.is_dir() {
                    bail!("output parent is not a directory: {}", cursor.display());
                }
            }
            Ok(metadata) if index + 1 < components.len() && !metadata.is_dir() => {
                bail!("output parent is not a directory: {}", cursor.display());
            }
            Ok(_) => {}
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("inspecting output path {}", cursor.display()));
            }
        }
    }
    let mut ancestor = output.as_path();
    while !ancestor.exists() {
        ancestor = ancestor
            .parent()
            .context("output has no existing ancestor")?;
    }
    let resolved_ancestor = fs::canonicalize(ancestor)
        .with_context(|| format!("resolving output ancestor {}", ancestor.display()))?;
    if !resolved_ancestor.starts_with(root) {
        bail!(
            "output escapes the project through a symlink: {}",
            output.display()
        );
    }
    let unresolved_suffix = output
        .strip_prefix(ancestor)
        .expect("existing ancestor is a prefix of the requested output");
    let resolved_output = if unresolved_suffix.as_os_str().is_empty() {
        resolved_ancestor
    } else {
        resolved_ancestor.join(unresolved_suffix)
    };
    let resolved_relative = resolved_output
        .strip_prefix(root)
        .expect("resolved output was checked to remain below the project root");
    if is_reserved_project_path(resolved_relative) {
        bail!(
            "output path resolves into reserved project metadata: {} -> {}",
            output.display(),
            resolved_output.display()
        );
    }
    Ok(resolved_output)
}

fn is_reserved_project_path(relative: &Path) -> bool {
    relative == Path::new(PROJECT_FILE)
        || [".aie", ".git", ".hg", ".svn"]
            .iter()
            .any(|directory| relative.starts_with(directory))
}

fn validate_project_name(name: &str) -> Result<()> {
    if name.trim().is_empty() || name.len() > 200 || name.chars().any(char::is_control) {
        bail!("project name must be 1-200 printable characters");
    }
    Ok(())
}

fn scientific_metadata_from_args(
    kind: ResourceKind,
    assembly: Option<&str>,
    annotation: Option<&str>,
) -> Result<(Option<String>, Option<ProjectAnnotationIdentity>)> {
    let (declared_assembly, annotation_identity) = match kind {
        ResourceKind::Annotation => match (assembly, annotation) {
            (None, None) => (None, None),
            (Some(_), None) | (None, Some(_)) => {
                bail!("annotation resources require --assembly and --annotation-label together")
            }
            (Some(assembly), Some(annotation)) => (
                None,
                Some(ProjectAnnotationIdentity {
                    assembly: assembly.to_owned(),
                    annotation: annotation.to_owned(),
                }),
            ),
        },
        ResourceKind::Archive | ResourceKind::Collection | ResourceKind::Genome => {
            if annotation.is_some() {
                bail!("--annotation-label may only be used for annotation resources");
            }
            (assembly.map(str::to_owned), None)
        }
        _ => {
            if assembly.is_some() || annotation.is_some() {
                bail!(
                    "assembly metadata may only be attached to archive, collection, genome, or annotation resources"
                );
            }
            (None, None)
        }
    };
    validate_scientific_metadata(
        kind,
        declared_assembly.as_deref(),
        annotation_identity.as_ref(),
    )?;
    Ok((declared_assembly, annotation_identity))
}

fn validate_scientific_metadata(
    kind: ResourceKind,
    assembly: Option<&str>,
    identity: Option<&ProjectAnnotationIdentity>,
) -> Result<()> {
    if let Some(assembly) = assembly {
        if !matches!(
            kind,
            ResourceKind::Archive | ResourceKind::Collection | ResourceKind::Genome
        ) {
            bail!(
                "assembly metadata may only be attached to archive, collection, or genome resources"
            );
        }
        validate_identity_label(assembly, "assembly")?;
    }
    if let Some(identity) = identity {
        if kind != ResourceKind::Annotation {
            bail!("annotation identity metadata may only be attached to an annotation resource");
        }
        if assembly.is_some() {
            bail!("annotation assembly must be nested in annotation_identity, not duplicated");
        }
        validate_identity_label(&identity.assembly, "annotation assembly")?;
        validate_identity_label(&identity.annotation, "annotation label")?;
    } else if kind == ResourceKind::Annotation && assembly.is_some() {
        bail!("annotation resources require an annotation label with their assembly");
    }
    Ok(())
}

fn validate_identity_label(value: &str, label: &str) -> Result<()> {
    if value.trim() != value
        || value.is_empty()
        || value.len() > 200
        || value.chars().any(char::is_control)
    {
        bail!("{label} must be 1-200 printable characters without surrounding whitespace");
    }
    Ok(())
}

fn is_false(value: &bool) -> bool {
    !*value
}

fn infer_resource_kind(path: &Path) -> ResourceKind {
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if name.ends_with(".aicollection") {
        ResourceKind::Collection
    } else if name.ends_with(".aie") {
        ResourceKind::Archive
    } else if name.ends_with(".gtf")
        || name.ends_with(".gtf.gz")
        || name.ends_with(".gff")
        || name.ends_with(".gff3")
        || name.ends_with(".aic")
    {
        ResourceKind::Annotation
    } else if name.ends_with(".fa")
        || name.ends_with(".fa.gz")
        || name.ends_with(".fasta")
        || name.ends_with(".fasta.gz")
        || name.ends_with(".fna")
        || name.ends_with(".fna.gz")
    {
        ResourceKind::Genome
    } else if name.ends_with(".bam") {
        ResourceKind::Bam
    } else if name.contains("design") {
        ResourceKind::Design
    } else if name.contains("barcode") {
        ResourceKind::Barcodes
    } else if name.contains("whitelist") {
        ResourceKind::Whitelist
    } else if name.contains("group") {
        ResourceKind::Groups
    } else if name.contains("cell") {
        ResourceKind::Cells
    } else if name.ends_with(".tsv") || name.ends_with(".csv") || name.ends_with(".json") {
        ResourceKind::Metadata
    } else {
        ResourceKind::File
    }
}

fn manifest_yaml(manifest: &ProjectManifest) -> Result<Vec<u8>> {
    let mut text = String::from(
        "# aie project manifest; paths are relative unless explicitly marked external\n",
    );
    text.push_str(&serde_yaml::to_string(manifest).context("serializing project manifest")?);
    Ok(text.into_bytes())
}

fn write_new_manifest(path: &Path, manifest: &ProjectManifest) -> Result<()> {
    let bytes = manifest_yaml(manifest)?;
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .with_context(|| format!("creating project manifest {}", path.display()))?;
    file.write_all(&bytes)?;
    file.sync_all()?;
    Ok(())
}

struct ManifestLock {
    file: fs::File,
}

impl Drop for ManifestLock {
    fn drop(&mut self) {
        let _ = self.file.unlock();
    }
}

fn lock_project_manifest(root: &Path) -> Result<ManifestLock> {
    let directory = ensure_project_directory(root, Path::new(".aie"))?;
    let path = directory.join("project-manifest.lock");
    #[cfg(unix)]
    let file = {
        use std::os::unix::fs::OpenOptionsExt;
        OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .mode(0o600)
            .custom_flags(libc::O_NOFOLLOW | libc::O_CLOEXEC)
            .open(&path)
    };
    #[cfg(not(unix))]
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(&path);
    let file = file.with_context(|| format!("opening project manifest lock {}", path.display()))?;
    if !file.metadata()?.is_file() {
        bail!(
            "project manifest lock is not a regular file: {}",
            path.display()
        );
    }
    file.lock()
        .with_context(|| format!("locking project manifest {}", path.display()))?;
    Ok(ManifestLock { file })
}

fn save_manifest(path: &Path, manifest: &ProjectManifest, expected_digest: &str) -> Result<()> {
    let bytes = manifest_yaml(manifest)?;
    let file_name = path
        .file_name()
        .context("project manifest has no file name")?
        .to_string_lossy();
    for attempt in 0u32..1_000 {
        let temporary =
            path.with_file_name(format!(".{file_name}.tmp.{}-{attempt}", std::process::id()));
        let mut file = match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
        {
            Ok(file) => file,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error).context("creating temporary project manifest"),
        };
        if let Err(error) = file.write_all(&bytes).and_then(|_| file.sync_all()) {
            let _ = fs::remove_file(&temporary);
            return Err(error).context("writing temporary project manifest");
        }
        drop(file);
        let current = fs::read(path)
            .with_context(|| format!("re-reading project manifest {}", path.display()))?;
        let current_digest = blake3::hash(&current).to_hex().to_string();
        if current_digest != expected_digest {
            let _ = fs::remove_file(&temporary);
            bail!(
                "project manifest changed while it was being updated; reload and retry: {}",
                path.display()
            );
        }
        if let Err(error) = fs::rename(&temporary, path) {
            let _ = fs::remove_file(&temporary);
            return Err(error).with_context(|| format!("installing {}", path.display()));
        }
        return Ok(());
    }
    bail!("could not allocate a temporary project manifest");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_unsafe_relative_paths() {
        assert!(clean_relative_path(Path::new("../outside"), "test").is_err());
        assert!(clean_relative_path(Path::new("/absolute"), "test").is_err());
        assert_eq!(
            clean_relative_path(Path::new("results/./answer.json"), "test").unwrap(),
            Path::new("results/answer.json")
        );
    }

    #[test]
    fn project_metadata_is_not_an_output_target() {
        let root = fs::canonicalize(".").unwrap();
        assert!(resolve_project_output(&root, Path::new(".aie/resolved-plans/fake.json")).is_err());
        assert!(resolve_project_output(&root, Path::new(PROJECT_FILE)).is_err());
    }

    #[test]
    fn resource_kind_inference_handles_compressed_references() {
        assert_eq!(
            infer_resource_kind(Path::new("GRCh38.primary.fa.gz")),
            ResourceKind::Genome
        );
        assert_eq!(
            infer_resource_kind(Path::new("gencode.v48.gtf.gz")),
            ResourceKind::Annotation
        );
        assert_eq!(
            infer_resource_kind(Path::new("sample.aie")),
            ResourceKind::Archive
        );
        assert_eq!(
            infer_resource_kind(Path::new("atlas.aicollection")),
            ResourceKind::Collection
        );
    }

    #[test]
    fn stale_manifest_snapshots_cannot_replace_newer_updates() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gravlax-project-manifest-cas-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        let path = root.join(PROJECT_FILE);
        let initial = ProjectManifest {
            schema_version: PROJECT_SCHEMA_VERSION,
            name: "initial".to_owned(),
            resources: BTreeMap::new(),
        };
        write_new_manifest(&path, &initial).unwrap();
        let stale = load_project(&path).unwrap();

        let mut concurrent = initial.clone();
        concurrent.name = "concurrent".to_owned();
        fs::write(&path, manifest_yaml(&concurrent).unwrap()).unwrap();
        let mut proposed = stale.manifest;
        proposed.name = "stale-writer".to_owned();
        let error = save_manifest(&path, &proposed, &stale.manifest_digest)
            .unwrap_err()
            .to_string();
        assert!(error.contains("changed while it was being updated"));
        assert_eq!(load_project(&path).unwrap().manifest.name, "concurrent");
        fs::remove_dir_all(root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn broken_output_symlinks_are_rejected() {
        use std::os::unix::fs::symlink;

        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gravlax-project-path-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(root.join("results")).unwrap();
        let root = fs::canonicalize(root).unwrap();
        symlink("missing-target", root.join("results/out.json")).unwrap();
        let error = resolve_project_output(&root, Path::new("results/out.json"))
            .unwrap_err()
            .to_string();
        assert!(error.contains("broken symlink"));
        fs::remove_dir_all(&root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn in_project_output_aliases_resolve_to_one_physical_path() {
        use std::os::unix::fs::symlink;

        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gravlax-project-alias-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(root.join("results")).unwrap();
        let root = fs::canonicalize(root).unwrap();
        symlink("results", root.join("alias")).unwrap();
        let direct = resolve_project_output(&root, Path::new("results/out.gtf")).unwrap();
        let aliased = resolve_project_output(&root, Path::new("alias/out.gtf")).unwrap();
        assert_eq!(direct, aliased);
        fs::remove_dir_all(&root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn output_aliases_into_reserved_metadata_are_rejected() {
        use std::os::unix::fs::symlink;

        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "gravlax-project-reserved-alias-test-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir_all(root.join(".aie")).unwrap();
        let root = fs::canonicalize(root).unwrap();
        symlink(".aie", root.join("alias")).unwrap();
        let error = resolve_project_output(&root, Path::new("alias/evil.aic"))
            .unwrap_err()
            .to_string();
        assert!(error.contains("reserved project metadata"));
        fs::remove_dir_all(&root).unwrap();
    }
}
