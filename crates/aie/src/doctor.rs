//! Installation and input diagnostics for the user-facing CLI.
//!
//! `doctor` is deliberately observational.  Apart from a create/remove probe in the selected
//! workspace, it does not mutate projects or archives.  Archive checks use the same guarded
//! section reader as ordinary commands, so a successful check means the directory/root and the
//! small metadata sections were genuinely read rather than guessed from a file extension.

use anyhow::{bail, Context, Result};
use clap::Parser;
use evidence_io::format::{SectionReader, SEEKABLE_VERSION, VERSION};
use serde_json::{json, Value};
use std::collections::{BTreeSet, HashSet};
use std::ffi::OsStr;
use std::fs::OpenOptions;
use std::path::{Path, PathBuf};

const PROJECT_MANIFEST: &str = "aie-project.yaml";

#[derive(Parser, Debug)]
pub struct Args {
    /// Files to validate (`.aie` archives and `.aic` compiled annotations are inspected).
    #[arg(value_name = "PATH")]
    pub paths: Vec<PathBuf>,

    /// Project directory or `aie-project.yaml`; defaults to discovery from the current directory.
    #[arg(long, value_name = "PATH")]
    pub project: Option<PathBuf>,

    /// Read and authenticate every compressed archive payload, not only selected metadata.
    #[arg(long)]
    pub verify_content: bool,

    /// Treat warnings (such as missing optional alignment tools) as an unsuccessful diagnosis.
    #[arg(long)]
    pub strict: bool,

    /// Emit one stable JSON report.
    #[arg(long)]
    pub json: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Status {
    Pass,
    Warn,
    Fail,
}

impl Status {
    fn as_str(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Warn => "warn",
            Self::Fail => "fail",
        }
    }

    fn marker(self) -> &'static str {
        match self {
            Self::Pass => "ok",
            Self::Warn => "warning",
            Self::Fail => "failed",
        }
    }
}

#[derive(Debug)]
struct Check {
    id: String,
    status: Status,
    summary: String,
    detail: Option<String>,
    remedy: Option<String>,
    data: Value,
}

impl Check {
    fn new(id: impl Into<String>, status: Status, summary: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            status,
            summary: summary.into(),
            detail: None,
            remedy: None,
            data: Value::Null,
        }
    }

    fn detail(mut self, detail: impl Into<String>) -> Self {
        self.detail = Some(detail.into());
        self
    }

    fn remedy(mut self, remedy: impl Into<String>) -> Self {
        self.remedy = Some(remedy.into());
        self
    }

    fn data(mut self, data: Value) -> Self {
        self.data = data;
        self
    }

    fn value(&self) -> Value {
        json!({
            "id": self.id,
            "status": self.status.as_str(),
            "summary": self.summary,
            "detail": self.detail,
            "remedy": self.remedy,
            "data": self.data,
        })
    }
}

#[derive(Debug)]
struct Report {
    checks: Vec<Check>,
}

impl Report {
    fn failures(&self) -> usize {
        self.checks
            .iter()
            .filter(|check| check.status == Status::Fail)
            .count()
    }

    fn warnings(&self) -> usize {
        self.checks
            .iter()
            .filter(|check| check.status == Status::Warn)
            .count()
    }

    fn value(&self, strict: bool) -> Value {
        json!({
            "schema": "gravlax.doctor.v1",
            "version": env!("CARGO_PKG_VERSION"),
            "ok": self.failures() == 0 && (!strict || self.warnings() == 0),
            "strict": strict,
            "summary": {
                "passed": self.checks.len() - self.failures() - self.warnings(),
                "warnings": self.warnings(),
                "failures": self.failures(),
            },
            "checks": self.checks.iter().map(Check::value).collect::<Vec<_>>(),
        })
    }
}

fn executable_in_path(name: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    let candidates: Vec<&OsStr> = if cfg!(windows) {
        vec![
            OsStr::new(""),
            OsStr::new(".exe"),
            OsStr::new(".cmd"),
            OsStr::new(".bat"),
        ]
    } else {
        vec![OsStr::new("")]
    };
    for directory in std::env::split_paths(&path) {
        for suffix in &candidates {
            let mut file_name = name.to_owned();
            file_name.push_str(&suffix.to_string_lossy());
            let candidate = directory.join(file_name);
            if is_executable_file(&candidate) {
                return Some(candidate);
            }
        }
    }
    None
}

fn is_executable_file(path: &Path) -> bool {
    let Ok(metadata) = std::fs::metadata(path) else {
        return false;
    };
    if !metadata.is_file() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        metadata.permissions().mode() & 0o111 != 0
    }
    #[cfg(not(unix))]
    {
        true
    }
}

fn discover_manifest(start: &Path) -> Option<PathBuf> {
    let mut directory = if start.is_file() {
        start.parent()?.to_path_buf()
    } else {
        start.to_path_buf()
    };
    loop {
        let candidate = directory.join(PROJECT_MANIFEST);
        if candidate.is_file() {
            return Some(candidate);
        }
        if !directory.pop() {
            return None;
        }
    }
}

fn workspace_probe(directory: &Path) -> Check {
    if !directory.is_dir() {
        return Check::new(
            "workspace",
            Status::Fail,
            format!("{} is not a directory", directory.display()),
        )
        .remedy("Choose an existing project directory with --project.");
    }
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |duration| duration.as_nanos());
    let probe = directory.join(format!(
        ".aie-doctor-write-probe-{}-{nonce}",
        std::process::id()
    ));
    match OpenOptions::new().write(true).create_new(true).open(&probe) {
        Ok(file) => {
            drop(file);
            match std::fs::remove_file(&probe) {
                Ok(()) => Check::new(
                    "workspace",
                    Status::Pass,
                    format!("{} is readable and writable", directory.display()),
                )
                .data(json!({"path": directory})),
                Err(error) => Check::new(
                    "workspace",
                    Status::Warn,
                    format!(
                        "{} is writable, but the probe could not be removed",
                        directory.display()
                    ),
                )
                .detail(error.to_string())
                .remedy(format!("Remove {} when convenient.", probe.display())),
            }
        }
        Err(error) => Check::new(
            "workspace",
            Status::Fail,
            format!("{} is not writable", directory.display()),
        )
        .detail(error.to_string())
        .remedy("Choose a writable project directory or correct its permissions."),
    }
}

fn check_project(explicit: Option<&Path>, current: &Path) -> (Check, PathBuf) {
    let selected = explicit.unwrap_or(current);
    if explicit.is_some() && !selected.exists() {
        return (
            Check::new(
                "project",
                Status::Fail,
                format!("project path {} does not exist", selected.display()),
            )
            .remedy("Pass a project directory or its aie-project.yaml manifest."),
            selected.to_path_buf(),
        );
    }
    let manifest = if explicit.is_some() && selected.is_file() {
        if selected
            .file_name()
            .is_some_and(|name| name == PROJECT_MANIFEST)
        {
            Some(selected.to_path_buf())
        } else {
            return (
                Check::new(
                    "project",
                    Status::Fail,
                    format!("project file must be named {PROJECT_MANIFEST}"),
                )
                .detail(format!("selected: {}", selected.display()))
                .remedy("Pass the containing project directory or its aie-project.yaml manifest."),
                selected.parent().unwrap_or(current).to_path_buf(),
            );
        }
    } else if selected
        .file_name()
        .is_some_and(|name| name == PROJECT_MANIFEST)
    {
        selected.is_file().then(|| selected.to_path_buf())
    } else {
        discover_manifest(selected)
    };
    if let Some(manifest) = manifest {
        let root = manifest
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .to_path_buf();
        match crate::projectcmd::load_project(&manifest) {
            Ok(project) => {
                let resources: Vec<Value> = project
                    .manifest
                    .resources
                    .iter()
                    .map(|(name, resource)| {
                        match project.resolve_resource(name, &[resource.kind]) {
                            Ok((_, path)) => json!({
                                "name": name,
                                "kind": resource.kind.as_str(),
                                "path": path,
                                "external": resource.external,
                                "status": "ok",
                            }),
                            Err(error) => json!({
                                "name": name,
                                "kind": resource.kind.as_str(),
                                "path": resource.path,
                                "external": resource.external,
                                "status": "invalid",
                                "error": format!("{error:#}"),
                            }),
                        }
                    })
                    .collect();
                let invalid = resources
                    .iter()
                    .filter(|resource| resource["status"] == "invalid")
                    .count();
                let status = if invalid == 0 {
                    Status::Pass
                } else {
                    Status::Fail
                };
                let detail = resources
                    .iter()
                    .find(|resource| resource["status"] == "invalid")
                    .and_then(|resource| {
                        Some(format!(
                            "resource `{}`: {}; manifest: {}",
                            resource["name"].as_str()?,
                            resource["error"].as_str()?,
                            project.manifest_path.display()
                        ))
                    })
                    .unwrap_or_else(|| format!("manifest: {}", project.manifest_path.display()));
                let mut check = Check::new(
                    "project",
                    status,
                    if invalid == 0 {
                        format!("project `{}` is valid", project.manifest.name)
                    } else {
                        format!(
                            "project `{}` has {invalid} unavailable resource(s)",
                            project.manifest.name
                        )
                    },
                )
                .detail(detail)
                .data(json!({
                    "manifest": project.manifest_path,
                    "root": project.root,
                    "name": project.manifest.name,
                    "schema_version": project.manifest.schema_version,
                    "resources": resources,
                }));
                if invalid != 0 {
                    check = check.remedy(
                        "Restore missing inputs or update them with `aie project add --replace`.",
                    );
                }
                (check, root)
            }
            Err(error) => (
                Check::new(
                    "project",
                    Status::Fail,
                    format!("project manifest {} is invalid", manifest.display()),
                )
                .detail(format!("{error:#}"))
                .remedy("Correct the manifest, then run `aie plan check --explain` again."),
                root,
            ),
        }
    } else {
        let root = if selected.is_dir() {
            selected.to_path_buf()
        } else {
            selected.parent().unwrap_or(current).to_path_buf()
        };
        (
            Check::new(
                "project",
                Status::Warn,
                format!(
                    "no {PROJECT_MANIFEST} was found from {}",
                    selected.display()
                ),
            )
            .remedy(
                "Create a project with `aie project init`, or pass standalone inputs to commands.",
            ),
            root,
        )
    }
}

fn check_archive(path: &Path, verify_content: bool) -> Check {
    let result = inspect_archive(path, verify_content);
    match result {
        Ok(data) => {
            let version = data["format_version"].as_u64().unwrap_or(0) as u32;
            let status = if version == SEEKABLE_VERSION {
                Status::Warn
            } else {
                Status::Pass
            };
            let mut check = Check::new(
                format!("archive:{}", path.display()),
                status,
                if verify_content {
                    format!(
                        "{} is a valid archive and all payloads verified",
                        path.display()
                    )
                } else {
                    format!("{} is a valid archive", path.display())
                },
            )
            .data(data);
            if version == SEEKABLE_VERSION {
                check = check.remedy(format!(
                    "Seal this legacy v1 archive with `aie seal-archive {} --out <rooted.aie>`.",
                    path.display()
                ));
            } else if version != VERSION {
                check.status = Status::Warn;
                check = check.detail(format!("reader expected current archive version {VERSION}"));
            }
            check
        }
        Err(error) => Check::new(
            format!("archive:{}", path.display()),
            Status::Fail,
            format!("{} could not be validated", path.display()),
        )
        .detail(format!("{error:#}"))
        .remedy("Confirm that the path names a complete Gravlax .aie archive."),
    }
}

fn meta_u64(meta: &Value, key: &str) -> Result<u64> {
    meta.get(key)
        .and_then(Value::as_u64)
        .with_context(|| format!("archive meta.{key} is missing or is not an unsigned integer"))
}

fn require_sections(names: &BTreeSet<String>, required: &[&str]) -> Result<()> {
    for name in required {
        if !names.contains(*name) {
            bail!("archive is missing required section {name}");
        }
    }
    Ok(())
}

fn indexed_section_number(name: &str, prefix: &str) -> Option<usize> {
    let suffix = name.strip_prefix(prefix)?;
    if suffix.is_empty() || !suffix.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    suffix.parse().ok()
}

fn validate_fast_archive_structure(
    archive: &mut crate::archivecmd::LazyArchive,
    meta: &Value,
    names: &BTreeSet<String>,
) -> Result<Vec<crate::archivecmd::ChunkInfo>> {
    require_sections(
        names,
        &[
            "meta",
            "chroms",
            "cells",
            "shapes",
            "patterns",
            "rans.tables",
            "index.chunks",
            "index.junctions",
            "index.jpost",
            "edges",
        ],
    )?;
    let molecules = meta_u64(meta, "mols")?;
    let classes = meta_u64(meta, "classes")?;
    let cells = meta_u64(meta, "cells")?;
    meta_u64(meta, "shapes")?;
    meta_u64(meta, "patterns")?;
    let chunk_bp = meta_u64(meta, "chunk_bp")?;
    if chunk_bp == 0 || chunk_bp > u32::MAX as u64 {
        bail!("archive meta.chunk_bp must be between 1 and {}", u32::MAX);
    }
    if classes > u32::MAX as u64 || cells > u32::MAX as u64 {
        bail!("archive class or cell count exceeds u32");
    }
    if classes > molecules {
        bail!("archive has {classes} classes but only {molecules} molecules");
    }
    let mut unique_chroms = HashSet::new();
    for (index, name) in archive.chrom_names.iter().enumerate() {
        if name.is_empty() {
            bail!("archive chromosome name {index} is empty");
        }
        if !unique_chroms.insert(name) {
            bail!("archive chromosome dictionary contains duplicate name {name}");
        }
    }

    let expected_cell_bytes = cells
        .checked_mul(4)
        .context("cell dictionary length overflow")?;
    let actual_cell_bytes = archive
        .reader()
        .entries()
        .iter()
        .find(|(name, _, _, _)| name == "cells")
        .map(|(_, _, raw_len, _)| *raw_len)
        .context("archive is missing required section cells")?;
    if actual_cell_bytes != expected_cell_bytes {
        bail!(
            "cells section has {actual_cell_bytes} raw bytes for meta.cells={cells}; expected {expected_cell_bytes}"
        );
    }

    let chunks = crate::archivecmd::read_chunk_index(archive.reader())?;
    let indexed_molecules = chunks.iter().try_fold(0u64, |sum, chunk| {
        sum.checked_add(chunk.n_mols as u64)
            .context("indexed molecule count overflow")
    })?;
    if indexed_molecules != molecules {
        bail!("molecule count mismatch: {indexed_molecules} indexed vs {molecules} in meta.mols");
    }
    if chunks.is_empty() && classes != 0 {
        bail!("empty chunk index cannot describe {classes} classes");
    }
    if chunks.first().is_some_and(|chunk| chunk.class_base != 0) {
        bail!("first chunk class base is not zero");
    }
    let chunk_bp = chunk_bp as u32;
    for (index, chunk) in chunks.iter().enumerate() {
        if chunk.n_mols == 0 {
            bail!("chunk {index} declares zero molecules");
        }
        if chunk.chrom as usize >= archive.chrom_names.len() {
            bail!(
                "chunk {index} references chromosome {} but the dictionary has {} entries",
                chunk.chrom,
                archive.chrom_names.len()
            );
        }
        if chunk.class_base as u64 > classes {
            bail!("chunk {index} class base exceeds meta.classes");
        }
        if chunk.n_cells > chunk.n_mols || chunk.n_cells as u64 > cells {
            bail!("chunk {index} has an impossible distinct-cell count");
        }
        let bin_end = chunk
            .bin_start
            .checked_add(chunk_bp)
            .context("chunk genomic bin overflows u32")?;
        if chunk.max_anchor < chunk.bin_start || chunk.max_anchor >= bin_end {
            bail!("chunk {index} max anchor lies outside its genomic bin");
        }
        if let Some(previous) = index.checked_sub(1).map(|previous| &chunks[previous]) {
            if (chunk.chrom, chunk.bin_start) <= (previous.chrom, previous.bin_start) {
                bail!("chunk index is not strictly sorted by chromosome and bin");
            }
            let maximum_next_base = previous
                .class_base
                .checked_add(previous.n_mols)
                .context("chunk class range overflow")?;
            if chunk.class_base < previous.class_base || chunk.class_base > maximum_next_base {
                bail!("chunk {index} has an impossible class base");
            }
        }
    }
    if let Some(last) = chunks.last() {
        let maximum_classes = last
            .class_base
            .checked_add(last.n_mols)
            .context("final chunk class range overflow")?;
        if (classes as u32) < last.class_base || classes as u32 > maximum_classes {
            bail!("meta.classes is outside the final chunk's possible class range");
        }
    }

    for index in 0..chunks.len() {
        if !names.contains(&format!("c{index}")) {
            bail!("chunk index references missing section c{index}");
        }
    }
    for name in names {
        if let Some(index) = indexed_section_number(name, "c") {
            if index >= chunks.len() || name != &format!("c{index}") {
                bail!("archive has unindexed chunk section {name}");
            }
        }
    }

    let coc_blocks = classes.div_ceil(crate::archivecmd::COC_BLOCK as u64);
    for block in 0..coc_blocks {
        if !names.contains(&format!("coc.{block}")) {
            bail!("archive is missing cell-of-class section coc.{block}");
        }
    }
    for name in names {
        if let Some(block) = indexed_section_number(name, "coc.") {
            if block as u64 >= coc_blocks || name != &format!("coc.{block}") {
                bail!("archive has unexpected cell-of-class section {name}");
            }
        }
    }
    Ok(chunks)
}

fn validate_shape_at(
    position: u32,
    shape: u32,
    shapes: &[evidence_io::archive::Shape],
    context: &str,
) -> Result<()> {
    let shape = shapes
        .get(shape as usize)
        .with_context(|| format!("{context} references missing shape {shape}"))?;
    if shape.blocks.is_empty() {
        bail!("{context} references an empty shape");
    }
    let mut previous_end = None;
    for &(offset, length) in &shape.blocks {
        if length == 0 {
            bail!("{context} references a zero-length shape block");
        }
        let end_offset = offset
            .checked_add(length)
            .with_context(|| format!("{context} shape block overflows u32"))?;
        position
            .checked_add(end_offset)
            .with_context(|| format!("{context} placement coordinate overflows u32"))?;
        if previous_end.is_some_and(|end| offset < end) {
            bail!("{context} references an overlapping or unsorted shape");
        }
        previous_end = Some(end_offset);
    }
    Ok(())
}

fn validate_archive_semantics(
    path: &Path,
    meta: &Value,
    chunks: &[crate::archivecmd::ChunkInfo],
) -> Result<(u64, usize)> {
    use crate::rows::SAME_SHAPE;

    let reader = SectionReader::open(path)?;
    let verified_bytes = reader.verify_all_payloads()?;
    let mut reader = SectionReader::open(path)?;
    let dictionaries = crate::archivecmd::read_dicts(&mut reader)?;
    let expected_cells =
        usize::try_from(meta_u64(meta, "cells")?).context("meta.cells exceeds usize")?;
    let expected_shapes =
        usize::try_from(meta_u64(meta, "shapes")?).context("meta.shapes exceeds usize")?;
    let expected_patterns =
        usize::try_from(meta_u64(meta, "patterns")?).context("meta.patterns exceeds usize")?;
    let expected_classes =
        u32::try_from(meta_u64(meta, "classes")?).context("meta.classes exceeds u32")?;
    let expected_molecules =
        usize::try_from(meta_u64(meta, "mols")?).context("meta.mols exceeds usize")?;
    if dictionaries.cells.len() != expected_cells
        || dictionaries.shapes.len() != expected_shapes
        || dictionaries.patterns.len() != expected_patterns
        || dictionaries.n_classes != expected_classes
        || dictionaries.n_mols != expected_molecules
    {
        bail!("decoded dictionary counts disagree with archive metadata");
    }
    for (class, &cell) in dictionaries.cell_of_class.iter().enumerate() {
        if cell as usize >= dictionaries.cells.len() {
            bail!("class {class} references missing cell {cell}");
        }
    }
    for (pattern_index, pattern) in dictionaries.patterns.iter().enumerate() {
        for (alternative_index, alternative) in pattern.iter().enumerate() {
            if alternative.chrom as usize >= dictionaries.chrom_names.len() {
                bail!(
                    "pattern {pattern_index} alternative {alternative_index} references missing chromosome {}",
                    alternative.chrom
                );
            }
            if alternative.shape != SAME_SHAPE
                && alternative.shape as usize >= dictionaries.shapes.len()
            {
                bail!(
                    "pattern {pattern_index} alternative {alternative_index} references missing shape {}",
                    alternative.shape
                );
            }
        }
    }

    let mut decoded_molecules = 0usize;
    let mut next_class = 0u32;
    for (chunk_index, info) in chunks.iter().enumerate() {
        if info.class_base != next_class {
            bail!(
                "chunk {chunk_index} class base {} does not follow decoded class count {next_class}",
                info.class_base
            );
        }
        let raw = reader.read(&format!("c{chunk_index}"))?;
        let molecules = crate::archivecmd::decode_chunk(
            &raw,
            info,
            Some(&dictionaries.cell_of_class),
            &dictionaries.rans_tables,
        )?;
        if molecules.len() != info.n_mols as usize {
            bail!("chunk {chunk_index} decoded molecule count disagrees with its index");
        }
        let mut cells = HashSet::new();
        let mut previous_anchor = None;
        for (molecule_index, molecule) in molecules.iter().enumerate() {
            let context = format!("chunk {chunk_index} molecule {molecule_index}");
            if molecule.cell as usize >= dictionaries.cells.len() {
                bail!("{context} references missing cell {}", molecule.cell);
            }
            if molecule.umi_class >= dictionaries.n_classes {
                bail!("{context} references missing class {}", molecule.umi_class);
            }
            if molecule.umi_class == next_class {
                next_class = next_class.checked_add(1).context("class count overflow")?;
            } else if molecule.umi_class > next_class {
                bail!(
                    "{context} references class {} before its introduction",
                    molecule.umi_class
                );
            }
            cells.insert(molecule.cell);
            if molecule.chains.is_empty() && molecule.mms.is_empty() {
                bail!("{context} has no evidence children");
            }
            for chain in &molecule.chains {
                if chain.weight == 0 || chain.reps.is_empty() {
                    bail!("{context} has an empty or zero-weight chain");
                }
                for &(position, shape) in &chain.reps {
                    validate_shape_at(position, shape, &dictionaries.shapes, &context)?;
                }
            }
            for &(position, shape, pattern, weight) in &molecule.mms {
                if weight == 0 {
                    bail!("{context} has a zero-weight multimapper");
                }
                validate_shape_at(position, shape, &dictionaries.shapes, &context)?;
                let alternatives = dictionaries
                    .patterns
                    .get(pattern as usize)
                    .with_context(|| format!("{context} references missing pattern {pattern}"))?;
                for alternative in alternatives {
                    let alternative_position = i64::from(position)
                        .checked_add(alternative.offset)
                        .context("alternative placement coordinate overflow")?;
                    let alternative_position = u32::try_from(alternative_position)
                        .with_context(|| format!("{context} alternative lies outside u32"))?;
                    let alternative_shape = if alternative.shape == SAME_SHAPE {
                        shape
                    } else {
                        alternative.shape
                    };
                    validate_shape_at(
                        alternative_position,
                        alternative_shape,
                        &dictionaries.shapes,
                        &context,
                    )?;
                }
            }
            let anchor = molecule.anchor();
            if anchor < info.bin_start || anchor > info.max_anchor {
                bail!("{context} anchor lies outside its indexed range");
            }
            if previous_anchor.is_some_and(|previous| anchor < previous) {
                bail!("chunk {chunk_index} molecule anchors are not sorted");
            }
            previous_anchor = Some(anchor);
        }
        if previous_anchor != Some(info.max_anchor) {
            bail!("chunk {chunk_index} max anchor disagrees with decoded molecules");
        }
        if cells.len() != info.n_cells as usize {
            bail!("chunk {chunk_index} distinct-cell count disagrees with its index");
        }
        decoded_molecules = decoded_molecules
            .checked_add(molecules.len())
            .context("decoded molecule count overflow")?;
    }
    if next_class != dictionaries.n_classes {
        bail!(
            "decoded chunks introduce {next_class} classes, expected {}",
            dictionaries.n_classes
        );
    }
    if decoded_molecules != dictionaries.n_mols {
        bail!(
            "molecule count mismatch: {decoded_molecules} decoded vs {} in meta.mols",
            dictionaries.n_mols
        );
    }
    Ok((verified_bytes, decoded_molecules))
}

fn inspect_archive(path: &Path, verify_content: bool) -> Result<Value> {
    // Use the same open path as ordinary queries. A doctor pass therefore cannot bless an archive
    // whose layout, chromosome dictionary, or rANS tables the query layer itself rejects.
    let archive_path = path.to_path_buf();
    let mut archive = crate::archivecmd::LazyArchive::open(&archive_path)
        .context("opening archive through the production query reader")?;
    let version = archive.reader().archive_version();
    let sections = archive.reader().entries().len();
    let archive_root = archive
        .reader()
        .content_commitment()
        .map(|root| root.to_hex());
    let names: BTreeSet<String> = archive.reader().names().map(str::to_owned).collect();
    let meta: Value = serde_json::from_slice(&archive.reader().read("meta")?)
        .context("archive metadata is not valid JSON")?;
    let chunks = validate_fast_archive_structure(&mut archive, &meta, &names)?;
    let molecules = meta_u64(&meta, "mols")?;
    let (verified_bytes, decoded_molecules) = if verify_content {
        let (bytes, decoded) = validate_archive_semantics(path, &meta, &chunks)?;
        (bytes, Some(decoded))
    } else {
        (0, None)
    };
    Ok(json!({
        "path": path,
        "format_version": version,
        "sections": sections,
        "chromosomes": archive.chrom_names.len(),
        "molecules": molecules,
        "chunks": chunks.len(),
        "archive_root": archive_root,
        "all_payloads_verified": verify_content,
        "semantic_content_verified": verify_content,
        "verified_compressed_bytes": verified_bytes,
        "decoded_molecules": decoded_molecules,
    }))
}

fn check_annotation(path: &Path) -> Check {
    match anno::Annotation::from_compiled(path) {
        Ok(annotation) => {
            let exons: usize = annotation
                .transcripts
                .iter()
                .map(|transcript| transcript.exons.len())
                .sum();
            Check::new(
                format!("annotation:{}", path.display()),
                Status::Pass,
                format!("{} is a valid compiled annotation", path.display()),
            )
            .data(json!({
                "path": path,
                "genes": annotation.gene_ids.len(),
                "transcripts": annotation.transcripts.len(),
                "exons": exons,
            }))
        }
        Err(error) => Check::new(
            format!("annotation:{}", path.display()),
            Status::Fail,
            format!("{} could not be validated", path.display()),
        )
        .detail(format!("{error:#}"))
        .remedy("Rebuild it with `aie compile-annotation <source.gtf> --out <annotation.aic>`."),
    }
}

fn check_path(path: &Path, verify_content: bool) -> Check {
    if !path.exists() {
        return Check::new(
            format!("input:{}", path.display()),
            Status::Fail,
            format!("{} does not exist", path.display()),
        )
        .remedy("Correct the path and run the diagnosis again.");
    }
    if path.is_dir() {
        return Check::new(
            format!("input:{}", path.display()),
            Status::Pass,
            format!("{} is an accessible directory", path.display()),
        );
    }
    match path.extension().and_then(OsStr::to_str) {
        Some("aie") => check_archive(path, verify_content),
        Some("aic") => check_annotation(path),
        _ => match std::fs::metadata(path) {
            Ok(metadata) => Check::new(
                format!("input:{}", path.display()),
                Status::Pass,
                format!("{} is readable", path.display()),
            )
            .detail("No format-specific validation is available for this file type.")
            .data(json!({"path": path, "bytes": metadata.len()})),
            Err(error) => Check::new(
                format!("input:{}", path.display()),
                Status::Fail,
                format!("{} is not readable", path.display()),
            )
            .detail(error.to_string()),
        },
    }
}

fn collect(args: &Args) -> Report {
    let mut checks = Vec::new();
    match std::env::current_exe() {
        Ok(executable) => checks.push(
            Check::new(
                "installation",
                Status::Pass,
                format!("aie {} is runnable", env!("CARGO_PKG_VERSION")),
            )
            .data(json!({"executable": executable, "version": env!("CARGO_PKG_VERSION")})),
        ),
        Err(error) => checks.push(
            Check::new(
                "installation",
                Status::Fail,
                "the aie executable path is unavailable",
            )
            .detail(error.to_string())
            .remedy("Reinstall Gravlax and ensure the executable is on PATH."),
        ),
    }

    match std::thread::available_parallelism() {
        Ok(threads) => {
            let configured = std::env::var("RAYON_NUM_THREADS").ok();
            let invalid = configured
                .as_deref()
                .is_some_and(|value| value.parse::<usize>().ok().is_none_or(|value| value == 0));
            let status = if invalid { Status::Warn } else { Status::Pass };
            let mut check = Check::new(
                "compute",
                status,
                format!("{} hardware thread(s) available", threads.get()),
            )
            .data(json!({
                "hardware_threads": threads.get(),
                "rayon_num_threads": configured,
                "default_thread_cap": 24,
            }));
            if invalid {
                check = check.remedy("Unset RAYON_NUM_THREADS or set it to a positive integer.");
            }
            checks.push(check);
        }
        Err(error) => checks.push(
            Check::new(
                "compute",
                Status::Warn,
                "hardware thread count is unavailable",
            )
            .detail(error.to_string()),
        ),
    }

    let current = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    let (project_check, workspace) = check_project(args.project.as_deref(), &current);
    checks.push(project_check);
    checks.push(workspace_probe(&workspace));

    for (id, executable, purpose) in [
        ("tool:STAR", "STAR", "one-time annotation-free alignment"),
        ("tool:samtools", "samtools", "optional BAM inspection"),
    ] {
        if let Some(path) = executable_in_path(executable) {
            checks.push(
                Check::new(id, Status::Pass, format!("{executable} is available"))
                    .detail(purpose)
                    .data(json!({"executable": path})),
            );
        } else {
            checks.push(
                Check::new(
                    id,
                    Status::Warn,
                    format!("{executable} was not found on PATH"),
                )
                .detail(format!(
                    "It is used for {purpose}, but is not needed for archive queries."
                ))
                .remedy(format!(
                    "Install {executable} only if this machine performs that step."
                )),
            );
        }
    }

    checks.extend(
        args.paths
            .iter()
            .map(|path| check_path(path, args.verify_content)),
    );
    Report { checks }
}

fn print_human(report: &Report) {
    println!("Gravlax doctor {}", env!("CARGO_PKG_VERSION"));
    for check in &report.checks {
        println!("  [{:7}] {}", check.status.marker(), check.summary);
        if let Some(detail) = &check.detail {
            println!("            {detail}");
        }
        if let Some(remedy) = &check.remedy {
            println!("            Next: {remedy}");
        }
    }
    println!(
        "\n{} passed; {} warning(s); {} failure(s)",
        report.checks.len() - report.failures() - report.warnings(),
        report.warnings(),
        report.failures()
    );
}

pub fn run(args: Args) -> Result<()> {
    let strict = args.strict;
    let report = collect(&args);
    if args.json {
        println!("{}", serde_json::to_string_pretty(&report.value(strict))?);
    } else {
        print_human(&report);
    }
    if report.failures() > 0 || (strict && report.warnings() > 0) {
        bail!(
            "doctor found {} failure(s) and {} warning(s)",
            report.failures(),
            report.warnings()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rows::{Extracted, MolChain, MolRec, PatAlt, SAME_SHAPE};
    use evidence_io::archive::{put_varint, Shape};
    use smallvec::smallvec;

    struct Scratch(PathBuf);

    impl Scratch {
        fn new() -> Self {
            static NEXT_SCRATCH: std::sync::atomic::AtomicU64 =
                std::sync::atomic::AtomicU64::new(0);
            let nonce = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let sequence = NEXT_SCRATCH.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "gravlax-doctor-test-{}-{nonce}-{sequence}",
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

    fn args(project: PathBuf) -> Args {
        Args {
            paths: Vec::new(),
            project: Some(project),
            verify_content: false,
            strict: false,
            json: true,
        }
    }

    fn valid_rans_tables() -> Vec<u8> {
        let mut counts = [0u64; evidence_io::rans::NSYM];
        counts[0] = 1;
        let mut encoded = Vec::new();
        for _ in 0..6 {
            evidence_io::rans::Table::from_counts(&counts)
                .unwrap()
                .serialize(&mut encoded);
        }
        encoded
    }

    fn write_minimal_archive(path: &Path, rans_tables: &[u8], chunk_index: &[u8]) {
        let mut writer = evidence_io::format::SectionWriter::create(path, 1).unwrap();
        writer
            .section(
                "meta",
                br#"{"mols":0,"edges":0,"cells":0,"shapes":0,"patterns":0,"classes":0,"chunk_bp":4000000,"chunk_streams":10,"coc_block":65536,"codec":"rans2"}"#,
            )
            .unwrap();
        for (name, bytes) in [
            ("chroms", b"chr1\n".as_slice()),
            ("cells", b"".as_slice()),
            ("shapes", b"".as_slice()),
            ("patterns", b"".as_slice()),
            ("rans.tables", rans_tables),
            ("index.chunks", chunk_index),
            ("index.junctions", b"".as_slice()),
            ("index.jpost", b"".as_slice()),
            ("edges", b"".as_slice()),
        ] {
            writer.section(name, bytes).unwrap();
        }
        writer.finish().unwrap();
    }

    #[test]
    fn discovers_a_project_from_a_nested_directory() {
        let scratch = Scratch::new();
        std::fs::write(
            scratch.0.join(PROJECT_MANIFEST),
            "schema_version: 1\nname: test-project\nresources: {}\n",
        )
        .unwrap();
        let nested = scratch.0.join("one/two");
        std::fs::create_dir_all(&nested).unwrap();
        let (check, root) = check_project(None, &nested);
        assert_eq!(check.status, Status::Pass);
        assert_eq!(root, scratch.0);
    }

    #[test]
    fn missing_explicit_project_is_a_failure() {
        let scratch = Scratch::new();
        let (check, _) = check_project(Some(&scratch.0.join("missing")), &scratch.0);
        assert_eq!(check.status, Status::Fail);
    }

    #[test]
    fn corrupt_archive_is_reported_without_panicking() {
        let scratch = Scratch::new();
        let path = scratch.0.join("broken.aie");
        std::fs::write(&path, b"not an archive").unwrap();
        let check = check_archive(&path, true);
        assert_eq!(check.status, Status::Fail);
        assert!(check.detail.unwrap().contains("truncated") || check.summary.contains("validated"));
    }

    #[test]
    fn report_json_has_a_stable_shape() {
        let scratch = Scratch::new();
        let report = collect(&args(scratch.0.clone()));
        let value = report.value(false);
        assert_eq!(value["schema"], "gravlax.doctor.v1");
        assert!(value["checks"].is_array());
        assert!(value["summary"]["passed"].is_number());
    }

    #[test]
    fn valid_rooted_archive_can_be_fully_verified() {
        let scratch = Scratch::new();
        let path = scratch.0.join("minimal.aie");
        write_minimal_archive(&path, &valid_rans_tables(), b"");
        let check = check_archive(&path, true);
        assert_eq!(check.status, Status::Pass, "{:?}", check.detail);
        assert_eq!(check.data["format_version"], VERSION);
        assert_eq!(check.data["all_payloads_verified"], true);
        assert_eq!(check.data["semantic_content_verified"], true);
        assert_eq!(check.data["molecules"], 0);
        assert_eq!(check.data["decoded_molecules"], 0);
    }

    #[test]
    fn writer_archive_can_be_fully_verified() {
        let scratch = Scratch::new();
        let path = scratch.0.join("writer.aie");
        let extracted = Extracted {
            mols: vec![MolRec {
                cell: 0,
                umi_class: 0,
                chrom: 0,
                strand_rev: false,
                chains: smallvec![MolChain {
                    weight: 1,
                    reps: smallvec![(100, 0), (110, 0)],
                }],
                mms: smallvec![(120, 0, 0, 1)],
            }],
            edges: Vec::new(),
            cells: vec![0],
            shapes: vec![Shape {
                blocks: vec![(0, 10)],
            }],
            patterns: vec![vec![PatAlt {
                chrom: 0,
                offset: 0,
                strand_flip: false,
                shape: SAME_SHAPE,
            }]],
            n_classes: 1,
            chrom_names: vec!["chr1".into()],
        };
        crate::archivecmd::write_archive(&extracted, &path, 1, 4_000_000, None).unwrap();

        let check = check_archive(&path, true);

        assert_eq!(check.status, Status::Pass, "{:?}", check.detail);
        assert_eq!(check.data["molecules"], 1);
        assert_eq!(check.data["decoded_molecules"], 1);
    }

    #[test]
    fn archive_with_invalid_rans_tables_fails_default_diagnosis() {
        let scratch = Scratch::new();
        let path = scratch.0.join("bad-rans.aie");
        write_minimal_archive(&path, b"", b"");
        let check = check_archive(&path, false);
        assert_eq!(check.status, Status::Fail);
        assert!(check
            .detail
            .as_deref()
            .is_some_and(|detail| detail.contains("rans") || detail.contains("varint")));
    }

    #[test]
    fn archive_with_invalid_chunk_index_fails_default_diagnosis() {
        let scratch = Scratch::new();
        let path = scratch.0.join("bad-index.aie");
        let mut invalid_index = Vec::new();
        put_varint(&mut invalid_index, 0);
        write_minimal_archive(&path, &valid_rans_tables(), &invalid_index);
        let check = check_archive(&path, false);
        assert_eq!(check.status, Status::Fail);
        assert!(check
            .detail
            .as_deref()
            .is_some_and(|detail| detail.contains("varint past end")));
    }

    #[test]
    fn workspace_probe_cleans_up_after_itself() {
        let scratch = Scratch::new();
        let check = workspace_probe(&scratch.0);
        assert_eq!(check.status, Status::Pass);
        let names: Vec<_> = std::fs::read_dir(&scratch.0).unwrap().collect();
        assert!(names.is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn path_candidate_must_have_an_executable_permission_bit() {
        use std::os::unix::fs::PermissionsExt;

        let scratch = Scratch::new();
        let candidate = scratch.0.join("STAR");
        std::fs::write(&candidate, b"#!/bin/sh\n").unwrap();
        std::fs::set_permissions(&candidate, std::fs::Permissions::from_mode(0o600)).unwrap();
        assert!(!is_executable_file(&candidate));

        std::fs::set_permissions(&candidate, std::fs::Permissions::from_mode(0o700)).unwrap();
        assert!(is_executable_file(&candidate));
    }
}
