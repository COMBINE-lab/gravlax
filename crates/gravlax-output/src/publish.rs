use crate::OutputError;
use serde::{Deserialize, Serialize};
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_STAGING_FILE: AtomicU64 = AtomicU64::new(0);

fn same_file_identity(expected: &fs::Metadata, observed: &fs::Metadata) -> bool {
    if !expected.is_file() || !observed.is_file() || expected.len() != observed.len() {
        return false;
    }
    same_file_object(expected, observed)
}

fn same_file_object(expected: &fs::Metadata, observed: &fs::Metadata) -> bool {
    if !expected.is_file() || !observed.is_file() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        expected.dev() == observed.dev() && expected.ino() == observed.ino()
    }
    #[cfg(not(unix))]
    {
        expected.modified().ok() == observed.modified().ok()
    }
}

#[cfg(target_os = "linux")]
fn link_open_file_no_replace(file: &File, destination: &Path) -> std::io::Result<()> {
    use std::ffi::CString;
    use std::os::fd::AsRawFd;
    use std::os::unix::ffi::OsStrExt;

    let empty = CString::new("").expect("empty C string is valid");
    let destination = CString::new(destination.as_os_str().as_bytes()).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "output path contains a NUL byte",
        )
    })?;
    let direct = unsafe {
        libc::linkat(
            file.as_raw_fd(),
            empty.as_ptr(),
            libc::AT_FDCWD,
            destination.as_ptr(),
            libc::AT_EMPTY_PATH,
        )
    };
    if direct == 0 {
        return Ok(());
    }
    let direct_error = std::io::Error::last_os_error();
    if !matches!(
        direct_error.raw_os_error(),
        Some(code)
            if code == libc::EPERM
                || code == libc::EINVAL
                || code == libc::ENOENT
                || code == libc::EOPNOTSUPP
    ) {
        return Err(direct_error);
    }

    // Some Linux configurations disallow AT_EMPTY_PATH without a capability. `/proc/self/fd`
    // with AT_SYMLINK_FOLLOW still links the held descriptor's inode, never a replaceable staging
    // pathname.
    let descriptor = CString::new(format!("/proc/self/fd/{}", file.as_raw_fd()))
        .expect("numeric descriptor path is a valid C string");
    let result = unsafe {
        libc::linkat(
            libc::AT_FDCWD,
            descriptor.as_ptr(),
            libc::AT_FDCWD,
            destination.as_ptr(),
            libc::AT_SYMLINK_FOLLOW,
        )
    };
    if result == 0 {
        Ok(())
    } else {
        Err(std::io::Error::last_os_error())
    }
}

#[cfg(not(target_os = "linux"))]
fn link_open_file_no_replace(
    _file: &File,
    staging_path: &Path,
    destination: &Path,
) -> std::io::Result<()> {
    fs::hard_link(staging_path, destination)
}

fn remove_staging_if_owned(path: &Path, expected: &fs::Metadata) -> std::io::Result<bool> {
    let observed = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(true),
        Err(error) => return Err(error),
    };
    if !same_file_object(expected, &observed) {
        return Ok(false);
    }
    fs::remove_file(path)?;
    Ok(true)
}

/// Internal plan-execution mapping from physical staging paths to the logical destinations that
/// should appear in result provenance. Plan execution overrides this variable for every child,
/// including an empty mapping, so an ambient value cannot leak into a planned result.
pub const LOGICAL_OUTPUT_MAP_ENV: &str = "GRAVLAX_LOGICAL_OUTPUT_MAP_V1";

/// Return a stable key for a not-yet-created destination by resolving its existing parent and
/// retaining its final file name. This treats spelling variants such as `out` and `./out`, and
/// paths through symlinked parent directories, as the same destination without following a final
/// component that a caller intends to create.
pub fn canonical_destination_key(path: &Path) -> Result<PathBuf, OutputError> {
    let file_name = path.file_name().ok_or_else(|| {
        OutputError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("output path must name a file: {}", path.display()),
        ))
    })?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let parent = fs::canonicalize(parent)?;
    if !parent.is_dir() {
        return Err(OutputError::Io(std::io::Error::new(
            std::io::ErrorKind::NotADirectory,
            format!("output parent is not a directory: {}", parent.display()),
        )));
    }
    Ok(parent.join(file_name))
}

/// Return the logical destination to report for a physical output path. Direct commands normally
/// see no mapping and report the supplied path unchanged. A malformed mapping fails closed rather
/// than allowing a result to claim an ambiguous destination.
pub fn reported_output_path(path: &Path) -> Result<String, OutputError> {
    let physical = path.to_str().ok_or_else(|| {
        OutputError::InvalidSchema("uniform reports require UTF-8 output paths".into())
    })?;
    let Some(raw) = std::env::var_os(LOGICAL_OUTPUT_MAP_ENV) else {
        return Ok(physical.to_owned());
    };
    let raw = raw.to_str().ok_or_else(|| {
        OutputError::InvalidSchema(format!("{LOGICAL_OUTPUT_MAP_ENV} must contain UTF-8 JSON"))
    })?;
    reported_output_path_from_map(physical, raw)
}

fn reported_output_path_from_map(physical: &str, raw: &str) -> Result<String, OutputError> {
    let value: serde_json::Value = serde_json::from_str(raw).map_err(|error| {
        OutputError::InvalidSchema(format!("invalid {LOGICAL_OUTPUT_MAP_ENV}: {error}"))
    })?;
    let object = value.as_object().ok_or_else(|| {
        OutputError::InvalidSchema(format!("{LOGICAL_OUTPUT_MAP_ENV} must be a JSON object"))
    })?;
    for (source, destination) in object {
        if destination.as_str().is_none() {
            return Err(OutputError::InvalidSchema(format!(
                "{LOGICAL_OUTPUT_MAP_ENV} value for {source:?} must be a string"
            )));
        }
    }
    Ok(object
        .get(physical)
        .and_then(serde_json::Value::as_str)
        .unwrap_or(physical)
        .to_owned())
}

/// Durability is deliberately separate from the uniform data contract. `Flush` makes bytes visible
/// to the operating system; `File` additionally synchronizes the staged inode before installation;
/// `FileAndDirectory` also synchronizes the destination directory after the no-clobber link.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Durability {
    Flush,
    File,
    FileAndDirectory,
}

/// Publication may have completed even when best-effort staging cleanup or a post-install
/// directory sync reports a problem. Callers therefore receive an explicit outcome instead of an
/// error that ambiguously suggests the destination was never installed.
#[derive(Clone, Debug, PartialEq, Eq)]
#[must_use = "publication warnings and achieved durability should be inspected"]
pub struct PublicationOutcome {
    pub durability: Durability,
    pub installed_directory_synced: bool,
    pub staging_cleanup_complete: bool,
    pub cleanup_directory_synced: bool,
    pub warnings: Vec<String>,
}

impl PublicationOutcome {
    pub fn requested_durability_achieved(&self) -> bool {
        self.durability != Durability::FileAndDirectory || self.installed_directory_synced
    }
}

/// Atomically install the exact inode referenced by `file` without replacing `destination`.
/// `staging_path` is used only for ownership-guarded cleanup (and as a guarded portability
/// fallback where the operating system cannot link an open descriptor directly). A concurrent
/// replacement of that pathname is never installed on Linux and is never deleted during cleanup.
pub fn install_open_file_no_clobber(
    file: &File,
    staging_path: &Path,
    destination: &Path,
    durability: Durability,
) -> Result<PublicationOutcome, OutputError> {
    let expected = file.metadata()?;
    if !expected.is_file() {
        return Err(OutputError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "held output is not a regular file",
        )));
    }
    if matches!(durability, Durability::File | Durability::FileAndDirectory) {
        if let Err(error) = file.sync_all() {
            let _ = remove_staging_if_owned(staging_path, &expected);
            return Err(OutputError::Io(error));
        }
    }

    #[cfg(not(target_os = "linux"))]
    {
        let staged = fs::symlink_metadata(staging_path)?;
        if !same_file_identity(&expected, &staged) {
            return Err(OutputError::Io(std::io::Error::other(format!(
                "staging path {} no longer names the produced file",
                staging_path.display()
            ))));
        }
    }

    #[cfg(target_os = "linux")]
    let installed = link_open_file_no_replace(file, destination);
    #[cfg(not(target_os = "linux"))]
    let installed = link_open_file_no_replace(file, staging_path, destination);
    if let Err(error) = installed {
        let _ = remove_staging_if_owned(staging_path, &expected);
        if error.kind() == std::io::ErrorKind::AlreadyExists {
            return Err(OutputError::Io(std::io::Error::new(
                std::io::ErrorKind::AlreadyExists,
                format!(
                    "refusing to replace existing output {}",
                    destination.display()
                ),
            )));
        }
        return Err(OutputError::Io(error));
    }

    let observed = match fs::symlink_metadata(destination) {
        Ok(metadata) => metadata,
        Err(error) => {
            let _ = remove_staging_if_owned(staging_path, &expected);
            return Err(OutputError::Io(error));
        }
    };
    if !same_file_identity(&expected, &observed) {
        let _ = remove_staging_if_owned(staging_path, &expected);
        return Err(OutputError::Io(std::io::Error::other(format!(
            "installed output is not the produced inode; preserving {} because a concurrent \
             process may have replaced that path",
            destination.display()
        ))));
    }

    // From here on the complete destination exists. Report post-install limitations in the
    // outcome instead of returning an error whose retry would encounter that destination.
    let parent = destination
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let mut outcome = PublicationOutcome {
        durability,
        installed_directory_synced: durability != Durability::FileAndDirectory,
        staging_cleanup_complete: false,
        cleanup_directory_synced: durability != Durability::FileAndDirectory,
        warnings: Vec::new(),
    };
    if durability == Durability::FileAndDirectory {
        match File::open(parent).and_then(|directory| directory.sync_all()) {
            Ok(()) => outcome.installed_directory_synced = true,
            Err(error) => outcome.warnings.push(format!(
                "output {} was installed, but its directory could not be synchronized: {error}",
                destination.display()
            )),
        }
    }

    match remove_staging_if_owned(staging_path, &expected) {
        Ok(true) => outcome.staging_cleanup_complete = true,
        Ok(false) => outcome.warnings.push(format!(
            "output {} was installed, but staging path {} now names a different file and was preserved",
            destination.display(),
            staging_path.display()
        )),
        Err(error) => outcome.warnings.push(format!(
            "output {} was installed, but staging link {} could not be removed: {error}",
            destination.display(),
            staging_path.display()
        )),
    }
    if durability == Durability::FileAndDirectory && outcome.staging_cleanup_complete {
        match File::open(parent).and_then(|directory| directory.sync_all()) {
            Ok(()) => outcome.cleanup_directory_synced = true,
            Err(error) => outcome.warnings.push(format!(
                "output {} is durable, but staging cleanup could not be directory-synchronized: {error}",
                destination.display()
            )),
        }
    }
    Ok(outcome)
}

/// Write a sibling staging file and atomically install it without replacing an existing path.
/// Competing writers may both complete their private files, but exactly one hard-link operation can
/// create the destination. A producer or pre-install synchronization error leaves no destination.
pub fn publish_file_no_clobber<F>(
    path: &Path,
    durability: Durability,
    produce: F,
) -> Result<PublicationOutcome, OutputError>
where
    F: FnOnce(&mut BufWriter<File>) -> Result<(), OutputError>,
{
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let metadata = fs::metadata(parent)?;
    if !metadata.is_dir() {
        return Err(OutputError::Io(std::io::Error::new(
            std::io::ErrorKind::NotADirectory,
            format!("output parent is not a directory: {}", parent.display()),
        )));
    }
    let file_name = path.file_name().ok_or_else(|| {
        OutputError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "output path must name a file",
        ))
    })?;

    let mut staged = None;
    for _ in 0..4_096 {
        let sequence = NEXT_STAGING_FILE.fetch_add(1, Ordering::Relaxed);
        let mut staging_name = std::ffi::OsString::from(".");
        staging_name.push(file_name);
        staging_name.push(format!(
            ".gravlax-stage.{}.{}",
            std::process::id(),
            sequence
        ));
        let staging_path = parent.join(staging_name);
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&staging_path)
        {
            Ok(file) => {
                staged = Some((staging_path, file));
                break;
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(OutputError::Io(error)),
        }
    }
    let (staging_path, file) = staged.ok_or_else(|| {
        OutputError::Io(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!(
                "could not allocate a private staging file beside {}",
                path.display()
            ),
        ))
    })?;

    let mut writer = BufWriter::new(file);
    if let Err(error) = produce(&mut writer).and_then(|()| {
        writer.flush()?;
        Ok(())
    }) {
        let expected = writer.get_ref().metadata().ok();
        drop(writer);
        if let Some(expected) = expected.as_ref() {
            let _ = remove_staging_if_owned(&staging_path, expected);
        }
        return Err(error);
    }
    let outcome = install_open_file_no_clobber(writer.get_ref(), &staging_path, path, durability);
    drop(writer);
    outcome
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        CellValueRef, DataType, Field, OutputFormat, ResultContext, StreamingTableWriter,
        TableSchema,
    };
    use std::sync::{Arc, Barrier};

    struct TempDir(std::path::PathBuf);

    impl TempDir {
        fn new(label: &str) -> Self {
            let sequence = NEXT_STAGING_FILE.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "gravlax-output-{label}-{}-{sequence}",
                std::process::id()
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }

        fn path(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn publishes_complete_file_and_never_overwrites() {
        let root = TempDir::new("publish");
        let output = root.path().join("result.json");
        let outcome = publish_file_no_clobber(&output, Durability::FileAndDirectory, |writer| {
            writer.write_all(b"complete\n")?;
            Ok(())
        })
        .unwrap();
        assert_eq!(fs::read(&output).unwrap(), b"complete\n");
        assert!(outcome.requested_durability_achieved());
        assert!(outcome.staging_cleanup_complete);

        let error = publish_file_no_clobber(&output, Durability::Flush, |writer| {
            writer.write_all(b"replacement\n")?;
            Ok(())
        })
        .unwrap_err();
        assert!(error.to_string().contains("refusing to replace"));
        assert_eq!(fs::read(&output).unwrap(), b"complete\n");
    }

    #[test]
    fn held_install_refuses_occupied_destination_and_cleans_only_owned_stage() {
        let root = TempDir::new("held-no-clobber");
        let staging = root.path().join("staging");
        let output = root.path().join("result.tsv");
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&staging)
            .unwrap();
        file.write_all(b"produced\n").unwrap();
        file.flush().unwrap();
        fs::write(&output, b"existing\n").unwrap();

        let error = install_open_file_no_clobber(&file, &staging, &output, Durability::Flush)
            .unwrap_err();
        assert!(error.to_string().contains("refusing to replace"));
        assert_eq!(fs::read(&output).unwrap(), b"existing\n");
        assert!(!staging.exists());
    }

    #[test]
    fn logical_output_map_is_exact_and_fails_closed() {
        let raw = r#"{"/stage/result.json":"/project/results/result.json"}"#;
        assert_eq!(
            reported_output_path_from_map("/stage/result.json", raw).unwrap(),
            "/project/results/result.json"
        );
        assert_eq!(
            reported_output_path_from_map("/stage/other.json", raw).unwrap(),
            "/stage/other.json"
        );
        assert!(reported_output_path_from_map("/stage/result.json", "[]").is_err());
        assert!(reported_output_path_from_map(
            "/stage/result.json",
            r#"{"/stage/result.json":3}"#
        )
        .is_err());
    }

    #[test]
    fn destination_key_resolves_spelling_and_parent_aliases() {
        let root = TempDir::new("destination-key");
        let direct = root.path().join("result.json");
        let dotted = root.path().join(".").join("result.json");
        assert_eq!(
            canonical_destination_key(&direct).unwrap(),
            canonical_destination_key(&dotted).unwrap()
        );

        #[cfg(unix)]
        {
            std::os::unix::fs::symlink(root.path(), root.path().join("alias")).unwrap();
            let through_alias = root.path().join("alias").join("result.json");
            assert_eq!(
                canonical_destination_key(&direct).unwrap(),
                canonical_destination_key(&through_alias).unwrap()
            );
        }
    }

    #[test]
    fn producer_failure_leaves_no_destination_or_staging_file() {
        let root = TempDir::new("failure");
        let output = root.path().join("result.tsv");
        let error = publish_file_no_clobber(&output, Durability::Flush, |writer| {
            writer.write_all(b"partial")?;
            Err(OutputError::Sink("intentional producer failure".into()))
        })
        .unwrap_err();
        assert!(error.to_string().contains("intentional producer failure"));
        assert!(!output.exists());
        assert_eq!(fs::read_dir(root.path()).unwrap().count(), 0);
    }

    #[test]
    fn competing_writers_have_exactly_one_winner() {
        let root = TempDir::new("race");
        let output = root.path().join("result.tsv");
        let barrier = Arc::new(Barrier::new(3));
        let mut handles = Vec::new();
        for payload in [b"first\n".as_slice(), b"second\n".as_slice()] {
            let output = output.clone();
            let barrier = Arc::clone(&barrier);
            handles.push(std::thread::spawn(move || {
                publish_file_no_clobber(&output, Durability::Flush, |writer| {
                    writer.write_all(payload)?;
                    barrier.wait();
                    Ok(())
                })
            }));
        }
        barrier.wait();
        let results: Vec<_> = handles
            .into_iter()
            .map(|handle| handle.join().unwrap())
            .collect();
        assert_eq!(results.iter().filter(|result| result.is_ok()).count(), 1);
        assert_eq!(results.iter().filter(|result| result.is_err()).count(), 1);
        let output = fs::read(&output).unwrap();
        assert!(output == b"first\n" || output == b"second\n");
        assert_eq!(fs::read_dir(root.path()).unwrap().count(), 1);
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn staging_path_replacement_cannot_change_installed_bytes_or_delete_replacement() {
        let root = TempDir::new("staging-replacement");
        let output = root.path().join("result.tsv");
        let displaced = root.path().join("held-produced-file");
        let mut replacement = None;
        let outcome = publish_file_no_clobber(&output, Durability::Flush, |writer| {
            writer.write_all(b"produced bytes\n")?;
            let staging = fs::read_dir(root.path())?
                .next()
                .expect("publisher created a staging path")?
                .path();
            fs::rename(&staging, &displaced)?;
            fs::write(&staging, b"concurrent replacement\n")?;
            replacement = Some(staging);
            Ok(())
        })
        .unwrap();

        let replacement = replacement.expect("test recorded replacement path");
        assert_eq!(fs::read(&output).unwrap(), b"produced bytes\n");
        assert_eq!(fs::read(&replacement).unwrap(), b"concurrent replacement\n");
        assert_eq!(fs::read(&displaced).unwrap(), b"produced bytes\n");
        assert!(!outcome.staging_cleanup_complete);
        assert!(outcome
            .warnings
            .iter()
            .any(|warning| warning.contains("now names a different file")));
    }

    fn table_schema() -> TableSchema {
        TableSchema::new(
            "gravlax.test.atomic-stream.v1",
            vec![
                Field::new("label", DataType::String),
                Field::new("count", DataType::UInt64),
            ],
        )
        .unwrap()
    }

    #[test]
    fn all_streaming_text_formats_publish_only_after_complete_render() {
        let root = TempDir::new("formats");
        for format in [OutputFormat::Text, OutputFormat::Tsv, OutputFormat::Json] {
            let output = root.path().join(format!("result-{format:?}"));
            let outcome = publish_file_no_clobber(&output, Durability::Flush, |writer| {
                let schema = table_schema();
                let mut stream = StreamingTableWriter::new(
                    &mut *writer,
                    &schema,
                    format,
                    &ResultContext::default(),
                    None,
                )?;
                stream
                    .write_row_iter([CellValueRef::String("complete"), CellValueRef::UInt64(1)])?;
                stream.finish()?;
                Ok(())
            })
            .unwrap();
            assert!(outcome.requested_durability_achieved());
            let bytes = fs::read(&output).unwrap();
            assert!(!bytes.is_empty());
            if format == OutputFormat::Json {
                serde_json::from_slice::<serde_json::Value>(&bytes).unwrap();
            }
        }
    }

    #[test]
    fn late_row_validation_failure_never_publishes_partial_render() {
        let root = TempDir::new("late-row-failure");
        for format in [OutputFormat::Text, OutputFormat::Tsv, OutputFormat::Json] {
            let output = root.path().join(format!("result-{format:?}"));
            let error = publish_file_no_clobber(&output, Durability::Flush, |writer| {
                let schema = table_schema();
                let mut stream = StreamingTableWriter::new(
                    &mut *writer,
                    &schema,
                    format,
                    &ResultContext::default(),
                    None,
                )?;
                stream.write_row_iter([CellValueRef::String("valid"), CellValueRef::UInt64(1)])?;
                // The error occurs after an envelope and one complete row have reached staging.
                stream
                    .write_row_iter([CellValueRef::String("invalid"), CellValueRef::Int64(-1)])?;
                stream.finish()?;
                Ok(())
            })
            .unwrap_err();
            assert!(error.to_string().contains("column 1 (count)"));
            assert!(!output.exists());
        }
        assert_eq!(fs::read_dir(root.path()).unwrap().count(), 0);
    }
}
