use evidence_io::format::{compress, SectionReader, MAGIC, SEEKABLE_VERSION, VERSION};
use std::io::{Seek, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-archive-root-{}-{nonce}",
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

fn write_v1(path: &Path, sections: &[(&str, &[u8])]) {
    let mut file = std::fs::File::create(path).unwrap();
    file.write_all(MAGIC).unwrap();
    file.write_all(&SEEKABLE_VERSION.to_le_bytes()).unwrap();
    let mut directory = Vec::new();
    for (name, raw) in sections {
        let compressed = compress(raw, 1).unwrap();
        let offset = file.stream_position().unwrap();
        file.write_all(&[name.len() as u8]).unwrap();
        file.write_all(name.as_bytes()).unwrap();
        file.write_all(&(raw.len() as u64).to_le_bytes()).unwrap();
        file.write_all(&(compressed.len() as u64).to_le_bytes())
            .unwrap();
        file.write_all(&compressed).unwrap();
        directory.push((
            name.to_string(),
            offset,
            raw.len() as u64,
            compressed.len() as u64,
        ));
    }
    file.write_all(&[0]).unwrap();
    let directory_offset = file.stream_position().unwrap();
    file.write_all(&(directory.len() as u32).to_le_bytes())
        .unwrap();
    for (name, offset, raw_len, compressed_len) in directory {
        file.write_all(&[name.len() as u8]).unwrap();
        file.write_all(name.as_bytes()).unwrap();
        file.write_all(&offset.to_le_bytes()).unwrap();
        file.write_all(&raw_len.to_le_bytes()).unwrap();
        file.write_all(&compressed_len.to_le_bytes()).unwrap();
    }
    file.write_all(&directory_offset.to_le_bytes()).unwrap();
    file.write_all(b"AIED").unwrap();
}

fn run(arguments: &[&std::ffi::OsStr]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(arguments)
        .output()
        .unwrap()
}

fn json(output: &Output) -> serde_json::Value {
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn seal_is_byte_preserving_deterministic_atomic_and_no_clobber() {
    let scratch = Scratch::new();
    let source = scratch.0.join("source.v1.aie");
    let sealed = scratch.0.join("sealed.v2.aie");
    let repeated = scratch.0.join("repeated.v2.aie");
    write_v1(
        &source,
        &[("alpha", b"one"), ("optional.future", b"two two")],
    );
    let source_reader = SectionReader::open(&source).unwrap();
    let source_identity = source_reader.scan_legacy_identities().unwrap();

    let output = run(&[
        "seal-archive".as_ref(),
        source.as_os_str(),
        "--out".as_ref(),
        sealed.as_os_str(),
        "--json".as_ref(),
    ]);
    let report = json(&output);
    assert_eq!(report["schema"], "gravlax.archive.seal.v1");
    assert_eq!(report["source_format_version"], SEEKABLE_VERSION);
    assert_eq!(report["output_format_version"], VERSION);
    assert_eq!(report["sections"], 2);
    assert_eq!(
        report["source_identity_content_bytes_read"],
        source.metadata().unwrap().len()
    );
    assert_eq!(
        report["source_full_file_blake3"],
        source_identity
            .full_file_blake3
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>()
    );

    let mut before = SectionReader::open(&source).unwrap();
    let mut after = SectionReader::open(&sealed).unwrap();
    assert_eq!(after.archive_version(), VERSION);
    assert!(after.content_commitment().is_some());
    assert_eq!(
        after.encoded_content_identity().unwrap().unwrap(),
        source_identity.encoded_sections_blake3
    );
    let before_entries = before.entries().to_vec();
    let after_entries = after.entries().to_vec();
    assert_eq!(before_entries, after_entries);
    for (name, _, raw_len, compressed_len) in before_entries {
        let (left, left_raw) = before.read_compressed(&name).unwrap();
        let (right, right_raw) = after.read_compressed(&name).unwrap();
        assert_eq!(left_raw as u64, raw_len);
        assert_eq!(right_raw as u64, raw_len);
        assert_eq!(left.len() as u64, compressed_len);
        assert_eq!(left, right);
    }

    let repeat_output = run(&[
        "seal-archive".as_ref(),
        source.as_os_str(),
        "--out".as_ref(),
        repeated.as_os_str(),
        "--json".as_ref(),
    ]);
    json(&repeat_output);
    assert_eq!(
        std::fs::read(&sealed).unwrap(),
        std::fs::read(&repeated).unwrap()
    );

    let sealed_before = std::fs::read(&sealed).unwrap();
    let refused = run(&[
        "seal-archive".as_ref(),
        source.as_os_str(),
        "--out".as_ref(),
        sealed.as_os_str(),
        "--json".as_ref(),
    ]);
    assert!(!refused.status.success());
    assert!(refused.stdout.is_empty());
    assert_eq!(std::fs::read(&sealed).unwrap(), sealed_before);
    assert!(std::fs::read_dir(&scratch.0).unwrap().all(|entry| !entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .contains("seal-tmp")));

    let v2_source_rejected = run(&[
        "seal-archive".as_ref(),
        sealed.as_os_str(),
        "--out".as_ref(),
        scratch.0.join("v2-reseal.aie").as_os_str(),
        "--json".as_ref(),
    ]);
    assert!(!v2_source_rejected.status.success());
    assert!(v2_source_rejected.stdout.is_empty());
}

#[test]
fn inspect_distinguishes_lazy_root_validation_from_full_payload_verification() {
    let scratch = Scratch::new();
    let source = scratch.0.join("source.v1.aie");
    let sealed = scratch.0.join("sealed.v2.aie");
    write_v1(&source, &[("selected", b"one"), ("unselected", b"two two")]);
    let legacy = json(&run(&[
        "inspect-archive".as_ref(),
        source.as_os_str(),
        "--json".as_ref(),
    ]));
    assert_eq!(legacy["format_version"], SEEKABLE_VERSION);
    assert_eq!(legacy["native_identity"]["scheme"], "full-file-blake3-v1");
    assert_eq!(legacy["verification"]["directory_and_root"], false);
    assert_eq!(legacy["verification"]["all_payloads"], false);
    assert_eq!(
        legacy["verification"]["identity_content_bytes_read"],
        source.metadata().unwrap().len()
    );
    let legacy_text = run(&[
        "inspect-archive".as_ref(),
        source.as_os_str(),
        "--verify-content".as_ref(),
    ]);
    assert!(legacy_text.status.success());
    let legacy_text = String::from_utf8(legacy_text.stdout).unwrap();
    let legacy_payload_bytes: u64 = SectionReader::open(&source)
        .unwrap()
        .section_metadata()
        .map(|entry| entry.compressed_len)
        .sum();
    assert!(legacy_text.contains(&format!(
        "{}-byte full-file identity scan",
        source.metadata().unwrap().len()
    )));
    assert!(legacy_text.contains(&format!(
        "decoded all {legacy_payload_bytes} compressed payload bytes"
    )));

    json(&run(&[
        "seal-archive".as_ref(),
        source.as_os_str(),
        "--out".as_ref(),
        sealed.as_os_str(),
        "--json".as_ref(),
    ]));

    let lazy = json(&run(&[
        "inspect-archive".as_ref(),
        sealed.as_os_str(),
        "--json".as_ref(),
    ]));
    assert_eq!(lazy["schema"], "gravlax.archive.identity.v1");
    assert_eq!(lazy["format_version"], VERSION);
    assert_eq!(lazy["verification"]["directory_and_root"], true);
    assert_eq!(lazy["verification"]["all_payloads"], false);
    assert_eq!(lazy["verification"]["identity_content_bytes_read"], 0);
    assert_eq!(
        lazy["verification"]
            .as_object()
            .unwrap()
            .keys()
            .cloned()
            .collect::<std::collections::BTreeSet<_>>(),
        [
            "all_payloads",
            "directory_and_root",
            "identity_content_bytes_read",
            "ordinary_reads_verify_selected_payloads_only",
        ]
        .into_iter()
        .map(str::to_owned)
        .collect()
    );

    let full = json(&run(&[
        "inspect-archive".as_ref(),
        sealed.as_os_str(),
        "--verify-content".as_ref(),
        "--json".as_ref(),
    ]));
    let payload_bytes: u64 = SectionReader::open(&sealed)
        .unwrap()
        .section_metadata()
        .map(|entry| entry.compressed_len)
        .sum();
    assert_eq!(full["verification"]["all_payloads"], true);
    assert_eq!(
        full["verification"]["identity_content_bytes_read"],
        payload_bytes
    );

    let reader = SectionReader::open(&sealed).unwrap();
    let entry = reader
        .section_metadata()
        .find(|entry| entry.name == "unselected")
        .unwrap();
    let payload_offset = entry.offset + 1 + entry.name.len() as u64 + 16;
    let mut bytes = std::fs::read(&sealed).unwrap();
    bytes[payload_offset as usize] ^= 1;
    std::fs::write(&sealed, bytes).unwrap();

    // The directory/root remains valid, so lazy inspection is intentionally successful.
    json(&run(&[
        "inspect-archive".as_ref(),
        sealed.as_os_str(),
        "--json".as_ref(),
    ]));
    let corrupted_full = run(&[
        "inspect-archive".as_ref(),
        sealed.as_os_str(),
        "--verify-content".as_ref(),
        "--json".as_ref(),
    ]);
    assert!(!corrupted_full.status.success());
    assert!(corrupted_full.stdout.is_empty());
}

#[cfg(unix)]
#[test]
fn full_inspection_rejects_concurrent_same_inode_metadata_change_without_output() {
    use std::os::unix::fs::PermissionsExt;
    use std::process::Stdio;

    let scratch = Scratch::new();
    let source = scratch.0.join("changing.v1.aie");
    // A large declared raw section keeps decompression active long enough for a deterministic
    // same-inode ctime mutation between the command's before/after metadata snapshots.
    let raw = vec![0u8; 64 * 1024 * 1024];
    write_v1(&source, &[("large", raw.as_slice())]);
    drop(raw);

    let mut child = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args([
            "inspect-archive".as_ref(),
            source.as_os_str(),
            "--verify-content".as_ref(),
            "--json".as_ref(),
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut mode = 0o600;
    while child.try_wait().unwrap().is_none() {
        std::fs::set_permissions(&source, std::fs::Permissions::from_mode(mode)).unwrap();
        mode = if mode == 0o600 { 0o640 } else { 0o600 };
        std::thread::yield_now();
    }
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr)
        .contains("archive changed while its complete content was being inspected"));
}

#[test]
fn failed_seal_never_installs_a_partial_output() {
    let scratch = Scratch::new();
    let source = scratch.0.join("corrupt.v1.aie");
    let output = scratch.0.join("must-not-exist.v2.aie");
    write_v1(&source, &[("alpha", b"one"), ("beta", b"two")]);
    let reader = SectionReader::open(&source).unwrap();
    let entry = reader.section_metadata().next().unwrap();
    let payload_offset = entry.offset + 1 + entry.name.len() as u64 + 16;
    let mut bytes = std::fs::read(&source).unwrap();
    bytes[payload_offset as usize] ^= 1;
    std::fs::write(&source, bytes).unwrap();

    let failed = run(&[
        "seal-archive".as_ref(),
        source.as_os_str(),
        "--out".as_ref(),
        output.as_os_str(),
        "--json".as_ref(),
    ]);
    assert!(!failed.status.success());
    assert!(failed.stdout.is_empty());
    assert!(!output.exists());
    assert!(std::fs::read_dir(&scratch.0).unwrap().all(|entry| !entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .contains("seal-tmp")));
}
