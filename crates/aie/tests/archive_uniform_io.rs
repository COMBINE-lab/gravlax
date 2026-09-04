use evidence_io::format::{compress, SectionReader, SectionWriter, MAGIC, SEEKABLE_VERSION};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam as sam;
use sam::alignment::{
    io::Write as _,
    record::{
        cigar::{op::Kind, Op},
        data::field::Tag,
        Flags,
    },
    record_buf::{data::field::Value, Cigar, Data, QualityScores, RecordBuf, Sequence},
};
use sam::header::record::value::{map::ReferenceSequence, Map};
use std::io::{Seek, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BARCODE: &str = "AAAAAAAAAAAAAAAA";

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-archive-uniform-{}-{nonce}",
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

fn run(command: &mut Command) -> Output {
    let debug = format!("{command:?}");
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "command failed: {debug}\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn run_failing(command: &mut Command) -> Output {
    let output = command.output().unwrap();
    assert!(
        !output.status.success(),
        "command unexpectedly succeeded: {command:?}"
    );
    output
}

fn hex(digest: [u8; 32]) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn write_v2(path: &Path) {
    let mut writer = SectionWriter::create(path, 1).unwrap();
    writer.section("meta", b"{}").unwrap();
    writer.section("chroms", b"chr1\n").unwrap();
    writer.section("payload", b"scientific payload\n").unwrap();
    writer.finish().unwrap();
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
    for (name, offset, raw, compressed) in directory {
        file.write_all(&[name.len() as u8]).unwrap();
        file.write_all(name.as_bytes()).unwrap();
        file.write_all(&offset.to_le_bytes()).unwrap();
        file.write_all(&raw.to_le_bytes()).unwrap();
        file.write_all(&compressed.to_le_bytes()).unwrap();
    }
    file.write_all(&directory_offset.to_le_bytes()).unwrap();
    file.write_all(b"AIED").unwrap();
}

fn one_read_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(10_000).unwrap()),
        )
        .build();
    let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
    let record = |name: &str, start: usize, umi: &str, flags: Flags, hits: u8| {
        let data: Data = [
            (Tag::new(b'C', b'R'), Value::from(BARCODE)),
            (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
            (Tag::new(b'U', b'R'), Value::from(umi)),
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(hits)),
        ]
        .into_iter()
        .collect();
        RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(cigar.clone())
            .set_sequence(Sequence::from(vec![b'A'; 50]))
            .set_quality_scores(QualityScores::from(vec![30; 50]))
            .set_data(data)
            .build()
    };
    let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for record in [
        record("read-1a", 101, "AAAAAAAAAAAA", Flags::empty(), 1),
        record("read-1b", 121, "AAAAAAAAAAAA", Flags::empty(), 1),
        record("read-edge", 101, "AAAAAAAAAAAC", Flags::empty(), 1),
        record("read-2", 151, "TTTTTTTTTTTT", Flags::empty(), 1),
        record("multi", 201, "CCCCCCCCCCCC", Flags::empty(), 2),
        record("multi", 401, "CCCCCCCCCCCC", Flags::SECONDARY, 2),
        record("read-far", 1_101, "AAAAAAAAAAAA", Flags::empty(), 1),
    ] {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

#[test]
fn inspect_and_seal_reports_are_typed_atomic_and_legacy_exact() {
    let scratch = Scratch::new();
    let bin = env!("CARGO_BIN_EXE_aie");
    let rooted = scratch.0.join("rooted.aie");
    write_v2(&rooted);

    let reader = SectionReader::open(&rooted).unwrap();
    let root = reader.content_commitment().unwrap().to_hex();
    let encoded = hex(reader.encoded_content_identity().unwrap().unwrap());
    let legacy = run(Command::new(bin).arg("inspect-archive").arg(&rooted));
    assert_eq!(
        String::from_utf8(legacy.stdout).unwrap(),
        format!(
            "{}: archive v2 with 3 sections and {} bytes\n\
native identity: aie-directory-root-v2:{root}\n\
encoded sections: aie-encoded-sections-v1:{encoded}\n\
molecular evidence schema: unavailable (legacy archive)\n\
alignment provenance: unavailable; junction discovery is unknown\n\
terminal tails: unavailable (extraction rule was not recorded as evaluated)\n\
genome reference binding: legacy/unattributed\n\
verified directory/root; payloads will be verified when selected\n",
            rooted.display(),
            std::fs::metadata(&rooted).unwrap().len(),
        )
    );

    let report = scratch.0.join("inspect.json");
    let uniform = run(Command::new(bin)
        .arg("inspect-archive")
        .arg(&rooted)
        .args(["--format", "json", "--output"])
        .arg(&report));
    assert!(uniform.stdout.is_empty());
    let report_bytes = std::fs::read(&report).unwrap();
    let value: serde_json::Value = serde_json::from_slice(&report_bytes).unwrap();
    assert_eq!(value["$schema"], "gravlax.result-envelope.v1");
    assert_eq!(value["result_schema"], "gravlax.archive.inspect-report.v1");
    assert_eq!(
        value["data"]["summary"]["archive"]["native_identity"]["blake3"],
        root
    );
    assert_eq!(
        value["data"]["tables"][0]["rows"].as_array().unwrap().len(),
        3
    );

    let failed = run_failing(
        Command::new(bin)
            .arg("inspect-archive")
            .arg(&rooted)
            .args(["--format", "json", "--output"])
            .arg(&report),
    );
    assert!(String::from_utf8_lossy(&failed.stderr).contains("refusing to replace"));
    assert_eq!(std::fs::read(&report).unwrap(), report_bytes);

    let malformed = scratch.0.join("malformed.aie");
    let absent = scratch.0.join("malformed.json");
    std::fs::write(&malformed, b"not an archive").unwrap();
    run_failing(
        Command::new(bin)
            .arg("inspect-archive")
            .arg(&malformed)
            .args(["--format", "json", "--output"])
            .arg(&absent),
    );
    assert!(!absent.exists());
    assert!(std::fs::read_dir(&scratch.0).unwrap().all(|entry| !entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .contains("gravlax-stage")));

    let legacy_archive = scratch.0.join("legacy.aie");
    write_v1(&legacy_archive, &[("one", b"alpha"), ("two", b"beta")]);
    let sealed = scratch.0.join("sealed.aie");
    let seal_report = scratch.0.join("seal.json");
    run(Command::new(bin)
        .arg("seal-archive")
        .arg(&legacy_archive)
        .arg("--out")
        .arg(&sealed)
        .args(["--report-format", "json", "--report-output"])
        .arg(&seal_report));
    let value: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&seal_report).unwrap()).unwrap();
    assert_eq!(value["result_schema"], "gravlax.archive.seal-report.v1");
    let legacy_identity = SectionReader::open(&legacy_archive)
        .unwrap()
        .scan_legacy_identities()
        .unwrap()
        .encoded_sections_blake3;
    assert_eq!(
        SectionReader::open(&sealed)
            .unwrap()
            .encoded_content_identity()
            .unwrap()
            .unwrap(),
        legacy_identity
    );

    let blocked_output = scratch.0.join("must-not-exist.aie");
    let existing_report = scratch.0.join("existing-report.json");
    std::fs::write(&existing_report, b"keep").unwrap();
    run_failing(
        Command::new(bin)
            .arg("seal-archive")
            .arg(&legacy_archive)
            .arg("--out")
            .arg(&blocked_output)
            .args(["--report-format", "json", "--report-output"])
            .arg(&existing_report),
    );
    assert!(
        !blocked_output.exists(),
        "report preflight must precede sealing"
    );
    assert_eq!(std::fs::read(&existing_report).unwrap(), b"keep");
}

#[test]
fn artifact_report_preflight_rejects_normalized_destination_aliases() {
    let scratch = Scratch::new();
    let bin = env!("CARGO_BIN_EXE_aie");
    let artifacts = scratch.0.join("artifacts");
    std::fs::create_dir(&artifacts).unwrap();
    let primary = artifacts.join("result.aie");

    let dotted_alias = artifacts.join(".").join("result.aie");
    let dotted = run_failing(
        Command::new(bin)
            .arg("ingest-archive")
            .arg(scratch.0.join("unused.bam"))
            .arg("--whitelist")
            .arg(scratch.0.join("unused-whitelist.txt"))
            .arg("--out")
            .arg(&primary)
            .args(["--report-format", "json", "--report-output"])
            .arg(&dotted_alias),
    );
    assert!(String::from_utf8_lossy(&dotted.stderr)
        .contains("operation report path must differ from the primary artifact path"));
    assert!(!primary.exists(), "alias rejection must precede ingest");

    #[cfg(unix)]
    {
        use std::os::unix::fs::symlink;

        let alias_parent = scratch.0.join("artifact-alias");
        symlink(&artifacts, &alias_parent).unwrap();
        let through_symlink = alias_parent.join("result.aie");
        let symlinked = run_failing(
            Command::new(bin)
                .arg("ingest-archive")
                .arg(scratch.0.join("unused.bam"))
                .arg("--whitelist")
                .arg(scratch.0.join("unused-whitelist.txt"))
                .arg("--out")
                .arg(&primary)
                .args(["--report-format", "json", "--report-output"])
                .arg(&through_symlink),
        );
        assert!(String::from_utf8_lossy(&symlinked.stderr)
            .contains("operation report path must differ from the primary artifact path"));
        assert!(!primary.exists(), "symlink alias rejection must precede ingest");
    }
}

#[test]
fn uniform_ingest_replay_and_stamp_preserve_primary_artifact_bytes() {
    let scratch = Scratch::new();
    let bin = env!("CARGO_BIN_EXE_aie");
    let bam = scratch.0.join("one.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let barcodes = scratch.0.join("barcodes.tsv");
    let gtf = scratch.0.join("genes.gtf");
    one_read_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    std::fs::write(&barcodes, format!("{BARCODE}\n")).unwrap();
    std::fs::write(
        &gtf,
        "chr1\ttest\texon\t1\t10000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
    )
    .unwrap();

    let legacy_archive = scratch.0.join("legacy.aie");
    let uniform_archive = scratch.0.join("uniform.aie");
    run(Command::new(bin)
        .arg("ingest-archive")
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&legacy_archive)
        .args(["--zstd-level", "1", "--chunk-mb", "1"]));
    let ingest_report = scratch.0.join("ingest.json");
    let ingest = run(Command::new(bin)
        .arg("ingest-archive")
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&uniform_archive)
        .args(["--zstd-level", "1", "--chunk-mb", "1"])
        .args(["--report-format", "json", "--report-output"])
        .arg(&ingest_report));
    assert!(ingest.stdout.is_empty());
    assert_eq!(
        std::fs::read(&legacy_archive).unwrap(),
        std::fs::read(&uniform_archive).unwrap(),
        "uniform ingest reporting changed archive bytes"
    );
    let ingest_value: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&ingest_report).unwrap()).unwrap();
    assert_eq!(
        ingest_value["result_schema"],
        "gravlax.archive.ingest-report.v1"
    );
    assert_eq!(
        ingest_value["data"]["summary"]["inputs"]["bam"]["blake3"],
        blake3::hash(&std::fs::read(&bam).unwrap())
            .to_hex()
            .to_string()
    );
    assert_eq!(
        ingest_value["data"]["summary"]["inputs"]["whitelist"]["blake3"],
        blake3::hash(&std::fs::read(&whitelist).unwrap())
            .to_hex()
            .to_string()
    );
    let uniform_root = SectionReader::open(&uniform_archive)
        .unwrap()
        .content_commitment()
        .unwrap()
        .to_hex();
    assert_eq!(
        ingest_value["data"]["summary"]["output_archive"]["native_identity"]["blake3"],
        uniform_root
    );

    let occupied_archive = scratch.0.join("occupied.aie");
    let blocked_ingest_report = scratch.0.join("blocked-ingest.json");
    std::fs::write(&occupied_archive, b"keep existing archive\n").unwrap();
    let blocked_ingest = run_failing(
        Command::new(bin)
            .arg("ingest-archive")
            .arg(&bam)
            .arg("--whitelist")
            .arg(&whitelist)
            .arg("--out")
            .arg(&occupied_archive)
            .args(["--zstd-level", "1", "--chunk-mb", "1"])
            .args(["--report-format", "json", "--report-output"])
            .arg(&blocked_ingest_report),
    );
    assert!(String::from_utf8_lossy(&blocked_ingest.stderr).contains("refusing to replace"));
    assert_eq!(
        std::fs::read(&occupied_archive).unwrap(),
        b"keep existing archive\n"
    );
    assert!(!blocked_ingest_report.exists());
    assert!(std::fs::read_dir(&scratch.0).unwrap().all(|entry| !entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .contains("ingest-tmp")));

    let legacy_raw_mex = scratch.0.join("legacy-raw-mex");
    run(Command::new(bin)
        .arg("replay-rows")
        .arg(&bam)
        .arg("--from-bam")
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&legacy_raw_mex));
    let raw_mex = scratch.0.join("uniform-raw-mex");
    let raw_report = scratch.0.join("raw-replay.json");
    run(Command::new(bin)
        .arg("replay-rows")
        .arg(&bam)
        .arg("--from-bam")
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&raw_mex)
        .args(["--report-format", "json", "--report-output"])
        .arg(&raw_report));
    let raw_replay: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&raw_report).unwrap()).unwrap();
    assert_eq!(raw_replay["data"]["summary"]["source_kind"], "bam");
    assert_eq!(
        raw_replay["data"]["summary"]["input_identity"]["blake3"],
        blake3::hash(&std::fs::read(&bam).unwrap())
            .to_hex()
            .to_string()
    );
    assert_eq!(
        raw_replay["data"]["summary"]["whitelist_identity"]["blake3"],
        blake3::hash(&std::fs::read(&whitelist).unwrap())
            .to_hex()
            .to_string()
    );
    for name in ["matrix.mtx", "features.tsv", "barcodes.tsv"] {
        assert_eq!(
            std::fs::read(legacy_raw_mex.join(name)).unwrap(),
            std::fs::read(raw_mex.join(name)).unwrap(),
            "uniform raw-BAM replay reporting changed {name}"
        );
    }

    let legacy_mex = scratch.0.join("legacy-mex");
    let uniform_mex = scratch.0.join("uniform-mex");
    run(Command::new(bin)
        .arg("replay-rows")
        .arg(&legacy_archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&legacy_mex));
    let replay_report = scratch.0.join("replay.json");
    run(Command::new(bin)
        .arg("replay-rows")
        .arg(&uniform_archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&uniform_mex)
        .args(["--report-format", "json", "--report-output"])
        .arg(&replay_report));
    for name in ["matrix.mtx", "features.tsv", "barcodes.tsv"] {
        assert_eq!(
            std::fs::read(legacy_mex.join(name)).unwrap(),
            std::fs::read(uniform_mex.join(name)).unwrap(),
            "uniform replay reporting changed {name}"
        );
    }
    let metadata: serde_json::Value =
        serde_json::from_slice(&std::fs::read(uniform_mex.join("metadata.json")).unwrap()).unwrap();
    assert_eq!(metadata["$schema"], "gravlax.result-envelope.v1");
    assert_eq!(metadata["result_schema"], "gravlax.replay.mex-artifact.v1");
    assert_eq!(metadata["data"]["format"], "matrix_market_coordinate");
    assert_eq!(metadata["data"]["matrix"], "matrix.mtx");
    assert_eq!(metadata["data"]["features"], "features.tsv");
    assert_eq!(metadata["data"]["barcodes"], "barcodes.tsv");
    assert_eq!(metadata["data"]["index_base"], 1);
    assert_eq!(metadata["data"]["value_type"], "integer");
    let replay: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&replay_report).unwrap()).unwrap();
    assert_eq!(replay["result_schema"], "gravlax.archive.replay-report.v1");
    assert_eq!(replay["data"]["summary"]["counted_umis"], 3);

    let legacy_velocity = scratch.0.join("legacy-velocity");
    let uniform_velocity = scratch.0.join("uniform-velocity");
    run(Command::new(bin)
        .arg("replay-rows")
        .arg(&legacy_archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&legacy_velocity)
        .arg("--velocity"));
    let velocity_report = scratch.0.join("velocity.json");
    run(Command::new(bin)
        .arg("replay-rows")
        .arg(&uniform_archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&uniform_velocity)
        .arg("--velocity")
        .args(["--report-format", "json", "--report-output"])
        .arg(&velocity_report));
    for name in [
        "spliced.mtx",
        "unspliced.mtx",
        "ambiguous.mtx",
        "entries.complete.tsv",
        "features.tsv",
        "barcodes.tsv",
    ] {
        assert_eq!(
            std::fs::read(legacy_velocity.join(name)).unwrap(),
            std::fs::read(uniform_velocity.join(name)).unwrap(),
            "uniform velocity reporting changed {name}"
        );
    }
    let velocity_metadata: serde_json::Value =
        serde_json::from_slice(&std::fs::read(uniform_velocity.join("metadata.json")).unwrap())
            .unwrap();
    assert_eq!(velocity_metadata["data"]["matrix"], "spliced.mtx");
    assert_eq!(
        velocity_metadata["data"]["matrices"]
            .as_array()
            .unwrap()
            .len(),
        3
    );

    let audit_report = scratch.0.join("audit.json");
    let audit = run(Command::new(bin)
        .arg("replay-rows")
        .arg(&uniform_archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(scratch.0.join("unused-audit-output"))
        .arg("--audit-multigene")
        .args(["--report-format", "json", "--report-output"])
        .arg(&audit_report));
    assert!(audit.stdout.is_empty());
    let audit: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&audit_report).unwrap()).unwrap();
    assert_eq!(audit["data"]["summary"]["multigene_audit"], true);
    assert!(audit["data"]["summary"]["audit"]["gene_evidence_classes"].is_number());

    let blocked_mex = scratch.0.join("blocked-mex");
    let existing_report = scratch.0.join("blocked-report.json");
    std::fs::write(&existing_report, b"keep").unwrap();
    run_failing(
        Command::new(bin)
            .arg("replay-rows")
            .arg(&uniform_archive)
            .arg("--gtf")
            .arg(&gtf)
            .arg("--barcodes")
            .arg(&barcodes)
            .arg("--out-dir")
            .arg(&blocked_mex)
            .args(["--report-format", "json", "--report-output"])
            .arg(&existing_report),
    );
    assert!(
        !blocked_mex.exists(),
        "report preflight must precede replay"
    );

    let colliding_mex = scratch.0.join("colliding-mex");
    std::fs::create_dir(&colliding_mex).unwrap();
    let collision = run_failing(
        Command::new(bin)
            .arg("replay-rows")
            .arg(&uniform_archive)
            .arg("--gtf")
            .arg(&gtf)
            .arg("--barcodes")
            .arg(&barcodes)
            .arg("--out-dir")
            .arg(&colliding_mex)
            .args(["--report-format", "json", "--report-output"])
            .arg(colliding_mex.join("matrix.mtx")),
    );
    assert!(String::from_utf8_lossy(&collision.stderr)
        .contains("must differ from every replay artifact component"));
    assert_eq!(std::fs::read_dir(&colliding_mex).unwrap().count(), 0);

    let source_a = scratch.0.join("stamp-source-a.aie");
    let source_b = scratch.0.join("stamp-source-b.aie");
    write_v2(&source_a);
    std::fs::copy(&source_a, &source_b).unwrap();
    let fasta = scratch.0.join("genome.fa");
    std::fs::write(&fasta, ">chr1\nACGTACGT\n").unwrap();
    let legacy_stamped = scratch.0.join("legacy-stamped.aie");
    let uniform_stamped = scratch.0.join("uniform-stamped.aie");
    run(Command::new(bin)
        .arg("stamp-genome")
        .arg(&source_a)
        .arg("--genome")
        .arg(&fasta)
        .arg("--out")
        .arg(&legacy_stamped));
    let stamp_report = scratch.0.join("stamp.json");
    run(Command::new(bin)
        .arg("stamp-genome")
        .arg(&source_b)
        .arg("--genome")
        .arg(&fasta)
        .arg("--out")
        .arg(&uniform_stamped)
        .args(["--report-format", "json", "--report-output"])
        .arg(&stamp_report));
    assert_eq!(
        std::fs::read(&legacy_stamped).unwrap(),
        std::fs::read(&uniform_stamped).unwrap(),
        "uniform stamp reporting changed archive bytes"
    );
    let stamp: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&stamp_report).unwrap()).unwrap();
    assert_eq!(
        stamp["result_schema"],
        "gravlax.archive.stamp-genome-report.v1"
    );
    assert_eq!(stamp["data"]["summary"]["changed"], true);
    assert_eq!(
        stamp["data"]["summary"]["genome_file"]["blake3"],
        blake3::hash(&std::fs::read(&fasta).unwrap())
            .to_hex()
            .to_string()
    );

    let in_place_stamped = scratch.0.join("in-place-stamped.aie");
    let in_place_report = scratch.0.join("in-place-stamp.json");
    write_v2(&in_place_stamped);
    run(Command::new(bin)
        .arg("stamp-genome")
        .arg(&in_place_stamped)
        .arg("--genome")
        .arg(&fasta)
        .args(["--report-format", "json", "--report-output"])
        .arg(&in_place_report));
    assert_eq!(
        std::fs::read(&in_place_stamped).unwrap(),
        std::fs::read(&uniform_stamped).unwrap(),
        "held-file in-place stamping changed archive bytes"
    );
    let in_place: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&in_place_report).unwrap()).unwrap();
    assert_eq!(in_place["data"]["summary"]["changed"], true);
    assert_eq!(
        in_place["data"]["summary"]["output_archive"]["native_identity"]["blake3"],
        SectionReader::open(&in_place_stamped)
            .unwrap()
            .content_commitment()
            .unwrap()
            .to_hex()
    );
    assert!(std::fs::read_dir(&scratch.0).unwrap().all(|entry| {
        let name = entry.unwrap().file_name();
        let name = name.to_string_lossy();
        !name.contains("stamp-tmp") && !name.contains("stamp-commit")
    }));

    let unchanged_copy = scratch.0.join("already-stamped-copy.aie");
    let unchanged_report = scratch.0.join("already-stamped-copy.json");
    run(Command::new(bin)
        .arg("stamp-genome")
        .arg(&uniform_stamped)
        .arg("--genome")
        .arg(&fasta)
        .arg("--out")
        .arg(&unchanged_copy)
        .args(["--report-format", "json", "--report-output"])
        .arg(&unchanged_report));
    assert_eq!(
        std::fs::read(&uniform_stamped).unwrap(),
        std::fs::read(&unchanged_copy).unwrap(),
        "a no-op stamp with --out must copy the exact source archive"
    );
    let unchanged: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&unchanged_report).unwrap()).unwrap();
    assert_eq!(unchanged["data"]["summary"]["changed"], false);
    assert_eq!(
        unchanged["data"]["summary"]["source_archive"],
        unchanged["data"]["summary"]["output_archive"]
    );
}
