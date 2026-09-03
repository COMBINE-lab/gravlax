use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam as sam;
use sam::alignment::{
    io::Write,
    record::{
        cigar::{op::Kind, Op},
        data::field::Tag,
        Flags,
    },
    record_buf::{data::field::Value, Cigar, Data, QualityScores, RecordBuf, Sequence},
};
use sam::header::record::value::{map::ReferenceSequence, Map};
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
            "gravlax-extend-uniform-{}-{nonce}",
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

fn run(mut command: Command) -> Output {
    let debug = format!("{command:?}");
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "command failed: {debug}\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    output
}

fn run_failure(mut command: Command) -> Output {
    let debug = format!("{command:?}");
    let output = command.output().unwrap();
    assert!(
        !output.status.success(),
        "command unexpectedly succeeded: {debug}\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    output
}

fn write_fixture_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1_000_000).unwrap()),
        )
        .build();
    let record = |name: &str, start: usize, umi: &str, flags: Flags, nh: u8| {
        let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
        let data: Data = [
            (Tag::new(b'C', b'R'), Value::from(BARCODE)),
            (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
            (Tag::new(b'U', b'R'), Value::from(umi)),
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(nh)),
        ]
        .into_iter()
        .collect();
        RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(vec![b'A'; 50]))
            .set_quality_scores(QualityScores::from(vec![30; 50]))
            .set_data(data)
            .build()
    };
    let records = [
        record("downstream", 151, "CCCCCCCCCCCC", Flags::empty(), 1),
        record("edge-a", 100_101, "AAAAAAAAAAAA", Flags::empty(), 1),
        record("edge-b", 100_121, "AAAAAAAAAAAA", Flags::empty(), 1),
        record("edge-c", 100_101, "AAAAAAAAAAAC", Flags::empty(), 1),
        record("mm-a", 500_201, "GGGGGGGGGGGG", Flags::empty(), 3),
        record("mm-a", 500_201, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        record("mm-a", 500_401, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        record("edge-far", 800_101, "AAAAAAAAAAAA", Flags::empty(), 1),
    ];
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in &records {
        writer.write_alignment_record(&header, record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn prepare_fixture(root: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let bam = root.join("fixture.bam");
    let whitelist = root.join("whitelist.txt");
    let archive = root.join("fixture.aie");
    let gtf = root.join("fixture.gtf");
    let genome = root.join("genome.fa");
    write_fixture_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    std::fs::write(
        &gtf,
        concat!(
            "chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\"; gene_name \"Gene1\";\n",
            "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
            "chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
        ),
    )
    .unwrap();
    std::fs::write(&genome, format!(">chr1\n{}\n", "C".repeat(1_000_000))).unwrap();
    let mut ingest = Command::new(env!("CARGO_BIN_EXE_aie"));
    ingest
        .arg("ingest-archive")
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .arg("--genome")
        .arg(&genome)
        .arg("--zstd-level")
        .arg("1")
        .arg("--chunk-mb")
        .arg("1");
    run(ingest);
    (archive, gtf, genome)
}

fn extend_command(archive: &Path, gtf: &Path, out_gtf: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_aie"));
    command
        .arg("extend")
        .arg(archive)
        .arg("--gtf")
        .arg(gtf)
        .arg("--out-gtf")
        .arg(out_gtf)
        .arg("--min-umis")
        .arg("1")
        .arg("--min-cells")
        .arg("1")
        .arg("--min-extension")
        .arg("1");
    command
}

fn table<'a>(result: &'a serde_json::Value, name: &str) -> &'a serde_json::Value {
    result["data"]["tables"]
        .as_array()
        .unwrap()
        .iter()
        .find(|table| table["name"] == name)
        .unwrap()
}

#[test]
fn uniform_extend_preserves_primary_artifacts_and_exposes_typed_bundle() {
    let scratch = Scratch::new();
    let (archive, gtf, genome) = prepare_fixture(&scratch.0);
    let legacy_gtf = scratch.0.join("legacy.gtf");
    let legacy_report = scratch.0.join("legacy.tsv");
    let mut legacy = extend_command(&archive, &gtf, &legacy_gtf);
    legacy
        .arg("--report")
        .arg(&legacy_report)
        .arg("--genome")
        .arg(&genome);
    let legacy_output = run(legacy);
    assert!(
        String::from_utf8_lossy(&legacy_output.stdout).starts_with("extend: 1 of 1 genes extended")
    );

    let uniform_gtf = scratch.0.join("uniform.gtf");
    let uniform_legacy_report = scratch.0.join("uniform-legacy.tsv");
    let mut uniform = extend_command(&archive, &gtf, &uniform_gtf);
    uniform
        .arg("--report")
        .arg(&uniform_legacy_report)
        .arg("--genome")
        .arg(&genome)
        .arg("--report-format")
        .arg("json");
    let uniform_output = run(uniform);
    let result: serde_json::Value = serde_json::from_slice(&uniform_output.stdout).unwrap();
    assert_eq!(result["$schema"], "gravlax.result-envelope.v1");
    assert_eq!(result["result_schema"], "gravlax.extend.result.v1");
    assert_eq!(result["data"]["summary"]["genes_extended"], 1);
    assert_eq!(result["data"]["summary"]["total_extension_bp"], 100);
    assert!(result["provenance"]["archives"][0]
        .as_str()
        .unwrap()
        .starts_with("aie-directory-root-v2:"));
    assert!(result["provenance"]["annotation_digest"]
        .as_str()
        .unwrap()
        .starts_with("blake3:"));
    assert!(result["provenance"]["parameters"]["genome_identity"]
        .as_str()
        .unwrap()
        .starts_with("aie-genome-blake3-v1:"));

    let artifacts = table(&result, "artifacts");
    assert_eq!(artifacts["schema"]["id"], "gravlax.extend.artifacts.v1");
    assert_eq!(artifacts["schema"]["semantics"]["row_semantics"], "set");
    assert_eq!(
        artifacts["schema"]["semantics"]["key"],
        serde_json::json!(["artifact_kind", "path"])
    );
    assert_eq!(artifacts["rows"].as_array().unwrap().len(), 2);

    let extensions = table(&result, "extensions");
    assert_eq!(extensions["schema"]["id"], "gravlax.extend.genes.v1");
    assert_eq!(extensions["schema"]["semantics"]["row_semantics"], "set");
    assert_eq!(
        extensions["schema"]["semantics"]["key"],
        serde_json::json!(["gene_index"])
    );
    assert!(extensions["schema"]["semantics"]
        .get("ordered_by")
        .is_none());
    let row = &extensions["rows"][0];
    assert_eq!(row[1], "G1");
    assert_eq!(row[5], 100);
    assert_eq!(row[6], 200);
    assert_eq!(row[7], 100);
    assert_eq!(row[8], 1);
    assert_eq!(row[9], 1);

    assert_eq!(
        std::fs::read(&legacy_gtf).unwrap(),
        std::fs::read(&uniform_gtf).unwrap()
    );
    assert_eq!(
        std::fs::read(&legacy_report).unwrap(),
        std::fs::read(&uniform_legacy_report).unwrap(),
    );
    assert_eq!(
        std::fs::read_to_string(&legacy_report).unwrap(),
        concat!(
            "#gene_id\tgene_name\tchrom\tstrand\told_end\tnew_end\text_bp\tsite_umis\tsite_cells\tn_sites\tip_dropped\tclip\n",
            "G1\tGene1\tchr1\t+\t100\t200\t100\t1\t1\t1\t0\tmax\n",
        ),
    );
    assert!(
        String::from_utf8_lossy(&uniform_output.stderr).contains("extend: 1 of 1 genes extended")
    );
}

#[test]
fn uniform_extend_preflights_no_clobber_conflicts_and_malformed_flags() {
    let scratch = Scratch::new();
    let (archive, gtf, genome) = prepare_fixture(&scratch.0);

    let unstamped = scratch.0.join("unstamped.aie");
    let mut ingest = Command::new(env!("CARGO_BIN_EXE_aie"));
    ingest
        .arg("ingest-archive")
        .arg(scratch.0.join("fixture.bam"))
        .arg("--whitelist")
        .arg(scratch.0.join("whitelist.txt"))
        .arg("--out")
        .arg(&unstamped)
        .arg("--zstd-level")
        .arg("1")
        .arg("--chunk-mb")
        .arg("1");
    run(ingest);
    let unstamped_output = scratch.0.join("unstamped-output.gtf");
    let mut unverified = extend_command(&unstamped, &gtf, &unstamped_output);
    unverified
        .arg("--genome")
        .arg(&genome)
        .arg("--report-format")
        .arg("json");
    let failure = run_failure(unverified);
    assert!(String::from_utf8_lossy(&failure.stderr).contains("requires a stamped archive"));
    assert!(!unstamped_output.exists());

    let occupied = scratch.0.join("occupied.json");
    std::fs::write(&occupied, "keep me\n").unwrap();
    let guarded_gtf = scratch.0.join("guarded.gtf");
    let guarded_legacy = scratch.0.join("guarded.tsv");
    let mut guarded = extend_command(&archive, &gtf, &guarded_gtf);
    guarded
        .arg("--report")
        .arg(&guarded_legacy)
        .arg("--report-format")
        .arg("json")
        .arg("--report-output")
        .arg(&occupied);
    let failure = run_failure(guarded);
    assert!(String::from_utf8_lossy(&failure.stderr).contains("refusing to replace"));
    assert_eq!(std::fs::read_to_string(&occupied).unwrap(), "keep me\n");
    assert!(!guarded_gtf.exists());
    assert!(!guarded_legacy.exists());
    assert!(!std::fs::read_dir(&scratch.0).unwrap().any(|entry| {
        entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .contains("gravlax-stage")
    }));

    let conflict = scratch.0.join("conflict.gtf");
    let mut same_path = extend_command(&archive, &gtf, &conflict);
    same_path
        .arg("--report-format")
        .arg("json")
        .arg("--report-output")
        .arg(&conflict);
    let failure = run_failure(same_path);
    assert!(String::from_utf8_lossy(&failure.stderr).contains("must name different files"));
    assert!(!conflict.exists());

    let report_conflict_gtf = scratch.0.join("report-conflict.gtf");
    let report_conflict = scratch.0.join("report-conflict.tsv");
    let mut same_report_path = extend_command(&archive, &gtf, &report_conflict_gtf);
    same_report_path
        .arg("--report")
        .arg(&report_conflict)
        .arg("--report-format")
        .arg("tsv")
        .arg("--report-output")
        .arg(&report_conflict);
    let failure = run_failure(same_report_path);
    assert!(String::from_utf8_lossy(&failure.stderr).contains("must name different files"));
    assert!(!report_conflict_gtf.exists());
    assert!(!report_conflict.exists());

    let legacy_artifact_conflict = scratch.0.join("legacy-artifact-conflict.gtf");
    let mut same_legacy_artifact = extend_command(&archive, &gtf, &legacy_artifact_conflict);
    same_legacy_artifact
        .arg("--report")
        .arg(&legacy_artifact_conflict);
    let failure = run_failure(same_legacy_artifact);
    assert!(String::from_utf8_lossy(&failure.stderr).contains("must name different files"));
    assert!(!legacy_artifact_conflict.exists());

    let missing_parent_gtf = scratch.0.join("missing-parent.gtf");
    let missing_report = scratch.0.join("absent").join("report.json");
    let mut missing_parent = extend_command(&archive, &gtf, &missing_parent_gtf);
    missing_parent
        .arg("--report-format")
        .arg("json")
        .arg("--report-output")
        .arg(&missing_report);
    run_failure(missing_parent);
    assert!(!missing_parent_gtf.exists());
    assert!(!missing_report.exists());

    let malformed_gtf = scratch.0.join("malformed-flags.gtf");
    let mut malformed = extend_command(&archive, &gtf, &malformed_gtf);
    malformed
        .arg("--report-output")
        .arg(scratch.0.join("orphan.json"));
    let failure = run_failure(malformed);
    assert!(String::from_utf8_lossy(&failure.stderr).contains("--report-format"));
    assert!(!malformed_gtf.exists());
}

#[test]
fn uniform_extend_file_reports_are_complete_in_every_format_and_stdout_is_empty() {
    let scratch = Scratch::new();
    let (archive, gtf, _genome) = prepare_fixture(&scratch.0);
    for format in ["json", "tsv", "text"] {
        let out_gtf = scratch.0.join(format!("file-output-{format}.gtf"));
        let report = scratch.0.join(format!("result.{format}"));
        let mut command = extend_command(&archive, &gtf, &out_gtf);
        command
            .arg("--report-format")
            .arg(format)
            .arg("--report-output")
            .arg(&report);
        let output = run(command);
        assert!(output.stdout.is_empty());
        let bytes = std::fs::read(&report).unwrap();
        assert!(bytes.ends_with(b"\n"));
        let text = String::from_utf8(bytes.clone()).unwrap();
        match format {
            "json" => {
                let result: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
                assert_eq!(result["result_schema"], "gravlax.extend.result.v1");
                assert_eq!(
                    table(&result, "extensions")["selection"]["truncated"],
                    false
                );
            }
            "tsv" => {
                assert!(text.contains("# result_schema=gravlax.extend.result.v1"));
                assert!(text.contains("# table_schema=gravlax.extend.genes.v1"));
            }
            "text" => {
                assert!(text.contains("result: gravlax.extend.result.v1"));
                assert!(text.contains("table: extensions"));
            }
            _ => unreachable!(),
        }
    }
}
