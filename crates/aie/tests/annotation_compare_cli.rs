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
use std::sync::atomic::{AtomicU64, Ordering};

const BARCODE: &str = "AAAAAAAAAAAAAAAA";
static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

struct Fixture {
    root: PathBuf,
    archive: PathBuf,
    annotation_a: PathBuf,
    annotation_b: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let id = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
        let root = std::env::temp_dir().join(format!(
            "gravlax-annotation-compare-cli-{}-{id}",
            std::process::id()
        ));
        std::fs::create_dir(&root).unwrap();
        let bam = root.join("evidence.bam");
        let whitelist = root.join("whitelist.txt");
        let archive = root.join("evidence.aie");
        let annotation_a = root.join("a.gtf");
        let annotation_b = root.join("b.gtf");
        write_fixture_bam(&bam);
        std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
        std::fs::write(
            &annotation_a,
            concat!(
                "chr1\ttest\texon\t1\t300\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA\"; gene_name \"GA\";\n",
                "chr1\ttest\texon\t301\t2000000\t.\t+\t.\tgene_id \"GX.1\"; transcript_id \"TXA\"; gene_name \"GX\";\n"
            ),
        )
        .unwrap();
        std::fs::write(
            &annotation_b,
            concat!(
                "chr1\ttest\texon\t1\t300\t.\t+\t.\tgene_id \"GB.2\"; transcript_id \"TB\"; gene_name \"GB\";\n",
                "chr1\ttest\texon\t301\t2000000\t.\t+\t.\tgene_id \"GX.2\"; transcript_id \"TXB\"; gene_name \"GX\";\n"
            ),
        )
        .unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_aie"))
            .arg("ingest-archive")
            .arg(&bam)
            .arg("--whitelist")
            .arg(&whitelist)
            .arg("--out")
            .arg(&archive)
            .arg("--zstd-level")
            .arg("1")
            .arg("--chunk-mb")
            .arg("1")
            .output()
            .unwrap();
        assert_success(&output);
        Self {
            root,
            archive,
            annotation_a,
            annotation_b,
        }
    }

    fn command(&self) -> Command {
        self.command_for(&self.annotation_a, &self.annotation_b)
    }

    fn command_for(&self, annotation_a: &Path, annotation_b: &Path) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_aie"));
        command
            .arg("compare-annotations")
            .arg(&self.archive)
            .arg("--annotation-a")
            .arg(annotation_a)
            .arg("--annotation-b")
            .arg(annotation_b)
            .arg("--assembly")
            .arg("GRCh38.test")
            .arg("--annotation-a-label")
            .arg("release-A")
            .arg("--annotation-b-label")
            .arg("release-B");
        command
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.root).ok();
    }
}

fn tagged_record(name: &str, start: usize, umi: &str, flags: Flags, nh: u8) -> RecordBuf {
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
}

fn write_fixture_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let records = [
        tagged_record("u1a", 101, "AAAAAAAAAAAA", Flags::empty(), 1),
        tagged_record("u1b", 121, "AAAAAAAAAAAA", Flags::empty(), 1),
        tagged_record("u1c", 101, "AAAAAAAAAAAC", Flags::empty(), 1),
        tagged_record("u2", 151, "CCCCCCCCCCCC", Flags::empty(), 1),
        tagged_record("mm1", 201, "GGGGGGGGGGGG", Flags::empty(), 3),
        tagged_record("mm1", 201, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        tagged_record("mm1", 401, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        tagged_record("far", 1_100_101, "TTTTTTTTTTTT", Flags::empty(), 1),
    ];
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in records {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn help_exposes_the_paired_comparison_contract() {
    let output = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["compare-annotations", "--help"])
        .output()
        .unwrap();
    assert_success(&output);
    let help = String::from_utf8(output.stdout).unwrap();
    for required in [
        "--annotation-a <A>",
        "--annotation-b <B>",
        "--annotation-a-digest <DIGEST>",
        "--gene-key <GENE_KEY>",
        "--max-molecule-witnesses",
        "--allow-identical",
        "count-deltas",
        "contributing-causes",
    ] {
        assert!(help.contains(required), "help omitted {required:?}\n{help}");
    }
}

#[test]
fn json_tsv_and_text_are_pure_and_preserve_exactness_semantics() {
    let fixture = Fixture::new();

    let mut json_command = fixture.command();
    json_command.args(["--format", "json"]);
    let json_output = json_command.output().unwrap();
    assert_success(&json_output);
    assert!(json_output.stderr.is_empty());
    assert_eq!(json_output.stdout.first(), Some(&b'{'));
    let value: serde_json::Value = serde_json::from_slice(&json_output.stdout).unwrap();
    assert_eq!(value["$schema"], "gravlax.result-envelope.v1");
    assert_eq!(value["result_schema"], "gravlax.annotation.compare.v1");
    assert_eq!(value["data"]["summary"]["archive_passes"], 1);
    assert_eq!(
        value["data"]["semantics"]["final_count_deltas_are_exact"],
        true
    );
    assert_eq!(
        value["data"]["semantics"]["contributing_causes_are_nonexclusive"],
        true
    );
    assert_eq!(
        value["data"]["semantics"]["contributing_causes_are_additive_attributions"],
        false
    );
    assert_eq!(
        value["data"]["semantics"]["annotation_order_tie_break_is_biological_change"],
        false
    );
    assert!(value["data"]["semantics"]["annotation_order_tie_break"]
        .as_str()
        .unwrap()
        .contains("annotation_order_tie_break_changed"));
    let annotations = value["provenance"]["annotations"].as_array().unwrap();
    assert_eq!(annotations[0]["role"], "before");
    assert_eq!(annotations[0]["annotation"], "release-A");
    assert!(annotations[0]["digest"]
        .as_str()
        .unwrap()
        .starts_with("blake3:"));
    assert_eq!(annotations[1]["role"], "after");
    assert_eq!(annotations[1]["annotation"], "release-B");
    assert!(value["provenance"]["archives"][0]
        .as_str()
        .unwrap()
        .starts_with("aie-directory-root-v2:"));
    for (name, schema) in [
        ("count_deltas", "gravlax.annotation.compare.count-deltas.v1"),
        (
            "class_transitions",
            "gravlax.annotation.compare.class-transitions.v1",
        ),
        (
            "contributing_causes",
            "gravlax.annotation.compare.contributing-causes.v1",
        ),
        ("witnesses", "gravlax.annotation.compare.witnesses.v1"),
    ] {
        assert_eq!(value["data"][name]["schema"]["id"], schema);
        assert!(value["data"][name]["rows"].is_array());
        assert!(value["data"][name]["schema"]["fields"]
            .as_array()
            .unwrap()
            .iter()
            .any(|field| field["name"] == "cell_barcode"));
    }
    assert_eq!(value["data"]["count_deltas"]["rows"][0][1], BARCODE);

    let mut tsv_command = fixture.command();
    tsv_command.args(["--format", "tsv", "--table", "count-deltas"]);
    let tsv_output = tsv_command.output().unwrap();
    assert_success(&tsv_output);
    assert!(tsv_output.stderr.is_empty());
    let tsv = String::from_utf8(tsv_output.stdout).unwrap();
    assert!(tsv.starts_with("# envelope_schema=gravlax.result-envelope.v1\n"));
    assert!(tsv.contains("# result_schema=gravlax.annotation.compare.count-deltas.v1\n"));
    assert!(tsv.contains("\ncell\tcell_barcode\tcomparison_gene_id\tannotation_a_gene_id\t"));
    assert!(!tsv.contains("Annotation comparison:"));

    let text_output = fixture.command().output().unwrap();
    assert_success(&text_output);
    assert!(text_output.stderr.is_empty());
    let text = String::from_utf8(text_output.stdout).unwrap();
    assert!(text.contains("final count deltas are exact after two independent final collapses"));
    assert!(text.contains("Causes are non-exclusive observed state changes"));
    assert!(text.contains("not additive count-delta attributions"));
    assert!(text.contains("replay-method artifacts, not biological structural changes"));
    assert!(text.contains("witnesses are bounded"));
    assert!(text.contains(BARCODE));

    let mut missing_table = fixture.command();
    missing_table.args(["--format", "tsv"]);
    let missing_table = missing_table.output().unwrap();
    assert!(!missing_table.status.success());
    assert!(missing_table.stdout.is_empty());
    assert!(
        String::from_utf8_lossy(&missing_table.stderr).contains("--format tsv requires --table")
    );

    let mut forbidden_table = fixture.command();
    forbidden_table.args(["--format", "json", "--table", "witnesses"]);
    let forbidden_table = forbidden_table.output().unwrap();
    assert!(!forbidden_table.status.success());
    assert!(forbidden_table.stdout.is_empty());
    assert!(String::from_utf8_lossy(&forbidden_table.stderr)
        .contains("--table is only valid with --format tsv"));
}

#[test]
fn identical_digest_gate_and_atomic_output_failures_leave_no_partial() {
    let fixture = Fixture::new();

    let mut identical = fixture.command_for(&fixture.annotation_a, &fixture.annotation_a);
    identical.args(["--format", "json"]);
    let identical = identical.output().unwrap();
    assert!(!identical.status.success());
    assert!(identical.stdout.is_empty());
    assert!(
        String::from_utf8_lossy(&identical.stderr).contains("identical observed content digest")
    );

    let mut allowed = fixture.command_for(&fixture.annotation_a, &fixture.annotation_a);
    allowed.args(["--allow-identical", "--format", "json"]);
    let allowed = allowed.output().unwrap();
    assert_success(&allowed);
    let allowed: serde_json::Value = serde_json::from_slice(&allowed.stdout).unwrap();
    assert_eq!(
        allowed["data"]["count_deltas"]["rows"],
        serde_json::json!([])
    );
    assert_eq!(
        allowed["data"]["class_transitions"]["rows"],
        serde_json::json!([])
    );

    let failed_output = fixture.root.join("digest-failure.json");
    let mut digest_failure = fixture.command();
    digest_failure
        .args([
            "--annotation-a-digest",
            "blake3:0000000000000000000000000000000000000000000000000000000000000000",
            "--format",
            "json",
            "-o",
        ])
        .arg(&failed_output);
    let digest_failure = digest_failure.output().unwrap();
    assert!(!digest_failure.status.success());
    assert!(digest_failure.stdout.is_empty());
    assert!(!failed_output.exists());
    assert!(String::from_utf8_lossy(&digest_failure.stderr).contains("digest does not match"));

    let occupied = fixture.root.join("occupied.json");
    std::fs::write(&occupied, b"sentinel\n").unwrap();
    let mut no_clobber = fixture.command();
    no_clobber.args(["--format", "json", "-o"]).arg(&occupied);
    let no_clobber = no_clobber.output().unwrap();
    assert!(!no_clobber.status.success());
    assert!(no_clobber.stdout.is_empty());
    assert_eq!(std::fs::read(&occupied).unwrap(), b"sentinel\n");
    assert!(String::from_utf8_lossy(&no_clobber.stderr).contains("refusing to overwrite"));

    let installed = fixture.root.join("comparison.json");
    let mut success = fixture.command();
    success.args(["--format", "json", "-o"]).arg(&installed);
    let success = success.output().unwrap();
    assert_success(&success);
    assert!(success.stdout.is_empty());
    let installed_value: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&installed).unwrap()).unwrap();
    assert_eq!(
        installed_value["result_schema"],
        "gravlax.annotation.compare.v1"
    );
    assert!(!std::fs::read_dir(&fixture.root).unwrap().any(|entry| entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .contains(".tmp.")));
}
