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
use std::process::Command;

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("gravlax-ux-cli-{}-{nonce}", std::process::id()));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

fn write_annotation(path: &Path) {
    std::fs::write(
        path,
        concat!(
            "chr1\tX\texon\t101\t200\t.\t+\t.\tgene_id \"G1.2\"; transcript_id \"T1.4\"; exon_id \"E1.1\"; gene_name \"ALPHA\";\n",
            "chr2\tX\texon\t301\t400\t.\t+\t.\tgene_id \"G2.1\"; transcript_id \"T2.1\"; exon_id \"E2.1\"; gene_name \"DUP\";\n",
            "chr3\tX\texon\t501\t600\t.\t-\t.\tgene_id \"G3.1\"; transcript_id \"T3.1\"; exon_id \"E3.1\"; gene_name \"DUP\";\n",
        ),
    )
    .unwrap();
}

#[test]
fn resolve_emits_a_typed_envelope_and_ambiguity_leaves_no_file() {
    let scratch = Scratch::new();
    let gtf = scratch.0.join("annotation.gtf");
    write_annotation(&gtf);

    let output = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("resolve")
        .arg(&gtf)
        .args(["ALPHA", "transcript:T1"])
        .args(["--assembly", "GRCh38.p14", "--annotation", "test-release"])
        .args(["--format", "json"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "resolve failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["$schema"], "gravlax.result-envelope.v1");
    assert_eq!(value["result_schema"], "gravlax.annotation.resolve.v1");
    assert_eq!(value["producer"]["name"], "aie");
    assert_eq!(value["provenance"]["assembly"], "GRCh38.p14");
    assert_eq!(value["provenance"]["annotation"], "test-release");
    let digest = value["provenance"]["annotation_digest"].as_str().unwrap();
    assert!(digest.starts_with("blake3:"));
    assert_eq!(digest.len(), 71);
    assert_eq!(value["data"]["rows"].as_array().unwrap().len(), 2);
    assert_eq!(value["data"]["rows"][0][7], "chr1");
    assert_eq!(value["data"]["rows"][0][8], 100);
    assert_eq!(value["data"]["rows"][0][9], 200);

    let relative = Command::new(env!("CARGO_BIN_EXE_aie"))
        .current_dir(&scratch.0)
        .arg("resolve")
        .arg(&gtf)
        .arg("ALPHA")
        .args(["--assembly", "GRCh38.p14", "--annotation", "test-release"])
        .args(["--format", "tsv", "--output", "relative.tsv"])
        .output()
        .unwrap();
    assert!(relative.status.success());
    assert!(scratch.0.join("relative.tsv").is_file());

    let destination = scratch.0.join("ambiguous.json");
    let ambiguous = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("resolve")
        .arg(&gtf)
        .arg("DUP")
        .args(["--assembly", "GRCh38.p14", "--annotation", "test-release"])
        .args(["--format", "json", "--output"])
        .arg(&destination)
        .output()
        .unwrap();
    assert!(!ambiguous.status.success());
    assert!(String::from_utf8_lossy(&ambiguous.stderr).contains("ambiguous"));
    assert!(!destination.exists());
}

#[test]
fn public_help_groups_research_instruments_under_dev() {
    let root = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(root.status.success());
    let root_help = String::from_utf8_lossy(&root.stdout);
    assert!(root_help.contains("  dev"));
    assert!(!root_help.contains("  gate-c"));
    assert!(!root_help.contains("  umi-graph"));
    assert!(!root_help.contains("  em"));

    let dev = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["dev", "--help"])
        .output()
        .unwrap();
    assert!(dev.status.success());
    let dev_help = String::from_utf8_lossy(&dev.stdout);
    assert!(dev_help.contains("gate-c"));
    assert!(dev_help.contains("umi-graph"));
    assert!(dev_help.contains("em"));

    // The old path stays parseable for scripts even though it is absent from public help.
    let legacy = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["gate-c", "--help"])
        .output()
        .unwrap();
    assert!(legacy.status.success());

    let ingest = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["ingest-archive", "--help"])
        .output()
        .unwrap();
    assert!(ingest.status.success());
    let ingest_help = String::from_utf8_lossy(&ingest.stdout);
    assert!(!ingest_help.contains("frozen"));
    assert!(!ingest_help.contains("open block-size decision"));
    assert!(ingest_help.contains("authenticated .aie v2 archive"));

    let apa = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["query", "archive.aie", "apa", "--help"])
        .output()
        .unwrap();
    assert!(apa.status.success());
    let apa_help = String::from_utf8_lossy(&apa.stdout);
    assert!(!apa_help.contains("frozen"));
    assert!(!apa_help.contains("gates doc"));
    assert!(!apa_help.contains("count matrix cannot"));
}

#[test]
fn direct_apa_cli_rejects_invalid_option_combinations_before_opening_an_archive() {
    for (arguments, expected) in [
        (
            vec!["query", "missing.aie", "apa", "chr1:1-2", "--site-gap", "0"],
            "--site-gap",
        ),
        (
            vec!["query", "missing.aie", "apa", "chr1:1-2", "--drop-ip"],
            "--genome",
        ),
        (
            vec!["query", "missing.aie", "apa", "chr1:1-2", "--permute", "10"],
            "--groups",
        ),
        (
            vec![
                "query",
                "missing.aie",
                "apa-test",
                "--gtf",
                "genes.gtf",
                "--groups",
                "groups.tsv",
                "--site-gap",
                "0",
            ],
            "--site-gap",
        ),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_aie"))
            .args(arguments)
            .output()
            .unwrap();
        assert!(!output.status.success());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(expected), "unexpected error: {stderr}");
        assert!(
            !stderr.contains("failed to open"),
            "archive opened before validation: {stderr}"
        );
    }
}

#[test]
fn chemistry_recipes_are_explicit_and_annotation_free() {
    for (chemistry, umi_len, whitelist) in [
        ("10x-3p-v2", "10", "737K-august-2016.txt"),
        ("10x-3p-v3", "12", "3M-february-2018.txt"),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_aie"))
            .args(["ingest", "recipe", "--chemistry", chemistry])
            .output()
            .unwrap();
        assert!(output.status.success());
        let recipe = String::from_utf8_lossy(&output.stdout);
        assert!(recipe.contains(&format!("--soloUMIlen {umi_len}")));
        assert!(recipe.contains(whitelist));
        assert!(recipe.contains("--outSAMattributes NH HI AS nM CR CY UR UY"));
        assert!(recipe.contains("--outSAMmultNmax 50"));
        assert!(!recipe.contains("--sjdbGTFfile"));
    }
}

fn write_ingest_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(10_000).unwrap()),
        )
        .build();
    let data: Data = [
        (Tag::new(b'C', b'R'), Value::from("AAAAAAAAAAAAAAAA")),
        (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
        (Tag::new(b'U', b'R'), Value::from("ACGTACGTACGT")),
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
    ]
    .into_iter()
    .collect();
    let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
    let record = RecordBuf::builder()
        .set_name("read-1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(101).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; 50]))
        .set_quality_scores(QualityScores::from(vec![30; 50]))
        .set_data(data)
        .build();
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    writer.write_alignment_record(&header, &record).unwrap();
    writer.try_finish().unwrap();
}

fn write_broken_multimap_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(10_000).unwrap()),
        )
        .build();
    let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
    let primary_data: Data = [
        (Tag::new(b'C', b'R'), Value::from("AAAAAAAAAAAAAAAA")),
        (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
        (Tag::new(b'U', b'R'), Value::from("ACGTACGTACGT")),
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(2u8)),
    ]
    .into_iter()
    .collect();
    let primary = RecordBuf::builder()
        .set_name("multi")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(101).unwrap())
        .set_cigar(cigar.clone())
        .set_sequence(Sequence::from(vec![b'A'; 50]))
        .set_quality_scores(QualityScores::from(vec![30; 50]))
        .set_data(primary_data)
        .build();
    let secondary = RecordBuf::builder()
        .set_name("multi")
        .set_flags(Flags::SECONDARY)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(201).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; 50]))
        .set_quality_scores(QualityScores::from(vec![30; 50]))
        // Missing NH makes rows.rs treat this alternative as unique and lose it.
        .build();
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    writer.write_alignment_record(&header, &primary).unwrap();
    writer.write_alignment_record(&header, &secondary).unwrap();
    writer.try_finish().unwrap();
}

#[test]
fn ingest_check_scans_tags_order_and_whitelist_with_parseable_json() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("input.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    write_ingest_bam(&bam);
    std::fs::write(&whitelist, "AAAAAAAAAAAAAAAA\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["ingest", "check"])
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .args(["--chemistry", "10x-3p-v3", "--format", "json"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "preflight failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schema"], "gravlax.ingest.preflight.v1");
    assert_eq!(value["ok"], true);
    assert_eq!(value["counts"]["records"], 1);
    assert!(value["checks"]
        .as_array()
        .unwrap()
        .iter()
        .any(|check| check["id"] == "raw_tags" && check["status"] == "pass"));

    std::fs::write(&whitelist, "not-a-barcode\n").unwrap();
    let invalid = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["ingest", "check"])
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .args(["--format", "json"])
        .output()
        .unwrap();
    assert!(!invalid.status.success());
    let report: serde_json::Value = serde_json::from_slice(&invalid.stdout).unwrap();
    assert_eq!(report["ok"], false);

    let broken_bam = scratch.0.join("broken-multimap.bam");
    write_broken_multimap_bam(&broken_bam);
    std::fs::write(&whitelist, "AAAAAAAAAAAAAAAA\n").unwrap();
    let broken = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["ingest", "check"])
        .arg(&broken_bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .args(["--format", "json"])
        .output()
        .unwrap();
    assert!(!broken.status.success());
    let report: serde_json::Value = serde_json::from_slice(&broken.stdout).unwrap();
    let secondary_check = report["checks"]
        .as_array()
        .unwrap()
        .iter()
        .find(|check| check["id"] == "secondary_alignments")
        .unwrap();
    assert_eq!(secondary_check["status"], "fail");
}
