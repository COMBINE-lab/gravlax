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
const SECOND_BARCODE: &str = "CCCCCCCCCCCCCCCC";

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("gravlax-cli-golden-{}-{nonce}", std::process::id()));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

fn tagged_record(name: &str, start: usize, umi: &str, flags: Flags, nh: u8) -> RecordBuf {
    tagged_record_for_barcode(name, start, umi, flags, nh, BARCODE)
}

fn tagged_record_for_barcode(
    name: &str,
    start: usize,
    umi: &str,
    flags: Flags,
    nh: u8,
    barcode: &str,
) -> RecordBuf {
    let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
    let data: Data = [
        (Tag::new(b'C', b'R'), Value::from(barcode)),
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

fn spliced_record(name: &str, start: usize, umi: &str) -> RecordBuf {
    let cigar: Cigar = [
        Op::new(Kind::Match, 25),
        Op::new(Kind::Skip, 100),
        Op::new(Kind::Match, 25),
    ]
    .into_iter()
    .collect();
    let data: Data = [
        (Tag::new(b'C', b'R'), Value::from(BARCODE)),
        (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
        (Tag::new(b'U', b'R'), Value::from(umi)),
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
    ]
    .into_iter()
    .collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(start).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; 50]))
        .set_quality_scores(QualityScores::from(vec![30; 50]))
        .set_data(data)
        .build()
}

fn spliced_record_with_skip(name: &str, start: usize, skip: usize, umi: &str) -> RecordBuf {
    let cigar: Cigar = [
        Op::new(Kind::Match, 25),
        Op::new(Kind::Skip, skip),
        Op::new(Kind::Match, 25),
    ]
    .into_iter()
    .collect();
    let data: Data = [
        (Tag::new(b'C', b'R'), Value::from(BARCODE)),
        (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
        (Tag::new(b'U', b'R'), Value::from(umi)),
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
    ]
    .into_iter()
    .collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(start).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; 50]))
        .set_quality_scores(QualityScores::from(vec![30; 50]))
        .set_data(data)
        .build()
}

fn double_spliced_record(name: &str, start: usize, umi: &str) -> RecordBuf {
    double_spliced_record_with_flags(name, start, umi, Flags::empty())
}

fn double_spliced_record_with_flags(
    name: &str,
    start: usize,
    umi: &str,
    flags: Flags,
) -> RecordBuf {
    let cigar: Cigar = [
        Op::new(Kind::Match, 25),
        Op::new(Kind::Skip, 100),
        Op::new(Kind::Match, 25),
        Op::new(Kind::Skip, 100),
        Op::new(Kind::Match, 25),
    ]
    .into_iter()
    .collect();
    let data: Data = [
        (Tag::new(b'C', b'R'), Value::from(BARCODE)),
        (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
        (Tag::new(b'U', b'R'), Value::from(umi)),
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
    ]
    .into_iter()
    .collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(start).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; 75]))
        .set_quality_scores(QualityScores::from(vec![30; 75]))
        .set_data(data)
        .build()
}

const GRAPH_TEST_UMIS: [&str; 18] = [
    "AAAAAAAAAAAA",
    "CCCCCCCCCCCC",
    "GGGGGGGGGGGG",
    "TTTTTTTTTTTT",
    "ACACACACACAC",
    "CACACACACACA",
    "AGAGAGAGAGAG",
    "GAGAGAGAGAGA",
    "ATATATATATAT",
    "TATATATATATA",
    "CGCGCGCGCGCG",
    "GCGCGCGCGCGC",
    "CTCTCTCTCTCT",
    "TCTCTCTCTCTC",
    "GTGTGTGTGTGT",
    "TGTGTGTGTGTG",
    "ACGTACGTACGT",
    "TGCATGCATGCA",
];

fn write_cohort_graph_bam(path: &Path, left_umis: usize, right_umis: usize) {
    assert!(left_umis <= 9 && right_umis <= 9);
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for (index, &umi) in GRAPH_TEST_UMIS.iter().take(left_umis).enumerate() {
        writer
            .write_alignment_record(
                &header,
                &spliced_record(&format!("left-{index}"), 101, umi),
            )
            .unwrap();
        if index == 0 {
            writer
                .write_alignment_record(
                    &header,
                    &spliced_record("left-repeated-representative", 101, umi),
                )
                .unwrap();
        }
    }
    for (index, &umi) in GRAPH_TEST_UMIS.iter().skip(9).take(right_umis).enumerate() {
        writer
            .write_alignment_record(
                &header,
                &spliced_record(&format!("right-{index}"), 226, umi),
            )
            .unwrap();
    }
    for record in [
        tagged_record("mm-a", 1_001, "AACCAACCAACC", Flags::empty(), 2),
        tagged_record("mm-b", 1_201, "AACCAACCAACC", Flags::SECONDARY, 2),
        tagged_record("edge-a", 1_401, "AAGGAAGGAAGG", Flags::empty(), 1),
        tagged_record("edge-b", 1_421, "AAGGAAGGAAGG", Flags::empty(), 1),
        tagged_record("edge-c", 1_401, "CCTTCCTTCCTT", Flags::empty(), 1),
    ] {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn write_fixture_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let records = [
        // Two positions for one UMI exercise the two-representative reduction.
        tagged_record("u1a", 101, "AAAAAAAAAAAA", Flags::empty(), 1),
        tagged_record("u1b", 121, "AAAAAAAAAAAA", Flags::empty(), 1),
        // A one-mismatch neighbour exercises preservation of the archive UMI graph.
        tagged_record("u1c", 101, "AAAAAAAAAAAC", Flags::empty(), 1),
        tagged_record("u2", 151, "CCCCCCCCCCCC", Flags::empty(), 1),
        // Primary + secondary exercise the multimapper pattern streams.
        tagged_record("mm1", 201, "GGGGGGGGGGGG", Flags::empty(), 3),
        // Duplicate anchor geometry occurs in real STAR output; the archive intentionally
        // preserves the duplicate alternative without assigning it a separate identity.
        tagged_record("mm1", 201, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        tagged_record("mm1", 401, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        // Same UMI beyond a physical chunk boundary: streaming must merge the class globally.
        tagged_record("u1far", 1_100_101, "AAAAAAAAAAAA", Flags::empty(), 1),
    ];

    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in records {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn write_spliced_fixture_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    let records = [
        spliced_record("s1", 101, "AAAAAAAAAAAA"),
        spliced_record("s2", 101, "CCCCCCCCCCCC"),
        tagged_record("mm1", 201, "GGGGGGGGGGGG", Flags::empty(), 3),
        tagged_record("mm1", 201, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        tagged_record("mm1", 401, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        tagged_record("edge1a", 601, "ACACACACACAC", Flags::empty(), 1),
        tagged_record("edge1b", 621, "ACACACACACAC", Flags::empty(), 1),
        tagged_record("edge1c", 601, "ACACACACACAT", Flags::empty(), 1),
        tagged_record("far", 1_100_101, "TTTTTTTTTTTT", Flags::empty(), 1),
    ];
    for record in &records {
        writer.write_alignment_record(&header, record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn write_event_fixture_bam(path: &Path) {
    write_event_fixture_bam_variant(path, None);
}

fn write_event_fixture_bam_variant(path: &Path, extra: Option<(usize, &str)>) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let records = [
        // Inclusion-only cassette class: both flanks, counted once.
        spliced_record_with_skip("include-left", 101, 100, "AAAAAAAAAAAA"),
        spliced_record_with_skip("both-left", 101, 100, "GGGGGGGGGGGG"),
        // Exclusion-only and both-side classes use the 125-350 skipping junction.
        spliced_record_with_skip("exclude", 101, 225, "CCCCCCCCCCCC"),
        spliced_record_with_skip("both-skip", 101, 225, "GGGGGGGGGGGG"),
        spliced_record_with_skip("include-right", 226, 100, "AAAAAAAAAAAA"),
        spliced_record_with_skip("both-right", 226, 100, "GGGGGGGGGGGG"),
        // Keep all archive streams nonempty; these placements are outside the event window.
        tagged_record("mm", 1_001, "TTTTTTTTTTTT", Flags::empty(), 2),
        tagged_record("mm", 1_201, "TTTTTTTTTTTT", Flags::SECONDARY, 2),
    ];
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in &records {
        writer.write_alignment_record(&header, record).unwrap();
    }
    if let Some((start, umi)) = extra {
        writer
            .write_alignment_record(
                &header,
                &spliced_record_with_skip("sample-specific", start, 100, umi),
            )
            .unwrap();
    }
    writer.try_finish().unwrap();
}

fn run(mut command: Command) -> Output {
    let command_debug = format!("{command:?}");
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "command failed: {command_debug}\nstatus: {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn uniform_table<'a>(result: &'a serde_json::Value, name: &str) -> &'a Vec<serde_json::Value> {
    result["data"]["tables"]
        .as_array()
        .unwrap()
        .iter()
        .find(|table| table["name"] == name)
        .unwrap_or_else(|| panic!("uniform result has no {name:?} table"))["rows"]
        .as_array()
        .unwrap()
}

fn assert_uniform_table_contract(result: &serde_json::Value, name: &str, schema: &str) {
    let table = result["data"]["tables"]
        .as_array()
        .unwrap()
        .iter()
        .find(|table| table["name"] == name)
        .unwrap_or_else(|| panic!("uniform result has no {name:?} table"));
    assert_eq!(table["schema"]["id"], schema);
    assert!(table["schema"]["semantics"]["row_semantics"].is_string());
    assert!(table["selection"]["available_rows"].is_number());
    assert!(table["selection"]["emitted_rows"].is_number());
    assert!(table["selection"]["truncated"].is_boolean());
}

#[test]
fn bam_archive_and_post_correction_bam_replays_match() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("fixture.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let barcodes = scratch.0.join("barcodes.tsv");
    let gtf = scratch.0.join("fixture.gtf");
    let reverse_gtf = scratch.0.join("fixture-reverse.gtf");
    let archive = scratch.0.join("fixture.aie");
    let archive_out = scratch.0.join("archive-matrix");
    let compiled_annotation = scratch.0.join("fixture.aic");
    let compiled_archive_out = scratch.0.join("compiled-archive-matrix");
    let eager_archive_out = scratch.0.join("eager-archive-matrix");
    let bam_out = scratch.0.join("bam-matrix");
    let molecule_bam = scratch.0.join("post-correction.bam");
    let uniform_molecule_bam = scratch.0.join("post-correction-uniform.bam");
    let molecule_bam_report = scratch.0.join("post-correction-uniform.json");
    let molecule_bam_out = scratch.0.join("post-correction-matrix");
    let reverse_archive_out = scratch.0.join("reverse-archive-matrix");
    let reverse_eager_out = scratch.0.join("reverse-eager-matrix");
    let reverse_bam_out = scratch.0.join("reverse-bam-matrix");
    let fai = scratch.0.join("fixture.fa.fai");

    write_fixture_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    std::fs::write(&barcodes, format!("{BARCODE}\n")).unwrap();
    std::fs::write(&fai, "chr1\t2000000\t0\t80\t81\n").unwrap();
    std::fs::write(
        &gtf,
        "chr1\ttest\texon\t1\t2000000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
    )
    .unwrap();
    std::fs::write(
        &reverse_gtf,
        "chr1\ttest\texon\t1\t2000000\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
    )
    .unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .arg("--zstd-level")
        .arg("1")
        .arg("--chunk-mb")
        .arg("1");
    run(ingest);

    // The packed streaming EM and the retained v1 eager implementation must agree on logical
    // class partitioning.  This tiny fixture has no recoverable multi-gene class, so its complete
    // stdout is an exact golden comparison rather than a floating-point tolerance check.
    let mut packed_em = Command::new(bin);
    packed_em
        .arg("em")
        .arg(&archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--mask")
        .arg("0");
    let packed_em = run(packed_em);
    let mut eager_em = Command::new(bin);
    eager_em
        .arg("em")
        .arg(&archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--mask")
        .arg("0")
        .arg("--eager");
    let eager_em = run(eager_em);
    assert_eq!(
        packed_em.stdout, eager_em.stdout,
        "packed and eager EM summaries differ"
    );

    let groups = scratch.0.join("groups.tsv");
    let metrics = scratch.0.join("hierarchical-metrics.json");
    let candidates = scratch.0.join("candidate-genes.txt");
    std::fs::write(&groups, format!("{BARCODE}\tfixture-group\n")).unwrap();
    let mut hierarchical_em = Command::new(bin);
    hierarchical_em
        .arg("em")
        .arg(&archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--mask")
        .arg("0.2")
        .arg("--groups")
        .arg(&groups)
        .arg("--metrics-json")
        .arg(&metrics)
        .arg("--candidate-genes-out")
        .arg(&candidates);
    let hierarchical_em = run(hierarchical_em);
    let hierarchical_stdout = String::from_utf8_lossy(&hierarchical_em.stdout);
    assert!(hierarchical_stdout.contains("mode group"));
    assert!(hierarchical_stdout.contains("mode hierarchical"));
    assert!(hierarchical_stdout.contains("mode convex"));
    assert!(hierarchical_stdout.contains("mode dirichlet-proxy"));
    assert!(hierarchical_stdout.contains("mode depth-hybrid"));
    let metrics: serde_json::Value =
        serde_json::from_slice(&std::fs::read(metrics).unwrap()).unwrap();
    assert_eq!(metrics["schema_version"], 2);
    assert_eq!(metrics["scored_cells"], 1);
    assert_eq!(metrics["modes"].as_array().unwrap().len(), 9);
    assert_eq!(metrics["convex_candidate_normalized"], true);
    assert_eq!(metrics["convex_leave_one_cell_out_group"], true);
    assert_eq!(metrics["dirichlet_proxy_candidate_normalized"], true);
    assert_eq!(metrics["dirichlet_proxy_leave_one_cell_out_group"], true);
    assert_eq!(metrics["depth_hybrid_scale"], 8.0);
    assert_eq!(metrics["depth_hybrid_power"], 8.0);
    assert_eq!(metrics["depth_hybrid_only"], false);
    assert_eq!(metrics["paired_comparisons"].as_array().unwrap().len(), 5);
    assert!(metrics["paired_comparisons"][0]["negative_log_loss_ci95"].is_array());
    for mode in metrics["modes"].as_array().unwrap() {
        let strata = mode["evidence_depth_strata"].as_array().unwrap();
        assert_eq!(strata.len(), 4);
        let depth_n: u64 = strata.iter().map(|row| row["n"].as_u64().unwrap()).sum();
        assert_eq!(depth_n, mode["n"].as_u64().unwrap());
    }
    assert!(candidates.is_file());

    let convex_metrics = scratch.0.join("convex-only-metrics.json");
    let mut convex_only = Command::new(bin);
    convex_only
        .arg("em")
        .arg(&archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--mask")
        .arg("0.2")
        .arg("--groups")
        .arg(&groups)
        .arg("--convex-only")
        .arg("--metrics-json")
        .arg(&convex_metrics);
    run(convex_only);
    let convex_metrics: serde_json::Value =
        serde_json::from_slice(&std::fs::read(convex_metrics).unwrap()).unwrap();
    assert_eq!(convex_metrics["convex_only"], true);
    assert_eq!(convex_metrics["modes"].as_array().unwrap().len(), 1);
    assert_eq!(convex_metrics["modes"][0]["name"], "convex");

    let dirichlet_metrics = scratch.0.join("dirichlet-only-metrics.json");
    let mut dirichlet_only = Command::new(bin);
    dirichlet_only
        .arg("em")
        .arg(&archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--mask")
        .arg("0.2")
        .arg("--groups")
        .arg(&groups)
        .arg("--dirichlet-only")
        .arg("--metrics-json")
        .arg(&dirichlet_metrics);
    run(dirichlet_only);
    let dirichlet_metrics: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dirichlet_metrics).unwrap()).unwrap();
    assert_eq!(dirichlet_metrics["dirichlet_proxy_only"], true);
    assert_eq!(dirichlet_metrics["modes"].as_array().unwrap().len(), 1);
    assert_eq!(dirichlet_metrics["modes"][0]["name"], "dirichlet-proxy");

    let hybrid_metrics = scratch.0.join("hybrid-only-metrics.json");
    let mut hybrid_only = Command::new(bin);
    hybrid_only
        .arg("em")
        .arg(&archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--mask")
        .arg("0.2")
        .arg("--groups")
        .arg(&groups)
        .arg("--hybrid-only")
        .arg("--hybrid-depth-scale")
        .arg("32")
        .arg("--hybrid-depth-power")
        .arg("4")
        .arg("--metrics-json")
        .arg(&hybrid_metrics);
    run(hybrid_only);
    let hybrid_metrics: serde_json::Value =
        serde_json::from_slice(&std::fs::read(hybrid_metrics).unwrap()).unwrap();
    assert_eq!(hybrid_metrics["depth_hybrid_only"], true);
    assert_eq!(hybrid_metrics["depth_hybrid_scale"], 32.0);
    assert_eq!(hybrid_metrics["depth_hybrid_power"], 4.0);
    assert!(hybrid_metrics["paired_comparisons"].as_array().unwrap().is_empty());
    assert_eq!(hybrid_metrics["depth_hybrid_fixed_cell_weight"], 0.2);
    assert_eq!(hybrid_metrics["depth_hybrid_fixed_group_weight"], 0.6);
    assert_eq!(hybrid_metrics["depth_hybrid_adaptive_cell_prior"], 64.0);
    assert_eq!(hybrid_metrics["modes"].as_array().unwrap().len(), 1);
    assert_eq!(hybrid_metrics["modes"][0]["name"], "depth-hybrid");

    for (args, label, message) in [
        (
            vec!["--mask", "1.1"],
            "--mask",
            "--mask must be finite and between 0 and 1",
        ),
        (
            vec!["--alpha=-1"],
            "--alpha",
            "--alpha must be finite and non-negative",
        ),
        (
            vec!["--group-alpha=-1"],
            "--group-alpha",
            "--group-alpha must be finite and non-negative",
        ),
        (
            vec!["--global-alpha=-1"],
            "--global-alpha",
            "--global-alpha must be finite and non-negative",
        ),
        (
            vec!["--convex-cell-weight=1.1"],
            "--convex-cell-weight",
            "--convex-cell-weight must be finite and between 0 and 1",
        ),
        (
            vec!["--convex-group-weight=-0.1"],
            "--convex-group-weight",
            "--convex-group-weight must be finite and between 0 and 1",
        ),
        (
            vec!["--convex-cell-weight=0.6", "--convex-group-weight=0.5"],
            "convex simplex",
            "--convex-cell-weight plus --convex-group-weight must not exceed 1",
        ),
        (
            vec!["--convex-group-prior=-1"],
            "--convex-group-prior",
            "--convex-group-prior must be finite and non-negative",
        ),
        (
            vec!["--dirichlet-cell-prior=-1"],
            "--dirichlet-cell-prior",
            "--dirichlet-cell-prior must be finite and non-negative",
        ),
        (
            vec!["--dirichlet-group-prior=-1"],
            "--dirichlet-group-prior",
            "--dirichlet-group-prior must be finite and non-negative",
        ),
        (
            vec!["--hybrid-depth-scale=0"],
            "--hybrid-depth-scale",
            "--hybrid-depth-scale must be finite and greater than zero",
        ),
        (
            vec!["--hybrid-depth-power=0"],
            "--hybrid-depth-power",
            "--hybrid-depth-power must be finite and greater than zero",
        ),
    ] {
        let mut command = Command::new(bin);
        command
            .arg("em")
            .arg(&archive)
            .arg("--gtf")
            .arg(&gtf)
            .args(args);
        let output = command.output().unwrap();
        assert!(
            !output.status.success(),
            "invalid {label} unexpectedly succeeded"
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(message),
            "invalid {label} did not report the locked validation error"
        );
    }

    let mut replay_archive = Command::new(bin);
    replay_archive
        .arg("replay-rows")
        .arg(&archive)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&archive_out);
    run(replay_archive);

    let mut compile_annotation = Command::new(bin);
    compile_annotation
        .arg("compile-annotation")
        .arg(&gtf)
        .arg("--out")
        .arg(&compiled_annotation);
    run(compile_annotation);
    let mut replay_compiled = Command::new(bin);
    replay_compiled
        .arg("replay-rows")
        .arg(&archive)
        .arg("--gtf")
        .arg(&compiled_annotation)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&compiled_archive_out);
    run(replay_compiled);

    let mut replay_archive_eager = Command::new(bin);
    replay_archive_eager
        .arg("replay-rows")
        .arg(&archive)
        .arg("--eager")
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&eager_archive_out);
    run(replay_archive_eager);

    let mut replay_bam = Command::new(bin);
    replay_bam
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
        .arg(&bam_out);
    run(replay_bam);

    let mut export = Command::new(bin);
    export
        .arg("export-molecule-bam")
        .arg(&archive)
        .arg("--fai")
        .arg(&fai)
        .arg("--out")
        .arg(&molecule_bam);
    run(export);

    let mut uniform_export = Command::new(bin);
    uniform_export
        .arg("export-molecule-bam")
        .arg(&archive)
        .arg("--fai")
        .arg(&fai)
        .arg("--out")
        .arg(&uniform_molecule_bam)
        .args(["--report-format", "json", "--report-output"])
        .arg(&molecule_bam_report);
    run(uniform_export);
    assert_eq!(
        std::fs::read(&molecule_bam).unwrap(),
        std::fs::read(&uniform_molecule_bam).unwrap(),
        "uniform reporting changed the molecule-BAM artifact"
    );
    let molecule_report: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&molecule_bam_report).unwrap()).unwrap();
    assert_eq!(
        molecule_report["$schema"],
        "gravlax.result-envelope.v1"
    );
    assert_eq!(
        molecule_report["result_schema"],
        "gravlax.molecule-bam.export.result.v1"
    );
    assert!(molecule_report["provenance"]["parameters"]["fasta_index_digest"]
        .as_str()
        .unwrap()
        .starts_with("blake3:"));
    assert_uniform_table_contract(
        &molecule_report,
        "artifacts",
        "gravlax.molecule-bam.export.artifacts.v1",
    );
    assert!(uniform_table(&molecule_report, "artifacts")[0][4]
        .as_str()
        .unwrap()
        .starts_with("blake3:"));

    let alias_dir = scratch.0.join("molecule-output-alias");
    std::fs::create_dir(&alias_dir).unwrap();
    let aliased_bam = scratch.0.join("aliased-molecule.bam");
    let mut aliased_report = Command::new(bin);
    aliased_report
        .arg("export-molecule-bam")
        .arg(&archive)
        .arg("--fai")
        .arg(&fai)
        .arg("--out")
        .arg(&aliased_bam)
        .args(["--report-format", "json", "--report-output"])
        .arg(alias_dir.join("..").join("aliased-molecule.bam"));
    let aliased_report = aliased_report.output().unwrap();
    assert!(!aliased_report.status.success());
    assert!(String::from_utf8_lossy(&aliased_report.stderr)
        .contains("BAM output and report output must differ"));
    assert!(!aliased_bam.exists());

    let mut replay_molecule_bam = Command::new(bin);
    replay_molecule_bam
        .arg("replay-rows")
        .arg(&molecule_bam)
        .arg("--from-molecule-bam")
        .arg("--gtf")
        .arg(&gtf)
        .arg("--barcodes")
        .arg(&barcodes)
        .arg("--out-dir")
        .arg(&molecule_bam_out);
    run(replay_molecule_bam);

    for name in ["matrix.mtx", "features.tsv", "barcodes.tsv"] {
        let from_archive = std::fs::read(archive_out.join(name)).unwrap();
        assert_eq!(from_archive, std::fs::read(compiled_archive_out.join(name)).unwrap(),
            "{name} differs with compiled annotation input");
        let from_eager_archive = std::fs::read(eager_archive_out.join(name)).unwrap();
        assert_eq!(from_archive, from_eager_archive,
            "{name} differs between streaming and eager archive replay");
        let from_bam = std::fs::read(bam_out.join(name)).unwrap();
        assert_eq!(from_archive, from_bam, "{name} differs by replay source");
        let from_molecule_bam = std::fs::read(molecule_bam_out.join(name)).unwrap();
        assert_eq!(
            from_archive, from_molecule_bam,
            "{name} differs after post-correction BAM round trip"
        );
    }

    let matrix = std::fs::read_to_string(archive_out.join("matrix.mtx")).unwrap();
    assert!(matrix.contains("1 1 1\n"));
    assert!(matrix.lines().last().unwrap().ends_with(" 3"));

    // The same forward alignments are sense to the + transcript above and antisense to this -
    // transcript. Reverse-strand replay must therefore reproduce the complete forward result,
    // independently of streaming/eager/archive/BAM input paths.
    for (input, from_bam, eager, out) in [
        (&archive, false, false, &reverse_archive_out),
        (&archive, false, true, &reverse_eager_out),
        (&bam, true, false, &reverse_bam_out),
    ] {
        let mut command = Command::new(bin);
        command
            .arg("replay-rows")
            .arg(input)
            .arg("--gtf")
            .arg(&reverse_gtf)
            .arg("--barcodes")
            .arg(&barcodes)
            .arg("--out-dir")
            .arg(out)
            .arg("--solo-strand")
            .arg("reverse");
        if from_bam {
            command.arg("--from-bam").arg("--whitelist").arg(&whitelist);
        }
        if eager {
            command.arg("--eager");
        }
        run(command);
    }
    for name in ["matrix.mtx", "features.tsv", "barcodes.tsv"] {
        let forward = std::fs::read(archive_out.join(name)).unwrap();
        for out in [&reverse_archive_out, &reverse_eager_out, &reverse_bam_out] {
            assert_eq!(forward, std::fs::read(out.join(name)).unwrap(),
                "{name} differs under reverse-strand replay through {}", out.display());
        }
    }
}

#[test]
fn junction_listing_json_matches_point_query() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("spliced.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("spliced.aie");
    write_spliced_fixture_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .arg("--zstd-level")
        .arg("1")
        .arg("--chunk-mb")
        .arg("1");
    run(ingest);

    let mut listing = Command::new(bin);
    listing
        .arg("query")
        .arg(&archive)
        .arg("junctions")
        .arg("chr1:100-226")
        .arg("--with-cells")
        .arg("--json");
    let listing: serde_json::Value = serde_json::from_slice(&run(listing).stdout).unwrap();
    let row = &listing["junctions"][0];
    assert_eq!(listing["junctions"].as_array().unwrap().len(), 1);
    assert_eq!(row["donor"], 125);
    assert_eq!(row["acceptor"], 225);
    assert_eq!(row["supporting_children"], 2);
    assert_eq!(row["umis"], 2);
    assert_eq!(row["cells"], 1);

    let mut point = Command::new(bin);
    point
        .arg("query")
        .arg(&archive)
        .arg("junction")
        .arg("chr1:125-225")
        .arg("--top")
        .arg("0")
        .arg("--json");
    let point: serde_json::Value = serde_json::from_slice(&run(point).stdout).unwrap();
    for field in ["donor", "acceptor", "supporting_children", "posting_chunks", "umis", "cells"] {
        assert_eq!(row[field], point[field], "point/listing mismatch for {field}");
    }
    assert_eq!(row["cell_counts"], point["cell_rows"]);

    let mut uniform_listing = Command::new(bin);
    uniform_listing
        .arg("query")
        .arg(&archive)
        .arg("junctions")
        .arg("chr1:100-226")
        .arg("--with-cells")
        .args(["--format", "json"]);
    let uniform_listing: serde_json::Value =
        serde_json::from_slice(&run(uniform_listing).stdout).unwrap();
    assert_eq!(
        uniform_listing["result_schema"],
        "gravlax.query.junctions.result.v1"
    );
    assert_uniform_table_contract(
        &uniform_listing,
        "junctions",
        "gravlax.query.junctions.junctions.v1",
    );
    assert_uniform_table_contract(
        &uniform_listing,
        "counts",
        "gravlax.query.junctions.counts.v1",
    );
    let uniform_junction = &uniform_table(&uniform_listing, "junctions")[0];
    assert_eq!(uniform_junction[1], row["donor"]);
    assert_eq!(uniform_junction[2], row["acceptor"]);
    assert_eq!(uniform_junction[3], row["supporting_children"]);
    assert_eq!(uniform_junction[4], row["posting_chunks"]);
    assert_eq!(uniform_junction[5], row["umis"]);
    assert_eq!(uniform_junction[6], row["cells"]);
    assert_eq!(uniform_table(&uniform_listing, "counts")[0][4], row["umis"]);

    let uniform_query = |kind: &str, locus: &str, extra: &[&str]| {
        let mut command = Command::new(bin);
        command
            .arg("query")
            .arg(&archive)
            .arg(kind)
            .arg(locus)
            .args(extra)
            .arg("--format")
            .arg("json");
        serde_json::from_slice::<serde_json::Value>(&run(command).stdout).unwrap()
    };
    let uniform_point = uniform_query("junction", "chr1:125-225", &["--top", "0"]);
    assert_eq!(uniform_point["$schema"], "gravlax.result-envelope.v1");
    assert_eq!(
        uniform_point["result_schema"],
        "gravlax.query.junction.result.v1"
    );
    assert_eq!(uniform_point["data"]["summary"]["umis"], point["umis"]);
    assert_eq!(
        uniform_point["data"]["summary"]["archive_supporting_children"],
        point["supporting_children"]
    );
    let counts_table = &uniform_point["data"]["tables"][0];
    assert_eq!(counts_table["name"], "counts");
    assert_eq!(
        counts_table["schema"]["id"],
        "gravlax.query.junction.counts.v1"
    );
    assert_eq!(counts_table["schema"]["semantics"]["row_semantics"], "set");
    assert_eq!(
        counts_table["schema"]["semantics"]["key"],
        serde_json::json!(["aggregation", "entity"])
    );
    assert!(counts_table["schema"]["semantics"]["ordered_by"].is_null());
    assert_eq!(counts_table["selection"]["available_rows"], 1);
    assert_eq!(counts_table["selection"]["emitted_rows"], 1);
    assert_eq!(counts_table["selection"]["truncated"], false);
    assert_eq!(
        counts_table["rows"],
        serde_json::json!([["cell", BARCODE, 2, null, null]])
    );
    assert_eq!(
        uniform_point["provenance"]["parameters"]["selection_policy"]["comparator"],
        "umis descending, entity ascending (barcode)"
    );
    assert!(uniform_point["provenance"]["parameters"]
        .get("query_summary")
        .is_none());

    let uniform_region = uniform_query("region", "chr1:100-226", &["--top", "20"]);
    let mut legacy_region = Command::new(bin);
    legacy_region
        .arg("query")
        .arg(&archive)
        .arg("region")
        .arg("chr1:100-226")
        .arg("--top")
        .arg("20")
        .arg("--json");
    let legacy_region: serde_json::Value =
        serde_json::from_slice(&run(legacy_region).stdout).unwrap();
    assert_eq!(
        uniform_region["result_schema"],
        "gravlax.query.region.result.v1"
    );
    for field in ["molecules", "umis", "cells", "chunks_decoded"] {
        assert_eq!(
            uniform_region["data"]["summary"][field], legacy_region[field],
            "legacy/uniform region mismatch for {field}"
        );
    }
    assert_eq!(
        uniform_region["data"]["tables"][0]["rows"],
        serde_json::json!([[
            "cell",
            BARCODE,
            legacy_region["cell_rows"][0]["umis"],
            null,
            null
        ]])
    );

    // New uniform output standardizes top=0 to all rows. The omitted-format legacy path retains
    // its historical unscoped region behavior (top=0 emits no cell rows).
    let uniform_region_all = uniform_query("region", "chr1:100-226", &["--top", "0"]);
    let mut legacy_region_zero = Command::new(bin);
    legacy_region_zero
        .arg("query")
        .arg(&archive)
        .arg("region")
        .arg("chr1:100-226")
        .arg("--top")
        .arg("0")
        .arg("--json");
    let legacy_region_zero: serde_json::Value =
        serde_json::from_slice(&run(legacy_region_zero).stdout).unwrap();
    assert!(legacy_region_zero["cell_rows"]
        .as_array()
        .unwrap()
        .is_empty());
    assert_eq!(
        uniform_region_all["data"]["tables"][0]["selection"]["emitted_rows"],
        1
    );

    let cells = scratch.0.join("query-cells.txt");
    let groups = scratch.0.join("query-groups.tsv");
    std::fs::write(&cells, format!("{BARCODE}\n")).unwrap();
    std::fs::write(&groups, format!("{BARCODE}\tgroup-a\n")).unwrap();
    let cells_path = cells.to_str().unwrap();
    let groups_path = groups.to_str().unwrap();
    let scoped = uniform_query(
        "junction",
        "chr1:125-225",
        &["--cells", cells_path, "--agg", "cell", "--top", "0"],
    );
    assert!(
        scoped["provenance"]["parameters"]["cell_scope"]["source_content_blake3"]
            .as_str()
            .unwrap()
            .starts_with("blake3:")
    );
    assert!(
        scoped["provenance"]["parameters"]["cell_scope"]["resolved_mapping_blake3"]
            .as_str()
            .unwrap()
            .starts_with("blake3:")
    );
    let grouped = uniform_query(
        "junction",
        "chr1:125-225",
        &["--groups", groups_path, "--agg", "group", "--top", "0"],
    );
    assert_eq!(
        grouped["data"]["tables"][0]["rows"],
        serde_json::json!([["group", "group-a", 2, 1, 1]])
    );
    let bulk = uniform_query("junction", "chr1:125-225", &["--agg", "bulk", "--top", "0"]);
    assert_eq!(
        bulk["data"]["tables"][0]["rows"],
        serde_json::json!([["bulk", "bulk", 2, 1, 1]])
    );

    let mut federate = Command::new(bin);
    federate
        .arg("federate")
        .arg(&archive)
        .arg(&archive)
        .arg("chr1:125-225")
        .arg("--top")
        .arg("0")
        .args(["--format", "json"]);
    let federate: serde_json::Value = serde_json::from_slice(&run(federate).stdout).unwrap();
    assert_eq!(
        federate["result_schema"],
        "gravlax.federate.junction.result.v1"
    );
    assert_uniform_table_contract(
        &federate,
        "archives",
        "gravlax.federate.junction.archives.v1",
    );
    assert_uniform_table_contract(&federate, "counts", "gravlax.federate.junction.counts.v1");
    assert_eq!(federate["data"]["summary"]["totals"]["umis"], 4);
    assert_eq!(federate["data"]["summary"]["totals"]["cells"], 2);
    assert_eq!(uniform_table(&federate, "archives").len(), 2);
    assert_eq!(uniform_table(&federate, "counts").len(), 2);
    for archive_row in uniform_table(&federate, "archives") {
        assert_eq!(archive_row[5], point["umis"]);
        assert_eq!(archive_row[6], point["cells"]);
    }

    for format in ["text", "tsv"] {
        let mut command = Command::new(bin);
        command
            .arg("query")
            .arg(&archive)
            .arg("junction")
            .arg("chr1:125-225")
            .arg("--format")
            .arg(format);
        let output = run(command);
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.contains("gravlax.query.junction.result.v1"));
        assert!(stdout.contains("aggregation"));
        assert!(!stdout.contains("open 0."));
    }

    let output_path = scratch.0.join("uniform-junction.json");
    let mut to_file = Command::new(bin);
    to_file
        .arg("query")
        .arg(&archive)
        .arg("junction")
        .arg("chr1:125-225")
        .arg("--format")
        .arg("json")
        .arg("--output")
        .arg(&output_path);
    let written = run(to_file);
    assert!(written.stdout.is_empty());
    let original = std::fs::read(&output_path).unwrap();
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&original).unwrap()["result_schema"],
        "gravlax.query.junction.result.v1"
    );
    let mut overwrite = Command::new(bin);
    overwrite
        .arg("query")
        .arg(&archive)
        .arg("junction")
        .arg("chr1:125-225")
        .arg("--format")
        .arg("json")
        .arg("--output")
        .arg(&output_path);
    let overwrite = overwrite.output().unwrap();
    assert!(!overwrite.status.success());
    assert!(String::from_utf8_lossy(&overwrite.stderr).contains("refusing to replace"));
    assert_eq!(std::fs::read(&output_path).unwrap(), original);

    let mut missing_parent = Command::new(bin);
    missing_parent
        .arg("query")
        .arg(&archive)
        .arg("junction")
        .arg("chr1:125-225")
        .arg("--format")
        .arg("json")
        .arg("--output")
        .arg(scratch.0.join("missing-parent").join("result.json"));
    let missing_parent = missing_parent.output().unwrap();
    assert!(!missing_parent.status.success());
    assert!(String::from_utf8_lossy(&missing_parent.stderr).contains("checking parent directory"));

    for invalid in [
        vec!["--format", "json", "--json"],
        vec!["--format", "tsv", "--tsv"],
        vec!["--output", output_path.to_str().unwrap()],
    ] {
        let mut command = Command::new(bin);
        command
            .arg("query")
            .arg(&archive)
            .arg("junction")
            .arg("chr1:125-225")
            .args(invalid);
        assert!(!command.output().unwrap().status.success());
    }
    let mut mixed_artifacts = Command::new(bin);
    mixed_artifacts
        .arg("query")
        .arg(&archive)
        .arg("region")
        .arg("chr1:100-226")
        .arg("--format")
        .arg("json")
        .arg("--plot")
        .arg(scratch.0.join("mixed.svg"));
    let mixed_artifacts = mixed_artifacts.output().unwrap();
    assert!(!mixed_artifacts.status.success());
    assert!(String::from_utf8_lossy(&mixed_artifacts.stderr)
        .contains("run side-artifact export separately"));

    let mut boundary = Command::new(bin);
    boundary
        .arg("query")
        .arg(&archive)
        .arg("junctions")
        .arg("chr1:100-225")
        .arg("--json");
    let boundary: serde_json::Value = serde_json::from_slice(&run(boundary).stdout).unwrap();
    assert!(boundary["junctions"].as_array().unwrap().is_empty());

    let plan = scratch.0.join("plan.tsv");
    std::fs::write(
        &plan,
        concat!(
            "id\tkind\tlocus\n",
            "region-one\tregion\tchr1:100-226\n",
            "junction-one\tjunction\tchr1:125-225\n",
            "junction-absent\tjunction\tchr1:999-1099\n",
        ),
    )
    .unwrap();
    let batch_command = || {
        let mut command = Command::new(bin);
        command
            .arg("query")
            .arg(&archive)
            .arg("batch")
            .arg("--plan")
            .arg(&plan)
            .arg("--top")
            .arg("0");
        command
    };
    let batch_output = run(batch_command());
    let repeated = run(batch_command());
    assert_eq!(batch_output.stdout, repeated.stdout);
    let batch: serde_json::Value = serde_json::from_slice(&batch_output.stdout).unwrap();
    assert_eq!(batch["schema"], "gravlax.query.batch.v1");
    assert_eq!(batch["plan_queries"], 3);
    assert_eq!(batch["queries"][0]["id"], "region-one");
    assert_eq!(batch["queries"][1]["id"], "junction-one");
    assert_eq!(batch["queries"][2]["id"], "junction-absent");
    assert_eq!(batch["queries"][1]["supporting_children"], point["supporting_children"]);
    assert_eq!(batch["queries"][1]["posting_chunks"], point["posting_chunks"]);
    assert_eq!(batch["queries"][1]["umis"], point["umis"]);
    assert_eq!(batch["queries"][1]["cells"], point["cells"]);
    assert_eq!(batch["queries"][1]["cell_rows"], point["cell_rows"]);
    assert_eq!(batch["queries"][2]["present"], false);
    assert_eq!(batch["queries"][2]["umis"], 0);

    let mut uniform_batch = batch_command();
    uniform_batch.args(["--format", "json"]);
    let uniform_batch: serde_json::Value =
        serde_json::from_slice(&run(uniform_batch).stdout).unwrap();
    assert_eq!(
        uniform_batch["result_schema"],
        "gravlax.query.batch.result.v1"
    );
    assert_uniform_table_contract(&uniform_batch, "queries", "gravlax.query.batch.queries.v1");
    assert_uniform_table_contract(&uniform_batch, "counts", "gravlax.query.batch.counts.v1");
    let uniform_queries = uniform_table(&uniform_batch, "queries");
    assert_eq!(
        uniform_queries.len(),
        batch["queries"].as_array().unwrap().len()
    );
    for (uniform, legacy) in uniform_queries
        .iter()
        .zip(batch["queries"].as_array().unwrap())
    {
        assert_eq!(uniform[1], legacy["id"]);
        assert_eq!(uniform[2], legacy["kind"]);
        if legacy["kind"] == "junction" {
            assert_eq!(uniform[6], legacy["present"]);
        } else {
            assert_eq!(uniform[6], true);
        }
        assert_eq!(uniform[11], legacy["umis"]);
        assert_eq!(uniform[12], legacy["cells"]);
    }
    let uniform_counts = uniform_table(&uniform_batch, "counts");
    for legacy in batch["queries"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|query| !query["cell_rows"].as_array().unwrap().is_empty())
    {
        let row = uniform_counts
            .iter()
            .find(|row| row[0] == legacy["id"] && row[2] == BARCODE)
            .unwrap();
        assert_eq!(row[3], legacy["cell_rows"][0]["umis"]);
    }

    let mut region = Command::new(bin);
    region.arg("query").arg(&archive).arg("region").arg("chr1:100-226").arg("--top").arg("20");
    let region = String::from_utf8(run(region).stdout).unwrap();
    let mut region_lines = region.lines();
    let summary = region_lines.next().unwrap().split_once(": ").unwrap().1;
    let fields: Vec<&str> = summary.split_whitespace().collect();
    assert_eq!(batch["queries"][0]["molecules"], fields[0].parse::<u64>().unwrap());
    assert_eq!(batch["queries"][0]["umis"], fields[2].parse::<u64>().unwrap());
    assert_eq!(batch["queries"][0]["cells"], fields[5].parse::<u64>().unwrap());
    let standalone_cells: Vec<serde_json::Value> = region_lines
        .map(|line| {
            let mut fields = line.split_whitespace();
            serde_json::json!({
                "barcode": fields.next().unwrap(),
                "umis": fields.next().unwrap().parse::<u64>().unwrap(),
            })
        })
        .collect();
    assert_eq!(batch["queries"][0]["cell_rows"], serde_json::json!(standalone_cells));

    let duplicate_plan = scratch.0.join("duplicate-plan.tsv");
    std::fs::write(
        &duplicate_plan,
        "id\tkind\tlocus\nsame\tregion\tchr1:1-2\nsame\tjunction\tchr1:125-225\n",
    )
    .unwrap();
    let mut duplicate = Command::new(bin);
    duplicate.arg("query").arg(&archive).arg("batch").arg("--plan").arg(&duplicate_plan);
    let duplicate = duplicate.output().unwrap();
    assert!(!duplicate.status.success());
    assert!(String::from_utf8_lossy(&duplicate.stderr).contains("duplicate query id same"));

    let incomplete = scratch.0.join("duplicate-plan-result.json");
    let mut uniform_duplicate = Command::new(bin);
    uniform_duplicate
        .arg("query")
        .arg(&archive)
        .arg("batch")
        .arg("--plan")
        .arg(&duplicate_plan)
        .args(["--format", "json"])
        .arg("--output")
        .arg(&incomplete);
    let uniform_duplicate = uniform_duplicate.output().unwrap();
    assert!(!uniform_duplicate.status.success());
    assert!(String::from_utf8_lossy(&uniform_duplicate.stderr).contains("duplicate query id same"));
    assert!(!incomplete.exists());
}

#[test]
fn batch_uniform_counts_selection_reports_top_truncation() {
    let scratch = Scratch::new();
    let bam_path = scratch.0.join("batch-top.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("batch-top.aie");
    let plan = scratch.0.join("plan.tsv");
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let file = std::fs::File::create(&bam_path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in [
        tagged_record_for_barcode(
            "cell-a",
            101,
            "AAAAAAAAAAAA",
            Flags::empty(),
            1,
            BARCODE,
        ),
        tagged_record_for_barcode(
            "cell-b",
            151,
            "CCCCCCCCCCCC",
            Flags::empty(),
            1,
            SECOND_BARCODE,
        ),
        tagged_record("mm", 1_001, "GGGGGGGGGGGG", Flags::empty(), 2),
        tagged_record(
            "mm",
            1_201,
            "GGGGGGGGGGGG",
            Flags::SECONDARY,
            2,
        ),
        tagged_record("edge-a", 1_401, "TATATATATATA", Flags::empty(), 1),
        tagged_record("edge-b", 1_421, "TATATATATATA", Flags::empty(), 1),
        tagged_record("edge-c", 1_401, "TATATATATATC", Flags::empty(), 1),
    ] {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    std::fs::write(&whitelist, format!("{BARCODE}\n{SECOND_BARCODE}\n")).unwrap();
    std::fs::write(&plan, "id\tkind\tlocus\ntwo-cells\tregion\tchr1:0-1000\n").unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam_path)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .args(["--zstd-level", "1", "--chunk-mb", "1"]);
    run(ingest);

    let mut query = Command::new(bin);
    query
        .arg("query")
        .arg(&archive)
        .arg("batch")
        .arg("--plan")
        .arg(&plan)
        .args(["--top", "1", "--format", "json"]);
    let result: serde_json::Value = serde_json::from_slice(&run(query).stdout).unwrap();
    let tables = result["data"]["tables"].as_array().unwrap();
    let counts = tables.iter().find(|table| table["name"] == "counts").unwrap();
    assert_eq!(
        counts["selection"],
        serde_json::json!({
            "available_rows": 2,
            "emitted_rows": 1,
            "truncated": true,
        })
    );
    assert_eq!(counts["rows"].as_array().unwrap().len(), 1);
    let query_row = &uniform_table(&result, "queries")[0];
    assert_eq!(query_row[13], 2);
    assert_eq!(query_row[14], 1);
    assert_eq!(query_row[15], true);
    let expected_plan_digest = format!(
        "blake3:{}",
        blake3::hash(&std::fs::read(&plan).unwrap()).to_hex()
    );
    assert_eq!(
        result["provenance"]["parameters"]["plan_content_blake3"],
        expected_plan_digest
    );
}

#[test]
fn grouped_apa_tsv_keeps_statistical_diagnostics_on_stderr() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("apa.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("apa.aie");
    let groups = scratch.0.join("groups.tsv");
    write_fixture_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    std::fs::write(&groups, format!("{BARCODE}\tgroup-a\n")).unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .args(["--zstd-level", "1", "--chunk-mb", "1"]);
    run(ingest);

    let mut query = Command::new(bin);
    query
        .arg("query")
        .arg(&archive)
        .arg("apa")
        .arg("chr1:0-1000")
        .arg("--tsv")
        .arg("--groups")
        .arg(&groups);
    let output = run(query);
    let stdout = String::from_utf8(output.stdout).unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stdout.starts_with("#chrom\tstart\tend\tcleavage\tstrand\tumis\tcells\tgroup-a\n"));
    assert!(!stdout.contains("G-test"));
    assert!(stderr.contains("site x group G-test"));
    assert!(stdout.lines().all(|line| line.split('\t').count() == 8));

    let mut uniform = Command::new(bin);
    uniform
        .arg("query")
        .arg(&archive)
        .arg("apa")
        .arg("chr1:0-1000")
        .arg("--groups")
        .arg(&groups)
        .args(["--format", "json"]);
    let uniform: serde_json::Value = serde_json::from_slice(&run(uniform).stdout).unwrap();
    assert_eq!(uniform["result_schema"], "gravlax.query.apa.result.v1");
    assert_uniform_table_contract(&uniform, "sites", "gravlax.query.apa.sites.v1");
    assert_uniform_table_contract(
        &uniform,
        "group_counts",
        "gravlax.query.apa.group-counts.v1",
    );
    assert_uniform_table_contract(&uniform, "group_test", "gravlax.query.apa.group-test.v1");
    let sites = uniform_table(&uniform, "sites");
    assert_eq!(sites.len() + 1, stdout.lines().count());
    let site_umis: u64 = sites.iter().map(|row| row[6].as_u64().unwrap()).sum();
    let grouped_umis: u64 = uniform_table(&uniform, "group_counts")
        .iter()
        .map(|row| row[2].as_u64().unwrap())
        .sum();
    assert_eq!(grouped_umis, site_umis);
    assert_eq!(
        uniform["provenance"]["parameters"]["cell_scope"]["source_content_blake3"],
        format!(
            "blake3:{}",
            blake3::hash(&std::fs::read(&groups).unwrap()).to_hex()
        )
    );

    let malformed_groups = scratch.0.join("malformed-apa-groups.tsv");
    let incomplete_output = scratch.0.join("must-not-exist.json");
    std::fs::write(&malformed_groups, format!("{BARCODE}\n")).unwrap();
    let mut malformed = Command::new(bin);
    malformed
        .arg("query")
        .arg(&archive)
        .arg("apa")
        .arg("chr1:0-1000")
        .arg("--groups")
        .arg(&malformed_groups)
        .args(["--format", "json"])
        .arg("--output")
        .arg(&incomplete_output);
    let malformed = malformed.output().unwrap();
    assert!(!malformed.status.success());
    assert!(String::from_utf8_lossy(&malformed.stderr).contains("barcode<TAB>group"));
    assert!(!incomplete_output.exists());

    let mut unbound_genome = Command::new(bin);
    unbound_genome
        .arg("query")
        .arg(&archive)
        .arg("apa")
        .arg("chr1:0-1000")
        .arg("--genome")
        .arg(scratch.0.join("not-consulted.fa"))
        .args(["--format", "json"]);
    let unbound_genome = unbound_genome.output().unwrap();
    assert!(!unbound_genome.status.success());
    assert!(String::from_utf8_lossy(&unbound_genome.stderr)
        .contains("uniform APA output with --genome requires an archive stamped"));
}

#[test]
fn apa_test_uniform_contract_matches_legacy_scientific_rows() {
    let scratch = Scratch::new();
    let bam_path = scratch.0.join("apa-test.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("apa-test.aie");
    let groups = scratch.0.join("groups.tsv");
    let gtf = scratch.0.join("genes.gtf");
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let file = std::fs::File::create(&bam_path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in [
        tagged_record_for_barcode(
            "a-left",
            101,
            "AAAAAAAAAAAA",
            Flags::empty(),
            1,
            BARCODE,
        ),
        tagged_record_for_barcode(
            "a-right",
            301,
            "CCCCCCCCCCCC",
            Flags::empty(),
            1,
            BARCODE,
        ),
        tagged_record_for_barcode(
            "b-left",
            101,
            "GGGGGGGGGGGG",
            Flags::empty(),
            1,
            SECOND_BARCODE,
        ),
        tagged_record_for_barcode(
            "b-right",
            301,
            "TTTTTTTTTTTT",
            Flags::empty(),
            1,
            SECOND_BARCODE,
        ),
        tagged_record_for_barcode(
            "mm-a",
            1_101,
            "ACACACACACAC",
            Flags::empty(),
            2,
            BARCODE,
        ),
        tagged_record_for_barcode(
            "mm-b",
            1_301,
            "ACACACACACAC",
            Flags::SECONDARY,
            2,
            BARCODE,
        ),
        tagged_record_for_barcode(
            "edge-a",
            1_501,
            "TATATATATATA",
            Flags::empty(),
            1,
            BARCODE,
        ),
        tagged_record_for_barcode(
            "edge-b",
            1_521,
            "TATATATATATA",
            Flags::empty(),
            1,
            BARCODE,
        ),
        tagged_record_for_barcode(
            "edge-c",
            1_501,
            "TATATATATATC",
            Flags::empty(),
            1,
            BARCODE,
        ),
    ] {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    std::fs::write(&whitelist, format!("{BARCODE}\n{SECOND_BARCODE}\n")).unwrap();
    std::fs::write(
        &groups,
        format!("{BARCODE}\tA\n{SECOND_BARCODE}\tB\n"),
    )
    .unwrap();
    std::fs::write(
        &gtf,
        "chr1\tt\texon\t1\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
    )
    .unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam_path)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .args(["--zstd-level", "1", "--chunk-mb", "1"]);
    run(ingest);

    let command = |uniform: bool| {
        let mut command = Command::new(bin);
        command
            .arg("query")
            .arg(&archive)
            .arg("apa-test")
            .arg("--gtf")
            .arg(&gtf)
            .arg("--groups")
            .arg(&groups)
            .arg("--min-site-umis")
            .arg("1")
            .arg("--min-gene-umis")
            .arg("2");
        if uniform {
            command.args(["--format", "json"]);
        }
        command
    };
    let legacy = String::from_utf8(run(command(false)).stdout).unwrap();
    let legacy_row: Vec<&str> = legacy.lines().nth(1).unwrap().split('\t').collect();
    let uniform: serde_json::Value =
        serde_json::from_slice(&run(command(true)).stdout).unwrap();
    assert_eq!(
        uniform["result_schema"],
        "gravlax.query.apa-test.result.v1"
    );
    assert_uniform_table_contract(
        &uniform,
        "genes",
        "gravlax.query.apa-test.genes.v1",
    );
    let row = &uniform_table(&uniform, "genes")[0];
    assert_eq!(row[0], legacy_row[0]);
    assert_eq!(row[1], legacy_row[1]);
    assert_eq!(row[2], legacy_row[2].parse::<u64>().unwrap());
    assert_eq!(row[3], legacy_row[3].parse::<u64>().unwrap());
    assert_eq!(row[5], legacy_row[5].parse::<u64>().unwrap());
    assert_eq!(uniform["data"]["summary"]["genes_tested"], 1);
    let expected_annotation_digest = format!(
        "blake3:{}",
        blake3::hash(&std::fs::read(&gtf).unwrap()).to_hex()
    );
    assert_eq!(
        uniform["provenance"]["parameters"]["annotation_content_blake3"],
        expected_annotation_digest
    );
    let expected_groups_digest = format!(
        "blake3:{}",
        blake3::hash(&std::fs::read(&groups).unwrap()).to_hex()
    );
    assert_eq!(
        uniform["provenance"]["parameters"]["cell_scope"]["source_content_blake3"],
        expected_groups_digest
    );

    let mut unbound_genome = command(true);
    unbound_genome
        .arg("--genome")
        .arg(scratch.0.join("not-consulted.fa"));
    let unbound_genome = unbound_genome.output().unwrap();
    assert!(!unbound_genome.status.success());
    assert!(String::from_utf8_lossy(&unbound_genome.stderr)
        .contains("uniform APA-test output with --genome requires an archive stamped"));
}

#[test]
fn scoped_junction_set_is_exact_conservative_and_shared_across_query_paths() {
    let scratch = Scratch::new();
    let bam_path = scratch.0.join("junction-set.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("junction-set.aie");
    let groups = scratch.0.join("groups.tsv");
    let cells = scratch.0.join("cells.txt");
    let plan = scratch.0.join("plan.tsv");
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let file = std::fs::File::create(&bam_path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in [
        spliced_record("include", 101, "AAAAAAAAAAAA"),
        double_spliced_record("both", 101, "CCCCCCCCCCCC"),
        spliced_record("exclude", 226, "GGGGGGGGGGGG"),
        tagged_record("mm", 501, "ACACACACACAC", Flags::empty(), 2),
        tagged_record("mm", 701, "ACACACACACAC", Flags::SECONDARY, 2),
        // Keep every archive stream nonempty so this fixture also exercises the UMI graph.
        tagged_record("edge-a", 901, "TTTTTTTTTTTT", Flags::empty(), 1),
        tagged_record("edge-b", 921, "TTTTTTTTTTTT", Flags::empty(), 1),
        tagged_record("edge-c", 901, "TTTTTTTTTTTA", Flags::empty(), 1),
    ] {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    std::fs::write(&groups, format!("{BARCODE}\tgroup-a\n")).unwrap();
    std::fs::write(&cells, format!("{BARCODE}\n")).unwrap();
    std::fs::write(
        &plan,
        concat!(
            "id\tkind\tlocus\n",
            "region\tregion\tchr1:100-351\n",
            "include\tjunction\tchr1:125-225\n",
            "exclude\tjunction\tchr1:250-350\n",
        ),
    )
    .unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam_path)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .arg("--zstd-level")
        .arg("1")
        .arg("--chunk-mb")
        .arg("1");
    run(ingest);

    let jset = |scope_path: &Path, group_scope: bool| {
        let mut command = Command::new(bin);
        command
            .arg("query")
            .arg(&archive)
            .arg("jset")
            .arg("--include")
            .arg("chr1:125-225")
            .arg("--include")
            .arg("chr1:999-1099")
            .arg("--exclude")
            .arg("chr1:250-350")
            .arg("--json");
        if group_scope {
            command.arg("--groups");
        } else {
            command
                .arg("--cells")
                .arg(scope_path)
                .arg("--agg")
                .arg("cell");
            return serde_json::from_slice(&run(command).stdout).unwrap();
        }
        command.arg(scope_path);
        serde_json::from_slice(&run(command).stdout).unwrap()
    };
    let grouped: serde_json::Value = jset(&groups, true);
    assert_eq!(grouped["schema"], "gravlax.query.jset.v1");
    assert_eq!(grouped["scope"]["aggregation"], "group");
    assert_eq!(grouped["totals"]["include_only"], 1);
    assert_eq!(grouped["totals"]["exclude_only"], 1);
    assert_eq!(grouped["totals"]["both"], 1);
    assert_eq!(grouped["totals"]["informative_umis"], 2);
    assert_eq!(grouped["totals"]["usage_fraction"], 0.5);
    assert_eq!(grouped["group_rows"][0]["group"], "group-a");
    assert_eq!(grouped["group_rows"][0]["include_only"], 1);
    assert_eq!(grouped["inclusion_junctions"][1]["present"], false);
    assert_eq!(grouped["planning"]["unique_chunk_decodes"], 1);
    let cell_scoped: serde_json::Value = jset(&cells, false);
    assert_eq!(cell_scoped["totals"], grouped["totals"]);
    assert_eq!(cell_scoped["cell_rows"][0]["include_only"], 1);
    assert_eq!(cell_scoped["cell_rows"][0]["exclude_only"], 1);
    assert_eq!(cell_scoped["cell_rows"][0]["both"], 1);

    let mut uniform_jset = Command::new(bin);
    uniform_jset
        .arg("query")
        .arg(&archive)
        .arg("jset")
        .arg("--include")
        .arg("chr1:125-225")
        .arg("--include")
        .arg("chr1:999-1099")
        .arg("--exclude")
        .arg("chr1:250-350")
        .arg("--groups")
        .arg(&groups)
        .args(["--format", "json"]);
    let uniform_jset: serde_json::Value =
        serde_json::from_slice(&run(uniform_jset).stdout).unwrap();
    assert_eq!(
        uniform_jset["result_schema"],
        "gravlax.query.jset.result.v1"
    );
    assert_eq!(uniform_jset["data"]["summary"]["totals"], grouped["totals"]);
    assert_uniform_table_contract(
        &uniform_jset,
        "junctions",
        "gravlax.query.jset.junctions.v1",
    );
    assert_uniform_table_contract(&uniform_jset, "counts", "gravlax.query.jset.counts.v1");
    let uniform_group = &uniform_table(&uniform_jset, "counts")[0];
    assert_eq!(uniform_group[0], "group");
    assert_eq!(uniform_group[1], "group-a");
    assert_eq!(uniform_group[2], grouped["group_rows"][0]["include_only"]);
    assert_eq!(uniform_group[3], grouped["group_rows"][0]["exclude_only"]);
    assert_eq!(uniform_group[4], grouped["group_rows"][0]["both"]);

    let point = |locus: &str| {
        let mut command = Command::new(bin);
        command
            .arg("query")
            .arg(&archive)
            .arg("junction")
            .arg(locus)
            .arg("--top")
            .arg("0")
            .arg("--json");
        serde_json::from_slice::<serde_json::Value>(&run(command).stdout).unwrap()
    };
    let include = point("chr1:125-225");
    let exclude = point("chr1:250-350");
    assert_eq!(
        include["umis"].as_u64().unwrap(),
        grouped["totals"]["include_only"].as_u64().unwrap()
            + grouped["totals"]["both"].as_u64().unwrap()
    );
    assert_eq!(
        exclude["umis"].as_u64().unwrap(),
        grouped["totals"]["exclude_only"].as_u64().unwrap()
            + grouped["totals"]["both"].as_u64().unwrap()
    );

    let mut grouped_point = Command::new(bin);
    grouped_point
        .arg("query")
        .arg(&archive)
        .arg("junction")
        .arg("chr1:125-225")
        .arg("--groups")
        .arg(&groups)
        .arg("--json");
    let grouped_point: serde_json::Value =
        serde_json::from_slice(&run(grouped_point).stdout).unwrap();
    assert_eq!(grouped_point["schema"], "gravlax.query.junction.v2");
    assert_eq!(grouped_point["group_rows"][0]["umis"], 2);

    let mut grouped_region = Command::new(bin);
    grouped_region
        .arg("query")
        .arg(&archive)
        .arg("region")
        .arg("chr1:100-351")
        .arg("--groups")
        .arg(&groups)
        .arg("--json");
    let grouped_region: serde_json::Value =
        serde_json::from_slice(&run(grouped_region).stdout).unwrap();
    assert_eq!(grouped_region["group_rows"][0]["umis"], 3);

    let mut grouped_listing = Command::new(bin);
    grouped_listing
        .arg("query")
        .arg(&archive)
        .arg("junctions")
        .arg("chr1:100-351")
        .arg("--groups")
        .arg(&groups)
        .arg("--json");
    let grouped_listing: serde_json::Value =
        serde_json::from_slice(&run(grouped_listing).stdout).unwrap();
    assert_eq!(grouped_listing["schema"], "gravlax.query.junctions.v2");
    assert_eq!(
        grouped_listing["junctions"][0]["group_counts"][0]["umis"],
        2
    );

    let mut grouped_batch = Command::new(bin);
    grouped_batch
        .arg("query")
        .arg(&archive)
        .arg("batch")
        .arg("--plan")
        .arg(&plan)
        .arg("--groups")
        .arg(&groups);
    let grouped_batch: serde_json::Value =
        serde_json::from_slice(&run(grouped_batch).stdout).unwrap();
    assert_eq!(grouped_batch["schema"], "gravlax.query.batch.v2");
    assert_eq!(grouped_batch["queries"][0]["group_rows"][0]["umis"], 3);
    assert_eq!(grouped_batch["queries"][1]["group_rows"][0]["umis"], 2);

    let unknown = scratch.0.join("unknown.txt");
    std::fs::write(&unknown, "TTTTTTTTTTTTTTTT\n").unwrap();
    let mut bad_scope = Command::new(bin);
    bad_scope
        .arg("query")
        .arg(&archive)
        .arg("junction")
        .arg("chr1:125-225")
        .arg("--cells")
        .arg(&unknown);
    let bad_scope = bad_scope.output().unwrap();
    assert!(!bad_scope.status.success());
    assert!(String::from_utf8_lossy(&bad_scope.stderr).contains("is not in the archive"));

    let duplicate_cells = scratch.0.join("duplicate-cells.txt");
    std::fs::write(&duplicate_cells, format!("{BARCODE}\n{BARCODE}\n")).unwrap();
    let mut duplicate_scope = Command::new(bin);
    duplicate_scope.arg("query").arg(&archive).arg("junction").arg("chr1:125-225")
        .arg("--cells").arg(&duplicate_cells);
    let duplicate_scope = duplicate_scope.output().unwrap();
    assert!(!duplicate_scope.status.success());
    assert!(String::from_utf8_lossy(&duplicate_scope.stderr).contains("duplicate barcode"));

    let malformed_groups = scratch.0.join("malformed-groups.tsv");
    std::fs::write(&malformed_groups, format!("{BARCODE}\n")).unwrap();
    let mut malformed_scope = Command::new(bin);
    malformed_scope.arg("query").arg(&archive).arg("junction").arg("chr1:125-225")
        .arg("--groups").arg(&malformed_groups);
    let malformed_scope = malformed_scope.output().unwrap();
    assert!(!malformed_scope.status.success());
    assert!(String::from_utf8_lossy(&malformed_scope.stderr).contains("barcode<TAB>group"));

    let empty_cells = scratch.0.join("empty-cells.txt");
    std::fs::write(&empty_cells, "# no selected cells\n").unwrap();
    let mut empty_scope = Command::new(bin);
    empty_scope.arg("query").arg(&archive).arg("junction").arg("chr1:125-225")
        .arg("--cells").arg(&empty_cells);
    let empty_scope = empty_scope.output().unwrap();
    assert!(!empty_scope.status.success());
    assert!(String::from_utf8_lossy(&empty_scope.stderr).contains("contains no archive barcodes"));

    let mut group_without_mapping = Command::new(bin);
    group_without_mapping.arg("query").arg(&archive).arg("junction").arg("chr1:125-225")
        .arg("--agg").arg("group");
    let group_without_mapping = group_without_mapping.output().unwrap();
    assert!(!group_without_mapping.status.success());
    assert!(String::from_utf8_lossy(&group_without_mapping.stderr)
        .contains("--agg group requires --groups"));

    let mut overlap = Command::new(bin);
    overlap
        .arg("query")
        .arg(&archive)
        .arg("jset")
        .arg("--include")
        .arg("chr1:125-225")
        .arg("--exclude")
        .arg("chr1:125-225");
    let overlap = overlap.output().unwrap();
    assert!(!overlap.status.success());
    assert!(String::from_utf8_lossy(&overlap.stderr).contains("appears on both"));
}

#[test]
fn molecular_splice_graph_preserves_exact_paths_strands_and_scopes() {
    let scratch = Scratch::new();
    let bam_path = scratch.0.join("splice-graph.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("splice-graph.aie");
    let groups = scratch.0.join("groups.tsv");
    let cells = scratch.0.join("cells.txt");
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let file = std::fs::File::create(&bam_path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in [
        double_spliced_record("plus-both-a", 101, "AAAAAAAAAAAA"),
        // A repeated representative for the same UMI class must not add another path count.
        double_spliced_record("plus-both-b", 101, "AAAAAAAAAAAA"),
        double_spliced_record_with_flags(
            "minus-both",
            101,
            "CCCCCCCCCCCC",
            Flags::REVERSE_COMPLEMENTED,
        ),
        spliced_record("plus-left", 101, "GGGGGGGGGGGG"),
        // Keep the multimapper and UMI-neighbour streams nonempty outside the graph locus.
        tagged_record("mm-a", 1_001, "ACACACACACAC", Flags::empty(), 2),
        tagged_record("mm-b", 1_201, "ACACACACACAC", Flags::SECONDARY, 2),
        tagged_record("edge-a", 1_401, "TTTTTTTTTTTT", Flags::empty(), 1),
        tagged_record("edge-b", 1_421, "TTTTTTTTTTTT", Flags::empty(), 1),
        tagged_record("edge-c", 1_401, "TTTTTTTTTTTA", Flags::empty(), 1),
    ] {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    std::fs::write(&groups, format!("{BARCODE}\tgroup-a\n")).unwrap();
    std::fs::write(&cells, format!("{BARCODE}\n")).unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam_path)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .arg("--zstd-level")
        .arg("1")
        .arg("--chunk-mb")
        .arg("1");
    run(ingest);

    let mut cell_query = Command::new(bin);
    cell_query
        .arg("query")
        .arg(&archive)
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--cells")
        .arg(&cells)
        .arg("--agg")
        .arg("cell")
        .arg("--json");
    let cell_graph: serde_json::Value =
        serde_json::from_slice(&run(cell_query).stdout).unwrap();
    assert_eq!(cell_graph["schema"], "gravlax.query.splice-graph.v1");
    assert_eq!(cell_graph["semantics"]["lower_bound"], true);
    assert_eq!(cell_graph["semantics"]["complete_transcript_claim"], false);
    assert_eq!(cell_graph["totals"]["nodes"], 8);
    assert_eq!(cell_graph["totals"]["edges"], 4);
    assert_eq!(cell_graph["totals"]["paths"], 3);
    assert_eq!(cell_graph["totals"]["strand_path_umis"], 3);
    assert_eq!(cell_graph["planning"]["candidate_paths"], 3);
    assert_eq!(cell_graph["planning"]["scoped_distinct_umi_classes"], 3);
    assert_eq!(cell_graph["planning"]["candidate_strand_path_umis"], 3);
    assert_eq!(cell_graph["planning"]["unique_chunk_decodes"], 1);
    assert!(cell_graph["paths"]
        .as_array()
        .unwrap()
        .iter()
        .any(|path| path["strand"] == "+" && path["edge_ids"].as_array().unwrap().len() == 2));
    assert!(cell_graph["paths"]
        .as_array()
        .unwrap()
        .iter()
        .any(|path| path["strand"] == "-" && path["edge_ids"].as_array().unwrap().len() == 2));

    let edges = cell_graph["edges"].as_array().unwrap();
    let plus_left = edges
        .iter()
        .find(|edge| edge["strand"] == "+" && edge["donor"] == 125 && edge["acceptor"] == 225)
        .unwrap();
    assert_eq!(plus_left["umis"], 2);
    let minus_left = edges
        .iter()
        .find(|edge| edge["strand"] == "-" && edge["donor"] == 125 && edge["acceptor"] == 225)
        .unwrap();
    assert_eq!(minus_left["umis"], 1);
    let nodes = cell_graph["nodes"].as_array().unwrap();
    let reverse_source = minus_left["source"].as_u64().unwrap();
    let reverse_target = minus_left["target"].as_u64().unwrap();
    assert_eq!(nodes[reverse_source as usize]["coordinate"], 225);
    assert_eq!(nodes[reverse_target as usize]["coordinate"], 125);

    for edge in edges {
        let edge_id = edge["id"].as_u64().unwrap();
        let from_paths: u64 = cell_graph["paths"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|path| {
                path["edge_ids"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|id| id.as_u64() == Some(edge_id))
            })
            .map(|path| path["umis"].as_u64().unwrap())
            .sum();
        assert_eq!(edge["umis"].as_u64().unwrap(), from_paths);
    }

    let mut group_query = Command::new(bin);
    group_query
        .arg("query")
        .arg(&archive)
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--groups")
        .arg(&groups)
        .arg("--json");
    let group_graph: serde_json::Value =
        serde_json::from_slice(&run(group_query).stdout).unwrap();
    assert_eq!(group_graph["totals"], cell_graph["totals"]);
    assert_eq!(group_graph["edges"][0]["group_counts"][0]["group"], "group-a");
    assert_eq!(group_graph["edges"][0]["group_counts"][0]["umis"], 2);

    let mut uniform_graph = Command::new(bin);
    uniform_graph
        .arg("query")
        .arg(&archive)
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--groups")
        .arg(&groups)
        .args(["--format", "json"]);
    let uniform_graph: serde_json::Value =
        serde_json::from_slice(&run(uniform_graph).stdout).unwrap();
    assert_eq!(
        uniform_graph["result_schema"],
        "gravlax.query.splice-graph.result.v1"
    );
    assert_eq!(
        uniform_graph["data"]["summary"]["totals"],
        group_graph["totals"]
    );
    for (name, schema) in [
        ("nodes", "gravlax.query.splice-graph.nodes.v1"),
        ("edges", "gravlax.query.splice-graph.edges.v1"),
        ("paths", "gravlax.query.splice-graph.paths.v1"),
        ("group_counts", "gravlax.query.splice-graph.group-counts.v1"),
    ] {
        assert_uniform_table_contract(&uniform_graph, name, schema);
    }
    let uniform_paths = uniform_table(&uniform_graph, "paths");
    assert_eq!(
        uniform_paths.len(),
        group_graph["paths"].as_array().unwrap().len()
    );
    for (uniform, legacy) in uniform_paths
        .iter()
        .zip(group_graph["paths"].as_array().unwrap())
    {
        assert_eq!(uniform[1], legacy["strand"]);
        assert_eq!(uniform[4], legacy["umis"]);
        assert_eq!(uniform[5], legacy["cells"]);
    }
    let uniform_group_counts = uniform_table(&uniform_graph, "group_counts");
    let group_umis: u64 = uniform_group_counts
        .iter()
        .filter(|row| row[0] == "path" && row[2] == "group-a")
        .map(|row| row[3].as_u64().unwrap())
        .sum();
    assert_eq!(group_umis, group_graph["totals"]["strand_path_umis"]);

    let mut bulk_query = Command::new(bin);
    bulk_query
        .arg("query")
        .arg(&archive)
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--cells")
        .arg(&cells)
        .arg("--agg")
        .arg("bulk")
        .arg("--json");
    let bulk_graph: serde_json::Value =
        serde_json::from_slice(&run(bulk_query).stdout).unwrap();
    assert_eq!(bulk_graph["totals"], cell_graph["totals"]);

    let mut overflow = Command::new(bin);
    overflow
        .arg("query")
        .arg(&archive)
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--max-paths")
        .arg("2");
    let overflow = overflow.output().unwrap();
    assert!(!overflow.status.success());
    assert!(String::from_utf8_lossy(&overflow.stderr).contains("exceeding --max-paths 2"));

    for (flag, value, expected) in [
        ("--min-support", "0", "--min-support must be at least 1"),
        (
            "--min-path-umis",
            "0",
            "--min-path-umis must be at least 1",
        ),
        ("--max-paths", "0", "--max-paths must be at least 1"),
    ] {
        let mut invalid = Command::new(bin);
        invalid
            .arg("query")
            .arg(&archive)
            .arg("splice-graph")
            .arg("chr1:100-351")
            .arg(flag)
            .arg(value);
        let invalid = invalid.output().unwrap();
        assert!(!invalid.status.success());
        assert!(String::from_utf8_lossy(&invalid.stderr).contains(expected));
    }

    let mut malformed = Command::new(bin);
    malformed
        .arg("query")
        .arg(&archive)
        .arg("splice-graph")
        .arg("chr1:351-100");
    let malformed = malformed.output().unwrap();
    assert!(!malformed.status.success());
    assert!(String::from_utf8_lossy(&malformed.stderr).contains("start must be smaller"));

    let invalid_cells = scratch.0.join("invalid-cells.txt");
    std::fs::write(&invalid_cells, "TTTTTTTTTTTTTTTT\n").unwrap();
    let mut invalid_scope = Command::new(bin);
    invalid_scope
        .arg("query")
        .arg(&archive)
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--cells")
        .arg(&invalid_cells);
    let invalid_scope = invalid_scope.output().unwrap();
    assert!(!invalid_scope.status.success());
    assert!(String::from_utf8_lossy(&invalid_scope.stderr).contains("is not in the archive"));
}

#[test]
fn cohort_splice_graph_preserves_sample_rows_and_uses_sample_level_inference() {
    let scratch = Scratch::new();
    let whitelist = scratch.0.join("whitelist.txt");
    let design = scratch.0.join("design.tsv");
    let counts_design = scratch.0.join("counts-design.tsv");
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    let bin = env!("CARGO_BIN_EXE_aie");
    let mut design_text = String::from("sample\tcondition\tarchive\tcells\n");
    let mut counts_design_text = String::from("sample\tcondition\tarchive\tcells\n");
    let mut archives = Vec::new();
    for index in 0..8 {
        let sample = if index < 4 {
            format!("A{index}")
        } else {
            format!("B{}", index - 4)
        };
        let condition = if index < 4 { "control" } else { "treated" };
        let bam = scratch.0.join(format!("{sample}.bam"));
        let archive = scratch.0.join(format!("{sample}.aie"));
        write_cohort_graph_bam(
            &bam,
            if index < 4 { 9 } else { 0 },
            if index >= 4 { 9 } else { 0 },
        );
        let mut ingest = Command::new(bin);
        ingest
            .arg("ingest-archive")
            .arg(&bam)
            .arg("--whitelist")
            .arg(&whitelist)
            .arg("--out")
            .arg(&archive)
            .arg("--zstd-level")
            .arg("1")
            .arg("--chunk-mb")
            .arg("1");
        run(ingest);
        let archive_name = archive.file_name().unwrap().to_string_lossy();
        design_text.push_str(&format!("{sample}\t{condition}\t{archive_name}\t.\n"));
        counts_design_text.push_str(&format!("{sample}\tcohort\t{archive_name}\t.\n"));
        archives.push(archive);
    }
    std::fs::write(&design, design_text).unwrap();
    std::fs::write(&counts_design, counts_design_text).unwrap();

    let infer_command = |minimum_sample_umis: usize| {
        let mut command = Command::new(bin);
        command
            .arg("cohort")
            .arg("splice-graph")
            .arg("chr1:100-351")
            .arg("--design")
            .arg(&design)
            .arg("--contrast")
            .arg("control:treated")
            .arg("--min-sample-umis")
            .arg(minimum_sample_umis.to_string())
            .arg("--min-path-umis")
            .arg("1")
            .arg("--min-path-samples")
            .arg("2")
            .arg("--json");
        command
    };
    let graph: serde_json::Value =
        serde_json::from_slice(&run(infer_command(1)).stdout).unwrap();
    assert_eq!(graph["schema"], "gravlax.cohort.splice-graph.v1");
    assert_eq!(graph["semantics"]["replicate_unit"], "one unique design sample/archive row");
    assert_eq!(graph["semantics"]["cells_and_molecules_are_replicates"], false);
    assert_eq!(graph["planning"]["common_coordinate_edges"], 2);
    assert_eq!(graph["planning"]["strand_edges"], 2);
    assert_eq!(graph["planning"]["union_paths"], 2);
    assert_eq!(graph["planning"]["tested_paths"], 2);
    assert_eq!(graph["samples"].as_array().unwrap().len(), 8);
    for (index, sample) in graph["samples"].as_array().unwrap().iter().enumerate() {
        assert_eq!(sample["path_rows"].as_array().unwrap().len(), 2);
        assert_eq!(sample["edge_rows"].as_array().unwrap().len(), 2);
        assert_eq!(sample["strand_totals"][0]["umis"], 9);
        assert_eq!(sample["strand_totals"][0]["eligible"], true);
        let path_counts: Vec<u64> = sample["path_rows"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| row["umis"].as_u64().unwrap())
            .collect();
        assert_eq!(path_counts, if index < 4 { vec![9, 0] } else { vec![0, 9] });
        let edge_counts: Vec<u64> = sample["edge_rows"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| row["umis"].as_u64().unwrap())
            .collect();
        assert_eq!(edge_counts, path_counts);
    }
    let tests = graph["inference"]["tests"].as_array().unwrap();
    assert_eq!(tests.len(), 2);
    assert!(tests[0]["effect_b_minus_a"].as_f64().unwrap() < -0.9);
    assert!(tests[1]["effect_b_minus_a"].as_f64().unwrap() > 0.9);
    for test in tests {
        assert_eq!(test["sample_rows"].as_array().unwrap().len(), 8);
        assert_eq!(test["condition_a"]["eligible_samples"], 4);
        assert_eq!(test["condition_b"]["eligible_samples"], 4);
        assert!(test["beta_binomial"]["p_value"].as_f64().unwrap() < 0.01);
        assert!(test["beta_binomial"]["q_value"].as_f64().unwrap() < 0.01);
    }

    let mut uniform_graph = Command::new(bin);
    uniform_graph
        .arg("cohort")
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--design")
        .arg(&design)
        .arg("--contrast")
        .arg("control:treated")
        .arg("--min-sample-umis")
        .arg("1")
        .arg("--min-path-umis")
        .arg("1")
        .arg("--min-path-samples")
        .arg("2")
        .args(["--format", "json"]);
    let uniform_graph: serde_json::Value =
        serde_json::from_slice(&run(uniform_graph).stdout).unwrap();
    assert_eq!(
        uniform_graph["result_schema"],
        "gravlax.cohort.splice-graph.result.v1"
    );
    assert_eq!(
        uniform_graph["provenance"]["parameters"]["design_content_blake3"],
        format!(
            "blake3:{}",
            blake3::hash(&std::fs::read(&design).unwrap()).to_hex()
        )
    );
    for (name, schema) in [
        ("samples", "gravlax.cohort.splice-graph.samples.v1"),
        ("nodes", "gravlax.cohort.splice-graph.nodes.v1"),
        ("edges", "gravlax.cohort.splice-graph.edges.v1"),
        ("paths", "gravlax.cohort.splice-graph.paths.v1"),
        ("path_counts", "gravlax.cohort.splice-graph.path-counts.v1"),
        ("edge_counts", "gravlax.cohort.splice-graph.edge-counts.v1"),
        ("tests", "gravlax.cohort.splice-graph.tests.v1"),
        (
            "skipped_tests",
            "gravlax.cohort.splice-graph.skipped-tests.v1",
        ),
    ] {
        assert_uniform_table_contract(&uniform_graph, name, schema);
    }
    let path_counts = uniform_table(&uniform_graph, "path_counts");
    assert_eq!(path_counts.len(), 16);
    for (sample_index, sample) in graph["samples"].as_array().unwrap().iter().enumerate() {
        let sample_id = sample["sample"].as_str().unwrap();
        let rows: Vec<&serde_json::Value> = path_counts
            .iter()
            .filter(|row| row[0] == sample_id)
            .collect();
        assert_eq!(rows.len(), 2);
        assert_eq!(
            rows.iter()
                .map(|row| row[2].as_u64().unwrap())
                .collect::<Vec<_>>(),
            if sample_index < 4 {
                vec![9, 0]
            } else {
                vec![0, 9]
            }
        );
    }
    assert_eq!(uniform_table(&uniform_graph, "tests").len(), 2);

    let ineligible: serde_json::Value =
        serde_json::from_slice(&run(infer_command(10)).stdout).unwrap();
    assert_eq!(ineligible["planning"]["tested_paths"], 0);
    assert!(ineligible["inference"]["tests"].as_array().unwrap().is_empty());
    assert_eq!(
        ineligible["inference"]["skipped_tests"][0]["reason"],
        "insufficient_eligible_replicates"
    );

    let mut counts_only = Command::new(bin);
    counts_only
        .arg("cohort")
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--design")
        .arg(&counts_design)
        .arg("--counts-only")
        .arg("--json");
    let counts_only: serde_json::Value =
        serde_json::from_slice(&run(counts_only).stdout).unwrap();
    assert_eq!(counts_only["inference"]["enabled"], false);
    assert!(counts_only["inference"]["tests"].as_array().unwrap().is_empty());
    assert_eq!(counts_only["samples"][0]["path_rows"][0]["umis"], 9);

    let mut standalone = Command::new(bin);
    standalone
        .arg("query")
        .arg(&archives[0])
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--json");
    let standalone: serde_json::Value =
        serde_json::from_slice(&run(standalone).stdout).unwrap();
    assert_eq!(standalone["totals"]["strand_path_umis"], 9);
    assert_eq!(standalone["paths"][0]["umis"], graph["samples"][0]["path_rows"][0]["umis"]);

    let duplicate_design = scratch.0.join("duplicate-design.tsv");
    let duplicate_archive = archives[0].file_name().unwrap().to_string_lossy();
    std::fs::write(
        &duplicate_design,
        format!(
            "sample\tcondition\tarchive\tcells\nD0\ta\t{duplicate_archive}\t.\nD1\tb\t{duplicate_archive}\t.\n"
        ),
    )
    .unwrap();
    let mut duplicate = Command::new(bin);
    duplicate
        .arg("cohort")
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--design")
        .arg(&duplicate_design)
        .arg("--contrast")
        .arg("a:b");
    let duplicate = duplicate.output().unwrap();
    assert!(!duplicate.status.success());
    assert!(String::from_utf8_lossy(&duplicate.stderr).contains("reuses resolved archive"));

    let malformed_design = scratch.0.join("malformed-design.tsv");
    std::fs::write(
        &malformed_design,
        "sample\tcondition\tarchive\nA\tcontrol\tA0.aie\n",
    )
    .unwrap();
    let mut malformed = Command::new(bin);
    malformed
        .arg("cohort")
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--design")
        .arg(&malformed_design)
        .arg("--counts-only");
    let malformed = malformed.output().unwrap();
    assert!(!malformed.status.success());
    assert!(String::from_utf8_lossy(&malformed.stderr).contains("header must be exactly"));

    let mut missing_contrast = Command::new(bin);
    missing_contrast
        .arg("cohort")
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--design")
        .arg(&design);
    let missing_contrast = missing_contrast.output().unwrap();
    assert!(!missing_contrast.status.success());
    assert!(String::from_utf8_lossy(&missing_contrast.stderr).contains("requires --contrast"));

    let mut mismatched_contrast = Command::new(bin);
    mismatched_contrast
        .arg("cohort")
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--design")
        .arg(&design)
        .arg("--contrast")
        .arg("control:other");
    let mismatched_contrast = mismatched_contrast.output().unwrap();
    assert!(!mismatched_contrast.status.success());
    assert!(String::from_utf8_lossy(&mismatched_contrast.stderr)
        .contains("do not exactly match contrast"));

    let bad_cells = scratch.0.join("bad-cells.txt");
    std::fs::write(&bad_cells, "TTTTTTTTTTTTTTTT\n").unwrap();
    let unknown_cell_design = scratch.0.join("unknown-cell-design.tsv");
    std::fs::write(
        &unknown_cell_design,
        format!(
            "sample\tcondition\tarchive\tcells\nA\tcohort\t{}\t{}\nB\tcohort\t{}\t.\n",
            archives[0].file_name().unwrap().to_string_lossy(),
            bad_cells.file_name().unwrap().to_string_lossy(),
            archives[1].file_name().unwrap().to_string_lossy(),
        ),
    )
    .unwrap();
    let mut unknown_cell = Command::new(bin);
    unknown_cell
        .arg("cohort")
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--design")
        .arg(&unknown_cell_design)
        .arg("--counts-only");
    let unknown_cell = unknown_cell.output().unwrap();
    assert!(!unknown_cell.status.success());
    assert!(String::from_utf8_lossy(&unknown_cell.stderr).contains("is not in the archive"));

    for flag in [
        "--min-support",
        "--min-edge-samples",
        "--min-sample-umis",
        "--min-replicates",
        "--min-path-umis",
        "--min-path-samples",
        "--max-paths",
    ] {
        let mut invalid = Command::new(bin);
        invalid
            .arg("cohort")
            .arg("splice-graph")
            .arg("chr1:100-351")
            .arg("--design")
            .arg(&design)
            .arg("--contrast")
            .arg("control:treated")
            .arg(flag)
            .arg("0");
        let invalid = invalid.output().unwrap();
        assert!(!invalid.status.success(), "accepted {flag}=0");
        assert!(String::from_utf8_lossy(&invalid.stderr).contains("must be at least 1"));
    }

    let mut overflow = infer_command(1);
    overflow.arg("--max-paths").arg("1");
    let overflow = overflow.output().unwrap();
    assert!(!overflow.status.success());
    assert!(String::from_utf8_lossy(&overflow.stderr).contains("exceeding --max-paths 1"));

    let genome_a = scratch.0.join("genome-a.fa");
    let genome_b = scratch.0.join("genome-b.fa");
    std::fs::write(&genome_a, format!(">chr1\n{}\n", "A".repeat(1_000))).unwrap();
    std::fs::write(&genome_b, format!(">chr1\n{}\n", "C".repeat(1_000))).unwrap();
    let stamped_a = scratch.0.join("stamped-a.aie");
    let stamped_b = scratch.0.join("stamped-b.aie");
    for (archive, genome, stamped) in [
        (&archives[0], &genome_a, &stamped_a),
        (&archives[4], &genome_b, &stamped_b),
    ] {
        let mut stamp = Command::new(bin);
        stamp
            .arg("stamp-genome")
            .arg(archive)
            .arg("--genome")
            .arg(genome)
            .arg("--out")
            .arg(stamped);
        run(stamp);
    }
    let mismatched_reference_design = scratch.0.join("mismatched-reference-design.tsv");
    std::fs::write(
        &mismatched_reference_design,
        "sample\tcondition\tarchive\tcells\nA\tcontrol\tstamped-a.aie\t.\nB\ttreated\tstamped-b.aie\t.\n",
    )
    .unwrap();
    let mut mismatched_reference = Command::new(bin);
    mismatched_reference
        .arg("cohort")
        .arg("splice-graph")
        .arg("chr1:100-351")
        .arg("--design")
        .arg(&mismatched_reference_design)
        .arg("--contrast")
        .arg("control:treated");
    let mismatched_reference = mismatched_reference.output().unwrap();
    assert!(!mismatched_reference.status.success());
    assert!(String::from_utf8_lossy(&mismatched_reference.stderr)
        .contains("incompatible genome signatures"));
}

#[test]
fn event_and_cohort_engines_reduce_classes_exactly_once() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("events.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("events.aie");
    let groups = scratch.0.join("groups.tsv");
    let gtf = scratch.0.join("events.gtf");
    write_event_fixture_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    std::fs::write(&groups, format!("{BARCODE}\tgroup-a\n")).unwrap();
    std::fs::write(
        &gtf,
        concat!(
            "chr1\tt\texon\t101\t125\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
            "chr1\tt\texon\t226\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
            "chr1\tt\texon\t351\t375\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
            "chr1\tt\texon\t101\t125\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\"; gene_name \"Gene1\";\n",
            "chr1\tt\texon\t351\t375\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\"; gene_name \"Gene1\";\n",
        ),
    )
    .unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .arg("--zstd-level")
        .arg("1")
        .arg("--chunk-mb")
        .arg("1");
    run(ingest);

    let mut events = Command::new(bin);
    events
        .arg("query")
        .arg(&archive)
        .arg("events")
        .arg("chr1:100-351")
        .arg("--event-type")
        .arg("cassette")
        .arg("--min-support")
        .arg("1")
        .arg("--groups")
        .arg(&groups)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--json");
    let events: serde_json::Value = serde_json::from_slice(&run(events).stdout).unwrap();
    assert_eq!(events["schema"], "gravlax.query.events.v1");
    assert_eq!(events["events"].as_array().unwrap().len(), 1);
    let event = &events["events"][0];
    assert_eq!(event["event_type"], "cassette");
    assert_eq!(event["totals"]["include_only"], 1);
    assert_eq!(event["totals"]["exclude_only"], 1);
    assert_eq!(event["totals"]["both"], 1);
    assert_eq!(event["totals"]["informative_umis"], 2);
    assert_eq!(event["totals"]["usage_fraction"], 0.5);
    assert_eq!(event["group_rows"][0]["include_only"], 1);
    assert_eq!(event["annotation"]["fully_annotated"], true);
    assert_eq!(event["annotation"]["strand"], "+");
    assert_eq!(event["annotation"]["genes"][0]["gene_id"], "G1");
    assert_eq!(events["planning"]["unique_chunk_decodes"], 1);
    assert_eq!(events["planning"]["independent_chunk_decodes"], 3);

    let mut uniform_events = Command::new(bin);
    uniform_events
        .arg("query")
        .arg(&archive)
        .arg("events")
        .arg("chr1:100-351")
        .arg("--event-type")
        .arg("cassette")
        .arg("--min-support")
        .arg("1")
        .arg("--groups")
        .arg(&groups)
        .arg("--gtf")
        .arg(&gtf)
        .args(["--format", "json"]);
    let uniform_events: serde_json::Value =
        serde_json::from_slice(&run(uniform_events).stdout).unwrap();
    assert_eq!(
        uniform_events["result_schema"],
        "gravlax.query.events.result.v1"
    );
    let annotation_digest = format!(
        "blake3:{}",
        blake3::hash(&std::fs::read(&gtf).unwrap()).to_hex()
    );
    assert_eq!(
        uniform_events["provenance"]["parameters"]["annotation_content_blake3"],
        annotation_digest
    );
    for (name, schema) in [
        ("events", "gravlax.query.events.events.v1"),
        ("components", "gravlax.query.events.components.v1"),
        ("counts", "gravlax.query.events.counts.v1"),
    ] {
        assert_uniform_table_contract(&uniform_events, name, schema);
    }
    let uniform_event = &uniform_table(&uniform_events, "events")[0];
    assert_eq!(uniform_event[1], event["event_type"]);
    assert_eq!(uniform_event[4], event["totals"]["include_only"]);
    assert_eq!(uniform_event[5], event["totals"]["exclude_only"]);
    assert_eq!(uniform_event[6], event["totals"]["both"]);
    assert_eq!(uniform_event[7], event["totals"]["informative_umis"]);
    let uniform_count = &uniform_table(&uniform_events, "counts")[0];
    assert_eq!(uniform_count[1], "group");
    assert_eq!(uniform_count[2], "group-a");
    assert_eq!(uniform_count[3], event["group_rows"][0]["include_only"]);
    assert_eq!(uniform_count[4], event["group_rows"][0]["exclude_only"]);
    assert_eq!(uniform_count[5], event["group_rows"][0]["both"]);

    let mut jset = Command::new(bin);
    jset.arg("query")
        .arg(&archive)
        .arg("jset")
        .arg("--include")
        .arg("chr1:125-225")
        .arg("--include")
        .arg("chr1:250-350")
        .arg("--exclude")
        .arg("chr1:125-350")
        .arg("--groups")
        .arg(&groups)
        .arg("--json");
    let jset: serde_json::Value = serde_json::from_slice(&run(jset).stdout).unwrap();
    assert_eq!(event["totals"], jset["totals"]);
    assert_eq!(event["group_rows"], jset["group_rows"]);

    let cohort_command = |min_row_informative: Option<usize>| {
        let mut command = Command::new(bin);
        command
            .arg("cohort")
            .arg("events")
            .arg("chr1:100-351")
            .arg("--sample")
            .arg(format!("A={}", archive.display()))
            .arg("--sample")
            .arg(format!("B={}", archive.display()))
            .arg("--groups")
            .arg(format!("A={}", groups.display()))
            .arg("--groups")
            .arg(format!("B={}", groups.display()))
            .arg("--event-type")
            .arg("cassette")
            .arg("--min-support")
            .arg("1")
            .arg("--json");
        if let Some(minimum) = min_row_informative {
            command
                .arg("--min-row-informative")
                .arg(minimum.to_string());
        }
        command
    };
    let cohort: serde_json::Value =
        serde_json::from_slice(&run(cohort_command(None)).stdout).unwrap();
    assert_eq!(cohort["schema"], "gravlax.cohort.events.v1");
    assert!(cohort.get("min_row_informative").is_none());
    assert_eq!(cohort["events"].as_array().unwrap().len(), 1);
    for sample in cohort["events"][0]["sample_rows"].as_array().unwrap() {
        assert_eq!(sample["present"], true);
        assert_eq!(sample["totals"], event["totals"]);
        assert_eq!(sample["group_rows"], event["group_rows"]);
    }

    let mut uniform_cohort = Command::new(bin);
    uniform_cohort
        .arg("cohort")
        .arg("events")
        .arg("chr1:100-351")
        .arg("--sample")
        .arg(format!("A={}", archive.display()))
        .arg("--sample")
        .arg(format!("B={}", archive.display()))
        .arg("--groups")
        .arg(format!("A={}", groups.display()))
        .arg("--groups")
        .arg(format!("B={}", groups.display()))
        .arg("--event-type")
        .arg("cassette")
        .arg("--min-support")
        .arg("1")
        .arg("--gtf")
        .arg(&gtf)
        .args(["--format", "json"]);
    let uniform_cohort: serde_json::Value =
        serde_json::from_slice(&run(uniform_cohort).stdout).unwrap();
    assert_eq!(
        uniform_cohort["result_schema"],
        "gravlax.cohort.events.result.v1"
    );
    assert_eq!(
        uniform_cohort["provenance"]["parameters"]["annotation_content_blake3"],
        annotation_digest
    );
    for (name, schema) in [
        ("samples", "gravlax.cohort.events.samples.v1"),
        ("events", "gravlax.cohort.events.events.v1"),
        ("components", "gravlax.cohort.events.components.v1"),
        ("counts", "gravlax.cohort.events.counts.v1"),
    ] {
        assert_uniform_table_contract(&uniform_cohort, name, schema);
    }
    let cohort_counts = uniform_table(&uniform_cohort, "counts");
    for sample in ["A", "B"] {
        let total = cohort_counts
            .iter()
            .find(|row| row[1] == sample && row[3] == "total")
            .unwrap();
        let grouped: Vec<&serde_json::Value> = cohort_counts
            .iter()
            .filter(|row| row[1] == sample && row[3] == "group")
            .collect();
        for column in 5..=8 {
            assert_eq!(
                total[column].as_u64().unwrap(),
                grouped
                    .iter()
                    .map(|row| row[column].as_u64().unwrap())
                    .sum::<u64>()
            );
        }
    }

    let sparse_dir = scratch.0.join("cohort-sparse");
    let mut sparse = Command::new(bin);
    sparse
        .arg("cohort")
        .arg("events")
        .arg("chr1:100-351")
        .arg("--sample")
        .arg(format!("A={}", archive.display()))
        .arg("--sample")
        .arg(format!("B={}", archive.display()))
        .arg("--groups")
        .arg(format!("A={}", groups.display()))
        .arg("--groups")
        .arg(format!("B={}", groups.display()))
        .arg("--event-type")
        .arg("cassette")
        .arg("--min-support")
        .arg("1")
        .arg("--gtf")
        .arg(&gtf)
        .arg("--sparse-dir")
        .arg(&sparse_dir);
    let sparse_stdout: serde_json::Value = serde_json::from_slice(&run(sparse).stdout).unwrap();
    let sparse_metadata: serde_json::Value =
        serde_json::from_slice(&std::fs::read(sparse_dir.join("metadata.json")).unwrap()).unwrap();
    assert_eq!(sparse_stdout, sparse_metadata);
    assert_eq!(sparse_metadata["schema"], "gravlax.cohort.events.sparse.v1");
    assert_eq!(sparse_metadata["dimensions"]["events"], 1);
    assert_eq!(sparse_metadata["dimensions"]["samples"].as_array().unwrap().len(), 2);
    assert_eq!(sparse_metadata["output"]["rows"]["events"], 1);
    assert_eq!(sparse_metadata["output"]["rows"]["presence"], 2);
    assert_eq!(sparse_metadata["output"]["rows"]["nonzero_counts"], 4);
    let events_table = String::from_utf8(
        zstd::stream::decode_all(std::fs::File::open(sparse_dir.join("events.tsv.zst")).unwrap())
            .unwrap(),
    )
    .unwrap();
    assert!(events_table.contains("annotation_genes_json"));
    assert!(events_table.contains("\"gene_id\":\"G1\""));
    assert!(events_table.contains("\"gene_name\":\"Gene1\""));
    assert!(events_table.contains("\t+\ttrue"));
    let presence_table = String::from_utf8(
        zstd::stream::decode_all(
            std::fs::File::open(sparse_dir.join("presence.tsv.zst")).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    assert!(presence_table.lines().any(|line| line.ends_with("\tA")));
    assert!(presence_table.lines().any(|line| line.ends_with("\tB")));
    let counts_table = String::from_utf8(
        zstd::stream::decode_all(std::fs::File::open(sparse_dir.join("counts.tsv.zst")).unwrap())
            .unwrap(),
    )
    .unwrap();
    assert_eq!(counts_table.lines().count(), 5);
    assert!(counts_table.lines().any(|line| line.contains("\tA\ttotal\ttotal\t1\t1\t1\t1\t1")));
    assert!(counts_table.lines().any(|line| line.contains("\tB\tgroup\tgroup-a\t1\t1\t1\t1\t1")));

    let gated: serde_json::Value =
        serde_json::from_slice(&run(cohort_command(Some(2))).stdout).unwrap();
    assert_eq!(gated["min_row_informative"], 2);
    assert_eq!(gated["planning"]["candidate_events"], 1);
    assert_eq!(gated["planning"]["post_row_candidate_events"], 1);
    assert_eq!(gated["events"], cohort["events"]);
    let gated_out: serde_json::Value =
        serde_json::from_slice(&run(cohort_command(Some(3))).stdout).unwrap();
    assert_eq!(gated_out["planning"]["candidate_events"], 1);
    assert_eq!(gated_out["planning"]["post_row_candidate_events"], 0);
    assert!(gated_out["events"].as_array().unwrap().is_empty());

    let mut overflow = Command::new(bin);
    overflow
        .arg("query")
        .arg(&archive)
        .arg("events")
        .arg("chr1:100-351")
        .arg("--min-support")
        .arg("1")
        .arg("--max-events")
        .arg("1");
    let overflow = overflow.output().unwrap();
    assert!(!overflow.status.success());
    assert!(String::from_utf8_lossy(&overflow.stderr).contains("exceeds --max-events"));
}

#[test]
fn discovery_claim_modes_preserve_default_and_bound_residual_sites() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("discovery.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("discovery.aie");
    let gtf = scratch.0.join("old.gtf");
    write_spliced_fixture_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    std::fs::write(
        &gtf,
        concat!(
            "chr1\tt\texon\t101\t125\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
            "chr1\tt\texon\t226\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
            "chr1\tt\texon\t551\t575\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T2\";\n",
            "chr1\tt\texon\t676\t700\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T2\";\n",
        ),
    ).unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest.arg("ingest-archive").arg(&bam).arg("--whitelist").arg(&whitelist)
        .arg("--out").arg(&archive).arg("--zstd-level").arg("1").arg("--chunk-mb").arg("1");
    run(ingest);

    let discovery = |mode: Option<&str>, residual_min: Option<&str>| {
        let mut command = Command::new(bin);
        command.arg("query").arg(&archive).arg("discover").arg("--gtf").arg(&gtf)
            .arg("--min-umis").arg("1").arg("--tsv");
        if let Some(mode) = mode {
            command.arg("--claim-mode").arg(mode);
        }
        if let Some(minimum) = residual_min {
            command.arg("--residual-min-umis").arg(minimum);
        }
        run(command).stdout
    };

    let default = discovery(None, None);
    assert_eq!(default, discovery(Some("span"), None));
    let residual = String::from_utf8(discovery(Some("residual-sites"), Some("1"))).unwrap();
    let rows: Vec<Vec<&str>> = residual
        .lines()
        .map(|line| line.split('\t').collect())
        .collect();
    assert!(
        rows.iter().any(
            |row| row[1].parse::<u32>().unwrap() <= 669 && row[2].parse::<u32>().unwrap() > 669
        ),
        "intronic residual terminal site missing:\n{residual}"
    );
    assert!(rows
        .iter()
        .all(|row| row[2].parse::<u32>().unwrap() - row[1].parse::<u32>().unwrap() <= 1_001));

    let mut uniform = Command::new(bin);
    uniform
        .arg("query")
        .arg(&archive)
        .arg("discover")
        .arg("--gtf")
        .arg(&gtf)
        .arg("--min-umis")
        .arg("1")
        .arg("--claim-mode")
        .arg("residual-sites")
        .arg("--residual-min-umis")
        .arg("1")
        .args(["--format", "json"]);
    let uniform: serde_json::Value = serde_json::from_slice(&run(uniform).stdout).unwrap();
    assert_eq!(uniform["result_schema"], "gravlax.query.discover.result.v1");
    assert_eq!(
        uniform["provenance"]["parameters"]["annotation_content_blake3"],
        format!(
            "blake3:{}",
            blake3::hash(&std::fs::read(&gtf).unwrap()).to_hex()
        )
    );
    assert_uniform_table_contract(
        &uniform,
        "candidates",
        "gravlax.query.discover.candidates.v1",
    );
    let uniform_rows = uniform_table(&uniform, "candidates");
    assert_eq!(uniform_rows.len(), rows.len());
    for (uniform, legacy) in uniform_rows.iter().zip(&rows) {
        assert_eq!(uniform[0], legacy[0]);
        assert_eq!(uniform[1], legacy[1].parse::<u64>().unwrap());
        assert_eq!(uniform[2], legacy[2].parse::<u64>().unwrap());
        assert_eq!(uniform[3], legacy[3]);
        assert_eq!(uniform[4], legacy[4].parse::<u64>().unwrap());
        assert_eq!(uniform[5], legacy[5].parse::<u64>().unwrap());
    }
}

#[test]
fn collection_routes_exact_queries_and_guards_archive_identity() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("collection.bam");
    let bam_b = scratch.0.join("collection-b.bam");
    let bam_c = scratch.0.join("collection-c.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("collection.aie");
    let archive_b = scratch.0.join("collection-b.aie");
    let archive_c = scratch.0.join("collection-c.aie");
    let collection = scratch.0.join("atlas.aicollection");
    let routed_collection = scratch.0.join("atlas-routed.aicollection");
    let reverse_collection = scratch.0.join("atlas-reverse.aicollection");
    let extended_collection = scratch.0.join("atlas-extended.aicollection");
    let fresh_three_collection = scratch.0.join("atlas-fresh-three.aicollection");
    write_event_fixture_bam(&bam);
    write_event_fixture_bam_variant(&bam_b, Some((101, "ACACACACACAC")));
    write_event_fixture_bam_variant(&bam_c, Some((1_701, "AGAGAGAGAGAG")));
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    let mut ingest = Command::new(bin);
    ingest
        .arg("ingest-archive")
        .arg(&bam)
        .arg("--whitelist")
        .arg(&whitelist)
        .arg("--out")
        .arg(&archive)
        .arg("--zstd-level")
        .arg("1")
        .arg("--chunk-mb")
        .arg("1");
    run(ingest);
    for (input, out) in [(&bam_b, &archive_b), (&bam_c, &archive_c)] {
        let mut ingest = Command::new(bin);
        ingest
            .arg("ingest-archive")
            .arg(input)
            .arg("--whitelist")
            .arg(&whitelist)
            .arg("--out")
            .arg(out)
            .arg("--zstd-level")
            .arg("1")
            .arg("--chunk-mb")
            .arg("1");
        run(ingest);
    }

    let mut duplicate_source = Command::new(bin);
    duplicate_source
        .arg("collection")
        .arg("build")
        .arg("--sample")
        .arg(format!("A={}", archive.display()))
        .arg("--sample")
        .arg(format!("B={}", archive.display()))
        .arg("--allow-unstamped")
        .arg("--out")
        .arg(scratch.0.join("duplicate-source.aicollection"));
    let duplicate_source = duplicate_source.output().unwrap();
    assert!(!duplicate_source.status.success());
    assert!(String::from_utf8_lossy(&duplicate_source.stderr).contains("reuse resolved archive"));

    let byte_copy = scratch.0.join("collection-byte-copy.aie");
    std::fs::copy(&archive, &byte_copy).unwrap();
    let mut duplicate_content = Command::new(bin);
    duplicate_content
        .arg("collection")
        .arg("build")
        .arg("--sample")
        .arg(format!("A={}", archive.display()))
        .arg("--sample")
        .arg(format!("B={}", byte_copy.display()))
        .arg("--allow-unstamped")
        .arg("--out")
        .arg(scratch.0.join("duplicate-content.aicollection"));
    let duplicate_content = duplicate_content.output().unwrap();
    assert!(!duplicate_content.status.success());
    assert!(String::from_utf8_lossy(&duplicate_content.stderr).contains("identical archive content"));

    let mut forged_digests = Command::new(bin);
    forged_digests
        .arg("collection")
        .arg("build")
        .arg("--sample")
        .arg(format!("A={}", archive.display()))
        .arg("--sample")
        .arg(format!("B={}", byte_copy.display()))
        .arg("--source-digest")
        .arg(format!("A={}", "0".repeat(64)))
        .arg("--source-digest")
        .arg(format!("B={}", "1".repeat(64)))
        .arg("--allow-unstamped")
        .arg("--out")
        .arg(scratch.0.join("forged-digests.aicollection"));
    let forged_digests = forged_digests.output().unwrap();
    assert!(!forged_digests.status.success());
    assert!(String::from_utf8_lossy(&forged_digests.stderr).contains("source digest mismatch"));

    #[cfg(unix)]
    {
        let hardlink = scratch.0.join("collection-hardlink.aie");
        std::fs::hard_link(&archive, &hardlink).unwrap();
        let mut duplicate_inode = Command::new(bin);
        duplicate_inode
            .arg("collection")
            .arg("build")
            .arg("--sample")
            .arg(format!("A={}", archive.display()))
            .arg("--sample")
            .arg(format!("B={}", hardlink.display()))
            .arg("--allow-unstamped")
            .arg("--out")
            .arg(scratch.0.join("duplicate-inode.aicollection"));
        let duplicate_inode = duplicate_inode.output().unwrap();
        assert!(!duplicate_inode.status.success());
        assert!(String::from_utf8_lossy(&duplicate_inode.stderr).contains("archive inode"));
    }

    let build = |out: &Path, reverse: bool| {
        let mut command = Command::new(bin);
        command
            .arg("collection")
            .arg("build")
            .arg("--allow-unstamped");
        for (id, path) in if reverse {
            [("B", &archive_b), ("A", &archive)]
        } else {
            [("A", &archive), ("B", &archive_b)]
        } {
            command.arg("--sample").arg(format!("{id}={}", path.display()));
        }
        command.arg("--out").arg(out);
        run(command);
    };
    build(&collection, false);
    build(&reverse_collection, true);
    let build_report = scratch.0.join("atlas-routed.build.json");
    let mut routed_build = Command::new(bin);
    routed_build
        .arg("collection")
        .arg("build")
        .arg("--sample")
        .arg(format!("A={}", archive.display()))
        .arg("--sample")
        .arg(format!("B={}", archive_b.display()))
        .arg("--allow-unstamped")
        .arg("--shape-routes")
        .arg("--out")
        .arg(&routed_collection)
        .arg("--format")
        .arg("json")
        .arg("--output")
        .arg(&build_report);
    let routed_build = run(routed_build);
    assert!(routed_build.stdout.is_empty());
    let build_report: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&build_report).unwrap()).unwrap();
    assert_eq!(
        build_report["result_schema"],
        "gravlax.collection.build.result.v1"
    );
    assert_eq!(build_report["data"]["summary"]["new_archives"], 2);
    assert_eq!(
        std::fs::read(&collection).unwrap(),
        std::fs::read(&reverse_collection).unwrap(),
        "collection build depends on --sample order"
    );
    let base_bytes = std::fs::read(&collection).unwrap();
    let mut extend = Command::new(bin);
    extend
        .arg("collection")
        .arg("build")
        .arg("--base")
        .arg(&collection)
        .arg("--sample")
        .arg(format!("0C={}", archive_c.display()))
        .arg("--allow-unstamped")
        .arg("--out")
        .arg(&extended_collection);
    run(extend);
    assert_eq!(
        std::fs::read(&collection).unwrap(),
        base_bytes,
        "incremental build rewrote its immutable base"
    );
    let mut fresh_three = Command::new(bin);
    fresh_three
        .arg("collection")
        .arg("build")
        .arg("--sample")
        .arg(format!("A={}", archive.display()))
        .arg("--sample")
        .arg(format!("B={}", archive_b.display()))
        .arg("--sample")
        .arg(format!("0C={}", archive_c.display()))
        .arg("--allow-unstamped")
        .arg("--out")
        .arg(&fresh_three_collection);
    run(fresh_three);
    let json_command = |arguments: &[&str]| -> serde_json::Value {
        let mut command = Command::new(bin);
        command.args(arguments);
        serde_json::from_slice(&run(command).stdout).unwrap()
    };
    let uniform_json_command = |arguments: &[&str]| -> serde_json::Value {
        let mut command = Command::new(bin);
        command.args(arguments).arg("--format").arg("json");
        serde_json::from_slice(&run(command).stdout).unwrap()
    };
    let uniform_table = |value: &serde_json::Value, name: &str| -> Vec<serde_json::Value> {
        let table = value["data"]["tables"]
            .as_array()
            .unwrap()
            .iter()
            .find(|table| table["name"] == name)
            .unwrap();
        let fields = table["schema"]["fields"]
            .as_array()
            .unwrap()
            .iter()
            .map(|field| field["name"].as_str().unwrap())
            .collect::<Vec<_>>();
        table["rows"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| {
                serde_json::Value::Object(
                    fields
                        .iter()
                        .zip(row.as_array().unwrap())
                        .map(|(field, value)| ((*field).to_owned(), value.clone()))
                        .collect(),
                )
            })
            .collect()
    };
    let collection_path = collection.to_str().unwrap();
    let extended = json_command(&[
        "collection",
        "junction",
        extended_collection.to_str().unwrap(),
        "chr1:125-225",
        "--json",
        "--top",
        "0",
    ]);
    let compact = json_command(&[
        "collection",
        "junction",
        fresh_three_collection.to_str().unwrap(),
        "chr1:125-225",
        "--json",
        "--top",
        "0",
    ]);
    assert_eq!(extended["samples"], compact["samples"]);
    assert_eq!(extended["totals"], compact["totals"]);
    assert_eq!(extended["support_upper_bound"], compact["support_upper_bound"]);
    for subcommand in [
        vec!["region", "chr1:100-351", "--json", "--top", "9999"],
        vec![
            "jset",
            "--include",
            "chr1:125-225",
            "--include",
            "chr1:250-350",
            "--exclude",
            "chr1:125-350",
            "--json",
            "--top",
            "0",
        ],
    ] {
        let mut extension_args = vec![
            "collection",
            subcommand[0],
            extended_collection.to_str().unwrap(),
        ];
        extension_args.extend_from_slice(&subcommand[1..]);
        let mut compact_args = vec![
            "collection",
            subcommand[0],
            fresh_three_collection.to_str().unwrap(),
        ];
        compact_args.extend_from_slice(&subcommand[1..]);
        let extension_value = json_command(&extension_args);
        let compact_value = json_command(&compact_args);
        assert_eq!(extension_value["samples"], compact_value["samples"]);
        assert_eq!(extension_value["totals"], compact_value["totals"]);
    }

    for subcommand in [
        vec!["junction", "chr1:125-225", "--top", "0"],
        vec!["region", "chr1:100-351", "--top", "9999"],
        vec![
            "jset",
            "--include",
            "chr1:125-225",
            "--include",
            "chr1:250-350",
            "--exclude",
            "chr1:125-350",
            "--top",
            "0",
        ],
    ] {
        let mut extension_args = vec![
            "collection",
            subcommand[0],
            extended_collection.to_str().unwrap(),
        ];
        extension_args.extend_from_slice(&subcommand[1..]);
        let mut compact_args = vec![
            "collection",
            subcommand[0],
            fresh_three_collection.to_str().unwrap(),
        ];
        compact_args.extend_from_slice(&subcommand[1..]);
        let extension_value = uniform_json_command(&extension_args);
        let compact_value = uniform_json_command(&compact_args);
        assert_eq!(
            extension_value["data"]["summary"].get("umis"),
            compact_value["data"]["summary"].get("umis")
        );
        assert_eq!(
            extension_value["data"]["summary"].get("molecules"),
            compact_value["data"]["summary"].get("molecules")
        );
        assert_eq!(
            extension_value["data"]["summary"].get("totals"),
            compact_value["data"]["summary"].get("totals")
        );
        assert_eq!(
            uniform_table(&extension_value, "samples"),
            uniform_table(&compact_value, "samples")
        );
        assert_eq!(
            uniform_table(&extension_value, "cells"),
            uniform_table(&compact_value, "cells")
        );
    }

    let indexed_junction = json_command(&[
        "collection",
        "junction",
        collection_path,
        "chr1:125-225",
        "--json",
        "--explain",
        "--top",
        "0",
    ]);
    assert_eq!(indexed_junction["planning"]["archives_opened"], 2);
    assert_eq!(indexed_junction["planning"]["archive_catalogue_sections_read"], 0);
    for sample in indexed_junction["samples"].as_array().unwrap() {
        let source = if sample["sample"] == "A" {
            archive.to_str().unwrap()
        } else {
            archive_b.to_str().unwrap()
        };
        let naive_junction = json_command(&[
            "query",
            source,
            "junction",
            "chr1:125-225",
            "--json",
            "--top",
            "0",
        ]);
        assert_eq!(sample["umis"], naive_junction["umis"]);
        assert_eq!(sample["cells"], naive_junction["cells"]);
        assert_eq!(sample["top_cells"], naive_junction["cell_rows"]);
    }
    let uniform_junction = uniform_json_command(&[
        "collection",
        "junction",
        collection_path,
        "chr1:125-225",
        "--explain",
        "--top",
        "0",
    ]);
    let routed_uniform_junction = uniform_json_command(&[
        "collection",
        "junction",
        routed_collection.to_str().unwrap(),
        "chr1:125-225",
        "--explain",
        "--top",
        "0",
    ]);
    assert_eq!(
        uniform_junction["data"]["summary"]["umis"],
        indexed_junction["totals"]["umis"]
    );
    assert_eq!(
        uniform_junction["data"]["summary"]["cells"],
        indexed_junction["totals"]["cells"]
    );
    let scientific_sample_rows = |result: &serde_json::Value| {
        uniform_table(result, "samples")
            .into_iter()
            .map(|row| {
                serde_json::json!({
                    "sample": row["sample"],
                    "present": row["present"],
                    "supporting_children": row["supporting_children"],
                    "umis": row["umis"],
                    "cells": row["cells"],
                })
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(
        scientific_sample_rows(&uniform_junction),
        scientific_sample_rows(&routed_uniform_junction)
    );
    assert_eq!(
        uniform_table(&uniform_junction, "cells"),
        uniform_table(&routed_uniform_junction, "cells")
    );
    assert_eq!(
        uniform_junction["data"]["summary"]["planning"]["fallback_archives"],
        2
    );
    assert_eq!(
        routed_uniform_junction["data"]["summary"]["planning"]["routed_archives"],
        2
    );
    let uniform_samples = uniform_table(&uniform_junction, "samples");
    let uniform_cells = uniform_table(&uniform_junction, "cells");
    for sample in indexed_junction["samples"].as_array().unwrap() {
        let uniform_sample = uniform_samples
            .iter()
            .find(|row| row["sample"] == sample["sample"])
            .unwrap();
        for field in ["present", "supporting_children", "umis", "cells"] {
            assert_eq!(uniform_sample[field], sample[field]);
        }
        let cells = uniform_cells
            .iter()
            .filter(|row| row["sample"] == sample["sample"])
            .map(|row| serde_json::json!({"barcode": row["barcode"], "umis": row["umis"]}))
            .collect::<Vec<_>>();
        assert_eq!(cells, sample["top_cells"].as_array().unwrap().clone());
    }
    let support_pruned = json_command(&[
        "collection",
        "junction",
        collection_path,
        "chr1:125-225",
        "--min-support",
        "1000000",
        "--json",
    ]);
    assert_eq!(support_pruned["planning"]["support_bound_pruned"], true);
    assert_eq!(support_pruned["planning"]["archives_opened"], 0);
    assert_eq!(support_pruned["totals"]["umis"], 0);

    let indexed_region = json_command(&[
        "collection",
        "region",
        collection_path,
        "chr1:100-351",
        "--json",
        "--top",
        "9999",
    ]);
    for sample in indexed_region["samples"].as_array().unwrap() {
        let source = if sample["sample"] == "A" {
            archive.to_str().unwrap()
        } else {
            archive_b.to_str().unwrap()
        };
        let naive_region = json_command(&[
            "query",
            source,
            "region",
            "chr1:100-351",
            "--json",
            "--top",
            "9999",
        ]);
        assert_eq!(sample["molecules"], naive_region["molecules"]);
        assert_eq!(sample["umis"], naive_region["umis"]);
        assert_eq!(sample["cells"], naive_region["cells"]);
        assert_eq!(sample["top_cells"], naive_region["cell_rows"]);
    }
    let uniform_region = uniform_json_command(&[
        "collection",
        "region",
        collection_path,
        "chr1:100-351",
        "--top",
        "9999",
    ]);
    assert_eq!(
        uniform_region["data"]["summary"]["molecules"],
        indexed_region["totals"]["molecules"]
    );
    assert_eq!(
        uniform_region["data"]["summary"]["umis"],
        indexed_region["totals"]["umis"]
    );
    assert_eq!(
        uniform_region["data"]["summary"]["cells"],
        indexed_region["totals"]["cells"]
    );
    let uniform_region_samples = uniform_table(&uniform_region, "samples");
    for sample in indexed_region["samples"].as_array().unwrap() {
        let row = uniform_region_samples
            .iter()
            .find(|row| row["sample"] == sample["sample"])
            .unwrap();
        for field in ["present", "molecules", "umis", "cells"] {
            assert_eq!(row[field], sample[field]);
        }
    }

    let indexed_jset = json_command(&[
        "collection",
        "jset",
        collection_path,
        "--include",
        "chr1:125-225",
        "--include",
        "chr1:250-350",
        "--exclude",
        "chr1:125-350",
        "--json",
        "--top",
        "0",
    ]);
    for sample in indexed_jset["samples"].as_array().unwrap() {
        let source = if sample["sample"] == "A" {
            archive.to_str().unwrap()
        } else {
            archive_b.to_str().unwrap()
        };
        let naive_jset = json_command(&[
            "query",
            source,
            "jset",
            "--include",
            "chr1:125-225",
            "--include",
            "chr1:250-350",
            "--exclude",
            "chr1:125-350",
            "--json",
            "--top",
            "0",
        ]);
        assert_eq!(sample["totals"], naive_jset["totals"]);
        assert_eq!(sample["top_cells"], naive_jset["cell_rows"]);
    }
    let uniform_jset = uniform_json_command(&[
        "collection",
        "jset",
        collection_path,
        "--include",
        "chr1:125-225",
        "--include",
        "chr1:250-350",
        "--exclude",
        "chr1:125-350",
        "--top",
        "0",
    ]);
    assert_eq!(
        uniform_jset["data"]["summary"]["totals"],
        indexed_jset["totals"]
    );
    let uniform_jset_samples = uniform_table(&uniform_jset, "samples");
    let uniform_jset_cells = uniform_table(&uniform_jset, "cells");
    for sample in indexed_jset["samples"].as_array().unwrap() {
        let row = uniform_jset_samples
            .iter()
            .find(|row| row["sample"] == sample["sample"])
            .unwrap();
        for field in [
            "include_only",
            "exclude_only",
            "both",
            "informative_umis",
            "usage_fraction",
        ] {
            assert_eq!(row[field], sample["totals"][field]);
        }
        let cells = uniform_jset_cells
            .iter()
            .filter(|row| row["sample"] == sample["sample"])
            .map(|row| {
                serde_json::json!({
                    "barcode": row["barcode"],
                    "include_only": row["include_only"],
                    "exclude_only": row["exclude_only"],
                    "both": row["both"],
                    "informative_umis": row["informative_umis"],
                    "usage_fraction": row["usage_fraction"],
                })
            })
            .collect::<Vec<_>>();
        assert_eq!(cells, sample["top_cells"].as_array().unwrap().clone());
    }

    let inspect_report = scratch.0.join("atlas.inspect.json");
    let mut inspect = Command::new(bin);
    inspect
        .arg("collection")
        .arg("inspect")
        .arg(&routed_collection)
        .arg("--format")
        .arg("json")
        .arg("--output")
        .arg(&inspect_report);
    let inspect = run(inspect);
    assert!(inspect.stdout.is_empty());
    let inspect: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&inspect_report).unwrap()).unwrap();
    assert_eq!(
        inspect["result_schema"],
        "gravlax.collection.inspect.result.v1"
    );
    assert_eq!(uniform_table(&inspect, "archives").len(), 2);
    assert!(!uniform_table(&inspect, "shape_route_blocks").is_empty());

    let occupied = scratch.0.join("occupied-result.json");
    std::fs::write(&occupied, b"keep\n").unwrap();
    let mut occupied_query = Command::new(bin);
    occupied_query
        .arg("collection")
        .arg("junction")
        .arg(&collection)
        .arg("chr1:125-225")
        .arg("--format")
        .arg("json")
        .arg("--output")
        .arg(&occupied);
    let occupied_query = occupied_query.output().unwrap();
    assert!(!occupied_query.status.success());
    assert!(String::from_utf8_lossy(&occupied_query.stderr).contains("refusing to replace"));
    assert_eq!(std::fs::read(&occupied).unwrap(), b"keep\n");

    let reject_collection = |path: &Path| {
        let mut inspect = Command::new(bin);
        inspect.arg("collection").arg("inspect").arg(path);
        assert!(!inspect.output().unwrap().status.success());
        let result = path.with_extension("uniform.json");
        let mut uniform = Command::new(bin);
        uniform
            .arg("collection")
            .arg("inspect")
            .arg(path)
            .arg("--format")
            .arg("json")
            .arg("--output")
            .arg(&result);
        assert!(!uniform.output().unwrap().status.success());
        assert!(
            !result.exists(),
            "malformed input installed a uniform result"
        );
    };
    let valid_bytes = std::fs::read(&collection).unwrap();
    let unknown_version = scratch.0.join("unknown-version.aicollection");
    let mut changed = valid_bytes.clone();
    changed[8..12].copy_from_slice(&999u32.to_le_bytes());
    std::fs::write(&unknown_version, changed).unwrap();
    reject_collection(&unknown_version);
    let truncated = scratch.0.join("truncated.aicollection");
    std::fs::write(&truncated, &valid_bytes[..valid_bytes.len() - 1]).unwrap();
    reject_collection(&truncated);
    let corrupt = scratch.0.join("corrupt.aicollection");
    let mut changed = valid_bytes;
    let sections = u32::from_le_bytes(changed[12..16].try_into().unwrap()) as usize;
    let mut payload = 48usize;
    for _ in 0..sections {
        let name_len = changed[payload] as usize;
        payload += 1 + name_len + 8 + 8 + 32;
    }
    changed[payload] ^= 0x80;
    std::fs::write(&corrupt, changed).unwrap();
    reject_collection(&corrupt);

    let corrupt_directory = scratch.0.join("corrupt-directory.aicollection");
    let mut changed = std::fs::read(&collection).unwrap();
    changed[49] ^= 0x01;
    std::fs::write(&corrupt_directory, changed).unwrap();
    reject_collection(&corrupt_directory);

    // Even a root-authenticated directory must use the one canonical textual segment name:
    // otherwise lookup and full inspection could disagree about which bin a payload represents.
    let noncanonical_directory = scratch.0.join("noncanonical-directory.aicollection");
    let mut changed = std::fs::read(&collection).unwrap();
    let sections = u32::from_le_bytes(changed[12..16].try_into().unwrap()) as usize;
    let mut cursor = 48usize;
    let mut junction_name = None;
    for _ in 0..sections {
        let name_len = changed[cursor] as usize;
        let name_start = cursor + 1;
        if changed[name_start..name_start + name_len].starts_with(b"j.")
            && junction_name.is_none()
        {
            junction_name = Some((cursor, name_start + 2));
        }
        cursor = name_start + name_len + 48;
    }
    let (length_offset, insertion_offset) = junction_name.unwrap();
    changed.insert(insertion_offset, b'0');
    changed[length_offset] += 1;
    let directory_end = cursor + 1;
    let mut hasher = blake3::Hasher::new();
    hasher.update(&changed[..16]);
    hasher.update(&changed[48..directory_end]);
    changed[16..48].copy_from_slice(hasher.finalize().as_bytes());
    std::fs::write(&noncanonical_directory, changed).unwrap();
    reject_collection(&noncanonical_directory);
    let mut noncanonical_query = Command::new(bin);
    noncanonical_query
        .arg("collection")
        .arg("junction")
        .arg(&noncanonical_directory)
        .arg("chr1:125-225")
        .arg("--json");
    assert!(!noncanonical_query.output().unwrap().status.success());

    // A same-length change with its mtime restored still changes Unix ctime, so the default
    // open-file identity guard fails before stale routes can be used.
    #[cfg(unix)]
    {
        use std::os::unix::fs::FileExt as _;
        let modified = std::fs::metadata(&archive).unwrap().modified().unwrap();
        let changed = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(&archive)
            .unwrap();
        let mut byte = [0u8; 1];
        changed.read_exact_at(&mut byte, 8).unwrap();
        byte[0] ^= 0x01;
        changed.write_all_at(&byte, 8).unwrap();
        changed
            .set_times(std::fs::FileTimes::new().set_modified(modified))
            .unwrap();
    }
    let mut stale = Command::new(bin);
    stale.arg("collection").arg("inspect").arg(&collection);
    let stale = stale.output().unwrap();
    assert!(!stale.status.success());
    assert!(String::from_utf8_lossy(&stale.stderr).contains("archive identity changed"));
}
