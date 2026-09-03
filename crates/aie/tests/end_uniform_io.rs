use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam as sam;
use sam::alignment::{
    io::Write,
    record::{cigar::op::Kind, cigar::Op, data::field::Tag, Flags},
    record_buf::{data::field::Value, Cigar, Data, QualityScores, RecordBuf, Sequence},
};
use sam::header::record::value::{map::ReferenceSequence, Map};
use std::collections::BTreeMap;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const GROUP_A_BARCODE: &str = "AAAAAAAAAAAAAAAA";
const GROUP_B_BARCODE: &str = "CCCCCCCCCCCCCCCC";
const UNUSED_BARCODE: &str = "GGGGGGGGGGGGGGGG";
const UMIS: [&str; 12] = [
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
];

struct Scratch(PathBuf);

impl Scratch {
    fn new(label: &str) -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-end-uniform-{label}-{}-{nonce}",
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
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn alignment_record(
    name: &str,
    start: usize,
    barcode: &str,
    umi: &str,
    flags: Flags,
    nh: u8,
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
        .set_sequence(Sequence::from(vec![b'C'; 50]))
        .set_quality_scores(QualityScores::from(vec![30; 50]))
        .set_data(data)
        .build()
}

fn write_sample_bam(path: &Path, counts: [[usize; 2]; 2]) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(5_000).unwrap()),
        )
        .build();
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    let mut umi = 0usize;
    for (site, start) in [101usize, 2_601].into_iter().enumerate() {
        for (group, barcode) in [GROUP_A_BARCODE, GROUP_B_BARCODE].into_iter().enumerate() {
            for replicate in 0..counts[site][group] {
                let record = alignment_record(
                    &format!("s{site}-g{group}-u{replicate}"),
                    start,
                    barcode,
                    UMIS[umi],
                    Flags::empty(),
                    1,
                );
                writer.write_alignment_record(&header, &record).unwrap();
                umi += 1;
            }
        }
    }
    for (name, start, flags) in [
        ("unused-mm-primary", 4_001, Flags::empty()),
        ("unused-mm-secondary", 4_201, Flags::SECONDARY),
    ] {
        let record = alignment_record(name, start, UNUSED_BARCODE, UMIS[11], flags, 2);
        writer.write_alignment_record(&header, &record).unwrap();
    }
    for (name, start, umi) in [
        ("unused-edge-a1", 4_401, "AACCAACCAACC"),
        ("unused-edge-b", 4_401, "AACCAACCAACA"),
        ("unused-edge-a2", 4_421, "AACCAACCAACC"),
    ] {
        let record = alignment_record(name, start, UNUSED_BARCODE, umi, Flags::empty(), 1);
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

struct Fixture {
    root: Scratch,
    design: PathBuf,
    annotation: PathBuf,
    genome: PathBuf,
    polyasite: PathBuf,
}

impl Fixture {
    fn new(label: &str) -> Self {
        let root = Scratch::new(label);
        let genome = root.0.join("genome.fa");
        let annotation = root.0.join("annotation.gtf");
        let polyasite = root.0.join("polyasite.bed");
        let whitelist = root.0.join("whitelist.txt");
        std::fs::write(&genome, format!(">chr1\n{}\n", "C".repeat(5_000))).unwrap();
        std::fs::write(
            &annotation,
            "chr1\ttest\texon\t1\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\";\n",
        )
        .unwrap();
        std::fs::write(
            &polyasite,
            "chr1\t250\t251\tP1\t.\t+\nchr1\t2750\t2751\tP2\t.\t+\n",
        )
        .unwrap();
        std::fs::write(
            &whitelist,
            format!("{GROUP_A_BARCODE}\n{GROUP_B_BARCODE}\n{UNUSED_BARCODE}\n"),
        )
        .unwrap();

        let mut design = String::from("sample\tcondition\tarchive\tgroups\n");
        for (index, counts) in [
            [[2usize, 1usize], [1usize, 2usize]],
            [[3usize, 1usize], [1usize, 3usize]],
        ]
        .into_iter()
        .enumerate()
        {
            let bam = root.0.join(format!("sample-{index}.bam"));
            let archive = root.0.join(format!("sample-{index}.aie"));
            let groups = root.0.join(format!("sample-{index}.groups.tsv"));
            write_sample_bam(&bam, counts);
            std::fs::write(
                &groups,
                format!("{GROUP_A_BARCODE}\tA\n{GROUP_B_BARCODE}\tB\n"),
            )
            .unwrap();
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
            design.push_str(&format!(
                "S{index}\tcohort\t{}\t{}\n",
                archive.display(),
                groups.display()
            ));
        }
        let design_path = root.0.join("design.tsv");
        std::fs::write(&design_path, design).unwrap();
        Self {
            root,
            design: design_path,
            annotation,
            genome,
            polyasite,
        }
    }

    fn command(&self, kind: &str, out_dir: &Path) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_aie"));
        command
            .arg("cohort")
            .arg(kind)
            .arg("--design")
            .arg(&self.design)
            .arg("--gtf")
            .arg(&self.annotation)
            .arg("--genome")
            .arg(&self.genome)
            .arg("--polyasite")
            .arg(&self.polyasite)
            .arg("--group-contrast")
            .arg("A:B")
            .arg("--out-dir")
            .arg(out_dir)
            .arg("--site-gap")
            .arg("24")
            .arg("--tail-extend")
            .arg("1")
            .arg("--min-site-umis")
            .arg("1")
            .arg("--min-site-samples")
            .arg("1")
            .arg("--min-group-gene-umis")
            .arg("1")
            .arg("--min-samples")
            .arg("2")
            .arg("--min-distal-umis")
            .arg("1");
        if kind == "transcript-ends" {
            command.arg("--motif-min-samples").arg("1");
        }
        command
    }
}

fn table<'a>(document: &'a serde_json::Value, name: &str) -> &'a serde_json::Value {
    document["data"]["tables"]
        .as_array()
        .unwrap()
        .iter()
        .find(|value| value["name"] == name)
        .unwrap_or_else(|| panic!("missing table {name}"))
}

fn assert_tables(document: &serde_json::Value, expected: &[(&str, &str)]) {
    let tables = document["data"]["tables"].as_array().unwrap();
    assert_eq!(tables.len(), expected.len());
    for &(name, schema) in expected {
        let table = table(document, name);
        assert_eq!(table["schema"]["id"], schema);
        assert!(table["schema"]["semantics"]["row_semantics"].is_string());
    }
}

fn compare_artifacts(left: &Path, right: &Path, names: &[&str]) {
    for name in names {
        assert_eq!(
            std::fs::read(left.join(name)).unwrap(),
            std::fs::read(right.join(name)).unwrap(),
            "scientific artifact {name} changed under uniform reporting"
        );
    }
}

fn tsv_site_counts(path: &Path) -> BTreeMap<(u64, String, String), f64> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let count_columns: Vec<(usize, String, String)> = header
        .iter()
        .enumerate()
        .skip_while(|(_, name)| !name.contains(':'))
        .map(|(column, name)| {
            let (sample, group) = name.split_once(':').unwrap();
            (column, sample.to_owned(), group.to_owned())
        })
        .collect();
    let mut counts = BTreeMap::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let site = fields[0].parse::<u64>().unwrap();
        for (column, sample, group) in &count_columns {
            counts.insert(
                (site, sample.clone(), group.clone()),
                fields[*column].parse::<f64>().unwrap(),
            );
        }
    }
    counts
}

fn uniform_site_counts(
    document: &serde_json::Value,
    table_name: &str,
) -> BTreeMap<(u64, String, String), f64> {
    table(document, table_name)["rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            let row = row.as_array().unwrap();
            (
                (
                    row[0].as_u64().unwrap(),
                    row[1].as_str().unwrap().to_owned(),
                    row[2].as_str().unwrap().to_owned(),
                ),
                row[3].as_f64().unwrap(),
            )
        })
        .collect()
}

#[test]
fn uniform_end_reports_preserve_artifacts_and_normalize_the_same_counts() {
    let fixture = Fixture::new("parity");
    let transcript_legacy = fixture.root.0.join("transcript-legacy");
    let transcript_uniform = fixture.root.0.join("transcript-uniform");
    let legacy = run(fixture.command("transcript-ends", &transcript_legacy));
    let legacy_summary: serde_json::Value = serde_json::from_slice(&legacy.stdout).unwrap();
    assert_eq!(
        legacy_summary["schema"],
        "gravlax.cohort.transcript-ends.v1"
    );
    let mut command = fixture.command("transcript-ends", &transcript_uniform);
    command.arg("--report-format").arg("json");
    let uniform = run(command);
    let uniform: serde_json::Value = serde_json::from_slice(&uniform.stdout).unwrap();
    assert_eq!(
        uniform["result_schema"],
        "gravlax.cohort.transcript-ends.result.v1"
    );
    assert_tables(
        &uniform,
        &[
            ("samples", "gravlax.cohort.transcript-ends.samples.v1"),
            ("sites", "gravlax.cohort.transcript-ends.sites.v1"),
            (
                "site_counts",
                "gravlax.cohort.transcript-ends.site-counts.v1",
            ),
            (
                "mixture_sites",
                "gravlax.cohort.transcript-ends.polyasite-mixture.sites.v1",
            ),
            (
                "mixture_site_counts",
                "gravlax.cohort.transcript-ends.polyasite-mixture.site-counts.v1",
            ),
            ("genes", "gravlax.cohort.transcript-ends.genes.v1"),
            (
                "gene_usages",
                "gravlax.cohort.transcript-ends.gene-usages.v1",
            ),
            (
                "fragment_kernel",
                "gravlax.cohort.transcript-ends.fragment-kernel.v1",
            ),
            ("artifacts", "gravlax.cohort.transcript-ends.artifacts.v1"),
        ],
    );
    compare_artifacts(
        &transcript_legacy,
        &transcript_uniform,
        &[
            "sites.tsv",
            "genes.tsv",
            "genes.polyasite.tsv",
            "polyasite-mixture-sites.tsv",
            "polyasite-mixture-genes.tsv",
            "fragment-kernel.tsv",
        ],
    );
    assert_eq!(
        tsv_site_counts(&transcript_uniform.join("sites.tsv")),
        uniform_site_counts(&uniform, "site_counts")
    );
    assert_eq!(
        tsv_site_counts(&transcript_uniform.join("polyasite-mixture-sites.tsv")),
        uniform_site_counts(&uniform, "mixture_site_counts")
    );
    assert_eq!(
        table(&uniform, "sites")["schema"]["semantics"]["row_semantics"],
        "set"
    );
    assert!(table(&uniform, "sites").get("selection").is_none());
    assert!(uniform["provenance"]["parameters"]["design_identity"]
        .as_str()
        .unwrap()
        .starts_with("full-file-blake3-v1:"));
    assert!(uniform["provenance"]["parameters"]["genome_identity"]
        .as_str()
        .unwrap()
        .starts_with("aie-genome-blake3-v1:"));
    assert!(
        uniform["provenance"]["parameters"]["archive_genome_signature"]
            .as_str()
            .unwrap()
            .starts_with("aie-genome-blake3-v1:")
    );

    let mixture_legacy = fixture.root.0.join("mixture-legacy");
    let mixture_uniform = fixture.root.0.join("mixture-uniform");
    run(fixture.command("polyasite-mixture", &mixture_legacy));
    let report = fixture.root.0.join("mixture-report.json");
    let mut command = fixture.command("polyasite-mixture", &mixture_uniform);
    command
        .arg("--report-format")
        .arg("json")
        .arg("--report-output")
        .arg(&report);
    let output = run(command);
    assert!(output.stdout.is_empty());
    let uniform: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&report).unwrap()).unwrap();
    assert_eq!(
        uniform["result_schema"],
        "gravlax.cohort.polyasite-mixture.result.v1"
    );
    assert_tables(
        &uniform,
        &[
            ("samples", "gravlax.cohort.polyasite-mixture.samples.v1"),
            ("sites", "gravlax.cohort.polyasite-mixture.sites.v1"),
            (
                "site_counts",
                "gravlax.cohort.polyasite-mixture.site-counts.v1",
            ),
            ("genes", "gravlax.cohort.polyasite-mixture.genes.v1"),
            (
                "gene_usages",
                "gravlax.cohort.polyasite-mixture.gene-usages.v1",
            ),
            (
                "fragment_kernel",
                "gravlax.cohort.polyasite-mixture.fragment-kernel.v1",
            ),
            (
                "heldout_kernel",
                "gravlax.cohort.polyasite-mixture.heldout-kernel.v1",
            ),
            ("artifacts", "gravlax.cohort.polyasite-mixture.artifacts.v1"),
        ],
    );
    compare_artifacts(
        &mixture_legacy,
        &mixture_uniform,
        &["sites.tsv", "genes.tsv", "fragment-kernel.tsv"],
    );
    assert_eq!(
        tsv_site_counts(&mixture_uniform.join("sites.tsv")),
        uniform_site_counts(&uniform, "site_counts")
    );
    assert_eq!(
        uniform["provenance"]["archives"].as_array().unwrap().len(),
        2
    );
    assert!(uniform["provenance"]["annotation_digest"]
        .as_str()
        .unwrap()
        .starts_with("blake3:"));
    assert_eq!(
        uniform["data"]["summary"]["artifact"]["directory_transactional"],
        false
    );

    let transcript_tsv = fixture.root.0.join("transcript-tsv");
    let mut command = fixture.command("transcript-ends", &transcript_tsv);
    command.arg("--report-format").arg("tsv");
    let output = run(command);
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.starts_with("# envelope_schema="));
    assert!(stdout.contains("# result_schema=gravlax.cohort.transcript-ends.result.v1\n"));
    assert!(stdout.contains("# table=site_counts\n"));

    let mixture_text = fixture.root.0.join("mixture-text");
    let mut command = fixture.command("polyasite-mixture", &mixture_text);
    command.arg("--report-format").arg("text");
    let output = run(command);
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.starts_with("result: gravlax.cohort.polyasite-mixture.result.v1\n"));
    assert!(stdout.contains("table: heldout_kernel\n"));
}

#[test]
fn malformed_caps_and_no_clobber_fail_without_publishing_results() {
    let fixture = Fixture::new("failures");

    let occupied_report = fixture.root.0.join("occupied.json");
    std::fs::write(&occupied_report, "keep me\n").unwrap();
    let occupied_out = fixture.root.0.join("occupied-out");
    let mut occupied = fixture.command("transcript-ends", &occupied_out);
    occupied
        .arg("--report-format")
        .arg("json")
        .arg("--report-output")
        .arg(&occupied_report);
    let output = occupied.output().unwrap();
    assert!(!output.status.success());
    assert_eq!(
        std::fs::read_to_string(&occupied_report).unwrap(),
        "keep me\n"
    );
    assert!(!occupied_out.exists());

    let capped_report = fixture.root.0.join("capped.json");
    let capped_out = fixture.root.0.join("capped-out");
    let mut capped = fixture.command("transcript-ends", &capped_out);
    capped
        .arg("--max-sites")
        .arg("1")
        .arg("--report-format")
        .arg("json")
        .arg("--report-output")
        .arg(&capped_report);
    let output = capped.output().unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("output was not truncated"));
    assert!(!capped_out.exists());
    assert!(!capped_report.exists());

    let malformed_out = fixture.root.0.join("malformed-out");
    let mut malformed = fixture.command("polyasite-mixture", &malformed_out);
    malformed.arg("--report-format").arg("csv");
    let output = malformed.output().unwrap();
    assert!(!output.status.success());
    assert!(!malformed_out.exists());
}
