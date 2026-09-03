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
use std::ffi::OsStr;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const CELL_A: &str = "AAAAAAAAAAAAAAAA";
const CELL_B: &str = "CCCCCCCCCCCCCCCC";

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-transcript-ec-cli-{}-{nonce}",
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

fn run(arguments: &[&OsStr]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(arguments)
        .output()
        .unwrap()
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn unspliced_record(
    name: &str,
    start: usize,
    cell: &str,
    umi: &str,
    flags: Flags,
    nh: u8,
) -> RecordBuf {
    let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
    let data: Data = [
        (Tag::new(b'C', b'R'), Value::from(cell)),
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

struct Fixture {
    _scratch: Scratch,
    archive: PathBuf,
    annotation: PathBuf,
    groups: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let scratch = Scratch::new();
        let bam_path = scratch.0.join("evidence.bam");
        let whitelist = scratch.0.join("whitelist.txt");
        let archive = scratch.0.join("evidence.aie");
        let annotation = scratch.0.join("annotation.gtf");
        let groups = scratch.0.join("groups.tsv");

        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
            )
            .build();
        let file = std::fs::File::create(&bam_path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        // Cell A, UMI 1: one multimapper record has candidates {T1,T2}; a retained unique
        // record in a later archive chunk has {T1}. The global intersection must be {T1}.
        // Cell A, UMI 2 and cell B, UMI 3 have only the multimapper record and remain {T1,T2}.
        // UMI 4 overlaps the requested window but extends beyond every selected exon, so it must
        // remain as an explicit no-compatible outcome. UMI 5 lies outside the selector entirely.
        for record in [
            unspliced_record("umi1-mm", 101, CELL_A, "AAAAAAAAAAAA", Flags::empty(), 2),
            unspliced_record("umi2-mm", 111, CELL_A, "CCCCCCCCCCCC", Flags::empty(), 2),
            unspliced_record("umi3-mm", 121, CELL_B, "GGGGGGGGGGGG", Flags::empty(), 2),
            unspliced_record(
                "umi4-overlap-discordant",
                171,
                CELL_A,
                "TTTTTTTTTTTT",
                Flags::empty(),
                1,
            ),
            unspliced_record("umi1-mm", 301, CELL_A, "AAAAAAAAAAAA", Flags::SECONDARY, 2),
            unspliced_record("umi2-mm", 311, CELL_A, "CCCCCCCCCCCC", Flags::SECONDARY, 2),
            unspliced_record("umi3-mm", 321, CELL_B, "GGGGGGGGGGGG", Flags::SECONDARY, 2),
            unspliced_record(
                "umi1-unique",
                1_200_001,
                CELL_A,
                "AAAAAAAAAAAA",
                Flags::empty(),
                1,
            ),
            unspliced_record(
                "umi5-unrelated",
                1_500_001,
                CELL_B,
                "ACGTACGTACGT",
                Flags::empty(),
                1,
            ),
            // Populate the archive's UMI-adjacency stream without adding selected evidence.
            unspliced_record(
                "edge-a",
                1_600_001,
                CELL_B,
                "AACCAACCAACC",
                Flags::empty(),
                1,
            ),
            // A second representative for this UMI/chain keeps the archive's representative-
            // offset value stream nonempty. Both placements remain outside the query selector.
            unspliced_record(
                "edge-a-extended",
                1_600_031,
                CELL_B,
                "AACCAACCAACC",
                Flags::empty(),
                1,
            ),
            unspliced_record(
                "edge-b",
                1_600_011,
                CELL_B,
                "AACCAACCAACA",
                Flags::empty(),
                1,
            ),
        ] {
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();

        std::fs::write(&whitelist, format!("{CELL_A}\n{CELL_B}\n")).unwrap();
        std::fs::write(
            &annotation,
            concat!(
                "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_id \"E1\"; gene_name \"DUP\";\n",
                "chr1\ttest\texon\t1200001\t1200100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_id \"E2\"; gene_name \"DUP\";\n",
                "chr1\ttest\texon\t301\t400\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T2\"; exon_id \"E3\"; gene_name \"DUP\";\n",
            ),
        )
        .unwrap();
        std::fs::write(&groups, format!("{CELL_A}\tcase\n{CELL_B}\tcontrol\n")).unwrap();

        let output = run(&[
            "ingest-archive".as_ref(),
            bam_path.as_os_str(),
            "--whitelist".as_ref(),
            whitelist.as_os_str(),
            "--out".as_ref(),
            archive.as_os_str(),
            "--zstd-level".as_ref(),
            "1".as_ref(),
            "--chunk-mb".as_ref(),
            "1".as_ref(),
        ]);
        assert_success(&output);

        Self {
            _scratch: scratch,
            archive,
            annotation,
            groups,
        }
    }

    fn query(&self, extra: &[&OsStr]) -> Output {
        let mut command = Command::new(env!("CARGO_BIN_EXE_aie"));
        command
            .arg("query")
            .arg(&self.archive)
            .arg("transcript-ecs")
            .arg("--annotation-file")
            .arg(&self.annotation)
            .args(["--assembly", "test-assembly"])
            .args(["--annotation-label", "test-release"])
            .args(extra);
        command.output().unwrap()
    }
}

fn json(output: &Output) -> serde_json::Value {
    assert_success(output);
    assert!(output.stderr.is_empty(), "machine query wrote diagnostics");
    serde_json::from_slice(&output.stdout).unwrap()
}

fn field_index(table: &serde_json::Value, name: &str) -> usize {
    table["schema"]["fields"]
        .as_array()
        .unwrap()
        .iter()
        .position(|field| field["name"] == name)
        .unwrap()
}

fn catalog_row<'a>(table: &'a serde_json::Value, transcripts: &[&str]) -> &'a serde_json::Value {
    let transcript_column = field_index(table, "transcript_ids");
    let expected = serde_json::json!(transcripts);
    table["rows"]
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row[transcript_column] == expected)
        .unwrap()
}

#[test]
fn command_unions_alternatives_intersects_records_and_uses_typed_envelope() {
    let fixture = Fixture::new();
    let arguments: Vec<&OsStr> = vec![
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--emit-membership".as_ref(),
        "--format".as_ref(),
        "json".as_ref(),
    ];
    let first = fixture.query(&arguments);
    let second = fixture.query(&arguments);
    assert_success(&first);
    assert_success(&second);
    assert_eq!(
        first.stdout, second.stdout,
        "JSON output is not deterministic"
    );

    let value = json(&first);
    assert_eq!(value["$schema"], "gravlax.result-envelope.v1");
    assert_eq!(value["result_schema"], "gravlax.query.transcript-ecs.v1");
    assert_eq!(
        value["data"]["catalog"]["schema"]["id"],
        "gravlax.query.transcript-ecs.catalog.v1"
    );
    assert_eq!(
        value["data"]["counts"]["schema"]["id"],
        "gravlax.query.transcript-ecs.counts.v1"
    );
    assert_eq!(
        value["data"]["membership"]["schema"]["id"],
        "gravlax.query.transcript-ecs.membership.v1"
    );
    assert_eq!(
        value["data"]["scope"]["archive_access"],
        "full_archive_scan"
    );
    assert_eq!(value["data"]["scope"]["chunk_pruning_applied"], false);
    assert!(
        value["data"]["summary"]["archive_umi_classes_scanned"]
            .as_u64()
            .unwrap()
            >= 5
    );
    assert_eq!(value["data"]["summary"]["selector_relevant_umi_classes"], 4);
    assert_eq!(value["data"]["summary"]["scoped_umi_classes"], 4);
    assert_eq!(value["data"]["summary"]["assigned_umi_classes"], 3);
    assert_eq!(value["data"]["summary"]["unassigned_umi_classes"], 1);
    assert_eq!(
        value["data"]["summary"]["no_compatible_transcript_umi_classes"],
        1
    );
    assert_eq!(value["data"]["summary"]["transcript_ecs"], 2);
    assert_eq!(value["data"]["summary"]["membership_rows"], 4);
    assert_eq!(
        value["data"]["semantics"]["compatibility"]["abundance_inferred"],
        false
    );
    assert_eq!(
        value["data"]["semantics"]["compatibility"]["full_isoform_phasing_claimed"],
        false
    );

    let catalog = &value["data"]["catalog"];
    let count_column = field_index(catalog, "archived_umi_class_count");
    let cell_column = field_index(catalog, "cell_count");
    assert_eq!(catalog_row(catalog, &["T1"])[count_column], 1);
    assert_eq!(catalog_row(catalog, &["T1"])[cell_column], 1);
    assert_eq!(catalog_row(catalog, &["T1", "T2"])[count_column], 2);
    assert_eq!(catalog_row(catalog, &["T1", "T2"])[cell_column], 2);

    let membership = &value["data"]["membership"];
    let ec_column = field_index(membership, "ec_id");
    let empty_column = field_index(membership, "no_compatible_transcript");
    let empty_rows: Vec<_> = membership["rows"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|row| row[empty_column] == true)
        .collect();
    assert_eq!(empty_rows.len(), 1);
    assert!(empty_rows[0][ec_column].is_null());

    let root = value["data"]["summary"]["archive_root_blake3"]
        .as_str()
        .unwrap();
    assert_eq!(root.len(), 64);
    assert_eq!(value["provenance"]["assembly"], "test-assembly");
    assert_eq!(value["provenance"]["annotation"], "test-release");
    assert!(value["provenance"]["annotation_digest"]
        .as_str()
        .unwrap()
        .starts_with("blake3:"));
    let archive_identity = value["provenance"]["archives"][0].as_str().unwrap();
    assert_eq!(archive_identity, format!("aie-directory-root-v2:{root}"));
}

#[test]
fn gene_resolution_and_group_aggregation_preserve_per_cell_class_counts() {
    let fixture = Fixture::new();
    let output = fixture.query(&[
        "--feature".as_ref(),
        "gene:G1".as_ref(),
        "--groups".as_ref(),
        fixture.groups.as_os_str(),
        "--format".as_ref(),
        "json".as_ref(),
    ]);
    let value = json(&output);
    assert_eq!(value["data"]["scope"]["selection"]["kind"], "gene");
    assert_eq!(value["data"]["scope"]["selection"]["stable_id"], "G1");
    assert_eq!(
        value["data"]["scope"]["selected_transcript_ids"],
        serde_json::json!(["T1"])
    );
    assert_eq!(value["data"]["summary"]["transcript_ecs"], 1);
    assert_eq!(value["data"]["summary"]["scoped_umi_classes"], 4);

    let counts = &value["data"]["counts"];
    let aggregation_column = field_index(counts, "aggregation");
    let key_column = field_index(counts, "key");
    let count_column = field_index(counts, "archived_umi_class_count");
    let rows = counts["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 3);
    assert!(rows.iter().all(|row| row[aggregation_column] == "group"));
    let case: u64 = rows
        .iter()
        .filter(|row| row[key_column] == "case")
        .map(|row| row[count_column].as_u64().unwrap())
        .sum();
    let control: u64 = rows
        .iter()
        .filter(|row| row[key_column] == "control")
        .map(|row| row[count_column].as_u64().unwrap())
        .sum();
    assert_eq!(case, 3);
    assert_eq!(control, 1);
    let cell_scope = &value["data"]["scope"]["cells"];
    assert_eq!(cell_scope["source"], "groups");
    assert_eq!(
        cell_scope["source_path"],
        fixture.groups.to_string_lossy().as_ref()
    );
    assert!(cell_scope["source_content_blake3"]
        .as_str()
        .unwrap()
        .starts_with("blake3:"));
    assert!(cell_scope["resolved_population_blake3"]
        .as_str()
        .unwrap()
        .starts_with("blake3:"));
    assert_eq!(
        value["provenance"]["parameters"]["cell_scope"],
        cell_scope.clone()
    );

    let swapped_groups = fixture._scratch.0.join("groups-swapped.tsv");
    std::fs::write(
        &swapped_groups,
        format!("{CELL_A}\tcontrol\n{CELL_B}\tcase\n"),
    )
    .unwrap();
    let swapped = json(&fixture.query(&[
        "--feature".as_ref(),
        "gene:G1".as_ref(),
        "--groups".as_ref(),
        swapped_groups.as_os_str(),
        "--format".as_ref(),
        "json".as_ref(),
    ]));
    assert_ne!(
        cell_scope["resolved_population_blake3"],
        swapped["data"]["scope"]["cells"]["resolved_population_blake3"]
    );

    let bulk = json(&fixture.query(&[
        "--feature".as_ref(),
        "gene:G1".as_ref(),
        "--agg".as_ref(),
        "bulk".as_ref(),
        "--format".as_ref(),
        "json".as_ref(),
    ]));
    let bulk_counts = &bulk["data"]["counts"];
    let bulk_aggregation = field_index(bulk_counts, "aggregation");
    let bulk_key = field_index(bulk_counts, "key");
    assert!(bulk_counts["rows"]
        .as_array()
        .unwrap()
        .iter()
        .all(|row| row[bulk_aggregation] == "bulk" && row[bulk_key] == "all"));

    let transcript = fixture.query(&[
        "--feature".as_ref(),
        "transcript:T1".as_ref(),
        "--format".as_ref(),
        "json".as_ref(),
    ]);
    assert!(!transcript.status.success());
    assert!(String::from_utf8_lossy(&transcript.stderr).contains("must resolve to a gene"));
    assert!(transcript.stdout.is_empty());

    let ambiguous = fixture.query(&[
        "--feature".as_ref(),
        "DUP".as_ref(),
        "--format".as_ref(),
        "json".as_ref(),
    ]);
    assert!(!ambiguous.status.success());
    assert!(String::from_utf8_lossy(&ambiguous.stderr).contains("ambiguous"));
    assert!(ambiguous.stdout.is_empty());
}

#[test]
fn tsv_is_typed_and_caps_or_no_clobber_never_leave_partial_output() {
    let fixture = Fixture::new();
    let tsv = fixture.query(&[
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--format".as_ref(),
        "tsv".as_ref(),
        "--table".as_ref(),
        "catalog".as_ref(),
    ]);
    assert_success(&tsv);
    assert!(tsv.stderr.is_empty());
    let text = String::from_utf8(tsv.stdout).unwrap();
    assert!(text.starts_with("# envelope_schema=gravlax.result-envelope.v1\n"));
    assert!(text.contains("# result_schema=gravlax.query.transcript-ecs.catalog.v1\n"));
    assert!(text.contains("ec_id\ttranscript_ids\tgene_ids\tambiguous"));
    assert!(!text.contains("transcript EC query ("));

    let capped = fixture._scratch.0.join("capped.json");
    let overflow = fixture.query(&[
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--format".as_ref(),
        "json".as_ref(),
        "--max-ecs".as_ref(),
        "1".as_ref(),
        "-o".as_ref(),
        capped.as_os_str(),
    ]);
    assert!(!overflow.status.success());
    assert!(overflow.stdout.is_empty());
    assert!(String::from_utf8_lossy(&overflow.stderr).contains("--max-ecs"));
    assert!(!capped.exists());

    let membership_capped = fixture._scratch.0.join("membership-capped.tsv");
    let overflow = fixture.query(&[
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--emit-membership".as_ref(),
        "--max-memberships".as_ref(),
        "2".as_ref(),
        "--format".as_ref(),
        "tsv".as_ref(),
        "--table".as_ref(),
        "membership".as_ref(),
        "-o".as_ref(),
        membership_capped.as_os_str(),
    ]);
    assert!(!overflow.status.success());
    assert!(overflow.stdout.is_empty());
    assert!(String::from_utf8_lossy(&overflow.stderr).contains("--max-memberships"));
    assert!(!membership_capped.exists());

    let existing = fixture._scratch.0.join("existing.json");
    std::fs::write(&existing, b"keep-me").unwrap();
    let refused = fixture.query(&[
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--format".as_ref(),
        "json".as_ref(),
        "-o".as_ref(),
        existing.as_os_str(),
    ]);
    assert!(!refused.status.success());
    assert!(refused.stdout.is_empty());
    assert_eq!(std::fs::read(&existing).unwrap(), b"keep-me");
    assert!(std::fs::read_dir(&fixture._scratch.0)
        .unwrap()
        .all(|entry| !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .contains("transcript-ecs-tmp")));
}

fn read_u32(bytes: &[u8], at: &mut usize) -> u32 {
    let end = *at + 4;
    let value = u32::from_le_bytes(bytes[*at..end].try_into().unwrap());
    *at = end;
    value
}

fn skip_string(bytes: &[u8], at: &mut usize) {
    let length = read_u32(bytes, at) as usize;
    *at += length;
    assert!(*at <= bytes.len());
}

fn write_aic_v1(v2: &Path, v1: &Path) {
    const HEADER: usize = 8 + 4 + 8 + 32;
    let bytes = std::fs::read(v2).unwrap();
    assert_eq!(&bytes[..8], b"GRVLXAIC");
    let payload = &bytes[HEADER..];
    let mut at = 0usize;
    let genes = read_u32(payload, &mut at) as usize;
    for _ in 0..genes * 2 {
        skip_string(payload, &mut at);
    }
    let chroms = read_u32(payload, &mut at) as usize;
    for _ in 0..chroms {
        skip_string(payload, &mut at);
    }
    let transcripts = read_u32(payload, &mut at) as usize;
    for _ in 0..transcripts {
        at += 4 + 4 + 1;
        let exons = read_u32(payload, &mut at) as usize;
        at += exons * 8;
    }
    let indexed_chroms = read_u32(payload, &mut at) as usize;
    for _ in 0..indexed_chroms {
        at += 4;
        let entries = read_u32(payload, &mut at) as usize;
        at += entries * 16;
    }
    let legacy_payload = &payload[..at];
    let mut legacy = Vec::new();
    legacy.extend_from_slice(b"GRVLXAIC");
    legacy.extend_from_slice(&1u32.to_le_bytes());
    legacy.extend_from_slice(&(legacy_payload.len() as u64).to_le_bytes());
    legacy.extend_from_slice(blake3::hash(legacy_payload).as_bytes());
    legacy.extend_from_slice(legacy_payload);
    std::fs::write(v1, legacy).unwrap();
}

#[test]
fn legacy_aic_and_invalid_output_contracts_fail_before_writing() {
    let fixture = Fixture::new();
    let compiled = fixture._scratch.0.join("annotation.aic");
    let compile = run(&[
        "compile-annotation".as_ref(),
        fixture.annotation.as_os_str(),
        "--out".as_ref(),
        compiled.as_os_str(),
    ]);
    assert_success(&compile);
    let legacy = fixture._scratch.0.join("annotation-v1.aic");
    write_aic_v1(&compiled, &legacy);
    let destination = fixture._scratch.0.join("legacy.json");
    let output = run(&[
        "query".as_ref(),
        fixture.archive.as_os_str(),
        "transcript-ecs".as_ref(),
        "--annotation-file".as_ref(),
        legacy.as_os_str(),
        "--assembly".as_ref(),
        "test-assembly".as_ref(),
        "--annotation-label".as_ref(),
        "legacy-release".as_ref(),
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--format".as_ref(),
        "json".as_ref(),
        "-o".as_ref(),
        destination.as_os_str(),
    ]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("AIC v1 omitted transcript IDs"));
    assert!(!destination.exists());

    let digest_mismatch = fixture._scratch.0.join("digest-mismatch.json");
    let output = fixture.query(&[
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--annotation-digest".as_ref(),
        "blake3:0000000000000000000000000000000000000000000000000000000000000000".as_ref(),
        "--format".as_ref(),
        "json".as_ref(),
        "-o".as_ref(),
        digest_mismatch.as_os_str(),
    ]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("digest"));
    assert!(!digest_mismatch.exists());

    let missing_table = fixture.query(&[
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--format".as_ref(),
        "tsv".as_ref(),
    ]);
    assert!(!missing_table.status.success());
    assert!(String::from_utf8_lossy(&missing_table.stderr).contains("requires --table"));
    assert!(missing_table.stdout.is_empty());

    let invalid_membership = fixture.query(&[
        "--locus".as_ref(),
        "chr1:0-1300000".as_ref(),
        "--format".as_ref(),
        "tsv".as_ref(),
        "--table".as_ref(),
        "membership".as_ref(),
    ]);
    assert!(!invalid_membership.status.success());
    assert!(
        String::from_utf8_lossy(&invalid_membership.stderr).contains("requires --emit-membership")
    );
    assert!(invalid_membership.stdout.is_empty());
}
