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

const TEST_BARCODE: &str = "AAAAAAAAAAAAAAAA";

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-project-plan-{}-{nonce}",
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

fn aie() -> Command {
    Command::new(env!("CARGO_BIN_EXE_aie"))
}

fn write_empty_archive(path: &Path) {
    evidence_io::format::SectionWriter::create(path, 1)
        .unwrap()
        .finish()
        .unwrap();
}

fn write_empty_collection(path: &Path) {
    let mut prefix = Vec::new();
    prefix.extend_from_slice(b"GRVLXCOL");
    prefix.extend_from_slice(&4u32.to_le_bytes());
    prefix.extend_from_slice(&1u32.to_le_bytes());
    let mut directory = Vec::new();
    directory.push(8);
    directory.extend_from_slice(b"manifest");
    directory.extend_from_slice(&0u64.to_le_bytes());
    directory.extend_from_slice(&0u64.to_le_bytes());
    directory.extend_from_slice(blake3::hash(b"").as_bytes());
    let mut hasher = blake3::Hasher::new();
    hasher.update(&prefix);
    hasher.update(&directory);
    let mut file = prefix;
    file.extend_from_slice(hasher.finalize().as_bytes());
    file.extend_from_slice(&directory);
    std::fs::write(path, file).unwrap();
}

fn write_legacy_aic_without_transcript_ids(path: &Path) {
    fn push_string(bytes: &mut Vec<u8>, value: &str) {
        bytes.extend_from_slice(&(value.len() as u32).to_le_bytes());
        bytes.extend_from_slice(value.as_bytes());
    }

    let mut payload = Vec::new();
    payload.extend_from_slice(&1u32.to_le_bytes());
    push_string(&mut payload, "g1");
    push_string(&mut payload, "G1");
    payload.extend_from_slice(&1u32.to_le_bytes());
    push_string(&mut payload, "chr1");
    payload.extend_from_slice(&1u32.to_le_bytes());
    payload.extend_from_slice(&0u32.to_le_bytes());
    payload.extend_from_slice(&0u32.to_le_bytes());
    payload.push(0);
    payload.extend_from_slice(&1u32.to_le_bytes());
    payload.extend_from_slice(&100u32.to_le_bytes());
    payload.extend_from_slice(&150u32.to_le_bytes());
    payload.extend_from_slice(&1u32.to_le_bytes());
    payload.extend_from_slice(&0u32.to_le_bytes());
    payload.extend_from_slice(&1u32.to_le_bytes());
    payload.extend_from_slice(&100u32.to_le_bytes());
    payload.extend_from_slice(&150u32.to_le_bytes());
    payload.extend_from_slice(&150u32.to_le_bytes());
    payload.extend_from_slice(&0u32.to_le_bytes());

    let mut file = Vec::new();
    file.extend_from_slice(b"GRVLXAIC");
    file.extend_from_slice(&1u32.to_le_bytes());
    file.extend_from_slice(&(payload.len() as u64).to_le_bytes());
    file.extend_from_slice(blake3::hash(&payload).as_bytes());
    file.extend_from_slice(&payload);
    std::fs::write(path, file).unwrap();
}

fn write_replay_fixture_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(2_000_000).unwrap()),
        )
        .build();
    let record = |name: &str, start: usize, umi: &str, flags: Flags, nh: u8| {
        let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
        let data: Data = [
            (Tag::new(b'C', b'R'), Value::from(TEST_BARCODE)),
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
        record("u1a", 101, "AAAAAAAAAAAA", Flags::empty(), 1),
        record("u1b", 121, "AAAAAAAAAAAA", Flags::empty(), 1),
        record("u1c", 101, "AAAAAAAAAAAC", Flags::empty(), 1),
        record("u2", 151, "CCCCCCCCCCCC", Flags::empty(), 1),
        record("mm1", 201, "GGGGGGGGGGGG", Flags::empty(), 3),
        record("mm1", 201, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        record("mm1", 401, "GGGGGGGGGGGG", Flags::SECONDARY, 3),
        record("u1far", 1_100_101, "AAAAAAAAAAAA", Flags::empty(), 1),
    ];
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for record in &records {
        writer.write_alignment_record(&header, record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn success(output: Output) -> Output {
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn init(root: &Path) {
    success(
        aie()
            .args(["project", "init"])
            .arg(root)
            .args(["--name", "test project"])
            .output()
            .unwrap(),
    );
}

fn add(root: &Path, name: &str, path: &Path, kind: &str) {
    success(
        aie()
            .args(["project", "add", name])
            .arg(path)
            .args(["--kind", kind, "--project"])
            .arg(root)
            .output()
            .unwrap(),
    );
}

fn add_external(root: &Path, name: &str, path: &Path, kind: &str) {
    success(
        aie()
            .args(["project", "add", name])
            .arg(path)
            .args(["--kind", kind, "--external", "--project"])
            .arg(root)
            .output()
            .unwrap(),
    );
}

fn add_annotation(root: &Path, name: &str, path: &Path, assembly: &str, annotation: &str) {
    success(
        aie()
            .args(["project", "add", name])
            .arg(path)
            .args([
                "--kind",
                "annotation",
                "--assembly",
                assembly,
                "--annotation-label",
                annotation,
                "--project",
            ])
            .arg(root)
            .output()
            .unwrap(),
    );
}

fn add_with_assembly(root: &Path, name: &str, path: &Path, kind: &str, assembly: &str) {
    success(
        aie()
            .args(["project", "add", name])
            .arg(path)
            .args(["--kind", kind, "--assembly", assembly, "--project"])
            .arg(root)
            .output()
            .unwrap(),
    );
}

fn regular_files_below(root: &Path) -> Vec<PathBuf> {
    let mut pending = vec![root.to_path_buf()];
    let mut files = Vec::new();
    while let Some(directory) = pending.pop() {
        for entry in std::fs::read_dir(directory).unwrap() {
            let entry = entry.unwrap();
            let file_type = entry.file_type().unwrap();
            if file_type.is_dir() {
                pending.push(entry.path());
            } else if file_type.is_file() {
                files.push(entry.path());
            }
        }
    }
    files.sort();
    files
}

#[test]
fn internal_resources_move_with_project_and_external_resources_stay_canonical() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let internal = root.join("data/annotation.gtf");
    std::fs::create_dir_all(internal.parent().unwrap()).unwrap();
    std::fs::write(&internal, b"# internal\n").unwrap();
    add(&root, "annotation", &internal, "annotation");

    let external = scratch.0.join("shared/large-external.aie");
    std::fs::create_dir_all(external.parent().unwrap()).unwrap();
    write_empty_archive(&external);
    let rejected = aie()
        .args(["project", "add", "external-archive"])
        .arg(&external)
        .args(["--kind", "archive", "--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("pass --external"));
    add_external(&root, "external-archive", &external, "archive");

    let plan = root.join("plans/external.yaml");
    std::fs::write(
        &plan,
        "schema_version: 1\nsteps:\n  - id: inspect\n    kind: inspect-archive\n    archive: external-archive\n    format: json\n",
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(resolved["resources"]["external-archive"]["external"], true);
    assert_eq!(
        resolved["resources"]["external-archive"]["identity"]["scheme"],
        "aie-directory-root-v2"
    );

    let moved = scratch.0.join("moved-workspace");
    std::fs::rename(&root, &moved).unwrap();
    let shown = success(
        aie()
            .args(["project", "show", "--project"])
            .arg(&moved)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let shown: serde_json::Value = serde_json::from_slice(&shown.stdout).unwrap();
    let resources = shown["resources"].as_array().unwrap();
    assert!(resources.iter().all(|resource| resource["status"] == "ok"));
    let external_row = resources
        .iter()
        .find(|resource| resource["name"] == "external-archive")
        .unwrap();
    assert_eq!(external_row["external"], true);
    assert_eq!(
        external_row["path"],
        std::fs::canonicalize(&external)
            .unwrap()
            .to_string_lossy()
            .as_ref()
    );
}

#[test]
fn concurrent_project_adds_preserve_both_manifest_updates() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let first = root.join("first.tsv");
    let second = root.join("second.tsv");
    std::fs::write(&first, b"first\n").unwrap();
    std::fs::write(&second, b"second\n").unwrap();

    let mut first_add = aie()
        .args(["project", "add", "first"])
        .arg(&first)
        .args(["--kind", "metadata", "--project"])
        .arg(&root)
        .spawn()
        .unwrap();
    let mut second_add = aie()
        .args(["project", "add", "second"])
        .arg(&second)
        .args(["--kind", "metadata", "--project"])
        .arg(&root)
        .spawn()
        .unwrap();
    assert!(first_add.wait().unwrap().success());
    assert!(second_add.wait().unwrap().success());

    let shown = success(
        aie()
            .args(["project", "show", "--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let shown: serde_json::Value = serde_json::from_slice(&shown.stdout).unwrap();
    let names = shown["resources"]
        .as_array()
        .unwrap()
        .iter()
        .map(|resource| resource["name"].as_str().unwrap())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(names, ["first", "second"].into_iter().collect());
}

#[test]
fn project_and_plan_check_resolve_typed_resources() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);

    let archive = root.join("data/sample.aie");
    let groups = root.join("metadata/groups.tsv");
    let collection = root.join("data/atlas.aicollection");
    std::fs::create_dir_all(archive.parent().unwrap()).unwrap();
    std::fs::create_dir_all(groups.parent().unwrap()).unwrap();
    write_empty_archive(&archive);
    write_empty_collection(&collection);
    std::fs::write(&groups, b"cell-a\tcase\n").unwrap();
    add(&root, "sample", &archive, "archive");
    add(&root, "atlas", &collection, "collection");
    add(&root, "groups", &groups, "groups");

    let plan = root.join("plans/queries.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
name: grouped queries
steps:
  - id: event-table
    kind: query-events
    archive: sample
    locus: chr1:100-500
    event_types: [cassette, alt-donor]
    groups: ignored
    min_support: 2
    uniform_output:
      format: json
      output: results/events.json
    scope:
      groups: groups
      aggregation: group
  - id: atlas-junctions
    kind: collection-jset
    collection: atlas
    include: ["chr1:150-200"]
    exclude: ["chr1:150-300"]
    uniform_output:
      format: json
      output: results/atlas-jset.json
"#,
    )
    .unwrap();

    // Unknown fields fail closed instead of being silently ignored.
    let rejected = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("unknown field `groups`"));

    let text = std::fs::read_to_string(&plan)
        .unwrap()
        .replace("    groups: ignored\n", "");
    std::fs::write(&plan, text).unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--explain")
            .output()
            .unwrap(),
    );
    let stdout = String::from_utf8_lossy(&checked.stdout);
    assert!(stdout.contains("groups `groups` ->"));
    assert!(stdout.contains("collection `atlas` ->"));
    assert!(stdout.contains("valid: 2 supported step(s)"));
    for schema in [
        "gravlax.query.events.result.v1",
        "gravlax.query.events.events.v1",
        "gravlax.query.events.components.v1",
        "gravlax.query.events.counts.v1",
        "gravlax.collection.jset.result.v1",
        "gravlax.collection.jset.requests.v1",
        "gravlax.collection.jset.samples.v1",
        "gravlax.collection.jset.cells.v1",
    ] {
        assert!(stdout.contains(schema), "missing schema {schema}: {stdout}");
    }

    // Discovery is relative to the plan as well as the shell's current directory.
    success(
        aie()
            .current_dir(&root)
            .args(["plan", "check", "plans/queries.yaml"])
            .output()
            .unwrap(),
    );

    let shown = success(
        aie()
            .args(["project", "show", "--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let shown: serde_json::Value = serde_json::from_slice(&shown.stdout).unwrap();
    assert_eq!(shown["name"], "test project");
    assert_eq!(shown["resources"].as_array().unwrap().len(), 3);
}

#[test]
fn annotation_identity_and_feature_resolution_are_explicit_and_reproducible() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let data = root.join("data");
    std::fs::create_dir_all(&data).unwrap();
    let archive = data.join("sample.aie");
    let gtf = data.join("genes.gtf");
    write_empty_archive(&archive);
    std::fs::write(
        &gtf,
        concat!(
            "chr17\ttest\tgene\t101\t200\t.\t+\t.\tgene_id \"ENSG000001.1\"; gene_name \"TP53\";\n",
            "chr17\ttest\ttranscript\t101\t200\t.\t+\t.\tgene_id \"ENSG000001.1\"; gene_name \"TP53\"; transcript_id \"ENST000001.2\";\n",
            "chr17\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"ENSG000001.1\"; gene_name \"TP53\"; transcript_id \"ENST000001.2\"; exon_id \"ENSE000001.3\";\n",
        ),
    )
    .unwrap();
    add_with_assembly(&root, "sample", &archive, "archive", "GRCh38.p14");
    add_annotation(&root, "genes", &gtf, "GRCh38.p14", "GENCODE 49");

    let shown = success(
        aie()
            .args(["project", "show", "--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let shown: serde_json::Value = serde_json::from_slice(&shown.stdout).unwrap();
    let annotation = shown["resources"]
        .as_array()
        .unwrap()
        .iter()
        .find(|resource| resource["name"] == "genes")
        .unwrap();
    assert_eq!(annotation["annotation_identity"]["assembly"], "GRCh38.p14");
    assert_eq!(
        annotation["annotation_identity"]["annotation"],
        "GENCODE 49"
    );

    let plan = root.join("plans/by-feature.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: tp53
    kind: query-region
    archive: sample
    feature:
      identifier: TP53
      assembly: GRCh38.p14
      annotation: GENCODE 49
    annotation: genes
    format: json
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(resolved["schema_version"], 6);
    assert_eq!(
        resolved["steps"][0]["biological_intent"]["stable_id"],
        "ENSG000001.1"
    );
    assert_eq!(
        resolved["steps"][0]["biological_intent"]["locus"],
        "chr17:100-200"
    );
    assert_eq!(
        resolved["steps"][0]["biological_intent"]["strand"],
        "forward"
    );
    assert!(
        resolved["steps"][0]["biological_intent"]["annotation_digest"]
            .as_str()
            .unwrap()
            .starts_with("blake3:")
    );
    assert_eq!(
        resolved["steps"][0]["output_schema_ids"][0],
        "gravlax.query.region.v1"
    );
    assert_eq!(
        resolved["steps"][0]["biological_intent"]["compatibility"][0]["status"],
        "verified"
    );
    assert_eq!(
        resolved["steps"][0]["io_estimate"]["read_bytes_upper_bound"],
        serde_json::Value::Null
    );
    assert_eq!(
        resolved["steps"][0]["io_estimate"]["bound"],
        "known-inputs-only"
    );
    assert!(resolved["steps"][0]["io_estimate"]["note"]
        .as_str()
        .unwrap()
        .contains("no execution read upper bound"));
    assert!(resolved["steps"][0]["args"]
        .as_array()
        .unwrap()
        .iter()
        .any(|argument| argument == "chr17:100-200"));

    let groups = data.join("groups.tsv");
    let genome = data.join("genome.fa");
    std::fs::write(&groups, b"cell-a\tcase\n").unwrap();
    std::fs::write(&genome, b">chr17\nACGT\n").unwrap();
    add(&root, "groups", &groups, "groups");
    add_with_assembly(&root, "genome", &genome, "genome", "GRCh38.p14");
    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: apa
    kind: query-apa
    archive: sample
    feature: TP53
    annotation: genes
    site_gap: 24
    strand: forward
    tsv: true
    groups: groups
    genome: genome
    drop_ip: true
    permute: 100
    seed: 7
    plot: results/tp53-apa.svg
    output: results/tp53-apa.tsv
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(
        resolved["steps"][0]["output_schema_ids"][0],
        "gravlax.query.apa.tsv.v1"
    );
    assert_eq!(resolved["steps"][0]["outputs"].as_array().unwrap().len(), 2);
    assert_eq!(
        resolved["steps"][0]["biological_intent"]["compatibility"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert!(resolved["steps"][0]["args"]
        .as_array()
        .unwrap()
        .iter()
        .any(|argument| argument == "+"));

    let aic = data.join("genes.aic");
    anno::Annotation::from_gtf(&gtf)
        .unwrap()
        .write_compiled(&aic)
        .unwrap();
    add_annotation(&root, "genes-aic", &aic, "GRCh38.p14", "GENCODE 49");
    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: transcript
    kind: query-region
    archive: sample
    feature: transcript:ENST000001
    annotation: genes-aic
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(
        resolved["steps"][0]["biological_intent"]["resolved_kind"],
        "transcript"
    );
    assert_eq!(
        resolved["steps"][0]["biological_intent"]["stable_id"],
        "ENST000001.2"
    );
}

#[test]
fn feature_resolution_fails_closed_on_bad_intent_or_missing_identity() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let archive = root.join("sample.aie");
    let gtf = root.join("ambiguous.gtf");
    write_empty_archive(&archive);
    std::fs::write(
        &gtf,
        concat!(
            "chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id \"g1\"; gene_name \"DUP\"; transcript_id \"t1\"; exon_id \"e1\";\n",
            "chr2\ttest\texon\t201\t250\t.\t-\t.\tgene_id \"g2\"; gene_name \"DUP\"; transcript_id \"t2\"; exon_id \"e2\";\n",
        ),
    )
    .unwrap();
    add(&root, "sample", &archive, "archive");
    add(&root, "unidentified", &gtf, "annotation");
    add_with_assembly(&root, "wrong-build", &archive, "archive", "GRCh37");
    add_annotation(&root, "genes", &gtf, "GRCh38.p14", "release-1");
    let plan = root.join("plans/invalid-feature.yaml");

    let cases = [
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-region
    archive: sample
    feature: g1
    annotation: unidentified
"#,
            "has no scientific identity",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-region
    archive: sample
    feature:
      identifier: g1
      assembly: GRCh37
    annotation: genes
"#,
            "expects assembly `GRCh37`",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-region
    archive: wrong-build
    feature: g1
    annotation: genes
"#,
            "assembly mismatch",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-region
    archive: sample
    feature: DUP
    annotation: genes
"#,
            "is ambiguous",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-region
    archive: sample
    locus: chr1:1-10
    feature: g1
    annotation: genes
"#,
            "exactly one of `locus` or `feature`",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-region
    archive: sample
"#,
            "exactly one of `locus` or `feature`",
        ),
    ];
    for (source, expected) in cases {
        std::fs::write(&plan, source).unwrap();
        let rejected = aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .output()
            .unwrap();
        assert!(
            !rejected.status.success(),
            "case unexpectedly succeeded: {source}"
        );
        assert!(
            String::from_utf8_lossy(&rejected.stderr).contains(expected),
            "stderr did not contain {expected:?}: {}",
            String::from_utf8_lossy(&rejected.stderr)
        );
    }

    let invalid_flags = aie()
        .args(["project", "add", "wrong-kind"])
        .arg(&archive)
        .args([
            "--kind",
            "archive",
            "--assembly",
            "GRCh38.p14",
            "--annotation-label",
            "release-1",
            "--project",
        ])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!invalid_flags.status.success());
    assert!(String::from_utf8_lossy(&invalid_flags.stderr)
        .contains("--annotation-label may only be used for annotation resources"));
}

#[test]
fn feature_resolution_can_consume_a_prior_compiled_annotation() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let archive = root.join("sample.aie");
    let gtf = root.join("genes.gtf");
    write_empty_archive(&archive);
    std::fs::write(
        &gtf,
        "chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id \"g1\"; gene_name \"G1\"; transcript_id \"t1\"; exon_id \"e1\";\n",
    )
    .unwrap();
    add(&root, "sample", &archive, "archive");
    add_annotation(&root, "genes", &gtf, "GRCh38.p14", "local-1");
    let plan = root.join("plans/compile-query.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: compile
    kind: compile-annotation
    annotation: genes
    output: results/genes.aic
  - id: query
    kind: query-region
    archive: sample
    feature: gene:g1
    annotation: step:compile:annotation
    format: json
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(resolved["steps"][1]["biological_intent"]["stable_id"], "g1");
    assert_eq!(
        resolved["steps"][1]["biological_intent"]["compatibility"][0]["status"],
        "unverified"
    );
    assert!(resolved["steps"][1]["biological_intent"]["annotation_path"]
        .as_str()
        .unwrap()
        .ends_with("results/genes.aic"));
    assert_eq!(
        resolved["steps"][1]["step_inputs"][0]["producer_step"],
        "compile"
    );
    assert_eq!(
        resolved["steps"][1]["io_estimate"]["read_bytes_upper_bound"],
        serde_json::Value::Null
    );
    assert_eq!(
        resolved["steps"][1]["io_estimate"]["bound"],
        "known-inputs-only"
    );
}

#[test]
fn annotation_comparison_and_transcript_ec_steps_disclose_scientific_contracts() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let data = root.join("data");
    std::fs::create_dir_all(&data).unwrap();
    let archive = data.join("sample.aie");
    let before = data.join("before.gtf");
    let after = data.join("after.gtf");
    write_empty_archive(&archive);
    std::fs::write(
        &before,
        "chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id \"g1.1\"; gene_name \"G1\"; transcript_id \"t1.1\"; exon_id \"e1\";\n",
    )
    .unwrap();
    std::fs::write(
        &after,
        "chr1\ttest\texon\t101\t175\t.\t+\t.\tgene_id \"g1.2\"; gene_name \"G1\"; transcript_id \"t1.2\"; exon_id \"e1b\";\n",
    )
    .unwrap();
    add_with_assembly(&root, "sample", &archive, "archive", "GRCh38");
    add_annotation(&root, "before", &before, "GRCh38", "release-before");
    add_annotation(&root, "after", &after, "GRCh38", "release-after");

    let plan = root.join("plans/scientific.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
name: scientific capabilities
steps:
  - id: compare
    kind: compare-annotations
    archive: sample
    annotation_a: before
    annotation_b: after
    gene_key: unversioned
    solo_strand: reverse
    max_molecule_witnesses: 0
    max_row_transitions_per_molecule: 0
    format: json
    output: results/comparison.json
  - id: ecs
    kind: query-transcript-ecs
    archive: sample
    feature:
      identifier: gene:g1
      assembly: GRCh38
      annotation: release-before
    annotation: before
    solo_strand: unstranded
    emit_membership: true
    max_ecs: 25
    max_memberships: 250
    format: json
    output: results/ecs.json
"#,
    )
    .unwrap();

    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(resolved["schema_version"], 6);
    let comparison = &resolved["steps"][0];
    assert_eq!(comparison["kind"], "compare-annotations");
    assert_eq!(comparison["annotation_inputs"][0]["role"], "a");
    assert_eq!(comparison["annotation_inputs"][1]["role"], "b");
    assert_eq!(
        comparison["annotation_inputs"][0]["compatibility"][0]["status"],
        "verified"
    );
    assert_eq!(comparison["annotation_comparison"]["assembly"], "GRCh38");
    assert!(
        comparison["annotation_comparison"]["transition_evidence_semantics"]
            .as_str()
            .unwrap()
            .contains("not additive")
    );
    assert_eq!(
        comparison["output_schema_ids"],
        serde_json::json!([
            "gravlax.annotation.compare.v1",
            "gravlax.annotation.compare.count-deltas.v1",
            "gravlax.annotation.compare.class-transitions.v1",
            "gravlax.annotation.compare.contributing-causes.v1",
            "gravlax.annotation.compare.witnesses.v1"
        ])
    );
    assert_eq!(
        comparison["io_estimate"]["read_bytes_upper_bound"],
        serde_json::Value::Null
    );
    let comparison_args = comparison["args"].as_array().unwrap();
    assert!(comparison_args
        .iter()
        .any(|argument| argument == "--annotation-a-digest"));
    assert!(comparison_args
        .iter()
        .any(|argument| argument == "--max-molecule-witnesses"));
    assert!(comparison_args.iter().any(|argument| argument == "0"));

    let ecs = &resolved["steps"][1];
    assert_eq!(ecs["kind"], "query-transcript-ecs");
    assert_eq!(ecs["biological_intent"]["resolved_kind"], "gene");
    assert_eq!(ecs["biological_intent"]["stable_id"], "g1.1");
    assert_eq!(ecs["biological_intent"]["locus"], "chr1:100-150");
    assert_eq!(ecs["annotation_inputs"][0]["role"], "query");
    assert_eq!(
        ecs["output_schema_ids"],
        serde_json::json!([
            "gravlax.query.transcript-ecs.v1",
            "gravlax.query.transcript-ecs.catalog.v1",
            "gravlax.query.transcript-ecs.counts.v1",
            "gravlax.query.transcript-ecs.membership.v1"
        ])
    );
    let ecs_args = ecs["args"].as_array().unwrap();
    let feature_index = ecs_args
        .iter()
        .position(|argument| argument == "--feature")
        .unwrap();
    assert_eq!(ecs_args[feature_index + 1], "gene:g1.1");
    assert!(ecs["explanation"]
        .as_array()
        .unwrap()
        .iter()
        .any(|line| line.as_str().unwrap().contains("not abundance estimates")));
}

#[test]
fn scientific_capability_steps_accept_prior_annotation_and_archive_outputs() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let data = root.join("data");
    std::fs::create_dir_all(&data).unwrap();
    let gtf = data.join("genes.gtf");
    let other = data.join("other.gtf");
    let bam = data.join("reads.bam");
    let whitelist = data.join("whitelist.txt");
    std::fs::write(
        &gtf,
        "chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id \"g1\"; gene_name \"G1\"; transcript_id \"t1\"; exon_id \"e1\";\n",
    )
    .unwrap();
    std::fs::write(
        &other,
        "chr1\ttest\texon\t101\t175\t.\t+\t.\tgene_id \"g1\"; gene_name \"G1\"; transcript_id \"t2\"; exon_id \"e2\";\n",
    )
    .unwrap();
    std::fs::write(&bam, b"plan check does not decode BAM\n").unwrap();
    std::fs::write(&whitelist, b"AAAAAAAAAAAAAAAA\n").unwrap();
    add_annotation(&root, "genes", &gtf, "GRCh38", "source");
    add_annotation(&root, "other", &other, "GRCh38", "other");
    add(&root, "reads", &bam, "bam");
    add(&root, "whitelist", &whitelist, "whitelist");

    let plan = root.join("plans/prior-science.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: ingest
    kind: ingest-archive
    bam: reads
    whitelist: whitelist
    output: results/sample.aie
  - id: compile
    kind: compile-annotation
    annotation: genes
    output: results/genes.aic
  - id: compare
    kind: compare-annotations
    archive: step:ingest
    annotation_a: step:compile:annotation
    annotation_b: other
    format: json
  - id: ecs
    kind: query-transcript-ecs
    archive: step:ingest
    locus: chr1:100-175
    annotation: step:compile:annotation
    format: json
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    for index in [2usize, 3usize] {
        let step = &resolved["steps"][index];
        assert!(step["step_inputs"].as_array().unwrap().iter().any(|input| {
            input["producer_step"] == "ingest" && input["resource_kind"] == "archive"
        }));
        assert_eq!(
            step["io_estimate"]["read_bytes_upper_bound"],
            serde_json::Value::Null
        );
    }
    let prior_annotation = &resolved["steps"][2]["annotation_inputs"][0];
    assert_eq!(prior_annotation["resource"], "step:compile:annotation");
    assert!(prior_annotation.get("expected_command_identity").is_none());
    assert!(resolved["steps"][2]["step_inputs"]
        .as_array()
        .unwrap()
        .iter()
        .any(|input| input["producer_step"] == "compile"));
    assert!(resolved["steps"][3]["annotation_inputs"][0]
        .get("expected_command_identity")
        .is_none());
    assert!(resolved["steps"][3]["biological_intent"].is_null());
    assert_eq!(
        resolved["steps"][3]["annotation_inputs"][0]["compatibility"][0]["status"],
        "unverified"
    );
}

#[test]
fn scientific_capability_steps_fail_closed_on_invalid_intent_and_identity() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let archive = root.join("sample.aie");
    let gtf = root.join("genes.gtf");
    let other = root.join("other.gtf");
    let collision = root.join("collision.gtf");
    let legacy = root.join("legacy.aic");
    write_empty_archive(&archive);
    std::fs::write(
        &gtf,
        "chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id \"g1\"; gene_name \"G1\"; transcript_id \"t1\"; exon_id \"e1\";\n",
    )
    .unwrap();
    std::fs::write(
        &other,
        "chr1\ttest\texon\t101\t175\t.\t+\t.\tgene_id \"g2\"; gene_name \"G2\"; transcript_id \"t2\"; exon_id \"e2\";\n",
    )
    .unwrap();
    std::fs::write(
        &collision,
        concat!(
            "chr1\ttest\texon\t101\t125\t.\t+\t.\tgene_id \"dup.1\"; transcript_id \"d1\";\n",
            "chr1\ttest\texon\t151\t175\t.\t+\t.\tgene_id \"dup.2\"; transcript_id \"d2\";\n",
        ),
    )
    .unwrap();
    write_legacy_aic_without_transcript_ids(&legacy);
    add_with_assembly(&root, "sample", &archive, "archive", "GRCh38");
    add_with_assembly(&root, "wrong-build", &archive, "archive", "GRCh37");
    add_annotation(&root, "genes", &gtf, "GRCh38", "release");
    add_annotation(&root, "other", &other, "GRCh37", "other");
    add_annotation(&root, "collision", &collision, "GRCh38", "collision");
    add_annotation(&root, "legacy", &legacy, "GRCh38", "legacy");
    add(&root, "unidentified", &gtf, "annotation");
    let plan = root.join("plans/invalid-science.yaml");
    let cases = [
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: compare-annotations
    archive: sample
    annotation_a: genes
    annotation_b: other
"#,
            "annotation assembly mismatch",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: compare-annotations
    archive: sample
    annotation_a: genes
    annotation_b: genes
"#,
            "identical source content",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: compare-annotations
    archive: sample
    annotation_a: unidentified
    annotation_b: genes
"#,
            "has no scientific identity",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: compare-annotations
    archive: sample
    annotation_a: collision
    annotation_b: genes
"#,
            "normalization collision",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: compare-annotations
    archive: sample
    annotation_a: sample
    annotation_b: genes
"#,
            "requires annotation",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: compare-annotations
    archive: sample
    annotation_a: genes
    annotation_b: collision
    gene_key: exact
    format: tsv
"#,
            "requires an explicit table selection",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-transcript-ecs
    archive: sample
    feature: transcript:t1
    annotation: genes
"#,
            "require a gene feature",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-transcript-ecs
    archive: sample
    locus: chr1:100-150
    feature: g1
    annotation: genes
"#,
            "exactly one of `locus` or `feature`",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-transcript-ecs
    archive: sample
    locus: chr1:100-150
    annotation: genes
    format: tsv
    table: membership
"#,
            "requires emit_membership: true",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-transcript-ecs
    archive: sample
    locus: chr1:100-150
    annotation: genes
    max_ecs: 0
"#,
            "max_ecs must be greater than zero",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-transcript-ecs
    archive: sample
    locus: chr1:100-150
    annotation: unidentified
"#,
            "has no scientific identity",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-transcript-ecs
    archive: sample
    locus: chr1:100-150
    annotation: genes
    format: json
    table: counts
"#,
            "only valid with format tsv",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-transcript-ecs
    archive: wrong-build
    locus: chr1:100-150
    annotation: genes
"#,
            "assembly mismatch",
        ),
        (
            r#"schema_version: 1
steps:
  - id: q
    kind: query-transcript-ecs
    archive: sample
    locus: chr1:100-150
    annotation: legacy
"#,
            "stable transcript IDs",
        ),
    ];
    for (source, expected) in cases {
        std::fs::write(&plan, source).unwrap();
        let rejected = aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .output()
            .unwrap();
        assert!(
            !rejected.status.success(),
            "case unexpectedly succeeded: {source}"
        );
        let stderr = String::from_utf8_lossy(&rejected.stderr);
        assert!(
            stderr.contains(expected),
            "stderr did not contain {expected:?}: {stderr}"
        );
    }
}

#[test]
fn complex_query_collection_and_cohort_steps_compile_to_the_cli() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let data = root.join("data");
    let metadata = root.join("metadata");
    std::fs::create_dir_all(&data).unwrap();
    std::fs::create_dir_all(&metadata).unwrap();
    let archive_a = data.join("a.aie");
    let archive_b = data.join("b.aie");
    let collection = data.join("atlas.aicollection");
    let cells_a = metadata.join("cells-a.tsv");
    let cells_b = metadata.join("cells-b.tsv");
    let groups_a = metadata.join("groups-a.tsv");
    let groups_b = metadata.join("groups-b.tsv");
    write_empty_archive(&archive_a);
    write_empty_archive(&archive_b);
    write_empty_collection(&collection);
    for path in [&cells_a, &cells_b] {
        std::fs::write(path, b"cell-a\n").unwrap();
    }
    for path in [&groups_a, &groups_b] {
        std::fs::write(path, b"cell-a\tcase\n").unwrap();
    }
    let design = metadata.join("design.tsv");
    std::fs::write(
        &design,
        concat!(
            "sample\tcondition\tarchive\tcells\n",
            "a\tcase\t../data/a.aie\tcells-a.tsv\n",
            "b\tcontrol\t../data/b.aie\tcells-b.tsv\n",
        ),
    )
    .unwrap();
    for (name, path, kind) in [
        ("a", &archive_a, "archive"),
        ("b", &archive_b, "archive"),
        ("atlas", &collection, "collection"),
        ("cells-a", &cells_a, "cells"),
        ("cells-b", &cells_b, "cells"),
        ("groups-a", &groups_a, "groups"),
        ("groups-b", &groups_b, "groups"),
        ("design", &design, "design"),
    ] {
        add(&root, name, path, kind);
    }

    let plan = root.join("plans/complex.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
name: complex workflow
steps:
  - id: local-jset
    kind: query-jset
    archive: a
    include: ["chr1:100-200"]
    exclude: ["chr1:100-300"]
    uniform_output:
      format: tsv
    scope:
      cells: cells-a
      aggregation: bulk
  - id: collection-region
    kind: collection-region
    collection: atlas
    locus: chr1:1-1000
    explain_routing: true
    uniform_output:
      format: json
  - id: collection-junction
    kind: collection-junction
    collection: atlas
    locus: chr1:100-200
    verify_content: true
    uniform_output:
      format: text
  - id: cohort-events
    kind: cohort-events
    samples:
      a: a
      b: b
    groups:
      a: groups-a
      b: groups-b
    locus: chr1:1-1000
    event_types: [cassette]
    uniform_output:
      format: json
  - id: cohort-graph
    kind: cohort-splice-graph
    locus: chr1:1-1000
    design: design
    counts_only: true
    min_replicates: 1
    uniform_output:
      format: json
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(resolved["steps"].as_array().unwrap().len(), 5);
    assert!(resolved["steps"]
        .as_array()
        .unwrap()
        .iter()
        .all(|step| step.get("uniform_io").is_some()));
    for (index, expected) in [
        (
            0usize,
            serde_json::json!([
                "gravlax.query.jset.result.v1",
                "gravlax.query.jset.junctions.v1",
                "gravlax.query.jset.counts.v1"
            ]),
        ),
        (
            1,
            serde_json::json!([
                "gravlax.collection.region.result.v1",
                "gravlax.collection.region.samples.v1",
                "gravlax.collection.region.cells.v1"
            ]),
        ),
        (
            2,
            serde_json::json!([
                "gravlax.collection.junction.result.v1",
                "gravlax.collection.junction.samples.v1",
                "gravlax.collection.junction.cells.v1"
            ]),
        ),
        (
            3,
            serde_json::json!([
                "gravlax.cohort.events.result.v1",
                "gravlax.cohort.events.samples.v1",
                "gravlax.cohort.events.events.v1",
                "gravlax.cohort.events.components.v1",
                "gravlax.cohort.events.counts.v1"
            ]),
        ),
        (
            4,
            serde_json::json!([
                "gravlax.cohort.splice-graph.result.v1",
                "gravlax.cohort.splice-graph.samples.v1",
                "gravlax.cohort.splice-graph.nodes.v1",
                "gravlax.cohort.splice-graph.edges.v1",
                "gravlax.cohort.splice-graph.paths.v1",
                "gravlax.cohort.splice-graph.path-counts.v1",
                "gravlax.cohort.splice-graph.edge-counts.v1"
            ]),
        ),
    ] {
        assert_eq!(resolved["steps"][index]["output_schema_ids"], expected);
    }
    assert!(resolved["resources"].get("design").is_some());
    assert_eq!(
        resolved["resources"]["a"]["identity"]["scheme"],
        "aie-directory-root-v2"
    );
    assert_eq!(
        resolved["resources"]["atlas"]["identity"]["scheme"],
        "aicollection-directory-root-v1"
    );
    assert_eq!(
        resolved["resources"]["design"]["identity"]["scheme"],
        "full-file-blake3-v1"
    );
    assert_eq!(resolved["embedded_resources"].as_object().unwrap().len(), 4);
    assert_eq!(
        resolved["embedded_resources"]["design/a/archive"]["identity"]["scheme"],
        "aie-directory-root-v2"
    );
    let prepared = &resolved["steps"][4]["prepared_inputs"][0];
    assert_eq!(prepared["identity"]["scheme"], "full-file-blake3-v1");
    assert!(prepared["content"]
        .as_str()
        .unwrap()
        .contains(archive_a.to_str().unwrap()));
    assert!(resolved["steps"][3]["args"]
        .as_array()
        .unwrap()
        .iter()
        .any(|argument| argument.as_str().unwrap().starts_with("a=")));
}

#[test]
fn archive_query_ingest_replay_and_extension_steps_compile_to_the_cli() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let data = root.join("data");
    let metadata = root.join("metadata");
    std::fs::create_dir_all(&data).unwrap();
    std::fs::create_dir_all(&metadata).unwrap();
    let archive_a = data.join("a.aie");
    let archive_b = data.join("b.aie");
    let annotation = data.join("genes.gtf");
    let bam = data.join("reads.bam");
    let genome = data.join("genome.fa");
    let whitelist = metadata.join("whitelist.txt");
    let barcodes = metadata.join("barcodes.tsv");
    let cells = metadata.join("cells.tsv");
    write_empty_archive(&archive_a);
    write_empty_archive(&archive_b);
    std::fs::write(&annotation, b"# minimal GTF fixture\n").unwrap();
    std::fs::write(&bam, b"not decoded during plan checking\n").unwrap();
    std::fs::write(&genome, b">chr1\nACGT\n").unwrap();
    std::fs::write(&whitelist, b"AAAA\n").unwrap();
    std::fs::write(&barcodes, b"AAAA\n").unwrap();
    std::fs::write(&cells, b"AAAA\n").unwrap();
    for (name, path, kind) in [
        ("a", &archive_a, "archive"),
        ("b", &archive_b, "archive"),
        ("genes", &annotation, "annotation"),
        ("reads", &bam, "bam"),
        ("genome", &genome, "genome"),
        ("whitelist", &whitelist, "whitelist"),
        ("barcodes", &barcodes, "barcodes"),
        ("cells", &cells, "cells"),
    ] {
        add(&root, name, path, kind);
    }

    let plan = root.join("plans/remaining.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
name: remaining command surfaces
steps:
  - id: inspect
    kind: inspect-archive
    archive: a
    verify_content: true
    uniform_output:
      format: json
      output: results/inspect.json
  - id: region
    kind: query-region
    archive: a
    locus: chr1:1-100
    annotation: genes
    scope:
      cells: cells
      aggregation: cell
    uniform_output:
      format: tsv
      output: results/region.tsv
  - id: junction
    kind: query-junction
    archive: a
    locus: chr1:10-20
    uniform_output:
      format: text
  - id: junctions
    kind: query-junctions
    archive: a
    locus: chr1:1-100
    either: true
    with_cells: true
    min_cells: 1
    annotation: genes
    uniform_output:
      format: json
  - id: federate
    kind: federate-junction
    archives: [a, b]
    locus: chr1:10-20
    uniform_output:
      format: json
      output: results/federated.json
  - id: ingest
    kind: ingest-archive
    bam: reads
    whitelist: whitelist
    genome: genome
    output: results/new.aie
    uniform_report:
      format: json
      output: results/ingest-report.json
  - id: replay
    kind: replay-rows
    archive: a
    annotation: genes
    barcodes: barcodes
    velocity: true
    solo_strand: reverse
    out_dir: results/replay
    uniform_report:
      format: json
      output: results/replay-report.json
  - id: extend
    kind: extend-annotation
    archive: a
    annotation: genes
    genome: genome
    out_gtf: results/extended.gtf
    report: results/extension-report.tsv
    clip_any_strand: true
    uniform_report:
      format: json
      output: results/extend-report.json
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(resolved["steps"].as_array().unwrap().len(), 8);
    assert_eq!(resolved["resources"].as_object().unwrap().len(), 8);
    assert!(resolved["steps"]
        .as_array()
        .unwrap()
        .iter()
        .all(|step| step.get("uniform_io").is_some()));
    for (index, expected) in [
        (
            0usize,
            serde_json::json!([
                "gravlax.archive.inspect-report.v1",
                "gravlax.archive.section-accounting.v1"
            ]),
        ),
        (
            1,
            serde_json::json!([
                "gravlax.query.region.result.v1",
                "gravlax.query.region.counts.v1"
            ]),
        ),
        (
            2,
            serde_json::json!([
                "gravlax.query.junction.result.v1",
                "gravlax.query.junction.counts.v1"
            ]),
        ),
        (
            3,
            serde_json::json!([
                "gravlax.query.junctions.result.v1",
                "gravlax.query.junctions.junctions.v1",
                "gravlax.query.junctions.counts.v1"
            ]),
        ),
        (
            4,
            serde_json::json!([
                "gravlax.federate.junction.result.v1",
                "gravlax.federate.junction.archives.v1",
                "gravlax.federate.junction.counts.v1"
            ]),
        ),
        (
            5,
            serde_json::json!([
                "gravlax.archive.v2",
                "gravlax.archive.ingest-report.v1",
                "gravlax.archive.section-accounting.v1"
            ]),
        ),
        (
            6,
            serde_json::json!([
                "gravlax.replay.gene-matrix.v1",
                "gravlax.replay.mex-artifact.v1",
                "gravlax.replay.velocity-matrices.v1",
                "gravlax.archive.replay-report.v1",
                "gravlax.archive.replay-artifact-files.v1"
            ]),
        ),
        (
            7,
            serde_json::json!([
                "gravlax.annotation.gtf.v1",
                "gravlax.annotation.extension-report.v1",
                "gravlax.extend.result.v1",
                "gravlax.extend.artifacts.v1",
                "gravlax.extend.genes.v1"
            ]),
        ),
    ] {
        assert_eq!(resolved["steps"][index]["output_schema_ids"], expected);
    }
    let replay_args = resolved["steps"][6]["args"].as_array().unwrap();
    assert!(replay_args
        .windows(2)
        .any(|pair| pair[0] == "--report-format" && pair[1] == "json"));
    assert!(replay_args
        .windows(2)
        .any(|pair| pair[0] == "--report-output"));
}

#[test]
fn plan_run_persists_resolution_and_uses_existing_command() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let gtf = root.join("data/genes.gtf");
    std::fs::create_dir_all(gtf.parent().unwrap()).unwrap();
    std::fs::write(
        &gtf,
        concat!(
            "chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\"; gene_name \"G1\";\n",
            "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
            "chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_id \"e1\";\n",
        ),
    )
    .unwrap();
    add(&root, "genes", &gtf, "annotation");
    let plan = root.join("plans/compile.json");
    std::fs::write(
        &plan,
        serde_json::to_vec_pretty(&serde_json::json!({
            "schema_version": 1,
            "name": "compile annotation",
            "steps": [{
                "id": "compile",
                "kind": "compile-annotation",
                "annotation": "genes",
                "output": "results/genes.aic",
                "uniform_report": {
                    "format": "json",
                    "output": "results/compile-report.json"
                }
            }]
        }))
        .unwrap(),
    )
    .unwrap();

    let dry_run = success(
        aie()
            .args(["plan", "run"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--dry-run")
            .output()
            .unwrap(),
    );
    assert!(String::from_utf8_lossy(&dry_run.stdout).contains("dry run"));
    assert!(!root.join("results/genes.aic").exists());
    assert_eq!(
        std::fs::read_dir(root.join(".aie/resolved-plans"))
            .unwrap()
            .count(),
        0
    );

    let first_run = success(
        aie()
            .args(["plan", "run"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .output()
            .unwrap(),
    );
    assert!(String::from_utf8_lossy(&first_run.stdout).contains("completion:"));
    assert!(root.join("results/genes.aic").is_file());
    let report_path = root.join("results/compile-report.json");
    assert!(report_path.is_file());
    let report: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&report_path).unwrap()).unwrap();
    let report_json = serde_json::to_string(&report).unwrap();
    assert!(
        report_json.contains(root.join("results/genes.aic").to_str().unwrap()),
        "{report_json}"
    );
    assert!(
        !report_json.contains(".aie-stage-"),
        "planned report leaked a staging path: {report_json}"
    );
    let snapshots = std::fs::read_dir(root.join(".aie/resolved-plans"))
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(snapshots.len(), 1);
    let snapshot: serde_json::Value =
        serde_json::from_slice(&std::fs::read(snapshots[0].path()).unwrap()).unwrap();
    assert_eq!(snapshot["schema_version"], 6);
    assert_eq!(snapshot["producer"]["name"], "aie");
    assert_eq!(
        snapshot["producer"]["executable_identity"]["scheme"],
        "full-file-blake3-v1"
    );
    assert_eq!(snapshot["steps"][0]["kind"], "compile-annotation");
    assert!(snapshot["steps"][0]["outputs"][0]["path"]
        .as_str()
        .unwrap()
        .ends_with("results/genes.aic"));
    let staging_path = snapshot["steps"][0]["outputs"][0]["staging_path"]
        .as_str()
        .unwrap();
    assert!(staging_path.ends_with("results/.genes.aie-stage-compile.aic"));
    assert!(snapshot["steps"][0]["args"]
        .as_array()
        .unwrap()
        .iter()
        .any(|argument| argument.as_str() == Some(staging_path)));
    assert_eq!(snapshot["steps"][0]["uniform_io"]["kind"], "report");
    assert_eq!(
        snapshot["steps"][0]["uniform_io"]["output"],
        report_path.display().to_string()
    );

    let completion_files = regular_files_below(&root.join(".aie/completions"));
    assert_eq!(completion_files.len(), 1);
    let completion_bytes = std::fs::read(&completion_files[0]).unwrap();
    let completion: serde_json::Value = serde_json::from_slice(&completion_bytes).unwrap();
    assert_eq!(completion["schema_version"], 2);
    assert_eq!(completion["step_id"], "compile");
    assert_eq!(
        completion["outputs"][0]["identity"]["scheme"],
        "full-file-blake3-v1"
    );

    let overwrite = aie()
        .args(["plan", "run"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!overwrite.status.success());
    let overwrite_stderr = String::from_utf8_lossy(&overwrite.stderr);
    assert!(
        overwrite_stderr.contains("refusing to overwrite output"),
        "{overwrite_stderr}"
    );

    let resumed = success(
        aie()
            .args(["plan", "run"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--resume")
            .output()
            .unwrap(),
    );
    let resumed_stdout = String::from_utf8_lossy(&resumed.stdout);
    assert!(resumed_stdout.contains("exact verified outputs"));
    assert!(resumed_stdout.contains("exact completion verified"));

    std::fs::remove_file(&completion_files[0]).unwrap();
    let unverified = aie()
        .args(["plan", "run"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .arg("--resume")
        .output()
        .unwrap();
    assert!(!unverified.status.success());
    assert!(
        String::from_utf8_lossy(&unverified.stderr).contains("without an exact completion record")
    );
    std::fs::write(&completion_files[0], completion_bytes).unwrap();

    std::fs::write(root.join("results/genes.aic"), b"tampered\n").unwrap();
    let tampered = aie()
        .args(["plan", "run"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .arg("--resume")
        .output()
        .unwrap();
    assert!(!tampered.status.success());
    assert!(String::from_utf8_lossy(&tampered.stderr).contains("output identities differ"));
}

#[test]
fn prior_step_outputs_drive_compile_and_replay_end_to_end() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let data = root.join("data");
    std::fs::create_dir_all(&data).unwrap();
    let bam = data.join("reads.bam");
    let whitelist = data.join("whitelist.txt");
    let barcodes = data.join("barcodes.tsv");
    let gtf = data.join("genes.gtf");
    write_replay_fixture_bam(&bam);
    std::fs::write(&whitelist, format!("{TEST_BARCODE}\n")).unwrap();
    std::fs::write(&barcodes, format!("{TEST_BARCODE}\n")).unwrap();
    std::fs::write(
        &gtf,
        concat!(
            "chr1\ttest\tgene\t101\t150\t.\t+\t.\tgene_id \"g1\"; gene_name \"G1\";\n",
            "chr1\ttest\ttranscript\t101\t150\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
            "chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_id \"e1\";\n",
        ),
    )
    .unwrap();
    for (name, path, kind) in [
        ("reads", &bam, "bam"),
        ("whitelist", &whitelist, "whitelist"),
        ("barcodes", &barcodes, "barcodes"),
        ("genes", &gtf, "annotation"),
    ] {
        add(&root, name, path, kind);
    }
    let plan = root.join("plans/dataflow.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
name: typed dataflow
steps:
  - id: ingest
    kind: ingest-archive
    bam: reads
    whitelist: whitelist
    output: results/sample.aie
    zstd_level: 1
    chunk_mb: 1
  - id: compile
    kind: compile-annotation
    annotation: genes
    output: results/genes.aic
  - id: replay
    kind: replay-rows
    archive: "step:ingest"
    annotation: "step:compile:annotation"
    barcodes: barcodes
    out_dir: results/matrix
    uniform_report:
      format: json
      output: results/replay-report.json
"#,
    )
    .unwrap();

    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(resolved["schema_version"], 6);
    assert_eq!(resolved["steps"][0]["outputs"][0]["name"], "archive");
    assert_eq!(
        resolved["steps"][1]["outputs"][0]["resource_kind"],
        "annotation"
    );
    assert_eq!(
        resolved["steps"][2]["step_inputs"][0]["producer_step"],
        "compile"
    );
    assert_eq!(
        resolved["steps"][2]["step_inputs"]
            .as_array()
            .unwrap()
            .len(),
        2
    );

    success(
        aie()
            .args(["plan", "run"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .output()
            .unwrap(),
    );
    assert!(root.join("results/sample.aie").is_file());
    assert!(root.join("results/genes.aic").is_file());
    assert!(root.join("results/matrix/matrix.mtx").is_file());
    let replay_report: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("results/replay-report.json")).unwrap())
            .unwrap();
    let replay_report_json = serde_json::to_string(&replay_report).unwrap();
    assert!(
        replay_report_json.contains(root.join("results/matrix").to_str().unwrap()),
        "{replay_report_json}"
    );
    assert!(
        !replay_report_json.contains(".aie-stage-"),
        "planned replay report leaked a staging path: {replay_report_json}"
    );

    let completion_files = regular_files_below(&root.join(".aie/completions"));
    assert_eq!(completion_files.len(), 3);
    let replay_completion = completion_files
        .iter()
        .find(|path| path.file_name().unwrap() == "replay.json")
        .unwrap();
    let completion: serde_json::Value =
        serde_json::from_slice(&std::fs::read(replay_completion).unwrap()).unwrap();
    assert_eq!(completion["schema_version"], 2);
    assert_eq!(completion["inputs"].as_array().unwrap().len(), 2);

    let resumed = success(
        aie()
            .args(["plan", "run"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--resume")
            .output()
            .unwrap(),
    );
    assert!(String::from_utf8_lossy(&resumed.stdout).contains("skipping step 3/3 `replay`"));
}

#[test]
fn step_output_references_reject_forward_cycles_ambiguity_and_wrong_types() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let gtf = root.join("genes.gtf");
    let archive = root.join("sample.aie");
    std::fs::write(&gtf, b"# fixture\n").unwrap();
    write_empty_archive(&archive);
    add(&root, "genes", &gtf, "annotation");
    add(&root, "sample", &archive, "archive");

    let plan = root.join("plans/invalid-dataflow.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: consume
    kind: inspect-archive
    archive: "step:compile"
  - id: compile
    kind: compile-annotation
    annotation: genes
    output: results/genes.aic
"#,
    )
    .unwrap();
    let forward = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!forward.status.success());
    assert!(String::from_utf8_lossy(&forward.stderr).contains("before its producer"));

    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: compile
    kind: compile-annotation
    annotation: genes
    output: results/genes.aic
  - id: consume
    kind: inspect-archive
    archive: "step:compile"
"#,
    )
    .unwrap();
    let wrong_type = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!wrong_type.status.success());
    assert!(String::from_utf8_lossy(&wrong_type.stderr)
        .contains("is annotation, but this input requires archive"));

    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: cycle
    kind: inspect-archive
    archive: "step:cycle"
"#,
    )
    .unwrap();
    let cycle = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!cycle.status.success());
    assert!(String::from_utf8_lossy(&cycle.stderr).contains("before its producer"));

    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: extend
    kind: extend-annotation
    archive: sample
    annotation: genes
    out_gtf: results/extended.gtf
    report: results/report.tsv
  - id: ambiguous
    kind: inspect-archive
    archive: "step:extend"
"#,
    )
    .unwrap();
    let ambiguous = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!ambiguous.status.success());
    assert!(String::from_utf8_lossy(&ambiguous.stderr).contains("is ambiguous"));
}

#[test]
fn multi_sample_steps_reject_archive_aliases_and_identity_cache_keeps_roles_distinct() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let archive = root.join("sample.aie");
    write_empty_archive(&archive);
    add(&root, "archive", &archive, "archive");
    add(&root, "alias", &archive, "archive");
    add(&root, "generic", &archive, "file");
    let plan = root.join("plans/aliases.yaml");

    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: native
    kind: inspect-archive
    archive: archive
  - id: generic
    kind: inspect-archive
    archive: generic
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(
        resolved["resources"]["archive"]["identity"]["scheme"],
        "aie-directory-root-v2"
    );
    assert_eq!(
        resolved["resources"]["generic"]["identity"]["scheme"],
        "full-file-blake3-v1"
    );

    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: duplicate
    kind: federate-junction
    archives: [archive, alias]
    locus: chr1:10-20
"#,
    )
    .unwrap();
    let duplicate_federate = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!duplicate_federate.status.success());
    assert!(
        String::from_utf8_lossy(&duplicate_federate.stderr).contains("resolve to the same archive")
    );

    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: duplicate
    kind: cohort-events
    samples:
      one: archive
      two: alias
    locus: chr1:1-100
"#,
    )
    .unwrap();
    let duplicate_cohort = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!duplicate_cohort.status.success());
    assert!(
        String::from_utf8_lossy(&duplicate_cohort.stderr).contains("resolve to the same archive")
    );
}

#[test]
fn unsafe_output_and_unknown_operation_fail_before_execution() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let gtf = root.join("genes.gtf");
    std::fs::write(&gtf, b"# fixture\n").unwrap();
    add(&root, "genes", &gtf, "annotation");
    let plan = root.join("plans/unsafe.yaml");
    std::fs::write(
        &plan,
        "schema_version: 1\nsteps:\n  - id: escape\n    kind: compile-annotation\n    annotation: genes\n    output: ../outside.aic\n",
    )
    .unwrap();
    let unsafe_result = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!unsafe_result.status.success());
    assert!(String::from_utf8_lossy(&unsafe_result.stderr).contains("must not contain '..'"));

    std::fs::write(
        &plan,
        "schema_version: 1\nsteps:\n  - id: future\n    kind: annotation-compare\n",
    )
    .unwrap();
    let unknown = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!unknown.status.success());
    let stderr = String::from_utf8_lossy(&unknown.stderr);
    assert!(stderr.contains("unsupported plan step kind `annotation-compare`"));
    assert!(stderr.contains("query-events"));

    std::fs::write(
        &plan,
        "schema_version: 2\nsteps:\n  - id: compile\n    kind: compile-annotation\n    annotation: genes\n    output: results/genes.aic\n",
    )
    .unwrap();
    let future_schema = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!future_schema.status.success());
    assert!(String::from_utf8_lossy(&future_schema.stderr)
        .contains("unsupported plan schema version 2"));
}

#[test]
fn uniform_plan_io_is_explicit_staged_and_legacy_argv_is_unchanged() {
    let scratch = Scratch::new();
    let root = scratch.0.join("workspace");
    init(&root);
    let archive = root.join("sample.aie");
    write_empty_archive(&archive);
    add(&root, "sample", &archive, "archive");
    let gtf = root.join("genes.gtf");
    std::fs::write(&gtf, b"# fixture\n").unwrap();
    add(&root, "genes", &gtf, "annotation");

    let plan = root.join("plans/io.yaml");
    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: legacy
    kind: query-region
    archive: sample
    locus: chr1:0-1
    format: json
    output: results/legacy.json
  - id: uniform
    kind: query-region
    archive: sample
    locus: chr1:0-1
    uniform_output:
      format: json
      output: results/uniform.json
  - id: inspect
    kind: inspect-archive
    archive: sample
    uniform_output:
      format: tsv
  - id: compile
    kind: compile-annotation
    annotation: genes
    output: results/genes.aic
    uniform_report:
      format: json
      output: results/compile-report.json
"#,
    )
    .unwrap();
    let checked = success(
        aie()
            .args(["plan", "check"])
            .arg(&plan)
            .args(["--project"])
            .arg(&root)
            .arg("--json")
            .output()
            .unwrap(),
    );
    let resolved: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(resolved["schema_version"], 6);
    assert_eq!(resolved["plan_schema_version"], 1);
    assert_eq!(
        resolved["producer"]["plan_engine"],
        "aie-declarative-plan-v2"
    );

    let legacy = &resolved["steps"][0];
    let legacy_args = legacy["args"].as_array().unwrap();
    assert!(legacy_args.iter().any(|value| value == "--json"));
    assert!(!legacy_args.iter().any(|value| value == "--format"));
    assert!(legacy.get("uniform_io").is_none());
    assert_eq!(
        legacy["stdout"],
        root.join("results/legacy.json").display().to_string()
    );

    let uniform = &resolved["steps"][1];
    assert!(uniform["stdout"].is_null());
    assert_eq!(uniform["uniform_io"]["kind"], "result");
    assert_eq!(uniform["uniform_io"]["format"], "json");
    assert_eq!(
        uniform["uniform_io"]["publication"],
        "atomic-no-clobber-file"
    );
    assert_eq!(
        uniform["uniform_io"]["output"],
        root.join("results/uniform.json").display().to_string()
    );
    assert_eq!(uniform["outputs"][0]["name"], "result");
    assert_eq!(
        uniform["outputs"][0]["staging_path"],
        root.join("results/.uniform.aie-stage-uniform.json")
            .display()
            .to_string()
    );
    assert_eq!(
        uniform["output_schema_ids"],
        serde_json::json!([
            "gravlax.query.region.result.v1",
            "gravlax.query.region.counts.v1"
        ])
    );
    let uniform_args = uniform["args"].as_array().unwrap();
    assert!(uniform_args
        .windows(2)
        .any(|pair| { pair[0] == "--format" && pair[1] == "json" }));
    assert!(uniform_args.windows(2).any(|pair| {
        pair[0] == "--output"
            && pair[1]
                == root
                    .join("results/.uniform.aie-stage-uniform.json")
                    .display()
                    .to_string()
    }));

    let inspect = &resolved["steps"][2];
    assert_eq!(inspect["uniform_io"]["publication"], "stdout");
    assert!(inspect["uniform_io"].get("output").is_none());
    assert_eq!(
        inspect["output_schema_ids"],
        serde_json::json!([
            "gravlax.archive.inspect-report.v1",
            "gravlax.archive.section-accounting.v1"
        ])
    );

    let compile = &resolved["steps"][3];
    assert_eq!(compile["uniform_io"]["kind"], "report");
    assert!(compile["args"]
        .as_array()
        .unwrap()
        .windows(2)
        .any(|pair| { pair[0] == "--report-format" && pair[1] == "json" }));
    assert!(compile["output_schema_ids"]
        .as_array()
        .unwrap()
        .iter()
        .any(|value| value == "gravlax.annotation.compile.result.v1"));

    std::fs::write(
        &plan,
        r#"schema_version: 1
steps:
  - id: conflict
    kind: query-region
    archive: sample
    locus: chr1:0-1
    format: json
    uniform_output:
      format: json
"#,
    )
    .unwrap();
    let conflict = aie()
        .args(["plan", "check"])
        .arg(&plan)
        .args(["--project"])
        .arg(&root)
        .output()
        .unwrap();
    assert!(!conflict.status.success());
    assert!(String::from_utf8_lossy(&conflict.stderr)
        .contains("uniform_output cannot be combined with legacy format or output fields"));
}
