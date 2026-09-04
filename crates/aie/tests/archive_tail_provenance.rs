use evidence_io::alignment_provenance::{
    AlignmentProvenanceManifest, GenomeBindingAction, GenomeReferenceBinding,
    JunctionCatalogueRole, JunctionDiscoveryMode, ProvenanceStatus, ALIGNMENT_PROVENANCE_SECTION,
    JUNCTION_CATALOGUE_SECTION, MOLECULAR_EVIDENCE_SCHEMA,
};
use evidence_io::format::{Cursor, SectionReader, SectionWriter};
use evidence_io::terminal_tail::{self, TerminalTailMetadata};
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
const UMI: &str = "ACGTACGTACGT";

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-tail-provenance-{}-{nonce}",
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

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "command failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn tail_record(name: &str, start: usize, clip: usize) -> RecordBuf {
    tail_record_signal(name, start, clip, clip, clip)
}

fn tail_record_signal(
    name: &str,
    start: usize,
    clip: usize,
    tail_bases: usize,
    terminal_run: usize,
) -> RecordBuf {
    let cigar: Cigar = [Op::new(Kind::Match, 20), Op::new(Kind::SoftClip, clip)]
        .into_iter()
        .collect();
    let data: Data = [
        (Tag::new(b'C', b'R'), Value::from(BARCODE)),
        (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
        (Tag::new(b'U', b'R'), Value::from(UMI)),
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
    ]
    .into_iter()
    .collect();
    assert!(terminal_run <= tail_bases && tail_bases <= clip);
    let mut clipped = vec![b'C'; clip];
    for base in &mut clipped[..tail_bases - terminal_run] {
        *base = b'A';
    }
    for base in &mut clipped[clip - terminal_run..] {
        *base = b'A';
    }
    let mut sequence = vec![b'C'; 20];
    sequence.extend(clipped);
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(start).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(sequence))
        .set_quality_scores(QualityScores::from(vec![30; 20 + clip]))
        .set_data(data)
        .build()
}

fn multimapper_record(name: &str, start: usize, flags: Flags) -> RecordBuf {
    // Deliberately tail-positive, but NH=2: without a placement association an exact terminal
    // anchor cannot be attributed among its alternatives, so the v1 tail rule must exclude it.
    let cigar: Cigar = [Op::new(Kind::Match, 20), Op::new(Kind::SoftClip, 6)]
        .into_iter()
        .collect();
    let data: Data = [
        (Tag::new(b'C', b'R'), Value::from(BARCODE)),
        (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
        (Tag::new(b'U', b'R'), Value::from("TGCATGCATGCA")),
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(2u8)),
    ]
    .into_iter()
    .collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(start).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from([vec![b'C'; 20], vec![b'A'; 6]].concat()))
        .set_quality_scores(QualityScores::from(vec![30; 26]))
        .set_data(data)
        .build()
}

fn write_fixture_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(10_000).unwrap()),
        )
        .build();
    // The 111 records are internal to the chain span. Their witness must survive reduction. The
    // duplicate anchor is adversarial for five-bit-first ranking: 29/31 is stronger than 36/40,
    // although saturating both counts to five bits first would make 36/40 appear to be 31/31.
    let records = [
        tail_record("left", 101, 6),
        tail_record_signal("middle-long", 111, 40, 36, 10),
        tail_record_signal("middle-pure", 111, 31, 29, 20),
        tail_record("right", 121, 6),
        multimapper_record("multimap", 501, Flags::empty()),
        multimapper_record("multimap", 701, Flags::SECONDARY),
    ];
    let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for record in records {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn write_orientation_fixture_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(10_000).unwrap()),
        )
        .build();
    let tagged_data = |umi: &'static str| -> Data {
        [
            (Tag::new(b'C', b'R'), Value::from(BARCODE)),
            (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
            (Tag::new(b'U', b'R'), Value::from(umi)),
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
        ]
        .into_iter()
        .collect()
    };
    let reverse = RecordBuf::builder()
        .set_name("reverse-leading-t")
        .set_flags(Flags::REVERSE_COMPLEMENTED)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(401).unwrap())
        .set_cigar(
            [Op::new(Kind::SoftClip, 6), Op::new(Kind::Match, 20)]
                .into_iter()
                .collect::<Cigar>(),
        )
        .set_sequence(Sequence::from([vec![b'T'; 6], vec![b'C'; 20]].concat()))
        .set_quality_scores(QualityScores::from(vec![30; 26]))
        .set_data(tagged_data("CCCCCCCCCCCC"))
        .build();
    let forward_wrong_edge = RecordBuf::builder()
        .set_name("forward-leading-a-wrong-edge")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(801).unwrap())
        .set_cigar(
            [Op::new(Kind::SoftClip, 6), Op::new(Kind::Match, 20)]
                .into_iter()
                .collect::<Cigar>(),
        )
        .set_sequence(Sequence::from([vec![b'A'; 6], vec![b'C'; 20]].concat()))
        .set_quality_scores(QualityScores::from(vec![30; 26]))
        .set_data(tagged_data("GGGGGGGGGGGG"))
        .build();
    let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    // A multi-read unique chain keeps every established core stream populated in this compact
    // end-to-end fixture; its three forward tails also provide an opposite-orientation control.
    writer
        .write_alignment_record(&header, &tail_record("control-left", 101, 6))
        .unwrap();
    writer
        .write_alignment_record(&header, &tail_record("control-middle", 111, 6))
        .unwrap();
    writer
        .write_alignment_record(&header, &tail_record("control-right", 121, 6))
        .unwrap();
    writer.write_alignment_record(&header, &reverse).unwrap();
    writer
        .write_alignment_record(&header, &forward_wrong_edge)
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &multimapper_record("orientation-multimap", 1201, Flags::empty()),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &multimapper_record("orientation-multimap", 1401, Flags::SECONDARY),
        )
        .unwrap();
    writer.try_finish().unwrap();
}

fn write_duplicate_origin_fixture_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(10_000).unwrap()),
        )
        .build();
    let data = || -> Data {
        [
            (Tag::new(b'C', b'R'), Value::from(BARCODE)),
            (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
            (Tag::new(b'U', b'R'), Value::from(UMI)),
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
        ]
        .into_iter()
        .collect()
    };
    // Both records have the same corrected cell, UMI, strand, and cleavage anchor (9000), but
    // their starts are farther apart than the 2 kb locus gap and therefore become distinct
    // MolRecs. The earlier serialized record carries a junction but has the weaker tail signal.
    let spliced_weak = RecordBuf::builder()
        .set_name("spliced-weak")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(4001).unwrap())
        .set_cigar(
            [
                Op::new(Kind::Match, 20),
                Op::new(Kind::Skip, 4960),
                Op::new(Kind::Match, 20),
                Op::new(Kind::SoftClip, 6),
            ]
            .into_iter()
            .collect::<Cigar>(),
        )
        .set_sequence(Sequence::from([vec![b'C'; 40], vec![b'A'; 6]].concat()))
        .set_quality_scores(QualityScores::from(vec![30; 46]))
        .set_data(data())
        .build();
    let unspliced_strong = RecordBuf::builder()
        .set_name("unspliced-strong")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(8981).unwrap())
        .set_cigar(
            [Op::new(Kind::Match, 20), Op::new(Kind::SoftClip, 10)]
                .into_iter()
                .collect::<Cigar>(),
        )
        .set_sequence(Sequence::from([vec![b'C'; 20], vec![b'A'; 10]].concat()))
        .set_quality_scores(QualityScores::from(vec![30; 30]))
        .set_data(data())
        .build();
    let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    // Keep all established core streams populated so this remains an archive-level regression,
    // not a special extractor-only fixture.
    for record in [
        tail_record("control-left", 101, 6),
        tail_record_signal("control-middle-long", 111, 40, 36, 10),
        tail_record_signal("control-middle-pure", 111, 31, 29, 20),
        tail_record("control-right", 121, 6),
        multimapper_record("control-multimap", 501, Flags::empty()),
        multimapper_record("control-multimap", 701, Flags::SECONDARY),
    ] {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer
        .write_alignment_record(&header, &spliced_weak)
        .unwrap();
    writer
        .write_alignment_record(&header, &unspliced_strong)
        .unwrap();
    writer.try_finish().unwrap();
}

fn ingest(root: &Path, tails: bool) -> PathBuf {
    let bam = root.join("reads.bam");
    let whitelist = root.join("whitelist.txt");
    let catalogue = root.join("SJ.out.tab");
    let reads_1 = root.join("R1.fastq.gz");
    let reads_2 = root.join("R2.fastq.gz");
    let annotation = root.join("alignment.gtf");
    let log = root.join("Log.final.out");
    if !bam.exists() {
        write_fixture_bam(&bam);
        std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
        std::fs::write(
            &catalogue,
            "# exact table retained by the caller\n\nchr1\t120\t140\t1\t0\t0\nchr1\t220\t240\t1\t0\t0\n",
        )
        .unwrap();
        std::fs::write(&reads_1, b"read-one-bytes").unwrap();
        std::fs::write(&reads_2, b"read-two-bytes").unwrap();
        std::fs::write(&annotation, b"annotation-bytes\n").unwrap();
        std::fs::write(&log, b"resolved defaults\n").unwrap();
    }
    let out = root.join(if tails { "with-tail.aie" } else { "core.aie" });
    let mut command = Command::new(env!("CARGO_BIN_EXE_aie"));
    command.args([
        "ingest-archive",
        bam.to_str().unwrap(),
        "--whitelist",
        whitelist.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--zstd-level",
        "1",
        "--junction-discovery",
        "per-library-two-pass",
        "--junction-catalogue",
        catalogue.to_str().unwrap(),
        "--alignment-annotation",
        annotation.to_str().unwrap(),
        "--alignment-index-identity",
        "star-index:fixture-v1",
        "--alignment-input",
        reads_1.to_str().unwrap(),
        "--alignment-input",
        reads_2.to_str().unwrap(),
        "--alignment-log",
        log.to_str().unwrap(),
        "--alignment-chemistry",
        "10x-3p-v3",
    ]);
    if tails {
        command.arg("--terminal-tails");
    }
    assert_success(&command.output().unwrap());
    out
}

fn chunk_molecule_counts(reader: &mut SectionReader) -> Vec<u32> {
    let raw = reader.read("index.chunks").unwrap();
    let mut cursor = Cursor::new(&raw);
    let mut counts = Vec::new();
    while !cursor.is_empty() {
        cursor.varint().unwrap(); // chrom
        cursor.varint().unwrap(); // bin start
        counts.push(u32::try_from(cursor.varint().unwrap()).unwrap());
        cursor.varint().unwrap(); // class base
        cursor.varint().unwrap(); // max anchor delta
        cursor.varint().unwrap(); // cell count
    }
    counts
}

#[test]
fn tail_capture_is_sparse_lossless_for_the_frozen_rule_and_core_bytes_are_unchanged() {
    let scratch = Scratch::new();
    let core_path = ingest(&scratch.0, false);
    let tail_path = ingest(&scratch.0, true);

    let terminal_query = |archive: &Path| {
        Command::new(env!("CARGO_BIN_EXE_aie"))
            .args([
                "query",
                archive.to_str().unwrap(),
                "cooccur",
                "--predicate",
                "tail=terminal:chr1:0-1000:+",
                "--where",
                "tail",
                "--universe",
                "tail",
                "--format",
                "json",
            ])
            .output()
            .unwrap()
    };
    let valid_query = terminal_query(&tail_path);
    assert_success(&valid_query);
    let valid_result: serde_json::Value = serde_json::from_slice(&valid_query.stdout).unwrap();
    let patterns = valid_result["data"]["tables"]
        .as_array()
        .unwrap()
        .iter()
        .find(|table| table["name"] == "patterns")
        .unwrap();
    let fields = patterns["schema"]["fields"]
        .as_array()
        .unwrap()
        .iter()
        .map(|field| field["name"].as_str().unwrap())
        .collect::<Vec<_>>();
    let row = patterns["rows"][0].as_array().unwrap();
    let field = |name: &str| row[fields.iter().position(|field| *field == name).unwrap()].clone();
    assert_eq!(field("pattern_mask"), "0x0000000000000001");
    assert_eq!(field("selection_state"), "true");
    assert_eq!(field("selected"), true);
    assert_eq!(field("evidence_units"), 1);

    let unavailable_query = terminal_query(&core_path);
    assert!(!unavailable_query.status.success());
    assert!(String::from_utf8_lossy(&unavailable_query.stderr)
        .contains("terminal predicates require the lossless terminal-tail capability"));

    let mut core = SectionReader::open(&core_path).unwrap();
    let mut tails = SectionReader::open(&tail_path).unwrap();

    let core_names: Vec<String> = core.names().map(str::to_owned).collect();
    for name in core_names {
        if name == "meta" || name == ALIGNMENT_PROVENANCE_SECTION {
            continue;
        }
        assert_eq!(
            core.read(&name).unwrap(),
            tails.read(&name).unwrap(),
            "core evidence section {name} changed when tail capture was enabled"
        );
    }

    let meta: serde_json::Value = serde_json::from_slice(&tails.read("meta").unwrap()).unwrap();
    assert_eq!(meta["evidence_schema"], MOLECULAR_EVIDENCE_SCHEMA);
    let capability: TerminalTailMetadata =
        serde_json::from_value(meta["terminal_tail"].clone()).unwrap();
    capability.validate().unwrap();
    assert_eq!(
        capability.alignment_scope,
        "mapped-primary-nonsupplementary-explicit-nh1"
    );
    assert_eq!(capability.selected_molecules, 1);
    assert_eq!(capability.events, 3);
    assert_eq!(capability.chunks, 1);

    let molecule_counts = chunk_molecule_counts(&mut tails);
    let routes = terminal_tail::decode_index(
        &tails
            .read(terminal_tail::TERMINAL_TAIL_INDEX_SECTION)
            .unwrap(),
        molecule_counts.len(),
    )
    .unwrap();
    assert_eq!(routes.len(), 1);
    let route = routes[0];
    let selected = terminal_tail::decode_chunk(
        &tails.read(&format!("tail.c{}", route.chunk)).unwrap(),
        molecule_counts[route.chunk as usize],
    )
    .unwrap();
    assert_eq!(selected.len(), 1);
    assert_eq!(
        selected[0]
            .events
            .iter()
            .map(|event| event.anchor_delta)
            .collect::<Vec<_>>(),
        vec![20, 30, 40]
    );
    assert_eq!(selected[0].events[1].signal.clip_len, 31);
    assert_eq!(selected[0].events[1].signal.tail_bases, 29);
    assert_eq!(selected[0].events[1].signal.terminal_run, 20);
}

#[test]
fn reverse_tail_uses_leading_t_and_zero_based_aligned_start() {
    let scratch = Scratch::new();
    // Populate the common whitelist/provenance fixtures, then replace only the BAM before the
    // tail-enabled ingest. The unused core archive remains a separate immutable artifact.
    let _ = ingest(&scratch.0, false);
    write_orientation_fixture_bam(&scratch.0.join("reads.bam"));
    let archive = ingest(&scratch.0, true);
    let mut reader = SectionReader::open(&archive).unwrap();
    let meta: serde_json::Value = serde_json::from_slice(&reader.read("meta").unwrap()).unwrap();
    let capability: TerminalTailMetadata =
        serde_json::from_value(meta["terminal_tail"].clone()).unwrap();
    assert_eq!(capability.selected_molecules, 2);
    assert_eq!(capability.events, 4);

    let molecule_counts = chunk_molecule_counts(&mut reader);
    let routes = terminal_tail::decode_index(
        &reader
            .read(terminal_tail::TERMINAL_TAIL_INDEX_SECTION)
            .unwrap(),
        molecule_counts.len(),
    )
    .unwrap();
    assert_eq!(routes.len(), 1);
    assert_eq!(routes[0].min_anchor, 120);
    assert_eq!(routes[0].max_anchor, 400);
    let selected = terminal_tail::decode_chunk(
        &reader.read(&format!("tail.c{}", routes[0].chunk)).unwrap(),
        molecule_counts[routes[0].chunk as usize],
    )
    .unwrap();
    assert_eq!(selected.len(), 2);
    let reverse_events = selected
        .iter()
        .flat_map(|molecule| &molecule.events)
        .filter(|event| event.reverse)
        .collect::<Vec<_>>();
    assert_eq!(reverse_events.len(), 1);
    assert_eq!(reverse_events[0].anchor_delta, 0);
    assert_eq!(reverse_events[0].signal.clip_len, 6);
}

#[test]
fn globally_deduplicated_tail_keeps_the_selected_witness_record_identity() {
    let scratch = Scratch::new();
    let _ = ingest(&scratch.0, false);
    write_duplicate_origin_fixture_bam(&scratch.0.join("reads.bam"));
    let archive = ingest(&scratch.0, true);

    let mut reader = SectionReader::open(&archive).unwrap();
    let molecule_counts = chunk_molecule_counts(&mut reader);
    let routes = terminal_tail::decode_index(
        &reader
            .read(terminal_tail::TERMINAL_TAIL_INDEX_SECTION)
            .unwrap(),
        molecule_counts.len(),
    )
    .unwrap();
    assert_eq!(routes.len(), 1);
    let selected = terminal_tail::decode_chunk(
        &reader.read(&format!("tail.c{}", routes[0].chunk)).unwrap(),
        molecule_counts[routes[0].chunk as usize],
    )
    .unwrap();
    assert_eq!(
        selected
            .iter()
            .flat_map(|molecule| &molecule.events)
            .filter(|event| event.signal.clip_len == 10)
            .count(),
        1
    );

    // The globally deduplicated key's strongest witness is the later, unspliced MolRec. Attaching
    // its signal independently to the earlier ordinal would incorrectly make this expression true.
    let query = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args([
            "query",
            archive.to_str().unwrap(),
            "cooccur",
            "--predicate",
            "tail=terminal:chr1:9000-9001:+",
            "--predicate",
            "splice=junction:chr1:4020-8980:+",
            "--where",
            "tail & splice",
            "--universe",
            "tail",
            "--format",
            "json",
        ])
        .output()
        .unwrap();
    assert_success(&query);
    let result: serde_json::Value = serde_json::from_slice(&query.stdout).unwrap();
    assert_eq!(result["data"]["summary"]["candidate_units"], 1);
    assert_eq!(result["data"]["summary"]["selected_units"], 0);
}

#[test]
fn provenance_is_canonical_explicit_and_inspection_fails_closed_on_partial_or_tampered_v2() {
    let scratch = Scratch::new();
    let archive = ingest(&scratch.0, true);
    let mut reader = SectionReader::open(&archive).unwrap();
    let manifest_bytes = reader.read(ALIGNMENT_PROVENANCE_SECTION).unwrap();
    let manifest = AlignmentProvenanceManifest::from_json(&manifest_bytes).unwrap();
    assert_eq!(
        manifest.alignment.junction_discovery,
        JunctionDiscoveryMode::PerLibraryTwoPass
    );
    assert_eq!(
        manifest.alignment.junction_catalogue.as_ref().unwrap().role,
        JunctionCatalogueRole::PerLibraryPass1
    );
    assert_eq!(
        manifest
            .alignment
            .junction_catalogue
            .as_ref()
            .unwrap()
            .data_rows,
        2
    );
    assert_eq!(manifest.alignment.ordered_inputs.len(), 2);
    assert!(manifest.alignment.ordered_inputs[0]
        .locator
        .ends_with("R1.fastq.gz"));
    assert!(manifest.alignment.ordered_inputs[1]
        .locator
        .ends_with("R2.fastq.gz"));
    assert_eq!(manifest.alignment.chemistry.as_deref(), Some("10x-3p-v3"));
    assert_eq!(
        manifest.ingest.terminal_tail_rule.as_deref(),
        Some(terminal_tail::TERMINAL_TAIL_RULE)
    );
    assert_eq!(
        reader.read(JUNCTION_CATALOGUE_SECTION).unwrap(),
        std::fs::read(scratch.0.join("SJ.out.tab")).unwrap()
    );

    let inspect = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["inspect-archive", archive.to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert_success(&inspect);
    let inspected: serde_json::Value = serde_json::from_slice(&inspect.stdout).unwrap();
    assert_eq!(
        inspected["molecular_evidence"]["schema"],
        MOLECULAR_EVIDENCE_SCHEMA
    );
    assert_eq!(
        inspected["molecular_evidence"]["alignment_provenance_status"],
        "available"
    );
    assert_eq!(
        inspected["molecular_evidence"]["terminal_tail"]["events"],
        3
    );

    // A valid outer root does not make contradictory logical provenance true. The manifest and
    // meta copies of recoverable physical-layout fields must agree after either is rewritten.
    let layout_mismatch = scratch.0.join("layout-provenance-mismatch.aie");
    let mut layout_manifest = manifest.clone();
    layout_manifest.ingest.chunk_bp += 1;
    let mut layout_writer = SectionWriter::create(&layout_mismatch, 1).unwrap();
    for name in reader.names().map(str::to_owned).collect::<Vec<_>>() {
        let raw = if name == ALIGNMENT_PROVENANCE_SECTION {
            layout_manifest.to_canonical_json().unwrap()
        } else {
            reader.read(&name).unwrap()
        };
        layout_writer.section(&name, &raw).unwrap();
    }
    layout_writer.finish().unwrap();
    let mismatch_inspect = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args([
            "inspect-archive",
            layout_mismatch.to_str().unwrap(),
            "--json",
        ])
        .output()
        .unwrap();
    assert!(!mismatch_inspect.status.success());
    assert!(String::from_utf8_lossy(&mismatch_inspect.stderr)
        .contains("provenance ingest layout disagrees with archive metadata"));

    let genome = scratch.0.join("genome.fa");
    let mut fasta = b">chr1\n".to_vec();
    fasta.extend(std::iter::repeat_n(b'A', 10_000));
    fasta.push(b'\n');
    std::fs::write(&genome, fasta).unwrap();
    let stamped = scratch.0.join("stamped.aie");
    let stamp = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args([
            "stamp-genome",
            archive.to_str().unwrap(),
            "--genome",
            genome.to_str().unwrap(),
            "--out",
            stamped.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert_success(&stamp);
    let mut stamped_reader = SectionReader::open(&stamped).unwrap();
    let stamped_meta: serde_json::Value =
        serde_json::from_slice(&stamped_reader.read("meta").unwrap()).unwrap();
    let stamped_manifest = AlignmentProvenanceManifest::from_json(
        &stamped_reader.read(ALIGNMENT_PROVENANCE_SECTION).unwrap(),
    )
    .unwrap();
    assert_eq!(stamped_manifest, manifest);
    assert_eq!(
        stamped_reader.read(ALIGNMENT_PROVENANCE_SECTION).unwrap(),
        manifest_bytes,
        "stamping a later reference must not rewrite original alignment provenance"
    );
    let binding: GenomeReferenceBinding =
        serde_json::from_value(stamped_meta["genome_reference_binding"].clone()).unwrap();
    assert_eq!(
        binding.signature.digest,
        stamped_meta["genome_sig"]["digest"]
    );
    assert_eq!(binding.bound_by, GenomeBindingAction::StampGenome);
    assert_eq!(
        binding.relationship_status,
        ProvenanceStatus::DeclaredByCaller
    );
    assert_eq!(
        binding.identity.blake3,
        blake3::hash(&std::fs::read(&genome).unwrap())
            .to_hex()
            .to_string()
    );
    let stamped_inspect = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["inspect-archive", stamped.to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert_success(&stamped_inspect);
    let stamped_inspection: serde_json::Value =
        serde_json::from_slice(&stamped_inspect.stdout).unwrap();
    assert_eq!(
        stamped_inspection["molecular_evidence"]["genome_reference_binding_status"],
        "available"
    );
    assert_eq!(
        stamped_inspection["molecular_evidence"]["genome_reference_binding"]["bound_by"],
        "stamp-genome"
    );

    // The manifest and typed side capability are one coordinated logical revision. Neither may
    // claim that the extraction rule ran without the other, even when an attacker recomputes a
    // valid outer-container root for the inconsistent section set.
    let missing_rule = scratch.0.join("tail-without-provenance-rule.aie");
    let mut missing_rule_manifest = manifest.clone();
    missing_rule_manifest.ingest.terminal_tail_rule = None;
    let mut missing_rule_writer = SectionWriter::create(&missing_rule, 1).unwrap();
    for name in reader.names().map(str::to_owned).collect::<Vec<_>>() {
        let raw = if name == ALIGNMENT_PROVENANCE_SECTION {
            missing_rule_manifest.to_canonical_json().unwrap()
        } else {
            reader.read(&name).unwrap()
        };
        missing_rule_writer.section(&name, &raw).unwrap();
    }
    missing_rule_writer.finish().unwrap();
    let missing_rule_inspect = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["inspect-archive", missing_rule.to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert!(!missing_rule_inspect.status.success());
    assert!(String::from_utf8_lossy(&missing_rule_inspect.stderr)
        .contains("capability and provenance extraction rule disagree"));

    let core_archive = ingest(&scratch.0, false);
    let mut core_reader = SectionReader::open(&core_archive).unwrap();
    let mut false_rule_manifest = AlignmentProvenanceManifest::from_json(
        &core_reader.read(ALIGNMENT_PROVENANCE_SECTION).unwrap(),
    )
    .unwrap();
    false_rule_manifest.ingest.terminal_tail_rule = Some(terminal_tail::TERMINAL_TAIL_RULE.into());
    let false_rule = scratch.0.join("provenance-rule-without-tail.aie");
    let mut false_rule_writer = SectionWriter::create(&false_rule, 1).unwrap();
    for name in core_reader.names().map(str::to_owned).collect::<Vec<_>>() {
        let raw = if name == ALIGNMENT_PROVENANCE_SECTION {
            false_rule_manifest.to_canonical_json().unwrap()
        } else {
            core_reader.read(&name).unwrap()
        };
        false_rule_writer.section(&name, &raw).unwrap();
    }
    false_rule_writer.finish().unwrap();
    let false_rule_inspect = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["inspect-archive", false_rule.to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert!(!false_rule_inspect.status.success());
    assert!(String::from_utf8_lossy(&false_rule_inspect.stderr)
        .contains("capability and provenance extraction rule disagree"));

    let partial = scratch.0.join("partial.aie");
    let names: Vec<String> = reader.names().map(str::to_owned).collect();
    let mut writer = SectionWriter::create(&partial, 1).unwrap();
    for name in names {
        if name != ALIGNMENT_PROVENANCE_SECTION {
            writer.section(&name, &reader.read(&name).unwrap()).unwrap();
        }
    }
    writer.finish().unwrap();
    let partial_inspect = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["inspect-archive", partial.to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert!(!partial_inspect.status.success());
    assert!(String::from_utf8_lossy(&partial_inspect.stderr).contains("lacks its declared"));

    // A rooted legacy-logical archive is valid, but inspection must label unavailable evidence
    // as unavailable rather than inventing one-pass alignment or a biological zero.
    let legacy = scratch.0.join("legacy-logical.aie");
    let mut legacy_writer = SectionWriter::create(&legacy, 1).unwrap();
    let names: Vec<String> = reader.names().map(str::to_owned).collect();
    for name in &names {
        if name == ALIGNMENT_PROVENANCE_SECTION
            || name == JUNCTION_CATALOGUE_SECTION
            || name == terminal_tail::TERMINAL_TAIL_INDEX_SECTION
            || name.starts_with("tail.c")
        {
            continue;
        }
        let mut raw = reader.read(name).unwrap();
        if name == "meta" {
            let mut meta: serde_json::Map<String, serde_json::Value> =
                serde_json::from_slice(&raw).unwrap();
            meta.remove("evidence_schema");
            meta.remove("alignment_provenance");
            meta.remove("terminal_tail");
            raw = serde_json::to_vec(&meta).unwrap();
        }
        legacy_writer.section(name, &raw).unwrap();
    }
    legacy_writer.finish().unwrap();
    let legacy_inspect = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["inspect-archive", legacy.to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert_success(&legacy_inspect);
    let legacy_value: serde_json::Value = serde_json::from_slice(&legacy_inspect.stdout).unwrap();
    assert_eq!(
        legacy_value["molecular_evidence"]["alignment_provenance_status"],
        "unavailable"
    );
    assert_eq!(
        legacy_value["molecular_evidence"]["terminal_tail_status"],
        "unavailable"
    );
    assert!(legacy_value["molecular_evidence"]["terminal_tail"].is_null());

    // The tail codec alone cannot know the attached core molecule's strand. Re-root a structurally
    // valid archive with one flipped side-event flag and prove the integrated query reader rejects
    // the disagreement rather than returning detached evidence.
    let molecule_counts = chunk_molecule_counts(&mut reader);
    let routes = terminal_tail::decode_index(
        &reader
            .read(terminal_tail::TERMINAL_TAIL_INDEX_SECTION)
            .unwrap(),
        molecule_counts.len(),
    )
    .unwrap();
    let route = routes[0];
    let tail_name = format!("tail.c{}", route.chunk);
    let mut selected = terminal_tail::decode_chunk(
        &reader.read(&tail_name).unwrap(),
        molecule_counts[route.chunk as usize],
    )
    .unwrap();
    selected[0].events[0].reverse = !selected[0].events[0].reverse;
    let mismatched_tail =
        terminal_tail::encode_chunk(molecule_counts[route.chunk as usize], &selected).unwrap();
    let strand_mismatch = scratch.0.join("strand-mismatch.aie");
    let mut mismatch_writer = SectionWriter::create(&strand_mismatch, 1).unwrap();
    for name in names {
        let raw = if name == tail_name {
            mismatched_tail.clone()
        } else {
            reader.read(&name).unwrap()
        };
        mismatch_writer.section(&name, &raw).unwrap();
    }
    mismatch_writer.finish().unwrap();
    let mismatch_query = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args([
            "query",
            strand_mismatch.to_str().unwrap(),
            "cooccur",
            "--predicate",
            "tail=terminal:chr1:0-1000:+",
            "--where",
            "tail",
            "--universe",
            "tail",
            "--format",
            "json",
        ])
        .output()
        .unwrap();
    assert!(!mismatch_query.status.success());
    assert!(String::from_utf8_lossy(&mismatch_query.stderr).contains("strand disagrees"));

    let tampered = scratch.0.join("tampered.aie");
    let mut bytes = std::fs::read(&archive).unwrap();
    let root_byte = bytes.len() - 20;
    bytes[root_byte] ^= 1;
    std::fs::write(&tampered, bytes).unwrap();
    let tampered_inspect = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["inspect-archive", tampered.to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert!(!tampered_inspect.status.success());
    assert!(String::from_utf8_lossy(&tampered_inspect.stderr).contains("root mismatch"));
}

#[test]
fn collection_tail_search_routes_exact_anchors_and_reports_mixed_capability() {
    let scratch = Scratch::new();
    let tail_root = scratch.0.join("tail-source");
    let core_root = scratch.0.join("core-source");
    std::fs::create_dir(&tail_root).unwrap();
    std::fs::create_dir(&core_root).unwrap();
    let tail_archive = ingest(&tail_root, true);
    let core_archive = ingest(&core_root, false);
    let collection = scratch.0.join("tail-atlas.aicollection");
    let bin = env!("CARGO_BIN_EXE_aie");

    let build = Command::new(bin)
        .args(["collection", "build", "--sample"])
        .arg(format!("tail={}", tail_archive.display()))
        .arg("--sample")
        .arg(format!("core={}", core_archive.display()))
        .arg("--allow-unstamped")
        .arg("--out")
        .arg(&collection)
        .output()
        .unwrap();
    assert_success(&build);

    let search = Command::new(bin)
        .args(["collection", "find-events"])
        .arg(&collection)
        .args([
            "--kind",
            "terminal-tail",
            "--terminal-cluster-bp",
            "25",
            "--format",
            "json",
        ])
        .output()
        .unwrap();
    assert_success(&search);
    let result: serde_json::Value = serde_json::from_slice(&search.stdout).unwrap();
    assert_eq!(
        result["result_schema"],
        "gravlax.collection.find-events.result.v1"
    );
    let table_rows = |name: &str| -> Vec<serde_json::Value> {
        let table = result["data"]["tables"]
            .as_array()
            .unwrap()
            .iter()
            .find(|table| table["name"] == name)
            .unwrap_or_else(|| panic!("missing result table {name}"));
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
    for name in [
        "capabilities",
        "entities",
        "components",
        "counts",
        "terminal_anchors",
        "terminal_counts",
    ] {
        assert!(result["data"]["tables"]
            .as_array()
            .unwrap()
            .iter()
            .any(|table| table["name"] == name));
    }
    let tail_capability = table_rows("capabilities")
        .into_iter()
        .find(|row| row["kind"] == "terminal_tail" && row["scope"] == "aggregate")
        .unwrap();
    assert_eq!(tail_capability["status"], "partially_available");
    assert_eq!(tail_capability["archives_available"], 1);
    assert_eq!(tail_capability["archives_unavailable"], 1);
    let archive_capabilities = table_rows("capabilities")
        .into_iter()
        .filter(|row| row["kind"] == "terminal_tail" && row["scope"] == "archive")
        .collect::<Vec<_>>();
    assert_eq!(archive_capabilities.len(), 2);
    assert!(archive_capabilities.iter().any(|row| {
        row["sample"] == "tail"
            && row["status"] == "available"
            && row["included_in_denominator"] == true
    }));
    assert!(archive_capabilities.iter().any(|row| {
        row["sample"] == "core"
            && row["status"] == "unavailable"
            && row["included_in_denominator"] == false
    }));

    let entities = table_rows("entities");
    assert_eq!(entities.len(), 1);
    assert_eq!(entities[0]["kind"], "terminal_tail");
    assert_eq!(entities[0]["start"], 120);
    assert_eq!(entities[0]["end"], 141);
    assert_eq!(entities[0]["summit"], 120);
    assert_eq!(entities[0]["exact_umi_classes"], 1);
    assert_eq!(entities[0]["exact_samples"], 1);
    assert_eq!(entities[0]["exact_donors"], 1);
    assert_eq!(table_rows("counts").len(), 1);
    assert_eq!(
        table_rows("terminal_anchors")
            .iter()
            .map(|row| row["anchor"].as_u64().unwrap())
            .collect::<Vec<_>>(),
        vec![120, 130, 140]
    );
    assert_eq!(table_rows("terminal_counts").len(), 3);

    let capped = Command::new(bin)
        .args(["collection", "find-events"])
        .arg(&collection)
        .args([
            "--kind",
            "terminal-tail",
            "--max-terminal-events",
            "2",
            "--format",
            "json",
        ])
        .output()
        .unwrap();
    assert!(!capped.status.success());
    assert!(String::from_utf8_lossy(&capped.stderr).contains("exceeding --max-terminal-events 2"));

    let wrong_strand = Command::new(bin)
        .args(["collection", "find-events"])
        .arg(&collection)
        .args(["--kind", "terminal-tail", "--solo-strand", "reverse"])
        .output()
        .unwrap();
    assert!(!wrong_strand.status.success());
    assert!(
        String::from_utf8_lossy(&wrong_strand.stderr).contains("requires --solo-strand forward")
    );
}
