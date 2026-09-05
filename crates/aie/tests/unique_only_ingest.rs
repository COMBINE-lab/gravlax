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
            "gravlax-unique-only-ingest-{}-{nonce}",
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
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn write_unique_spliced_bam(path: &Path) {
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(1_000).unwrap()),
        )
        .build();
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
        (Tag::new(b'U', b'R'), Value::from("ACGTACGTACGT")),
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
    ]
    .into_iter()
    .collect();
    let record = RecordBuf::builder()
        .set_name("unique-spliced")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from(101).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; 50]))
        .set_quality_scores(QualityScores::from(vec![30; 50]))
        .set_data(data)
        .build();

    let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    writer.write_alignment_record(&header, &record).unwrap();
    writer.try_finish().unwrap();
}

#[test]
fn unique_only_bam_ingests_queries_and_verifies() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("unique-only.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    let archive = scratch.0.join("unique-only.aie");
    write_unique_spliced_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();

    let bin = env!("CARGO_BIN_EXE_aie");
    run({
        let mut command = Command::new(bin);
        command
            .arg("ingest-archive")
            .arg(&bam)
            .arg("--whitelist")
            .arg(&whitelist)
            .arg("--out")
            .arg(&archive)
            .args(["--zstd-level", "1", "--chunk-mb", "1"]);
        command
    });

    let junction: serde_json::Value = serde_json::from_slice(
        &run({
            let mut command = Command::new(bin);
            command
                .arg("query")
                .arg(&archive)
                .arg("junction")
                .arg("chr1:125-225")
                .args(["--top", "0", "--json"]);
            command
        })
        .stdout,
    )
    .unwrap();
    assert_eq!(junction["donor"], 125);
    assert_eq!(junction["acceptor"], 225);
    assert_eq!(junction["umis"], 1);
    assert_eq!(junction["cells"], 1);

    let diagnosis: serde_json::Value = serde_json::from_slice(
        &run({
            let mut command = Command::new(bin);
            command
                .arg("doctor")
                .arg(&archive)
                .args(["--verify-content", "--json"]);
            command
        })
        .stdout,
    )
    .unwrap();
    assert_eq!(diagnosis["ok"], true);
    let archive_check = diagnosis["checks"]
        .as_array()
        .unwrap()
        .iter()
        .find(|check| {
            check["id"]
                .as_str()
                .is_some_and(|id| id.starts_with("archive:"))
        })
        .unwrap();
    assert_eq!(archive_check["status"], "pass");
    assert_eq!(archive_check["data"]["all_payloads_verified"], true);
    assert_eq!(archive_check["data"]["semantic_content_verified"], true);
    assert_eq!(archive_check["data"]["decoded_molecules"], 1);
}

#[test]
fn optional_access_index_matches_fallback_without_full_scan_permission() {
    let scratch = Scratch::new();
    let bam = scratch.0.join("input.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    write_unique_spliced_bam(&bam);
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    let mut results = Vec::new();
    for indexed in [false, true] {
        let archive = scratch.0.join(format!("{indexed}.aie"));
        run({
            let mut c = Command::new(env!("CARGO_BIN_EXE_aie"));
            c.arg("ingest-archive")
                .arg(&bam)
                .arg("--whitelist")
                .arg(&whitelist)
                .arg("--out")
                .arg(&archive)
                .args(["--zstd-level", "1"]);
            if indexed {
                c.arg("--access-index");
            }
            c
        });
        let value: serde_json::Value = serde_json::from_slice(
            &run({
                let mut c = Command::new(env!("CARGO_BIN_EXE_aie"));
                c.arg("query").arg(&archive).args([
                    "cooccur",
                    "--predicate",
                    "a=region:chr1:101-110",
                    "--predicate",
                    "j=junction:chr1:125-225",
                    "--where",
                    "a & j",
                    "--universe",
                    "a",
                    "--unit",
                    "umi-class",
                    "--region-match",
                    "aligned-block",
                    "--format",
                    "json",
                    "--emit-membership",
                ]);
                if !indexed {
                    c.arg("--allow-full-scan");
                }
                c
            })
            .stdout,
        )
        .unwrap();
        results.push(value["data"]["tables"].clone());
        run({
            let mut c = Command::new(env!("CARGO_BIN_EXE_aie"));
            c.arg("doctor")
                .arg(&archive)
                .args(["--verify-content", "--json"]);
            c
        });
    }
    assert!(results[0].is_array());
    assert_eq!(results[0], results[1]);
}

#[test]
fn fidelity_recovers_omitted_middle_geometry_without_joining_loci() {
    let scratch = Scratch::new();
    let input = scratch.0.join("geometry.bam");
    let whitelist = scratch.0.join("whitelist.txt");
    std::fs::write(&whitelist, format!("{BARCODE}\n")).unwrap();
    let header = sam::Header::builder()
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZero::new(10_000).unwrap()),
        )
        .build();
    let mut writer = bam::io::Writer::new(std::fs::File::create(&input).unwrap());
    writer.write_header(&header).unwrap();
    // The middle placement is absent from the legacy extremes, but witnessed twice in BAM.
    // The final placement shares the UMI yet belongs to a separate locus.
    for (i, pos) in [101, 121, 121, 141, 5001].into_iter().enumerate() {
        let data: Data = [
            (Tag::new(b'C', b'R'), Value::from(BARCODE)),
            (Tag::new(b'C', b'Y'), Value::from("IIIIIIIIIIIIIIII")),
            (Tag::new(b'U', b'R'), Value::from("ACGTACGTACGT")),
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(1u8)),
        ]
        .into_iter()
        .collect();
        let record = RecordBuf::builder()
            .set_name(format!("geometry-{i}"))
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(pos).unwrap())
            .set_cigar([Op::new(Kind::Match, 20)].into_iter().collect::<Cigar>())
            .set_sequence(Sequence::from(vec![b'A'; 20]))
            .set_quality_scores(QualityScores::from(vec![30; 20]))
            .set_data(data)
            .build();
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    let gtf = scratch.0.join("geometry.gtf");
    std::fs::write(&gtf, "chr1\ttest\texon\t101\t120\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\nchr1\ttest\texon\t141\t160\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n").unwrap();
    for (fidelity, compression) in [(false, false), (true, false), (true, true)] {
        let archive = scratch
            .0
            .join(format!("geometry-{fidelity}-{compression}.aie"));
        run({
            let mut c = Command::new(env!("CARGO_BIN_EXE_aie"));
            c.arg("ingest-archive")
                .arg(&input)
                .arg("--whitelist")
                .arg(&whitelist)
                .arg("--out")
                .arg(&archive)
                .args([
                    "--zstd-level",
                    "1",
                    "--access-index",
                    "--chunk-records",
                    "1",
                ]);
            if fidelity {
                c.arg("--geometry-fidelity");
            }
            if compression {
                c.arg("--compression-tuning");
            }
            c
        });
        let query: serde_json::Value = serde_json::from_slice(
            &run({
                let mut c = Command::new(env!("CARGO_BIN_EXE_aie"));
                c.arg("query").arg(&archive).args([
                    "cooccur",
                    "--predicate",
                    "a=region:chr1:100-160",
                    "--predicate",
                    "b=region:chr1:125-130",
                    "--where",
                    "b",
                    "--universe",
                    "a",
                    "--region-match",
                    "aligned-block",
                    "--format",
                    "json",
                ]);
                c
            })
            .stdout,
        )
        .unwrap();
        assert_eq!(query["data"]["summary"]["candidate_units"], 1);
        assert_eq!(
            query["data"]["summary"]["selected_units"],
            if fidelity { 1 } else { 0 }
        );
        assert_eq!(
            query["data"]["summary"]["indeterminate_units"],
            if fidelity { 0 } else { 1 }
        );
        for velocity in [false,true] {
            let mut outputs = Vec::new();
            for bam_reference in [true,false] {
                let out = scratch.0.join(format!("replay-{fidelity}-{compression}-{velocity}-{bam_reference}"));
                run({ let mut c = Command::new(env!("CARGO_BIN_EXE_aie"));
                    c.arg("replay-rows").arg(if bam_reference {&input} else {&archive})
                        .arg("--gtf").arg(&gtf).arg("--barcodes").arg(&whitelist).arg("--out-dir").arg(&out);
                    if bam_reference { c.arg("--from-bam").arg("--whitelist").arg(&whitelist);
                        if fidelity { c.arg("--geometry-fidelity"); }
                    }
                    if velocity { c.arg("--velocity"); } c });
                outputs.push(out);
            }
            let matrices: &[&str] = if velocity { &["spliced.mtx","unspliced.mtx","ambiguous.mtx"] } else { &["matrix.mtx"] };
            for name in matrices.iter().chain(["features.tsv","barcodes.tsv"].iter()) {
                assert_eq!(std::fs::read(outputs[0].join(name)).unwrap(), std::fs::read(outputs[1].join(name)).unwrap(), "BAM/archive mismatch: {name}");
            }
        }
        run({
            let mut c = Command::new(env!("CARGO_BIN_EXE_aie"));
            c.arg("doctor")
                .arg(&archive)
                .args(["--verify-content", "--json"]);
            c
        });
    }
}
