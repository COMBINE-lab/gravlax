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
