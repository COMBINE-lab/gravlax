use std::path::{Path, PathBuf};
use std::process::Command;

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-uniform-artifact-{}-{nonce}",
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

fn write_gtf(path: &Path) {
    std::fs::write(
        path,
        concat!(
            "chr1\ttest\tgene\t101\t250\t.\t+\t.\tgene_id \"G1.1\"; gene_name \"ONE\";\n",
            "chr1\ttest\ttranscript\t101\t250\t.\t+\t.\tgene_id \"G1.1\"; transcript_id \"T1.1\";\n",
            "chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id \"G1.1\"; transcript_id \"T1.1\"; exon_id \"E1.1\";\n",
            "chr1\ttest\texon\t201\t250\t.\t+\t.\tgene_id \"G1.1\"; transcript_id \"T1.1\"; exon_id \"E2.1\";\n",
        ),
    )
    .unwrap();
}

#[test]
fn compile_annotation_preserves_legacy_and_emits_atomic_uniform_reports() {
    let scratch = Scratch::new();
    let gtf = scratch.0.join("genes.gtf");
    write_gtf(&gtf);

    let legacy_aic = scratch.0.join("legacy.aic");
    let legacy = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("compile-annotation")
        .arg(&gtf)
        .arg("--out")
        .arg(&legacy_aic)
        .output()
        .unwrap();
    assert!(
        legacy.status.success(),
        "legacy compile failed: {}",
        String::from_utf8_lossy(&legacy.stderr)
    );
    assert!(String::from_utf8_lossy(&legacy.stdout).starts_with("compiled 1 genes, 1 transcripts"));
    assert!(legacy_aic.is_file());

    let uniform_aic = scratch.0.join("uniform.aic");
    let uniform = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("compile-annotation")
        .arg(&gtf)
        .arg("--out")
        .arg(&uniform_aic)
        .args(["--report-format", "json"])
        .output()
        .unwrap();
    assert!(
        uniform.status.success(),
        "uniform compile failed: {}",
        String::from_utf8_lossy(&uniform.stderr)
    );
    assert_eq!(
        std::fs::read(&legacy_aic).unwrap(),
        std::fs::read(&uniform_aic).unwrap()
    );
    let result: serde_json::Value = serde_json::from_slice(&uniform.stdout).unwrap();
    assert_eq!(result["$schema"], "gravlax.result-envelope.v1");
    assert_eq!(
        result["result_schema"],
        "gravlax.annotation.compile.result.v1"
    );
    assert_eq!(result["data"]["summary"]["genes"], 1);
    assert_eq!(result["data"]["summary"]["transcripts"], 1);
    assert_eq!(result["data"]["summary"]["exons"], 2);
    assert_eq!(
        result["provenance"]["annotation_digest"],
        format!(
            "blake3:{}",
            blake3::hash(&std::fs::read(&gtf).unwrap()).to_hex()
        )
    );
    assert_eq!(result["data"]["tables"][0]["name"], "artifacts");
    assert_eq!(
        result["data"]["tables"][0]["schema"]["semantics"]["row_semantics"],
        "set"
    );
    assert_eq!(
        result["data"]["tables"][0]["rows"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert!(result["data"]["tables"][0]["rows"][0][3]
        .as_str()
        .unwrap()
        .starts_with("aic-v2-payload-blake3:"));

    let report = scratch.0.join("compile.tsv");
    let file_aic = scratch.0.join("file.aic");
    let file_output = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("compile-annotation")
        .arg(&gtf)
        .arg("--out")
        .arg(&file_aic)
        .args(["--report-format", "tsv", "--report-output"])
        .arg(&report)
        .output()
        .unwrap();
    assert!(file_output.status.success());
    assert!(file_output.stdout.is_empty());
    assert!(std::fs::read_to_string(&report)
        .unwrap()
        .contains("artifact_kind\tpath\tbytes"));

    let blocked_aic = scratch.0.join("blocked.aic");
    let blocked = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("compile-annotation")
        .arg(&gtf)
        .arg("--out")
        .arg(&blocked_aic)
        .args(["--report-format", "json", "--report-output"])
        .arg(&report)
        .output()
        .unwrap();
    assert!(!blocked.status.success());
    assert!(!blocked_aic.exists());
    assert!(String::from_utf8_lossy(&blocked.stderr).contains("refusing to overwrite"));

    let alias_dir = scratch.0.join("alias-component");
    std::fs::create_dir(&alias_dir).unwrap();
    let aliased_aic = scratch.0.join("aliased.aic");
    let aliased_report = alias_dir.join("..").join("aliased.aic");
    let aliased = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("compile-annotation")
        .arg(&gtf)
        .arg("--out")
        .arg(&aliased_aic)
        .args(["--report-format", "json", "--report-output"])
        .arg(&aliased_report)
        .output()
        .unwrap();
    assert!(!aliased.status.success());
    assert!(String::from_utf8_lossy(&aliased.stderr)
        .contains("compiled annotation output and report output must differ"));
    assert!(!aliased_aic.exists());
}
