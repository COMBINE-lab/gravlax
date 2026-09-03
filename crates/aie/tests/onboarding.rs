use std::path::PathBuf;
use std::process::Command;

struct Scratch(PathBuf);

impl Scratch {
    fn project() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-onboarding-test-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&path).unwrap();
        std::fs::write(
            path.join("aie-project.yaml"),
            "schema_version: 1\nname: onboarding-test\nresources: {}\n",
        )
        .unwrap();
        Self(path)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

#[test]
fn doctor_json_reports_the_loaded_project() {
    let project = Scratch::project();
    let output = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("doctor")
        .arg("--project")
        .arg(&project.0)
        .arg("--json")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "doctor failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schema"], "gravlax.doctor.v1");
    assert_eq!(value["ok"], true);
    let project_check = value["checks"]
        .as_array()
        .unwrap()
        .iter()
        .find(|check| check["id"] == "project")
        .unwrap();
    assert_eq!(project_check["status"], "pass");
    assert_eq!(project_check["data"]["name"], "onboarding-test");
}

#[test]
fn doctor_failure_keeps_a_parseable_json_report() {
    let project = Scratch::project();
    let output = Command::new(env!("CARGO_BIN_EXE_aie"))
        .arg("doctor")
        .arg(project.0.join("missing.aie"))
        .arg("--project")
        .arg(&project.0)
        .arg("--json")
        .output()
        .unwrap();
    assert!(!output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["ok"], false);
    assert_eq!(value["summary"]["failures"], 1);
}

#[test]
fn every_documented_completion_shell_generates_a_script() {
    for (shell, marker) in [
        ("bash", "_aie"),
        ("zsh", "#compdef aie"),
        ("fish", "complete -c aie"),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_aie"))
            .args(["completions", shell])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{shell} completion generation failed"
        );
        assert!(
            String::from_utf8_lossy(&output.stdout).contains(marker),
            "{shell} completion script lacks {marker}"
        );
    }
}

#[test]
fn explorer_help_discloses_the_loopback_port_interface() {
    let output = Command::new(env!("CARGO_BIN_EXE_aie"))
        .args(["explore", "--help"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let help = String::from_utf8_lossy(&output.stdout);
    assert!(help.contains("Loopback TCP port"));
    assert!(help.contains("--project"));
    assert!(!help.contains("--host"));
}
