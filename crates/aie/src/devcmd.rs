//! Research and validation instruments kept out of the everyday command surface.

use anyhow::Result;
use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Compare per-read evidence between two mapping configurations.
    GateC(crate::gatec::Args),
    /// Bound annotation-independent UMI-grouping error against STARsolo.
    GateD(crate::gated::Args),
    /// Measure evidence-level archive size and stream composition.
    Build(crate::build::Args),
    /// Measure UMI adjacency structure and replay preservation.
    UmiGraph(crate::graph::Args),
    /// Run the original BAM-direct replay diagnostic.
    Replay(crate::replaycmd::Args),
    /// Compare STARsolo GX assignment with the local assignment port.
    AssignDiff(crate::assigndiff::Args),
    /// Measure signature, paralog-placement, and stream coding statistics.
    SigStats(crate::sigstats::Args),
    /// Audit encoded archive sections and candidate coding factorizations.
    Debug(crate::debugcmd::Args),
    /// Run multimapper-recovery modes and masked-evidence scoring.
    Em(crate::archivecmd::EmArgs),
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::GateC(args) => crate::gatec::run(args),
        Command::GateD(args) => crate::gated::run(args),
        Command::Build(args) => crate::build::run(args),
        Command::UmiGraph(args) => crate::graph::run(args),
        Command::Replay(args) => crate::replaycmd::run(args),
        Command::AssignDiff(args) => crate::assigndiff::run(args),
        Command::SigStats(args) => crate::sigstats::run(args),
        Command::Debug(args) => crate::debugcmd::run(args),
        Command::Em(args) => crate::archivecmd::run_em(args),
    }
}

pub fn warn_deprecated_alias(alias: &str) {
    eprintln!("warning: `aie {alias}` is deprecated; use `aie dev {alias}`");
}
