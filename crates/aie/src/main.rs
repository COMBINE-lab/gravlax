//! `aie` — annotation-independent evidence toolkit.

mod apastats;
mod archivecmd;
mod assigndiff;
mod build;
mod collectioncmd;
mod comparecmd;
mod completion;
mod debugcmd;
mod devcmd;
mod doctor;
mod endcmd;
mod explorer;
mod extendcmd;
mod gatec;
mod gated;
mod graph;
mod ingestcmd;
mod moleculebam;
mod plancmd;
mod plots;
mod projectcmd;
mod querycmd;
mod replaycmd;
mod resolvecmd;
mod rows;
mod shaperoute;
mod sigstats;
mod sparseout;
mod viz;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use gravlax_output::{
    canonical_destination_key, publish_file_no_clobber, reported_output_path, DataType,
    Durability, Field, OutputError, OutputFormat, Producer, Provenance, ResultContext,
    RowSemantics, SelectionSummary, StreamingBundleWriter, TableSchema, TableSemantics,
};
use serde::Serialize;
use serde_json::json;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "aie", version, about = "Annotation-independent molecular evidence for scRNA-seq")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Create and manage a portable analysis project with named input resources.
    Project(projectcmd::Args),
    /// Validate or run a versioned YAML/JSON analysis plan.
    Plan(plancmd::Args),
    /// Diagnose this installation, project, and selected Gravlax artifacts.
    Doctor(doctor::Args),
    /// Browse project plans and exact result artifacts in a local, read-only web UI.
    Explore(explorer::Args),
    /// Generate a shell completion script from the current command interface.
    Completions(completion::Args),
    /// Resolve gene symbols and stable feature identifiers against an explicit annotation.
    Resolve(resolvecmd::Args),
    /// Check ingest inputs or print a chemistry-specific annotation-free STAR recipe.
    Ingest(ingestcmd::Args),
    /// Research, validation, and archive-format development instruments.
    Dev(devcmd::Args),
    /// Compare per-read alignment evidence between two mapping configurations.
    ///
    /// Reports the fraction of reads whose annotation-independent evidence tuple is identical,
    /// plus a taxonomy of how the rest differ. Both BAMs must come from the same reads.
    #[command(hide = true)]
    GateC(gatec::Args),
    /// Compare annotation-independent UMI grouping with STARsolo's per-(cell,gene) collapse on
    /// the same alignments.
    #[command(hide = true)]
    GateD(gated::Args),
    /// Build the annotation-independent molecule archive from an annotation-free alignment and
    /// report where its bytes go, at evidence levels E0 and E1, with and without the UMI stream.
    #[command(hide = true)]
    Build(build::Args),
    /// Measure UMI collision structure, compare collapse with STARsolo output, and report storage
    /// cost.
    #[command(hide = true)]
    UmiGraph(graph::Args),
    /// Replay an annotation against annotation-free molecules extracted from fixed alignments and
    /// compare the result with a fresh STARsolo run.
    #[command(hide = true)]
    Replay(replaycmd::Args),
    /// Read-level confusion between STARsolo's GX assignment and our alignToTranscript port, on
    /// identical alignments — pinpoints which assignment rule differs.
    #[command(hide = true)]
    AssignDiff(assigndiff::Args),
    /// Measure signature multiplicity, paralog-placement patterns, and candidate position/cell
    /// stream encodings.
    #[command(hide = true)]
    SigStats(sigstats::Args),
    /// Build an authenticated .aie v2 archive from a tagged, coordinate-sorted ingest BAM.
    IngestArchive(archivecmd::IngestArgs),
    /// Compile a GTF into a deterministic, guarded annotation artifact for fast reuse.
    CompileAnnotation(CompileAnnotationArgs),
    /// Export the exact post-correction molecule abstraction as sequence-free BAM with explicit
    /// local tags, for interchange and a capability-matched BAM/CRAM storage baseline.
    ExportMoleculeBam(moleculebam::ExportArgs),
    /// Quantify a compatible GTF from an .aie archive (or from a BAM via --from-bam, the regression
    /// reference: both paths share one row abstraction and must produce identical matrices).
    ReplayRows(archivecmd::ReplayRowsArgs),
    /// Compare two annotations on one evidence archive and explain their exact count consequences.
    CompareAnnotations(comparecmd::Args),
    /// Indexed queries against an .aie archive: regions, junction enumeration, APA, and discovery.
    Query(querycmd::Args),
    /// Coding-stack measurements on an archive: per-section/stream byte accounting, per-stream
    /// value-entropy bounds, and sharing ratios for candidate dictionary factorizations.
    #[command(hide = true)]
    Debug(debugcmd::Args),
    /// One junction query across N archives — the miniature-atlas access pattern.
    Federate(querycmd::FederateArgs),
    /// Event-centric queries across named archives with exact per-sample/group reductions.
    Cohort(querycmd::CohortArgs),
    /// Build and query a guarded, derived index over a collection of independent .aie archives.
    Collection(collectioncmd::Args),
    /// EM multimapper-recovery experiment: masked-evidence recovery scoring for per-cell, pooled,
    /// and blended sharing against a uniform baseline.
    #[command(hide = true)]
    Em(archivecmd::EmArgs),
    /// Extend per-gene 3' annotation boundaries from archived evidence (the Pool et al.
    /// reference-optimization workflow): discover -> extend -> replay, without realignment.
    Extend(extendcmd::Args),
    /// Stamp (or re-stamp) the reference genome's blake3 signature into an archive so
    /// sequence-consulting analyses can verify they see the genome the reads were aligned to.
    StampGenome(archivecmd::StampGenomeArgs),
    /// Atomically copy a legacy v1 archive into the authenticated v2 container without
    /// recompressing or changing any encoded section payload.
    SealArchive(archivecmd::SealArchiveArgs),
    /// Report an archive's native and scheme-independent encoded identities, optionally verifying
    /// every compressed payload.
    InspectArchive(archivecmd::InspectArchiveArgs),
}

#[derive(clap::Args)]
struct CompileAnnotationArgs {
    /// Source annotation in uncompressed GTF format.
    gtf: PathBuf,
    /// Destination compiled annotation (`.aic`). Refuses to overwrite an existing file.
    #[arg(long)]
    out: PathBuf,
    /// Emit a versioned report while preserving the compiled artifact at --out.
    #[arg(long, value_enum)]
    report_format: Option<CompileAnnotationReportFormat>,
    /// Publish the uniform report atomically without replacing an existing file.
    #[arg(long, requires = "report_format")]
    report_output: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum CompileAnnotationReportFormat {
    Text,
    Tsv,
    Json,
}

impl From<CompileAnnotationReportFormat> for OutputFormat {
    fn from(value: CompileAnnotationReportFormat) -> Self {
        match value {
            CompileAnnotationReportFormat::Text => Self::Text,
            CompileAnnotationReportFormat::Tsv => Self::Tsv,
            CompileAnnotationReportFormat::Json => Self::Json,
        }
    }
}

#[derive(Serialize)]
struct CompileAnnotationSummary {
    genes: u64,
    transcripts: u64,
    exons: u64,
    output_bytes: u64,
}

fn compile_annotation_artifact_schema() -> Result<TableSchema, OutputError> {
    TableSchema::new(
        "gravlax.annotation.compile.artifacts.v1",
        vec![
            Field::new("artifact_kind", DataType::String),
            Field::new("path", DataType::String),
            Field::new("bytes", DataType::UInt64),
            Field::new("identity", DataType::String),
        ],
    )?
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["artifact_kind", "path"]))
}

fn write_compile_annotation_report<W: Write>(
    writer: W,
    format: OutputFormat,
    input: &str,
    input_digest: &str,
    output: &str,
    output_identity: &str,
    summary: &CompileAnnotationSummary,
) -> Result<(), OutputError> {
    let mut parameters = BTreeMap::new();
    parameters.insert("input_annotation".into(), json!(input));
    let context = ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            annotation: Some(input.to_owned()),
            annotation_digest: Some(input_digest.to_owned()),
            parameters,
            ..Default::default()
        },
        warnings: Vec::new(),
    };
    let schema = compile_annotation_artifact_schema()?;
    let selection = SelectionSummary::complete(1);
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        "gravlax.annotation.compile.result.v1",
        format,
        &context,
        summary,
    )?;
    bundle.write_table("artifacts", &schema, Some(&selection), |rows| {
        rows.write_row_with(|row| {
            row.string("compiled-annotation")?;
            row.string(output)?;
            row.uint64(summary.output_bytes)?;
            row.string(output_identity)?;
            Ok(())
        })
    })?;
    bundle.finish()?;
    Ok(())
}

fn run_compile_annotation(args: CompileAnnotationArgs) -> Result<()> {
    if args.out.exists() {
        bail!("refusing to overwrite {}", args.out.display());
    }
    if let Some(path) = &args.report_output {
        if std::fs::symlink_metadata(path).is_ok() {
            bail!("refusing to overwrite {}", path.display());
        }
    }
    if args.gtf == args.out {
        bail!("compiled annotation output must differ from its input");
    }
    if let Some(report) = args.report_output.as_deref() {
        if canonical_destination_key(report)? == canonical_destination_key(&args.out)? {
            bail!("compiled annotation output and report output must differ");
        }
    }
    let uniform_paths = args
        .report_format
        .map(|_| {
            Ok::<_, anyhow::Error>((
                args.gtf
                    .to_str()
                    .context("uniform report requires a UTF-8 input annotation path")?
                    .to_owned(),
                reported_output_path(&args.out)?,
            ))
        })
        .transpose()?;
    let t0 = std::time::Instant::now();
    let (annotation, input_digest) = if args.report_format.is_some() {
        let (annotation, digest) = anno::Annotation::from_gtf_with_digest(&args.gtf)?;
        (annotation, Some(digest))
    } else {
        (anno::Annotation::from_gtf(&args.gtf)?, None)
    };
    let parsed = t0.elapsed();
    let mut compiled_artifact = None;
    let publication = publish_file_no_clobber(&args.out, Durability::Flush, |writer| {
        let artifact = annotation
            .write_compiled_to(writer)
            .map_err(|error| OutputError::Sink(format!("writing compiled annotation: {error:#}")))?;
        compiled_artifact = Some(artifact);
        Ok(())
    })?;
    for warning in publication.warnings {
        eprintln!("warning: {warning}");
    }
    let (produced_identity, bytes) =
        compiled_artifact.expect("successful publication wrote a compiled annotation");
    let output_identity = args.report_format.map(|_| produced_identity);
    let exons: usize = annotation
        .transcripts
        .iter()
        .map(|transcript| transcript.exons.len())
        .sum();
    if let (Some(format), Some((input, output))) = (args.report_format, uniform_paths) {
        let summary = CompileAnnotationSummary {
            genes: annotation.gene_ids.len() as u64,
            transcripts: annotation.transcripts.len() as u64,
            exons: exons as u64,
            output_bytes: bytes,
        };
        match &args.report_output {
            Some(path) => {
                let outcome = publish_file_no_clobber(path, Durability::Flush, |writer| {
                    write_compile_annotation_report(
                        writer,
                        format.into(),
                        &input,
                        input_digest.as_deref().expect("uniform input identity"),
                        &output,
                        output_identity.as_deref().expect("uniform output identity"),
                        &summary,
                    )
                })?;
                for warning in outcome.warnings {
                    eprintln!("warning: {warning}");
                }
            }
            None => write_compile_annotation_report(
                std::io::stdout().lock(),
                format.into(),
                &input,
                input_digest.as_deref().expect("uniform input identity"),
                &output,
                output_identity.as_deref().expect("uniform output identity"),
                &summary,
            )?,
        }
        eprintln!(
            "compiled annotation in {:.2}s (parse {:.2}s)",
            t0.elapsed().as_secs_f32(),
            parsed.as_secs_f32(),
        );
    } else {
        println!(
            "compiled {} genes, {} transcripts, {exons} exons into {} bytes at {} (parse {:.2}s, total {:.2}s)",
            annotation.gene_ids.len(),
            annotation.transcripts.len(),
            bytes,
            args.out.display(),
            parsed.as_secs_f32(),
            t0.elapsed().as_secs_f32()
        );
    }
    Ok(())
}

fn main() -> Result<()> {
    // Standard Unix pipeline behavior: die silently when the reader closes the pipe
    // (`aie query … | head` used to panic with a broken-pipe message after head exited).
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
    // Shared interactive host: default the rayon pool to 24 threads (the project's standing
    // thread-count discipline) unless the user overrides via RAYON_NUM_THREADS.
    if std::env::var_os("RAYON_NUM_THREADS").is_none() {
        let n = std::thread::available_parallelism().map(|p| p.get()).unwrap_or(8).min(24);
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }
    match Cli::parse().cmd {
        Cmd::Project(a) => projectcmd::run(a),
        Cmd::Plan(a) => plancmd::run(a),
        Cmd::Doctor(a) => doctor::run(a),
        Cmd::Explore(a) => explorer::run(a),
        Cmd::Completions(a) => completion::run::<Cli>(a),
        Cmd::Resolve(a) => resolvecmd::run(a),
        Cmd::Ingest(a) => ingestcmd::run(a),
        Cmd::Dev(a) => devcmd::run(a),
        Cmd::GateC(a) => {
            devcmd::warn_deprecated_alias("gate-c");
            gatec::run(a)
        }
        Cmd::GateD(a) => {
            devcmd::warn_deprecated_alias("gate-d");
            gated::run(a)
        }
        Cmd::Build(a) => {
            devcmd::warn_deprecated_alias("build");
            build::run(a)
        }
        Cmd::UmiGraph(a) => {
            devcmd::warn_deprecated_alias("umi-graph");
            graph::run(a)
        }
        Cmd::Replay(a) => {
            devcmd::warn_deprecated_alias("replay");
            replaycmd::run(a)
        }
        Cmd::AssignDiff(a) => {
            devcmd::warn_deprecated_alias("assign-diff");
            assigndiff::run(a)
        }
        Cmd::SigStats(a) => {
            devcmd::warn_deprecated_alias("sig-stats");
            sigstats::run(a)
        }
        Cmd::IngestArchive(a) => archivecmd::run_ingest(a),
        Cmd::CompileAnnotation(a) => run_compile_annotation(a),
        Cmd::ExportMoleculeBam(a) => moleculebam::run_export(a),
        Cmd::ReplayRows(a) => archivecmd::run_replay_rows(a),
        Cmd::CompareAnnotations(a) => comparecmd::run(a),
        Cmd::Query(a) => querycmd::run(a),
        Cmd::Debug(a) => {
            devcmd::warn_deprecated_alias("debug");
            debugcmd::run(a)
        }
        Cmd::Federate(a) => querycmd::run_federate(a),
        Cmd::Cohort(a) => querycmd::run_cohort(a),
        Cmd::Collection(a) => collectioncmd::run(a),
        Cmd::Em(a) => {
            devcmd::warn_deprecated_alias("em");
            archivecmd::run_em(a)
        }
        Cmd::Extend(a) => extendcmd::run(a),
        Cmd::StampGenome(a) => archivecmd::run_stamp_genome(a),
        Cmd::SealArchive(a) => archivecmd::run_seal_archive(a),
        Cmd::InspectArchive(a) => archivecmd::run_inspect_archive(a),
    }
}
