//! GQ textual frontend and native evidence execution.
mod ast;
mod check;
mod eval;
mod execute;
mod model;
mod native;
mod parser;
mod physical;
mod resources;
#[cfg(test)]
mod tests;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use std::{collections::BTreeMap, path::PathBuf};

#[derive(Parser)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    Validate(Input),
    Explain(Input),
    Run(Input),
}
#[derive(Clone, clap::Args)]
struct Input {
    query: PathBuf,
    /// Bind a source alias to an authenticated archive or a federation manifest.
    #[arg(long)]
    bind: Vec<String>,
    #[arg(long)]
    project: Option<PathBuf>,
    /// Typed sample/cell metadata JSON; missing columns must be declared optional.
    #[arg(long)]
    metadata: Option<PathBuf>,
    #[arg(long)]
    allow_full_scan: bool,
    #[arg(long, default_value_t = 4096)]
    max_chunks: usize,
    #[arg(long, default_value_t = 1_000_000)]
    max_records: usize,
    #[arg(long, default_value_t = 100_000_000)]
    max_steps: u64,
    #[arg(long, default_value_t = 100_000)]
    max_rows: usize,
    #[arg(long, default_value_t = 1_000_000)]
    max_terminal_events: u64,
    #[arg(long)]
    output: Option<PathBuf>,
    /// Execution strategy; reference retains the general interpreter for comparison.
    #[arg(long, value_enum, default_value_t = Engine::Auto)]
    engine: Engine,
    /// Opt in to bounded parallel chunk decoding; may use more temporary memory.
    #[arg(long)]
    parallel_decode: bool,
    /// Report execution phase timings (parsing is excluded) to stderr.
    #[arg(long)]
    profile: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum Engine {
    Auto,
    Reference,
}
pub fn run(args: Args) -> Result<()> {
    let input = match &args.command {
        Command::Validate(i) | Command::Explain(i) | Command::Run(i) => i,
    };
    let source = std::fs::read_to_string(&input.query)
        .with_context(|| format!("reading {}", input.query.display()))?;
    let doc =
        parser::parse(&source).with_context(|| format!("parsing {}", input.query.display()))?;
    let post_parse = std::time::Instant::now();
    let assembly = doc
        .header
        .iter()
        .find_map(|(n, e)| {
            if n == "assembly" {
                if let ast::Kind::String(s) = &e.kind {
                    Some(s.as_str())
                } else {
                    None
                }
            } else {
                None
            }
        })
        .context("header assembly must be an explicit string")?;
    let resources = resources::Resources::open(
        &input.bind,
        input.metadata.as_deref(),
        input.project.as_deref(),
        assembly,
    )?;
    let plan = check::Compiler::new(&doc, resources.fields()?)?
        .with_resources(&resources)
        .compile()?;
    if input.max_chunks == 0
        || input.max_records == 0
        || input.max_rows == 0
        || input.max_steps == 0
    {
        bail!("execution budgets must be positive");
    }
    let mut bytes = Vec::new();
    let preparation_ms = post_parse.elapsed().as_secs_f64() * 1000.0;
    let execution_start = std::time::Instant::now();
    let mut profile =
        serde_json::json!({"preparation_ms": preparation_ms,"allocator":crate::allocator::name()});
    match &args.command {
        Command::Validate(_) => {
            // Unbound validation is syntax/type only. Bound validation additionally verifies
            // resource capabilities and physical-plan feasibility without decoding records.
            let bound = resources.sources.contains_key(&plan.source);
            if bound {
                execute::explain(&plan, &resources, input)?;
            }
            serde_json::to_writer_pretty(
                &mut bytes,
                &serde_json::json!({"valid":true,"language_version":1,"unit":plan.unit,"resources_checked":bound}),
            )?;
        }
        Command::Explain(_) => {
            serde_json::to_writer_pretty(&mut bytes, &execute::explain(&plan, &resources, input)?)?
        }
        Command::Run(_) => {
            let table = execute::execute(&plan, &resources, input)?;
            profile["execution_ms"] =
                serde_json::json!(execution_start.elapsed().as_secs_f64() * 1000.0);
            profile["phases_ms"] = table.summary["timings_ms"].clone();
            let serialization_start = std::time::Instant::now();
            let mut parameters = BTreeMap::new();
            parameters.insert(
                "query_digest".into(),
                serde_json::json!(blake3::hash(source.as_bytes()).to_hex().to_string()),
            );
            parameters.insert(
                "logical_plan_digest".into(),
                serde_json::json!(logical_digest(&plan)?),
            );
            parameters.insert(
                "resource_digests".into(),
                serde_json::json!(resources.digests),
            );
            parameters.insert("assembly".into(), serde_json::json!(plan.assembly));
            parameters.insert("library_rule".into(), serde_json::json!(plan.library));
            parameters.insert("execution".into(), table.summary.clone());
            let context = gravlax_output::ResultContext {
                provenance: gravlax_output::Provenance {
                    parameters,
                    ..Default::default()
                },
                ..Default::default()
            };
            write_result(&mut bytes, &table, &context)?;
            profile["serialization_ms"] =
                serde_json::json!(serialization_start.elapsed().as_secs_f64() * 1000.0);
        }
    };
    profile["post_parse_compute_ms"] =
        serde_json::json!(post_parse.elapsed().as_secs_f64() * 1000.0);
    let publication_start = std::time::Instant::now();
    use std::io::Write;
    if let Some(path) = &input.output {
        let outcome = gravlax_output::publish_file_no_clobber(
            path,
            gravlax_output::Durability::Flush,
            |writer| {
                writer.write_all(&bytes)?;
                Ok(())
            },
        )?;
        for warning in outcome.warnings {
            eprintln!("warning: {warning}");
        }
    } else {
        std::io::stdout().lock().write_all(&bytes)?;
    }
    if input.profile {
        profile["publication_ms"] =
            serde_json::json!(publication_start.elapsed().as_secs_f64() * 1000.0);
        eprintln!("gq_profile={profile}");
    }
    Ok(())
}

fn write_result(
    writer: impl std::io::Write,
    table: &execute::Table,
    context: &gravlax_output::ResultContext,
) -> Result<()> {
    use gravlax_output::*;
    use model::{Type, Value};
    let fields = table
        .columns
        .iter()
        .map(|(name, ty)| {
            Field::new(
                name,
                match ty {
                    Type::Count | Type::Bases => DataType::UInt64,
                    Type::Number => DataType::Float64,
                    Type::Truth | Type::String => DataType::String,
                    _ => DataType::Json,
                },
            )
            .nullable()
        })
        .collect();
    let schema = TableSchema::new("gravlax.gq.table.v1", fields)?
        .with_semantics(TableSemantics::new(RowSemantics::Multiset))?;
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        "gravlax.gq.result.v1",
        OutputFormat::Json,
        context,
        &table.summary,
    )?;
    let available = table
        .summary
        .get("available_rows")
        .and_then(|v| v.as_u64())
        .unwrap_or(table.rows.len() as u64);
    let selection = SelectionSummary::selected(available, table.rows.len() as u64)?;
    bundle.write_table("results", &schema, Some(&selection), |rows| {
        for data in &table.rows {
            rows.write_row_with(|row| {
                for (name, _) in &table.columns {
                    match data.get(name).unwrap_or(&Value::Null) {
                        Value::Null => row.null()?,
                        Value::Count(n) => row.uint64(*n)?,
                        Value::Bases(n) => row.uint64(u64::from(*n))?,
                        Value::Number(n) => row.float64(*n)?,
                        Value::String(s) => row.string(s)?,
                        Value::Truth(t) => row.string(t.name())?,
                        v => row.json(&v.json())?,
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    })?;
    bundle.finish()?;
    Ok(())
}

fn logical_digest(plan: &model::Plan) -> Result<String> {
    fn normalize(value: &mut serde_json::Value) {
        match value {
            serde_json::Value::Object(map) => {
                map.remove("at");
                for v in map.values_mut() {
                    normalize(v);
                }
            }
            serde_json::Value::Array(values) => {
                for v in values {
                    normalize(v);
                }
            }
            _ => {}
        }
    }
    let mut value = serde_json::to_value(plan)?;
    normalize(&mut value);
    Ok(blake3::hash(&serde_json::to_vec(&value)?)
        .to_hex()
        .to_string())
}
