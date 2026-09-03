//! Biological-identifier resolution at the CLI boundary.
//!
//! Resolution is deliberately completed before an output sink is opened. A batch containing an
//! unknown or ambiguous identifier therefore cannot leave a plausible-looking partial result.

use anno::intent::{
    AnnotationIdentity, FeatureKind, IntentResolver, MatchBasis, ResolvedFeature, Strand,
};
use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use gravlax_output::{
    publish_file_no_clobber, write_table, DataType, Durability, Field, OutputFormat, Producer,
    Provenance, ResultContext, ScalarValue, TableSchema, TypedTable,
};
use serde_json::json;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const RESULT_SCHEMA: &str = "gravlax.annotation.resolve.v1";

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Format {
    Text,
    Tsv,
    Json,
}

impl From<Format> for OutputFormat {
    fn from(value: Format) -> Self {
        match value {
            Format::Text => Self::Text,
            Format::Tsv => Self::Tsv,
            Format::Json => Self::Json,
        }
    }
}

#[derive(Parser, Debug)]
pub struct Args {
    /// GTF or compiled `.aic` annotation used to resolve the identifiers.
    #[arg(value_name = "ANNOTATION_FILE")]
    pub annotation_file: PathBuf,

    /// Gene symbol or stable identifier. Prefix with `gene:`, `transcript:`, or `exon:` to
    /// constrain its kind. Every identifier must resolve unambiguously or no result is written.
    #[arg(value_name = "IDENTIFIER", required = true)]
    pub identifiers: Vec<String>,

    /// Exact reference assembly identity (for example `GRCh38.p14`).
    #[arg(long, value_name = "ASSEMBLY")]
    pub assembly: String,

    /// Exact annotation release identity (for example `GENCODE 49`).
    #[arg(long, value_name = "RELEASE")]
    pub annotation: String,

    /// Optional expected annotation content identity, as `blake3:<64 lowercase hex>`.
    #[arg(long, alias = "digest", value_name = "BLAKE3")]
    pub annotation_digest: Option<String>,

    /// Result representation. JSON is a `gravlax.result-envelope.v1` typed-table result.
    #[arg(long, value_enum, default_value_t = Format::Text)]
    pub format: Format,

    /// Write to a new file instead of standard output. Existing files are never overwritten.
    #[arg(long, short = 'o', value_name = "FILE")]
    pub output: Option<PathBuf>,
}

fn schema() -> Result<TableSchema> {
    let mut fields = vec![
        Field::new("requested", DataType::String),
        Field::new("kind", DataType::String),
        Field::new("stable_id", DataType::String),
        Field::new("display_name", DataType::String).nullable(),
        Field::new("matched_by", DataType::String),
        Field::new("gene_ids", DataType::Json),
        Field::new("transcript_ids", DataType::Json),
        Field::new("contig", DataType::String).nullable(),
        Field::new("start", DataType::UInt64).nullable(),
        Field::new("end", DataType::UInt64).nullable(),
        Field::new("strand", DataType::String).nullable(),
    ];
    fields[7].description =
        Some("Reference contig; null only when the annotation feature has no locus".into());
    fields[8].description = Some("Zero-based inclusive genomic start".into());
    fields[9].description = Some("Zero-based exclusive genomic end".into());
    Ok(TableSchema::new(RESULT_SCHEMA, fields)?)
}

fn feature_kind(value: FeatureKind) -> &'static str {
    match value {
        FeatureKind::Gene => "gene",
        FeatureKind::Transcript => "transcript",
        FeatureKind::Exon => "exon",
    }
}

fn match_basis(value: MatchBasis) -> &'static str {
    match value {
        MatchBasis::StableId => "stable_id",
        MatchBasis::StableIdWithoutVersion => "stable_id_without_version",
        MatchBasis::GeneSymbol => "gene_symbol",
    }
}

fn strand(value: Strand) -> &'static str {
    match value {
        Strand::Forward => "+",
        Strand::Reverse => "-",
    }
}

fn common_cells(feature: &ResolvedFeature) -> Vec<ScalarValue> {
    vec![
        feature.requested.clone().into(),
        feature_kind(feature.kind).into(),
        feature.stable_id.clone().into(),
        feature
            .display_name
            .clone()
            .map_or(ScalarValue::Null, ScalarValue::String),
        match_basis(feature.matched_by).into(),
        ScalarValue::Json(json!(feature.gene_ids)),
        ScalarValue::Json(json!(feature.transcript_ids)),
    ]
}

fn rows(features: &[ResolvedFeature]) -> Vec<Vec<ScalarValue>> {
    let mut rows = Vec::new();
    for feature in features {
        let common = common_cells(feature);
        if feature.loci.is_empty() {
            let mut row = common;
            row.extend([
                ScalarValue::Null,
                ScalarValue::Null,
                ScalarValue::Null,
                ScalarValue::Null,
            ]);
            rows.push(row);
            continue;
        }
        for locus in &feature.loci {
            let mut row = common.clone();
            row.extend([
                locus.contig.clone().into(),
                u64::from(locus.start).into(),
                u64::from(locus.end).into(),
                strand(locus.strand).into(),
            ]);
            rows.push(row);
        }
    }
    rows
}

fn write_new_file(
    path: &Path,
    table: &TypedTable,
    format: OutputFormat,
    context: &ResultContext,
) -> Result<()> {
    if path.exists() {
        bail!("refusing to overwrite {}", path.display());
    }
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    if !parent.is_dir() {
        bail!("output directory does not exist: {}", parent.display());
    }
    path.file_name().context("output must name a file")?;
    let outcome = publish_file_no_clobber(path, Durability::File, |writer| {
        write_table(
            writer,
            &table.schema,
            table.rows.clone(),
            format,
            context,
        )
    })
    .with_context(|| format!("installing new output {}", path.display()))?;
    for warning in outcome.warnings {
        eprintln!("warning: {warning}");
    }
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    let annotation_source = args
        .annotation_file
        .to_str()
        .context("annotation file path is not valid UTF-8 and cannot be recorded exactly in result provenance")?
        .to_owned();
    let mut identity = AnnotationIdentity::new(&args.assembly, &args.annotation)?;
    if let Some(digest) = &args.annotation_digest {
        identity = identity.with_digest(digest)?;
    }
    let resolver = IntentResolver::from_path(&args.annotation_file, identity)?;

    // Resolve the whole request before creating a destination. This is the fail-closed boundary.
    let features: Vec<ResolvedFeature> = args
        .identifiers
        .iter()
        .map(|identifier| {
            resolver
                .resolve_str(identifier)
                .with_context(|| format!("resolving identifier {identifier:?}"))
        })
        .collect::<Result<_>>()?;

    let schema = schema()?;
    let table = TypedTable::new(schema, rows(&features))?;
    let mut parameters = std::collections::BTreeMap::new();
    parameters.insert("identifiers".into(), json!(args.identifiers));
    parameters.insert("annotation_file".into(), json!(annotation_source));
    let context = ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            assembly: Some(resolver.identity().assembly.clone()),
            annotation: Some(resolver.identity().annotation.clone()),
            annotation_digest: resolver.identity().digest.clone(),
            parameters,
            ..Default::default()
        },
        warnings: Vec::new(),
    };
    let format = OutputFormat::from(args.format);
    if let Some(output) = &args.output {
        write_new_file(output, &table, format, &context)
    } else {
        let stdout = std::io::stdout();
        let mut writer = BufWriter::new(stdout.lock());
        write_table(&mut writer, &table.schema, table.rows, format, &context)?;
        writer.flush()?;
        Ok(())
    }
}
