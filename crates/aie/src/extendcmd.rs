//! `aie extend` — evidence-supported per-gene 3' extension: the Pool et al. (Nat Methods 2023)
//! reference-optimization workflow as an index query.
//!
//! Pool et al. redefined gene 3' ends by ranking genes on intergenic read mass within 10 kb of
//! the annotated end and manually reviewing coverage in IGV — under a thousand genes, one genome
//! browser session at a time. The archive automates the evidence side of that review: for each
//! gene, molecules ending downstream of the annotated 3' end are clustered into candidate
//! cleavage sites; a site extends the gene only if it is (a) reachable from the annotated end
//! without a coverage gap larger than `--evidence-gap`, (b) supported by enough molecules and
//! cells, (c) clipped before the next downstream gene, and (d) — when the genome is supplied —
//! not explainable as oligo(dT) internal priming on an A-rich tract. The output GTF feeds
//! straight back into `replay-rows`: discover → extend → replay, no realignment.

use crate::archivecmd::{decode_chunk, read_chunk_index, LazyArchive};
use crate::querycmd::{call_sites, SiteCall};
use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use gravlax_output::{
    publish_file_no_clobber, reported_output_path, DataType, Durability, Field, OutputError,
    OutputFormat, Producer, Provenance, ResultContext, RowSemantics, SelectionSummary,
    StreamingBundleWriter, TableSchema, TableSemantics,
};
use rustc_hash::FxHashMap;
use serde::Serialize;
use serde_json::json;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

const EXTEND_RESULT_SCHEMA: &str = "gravlax.extend.result.v1";
const EXTEND_ARTIFACTS_SCHEMA: &str = "gravlax.extend.artifacts.v1";
const EXTEND_GENES_SCHEMA: &str = "gravlax.extend.genes.v1";

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum UniformExtendReportFormat {
    Text,
    Tsv,
    Json,
}

impl From<UniformExtendReportFormat> for OutputFormat {
    fn from(value: UniformExtendReportFormat) -> Self {
        match value {
            UniformExtendReportFormat::Text => Self::Text,
            UniformExtendReportFormat::Tsv => Self::Tsv,
            UniformExtendReportFormat::Json => Self::Json,
        }
    }
}

#[derive(Parser)]
pub struct Args {
    /// `.aie` archive.
    pub archive: PathBuf,
    /// Annotation to extend.
    #[arg(long)]
    pub gtf: PathBuf,
    /// Extended GTF out.
    #[arg(long)]
    pub out_gtf: PathBuf,
    /// Per-gene extension report TSV (default: summary to stdout only).
    #[arg(long)]
    pub report: Option<PathBuf>,
    /// Emit a versioned uniform report while preserving --out-gtf and the legacy --report TSV.
    #[arg(long, value_enum)]
    pub report_format: Option<UniformExtendReportFormat>,
    /// Atomically publish the uniform report without replacing an existing file.
    #[arg(long, requires = "report_format")]
    pub report_output: Option<PathBuf>,
    /// Reference genome FASTA: drop candidate sites that look like internal priming. Verified
    /// against the archive's stamped genome signature.
    #[arg(long)]
    pub genome: Option<PathBuf>,
    /// Furthest a gene may extend past its annotated 3' end (Pool et al. searched 10 kb).
    #[arg(long, default_value_t = 10_000)]
    pub max_extend: u32,
    /// Largest tolerated gap in molecule-span coverage between the annotated end and the new end.
    #[arg(long, default_value_t = 2_000)]
    pub evidence_gap: u32,
    /// Minimum molecules at the accepting 3' site.
    #[arg(long, default_value_t = 5)]
    pub min_umis: usize,
    /// Minimum distinct cells at the accepting 3' site.
    #[arg(long, default_value_t = 3)]
    pub min_cells: usize,
    #[arg(long, default_value_t = 24)]
    pub site_gap: u32,
    /// Extensions shorter than this are not worth an annotation edit.
    #[arg(long, default_value_t = 50)]
    pub min_extension: u32,
    /// Clip at the next gene on either strand instead of the next same-strand gene (3' scRNA-seq
    /// reads are stranded, so same-strand clipping is the default; this is the conservative mode).
    #[arg(long)]
    pub clip_any_strand: bool,
}

struct UniformPaths {
    archive: String,
    gtf: String,
    out_gtf: String,
    legacy_report: Option<String>,
    genome: Option<String>,
}

#[derive(Serialize)]
struct ExtendSummary {
    coordinates: &'static str,
    genes_in_annotation: u64,
    transcripts_in_annotation: u64,
    genes_extended: u64,
    total_extension_bp: u64,
    qualifying_sites_for_extended_genes: u64,
    internal_priming_sites_dropped_for_extended_genes: u64,
    extensions_clipped_by_neighbor: u64,
    gtf_lines_written: u64,
    genome_filter_enabled: bool,
}

struct ExtendArtifactReport<'a> {
    gtf_path: &'a str,
    gtf_bytes: u64,
    gtf_lines: u64,
    legacy_report_path: Option<&'a str>,
    legacy_report_bytes: Option<u64>,
    legacy_report_rows: u64,
}

fn existing_parent_key(path: &Path, label: &str) -> Result<PathBuf> {
    let file_name = path
        .file_name()
        .with_context(|| format!("{label} path must name a file: {}", path.display()))?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let parent = std::fs::canonicalize(parent)
        .with_context(|| format!("checking parent directory for {label} {}", path.display()))?;
    if !parent.is_dir() {
        bail!("{label} parent is not a directory: {}", parent.display());
    }
    Ok(parent.join(file_name))
}

/// Validate every destination before opening the archive. The existence check is an early UX
/// guard; `publish_file_no_clobber` remains the race-safe commit point for the uniform report.
fn preflight_destinations(args: &Args) -> Result<Option<UniformPaths>> {
    let out_gtf_key = existing_parent_key(&args.out_gtf, "--out-gtf")?;
    let legacy_report_key = args.report.as_deref()
        .map(|path| existing_parent_key(path, "--report"))
        .transpose()?;
    if legacy_report_key.as_ref().is_some_and(|key| key == &out_gtf_key) {
        bail!("--out-gtf and legacy --report must name different files");
    }
    if args.report_format.is_none() {
        return Ok(None);
    }

    let report_output_key = args.report_output.as_deref()
        .map(|path| existing_parent_key(path, "--report-output"))
        .transpose()?;
    if let (Some(path), Some(key)) = (args.report_output.as_deref(), report_output_key.as_ref()) {
        match std::fs::symlink_metadata(path) {
            Ok(_) => bail!("refusing to replace existing report output {}", path.display()),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(error).with_context(|| {
                format!("checking uniform report destination {}", path.display())
            }),
        }
        if key == &out_gtf_key {
            bail!("--report-output and --out-gtf must name different files");
        }
        if legacy_report_key.as_ref().is_some_and(|legacy| legacy == key) {
            bail!("--report-output and legacy --report must name different files");
        }
    }

    let utf8 = |path: &Path, label: &str| -> Result<String> {
        path.to_str().map(str::to_owned)
            .with_context(|| format!("uniform report requires a UTF-8 {label} path"))
    };
    Ok(Some(UniformPaths {
        archive: utf8(&args.archive, "archive")?,
        gtf: utf8(&args.gtf, "GTF input")?,
        out_gtf: reported_output_path(&args.out_gtf)?,
        legacy_report: args.report.as_deref()
            .map(reported_output_path)
            .transpose()?,
        genome: args.genome.as_deref()
            .map(|path| utf8(path, "genome"))
            .transpose()?,
    }))
}

struct Gene {
    id: String,
    name: String,
    chrom: String,
    rev: bool,
    start0: u32,
    end0: u32, // exclusive
    // corridor downstream of the 3' end, in genome coordinates
    cor_lo: u32,
    cor_hi: u32,
    clip_reason: &'static str,
    pts: Vec<(u32, bool, u32, u32)>,   // (tp, rev, cell, class) with tp in corridor
    cov: Vec<(u32, u32)>,              // strand-matched molecule spans clipped to the corridor
}

struct Ext {
    gene: usize,
    new_end0: u32, // exclusive end on +, inclusive start on −
    site_umis: usize,
    site_cells: usize,
    n_sites: usize,
    ip_dropped: usize,
}

fn attr(attrs: &str, key: &str) -> Option<String> {
    attrs.split(';').find_map(|kv| {
        let kv = kv.trim();
        kv.strip_prefix(key)
            .map(|v| v.trim().trim_matches('"').to_string())
            .filter(|v| !v.is_empty())
    })
}

fn artifact_schema() -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(
        EXTEND_ARTIFACTS_SCHEMA,
        vec![
            Field::new("artifact_kind", DataType::String),
            Field::new("path", DataType::String),
            Field::new("bytes", DataType::UInt64),
            Field::new("records", DataType::UInt64),
            Field::new("record_unit", DataType::String),
        ],
    )?
    .with_semantics(
        TableSemantics::new(RowSemantics::Set).with_key(["artifact_kind", "path"]),
    )
}

fn extension_schema() -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(
        EXTEND_GENES_SCHEMA,
        vec![
            Field::new("gene_index", DataType::UInt64),
            Field::new("gene_id", DataType::String),
            Field::new("gene_name", DataType::String),
            Field::new("chrom", DataType::String),
            Field::new("strand", DataType::String),
            Field::new("old_end", DataType::UInt64),
            Field::new("new_end", DataType::UInt64),
            Field::new("extension_bp", DataType::UInt64),
            Field::new("site_umis", DataType::UInt64),
            Field::new("site_cells", DataType::UInt64),
            Field::new("qualifying_sites", DataType::UInt64),
            Field::new("internal_priming_sites_dropped", DataType::UInt64),
            Field::new("corridor_clip_reason", DataType::String),
        ],
    )?
    // The existing result vector follows annotation traversal, but that incidental physical order
    // is not part of the contract. `gene_index` is the stable per-input key; no sort is added.
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["gene_index"]))
}

fn uniform_context(
    args: &Args,
    paths: &UniformPaths,
    archive_version: u32,
    archive_root: Option<String>,
    annotation_digest: String,
    genome_identity: Option<&str>,
) -> ResultContext {
    let mut parameters = BTreeMap::new();
    parameters.insert("archive_path".into(), json!(paths.archive));
    parameters.insert("archive_version".into(), json!(archive_version));
    parameters.insert("out_gtf".into(), json!(paths.out_gtf));
    parameters.insert("legacy_report".into(), json!(paths.legacy_report));
    parameters.insert("genome".into(), json!(paths.genome));
    parameters.insert("genome_identity".into(), json!(genome_identity));
    parameters.insert("max_extend".into(), json!(args.max_extend));
    parameters.insert("evidence_gap".into(), json!(args.evidence_gap));
    parameters.insert("min_umis".into(), json!(args.min_umis));
    parameters.insert("min_cells".into(), json!(args.min_cells));
    parameters.insert("site_gap".into(), json!(args.site_gap));
    parameters.insert("min_extension".into(), json!(args.min_extension));
    parameters.insert("clip_any_strand".into(), json!(args.clip_any_strand));
    let mut warnings = Vec::new();
    let archives = match archive_root {
        Some(root) => vec![format!("aie-directory-root-v2:{root}")],
        None => {
            warnings.push(
                "legacy v1 archive has no rooted content commitment; its path locator is not a portable content identity"
                    .into(),
            );
            vec![format!("aie-v{archive_version}-unrooted:{}", paths.archive)]
        }
    };
    ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            archives,
            annotation: Some(paths.gtf.clone()),
            annotation_digest: Some(annotation_digest),
            parameters,
            ..Default::default()
        },
        warnings,
    }
}

fn write_uniform_report<W: Write>(
    writer: W,
    format: OutputFormat,
    context: &ResultContext,
    summary: &ExtendSummary,
    artifacts: &ExtendArtifactReport<'_>,
    exts: &[Ext],
    genes: &[Gene],
) -> std::result::Result<W, OutputError> {
    let artifacts_schema = artifact_schema()?;
    let extensions_schema = extension_schema()?;
    let artifact_selection = SelectionSummary::complete(
        1 + u64::from(artifacts.legacy_report_path.is_some()),
    );
    let extension_selection = SelectionSummary::complete(exts.len() as u64);
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer, EXTEND_RESULT_SCHEMA, format, context, summary,
    )?;
    bundle.write_table("artifacts", &artifacts_schema, Some(&artifact_selection), |rows| {
        rows.write_row_with(|row| {
            row.string("extended_gtf")?;
            row.string(artifacts.gtf_path)?;
            row.uint64(artifacts.gtf_bytes)?;
            row.uint64(artifacts.gtf_lines)?;
            row.string("lines")?;
            Ok(())
        })?;
        if let Some(path) = artifacts.legacy_report_path {
            rows.write_row_with(|row| {
                row.string("legacy_per_gene_tsv")?;
                row.string(path)?;
                row.uint64(artifacts.legacy_report_bytes.unwrap_or(0))?;
                row.uint64(artifacts.legacy_report_rows)?;
                row.string("data_rows")?;
                Ok(())
            })?;
        }
        Ok(())
    })?;
    bundle.write_table("extensions", &extensions_schema, Some(&extension_selection), |rows| {
        for extension in exts {
            let gene = &genes[extension.gene];
            let (old_end, extension_bp) = if gene.rev {
                (gene.start0, gene.start0 - extension.new_end0)
            } else {
                (gene.end0, extension.new_end0 - gene.end0)
            };
            rows.write_row_with(|row| {
                row.uint64(extension.gene as u64)?;
                row.string(&gene.id)?;
                row.string(&gene.name)?;
                row.string(&gene.chrom)?;
                row.string(if gene.rev { "-" } else { "+" })?;
                row.uint64(old_end as u64)?;
                row.uint64(extension.new_end0 as u64)?;
                row.uint64(extension_bp as u64)?;
                row.uint64(extension.site_umis as u64)?;
                row.uint64(extension.site_cells as u64)?;
                row.uint64(extension.n_sites as u64)?;
                row.uint64(extension.ip_dropped as u64)?;
                row.string(gene.clip_reason)?;
                Ok(())
            })?;
        }
        Ok(())
    })?;
    bundle.finish()
}

fn write_legacy_report<W: Write>(mut writer: W, exts: &[Ext], genes: &[Gene]) -> std::io::Result<W> {
    writeln!(writer, "#gene_id\tgene_name\tchrom\tstrand\told_end\tnew_end\text_bp\tsite_umis\tsite_cells\tn_sites\tip_dropped\tclip")?;
    let mut sorted: Vec<&Ext> = exts.iter().collect();
    sorted.sort_by_key(|e| {
        let g = &genes[e.gene];
        std::cmp::Reverse(if !g.rev { e.new_end0 - g.end0 } else { g.start0 - e.new_end0 })
    });
    for e in sorted {
        let g = &genes[e.gene];
        let (old, ext) = if !g.rev { (g.end0, e.new_end0 - g.end0) } else { (g.start0, g.start0 - e.new_end0) };
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.id, g.name, g.chrom, if g.rev { '-' } else { '+' }, old, e.new_end0, ext,
            e.site_umis, e.site_cells, e.n_sites, e.ip_dropped, g.clip_reason
        )?;
    }
    writer.flush()?;
    Ok(writer)
}

pub fn run(args: Args) -> Result<()> {
    let t0 = std::time::Instant::now();
    let uniform_paths = preflight_destinations(&args)?;
    let gtf_text = if args.gtf.extension().is_some_and(|e| e == "gz") {
        bail!("--gtf must be uncompressed (the rewrite preserves the file line-for-line)");
    } else {
        std::fs::read_to_string(&args.gtf).with_context(|| format!("reading {}", args.gtf.display()))?
    };
    let (annotation_digest, lines): (Option<String>, Vec<&str>) = if uniform_paths.is_some() {
        // Exact identity is uniform-only. Hash in parallel with the already-required line-index
        // pass so a multi-gigabyte GTF does not incur another serial, result-wide traversal.
        let (digest, lines) = rayon::join(
            || blake3::hash(gtf_text.as_bytes()),
            || gtf_text.lines().collect(),
        );
        (Some(format!("blake3:{}", digest.to_hex())), lines)
    } else {
        (None, gtf_text.lines().collect())
    };

    // Pass 1 over the GTF: genes (id, span, strand) and, per transcript, the line indices of its
    // transcript record and terminal exon — everything the rewrite needs.
    let mut genes: Vec<Gene> = Vec::new();
    let mut gene_idx: FxHashMap<String, usize> = FxHashMap::default();
    // transcript_id -> (gene_id, tx line idx, tx end0, tx start0, terminal exon line idx + coord)
    struct Tx {
        gene: String,
        line: usize,
        start0: u32,
        end0: u32,
        term_exon_line: usize,
        term_exon_coord: u32, // end0 on +, start0 on −
        rev: bool,
    }
    let mut txs: FxHashMap<String, Tx> = FxHashMap::default();
    for (li, line) in lines.iter().enumerate() {
        if line.starts_with('#') {
            continue;
        }
        let mut f = line.split('\t');
        let (Some(chrom), Some(_src), Some(feat), Some(s1), Some(e1), Some(_score), Some(strand), Some(_frame), Some(attrs)) =
            (f.next(), f.next(), f.next(), f.next(), f.next(), f.next(), f.next(), f.next(), f.next())
        else {
            continue;
        };
        let (Ok(s1), Ok(e1)) = (s1.parse::<u32>(), e1.parse::<u32>()) else { continue };
        let (start0, end0) = (s1 - 1, e1); // GTF is 1-based inclusive
        let rev = strand == "-";
        match feat {
            "gene" => {
                let Some(id) = attr(attrs, "gene_id") else { continue };
                gene_idx.insert(id.clone(), genes.len());
                genes.push(Gene {
                    id,
                    name: attr(attrs, "gene_name").unwrap_or_default(),
                    chrom: chrom.to_string(),
                    rev,
                    start0,
                    end0,
                    cor_lo: 0,
                    cor_hi: 0,
                    clip_reason: "max",
                    pts: Vec::new(),
                    cov: Vec::new(),
                });
            }
            "transcript" => {
                let (Some(gid), Some(tid)) = (attr(attrs, "gene_id"), attr(attrs, "transcript_id")) else {
                    continue;
                };
                txs.insert(tid, Tx {
                    gene: gid,
                    line: li,
                    start0,
                    end0,
                    term_exon_line: usize::MAX,
                    term_exon_coord: if rev { u32::MAX } else { 0 },
                    rev,
                });
            }
            "exon" => {
                let Some(tid) = attr(attrs, "transcript_id") else { continue };
                if let Some(tx) = txs.get_mut(&tid) {
                    let better = if tx.rev { start0 < tx.term_exon_coord } else { end0 > tx.term_exon_coord };
                    if better || tx.term_exon_line == usize::MAX {
                        tx.term_exon_line = li;
                        tx.term_exon_coord = if tx.rev { start0 } else { end0 };
                    }
                }
            }
            _ => {}
        }
    }
    eprintln!("gtf: {} genes, {} transcripts ({:.1}s)", genes.len(), txs.len(), t0.elapsed().as_secs_f32());

    // Corridors, clipped at the next downstream gene (same strand unless --clip-any-strand).
    let mut by_chrom: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    for (i, g) in genes.iter().enumerate() {
        by_chrom.entry(g.chrom.clone()).or_default().push(i);
    }
    for idxs in by_chrom.values() {
        for &i in idxs {
            let (rev, e0, s0, strand_rev) = {
                let g = &genes[i];
                (g.rev, g.end0, g.start0, g.rev)
            };
            // Clip at any other gene occupying the corridor: a gene starting inside it, or one
            // already straddling our 3' end (which collapses the corridor entirely). Extending
            // into an occupied region converts uniquely-assigned molecules into discarded
            // multi-gene ones — the extension must stop where another gene's claim begins.
            let mut clip: Option<u32> = None;
            for &j in idxs {
                if i == j {
                    continue;
                }
                let o = &genes[j];
                if !args.clip_any_strand && o.rev != strand_rev {
                    continue;
                }
                if !rev {
                    if o.end0 > e0 {
                        let b = o.start0.max(e0);
                        clip = Some(clip.map_or(b, |c: u32| c.min(b)));
                    }
                } else if o.start0 < s0 {
                    let b = o.end0.min(s0);
                    clip = Some(clip.map_or(b, |c: u32| c.max(b)));
                }
            }
            let g = &mut genes[i];
            if !rev {
                g.cor_lo = e0;
                g.cor_hi = e0 + args.max_extend;
                if let Some(c) = clip {
                    if c < g.cor_hi {
                        g.cor_hi = c;
                        g.clip_reason = "neighbor";
                    }
                }
            } else {
                g.cor_hi = s0;
                g.cor_lo = s0.saturating_sub(args.max_extend);
                if let Some(c) = clip {
                    if c > g.cor_lo {
                        g.cor_lo = c;
                        g.clip_reason = "neighbor";
                    }
                }
            }
        }
    }

    // One pass over the chunks: molecules' 3' points into corridors; spans into coverage.
    let mut la = LazyArchive::open(&args.archive)?;
    if uniform_paths.is_some() && args.genome.is_some() && la.genome_sig.is_none() {
        bail!(
            "uniform extension with --genome requires a stamped archive; run `aie stamp-genome` \
             first so the scientific FASTA input can be verified"
        );
    }
    let chunks = read_chunk_index(la.reader())?;
    let chrom_names = la.chrom_names.clone();
    let shapes = la.shapes()?.to_vec();
    // Per archive chrom: corridor list sorted by cor_lo with running max cor_hi.
    let mut cors: Vec<Vec<usize>> = vec![Vec::new(); chrom_names.len()];
    for (i, g) in genes.iter().enumerate() {
        if g.cor_hi > g.cor_lo {
            if let Some(cid) = chrom_names.iter().position(|n| *n == g.chrom) {
                cors[cid].push(i);
            }
        }
    }
    let mut maxhi: Vec<Vec<u32>> = Vec::with_capacity(cors.len());
    for list in &mut cors {
        list.sort_unstable_by_key(|&i| genes[i].cor_lo);
        let mut run = 0u32;
        maxhi.push(list.iter().map(|&i| { run = run.max(genes[i].cor_hi); run }).collect());
    }
    for (ci, c) in chunks.iter().enumerate() {
        let list = &cors[c.chrom as usize];
        if list.is_empty() {
            continue;
        }
        let mh = &maxhi[c.chrom as usize];
        let raw = la.reader().read(&format!("c{ci}"))?;
        for m in decode_chunk(&raw, c, None, &la.rans_tables)? {
            // Evidence span of the molecule (reps and mm anchors), as discover computes it.
            let (mut lo, mut hi) = (u32::MAX, 0u32);
            for (pos, sh) in m.chains.iter().flat_map(|ch| ch.reps.iter()) {
                lo = lo.min(*pos);
                hi = hi.max(*pos + shapes[*sh as usize].blocks.last().map(|b| b.0 + b.1).unwrap_or(0));
            }
            for (pos, sh, _, _) in &m.mms {
                lo = lo.min(*pos);
                hi = hi.max(*pos + shapes[*sh as usize].blocks.last().map(|b| b.0 + b.1).unwrap_or(0));
            }
            if lo == u32::MAX {
                continue;
            }
            let tp = if m.strand_rev { lo } else { hi };
            // Corridors overlapping the span (coverage) and containing tp (site points).
            let mut k = list.partition_point(|&gi| genes[gi].cor_lo <= hi.min(u32::MAX - 1));
            let mut cell: Option<u32> = None;
            let mut tp_best: Option<usize> = None;
            while k > 0 {
                k -= 1;
                if mh[k] <= lo {
                    break;
                }
                let gi = list[k];
                let g = &genes[gi];
                if g.rev != m.strand_rev || g.cor_lo >= hi || g.cor_hi <= lo {
                    continue;
                }
                let (clo, chi) = (genes[gi].cor_lo, genes[gi].cor_hi);
                genes[gi].cov.push((lo.max(clo), hi.min(chi)));
                let g = &genes[gi];
                if tp >= g.cor_lo && tp < g.cor_hi {
                    // Nearest annotated 3' end wins when corridors overlap.
                    let better = match tp_best {
                        None => true,
                        Some(pb) => {
                            let (pg, ng) = (&genes[pb], g);
                            if ng.rev { ng.start0 < pg.start0 } else { ng.end0 > pg.end0 }
                        }
                    };
                    if better {
                        tp_best = Some(gi);
                    }
                }
            }
            if let Some(gi) = tp_best {
                let cell = *cell.get_or_insert(la.cell_of(m.umi_class)?);
                genes[gi].pts.push((tp, m.strand_rev, cell, m.umi_class));
            }
        }
    }

    // Per gene: sites, internal-priming filter, coverage walk, verdict.
    let mut exts: Vec<Ext> = Vec::new();
    let no_groups: FxHashMap<u32, u32> = FxHashMap::default();
    let decide = |gi: usize, g: &mut Gene, seq: Option<&[u8]>, exts: &mut Vec<Ext>| {
        if g.pts.is_empty() {
            return;
        }
        let mut pts = std::mem::take(&mut g.pts);
        let sites = call_sites(&mut pts, args.site_gap, &no_groups, 0, seq);
        let ip_dropped = sites.iter().filter(|s| s.ip).count();
        let qual: Vec<&SiteCall> = sites
            .iter()
            .filter(|s| !s.ip && s.umis >= args.min_umis && s.cells >= args.min_cells)
            .collect();
        if qual.is_empty() {
            return;
        }
        // Coverage walk from the annotated end, tolerating gaps up to evidence_gap.
        let mut cov = std::mem::take(&mut g.cov);
        let frontier = if !g.rev {
            cov.sort_unstable();
            let mut fr = g.end0;
            for (s, e) in cov {
                if s > fr.saturating_add(args.evidence_gap) {
                    break;
                }
                fr = fr.max(e);
            }
            fr
        } else {
            cov.sort_unstable_by_key(|&(_, e)| std::cmp::Reverse(e));
            let mut fr = g.start0;
            for (s, e) in cov {
                if e.saturating_add(args.evidence_gap) < fr {
                    break;
                }
                fr = fr.min(s);
            }
            fr
        };
        let pick = if !g.rev {
            qual.iter().filter(|s| s.cp <= frontier).max_by_key(|s| s.cp)
        } else {
            qual.iter().filter(|s| s.cp >= frontier).min_by_key(|s| s.cp)
        };
        let Some(site) = pick else { return };
        let ext_len = if !g.rev { site.cp.saturating_sub(g.end0) } else { g.start0.saturating_sub(site.cp) };
        if ext_len < args.min_extension {
            return;
        }
        exts.push(Ext {
            gene: gi,
            new_end0: site.cp,
            site_umis: site.umis,
            site_cells: site.cells,
            n_sites: qual.len(),
            ip_dropped,
        });
    };
    let mut genome_identity = None;
    if let Some(fasta) = &args.genome {
        if la.genome_sig.is_none() {
            eprintln!(
                "warning: archive carries no genome signature; --genome cannot be verified \
                 (stamp it with `aie stamp-genome`)"
            );
        }
        let sig = la.genome_sig.clone();
        let mut genes_by_chrom: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        for (i, g) in genes.iter().enumerate() {
            if !g.pts.is_empty() {
                genes_by_chrom.entry(g.chrom.clone()).or_default().push(i);
            }
        }
        let mut verr: Option<anyhow::Error> = None;
        let mut observed_contigs = Vec::new();
        let capture_genome_identity = uniform_paths.is_some();
        evidence_io::genome::for_each_contig(fasta, |name, seq| {
            if capture_genome_identity {
                observed_contigs.push(evidence_io::genome::ContigSig {
                    name: name.to_owned(),
                    len: seq.len() as u64,
                    blake3: blake3::hash(seq).to_hex().to_string(),
                });
            }
            let Some(idxs) = genes_by_chrom.remove(name) else { return true };
            if let Some(sig) = &sig {
                if let Err(e) = evidence_io::genome::verify_contig(sig, name, seq) {
                    verr = Some(e);
                    return false;
                }
            }
            for gi in idxs {
                let mut g = std::mem::replace(&mut genes[gi], Gene {
                    id: String::new(), name: String::new(), chrom: String::new(), rev: false,
                    start0: 0, end0: 0, cor_lo: 0, cor_hi: 0, clip_reason: "max",
                    pts: Vec::new(), cov: Vec::new(),
                });
                decide(gi, &mut g, Some(seq), &mut exts);
                genes[gi] = g;
            }
            true
        })?;
        if let Some(e) = verr {
            return Err(e);
        }
        if capture_genome_identity {
            genome_identity = Some(format!(
                "{}:{}",
                evidence_io::genome::GENOME_SIG_ALGO,
                evidence_io::genome::GenomeSig::combined_digest(&observed_contigs)
            ));
        }
        if let Some(missed) = genes_by_chrom.keys().next() {
            bail!("contig {missed} with candidate genes not found in {}", fasta.display());
        }
    } else {
        for gi in 0..genes.len() {
            if genes[gi].pts.is_empty() {
                continue;
            }
            let mut g = std::mem::replace(&mut genes[gi], Gene {
                id: String::new(), name: String::new(), chrom: String::new(), rev: false,
                start0: 0, end0: 0, cor_lo: 0, cor_hi: 0, clip_reason: "max",
                pts: Vec::new(), cov: Vec::new(),
            });
            decide(gi, &mut g, None, &mut exts);
            genes[gi] = g;
        }
    }

    // Rewrite the GTF: gene bounds, plus every transcript (and its terminal exon) that ends at
    // the gene's annotated 3' end.
    let mut edits: FxHashMap<usize, (Option<u32>, Option<u32>)> = FxHashMap::default(); // line -> (new start1, new end1)
    let mut gene_line_of: FxHashMap<String, usize> = FxHashMap::default();
    for (li, line) in lines.iter().enumerate() {
        if line.starts_with('#') {
            continue;
        }
        let mut f = line.split('\t');
        let (Some(_), Some(_), Some(feat)) = (f.next(), f.next(), f.next()) else { continue };
        if feat == "gene" {
            if let Some(id) = f.nth(5).and_then(|a| attr(a, "gene_id")) {
                if gene_idx.contains_key(id.as_str()) {
                    gene_line_of.insert(id, li);
                }
            }
        }
    }
    for e in &exts {
        let g = &genes[e.gene];
        if let Some(&li) = gene_line_of.get(g.id.as_str()) {
            if !g.rev {
                edits.insert(li, (None, Some(e.new_end0)));
            } else {
                edits.insert(li, (Some(e.new_end0 + 1), None));
            }
        }
    }
    for tx in txs.values() {
        let Some(&gi) = gene_idx.get(&tx.gene) else { continue };
        let Some(e) = exts.iter().find(|e| e.gene == gi) else { continue };
        let g = &genes[gi];
        let at_end = if !g.rev { tx.end0 == g.end0 } else { tx.start0 == g.start0 };
        if !at_end {
            continue;
        }
        if !g.rev {
            edits.insert(tx.line, (None, Some(e.new_end0)));
            if tx.term_exon_line != usize::MAX {
                edits.insert(tx.term_exon_line, (None, Some(e.new_end0)));
            }
        } else {
            edits.insert(tx.line, (Some(e.new_end0 + 1), None));
            if tx.term_exon_line != usize::MAX {
                edits.insert(tx.term_exon_line, (Some(e.new_end0 + 1), None));
            }
        }
    }
    let mut w = std::io::BufWriter::new(std::fs::File::create(&args.out_gtf)?);
    for (li, line) in lines.iter().enumerate() {
        match edits.get(&li) {
            None => writeln!(w, "{line}")?,
            Some((ns, ne)) => {
                let mut fields: Vec<&str> = line.split('\t').collect();
                let s_new = ns.map(|v| v.to_string());
                let e_new = ne.map(|v| v.to_string());
                if let Some(sv) = &s_new {
                    fields[3] = sv;
                }
                if let Some(ev) = &e_new {
                    fields[4] = ev;
                }
                writeln!(w, "{}", fields.join("\t"))?;
            }
        }
    }
    w.flush()?;

    if let Some(rp) = &args.report {
        let writer = std::io::BufWriter::new(std::fs::File::create(rp)?);
        let _writer = write_legacy_report(writer, &exts, &genes)?;
    }
    let total_ext: u64 = exts.iter().map(|e| {
        let g = &genes[e.gene];
        (if !g.rev { e.new_end0 - g.end0 } else { g.start0 - e.new_end0 }) as u64
    }).sum();
    if let Some(paths) = uniform_paths.as_ref() {
        let archive_version = la.reader().archive_version();
        let archive_root = la.reader().content_commitment().map(|commitment| commitment.to_hex());
        let mut context = uniform_context(
            &args,
            paths,
            archive_version,
            archive_root,
            annotation_digest.expect("uniform annotation identity was captured"),
            genome_identity.as_deref(),
        );
        if args.genome.is_some() && la.genome_sig.is_none() {
            context.warnings.push(
                "archive carries no genome signature; the supplied genome could not be checked against the archive"
                    .into(),
            );
        }
        let summary = ExtendSummary {
            coordinates: "strand-directed 3-prime boundary: plus is the 0-based exclusive end (equal to the GTF end); minus is the 0-based start",
            genes_in_annotation: genes.len() as u64,
            transcripts_in_annotation: txs.len() as u64,
            genes_extended: exts.len() as u64,
            total_extension_bp: total_ext,
            qualifying_sites_for_extended_genes: exts.iter().map(|ext| ext.n_sites as u64).sum(),
            internal_priming_sites_dropped_for_extended_genes: exts.iter()
                .map(|ext| ext.ip_dropped as u64).sum(),
            extensions_clipped_by_neighbor: exts.iter()
                .filter(|ext| genes[ext.gene].clip_reason == "neighbor").count() as u64,
            gtf_lines_written: lines.len() as u64,
            genome_filter_enabled: args.genome.is_some(),
        };
        let artifacts = ExtendArtifactReport {
            gtf_path: &paths.out_gtf,
            gtf_bytes: std::fs::metadata(&args.out_gtf)?.len(),
            gtf_lines: lines.len() as u64,
            legacy_report_path: paths.legacy_report.as_deref(),
            legacy_report_bytes: args.report.as_deref()
                .map(std::fs::metadata).transpose()?.map(|metadata| metadata.len()),
            legacy_report_rows: exts.len() as u64,
        };
        let format = args.report_format
            .expect("uniform paths require an explicit report format");
        if let Some(path) = args.report_output.as_deref() {
            let outcome = publish_file_no_clobber(path, Durability::Flush, |writer| {
                write_uniform_report(
                    writer, format.into(), &context, &summary, &artifacts, &exts, &genes,
                )?;
                Ok(())
            })?;
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
        } else {
            let stdout = std::io::stdout();
            let writer = std::io::BufWriter::new(stdout.lock());
            let _writer = write_uniform_report(
                writer, format.into(), &context, &summary, &artifacts, &exts, &genes,
            )?;
        }
        eprintln!(
            "extend: {} of {} genes extended, {} bp total, wrote {} ({:.1}s)",
            exts.len(), genes.len(), total_ext, args.out_gtf.display(), t0.elapsed().as_secs_f32()
        );
    } else {
        println!(
            "extend: {} of {} genes extended, {} bp total, wrote {} ({:.1}s)",
            exts.len(), genes.len(), total_ext, args.out_gtf.display(), t0.elapsed().as_secs_f32()
        );
    }
    Ok(())
}
