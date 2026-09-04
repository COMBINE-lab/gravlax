//! `aie query` — indexed access to the archive without decoding it whole.
//!
//! The format's dictionaries also serve as indexes for these queries:
//!
//! * `region`  — molecules whose archive anchor lies in a genomic window, per-cell UMI counts;
//!   reads only the chunks the window touches (the range index). Optional coverage/junction
//!   tracks reconstruct overlapping geometry.
//! * `junction` — per-cell molecule counts supporting one splice junction, via the junction
//!   catalogue + postings: Snaptron-at-molecule-resolution, the head-to-head with Malva
//!   (which answers "which cells contain the sequence", not "how many molecules per cell").
//! * `junctions` — enumerate every junction in a 0-based interval from the catalogue alone;
//!   optional postings decode adds exact per-cell molecule counts and annotation endpoint flags.
//!
//! Counting semantics: molecules are deduplicated per (cell, class) so the two span extremes of a
//! chain never count twice; no gene model is consulted — these are annotation-free queries.

use crate::apastats;
use crate::archivecmd::transcriptec::{
    derive_transcript_equivalence_classes_with_annotation, TranscriptEquivalenceOptions,
    TranscriptEquivalenceReport,
};
use crate::archivecmd::{decode_chunk, read_chunk_index, LazyArchive};
use crate::rows::{placement_from_parts_into, MolRec, PatAlt, SAME_SHAPE};
use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use evidence_io::format::Cursor;
use evidence_io::umi;
use gravlax_output::{
    publish_file_no_clobber, write_table, DataType, Durability, Field, OrderKey, OutputError,
    OutputFormat, Producer, Provenance, ResultContext, ResultEnvelope, RowSemantics, ScalarValue,
    SelectionSummary, StreamingBundleWriter, TableSchema, TableSemantics, TypedTable,
};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
pub struct Args {
    pub archive: PathBuf,
    #[command(subcommand)]
    pub what: What,
}

/// What makes archived evidence already explained by the supplied annotation during discovery.
/// `Span` is the historical conservative behavior and remains the default for compatibility.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum DiscoveryClaimMode {
    Span,
    StrandSpan,
    Compatible,
    ResidualSites,
}

impl DiscoveryClaimMode {
    fn name(self) -> &'static str {
        match self {
            Self::Span => "span",
            Self::StrandSpan => "strand_span",
            Self::Compatible => "compatible",
            Self::ResidualSites => "residual_sites",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum QueryAggregationArg {
    /// Per-cell without --groups; per-group with --groups.
    Auto,
    Cell,
    Group,
    Bulk,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, ValueEnum)]
pub enum EventTypeArg {
    AltAcceptor,
    AltDonor,
    Cassette,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum TranscriptEcFormat {
    Text,
    Json,
    Tsv,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum TranscriptEcTable {
    Catalog,
    Counts,
    Membership,
}

/// Opt-in uniform output for the first legacy-query migration.  This deliberately has no
/// default: omitting `--format` keeps the historical stdout byte contract intact.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum UniformQueryFormat {
    Text,
    Tsv,
    Json,
}

impl From<UniformQueryFormat> for OutputFormat {
    fn from(value: UniformQueryFormat) -> Self {
        match value {
            UniformQueryFormat::Text => Self::Text,
            UniformQueryFormat::Tsv => Self::Tsv,
            UniformQueryFormat::Json => Self::Json,
        }
    }
}

#[derive(clap::Args, Clone, Debug)]
pub struct UniformQueryOutputArgs {
    /// Use the versioned uniform result contract instead of the legacy presentation.
    #[arg(long, value_enum)]
    format: Option<UniformQueryFormat>,
    /// Publish atomically without replacing an existing file; requires --format.
    #[arg(short = 'o', long, requires = "format")]
    output: Option<PathBuf>,
}

#[derive(clap::Args, Clone, Debug)]
pub struct QueryScopeArgs {
    /// Headerless file with one archive barcode per line; only these cells contribute.
    #[arg(long, conflicts_with = "groups")]
    cells: Option<PathBuf>,
    /// Headerless barcode<TAB>group file; listed archive cells define the scope.
    #[arg(long, conflicts_with = "cells")]
    groups: Option<PathBuf>,
    /// Output aggregation. `auto` selects groups when --groups is present, cells otherwise.
    #[arg(long, value_enum, default_value_t = QueryAggregationArg::Auto)]
    agg: QueryAggregationArg,
}

#[derive(clap::Subcommand)]
pub enum What {
    /// Annotation-conditional compatible-transcript sets for archived UMI classes.
    TranscriptEcs {
        /// GTF or compiled AIC whose transcripts define compatibility.
        #[arg(long)]
        annotation_file: PathBuf,
        /// Exact reference assembly label attached to annotation resolution.
        #[arg(long)]
        assembly: String,
        /// Exact annotation release or immutable label.
        #[arg(long)]
        annotation_label: String,
        /// Optional expected annotation identity as blake3:<64 lowercase hex characters>.
        #[arg(long)]
        annotation_digest: Option<String>,
        /// Gene identifier or symbol resolved against the bound annotation.
        #[arg(long, required_unless_present = "locus", conflicts_with = "locus")]
        feature: Option<String>,
        /// 0-based half-open window; selects transcripts whose spans overlap it.
        #[arg(long, required_unless_present = "feature", conflicts_with = "feature")]
        locus: Option<String>,
        /// STARsolo alignment/transcript strand relationship.
        #[arg(long, value_enum, default_value_t = crate::archivecmd::SoloStrandArg::Forward)]
        solo_strand: crate::archivecmd::SoloStrandArg,
        #[command(flatten)]
        scope: QueryScopeArgs,
        /// Include one row per scoped archived UMI class.
        #[arg(long)]
        emit_membership: bool,
        /// Hard output limit for distinct compatible transcript sets; never truncates.
        #[arg(long, default_value_t = 100_000)]
        max_ecs: usize,
        /// Hard output limit for membership rows; never truncates.
        #[arg(long, default_value_t = 1_000_000)]
        max_memberships: usize,
        #[arg(long, value_enum, default_value_t = TranscriptEcFormat::Text)]
        format: TranscriptEcFormat,
        /// Required for TSV; unavailable for JSON and text.
        #[arg(long, value_enum)]
        table: Option<TranscriptEcTable>,
        /// Write atomically without replacing an existing path.
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },
    /// Execute many anchor-region and exact-junction predicates with shared archive/chunk work.
    Batch {
        /// Strict TSV with header: id, kind, locus (kind is region or junction).
        #[arg(long)]
        plan: PathBuf,
        /// Number of per-cell rows per result (0 = all).
        #[arg(long, default_value_t = 20)]
        top: usize,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
        #[command(flatten)]
        scope: QueryScopeArgs,
    },
    /// Molecules with archive anchors in chrom:start-end; per-cell UMI counts to stdout.
    Region {
        /// e.g. chr6:73489308-73525587
        locus: String,
        /// Number of cell rows. With --format, 0 = all; legacy unscoped output retains its
        /// historical 0-row behavior.
        #[arg(long, default_value_t = 20)]
        top: usize,
        /// Write a sashimi-style SVG portrait of the window: per-strand molecule coverage with
        /// junction arcs weighted by support (window limited to 5 Mb).
        #[arg(long)]
        plot: Option<PathBuf>,
        /// Write IGV-ready tracks: <prefix>.plus.bedgraph, <prefix>.minus.bedgraph (molecule
        /// coverage per strand) and <prefix>.junctions.bed (BED12, TopHat-style).
        #[arg(long)]
        export_prefix: Option<PathBuf>,
        /// GTF for a gene underlay in --plot.
        #[arg(long)]
        gtf: Option<PathBuf>,
        /// Emit a header and scoped count rows; diagnostics go to stderr.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit one JSON object; diagnostics go to stderr.
        #[arg(long, conflicts_with = "tsv")]
        json: bool,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
        #[command(flatten)]
        scope: QueryScopeArgs,
    },
    /// Per-cell molecule counts supporting the junction chrom:donor-acceptor (0-based, exact).
    Junction {
        /// e.g. chr1:155234452-155235327
        locus: String,
        /// Number of per-cell rows (0 = all).
        #[arg(long, default_value_t = 20)]
        top: usize,
        /// Emit a header and per-cell rows only; diagnostics go to stderr.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit one JSON object only; diagnostics go to stderr.
        #[arg(long, conflicts_with = "tsv")]
        json: bool,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
        #[command(flatten)]
        scope: QueryScopeArgs,
    },
    /// Enumerate splice junctions whose endpoints lie in a 0-based half-open window.
    Junctions {
        /// e.g. chr11:35138870-35232402 (0-based, half-open)
        locus: String,
        /// Include a junction when either endpoint is in the window; default requires both.
        #[arg(long)]
        either: bool,
        /// Minimum index supporting-child count.
        #[arg(long, default_value_t = 1)]
        min_support: u64,
        /// Decode selected postings and report exact class-deduplicated UMI and cell counts.
        #[arg(long)]
        with_cells: bool,
        /// Minimum exact cell count; implies --with-cells when nonzero.
        #[arg(long, default_value_t = 0)]
        min_cells: usize,
        /// Mark exact and endpoint membership in this annotation (either strand in format v1).
        #[arg(long)]
        gtf: Option<PathBuf>,
        /// Emit a header and rows only; diagnostics go to stderr.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit one JSON object only; diagnostics go to stderr.
        #[arg(long, conflicts_with = "tsv")]
        json: bool,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
        #[command(flatten)]
        scope: QueryScopeArgs,
    },
    /// Class-level inclusion/exclusion junction-set usage with an explicit both-side category.
    Jset {
        /// Inclusion junction, chrom:donor-acceptor; repeat for a junction set.
        #[arg(long = "include", required = true, num_args = 1..)]
        include: Vec<String>,
        /// Exclusion junction, chrom:donor-acceptor; repeat for a junction set.
        #[arg(long = "exclude", required = true, num_args = 1..)]
        exclude: Vec<String>,
        /// Number of cell rows (0 = all); ignored for group/bulk aggregation.
        #[arg(long, default_value_t = 20)]
        top: usize,
        /// Emit a header and count rows; diagnostics go to stderr.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit one JSON object; diagnostics go to stderr.
        #[arg(long, conflicts_with = "tsv")]
        json: bool,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
        #[command(flatten)]
        scope: QueryScopeArgs,
    },
    /// Discover coordinate-defined splice events and reduce all of them in one union decode.
    Events {
        /// Genomic window containing all component junctions, e.g. chr1:45550000-45571000.
        locus: String,
        /// Event type; repeat to select several. The default is all supported types.
        #[arg(long = "event-type", value_enum)]
        event_types: Vec<EventTypeArg>,
        /// Minimum catalogue supporting-child count for every event component.
        #[arg(long, default_value_t = 2)]
        min_support: u64,
        /// Omit events with fewer scoped informative molecules after exact reduction.
        #[arg(long, default_value_t = 1)]
        min_informative: usize,
        /// Hard candidate limit. Exceeding it is an error; results are never truncated.
        #[arg(long, default_value_t = 100_000)]
        max_events: usize,
        /// Number of cell rows per event (0 = all); ignored for group/bulk aggregation.
        #[arg(long, default_value_t = 20)]
        top: usize,
        /// Optional GTF or compiled AIC used only to label coordinate-defined events.
        #[arg(long)]
        gtf: Option<PathBuf>,
        /// Emit a long TSV table; diagnostics go to stderr.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit one versioned JSON object; diagnostics go to stderr.
        #[arg(long, conflicts_with = "tsv")]
        json: bool,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
        #[command(flatten)]
        scope: QueryScopeArgs,
    },
    /// Strand-aware junction graph with molecule-class path-fragment hyperedges.
    SpliceGraph {
        /// Genomic window containing the graph, e.g. chr1:45550000-45571000.
        locus: String,
        /// Minimum catalogue supporting-child count for a junction to enter the graph.
        #[arg(long, default_value_t = 1)]
        min_support: u64,
        /// Retain exact path fragments with at least this many scoped UMI classes.
        #[arg(long, default_value_t = 1)]
        min_path_umis: usize,
        /// Hard candidate-path limit. Exceeding it is an error; paths are never truncated.
        #[arg(long, default_value_t = 100_000)]
        max_paths: usize,
        /// Emit one versioned JSON object; diagnostics go to stderr.
        #[arg(long)]
        json: bool,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
        #[command(flatten)]
        scope: QueryScopeArgs,
    },
    /// Strand-aware 3'-end site usage in a window, with per-site UMI and cell counts.
    Apa {
        locus: String,
        /// Maximum distance in bp for clustering neighboring cleavage coordinates (minimum 1).
        #[arg(long, default_value_t = 24, value_parser = clap::value_parser!(u32).range(1..))]
        site_gap: u32,
        /// Restrict to one strand: + or -.
        #[arg(long)]
        strand: Option<String>,
        /// Emit TSV rows instead of a summary.
        #[arg(long)]
        tsv: bool,
        /// Two-column TSV (barcode, group): emit per-site per-group UMI counts — differential
        /// 3' usage between cell populations, straight from the archive — plus a site×group
        /// G-test over the window.
        #[arg(long)]
        groups: Option<PathBuf>,
        /// Reference genome FASTA (plain or gzipped): enables the internal-priming filter
        /// (A-rich sequence downstream of the putative cleavage site). The contig is verified
        /// against the archive's stamped genome signature before any sequence is consulted.
        #[arg(long)]
        genome: Option<PathBuf>,
        /// With --genome: drop flagged sites from output and tests instead of only flagging them.
        #[arg(long, requires = "genome")]
        drop_ip: bool,
        /// Positive number of cell-label permutations for the G-test; omitted uses chi-square only.
        #[arg(long, requires = "groups")]
        permute: Option<std::num::NonZeroUsize>,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Write an SVG lollipop portrait of the window's 3'-sites (per-group usage band with
        /// --groups; internal-priming sites hollow with --genome).
        #[arg(long)]
        plot: Option<PathBuf>,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
    },
    /// Differential 3'-site usage per gene, genome-wide: an annotation-free site×group table per
    /// gene, a multinomial G-test, and BH-FDR across genes. The GTF supplies gene spans and ids
    /// only — sites come from the archive.
    ApaTest {
        #[arg(long)]
        gtf: PathBuf,
        /// Two-column TSV (barcode, group).
        #[arg(long)]
        groups: PathBuf,
        /// Reference genome FASTA for the internal-priming filter (verified against the stamp).
        #[arg(long)]
        genome: Option<PathBuf>,
        #[arg(long, default_value_t = 24, value_parser = clap::value_parser!(u32).range(1..))]
        site_gap: u32,
        /// Minimum grouped UMIs for a site to enter the table.
        #[arg(long, default_value_t = 5)]
        min_site_umis: usize,
        /// Minimum grouped UMIs across a gene's sites for the gene to be tested.
        #[arg(long, default_value_t = 20)]
        min_gene_umis: usize,
        /// Count molecules ending up to this far past the annotated gene 3' end.
        #[arg(long, default_value_t = 2_000)]
        tail_extend: u32,
        /// Cell-label permutations per gene (0 = chi-square approximation only).
        #[arg(long, default_value_t = 0)]
        permute: usize,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
    },
    /// Unannotated transcription: cluster molecules unclaimed by --gtf into candidate loci.
    Discover {
        #[arg(long)]
        gtf: PathBuf,
        /// Maximum gap in bp when clustering molecules into candidate loci.
        #[arg(long, default_value_t = 1_000)]
        merge_gap: u32,
        /// Minimum class-deduplicated UMIs required for a candidate locus.
        #[arg(long, default_value_t = 10)]
        min_umis: usize,
        /// Minimum UMI classes for the residual channel in `residual-sites` mode. The unchanged
        /// span channel continues to use --min-umis.
        #[arg(long, default_value_t = 10)]
        residual_min_umis: usize,
        /// Rule used to suppress evidence already explained by the supplied annotation.
        #[arg(long, value_enum, default_value_t = DiscoveryClaimMode::Span)]
        claim_mode: DiscoveryClaimMode,
        /// STARsolo alignment/transcript strand relationship. 10x 3' uses forward; R2-only 10x
        /// 5' uses reverse. Ignored by the historical `span` claiming mode.
        #[arg(long, value_enum, default_value_t = crate::archivecmd::SoloStrandArg::Forward)]
        solo_strand: crate::archivecmd::SoloStrandArg,
        /// Emit TSV rows (chrom, start, end, strand, umis, cells) instead of a summary.
        #[arg(long)]
        tsv: bool,
        /// Also write candidates as a GTF (single-exon genes, AIENOVEL ids) — the
        /// discover→refine→replay loop: this file feeds straight back into `replay-rows`.
        #[arg(long)]
        emit_gtf: Option<PathBuf>,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
    },
}

pub(crate) fn parse_locus(s: &str) -> Result<(String, u32, u32)> {
    let (chrom, range) = s.split_once(':').context("locus must be chrom:start-end")?;
    let (a, b) = range.split_once('-').context("locus must be chrom:start-end")?;
    let start = a.replace(',', "").parse()?;
    let end = b.replace(',', "").parse()?;
    if start >= end {
        bail!("locus start must be smaller than end (0-based half-open coordinates)");
    }
    Ok((chrom.to_string(), start, end))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum QueryAggregation {
    Cell,
    Group,
    Bulk,
}

struct QueryScope {
    selected: Option<FxHashSet<u32>>,
    group_names: Vec<String>,
    group_of: FxHashMap<u32, u32>,
    selected_per_group: Vec<usize>,
    selected_cells: usize,
    aggregation: QueryAggregation,
    source: &'static str,
    source_path: Option<PathBuf>,
    source_content_blake3: Option<String>,
    resolved_mapping_blake3: Option<String>,
    archive_cells: usize,
    active: bool,
}

impl QueryScope {
    fn includes(&self, cell: u32) -> bool {
        self.selected
            .as_ref()
            .is_none_or(|selected| selected.contains(&cell))
    }

    fn json(&self) -> serde_json::Value {
        let aggregation = self.aggregation_name();
        let mut value = json!({
            "source": self.source,
            "aggregation": aggregation,
            "selected_cells": self.selected_cells,
        });
        if !self.group_names.is_empty() {
            value.as_object_mut().unwrap().insert(
                "groups".into(),
                json!(self
                    .group_names
                    .iter()
                    .enumerate()
                    .map(|(index, name)| json!({
                        "name": name,
                        "selected_cells": self.selected_per_group[index],
                    }))
                    .collect::<Vec<_>>()),
            );
        }
        value
    }

    fn aggregation_name(&self) -> &'static str {
        match self.aggregation {
            QueryAggregation::Cell => "cell",
            QueryAggregation::Group => "group",
            QueryAggregation::Bulk => "bulk",
        }
    }

    /// Complete scope provenance for uniform output. The raw file digest binds the exact user
    /// input while the canonical resolution digest binds the archive cell/group population it
    /// selected, independent of input line order.
    fn provenance_json(&self) -> serde_json::Value {
        let mut value = self.json();
        let object = value.as_object_mut().expect("scope JSON is an object");
        if let Some(path) = &self.source_path {
            object.insert("source_path".into(), json!(path));
        }
        if let Some(digest) = &self.source_content_blake3 {
            object.insert("source_content_blake3".into(), json!(digest));
        }
        if let Some(digest) = &self.resolved_mapping_blake3 {
            object.insert("resolved_mapping_blake3".into(), json!(digest));
        }
        value
    }

    /// Compute canonical resolved-population provenance only for the opt-in uniform path. Legacy
    /// scoped commands therefore do not acquire a second dense dictionary traversal.
    fn ensure_resolved_mapping_digest(&mut self) -> Result<()> {
        if !self.active || self.resolved_mapping_blake3.is_some() {
            return Ok(());
        }
        let dictionary_len = u32::try_from(self.archive_cells)
            .context("archive cell dictionary exceeds the u32 cell-id domain")?;
        let Some(selected) = self.selected.as_ref() else {
            // `--agg cell|bulk` without a scope file activates aggregation metadata but still
            // denotes the archive's complete population, already bound by archive identity.
            return Ok(());
        };
        let mut hasher = blake3::Hasher::new();
        hasher.update(b"gravlax-query-resolved-scope-v1\0");
        if selected.len().saturating_mul(4) < self.archive_cells {
            // Sparse scope: sorting k compact ids avoids a pathological scan of a huge archive
            // dictionary. This allocation exists only while uniform provenance is constructed.
            let mut cells: Vec<u32> = selected.iter().copied().collect();
            cells.sort_unstable();
            for cell in cells {
                update_resolved_scope_hasher(&mut hasher, cell, &self.group_of, &self.group_names);
            }
        } else {
            // Dense scope: dictionary-order traversal is already canonical and avoids cloning a
            // result-sized set merely to sort it.
            for cell in 0..dictionary_len {
                if selected.contains(&cell) {
                    update_resolved_scope_hasher(
                        &mut hasher,
                        cell,
                        &self.group_of,
                        &self.group_names,
                    );
                }
            }
        }
        self.resolved_mapping_blake3 = Some(format!("blake3:{}", hasher.finalize().to_hex()));
        Ok(())
    }
}

fn update_resolved_scope_hasher(
    hasher: &mut blake3::Hasher,
    cell: u32,
    group_of: &FxHashMap<u32, u32>,
    group_names: &[String],
) {
    hasher.update(&cell.to_le_bytes());
    match group_of.get(&cell) {
        Some(group) => {
            let group = &group_names[*group as usize];
            hasher.update(&(group.len() as u64).to_le_bytes());
            hasher.update(group.as_bytes());
        }
        None => {
            hasher.update(&0_u64.to_le_bytes());
        }
    }
}

fn archive_cell_ids(la: &mut LazyArchive) -> Result<(Vec<u32>, FxHashMap<u32, u32>)> {
    let dictionary = la.cells()?.to_vec();
    let ids = dictionary
        .iter()
        .enumerate()
        .map(|(index, &packed)| (packed, index as u32))
        .collect();
    Ok((dictionary, ids))
}

fn load_query_scope(la: &mut LazyArchive, args: &QueryScopeArgs) -> Result<QueryScope> {
    let (dictionary, packed_to_id) = archive_cell_ids(la)?;
    load_query_scope_from_dictionary(dictionary.len(), &packed_to_id, args)
}

fn load_transcript_ec_scope(
    cells: &[crate::archivecmd::transcriptec::TranscriptCellSummary],
    args: &QueryScopeArgs,
) -> Result<QueryScope> {
    let mut packed_to_id = FxHashMap::default();
    for cell in cells {
        let packed = umi::pack(cell.barcode.as_bytes()).with_context(|| {
            format!(
                "archive cell {} has an invalid stored barcode {}",
                cell.cell_id, cell.barcode
            )
        })?;
        if packed_to_id.insert(packed, cell.cell_id).is_some() {
            bail!(
                "archive cell dictionary contains duplicate barcode {}",
                cell.barcode
            );
        }
    }
    load_query_scope_from_dictionary(cells.len(), &packed_to_id, args)
}

fn load_query_scope_from_dictionary(
    dictionary_len: usize,
    packed_to_id: &FxHashMap<u32, u32>,
    args: &QueryScopeArgs,
) -> Result<QueryScope> {
    let mut selected: Option<FxHashSet<u32>> = None;
    let mut group_names = Vec::new();
    let mut group_of = FxHashMap::default();
    let mut source_path = None;
    let mut source_content_blake3 = None;
    let source;

    if let Some(path) = &args.cells {
        source = "cells";
        let (text, digest) = read_query_scope_file(path, "cell")?;
        source_path = Some(path.clone());
        source_content_blake3 = Some(digest);
        let mut ids = FxHashSet::default();
        for (line_index, raw) in text.lines().enumerate() {
            let line_no = line_index + 1;
            let line = raw.trim_end_matches('\r');
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if line.trim() != line || line.chars().any(char::is_whitespace) {
                bail!("cell scope line {line_no} must contain exactly one barcode");
            }
            let packed = umi::pack(line.as_bytes())
                .with_context(|| format!("cell scope line {line_no} has an invalid barcode"))?;
            let cell = packed_to_id.get(&packed).copied().with_context(|| {
                format!("cell scope line {line_no}: barcode {line} is not in the archive")
            })?;
            if !ids.insert(cell) {
                bail!("cell scope contains duplicate barcode {line}");
            }
        }
        if ids.is_empty() {
            bail!("cell scope contains no archive barcodes");
        }
        selected = Some(ids);
    } else if let Some(path) = &args.groups {
        source = "groups";
        let (text, digest) = read_query_scope_file(path, "group")?;
        source_path = Some(path.clone());
        source_content_blake3 = Some(digest);
        let mut ids = FxHashSet::default();
        for (line_index, raw) in text.lines().enumerate() {
            let line_no = line_index + 1;
            let line = raw.trim_end_matches('\r');
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 2 || fields[0].is_empty() || fields[1].is_empty() {
                bail!("group scope line {line_no} must be barcode<TAB>group");
            }
            if fields[0].trim() != fields[0]
                || fields[0].chars().any(char::is_whitespace)
                || fields[1].trim() != fields[1]
            {
                bail!("group scope line {line_no} has invalid surrounding whitespace");
            }
            let packed = umi::pack(fields[0].as_bytes())
                .with_context(|| format!("group scope line {line_no} has an invalid barcode"))?;
            let cell = packed_to_id.get(&packed).copied().with_context(|| {
                format!(
                    "group scope line {line_no}: barcode {} is not in the archive",
                    fields[0]
                )
            })?;
            if !ids.insert(cell) {
                bail!("group scope contains duplicate barcode {}", fields[0]);
            }
            let group = match group_names.iter().position(|name| name == fields[1]) {
                Some(index) => index as u32,
                None => {
                    group_names.push(fields[1].to_owned());
                    (group_names.len() - 1) as u32
                }
            };
            group_of.insert(cell, group);
        }
        if ids.is_empty() {
            bail!("group scope contains no archive barcodes");
        }
        selected = Some(ids);
    } else {
        source = "all";
    }

    let aggregation = match args.agg {
        QueryAggregationArg::Auto => {
            if args.groups.is_some() {
                QueryAggregation::Group
            } else {
                QueryAggregation::Cell
            }
        }
        QueryAggregationArg::Cell => QueryAggregation::Cell,
        QueryAggregationArg::Group => QueryAggregation::Group,
        QueryAggregationArg::Bulk => QueryAggregation::Bulk,
    };
    if aggregation == QueryAggregation::Group && args.groups.is_none() {
        bail!("--agg group requires --groups");
    }
    let selected_cells = selected.as_ref().map_or(dictionary_len, FxHashSet::len);
    let mut selected_per_group = vec![0usize; group_names.len()];
    for &group in group_of.values() {
        selected_per_group[group as usize] += 1;
    }
    Ok(QueryScope {
        selected,
        group_names,
        group_of,
        selected_per_group,
        selected_cells,
        aggregation,
        source,
        source_path,
        source_content_blake3,
        resolved_mapping_blake3: None,
        archive_cells: dictionary_len,
        active: args.cells.is_some()
            || args.groups.is_some()
            || args.agg != QueryAggregationArg::Auto,
    })
}

fn read_query_scope_file(path: &Path, kind: &str) -> Result<(String, String)> {
    let bytes = std::fs::read(path)
        .with_context(|| format!("reading {kind} scope {}", path.display()))?;
    let digest = format!("blake3:{}", blake3::hash(&bytes).to_hex());
    let text = String::from_utf8(bytes)
        .with_context(|| format!("{kind} scope {} is not UTF-8", path.display()))?;
    Ok((text, digest))
}

#[cfg(unix)]
fn same_query_input_snapshot(before: &std::fs::Metadata, after: &std::fs::Metadata) -> bool {
    use std::os::unix::fs::MetadataExt;

    before.dev() == after.dev()
        && before.ino() == after.ino()
        && before.len() == after.len()
        && before.mtime() == after.mtime()
        && before.mtime_nsec() == after.mtime_nsec()
        && before.ctime() == after.ctime()
        && before.ctime_nsec() == after.ctime_nsec()
}

#[cfg(not(unix))]
fn same_query_input_snapshot(before: &std::fs::Metadata, after: &std::fs::Metadata) -> bool {
    before.len() == after.len() && before.modified().ok() == after.modified().ok()
}

/// UTF-8 input and its canonical digest, both captured from one stable open file description.
/// This is reserved for uniform provenance paths so legacy parsing and performance stay intact.
struct BoundQueryText {
    text: String,
    content_blake3: String,
}

fn read_bound_query_text(path: &Path, label: &str) -> Result<BoundQueryText> {
    let mut file = std::fs::File::open(path)
        .with_context(|| format!("opening {label} {}", path.display()))?;
    let before = file
        .metadata()
        .with_context(|| format!("inspecting {label} {}", path.display()))?;
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)
        .with_context(|| format!("reading {label} {}", path.display()))?;
    let after = file
        .metadata()
        .with_context(|| format!("re-inspecting {label} {}", path.display()))?;
    if !same_query_input_snapshot(&before, &after) {
        bail!(
            "{label} {} changed while its content and provenance were being loaded",
            path.display()
        );
    }
    let content_blake3 = format!("blake3:{}", blake3::hash(&bytes).to_hex());
    let text = String::from_utf8(bytes)
        .with_context(|| format!("{label} {} is not UTF-8", path.display()))?;
    Ok(BoundQueryText {
        text,
        content_blake3,
    })
}

/// Parse an annotation and obtain its digest from the same guarded descriptor for uniform
/// provenance. The labels here identify content only; they do not assert an assembly/release.
fn load_query_annotation(
    path: &Path,
    label: &str,
    bind_provenance: bool,
) -> Result<(anno::Annotation, Option<String>)> {
    if !bind_provenance {
        return Ok((anno::Annotation::from_path(path)?, None));
    }
    let identity = anno::intent::AnnotationIdentity::new("unspecified", label)
        .context("constructing query annotation content identity")?;
    let bound = anno::intent::BoundAnnotation::from_path(path, identity)
        .with_context(|| format!("loading {label} {}", path.display()))?;
    let (annotation, identity) = bound.into_parts();
    let digest = identity
        .digest
        .context("bound query annotation did not report its canonical content digest")?;
    Ok((annotation, Some(digest)))
}

struct ScopedCounts {
    cells: Vec<(u32, usize)>,
    groups: Vec<(usize, usize)>,
    total_umis: usize,
}

#[derive(Serialize)]
struct RegionUniformSummary<'a> {
    coordinates: &'static str,
    anchor_semantics: bool,
    chrom: &'a str,
    start: u32,
    end: u32,
    molecules: u64,
    umis: u64,
    cells: u64,
    chunks_decoded: u64,
}

#[derive(Serialize)]
struct JunctionUniformSummary<'a> {
    coordinates: &'static str,
    chrom: &'a str,
    donor: u32,
    acceptor: u32,
    archive_supporting_children: u64,
    archive_posting_chunks: u64,
    umis: u64,
    cells: u64,
}

#[derive(Serialize)]
struct UniformSelectionPolicy {
    requested_top: usize,
    top_zero_means_all: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    comparator: Option<&'static str>,
}

const REGION_UNIFORM_RESULT_SCHEMA: &str = "gravlax.query.region.result.v1";
const REGION_UNIFORM_COUNTS_SCHEMA: &str = "gravlax.query.region.counts.v1";
const JUNCTION_UNIFORM_RESULT_SCHEMA: &str = "gravlax.query.junction.result.v1";
const JUNCTION_UNIFORM_COUNTS_SCHEMA: &str = "gravlax.query.junction.counts.v1";

fn uniform_count_schema(id: &'static str) -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(
        id,
        vec![
            Field::new("aggregation", DataType::String),
            Field::new("entity", DataType::String),
            Field::new("umis", DataType::UInt64),
            Field::new("cells", DataType::UInt64).nullable(),
            Field::new("selected_cells", DataType::UInt64).nullable(),
        ],
    )?
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["aggregation", "entity"]))
}

fn unpack_cell_bytes(packed: u32) -> [u8; 16] {
    let mut packed = packed;
    let mut barcode = [0_u8; 16];
    for index in (0..barcode.len()).rev() {
        barcode[index] = b"ACGT"[(packed & 0b11) as usize];
        packed >>= 2;
    }
    barcode
}

#[allow(clippy::too_many_arguments)]
fn stream_uniform_counts<W, S>(
    writer: W,
    result_schema: &'static str,
    table_schema: &'static str,
    format: UniformQueryFormat,
    context: &ResultContext,
    summary: &S,
    counts: &ScopedCounts,
    scope: &QueryScope,
    cells: &[u32],
    emitted_rows: usize,
) -> std::result::Result<W, OutputError>
where
    W: Write,
    S: Serialize,
{
    let schema = uniform_count_schema(table_schema)?;
    let (available_rows, expected_emitted) = uniform_row_limit(counts, scope, usize::MAX);
    let available_rows = u64::try_from(available_rows)
        .map_err(|_| OutputError::InvalidSchema("available row count exceeds u64".into()))?;
    let emitted_rows_u64 = u64::try_from(emitted_rows)
        .map_err(|_| OutputError::InvalidSchema("emitted row count exceeds u64".into()))?;
    let selection = SelectionSummary::selected(available_rows, emitted_rows_u64)?;
    if expected_emitted < emitted_rows {
        return Err(OutputError::InvalidSchema(
            "uniform query row limit exceeds available rows".into(),
        ));
    }
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        result_schema,
        OutputFormat::from(format),
        context,
        summary,
    )?;
    bundle.write_table("counts", &schema, Some(&selection), |rows| {
        match scope.aggregation {
            QueryAggregation::Cell => {
                for (cell, umis) in counts.cells.iter().take(emitted_rows) {
                    let barcode = unpack_cell_bytes(cells[*cell as usize]);
                    let barcode = std::str::from_utf8(&barcode)
                        .expect("packed archive barcode decodes to ASCII");
                    rows.write_row_with(|row| {
                        row.string("cell")?;
                        row.string(barcode)?;
                        row.uint64(*umis as u64)?;
                        row.null()?;
                        row.null()?;
                        Ok(())
                    })?;
                }
            }
            QueryAggregation::Group => {
                for (group, (umis, contributing_cells)) in counts.groups.iter().enumerate() {
                    rows.write_row_with(|row| {
                        row.string("group")?;
                        row.string(&scope.group_names[group])?;
                        row.uint64(*umis as u64)?;
                        row.uint64(*contributing_cells as u64)?;
                        row.uint64(scope.selected_per_group[group] as u64)?;
                        Ok(())
                    })?;
                }
            }
            QueryAggregation::Bulk => {
                rows.write_row_with(|row| {
                    row.string("bulk")?;
                    row.string("bulk")?;
                    row.uint64(counts.total_umis as u64)?;
                    row.uint64(counts.cells.len() as u64)?;
                    row.uint64(scope.selected_cells as u64)?;
                    Ok(())
                })?;
            }
        }
        Ok(())
    })?;
    bundle.finish()
}

#[allow(clippy::too_many_arguments)]
fn write_uniform_count_output<S: Serialize>(
    output: &UniformQueryOutputArgs,
    result_schema: &'static str,
    table_schema: &'static str,
    context: &ResultContext,
    summary: &S,
    counts: &ScopedCounts,
    scope: &QueryScope,
    cells: &[u32],
    emitted_rows: usize,
) -> Result<()> {
    let format = output
        .format
        .context("uniform output requires an explicit --format")?;
    if let Some(path) = output.output.as_deref() {
        let outcome = publish_file_no_clobber(path, Durability::Flush, |writer| {
            stream_uniform_counts(
                &mut *writer,
                result_schema,
                table_schema,
                format,
                context,
                summary,
                counts,
                scope,
                cells,
                emitted_rows,
            )?;
            Ok(())
        })?;
        for warning in outcome.warnings {
            eprintln!("warning: {warning}");
        }
    } else {
        let stdout = std::io::stdout();
        let writer = std::io::BufWriter::new(stdout.lock());
        stream_uniform_counts(
            writer,
            result_schema,
            table_schema,
            format,
            context,
            summary,
            counts,
            scope,
            cells,
            emitted_rows,
        )?;
    }
    Ok(())
}

fn uniform_row_limit(counts: &ScopedCounts, scope: &QueryScope, top: usize) -> (usize, usize) {
    match scope.aggregation {
        QueryAggregation::Cell => {
            let available = counts.cells.len();
            (
                available,
                if top == 0 {
                    available
                } else {
                    top.min(available)
                },
            )
        }
        QueryAggregation::Group => (counts.groups.len(), counts.groups.len()),
        QueryAggregation::Bulk => (1, 1),
    }
}

fn uniform_query_parameters(
    scope: &QueryScope,
    top: usize,
    archive_access: &'static str,
) -> Result<BTreeMap<String, serde_json::Value>> {
    let mut parameters = BTreeMap::new();
    parameters.insert("cell_scope".into(), scope.provenance_json());
    parameters.insert("aggregation".into(), json!(scope.aggregation_name()));
    parameters.insert("archive_access".into(), json!(archive_access));
    parameters.insert(
        "selection_policy".into(),
        serde_json::to_value(UniformSelectionPolicy {
            requested_top: top,
            top_zero_means_all: true,
            comparator: matches!(scope.aggregation, QueryAggregation::Cell)
                .then_some("umis descending, entity ascending (barcode)"),
        })?,
    );
    Ok(parameters)
}

fn scoped_counts(per_cell: &FxHashMap<u32, FxHashSet<u32>>, scope: &QueryScope) -> ScopedCounts {
    scoped_counts_with_cell_order(per_cell, scope, None)
}

fn scoped_counts_with_cell_order(
    per_cell: &FxHashMap<u32, FxHashSet<u32>>,
    scope: &QueryScope,
    packed_barcodes: Option<&[u32]>,
) -> ScopedCounts {
    let mut cells: Vec<(u32, usize)> = per_cell
        .iter()
        .filter(|(cell, _)| scope.includes(**cell))
        .map(|(cell, classes)| (*cell, classes.len()))
        .collect();
    sort_scoped_cell_counts(&mut cells, packed_barcodes);
    let total_umis = cells.iter().map(|(_, count)| count).sum();
    let mut groups = vec![(0usize, 0usize); scope.group_names.len()];
    for (cell, count) in &cells {
        if let Some(&group) = scope.group_of.get(cell) {
            groups[group as usize].0 += *count;
            groups[group as usize].1 += 1;
        }
    }
    ScopedCounts {
        cells,
        groups,
        total_umis,
    }
}

fn scoped_counts_unordered(
    per_cell: &FxHashMap<u32, FxHashSet<u32>>,
    scope: &QueryScope,
) -> ScopedCounts {
    let cells: Vec<(u32, usize)> = per_cell
        .iter()
        .filter(|(cell, _)| scope.includes(**cell))
        .map(|(cell, classes)| (*cell, classes.len()))
        .collect();
    let total_umis = cells.iter().map(|(_, count)| count).sum();
    let mut groups = vec![(0usize, 0usize); scope.group_names.len()];
    for (cell, count) in &cells {
        if let Some(&group) = scope.group_of.get(cell) {
            groups[group as usize].0 += count;
            groups[group as usize].1 += 1;
        }
    }
    ScopedCounts {
        cells,
        groups,
        total_umis,
    }
}

fn sort_scoped_cell_counts(cells: &mut [(u32, usize)], packed_barcodes: Option<&[u32]>) {
    match packed_barcodes {
        Some(barcodes) => cells.sort_unstable_by_key(|(cell, count)| {
            (std::cmp::Reverse(*count), barcodes[*cell as usize])
        }),
        None => cells.sort_unstable_by_key(|(cell, count)| (std::cmp::Reverse(*count), *cell)),
    }
}

fn uniform_query_context(
    archive: &Path,
    archive_version: u32,
    archive_root: Option<String>,
    mut parameters: BTreeMap<String, serde_json::Value>,
) -> ResultContext {
    parameters.insert("archive_path".into(), json!(archive));
    parameters.insert("archive_version".into(), json!(archive_version));
    let mut warnings = Vec::new();
    let archives = match archive_root {
        Some(root) => vec![format!("aie-directory-root-v2:{root}")],
        None => {
            warnings.push(
                "legacy v1 archive has no rooted content commitment; its path locator is not a portable content identity"
                    .into(),
            );
            vec![format!(
                "aie-v{archive_version}-unrooted:{}",
                archive.display()
            )]
        }
    };
    ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            archives,
            parameters,
            ..Default::default()
        },
        warnings,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BatchKind {
    Region,
    Junction,
}

struct BatchSpec {
    id: String,
    kind: BatchKind,
    chrom: String,
    chrom_id: u32,
    start: u32,
    end: u32,
}

fn parse_batch_plan(text: &str, chrom_names: &[String]) -> Result<Vec<BatchSpec>> {
    let mut lines = text.lines();
    let header = lines.next().unwrap_or_default().trim_end_matches('\r');
    if header != "id\tkind\tlocus" {
        bail!("batch plan header must be exactly: id<TAB>kind<TAB>locus");
    }
    let mut ids = FxHashSet::default();
    let mut out = Vec::new();
    for (line_index, raw) in lines.enumerate() {
        let line_no = line_index + 2;
        let line = raw.trim_end_matches('\r');
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 3 {
            bail!("batch plan line {line_no} must contain exactly three tab-separated fields");
        }
        let id = fields[0];
        if id.is_empty() || id.chars().any(char::is_whitespace) {
            bail!("batch plan line {line_no} has an empty or whitespace-containing id");
        }
        if !ids.insert(id.to_owned()) {
            bail!("batch plan contains duplicate query id {id}");
        }
        let kind = match fields[1] {
            "region" => BatchKind::Region,
            "junction" => BatchKind::Junction,
            other => bail!("batch plan line {line_no} has unsupported query kind {other}"),
        };
        let (chrom, start, end) =
            parse_locus(fields[2]).with_context(|| format!("batch plan line {line_no}"))?;
        let chrom_id = chrom_names
            .iter()
            .position(|name| name == &chrom)
            .map(|id| id as u32)
            .with_context(|| format!("batch plan line {line_no}: unknown chromosome {chrom}"))?;
        out.push(BatchSpec {
            id: id.to_owned(),
            kind,
            chrom,
            chrom_id,
            start,
            end,
        });
        if out.len() > 100_000 {
            bail!("batch plan exceeds the 100000-query safety limit");
        }
    }
    if out.is_empty() {
        bail!("batch plan contains no queries");
    }
    Ok(out)
}

fn read_batch_plan(path: &std::path::Path, chrom_names: &[String]) -> Result<Vec<BatchSpec>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading batch query plan {}", path.display()))?;
    parse_batch_plan(&text, chrom_names)
}

pub(crate) fn region_selects_chunk(
    chunk: &crate::archivecmd::ChunkInfo,
    chrom: u32,
    start: u32,
    end: u32,
) -> bool {
    chunk.chrom == chrom
        && chunk.bin_start < end
        && !((chunk.n_mols == 0 || chunk.max_anchor < start)
            && chunk.bin_start.saturating_add(8_000_000) < start)
}

struct BatchRun<'a> {
    archive: &'a Path,
    la: &'a mut LazyArchive,
    chunks: &'a [crate::archivecmd::ChunkInfo],
    chrom_names: &'a [String],
    plan_path: &'a Path,
    top: usize,
    uniform_output: &'a UniformQueryOutputArgs,
    scope_args: &'a QueryScopeArgs,
    t0: std::time::Instant,
    t_open: f32,
}

fn run_batch(args: BatchRun<'_>) -> Result<()> {
    let BatchRun {
        archive,
        la,
        chunks,
        chrom_names,
        plan_path,
        top,
        uniform_output,
        scope_args,
        t0,
        t_open,
    } = args;
    let mut scope = load_query_scope(la, scope_args)?;
    if uniform_output.format.is_some() {
        scope.ensure_resolved_mapping_digest()?;
    }
    let (specs, plan_content_blake3) = if uniform_output.format.is_some() {
        let bound = read_bound_query_text(plan_path, "batch plan")?;
        (
            parse_batch_plan(&bound.text, chrom_names)?,
            Some(bound.content_blake3),
        )
    } else {
        (read_batch_plan(plan_path, chrom_names)?, None)
    };
    let n_regions = specs
        .iter()
        .filter(|query| query.kind == BatchKind::Region)
        .count();
    let n_junctions = specs.len() - n_regions;
    let metadata = if n_junctions > 0 {
        read_junction_metadata(la)?
    } else {
        Vec::new()
    };
    let metadata_of: FxHashMap<(u32, u32, u32), usize> = metadata
        .iter()
        .enumerate()
        .map(|(index, row)| ((row.chrom, row.donor, row.acceptor), index))
        .collect();
    let mut junction_info: Vec<Option<(u64, usize)>> = vec![None; specs.len()];
    let mut chunk_queries: Vec<Vec<usize>> = (0..chunks.len()).map(|_| Vec::new()).collect();
    let mut independent_chunk_decodes = 0usize;
    for (query_index, query) in specs.iter().enumerate() {
        match query.kind {
            BatchKind::Region => {
                for (chunk_index, chunk) in chunks.iter().enumerate() {
                    if region_selects_chunk(chunk, query.chrom_id, query.start, query.end) {
                        chunk_queries[chunk_index].push(query_index);
                        independent_chunk_decodes += 1;
                    }
                }
            }
            BatchKind::Junction => {
                if let Some(&metadata_index) =
                    metadata_of.get(&(query.chrom_id, query.start, query.end))
                {
                    let row = &metadata[metadata_index];
                    junction_info[query_index] = Some((row.supporting_children, row.posts.len()));
                    for &post in &row.posts {
                        let tasks = chunk_queries.get_mut(post as usize).with_context(|| {
                            format!("junction {} references missing chunk {post}", query.id)
                        })?;
                        tasks.push(query_index);
                        independent_chunk_decodes += 1;
                    }
                }
            }
        }
    }
    for tasks in &mut chunk_queries {
        tasks.sort_unstable();
        tasks.dedup();
    }
    let selected: Vec<usize> = chunk_queries
        .iter()
        .enumerate()
        .filter_map(|(index, tasks)| (!tasks.is_empty()).then_some(index))
        .collect();
    let unique_chunk_decodes = selected.len();
    let shapes = (n_junctions > 0).then(|| la.shapes()).transpose()?;

    struct ChunkBatchHits {
        hits: Vec<(usize, u32)>,
        region_molecules: Vec<(usize, u64)>,
        scoped_region_molecules: Vec<(usize, u32)>,
    }
    let chunk_hits: Vec<ChunkBatchHits> = {
        let (reader, tables) = la.reader_and_tables();
        let reader = &*reader;
        selected
            .par_iter()
            .map(|&chunk_index| -> Result<ChunkBatchHits> {
                let (compressed, raw_len) =
                    reader.read_compressed_at(&format!("c{chunk_index}"))?;
                let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                let molecules = decode_chunk(&raw, &chunks[chunk_index], None, tables)?;
                let tasks = &chunk_queries[chunk_index];
                let region_tasks: Vec<usize> = tasks
                    .iter()
                    .copied()
                    .filter(|&query| specs[query].kind == BatchKind::Region)
                    .collect();
                let mut wanted_junctions: FxHashMap<(u32, u32), Vec<usize>> =
                    FxHashMap::default();
                for &query in tasks {
                    if specs[query].kind == BatchKind::Junction
                        && junction_info[query].is_some()
                    {
                        wanted_junctions
                            .entry((specs[query].start, specs[query].end))
                            .or_default()
                            .push(query);
                    }
                }
                let mut out = ChunkBatchHits {
                    hits: Vec::new(),
                    region_molecules: region_tasks
                        .iter()
                        .copied()
                        .map(|query| (query, 0))
                        .collect(),
                    scoped_region_molecules: Vec::new(),
                };
                for molecule in &molecules {
                    let anchor = molecule.anchor();
                    for (local_query, &query) in region_tasks.iter().enumerate() {
                        let spec = &specs[query];
                        if anchor >= spec.start && anchor < spec.end {
                            out.region_molecules[local_query].1 += 1;
                            out.hits.push((query, molecule.umi_class));
                            if scope.active {
                                out.scoped_region_molecules
                                    .push((query, molecule.umi_class));
                            }
                        }
                    }
                    if wanted_junctions.is_empty() {
                        continue;
                    }
                    let shapes = shapes.as_ref().unwrap();
                    let mut seen: FxHashSet<usize> = FxHashSet::default();
                    let mut inspect = |position: u32, shape: u32| {
                        for blocks in shapes[shape as usize].blocks.windows(2) {
                            let donor = position + blocks[0].0 + blocks[0].1;
                            let acceptor = position + blocks[1].0;
                            if let Some(queries) =
                                wanted_junctions.get(&(donor, acceptor))
                            {
                                seen.extend(queries);
                            }
                        }
                    };
                    for chain in &molecule.chains {
                        for &(position, shape) in &chain.reps {
                            inspect(position, shape);
                        }
                    }
                    for &(position, shape, _, _) in &molecule.mms {
                        inspect(position, shape);
                    }
                    out.hits.extend(seen.into_iter().map(|query| (query, molecule.umi_class)));
                }
                out.hits.sort_unstable();
                out.hits.dedup();
                Ok(out)
            })
            .collect::<Result<_>>()?
    };

    la.prefetch_coc(chunk_hits.iter().flat_map(|chunk| {
        chunk
            .hits
            .iter()
            .map(|row| row.1)
            .chain(chunk.scoped_region_molecules.iter().map(|row| row.1))
    }))?;
    let mut per_query: Vec<FxHashMap<u32, FxHashSet<u32>>> =
        (0..specs.len()).map(|_| FxHashMap::default()).collect();
    let mut region_molecules = vec![0u64; specs.len()];
    let mut scoped_region_molecules = vec![0u64; specs.len()];
    for chunk in chunk_hits {
        for (query, count) in chunk.region_molecules {
            region_molecules[query] += count;
        }
        for (query, class) in chunk.scoped_region_molecules {
            if scope.includes(la.cell_of(class)?) {
                scoped_region_molecules[query] += 1;
            }
        }
        for (query, class) in chunk.hits {
            let cell = la.cell_of(class)?;
            if scope.includes(cell) {
                per_query[query].entry(cell).or_default().insert(class);
            }
        }
    }
    let cell_dictionary = la.cells()?.to_vec();
    let uniform = uniform_output.format.is_some();
    let mut results = Vec::with_capacity(specs.len());
    let mut uniform_counts = Vec::with_capacity(specs.len());
    for (query_index, query) in specs.iter().enumerate() {
        let counts = scoped_counts_with_cell_order(
            &per_query[query_index],
            &scope,
            uniform.then_some(cell_dictionary.as_slice()),
        );
        if uniform {
            uniform_counts.push(counts);
            continue;
        }
        let limit = if top == 0 {
            counts.cells.len()
        } else {
            top.min(counts.cells.len())
        };
        let mut result = json!({
            "id": query.id,
            "kind": match query.kind { BatchKind::Region => "region", BatchKind::Junction => "junction" },
            "chrom": query.chrom,
            "start": query.start,
            "end": query.end,
            "umis": counts.total_umis,
            "cells": counts.cells.len(),
        });
        let object = result.as_object_mut().unwrap();
        if !scope.active {
            object.insert(
                "cell_rows".into(),
                json!(counts
                    .cells
                    .iter()
                    .take(limit)
                    .map(|(cell, count)| json!({
                        "barcode": unpack_cell(&cell_dictionary, *cell),
                        "umis": count,
                    }))
                    .collect::<Vec<_>>()),
            );
            object.insert(
                "cell_rows_truncated".into(),
                json!(limit < counts.cells.len()),
            );
        } else {
            object.insert("scope".into(), scope.json());
            match scope.aggregation {
                QueryAggregation::Cell => {
                    object.insert(
                        "cell_rows".into(),
                        json!(counts
                            .cells
                            .iter()
                            .take(limit)
                            .map(|(cell, count)| json!({
                                "barcode": unpack_cell(&cell_dictionary, *cell),
                                "umis": count,
                            }))
                            .collect::<Vec<_>>()),
                    );
                    object.insert(
                        "cell_rows_truncated".into(),
                        json!(limit < counts.cells.len()),
                    );
                }
                QueryAggregation::Group => {
                    object.insert(
                        "group_rows".into(),
                        json!(counts
                            .groups
                            .iter()
                            .enumerate()
                            .map(|(group, (umis, cells))| json!({
                                "group": scope.group_names[group],
                                "umis": umis,
                                "cells": cells,
                                "selected_cells": scope.selected_per_group[group],
                            }))
                            .collect::<Vec<_>>()),
                    );
                }
                QueryAggregation::Bulk => {
                    object.insert(
                        "bulk".into(),
                        json!({"umis": counts.total_umis, "cells": counts.cells.len()}),
                    );
                }
            }
        }
        match query.kind {
            BatchKind::Region => {
                object.insert(
                    "molecules".into(),
                    json!(if scope.active {
                        scoped_region_molecules[query_index]
                    } else {
                        region_molecules[query_index]
                    }),
                );
                object.insert("anchor_semantics".into(), json!(true));
            }
            BatchKind::Junction => {
                object.insert(
                    "present".into(),
                    json!(junction_info[query_index].is_some()),
                );
                if let Some((support, postings)) = junction_info[query_index] {
                    object.insert("supporting_children".into(), json!(support));
                    object.insert("posting_chunks".into(), json!(postings));
                }
            }
        }
        results.push(result);
    }
    let reduction = if independent_chunk_decodes == 0 {
        0.0
    } else {
        1.0 - unique_chunk_decodes as f64 / independent_chunk_decodes as f64
    };
    if uniform {
        let summary = json!({
            "coordinates": "0-based half-open",
            "query_order": "plan",
            "scope": scope.json(),
            "plan_queries": specs.len(),
            "region_queries": n_regions,
            "junction_queries": n_junctions,
            "planning": {
                "independent_chunk_decodes": independent_chunk_decodes,
                "unique_chunk_decodes": unique_chunk_decodes,
                "chunk_decode_reduction_fraction": reduction,
            },
        });
        let query_schema = TableSchema::new(
            "gravlax.query.batch.queries.v1",
            vec![
                Field::new("plan_index", DataType::UInt64),
                Field::new("id", DataType::String),
                Field::new("kind", DataType::String),
                Field::new("chrom", DataType::String),
                Field::new("start", DataType::UInt64),
                Field::new("end", DataType::UInt64),
                Field::new("present", DataType::Boolean),
                Field::new("archive_supporting_children", DataType::UInt64).nullable(),
                Field::new("archive_posting_chunks", DataType::UInt64).nullable(),
                Field::new("molecules", DataType::UInt64).nullable(),
                Field::new("anchor_semantics", DataType::Boolean).nullable(),
                Field::new("umis", DataType::UInt64),
                Field::new("cells", DataType::UInt64),
                Field::new("available_count_rows", DataType::UInt64),
                Field::new("emitted_count_rows", DataType::UInt64),
                Field::new("count_rows_truncated", DataType::Boolean),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["id"])
                .ordered_by([OrderKey::ascending("plan_index")]),
        )?;
        let count_schema = TableSchema::new(
            "gravlax.query.batch.counts.v1",
            vec![
                Field::new("query_id", DataType::String),
                Field::new("aggregation", DataType::String),
                Field::new("entity", DataType::String),
                Field::new("umis", DataType::UInt64),
                Field::new("cells", DataType::UInt64).nullable(),
                Field::new("selected_cells", DataType::UInt64).nullable(),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key([
            "query_id",
            "aggregation",
            "entity",
        ]))?;
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "archive_access".into(),
            json!("union of plan-selected archive chunks"),
        );
        parameters.insert("plan_path".into(), json!(plan_path));
        parameters.insert(
            "plan_content_blake3".into(),
            json!(plan_content_blake3
                .as_deref()
                .context("uniform batch plan is missing its bound content digest")?),
        );
        parameters.insert("cell_scope".into(), scope.provenance_json());
        parameters.insert("aggregation".into(), json!(scope.aggregation_name()));
        parameters.insert(
            "selection_policy".into(),
            serde_json::to_value(UniformSelectionPolicy {
                requested_top: top,
                top_zero_means_all: true,
                comparator: matches!(scope.aggregation, QueryAggregation::Cell)
                    .then_some("umis descending, entity ascending (barcode)"),
            })?,
        );
        let context = uniform_query_context(
            archive,
            la.reader().archive_version(),
            la.reader()
                .content_commitment()
                .map(|commitment| commitment.to_hex()),
            parameters,
        );
        let query_selection = SelectionSummary::complete(specs.len() as u64);
        let (available_count_rows, emitted_count_rows) = uniform_counts.iter().try_fold(
            (0u64, 0u64),
            |(available_total, emitted_total), counts| -> Result<_> {
                let (available, emitted) = uniform_row_limit(counts, &scope, top);
                Ok((
                    available_total
                        .checked_add(u64::try_from(available)?)
                        .context("batch available-row count overflow")?,
                    emitted_total
                        .checked_add(u64::try_from(emitted)?)
                        .context("batch emitted-row count overflow")?,
                ))
            },
        )?;
        let count_selection =
            SelectionSummary::selected(available_count_rows, emitted_count_rows)?;
        write_uniform_bundle_output(uniform_output, |writer, format| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                "gravlax.query.batch.result.v1",
                OutputFormat::from(format),
                &context,
                &summary,
            )?;
            bundle.write_table("queries", &query_schema, Some(&query_selection), |rows| {
                for (query_index, query) in specs.iter().enumerate() {
                    let counts = &uniform_counts[query_index];
                    let (available, emitted) = uniform_row_limit(counts, &scope, top);
                    rows.write_row_with(|row| {
                        row.uint64(query_index as u64)?;
                        row.string(&query.id)?;
                        row.string(match query.kind {
                            BatchKind::Region => "region",
                            BatchKind::Junction => "junction",
                        })?;
                        row.string(&query.chrom)?;
                        row.uint64(query.start as u64)?;
                        row.uint64(query.end as u64)?;
                        row.boolean(
                            query.kind == BatchKind::Region || junction_info[query_index].is_some(),
                        )?;
                        if let Some((support, _)) = junction_info[query_index] {
                            row.uint64(support)?;
                        } else {
                            row.null()?;
                        }
                        if let Some((_, postings)) = junction_info[query_index] {
                            row.uint64(postings as u64)?;
                        } else {
                            row.null()?;
                        }
                        if query.kind == BatchKind::Region {
                            row.uint64(if scope.active {
                                scoped_region_molecules[query_index]
                            } else {
                                region_molecules[query_index]
                            })?;
                            row.boolean(true)?;
                        } else {
                            row.null()?;
                            row.null()?;
                        }
                        row.uint64(counts.total_umis as u64)?;
                        row.uint64(counts.cells.len() as u64)?;
                        row.uint64(available as u64)?;
                        row.uint64(emitted as u64)?;
                        row.boolean(emitted < available)?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table("counts", &count_schema, Some(&count_selection), |rows| {
                for (query, counts) in specs.iter().zip(&uniform_counts) {
                    let (_, emitted) = uniform_row_limit(counts, &scope, top);
                    match scope.aggregation {
                        QueryAggregation::Cell => {
                            for (cell, umis) in counts.cells.iter().take(emitted) {
                                let barcode = unpack_cell_bytes(cell_dictionary[*cell as usize]);
                                let barcode = std::str::from_utf8(&barcode)
                                    .expect("packed archive barcode decodes to ASCII");
                                rows.write_row_with(|row| {
                                    row.string(&query.id)?;
                                    row.string("cell")?;
                                    row.string(barcode)?;
                                    row.uint64(*umis as u64)?;
                                    row.null()?;
                                    row.null()?;
                                    Ok(())
                                })?;
                            }
                        }
                        QueryAggregation::Group => {
                            for (group, (umis, cells)) in counts.groups.iter().enumerate() {
                                rows.write_row_with(|row| {
                                    row.string(&query.id)?;
                                    row.string("group")?;
                                    row.string(&scope.group_names[group])?;
                                    row.uint64(*umis as u64)?;
                                    row.uint64(*cells as u64)?;
                                    row.uint64(scope.selected_per_group[group] as u64)?;
                                    Ok(())
                                })?;
                            }
                        }
                        QueryAggregation::Bulk => {
                            rows.write_row_with(|row| {
                                row.string(&query.id)?;
                                row.string("bulk")?;
                                row.string("bulk")?;
                                row.uint64(counts.total_umis as u64)?;
                                row.uint64(counts.cells.len() as u64)?;
                                row.uint64(scope.selected_cells as u64)?;
                                Ok(())
                            })?;
                        }
                    }
                }
                Ok(())
            })?;
            bundle.finish()?;
            Ok(())
        })?;
    } else {
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "schema": if scope.active { "gravlax.query.batch.v2" } else { "gravlax.query.batch.v1" },
                "coordinates": "0-based half-open",
                "query_order": "plan",
                "plan_queries": specs.len(),
                "region_queries": n_regions,
                "junction_queries": n_junctions,
                "planning": {
                    "independent_chunk_decodes": independent_chunk_decodes,
                    "unique_chunk_decodes": unique_chunk_decodes,
                    "chunk_decode_reduction_fraction": reduction,
                },
                "queries": results,
            }))?
        );
    }
    eprintln!(
        "batch {} queries: {} unique / {} independent chunk decodes (open {t_open:.2}s, total {:.2}s)",
        specs.len(),
        unique_chunk_decodes,
        independent_chunk_decodes,
        t0.elapsed().as_secs_f32()
    );
    Ok(())
}

fn discovery_molecule_span(
    molecule: &MolRec,
    shapes: &[evidence_io::archive::Shape],
) -> Option<(u32, u32)> {
    let mut lo = u32::MAX;
    let mut hi = 0u32;
    for (position, shape) in molecule.chains.iter().flat_map(|chain| chain.reps.iter()) {
        lo = lo.min(*position);
        hi = hi.max(
            *position
                + shapes[*shape as usize]
                    .blocks
                    .last()
                    .map(|block| block.0 + block.1)
                    .unwrap_or(0),
        );
    }
    for (position, shape, _, _) in &molecule.mms {
        lo = lo.min(*position);
        hi = hi.max(
            *position
                + shapes[*shape as usize]
                    .blocks
                    .last()
                    .map(|block| block.0 + block.1)
                    .unwrap_or(0),
        );
    }
    (lo != u32::MAX).then_some((lo, hi))
}

#[allow(clippy::too_many_arguments)]
fn discovery_placement_compatible(
    chrom: u32,
    position: u32,
    strand_rev: bool,
    shape: u32,
    nh: u16,
    shapes: &[evidence_io::archive::Shape],
    annotation: &anno::Annotation,
    anno_of: &[Option<u32>],
    solo_strand: anno::assign::SoloStrand,
    placement: &mut evidence_io::Placement,
    txbuf: &mut Vec<u32>,
    genes: &mut Vec<u32>,
) -> bool {
    let Some(ac) = anno_of.get(chrom as usize).copied().flatten() else {
        return false;
    };
    placement_from_parts_into(
        placement,
        chrom,
        position,
        strand_rev,
        &shapes[shape as usize],
        nh,
    );
    anno::assign::concordant_genes_stranded_into(
        placement,
        annotation,
        ac,
        solo_strand,
        txbuf,
        genes,
    );
    !genes.is_empty()
}

#[allow(clippy::too_many_arguments)]
fn discovery_molecule_compatible(
    molecule: &MolRec,
    shapes: &[evidence_io::archive::Shape],
    patterns: &[Vec<PatAlt>],
    annotation: &anno::Annotation,
    anno_of: &[Option<u32>],
    solo_strand: anno::assign::SoloStrand,
    placement: &mut evidence_io::Placement,
    txbuf: &mut Vec<u32>,
    genes: &mut Vec<u32>,
) -> bool {
    for chain in &molecule.chains {
        for (position, shape) in &chain.reps {
            if discovery_placement_compatible(
                molecule.chrom, *position, molecule.strand_rev, *shape, 1, shapes,
                annotation, anno_of, solo_strand, placement, txbuf, genes,
            ) {
                return true;
            }
        }
    }
    for (position, shape, pattern, _) in &molecule.mms {
        for alt in &patterns[*pattern as usize] {
            let alt_shape = if alt.shape == SAME_SHAPE { *shape } else { alt.shape };
            if discovery_placement_compatible(
                alt.chrom, (*position as i64 + alt.offset) as u32,
                molecule.strand_rev != alt.strand_flip, alt_shape, 2, shapes,
                annotation, anno_of, solo_strand, placement, txbuf, genes,
            ) {
                return true;
            }
        }
    }
    false
}

/// Return the unclaimed molecule summaries for one archive chunk. Scratch is worker-local, so
/// archive chunks can be classified in parallel without changing output order.
fn discovery_unclaimed(
    molecules: &[MolRec],
    shapes: &[evidence_io::archive::Shape],
    patterns: Option<&Vec<Vec<PatAlt>>>,
    annotation: &anno::Annotation,
    anno_of: &[Option<u32>],
    mode: DiscoveryClaimMode,
    solo_strand: anno::assign::SoloStrand,
) -> Vec<DiscoveryUnclaimed> {
    let mut out = Vec::new();
    let mut txbuf = Vec::new();
    let mut genes = Vec::new();
    let mut placement = evidence_io::Placement {
        chrom: 0,
        strand: evidence_io::Strand::Forward,
        blocks: Vec::new(),
        junctions: Vec::new(),
        nm: 0,
        score: 0,
        nh: 1,
        clip: (0, 0),
    };
    for molecule in molecules {
        let Some((lo, hi)) = discovery_molecule_span(molecule, shapes) else {
            continue;
        };
        let span_claimed = anno_of
            .get(molecule.chrom as usize)
            .copied()
            .flatten()
            .is_some_and(|ac| {
                annotation.overlapping_into(ac, lo, hi, &mut txbuf);
                !txbuf.is_empty()
            });
        let channel = match mode {
            DiscoveryClaimMode::Span => (!span_claimed).then_some(None),
            DiscoveryClaimMode::StrandSpan => {
                let claimed = anno_of
                    .get(molecule.chrom as usize)
                    .copied()
                    .flatten()
                    .is_some_and(|ac| {
                        annotation.overlapping_into(ac, lo, hi, &mut txbuf);
                        txbuf.iter().any(|&tx| {
                            solo_strand.accepts(
                                molecule.strand_rev,
                                annotation.transcripts[tx as usize].strand_rev,
                            )
                        })
                    });
                (!claimed).then_some(None)
            }
            DiscoveryClaimMode::Compatible => {
                let patterns = patterns.expect("compatible discovery loads patterns");
                (!discovery_molecule_compatible(
                    molecule,
                    shapes,
                    patterns,
                    annotation,
                    anno_of,
                    solo_strand,
                    &mut placement,
                    &mut txbuf,
                    &mut genes,
                ))
                .then_some(None)
            }
            DiscoveryClaimMode::ResidualSites => {
                if !span_claimed {
                    Some(None)
                } else {
                    let patterns = patterns.expect("residual-site discovery loads patterns");
                    (!discovery_molecule_compatible(
                        molecule, shapes, patterns, annotation, anno_of, solo_strand,
                        &mut placement, &mut txbuf, &mut genes,
                    )).then_some(Some(if molecule.strand_rev {
                        lo
                    } else {
                        hi.saturating_sub(1)
                    }))
                }
            }
        };
        if let Some(site) = channel {
            out.push((
                molecule.chrom,
                lo,
                hi,
                molecule.strand_rev,
                molecule.umi_class,
                site,
            ));
        }
    }
    out
}

type DiscoveryUnclaimed = (u32, u32, u32, bool, u32, Option<u32>);

#[derive(Clone, Debug, PartialEq, Eq)]
struct JunctionMeta {
    chrom: u32,
    donor: u32,
    acceptor: u32,
    supporting_children: u64,
    posts: Vec<u32>,
}

fn parse_junction_catalogue(raw: &[u8]) -> Result<Vec<(u32, u32, u32)>> {
    let mut cursor = Cursor::new(raw);
    let mut rows = Vec::new();
    let (mut last_chrom, mut last_donor) = (u32::MAX, 0u32);
    while !cursor.is_empty() {
        let chrom = u32::try_from(cursor.varint()?).context("junction chromosome exceeds u32")?;
        if chrom != last_chrom {
            last_chrom = chrom;
            last_donor = 0;
        }
        let donor_delta =
            u32::try_from(cursor.varint()?).context("junction donor delta exceeds u32")?;
        let donor = last_donor
            .checked_add(donor_delta)
            .context("junction donor delta overflow")?;
        last_donor = donor;
        let intron = u32::try_from(cursor.varint()?).context("junction span exceeds u32")?;
        let acceptor = donor.checked_add(intron).context("junction acceptor overflow")?;
        rows.push((chrom, donor, acceptor));
    }
    Ok(rows)
}

fn parse_junction_postings(raw: &[u8], expected: usize) -> Result<Vec<(u64, Vec<u32>)>> {
    let mut cursor = Cursor::new(raw);
    let mut rows = Vec::with_capacity(expected);
    while !cursor.is_empty() {
        let support = cursor.varint()?;
        let n = usize::try_from(cursor.varint()?).context("junction posting count exceeds usize")?;
        let mut posts = Vec::with_capacity(n);
        let mut last = 0u32;
        for _ in 0..n {
            let delta = u32::try_from(cursor.varint()?).context("junction posting delta exceeds u32")?;
            last = last.checked_add(delta).context("junction posting delta overflow")?;
            posts.push(last);
        }
        rows.push((support, posts));
    }
    if rows.len() != expected {
        bail!(
            "junction catalogue/postings length mismatch: {} coordinates, {} postings",
            expected,
            rows.len()
        );
    }
    Ok(rows)
}

fn read_junction_metadata(la: &mut LazyArchive) -> Result<Vec<JunctionMeta>> {
    let catalogue = parse_junction_catalogue(&la.reader().read("index.junctions")?)?;
    let postings = parse_junction_postings(&la.reader().read("index.jpost")?, catalogue.len())?;
    Ok(catalogue
        .into_iter()
        .zip(postings)
        .map(
            |((chrom, donor, acceptor), (supporting_children, posts))| JunctionMeta {
                chrom,
                donor,
                acceptor,
                supporting_children,
                posts,
            },
        )
        .collect())
}

fn junction_in_window(row: &JunctionMeta, start: u32, end: u32, either: bool) -> bool {
    let donor = row.donor >= start && row.donor < end;
    let acceptor = row.acceptor >= start && row.acceptor < end;
    if either {
        donor || acceptor
    } else {
        donor && acceptor
    }
}

#[derive(Default)]
struct JunctionAnnotation {
    exact: FxHashSet<(u32, u32)>,
    donors: FxHashSet<u32>,
    acceptors: FxHashSet<u32>,
}

fn load_junction_annotation(
    gtf: &std::path::Path,
    chrom: &str,
    bind_provenance: bool,
) -> Result<(JunctionAnnotation, Option<String>)> {
    let (annotation, content_blake3) =
        load_query_annotation(gtf, "junction annotation", bind_provenance)?;
    let Some(&chrom_id) = annotation.chrom_ids.get(chrom) else {
        return Ok((JunctionAnnotation::default(), content_blake3));
    };
    let mut out = JunctionAnnotation::default();
    for transcript in annotation.transcripts.iter().filter(|t| t.chrom == chrom_id) {
        for pair in transcript.exons.windows(2) {
            out.exact.insert((pair[0].end, pair[1].start));
            out.donors.insert(pair[0].end);
            out.acceptors.insert(pair[1].start);
        }
    }
    Ok((out, content_blake3))
}

fn junction_counts_many(
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    rows: &[JunctionMeta],
) -> Result<Vec<FxHashMap<u32, FxHashSet<u32>>>> {
    if rows.is_empty() {
        return Ok(Vec::new());
    }
    let shapes = la.shapes()?;
    let wanted: FxHashMap<(u32, u32), usize> = rows
        .iter()
        .enumerate()
        .map(|(i, row)| ((row.donor, row.acceptor), i))
        .collect();
    let mut posts: Vec<u32> = rows.iter().flat_map(|row| row.posts.iter().copied()).collect();
    posts.sort_unstable();
    posts.dedup();
    if let Some(&bad) = posts.iter().find(|&&post| post as usize >= chunks.len()) {
        bail!("junction posting references missing chunk {bad}");
    }

    let hits: Vec<Vec<(usize, u32)>> = {
        let (reader, tables) = la.reader_and_tables();
        let reader = &*reader;
        posts
            .par_iter()
            .map(|post| {
                let (compressed, raw_len) = reader.read_compressed_at(&format!("c{post}"))?;
                let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                let molecules = decode_chunk(&raw, &chunks[*post as usize], None, tables)?;
                let mut chunk_hits = Vec::new();
                for molecule in &molecules {
                    let mut seen: FxHashSet<usize> = FxHashSet::default();
                    for chain in &molecule.chains {
                        for (position, shape) in &chain.reps {
                            for blocks in shapes[*shape as usize].blocks.windows(2) {
                                let donor = position + blocks[0].0 + blocks[0].1;
                                let acceptor = position + blocks[1].0;
                                if let Some(&row) = wanted.get(&(donor, acceptor)) {
                                    seen.insert(row);
                                }
                            }
                        }
                    }
                    for (position, shape, _, _) in &molecule.mms {
                        for blocks in shapes[*shape as usize].blocks.windows(2) {
                            let donor = position + blocks[0].0 + blocks[0].1;
                            let acceptor = position + blocks[1].0;
                            if let Some(&row) = wanted.get(&(donor, acceptor)) {
                                seen.insert(row);
                            }
                        }
                    }
                    chunk_hits.extend(seen.into_iter().map(|row| (row, molecule.umi_class)));
                }
                Ok(chunk_hits)
            })
            .collect::<Result<_>>()?
    };
    la.prefetch_coc(hits.iter().flatten().map(|(_, class)| *class))?;
    let mut counts: Vec<FxHashMap<u32, FxHashSet<u32>>> =
        (0..rows.len()).map(|_| FxHashMap::default()).collect();
    for (row, class) in hits.into_iter().flatten() {
        counts[row].entry(la.cell_of(class)?).or_default().insert(class);
    }
    Ok(counts)
}

#[derive(Clone)]
struct JunctionSetRequest {
    locus: String,
    chrom: String,
    chrom_id: u32,
    donor: u32,
    acceptor: u32,
    side: u8,
    metadata: Option<JunctionMeta>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct JunctionSetCounts {
    include_only: usize,
    exclude_only: usize,
    both: usize,
}

impl JunctionSetCounts {
    fn add_mask(&mut self, mask: u8) {
        match mask {
            1 => self.include_only += 1,
            2 => self.exclude_only += 1,
            3 => self.both += 1,
            _ => {}
        }
    }

    fn add(&mut self, other: Self) {
        self.include_only += other.include_only;
        self.exclude_only += other.exclude_only;
        self.both += other.both;
    }

    fn informative(&self) -> usize {
        self.include_only + self.exclude_only
    }

    fn total(&self) -> usize {
        self.informative() + self.both
    }

    fn usage(&self) -> Option<f64> {
        (self.informative() > 0).then_some(self.include_only as f64 / self.informative() as f64)
    }

    fn json(&self) -> serde_json::Value {
        json!({
            "include_only": self.include_only,
            "exclude_only": self.exclude_only,
            "both": self.both,
            "informative_umis": self.informative(),
            "usage_fraction": self.usage(),
        })
    }
}

fn prepare_junction_set(
    includes: &[String],
    excludes: &[String],
    chrom_names: &[String],
    metadata: &[JunctionMeta],
) -> Result<Vec<JunctionSetRequest>> {
    let metadata_of: FxHashMap<(u32, u32, u32), &JunctionMeta> = metadata
        .iter()
        .map(|row| ((row.chrom, row.donor, row.acceptor), row))
        .collect();
    let mut seen = FxHashMap::default();
    let mut requests = Vec::new();
    for (side_name, side, loci) in [("include", 1u8, includes), ("exclude", 2u8, excludes)] {
        for locus in loci {
            let (chrom, donor, acceptor) = parse_locus(locus)
                .with_context(|| format!("invalid --{side_name} junction {locus}"))?;
            let chrom_id = chrom_names
                .iter()
                .position(|name| name == &chrom)
                .map(|index| index as u32)
                .with_context(|| format!("unknown chromosome {chrom} in --{side_name}"))?;
            let key = (chrom_id, donor, acceptor);
            if let Some(previous) = seen.insert(key, side) {
                if previous == side {
                    bail!("duplicate --{side_name} junction {locus}");
                }
                bail!("junction {locus} appears on both inclusion and exclusion sides");
            }
            requests.push(JunctionSetRequest {
                locus: locus.clone(),
                chrom,
                chrom_id,
                donor,
                acceptor,
                side,
                metadata: metadata_of.get(&key).map(|row| (*row).clone()),
            });
        }
    }
    Ok(requests)
}

fn junction_set_class_masks(
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    requests: &[JunctionSetRequest],
) -> Result<(FxHashMap<u32, u8>, usize, usize)> {
    let shapes = la.shapes()?;
    let mut chunk_wanted: Vec<FxHashMap<(u32, u32), u8>> =
        (0..chunks.len()).map(|_| FxHashMap::default()).collect();
    let mut independent_chunk_decodes = 0usize;
    for request in requests {
        let Some(metadata) = &request.metadata else {
            continue;
        };
        for &post in &metadata.posts {
            let chunk = chunks.get(post as usize).with_context(|| {
                format!("junction {} references missing chunk {post}", request.locus)
            })?;
            if chunk.chrom != request.chrom_id {
                bail!(
                    "junction {} posting references a different chromosome",
                    request.locus
                );
            }
            *chunk_wanted[post as usize]
                .entry((request.donor, request.acceptor))
                .or_insert(0) |= request.side;
            independent_chunk_decodes += 1;
        }
    }
    let selected: Vec<usize> = chunk_wanted
        .iter()
        .enumerate()
        .filter_map(|(index, wanted)| (!wanted.is_empty()).then_some(index))
        .collect();

    let hits: Vec<Vec<(u32, u8)>> = {
        let (reader, tables) = la.reader_and_tables();
        let reader = &*reader;
        selected
            .par_iter()
            .map(|&chunk_index| -> Result<Vec<(u32, u8)>> {
                let (compressed, raw_len) =
                    reader.read_compressed_at(&format!("c{chunk_index}"))?;
                let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                let molecules = decode_chunk(&raw, &chunks[chunk_index], None, tables)?;
                let wanted = &chunk_wanted[chunk_index];
                let mut chunk_hits = Vec::new();
                for molecule in &molecules {
                    let mut mask = 0u8;
                    let mut inspect = |position: u32, shape: u32| {
                        for blocks in shapes[shape as usize].blocks.windows(2) {
                            let donor = position + blocks[0].0 + blocks[0].1;
                            let acceptor = position + blocks[1].0;
                            if let Some(side) = wanted.get(&(donor, acceptor)) {
                                mask |= *side;
                            }
                        }
                    };
                    for chain in &molecule.chains {
                        for &(position, shape) in &chain.reps {
                            inspect(position, shape);
                        }
                    }
                    for &(position, shape, _, _) in &molecule.mms {
                        inspect(position, shape);
                    }
                    if mask != 0 {
                        chunk_hits.push((molecule.umi_class, mask));
                    }
                }
                chunk_hits.sort_unstable();
                let mut reduced: Vec<(u32, u8)> = Vec::with_capacity(chunk_hits.len());
                for (class, mask) in chunk_hits {
                    match reduced.last_mut() {
                        Some((previous, combined)) if *previous == class => *combined |= mask,
                        _ => reduced.push((class, mask)),
                    }
                }
                Ok(reduced)
            })
            .collect::<Result<_>>()?
    };
    la.prefetch_coc(hits.iter().flatten().map(|(class, _)| *class))?;
    let mut class_masks = FxHashMap::default();
    for (class, mask) in hits.into_iter().flatten() {
        *class_masks.entry(class).or_insert(0) |= mask;
    }
    Ok((class_masks, selected.len(), independent_chunk_decodes))
}

#[allow(clippy::too_many_arguments)]
fn run_junction_set(
    archive: &Path,
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    chrom_names: &[String],
    includes: &[String],
    excludes: &[String],
    top: usize,
    tsv: bool,
    json_output: bool,
    uniform_output: &UniformQueryOutputArgs,
    scope_args: &QueryScopeArgs,
    t0: std::time::Instant,
    t_open: f32,
) -> Result<()> {
    let mut scope = load_query_scope(la, scope_args)?;
    if uniform_output.format.is_some() {
        scope.ensure_resolved_mapping_digest()?;
    }
    let metadata = read_junction_metadata(la)?;
    let requests = prepare_junction_set(includes, excludes, chrom_names, &metadata)?;
    let (class_masks, unique_chunk_decodes, independent_chunk_decodes) =
        junction_set_class_masks(la, chunks, &requests)?;
    let mut per_cell: FxHashMap<u32, JunctionSetCounts> = FxHashMap::default();
    for (class, mask) in class_masks {
        let cell = la.cell_of(class)?;
        if scope.includes(cell) {
            per_cell.entry(cell).or_default().add_mask(mask);
        }
    }
    let mut cells: Vec<(u32, JunctionSetCounts)> = per_cell.into_iter().collect();
    let dictionary = la.cells()?.to_vec();
    if uniform_output.format.is_some() {
        cells.sort_unstable_by_key(|(cell, counts)| {
            (
                std::cmp::Reverse(counts.total()),
                dictionary[*cell as usize],
            )
        });
    } else {
        cells.sort_unstable_by_key(|(cell, counts)| (std::cmp::Reverse(counts.total()), *cell));
    }
    let mut totals = JunctionSetCounts::default();
    for (_, counts) in &cells {
        totals.add(*counts);
    }
    let mut groups = vec![(JunctionSetCounts::default(), 0usize); scope.group_names.len()];
    for (cell, counts) in &cells {
        if let Some(&group) = scope.group_of.get(cell) {
            groups[group as usize].0.add(*counts);
            groups[group as usize].1 += 1;
        }
    }
    let limit = if top == 0 {
        cells.len()
    } else {
        top.min(cells.len())
    };
    let request_json = |request: &JunctionSetRequest| {
        let mut value = json!({
            "locus": request.locus,
            "chrom": request.chrom,
            "donor": request.donor,
            "acceptor": request.acceptor,
            "present": request.metadata.is_some(),
        });
        if let Some(metadata) = &request.metadata {
            let object = value.as_object_mut().unwrap();
            object.insert(
                "supporting_children".into(),
                json!(metadata.supporting_children),
            );
            object.insert("posting_chunks".into(), json!(metadata.posts.len()));
        }
        value
    };

    if uniform_output.format.is_some() {
        let summary = json!({
            "coordinates": "0-based half-open junction boundaries",
            "scope": scope.json(),
            "semantics": {
                "class_categories": ["include_only", "exclude_only", "both"],
                "informative_umis": "include_only + exclude_only",
                "usage_fraction": "include_only / informative_umis",
                "both_in_usage_denominator": false,
                "zero_denominator": null,
            },
            "totals": totals.json(),
            "supporting_cells": cells.len(),
            "planning": {
                "independent_chunk_decodes": independent_chunk_decodes,
                "unique_chunk_decodes": unique_chunk_decodes,
                "chunk_decode_reduction_fraction": if independent_chunk_decodes == 0 {
                    0.0
                } else {
                    1.0 - unique_chunk_decodes as f64 / independent_chunk_decodes as f64
                },
            },
        });
        let request_schema = TableSchema::new(
            "gravlax.query.jset.junctions.v1",
            vec![
                Field::new("request_index", DataType::UInt64),
                Field::new("side", DataType::String),
                Field::new("side_index", DataType::UInt64),
                Field::new("locus", DataType::String),
                Field::new("chrom", DataType::String),
                Field::new("donor", DataType::UInt64),
                Field::new("acceptor", DataType::UInt64),
                Field::new("present", DataType::Boolean),
                Field::new("archive_supporting_children", DataType::UInt64).nullable(),
                Field::new("archive_posting_chunks", DataType::UInt64).nullable(),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["request_index"])
                .ordered_by([OrderKey::ascending("request_index")]),
        )?;
        let count_schema = TableSchema::new(
            "gravlax.query.jset.counts.v1",
            vec![
                Field::new("aggregation", DataType::String),
                Field::new("entity", DataType::String),
                Field::new("include_only", DataType::UInt64),
                Field::new("exclude_only", DataType::UInt64),
                Field::new("both", DataType::UInt64),
                Field::new("informative_umis", DataType::UInt64),
                Field::new("usage_fraction", DataType::Float64).nullable(),
                Field::new("cells", DataType::UInt64).nullable(),
                Field::new("selected_cells", DataType::UInt64).nullable(),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Set).with_key(["aggregation", "entity"]),
        )?;
        let (available, emitted) = match scope.aggregation {
            QueryAggregation::Cell => (cells.len(), limit),
            QueryAggregation::Group => (groups.len(), groups.len()),
            QueryAggregation::Bulk => (1, 1),
        };
        let request_selection = SelectionSummary::complete(requests.len() as u64);
        let count_selection = SelectionSummary::selected(available as u64, emitted as u64)?;
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "archive_access".into(),
            json!("junction catalogue and union of requested postings"),
        );
        parameters.insert("inclusion_junctions".into(), json!(includes));
        parameters.insert("exclusion_junctions".into(), json!(excludes));
        parameters.insert("cell_scope".into(), scope.provenance_json());
        parameters.insert("aggregation".into(), json!(scope.aggregation_name()));
        parameters.insert(
            "selection_policy".into(),
            serde_json::to_value(UniformSelectionPolicy {
                requested_top: top,
                top_zero_means_all: true,
                comparator: matches!(scope.aggregation, QueryAggregation::Cell)
                    .then_some("total categorized UMIs descending, entity ascending (barcode)"),
            })?,
        );
        let context = uniform_query_context(
            archive,
            la.reader().archive_version(),
            la.reader()
                .content_commitment()
                .map(|commitment| commitment.to_hex()),
            parameters,
        );
        write_uniform_bundle_output(uniform_output, |writer, format| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                "gravlax.query.jset.result.v1",
                OutputFormat::from(format),
                &context,
                &summary,
            )?;
            bundle.write_table(
                "junctions",
                &request_schema,
                Some(&request_selection),
                |rows| {
                    let mut inclusion_index = 0_u64;
                    let mut exclusion_index = 0_u64;
                    for (request_index, request) in requests.iter().enumerate() {
                        let (side, side_index) = if request.side == 1 {
                            let index = inclusion_index;
                            inclusion_index += 1;
                            ("include", index)
                        } else {
                            let index = exclusion_index;
                            exclusion_index += 1;
                            ("exclude", index)
                        };
                        rows.write_row_with(|row| {
                            row.uint64(request_index as u64)?;
                            row.string(side)?;
                            row.uint64(side_index)?;
                            row.string(&request.locus)?;
                            row.string(&request.chrom)?;
                            row.uint64(request.donor as u64)?;
                            row.uint64(request.acceptor as u64)?;
                            row.boolean(request.metadata.is_some())?;
                            if let Some(metadata) = &request.metadata {
                                row.uint64(metadata.supporting_children)?;
                                row.uint64(metadata.posts.len() as u64)?;
                            } else {
                                row.null()?;
                                row.null()?;
                            }
                            Ok(())
                        })?;
                    }
                    Ok(())
                },
            )?;
            bundle.write_table("counts", &count_schema, Some(&count_selection), |rows| {
                let mut emit = |aggregation: &str,
                                entity: &str,
                                counts: &JunctionSetCounts,
                                support_cells: Option<usize>,
                                selected_cells: Option<usize>|
                 -> std::result::Result<(), OutputError> {
                    rows.write_row_with(|row| {
                        row.string(aggregation)?;
                        row.string(entity)?;
                        row.uint64(counts.include_only as u64)?;
                        row.uint64(counts.exclude_only as u64)?;
                        row.uint64(counts.both as u64)?;
                        row.uint64(counts.informative() as u64)?;
                        if let Some(usage) = counts.usage() {
                            row.float64(usage)?;
                        } else {
                            row.null()?;
                        }
                        if let Some(cells) = support_cells {
                            row.uint64(cells as u64)?;
                        } else {
                            row.null()?;
                        }
                        if let Some(cells) = selected_cells {
                            row.uint64(cells as u64)?;
                        } else {
                            row.null()?;
                        }
                        Ok(())
                    })
                };
                match scope.aggregation {
                    QueryAggregation::Cell => {
                        for (cell, counts) in cells.iter().take(limit) {
                            let barcode = unpack_cell_bytes(dictionary[*cell as usize]);
                            let barcode = std::str::from_utf8(&barcode)
                                .expect("packed archive barcode decodes to ASCII");
                            emit("cell", barcode, counts, None, None)?;
                        }
                    }
                    QueryAggregation::Group => {
                        for (group, (counts, support_cells)) in groups.iter().enumerate() {
                            emit(
                                "group",
                                &scope.group_names[group],
                                counts,
                                Some(*support_cells),
                                Some(scope.selected_per_group[group]),
                            )?;
                        }
                    }
                    QueryAggregation::Bulk => emit(
                        "bulk",
                        "bulk",
                        &totals,
                        Some(cells.len()),
                        Some(scope.selected_cells),
                    )?,
                }
                Ok(())
            })?;
            bundle.finish()?;
            Ok(())
        })?;
    } else if tsv {
        let usage = |counts: &JunctionSetCounts| {
            counts
                .usage()
                .map_or_else(|| "NA".to_owned(), |value| format!("{value:.9}"))
        };
        match scope.aggregation {
            QueryAggregation::Cell => {
                println!(
                    "barcode\tinclude_only\texclude_only\tboth\tinformative_umis\tusage_fraction"
                );
                for (cell, counts) in cells.iter().take(limit) {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        unpack_cell(&dictionary, *cell),
                        counts.include_only,
                        counts.exclude_only,
                        counts.both,
                        counts.informative(),
                        usage(counts)
                    );
                }
            }
            QueryAggregation::Group => {
                println!("group\tinclude_only\texclude_only\tboth\tinformative_umis\tusage_fraction\tcells\tselected_cells");
                for (group, (counts, support_cells)) in groups.iter().enumerate() {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        scope.group_names[group],
                        counts.include_only,
                        counts.exclude_only,
                        counts.both,
                        counts.informative(),
                        usage(counts),
                        support_cells,
                        scope.selected_per_group[group]
                    );
                }
            }
            QueryAggregation::Bulk => {
                println!("scope\tinclude_only\texclude_only\tboth\tinformative_umis\tusage_fraction\tcells\tselected_cells");
                println!(
                    "bulk\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    totals.include_only,
                    totals.exclude_only,
                    totals.both,
                    totals.informative(),
                    usage(&totals),
                    cells.len(),
                    scope.selected_cells
                );
            }
        }
    } else if json_output {
        let mut value = json!({
            "schema": "gravlax.query.jset.v1",
            "coordinates": "0-based half-open junction boundaries",
            "semantics": {
                "class_categories": ["include_only", "exclude_only", "both"],
                "informative_umis": "include_only + exclude_only",
                "usage_fraction": "include_only / informative_umis",
                "both_in_usage_denominator": false,
                "zero_denominator": null,
            },
            "scope": scope.json(),
            "inclusion_junctions": requests.iter().filter(|request| request.side == 1).map(request_json).collect::<Vec<_>>(),
            "exclusion_junctions": requests.iter().filter(|request| request.side == 2).map(request_json).collect::<Vec<_>>(),
            "planning": {
                "independent_chunk_decodes": independent_chunk_decodes,
                "unique_chunk_decodes": unique_chunk_decodes,
                "chunk_decode_reduction_fraction": if independent_chunk_decodes == 0 {
                    0.0
                } else {
                    1.0 - unique_chunk_decodes as f64 / independent_chunk_decodes as f64
                },
            },
            "totals": totals.json(),
            "cells": cells.len(),
        });
        let object = value.as_object_mut().unwrap();
        match scope.aggregation {
            QueryAggregation::Cell => {
                object.insert(
                    "cell_rows".into(),
                    json!(cells
                        .iter()
                        .take(limit)
                        .map(|(cell, counts)| {
                            let mut row = counts.json();
                            row.as_object_mut()
                                .unwrap()
                                .insert("barcode".into(), json!(unpack_cell(&dictionary, *cell)));
                            row
                        })
                        .collect::<Vec<_>>()),
                );
                object.insert("cell_rows_truncated".into(), json!(limit < cells.len()));
            }
            QueryAggregation::Group => {
                object.insert(
                    "group_rows".into(),
                    json!(groups
                        .iter()
                        .enumerate()
                        .map(|(group, (counts, support_cells))| {
                            let mut row = counts.json();
                            let object = row.as_object_mut().unwrap();
                            object.insert("group".into(), json!(scope.group_names[group]));
                            object.insert("cells".into(), json!(support_cells));
                            object.insert(
                                "selected_cells".into(),
                                json!(scope.selected_per_group[group]),
                            );
                            row
                        })
                        .collect::<Vec<_>>()),
                );
            }
            QueryAggregation::Bulk => {
                object.insert("bulk".into(), totals.json());
            }
        }
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        println!(
            "jset: {} include-only / {} exclude-only / {} both classes across {} cells; usage={} ({} unique / {} independent chunks; open {t_open:.2}s, total {:.2}s)",
            totals.include_only,
            totals.exclude_only,
            totals.both,
            cells.len(),
            totals.usage().map_or_else(|| "NA".to_owned(), |value| format!("{value:.4}")),
            unique_chunk_decodes,
            independent_chunk_decodes,
            t0.elapsed().as_secs_f32()
        );
        match scope.aggregation {
            QueryAggregation::Cell => {
                for (cell, counts) in cells.iter().take(limit) {
                    println!(
                        "  {}\t{}\t{}\t{}\t{}",
                        unpack_cell(&dictionary, *cell),
                        counts.include_only,
                        counts.exclude_only,
                        counts.both,
                        counts
                            .usage()
                            .map_or_else(|| "NA".to_owned(), |value| format!("{value:.4}"))
                    );
                }
            }
            QueryAggregation::Group => {
                for (group, (counts, _)) in groups.iter().enumerate() {
                    println!(
                        "  {}\t{}\t{}\t{}\t{}",
                        scope.group_names[group],
                        counts.include_only,
                        counts.exclude_only,
                        counts.both,
                        counts
                            .usage()
                            .map_or_else(|| "NA".to_owned(), |value| format!("{value:.4}"))
                    );
                }
            }
            QueryAggregation::Bulk => {}
        }
    }
    eprintln!(
        "jset: {} unique / {} independent chunk decodes (open {t_open:.2}s, total {:.2}s)",
        unique_chunk_decodes,
        independent_chunk_decodes,
        t0.elapsed().as_secs_f32()
    );
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct EventKey {
    kind: EventTypeArg,
    chrom: String,
    includes: Vec<(u32, u32)>,
    excludes: Vec<(u32, u32)>,
}

impl EventKey {
    fn id(&self) -> String {
        let coordinates = self
            .includes
            .iter()
            .map(|(donor, acceptor)| format!("{donor}-{acceptor}"))
            .chain(
                self.excludes
                    .iter()
                    .map(|(donor, acceptor)| format!("{donor}-{acceptor}")),
            )
            .collect::<Vec<_>>()
            .join(":");
        format!("{}:{}:{coordinates}", self.kind.name(), self.chrom)
    }

    fn components(&self) -> impl Iterator<Item = ((u32, u32), u8)> + '_ {
        self.includes
            .iter()
            .copied()
            .map(|junction| (junction, 1u8))
            .chain(
                self.excludes
                    .iter()
                    .copied()
                    .map(|junction| (junction, 2u8)),
            )
    }
}

impl EventTypeArg {
    fn name(self) -> &'static str {
        match self {
            Self::AltAcceptor => "alt_acceptor",
            Self::AltDonor => "alt_donor",
            Self::Cassette => "cassette",
        }
    }
}

fn selected_event_types(requested: &[EventTypeArg]) -> BTreeSet<EventTypeArg> {
    if requested.is_empty() {
        [
            EventTypeArg::AltAcceptor,
            EventTypeArg::AltDonor,
            EventTypeArg::Cassette,
        ]
        .into_iter()
        .collect()
    } else {
        requested.iter().copied().collect()
    }
}

#[expect(
    clippy::too_many_arguments,
    reason = "event discovery keeps genomic bounds, catalogue input, and hard scientific thresholds explicit"
)]
fn discover_event_keys(
    chrom: &str,
    chrom_id: u32,
    start: u32,
    end: u32,
    metadata: &[JunctionMeta],
    requested: &[EventTypeArg],
    min_support: u64,
    max_events: usize,
) -> Result<Vec<EventKey>> {
    if min_support == 0 {
        bail!("--min-support must be at least 1");
    }
    if max_events == 0 {
        bail!("--max-events must be at least 1");
    }
    if max_events > MAX_PACKED_EVENTS {
        bail!("--max-events must not exceed {MAX_PACKED_EVENTS}");
    }
    let types = selected_event_types(requested);
    let mut rows: Vec<&JunctionMeta> = metadata
        .iter()
        .filter(|row| {
            row.chrom == chrom_id
                && row.supporting_children >= min_support
                && row.donor >= start
                && row.acceptor < end
        })
        .collect();
    rows.sort_unstable_by_key(|row| (row.donor, row.acceptor));
    let coordinate_set: FxHashSet<(u32, u32)> =
        rows.iter().map(|row| (row.donor, row.acceptor)).collect();
    let mut by_donor: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    let mut by_acceptor: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for row in &rows {
        by_donor.entry(row.donor).or_default().push(row.acceptor);
        by_acceptor.entry(row.acceptor).or_default().push(row.donor);
    }
    for acceptors in by_donor.values_mut() {
        acceptors.sort_unstable();
        acceptors.dedup();
    }
    for donors in by_acceptor.values_mut() {
        donors.sort_unstable();
        donors.dedup();
    }

    let mut events = BTreeSet::new();
    let mut insert = |event: EventKey| -> Result<()> {
        events.insert(event);
        if events.len() > max_events {
            bail!(
                "event catalogue exceeds --max-events {max_events}; narrow the locus, increase --min-support, or raise the explicit limit"
            );
        }
        Ok(())
    };

    if types.contains(&EventTypeArg::AltAcceptor) {
        for (&donor, acceptors) in &by_donor {
            for lower in 0..acceptors.len() {
                for higher in lower + 1..acceptors.len() {
                    insert(EventKey {
                        kind: EventTypeArg::AltAcceptor,
                        chrom: chrom.to_owned(),
                        includes: vec![(donor, acceptors[lower])],
                        excludes: vec![(donor, acceptors[higher])],
                    })?;
                }
            }
        }
    }
    if types.contains(&EventTypeArg::AltDonor) {
        for (&acceptor, donors) in &by_acceptor {
            for lower in 0..donors.len() {
                for higher in lower + 1..donors.len() {
                    insert(EventKey {
                        kind: EventTypeArg::AltDonor,
                        chrom: chrom.to_owned(),
                        includes: vec![(donors[lower], acceptor)],
                        excludes: vec![(donors[higher], acceptor)],
                    })?;
                }
            }
        }
    }
    if types.contains(&EventTypeArg::Cassette) {
        for skip in &rows {
            let Some(left_acceptors) = by_donor.get(&skip.donor) else {
                continue;
            };
            let Some(right_donors) = by_acceptor.get(&skip.acceptor) else {
                continue;
            };
            for &inner_acceptor in left_acceptors {
                if inner_acceptor >= skip.acceptor {
                    break;
                }
                for &inner_donor in right_donors {
                    if inner_donor < inner_acceptor || inner_donor <= skip.donor {
                        continue;
                    }
                    if inner_donor >= skip.acceptor {
                        break;
                    }
                    if coordinate_set.contains(&(skip.donor, inner_acceptor))
                        && coordinate_set.contains(&(inner_donor, skip.acceptor))
                    {
                        insert(EventKey {
                            kind: EventTypeArg::Cassette,
                            chrom: chrom.to_owned(),
                            includes: vec![
                                (skip.donor, inner_acceptor),
                                (inner_donor, skip.acceptor),
                            ],
                            excludes: vec![(skip.donor, skip.acceptor)],
                        })?;
                    }
                }
            }
        }
    }
    Ok(events.into_iter().collect())
}

#[derive(Clone)]
struct EventComponent {
    donor: u32,
    acceptor: u32,
    side: u8,
    metadata: Option<usize>,
}

#[derive(Clone)]
struct EventDefinition {
    key: EventKey,
    chrom_id: u32,
    components: Vec<EventComponent>,
    catalogue_present: bool,
}

type EventComponentTargets = FxHashMap<(u32, u32), Vec<(u32, u8)>>;

fn prepare_event_definitions(
    keys: &[EventKey],
    chrom_names: &[String],
    metadata: &[JunctionMeta],
) -> Result<Vec<EventDefinition>> {
    let metadata_of: FxHashMap<(u32, u32, u32), usize> = metadata
        .iter()
        .enumerate()
        .map(|(index, row)| ((row.chrom, row.donor, row.acceptor), index))
        .collect();
    keys.iter()
        .map(|key| {
            let chrom_id = chrom_names
                .iter()
                .position(|name| name == &key.chrom)
                .map(|index| index as u32)
                .with_context(|| format!("unknown chromosome {}", key.chrom))?;
            let components: Vec<EventComponent> = key
                .components()
                .map(|((donor, acceptor), side)| EventComponent {
                    donor,
                    acceptor,
                    side,
                    metadata: metadata_of.get(&(chrom_id, donor, acceptor)).copied(),
                })
                .collect();
            Ok(EventDefinition {
                key: key.clone(),
                chrom_id,
                catalogue_present: components
                    .iter()
                    .all(|component| component.metadata.is_some()),
                components,
            })
        })
        .collect()
}

const EVENT_HIT_MASK_BITS: u32 = 2;
const EVENT_HIT_CLASS_BITS: u32 = 32;
const EVENT_HIT_EVENT_SHIFT: u32 = EVENT_HIT_CLASS_BITS + EVENT_HIT_MASK_BITS;
const MAX_PACKED_EVENTS: usize = 1usize << (64 - EVENT_HIT_EVENT_SHIFT);

fn pack_event_hit(event: u32, class: u32, mask: u8) -> u64 {
    debug_assert!(mask & !3 == 0);
    debug_assert!((event as usize) < MAX_PACKED_EVENTS);
    ((event as u64) << EVENT_HIT_EVENT_SHIFT)
        | ((class as u64) << EVENT_HIT_MASK_BITS)
        | u64::from(mask)
}

fn unpack_event_hit(hit: u64) -> (u32, u32, u8) {
    (
        (hit >> EVENT_HIT_EVENT_SHIFT) as u32,
        ((hit >> EVENT_HIT_MASK_BITS) & u64::from(u32::MAX)) as u32,
        (hit & 3) as u8,
    )
}

fn reduce_sorted_packed_hits(hits: Vec<u64>) -> Vec<u64> {
    let mut reduced = Vec::with_capacity(hits.len());
    for hit in hits {
        let key = hit & !3;
        match reduced.last_mut() {
            Some(previous) if *previous & !3 == key => *previous |= hit & 3,
            _ => reduced.push(hit),
        }
    }
    reduced
}

fn event_packed_hits(
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    metadata: &[JunctionMeta],
    events: &[EventDefinition],
) -> Result<(Vec<u64>, usize, usize)> {
    let shapes = la.shapes()?;
    let mut chunk_wanted: Vec<EventComponentTargets> =
        (0..chunks.len()).map(|_| FxHashMap::default()).collect();
    let mut independent_chunk_decodes = 0usize;
    for (event_index, event) in events.iter().enumerate() {
        for component in &event.components {
            let Some(metadata_index) = component.metadata else {
                continue;
            };
            let row = &metadata[metadata_index];
            for &post in &row.posts {
                let chunk = chunks.get(post as usize).with_context(|| {
                    format!("event {} references missing chunk {post}", event.key.id())
                })?;
                if chunk.chrom != event.chrom_id {
                    bail!(
                        "event {} posting references a different chromosome",
                        event.key.id()
                    );
                }
                chunk_wanted[post as usize]
                    .entry((component.donor, component.acceptor))
                    .or_default()
                    .push((event_index as u32, component.side));
                independent_chunk_decodes += 1;
            }
        }
    }
    for wanted in &mut chunk_wanted {
        for targets in wanted.values_mut() {
            targets.sort_unstable();
            targets.dedup();
        }
    }
    let selected: Vec<usize> = chunk_wanted
        .iter()
        .enumerate()
        .filter_map(|(index, wanted)| (!wanted.is_empty()).then_some(index))
        .collect();

    let chunk_hits: Vec<Vec<u64>> = {
        let (reader, tables) = la.reader_and_tables();
        let reader = &*reader;
        selected
            .par_iter()
            .map(|&chunk_index| -> Result<Vec<u64>> {
                let (compressed, raw_len) =
                    reader.read_compressed_at(&format!("c{chunk_index}"))?;
                let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                let molecules = decode_chunk(&raw, &chunks[chunk_index], None, tables)?;
                let wanted = &chunk_wanted[chunk_index];
                let mut hits = Vec::new();
                for molecule in &molecules {
                    let mut event_sides = Vec::new();
                    let mut inspect = |position: u32, shape: u32| {
                        for blocks in shapes[shape as usize].blocks.windows(2) {
                            let donor = position + blocks[0].0 + blocks[0].1;
                            let acceptor = position + blocks[1].0;
                            if let Some(targets) = wanted.get(&(donor, acceptor)) {
                                event_sides.extend_from_slice(targets);
                            }
                        }
                    };
                    for chain in &molecule.chains {
                        for &(position, shape) in &chain.reps {
                            inspect(position, shape);
                        }
                    }
                    for &(position, shape, _, _) in &molecule.mms {
                        inspect(position, shape);
                    }
                    event_sides.sort_unstable();
                    let mut reduced: Vec<(u32, u8)> = Vec::with_capacity(event_sides.len());
                    for (event, side) in event_sides {
                        match reduced.last_mut() {
                            Some((previous, mask)) if *previous == event => *mask |= side,
                            _ => reduced.push((event, side)),
                        }
                    }
                    hits.extend(
                        reduced
                            .into_iter()
                            .map(|(event, mask)| pack_event_hit(event, molecule.umi_class, mask)),
                    );
                }
                hits.sort_unstable();
                Ok(reduce_sorted_packed_hits(hits))
            })
            .collect::<Result<_>>()?
    };
    la.prefetch_coc(
        chunk_hits
            .iter()
            .flatten()
            .map(|hit| unpack_event_hit(*hit).1),
    )?;
    let mut all_hits: Vec<u64> = chunk_hits.into_iter().flatten().collect();
    all_hits.sort_unstable();
    Ok((
        reduce_sorted_packed_hits(all_hits),
        selected.len(),
        independent_chunk_decodes,
    ))
}

fn packed_hits_by_event(hits: Vec<u64>, event_count: usize) -> Vec<Vec<(u32, u8)>> {
    let mut by_event = vec![Vec::new(); event_count];
    for hit in hits {
        let (event, class, mask) = unpack_event_hit(hit);
        by_event[event as usize].push((class, mask));
    }
    by_event
}

#[derive(Debug, PartialEq, Eq)]
struct EventResult {
    totals: JunctionSetCounts,
    cells: Vec<(u32, JunctionSetCounts)>,
    support_cells: usize,
    groups: Vec<(JunctionSetCounts, usize)>,
}

#[cfg(test)]
fn reduce_event_results_with<F>(
    class_masks: Vec<Vec<(u32, u8)>>,
    scope: &QueryScope,
    mut cell_of: F,
) -> Result<Vec<EventResult>>
where
    F: FnMut(u32) -> Result<u32>,
{
    reduce_event_results_with_order(class_masks, scope, None, &mut cell_of)
}

fn reduce_event_results_with_order<F>(
    class_masks: Vec<Vec<(u32, u8)>>,
    scope: &QueryScope,
    packed_barcodes: Option<&[u32]>,
    mut cell_of: F,
) -> Result<Vec<EventResult>>
where
    F: FnMut(u32) -> Result<u32>,
{
    class_masks
        .into_iter()
        .map(|hits| {
            let mut per_cell: FxHashMap<u32, JunctionSetCounts> = FxHashMap::default();
            for (class, mask) in hits {
                let cell = cell_of(class)?;
                if scope.includes(cell) {
                    per_cell.entry(cell).or_default().add_mask(mask);
                }
            }
            let mut cells: Vec<(u32, JunctionSetCounts)> = per_cell.into_iter().collect();
            match packed_barcodes {
                Some(barcodes) => cells.sort_unstable_by_key(|(cell, counts)| {
                    (std::cmp::Reverse(counts.total()), barcodes[*cell as usize])
                }),
                None => cells.sort_unstable_by_key(|(cell, counts)| {
                    (std::cmp::Reverse(counts.total()), *cell)
                }),
            }
            let mut totals = JunctionSetCounts::default();
            let mut groups = vec![(JunctionSetCounts::default(), 0usize); scope.group_names.len()];
            for (cell, counts) in &cells {
                totals.add(*counts);
                if let Some(&group) = scope.group_of.get(cell) {
                    groups[group as usize].0.add(*counts);
                    groups[group as usize].1 += 1;
                }
            }
            Ok(EventResult {
                totals,
                support_cells: cells.len(),
                cells,
                groups,
            })
        })
        .collect()
}

fn reduce_event_results(
    la: &mut LazyArchive,
    class_masks: Vec<Vec<(u32, u8)>>,
    scope: &QueryScope,
    packed_barcodes: Option<&[u32]>,
) -> Result<Vec<EventResult>> {
    reduce_event_results_with_order(class_masks, scope, packed_barcodes, |class| {
        la.cell_of(class)
    })
}

fn reduce_packed_event_results_with<F>(
    hits: &[u64],
    event_count: usize,
    scope: &QueryScope,
    mut cell_of: F,
) -> Result<Vec<EventResult>>
where
    F: FnMut(u32) -> Result<u32>,
{
    let mut results: Vec<EventResult> = (0..event_count)
        .map(|_| EventResult {
            totals: JunctionSetCounts::default(),
            cells: Vec::new(),
            support_cells: 0,
            groups: vec![(JunctionSetCounts::default(), 0); scope.group_names.len()],
        })
        .collect();
    let mut active_event = None;
    let mut support_cells = FxHashSet::default();
    for &hit in hits {
        let (event, class, mask) = unpack_event_hit(hit);
        let event = event as usize;
        if event >= event_count {
            bail!("packed event hit references out-of-range event {event}");
        }
        if active_event != Some(event) {
            support_cells.clear();
            active_event = Some(event);
        }
        let cell = cell_of(class)?;
        if !scope.includes(cell) {
            continue;
        }
        let result = &mut results[event];
        result.totals.add_mask(mask);
        if let Some(&group) = scope.group_of.get(&cell) {
            result.groups[group as usize].0.add_mask(mask);
        }
        if support_cells.insert(cell) {
            result.support_cells += 1;
            if let Some(&group) = scope.group_of.get(&cell) {
                result.groups[group as usize].1 += 1;
            }
        }
    }
    Ok(results)
}

fn reduce_packed_event_results(
    la: &mut LazyArchive,
    hits: &[u64],
    event_count: usize,
    scope: &QueryScope,
) -> Result<Vec<EventResult>> {
    reduce_packed_event_results_with(hits, event_count, scope, |class| la.cell_of(class))
}

struct EventAnnotationIndex<'a> {
    annotation: &'a anno::Annotation,
    junctions: FxHashMap<(u32, u32), Vec<(u32, bool)>>,
}

fn build_event_annotation_index<'a>(
    annotation: &'a anno::Annotation,
    events: &[EventDefinition],
) -> EventAnnotationIndex<'a> {
    let Some(first) = events.first() else {
        return EventAnnotationIndex {
            annotation,
            junctions: FxHashMap::default(),
        };
    };
    let wanted: FxHashSet<(u32, u32)> = events
        .iter()
        .flat_map(|event| {
            event
                .components
                .iter()
                .map(|component| (component.donor, component.acceptor))
        })
        .collect();
    let mut junctions: FxHashMap<(u32, u32), Vec<(u32, bool)>> = FxHashMap::default();
    if let Some(&annotation_chrom) = annotation.chrom_ids.get(&first.key.chrom) {
        for transcript in annotation
            .transcripts
            .iter()
            .filter(|transcript| transcript.chrom == annotation_chrom)
        {
            for exons in transcript.exons.windows(2) {
                let junction = (exons[0].end, exons[1].start);
                if wanted.contains(&junction) {
                    junctions
                        .entry(junction)
                        .or_default()
                        .push((transcript.gene, transcript.strand_rev));
                }
            }
        }
    }
    for labels in junctions.values_mut() {
        labels.sort_unstable();
        labels.dedup();
    }
    EventAnnotationIndex {
        annotation,
        junctions,
    }
}

fn event_annotation_json(
    index: &EventAnnotationIndex<'_>,
    event: &EventDefinition,
) -> serde_json::Value {
    let component_coordinates: Vec<(u32, u32)> = event
        .components
        .iter()
        .map(|component| (component.donor, component.acceptor))
        .collect();
    let mut genes: BTreeMap<u32, BTreeSet<bool>> = BTreeMap::new();
    for junction in &component_coordinates {
        if let Some(labels) = index.junctions.get(junction) {
            for &(gene, strand_rev) in labels {
                genes.entry(gene).or_default().insert(strand_rev);
            }
        }
    }
    let strands: BTreeSet<bool> = genes
        .values()
        .flat_map(|values| values.iter().copied())
        .collect();
    let strand = if strands.len() == 1 {
        if strands.contains(&true) { Some("-") } else { Some("+") }
    } else {
        None
    };
    json!({
        "genes": genes.keys().map(|gene| json!({
            "gene_id": index.annotation.gene_ids[*gene as usize],
            "gene_name": index.annotation.gene_names[*gene as usize],
        })).collect::<Vec<_>>(),
        "strand": strand,
        "fully_annotated": component_coordinates.iter().all(|junction| index.junctions.contains_key(junction)),
    })
}

fn event_annotation_tsv(
    index: Option<&EventAnnotationIndex<'_>>,
    event: &EventDefinition,
) -> (String, String, String) {
    let Some(index) = index else {
        return ("NA".to_owned(), "NA".to_owned(), "NA".to_owned());
    };
    let value = event_annotation_json(index, event);
    let genes = value["genes"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|gene| gene["gene_id"].as_str())
        .collect::<Vec<_>>()
        .join(",");
    let strand = value["strand"].as_str().unwrap_or("NA").to_owned();
    let fully_annotated = value["fully_annotated"].as_bool().unwrap().to_string();
    (genes, strand, fully_annotated)
}

fn event_component_json(
    event: &EventDefinition,
    component: &EventComponent,
    metadata: &[JunctionMeta],
) -> serde_json::Value {
    let mut value = json!({
        "locus": format!("{}:{}-{}", event.key.chrom, component.donor, component.acceptor),
        "donor": component.donor,
        "acceptor": component.acceptor,
        "present": component.metadata.is_some(),
    });
    if let Some(index) = component.metadata {
        let row = &metadata[index];
        let object = value.as_object_mut().unwrap();
        object.insert("supporting_children".into(), json!(row.supporting_children));
        object.insert("posting_chunks".into(), json!(row.posts.len()));
    }
    value
}

fn event_json(
    event: &EventDefinition,
    result: &EventResult,
    metadata: &[JunctionMeta],
    annotation: Option<&EventAnnotationIndex<'_>>,
    scope: &QueryScope,
    dictionary: &[u32],
    top: usize,
) -> serde_json::Value {
    let limit = if top == 0 {
        result.cells.len()
    } else {
        top.min(result.cells.len())
    };
    let mut value = json!({
        "id": event.key.id(),
        "event_type": event.key.kind.name(),
        "chrom": event.key.chrom,
        "present": event.catalogue_present,
        "inclusion_junctions": event.components.iter().filter(|component| component.side == 1).map(|component| event_component_json(event, component, metadata)).collect::<Vec<_>>(),
        "exclusion_junctions": event.components.iter().filter(|component| component.side == 2).map(|component| event_component_json(event, component, metadata)).collect::<Vec<_>>(),
        "totals": result.totals.json(),
        "cells": result.support_cells,
    });
    let object = value.as_object_mut().unwrap();
    if let Some(annotation) = annotation {
        object.insert(
            "annotation".into(),
            event_annotation_json(annotation, event),
        );
    }
    match scope.aggregation {
        QueryAggregation::Cell => {
            object.insert(
                "cell_rows".into(),
                json!(result
                    .cells
                    .iter()
                    .take(limit)
                    .map(|(cell, counts)| {
                        let mut row = counts.json();
                        row.as_object_mut()
                            .unwrap()
                            .insert("barcode".into(), json!(unpack_cell(dictionary, *cell)));
                        row
                    })
                    .collect::<Vec<_>>()),
            );
            object.insert(
                "cell_rows_truncated".into(),
                json!(limit < result.cells.len()),
            );
        }
        QueryAggregation::Group => {
            object.insert(
                "group_rows".into(),
                json!(result
                    .groups
                    .iter()
                    .enumerate()
                    .map(|(group, (counts, support_cells))| {
                        let mut row = counts.json();
                        let row_object = row.as_object_mut().unwrap();
                        row_object.insert("group".into(), json!(scope.group_names[group]));
                        row_object.insert("cells".into(), json!(support_cells));
                        row_object.insert(
                            "selected_cells".into(),
                            json!(scope.selected_per_group[group]),
                        );
                        row
                    })
                    .collect::<Vec<_>>()),
            );
        }
        QueryAggregation::Bulk => {
            object.insert("bulk".into(), result.totals.json());
        }
    }
    value
}

#[allow(clippy::too_many_arguments)]
fn run_events(
    archive: &Path,
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    chrom_names: &[String],
    locus: &str,
    event_types: &[EventTypeArg],
    min_support: u64,
    min_informative: usize,
    max_events: usize,
    top: usize,
    gtf: Option<&std::path::Path>,
    tsv: bool,
    json_output: bool,
    uniform_output: &UniformQueryOutputArgs,
    scope_args: &QueryScopeArgs,
    t0: std::time::Instant,
    t_open: f32,
) -> Result<()> {
    let (chrom, start, end) = parse_locus(locus)?;
    let chrom_id = chrom_names
        .iter()
        .position(|name| name == &chrom)
        .map(|index| index as u32)
        .with_context(|| format!("unknown chromosome {chrom}"))?;
    let mut scope = load_query_scope(la, scope_args)?;
    if uniform_output.format.is_some() {
        scope.ensure_resolved_mapping_digest()?;
    }
    let metadata = read_junction_metadata(la)?;
    let keys = discover_event_keys(
        &chrom,
        chrom_id,
        start,
        end,
        &metadata,
        event_types,
        min_support,
        max_events,
    )?;
    let events = prepare_event_definitions(&keys, chrom_names, &metadata)?;
    let (packed_hits, unique_chunk_decodes, independent_chunk_decodes) =
        event_packed_hits(la, chunks, &metadata, &events)?;
    let dictionary = la.cells()?.to_vec();
    let results = match scope.aggregation {
        QueryAggregation::Cell => reduce_event_results(
            la,
            packed_hits_by_event(packed_hits, events.len()),
            &scope,
            uniform_output.format.map(|_| dictionary.as_slice()),
        )?,
        QueryAggregation::Group | QueryAggregation::Bulk => {
            reduce_packed_event_results(la, &packed_hits, events.len(), &scope)?
        }
    };
    let retained: Vec<usize> = results
        .iter()
        .enumerate()
        .filter_map(|(index, result)| {
            (result.totals.informative() >= min_informative).then_some(index)
        })
        .collect();
    let (annotation, annotation_content_blake3) = match gtf {
        Some(path) => {
            let (annotation, digest) = load_query_annotation(
                path,
                "event annotation",
                uniform_output.format.is_some(),
            )?;
            (Some(annotation), digest)
        }
        None => (None, None),
    };
    let annotation_index = annotation
        .as_ref()
        .map(|annotation| build_event_annotation_index(annotation, &events));

    if uniform_output.format.is_some() {
        let selected_types: Vec<&str> = selected_event_types(event_types)
            .into_iter()
            .map(EventTypeArg::name)
            .collect();
        let summary = json!({
            "coordinates": "0-based half-open junction boundaries",
            "locus": locus,
            "chrom": chrom,
            "start": start,
            "end": end,
            "event_types": selected_types,
            "thresholds": {
                "min_archive_supporting_children": min_support,
                "min_scoped_informative_umis": min_informative,
                "max_candidate_events": max_events,
            },
            "semantics": {
                "candidate_source": "archive junction coordinates and support only",
                "class_categories": ["include_only", "exclude_only", "both"],
                "informative_umis": "include_only + exclude_only",
                "usage_fraction": "include_only / informative_umis",
                "both_in_usage_denominator": false,
                "alternative_site_side_order": "lower genomic coordinate is inclusion",
            },
            "scope": scope.json(),
            "planning": {
                "candidate_events": events.len(),
                "retained_events": retained.len(),
                "independent_chunk_decodes": independent_chunk_decodes,
                "unique_chunk_decodes": unique_chunk_decodes,
                "chunk_decode_reduction_fraction": if independent_chunk_decodes == 0 {
                    0.0
                } else {
                    1.0 - unique_chunk_decodes as f64 / independent_chunk_decodes as f64
                },
            },
        });
        let event_schema = TableSchema::new(
            "gravlax.query.events.events.v1",
            vec![
                Field::new("event_id", DataType::String),
                Field::new("event_type", DataType::String),
                Field::new("chrom", DataType::String),
                Field::new("catalogue_present", DataType::Boolean),
                Field::new("include_only", DataType::UInt64),
                Field::new("exclude_only", DataType::UInt64),
                Field::new("both", DataType::UInt64),
                Field::new("informative_umis", DataType::UInt64),
                Field::new("usage_fraction", DataType::Float64).nullable(),
                Field::new("supporting_cells", DataType::UInt64),
                Field::new("genes", DataType::Json).nullable(),
                Field::new("strand", DataType::String).nullable(),
                Field::new("fully_annotated", DataType::Boolean).nullable(),
                Field::new("available_count_rows", DataType::UInt64),
                Field::new("emitted_count_rows", DataType::UInt64),
                Field::new("count_rows_truncated", DataType::Boolean),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["event_id"]))?;
        let component_schema = TableSchema::new(
            "gravlax.query.events.components.v1",
            vec![
                Field::new("event_id", DataType::String),
                Field::new("side", DataType::String),
                Field::new("side_index", DataType::UInt64),
                Field::new("donor", DataType::UInt64),
                Field::new("acceptor", DataType::UInt64),
                Field::new("present", DataType::Boolean),
                Field::new("archive_supporting_children", DataType::UInt64).nullable(),
                Field::new("archive_posting_chunks", DataType::UInt64).nullable(),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key([
            "event_id",
            "side",
            "side_index",
        ]))?;
        let count_schema = TableSchema::new(
            "gravlax.query.events.counts.v1",
            vec![
                Field::new("event_id", DataType::String),
                Field::new("aggregation", DataType::String),
                Field::new("entity", DataType::String),
                Field::new("include_only", DataType::UInt64),
                Field::new("exclude_only", DataType::UInt64),
                Field::new("both", DataType::UInt64),
                Field::new("informative_umis", DataType::UInt64),
                Field::new("usage_fraction", DataType::Float64).nullable(),
                Field::new("cells", DataType::UInt64).nullable(),
                Field::new("selected_cells", DataType::UInt64).nullable(),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key([
            "event_id",
            "aggregation",
            "entity",
        ]))?;
        let event_selection = SelectionSummary::complete(retained.len() as u64);
        let component_rows: usize = retained
            .iter()
            .map(|&index| events[index].components.len())
            .sum();
        let component_selection = SelectionSummary::complete(component_rows as u64);
        let available_count_rows: usize = retained
            .iter()
            .map(|&index| match scope.aggregation {
                QueryAggregation::Cell => results[index].cells.len(),
                QueryAggregation::Group => results[index].groups.len(),
                QueryAggregation::Bulk => 1,
            })
            .sum();
        let emitted_count_rows: usize = retained
            .iter()
            .map(|&index| match scope.aggregation {
                QueryAggregation::Cell => {
                    if top == 0 {
                        results[index].cells.len()
                    } else {
                        top.min(results[index].cells.len())
                    }
                }
                QueryAggregation::Group => results[index].groups.len(),
                QueryAggregation::Bulk => 1,
            })
            .sum();
        let count_selection =
            SelectionSummary::selected(available_count_rows as u64, emitted_count_rows as u64)?;
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "archive_access".into(),
            json!("junction catalogue and union of candidate-event postings"),
        );
        parameters.insert("locus".into(), json!(locus));
        parameters.insert("event_types".into(), json!(selected_types));
        parameters.insert("min_archive_supporting_children".into(), json!(min_support));
        parameters.insert("min_scoped_informative_umis".into(), json!(min_informative));
        parameters.insert("max_candidate_events".into(), json!(max_events));
        parameters.insert("cell_scope".into(), scope.provenance_json());
        parameters.insert("aggregation".into(), json!(scope.aggregation_name()));
        parameters.insert(
            "selection_policy".into(),
            serde_json::to_value(UniformSelectionPolicy {
                requested_top: top,
                top_zero_means_all: true,
                comparator: matches!(scope.aggregation, QueryAggregation::Cell).then_some(
                    "total categorized UMIs descending, entity ascending (barcode), independently per event",
                ),
            })?,
        );
        if let Some(path) = gtf {
            parameters.insert("annotation_path".into(), json!(path));
            parameters.insert(
                "annotation_content_blake3".into(),
                json!(annotation_content_blake3
                    .as_deref()
                    .context("uniform event annotation is missing its bound content digest")?),
            );
        }
        let context = uniform_query_context(
            archive,
            la.reader().archive_version(),
            la.reader()
                .content_commitment()
                .map(|commitment| commitment.to_hex()),
            parameters,
        );
        write_uniform_bundle_output(uniform_output, |writer, format| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                "gravlax.query.events.result.v1",
                OutputFormat::from(format),
                &context,
                &summary,
            )?;
            bundle.write_table("events", &event_schema, Some(&event_selection), |rows| {
                for &index in &retained {
                    let event = &events[index];
                    let result = &results[index];
                    let (available, emitted) = match scope.aggregation {
                        QueryAggregation::Cell => {
                            let available = result.cells.len();
                            let emitted = if top == 0 {
                                available
                            } else {
                                top.min(available)
                            };
                            (available, emitted)
                        }
                        QueryAggregation::Group => (result.groups.len(), result.groups.len()),
                        QueryAggregation::Bulk => (1, 1),
                    };
                    let annotation_value = annotation_index
                        .as_ref()
                        .map(|annotation| event_annotation_json(annotation, event));
                    rows.write_row_with(|row| {
                        row.string(&event.key.id())?;
                        row.string(event.key.kind.name())?;
                        row.string(&event.key.chrom)?;
                        row.boolean(event.catalogue_present)?;
                        row.uint64(result.totals.include_only as u64)?;
                        row.uint64(result.totals.exclude_only as u64)?;
                        row.uint64(result.totals.both as u64)?;
                        row.uint64(result.totals.informative() as u64)?;
                        if let Some(usage) = result.totals.usage() {
                            row.float64(usage)?;
                        } else {
                            row.null()?;
                        }
                        row.uint64(result.support_cells as u64)?;
                        if let Some(annotation) = &annotation_value {
                            row.json(&annotation["genes"])?;
                            if let Some(strand) = annotation["strand"].as_str() {
                                row.string(strand)?;
                            } else {
                                row.null()?;
                            }
                            row.boolean(annotation["fully_annotated"].as_bool().unwrap_or(false))?;
                        } else {
                            row.null()?;
                            row.null()?;
                            row.null()?;
                        }
                        row.uint64(available as u64)?;
                        row.uint64(emitted as u64)?;
                        row.boolean(emitted < available)?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table(
                "components",
                &component_schema,
                Some(&component_selection),
                |rows| {
                    for &index in &retained {
                        let event = &events[index];
                        let mut inclusion_index = 0_u64;
                        let mut exclusion_index = 0_u64;
                        for component in &event.components {
                            let (side, side_index) = if component.side == 1 {
                                let value = inclusion_index;
                                inclusion_index += 1;
                                ("include", value)
                            } else {
                                let value = exclusion_index;
                                exclusion_index += 1;
                                ("exclude", value)
                            };
                            rows.write_row_with(|row| {
                                row.string(&event.key.id())?;
                                row.string(side)?;
                                row.uint64(side_index)?;
                                row.uint64(component.donor as u64)?;
                                row.uint64(component.acceptor as u64)?;
                                row.boolean(component.metadata.is_some())?;
                                if let Some(metadata_index) = component.metadata {
                                    let metadata = &metadata[metadata_index];
                                    row.uint64(metadata.supporting_children)?;
                                    row.uint64(metadata.posts.len() as u64)?;
                                } else {
                                    row.null()?;
                                    row.null()?;
                                }
                                Ok(())
                            })?;
                        }
                    }
                    Ok(())
                },
            )?;
            bundle.write_table("counts", &count_schema, Some(&count_selection), |rows| {
                let mut emit = |event_id: &str,
                                aggregation: &str,
                                entity: &str,
                                counts: &JunctionSetCounts,
                                cells: Option<usize>,
                                selected_cells: Option<usize>|
                 -> std::result::Result<(), OutputError> {
                    rows.write_row_with(|row| {
                        row.string(event_id)?;
                        row.string(aggregation)?;
                        row.string(entity)?;
                        row.uint64(counts.include_only as u64)?;
                        row.uint64(counts.exclude_only as u64)?;
                        row.uint64(counts.both as u64)?;
                        row.uint64(counts.informative() as u64)?;
                        if let Some(usage) = counts.usage() {
                            row.float64(usage)?;
                        } else {
                            row.null()?;
                        }
                        if let Some(cells) = cells {
                            row.uint64(cells as u64)?;
                        } else {
                            row.null()?;
                        }
                        if let Some(selected_cells) = selected_cells {
                            row.uint64(selected_cells as u64)?;
                        } else {
                            row.null()?;
                        }
                        Ok(())
                    })
                };
                for &index in &retained {
                    let event_id = events[index].key.id();
                    let result = &results[index];
                    match scope.aggregation {
                        QueryAggregation::Cell => {
                            let limit = if top == 0 {
                                result.cells.len()
                            } else {
                                top.min(result.cells.len())
                            };
                            for (cell, counts) in result.cells.iter().take(limit) {
                                let barcode = unpack_cell_bytes(dictionary[*cell as usize]);
                                let barcode = std::str::from_utf8(&barcode)
                                    .expect("packed archive barcode decodes to ASCII");
                                emit(&event_id, "cell", barcode, counts, None, None)?;
                            }
                        }
                        QueryAggregation::Group => {
                            for (group, (counts, cells)) in result.groups.iter().enumerate() {
                                emit(
                                    &event_id,
                                    "group",
                                    &scope.group_names[group],
                                    counts,
                                    Some(*cells),
                                    Some(scope.selected_per_group[group]),
                                )?;
                            }
                        }
                        QueryAggregation::Bulk => emit(
                            &event_id,
                            "bulk",
                            "bulk",
                            &result.totals,
                            Some(result.support_cells),
                            Some(scope.selected_cells),
                        )?,
                    }
                }
                Ok(())
            })?;
            bundle.finish()?;
            Ok(())
        })?;
    } else if tsv {
        println!("event_id\tevent_type\tchrom\tinclude_junctions\texclude_junctions\tgenes\tstrand\tfully_annotated\taggregation\tlabel\tinclude_only\texclude_only\tboth\tinformative_umis\tusage_fraction\tcells\tselected_cells");
        for &index in &retained {
            let event = &events[index];
            let result = &results[index];
            let includes = event
                .key
                .includes
                .iter()
                .map(|(donor, acceptor)| format!("{chrom}:{donor}-{acceptor}"))
                .collect::<Vec<_>>()
                .join(",");
            let excludes = event
                .key
                .excludes
                .iter()
                .map(|(donor, acceptor)| format!("{chrom}:{donor}-{acceptor}"))
                .collect::<Vec<_>>()
                .join(",");
            let (genes, strand, fully_annotated) =
                event_annotation_tsv(annotation_index.as_ref(), event);
            let emit = |aggregation: &str,
                        label: &str,
                        counts: &JunctionSetCounts,
                        cells: usize,
                        selected_cells: usize| {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    event.key.id(),
                    event.key.kind.name(),
                    chrom,
                    includes,
                    excludes,
                    genes,
                    strand,
                    fully_annotated,
                    aggregation,
                    label,
                    counts.include_only,
                    counts.exclude_only,
                    counts.both,
                    counts.informative(),
                    counts
                        .usage()
                        .map_or_else(|| "NA".to_owned(), |usage| format!("{usage:.9}")),
                    cells,
                    selected_cells,
                );
            };
            match scope.aggregation {
                QueryAggregation::Cell => {
                    let limit = if top == 0 {
                        result.cells.len()
                    } else {
                        top.min(result.cells.len())
                    };
                    for (cell, counts) in result.cells.iter().take(limit) {
                        emit(
                            "cell",
                            &unpack_cell(&dictionary, *cell),
                            counts,
                            usize::from(counts.total() > 0),
                            1,
                        );
                    }
                }
                QueryAggregation::Group => {
                    for (group, (counts, cells)) in result.groups.iter().enumerate() {
                        emit(
                            "group",
                            &scope.group_names[group],
                            counts,
                            *cells,
                            scope.selected_per_group[group],
                        );
                    }
                }
                QueryAggregation::Bulk => emit(
                    "bulk",
                    "bulk",
                    &result.totals,
                    result.support_cells,
                    scope.selected_cells,
                ),
            }
        }
    } else if json_output {
        let value = json!({
            "schema": "gravlax.query.events.v1",
            "coordinates": "0-based half-open junction boundaries",
            "semantics": {
                "candidate_source": "archive junction coordinates and support only",
                "class_categories": ["include_only", "exclude_only", "both"],
                "informative_umis": "include_only + exclude_only",
                "usage_fraction": "include_only / informative_umis",
                "both_in_usage_denominator": false,
                "alternative_site_side_order": "lower genomic coordinate is inclusion",
            },
            "locus": locus,
            "min_support": min_support,
            "min_informative": min_informative,
            "scope": scope.json(),
            "planning": {
                "candidate_events": events.len(),
                "retained_events": retained.len(),
                "independent_chunk_decodes": independent_chunk_decodes,
                "unique_chunk_decodes": unique_chunk_decodes,
                "chunk_decode_reduction_fraction": if independent_chunk_decodes == 0 { 0.0 } else { 1.0 - unique_chunk_decodes as f64 / independent_chunk_decodes as f64 },
            },
            "events": retained.iter().map(|&index| event_json(&events[index], &results[index], &metadata, annotation_index.as_ref(), &scope, &dictionary, top)).collect::<Vec<_>>(),
        });
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        println!(
            "events {locus}: {} retained / {} candidates ({} unique / {} independent chunk decodes; open {t_open:.2}s, total {:.2}s)",
            retained.len(), events.len(), unique_chunk_decodes, independent_chunk_decodes,
            t0.elapsed().as_secs_f32(),
        );
        for &index in retained.iter().take(20) {
            let result = &results[index];
            println!(
                "  {}\t{} informative\tusage {}\t{} cells",
                events[index].key.id(), result.totals.informative(),
                result.totals.usage().map_or_else(|| "NA".to_owned(), |usage| format!("{usage:.4}")),
                result.support_cells,
            );
        }
    }
    eprintln!(
        "event engine: {} events retained from {} candidates; open {t_open:.2}s, total {:.2}s",
        retained.len(), events.len(), t0.elapsed().as_secs_f32(),
    );
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct GraphPathKey {
    strand_rev: bool,
    junctions: Vec<(u32, u32)>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct GraphEdgeKey {
    strand_rev: bool,
    donor: u32,
    acceptor: u32,
}

struct GraphAggregate {
    umis: usize,
    cells: FxHashSet<u32>,
    group_umis: Vec<usize>,
    group_cells: Vec<FxHashSet<u32>>,
}

impl GraphAggregate {
    fn new(group_count: usize) -> Self {
        Self {
            umis: 0,
            cells: FxHashSet::default(),
            group_umis: vec![0; group_count],
            group_cells: (0..group_count).map(|_| FxHashSet::default()).collect(),
        }
    }

    fn add_class(&mut self, cell: u32, group: Option<u32>) {
        self.umis += 1;
        self.cells.insert(cell);
        if let Some(group) = group {
            let group = group as usize;
            self.group_umis[group] += 1;
            self.group_cells[group].insert(cell);
        }
    }

    fn merge(&mut self, other: &Self) {
        self.umis += other.umis;
        self.cells.extend(other.cells.iter().copied());
        for group in 0..self.group_umis.len() {
            self.group_umis[group] += other.group_umis[group];
            self.group_cells[group].extend(other.group_cells[group].iter().copied());
        }
    }

    fn group_json(&self, scope: &QueryScope) -> Vec<serde_json::Value> {
        scope
            .group_names
            .iter()
            .enumerate()
            .map(|(group, name)| {
                json!({
                    "group": name,
                    "umis": self.group_umis[group],
                    "cells": self.group_cells[group].len(),
                    "selected_cells": scope.selected_per_group[group],
                })
            })
            .collect()
    }
}

struct GraphPathReduction {
    paths: BTreeMap<GraphPathKey, GraphAggregate>,
    scoped_distinct_umi_classes: usize,
    unique_chunk_decodes: usize,
    independent_chunk_decodes: usize,
}

type GraphPathHit = (u32, bool, Vec<(u32, u32)>);

fn reduce_graph_paths(
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    chrom_id: u32,
    selected_metadata: &[&JunctionMeta],
    scope: &QueryScope,
    max_paths: usize,
) -> Result<GraphPathReduction> {
    let wanted: FxHashSet<(u32, u32)> = selected_metadata
        .iter()
        .map(|row| (row.donor, row.acceptor))
        .collect();
    let independent_chunk_decodes: usize = selected_metadata.iter().map(|row| row.posts.len()).sum();
    let mut selected_chunks = BTreeSet::new();
    for row in selected_metadata {
        for &post in &row.posts {
            let chunk = chunks.get(post as usize).with_context(|| {
                format!(
                    "junction {}-{} references missing chunk {post}",
                    row.donor, row.acceptor
                )
            })?;
            if chunk.chrom != chrom_id {
                bail!(
                    "junction {}-{} posting references a different chromosome",
                    row.donor,
                    row.acceptor
                );
            }
            selected_chunks.insert(post as usize);
        }
    }
    let selected_chunks: Vec<usize> = selected_chunks.into_iter().collect();
    let shapes = la.shapes()?;
    let chunk_hits: Vec<Vec<GraphPathHit>> = {
        let (reader, tables) = la.reader_and_tables();
        let reader = &*reader;
        selected_chunks
            .par_iter()
            .map(
                |&chunk_index| -> Result<Vec<GraphPathHit>> {
                    let (compressed, raw_len) =
                        reader.read_compressed_at(&format!("c{chunk_index}"))?;
                    let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                    let molecules = decode_chunk(&raw, &chunks[chunk_index], None, tables)?;
                    let mut hits = Vec::new();
                    for molecule in molecules {
                        let mut junctions = Vec::new();
                        let mut inspect = |position: u32, shape: u32| {
                            for blocks in shapes[shape as usize].blocks.windows(2) {
                                let donor = position + blocks[0].0 + blocks[0].1;
                                let acceptor = position + blocks[1].0;
                                if wanted.contains(&(donor, acceptor)) {
                                    junctions.push((donor, acceptor));
                                }
                            }
                        };
                        for chain in &molecule.chains {
                            for &(position, shape) in &chain.reps {
                                inspect(position, shape);
                            }
                        }
                        for &(position, shape, _, _) in &molecule.mms {
                            inspect(position, shape);
                        }
                        if !junctions.is_empty() {
                            junctions.sort_unstable();
                            junctions.dedup();
                            hits.push((molecule.umi_class, molecule.strand_rev, junctions));
                        }
                    }
                    Ok(hits)
                },
            )
            .collect::<Result<_>>()?
    };

    let mut hits: Vec<GraphPathHit> = chunk_hits.into_iter().flatten().collect();
    hits.sort_unstable_by_key(|(class, strand_rev, _)| (*class, *strand_rev));
    let mut merged_hits: Vec<GraphPathHit> = Vec::with_capacity(hits.len());
    for (class, strand_rev, junctions) in hits {
        match merged_hits.last_mut() {
            Some((previous_class, previous_strand, previous_junctions))
                if *previous_class == class && *previous_strand == strand_rev =>
            {
                previous_junctions.extend(junctions);
                previous_junctions.sort_unstable();
                previous_junctions.dedup();
            }
            _ => merged_hits.push((class, strand_rev, junctions)),
        }
    }
    la.prefetch_coc(merged_hits.iter().map(|(class, _, _)| *class))?;

    let mut paths: BTreeMap<GraphPathKey, GraphAggregate> = BTreeMap::new();
    let mut distinct_classes = FxHashSet::default();
    for (class, strand_rev, junctions) in merged_hits {
        let cell = la.cell_of(class)?;
        if !scope.includes(cell) {
            continue;
        }
        distinct_classes.insert(class);
        let path = GraphPathKey {
            strand_rev,
            junctions,
        };
        if !paths.contains_key(&path) && paths.len() == max_paths {
            bail!(
                "splice graph has more than {max_paths} candidate paths, exceeding --max-paths {max_paths}; narrow the locus or increase support thresholds"
            );
        }
        paths
            .entry(path)
            .or_insert_with(|| GraphAggregate::new(scope.group_names.len()))
            .add_class(cell, scope.group_of.get(&cell).copied());
    }
    Ok(GraphPathReduction {
        paths,
        scoped_distinct_umi_classes: distinct_classes.len(),
        unique_chunk_decodes: selected_chunks.len(),
        independent_chunk_decodes,
    })
}

fn aggregate_graph_edges<'a, I>(
    paths: I,
    group_count: usize,
) -> BTreeMap<GraphEdgeKey, GraphAggregate>
where
    I: IntoIterator<Item = (&'a GraphPathKey, &'a GraphAggregate)>,
{
    let mut edge_counts = BTreeMap::new();
    for (path, counts) in paths {
        for &(donor, acceptor) in &path.junctions {
            edge_counts
                .entry(GraphEdgeKey {
                    strand_rev: path.strand_rev,
                    donor,
                    acceptor,
                })
                .or_insert_with(|| GraphAggregate::new(group_count))
                .merge(counts);
        }
    }
    edge_counts
}

fn ordered_graph_junctions(path: &GraphPathKey) -> Vec<(u32, u32)> {
    let mut junctions = path.junctions.clone();
    if path.strand_rev {
        junctions.sort_unstable_by(|left, right| right.cmp(left));
    }
    junctions
}

#[allow(clippy::too_many_arguments)]
fn run_splice_graph(
    archive: &Path,
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    chrom_names: &[String],
    locus: &str,
    min_support: u64,
    min_path_umis: usize,
    max_paths: usize,
    json_output: bool,
    uniform_output: &UniformQueryOutputArgs,
    scope_args: &QueryScopeArgs,
    t0: std::time::Instant,
    t_open: f32,
) -> Result<()> {
    if min_support == 0 {
        bail!("--min-support must be at least 1");
    }
    if min_path_umis == 0 {
        bail!("--min-path-umis must be at least 1");
    }
    if max_paths == 0 {
        bail!("--max-paths must be at least 1");
    }
    let mut scope = load_query_scope(la, scope_args)?;
    if uniform_output.format.is_some() {
        scope.ensure_resolved_mapping_digest()?;
    }
    let (chrom, start, end) = parse_locus(locus)?;
    let chrom_id = chrom_names
        .iter()
        .position(|name| name == &chrom)
        .map(|index| index as u32)
        .with_context(|| format!("unknown chromosome {chrom}"))?;
    let metadata = read_junction_metadata(la)?;
    let selected_metadata: Vec<&JunctionMeta> = metadata
        .iter()
        .filter(|row| {
            row.chrom == chrom_id
                && row.donor >= start
                && row.acceptor < end
                && row.supporting_children >= min_support
        })
        .collect();
    let coordinate_metadata: FxHashMap<(u32, u32), (u64, usize)> = selected_metadata
        .iter()
        .map(|row| {
            (
                (row.donor, row.acceptor),
                (row.supporting_children, row.posts.len()),
            )
        })
        .collect();
    let reduction =
        reduce_graph_paths(la, chunks, chrom_id, &selected_metadata, &scope, max_paths)?;
    let candidate_paths = reduction.paths;
    let candidate_path_count = candidate_paths.len();
    let candidate_strand_path_umis: usize =
        candidate_paths.values().map(|counts| counts.umis).sum();
    let mut paths: Vec<(GraphPathKey, GraphAggregate)> = candidate_paths
        .into_iter()
        .filter(|(_, counts)| counts.umis >= min_path_umis)
        .collect();
    paths.sort_by(|(left_key, left_counts), (right_key, right_counts)| {
        std::cmp::Reverse(left_counts.umis)
            .cmp(&std::cmp::Reverse(right_counts.umis))
            .then_with(|| left_key.cmp(right_key))
    });

    let edge_counts = aggregate_graph_edges(
        paths.iter().map(|(path, counts)| (path, counts)),
        scope.group_names.len(),
    );
    let mut node_keys = BTreeSet::new();
    for edge in edge_counts.keys() {
        node_keys.insert((edge.strand_rev, edge.donor));
        node_keys.insert((edge.strand_rev, edge.acceptor));
    }
    let node_ids: BTreeMap<(bool, u32), usize> = node_keys
        .into_iter()
        .enumerate()
        .map(|(id, key)| (key, id))
        .collect();
    let edge_ids: BTreeMap<GraphEdgeKey, usize> = edge_counts
        .keys()
        .copied()
        .enumerate()
        .map(|(id, key)| (key, id))
        .collect();
    let strand_path_umis: usize = paths.iter().map(|(_, counts)| counts.umis).sum();

    if uniform_output.format.is_some() {
        let summary = json!({
            "coordinates": "0-based half-open",
            "chrom": chrom,
            "start": start,
            "end": end,
            "scope": scope.json(),
            "semantics": {
                "path_fragment": "selected junction set co-supported by one archive UMI class and strand",
                "lower_bound": true,
                "complete_transcript_claim": false,
                "multimapper_policy": "archived anchor placement; alternatives are not resolved",
                "edge_count_basis": "retained path fragments",
                "archive_catalogue_support_is_strand_combined": true,
            },
            "thresholds": {
                "min_archive_supporting_children": min_support,
                "min_scoped_path_umis": min_path_umis,
                "max_candidate_paths": max_paths,
            },
            "totals": {
                "nodes": node_ids.len(),
                "edges": edge_counts.len(),
                "paths": paths.len(),
                "strand_path_umis": strand_path_umis,
            },
            "planning": {
                "catalogue_junctions": selected_metadata.len(),
                "candidate_paths": candidate_path_count,
                "scoped_distinct_umi_classes": reduction.scoped_distinct_umi_classes,
                "candidate_strand_path_umis": candidate_strand_path_umis,
                "retained_paths": paths.len(),
                "unique_chunk_decodes": reduction.unique_chunk_decodes,
                "independent_chunk_decodes": reduction.independent_chunk_decodes,
            },
        });
        let node_schema = TableSchema::new(
            "gravlax.query.splice-graph.nodes.v1",
            vec![
                Field::new("node_id", DataType::UInt64),
                Field::new("coordinate", DataType::UInt64),
                Field::new("strand", DataType::String),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["node_id"])
                .ordered_by([OrderKey::ascending("node_id")]),
        )?;
        let edge_schema = TableSchema::new(
            "gravlax.query.splice-graph.edges.v1",
            vec![
                Field::new("edge_id", DataType::UInt64),
                Field::new("strand", DataType::String),
                Field::new("donor", DataType::UInt64),
                Field::new("acceptor", DataType::UInt64),
                Field::new("source_node_id", DataType::UInt64),
                Field::new("target_node_id", DataType::UInt64),
                Field::new("umis", DataType::UInt64),
                Field::new("cells", DataType::UInt64),
                Field::new("archive_supporting_children", DataType::UInt64),
                Field::new("archive_posting_chunks", DataType::UInt64),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["edge_id"])
                .ordered_by([OrderKey::ascending("edge_id")]),
        )?;
        let path_schema = TableSchema::new(
            "gravlax.query.splice-graph.paths.v1",
            vec![
                Field::new("path_id", DataType::UInt64),
                Field::new("strand", DataType::String),
                Field::new("edge_ids", DataType::Json),
                Field::new("junctions", DataType::Json),
                Field::new("umis", DataType::UInt64),
                Field::new("cells", DataType::UInt64),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["path_id"])
                .ordered_by([OrderKey::ascending("path_id")]),
        )?;
        let group_schema = TableSchema::new(
            "gravlax.query.splice-graph.group-counts.v1",
            vec![
                Field::new("object_kind", DataType::String),
                Field::new("object_id", DataType::UInt64),
                Field::new("group", DataType::String),
                Field::new("umis", DataType::UInt64),
                Field::new("cells", DataType::UInt64),
                Field::new("selected_cells", DataType::UInt64),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key([
            "object_kind",
            "object_id",
            "group",
        ]))?;
        let node_selection = SelectionSummary::complete(node_ids.len() as u64);
        let edge_selection = SelectionSummary::complete(edge_counts.len() as u64);
        let path_selection = SelectionSummary::complete(paths.len() as u64);
        let group_rows = scope
            .group_names
            .len()
            .saturating_mul(edge_counts.len().saturating_add(paths.len()));
        let group_selection = SelectionSummary::complete(group_rows as u64);
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "archive_access".into(),
            json!("junction catalogue and union of selected-junction postings"),
        );
        parameters.insert("locus".into(), json!(locus));
        parameters.insert("min_archive_supporting_children".into(), json!(min_support));
        parameters.insert("min_scoped_path_umis".into(), json!(min_path_umis));
        parameters.insert("max_candidate_paths".into(), json!(max_paths));
        parameters.insert("cell_scope".into(), scope.provenance_json());
        parameters.insert("aggregation".into(), json!(scope.aggregation_name()));
        parameters.insert(
            "path_order".into(),
            json!("umis descending, then strand and junction coordinates ascending"),
        );
        let context = uniform_query_context(
            archive,
            la.reader().archive_version(),
            la.reader()
                .content_commitment()
                .map(|commitment| commitment.to_hex()),
            parameters,
        );
        write_uniform_bundle_output(uniform_output, |writer, format| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                "gravlax.query.splice-graph.result.v1",
                OutputFormat::from(format),
                &context,
                &summary,
            )?;
            bundle.write_table("nodes", &node_schema, Some(&node_selection), |rows| {
                for (&(strand_rev, coordinate), &id) in &node_ids {
                    rows.write_row_with(|row| {
                        row.uint64(id as u64)?;
                        row.uint64(coordinate as u64)?;
                        row.string(if strand_rev { "-" } else { "+" })?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table("edges", &edge_schema, Some(&edge_selection), |rows| {
                for (edge, counts) in &edge_counts {
                    let id = edge_ids[edge];
                    let source_coordinate = if edge.strand_rev {
                        edge.acceptor
                    } else {
                        edge.donor
                    };
                    let target_coordinate = if edge.strand_rev {
                        edge.donor
                    } else {
                        edge.acceptor
                    };
                    let (supporting_children, posting_chunks) = coordinate_metadata
                        .get(&(edge.donor, edge.acceptor))
                        .copied()
                        .unwrap_or((0, 0));
                    rows.write_row_with(|row| {
                        row.uint64(id as u64)?;
                        row.string(if edge.strand_rev { "-" } else { "+" })?;
                        row.uint64(edge.donor as u64)?;
                        row.uint64(edge.acceptor as u64)?;
                        row.uint64(node_ids[&(edge.strand_rev, source_coordinate)] as u64)?;
                        row.uint64(node_ids[&(edge.strand_rev, target_coordinate)] as u64)?;
                        row.uint64(counts.umis as u64)?;
                        row.uint64(counts.cells.len() as u64)?;
                        row.uint64(supporting_children)?;
                        row.uint64(posting_chunks as u64)?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table("paths", &path_schema, Some(&path_selection), |rows| {
                for (id, (path, counts)) in paths.iter().enumerate() {
                    let junctions = ordered_graph_junctions(path);
                    let path_edges = json!(junctions
                        .iter()
                        .map(|&(donor, acceptor)| edge_ids[&GraphEdgeKey {
                            strand_rev: path.strand_rev,
                            donor,
                            acceptor,
                        }])
                        .collect::<Vec<_>>());
                    let junction_values = json!(junctions
                        .iter()
                        .map(|(donor, acceptor)| json!({
                            "donor": donor,
                            "acceptor": acceptor,
                        }))
                        .collect::<Vec<_>>());
                    rows.write_row_with(|row| {
                        row.uint64(id as u64)?;
                        row.string(if path.strand_rev { "-" } else { "+" })?;
                        row.json(&path_edges)?;
                        row.json(&junction_values)?;
                        row.uint64(counts.umis as u64)?;
                        row.uint64(counts.cells.len() as u64)?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            if !scope.group_names.is_empty() {
                bundle.write_table(
                    "group_counts",
                    &group_schema,
                    Some(&group_selection),
                    |rows| {
                        for (edge, counts) in &edge_counts {
                            let id = edge_ids[edge];
                            for (group, name) in scope.group_names.iter().enumerate() {
                                rows.write_row_with(|row| {
                                    row.string("edge")?;
                                    row.uint64(id as u64)?;
                                    row.string(name)?;
                                    row.uint64(counts.group_umis[group] as u64)?;
                                    row.uint64(counts.group_cells[group].len() as u64)?;
                                    row.uint64(scope.selected_per_group[group] as u64)?;
                                    Ok(())
                                })?;
                            }
                        }
                        for (id, (_, counts)) in paths.iter().enumerate() {
                            for (group, name) in scope.group_names.iter().enumerate() {
                                rows.write_row_with(|row| {
                                    row.string("path")?;
                                    row.uint64(id as u64)?;
                                    row.string(name)?;
                                    row.uint64(counts.group_umis[group] as u64)?;
                                    row.uint64(counts.group_cells[group].len() as u64)?;
                                    row.uint64(scope.selected_per_group[group] as u64)?;
                                    Ok(())
                                })?;
                            }
                        }
                        Ok(())
                    },
                )?;
            }
            bundle.finish()?;
            Ok(())
        })?;
        eprintln!(
            "splice graph: {} retained / {} candidate paths; open {t_open:.2}s, total {:.2}s",
            paths.len(),
            candidate_path_count,
            t0.elapsed().as_secs_f32(),
        );
        return Ok(());
    }

    let nodes: Vec<serde_json::Value> = node_ids
        .iter()
        .map(|(&(strand_rev, coordinate), &id)| {
            json!({
                "id": id,
                "coordinate": coordinate,
                "strand": if strand_rev { "-" } else { "+" },
            })
        })
        .collect();
    let edges: Vec<serde_json::Value> = edge_counts
        .iter()
        .map(|(edge, counts)| {
            let id = edge_ids[edge];
            let source_coordinate = if edge.strand_rev {
                edge.acceptor
            } else {
                edge.donor
            };
            let target_coordinate = if edge.strand_rev {
                edge.donor
            } else {
                edge.acceptor
            };
            let (supporting_children, posting_chunks) = coordinate_metadata
                .get(&(edge.donor, edge.acceptor))
                .copied()
                .unwrap_or((0, 0));
            let mut value = json!({
                "id": id,
                "strand": if edge.strand_rev { "-" } else { "+" },
                "donor": edge.donor,
                "acceptor": edge.acceptor,
                "source": node_ids[&(edge.strand_rev, source_coordinate)],
                "target": node_ids[&(edge.strand_rev, target_coordinate)],
                "umis": counts.umis,
                "cells": counts.cells.len(),
                "catalogue_supporting_children": supporting_children,
                "catalogue_posting_chunks": posting_chunks,
            });
            if !scope.group_names.is_empty() {
                value
                    .as_object_mut()
                    .unwrap()
                    .insert("group_counts".into(), json!(counts.group_json(&scope)));
            }
            value
        })
        .collect();
    let path_rows: Vec<serde_json::Value> = paths
        .iter()
        .enumerate()
        .map(|(id, (path, counts))| {
            let ordered = ordered_graph_junctions(path);
            let path_edges: Vec<usize> = ordered
                .iter()
                .map(|&(donor, acceptor)| {
                    edge_ids[&GraphEdgeKey {
                        strand_rev: path.strand_rev,
                        donor,
                        acceptor,
                    }]
                })
                .collect();
            let mut value = json!({
                "id": id,
                "strand": if path.strand_rev { "-" } else { "+" },
                "edge_ids": path_edges,
                "junctions": ordered.iter().map(|(donor, acceptor)| json!({
                    "donor": donor,
                    "acceptor": acceptor,
                })).collect::<Vec<_>>(),
                "umis": counts.umis,
                "cells": counts.cells.len(),
            });
            if !scope.group_names.is_empty() {
                value
                    .as_object_mut()
                    .unwrap()
                    .insert("group_counts".into(), json!(counts.group_json(&scope)));
            }
            value
        })
        .collect();
    let value = json!({
        "schema": "gravlax.query.splice-graph.v1",
        "coordinates": "0-based half-open",
        "chrom": chrom,
        "start": start,
        "end": end,
        "scope": scope.json(),
        "semantics": {
            "path_fragment": "selected junction set co-supported by one archive UMI class and strand",
            "lower_bound": true,
            "complete_transcript_claim": false,
            "multimapper_policy": "archived anchor placement; alternatives are not resolved",
            "edge_count_basis": "retained path fragments",
            "catalogue_support_is_strand_combined": true,
        },
        "thresholds": {
            "min_support": min_support,
            "min_path_umis": min_path_umis,
            "max_paths": max_paths,
        },
        "totals": {
            "nodes": nodes.len(),
            "edges": edges.len(),
            "paths": path_rows.len(),
            "strand_path_umis": strand_path_umis,
        },
        "planning": {
            "catalogue_junctions": selected_metadata.len(),
            "candidate_paths": candidate_path_count,
            "scoped_distinct_umi_classes": reduction.scoped_distinct_umi_classes,
            "candidate_strand_path_umis": candidate_strand_path_umis,
            "retained_paths": path_rows.len(),
            "unique_chunk_decodes": reduction.unique_chunk_decodes,
            "independent_chunk_decodes": reduction.independent_chunk_decodes,
        },
        "nodes": nodes,
        "edges": edges,
        "paths": path_rows,
    });
    if json_output {
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        println!(
            "splice-graph {locus}: {} nodes, {} edges, {} paths, {} retained strand-path UMI counts ({} unique / {} independent chunk decodes; open {t_open:.2}s, total {:.2}s)",
            value["totals"]["nodes"],
            value["totals"]["edges"],
            value["totals"]["paths"],
            value["totals"]["strand_path_umis"],
            reduction.unique_chunk_decodes,
            reduction.independent_chunk_decodes,
            t0.elapsed().as_secs_f32(),
        );
        for path in value["paths"].as_array().into_iter().flatten().take(20) {
            println!(
                "  path {}	{}	{} UMIs	{} cells",
                path["id"],
                path["strand"].as_str().unwrap_or("?"),
                path["umis"],
                path["cells"]
            );
        }
    }
    eprintln!(
        "splice graph: {} retained / {} candidate paths; open {t_open:.2}s, total {:.2}s",
        value["planning"]["retained_paths"],
        value["planning"]["candidate_paths"],
        t0.elapsed().as_secs_f32(),
    );
    Ok(())
}

const TRANSCRIPT_EC_QUERY_SCHEMA: &str = "gravlax.query.transcript-ecs.v1";
const TRANSCRIPT_EC_CATALOG_SCHEMA: &str = "gravlax.query.transcript-ecs.catalog.v1";
const TRANSCRIPT_EC_COUNTS_SCHEMA: &str = "gravlax.query.transcript-ecs.counts.v1";
const TRANSCRIPT_EC_MEMBERSHIP_SCHEMA: &str = "gravlax.query.transcript-ecs.membership.v1";
const TRANSCRIPT_EC_MAX_COUNT_ROWS: usize = 1_000_000;

#[derive(Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum TranscriptEcSelection {
    Gene {
        requested: String,
        stable_id: String,
        display_name: Option<String>,
        matched_by: anno::intent::MatchBasis,
        loci: Vec<anno::intent::GenomicLocus>,
    },
    Locus {
        requested: String,
        contig: String,
        start: u32,
        end: u32,
    },
}

#[derive(Serialize)]
struct TranscriptEcQueryScope {
    selection: TranscriptEcSelection,
    selected_transcript_count: usize,
    selected_transcript_ids: Vec<String>,
    cells: serde_json::Value,
    archive_access: &'static str,
    chunk_pruning_applied: bool,
}

#[derive(Serialize)]
struct TranscriptEcQuerySemantics {
    compatibility: crate::archivecmd::transcriptec::TranscriptEquivalenceSemantics,
    count_unit: &'static str,
    aggregation: &'static str,
    inference: &'static str,
    phasing: &'static str,
}

struct TranscriptEcCatalogRow {
    ec_id: String,
    transcript_ids: Vec<String>,
    gene_ids: Vec<String>,
    ambiguous: bool,
    archived_umi_class_count: u64,
    cell_count: u64,
    complete_umi_class_count: u64,
}

struct TranscriptEcCountRow {
    aggregation: &'static str,
    key: String,
    cell_id: Option<u32>,
    ec_id: Option<String>,
    archived_umi_class_count: u64,
    ambiguous_umi_class_count: u64,
    no_compatible_transcript_umi_class_count: u64,
    conflicting_umi_class_count: u64,
    complete_umi_class_count: u64,
    incomplete_umi_class_count: u64,
    retained_record_count: u64,
    represented_alignment_count: u64,
}

struct TranscriptEcMembershipRow {
    umi_class: u32,
    cell_id: u32,
    barcode: String,
    aggregation: &'static str,
    key: String,
    ec_id: Option<String>,
    retained_record_count: u64,
    represented_alignment_count: u64,
    compatible_record_count: u64,
    unmatched_record_count: u64,
    ambiguous: bool,
    no_compatible_transcript: bool,
    conflict: bool,
    complete_within_archive_quotient: bool,
    retained_representatives_complete: bool,
}

#[derive(Serialize)]
struct TranscriptEcQuerySummary {
    archive_version: u32,
    archive_root_blake3: Option<String>,
    strand_policy: crate::archivecmd::transcriptec::TranscriptStrandPolicy,
    archive_umi_classes_scanned: u64,
    selector_relevant_umi_classes: u64,
    scoped_umi_classes: u64,
    scoped_cells_with_classes: u64,
    transcript_ecs: u64,
    count_rows: u64,
    membership_rows: u64,
    assigned_umi_classes: u64,
    unassigned_umi_classes: u64,
    ambiguous_umi_classes: u64,
    no_compatible_transcript_umi_classes: u64,
    conflicting_umi_classes: u64,
    complete_umi_classes: u64,
    incomplete_umi_classes: u64,
}

#[derive(Serialize)]
struct TranscriptEcResultData {
    scope: TranscriptEcQueryScope,
    semantics: TranscriptEcQuerySemantics,
    summary: TranscriptEcQuerySummary,
    catalog: TypedTable,
    counts: TypedTable,
    #[serde(skip_serializing_if = "Option::is_none")]
    membership: Option<TypedTable>,
}

#[derive(Default)]
struct TranscriptEcCatalogAccumulator {
    cells: BTreeSet<u32>,
    class_count: u64,
    complete_class_count: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct TranscriptEcCountKey {
    key: String,
    cell_id: Option<u32>,
    ec_id: Option<String>,
}

#[derive(Default)]
struct TranscriptEcCountAccumulator {
    class_count: u64,
    ambiguous_class_count: u64,
    no_compatible_class_count: u64,
    conflicting_class_count: u64,
    complete_class_count: u64,
    incomplete_class_count: u64,
    retained_record_count: u64,
    represented_alignment_count: u64,
}

fn validate_transcript_ec_output_args(
    format: TranscriptEcFormat,
    table: Option<TranscriptEcTable>,
    emit_membership: bool,
) -> Result<()> {
    match (format, table) {
        (TranscriptEcFormat::Tsv, None) => {
            bail!("--format tsv requires --table catalog, counts, or membership")
        }
        (TranscriptEcFormat::Json | TranscriptEcFormat::Text, Some(_)) => {
            bail!("--table is only valid with --format tsv")
        }
        _ => {}
    }
    if table == Some(TranscriptEcTable::Membership) && !emit_membership {
        bail!("--table membership requires --emit-membership");
    }
    Ok(())
}

fn transcript_ec_output_parent(path: &Path) -> &Path {
    match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    }
}

fn validate_transcript_ec_output_path(path: &Path) -> Result<()> {
    match std::fs::symlink_metadata(path) {
        Ok(_) => bail!("refusing to replace existing output {}", path.display()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error)
                .with_context(|| format!("inspecting output path {}", path.display()));
        }
    }
    let parent = transcript_ec_output_parent(path);
    let metadata = std::fs::metadata(parent)
        .with_context(|| format!("inspecting output directory {}", parent.display()))?;
    if !metadata.is_dir() {
        bail!("output parent is not a directory: {}", parent.display());
    }
    path.file_name().context("output path must name a file")?;
    Ok(())
}

fn write_transcript_ec_output(path: Option<&Path>, bytes: &[u8]) -> Result<()> {
    let Some(path) = path else {
        let stdout = std::io::stdout();
        let mut lock = stdout.lock();
        lock.write_all(bytes)
            .context("writing transcript-EC output")?;
        lock.flush().context("flushing transcript-EC output")?;
        return Ok(());
    };
    validate_transcript_ec_output_path(path)?;
    let outcome = publish_file_no_clobber(path, Durability::File, |writer| {
        writer.write_all(bytes)?;
        Ok(())
    })
    .with_context(|| format!("installing staged transcript-EC output {}", path.display()))?;
    for warning in outcome.warnings {
        eprintln!("warning: {warning}");
    }
    Ok(())
}

fn annotation_digest_bytes(identity: &anno::intent::AnnotationIdentity) -> Result<[u8; 32]> {
    let digest = identity
        .digest
        .as_deref()
        .context("bound annotation did not report its content digest")?;
    let encoded = digest
        .strip_prefix("blake3:")
        .context("bound annotation digest is not canonical BLAKE3")?;
    let hash = blake3::Hash::from_hex(encoded)
        .with_context(|| format!("bound annotation reported invalid digest {digest}"))?;
    Ok(*hash.as_bytes())
}

fn resolve_transcript_ec_selection(
    bound: &anno::intent::BoundAnnotation,
    feature: Option<&str>,
    locus: Option<&str>,
) -> Result<(
    BTreeSet<u32>,
    Vec<crate::archivecmd::transcriptec::TranscriptSelectorLocus>,
    TranscriptEcSelection,
)> {
    match (feature, locus) {
        (Some(requested), None) => {
            let resolver = anno::intent::IntentResolver::from_bound_annotation(bound)
                .context("building the bound annotation identifier resolver")?;
            let query = anno::intent::IdentifierQuery::parse(requested)
                .with_context(|| format!("invalid gene feature {requested:?}"))?;
            let resolved = resolver
                .resolve(&query)
                .with_context(|| format!("resolving gene feature {requested:?}"))?;
            if resolved.kind != anno::intent::FeatureKind::Gene {
                bail!(
                    "--feature must resolve to a gene, but {requested:?} resolved to {:?} {}",
                    resolved.kind,
                    resolved.stable_id
                );
            }
            let gene_matches: Vec<usize> = bound
                .annotation()
                .gene_ids
                .iter()
                .enumerate()
                .filter_map(|(index, id)| (id == &resolved.stable_id).then_some(index))
                .collect();
            if gene_matches.len() != 1 {
                bail!(
                    "resolved gene {} maps to {} annotation gene records",
                    resolved.stable_id,
                    gene_matches.len()
                );
            }
            let gene_index = u32::try_from(gene_matches[0])
                .context("resolved annotation gene index exceeds u32")?;
            let selected = bound
                .annotation()
                .transcripts
                .iter()
                .enumerate()
                .filter(|(_, transcript)| transcript.gene == gene_index)
                .map(|(index, _)| u32::try_from(index))
                .collect::<std::result::Result<BTreeSet<_>, _>>()
                .context("selected transcript index exceeds u32")?;
            let selector_loci = resolved
                .loci
                .iter()
                .map(|locus| crate::archivecmd::transcriptec::TranscriptSelectorLocus {
                    contig: locus.contig.clone(),
                    start: locus.start,
                    end: locus.end,
                })
                .collect();
            Ok((
                selected,
                selector_loci,
                TranscriptEcSelection::Gene {
                    requested: requested.to_owned(),
                    stable_id: resolved.stable_id,
                    display_name: resolved.display_name,
                    matched_by: resolved.matched_by,
                    loci: resolved.loci,
                },
            ))
        }
        (None, Some(requested)) => {
            let (contig, start, end) = parse_locus(requested)?;
            let chrom = bound
                .annotation()
                .chrom_ids
                .get(&contig)
                .copied()
                .with_context(|| format!("annotation has no contig {contig}"))?;
            let selected = bound
                .annotation()
                .overlapping(chrom, start, end)
                .into_iter()
                .collect();
            Ok((
                selected,
                vec![crate::archivecmd::transcriptec::TranscriptSelectorLocus {
                    contig: contig.clone(),
                    start,
                    end,
                }],
                TranscriptEcSelection::Locus {
                    requested: requested.to_owned(),
                    contig,
                    start,
                    end,
                },
            ))
        }
        (Some(_), Some(_)) => {
            bail!("transcript-ecs requires exactly one of --feature or --locus, not both")
        }
        (None, None) => {
            bail!("transcript-ecs requires exactly one of --feature or --locus")
        }
    }
}

fn checked_transcript_ec_add(target: &mut u64, value: u64, label: &str) -> Result<()> {
    *target = target
        .checked_add(value)
        .with_context(|| format!("{label} overflow"))?;
    Ok(())
}

fn transcript_ec_aggregation_name(aggregation: QueryAggregation) -> &'static str {
    match aggregation {
        QueryAggregation::Cell => "cell",
        QueryAggregation::Group => "group",
        QueryAggregation::Bulk => "bulk",
    }
}

fn transcript_ec_aggregation_target(
    scope: &QueryScope,
    cell_id: u32,
    barcode: &str,
) -> Result<(String, Option<u32>)> {
    match scope.aggregation {
        QueryAggregation::Cell => Ok((barcode.to_owned(), Some(cell_id))),
        QueryAggregation::Group => {
            let group = scope
                .group_of
                .get(&cell_id)
                .copied()
                .with_context(|| format!("scoped cell {cell_id} has no group assignment"))?;
            let name = scope
                .group_names
                .get(group as usize)
                .with_context(|| format!("cell {cell_id} references missing group {group}"))?;
            Ok((name.clone(), None))
        }
        QueryAggregation::Bulk => Ok(("all".to_owned(), None)),
    }
}

fn transcript_ec_hash_frame(hasher: &mut blake3::Hasher, value: &[u8]) -> Result<()> {
    let length = u64::try_from(value.len()).context("scope identity field exceeds u64")?;
    hasher.update(&length.to_le_bytes());
    hasher.update(value);
    Ok(())
}

fn transcript_ec_scope_json(
    scope: &QueryScope,
    cell_barcodes: &BTreeMap<u32, String>,
) -> Result<serde_json::Value> {
    let mut members = Vec::with_capacity(scope.selected_cells);
    for (&cell_id, barcode) in cell_barcodes {
        if !scope.includes(cell_id) {
            continue;
        }
        let group = if scope.source == "groups" {
            let group_id = scope.group_of.get(&cell_id).copied().with_context(|| {
                format!("selected cell {cell_id} has no resolved group assignment")
            })?;
            Some(
                scope
                    .group_names
                    .get(group_id as usize)
                    .with_context(|| format!("cell {cell_id} references missing group {group_id}"))?
                    .clone(),
            )
        } else {
            None
        };
        members.push((barcode.clone(), group));
    }
    members.sort();
    if members.len() != scope.selected_cells {
        bail!(
            "resolved cell scope has {} barcodes, expected {}",
            members.len(),
            scope.selected_cells
        );
    }

    let mut hasher = blake3::Hasher::new();
    hasher.update(b"gravlax.query.resolved-cell-scope.v1\0");
    transcript_ec_hash_frame(&mut hasher, scope.source.as_bytes())?;
    let count = u64::try_from(members.len()).context("resolved cell scope exceeds u64")?;
    hasher.update(&count.to_le_bytes());
    for (barcode, group) in &members {
        transcript_ec_hash_frame(&mut hasher, barcode.as_bytes())?;
        match group {
            Some(group) => {
                hasher.update(&[1]);
                transcript_ec_hash_frame(&mut hasher, group.as_bytes())?;
            }
            None => {
                hasher.update(&[0]);
            }
        }
    }

    let mut value = scope.json();
    let object = value
        .as_object_mut()
        .context("query scope JSON is not an object")?;
    object.insert(
        "resolved_population_blake3".into(),
        json!(format!("blake3:{}", hasher.finalize().to_hex())),
    );
    if let Some(path) = &scope.source_path {
        object.insert("source_path".into(), json!(path));
    }
    if let Some(digest) = &scope.source_content_blake3 {
        object.insert("source_content_blake3".into(), json!(digest));
    }
    Ok(value)
}

fn transcript_ec_catalog_table(rows: Vec<TranscriptEcCatalogRow>) -> Result<TypedTable> {
    let schema = TableSchema::new(
        TRANSCRIPT_EC_CATALOG_SCHEMA,
        vec![
            Field::new("ec_id", DataType::String),
            Field::new("transcript_ids", DataType::Json),
            Field::new("gene_ids", DataType::Json),
            Field::new("ambiguous", DataType::Boolean),
            Field::new("archived_umi_class_count", DataType::UInt64),
            Field::new("cell_count", DataType::UInt64),
            Field::new("complete_umi_class_count", DataType::UInt64),
        ],
    )?;
    let typed_rows = rows
        .into_iter()
        .map(|row| {
            vec![
                ScalarValue::String(row.ec_id),
                ScalarValue::Json(json!(row.transcript_ids)),
                ScalarValue::Json(json!(row.gene_ids)),
                ScalarValue::Boolean(row.ambiguous),
                ScalarValue::UInt64(row.archived_umi_class_count),
                ScalarValue::UInt64(row.cell_count),
                ScalarValue::UInt64(row.complete_umi_class_count),
            ]
        })
        .collect();
    Ok(TypedTable::new(schema, typed_rows)?)
}

fn transcript_ec_counts_table(rows: Vec<TranscriptEcCountRow>) -> Result<TypedTable> {
    let schema = TableSchema::new(
        TRANSCRIPT_EC_COUNTS_SCHEMA,
        vec![
            Field::new("aggregation", DataType::String),
            Field::new("key", DataType::String),
            Field::new("cell_id", DataType::UInt64).nullable(),
            Field::new("ec_id", DataType::String).nullable(),
            Field::new("archived_umi_class_count", DataType::UInt64),
            Field::new("ambiguous_umi_class_count", DataType::UInt64),
            Field::new("no_compatible_transcript_umi_class_count", DataType::UInt64),
            Field::new("conflicting_umi_class_count", DataType::UInt64),
            Field::new("complete_umi_class_count", DataType::UInt64),
            Field::new("incomplete_umi_class_count", DataType::UInt64),
            Field::new("retained_record_count", DataType::UInt64),
            Field::new("represented_alignment_count", DataType::UInt64),
        ],
    )?;
    let typed_rows = rows
        .into_iter()
        .map(|row| {
            vec![
                ScalarValue::String(row.aggregation.into()),
                ScalarValue::String(row.key),
                row.cell_id
                    .map_or(ScalarValue::Null, |value| ScalarValue::UInt64(value.into())),
                row.ec_id.map_or(ScalarValue::Null, ScalarValue::String),
                ScalarValue::UInt64(row.archived_umi_class_count),
                ScalarValue::UInt64(row.ambiguous_umi_class_count),
                ScalarValue::UInt64(row.no_compatible_transcript_umi_class_count),
                ScalarValue::UInt64(row.conflicting_umi_class_count),
                ScalarValue::UInt64(row.complete_umi_class_count),
                ScalarValue::UInt64(row.incomplete_umi_class_count),
                ScalarValue::UInt64(row.retained_record_count),
                ScalarValue::UInt64(row.represented_alignment_count),
            ]
        })
        .collect();
    Ok(TypedTable::new(schema, typed_rows)?)
}

fn transcript_ec_membership_table(rows: Vec<TranscriptEcMembershipRow>) -> Result<TypedTable> {
    let schema = TableSchema::new(
        TRANSCRIPT_EC_MEMBERSHIP_SCHEMA,
        vec![
            Field::new("umi_class", DataType::UInt64),
            Field::new("cell_id", DataType::UInt64),
            Field::new("barcode", DataType::String),
            Field::new("aggregation", DataType::String),
            Field::new("key", DataType::String),
            Field::new("ec_id", DataType::String).nullable(),
            Field::new("retained_record_count", DataType::UInt64),
            Field::new("represented_alignment_count", DataType::UInt64),
            Field::new("compatible_record_count", DataType::UInt64),
            Field::new("unmatched_record_count", DataType::UInt64),
            Field::new("ambiguous", DataType::Boolean),
            Field::new("no_compatible_transcript", DataType::Boolean),
            Field::new("conflict", DataType::Boolean),
            Field::new("complete_within_archive_quotient", DataType::Boolean),
            Field::new("retained_representatives_complete", DataType::Boolean),
        ],
    )?;
    let typed_rows = rows
        .into_iter()
        .map(|row| {
            vec![
                ScalarValue::UInt64(row.umi_class.into()),
                ScalarValue::UInt64(row.cell_id.into()),
                ScalarValue::String(row.barcode),
                ScalarValue::String(row.aggregation.into()),
                ScalarValue::String(row.key),
                row.ec_id.map_or(ScalarValue::Null, ScalarValue::String),
                ScalarValue::UInt64(row.retained_record_count),
                ScalarValue::UInt64(row.represented_alignment_count),
                ScalarValue::UInt64(row.compatible_record_count),
                ScalarValue::UInt64(row.unmatched_record_count),
                ScalarValue::Boolean(row.ambiguous),
                ScalarValue::Boolean(row.no_compatible_transcript),
                ScalarValue::Boolean(row.conflict),
                ScalarValue::Boolean(row.complete_within_archive_quotient),
                ScalarValue::Boolean(row.retained_representatives_complete),
            ]
        })
        .collect();
    Ok(TypedTable::new(schema, typed_rows)?)
}

#[expect(
    clippy::too_many_arguments,
    reason = "the envelope constructor keeps archive provenance, annotation identity, scope, and output limits explicit"
)]
fn build_transcript_ec_envelope(
    report: TranscriptEquivalenceReport,
    archive_path: &Path,
    annotation_identity: &anno::intent::AnnotationIdentity,
    selection: TranscriptEcSelection,
    scope: &QueryScope,
    emit_membership: bool,
    max_ecs: usize,
    max_memberships: usize,
) -> Result<(TranscriptEcResultData, ResultContext)> {
    let TranscriptEquivalenceReport {
        archive_version,
        archive_root_blake3,
        strand_policy,
        scope: derivation_scope,
        semantics: compatibility_semantics,
        transcript_sets,
        classes,
        cells,
        totals,
        ..
    } = report;
    let selected_transcript_ids = derivation_scope
        .selected_transcript_ids
        .context("transcript-EC query did not retain its selected transcript universe")?;
    let cell_barcodes: BTreeMap<u32, String> = cells
        .into_iter()
        .map(|cell| (cell.cell_id, cell.barcode))
        .collect();
    let set_metadata: BTreeMap<String, crate::archivecmd::transcriptec::TranscriptSetSummary> =
        transcript_sets
            .into_iter()
            .map(|set| (set.ec_id.clone(), set))
            .collect();
    let selector_relevant_umi_classes = classes
        .iter()
        .filter(|class| class.selector_relevant)
        .count();
    let scoped_classes: Vec<_> = classes
        .iter()
        .filter(|class| class.selector_relevant && scope.includes(class.cell_id))
        .collect();
    if emit_membership && scoped_classes.len() > max_memberships {
        bail!(
            "transcript-EC membership has {} rows, exceeding --max-memberships {max_memberships}; output was not truncated",
            scoped_classes.len()
        );
    }

    let aggregation = transcript_ec_aggregation_name(scope.aggregation);
    let mut catalog_accumulators: BTreeMap<String, TranscriptEcCatalogAccumulator> =
        BTreeMap::new();
    let mut count_accumulators: BTreeMap<TranscriptEcCountKey, TranscriptEcCountAccumulator> =
        BTreeMap::new();
    let mut memberships = Vec::with_capacity(if emit_membership {
        scoped_classes.len()
    } else {
        0
    });
    let mut scoped_cells = BTreeSet::new();
    let mut assigned = 0u64;
    let mut ambiguous = 0u64;
    let mut no_compatible = 0u64;
    let mut conflicting = 0u64;
    let mut complete = 0u64;

    for class in scoped_classes {
        let barcode = cell_barcodes.get(&class.cell_id).with_context(|| {
            format!(
                "transcript-EC class {} references absent cell {}",
                class.umi_class, class.cell_id
            )
        })?;
        scoped_cells.insert(class.cell_id);
        let (key, cell_id) = transcript_ec_aggregation_target(scope, class.cell_id, barcode)?;
        if let Some(ec_id) = &class.ec_id {
            let catalog = catalog_accumulators.entry(ec_id.clone()).or_default();
            checked_transcript_ec_add(&mut catalog.class_count, 1, "catalog class count")?;
            checked_transcript_ec_add(
                &mut catalog.complete_class_count,
                u64::from(class.complete_within_archive_quotient),
                "catalog complete-class count",
            )?;
            catalog.cells.insert(class.cell_id);
            checked_transcript_ec_add(&mut assigned, 1, "assigned class count")?;
        }
        checked_transcript_ec_add(
            &mut ambiguous,
            u64::from(class.ambiguous),
            "ambiguous class count",
        )?;
        checked_transcript_ec_add(
            &mut no_compatible,
            u64::from(class.no_compatible_transcript),
            "no-compatible class count",
        )?;
        checked_transcript_ec_add(
            &mut conflicting,
            u64::from(class.conflict),
            "conflicting class count",
        )?;
        checked_transcript_ec_add(
            &mut complete,
            u64::from(class.complete_within_archive_quotient),
            "complete class count",
        )?;

        let count_key = TranscriptEcCountKey {
            key: key.clone(),
            cell_id,
            ec_id: class.ec_id.clone(),
        };
        let counts = count_accumulators.entry(count_key).or_default();
        checked_transcript_ec_add(&mut counts.class_count, 1, "count-row class count")?;
        checked_transcript_ec_add(
            &mut counts.ambiguous_class_count,
            u64::from(class.ambiguous),
            "count-row ambiguous class count",
        )?;
        checked_transcript_ec_add(
            &mut counts.no_compatible_class_count,
            u64::from(class.no_compatible_transcript),
            "count-row no-compatible class count",
        )?;
        checked_transcript_ec_add(
            &mut counts.conflicting_class_count,
            u64::from(class.conflict),
            "count-row conflicting class count",
        )?;
        checked_transcript_ec_add(
            &mut counts.complete_class_count,
            u64::from(class.complete_within_archive_quotient),
            "count-row complete class count",
        )?;
        checked_transcript_ec_add(
            &mut counts.incomplete_class_count,
            u64::from(!class.complete_within_archive_quotient),
            "count-row incomplete class count",
        )?;
        checked_transcript_ec_add(
            &mut counts.retained_record_count,
            class.retained_record_count,
            "count-row retained-record count",
        )?;
        checked_transcript_ec_add(
            &mut counts.represented_alignment_count,
            class.represented_alignment_count,
            "count-row represented-alignment count",
        )?;

        if emit_membership {
            memberships.push(TranscriptEcMembershipRow {
                umi_class: class.umi_class,
                cell_id: class.cell_id,
                barcode: barcode.clone(),
                aggregation,
                key,
                ec_id: class.ec_id.clone(),
                retained_record_count: class.retained_record_count,
                represented_alignment_count: class.represented_alignment_count,
                compatible_record_count: class.compatible_record_count,
                unmatched_record_count: class.unmatched_record_count,
                ambiguous: class.ambiguous,
                no_compatible_transcript: class.no_compatible_transcript,
                conflict: class.conflict,
                complete_within_archive_quotient: class.complete_within_archive_quotient,
                retained_representatives_complete: class.retained_representatives_complete,
            });
        }
    }
    memberships.sort_by_key(|row| row.umi_class);

    if catalog_accumulators.len() > max_ecs {
        bail!(
            "transcript-EC catalog has {} rows, exceeding --max-ecs {max_ecs}; output was not truncated",
            catalog_accumulators.len()
        );
    }
    if count_accumulators.len() > TRANSCRIPT_EC_MAX_COUNT_ROWS {
        bail!(
            "transcript-EC counts table has {} rows, exceeding the hard safety limit {}; output was not truncated",
            count_accumulators.len(),
            TRANSCRIPT_EC_MAX_COUNT_ROWS
        );
    }
    let mut catalog_rows = Vec::with_capacity(catalog_accumulators.len());
    for (ec_id, aggregate) in catalog_accumulators {
        let metadata = set_metadata
            .get(&ec_id)
            .with_context(|| format!("transcript-EC {ec_id} is absent from the core catalog"))?;
        catalog_rows.push(TranscriptEcCatalogRow {
            ec_id,
            transcript_ids: metadata.transcript_ids.clone(),
            gene_ids: metadata.gene_ids.clone(),
            ambiguous: metadata.ambiguous,
            archived_umi_class_count: aggregate.class_count,
            cell_count: u64::try_from(aggregate.cells.len())
                .context("catalog cell count exceeds u64")?,
            complete_umi_class_count: aggregate.complete_class_count,
        });
    }

    let count_rows: Vec<_> = count_accumulators
        .into_iter()
        .map(|(key, counts)| TranscriptEcCountRow {
            aggregation,
            key: key.key,
            cell_id: key.cell_id,
            ec_id: key.ec_id,
            archived_umi_class_count: counts.class_count,
            ambiguous_umi_class_count: counts.ambiguous_class_count,
            no_compatible_transcript_umi_class_count: counts.no_compatible_class_count,
            conflicting_umi_class_count: counts.conflicting_class_count,
            complete_umi_class_count: counts.complete_class_count,
            incomplete_umi_class_count: counts.incomplete_class_count,
            retained_record_count: counts.retained_record_count,
            represented_alignment_count: counts.represented_alignment_count,
        })
        .collect();
    let scoped_umi_classes = count_rows
        .iter()
        .try_fold(0u64, |sum, row| {
            sum.checked_add(row.archived_umi_class_count)
        })
        .context("scoped class count overflow")?;
    let incomplete = scoped_umi_classes
        .checked_sub(complete)
        .context("complete class count exceeds scoped class count")?;
    let unassigned = scoped_umi_classes
        .checked_sub(assigned)
        .context("assigned class count exceeds scoped class count")?;
    let membership_rows =
        u64::try_from(memberships.len()).context("membership row count exceeds u64")?;
    let summary = TranscriptEcQuerySummary {
        archive_version,
        archive_root_blake3: archive_root_blake3.clone(),
        strand_policy,
        archive_umi_classes_scanned: totals.archive_classes_scanned,
        selector_relevant_umi_classes: u64::try_from(selector_relevant_umi_classes)
            .context("selector-relevant class count exceeds u64")?,
        scoped_umi_classes,
        scoped_cells_with_classes: u64::try_from(scoped_cells.len())
            .context("scoped cell count exceeds u64")?,
        transcript_ecs: u64::try_from(catalog_rows.len())
            .context("transcript-EC count exceeds u64")?,
        count_rows: u64::try_from(count_rows.len()).context("count row count exceeds u64")?,
        membership_rows,
        assigned_umi_classes: assigned,
        unassigned_umi_classes: unassigned,
        ambiguous_umi_classes: ambiguous,
        no_compatible_transcript_umi_classes: no_compatible,
        conflicting_umi_classes: conflicting,
        complete_umi_classes: complete,
        incomplete_umi_classes: incomplete,
    };
    if summary.assigned_umi_classes + summary.unassigned_umi_classes != summary.scoped_umi_classes
        || summary.complete_umi_classes + summary.incomplete_umi_classes
            != summary.scoped_umi_classes
    {
        bail!("scoped transcript-EC conservation invariant failed");
    }

    let annotation_digest = annotation_identity
        .digest
        .clone()
        .context("bound annotation identity has no digest")?;
    let cell_scope = transcript_ec_scope_json(scope, &cell_barcodes)?;
    let selection_value = serde_json::to_value(&selection)?;
    let selected_transcript_count = selected_transcript_ids.len();
    let data = TranscriptEcResultData {
        scope: TranscriptEcQueryScope {
            selection,
            selected_transcript_count,
            selected_transcript_ids,
            cells: cell_scope.clone(),
            archive_access: "full_archive_scan",
            chunk_pruning_applied: false,
        },
        semantics: TranscriptEcQuerySemantics {
            compatibility: compatibility_semantics,
            count_unit: "archived UMI class",
            aggregation: "cell scope and aggregation are applied after per-cell transcript-EC assignment; classes are never collapsed across cells",
            inference: "counts are archived UMI-class counts, not gene-collapse counts or transcript abundance estimates",
            phasing: "compatibility does not claim full-isoform phasing",
        },
        summary,
        catalog: transcript_ec_catalog_table(catalog_rows)?,
        counts: transcript_ec_counts_table(count_rows)?,
        membership: if emit_membership {
            Some(transcript_ec_membership_table(memberships)?)
        } else {
            None
        },
    };
    let mut parameters = BTreeMap::new();
    parameters.insert("selection".into(), selection_value);
    parameters.insert(
        "selected_transcript_count".into(),
        json!(selected_transcript_count),
    );
    parameters.insert("cell_scope".into(), cell_scope);
    parameters.insert("archive_access".into(), json!("full_archive_scan"));
    parameters.insert("chunk_pruning_applied".into(), json!(false));
    parameters.insert("emit_membership".into(), json!(emit_membership));
    parameters.insert("max_ecs".into(), json!(max_ecs));
    parameters.insert("max_memberships".into(), json!(max_memberships));
    parameters.insert(
        "max_count_rows".into(),
        json!(TRANSCRIPT_EC_MAX_COUNT_ROWS),
    );
    parameters.insert("archive_path".into(), json!(archive_path));
    parameters.insert("strand_policy".into(), serde_json::to_value(strand_policy)?);
    let mut warnings = Vec::new();
    let archives = match archive_root_blake3 {
        Some(root) => vec![format!("aie-directory-root-v2:{root}")],
        None => {
            warnings.push(
                "legacy v1 archive has no rooted content commitment; its path locator is not a portable content identity"
                    .into(),
            );
            vec![format!(
                "aie-v{archive_version}-unrooted:{}",
                archive_path.display()
            )]
        }
    };
    let context = ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            archives,
            assembly: Some(annotation_identity.assembly.clone()),
            annotation: Some(annotation_identity.annotation.clone()),
            annotation_digest: Some(annotation_digest),
            parameters,
            ..Default::default()
        },
        warnings,
    };
    Ok((data, context))
}

fn selected_transcript_ec_table(
    data: &TranscriptEcResultData,
    table: TranscriptEcTable,
) -> Result<&TypedTable> {
    match table {
        TranscriptEcTable::Catalog => Ok(&data.catalog),
        TranscriptEcTable::Counts => Ok(&data.counts),
        TranscriptEcTable::Membership => data
            .membership
            .as_ref()
            .context("membership table was not requested"),
    }
}

fn render_transcript_ec_output(
    data: TranscriptEcResultData,
    context: &ResultContext,
    format: TranscriptEcFormat,
    table: Option<TranscriptEcTable>,
) -> Result<Vec<u8>> {
    match format {
        TranscriptEcFormat::Json => {
            let envelope = ResultEnvelope::new(TRANSCRIPT_EC_QUERY_SCHEMA, context.clone(), data)?;
            let mut output = serde_json::to_vec(&envelope)
                .context("serializing transcript-EC result envelope")?;
            output.push(b'\n');
            Ok(output)
        }
        TranscriptEcFormat::Tsv => {
            let mut output = Vec::new();
            let selected = selected_transcript_ec_table(
                &data,
                table.context("--format tsv requires --table")?,
            )?;
            write_table(
                &mut output,
                &selected.schema,
                selected.rows.clone(),
                OutputFormat::Tsv,
                context,
            )?;
            Ok(output)
        }
        TranscriptEcFormat::Text => {
            let mut output = Vec::new();
            writeln!(output, "transcript EC query ({TRANSCRIPT_EC_QUERY_SCHEMA})")?;
            writeln!(
                output,
                "archive: v{} root {}",
                data.summary.archive_version,
                data.summary
                    .archive_root_blake3
                    .as_deref()
                    .unwrap_or("unavailable (legacy archive)")
            )?;
            writeln!(
                output,
                "annotation: {} / {} / {}",
                context.provenance.assembly.as_deref().unwrap_or("unknown"),
                context
                    .provenance
                    .annotation
                    .as_deref()
                    .unwrap_or("unknown"),
                context
                    .provenance
                    .annotation_digest
                    .as_deref()
                    .unwrap_or("unknown")
            )?;
            writeln!(
                output,
                "scope: {}; {} selected transcripts; full archive scan",
                serde_json::to_string(&data.scope.selection)?,
                data.scope.selected_transcript_count
            )?;
            writeln!(
                output,
                "result: {} scoped archived UMI classes, {} compatible transcript sets, {} scoped cells",
                data.summary.scoped_umi_classes,
                data.summary.transcript_ecs,
                data.summary.scoped_cells_with_classes
            )?;
            writeln!(output)?;
            write_table(
                &mut output,
                &data.catalog.schema,
                data.catalog.rows.clone(),
                OutputFormat::Text,
                context,
            )?;
            writeln!(output)?;
            write_table(
                &mut output,
                &data.counts.schema,
                data.counts.rows.clone(),
                OutputFormat::Text,
                context,
            )?;
            if let Some(membership) = &data.membership {
                writeln!(output)?;
                write_table(
                    &mut output,
                    &membership.schema,
                    membership.rows.clone(),
                    OutputFormat::Text,
                    context,
                )?;
            }
            Ok(output)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_transcript_ecs(
    archive: &Path,
    annotation_file: &Path,
    assembly: String,
    annotation_label: String,
    annotation_digest: Option<String>,
    feature: Option<String>,
    locus: Option<String>,
    solo_strand: crate::archivecmd::SoloStrandArg,
    scope_args: QueryScopeArgs,
    emit_membership: bool,
    max_ecs: usize,
    max_memberships: usize,
    format: TranscriptEcFormat,
    table: Option<TranscriptEcTable>,
    output: Option<PathBuf>,
) -> Result<()> {
    validate_transcript_ec_output_args(format, table, emit_membership)?;
    if max_ecs == 0 {
        bail!("--max-ecs must be at least 1");
    }
    if max_memberships == 0 {
        bail!("--max-memberships must be at least 1");
    }
    if let Some(path) = output.as_deref() {
        validate_transcript_ec_output_path(path)?;
    }
    let mut identity = anno::intent::AnnotationIdentity::new(assembly, annotation_label)
        .context("validating annotation identity")?;
    if let Some(digest) = annotation_digest {
        identity = identity
            .with_digest(digest)
            .context("validating expected annotation digest")?;
    }
    let bound = anno::intent::BoundAnnotation::from_path(annotation_file, identity)
        .with_context(|| format!("loading bound annotation {}", annotation_file.display()))?;
    let (selected_transcripts, selector_loci, selection) =
        resolve_transcript_ec_selection(&bound, feature.as_deref(), locus.as_deref())?;
    let options = TranscriptEquivalenceOptions::new(
        solo_strand.into(),
        annotation_digest_bytes(bound.identity())?,
    )
    .with_transcript_universe(selected_transcripts)
    .with_selector_loci(selector_loci);
    // Keep all archive classes in the core report so scope files can name any archive cell.
    // The query result below retains classes with any retained placement block overlapping the
    // selector, independently of compatibility. An overlapping discordant class therefore stays
    // visible as an explicit no-compatible-transcript outcome.
    let report = derive_transcript_equivalence_classes_with_annotation(
        archive,
        bound.annotation(),
        &options,
    )?;
    let scope = load_transcript_ec_scope(&report.cells, &scope_args)?;
    let (data, context) = build_transcript_ec_envelope(
        report,
        archive,
        bound.identity(),
        selection,
        &scope,
        emit_membership,
        max_ecs,
        max_memberships,
    )?;
    let bytes = render_transcript_ec_output(data, &context, format, table)?;
    write_transcript_ec_output(output.as_deref(), &bytes)
}

fn validate_query_args(what: &What) -> Result<()> {
    if let What::Apa {
        site_gap,
        groups,
        genome,
        drop_ip,
        permute,
        ..
    } = what
    {
        if *site_gap == 0 {
            bail!("--site-gap must be at least 1");
        }
        if *drop_ip && genome.is_none() {
            bail!("--drop-ip requires --genome");
        }
        if permute.is_some() && groups.is_none() {
            bail!("--permute requires --groups");
        }
    }
    if let What::ApaTest { site_gap, .. } = what {
        if *site_gap == 0 {
            bail!("--site-gap must be at least 1");
        }
    }
    match what {
        What::Batch { uniform_output, .. } => {
            validate_uniform_output_flags(uniform_output, false, false)?;
        }
        What::Region {
            plot,
            export_prefix,
            tsv,
            json,
            uniform_output,
            ..
        } => {
            validate_uniform_output_flags(uniform_output, *tsv, *json)?;
            if uniform_output.format.is_some() && (plot.is_some() || export_prefix.is_some()) {
                bail!(
                    "uniform region output cannot be combined with --plot or --export-prefix; run side-artifact export separately"
                );
            }
        }
        What::Junction {
            tsv,
            json,
            uniform_output,
            ..
        }
        | What::Junctions {
            tsv,
            json,
            uniform_output,
            ..
        }
        | What::Jset {
            tsv,
            json,
            uniform_output,
            ..
        }
        | What::Events {
            tsv,
            json,
            uniform_output,
            ..
        } => validate_uniform_output_flags(uniform_output, *tsv, *json)?,
        What::SpliceGraph {
            json,
            uniform_output,
            ..
        } => validate_uniform_output_flags(uniform_output, false, *json)?,
        What::Apa {
            tsv,
            plot,
            uniform_output,
            ..
        } => {
            validate_uniform_output_flags(uniform_output, *tsv, false)?;
            if uniform_output.format.is_some() && plot.is_some() {
                bail!(
                    "uniform APA output cannot be combined with --plot; run side-artifact export separately"
                );
            }
        }
        What::ApaTest { uniform_output, .. } => {
            validate_uniform_output_flags(uniform_output, false, false)?;
        }
        What::Discover {
            tsv,
            emit_gtf,
            uniform_output,
            ..
        } => {
            validate_uniform_output_flags(uniform_output, *tsv, false)?;
            if uniform_output.format.is_some() && emit_gtf.is_some() {
                bail!(
                    "uniform discovery output cannot be combined with --emit-gtf; run side-artifact export separately"
                );
            }
        }
        What::TranscriptEcs { .. } => {}
    }
    Ok(())
}

fn validate_uniform_output_flags(
    output: &UniformQueryOutputArgs,
    legacy_tsv: bool,
    legacy_json: bool,
) -> Result<()> {
    if output.format.is_some() && (legacy_tsv || legacy_json) {
        bail!("--format cannot be combined with legacy --tsv or --json");
    }
    preflight_uniform_output(output)
}

/// Avoid an expensive archive query when the destination is already occupied. This check is only
/// an early UX guard: the publisher's no-clobber hard-link remains the authoritative race check.
fn preflight_uniform_output(output: &UniformQueryOutputArgs) -> Result<()> {
    let Some(path) = output.output.as_deref() else {
        return Ok(());
    };
    if path.file_name().is_none() {
        bail!("uniform output path must name a file: {}", path.display());
    }
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| std::path::Path::new("."));
    let parent_metadata = std::fs::metadata(parent).with_context(|| {
        format!(
            "checking parent directory {} for uniform output {}",
            parent.display(),
            path.display()
        )
    })?;
    if !parent_metadata.is_dir() {
        bail!(
            "uniform output parent is not a directory: {}",
            parent.display()
        );
    }
    match std::fs::symlink_metadata(path) {
        Ok(_) => bail!("refusing to replace existing output {}", path.display()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error)
            .with_context(|| format!("checking uniform output destination {}", path.display())),
    }
}

fn write_uniform_bundle_output<F>(output: &UniformQueryOutputArgs, render: F) -> Result<()>
where
    F: Fn(&mut dyn Write, UniformQueryFormat) -> std::result::Result<(), OutputError>,
{
    let format = output
        .format
        .context("uniform output requires an explicit --format")?;
    if let Some(path) = output.output.as_deref() {
        let outcome =
            publish_file_no_clobber(path, Durability::Flush, |writer| render(writer, format))?;
        for warning in outcome.warnings {
            eprintln!("warning: {warning}");
        }
    } else {
        let stdout = std::io::stdout();
        let mut writer = std::io::BufWriter::new(stdout.lock());
        render(&mut writer, format)?;
    }
    Ok(())
}

fn uniform_group_scope_provenance(
    path: &Path,
    source_content_blake3: &str,
    group_names: &[String],
    cell_group: &FxHashMap<u32, u32>,
    archive_cells: usize,
) -> Result<serde_json::Value> {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"gravlax-query-resolved-groups-v1\0");
    let mut update = |cell: u32, group: u32| -> Result<()> {
        let name = group_names
            .get(group as usize)
            .context("resolved group index is out of range")?;
        hasher.update(&cell.to_le_bytes());
        hasher.update(&(name.len() as u64).to_le_bytes());
        hasher.update(name.as_bytes());
        Ok(())
    };
    if cell_group.len().saturating_mul(4) < archive_cells {
        let mut cells: Vec<u32> = cell_group.keys().copied().collect();
        cells.sort_unstable();
        for cell in cells {
            update(cell, cell_group[&cell])?;
        }
    } else {
        for cell in 0..u32::try_from(archive_cells)
            .context("archive cell dictionary exceeds the u32 cell-id domain")?
        {
            if let Some(&group) = cell_group.get(&cell) {
                update(cell, group)?;
            }
        }
    }
    let mut selected_per_group = vec![0_u64; group_names.len()];
    for &group in cell_group.values() {
        selected_per_group[group as usize] += 1;
    }
    Ok(json!({
        "source": "groups",
        "source_path": path,
        "source_content_blake3": source_content_blake3,
        "resolved_mapping_blake3": format!("blake3:{}", hasher.finalize().to_hex()),
        "archive_cells": archive_cells,
        "selected_cells": cell_group.len(),
        "groups": group_names.iter().enumerate().map(|(index, name)| json!({
            "name": name,
            "selected_cells": selected_per_group[index],
        })).collect::<Vec<_>>(),
    }))
}

pub fn run(args: Args) -> Result<()> {
    validate_query_args(&args.what)?;
    let Args { archive, what } = args;
    let what = match what {
        What::TranscriptEcs {
            annotation_file,
            assembly,
            annotation_label,
            annotation_digest,
            feature,
            locus,
            solo_strand,
            scope,
            emit_membership,
            max_ecs,
            max_memberships,
            format,
            table,
            output,
        } => {
            return run_transcript_ecs(
                &archive,
                &annotation_file,
                assembly,
                annotation_label,
                annotation_digest,
                feature,
                locus,
                solo_strand,
                scope,
                emit_membership,
                max_ecs,
                max_memberships,
                format,
                table,
                output,
            );
        }
        other => other,
    };
    let t0 = std::time::Instant::now();
    let mut la = LazyArchive::open(&archive)?;
    let chunks = read_chunk_index(la.reader())?;
    let chrom_names = la.chrom_names.clone();
    let chrom_id = |name: &str| -> Result<u32> {
        chrom_names
            .iter()
            .position(|n| n == name)
            .map(|i| i as u32)
            .with_context(|| format!("unknown chromosome {name}"))
    };
    let t_open = t0.elapsed().as_secs_f32();

    match what {
        What::TranscriptEcs { .. } => unreachable!("transcript-EC query was dispatched above"),
        What::Batch {
            plan,
            top,
            uniform_output,
            scope,
        } => {
            run_batch(BatchRun {
                archive: &archive,
                la: &mut la,
                chunks: &chunks,
                chrom_names: &chrom_names,
                plan_path: &plan,
                top,
                uniform_output: &uniform_output,
                scope_args: &scope,
                t0,
                t_open,
            })?;
        }
        What::Jset {
            include,
            exclude,
            top,
            tsv,
            json,
            uniform_output,
            scope,
        } => {
            run_junction_set(
                &archive,
                &mut la,
                &chunks,
                &chrom_names,
                &include,
                &exclude,
                top,
                tsv,
                json,
                &uniform_output,
                &scope,
                t0,
                t_open,
            )?;
        }
        What::Events {
            locus,
            event_types,
            min_support,
            min_informative,
            max_events,
            top,
            gtf,
            tsv,
            json,
            uniform_output,
            scope,
        } => {
            run_events(
                &archive,
                &mut la,
                &chunks,
                &chrom_names,
                &locus,
                &event_types,
                min_support,
                min_informative,
                max_events,
                top,
                gtf.as_deref(),
                tsv,
                json,
                &uniform_output,
                &scope,
                t0,
                t_open,
            )?;
        }
        What::SpliceGraph {
            locus,
            min_support,
            min_path_umis,
            max_paths,
            json,
            uniform_output,
            scope,
        } => {
            run_splice_graph(
                &archive,
                &mut la,
                &chunks,
                &chrom_names,
                &locus,
                min_support,
                min_path_umis,
                max_paths,
                json,
                &uniform_output,
                &scope,
                t0,
                t_open,
            )?;
        }
        What::Region {
            locus,
            top,
            plot,
            export_prefix,
            gtf,
            tsv,
            json: json_output,
            uniform_output,
            scope: scope_args,
        } => {
            let mut scope = load_query_scope(&mut la, &scope_args)?;
            if uniform_output.format.is_some() {
                scope.ensure_resolved_mapping_digest()?;
            }
            let (chrom, start, end) = parse_locus(&locus)?;
            let cid = chrom_id(&chrom)?;
            let tracking = plot.is_some() || export_prefix.is_some();
            if tracking && end.saturating_sub(start) > 5_000_000 {
                bail!("--plot/--export-prefix support windows up to 5 Mb");
            }
            let shapes = if tracking { Some(la.shapes()?) } else { None };
            let wlen = (end - start) as usize;
            let mut covdiff = [vec![0i64; wlen + 1], vec![0i64; wlen + 1]];
            let mut jcount: [FxHashMap<(u32, u32), u64>; 2] =
                [FxHashMap::default(), FxHashMap::default()];
            // Range index: chunks are (chrom, bin) sorted; select the touched ones. The chunk
            // index's max_anchor answers "is this chunk entirely left of the window" without a
            // decode — the exact test the old code ran on the decoded molecules (its last anchor
            // IS max_anchor), so the selected set and the reported chunk count are unchanged.
            let selected: Vec<usize> = chunks
                .iter()
                .enumerate()
                .filter(|(_, c)| region_selects_chunk(c, cid, start, end))
                .map(|(i, _)| i)
                .collect();
            // Fused parallel pread + decompress + decode per selected chunk.
            let decoded: Vec<Vec<crate::rows::MolRec>> = {
                let (reader, tables) = la.reader_and_tables();
                let reader = &*reader;
                selected
                    .par_iter()
                    .map(|&i| {
                        let (c, raw_len) = reader.read_compressed_at(&format!("c{i}"))?;
                        let raw = evidence_io::format::decompress(&c, raw_len)?;
                        decode_chunk(&raw, &chunks[i], None, tables)
                    })
                    .collect::<Result<_>>()?
            };
            la.prefetch_coc(
                decoded
                    .iter()
                    .flatten()
                    .filter(|m| {
                        let anchor_hit = m.anchor() >= start && m.anchor() < end;
                        anchor_hit || (tracking && scope.active)
                    })
                    .map(|m| m.umi_class),
            )?;
            let mut per_cell: FxHashMap<u32, FxHashSet<u32>> = FxHashMap::default();
            let mut n_mols = 0u64;
            let n_chunks = selected.len() as u64;
            for mols in &decoded {
                for m in mols {
                    let a = m.anchor();
                    let anchor_hit = a >= start && a < end;
                    let cell = if anchor_hit || (tracking && scope.active) {
                        Some(la.cell_of(m.umi_class)?)
                    } else {
                        None
                    };
                    let in_scope = !scope.active || cell.is_some_and(|cell| scope.includes(cell));
                    if anchor_hit && in_scope {
                        n_mols += 1;
                        per_cell
                            .entry(cell.unwrap())
                            .or_default()
                            .insert(m.umi_class);
                    }
                    if tracking && in_scope {
                        let shapes = shapes.as_ref().unwrap();
                        let si = m.strand_rev as usize;
                        // Molecule coverage: union of the unique reps' aligned blocks.
                        let mut ivs: Vec<(u32, u32)> = Vec::new();
                        for ch in &m.chains {
                            for (pos, sh) in &ch.reps {
                                for (off, len) in &shapes[*sh as usize].blocks {
                                    ivs.push((pos + off, pos + off + len));
                                }
                                for w2 in shapes[*sh as usize].blocks.windows(2) {
                                    let dn = pos + w2[0].0 + w2[0].1;
                                    let ac = pos + w2[1].0;
                                    if dn >= start && ac < end {
                                        jcount[si].entry((dn, ac)).or_insert(0);
                                    }
                                }
                            }
                        }
                        ivs.sort_unstable();
                        let mut merged: Vec<(u32, u32)> = Vec::new();
                        for (s2, e2) in ivs {
                            match merged.last_mut() {
                                Some((_, le)) if s2 <= *le => *le = (*le).max(e2),
                                _ => merged.push((s2, e2)),
                            }
                        }
                        for (s2, e2) in merged {
                            let (cs, ce) = (s2.max(start), e2.min(end));
                            if cs < ce {
                                covdiff[si][(cs - start) as usize] += 1;
                                covdiff[si][(ce - start) as usize] -= 1;
                            }
                        }
                        // Per-molecule junction dedup: the entry() above created keys; count once.
                        let mut seen: FxHashSet<(u32, u32)> = FxHashSet::default();
                        for ch in &m.chains {
                            for (pos, sh) in &ch.reps {
                                for w2 in shapes[*sh as usize].blocks.windows(2) {
                                    let dn = pos + w2[0].0 + w2[0].1;
                                    let ac = pos + w2[1].0;
                                    if dn >= start && ac < end && seen.insert((dn, ac)) {
                                        *jcount[si].get_mut(&(dn, ac)).unwrap() += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if tracking {
                let cov: Vec<Vec<u32>> = covdiff
                    .iter()
                    .map(|d| {
                        let mut run = 0i64;
                        d[..wlen].iter().map(|v| { run += v; run.max(0) as u32 }).collect()
                    })
                    .collect();
                let juncs: Vec<Vec<(u32, u32, u64)>> = jcount
                    .iter()
                    .map(|m2| m2.iter().map(|(&(d2, a2), &n2)| (d2, a2, n2)).collect())
                    .collect();
                if let Some(prefix) = &export_prefix {
                    export_igv(prefix, &chrom, start, &cov, &juncs)?;
                }
                if let Some(plot_path) = &plot {
                    let genes = match &gtf {
                        Some(g) => gene_underlay(g, &chrom, start, end)?,
                        None => Vec::new(),
                    };
                    crate::plots::region_plot(
                        plot_path, &chrom, start, end,
                        [&cov[0], &cov[1]],
                        [&juncs[0], &juncs[1]],
                        &genes,
                        &format!("molecule evidence, {locus}"),
                    )?;
                    eprintln!("wrote {}", plot_path.display());
                }
            }
            let archive_identity = uniform_output.format.map(|_| {
                (
                    la.reader().archive_version(),
                    la.reader()
                        .content_commitment()
                        .map(|commitment| commitment.to_hex()),
                )
            });
            let dict = la.cells()?;
            let counts = scoped_counts_with_cell_order(
                &per_cell,
                &scope,
                uniform_output.format.map(|_| dict),
            );
            // Historical unscoped `region --top 0` printed no cell rows. Preserve that
            // interface exactly; scoped/tabular paths use the shared 0 = all convention.
            let legacy_limit = if top == 0 && scope.active {
                counts.cells.len()
            } else {
                top.min(counts.cells.len())
            };
            let (_, uniform_limit) = uniform_row_limit(&counts, &scope, top);
            if let Some(format) = uniform_output.format {
                let summary = RegionUniformSummary {
                    coordinates: "0-based half-open",
                    anchor_semantics: true,
                    chrom: &chrom,
                    start,
                    end,
                    molecules: n_mols,
                    umis: counts.total_umis as u64,
                    cells: counts.cells.len() as u64,
                    chunks_decoded: n_chunks,
                };
                let parameters =
                    uniform_query_parameters(&scope, top, "range-index-selected archive chunks")?;
                let (archive_version, archive_root) =
                    archive_identity.expect("uniform output captured archive identity");
                let context =
                    uniform_query_context(&archive, archive_version, archive_root, parameters);
                write_uniform_count_output(
                    &uniform_output,
                    REGION_UNIFORM_RESULT_SCHEMA,
                    REGION_UNIFORM_COUNTS_SCHEMA,
                    &context,
                    &summary,
                    &counts,
                    &scope,
                    dict,
                    uniform_limit,
                )?;
                eprintln!(
                    "region {locus}: {n_mols} molecules, {} UMIs across {} cells ({} chunks; uniform {format:?}; open {t_open:.2}s, total {:.2}s)",
                    counts.total_umis,
                    counts.cells.len(),
                    n_chunks,
                    t0.elapsed().as_secs_f32()
                );
            } else if tsv {
                match scope.aggregation {
                    QueryAggregation::Cell => {
                        println!("barcode\tumis");
                        for (cell, umis) in counts.cells.iter().take(legacy_limit) {
                            println!("{}\t{umis}", unpack_cell(dict, *cell));
                        }
                    }
                    QueryAggregation::Group => {
                        println!("group\tumis\tcells\tselected_cells");
                        for (group, (umis, cells)) in counts.groups.iter().enumerate() {
                            println!(
                                "{}\t{umis}\t{cells}\t{}",
                                scope.group_names[group], scope.selected_per_group[group]
                            );
                        }
                    }
                    QueryAggregation::Bulk => {
                        println!("scope\tumis\tcells\tselected_cells");
                        println!(
                            "bulk\t{}\t{}\t{}",
                            counts.total_umis,
                            counts.cells.len(),
                            scope.selected_cells
                        );
                    }
                }
                eprintln!(
                    "region {locus}: {n_mols} molecules, {} UMIs across {} cells ({} chunks; open {t_open:.2}s, total {:.2}s)",
                    counts.total_umis,
                    counts.cells.len(),
                    n_chunks,
                    t0.elapsed().as_secs_f32()
                );
            } else if json_output {
                let mut value = json!({
                    "schema": "gravlax.query.region.v1",
                    "coordinates": "0-based half-open",
                    "anchor_semantics": true,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "molecules": n_mols,
                    "umis": counts.total_umis,
                    "cells": counts.cells.len(),
                    "scope": scope.json(),
                    "chunks_decoded": n_chunks,
                });
                let object = value.as_object_mut().unwrap();
                match scope.aggregation {
                    QueryAggregation::Cell => {
                        object.insert(
                            "cell_rows".into(),
                            json!(counts
                                .cells
                                .iter()
                                .take(legacy_limit)
                                .map(|(cell, umis)| json!({
                                    "barcode": unpack_cell(dict, *cell),
                                    "umis": umis,
                                }))
                                .collect::<Vec<_>>()),
                        );
                        object.insert(
                            "cell_rows_truncated".into(),
                            json!(legacy_limit < counts.cells.len()),
                        );
                    }
                    QueryAggregation::Group => {
                        object.insert(
                            "group_rows".into(),
                            json!(counts
                                .groups
                                .iter()
                                .enumerate()
                                .map(|(group, (umis, cells))| json!({
                                    "group": scope.group_names[group],
                                    "umis": umis,
                                    "cells": cells,
                                    "selected_cells": scope.selected_per_group[group],
                                }))
                                .collect::<Vec<_>>()),
                        );
                    }
                    QueryAggregation::Bulk => {
                        object.insert(
                            "bulk".into(),
                            json!({"umis": counts.total_umis, "cells": counts.cells.len()}),
                        );
                    }
                }
                println!("{}", serde_json::to_string_pretty(&value)?);
                eprintln!(
                    "region {locus}: open {t_open:.2}s, total {:.2}s",
                    t0.elapsed().as_secs_f32()
                );
            } else {
                println!(
                    "region {locus}: {} molecules, {} UMIs across {} cells ({} chunks decoded; open {t_open:.2}s, total {:.2}s)",
                    n_mols,
                    counts.total_umis,
                    counts.cells.len(),
                    n_chunks,
                    t0.elapsed().as_secs_f32()
                );
                match scope.aggregation {
                    QueryAggregation::Cell => {
                        for (cell, umis) in counts.cells.iter().take(legacy_limit) {
                            println!("  {}\t{}", unpack_cell(dict, *cell), umis);
                        }
                    }
                    QueryAggregation::Group => {
                        for (group, (umis, cells)) in counts.groups.iter().enumerate() {
                            println!("  {}\t{umis}\t{cells}", scope.group_names[group]);
                        }
                    }
                    QueryAggregation::Bulk => {}
                }
            }
        }
        What::Apa {
            locus,
            site_gap,
            strand,
            tsv,
            groups,
            genome,
            drop_ip,
            permute,
            seed,
            plot,
            uniform_output,
        } => {
            let permute = permute.map(std::num::NonZeroUsize::get).unwrap_or(0);
            require_uniform_genome_binding(
                &la,
                genome.as_deref(),
                uniform_output.format.is_some(),
                "APA",
            )?;
            let (group_names, cell_group, group_source_content_blake3) = match &groups {
                Some(gpath) if uniform_output.format.is_some() => {
                    let (names, mapping, digest) = load_groups_strict(&mut la, gpath)?;
                    (names, mapping, Some(digest))
                }
                Some(gpath) => {
                    let (names, mapping) = load_groups(&mut la, gpath)?;
                    (names, mapping, None)
                }
                None => (Vec::new(), FxHashMap::default(), None),
            };
            let (chrom, start, end) = parse_locus(&locus)?;
            let cid = chrom_id(&chrom)?;
            let seq = load_query_genome(&la, &genome, &chrom)?;
            let want_rev = match strand.as_deref() {
                Some("+") => Some(false),
                Some("-") => Some(true),
                None => None,
                _ => bail!("--strand must be + or -"),
            };
            let shapes = la.shapes()?;
            // (3' coordinate, strand, cell, class) for deduped molecules in the window.
            let mut pts: Vec<(u32, bool, u32, u32)> = Vec::new();
            for (i, c) in chunks.iter().enumerate() {
                if c.chrom != cid || c.bin_start >= end {
                    continue;
                }
                let raw = la.reader().read(&format!("c{i}"))?;
                for m in decode_chunk(&raw, c, None, &la.rans_tables)? {
                    let a = m.anchor();
                    if a < start.saturating_sub(2_000) || a >= end {
                        continue;
                    }
                    if let Some(wr) = want_rev {
                        if m.strand_rev != wr {
                            continue;
                        }
                    }
                    let tp = three_prime(&m, &shapes);
                    if tp >= start && tp < end {
                        pts.push((tp, m.strand_rev, la.cell_of(m.umi_class)?, m.umi_class));
                    }
                }
            }
            let sites = call_sites(
                &mut pts,
                site_gap,
                &cell_group,
                group_names.len(),
                seq.as_deref(),
            );
            let n_flagged = sites.iter().filter(|st| st.ip).count();
            let kept: Vec<&SiteCall> = sites.iter().filter(|st| !(drop_ip && st.ip)).collect();
            let mut order: Vec<&SiteCall> = kept.clone();
            order.sort_by_key(|st| (std::cmp::Reverse(st.umis), st.lo));
            let ip_cols = |st: &SiteCall| -> String {
                if seq.is_some() {
                    format!("\t{}\t{}\t{}", st.ip as u8, st.a20, st.arun)
                } else {
                    String::new()
                }
            };
            let gsuffix = |gc: &[u64]| gc.iter().map(|n| format!("\t{n}")).collect::<String>();
            let group_test = if group_names.is_empty() {
                None
            } else {
                let table: Vec<Vec<u64>> = kept.iter().map(|st| st.gc.clone()).collect();
                let (g, df) = apastats::g_statistic(&table);
                let p = apastats::chi2_sf(g, df);
                let permutation_p = if permute > 0 {
                    let (csc, glab) = per_cell_site_counts(&kept, &cell_group);
                    Some(apastats::permutation_p(
                        &csc,
                        &glab,
                        kept.len(),
                        group_names.len(),
                        g,
                        permute,
                        seed,
                    ))
                } else {
                    None
                };
                Some((g, df, p, permutation_p))
            };
            if uniform_output.format.is_some() {
                let total_umis: usize = order.iter().map(|site| site.umis).sum();
                let grouped_umis: u64 = order.iter().flat_map(|site| &site.gc).sum();
                let summary = json!({
                    "coordinates": "0-based point coordinates; start and end are inclusive observed cluster extrema",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand_filter": strand,
                    "site_gap": site_gap,
                    "semantics": {
                        "site_clustering": "single-linkage by three-prime coordinate within strand",
                        "cleavage_coordinate": "modal three-prime coordinate; ties resolve transcript-downstream",
                        "count_unit": "archive UMI class",
                        "group_counts_include_only_mapped_group_cells": !group_names.is_empty(),
                    },
                    "internal_priming": {
                        "reference_consulted": seq.is_some(),
                        "flagged_sites": n_flagged,
                        "drop_flagged": drop_ip,
                    },
                    "totals": {
                        "sites": order.len(),
                        "umis": total_umis,
                        "grouped_umis": grouped_umis,
                    },
                });
                let site_schema = TableSchema::new(
                    "gravlax.query.apa.sites.v1",
                    vec![
                        Field::new("rank", DataType::UInt64),
                        Field::new("chrom", DataType::String),
                        Field::new("start", DataType::UInt64),
                        Field::new("end", DataType::UInt64),
                        Field::new("cleavage", DataType::UInt64),
                        Field::new("strand", DataType::String),
                        Field::new("umis", DataType::UInt64),
                        Field::new("cells", DataType::UInt64),
                        Field::new("grouped_umis", DataType::UInt64),
                        Field::new("internal_priming", DataType::Boolean).nullable(),
                        Field::new("downstream_a_in_20", DataType::UInt64).nullable(),
                        Field::new("downstream_a_run", DataType::UInt64).nullable(),
                    ],
                )?
                .with_semantics(
                    TableSemantics::new(RowSemantics::Sequence)
                        .with_key(["rank"])
                        .ordered_by([OrderKey::ascending("rank")]),
                )?;
                let group_schema = TableSchema::new(
                    "gravlax.query.apa.group-counts.v1",
                    vec![
                        Field::new("site_rank", DataType::UInt64),
                        Field::new("group", DataType::String),
                        Field::new("umis", DataType::UInt64),
                        Field::new("selected_cells", DataType::UInt64),
                    ],
                )?
                .with_semantics(
                    TableSemantics::new(RowSemantics::Set).with_key(["site_rank", "group"]),
                )?;
                let test_schema = TableSchema::new(
                    "gravlax.query.apa.group-test.v1",
                    vec![
                        Field::new("test", DataType::String),
                        Field::new("sites", DataType::UInt64),
                        Field::new("groups", DataType::UInt64),
                        Field::new("g_statistic", DataType::Float64),
                        Field::new("degrees_of_freedom", DataType::UInt64),
                        Field::new("p_value", DataType::Float64),
                        Field::new("permutation_p_value", DataType::Float64).nullable(),
                        Field::new("permutations", DataType::UInt64),
                        Field::new("seed", DataType::UInt64),
                    ],
                )?
                .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["test"]))?;
                let site_selection = SelectionSummary::complete(order.len() as u64);
                let group_rows = order.len().saturating_mul(group_names.len());
                let group_selection = SelectionSummary::complete(group_rows as u64);
                let test_selection = SelectionSummary::complete(u64::from(group_test.is_some()));
                let mut parameters = BTreeMap::new();
                parameters.insert(
                    "archive_access".into(),
                    json!("chunks overlapping the locus and its upstream anchor margin"),
                );
                parameters.insert("locus".into(), json!(locus));
                parameters.insert("site_gap".into(), json!(site_gap));
                parameters.insert("strand".into(), json!(strand));
                parameters.insert("drop_internal_priming".into(), json!(drop_ip));
                parameters.insert("permutations".into(), json!(permute));
                parameters.insert("seed".into(), json!(seed));
                if let Some(path) = groups.as_deref() {
                    parameters.insert(
                        "cell_scope".into(),
                        uniform_group_scope_provenance(
                            path,
                            group_source_content_blake3.as_deref().context(
                                "uniform APA group scope is missing its bound content digest",
                            )?,
                            &group_names,
                            &cell_group,
                            la.cells()?.len(),
                        )?,
                    );
                }
                if let Some(path) = genome.as_deref() {
                    parameters.insert("reference_path".into(), json!(path));
                    parameters.insert(
                        "archive_genome_signature".into(),
                        serde_json::to_value(&la.genome_sig)?,
                    );
                }
                parameters.insert(
                    "site_order".into(),
                    json!("umis descending, then start ascending"),
                );
                let context = uniform_query_context(
                    &archive,
                    la.reader().archive_version(),
                    la.reader()
                        .content_commitment()
                        .map(|commitment| commitment.to_hex()),
                    parameters,
                );
                let selected_per_group: Vec<usize> = (0..group_names.len())
                    .map(|group| {
                        cell_group
                            .values()
                            .filter(|&&value| value as usize == group)
                            .count()
                    })
                    .collect();
                write_uniform_bundle_output(&uniform_output, |writer, format| {
                    let mut bundle = StreamingBundleWriter::new_with_summary(
                        writer,
                        "gravlax.query.apa.result.v1",
                        OutputFormat::from(format),
                        &context,
                        &summary,
                    )?;
                    bundle.write_table("sites", &site_schema, Some(&site_selection), |rows| {
                        for (rank, site) in order.iter().enumerate() {
                            rows.write_row_with(|row| {
                                row.uint64(rank as u64)?;
                                row.string(&chrom)?;
                                row.uint64(site.lo as u64)?;
                                row.uint64(site.hi as u64)?;
                                row.uint64(site.cp as u64)?;
                                row.string(if site.rev { "-" } else { "+" })?;
                                row.uint64(site.umis as u64)?;
                                row.uint64(site.cells as u64)?;
                                row.uint64(site.gc.iter().sum())?;
                                if seq.is_some() {
                                    row.boolean(site.ip)?;
                                    row.uint64(site.a20 as u64)?;
                                    row.uint64(site.arun as u64)?;
                                } else {
                                    row.null()?;
                                    row.null()?;
                                    row.null()?;
                                }
                                Ok(())
                            })?;
                        }
                        Ok(())
                    })?;
                    if !group_names.is_empty() {
                        bundle.write_table(
                            "group_counts",
                            &group_schema,
                            Some(&group_selection),
                            |rows| {
                                for (rank, site) in order.iter().enumerate() {
                                    for (group, name) in group_names.iter().enumerate() {
                                        rows.write_row_with(|row| {
                                            row.uint64(rank as u64)?;
                                            row.string(name)?;
                                            row.uint64(site.gc[group])?;
                                            row.uint64(selected_per_group[group] as u64)?;
                                            Ok(())
                                        })?;
                                    }
                                }
                                Ok(())
                            },
                        )?;
                        bundle.write_table(
                            "group_test",
                            &test_schema,
                            Some(&test_selection),
                            |rows| {
                                if let Some((g, df, p_value, permutation_p)) = group_test {
                                    rows.write_row_with(|row| {
                                        row.string("site_x_group_g_test")?;
                                        row.uint64(kept.len() as u64)?;
                                        row.uint64(group_names.len() as u64)?;
                                        row.float64(g)?;
                                        row.uint64(df)?;
                                        row.float64(p_value)?;
                                        if let Some(value) = permutation_p {
                                            row.float64(value)?;
                                        } else {
                                            row.null()?;
                                        }
                                        row.uint64(permute as u64)?;
                                        row.uint64(seed)?;
                                        Ok(())
                                    })?;
                                }
                                Ok(())
                            },
                        )?;
                    }
                    bundle.finish()?;
                    Ok(())
                })?;
            } else if tsv {
                let ghdr: String = group_names.iter().map(|g| format!("\t{g}")).collect();
                let iphdr = if seq.is_some() { "\tip\ta20\tarun" } else { "" };
                println!("#chrom\tstart\tend\tcleavage\tstrand\tumis\tcells{ghdr}{iphdr}");
                for st in &order {
                    println!(
                        "{chrom}\t{}\t{}\t{}\t{}\t{}\t{}{}{}",
                        st.lo, st.hi, st.cp, if st.rev { '-' } else { '+' }, st.umis, st.cells,
                        gsuffix(&st.gc), ip_cols(st)
                    );
                }
            } else {
                let total: usize = order.iter().map(|st| st.umis).sum();
                let ipnote = if seq.is_some() {
                    format!(
                        ", {n_flagged} flagged internal-priming{}",
                        if drop_ip { " (dropped)" } else { "" }
                    )
                } else {
                    String::new()
                };
                println!(
                    "apa {locus}: {} sites, {total} UMIs{ipnote} (open {t_open:.2}s, total {:.2}s)",
                    order.len(), t0.elapsed().as_secs_f32()
                );
                for st in order.iter().take(12) {
                    println!(
                        "  {chrom}:{}-{} ({}) {} UMIs / {} cells{}{}",
                        st.lo, st.hi, if st.rev { '-' } else { '+' }, st.umis, st.cells,
                        gsuffix(&st.gc), ip_cols(st)
                    );
                }
            }
            if let Some((g, df, p, permutation_p)) = group_test {
                let perm = permutation_p.map_or_else(String::new, |pp| {
                    format!(", permutation p = {pp:.4} ({permute} perms)")
                });
                let diagnostic = format!(
                    "site x group G-test ({} sites x {} groups): G = {g:.2}, df = {df}, p = {p:.3e}{perm}",
                    kept.len(), group_names.len()
                );
                if tsv || uniform_output.format.is_some() {
                    eprintln!("{diagnostic}");
                } else {
                    println!("{diagnostic}");
                }
            }
            if let Some(plot_path) = &plot {
                // Readability: plot only supported sites, heaviest 250.
                let mut by_umis: Vec<&SiteCall> = sites.iter().filter(|st| st.umis >= 5).collect();
                by_umis.sort_by_key(|st| std::cmp::Reverse(st.umis));
                by_umis.truncate(250);
                let dots: Vec<crate::plots::SiteDot> = by_umis
                    .iter()
                    .map(|st| crate::plots::SiteDot {
                        cp: st.cp, rev: st.rev, umis: st.umis, ip: st.ip, gc: st.gc.clone(),
                    })
                    .collect();
                crate::plots::apa_plot(
                    plot_path, &chrom, start, end, &dots, &group_names,
                    &format!("3'-site usage, {locus}"),
                )?;
                eprintln!("wrote {}", plot_path.display());
            }
        }
        What::ApaTest {
            gtf,
            groups,
            genome,
            site_gap,
            min_site_umis,
            min_gene_umis,
            tail_extend,
            permute,
            seed,
            uniform_output,
        } => {
            apa_test(
                &archive,
                &mut la,
                &chunks,
                &chrom_names,
                ApaTestParams {
                    gtf, groups, genome, site_gap, min_site_umis, min_gene_umis, tail_extend,
                    permute, seed, t0, t_open,
                },
                &uniform_output,
            )?;
        }
        What::Discover {
            gtf,
            merge_gap,
            min_umis,
            residual_min_umis,
            claim_mode,
            solo_strand,
            tsv,
            emit_gtf,
            uniform_output,
        } => {
            let (anno, annotation_content_blake3) = load_query_annotation(
                &gtf,
                "discovery annotation",
                uniform_output.format.is_some(),
            )?;
            let shapes = la.shapes()?;
            let patterns = if matches!(
                claim_mode,
                DiscoveryClaimMode::Compatible | DiscoveryClaimMode::ResidualSites
            ) {
                Some(la.patterns()?)
            } else {
                None
            };
            let anno_of: Vec<Option<u32>> = chrom_names
                .iter()
                .map(|n| anno.chrom_ids.get(n).copied())
                .collect();
            let solo_strand_name = match solo_strand {
                crate::archivecmd::SoloStrandArg::Forward => "forward",
                crate::archivecmd::SoloStrandArg::Reverse => "reverse",
                crate::archivecmd::SoloStrandArg::Unstranded => "unstranded",
            };
            let solo_strand = solo_strand.into();
            let mut candidates: Vec<(u32, u32, u32, bool, usize, usize)> = Vec::new(); // chrom, s, e, rev, umis, cells
            let mut n_unclaimed = 0u64;
            // Genome scan; bounded chunk batches decode/classify in parallel and are consumed in
            // input order, so no whole-archive molecule table is retained.
            let mut runbuf: Vec<(u32, u32, u32, bool, u32, u32)> = Vec::new(); // chrom,s,e,rev,cell,class
            let mut sitebuf: Vec<(u32, u32, u32, bool, u32, u32, u32)> = Vec::new(); // + site
            let flush_runs =
                |buf: &mut Vec<(u32, u32, u32, bool, u32, u32)>,
                 candidates: &mut Vec<(u32, u32, u32, bool, usize, usize)>| {
                    for rev in [false, true] {
                        let mut sp: Vec<&(u32, u32, u32, bool, u32, u32)> =
                            buf.iter().filter(|p| p.3 == rev).collect();
                        sp.sort_unstable_by_key(|p| p.1);
                        let mut i2 = 0usize;
                        while i2 < sp.len() {
                            let mut j2 = i2;
                            let mut span_end = sp[i2].2;
                            while j2 + 1 < sp.len() && sp[j2 + 1].1 <= span_end + merge_gap {
                                j2 += 1;
                                span_end = span_end.max(sp[j2].2);
                            }
                            let members = &sp[i2..=j2];
                            let classes: FxHashSet<u32> = members.iter().map(|p| p.5).collect();
                            if classes.len() >= min_umis {
                                let cells2: FxHashSet<u32> = members.iter().map(|p| p.4).collect();
                                candidates.push((
                                    members[0].0,
                                    members.iter().map(|p| p.1).min().unwrap(),
                                    span_end,
                                    rev,
                                    classes.len(),
                                    cells2.len(),
                                ));
                            }
                            i2 = j2 + 1;
                        }
                    }
                    buf.clear();
                };
            let flush_sites =
                |buf: &mut Vec<(u32, u32, u32, bool, u32, u32, u32)>,
                 candidates: &mut Vec<(u32, u32, u32, bool, usize, usize)>| {
                    for rev in [false, true] {
                        let mut sp: Vec<&(u32, u32, u32, bool, u32, u32, u32)> =
                            buf.iter().filter(|p| p.3 == rev).collect();
                        sp.sort_unstable_by_key(|p| p.6);
                        let mut i2 = 0usize;
                        while i2 < sp.len() {
                            let first_site = sp[i2].6;
                            let mut j2 = i2;
                            while j2 + 1 < sp.len()
                                && sp[j2 + 1].6 <= first_site.saturating_add(merge_gap)
                            {
                                j2 += 1;
                            }
                            let members = &sp[i2..=j2];
                            let classes: FxHashSet<u32> = members.iter().map(|p| p.5).collect();
                            if classes.len() >= residual_min_umis {
                                let cells2: FxHashSet<u32> = members.iter().map(|p| p.4).collect();
                                candidates.push((
                                    members[0].0,
                                    members.iter().map(|p| p.6).min().unwrap(),
                                    members.iter().map(|p| p.6).max().unwrap().saturating_add(1),
                                    rev,
                                    classes.len(),
                                    cells2.len(),
                                ));
                            }
                            i2 = j2 + 1;
                        }
                    }
                    buf.clear();
                };
            let mut last_chrom = u32::MAX;
            // Pair sparse access units to overlap decode/classification, but process dense units
            // alone: the molecule cap is the memory guard, independent of thread count.
            const MAX_BATCH_MOLECULES: u64 = 500_000;
            let mut first = 0usize;
            while first < chunks.len() {
                let mut end = first + 1;
                if end < chunks.len()
                    && chunks[first].n_mols as u64 + chunks[end].n_mols as u64
                        <= MAX_BATCH_MOLECULES
                {
                    end += 1;
                }
                let batch = &chunks[first..end];
                let unclaimed: Vec<Vec<DiscoveryUnclaimed>> = {
                    let (reader, tables) = la.reader_and_tables();
                    let reader = &*reader;
                    batch.par_iter().enumerate().map(|(j, info)| {
                        let i = first + j;
                        let (compressed, raw_len) = reader.read_compressed_at(&format!("c{i}"))?;
                        let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                        let molecules = decode_chunk(&raw, info, None, tables)?;
                        Ok(discovery_unclaimed(
                            &molecules, &shapes, patterns.as_deref(), &anno, &anno_of,
                            claim_mode, solo_strand,
                        ))
                    }).collect::<Result<_>>()?
                };
                la.prefetch_coc(unclaimed.iter().flatten().map(|row| row.4))?;
                for (j, rows) in unclaimed.into_iter().enumerate() {
                    let chrom = batch[j].chrom;
                    if chrom != last_chrom {
                        flush_runs(&mut runbuf, &mut candidates);
                        flush_sites(&mut sitebuf, &mut candidates);
                        last_chrom = chrom;
                    }
                    n_unclaimed += rows.len() as u64;
                    for (chrom, lo, hi, rev, class, site) in rows {
                        let cell = la.cell_of(class)?;
                        match site {
                            Some(site) => sitebuf.push((chrom, lo, hi, rev, cell, class, site)),
                            None => runbuf.push((chrom, lo, hi, rev, cell, class)),
                        }
                    }
                }
                first = end;
            }
            flush_runs(&mut runbuf, &mut candidates);
            flush_sites(&mut sitebuf, &mut candidates);
            candidates.sort_unstable_by_key(|c2| (std::cmp::Reverse(c2.4), c2.0, c2.1));
            if let Some(gtf_out) = &emit_gtf {
                use std::io::Write as _;
                let mut w = std::io::BufWriter::new(std::fs::File::create(gtf_out)?);
                let mut by_pos: Vec<&(u32, u32, u32, bool, usize, usize)> = candidates.iter().collect();
                by_pos.sort_unstable_by_key(|c2| (c2.0, c2.1));
                for (k, (ch, s2, e2, rev, _u, _c)) in by_pos.iter().enumerate() {
                    let (cn, st) = (&chrom_names[*ch as usize], if *rev { '-' } else { '+' });
                    let gid = format!("AIENOVEL{k:06}");
                    for feat in ["gene", "transcript", "exon"] {
                        writeln!(
                            w,
                            "{cn}\taie\t{feat}\t{}\t{e2}\t.\t{st}\t.\tgene_id \"{gid}\"; transcript_id \"{gid}.1\"; gene_name \"{gid}\";",
                            s2 + 1
                        )?;
                    }
                }
                eprintln!("wrote {} candidate loci to {}", candidates.len(), gtf_out.display());
            }
            if uniform_output.format.is_some() {
                let summary = json!({
                    "coordinates": "0-based half-open candidate extents",
                    "semantics": {
                        "candidate_source": "archive molecules unclaimed by the supplied annotation",
                        "claim_mode": claim_mode.name(),
                        "solo_strand": solo_strand_name,
                        "count_unit": "distinct archive UMI class within each candidate",
                        "candidate_rows_may_share_coordinates": true,
                    },
                    "thresholds": {
                        "merge_gap": merge_gap,
                        "span_min_umis": min_umis,
                        "residual_site_min_umis": residual_min_umis,
                    },
                    "planning": {
                        "archive_chunks": chunks.len(),
                        "max_decode_batch_molecules": MAX_BATCH_MOLECULES,
                    },
                    "unclaimed_molecules": n_unclaimed,
                    "candidate_loci": candidates.len(),
                });
                let schema = TableSchema::new(
                    "gravlax.query.discover.candidates.v1",
                    vec![
                        Field::new("chrom", DataType::String),
                        Field::new("start", DataType::UInt64),
                        Field::new("end", DataType::UInt64),
                        Field::new("strand", DataType::String),
                        Field::new("umis", DataType::UInt64),
                        Field::new("cells", DataType::UInt64),
                    ],
                )?
                .with_semantics(TableSemantics::new(RowSemantics::Multiset))?;
                let selection = SelectionSummary::complete(candidates.len() as u64);
                let mut parameters = BTreeMap::new();
                parameters.insert(
                    "archive_access".into(),
                    json!("bounded complete archive scan"),
                );
                parameters.insert("annotation_path".into(), json!(gtf));
                parameters.insert(
                    "annotation_content_blake3".into(),
                    json!(annotation_content_blake3
                        .as_deref()
                        .context("uniform discovery annotation is missing its bound content digest")?),
                );
                parameters.insert("merge_gap".into(), json!(merge_gap));
                parameters.insert("span_min_umis".into(), json!(min_umis));
                parameters.insert("residual_site_min_umis".into(), json!(residual_min_umis));
                parameters.insert("claim_mode".into(), json!(claim_mode.name()));
                parameters.insert("solo_strand".into(), json!(solo_strand_name));
                parameters.insert(
                    "presentation_order".into(),
                    json!(
                        "umis descending, chromosome id and start ascending; not a selection rule"
                    ),
                );
                let context = uniform_query_context(
                    &archive,
                    la.reader().archive_version(),
                    la.reader()
                        .content_commitment()
                        .map(|commitment| commitment.to_hex()),
                    parameters,
                );
                write_uniform_bundle_output(&uniform_output, |writer, format| {
                    let mut bundle = StreamingBundleWriter::new_with_summary(
                        writer,
                        "gravlax.query.discover.result.v1",
                        OutputFormat::from(format),
                        &context,
                        &summary,
                    )?;
                    bundle.write_table("candidates", &schema, Some(&selection), |rows| {
                        for (chrom_id, start, end, reverse, umis, cells) in &candidates {
                            rows.write_row_with(|row| {
                                row.string(&chrom_names[*chrom_id as usize])?;
                                row.uint64(*start as u64)?;
                                row.uint64(*end as u64)?;
                                row.string(if *reverse { "-" } else { "+" })?;
                                row.uint64(*umis as u64)?;
                                row.uint64(*cells as u64)?;
                                Ok(())
                            })?;
                        }
                        Ok(())
                    })?;
                    bundle.finish()?;
                    Ok(())
                })?;
                eprintln!(
                    "discover: {} unclaimed molecules -> {} candidate loci (uniform; open {t_open:.2}s, total {:.2}s)",
                    n_unclaimed,
                    candidates.len(),
                    t0.elapsed().as_secs_f32()
                );
            } else if tsv {
                for (ch, s2, e2, rev, u, c2) in &candidates {
                    println!(
                        "{}	{s2}	{e2}	{}	{u}	{c2}",
                        chrom_names[*ch as usize],
                        if *rev { '-' } else { '+' }
                    );
                }
            } else {
                let threshold_note = if claim_mode == DiscoveryClaimMode::ResidualSites {
                    format!("span>={min_umis}, residual>={residual_min_umis} UMIs")
                } else {
                    format!(">={min_umis} UMIs")
                };
                println!(
                    "discover: {} unclaimed molecules -> {} candidate loci ({threshold_note}; claim={claim_mode:?}) (open {t_open:.2}s, total {:.2}s)",
                    n_unclaimed, candidates.len(), t0.elapsed().as_secs_f32()
                );
                for (ch, s2, e2, rev, u, c2) in candidates.iter().take(12) {
                    println!(
                        "  {}:{s2}-{e2} ({}) {u} UMIs / {c2} cells",
                        chrom_names[*ch as usize],
                        if *rev { '-' } else { '+' }
                    );
                }
            }
        }
        What::Junction {
            locus,
            top,
            tsv,
            json: json_output,
            uniform_output,
            scope: scope_args,
        } => {
            let mut scope = load_query_scope(&mut la, &scope_args)?;
            if uniform_output.format.is_some() {
                scope.ensure_resolved_mapping_digest()?;
            }
            let (chrom, donor, acceptor) = parse_locus(&locus)?;
            let cid = chrom_id(&chrom)?;
            let Some((total_support, n_chunks2, per_cell)) =
                junction_counts(&mut la, &chunks, cid, donor, acceptor)?
            else {
                bail!("junction {locus} not present in the archive");
            };
            let archive_identity = uniform_output.format.map(|_| {
                (
                    la.reader().archive_version(),
                    la.reader()
                        .content_commitment()
                        .map(|commitment| commitment.to_hex()),
                )
            });
            let dict = la.cells()?;
            let counts = scoped_counts_with_cell_order(
                &per_cell,
                &scope,
                uniform_output.format.map(|_| dict),
            );
            let legacy_limit = if top == 0 {
                counts.cells.len()
            } else {
                top.min(counts.cells.len())
            };
            let (_, uniform_limit) = uniform_row_limit(&counts, &scope, top);
            if let Some(format) = uniform_output.format {
                let summary = JunctionUniformSummary {
                    coordinates: "0-based half-open junction boundaries",
                    chrom: &chrom,
                    donor,
                    acceptor,
                    archive_supporting_children: total_support,
                    archive_posting_chunks: n_chunks2 as u64,
                    umis: counts.total_umis as u64,
                    cells: counts.cells.len() as u64,
                };
                let parameters = uniform_query_parameters(
                    &scope,
                    top,
                    "junction catalogue and postings-selected archive chunks",
                )?;
                let (archive_version, archive_root) =
                    archive_identity.expect("uniform output captured archive identity");
                let context =
                    uniform_query_context(&archive, archive_version, archive_root, parameters);
                write_uniform_count_output(
                    &uniform_output,
                    JUNCTION_UNIFORM_RESULT_SCHEMA,
                    JUNCTION_UNIFORM_COUNTS_SCHEMA,
                    &context,
                    &summary,
                    &counts,
                    &scope,
                    dict,
                    uniform_limit,
                )?;
                eprintln!(
                    "junction {locus}: {} UMIs across {} cells; {total_support} supporting children, {n_chunks2} chunks (uniform {format:?}; open {t_open:.2}s, total {:.2}s)",
                    counts.total_umis,
                    counts.cells.len(),
                    t0.elapsed().as_secs_f32()
                );
            } else if tsv {
                match scope.aggregation {
                    QueryAggregation::Cell => {
                        println!("barcode\tumis");
                        for (cell, umis) in counts.cells.iter().take(legacy_limit) {
                            println!("{}\t{umis}", unpack_cell(dict, *cell));
                        }
                    }
                    QueryAggregation::Group => {
                        println!("group\tumis\tcells\tselected_cells");
                        for (group, (umis, cells)) in counts.groups.iter().enumerate() {
                            println!(
                                "{}\t{umis}\t{cells}\t{}",
                                scope.group_names[group], scope.selected_per_group[group]
                            );
                        }
                    }
                    QueryAggregation::Bulk => {
                        println!("scope\tumis\tcells\tselected_cells");
                        println!(
                            "bulk\t{}\t{}\t{}",
                            counts.total_umis,
                            counts.cells.len(),
                            scope.selected_cells
                        );
                    }
                }
                eprintln!(
                    "junction {locus}: {} UMIs across {} cells; {total_support} supporting children, {n_chunks2} chunks; open {t_open:.2}s, total {:.2}s",
                    counts.total_umis,
                    counts.cells.len(),
                    t0.elapsed().as_secs_f32()
                );
            } else if json_output {
                let mut value = json!({
                    "schema": if scope.active { "gravlax.query.junction.v2" } else { "gravlax.query.junction.v1" },
                    "coordinates": "0-based half-open junction boundaries",
                    "chrom": chrom,
                    "donor": donor,
                    "acceptor": acceptor,
                    "supporting_children": total_support,
                    "posting_chunks": n_chunks2,
                    "umis": counts.total_umis,
                    "cells": counts.cells.len(),
                });
                let object = value.as_object_mut().unwrap();
                if scope.active {
                    object.insert("scope".into(), scope.json());
                }
                match scope.aggregation {
                    QueryAggregation::Cell => {
                        object.insert(
                            "cell_rows".into(),
                            json!(counts
                                .cells
                                .iter()
                                .take(legacy_limit)
                                .map(|(cell, umis)| json!({
                                    "barcode": unpack_cell(dict, *cell),
                                    "umis": umis,
                                }))
                                .collect::<Vec<_>>()),
                        );
                        object.insert(
                            "cell_rows_truncated".into(),
                            json!(legacy_limit < counts.cells.len()),
                        );
                    }
                    QueryAggregation::Group => {
                        object.insert(
                            "group_rows".into(),
                            json!(counts
                                .groups
                                .iter()
                                .enumerate()
                                .map(|(group, (umis, cells))| json!({
                                    "group": scope.group_names[group],
                                    "umis": umis,
                                    "cells": cells,
                                    "selected_cells": scope.selected_per_group[group],
                                }))
                                .collect::<Vec<_>>()),
                        );
                    }
                    QueryAggregation::Bulk => {
                        object.insert(
                            "bulk".into(),
                            json!({"umis": counts.total_umis, "cells": counts.cells.len()}),
                        );
                    }
                }
                println!("{}", serde_json::to_string_pretty(&value)?);
                eprintln!(
                    "junction {locus}: open {t_open:.2}s, total {:.2}s",
                    t0.elapsed().as_secs_f32()
                );
            } else {
                println!(
                    "junction {locus}: {} UMIs across {} cells ({} supporting children per index, {} chunks decoded; open {t_open:.2}s, total {:.2}s)",
                    counts.total_umis,
                    counts.cells.len(),
                    total_support,
                    n_chunks2,
                    t0.elapsed().as_secs_f32()
                );
                match scope.aggregation {
                    QueryAggregation::Cell => {
                        for (cell, umis) in counts.cells.iter().take(legacy_limit) {
                            println!("  {}\t{}", unpack_cell(dict, *cell), umis);
                        }
                    }
                    QueryAggregation::Group => {
                        for (group, (umis, cells)) in counts.groups.iter().enumerate() {
                            println!("  {}\t{umis}\t{cells}", scope.group_names[group]);
                        }
                    }
                    QueryAggregation::Bulk => {}
                }
            }
        }
        What::Junctions {
            locus,
            either,
            min_support,
            with_cells,
            min_cells,
            gtf,
            tsv,
            json: json_output,
            uniform_output,
            scope: scope_args,
        } => {
            let mut scope = load_query_scope(&mut la, &scope_args)?;
            if uniform_output.format.is_some() {
                scope.ensure_resolved_mapping_digest()?;
            }
            let (chrom, start, end) = parse_locus(&locus)?;
            let cid = chrom_id(&chrom)?;
            let metadata = read_junction_metadata(&mut la)?;
            let mut selected: Vec<JunctionMeta> = metadata
                .into_iter()
                .filter(|row| {
                    row.chrom == cid
                        && row.supporting_children >= min_support
                        && junction_in_window(row, start, end, either)
                })
                .collect();
            let decode_cells = with_cells || min_cells > 0 || scope.active;
            let mut counts = if decode_cells {
                junction_counts_many(&mut la, &chunks, &selected)?
            } else {
                Vec::new()
            };
            if scope.active {
                for per_cell in &mut counts {
                    per_cell.retain(|cell, _| scope.includes(*cell));
                }
            }
            if decode_cells && min_cells > 0 {
                let keep: Vec<bool> = counts.iter().map(|per_cell| per_cell.len() >= min_cells).collect();
                selected = selected
                    .into_iter()
                    .zip(&keep)
                    .filter_map(|(row, keep)| (*keep).then_some(row))
                    .collect();
            }
            // Counts were computed before the optional filter; retain the matching entries.
            let selected_counts: Vec<&FxHashMap<u32, FxHashSet<u32>>> = if decode_cells {
                counts
                    .iter()
                    .filter(|per_cell| min_cells == 0 || per_cell.len() >= min_cells)
                    .collect()
            } else {
                Vec::new()
            };
            let (annotation, annotation_content_blake3) = match &gtf {
                Some(path) => {
                    let (annotation, digest) = load_junction_annotation(
                        path,
                        &chrom,
                        uniform_output.format.is_some(),
                    )?;
                    (Some(annotation), digest)
                }
                None => (None, None),
            };
            let cell_dictionary =
                if decode_cells && (json_output || uniform_output.format.is_some()) {
                    Some(la.cells()?.to_vec())
                } else {
                    None
                };

            if uniform_output.format.is_some() {
                let summary = json!({
                    "coordinates": "0-based half-open",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "interval_semantics": if either {
                        "either_endpoint"
                    } else {
                        "both_endpoints_contained"
                    },
                    "min_archive_supporting_children": min_support,
                    "min_scoped_cells": min_cells,
                    "exact_cell_counts": decode_cells,
                    "annotation_strand": if annotation.is_some() {
                        "either"
                    } else {
                        "not_supplied"
                    },
                    "junctions": selected.len(),
                    "scope": scope.json(),
                });
                let junction_schema = TableSchema::new(
                    "gravlax.query.junctions.junctions.v1",
                    vec![
                        Field::new("chrom", DataType::String),
                        Field::new("donor", DataType::UInt64),
                        Field::new("acceptor", DataType::UInt64),
                        Field::new("archive_supporting_children", DataType::UInt64),
                        Field::new("archive_posting_chunks", DataType::UInt64),
                        Field::new("umis", DataType::UInt64).nullable(),
                        Field::new("cells", DataType::UInt64).nullable(),
                        Field::new("annotated", DataType::Boolean).nullable(),
                        Field::new("donor_annotated", DataType::Boolean).nullable(),
                        Field::new("acceptor_annotated", DataType::Boolean).nullable(),
                    ],
                )?
                .with_semantics(
                    TableSemantics::new(RowSemantics::Set).with_key(["chrom", "donor", "acceptor"]),
                )?;
                let count_schema = TableSchema::new(
                    "gravlax.query.junctions.counts.v1",
                    vec![
                        Field::new("donor", DataType::UInt64),
                        Field::new("acceptor", DataType::UInt64),
                        Field::new("aggregation", DataType::String),
                        Field::new("entity", DataType::String),
                        Field::new("umis", DataType::UInt64),
                        Field::new("cells", DataType::UInt64).nullable(),
                        Field::new("selected_cells", DataType::UInt64).nullable(),
                    ],
                )?
                .with_semantics(
                    TableSemantics::new(RowSemantics::Set).with_key([
                        "donor",
                        "acceptor",
                        "aggregation",
                        "entity",
                    ]),
                )?;
                let junction_selection = SelectionSummary::complete(selected.len() as u64);
                let count_rows = if decode_cells {
                    match scope.aggregation {
                        QueryAggregation::Cell => selected_counts
                            .iter()
                            .map(|per_cell| per_cell.len())
                            .sum::<usize>(),
                        QueryAggregation::Group => {
                            selected.len().saturating_mul(scope.group_names.len())
                        }
                        QueryAggregation::Bulk => selected.len(),
                    }
                } else {
                    0
                };
                let count_selection = SelectionSummary::complete(count_rows as u64);
                let mut parameters = BTreeMap::new();
                parameters.insert(
                    "archive_access".into(),
                    json!(if decode_cells {
                        "junction catalogue and selected postings"
                    } else {
                        "junction catalogue only"
                    }),
                );
                parameters.insert("locus".into(), json!(locus));
                parameters.insert("interval_either_endpoint".into(), json!(either));
                parameters.insert("min_archive_supporting_children".into(), json!(min_support));
                parameters.insert("with_cells".into(), json!(with_cells));
                parameters.insert("min_scoped_cells".into(), json!(min_cells));
                parameters.insert("cell_scope".into(), scope.provenance_json());
                parameters.insert("aggregation".into(), json!(scope.aggregation_name()));
                if let Some(path) = gtf.as_deref() {
                    parameters.insert("annotation_path".into(), json!(path));
                    parameters.insert(
                        "annotation_content_blake3".into(),
                        json!(annotation_content_blake3
                            .as_deref()
                            .context("uniform junction annotation is missing its bound content digest")?),
                    );
                }
                let context = uniform_query_context(
                    &archive,
                    la.reader().archive_version(),
                    la.reader()
                        .content_commitment()
                        .map(|commitment| commitment.to_hex()),
                    parameters,
                );
                write_uniform_bundle_output(&uniform_output, |writer, format| {
                    let mut bundle = StreamingBundleWriter::new_with_summary(
                        writer,
                        "gravlax.query.junctions.result.v1",
                        OutputFormat::from(format),
                        &context,
                        &summary,
                    )?;
                    bundle.write_table(
                        "junctions",
                        &junction_schema,
                        Some(&junction_selection),
                        |rows| {
                            for (index, junction) in selected.iter().enumerate() {
                                let exact = decode_cells.then(|| {
                                    let per_cell = selected_counts[index];
                                    (
                                        per_cell.values().map(FxHashSet::len).sum::<usize>(),
                                        per_cell.len(),
                                    )
                                });
                                rows.write_row_with(|row| {
                                    row.string(&chrom)?;
                                    row.uint64(junction.donor as u64)?;
                                    row.uint64(junction.acceptor as u64)?;
                                    row.uint64(junction.supporting_children)?;
                                    row.uint64(junction.posts.len() as u64)?;
                                    if let Some((umis, cells)) = exact {
                                        row.uint64(umis as u64)?;
                                        row.uint64(cells as u64)?;
                                    } else {
                                        row.null()?;
                                        row.null()?;
                                    }
                                    if let Some(annotation) = &annotation {
                                        row.boolean(
                                            annotation
                                                .exact
                                                .contains(&(junction.donor, junction.acceptor)),
                                        )?;
                                        row.boolean(annotation.donors.contains(&junction.donor))?;
                                        row.boolean(
                                            annotation.acceptors.contains(&junction.acceptor),
                                        )?;
                                    } else {
                                        row.null()?;
                                        row.null()?;
                                        row.null()?;
                                    }
                                    Ok(())
                                })?;
                            }
                            Ok(())
                        },
                    )?;
                    if decode_cells {
                        let dictionary = cell_dictionary
                            .as_deref()
                            .expect("uniform cell rows loaded the archive dictionary");
                        bundle.write_table(
                            "counts",
                            &count_schema,
                            Some(&count_selection),
                            |rows| {
                                for (index, junction) in selected.iter().enumerate() {
                                    let counts =
                                        scoped_counts_unordered(selected_counts[index], &scope);
                                    match scope.aggregation {
                                        QueryAggregation::Cell => {
                                            for (cell, umis) in &counts.cells {
                                                let barcode =
                                                    unpack_cell_bytes(dictionary[*cell as usize]);
                                                let barcode = std::str::from_utf8(&barcode)
                                                    .expect("packed barcode is ASCII");
                                                rows.write_row_with(|row| {
                                                    row.uint64(junction.donor as u64)?;
                                                    row.uint64(junction.acceptor as u64)?;
                                                    row.string("cell")?;
                                                    row.string(barcode)?;
                                                    row.uint64(*umis as u64)?;
                                                    row.null()?;
                                                    row.null()?;
                                                    Ok(())
                                                })?;
                                            }
                                        }
                                        QueryAggregation::Group => {
                                            for (group, (umis, cells)) in
                                                counts.groups.iter().enumerate()
                                            {
                                                rows.write_row_with(|row| {
                                                    row.uint64(junction.donor as u64)?;
                                                    row.uint64(junction.acceptor as u64)?;
                                                    row.string("group")?;
                                                    row.string(&scope.group_names[group])?;
                                                    row.uint64(*umis as u64)?;
                                                    row.uint64(*cells as u64)?;
                                                    row.uint64(
                                                        scope.selected_per_group[group] as u64,
                                                    )?;
                                                    Ok(())
                                                })?;
                                            }
                                        }
                                        QueryAggregation::Bulk => {
                                            rows.write_row_with(|row| {
                                                row.uint64(junction.donor as u64)?;
                                                row.uint64(junction.acceptor as u64)?;
                                                row.string("bulk")?;
                                                row.string("bulk")?;
                                                row.uint64(counts.total_umis as u64)?;
                                                row.uint64(counts.cells.len() as u64)?;
                                                row.uint64(scope.selected_cells as u64)?;
                                                Ok(())
                                            })?;
                                        }
                                    }
                                }
                                Ok(())
                            },
                        )?;
                    }
                    bundle.finish()?;
                    Ok(())
                })?;
                eprintln!(
                    "junctions {locus}: {} rows ({} path; uniform; open {t_open:.2}s, total {:.2}s)",
                    selected.len(),
                    if decode_cells { "postings" } else { "index-only" },
                    t0.elapsed().as_secs_f32()
                );
            } else if tsv {
                print!("chrom\tdonor\tacceptor\tsupporting_children\tposting_chunks");
                if decode_cells {
                    print!("\tumis\tcells");
                }
                if decode_cells && scope.aggregation == QueryAggregation::Group {
                    print!("\tgroup\tgroup_umis\tgroup_cells\tgroup_selected_cells");
                }
                if annotation.is_some() {
                    print!("\tannotated\tdonor_annotated\tacceptor_annotated");
                }
                println!();
                for (index, row) in selected.iter().enumerate() {
                    let summary =
                        decode_cells.then(|| scoped_counts(selected_counts[index], &scope));
                    let group_repetitions =
                        if decode_cells && scope.aggregation == QueryAggregation::Group {
                            scope.group_names.len()
                        } else {
                            1
                        };
                    for group in 0..group_repetitions {
                        print!(
                            "{chrom}\t{}\t{}\t{}\t{}",
                            row.donor,
                            row.acceptor,
                            row.supporting_children,
                            row.posts.len()
                        );
                        if let Some(summary) = &summary {
                            print!("\t{}\t{}", summary.total_umis, summary.cells.len());
                            if scope.aggregation == QueryAggregation::Group {
                                let (umis, cells) = summary.groups[group];
                                print!(
                                    "\t{}\t{umis}\t{cells}\t{}",
                                    scope.group_names[group], scope.selected_per_group[group]
                                );
                            }
                        }
                        if let Some(a) = &annotation {
                            print!(
                                "\t{}\t{}\t{}",
                                u8::from(a.exact.contains(&(row.donor, row.acceptor))),
                                u8::from(a.donors.contains(&row.donor)),
                                u8::from(a.acceptors.contains(&row.acceptor))
                            );
                        }
                        println!();
                    }
                }
                eprintln!(
                    "junctions {locus}: {} rows ({} path; open {t_open:.2}s, total {:.2}s)",
                    selected.len(),
                    if decode_cells { "postings" } else { "index-only" },
                    t0.elapsed().as_secs_f32()
                );
            } else if json_output {
                let rows: Vec<_> = selected
                    .iter()
                    .enumerate()
                    .map(|(index, row)| {
                        let mut value = json!({
                            "chrom": chrom,
                            "donor": row.donor,
                            "acceptor": row.acceptor,
                            "supporting_children": row.supporting_children,
                            "posting_chunks": row.posts.len(),
                        });
                        let object = value.as_object_mut().unwrap();
                        if decode_cells {
                            let per_cell = selected_counts[index];
                            let summary = scoped_counts(per_cell, &scope);
                            object.insert("umis".into(), json!(summary.total_umis));
                            object.insert("cells".into(), json!(summary.cells.len()));
                            let dictionary = cell_dictionary.as_ref().unwrap();
                            match scope.aggregation {
                                QueryAggregation::Cell => {
                                    object.insert(
                                        "cell_counts".into(),
                                        json!(summary
                                            .cells
                                            .into_iter()
                                            .map(|(cell, umis)| json!({
                                                "barcode": unpack_cell(dictionary, cell),
                                                "umis": umis,
                                            }))
                                            .collect::<Vec<_>>()),
                                    );
                                }
                                QueryAggregation::Group => {
                                    object.insert(
                                        "group_counts".into(),
                                        json!(summary
                                            .groups
                                            .iter()
                                            .enumerate()
                                            .map(|(group, (umis, cells))| json!({
                                                "group": scope.group_names[group],
                                                "umis": umis,
                                                "cells": cells,
                                                "selected_cells": scope.selected_per_group[group],
                                            }))
                                            .collect::<Vec<_>>()),
                                    );
                                }
                                QueryAggregation::Bulk => {}
                            }
                        }
                        if let Some(a) = &annotation {
                            object.insert(
                                "annotated".into(),
                                json!(a.exact.contains(&(row.donor, row.acceptor))),
                            );
                            object.insert(
                                "donor_annotated".into(),
                                json!(a.donors.contains(&row.donor)),
                            );
                            object.insert(
                                "acceptor_annotated".into(),
                                json!(a.acceptors.contains(&row.acceptor)),
                            );
                        }
                        value
                    })
                    .collect();
                let mut value = json!({
                    "schema": if scope.active { "gravlax.query.junctions.v2" } else { "gravlax.query.junctions.v1" },
                    "coordinates": "0-based half-open",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "interval_semantics": if either { "either_endpoint" } else { "both_endpoints_contained" },
                    "min_supporting_children": min_support,
                    "min_cells": min_cells,
                    "exact_cell_counts": decode_cells,
                    "annotation_strand": if annotation.is_some() { "either" } else { "not_supplied" },
                    "junctions": rows,
                });
                if scope.active {
                    value
                        .as_object_mut()
                        .unwrap()
                        .insert("scope".into(), scope.json());
                }
                println!("{}", serde_json::to_string_pretty(&value)?);
                eprintln!(
                    "junctions {locus}: {} rows ({} path; open {t_open:.2}s, total {:.2}s)",
                    selected.len(),
                    if decode_cells { "postings" } else { "index-only" },
                    t0.elapsed().as_secs_f32()
                );
            } else {
                println!(
                    "junctions {locus}: {} rows ({} path; open {t_open:.2}s, total {:.2}s)",
                    selected.len(),
                    if decode_cells { "postings" } else { "index-only" },
                    t0.elapsed().as_secs_f32()
                );
                for (index, row) in selected.iter().enumerate() {
                    print!(
                        "  {chrom}:{}-{} support={} chunks={}",
                        row.donor, row.acceptor, row.supporting_children, row.posts.len()
                    );
                    if decode_cells {
                        let per_cell = selected_counts[index];
                        let umis: usize = per_cell.values().map(FxHashSet::len).sum();
                        print!(" umis={umis} cells={}", per_cell.len());
                    }
                    if let Some(a) = &annotation {
                        print!(
                            " annotated={}",
                            a.exact.contains(&(row.donor, row.acceptor))
                        );
                    }
                    println!();
                }
            }
        }
    }
    Ok(())
}

/// Exact junction support, decoded chunk count, and per-cell molecule classes.
pub type JunctionCountResult = (u64, usize, FxHashMap<u32, FxHashSet<u32>>);

/// The junction query core, shared by `query junction` and `federate`: catalogue lookup,
/// postings walk, chunk decode, per-cell class-deduplicated counts. None if the junction is not
/// in this archive's catalogue.
pub fn junction_counts(
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    cid: u32,
    donor: u32,
    acceptor: u32,
) -> Result<Option<JunctionCountResult>> {
    let cat = la.reader().read("index.junctions")?;
    let mut c = Cursor::new(&cat);
    let mut jid = 0u32;
    let mut found: Option<u32> = None;
    let (mut lch, mut ld) = (u32::MAX, 0u32);
    while !c.is_empty() {
        let ch = c.varint()? as u32;
        if ch != lch {
            lch = ch;
            ld = 0;
        }
        let dn = ld + c.varint()? as u32;
        ld = dn;
        let ac = dn + c.varint()? as u32;
        if ch == cid && dn == donor && ac == acceptor {
            found = Some(jid);
            break;
        }
        jid += 1;
    }
    let Some(jid) = found else { return Ok(None) };
    let jp = la.reader().read("index.jpost")?;
    let mut pc = Cursor::new(&jp);
    let mut posts: Vec<u32> = Vec::new();
    let mut total_support = 0u64;
    for k in 0..=jid {
        let support = pc.varint()?;
        if k == jid {
            total_support = support;
        }
        let n = pc.varint()? as usize;
        let mut last = 0u32;
        let mut cur = Vec::with_capacity(n);
        for _ in 0..n {
            last += pc.varint()? as u32;
            cur.push(last);
        }
        if k == jid {
            posts = cur;
        }
    }
    junction_counts_routed(la, chunks, cid, donor, acceptor, total_support, &posts).map(Some)
}

/// Exact point-junction reduction when a collection planner has already resolved the archive-local
/// catalogue row. This deliberately skips `index.junctions` and `index.jpost`: the collection
/// sidecar is authoritative only for routing, while molecule-class membership is still recomputed
/// from the immutable archive chunks.
pub(crate) fn junction_counts_routed(
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    chrom: u32,
    donor: u32,
    acceptor: u32,
    total_support: u64,
    posts: &[u32],
) -> Result<JunctionCountResult> {
    junction_counts_routed_with_shape_route(
        la,
        chunks,
        chrom,
        donor,
        acceptor,
        total_support,
        posts,
        None,
    )
}

/// Route-aware form of [`junction_counts_routed`]. A collection-local span route replaces the
/// source archive's complete shape dictionary, but never replaces the source molecule records or
/// their class-to-cell mapping. Every candidate placement is checked against both requested
/// genomic boundaries by [`crate::shaperoute::SpanRoute::matches`].
#[expect(
    clippy::too_many_arguments,
    reason = "the routed junction kernel keeps exact coordinates, support, postings, and optional shape routing explicit"
)]
pub(crate) fn junction_counts_routed_with_shape_route(
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    chrom: u32,
    donor: u32,
    acceptor: u32,
    total_support: u64,
    posts: &[u32],
    shape_route: Option<(&crate::shaperoute::SpanRoute, u32)>,
) -> Result<JunctionCountResult> {
    if let Some(&bad) = posts.iter().find(|&&post| post as usize >= chunks.len()) {
        bail!("junction posting references missing chunk {bad}");
    }
    if let Some(&bad) = posts
        .iter()
        .find(|&&post| chunks[post as usize].chrom != chrom)
    {
        bail!(
            "junction posting references chromosome {} chunk {bad} for chromosome {chrom}",
            chunks[bad as usize].chrom
        );
    }
    let shapes = shape_route.is_none().then(|| la.shapes()).transpose()?;
    // Fused parallel pread + decompress + decode + junction-scan; only the per-cell aggregation
    // (which touches the coc cache) stays serial.
    let hits: Vec<Vec<u32>> = {
        let (reader, tables) = la.reader_and_tables();
        let reader = &*reader;
        posts
        .par_iter()
        .map(|i| {
            let (c2, raw_len) = reader.read_compressed_at(&format!("c{i}"))?;
            let raw = evidence_io::format::decompress(&c2, raw_len)?;
            let mols = decode_chunk(&raw, &chunks[*i as usize], None, tables)?;
            let mut classes = Vec::new();
            for m in &mols {
                let mut hit = false;
                for ch2 in &m.chains {
                    for (pos, shape) in &ch2.reps {
                        let placement_hit = match shape_route {
                            Some((route, n_shapes)) => {
                                if *shape >= n_shapes {
                                    bail!("molecule references shape {shape} outside the routed {n_shapes}-shape dictionary");
                                }
                                route.matches(*shape, *pos, donor, acceptor)?
                            }
                            None => {
                                let shape = shapes.as_ref().unwrap().get(*shape as usize)
                                    .with_context(|| format!("molecule references missing shape {shape}"))?;
                                has_junction(shape, *pos, donor, acceptor)?
                            }
                        };
                        if placement_hit {
                            hit = true;
                        }
                    }
                }
                for (pos, shape, _, _) in &m.mms {
                    let placement_hit = match shape_route {
                        Some((route, n_shapes)) => {
                            if *shape >= n_shapes {
                                bail!("multimapper references shape {shape} outside the routed {n_shapes}-shape dictionary");
                            }
                            route.matches(*shape, *pos, donor, acceptor)?
                        }
                        None => {
                            let shape = shapes.as_ref().unwrap().get(*shape as usize)
                                .with_context(|| format!("multimapper references missing shape {shape}"))?;
                            has_junction(shape, *pos, donor, acceptor)?
                        }
                    };
                    if placement_hit {
                        hit = true;
                    }
                }
                if hit {
                    classes.push(m.umi_class);
                }
            }
            Ok(classes)
        })
        .collect::<Result<_>>()?
    };
    la.prefetch_coc(hits.iter().flatten().copied())?;
    let mut per_cell: FxHashMap<u32, FxHashSet<u32>> = FxHashMap::default();
    for chunk_hits in &hits {
        for &cls in chunk_hits {
            per_cell.entry(la.cell_of(cls)?).or_default().insert(cls);
        }
    }
    Ok((total_support, posts.len(), per_cell))
}

#[derive(Parser)]
pub struct CohortArgs {
    #[command(subcommand)]
    pub what: CohortWhat,
}

#[derive(clap::Subcommand)]
pub enum CohortWhat {
    /// Discover and reduce the same coordinate-defined splice events across named archives.
    Events {
        /// Genomic window containing all event component junctions.
        locus: String,
        /// Named archive as ID=PATH; repeat at least twice.
        #[arg(long = "sample", required = true)]
        samples: Vec<String>,
        /// Optional named barcode/group map as ID=PATH; repeat for grouped samples.
        #[arg(long = "groups")]
        groups: Vec<String>,
        /// Event type; repeat to select several. The default is all supported types.
        #[arg(long = "event-type", value_enum)]
        event_types: Vec<EventTypeArg>,
        /// Minimum catalogue supporting-child count for every component in a sample.
        #[arg(long, default_value_t = 2)]
        min_support: u64,
        /// Minimum number of sample catalogues containing every event component.
        #[arg(long, default_value_t = 2)]
        min_samples: usize,
        /// Omit events with fewer informative molecules summed across samples.
        #[arg(long, default_value_t = 1)]
        min_informative: usize,
        /// Require every emitted group row (or bulk row) in every sample to have at least this
        /// many conservative informative molecules. Zero disables this row-level filter.
        #[arg(long, default_value_t = 0)]
        min_row_informative: usize,
        /// Hard union-catalogue limit; exceeding it fails without truncation.
        #[arg(long, default_value_t = 100_000)]
        max_events: usize,
        /// Optional GTF or compiled AIC used only for event labels.
        #[arg(long)]
        gtf: Option<PathBuf>,
        /// Emit a long sample/group TSV table.
        #[arg(long, conflicts_with = "json")]
        tsv: bool,
        /// Emit one versioned JSON object.
        #[arg(long, conflicts_with = "tsv")]
        json: bool,
        /// Stream a versioned sparse, zero-reconstructible table bundle into a new directory.
        #[arg(long, value_name = "DIR", conflicts_with_all = ["tsv", "json"])]
        sparse_dir: Option<PathBuf>,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
    },
    /// Exact per-sample molecular splice graphs with an optional replicate-aware contrast.
    SpliceGraph {
        /// Genomic window containing the common graph.
        locus: String,
        /// Strict sample<TAB>condition<TAB>archive<TAB>cells design TSV.
        #[arg(long)]
        design: PathBuf,
        /// Ordered CONDITION_A:CONDITION_B contrast; the reported effect is B minus A.
        #[arg(long)]
        contrast: Option<String>,
        /// Emit exact sample rows without fitting a contrast.
        #[arg(long, conflicts_with = "contrast")]
        counts_only: bool,
        /// Minimum local catalogue support for an edge to count toward recurrence.
        #[arg(long, default_value_t = 1)]
        min_support: u64,
        /// Minimum sample catalogues meeting --min-support for a common edge.
        #[arg(long, default_value_t = 2)]
        min_edge_samples: usize,
        /// Minimum same-strand path-fragment UMIs for a sample to enter inference.
        #[arg(long, default_value_t = 10)]
        min_sample_umis: usize,
        /// Minimum eligible biological samples in each contrasted condition.
        #[arg(long, default_value_t = 2)]
        min_replicates: usize,
        /// Minimum path UMIs across eligible samples before testing.
        #[arg(long, default_value_t = 5)]
        min_path_umis: usize,
        /// Minimum eligible samples with nonzero support for a path before testing.
        #[arg(long, default_value_t = 2)]
        min_path_samples: usize,
        /// Hard union-path limit. Exceeding it fails; paths are never truncated.
        #[arg(long, default_value_t = 100_000)]
        max_paths: usize,
        /// Emit one versioned JSON object.
        #[arg(long)]
        json: bool,
        #[command(flatten)]
        uniform_output: UniformQueryOutputArgs,
    },
    /// Donor-aware transcript-end atlas: recurrent 3'-endpoint sites, exact group counts,
    /// reference-supported confidence tiers, and paired biological-sample inference.
    TranscriptEnds(crate::endcmd::Args),
    /// Protocol-aware 3'-tag analysis: learn the fragment-distance kernel, deconvolve
    /// externally catalogued poly(A) sites, and test paired biological samples.
    PolyasiteMixture(crate::endcmd::MixtureArgs),
}

fn parse_named_paths(values: &[String], label: &str) -> Result<Vec<(String, PathBuf)>> {
    let mut seen = BTreeSet::new();
    values
        .iter()
        .map(|value| {
            let (id, path) = value
                .split_once('=')
                .with_context(|| format!("--{label} must be ID=PATH, got {value}"))?;
            if id.is_empty()
                || !id
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || b"_.-".contains(&byte))
            {
                bail!(
                    "--{label} ID {id:?} must contain only ASCII letters, digits, '.', '_', or '-'"
                );
            }
            if path.is_empty() {
                bail!("--{label} {id} has an empty path");
            }
            if !seen.insert(id.to_owned()) {
                bail!("duplicate --{label} ID {id}");
            }
            Ok((id.to_owned(), PathBuf::from(path)))
        })
        .collect()
}

#[derive(Clone)]
struct GraphDesignRow {
    sample: String,
    condition: String,
    archive: PathBuf,
    archive_label: String,
    cells: Option<PathBuf>,
    cells_label: String,
}

fn valid_cohort_identifier(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"_.-".contains(&byte))
}

fn resolve_design_path(base: &std::path::Path, value: &str) -> PathBuf {
    let path = PathBuf::from(value);
    if path.is_absolute() {
        path
    } else {
        base.join(path)
    }
}

fn parse_graph_design_text(path: &std::path::Path, text: &str) -> Result<Vec<GraphDesignRow>> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .map(|line| line.trim_end_matches('\r'))
        .context("graph design is empty")?;
    if header != "sample\tcondition\tarchive\tcells" {
        bail!("graph design header must be exactly: sample<TAB>condition<TAB>archive<TAB>cells");
    }
    let base = path.parent().unwrap_or_else(|| std::path::Path::new("."));
    let mut sample_ids = BTreeSet::new();
    let mut archive_paths = BTreeSet::new();
    let mut rows = Vec::new();
    for (line_index, raw) in lines.enumerate() {
        let line_no = line_index + 2;
        let line = raw.trim_end_matches('\r');
        if line.is_empty() {
            bail!("graph design line {line_no} is empty");
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 4 || fields.iter().any(|field| field.is_empty()) {
            bail!("graph design line {line_no} must have four nonempty tab-separated fields");
        }
        if !valid_cohort_identifier(fields[0]) {
            bail!("graph design line {line_no} has invalid sample ID {:?}", fields[0]);
        }
        if !valid_cohort_identifier(fields[1]) {
            bail!(
                "graph design line {line_no} has invalid condition ID {:?}",
                fields[1]
            );
        }
        if !sample_ids.insert(fields[0].to_owned()) {
            bail!("graph design contains duplicate sample ID {}", fields[0]);
        }
        let archive = resolve_design_path(base, fields[2]);
        let canonical_archive = std::fs::canonicalize(&archive).with_context(|| {
            format!(
                "resolving graph design archive on line {line_no}: {}",
                archive.display()
            )
        })?;
        if !archive_paths.insert(canonical_archive.clone()) {
            bail!(
                "graph design reuses resolved archive {} as more than one biological sample",
                canonical_archive.display()
            );
        }
        let cells = (fields[3] != ".").then(|| resolve_design_path(base, fields[3]));
        rows.push(GraphDesignRow {
            sample: fields[0].to_owned(),
            condition: fields[1].to_owned(),
            archive: canonical_archive,
            archive_label: fields[2].to_owned(),
            cells,
            cells_label: fields[3].to_owned(),
        });
    }
    if rows.len() < 2 {
        bail!("graph design requires at least two biological samples");
    }
    Ok(rows)
}

fn parse_graph_design(path: &std::path::Path) -> Result<Vec<GraphDesignRow>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading graph design {}", path.display()))?;
    parse_graph_design_text(path, &text)
}

fn parse_graph_contrast(value: &str) -> Result<(String, String)> {
    let fields: Vec<&str> = value.split(':').collect();
    if fields.len() != 2
        || !valid_cohort_identifier(fields[0])
        || !valid_cohort_identifier(fields[1])
        || fields[0] == fields[1]
    {
        bail!("--contrast must be two distinct condition IDs formatted CONDITION_A:CONDITION_B");
    }
    Ok((fields[0].to_owned(), fields[1].to_owned()))
}

struct CohortSampleWork {
    caller_index: usize,
    id: String,
    archive: PathBuf,
    archive_data: LazyArchive,
    chunks: Vec<crate::archivecmd::ChunkInfo>,
    chrom_names: Vec<String>,
    metadata: Vec<JunctionMeta>,
    catalogue: BTreeSet<EventKey>,
    scope: QueryScope,
    archive_bytes: u64,
}

struct CohortSampleResult {
    caller_index: usize,
    id: String,
    archive: PathBuf,
    archive_version: u32,
    archive_root: Option<String>,
    definitions: Vec<EventDefinition>,
    results: Vec<EventResult>,
    scope: QueryScope,
    unique_chunk_decodes: usize,
    independent_chunk_decodes: usize,
}

fn event_passes_row_gate(
    result: &EventResult,
    scope: &QueryScope,
    min_row_informative: usize,
) -> bool {
    if min_row_informative == 0 {
        return true;
    }
    match scope.aggregation {
        QueryAggregation::Group => result
            .groups
            .iter()
            .all(|(counts, _)| counts.informative() >= min_row_informative),
        QueryAggregation::Cell => result
            .cells
            .iter()
            .all(|(_, counts)| counts.informative() >= min_row_informative),
        QueryAggregation::Bulk => result.totals.informative() >= min_row_informative,
    }
}

fn row_gate_sort_key(sample: &CohortSampleWork) -> (usize, usize, u64, usize) {
    let shallowest_row = match sample.scope.aggregation {
        QueryAggregation::Group => sample
            .scope
            .selected_per_group
            .iter()
            .copied()
            .min()
            .unwrap_or(0),
        QueryAggregation::Cell | QueryAggregation::Bulk => sample.scope.selected_cells,
    };
    (
        shallowest_row,
        sample.scope.selected_cells,
        sample.archive_bytes,
        sample.caller_index,
    )
}

fn retain_by_mask<T>(values: Vec<T>, keep: &[bool]) -> Vec<T> {
    debug_assert_eq!(values.len(), keep.len());
    values
        .into_iter()
        .zip(keep)
        .filter_map(|(value, &keep)| keep.then_some(value))
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn run_cohort_events(
    locus: &str,
    sample_values: &[String],
    group_values: &[String],
    event_types: &[EventTypeArg],
    min_support: u64,
    min_samples: usize,
    min_informative: usize,
    min_row_informative: usize,
    max_events: usize,
    gtf: Option<&std::path::Path>,
    tsv: bool,
    json_output: bool,
    sparse_dir: Option<&Path>,
    uniform_output: &UniformQueryOutputArgs,
) -> Result<()> {
    validate_uniform_output_flags(uniform_output, tsv, json_output)?;
    if uniform_output.format.is_some() && sparse_dir.is_some() {
        bail!("--format cannot be combined with --sparse-dir; the sparse directory is a separate multi-file artifact");
    }
    let t0 = std::time::Instant::now();
    let samples = parse_named_paths(sample_values, "sample")?;
    if samples.len() < 2 {
        bail!("cohort events requires at least two --sample ID=PATH arguments");
    }
    if min_samples == 0 || min_samples > samples.len() {
        bail!(
            "--min-samples must be between 1 and the {} supplied samples",
            samples.len()
        );
    }
    if max_events == 0 {
        bail!("--max-events must be at least 1");
    }
    let groups = parse_named_paths(group_values, "groups")?;
    let sample_ids: BTreeSet<&str> = samples.iter().map(|(id, _)| id.as_str()).collect();
    for (id, _) in &groups {
        if !sample_ids.contains(id.as_str()) {
            bail!("--groups ID {id} has no matching --sample");
        }
    }
    let group_of: BTreeMap<String, PathBuf> = groups.into_iter().collect();
    let (chrom, start, end) = parse_locus(locus)?;
    let mut work = Vec::with_capacity(samples.len());
    let mut reference_digest: Option<String> = None;
    for (caller_index, (id, archive)) in samples.into_iter().enumerate() {
        let mut archive_data = LazyArchive::open(&archive)
            .with_context(|| format!("opening cohort sample {id} at {}", archive.display()))?;
        if let Some(signature) = &archive_data.genome_sig {
            if let Some(expected) = &reference_digest {
                if expected != &signature.digest {
                    bail!(
                        "cohort samples use incompatible genome signatures: {id} has {} but the first stamped sample has {expected}",
                        signature.digest
                    );
                }
            } else {
                reference_digest = Some(signature.digest.clone());
            }
        } else if reference_digest.is_some() {
            bail!("cohort sample {id} lacks the genome signature present on earlier samples");
        }
        let chunks = read_chunk_index(archive_data.reader())?;
        let chrom_names = archive_data.chrom_names.clone();
        let chrom_id = chrom_names
            .iter()
            .position(|name| name == &chrom)
            .map(|index| index as u32)
            .with_context(|| format!("cohort sample {id} lacks chromosome {chrom}"))?;
        let metadata = read_junction_metadata(&mut archive_data)?;
        let catalogue: BTreeSet<EventKey> = discover_event_keys(
            &chrom,
            chrom_id,
            start,
            end,
            &metadata,
            event_types,
            min_support,
            max_events,
        )?
        .into_iter()
        .collect();
        let group_path = group_of.get(&id).cloned();
        let scope_args = QueryScopeArgs {
            cells: None,
            groups: group_path.clone(),
            agg: if group_path.is_some() {
                QueryAggregationArg::Group
            } else {
                QueryAggregationArg::Bulk
            },
        };
        let mut scope = load_query_scope(&mut archive_data, &scope_args)?;
        if uniform_output.format.is_some() {
            scope.ensure_resolved_mapping_digest()?;
        }
        let archive_bytes = std::fs::metadata(&archive).map_or(0, |metadata| metadata.len());
        work.push(CohortSampleWork {
            caller_index,
            id: id.clone(),
            archive,
            archive_data,
            chunks,
            chrom_names,
            metadata,
            catalogue,
            scope,
            archive_bytes,
        });
    }
    if work.iter().any(|sample| sample.archive_data.genome_sig.is_none())
        && work.iter().any(|sample| sample.archive_data.genome_sig.is_some())
    {
        bail!("cohort samples mix stamped and unstamped genome identities");
    }

    let mut occurrences: BTreeMap<EventKey, usize> = BTreeMap::new();
    for sample in &work {
        for key in &sample.catalogue {
            *occurrences.entry(key.clone()).or_default() += 1;
        }
    }
    let mut keys: Vec<EventKey> = occurrences
        .into_iter()
        .filter_map(|(key, count)| (count >= min_samples).then_some(key))
        .collect();
    if keys.len() > max_events {
        bail!(
            "cohort union contains {} events, exceeding --max-events {max_events}",
            keys.len()
        );
    }
    let pre_row_candidate_events = keys.len();
    if min_row_informative > 0 {
        work.sort_by_key(row_gate_sort_key);
    }
    let row_gate_evaluation_order: Vec<String> =
        work.iter().map(|sample| sample.id.clone()).collect();

    let mut sample_results: Vec<CohortSampleResult> = Vec::with_capacity(work.len());
    for mut sample in work {
        let mut definitions =
            prepare_event_definitions(&keys, &sample.chrom_names, &sample.metadata)?;
        for definition in &mut definitions {
            definition.catalogue_present = sample.catalogue.contains(&definition.key);
        }
        let (packed_hits, unique_chunk_decodes, independent_chunk_decodes) = event_packed_hits(
            &mut sample.archive_data,
            &sample.chunks,
            &sample.metadata,
            &definitions,
        )?;
        let mut results = reduce_packed_event_results(
            &mut sample.archive_data,
            &packed_hits,
            definitions.len(),
            &sample.scope,
        )?;
        if min_row_informative > 0 {
            let keep: Vec<bool> = definitions
                .iter()
                .zip(&results)
                .map(|(_, result)| {
                    event_passes_row_gate(result, &sample.scope, min_row_informative)
                })
                .collect();
            keys = retain_by_mask(keys, &keep);
            definitions = retain_by_mask(definitions, &keep);
            results = retain_by_mask(results, &keep);
            for previous in &mut sample_results {
                previous.definitions =
                    retain_by_mask(std::mem::take(&mut previous.definitions), &keep);
                previous.results = retain_by_mask(std::mem::take(&mut previous.results), &keep);
            }
        }
        let archive_version = sample.archive_data.reader().archive_version();
        let archive_root = sample
            .archive_data
            .reader()
            .content_commitment()
            .map(|commitment| commitment.to_hex());
        sample_results.push(CohortSampleResult {
            caller_index: sample.caller_index,
            id: sample.id,
            archive: sample.archive,
            archive_version,
            archive_root,
            definitions,
            results,
            scope: sample.scope,
            unique_chunk_decodes,
            independent_chunk_decodes,
        });
    }
    let post_row_candidate_events = keys.len();
    sample_results.sort_by_key(|sample| sample.caller_index);
    let retained: Vec<usize> = (0..keys.len())
        .filter(|&event| {
            sample_results
                .iter()
                .map(|sample| sample.results[event].totals.informative())
                .sum::<usize>()
                >= min_informative
        })
        .collect();
    let (annotation, annotation_content_blake3) = match gtf {
        Some(path) => {
            let (annotation, digest) = load_query_annotation(
                path,
                "cohort event annotation",
                uniform_output.format.is_some(),
            )?;
            (Some(annotation), digest)
        }
        None => (None, None),
    };
    let annotation_index = annotation
        .as_ref()
        .map(|annotation| build_event_annotation_index(annotation, &sample_results[0].definitions));
    let mut planning = json!({
        "candidate_events": pre_row_candidate_events,
        "retained_events": retained.len(),
        "samples": sample_results.iter().map(|sample| json!({
            "sample": sample.id,
            "unique_chunk_decodes": sample.unique_chunk_decodes,
            "independent_chunk_decodes": sample.independent_chunk_decodes,
        })).collect::<Vec<_>>(),
    });
    if min_row_informative > 0 {
        let object = planning.as_object_mut().unwrap();
        object.insert(
            "post_row_candidate_events".into(),
            json!(post_row_candidate_events),
        );
        object.insert(
            "row_gate_evaluation_order".into(),
            json!(row_gate_evaluation_order),
        );
    }

    if uniform_output.format.is_some() {
        let summary = json!({
            "coordinates": "0-based half-open junction boundaries",
            "locus": locus,
            "chrom": chrom,
            "start": start,
            "end": end,
            "semantics": {
                "candidate_union": "coordinate-defined per-archive catalogues",
                "missing_event": "present=false; evidence is not imputed",
                "statistics": "descriptive counts and usage only",
                "both_in_usage_denominator": false,
                "group_totals": "total rows equal the exact sum of the emitted group rows",
            },
            "thresholds": {
                "min_archive_supporting_children": min_support,
                "min_catalogue_samples": min_samples,
                "min_cohort_informative_umis": min_informative,
                "min_row_informative_umis": min_row_informative,
                "max_union_events": max_events,
            },
            "reference_digest": reference_digest,
            "samples": sample_results.len(),
            "planning": planning,
        });
        let sample_schema = TableSchema::new(
            "gravlax.cohort.events.samples.v1",
            vec![
                Field::new("sample_index", DataType::UInt64),
                Field::new("sample", DataType::String),
                Field::new("archive", DataType::String),
                Field::new("aggregation", DataType::String),
                Field::new("selected_cells", DataType::UInt64),
                Field::new("scope", DataType::Json),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["sample"])
                .ordered_by([OrderKey::ascending("sample_index")]),
        )?;
        let event_schema = TableSchema::new(
            "gravlax.cohort.events.events.v1",
            vec![
                Field::new("event_id", DataType::String),
                Field::new("event_type", DataType::String),
                Field::new("chrom", DataType::String),
                Field::new("include_only", DataType::UInt64),
                Field::new("exclude_only", DataType::UInt64),
                Field::new("both", DataType::UInt64),
                Field::new("informative_umis", DataType::UInt64),
                Field::new("present_samples", DataType::UInt64),
                Field::new("genes", DataType::Json).nullable(),
                Field::new("strand", DataType::String).nullable(),
                Field::new("fully_annotated", DataType::Boolean).nullable(),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["event_id"]))?;
        let component_schema = TableSchema::new(
            "gravlax.cohort.events.components.v1",
            vec![
                Field::new("event_id", DataType::String),
                Field::new("side", DataType::String),
                Field::new("side_index", DataType::UInt64),
                Field::new("donor", DataType::UInt64),
                Field::new("acceptor", DataType::UInt64),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key([
            "event_id",
            "side",
            "side_index",
        ]))?;
        let count_schema = TableSchema::new(
            "gravlax.cohort.events.counts.v1",
            vec![
                Field::new("event_id", DataType::String),
                Field::new("sample", DataType::String),
                Field::new("present", DataType::Boolean),
                Field::new("aggregation", DataType::String),
                Field::new("entity", DataType::String),
                Field::new("include_only", DataType::UInt64),
                Field::new("exclude_only", DataType::UInt64),
                Field::new("both", DataType::UInt64),
                Field::new("informative_umis", DataType::UInt64),
                Field::new("usage_fraction", DataType::Float64).nullable(),
                Field::new("cells", DataType::UInt64),
                Field::new("selected_cells", DataType::UInt64),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key([
            "event_id",
            "sample",
            "aggregation",
            "entity",
        ]))?;
        let sample_selection = SelectionSummary::complete(sample_results.len() as u64);
        let event_selection = SelectionSummary::complete(retained.len() as u64);
        let component_rows: usize = retained
            .iter()
            .map(|&event| keys[event].includes.len() + keys[event].excludes.len())
            .sum();
        let component_selection = SelectionSummary::complete(component_rows as u64);
        let rows_per_event: usize = sample_results
            .iter()
            .map(|sample| match sample.scope.aggregation {
                QueryAggregation::Group => sample.scope.group_names.len() + 1,
                QueryAggregation::Cell | QueryAggregation::Bulk => 1,
            })
            .sum();
        let count_selection =
            SelectionSummary::complete(retained.len().saturating_mul(rows_per_event) as u64);
        let mut archive_identities = Vec::with_capacity(sample_results.len());
        let mut warnings = Vec::new();
        for sample in &sample_results {
            match &sample.archive_root {
                Some(root) => archive_identities
                    .push(format!("aie-directory-root-v2:{root};sample={}", sample.id)),
                None => {
                    archive_identities.push(format!(
                        "aie-v{}-unrooted:{};sample={}",
                        sample.archive_version,
                        sample.archive.display(),
                        sample.id
                    ));
                    warnings.push(format!(
                        "legacy v{} archive {} has no rooted content commitment",
                        sample.archive_version,
                        sample.archive.display()
                    ));
                }
            }
        }
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "archive_access".into(),
            json!("per-sample junction catalogues and union-event postings"),
        );
        parameters.insert("locus".into(), json!(locus));
        parameters.insert("samples".into(), json!(sample_values));
        parameters.insert(
            "event_types".into(),
            json!(event_types
                .iter()
                .map(|kind| kind.name())
                .collect::<Vec<_>>()),
        );
        parameters.insert("min_archive_supporting_children".into(), json!(min_support));
        parameters.insert("min_catalogue_samples".into(), json!(min_samples));
        parameters.insert("min_cohort_informative_umis".into(), json!(min_informative));
        parameters.insert(
            "min_row_informative_umis".into(),
            json!(min_row_informative),
        );
        parameters.insert("max_union_events".into(), json!(max_events));
        parameters.insert(
            "sample_scopes".into(),
            json!(sample_results
                .iter()
                .map(|sample| json!({
                    "sample": sample.id,
                    "archive_path": sample.archive,
                    "scope": sample.scope.provenance_json(),
                }))
                .collect::<Vec<_>>()),
        );
        if let Some(path) = gtf {
            parameters.insert("annotation_path".into(), json!(path));
            parameters.insert(
                "annotation_content_blake3".into(),
                json!(annotation_content_blake3
                    .as_deref()
                    .context("uniform cohort event annotation is missing its bound content digest")?),
            );
        }
        let context = ResultContext {
            producer: Producer {
                name: "aie".into(),
                version: env!("CARGO_PKG_VERSION").into(),
            },
            provenance: Provenance {
                archives: archive_identities,
                parameters,
                ..Default::default()
            },
            warnings,
        };
        write_uniform_bundle_output(uniform_output, |writer, format| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                "gravlax.cohort.events.result.v1",
                OutputFormat::from(format),
                &context,
                &summary,
            )?;
            bundle.write_table("samples", &sample_schema, Some(&sample_selection), |rows| {
                for (index, sample) in sample_results.iter().enumerate() {
                    let scope = sample.scope.json();
                    rows.write_row_with(|row| {
                        row.uint64(index as u64)?;
                        row.string(&sample.id)?;
                        row.string(&sample.archive.display().to_string())?;
                        row.string(sample.scope.aggregation_name())?;
                        row.uint64(sample.scope.selected_cells as u64)?;
                        row.json(&scope)?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table("events", &event_schema, Some(&event_selection), |rows| {
                for &event in &retained {
                    let key = &keys[event];
                    let mut totals = JunctionSetCounts::default();
                    let mut present_samples = 0_usize;
                    for sample in &sample_results {
                        totals.add(sample.results[event].totals);
                        present_samples += usize::from(sample.definitions[event].catalogue_present);
                    }
                    let annotation_value = annotation_index.as_ref().map(|annotation| {
                        event_annotation_json(annotation, &sample_results[0].definitions[event])
                    });
                    rows.write_row_with(|row| {
                        row.string(&key.id())?;
                        row.string(key.kind.name())?;
                        row.string(&key.chrom)?;
                        row.uint64(totals.include_only as u64)?;
                        row.uint64(totals.exclude_only as u64)?;
                        row.uint64(totals.both as u64)?;
                        row.uint64(totals.informative() as u64)?;
                        row.uint64(present_samples as u64)?;
                        if let Some(annotation) = &annotation_value {
                            row.json(&annotation["genes"])?;
                            if let Some(strand) = annotation["strand"].as_str() {
                                row.string(strand)?;
                            } else {
                                row.null()?;
                            }
                            row.boolean(annotation["fully_annotated"].as_bool().unwrap_or(false))?;
                        } else {
                            row.null()?;
                            row.null()?;
                            row.null()?;
                        }
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table(
                "components",
                &component_schema,
                Some(&component_selection),
                |rows| {
                    for &event in &retained {
                        let key = &keys[event];
                        for (side_index, &(donor, acceptor)) in key.includes.iter().enumerate() {
                            rows.write_row_with(|row| {
                                row.string(&key.id())?;
                                row.string("include")?;
                                row.uint64(side_index as u64)?;
                                row.uint64(donor as u64)?;
                                row.uint64(acceptor as u64)?;
                                Ok(())
                            })?;
                        }
                        for (side_index, &(donor, acceptor)) in key.excludes.iter().enumerate() {
                            rows.write_row_with(|row| {
                                row.string(&key.id())?;
                                row.string("exclude")?;
                                row.uint64(side_index as u64)?;
                                row.uint64(donor as u64)?;
                                row.uint64(acceptor as u64)?;
                                Ok(())
                            })?;
                        }
                    }
                    Ok(())
                },
            )?;
            bundle.write_table("counts", &count_schema, Some(&count_selection), |rows| {
                let mut emit = |event_id: &str,
                                sample: &CohortSampleResult,
                                present: bool,
                                aggregation: &str,
                                entity: &str,
                                counts: &JunctionSetCounts,
                                cells: usize,
                                selected_cells: usize|
                 -> std::result::Result<(), OutputError> {
                    rows.write_row_with(|row| {
                        row.string(event_id)?;
                        row.string(&sample.id)?;
                        row.boolean(present)?;
                        row.string(aggregation)?;
                        row.string(entity)?;
                        row.uint64(counts.include_only as u64)?;
                        row.uint64(counts.exclude_only as u64)?;
                        row.uint64(counts.both as u64)?;
                        row.uint64(counts.informative() as u64)?;
                        if let Some(usage) = counts.usage() {
                            row.float64(usage)?;
                        } else {
                            row.null()?;
                        }
                        row.uint64(cells as u64)?;
                        row.uint64(selected_cells as u64)?;
                        Ok(())
                    })
                };
                for &event in &retained {
                    let event_id = keys[event].id();
                    for sample in &sample_results {
                        let result = &sample.results[event];
                        let present = sample.definitions[event].catalogue_present;
                        match sample.scope.aggregation {
                            QueryAggregation::Group => {
                                emit(
                                    &event_id,
                                    sample,
                                    present,
                                    "total",
                                    "total",
                                    &result.totals,
                                    result.support_cells,
                                    sample.scope.selected_cells,
                                )?;
                                for (group, (counts, cells)) in result.groups.iter().enumerate() {
                                    emit(
                                        &event_id,
                                        sample,
                                        present,
                                        "group",
                                        &sample.scope.group_names[group],
                                        counts,
                                        *cells,
                                        sample.scope.selected_per_group[group],
                                    )?;
                                }
                            }
                            QueryAggregation::Cell | QueryAggregation::Bulk => emit(
                                &event_id,
                                sample,
                                present,
                                "bulk",
                                "bulk",
                                &result.totals,
                                result.support_cells,
                                sample.scope.selected_cells,
                            )?,
                        }
                    }
                }
                Ok(())
            })?;
            bundle.finish()?;
            Ok(())
        })?;
    } else if let Some(out_dir) = sparse_dir {
        let mut writer = crate::sparseout::SparseCohortWriter::create(out_dir)?;
        for &event in &retained {
            let key = &keys[event];
            let event_id = key.id();
            let includes = key
                .includes
                .iter()
                .map(|(donor, acceptor)| format!("{chrom}:{donor}-{acceptor}"))
                .collect::<Vec<_>>()
                .join(",");
            let excludes = key
                .excludes
                .iter()
                .map(|(donor, acceptor)| format!("{chrom}:{donor}-{acceptor}"))
                .collect::<Vec<_>>()
                .join(",");
            let (annotation_genes_json, strand, fully_annotated) =
                if let Some(annotation) = &annotation_index {
                    let value =
                        event_annotation_json(annotation, &sample_results[0].definitions[event]);
                    (
                        serde_json::to_string(&value["genes"]).unwrap(),
                        value["strand"].as_str().unwrap_or("NA").to_owned(),
                        value["fully_annotated"].as_bool().unwrap().to_string(),
                    )
                } else {
                    ("null".to_owned(), "NA".to_owned(), "NA".to_owned())
                };
            writer.event(
                &event_id,
                key.kind.name(),
                &key.chrom,
                &includes,
                &excludes,
                &annotation_genes_json,
                &strand,
                &fully_annotated,
            )?;
            for sample in &sample_results {
                let result = &sample.results[event];
                if sample.definitions[event].catalogue_present {
                    writer.present(&event_id, &sample.id)?;
                }
                if sample.scope.aggregation == QueryAggregation::Group {
                    writer.count(
                        &event_id,
                        &sample.id,
                        "total",
                        "total",
                        result.totals.include_only,
                        result.totals.exclude_only,
                        result.totals.both,
                        result.support_cells,
                        sample.scope.selected_cells,
                    )?;
                    for (group, (counts, cells)) in result.groups.iter().enumerate() {
                        writer.count(
                            &event_id,
                            &sample.id,
                            "group",
                            &sample.scope.group_names[group],
                            counts.include_only,
                            counts.exclude_only,
                            counts.both,
                            *cells,
                            sample.scope.selected_per_group[group],
                        )?;
                    }
                } else {
                    writer.count(
                        &event_id,
                        &sample.id,
                        "bulk",
                        "bulk",
                        result.totals.include_only,
                        result.totals.exclude_only,
                        result.totals.both,
                        result.support_cells,
                        sample.scope.selected_cells,
                    )?;
                }
            }
        }
        let mut metadata = json!({
            "schema": "gravlax.cohort.events.sparse.v1",
            "source_schema": "gravlax.cohort.events.v1",
            "coordinates": "0-based half-open junction boundaries",
            "semantics": {
                "candidate_union": "coordinate-defined per-archive catalogues",
                "missing_presence_row": "catalogue present=false; evidence is not imputed",
                "missing_count_row": "exact logical zero for include_only, exclude_only, both, and cells over the explicit event/sample/scope dimensions",
                "statistics": "descriptive counts and usage only",
                "both_in_usage_denominator": false,
                "derived_fields": {
                    "informative_umis": "include_only + exclude_only",
                    "usage_fraction": "include_only / informative_umis, or null when informative_umis is zero",
                },
            },
            "locus": locus,
            "min_support": min_support,
            "min_samples": min_samples,
            "min_informative": min_informative,
            "dimensions": {
                "events": retained.len(),
                "samples": sample_results.iter().map(|sample| json!({
                    "sample": sample.id,
                    "archive": sample.archive,
                    "scope": sample.scope.json(),
                })).collect::<Vec<_>>(),
            },
            "planning": planning,
        });
        if min_row_informative > 0 {
            metadata
                .as_object_mut()
                .unwrap()
                .insert("min_row_informative".into(), json!(min_row_informative));
        }
        let metadata = writer.finish(metadata)?;
        println!("{}", serde_json::to_string_pretty(&metadata)?);
    } else if tsv {
        println!("event_id\tevent_type\tgenes\tstrand\tfully_annotated\tsample\tarchive\tpresent\taggregation\tgroup\tinclude_only\texclude_only\tboth\tinformative_umis\tusage_fraction\tcells\tselected_cells");
        for &event in &retained {
            let (genes, strand, fully_annotated) = event_annotation_tsv(
                annotation_index.as_ref(),
                &sample_results[0].definitions[event],
            );
            for sample in &sample_results {
                let result = &sample.results[event];
                let present = sample.definitions[event].catalogue_present;
                let emit = |aggregation: &str,
                            group: &str,
                            counts: &JunctionSetCounts,
                            cells: usize,
                            selected_cells: usize| {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        keys[event].id(),
                        keys[event].kind.name(),
                        genes,
                        strand,
                        fully_annotated,
                        sample.id,
                        sample.archive.display(),
                        present,
                        aggregation,
                        group,
                        counts.include_only,
                        counts.exclude_only,
                        counts.both,
                        counts.informative(),
                        counts
                            .usage()
                            .map_or_else(|| "NA".to_owned(), |usage| format!("{usage:.9}")),
                        cells,
                        selected_cells,
                    );
                };
                if sample.scope.aggregation == QueryAggregation::Group {
                    for (group, (counts, cells)) in result.groups.iter().enumerate() {
                        emit(
                            "group",
                            &sample.scope.group_names[group],
                            counts,
                            *cells,
                            sample.scope.selected_per_group[group],
                        );
                    }
                } else {
                    emit(
                        "bulk",
                        "bulk",
                        &result.totals,
                        result.support_cells,
                        sample.scope.selected_cells,
                    );
                }
            }
        }
    } else if json_output {
        let events = retained.iter().map(|&event| {
            let mut value = json!({
                "id": keys[event].id(),
                "event_type": keys[event].kind.name(),
                "chrom": keys[event].chrom,
                "inclusion_junctions": keys[event].includes.iter().map(|(donor, acceptor)| format!("{chrom}:{donor}-{acceptor}")).collect::<Vec<_>>(),
                "exclusion_junctions": keys[event].excludes.iter().map(|(donor, acceptor)| format!("{chrom}:{donor}-{acceptor}")).collect::<Vec<_>>(),
                "sample_rows": sample_results.iter().map(|sample| {
                    let result = &sample.results[event];
                    let present = sample.definitions[event].catalogue_present;
                    let mut row = json!({
                        "sample": sample.id,
                        "archive": sample.archive,
                        "present": present,
                        "totals": result.totals.json(),
                        "cells": result.support_cells,
                        "scope": sample.scope.json(),
                    });
                    if sample.scope.aggregation == QueryAggregation::Group {
                        row.as_object_mut().unwrap().insert("group_rows".into(), json!(result.groups.iter().enumerate().map(|(group, (counts, cells))| {
                            let mut group_row = counts.json();
                            let object = group_row.as_object_mut().unwrap();
                            object.insert("group".into(), json!(sample.scope.group_names[group]));
                            object.insert("cells".into(), json!(cells));
                            object.insert("selected_cells".into(), json!(sample.scope.selected_per_group[group]));
                            group_row
                        }).collect::<Vec<_>>()));
                    }
                    row
                }).collect::<Vec<_>>(),
            });
            if let Some(annotation) = &annotation_index {
                value.as_object_mut().unwrap().insert(
                    "annotation".into(),
                    event_annotation_json(annotation, &sample_results[0].definitions[event]),
                );
            }
            value
        }).collect::<Vec<_>>();
        let mut value = json!({
            "schema": "gravlax.cohort.events.v1",
            "coordinates": "0-based half-open junction boundaries",
            "semantics": {
                "candidate_union": "coordinate-defined per-archive catalogues",
                "missing_event": "present=false; evidence is not imputed",
                "statistics": "descriptive counts and usage only",
                "both_in_usage_denominator": false,
            },
            "locus": locus,
            "min_support": min_support,
            "min_samples": min_samples,
            "min_informative": min_informative,
            "planning": planning,
            "events": events,
        });
        if min_row_informative > 0 {
            value
                .as_object_mut()
                .unwrap()
                .insert("min_row_informative".into(), json!(min_row_informative));
        }
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        println!(
            "cohort events {locus}: {} retained events across {} samples ({:.2}s)",
            retained.len(), sample_results.len(), t0.elapsed().as_secs_f32(),
        );
        for &event in retained.iter().take(20) {
            println!("  {}", keys[event].id());
            for sample in &sample_results {
                let counts = sample.results[event].totals;
                println!(
                    "    {}\t{} informative\tusage {}",
                    sample.id, counts.informative(),
                    counts.usage().map_or_else(|| "NA".to_owned(), |usage| format!("{usage:.4}")),
                );
            }
        }
    }
    eprintln!(
        "cohort event engine: {} events retained across {} samples in {:.2}s",
        retained.len(), sample_results.len(), t0.elapsed().as_secs_f32(),
    );
    Ok(())
}

struct CohortGraphSampleWork {
    row: GraphDesignRow,
    archive_data: LazyArchive,
    chunks: Vec<crate::archivecmd::ChunkInfo>,
    chrom_id: u32,
    metadata: Vec<JunctionMeta>,
    qualifying_edges: BTreeSet<(u32, u32)>,
    scope: QueryScope,
}

struct CohortGraphSampleResult {
    row: GraphDesignRow,
    archive_version: u32,
    archive_root: Option<String>,
    scope: QueryScope,
    paths: BTreeMap<GraphPathKey, GraphAggregate>,
    edges: BTreeMap<GraphEdgeKey, GraphAggregate>,
    strand_umis: [usize; 2],
    scoped_distinct_umi_classes: usize,
    unique_chunk_decodes: usize,
    independent_chunk_decodes: usize,
}

struct PendingGraphTest {
    path_id: usize,
    fit: apastats::BetaBinomialContrast,
    sample_rows: Vec<serde_json::Value>,
    total_umis: usize,
    supporting_samples: usize,
    condition_a_samples: usize,
    condition_b_samples: usize,
    condition_a_umis: u64,
    condition_a_total: u64,
    condition_b_umis: u64,
    condition_b_total: u64,
}

#[allow(clippy::too_many_arguments)]
fn run_cohort_splice_graph(
    locus: &str,
    design_path: &std::path::Path,
    contrast_value: Option<&str>,
    counts_only: bool,
    min_support: u64,
    min_edge_samples: usize,
    min_sample_umis: usize,
    min_replicates: usize,
    min_path_umis: usize,
    min_path_samples: usize,
    max_paths: usize,
    json_output: bool,
    uniform_output: &UniformQueryOutputArgs,
) -> Result<()> {
    validate_uniform_output_flags(uniform_output, false, json_output)?;
    let t0 = std::time::Instant::now();
    if min_support == 0 {
        bail!("--min-support must be at least 1");
    }
    for (name, value) in [
        ("--min-edge-samples", min_edge_samples),
        ("--min-sample-umis", min_sample_umis),
        ("--min-replicates", min_replicates),
        ("--min-path-umis", min_path_umis),
        ("--min-path-samples", min_path_samples),
        ("--max-paths", max_paths),
    ] {
        if value == 0 {
            bail!("{name} must be at least 1");
        }
    }
    if counts_only && contrast_value.is_some() {
        bail!("--counts-only cannot be combined with --contrast");
    }
    if !counts_only && contrast_value.is_none() {
        bail!("replicate-aware splice-graph inference requires --contrast A:B or --counts-only");
    }
    let (design, design_content_blake3) = if uniform_output.format.is_some() {
        let bound = read_bound_query_text(design_path, "cohort graph design")?;
        (
            parse_graph_design_text(design_path, &bound.text)?,
            Some(bound.content_blake3),
        )
    } else {
        (parse_graph_design(design_path)?, None)
    };
    if min_edge_samples > design.len() {
        bail!(
            "--min-edge-samples cannot exceed the {} design samples",
            design.len()
        );
    }
    let contrast = contrast_value.map(parse_graph_contrast).transpose()?;
    if let Some((condition_a, condition_b)) = &contrast {
        let observed: BTreeSet<&str> = design.iter().map(|row| row.condition.as_str()).collect();
        let expected = BTreeSet::from([condition_a.as_str(), condition_b.as_str()]);
        if observed != expected {
            bail!(
                "graph design conditions {:?} do not exactly match contrast {}:{}",
                observed,
                condition_a,
                condition_b
            );
        }
    }
    let (chrom, start, end) = parse_locus(locus)?;

    let mut work = Vec::with_capacity(design.len());
    let mut reference_digest: Option<String> = None;
    for row in design {
        let mut archive_data = LazyArchive::open(&row.archive).with_context(|| {
            format!(
                "opening graph sample {} at {}",
                row.sample,
                row.archive.display()
            )
        })?;
        if let Some(signature) = &archive_data.genome_sig {
            if let Some(expected) = &reference_digest {
                if expected != &signature.digest {
                    bail!(
                        "graph samples use incompatible genome signatures: {} has {} but expected {expected}",
                        row.sample,
                        signature.digest
                    );
                }
            } else {
                reference_digest = Some(signature.digest.clone());
            }
        }
        let chunks = read_chunk_index(archive_data.reader())?;
        let chrom_id = archive_data
            .chrom_names
            .iter()
            .position(|name| name == &chrom)
            .map(|index| index as u32)
            .with_context(|| format!("graph sample {} lacks chromosome {chrom}", row.sample))?;
        let metadata = read_junction_metadata(&mut archive_data)?;
        let qualifying_edges = metadata
            .iter()
            .filter(|junction| {
                junction.chrom == chrom_id
                    && junction.donor >= start
                    && junction.acceptor < end
                    && junction.supporting_children >= min_support
            })
            .map(|junction| (junction.donor, junction.acceptor))
            .collect();
        let scope_args = QueryScopeArgs {
            cells: row.cells.clone(),
            groups: None,
            agg: QueryAggregationArg::Bulk,
        };
        let mut scope = load_query_scope(&mut archive_data, &scope_args)
            .with_context(|| format!("loading cell scope for graph sample {}", row.sample))?;
        if uniform_output.format.is_some() {
            scope.ensure_resolved_mapping_digest()?;
        }
        work.push(CohortGraphSampleWork {
            row,
            archive_data,
            chunks,
            chrom_id,
            metadata,
            qualifying_edges,
            scope,
        });
    }
    if work
        .iter()
        .any(|sample| sample.archive_data.genome_sig.is_some())
        && work
            .iter()
            .any(|sample| sample.archive_data.genome_sig.is_none())
    {
        bail!("graph samples mix stamped and unstamped genome identities");
    }

    let mut edge_occurrences: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for sample in &work {
        for &edge in &sample.qualifying_edges {
            *edge_occurrences.entry(edge).or_default() += 1;
        }
    }
    let common_edges: BTreeSet<(u32, u32)> = edge_occurrences
        .iter()
        .filter_map(|(&edge, &samples)| (samples >= min_edge_samples).then_some(edge))
        .collect();

    let mut union_paths = BTreeSet::new();
    let mut sample_results = Vec::with_capacity(work.len());
    for mut sample in work {
        let selected_metadata: Vec<&JunctionMeta> = sample
            .metadata
            .iter()
            .filter(|junction| {
                junction.chrom == sample.chrom_id
                    && common_edges.contains(&(junction.donor, junction.acceptor))
            })
            .collect();
        let reduction = reduce_graph_paths(
            &mut sample.archive_data,
            &sample.chunks,
            sample.chrom_id,
            &selected_metadata,
            &sample.scope,
            max_paths,
        )?;
        for path in reduction.paths.keys() {
            if !union_paths.contains(path) && union_paths.len() == max_paths {
                bail!(
                    "cohort graph has more than {max_paths} union paths, exceeding --max-paths {max_paths}"
                );
            }
            union_paths.insert(path.clone());
        }
        let edges = aggregate_graph_edges(reduction.paths.iter(), 0);
        let mut strand_umis = [0usize; 2];
        for (path, counts) in &reduction.paths {
            strand_umis[path.strand_rev as usize] += counts.umis;
        }
        let archive_version = sample.archive_data.reader().archive_version();
        let archive_root = sample
            .archive_data
            .reader()
            .content_commitment()
            .map(|commitment| commitment.to_hex());
        sample_results.push(CohortGraphSampleResult {
            row: sample.row,
            archive_version,
            archive_root,
            scope: sample.scope,
            paths: reduction.paths,
            edges,
            strand_umis,
            scoped_distinct_umi_classes: reduction.scoped_distinct_umi_classes,
            unique_chunk_decodes: reduction.unique_chunk_decodes,
            independent_chunk_decodes: reduction.independent_chunk_decodes,
        });
    }
    let paths: Vec<GraphPathKey> = union_paths.into_iter().collect();
    let mut union_edges = BTreeSet::new();
    for path in &paths {
        for &(donor, acceptor) in &path.junctions {
            union_edges.insert(GraphEdgeKey {
                strand_rev: path.strand_rev,
                donor,
                acceptor,
            });
        }
    }
    let edges: Vec<GraphEdgeKey> = union_edges.into_iter().collect();
    let edge_ids: BTreeMap<GraphEdgeKey, usize> = edges
        .iter()
        .copied()
        .enumerate()
        .map(|(id, edge)| (edge, id))
        .collect();
    let mut node_keys = BTreeSet::new();
    for edge in &edges {
        node_keys.insert((edge.strand_rev, edge.donor));
        node_keys.insert((edge.strand_rev, edge.acceptor));
    }
    let node_ids: BTreeMap<(bool, u32), usize> = node_keys
        .into_iter()
        .enumerate()
        .map(|(id, node)| (node, id))
        .collect();

    let nodes_json = if uniform_output.format.is_some() {
        Vec::new()
    } else {
        node_ids
            .iter()
            .map(|(&(strand_rev, coordinate), &id)| {
                json!({
                    "id": id,
                    "coordinate": coordinate,
                    "strand": if strand_rev { "-" } else { "+" },
                })
            })
            .collect::<Vec<_>>()
    };
    let edges_json = if uniform_output.format.is_some() {
        Vec::new()
    } else {
        edges
            .iter()
            .enumerate()
            .map(|(id, edge)| {
            let source_coordinate = if edge.strand_rev {
                edge.acceptor
            } else {
                edge.donor
            };
            let target_coordinate = if edge.strand_rev {
                edge.donor
            } else {
                edge.acceptor
            };
            json!({
                "id": id,
                "strand": if edge.strand_rev { "-" } else { "+" },
                "donor": edge.donor,
                "acceptor": edge.acceptor,
                "source": node_ids[&(edge.strand_rev, source_coordinate)],
                "target": node_ids[&(edge.strand_rev, target_coordinate)],
                "catalogue_samples": edge_occurrences.get(&(edge.donor, edge.acceptor)).copied().unwrap_or(0),
            })
            })
            .collect::<Vec<_>>()
    };
    let paths_json = if uniform_output.format.is_some() {
        Vec::new()
    } else {
        paths
            .iter()
            .enumerate()
            .map(|(id, path)| {
                let ordered = ordered_graph_junctions(path);
                let aggregate_umis: usize = sample_results
                    .iter()
                    .map(|sample| sample.paths.get(path).map_or(0, |counts| counts.umis))
                    .sum();
                let supporting_samples = sample_results
                    .iter()
                    .filter(|sample| sample.paths.get(path).is_some_and(|counts| counts.umis > 0))
                    .count();
                json!({
                    "id": id,
                    "strand": if path.strand_rev { "-" } else { "+" },
                    "edge_ids": ordered.iter().map(|&(donor, acceptor)| edge_ids[&GraphEdgeKey {
                        strand_rev: path.strand_rev,
                        donor,
                        acceptor,
                    }]).collect::<Vec<_>>(),
                    "junctions": ordered.iter().map(|&(donor, acceptor)| json!({
                        "donor": donor,
                        "acceptor": acceptor,
                    })).collect::<Vec<_>>(),
                    "aggregate_umis": aggregate_umis,
                    "supporting_samples": supporting_samples,
                })
            })
            .collect::<Vec<_>>()
    };

    let samples_json = if uniform_output.format.is_some() {
        Vec::new()
    } else {
        sample_results
            .iter()
            .map(|sample| {
                json!({
                    "sample": sample.row.sample,
                    "condition": sample.row.condition,
                    "archive": sample.row.archive_label,
                    "cells": sample.row.cells_label,
                    "scope": sample.scope.json(),
                    "strand_totals": [
                        {
                            "strand": "+",
                            "umis": sample.strand_umis[0],
                            "eligible": sample.strand_umis[0] >= min_sample_umis,
                        },
                        {
                            "strand": "-",
                            "umis": sample.strand_umis[1],
                            "eligible": sample.strand_umis[1] >= min_sample_umis,
                        },
                    ],
                    "path_rows": paths.iter().enumerate().map(|(path_id, path)| {
                        let counts = sample.paths.get(path);
                        json!({
                            "path_id": path_id,
                            "umis": counts.map_or(0, |counts| counts.umis),
                            "cells": counts.map_or(0, |counts| counts.cells.len()),
                        })
                    }).collect::<Vec<_>>(),
                    "edge_rows": edges.iter().enumerate().map(|(edge_id, edge)| {
                        let counts = sample.edges.get(edge);
                        json!({
                            "edge_id": edge_id,
                            "umis": counts.map_or(0, |counts| counts.umis),
                            "cells": counts.map_or(0, |counts| counts.cells.len()),
                        })
                    }).collect::<Vec<_>>(),
                    "planning": {
                        "scoped_distinct_umi_classes": sample.scoped_distinct_umi_classes,
                        "unique_chunk_decodes": sample.unique_chunk_decodes,
                        "independent_chunk_decodes": sample.independent_chunk_decodes,
                    },
                })
            })
            .collect::<Vec<_>>()
    };

    let mut pending_tests = Vec::new();
    let mut skipped_tests = Vec::new();
    if let Some((condition_a, condition_b)) = &contrast {
        for (path_id, path) in paths.iter().enumerate() {
            let strand = path.strand_rev as usize;
            let sample_rows = if uniform_output.format.is_some() {
                Vec::new()
            } else {
                sample_results
                    .iter()
                    .map(|sample| {
                        let path_umis = sample.paths.get(path).map_or(0, |counts| counts.umis);
                        let strand_umis = sample.strand_umis[strand];
                        json!({
                            "sample": sample.row.sample,
                            "condition": sample.row.condition,
                            "path_umis": path_umis,
                            "strand_umis": strand_umis,
                            "eligible": strand_umis >= min_sample_umis,
                        })
                    })
                    .collect::<Vec<_>>()
            };
            let eligible: Vec<(usize, &CohortGraphSampleResult)> = sample_results
                .iter()
                .enumerate()
                .filter(|(_, sample)| sample.strand_umis[strand] >= min_sample_umis)
                .collect();
            let condition_a_rows: Vec<(u64, u64)> = eligible
                .iter()
                .filter(|(_, sample)| sample.row.condition == *condition_a)
                .map(|(_, sample)| {
                    (
                        sample.paths.get(path).map_or(0, |counts| counts.umis) as u64,
                        sample.strand_umis[strand] as u64,
                    )
                })
                .collect();
            let condition_b_rows: Vec<(u64, u64)> = eligible
                .iter()
                .filter(|(_, sample)| sample.row.condition == *condition_b)
                .map(|(_, sample)| {
                    (
                        sample.paths.get(path).map_or(0, |counts| counts.umis) as u64,
                        sample.strand_umis[strand] as u64,
                    )
                })
                .collect();
            let total_umis: usize = eligible
                .iter()
                .map(|(_, sample)| sample.paths.get(path).map_or(0, |counts| counts.umis))
                .sum();
            let supporting_samples = eligible
                .iter()
                .filter(|(_, sample)| sample.paths.get(path).is_some_and(|counts| counts.umis > 0))
                .count();
            let comparator_umis: usize = eligible
                .iter()
                .map(|(_, sample)| {
                    sample.strand_umis[strand]
                        - sample.paths.get(path).map_or(0, |counts| counts.umis)
                })
                .sum();
            let reason = if condition_a_rows.len() < min_replicates
                || condition_b_rows.len() < min_replicates
            {
                Some("insufficient_eligible_replicates")
            } else if total_umis < min_path_umis {
                Some("insufficient_path_umis")
            } else if supporting_samples < min_path_samples {
                Some("insufficient_supporting_samples")
            } else if comparator_umis == 0 {
                Some("no_same_strand_comparator")
            } else {
                None
            };
            if let Some(reason) = reason {
                skipped_tests.push(json!({
                    "path_id": path_id,
                    "reason": reason,
                    "eligible_condition_a": condition_a_rows.len(),
                    "eligible_condition_b": condition_b_rows.len(),
                    "total_path_umis": total_umis,
                    "supporting_samples": supporting_samples,
                }));
                continue;
            }
            let fit = apastats::beta_binomial_contrast(&condition_a_rows, &condition_b_rows)
                .context("eligible beta-binomial path rows unexpectedly failed validation")?;
            pending_tests.push(PendingGraphTest {
                path_id,
                fit,
                sample_rows,
                total_umis,
                supporting_samples,
                condition_a_samples: condition_a_rows.len(),
                condition_b_samples: condition_b_rows.len(),
                condition_a_umis: condition_a_rows.iter().map(|row| row.0).sum(),
                condition_a_total: condition_a_rows.iter().map(|row| row.1).sum(),
                condition_b_umis: condition_b_rows.iter().map(|row| row.0).sum(),
                condition_b_total: condition_b_rows.iter().map(|row| row.1).sum(),
            });
        }
    }
    let q_values = apastats::bh_fdr(
        &pending_tests
            .iter()
            .map(|test| test.fit.p_value)
            .collect::<Vec<_>>(),
    );
    if uniform_output.format.is_some() {
        let conditions: BTreeSet<&str> = sample_results
            .iter()
            .map(|sample| sample.row.condition.as_str())
            .collect();
        let summary = json!({
            "coordinates": "0-based half-open junction boundaries",
            "locus": locus,
            "chrom": chrom,
            "start": start,
            "end": end,
            "semantics": {
                "replicate_unit": "one unique design sample/archive row",
                "cells_and_molecules_are_replicates": false,
                "path_fragment": "selected junction set co-supported by one archive UMI class and strand",
                "complete_transcript_claim": false,
                "missing_count": "explicit zero; evidence is not imputed",
                "ineligible_sample": "reported but excluded from inference, never converted to a zero replicate",
                "archive_catalogue_support_is_strand_combined": true,
                "multimapper_policy": "archived anchor placement; alternatives are not resolved",
            },
            "design": {
                "samples": sample_results.len(),
                "conditions": conditions,
            },
            "contrast": contrast,
            "counts_only": counts_only,
            "reference_digest": reference_digest,
            "thresholds": {
                "min_archive_supporting_children": min_support,
                "min_edge_samples": min_edge_samples,
                "min_sample_umis": min_sample_umis,
                "min_replicates": min_replicates,
                "min_path_umis": min_path_umis,
                "min_path_samples": min_path_samples,
                "max_union_paths": max_paths,
            },
            "planning": {
                "catalogue_union_edges": edge_occurrences.len(),
                "common_coordinate_edges": common_edges.len(),
                "strand_edges": edges.len(),
                "nodes": node_ids.len(),
                "union_paths": paths.len(),
                "tested_paths": pending_tests.len(),
                "skipped_paths": skipped_tests.len(),
            },
            "inference": {
                "model": "path-versus-rest beta-binomial; condition-specific means and one shared alternative concentration",
                "test": "one-degree-of-freedom asymptotic likelihood-ratio test",
                "multiplicity": "Benjamini-Hochberg across tested locus paths",
            },
        });
        let sample_schema = TableSchema::new(
            "gravlax.cohort.splice-graph.samples.v1",
            vec![
                Field::new("sample_index", DataType::UInt64),
                Field::new("sample", DataType::String),
                Field::new("condition", DataType::String),
                Field::new("archive", DataType::String),
                Field::new("cells_scope", DataType::String),
                Field::new("selected_cells", DataType::UInt64),
                Field::new("plus_strand_umis", DataType::UInt64),
                Field::new("minus_strand_umis", DataType::UInt64),
                Field::new("plus_eligible", DataType::Boolean),
                Field::new("minus_eligible", DataType::Boolean),
                Field::new("scoped_distinct_umi_classes", DataType::UInt64),
                Field::new("unique_chunk_decodes", DataType::UInt64),
                Field::new("independent_chunk_decodes", DataType::UInt64),
                Field::new("scope", DataType::Json),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["sample"])
                .ordered_by([OrderKey::ascending("sample_index")]),
        )?;
        let node_schema = TableSchema::new(
            "gravlax.cohort.splice-graph.nodes.v1",
            vec![
                Field::new("node_id", DataType::UInt64),
                Field::new("coordinate", DataType::UInt64),
                Field::new("strand", DataType::String),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["node_id"])
                .ordered_by([OrderKey::ascending("node_id")]),
        )?;
        let edge_schema = TableSchema::new(
            "gravlax.cohort.splice-graph.edges.v1",
            vec![
                Field::new("edge_id", DataType::UInt64),
                Field::new("strand", DataType::String),
                Field::new("donor", DataType::UInt64),
                Field::new("acceptor", DataType::UInt64),
                Field::new("source_node_id", DataType::UInt64),
                Field::new("target_node_id", DataType::UInt64),
                Field::new("catalogue_samples", DataType::UInt64),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["edge_id"])
                .ordered_by([OrderKey::ascending("edge_id")]),
        )?;
        let path_schema = TableSchema::new(
            "gravlax.cohort.splice-graph.paths.v1",
            vec![
                Field::new("path_id", DataType::UInt64),
                Field::new("strand", DataType::String),
                Field::new("edge_ids", DataType::Json),
                Field::new("junctions", DataType::Json),
                Field::new("aggregate_umis", DataType::UInt64),
                Field::new("supporting_samples", DataType::UInt64),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["path_id"])
                .ordered_by([OrderKey::ascending("path_id")]),
        )?;
        let path_count_schema = TableSchema::new(
            "gravlax.cohort.splice-graph.path-counts.v1",
            vec![
                Field::new("sample", DataType::String),
                Field::new("path_id", DataType::UInt64),
                Field::new("umis", DataType::UInt64),
                Field::new("cells", DataType::UInt64),
                Field::new("strand_umis", DataType::UInt64),
                Field::new("eligible", DataType::Boolean),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["sample", "path_id"]))?;
        let edge_count_schema = TableSchema::new(
            "gravlax.cohort.splice-graph.edge-counts.v1",
            vec![
                Field::new("sample", DataType::String),
                Field::new("edge_id", DataType::UInt64),
                Field::new("umis", DataType::UInt64),
                Field::new("cells", DataType::UInt64),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["sample", "edge_id"]))?;
        let test_schema = TableSchema::new(
            "gravlax.cohort.splice-graph.tests.v1",
            vec![
                Field::new("path_id", DataType::UInt64),
                Field::new("strand", DataType::String),
                Field::new("total_path_umis", DataType::UInt64),
                Field::new("supporting_samples", DataType::UInt64),
                Field::new("condition_a", DataType::String),
                Field::new("condition_a_samples", DataType::UInt64),
                Field::new("condition_a_path_umis", DataType::UInt64),
                Field::new("condition_a_strand_umis", DataType::UInt64),
                Field::new("condition_a_observed_usage", DataType::Float64),
                Field::new("condition_a_fitted_usage", DataType::Float64),
                Field::new("condition_b", DataType::String),
                Field::new("condition_b_samples", DataType::UInt64),
                Field::new("condition_b_path_umis", DataType::UInt64),
                Field::new("condition_b_strand_umis", DataType::UInt64),
                Field::new("condition_b_observed_usage", DataType::Float64),
                Field::new("condition_b_fitted_usage", DataType::Float64),
                Field::new("effect_b_minus_a", DataType::Float64),
                Field::new("null_mean", DataType::Float64),
                Field::new("null_concentration", DataType::Float64),
                Field::new("alternative_concentration", DataType::Float64),
                Field::new("null_log_likelihood", DataType::Float64),
                Field::new("alternative_log_likelihood", DataType::Float64),
                Field::new("likelihood_ratio", DataType::Float64),
                Field::new("p_value", DataType::Float64),
                Field::new("q_value", DataType::Float64),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["path_id"]))?;
        let skipped_schema = TableSchema::new(
            "gravlax.cohort.splice-graph.skipped-tests.v1",
            vec![
                Field::new("path_id", DataType::UInt64),
                Field::new("reason", DataType::String),
                Field::new("eligible_condition_a", DataType::UInt64),
                Field::new("eligible_condition_b", DataType::UInt64),
                Field::new("total_path_umis", DataType::UInt64),
                Field::new("supporting_samples", DataType::UInt64),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["path_id"]))?;
        let sample_selection = SelectionSummary::complete(sample_results.len() as u64);
        let node_selection = SelectionSummary::complete(node_ids.len() as u64);
        let edge_selection = SelectionSummary::complete(edges.len() as u64);
        let path_selection = SelectionSummary::complete(paths.len() as u64);
        let path_count_selection =
            SelectionSummary::complete(sample_results.len().saturating_mul(paths.len()) as u64);
        let edge_count_selection =
            SelectionSummary::complete(sample_results.len().saturating_mul(edges.len()) as u64);
        let test_selection = SelectionSummary::complete(pending_tests.len() as u64);
        let skipped_selection = SelectionSummary::complete(skipped_tests.len() as u64);
        let mut archive_identities = Vec::with_capacity(sample_results.len());
        let mut warnings = Vec::new();
        for sample in &sample_results {
            match &sample.archive_root {
                Some(root) => archive_identities.push(format!(
                    "aie-directory-root-v2:{root};sample={}",
                    sample.row.sample
                )),
                None => {
                    archive_identities.push(format!(
                        "aie-v{}-unrooted:{};sample={}",
                        sample.archive_version,
                        sample.row.archive.display(),
                        sample.row.sample
                    ));
                    warnings.push(format!(
                        "legacy v{} archive {} has no rooted content commitment",
                        sample.archive_version,
                        sample.row.archive.display()
                    ));
                }
            }
        }
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "archive_access".into(),
            json!("per-sample junction catalogues and common-edge postings"),
        );
        parameters.insert("locus".into(), json!(locus));
        parameters.insert("design_path".into(), json!(design_path));
        parameters.insert(
            "design_content_blake3".into(),
            json!(design_content_blake3
                .as_deref()
                .context("uniform cohort graph design is missing its bound content digest")?),
        );
        parameters.insert("contrast".into(), json!(contrast_value));
        parameters.insert("counts_only".into(), json!(counts_only));
        parameters.insert("min_archive_supporting_children".into(), json!(min_support));
        parameters.insert("min_edge_samples".into(), json!(min_edge_samples));
        parameters.insert("min_sample_umis".into(), json!(min_sample_umis));
        parameters.insert("min_replicates".into(), json!(min_replicates));
        parameters.insert("min_path_umis".into(), json!(min_path_umis));
        parameters.insert("min_path_samples".into(), json!(min_path_samples));
        parameters.insert("max_union_paths".into(), json!(max_paths));
        parameters.insert(
            "sample_scopes".into(),
            json!(sample_results
                .iter()
                .map(|sample| json!({
                    "sample": sample.row.sample,
                    "archive_path": sample.row.archive,
                    "scope": sample.scope.provenance_json(),
                }))
                .collect::<Vec<_>>()),
        );
        let context = ResultContext {
            producer: Producer {
                name: "aie".into(),
                version: env!("CARGO_PKG_VERSION").into(),
            },
            provenance: Provenance {
                archives: archive_identities,
                parameters,
                ..Default::default()
            },
            warnings,
        };
        write_uniform_bundle_output(uniform_output, |writer, format| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                "gravlax.cohort.splice-graph.result.v1",
                OutputFormat::from(format),
                &context,
                &summary,
            )?;
            bundle.write_table("samples", &sample_schema, Some(&sample_selection), |rows| {
                for (index, sample) in sample_results.iter().enumerate() {
                    let scope = sample.scope.json();
                    rows.write_row_with(|row| {
                        row.uint64(index as u64)?;
                        row.string(&sample.row.sample)?;
                        row.string(&sample.row.condition)?;
                        row.string(&sample.row.archive_label)?;
                        row.string(&sample.row.cells_label)?;
                        row.uint64(sample.scope.selected_cells as u64)?;
                        row.uint64(sample.strand_umis[0] as u64)?;
                        row.uint64(sample.strand_umis[1] as u64)?;
                        row.boolean(sample.strand_umis[0] >= min_sample_umis)?;
                        row.boolean(sample.strand_umis[1] >= min_sample_umis)?;
                        row.uint64(sample.scoped_distinct_umi_classes as u64)?;
                        row.uint64(sample.unique_chunk_decodes as u64)?;
                        row.uint64(sample.independent_chunk_decodes as u64)?;
                        row.json(&scope)?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table("nodes", &node_schema, Some(&node_selection), |rows| {
                for (&(strand_rev, coordinate), &node_id) in &node_ids {
                    rows.write_row_with(|row| {
                        row.uint64(node_id as u64)?;
                        row.uint64(coordinate as u64)?;
                        row.string(if strand_rev { "-" } else { "+" })?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table("edges", &edge_schema, Some(&edge_selection), |rows| {
                for (edge_id, edge) in edges.iter().enumerate() {
                    let source = if edge.strand_rev {
                        edge.acceptor
                    } else {
                        edge.donor
                    };
                    let target = if edge.strand_rev {
                        edge.donor
                    } else {
                        edge.acceptor
                    };
                    rows.write_row_with(|row| {
                        row.uint64(edge_id as u64)?;
                        row.string(if edge.strand_rev { "-" } else { "+" })?;
                        row.uint64(edge.donor as u64)?;
                        row.uint64(edge.acceptor as u64)?;
                        row.uint64(node_ids[&(edge.strand_rev, source)] as u64)?;
                        row.uint64(node_ids[&(edge.strand_rev, target)] as u64)?;
                        row.uint64(
                            edge_occurrences
                                .get(&(edge.donor, edge.acceptor))
                                .copied()
                                .unwrap_or(0) as u64,
                        )?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table("paths", &path_schema, Some(&path_selection), |rows| {
                for (path_id, path) in paths.iter().enumerate() {
                    let ordered = ordered_graph_junctions(path);
                    let edge_values = json!(ordered
                        .iter()
                        .map(|&(donor, acceptor)| {
                            edge_ids[&GraphEdgeKey {
                                strand_rev: path.strand_rev,
                                donor,
                                acceptor,
                            }]
                        })
                        .collect::<Vec<_>>());
                    let junction_values = json!(ordered
                        .iter()
                        .map(|&(donor, acceptor)| json!({
                            "donor": donor,
                            "acceptor": acceptor,
                        }))
                        .collect::<Vec<_>>());
                    let aggregate_umis: usize = sample_results
                        .iter()
                        .map(|sample| sample.paths.get(path).map_or(0, |counts| counts.umis))
                        .sum();
                    let supporting_samples = sample_results
                        .iter()
                        .filter(|sample| {
                            sample.paths.get(path).is_some_and(|counts| counts.umis > 0)
                        })
                        .count();
                    rows.write_row_with(|row| {
                        row.uint64(path_id as u64)?;
                        row.string(if path.strand_rev { "-" } else { "+" })?;
                        row.json(&edge_values)?;
                        row.json(&junction_values)?;
                        row.uint64(aggregate_umis as u64)?;
                        row.uint64(supporting_samples as u64)?;
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.write_table(
                "path_counts",
                &path_count_schema,
                Some(&path_count_selection),
                |rows| {
                    for sample in &sample_results {
                        for (path_id, path) in paths.iter().enumerate() {
                            let counts = sample.paths.get(path);
                            let strand_umis = sample.strand_umis[path.strand_rev as usize];
                            rows.write_row_with(|row| {
                                row.string(&sample.row.sample)?;
                                row.uint64(path_id as u64)?;
                                row.uint64(counts.map_or(0, |counts| counts.umis) as u64)?;
                                row.uint64(counts.map_or(0, |counts| counts.cells.len()) as u64)?;
                                row.uint64(strand_umis as u64)?;
                                row.boolean(strand_umis >= min_sample_umis)?;
                                Ok(())
                            })?;
                        }
                    }
                    Ok(())
                },
            )?;
            bundle.write_table(
                "edge_counts",
                &edge_count_schema,
                Some(&edge_count_selection),
                |rows| {
                    for sample in &sample_results {
                        for (edge_id, edge) in edges.iter().enumerate() {
                            let counts = sample.edges.get(edge);
                            rows.write_row_with(|row| {
                                row.string(&sample.row.sample)?;
                                row.uint64(edge_id as u64)?;
                                row.uint64(counts.map_or(0, |counts| counts.umis) as u64)?;
                                row.uint64(counts.map_or(0, |counts| counts.cells.len()) as u64)?;
                                Ok(())
                            })?;
                        }
                    }
                    Ok(())
                },
            )?;
            if let Some((condition_a, condition_b)) = &contrast {
                bundle.write_table("tests", &test_schema, Some(&test_selection), |rows| {
                    for (test, &q_value) in pending_tests.iter().zip(&q_values) {
                        rows.write_row_with(|row| {
                            row.uint64(test.path_id as u64)?;
                            row.string(if paths[test.path_id].strand_rev {
                                "-"
                            } else {
                                "+"
                            })?;
                            row.uint64(test.total_umis as u64)?;
                            row.uint64(test.supporting_samples as u64)?;
                            row.string(condition_a)?;
                            row.uint64(test.condition_a_samples as u64)?;
                            row.uint64(test.condition_a_umis)?;
                            row.uint64(test.condition_a_total)?;
                            row.float64(
                                test.condition_a_umis as f64 / test.condition_a_total as f64,
                            )?;
                            row.float64(test.fit.condition_a_mean)?;
                            row.string(condition_b)?;
                            row.uint64(test.condition_b_samples as u64)?;
                            row.uint64(test.condition_b_umis)?;
                            row.uint64(test.condition_b_total)?;
                            row.float64(
                                test.condition_b_umis as f64 / test.condition_b_total as f64,
                            )?;
                            row.float64(test.fit.condition_b_mean)?;
                            row.float64(test.fit.condition_b_mean - test.fit.condition_a_mean)?;
                            row.float64(test.fit.null_mean)?;
                            row.float64(test.fit.null_concentration)?;
                            row.float64(test.fit.alternative_concentration)?;
                            row.float64(test.fit.null_log_likelihood)?;
                            row.float64(test.fit.alternative_log_likelihood)?;
                            row.float64(test.fit.likelihood_ratio)?;
                            row.float64(test.fit.p_value)?;
                            row.float64(q_value)?;
                            Ok(())
                        })?;
                    }
                    Ok(())
                })?;
                bundle.write_table(
                    "skipped_tests",
                    &skipped_schema,
                    Some(&skipped_selection),
                    |rows| {
                        for skipped in &skipped_tests {
                            rows.write_row_with(|row| {
                                row.uint64(skipped["path_id"].as_u64().unwrap_or(0))?;
                                row.string(skipped["reason"].as_str().unwrap_or("unknown"))?;
                                row.uint64(skipped["eligible_condition_a"].as_u64().unwrap_or(0))?;
                                row.uint64(skipped["eligible_condition_b"].as_u64().unwrap_or(0))?;
                                row.uint64(skipped["total_path_umis"].as_u64().unwrap_or(0))?;
                                row.uint64(skipped["supporting_samples"].as_u64().unwrap_or(0))?;
                                Ok(())
                            })?;
                        }
                        Ok(())
                    },
                )?;
            }
            bundle.finish()?;
            Ok(())
        })?;
        eprintln!(
            "cohort splice graph: {} paths across {} samples, {} tests in {:.2}s",
            paths.len(),
            sample_results.len(),
            pending_tests.len(),
            t0.elapsed().as_secs_f32(),
        );
        return Ok(());
    }
    let tests_json = pending_tests
        .into_iter()
        .zip(q_values)
        .map(|(test, q_value)| {
            let (condition_a, condition_b) = contrast.as_ref().unwrap();
            json!({
                "path_id": test.path_id,
                "strand": if paths[test.path_id].strand_rev { "-" } else { "+" },
                "sample_rows": test.sample_rows,
                "total_path_umis": test.total_umis,
                "supporting_samples": test.supporting_samples,
                "condition_a": {
                    "condition": condition_a,
                    "eligible_samples": test.condition_a_samples,
                    "path_umis": test.condition_a_umis,
                    "strand_umis": test.condition_a_total,
                    "observed_usage": test.condition_a_umis as f64 / test.condition_a_total as f64,
                    "fitted_usage": test.fit.condition_a_mean,
                },
                "condition_b": {
                    "condition": condition_b,
                    "eligible_samples": test.condition_b_samples,
                    "path_umis": test.condition_b_umis,
                    "strand_umis": test.condition_b_total,
                    "observed_usage": test.condition_b_umis as f64 / test.condition_b_total as f64,
                    "fitted_usage": test.fit.condition_b_mean,
                },
                "effect_b_minus_a": test.fit.condition_b_mean - test.fit.condition_a_mean,
                "beta_binomial": {
                    "null_mean": test.fit.null_mean,
                    "null_concentration": test.fit.null_concentration,
                    "alternative_concentration": test.fit.alternative_concentration,
                    "null_log_likelihood": test.fit.null_log_likelihood,
                    "alternative_log_likelihood": test.fit.alternative_log_likelihood,
                    "likelihood_ratio": test.fit.likelihood_ratio,
                    "degrees_of_freedom": 1,
                    "p_value": test.fit.p_value,
                    "q_value": q_value,
                    "calibration": "asymptotic chi-square likelihood-ratio tail",
                },
            })
        })
        .collect::<Vec<_>>();

    let inference = if let Some((condition_a, condition_b)) = &contrast {
        json!({
            "enabled": true,
            "unit": "biological sample from one unique archive row",
            "contrast": {
                "condition_a": condition_a,
                "condition_b": condition_b,
                "effect": "condition_b fitted usage minus condition_a fitted usage",
            },
            "model": "path-versus-rest beta-binomial; condition-specific means and one shared alternative concentration",
            "test": "one-degree-of-freedom asymptotic likelihood-ratio test",
            "multiplicity": "Benjamini-Hochberg across tested locus paths",
            "tests": tests_json,
            "skipped_tests": skipped_tests,
        })
    } else {
        json!({
            "enabled": false,
            "reason": "counts-only",
            "tests": [],
        })
    };
    let value = json!({
        "schema": "gravlax.cohort.splice-graph.v1",
        "coordinates": "0-based half-open junction boundaries",
        "semantics": {
            "replicate_unit": "one unique design sample/archive row",
            "cells_and_molecules_are_replicates": false,
            "path_fragment": "selected junction set co-supported by one archive UMI class and strand",
            "complete_transcript_claim": false,
            "missing_count": "explicit zero; evidence is not imputed",
            "ineligible_sample": "reported but excluded from inference, never converted to a zero replicate",
            "catalogue_support_is_strand_combined": true,
            "multimapper_policy": "archived anchor placement; alternatives are not resolved",
        },
        "locus": locus,
        "design": {
            "path": design_path,
            "samples": sample_results.len(),
            "conditions": sample_results.iter().map(|sample| sample.row.condition.clone()).collect::<BTreeSet<_>>(),
        },
        "reference_digest": reference_digest,
        "thresholds": {
            "min_support": min_support,
            "min_edge_samples": min_edge_samples,
            "min_sample_umis": min_sample_umis,
            "min_replicates": min_replicates,
            "min_path_umis": min_path_umis,
            "min_path_samples": min_path_samples,
            "max_paths": max_paths,
        },
        "planning": {
            "catalogue_union_edges": edge_occurrences.len(),
            "common_coordinate_edges": common_edges.len(),
            "strand_edges": edges.len(),
            "nodes": node_ids.len(),
            "union_paths": paths.len(),
            "tested_paths": inference["tests"].as_array().map_or(0, Vec::len),
        },
        "nodes": nodes_json,
        "edges": edges_json,
        "paths": paths_json,
        "samples": samples_json,
        "inference": inference,
    });
    if json_output {
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        println!(
            "cohort splice-graph {locus}: {} samples, {} strand edges, {} paths, {} tests ({:.2}s)",
            sample_results.len(),
            edges.len(),
            paths.len(),
            value["planning"]["tested_paths"],
            t0.elapsed().as_secs_f32(),
        );
        for sample in &sample_results {
            println!(
                "  {}\t{}\t{} plus / {} minus strand-path UMIs",
                sample.row.sample,
                sample.row.condition,
                sample.strand_umis[0],
                sample.strand_umis[1]
            );
        }
    }
    eprintln!(
        "cohort splice graph: {} paths across {} samples, {} tests in {:.2}s",
        paths.len(),
        sample_results.len(),
        value["planning"]["tested_paths"],
        t0.elapsed().as_secs_f32(),
    );
    Ok(())
}

pub fn run_cohort(args: CohortArgs) -> Result<()> {
    match args.what {
        CohortWhat::Events {
            locus,
            samples,
            groups,
            event_types,
            min_support,
            min_samples,
            min_informative,
            min_row_informative,
            max_events,
            gtf,
            tsv,
            json,
            sparse_dir,
            uniform_output,
        } => run_cohort_events(
            &locus,
            &samples,
            &groups,
            &event_types,
            min_support,
            min_samples,
            min_informative,
            min_row_informative,
            max_events,
            gtf.as_deref(),
            tsv,
            json,
            sparse_dir.as_deref(),
            &uniform_output,
        ),
        CohortWhat::SpliceGraph {
            locus,
            design,
            contrast,
            counts_only,
            min_support,
            min_edge_samples,
            min_sample_umis,
            min_replicates,
            min_path_umis,
            min_path_samples,
            max_paths,
            json,
            uniform_output,
        } => run_cohort_splice_graph(
            &locus,
            &design,
            contrast.as_deref(),
            counts_only,
            min_support,
            min_edge_samples,
            min_sample_umis,
            min_replicates,
            min_path_umis,
            min_path_samples,
            max_paths,
            json,
            &uniform_output,
        ),
        CohortWhat::TranscriptEnds(args) => crate::endcmd::run(args),
        CohortWhat::PolyasiteMixture(args) => crate::endcmd::run_polyasite_mixture(args),
    }
}

/// `aie federate` — one junction query across N archives: the miniature-atlas access pattern.
/// Per archive: its own barcode namespace, its own catalogue; the merged view is per-sample rows.
#[derive(Parser)]
pub struct FederateArgs {
    /// Two or more .aie archives.
    #[arg(required = true, num_args = 2..)]
    pub archives: Vec<PathBuf>,
    /// junction chrom:donor-acceptor
    pub locus: String,
    #[arg(long, default_value_t = 5)]
    pub top: usize,
    #[command(flatten)]
    uniform_output: UniformQueryOutputArgs,
}

pub fn run_federate(args: FederateArgs) -> Result<()> {
    validate_uniform_output_flags(&args.uniform_output, false, false)?;
    if args.uniform_output.format.is_some() {
        struct ArchiveResult {
            path: PathBuf,
            archive_version: u32,
            archive_root: Option<String>,
            status: &'static str,
            archive_supporting_children: Option<u64>,
            archive_posting_chunks: Option<usize>,
            total_umis: usize,
            cells: Vec<(u32, usize)>,
            dictionary: Vec<u32>,
        }

        let t0 = std::time::Instant::now();
        let (chrom, donor, acceptor) = parse_locus(&args.locus)?;
        let results: Vec<ArchiveResult> = args
            .archives
            .par_iter()
            .map(|path| -> Result<ArchiveResult> {
                let mut archive = LazyArchive::open(path)?;
                let archive_version = archive.reader().archive_version();
                let archive_root = archive
                    .reader()
                    .content_commitment()
                    .map(|commitment| commitment.to_hex());
                let chunks = read_chunk_index(archive.reader())?;
                let Some(chrom_id) = archive.chrom_names.iter().position(|name| *name == chrom)
                else {
                    return Ok(ArchiveResult {
                        path: path.clone(),
                        archive_version,
                        archive_root,
                        status: "chromosome_absent",
                        archive_supporting_children: None,
                        archive_posting_chunks: None,
                        total_umis: 0,
                        cells: Vec::new(),
                        dictionary: Vec::new(),
                    });
                };
                let Some((support, posting_chunks, per_cell)) =
                    junction_counts(&mut archive, &chunks, chrom_id as u32, donor, acceptor)?
                else {
                    return Ok(ArchiveResult {
                        path: path.clone(),
                        archive_version,
                        archive_root,
                        status: "junction_absent",
                        archive_supporting_children: None,
                        archive_posting_chunks: None,
                        total_umis: 0,
                        cells: Vec::new(),
                        dictionary: Vec::new(),
                    });
                };
                let dictionary = archive.cells()?.to_vec();
                let mut cells: Vec<(u32, usize)> = per_cell
                    .iter()
                    .map(|(&cell, classes)| (cell, classes.len()))
                    .collect();
                cells.sort_unstable_by_key(|(cell, umis)| {
                    (std::cmp::Reverse(*umis), dictionary[*cell as usize])
                });
                Ok(ArchiveResult {
                    path: path.clone(),
                    archive_version,
                    archive_root,
                    status: "present",
                    archive_supporting_children: Some(support),
                    archive_posting_chunks: Some(posting_chunks),
                    total_umis: cells.iter().map(|(_, umis)| *umis).sum(),
                    cells,
                    dictionary,
                })
            })
            .collect::<Result<_>>()?;
        let emitted = |result: &ArchiveResult| {
            if args.top == 0 {
                result.cells.len()
            } else {
                args.top.min(result.cells.len())
            }
        };
        let grand_umis: usize = results.iter().map(|result| result.total_umis).sum();
        let grand_cells: usize = results.iter().map(|result| result.cells.len()).sum();
        let count_rows: usize = results.iter().map(emitted).sum();
        let available_count_rows: usize = results.iter().map(|result| result.cells.len()).sum();
        let summary = json!({
            "coordinates": "0-based half-open junction boundaries",
            "chrom": chrom,
            "donor": donor,
            "acceptor": acceptor,
            "archives": results.len(),
            "totals": {
                "umis": grand_umis,
                "cells": grand_cells,
            },
        });
        let archive_schema = TableSchema::new(
            "gravlax.federate.junction.archives.v1",
            vec![
                Field::new("archive_index", DataType::UInt64),
                Field::new("archive", DataType::String),
                Field::new("status", DataType::String),
                Field::new("archive_supporting_children", DataType::UInt64).nullable(),
                Field::new("archive_posting_chunks", DataType::UInt64).nullable(),
                Field::new("umis", DataType::UInt64),
                Field::new("cells", DataType::UInt64),
                Field::new("available_count_rows", DataType::UInt64),
                Field::new("emitted_count_rows", DataType::UInt64),
                Field::new("count_rows_truncated", DataType::Boolean),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["archive_index"])
                .ordered_by([OrderKey::ascending("archive_index")]),
        )?;
        let count_schema = TableSchema::new(
            "gravlax.federate.junction.counts.v1",
            vec![
                Field::new("archive_index", DataType::UInt64),
                Field::new("rank", DataType::UInt64),
                Field::new("barcode", DataType::String),
                Field::new("umis", DataType::UInt64),
            ],
        )?
        .with_semantics(
            TableSemantics::new(RowSemantics::Sequence)
                .with_key(["archive_index", "barcode"])
                .ordered_by([
                    OrderKey::ascending("archive_index"),
                    OrderKey::ascending("rank"),
                ]),
        )?;
        let archive_selection = SelectionSummary::complete(results.len() as u64);
        let count_selection =
            SelectionSummary::selected(available_count_rows as u64, count_rows as u64)?;
        let mut archives = Vec::with_capacity(results.len());
        let mut warnings = Vec::new();
        for (index, result) in results.iter().enumerate() {
            match &result.archive_root {
                Some(root) => archives.push(format!(
                    "aie-directory-root-v2:{root};archive_index={index}"
                )),
                None => {
                    archives.push(format!(
                        "aie-v{}-unrooted:{};archive_index={index}",
                        result.archive_version,
                        result.path.display()
                    ));
                    warnings.push(format!(
                        "legacy v{} archive {} has no rooted content commitment",
                        result.archive_version,
                        result.path.display()
                    ));
                }
            }
        }
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "archive_access".into(),
            json!("per-archive junction catalogue and postings-selected chunks"),
        );
        parameters.insert("locus".into(), json!(args.locus));
        parameters.insert(
            "selection_policy".into(),
            serde_json::to_value(UniformSelectionPolicy {
                requested_top: args.top,
                top_zero_means_all: true,
                comparator: Some("umis descending, barcode ascending, independently per archive"),
            })?,
        );
        let context = ResultContext {
            producer: Producer {
                name: "aie".into(),
                version: env!("CARGO_PKG_VERSION").into(),
            },
            provenance: Provenance {
                archives,
                parameters,
                ..Default::default()
            },
            warnings,
        };
        write_uniform_bundle_output(&args.uniform_output, |writer, format| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                "gravlax.federate.junction.result.v1",
                OutputFormat::from(format),
                &context,
                &summary,
            )?;
            bundle.write_table(
                "archives",
                &archive_schema,
                Some(&archive_selection),
                |rows| {
                    for (index, result) in results.iter().enumerate() {
                        let emitted = emitted(result);
                        rows.write_row_with(|row| {
                            row.uint64(index as u64)?;
                            row.string(&result.path.display().to_string())?;
                            row.string(result.status)?;
                            if let Some(support) = result.archive_supporting_children {
                                row.uint64(support)?;
                            } else {
                                row.null()?;
                            }
                            if let Some(chunks) = result.archive_posting_chunks {
                                row.uint64(chunks as u64)?;
                            } else {
                                row.null()?;
                            }
                            row.uint64(result.total_umis as u64)?;
                            row.uint64(result.cells.len() as u64)?;
                            row.uint64(result.cells.len() as u64)?;
                            row.uint64(emitted as u64)?;
                            row.boolean(emitted < result.cells.len())?;
                            Ok(())
                        })?;
                    }
                    Ok(())
                },
            )?;
            bundle.write_table("counts", &count_schema, Some(&count_selection), |rows| {
                for (archive_index, result) in results.iter().enumerate() {
                    for (rank, (cell, umis)) in
                        result.cells.iter().take(emitted(result)).enumerate()
                    {
                        let barcode = unpack_cell_bytes(result.dictionary[*cell as usize]);
                        let barcode = std::str::from_utf8(&barcode)
                            .expect("packed archive barcode decodes to ASCII");
                        rows.write_row_with(|row| {
                            row.uint64(archive_index as u64)?;
                            row.uint64(rank as u64)?;
                            row.string(barcode)?;
                            row.uint64(*umis as u64)?;
                            Ok(())
                        })?;
                    }
                }
                Ok(())
            })?;
            bundle.finish()?;
            Ok(())
        })?;
        eprintln!(
            "federated {} archives: {grand_umis} UMIs across {grand_cells} cells total (uniform; {:.2}s)",
            args.archives.len(),
            t0.elapsed().as_secs_f32()
        );
        return Ok(());
    }

    let t0 = std::time::Instant::now();
    let (chrom, donor, acceptor) = parse_locus(&args.locus)?;
    // Archives are independent: query them in parallel, then print per-archive results in the
    // original argument order so the output is byte-identical to the old serial loop.
    let per_archive: Vec<(String, usize, usize)> = args
        .archives
        .par_iter()
        .map(|path| -> Result<(String, usize, usize)> {
            let mut la = LazyArchive::open(path)?;
            let chunks = read_chunk_index(la.reader())?;
            let Some(cid) = la.chrom_names.iter().position(|n| *n == chrom) else {
                return Ok((
                    format!("{}: chromosome {chrom} absent", path.display()),
                    0,
                    0,
                ));
            };
            match junction_counts(&mut la, &chunks, cid as u32, donor, acceptor)? {
                None => Ok((format!("{}: junction absent", path.display()), 0, 0)),
                Some((support, _nch, per_cell)) => {
                    let umis: usize = per_cell.values().map(|s| s.len()).sum();
                    let mut cells: Vec<(u32, usize)> =
                        per_cell.iter().map(|(c2, s)| (*c2, s.len())).collect();
                    cells.sort_unstable_by_key(|(c2, n)| (std::cmp::Reverse(*n), *c2));
                    let dict = la.cells()?;
                    let tops: Vec<String> = cells
                        .iter()
                        .take(args.top)
                        .map(|(c2, n)| format!("{}:{}", unpack_cell(dict, *c2), n))
                        .collect();
                    Ok((
                        format!(
                            "{}: {umis} UMIs / {} cells (index support {support}) top [{}]",
                            path.display(),
                            per_cell.len(),
                            tops.join(", ")
                        ),
                        umis,
                        per_cell.len(),
                    ))
                }
            }
        })
        .collect::<Result<_>>()?;
    let mut grand_umis = 0usize;
    let mut grand_cells = 0usize;
    for (line, umis, cells) in &per_archive {
        println!("{line}");
        grand_umis += umis;
        grand_cells += cells;
    }
    println!(
        "federated {} archives: {grand_umis} UMIs across {grand_cells} cells total ({:.2}s)",
        args.archives.len(), t0.elapsed().as_secs_f32()
    );
    Ok(())
}

/// A molecule's 3'-most coordinate: on + the maximum end over its unique reps, on − the minimum
/// start. Multimapper-only molecules use their anchors the same way.
fn three_prime(m: &crate::rows::MolRec, shapes: &[evidence_io::archive::Shape]) -> u32 {
    let ends = m
        .chains
        .iter()
        .flat_map(|c| c.reps.iter())
        .map(|(pos, sh)| {
            (
                *pos,
                *pos + shapes[*sh as usize]
                    .blocks
                    .last()
                    .map(|b| b.0 + b.1)
                    .unwrap_or(0),
            )
        })
        .chain(m.mms.iter().map(|(pos, sh, _, _)| {
            (
                *pos,
                *pos + shapes[*sh as usize]
                    .blocks
                    .last()
                    .map(|b| b.0 + b.1)
                    .unwrap_or(0),
            )
        }));
    if m.strand_rev {
        ends.map(|(s, _)| s).min().unwrap_or(0)
    } else {
        ends.map(|(_, e)| e).max().unwrap_or(0)
    }
}

fn has_junction(
    shape: &evidence_io::archive::Shape,
    pos: u32,
    donor: u32,
    acceptor: u32,
) -> Result<bool> {
    for blocks in shape.blocks.windows(2) {
        let donor_offset = blocks[0].0.checked_add(blocks[0].1)
            .context("junction donor offset overflow")?;
        let observed_donor = pos.checked_add(donor_offset)
            .context("junction genomic donor overflow")?;
        let observed_acceptor = pos.checked_add(blocks[1].0)
            .context("junction genomic acceptor overflow")?;
        if observed_donor == donor && observed_acceptor == acceptor {
            return Ok(true);
        }
    }
    Ok(false)
}

fn unpack_cell(cells: &[u32], id: u32) -> String {
    String::from_utf8(umi::unpack(cells[id as usize], 16)).unwrap_or_default()
}

/// One called 3'-site: cluster extent, usage, per-group counts, the modal cleavage position, and
/// (when sequence was consulted) the internal-priming stats at that position.
pub struct SiteCall {
    pub lo: u32,
    pub hi: u32,
    pub rev: bool,
    pub umis: usize,
    pub cells: usize,
    pub gc: Vec<u64>,
    pub cp: u32,
    pub ip: bool,
    pub a20: u32,
    pub arun: u32,
    /// Grouped cells only: cell id -> UMIs at this site (feeds the permutation test).
    pub cell_counts: FxHashMap<u32, u64>,
}

/// Cluster molecule 3' points into sites (single-linkage within `site_gap`, per strand). The
/// cleavage position is the modal 3' coordinate, ties resolved downstream in transcript
/// orientation; internal-priming stats are computed there when `seq` is present.
pub fn call_sites(
    pts: &mut Vec<(u32, bool, u32, u32)>, // (3' coord, strand_rev, cell, class)
    site_gap: u32,
    cell_group: &FxHashMap<u32, u32>,
    n_groups: usize,
    seq: Option<&[u8]>,
) -> Vec<SiteCall> {
    pts.sort_unstable();
    pts.dedup();
    let mut out: Vec<SiteCall> = Vec::new();
    for rev in [false, true] {
        let sp: Vec<&(u32, bool, u32, u32)> = pts.iter().filter(|p| p.1 == rev).collect();
        let mut i = 0usize;
        while i < sp.len() {
            let mut j = i;
            while j + 1 < sp.len() && sp[j + 1].0 - sp[j].0 <= site_gap {
                j += 1;
            }
            let site = &sp[i..=j];
            let cells: FxHashSet<u32> = site.iter().map(|p| p.2).collect();
            let mut gc = vec![0u64; n_groups];
            let mut cell_counts: FxHashMap<u32, u64> = FxHashMap::default();
            for p in site {
                if let Some(&g) = cell_group.get(&p.2) {
                    gc[g as usize] += 1;
                    *cell_counts.entry(p.2).or_insert(0) += 1;
                }
            }
            // Modal 3' coordinate; ties toward transcript-downstream (larger on +, smaller on -).
            let mut freq: FxHashMap<u32, u32> = FxHashMap::default();
            for p in site {
                *freq.entry(p.0).or_insert(0) += 1;
            }
            let cp = freq
                .iter()
                .map(|(&tp, &n)| (n, if rev { u32::MAX - tp } else { tp }, tp))
                .max()
                .map(|(_, _, tp)| tp)
                .unwrap_or(site[0].0);
            let (a20, arun) = match seq {
                Some(sq) => apastats::ip_stats(sq, cp, rev),
                None => (0, 0),
            };
            out.push(SiteCall {
                lo: site[0].0,
                hi: site[j - i].0,
                rev,
                umis: site.len(),
                cells: cells.len(),
                gc,
                cp,
                ip: seq.is_some() && apastats::is_internal_priming(a20, arun),
                a20,
                arun,
                cell_counts,
            });
            i = j + 1;
        }
    }
    out
}

/// Parse a two-column (barcode, group) TSV against the archive's cell dictionary.
fn load_groups(
    la: &mut LazyArchive,
    path: &std::path::Path,
) -> Result<(Vec<String>, FxHashMap<u32, u32>)> {
    let cells_dict = la.cells()?.to_vec();
    let mut packed_to_id: FxHashMap<u32, u32> = FxHashMap::default();
    for (i, p2) in cells_dict.iter().enumerate() {
        packed_to_id.insert(*p2, i as u32);
    }
    let mut group_names: Vec<String> = Vec::new();
    let mut cell_group: FxHashMap<u32, u32> = FxHashMap::default();
    for line in std::fs::read_to_string(path)?.lines() {
        let Some((bc, grp)) = line.trim().split_once('\t') else {
            continue;
        };
        let Some(pk) = umi::pack(bc.as_bytes()) else {
            continue;
        };
        let Some(&cid) = packed_to_id.get(&pk) else {
            continue;
        };
        let gi = match group_names.iter().position(|g| g == grp) {
            Some(i) => i as u32,
            None => {
                group_names.push(grp.to_string());
                (group_names.len() - 1) as u32
            }
        };
        cell_group.insert(cid, gi);
    }
    Ok((group_names, cell_group))
}

/// Uniform results bind the resolved population, so silently skipping malformed, duplicate, or
/// unknown rows would make the requested scientific scope ambiguous. The legacy APA parser stays
/// permissive to preserve its historical behavior and bytes.
fn load_groups_strict(
    la: &mut LazyArchive,
    path: &Path,
) -> Result<(Vec<String>, FxHashMap<u32, u32>, String)> {
    let scope = load_query_scope(
        la,
        &QueryScopeArgs {
            cells: None,
            groups: Some(path.to_path_buf()),
            agg: QueryAggregationArg::Group,
        },
    )?;
    let source_content_blake3 = scope
        .source_content_blake3
        .context("strict group scope did not retain its source content digest")?;
    Ok((scope.group_names, scope.group_of, source_content_blake3))
}

/// Uniform results may only claim sequence-derived answers when the archive authenticates the
/// supplied FASTA. Legacy APA keeps its historical warning-and-continue behavior.
fn require_uniform_genome_binding(
    la: &LazyArchive,
    genome: Option<&Path>,
    uniform: bool,
    analysis: &str,
) -> Result<()> {
    if uniform && genome.is_some() && la.genome_sig.is_none() {
        bail!(
            "uniform {analysis} output with --genome requires an archive stamped with that genome identity \
             (stamp it with `aie stamp-genome`)"
        );
    }
    Ok(())
}

/// Load and verify one contig for a windowed query. Unsigned archives get a loud warning, not a
/// refusal — `aie stamp-genome` is the fix.
fn load_query_genome(
    la: &LazyArchive,
    genome: &Option<PathBuf>,
    chrom: &str,
) -> Result<Option<Vec<u8>>> {
    let Some(g) = genome else { return Ok(None) };
    match &la.genome_sig {
        Some(sig) => Ok(Some(evidence_io::genome::load_contig(g, chrom, Some(sig))?)),
        None => {
            eprintln!(
                "warning: archive carries no genome signature; --genome cannot be verified \
                 (stamp it with `aie stamp-genome`)"
            );
            Ok(Some(evidence_io::genome::load_contig(g, chrom, None)?))
        }
    }
}

/// Dense per-cell site-count vectors for the permutation test, grouped cells only.
fn per_cell_site_counts(
    kept: &[&SiteCall],
    cell_group: &FxHashMap<u32, u32>,
) -> (Vec<Vec<(usize, u64)>>, Vec<usize>) {
    let mut dense: FxHashMap<u32, usize> = FxHashMap::default();
    let mut labels: Vec<usize> = Vec::new();
    let mut csc: Vec<Vec<(usize, u64)>> = Vec::new();
    for (si, st) in kept.iter().enumerate() {
        for (&cell, &n) in &st.cell_counts {
            let idx = *dense.entry(cell).or_insert_with(|| {
                labels.push(cell_group[&cell] as usize);
                csc.push(Vec::new());
                labels.len() - 1
            });
            csc[idx].push((si, n));
        }
    }
    (csc, labels)
}

struct ApaTestParams {
    gtf: PathBuf,
    groups: PathBuf,
    genome: Option<PathBuf>,
    site_gap: u32,
    min_site_umis: usize,
    min_gene_umis: usize,
    tail_extend: u32,
    permute: usize,
    seed: u64,
    t0: std::time::Instant,
    t_open: f32,
}

/// Genome-wide per-gene differential 3'-site usage. One sequential pass over the chunks collects
/// grouped molecules' 3' points per gene window; sites are called annotation-free per gene; the
/// G-test runs per gene and BH-FDR is applied across all tested genes.
fn apa_test(
    archive: &Path,
    la: &mut LazyArchive,
    chunks: &[crate::archivecmd::ChunkInfo],
    chrom_names: &[String],
    p: ApaTestParams,
    uniform_output: &UniformQueryOutputArgs,
) -> Result<()> {
    require_uniform_genome_binding(
        la,
        p.genome.as_deref(),
        uniform_output.format.is_some(),
        "APA-test",
    )?;
    let (group_names, cell_group, group_source_content_blake3) = if uniform_output.format.is_some() {
        let (names, mapping, digest) = load_groups_strict(la, &p.groups)?;
        (names, mapping, Some(digest))
    } else {
        let (names, mapping) = load_groups(la, &p.groups)?;
        (names, mapping, None)
    };
    if group_names.len() < 2 {
        bail!("--groups must define at least two groups over archive barcodes");
    }
    let (anno, annotation_content_blake3) = load_query_annotation(
        &p.gtf,
        "APA-test annotation",
        uniform_output.format.is_some(),
    )?;
    // Gene windows: span over transcripts, extended tail_extend on the 3' side.
    struct GeneWin {
        gene: u32,
        rev: bool,
        lo: u32,
        hi: u32,
    }
    let mut span: FxHashMap<u32, (u32, bool, u32, u32)> = FxHashMap::default(); // gene -> chrom, rev, lo, hi
    for t in &anno.transcripts {
        let (s, e) = t.span();
        let ent = span.entry(t.gene).or_insert((t.chrom, t.strand_rev, s, e));
        ent.2 = ent.2.min(s);
        ent.3 = ent.3.max(e);
    }
    // Archive chrom id -> gene windows sorted by lo, with a running-max hi for interval stabbing.
    let anno_of: Vec<Option<u32>> =
        chrom_names.iter().map(|n| anno.chrom_ids.get(n).copied()).collect();
    let mut per_chrom: Vec<Vec<GeneWin>> = (0..chrom_names.len()).map(|_| Vec::new()).collect();
    for (gene, (ac, rev, lo, hi)) in &span {
        if let Some(cid) = anno_of.iter().position(|a| *a == Some(*ac)) {
            let (wlo, whi) = if *rev {
                (lo.saturating_sub(p.tail_extend), *hi)
            } else {
                (*lo, hi + p.tail_extend)
            };
            per_chrom[cid].push(GeneWin {
                gene: *gene,
                rev: *rev,
                lo: wlo,
                hi: whi,
            });
        }
    }
    let mut maxhi: Vec<Vec<u32>> = Vec::with_capacity(per_chrom.len());
    for wins in &mut per_chrom {
        wins.sort_unstable_by_key(|w| w.lo);
        let mut run = 0u32;
        maxhi.push(
            wins.iter()
                .map(|w| {
                    run = run.max(w.hi);
                    run
                })
                .collect(),
        );
    }
    // Pass over chunks: grouped molecules' (tp, cell, class) into overlapping gene buckets.
    let shapes = la.shapes()?;
    let mut pts_of: FxHashMap<u32, Vec<(u32, bool, u32, u32)>> = FxHashMap::default();
    for (i, c) in chunks.iter().enumerate() {
        let wins = &per_chrom[c.chrom as usize];
        if wins.is_empty() {
            continue;
        }
        let mh = &maxhi[c.chrom as usize];
        let raw = la.reader().read(&format!("c{i}"))?;
        for m in decode_chunk(&raw, c, None, &la.rans_tables)? {
            let cell = la.cell_of(m.umi_class)?;
            if !cell_group.contains_key(&cell) {
                continue;
            }
            let tp = three_prime(&m, &shapes);
            // Stab wins for lo <= tp < hi, strand-matched.
            let mut k = wins.partition_point(|w| w.lo <= tp);
            while k > 0 {
                k -= 1;
                if mh[k] <= tp {
                    break;
                }
                let w = &wins[k];
                if w.rev == m.strand_rev && tp < w.hi && tp >= w.lo {
                    pts_of.entry(w.gene).or_default().push((tp, m.strand_rev, cell, m.umi_class));
                }
            }
        }
    }
    // Per gene: call sites (with sequence when supplied), test.
    struct Row {
        gene: u32,
        n_sites: usize,
        umis: u64,
        ip_dropped: usize,
        g: f64,
        df: u64,
        pval: f64,
        p_perm: Option<f64>,
    }
    let mut rows: Vec<Row> = Vec::new();
    let process = |gene: u32,
                       pts: &mut Vec<(u32, bool, u32, u32)>,
                       seq: Option<&[u8]>,
                       rows: &mut Vec<Row>| {
        let sites = call_sites(pts, p.site_gap, &cell_group, group_names.len(), seq);
        let ip_dropped = sites.iter().filter(|st| st.ip).count();
        let kept: Vec<&SiteCall> = sites
            .iter()
            .filter(|st| !st.ip && st.gc.iter().sum::<u64>() >= p.min_site_umis as u64)
            .collect();
        let umis: u64 = kept.iter().map(|st| st.gc.iter().sum::<u64>()).sum();
        if kept.len() < 2 || umis < p.min_gene_umis as u64 {
            return;
        }
        let table: Vec<Vec<u64>> = kept.iter().map(|st| st.gc.clone()).collect();
        let (g, df) = apastats::g_statistic(&table);
        if df == 0 {
            return;
        }
        let pval = apastats::chi2_sf(g, df);
        let p_perm = (p.permute > 0).then(|| {
            let (csc, glab) = per_cell_site_counts(&kept, &cell_group);
            apastats::permutation_p(
                &csc,
                &glab,
                kept.len(),
                group_names.len(),
                g,
                p.permute,
                p.seed,
            )
        });
        rows.push(Row {
            gene,
            n_sites: kept.len(),
            umis,
            ip_dropped,
            g,
            df,
            pval,
            p_perm,
        });
    };
    if let Some(fasta) = &p.genome {
        if la.genome_sig.is_none() {
            eprintln!(
                "warning: archive carries no genome signature; --genome cannot be verified \
                 (stamp it with `aie stamp-genome`)"
            );
        }
        let mut genes_by_chrom: FxHashMap<String, Vec<u32>> = FxHashMap::default();
        for (cid2, wins) in per_chrom.iter().enumerate() {
            for w in wins {
                if pts_of.contains_key(&w.gene) {
                    genes_by_chrom
                        .entry(chrom_names[cid2].clone())
                        .or_default()
                        .push(w.gene);
                }
            }
        }
        let sig = la.genome_sig.clone();
        let mut verr: Option<anyhow::Error> = None;
        evidence_io::genome::for_each_contig(fasta, |name, seq| {
            let Some(genes) = genes_by_chrom.remove(name) else {
                return true;
            };
            if let Some(sig) = &sig {
                if let Err(e) = evidence_io::genome::verify_contig(sig, name, seq) {
                    verr = Some(e);
                    return false;
                }
            }
            for gene in genes {
                if let Some(mut pts) = pts_of.remove(&gene) {
                    process(gene, &mut pts, Some(seq), &mut rows);
                }
            }
            true
        })?;
        if let Some(e) = verr {
            return Err(e);
        }
        if let Some(missed) = genes_by_chrom.keys().next() {
            bail!("contig {missed} with tested genes not found in {}", fasta.display());
        }
    } else {
        let genes: Vec<u32> = pts_of.keys().copied().collect();
        for gene in genes {
            let mut pts = pts_of.remove(&gene).unwrap();
            process(gene, &mut pts, None, &mut rows);
        }
    }
    let qs = apastats::bh_fdr(&rows.iter().map(|r| r.pval).collect::<Vec<_>>());
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        rows[a]
            .pval
            .partial_cmp(&rows[b].pval)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    if uniform_output.format.is_some() {
        let summary = json!({
            "coordinates": "0-based half-open gene windows and three-prime coordinates",
            "semantics": {
                "candidate_sites": "annotation-free calls within annotation-defined gene windows",
                "test": "multinomial site-by-group G-test per gene",
                "multiple_testing": "Benjamini-Hochberg over tested genes",
                "count_unit": "archive UMI class assigned to a supplied group",
            },
            "thresholds": {
                "site_gap": p.site_gap,
                "min_site_umis": p.min_site_umis,
                "min_gene_umis": p.min_gene_umis,
                "tail_extend": p.tail_extend,
            },
            "permutations": p.permute,
            "seed": p.seed,
            "groups": group_names,
            "genes_tested": rows.len(),
            "reference_consulted": p.genome.is_some(),
        });
        let schema = TableSchema::new(
            "gravlax.query.apa-test.genes.v1",
            vec![
                Field::new("gene_id", DataType::String),
                Field::new("gene_name", DataType::String),
                Field::new("sites", DataType::UInt64),
                Field::new("umis", DataType::UInt64),
                Field::new("g_statistic", DataType::Float64),
                Field::new("degrees_of_freedom", DataType::UInt64),
                Field::new("p_value", DataType::Float64),
                Field::new("q_value", DataType::Float64),
                Field::new("permutation_p_value", DataType::Float64).nullable(),
                Field::new("internal_priming_sites_dropped", DataType::UInt64).nullable(),
            ],
        )?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["gene_id"]))?;
        let selection = SelectionSummary::complete(rows.len() as u64);
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "archive_access".into(),
            json!("single sequential archive scan"),
        );
        parameters.insert("annotation_path".into(), json!(p.gtf));
        parameters.insert(
            "annotation_content_blake3".into(),
            json!(annotation_content_blake3
                .as_deref()
                .context("uniform APA-test annotation is missing its bound content digest")?),
        );
        parameters.insert(
            "cell_scope".into(),
            uniform_group_scope_provenance(
                &p.groups,
                group_source_content_blake3
                    .as_deref()
                    .context("uniform APA-test group scope is missing its bound content digest")?,
                &group_names,
                &cell_group,
                la.cells()?.len(),
            )?,
        );
        parameters.insert("site_gap".into(), json!(p.site_gap));
        parameters.insert("min_site_umis".into(), json!(p.min_site_umis));
        parameters.insert("min_gene_umis".into(), json!(p.min_gene_umis));
        parameters.insert("tail_extend".into(), json!(p.tail_extend));
        parameters.insert("permutations".into(), json!(p.permute));
        parameters.insert("seed".into(), json!(p.seed));
        parameters.insert(
            "presentation_order".into(),
            json!("p-value ascending; ties are unspecified and are not a selection rule"),
        );
        if let Some(path) = p.genome.as_deref() {
            parameters.insert("reference_path".into(), json!(path));
            parameters.insert(
                "archive_genome_signature".into(),
                serde_json::to_value(&la.genome_sig)?,
            );
        }
        let context = uniform_query_context(
            archive,
            la.reader().archive_version(),
            la.reader()
                .content_commitment()
                .map(|commitment| commitment.to_hex()),
            parameters,
        );
        write_uniform_bundle_output(uniform_output, |writer, format| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                "gravlax.query.apa-test.result.v1",
                OutputFormat::from(format),
                &context,
                &summary,
            )?;
            bundle.write_table("genes", &schema, Some(&selection), |table| {
                for &index in &order {
                    let result = &rows[index];
                    table.write_row_with(|row| {
                        row.string(&anno.gene_ids[result.gene as usize])?;
                        row.string(&anno.gene_names[result.gene as usize])?;
                        row.uint64(result.n_sites as u64)?;
                        row.uint64(result.umis)?;
                        row.float64(result.g)?;
                        row.uint64(result.df)?;
                        row.float64(result.pval)?;
                        row.float64(qs[index])?;
                        if let Some(value) = result.p_perm {
                            row.float64(value)?;
                        } else {
                            row.null()?;
                        }
                        if p.genome.is_some() {
                            row.uint64(result.ip_dropped as u64)?;
                        } else {
                            row.null()?;
                        }
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
            bundle.finish()?;
            Ok(())
        })?;
    } else {
        let permhdr = if p.permute > 0 { "\tp_perm" } else { "" };
        let iphdr = if p.genome.is_some() {
            "\tip_dropped"
        } else {
            ""
        };
        println!("#gene_id\tgene_name\tn_sites\tumis\tG\tdf\tp\tq{permhdr}{iphdr}");
        for &i in &order {
            let r = &rows[i];
            let perm = r.p_perm.map(|v| format!("\t{v:.5}")).unwrap_or_default();
            let ip = if p.genome.is_some() {
                format!("\t{}", r.ip_dropped)
            } else {
                String::new()
            };
            println!(
                "{}\t{}\t{}\t{}\t{:.3}\t{}\t{:.4e}\t{:.4e}{perm}{ip}",
                anno.gene_ids[r.gene as usize],
                anno.gene_names[r.gene as usize],
                r.n_sites,
                r.umis,
                r.g,
                r.df,
                r.pval,
                qs[i]
            );
        }
    }
    eprintln!(
        "apa-test: {} genes tested across {} groups (open {:.2}s, total {:.2}s)",
        rows.len(),
        group_names.len(),
        p.t_open,
        p.t0.elapsed().as_secs_f32()
    );
    Ok(())
}

/// Gene spans overlapping the window, for the --plot underlay.
fn gene_underlay(
    gtf: &Path,
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<Vec<crate::plots::GeneBox>> {
    let anno = anno::Annotation::from_path(gtf)?;
    let Some(&ac) = anno.chrom_ids.get(chrom) else {
        return Ok(Vec::new());
    };
    let mut span: FxHashMap<u32, (bool, u32, u32)> = FxHashMap::default();
    for t in &anno.transcripts {
        if t.chrom != ac {
            continue;
        }
        let (s, e) = t.span();
        let ent = span.entry(t.gene).or_insert((t.strand_rev, s, e));
        ent.1 = ent.1.min(s);
        ent.2 = ent.2.max(e);
    }
    let mut out: Vec<crate::plots::GeneBox> = span
        .into_iter()
        .filter(|(_, (_, s, e))| *s < end && *e > start)
        .map(|(g, (rev, s, e))| crate::plots::GeneBox {
            name: anno.gene_names[g as usize].clone(),
            start: s as f64,
            end: e as f64,
            rev,
        })
        .collect();
    out.sort_by(|a, b| a.start.partial_cmp(&b.start).unwrap());
    out.truncate(12);
    Ok(out)
}

/// IGV-ready exports: per-strand molecule-coverage bedGraph + TopHat-style BED12 junctions.
fn export_igv(
    prefix: &Path,
    chrom: &str,
    start: u32,
    cov: &[Vec<u32>],
    juncs: &[Vec<(u32, u32, u64)>],
) -> Result<()> {
    use std::io::Write as _;
    for (si, name) in [(0usize, "plus"), (1, "minus")] {
        let path = PathBuf::from(format!("{}.{}.bedgraph", prefix.display(), name));
        let mut w = std::io::BufWriter::new(std::fs::File::create(&path)?);
        writeln!(w, "track type=bedGraph name=\"molecule coverage ({})\"", if si == 0 { "+" } else { "-" })?;
        let cv = &cov[si];
        let mut i = 0usize;
        while i < cv.len() {
            let v = cv[i];
            let mut j = i;
            while j + 1 < cv.len() && cv[j + 1] == v {
                j += 1;
            }
            if v > 0 {
                writeln!(w, "{chrom}\t{}\t{}\t{v}", start as usize + i, start as usize + j + 1)?;
            }
            i = j + 1;
        }
        eprintln!("wrote {}", path.display());
    }
    let path = PathBuf::from(format!("{}.junctions.bed", prefix.display()));
    let mut w = std::io::BufWriter::new(std::fs::File::create(&path)?);
    writeln!(w, "track name=\"molecule junctions\" graphType=junctions")?;
    let mut k = 0usize;
    for (si, js) in juncs.iter().enumerate() {
        let strand = if si == 0 { '+' } else { '-' };
        let mut js2: Vec<&(u32, u32, u64)> = js.iter().collect();
        js2.sort_unstable();
        for (dn, ac, n) in js2 {
            let (bs, be) = (dn.saturating_sub(20), ac + 20);
            writeln!(
                w,
                "{chrom}\t{bs}\t{be}\tJUNC{k:05}\t{}\t{strand}\t{bs}\t{be}\t0,0,0\t2\t20,20\t0,{}",
                n.min(&1000), ac - bs
            )?;
            k += 1;
        }
    }
    eprintln!("wrote {}", path.display());
    Ok(())
}

#[cfg(test)]
mod junction_listing_tests {
    use super::*;
    use evidence_io::archive::put_varint;

    fn catalogue_fixture() -> Vec<u8> {
        let mut raw = Vec::new();
        // chr0:100-200, chr0:150-250, chr1:10-20.
        for value in [0, 100, 100, 0, 50, 100, 1, 10, 10] {
            put_varint(&mut raw, value);
        }
        raw
    }

    fn postings_fixture() -> Vec<u8> {
        let mut raw = Vec::new();
        // support, n_chunks, delta-coded chunks.
        for value in [4, 2, 1, 2, 8, 1, 2, 3, 0] {
            put_varint(&mut raw, value);
        }
        raw
    }

    #[test]
    fn decodes_coordinate_and_posting_delta_streams() {
        assert_eq!(
            parse_junction_catalogue(&catalogue_fixture()).unwrap(),
            vec![(0, 100, 200), (0, 150, 250), (1, 10, 20)]
        );
        assert_eq!(
            parse_junction_postings(&postings_fixture(), 3).unwrap(),
            vec![(4, vec![1, 3]), (8, vec![2]), (3, vec![])]
        );
        assert!(parse_junction_postings(&postings_fixture(), 2).is_err());
    }

    #[test]
    fn uniform_top_selection_uses_visible_barcode_ties_without_ordering_the_table() {
        let cell_dictionary = vec![
            umi::pack(b"TTTTTTTTTTTTTTTT").unwrap(),
            umi::pack(b"AAAAAAAAAAAAAAAA").unwrap(),
            umi::pack(b"CCCCCCCCCCCCCCCC").unwrap(),
        ];
        let mut counts = ScopedCounts {
            // Historical selection used archive cell id as its secondary key.
            cells: vec![(0, 7), (1, 7), (2, 7)],
            groups: Vec::new(),
            total_umis: 21,
        };
        assert_eq!(
            counts
                .cells
                .iter()
                .take(2)
                .map(|row| row.0)
                .collect::<Vec<_>>(),
            vec![0, 1],
            "legacy cell-id tie behavior changed"
        );
        sort_scoped_cell_counts(&mut counts.cells, Some(&cell_dictionary));
        assert_eq!(
            counts
                .cells
                .iter()
                .take(2)
                .map(|row| row.0)
                .collect::<Vec<_>>(),
            vec![1, 2],
            "uniform top-N did not select barcodes in ascending order across the boundary tie"
        );

        let scope = QueryScope {
            selected: None,
            group_names: Vec::new(),
            group_of: FxHashMap::default(),
            selected_per_group: Vec::new(),
            selected_cells: 3,
            aggregation: QueryAggregation::Cell,
            source: "all",
            source_path: None,
            source_content_blake3: None,
            resolved_mapping_blake3: None,
            archive_cells: 3,
            active: false,
        };
        let output = stream_uniform_counts(
            Vec::new(),
            JUNCTION_UNIFORM_RESULT_SCHEMA,
            JUNCTION_UNIFORM_COUNTS_SCHEMA,
            UniformQueryFormat::Json,
            &ResultContext::default(),
            &json!({"fixture": "boundary-tie"}),
            &counts,
            &scope,
            &cell_dictionary,
            2,
        )
        .unwrap();
        let value: serde_json::Value = serde_json::from_slice(&output).unwrap();
        let table = &value["data"]["tables"][0];
        assert_eq!(
            table["rows"],
            json!([
                ["cell", "AAAAAAAAAAAAAAAA", 7, null, null],
                ["cell", "CCCCCCCCCCCCCCCC", 7, null, null]
            ])
        );
        assert_eq!(
            table["selection"],
            json!({"available_rows": 3, "emitted_rows": 2, "truncated": true})
        );
        assert!(table["schema"]["semantics"]["ordered_by"].is_null());
    }

    #[test]
    fn uniform_scope_digest_is_lazy_and_exact_across_sparse_and_dense_paths() {
        let scope = |archive_cells| QueryScope {
            selected: Some([0, 1].into_iter().collect()),
            group_names: vec!["A".to_owned(), "B".to_owned()],
            group_of: [(0, 0), (1, 1)].into_iter().collect(),
            selected_per_group: vec![1, 1],
            selected_cells: 2,
            aggregation: QueryAggregation::Group,
            source: "groups",
            source_path: None,
            source_content_blake3: None,
            resolved_mapping_blake3: None,
            archive_cells,
            active: true,
        };
        let mut sparse = scope(10_000);
        let mut dense = scope(2);
        assert!(sparse.resolved_mapping_blake3.is_none());
        assert!(dense.resolved_mapping_blake3.is_none());
        sparse.ensure_resolved_mapping_digest().unwrap();
        dense.ensure_resolved_mapping_digest().unwrap();
        assert_eq!(
            sparse.resolved_mapping_blake3,
            dense.resolved_mapping_blake3
        );
        assert!(sparse
            .resolved_mapping_blake3
            .as_deref()
            .unwrap()
            .starts_with("blake3:"));

        let mut all_cells = QueryScope {
            selected: None,
            group_names: Vec::new(),
            group_of: FxHashMap::default(),
            selected_per_group: Vec::new(),
            selected_cells: 10_000,
            aggregation: QueryAggregation::Bulk,
            source: "all",
            source_path: None,
            source_content_blake3: None,
            resolved_mapping_blake3: None,
            archive_cells: 10_000,
            active: true,
        };
        all_cells.ensure_resolved_mapping_digest().unwrap();
        assert!(all_cells.resolved_mapping_blake3.is_none());
    }

    #[test]
    fn uniform_text_input_keeps_parsing_and_digest_on_one_snapshot() {
        let path = std::env::temp_dir().join(format!(
            "gravlax-bound-query-plan-{}-{}.tsv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let original = b"id\tkind\tlocus\noriginal\tregion\tchr1:10-20\n";
        std::fs::write(&path, original).unwrap();
        let bound = read_bound_query_text(&path, "test batch plan").unwrap();

        // A path mutation after capture cannot make reported provenance describe the new bytes
        // while execution still uses the old parsed snapshot.
        std::fs::write(&path, b"id\tkind\tlocus\nreplacement\tregion\tchr1:30-40\n").unwrap();
        let specs = parse_batch_plan(&bound.text, &["chr1".to_owned()]).unwrap();
        assert_eq!(specs[0].id, "original");
        assert_eq!(
            bound.content_blake3,
            format!("blake3:{}", blake3::hash(original).to_hex())
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn uniform_annotation_keeps_model_and_digest_on_one_snapshot() {
        let path = std::env::temp_dir().join(format!(
            "gravlax-bound-query-annotation-{}-{}.gtf",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let original = b"chr1\tt\texon\t1\t10\t.\t+\t.\tgene_id \"G_OLD\"; transcript_id \"T_OLD\";\n";
        std::fs::write(&path, original).unwrap();
        let (annotation, digest) =
            load_query_annotation(&path, "test annotation", true).unwrap();

        std::fs::write(
            &path,
            b"chr1\tt\texon\t20\t30\t.\t+\t.\tgene_id \"G_NEW\"; transcript_id \"T_NEW\";\n",
        )
        .unwrap();
        assert_eq!(annotation.gene_ids, vec!["G_OLD".to_owned()]);
        let expected = format!("blake3:{}", blake3::hash(original).to_hex());
        assert_eq!(digest.as_deref(), Some(expected.as_str()));
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn interval_semantics_are_half_open_and_explicit() {
        let row = JunctionMeta {
            chrom: 0,
            donor: 100,
            acceptor: 200,
            supporting_children: 1,
            posts: Vec::new(),
        };
        assert!(junction_in_window(&row, 100, 201, false));
        assert!(!junction_in_window(&row, 100, 200, false));
        assert!(junction_in_window(&row, 100, 200, true));
        assert!(!junction_in_window(&row, 101, 200, true));
    }

    #[test]
    fn annotation_flags_distinguish_exact_and_endpoint_membership() {
        let mut annotation = JunctionAnnotation::default();
        annotation.exact.insert((100, 200));
        annotation.donors.insert(100);
        annotation.acceptors.insert(200);
        annotation.acceptors.insert(300);
        assert!(annotation.exact.contains(&(100, 200)));
        assert!(!annotation.exact.contains(&(100, 300)));
        assert!(annotation.donors.contains(&100));
        assert!(annotation.acceptors.contains(&300));
    }

    #[test]
    fn junction_set_categories_are_exclusive_and_both_is_not_informative() {
        let mut counts = JunctionSetCounts::default();
        counts.add_mask(1);
        counts.add_mask(2);
        counts.add_mask(3);
        counts.add_mask(3);
        assert_eq!(
            counts,
            JunctionSetCounts {
                include_only: 1,
                exclude_only: 1,
                both: 2
            }
        );
        assert_eq!(counts.informative(), 2);
        assert_eq!(counts.total(), 4);
        assert_eq!(counts.usage(), Some(0.5));

        let only_both = JunctionSetCounts {
            include_only: 0,
            exclude_only: 0,
            both: 3,
        };
        assert_eq!(only_both.informative(), 0);
        assert_eq!(only_both.usage(), None);
        assert!(only_both.json()["usage_fraction"].is_null());
    }

    #[test]
    fn junction_set_rejects_duplicates_and_represents_absence() {
        let chroms = vec!["chr1".to_owned()];
        let metadata = vec![JunctionMeta {
            chrom: 0,
            donor: 100,
            acceptor: 200,
            supporting_children: 4,
            posts: vec![1],
        }];
        let requests = prepare_junction_set(
            &["chr1:100-200".to_owned()],
            &["chr1:300-400".to_owned()],
            &chroms,
            &metadata,
        )
        .unwrap();
        assert!(requests[0].metadata.is_some());
        assert!(requests[1].metadata.is_none());

        assert!(prepare_junction_set(
            &["chr1:100-200".to_owned(), "chr1:100-200".to_owned()],
            &["chr1:300-400".to_owned()],
            &chroms,
            &metadata,
        )
        .is_err());
        assert!(prepare_junction_set(
            &["chr1:100-200".to_owned()],
            &["chr1:100-200".to_owned()],
            &chroms,
            &metadata,
        )
        .is_err());
    }

    #[test]
    fn event_catalogue_discovers_all_coordinate_defined_types() {
        let row = |donor, acceptor| JunctionMeta {
            chrom: 0,
            donor,
            acceptor,
            supporting_children: 5,
            posts: vec![0],
        };
        let metadata = vec![row(100, 200), row(100, 300), row(150, 300), row(250, 300)];
        let events = discover_event_keys("chr1", 0, 0, 1_000, &metadata, &[], 2, 100).unwrap();
        assert!(events.iter().any(|event| {
            event.kind == EventTypeArg::AltAcceptor
                && event.includes == vec![(100, 200)]
                && event.excludes == vec![(100, 300)]
        }));
        assert!(events.iter().any(|event| {
            event.kind == EventTypeArg::AltDonor
                && event.includes == vec![(100, 300)]
                && event.excludes == vec![(150, 300)]
        }));
        assert!(events.iter().any(|event| {
            event.kind == EventTypeArg::Cassette
                && event.includes == vec![(100, 200), (250, 300)]
                && event.excludes == vec![(100, 300)]
        }));
    }

    #[test]
    fn event_catalogue_limit_is_a_hard_failure() {
        let metadata = vec![
            JunctionMeta {
                chrom: 0,
                donor: 100,
                acceptor: 200,
                supporting_children: 5,
                posts: vec![],
            },
            JunctionMeta {
                chrom: 0,
                donor: 100,
                acceptor: 300,
                supporting_children: 5,
                posts: vec![],
            },
            JunctionMeta {
                chrom: 0,
                donor: 100,
                acceptor: 400,
                supporting_children: 5,
                posts: vec![],
            },
        ];
        assert!(discover_event_keys(
            "chr1",
            0,
            0,
            1_000,
            &metadata,
            &[EventTypeArg::AltAcceptor],
            1,
            1,
        )
        .is_err());
        assert!(discover_event_keys(
            "chr1",
            0,
            0,
            1_000,
            &metadata,
            &[EventTypeArg::AltAcceptor],
            1,
            MAX_PACKED_EVENTS + 1,
        )
        .unwrap_err()
        .to_string()
        .contains("--max-events"));
    }

    fn grouped_test_scope() -> QueryScope {
        QueryScope {
            selected: Some([10, 11].into_iter().collect()),
            group_names: vec!["A".to_owned(), "B".to_owned()],
            group_of: [(10, 0), (11, 1)].into_iter().collect(),
            selected_per_group: vec![1, 1],
            selected_cells: 2,
            aggregation: QueryAggregation::Group,
            source: "groups",
            source_path: None,
            source_content_blake3: None,
            resolved_mapping_blake3: None,
            archive_cells: 2,
            active: true,
        }
    }

    fn test_event_definition(present: bool) -> EventDefinition {
        EventDefinition {
            key: EventKey {
                kind: EventTypeArg::AltDonor,
                chrom: "chr1".to_owned(),
                includes: vec![(100, 300)],
                excludes: vec![(200, 300)],
            },
            chrom_id: 0,
            components: Vec::new(),
            catalogue_present: present,
        }
    }

    #[test]
    fn packed_group_reduction_matches_cell_reference_without_retaining_cells() {
        let scope = grouped_test_scope();
        let mut hits = vec![
            pack_event_hit(0, 1, 1),
            pack_event_hit(0, 1, 2),
            pack_event_hit(0, 2, 1),
            pack_event_hit(0, 3, 2),
            pack_event_hit(0, 4, 1),
        ];
        hits.sort_unstable();
        let hits = reduce_sorted_packed_hits(hits);
        let cell_of = |class| match class {
            1 | 2 => Ok(10),
            3 => Ok(11),
            4 => Ok(12),
            _ => bail!("unexpected class"),
        };
        let reference =
            reduce_event_results_with(packed_hits_by_event(hits.clone(), 1), &scope, cell_of)
                .unwrap();
        let packed = reduce_packed_event_results_with(&hits, 1, &scope, cell_of).unwrap();

        assert_eq!(packed[0].totals, reference[0].totals);
        assert_eq!(packed[0].groups, reference[0].groups);
        assert_eq!(packed[0].support_cells, reference[0].support_cells);
        assert!(packed[0].cells.is_empty());
        assert_eq!(packed[0].totals.include_only, 1);
        assert_eq!(packed[0].totals.exclude_only, 1);
        assert_eq!(packed[0].totals.both, 1);

        let event = test_event_definition(true);
        assert_eq!(
            event_json(&event, &packed[0], &[], None, &scope, &[], 20),
            event_json(&event, &reference[0], &[], None, &scope, &[], 20),
        );

        let bulk_scope = QueryScope {
            selected: None,
            group_names: Vec::new(),
            group_of: FxHashMap::default(),
            selected_per_group: Vec::new(),
            selected_cells: 3,
            aggregation: QueryAggregation::Bulk,
            source: "all",
            source_path: None,
            source_content_blake3: None,
            resolved_mapping_blake3: None,
            archive_cells: 3,
            active: false,
        };
        let bulk_reference =
            reduce_event_results_with(packed_hits_by_event(hits.clone(), 1), &bulk_scope, cell_of)
                .unwrap();
        let bulk_packed = reduce_packed_event_results_with(&hits, 1, &bulk_scope, cell_of).unwrap();
        assert_eq!(bulk_packed[0].totals, bulk_reference[0].totals);
        assert_eq!(
            bulk_packed[0].support_cells,
            bulk_reference[0].support_cells
        );
        assert!(bulk_packed[0].cells.is_empty());
        assert_eq!(
            event_json(
                &event,
                &bulk_packed[0],
                &[],
                None,
                &bulk_scope,
                &[],
                20,
            ),
            event_json(
                &event,
                &bulk_reference[0],
                &[],
                None,
                &bulk_scope,
                &[],
                20,
            ),
        );
    }

    #[test]
    fn row_gate_requires_every_conservative_group_denominator() {
        let scope = grouped_test_scope();
        let result = EventResult {
            totals: JunctionSetCounts {
                include_only: 2,
                exclude_only: 1,
                both: 10,
            },
            cells: Vec::new(),
            support_cells: 2,
            groups: vec![
                (
                    JunctionSetCounts {
                        include_only: 1,
                        exclude_only: 1,
                        both: 8,
                    },
                    1,
                ),
                (
                    JunctionSetCounts {
                        include_only: 1,
                        exclude_only: 0,
                        both: 2,
                    },
                    1,
                ),
            ],
        };
        assert!(event_passes_row_gate(&result, &scope, 1));
        assert!(!event_passes_row_gate(&result, &scope, 2));
        assert!(event_passes_row_gate(&result, &scope, 0));
    }

    #[test]
    fn cohort_row_gate_cli_defaults_off_and_rejects_negative_values() {
        let args = CohortArgs::try_parse_from([
            "cohort",
            "events",
            "chr1:1-2",
            "--sample",
            "D0=a.aie",
            "--sample",
            "D1=b.aie",
        ])
        .unwrap();
        let min_row_informative = match args.what {
            CohortWhat::Events {
                min_row_informative,
                ..
            } => min_row_informative,
            _ => panic!("expected events arguments"),
        };
        assert_eq!(min_row_informative, 0);
        assert!(CohortArgs::try_parse_from([
            "cohort",
            "events",
            "chr1:1-2",
            "--sample",
            "D0=a.aie",
            "--sample",
            "D1=b.aie",
            "--min-row-informative",
            "-1",
        ])
        .is_err());
    }

    #[test]
    fn apa_cli_and_runtime_reject_invalid_scientific_options() {
        for arguments in [
            vec!["query", "archive.aie", "apa", "chr1:1-2", "--site-gap", "0"],
            vec!["query", "archive.aie", "apa", "chr1:1-2", "--drop-ip"],
            vec!["query", "archive.aie", "apa", "chr1:1-2", "--permute", "10"],
            vec![
                "query",
                "archive.aie",
                "apa-test",
                "--gtf",
                "genes.gtf",
                "--groups",
                "groups.tsv",
                "--site-gap",
                "0",
            ],
        ] {
            assert!(Args::try_parse_from(arguments).is_err());
        }

        let mut args = Args::try_parse_from([
            "query",
            "archive.aie",
            "apa",
            "chr1:1-2",
            "--genome",
            "genome.fa",
            "--drop-ip",
            "--groups",
            "groups.tsv",
            "--permute",
            "10",
        ])
        .unwrap();
        assert!(validate_query_args(&args.what).is_ok());

        if let What::Apa { site_gap, .. } = &mut args.what {
            *site_gap = 0;
        }
        assert!(validate_query_args(&args.what)
            .unwrap_err()
            .to_string()
            .contains("--site-gap must be at least 1"));

        if let What::Apa {
            site_gap,
            genome,
            groups,
            ..
        } = &mut args.what
        {
            *site_gap = 24;
            *genome = None;
            *groups = None;
        }
        let error = validate_query_args(&args.what).unwrap_err().to_string();
        assert!(error.contains("--drop-ip requires --genome"));

        if let What::Apa { drop_ip, .. } = &mut args.what {
            *drop_ip = false;
        }
        assert!(validate_query_args(&args.what)
            .unwrap_err()
            .to_string()
            .contains("--permute requires --groups"));
    }

    #[test]
    fn cohort_splice_graph_cli_has_registered_defaults_and_strict_thresholds() {
        let args = CohortArgs::try_parse_from([
            "cohort",
            "splice-graph",
            "chr1:1-2",
            "--design",
            "design.tsv",
            "--counts-only",
        ])
        .unwrap();
        match args.what {
            CohortWhat::SpliceGraph {
                locus,
                design,
                contrast,
                counts_only,
                min_support,
                min_edge_samples,
                min_sample_umis,
                min_replicates,
                min_path_umis,
                min_path_samples,
                max_paths,
                json,
                uniform_output,
            } => {
                assert_eq!(locus, "chr1:1-2");
                assert_eq!(design, PathBuf::from("design.tsv"));
                assert_eq!(contrast, None);
                assert!(counts_only);
                assert_eq!(min_support, 1);
                assert_eq!(min_edge_samples, 2);
                assert_eq!(min_sample_umis, 10);
                assert_eq!(min_replicates, 2);
                assert_eq!(min_path_umis, 5);
                assert_eq!(min_path_samples, 2);
                assert_eq!(max_paths, 100_000);
                assert!(!json);
                assert!(uniform_output.format.is_none());
            }
            _ => panic!("expected splice-graph arguments"),
        }
        assert!(CohortArgs::try_parse_from([
            "cohort",
            "splice-graph",
            "chr1:1-2",
            "--design",
            "design.tsv",
            "--min-edge-samples",
            "-1",
        ])
        .is_err());
        assert!(CohortArgs::try_parse_from([
            "cohort",
            "splice-graph",
            "chr1:1-2",
            "--design",
            "design.tsv",
            "--counts-only",
            "--contrast",
            "A:B",
        ])
        .is_err());
    }

    #[test]
    fn cohort_splice_graph_contrasts_are_ordered_and_strict() {
        assert_eq!(
            parse_graph_contrast("control:treated").unwrap(),
            ("control".to_owned(), "treated".to_owned())
        );
        for invalid in ["control", "control:", ":treated", "a:a", "a:b:c", "a/b:c"] {
            assert!(parse_graph_contrast(invalid).is_err(), "accepted {invalid}");
        }
    }

    #[test]
    fn cohort_named_paths_are_strict_and_deterministic() {
        assert_eq!(
            parse_named_paths(&["D0=/tmp/d0.aie".to_owned()], "sample")
                .unwrap(),
            vec![("D0".to_owned(), PathBuf::from("/tmp/d0.aie"))]
        );
        assert!(parse_named_paths(
            &["D0=/tmp/a".to_owned(), "D0=/tmp/b".to_owned()],
            "sample",
        )
        .is_err());
        assert!(parse_named_paths(&["bad/id=/tmp/a".to_owned()], "sample").is_err());
        assert!(parse_named_paths(&["D0".to_owned()], "sample").is_err());
    }

    #[test]
    fn fallback_junction_geometry_fails_closed_on_coordinate_overflow() {
        let shape = evidence_io::archive::Shape {
            blocks: vec![(0, 2), (10, 2)],
        };
        assert!(has_junction(&shape, 100, 102, 110).unwrap());
        assert!(has_junction(&shape, u32::MAX - 1, 0, 0).is_err());
    }
}
