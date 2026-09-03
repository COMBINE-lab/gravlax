//! Replicate-aware transcript-end analysis for 3'-tagged cohorts.
//!
//! The hot path reduces one archive to exact `(gene, endpoint, group)` counts after deduplicating
//! `(gene, UMI-class)`. Site construction happens once over the cohort. Cells and molecules are
//! measurements; the paired biological sample is the only inferential unit.

use crate::apastats;
use crate::archivecmd::{decode_chunk, read_chunk_index, LazyArchive};
use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use flate2::bufread::MultiGzDecoder;
use gravlax_output::{
    publish_file_no_clobber, reported_output_path, DataType, Durability, Field, OrderKey,
    OutputError, OutputFormat, Producer, Provenance, ResultContext, RowSemantics, SortDirection,
    StreamingBundleWriter, TableSchema, TableSemantics,
};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Serialize;
use serde_json::json;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

const BIN_BP: u32 = 65_536;
const KERNEL_MIN_BP: i32 = -50;
const KERNEL_MAX_BP: i32 = 2_000;
const KERNEL_BIN_BP: i32 = 10;
const PAS_MOTIFS: [&[u8; 6]; 12] = [
    b"AATAAA", b"ATTAAA", b"TATAAA", b"AGTAAA", b"AAGAAA", b"AATACA", b"CATAAA", b"GATAAA",
    b"AATATA", b"AATAGA", b"ACTAAA", b"AATGAA",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum UniformEndReportFormat {
    Text,
    Tsv,
    Json,
}

impl From<UniformEndReportFormat> for OutputFormat {
    fn from(value: UniformEndReportFormat) -> Self {
        match value {
            UniformEndReportFormat::Text => Self::Text,
            UniformEndReportFormat::Tsv => Self::Tsv,
            UniformEndReportFormat::Json => Self::Json,
        }
    }
}

#[derive(Parser, Clone)]
pub struct Args {
    /// Strict sample<TAB>condition<TAB>archive<TAB>groups design TSV.
    #[arg(long)]
    pub design: PathBuf,
    /// GTF or compiled AIC. It supplies gene windows and labels, never observed site positions.
    #[arg(long)]
    pub gtf: PathBuf,
    /// Stamped reference FASTA used for internal-priming and PAS-motif evidence.
    #[arg(long)]
    pub genome: PathBuf,
    /// Same-assembly PolyASite BED or BED.gz.
    #[arg(long)]
    pub polyasite: PathBuf,
    /// Ordered GROUP_A:GROUP_B paired contrast; effects are B minus A.
    #[arg(long)]
    pub group_contrast: String,
    /// Destination directory. It must not already exist.
    #[arg(long)]
    pub out_dir: PathBuf,
    /// Adjacent endpoint coordinates at most this far apart form one site.
    #[arg(long, default_value_t = 24)]
    pub site_gap: u32,
    /// Extend gene windows this far in transcript 3' direction.
    #[arg(long, default_value_t = 2_000)]
    pub tail_extend: u32,
    /// A donor supports site recurrence with at least this many primary-group UMIs.
    #[arg(long, default_value_t = 10)]
    pub min_site_umis: u64,
    /// Minimum recurrent donors for a site to enter the output catalogue.
    #[arg(long, default_value_t = 3)]
    pub min_site_samples: usize,
    /// Motif-only high-confidence sites require at least this many recurrent donors.
    #[arg(long, default_value_t = 4)]
    pub motif_min_samples: usize,
    /// Per donor, each contrasted group needs this many gene-site UMIs.
    #[arg(long, default_value_t = 20)]
    pub min_group_gene_umis: u64,
    /// Minimum eligible paired biological samples for a gene test.
    #[arg(long, default_value_t = 6)]
    pub min_samples: usize,
    /// Minimum most-distal-site UMIs across eligible sample/group rows.
    #[arg(long, default_value_t = 20)]
    pub min_distal_umis: u64,
    /// Hard site-catalogue bound. Exceeding it is an error; output is never truncated.
    #[arg(long, default_value_t = 1_000_000)]
    pub max_sites: usize,
    /// Deterministically shuffle the two group labels within each donor (negative control).
    #[arg(long)]
    pub shuffle_seed: Option<u64>,
    /// Emit a versioned uniform report while preserving the scientific output directory.
    #[arg(long, value_enum)]
    pub report_format: Option<UniformEndReportFormat>,
    /// Atomically publish the uniform report here without replacing an existing file.
    #[arg(long, requires = "report_format")]
    pub report_output: Option<PathBuf>,
}

#[derive(Parser, Clone)]
pub struct MixtureArgs {
    /// Strict sample<TAB>condition<TAB>archive<TAB>groups design TSV.
    #[arg(long)]
    pub design: PathBuf,
    /// GTF or compiled AIC supplying terminal-exon regions and labels.
    #[arg(long)]
    pub gtf: PathBuf,
    /// Stamped reference FASTA used to reject internal-priming candidates.
    #[arg(long)]
    pub genome: PathBuf,
    /// Same-assembly PolyASite BED or BED.gz supplying candidate site identities.
    #[arg(long)]
    pub polyasite: PathBuf,
    /// Ordered GROUP_A:GROUP_B paired contrast; effects are B minus A.
    #[arg(long)]
    pub group_contrast: String,
    /// Destination directory. It must not already exist.
    #[arg(long)]
    pub out_dir: PathBuf,
    /// Merge adjacent external catalogue coordinates within this distance.
    #[arg(long, default_value_t = 24)]
    pub site_gap: u32,
    /// Extend terminal-exon regions this far in transcript 3' direction.
    #[arg(long, default_value_t = 2_000)]
    pub tail_extend: u32,
    /// A donor supports a deconvolved site with at least this many expected UMIs.
    #[arg(long, default_value_t = 10)]
    pub min_site_umis: u64,
    /// Minimum recurrent donors for a site to enter the fitted catalogue.
    #[arg(long, default_value_t = 3)]
    pub min_site_samples: usize,
    /// Per donor, each contrasted group needs this many fitted gene-site UMIs.
    #[arg(long, default_value_t = 20)]
    pub min_group_gene_umis: u64,
    /// Minimum eligible paired biological samples for a gene test.
    #[arg(long, default_value_t = 6)]
    pub min_samples: usize,
    /// Minimum fitted most-distal-site UMIs across eligible sample/group rows.
    #[arg(long, default_value_t = 20)]
    pub min_distal_umis: u64,
    /// Hard candidate/recurrent-site bound. Exceeding it is an error; output is never truncated.
    #[arg(long, default_value_t = 1_000_000)]
    pub max_sites: usize,
    /// Deterministically shuffle the two group labels within each donor (negative control).
    #[arg(long)]
    pub shuffle_seed: Option<u64>,
    /// Emit a versioned uniform report while preserving the scientific output directory.
    #[arg(long, value_enum)]
    pub report_format: Option<UniformEndReportFormat>,
    /// Atomically publish the uniform report here without replacing an existing file.
    #[arg(long, requires = "report_format")]
    pub report_output: Option<PathBuf>,
}

impl From<MixtureArgs> for Args {
    fn from(value: MixtureArgs) -> Self {
        Self {
            design: value.design,
            gtf: value.gtf,
            genome: value.genome,
            polyasite: value.polyasite,
            group_contrast: value.group_contrast,
            out_dir: value.out_dir,
            site_gap: value.site_gap,
            tail_extend: value.tail_extend,
            min_site_umis: value.min_site_umis,
            min_site_samples: value.min_site_samples,
            motif_min_samples: 4,
            min_group_gene_umis: value.min_group_gene_umis,
            min_samples: value.min_samples,
            min_distal_umis: value.min_distal_umis,
            max_sites: value.max_sites,
            shuffle_seed: value.shuffle_seed,
            report_format: value.report_format,
            report_output: value.report_output,
        }
    }
}

#[derive(Clone)]
struct DesignRow {
    sample: String,
    condition: String,
    archive: PathBuf,
    archive_label: String,
    groups: PathBuf,
    groups_label: String,
}

#[derive(Clone, Copy)]
struct GeneWindow {
    gene: u32,
    rev: bool,
    lo: u32,
    hi: u32,
}

struct ChromWindows {
    windows: Vec<GeneWindow>,
    bins: Vec<Vec<u32>>,
}

struct GeneIndex {
    by_name: BTreeMap<String, Arc<ChromWindows>>,
    rev: Vec<bool>,
    chrom: Vec<String>,
}

#[derive(Clone, Copy)]
struct EndpointRecord {
    gene: u32,
    class: u32,
    endpoint: u32,
    group: u8,
}

struct SampleCounts {
    row: DesignRow,
    selected_cells: [usize; 2],
    absent_design_cells: [usize; 2],
    endpoint_records: usize,
    deduplicated_gene_umis: usize,
    outside_terminal_region_records: usize,
    ambiguous_terminal_region_records: usize,
    /// Sorted `(gene << 32 | endpoint, [group0, group1])` rows. A flat array is both smaller
    /// than a hash table and lets downstream gene fits borrow contiguous ranges without building
    /// a second sample-by-gene hierarchy.
    counts: Vec<(u64, [u32; 2])>,
    reference_signature: Option<evidence_io::genome::GenomeSig>,
    chunks: usize,
    archive_format_version: Option<u32>,
    archive_identity: Option<String>,
    groups_identity: Option<String>,
}

fn gene_count_rows(sample: &SampleCounts, gene: usize) -> &[(u64, [u32; 2])] {
    let lo = (gene as u64) << 32;
    let hi = ((gene as u64) + 1) << 32;
    let start = sample.counts.partition_point(|row| row.0 < lo);
    let end = sample.counts.partition_point(|row| row.0 < hi);
    &sample.counts[start..end]
}

fn collapse_endpoint_records(
    records: &mut Vec<EndpointRecord>,
    gene_rev: &[bool],
) -> Result<(Vec<(u64, [u32; 2])>, usize)> {
    records.par_sort_unstable_by_key(|record| (record.gene, record.class, record.endpoint));
    let mut deduplicated_gene_umis = 0usize;
    let mut write = 0usize;
    let mut start = 0usize;
    while start < records.len() {
        let mut end = start + 1;
        while end < records.len()
            && records[end].gene == records[start].gene
            && records[end].class == records[start].class
        {
            end += 1;
        }
        let gene = records[start].gene;
        let chosen = records[start..end].iter().map(|record| record.endpoint);
        let endpoint = if gene_rev[gene as usize] {
            chosen.min().unwrap()
        } else {
            chosen.max().unwrap()
        };
        records[write] = EndpointRecord {
            gene,
            class: 0,
            endpoint,
            group: records[start].group,
        };
        write += 1;
        deduplicated_gene_umis += 1;
        start = end;
    }
    records.truncate(write);
    records.par_sort_unstable_by_key(|record| (record.gene, record.endpoint, record.group));
    let mut counts = Vec::new();
    let mut start = 0usize;
    while start < records.len() {
        let gene = records[start].gene;
        let endpoint = records[start].endpoint;
        let mut end = start + 1;
        while end < records.len()
            && records[end].gene == gene
            && records[end].endpoint == endpoint
        {
            end += 1;
        }
        let mut groups = [0u32; 2];
        for record in &records[start..end] {
            groups[record.group as usize] = groups[record.group as usize]
                .checked_add(1)
                .context("transcript-end coordinate count overflow")?;
        }
        counts.push((((gene as u64) << 32) | endpoint as u64, groups));
        start = end;
    }
    Ok((counts, deduplicated_gene_umis))
}

struct Site {
    gene: u32,
    lo: u32,
    hi: u32,
    representative: u32,
    recurrent_samples: usize,
    counts: Vec<[u64; 2]>,
    ip: bool,
    a20: u32,
    arun: u32,
    motif: Option<String>,
    polyasite_distance: Option<u32>,
    high_confidence: bool,
}

struct GeneTest {
    gene: u32,
    site_ids: Vec<usize>,
    eligible_samples: Vec<usize>,
    usages: Vec<Option<[f64; 2]>>,
    effect: f64,
    t: f64,
    p: f64,
    p_flip: f64,
    q: f64,
    concordant: usize,
    lodo_stable: bool,
    distal_umis: u64,
    reported: bool,
}

fn valid_id(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"_.-".contains(&byte))
}

fn resolve(base: &Path, value: &str) -> PathBuf {
    let path = PathBuf::from(value);
    if path.is_absolute() {
        path
    } else {
        base.join(path)
    }
}

fn digest_bytes(scheme: &str, bytes: &[u8]) -> String {
    format!("{scheme}:{}", blake3::hash(bytes).to_hex())
}

fn same_file_snapshot(before: &std::fs::Metadata, after: &std::fs::Metadata) -> bool {
    same_file_object(before, after)
        && before.len() == after.len()
        && before.modified().ok() == after.modified().ok()
}

#[cfg(unix)]
fn same_file_object(left: &std::fs::Metadata, right: &std::fs::Metadata) -> bool {
    use std::os::unix::fs::MetadataExt;
    left.dev() == right.dev() && left.ino() == right.ino()
}

#[cfg(not(unix))]
fn same_file_object(left: &std::fs::Metadata, right: &std::fs::Metadata) -> bool {
    left.len() == right.len() && left.modified().ok() == right.modified().ok()
}

fn parse_design(path: &Path, capture_identity: bool) -> Result<(Vec<DesignRow>, Option<String>)> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading transcript-end design {}", path.display()))?;
    let identity = capture_identity.then(|| digest_bytes("full-file-blake3-v1", text.as_bytes()));
    let mut lines = text.lines();
    if lines.next().map(|line| line.trim_end_matches('\r'))
        != Some("sample\tcondition\tarchive\tgroups")
    {
        bail!(
            "transcript-end design header must be exactly: sample<TAB>condition<TAB>archive<TAB>groups"
        );
    }
    let base = path.parent().unwrap_or_else(|| Path::new("."));
    let mut samples = BTreeSet::new();
    let mut archives = BTreeSet::new();
    let mut out = Vec::new();
    for (line_index, raw) in lines.enumerate() {
        let line_no = line_index + 2;
        let line = raw.trim_end_matches('\r');
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 4 || fields.iter().any(|field| field.is_empty()) {
            bail!(
                "transcript-end design line {line_no} must have four nonempty tab-separated fields"
            );
        }
        if !valid_id(fields[0]) || !valid_id(fields[1]) {
            bail!(
                "transcript-end design line {line_no} has an invalid sample or condition identifier"
            );
        }
        if !samples.insert(fields[0].to_owned()) {
            bail!("duplicate transcript-end sample {}", fields[0]);
        }
        let archive = std::fs::canonicalize(resolve(base, fields[2])).with_context(|| {
            format!(
                "resolving transcript-end archive on line {line_no}: {}",
                fields[2]
            )
        })?;
        if !archives.insert(archive.clone()) {
            bail!("transcript-end design reuses archive {}", archive.display());
        }
        let groups = resolve(base, fields[3]);
        if !groups.is_file() {
            bail!(
                "transcript-end groups file on line {line_no} is absent: {}",
                groups.display()
            );
        }
        out.push(DesignRow {
            sample: fields[0].to_owned(),
            condition: fields[1].to_owned(),
            archive,
            archive_label: fields[2].to_owned(),
            groups,
            groups_label: fields[3].to_owned(),
        });
    }
    if out.len() < 2 {
        bail!("transcript-end design requires at least two biological samples");
    }
    Ok((out, identity))
}

fn parse_contrast(value: &str) -> Result<[String; 2]> {
    let fields: Vec<&str> = value.split(':').collect();
    if fields.len() != 2 || !valid_id(fields[0]) || !valid_id(fields[1]) || fields[0] == fields[1] {
        bail!("--group-contrast must be two distinct group IDs formatted GROUP_A:GROUP_B");
    }
    Ok([fields[0].to_owned(), fields[1].to_owned()])
}

fn build_gene_index(annotation: &anno::Annotation, tail_extend: u32) -> Result<GeneIndex> {
    let mut chrom_names = vec![String::new(); annotation.chrom_ids.len()];
    for (name, &id) in &annotation.chrom_ids {
        chrom_names[id as usize] = name.clone();
    }
    let mut raw: Vec<Vec<GeneWindow>> = (0..chrom_names.len()).map(|_| Vec::new()).collect();
    let mut identities: FxHashMap<u32, (u32, bool)> = FxHashMap::default();
    for transcript in &annotation.transcripts {
        let identity = identities
            .entry(transcript.gene)
            .or_insert((transcript.chrom, transcript.strand_rev));
        if *identity != (transcript.chrom, transcript.strand_rev) {
            bail!(
                "gene {} spans inconsistent chromosome/strand records",
                transcript.gene
            );
        }
        let Some(terminal) = (if transcript.strand_rev {
            transcript.exons.first()
        } else {
            transcript.exons.last()
        }) else {
            continue;
        };
        let (lo, hi) = if transcript.strand_rev {
            (terminal.start.saturating_sub(tail_extend), terminal.end)
        } else {
            (
                terminal.start,
                terminal
                    .end
                    .checked_add(tail_extend)
                    .context("terminal-exon tail-extension overflow")?,
            )
        };
        raw[transcript.chrom as usize].push(GeneWindow {
            gene: transcript.gene,
            rev: transcript.strand_rev,
            lo,
            hi,
        });
    }
    let mut rev = vec![false; annotation.gene_ids.len()];
    let mut gene_chrom = vec![String::new(); annotation.gene_ids.len()];
    for (gene, (chrom, strand_rev)) in identities {
        rev[gene as usize] = strand_rev;
        gene_chrom[gene as usize] = chrom_names[chrom as usize].clone();
    }
    let mut by_name = BTreeMap::new();
    for (chrom, mut windows) in raw.into_iter().enumerate() {
        if windows.is_empty() {
            continue;
        }
        windows.sort_unstable_by_key(|window| (window.lo, window.hi, window.gene));
        windows.dedup_by_key(|window| (window.lo, window.hi, window.gene));
        let max_end = windows.iter().map(|window| window.hi).max().unwrap_or(0);
        let mut bins = vec![Vec::new(); (max_end / BIN_BP + 1) as usize];
        for (index, window) in windows.iter().enumerate() {
            if window.hi <= window.lo {
                continue;
            }
            let first = window.lo / BIN_BP;
            let last = (window.hi - 1) / BIN_BP;
            for bin in first..=last {
                bins[bin as usize].push(index as u32);
            }
        }
        by_name.insert(
            chrom_names[chrom].clone(),
            Arc::new(ChromWindows { windows, bins }),
        );
    }
    Ok(GeneIndex {
        by_name,
        rev,
        chrom: gene_chrom,
    })
}

enum GeneMatch {
    Unique(u32),
    None,
    Ambiguous,
}

fn unique_gene(index: &ChromWindows, endpoint: u32, rev: bool) -> GeneMatch {
    let Some(bin) = index.bins.get((endpoint / BIN_BP) as usize) else {
        return GeneMatch::None;
    };
    let mut found = None;
    for &window_index in bin {
        let window = index.windows[window_index as usize];
        if window.rev == rev && endpoint >= window.lo && endpoint < window.hi {
            if found.is_some_and(|gene| gene != window.gene) {
                return GeneMatch::Ambiguous;
            }
            found = Some(window.gene);
        }
    }
    found.map_or(GeneMatch::None, GeneMatch::Unique)
}

fn molecule_endpoint(
    molecule: &crate::rows::MolRec,
    shapes: &[evidence_io::archive::Shape],
) -> u32 {
    let ends = molecule
        .chains
        .iter()
        .flat_map(|chain| chain.reps.iter())
        .map(|(position, shape)| {
            let end = *position
                + shapes[*shape as usize]
                    .blocks
                    .last()
                    .map(|block| block.0 + block.1)
                    .unwrap_or(0);
            (*position, end)
        })
        .chain(molecule.mms.iter().map(|(position, shape, _, _)| {
            let end = *position
                + shapes[*shape as usize]
                    .blocks
                    .last()
                    .map(|block| block.0 + block.1)
                    .unwrap_or(0);
            (*position, end)
        }));
    if molecule.strand_rev {
        ends.map(|(start, _)| start).min().unwrap_or(0)
    } else {
        ends.map(|(_, end)| end).max().unwrap_or(0)
    }
}

fn splitmix(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

fn load_group_map(
    archive: &mut LazyArchive,
    path: &Path,
    contrast: &[String; 2],
    shuffle_seed: Option<u64>,
    capture_identity: bool,
) -> Result<(FxHashMap<u32, u8>, [usize; 2], [usize; 2], Option<String>)> {
    let packed_cells = archive.cells()?.to_vec();
    let packed_to_id: FxHashMap<u32, u32> = packed_cells
        .iter()
        .enumerate()
        .map(|(index, &packed)| (packed, index as u32))
        .collect();
    let mut seen = BTreeSet::new();
    let mut assignments = Vec::new();
    let mut absent = [0usize; 2];
    let text = std::fs::read_to_string(path)?;
    let identity = capture_identity.then(|| digest_bytes("full-file-blake3-v1", text.as_bytes()));
    for (line_index, raw) in text.lines().enumerate() {
        let fields: Vec<&str> = raw.trim_end_matches('\r').split('\t').collect();
        if fields.len() != 2 || fields.iter().any(|field| field.is_empty()) {
            bail!(
                "{} line {} must contain barcode<TAB>group",
                path.display(),
                line_index + 1
            );
        }
        let Some(group) = contrast.iter().position(|value| value == fields[1]) else {
            continue;
        };
        let packed = evidence_io::umi::pack(fields[0].as_bytes()).with_context(|| {
            format!(
                "invalid barcode on {} line {}",
                path.display(),
                line_index + 1
            )
        })?;
        if !seen.insert(packed) {
            bail!(
                "duplicate primary-group barcode on {} line {}",
                path.display(),
                line_index + 1
            );
        }
        if let Some(&cell) = packed_to_id.get(&packed) {
            assignments.push((cell, group as u8));
        } else {
            absent[group] += 1;
        }
    }
    if assignments.iter().all(|assignment| assignment.1 != 0)
        || assignments.iter().all(|assignment| assignment.1 != 1)
    {
        bail!(
            "{} does not select archive cells in both contrasted groups",
            path.display()
        );
    }
    if let Some(seed) = shuffle_seed {
        let mut labels: Vec<u8> = assignments.iter().map(|assignment| assignment.1).collect();
        for index in (1..labels.len()).rev() {
            let other = (splitmix(seed ^ index as u64) % (index as u64 + 1)) as usize;
            labels.swap(index, other);
        }
        for ((_, group), shuffled) in assignments.iter_mut().zip(labels) {
            *group = shuffled;
        }
    }
    let mut selected = [0usize; 2];
    let mut map = FxHashMap::default();
    for (cell, group) in assignments {
        selected[group as usize] += 1;
        map.insert(cell, group);
    }
    Ok((map, selected, absent, identity))
}

fn reduce_sample(
    row: &DesignRow,
    index: &GeneIndex,
    contrast: &[String; 2],
    shuffle_seed: Option<u64>,
    capture_identity: bool,
) -> Result<SampleCounts> {
    let mut archive = LazyArchive::open(&row.archive)
        .with_context(|| format!("opening transcript-end sample {}", row.sample))?;
    let archive_before = capture_identity
        .then(|| archive.reader().file_metadata())
        .transpose()?;
    let (archive_format_version, archive_identity) = if capture_identity {
        let reader = archive.reader();
        let version = reader.archive_version();
        let identity = if let Some(commitment) = reader.content_commitment() {
            format!("aie-directory-root-v2:{}", commitment.to_hex())
        } else {
            let scan = reader.scan_legacy_identities()?;
            format!(
                "full-file-blake3-v1:{}",
                scan.full_file_blake3
                    .iter()
                    .map(|byte| format!("{byte:02x}"))
                    .collect::<String>()
            )
        };
        (Some(version), Some(identity))
    } else {
        (None, None)
    };
    let chunks = read_chunk_index(archive.reader())?;
    if chunks.iter().any(|chunk| chunk.chrom as usize >= archive.chrom_names.len()) {
        bail!("transcript-end sample {} has an out-of-range chunk chromosome", row.sample);
    }
    if chunks.windows(2).any(|pair| pair[1].chrom < pair[0].chrom) {
        bail!(
            "transcript-end sample {} has a non-monotone chromosome chunk index",
            row.sample
        );
    }
    let shapes = archive.shapes()?;
    let (group_of, selected_cells, absent_design_cells, groups_identity) = load_group_map(
        &mut archive,
        &row.groups,
        contrast,
        shuffle_seed,
        capture_identity,
    )?;
    let chrom_windows: Vec<Option<Arc<ChromWindows>>> = archive
        .chrom_names
        .iter()
        .map(|name| index.by_name.get(name).cloned())
        .collect();
    let mut counts = Vec::new();
    let mut endpoint_records = 0usize;
    let mut deduplicated_gene_umis = 0usize;
    let mut outside_terminal_region_records = 0usize;
    let mut ambiguous_terminal_region_records = 0usize;
    const CHUNK_BATCH: usize = 8;
    let mut chromosome_first = 0usize;
    while chromosome_first < chunks.len() {
        let chromosome = chunks[chromosome_first].chrom;
        let mut chromosome_last = chromosome_first + 1;
        while chromosome_last < chunks.len() && chunks[chromosome_last].chrom == chromosome {
            chromosome_last += 1;
        }
        let mut records = Vec::new();
        for first in (chromosome_first..chromosome_last).step_by(CHUNK_BATCH) {
            let last = (first + CHUNK_BATCH).min(chromosome_last);
            let decoded = {
                let (reader, tables) = archive.reader_and_tables();
                (first..last)
                    .into_par_iter()
                    .map(|chunk_index| -> Result<Vec<crate::rows::MolRec>> {
                        let (compressed, raw_len) =
                            reader.read_compressed_at(&format!("c{chunk_index}"))?;
                        let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                        decode_chunk(&raw, &chunks[chunk_index], None, tables)
                    })
                    .collect::<Result<Vec<_>>>()?
            };
            archive.prefetch_coc(
                decoded
                    .iter()
                    .flat_map(|molecules| molecules.iter().map(|m| m.umi_class)),
            )?;
            let classified = decoded
                .into_par_iter()
                .enumerate()
                .map(
                    |(offset, molecules)| -> Result<(Vec<EndpointRecord>, usize, usize)> {
                        let chunk = &chunks[first + offset];
                        let Some(windows) = chrom_windows[chunk.chrom as usize].as_deref() else {
                            return Ok((Vec::new(), 0, 0));
                        };
                        let mut local = Vec::new();
                        let mut outside = 0usize;
                        let mut ambiguous = 0usize;
                        for molecule in molecules {
                            let cell = archive.cell_of_cached(molecule.umi_class)?;
                            let Some(&group) = group_of.get(&cell) else {
                                continue;
                            };
                            let endpoint = molecule_endpoint(&molecule, &shapes);
                            match unique_gene(windows, endpoint, molecule.strand_rev) {
                                GeneMatch::Unique(gene) => local.push(EndpointRecord {
                                    gene,
                                    class: molecule.umi_class,
                                    endpoint,
                                    group,
                                }),
                                GeneMatch::None => outside += 1,
                                GeneMatch::Ambiguous => ambiguous += 1,
                            }
                        }
                        Ok((local, outside, ambiguous))
                    },
                )
                .collect::<Result<Vec<_>>>()?;
            for (mut local, outside, ambiguous) in classified {
                records.append(&mut local);
                outside_terminal_region_records += outside;
                ambiguous_terminal_region_records += ambiguous;
            }
        }
        endpoint_records += records.len();
        let (mut chromosome_counts, chromosome_umis) =
            collapse_endpoint_records(&mut records, &index.rev)?;
        deduplicated_gene_umis += chromosome_umis;
        counts.append(&mut chromosome_counts);
        archive.clear_coc_cache();
        chromosome_first = chromosome_last;
    }
    counts.par_sort_unstable_by_key(|row| row.0);
    let mut write = 0usize;
    for read in 0..counts.len() {
        let row = counts[read];
        if write > 0 && counts[write - 1].0 == row.0 {
            for group in 0..2 {
                counts[write - 1].1[group] = counts[write - 1].1[group]
                    .checked_add(row.1[group])
                    .context("transcript-end coordinate count overflow")?;
            }
        } else {
            counts[write] = row;
            write += 1;
        }
    }
    counts.truncate(write);
    if let Some(before) = archive_before.as_ref() {
        let after = archive.reader().file_metadata()?;
        let path_after = std::fs::metadata(&row.archive)?;
        if !same_file_snapshot(before, &after) || !same_file_object(&after, &path_after) {
            bail!(
                "transcript-end archive {} changed while it was analyzed",
                row.archive.display()
            );
        }
    }
    Ok(SampleCounts {
        row: row.clone(),
        selected_cells,
        absent_design_cells,
        endpoint_records,
        deduplicated_gene_umis,
        outside_terminal_region_records,
        ambiguous_terminal_region_records,
        counts,
        reference_signature: archive.genome_sig.clone(),
        chunks: chunks.len(),
        archive_format_version,
        archive_identity,
        groups_identity,
    })
}

fn reduce_design_bounded(
    design: &[DesignRow],
    index: &GeneIndex,
    contrast: &[String; 2],
    shuffle_seed: Option<u64>,
    capture_identity: bool,
) -> Result<Vec<SampleCounts>> {
    const SAMPLE_BATCH: usize = 2;
    const MAX_BATCH_ARCHIVE_BYTES: u64 = 256 * 1024 * 1024;
    let mut samples = Vec::with_capacity(design.len());
    let mut first = 0usize;
    while first < design.len() {
        let mut last = first + 1;
        let mut archive_bytes = std::fs::metadata(&design[first].archive)
            .map(|metadata| metadata.len())
            .unwrap_or(MAX_BATCH_ARCHIVE_BYTES);
        while last < design.len() && last - first < SAMPLE_BATCH {
            let next_bytes = std::fs::metadata(&design[last].archive)
                .map(|metadata| metadata.len())
                .unwrap_or(MAX_BATCH_ARCHIVE_BYTES);
            if archive_bytes.saturating_add(next_bytes) > MAX_BATCH_ARCHIVE_BYTES {
                break;
            }
            archive_bytes += next_bytes;
            last += 1;
        }
        let mut reduced = design[first..last]
            .par_iter()
            .map(|row| reduce_sample(row, index, contrast, shuffle_seed, capture_identity))
            .collect::<Result<Vec<_>>>()?;
        samples.append(&mut reduced);
        first = last;
    }
    Ok(samples)
}

fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(
            BufReader::new(file),
        ))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn normalize_chrom(value: &str) -> String {
    if value.starts_with("chr") {
        value.to_owned()
    } else {
        format!("chr{value}")
    }
}

fn polyasite_catalogue_identity(sites: &BTreeMap<(String, bool), Vec<u32>>) -> String {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"gravlax-polyasite-catalogue-v1\0");
    for ((chrom, reverse), positions) in sites {
        hasher.update(&(chrom.len() as u64).to_le_bytes());
        hasher.update(chrom.as_bytes());
        hasher.update(&[*reverse as u8]);
        hasher.update(&(positions.len() as u64).to_le_bytes());
        for position in positions {
            hasher.update(&position.to_le_bytes());
        }
    }
    format!(
        "gravlax-polyasite-catalogue-v1:{}",
        hasher.finalize().to_hex()
    )
}

fn load_polyasite(
    path: &Path,
    capture_identity: bool,
) -> Result<(BTreeMap<(String, bool), Vec<u32>>, Option<String>)> {
    let mut sites: BTreeMap<(String, bool), Vec<u32>> = BTreeMap::new();
    for (line_index, line) in open_maybe_gzip(path)?.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 || !matches!(fields[5], "+" | "-") {
            continue;
        }
        let start: u32 = fields[1]
            .parse()
            .with_context(|| format!("PolyASite line {} start", line_index + 1))?;
        let end: u32 = fields[2]
            .parse()
            .with_context(|| format!("PolyASite line {} end", line_index + 1))?;
        sites
            .entry((normalize_chrom(fields[0]), fields[5] == "-"))
            .or_default()
            .push(start + (end - start) / 2);
    }
    for values in sites.values_mut() {
        values.sort_unstable();
        values.dedup();
    }
    let identity = capture_identity.then(|| polyasite_catalogue_identity(&sites));
    Ok((sites, identity))
}

fn external_catalogue_by_gene(
    index: &GeneIndex,
    polyasite: &BTreeMap<(String, bool), Vec<u32>>,
    gene_count: usize,
    merge_gap: u32,
) -> Vec<Vec<u32>> {
    let mut by_gene: Vec<Vec<u32>> = (0..gene_count).map(|_| Vec::new()).collect();
    for ((chrom, rev), positions) in polyasite {
        let Some(windows) = index.by_name.get(chrom) else {
            continue;
        };
        for &position in positions {
            if let GeneMatch::Unique(gene) = unique_gene(windows, position, *rev) {
                by_gene[gene as usize].push(position);
            }
        }
    }
    for (gene, positions) in by_gene.iter_mut().enumerate() {
        positions.sort_unstable();
        positions.dedup();
        let rev = index.rev[gene];
        let mut merged = Vec::new();
        let mut start = 0usize;
        while start < positions.len() {
            let mut end = start + 1;
            while end < positions.len() && positions[end] - positions[end - 1] <= merge_gap {
                end += 1;
            }
            let cluster = &positions[start..end];
            merged.push(if rev {
                cluster[cluster.len() / 2]
            } else {
                cluster[(cluster.len() - 1) / 2]
            });
            start = end;
        }
        *positions = merged;
    }
    by_gene
}

fn upstream_distance(endpoint: u32, site: u32, rev: bool) -> i64 {
    if rev {
        endpoint as i64 - site as i64
    } else {
        site as i64 - endpoint as i64
    }
}

#[derive(Clone)]
struct KernelAudit {
    bins: Vec<u64>,
    no_candidate_umis: u64,
    unique_candidate_umis: u64,
    multiple_candidate_umis: u64,
}

fn learn_fragment_kernel(
    samples: &[SampleCounts],
    candidates: &[Vec<u32>],
    gene_rev: &[bool],
) -> KernelAudit {
    let n_bins = ((KERNEL_MAX_BP - KERNEL_MIN_BP) / KERNEL_BIN_BP + 1) as usize;
    let mut audit = KernelAudit {
        bins: vec![0; n_bins],
        no_candidate_umis: 0,
        unique_candidate_umis: 0,
        multiple_candidate_umis: 0,
    };
    for sample in samples {
        for &(key, groups) in &sample.counts {
            let gene = (key >> 32) as usize;
            let endpoint = key as u32;
            let weight = groups[0] as u64 + groups[1] as u64;
            let (lo, hi) = compatible_site_range(endpoint, gene_rev[gene]);
            let sites = &candidates[gene];
            let start = sites.partition_point(|&site| site < lo);
            let end = sites.partition_point(|&site| site <= hi);
            match end - start {
                0 => audit.no_candidate_umis += weight,
                1 => {
                    audit.unique_candidate_umis += weight;
                    let distance = upstream_distance(endpoint, sites[start], gene_rev[gene]);
                    let bin = ((distance as i32 - KERNEL_MIN_BP) / KERNEL_BIN_BP) as usize;
                    audit.bins[bin] += weight;
                }
                _ => audit.multiple_candidate_umis += weight,
            }
        }
    }
    audit
}

fn remove_internal_priming_candidates(
    candidates: &mut [Vec<u32>],
    gene_index: &GeneIndex,
    genome: &Path,
    reference_signature: Option<&evidence_io::genome::GenomeSig>,
    capture_identity: bool,
) -> Result<(usize, Option<String>)> {
    let mut by_chrom: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    for (gene, positions) in candidates.iter().enumerate() {
        for candidate in 0..positions.len() {
            by_chrom
                .entry(gene_index.chrom[gene].clone())
                .or_default()
                .push((gene, candidate));
        }
    }
    let mut flagged: Vec<Vec<bool>> = candidates
        .iter()
        .map(|positions| vec![false; positions.len()])
        .collect();
    let mut seen = BTreeSet::new();
    let mut verification_error = None;
    let mut observed_contigs = Vec::new();
    evidence_io::genome::for_each_contig(genome, |name, sequence| {
        if capture_identity {
            observed_contigs.push(evidence_io::genome::ContigSig {
                name: name.to_owned(),
                len: sequence.len() as u64,
                blake3: blake3::hash(sequence).to_hex().to_string(),
            });
        }
        let Some(entries) = by_chrom.get(name) else {
            return true;
        };
        if let Some(signature) = reference_signature {
            if let Err(error) = evidence_io::genome::verify_contig(signature, name, sequence) {
                verification_error = Some(error);
                return false;
            }
        }
        seen.insert(name.to_owned());
        for &(gene, candidate) in entries {
            let (a20, arun) =
                apastats::ip_stats(sequence, candidates[gene][candidate], gene_index.rev[gene]);
            flagged[gene][candidate] = apastats::is_internal_priming(a20, arun);
        }
        true
    })?;
    if let Some(error) = verification_error {
        return Err(error);
    }
    if let Some(missing) = by_chrom.keys().find(|name| !seen.contains(*name)) {
        bail!("external-site chromosome {missing} is absent from --genome");
    }
    let mut dropped = 0usize;
    for (positions, flags) in candidates.iter_mut().zip(flagged) {
        let mut index = 0usize;
        positions.retain(|_| {
            let keep = !flags[index];
            dropped += (!keep) as usize;
            index += 1;
            keep
        });
    }
    let identity = if capture_identity {
        if observed_contigs.is_empty() {
            bail!("{} contains no FASTA records", genome.display());
        }
        Some(format!(
            "{}:{}",
            evidence_io::genome::GENOME_SIG_ALGO,
            evidence_io::genome::GenomeSig::combined_digest(&observed_contigs)
        ))
    } else {
        None
    };
    Ok((dropped, identity))
}

fn smoothed_kernel(audit: &KernelAudit) -> Vec<f64> {
    let mut probabilities = vec![0.0; audit.bins.len()];
    for (index, value) in probabilities.iter_mut().enumerate() {
        let lo = index.saturating_sub(2);
        let hi = (index + 3).min(audit.bins.len());
        *value = 1.0 + audit.bins[lo..hi].iter().sum::<u64>() as f64;
    }
    let total = probabilities.iter().sum::<f64>();
    for value in &mut probabilities {
        *value /= total;
    }
    probabilities
}

fn cross_fitted_kernels(
    samples: &[SampleCounts],
    candidates: &[Vec<u32>],
    gene_rev: &[bool],
) -> (KernelAudit, Vec<KernelAudit>, Vec<Vec<f64>>, Vec<f64>) {
    let per_sample: Vec<KernelAudit> = samples
        .par_iter()
        .map(|sample| learn_fragment_kernel(std::slice::from_ref(sample), candidates, gene_rev))
        .collect();
    let mut pooled = KernelAudit {
        bins: vec![0; ((KERNEL_MAX_BP - KERNEL_MIN_BP) / KERNEL_BIN_BP + 1) as usize],
        no_candidate_umis: 0,
        unique_candidate_umis: 0,
        multiple_candidate_umis: 0,
    };
    for sample in &per_sample {
        for (total, &count) in pooled.bins.iter_mut().zip(&sample.bins) {
            *total += count;
        }
        pooled.no_candidate_umis += sample.no_candidate_umis;
        pooled.unique_candidate_umis += sample.unique_candidate_umis;
        pooled.multiple_candidate_umis += sample.multiple_candidate_umis;
    }
    let mut kernels = Vec::with_capacity(samples.len());
    let mut heldout_bits_per_umi = Vec::with_capacity(samples.len());
    for heldout in &per_sample {
        let training = KernelAudit {
            bins: pooled
                .bins
                .iter()
                .zip(&heldout.bins)
                .map(|(&total, &sample)| total - sample)
                .collect(),
            no_candidate_umis: pooled.no_candidate_umis - heldout.no_candidate_umis,
            unique_candidate_umis: pooled.unique_candidate_umis - heldout.unique_candidate_umis,
            multiple_candidate_umis: pooled.multiple_candidate_umis
                - heldout.multiple_candidate_umis,
        };
        let probabilities = smoothed_kernel(&training);
        let observations = heldout.bins.iter().sum::<u64>();
        let uniform_bits = (probabilities.len() as f64).log2();
        let model_bits = if observations == 0 {
            f64::NAN
        } else {
            -heldout
                .bins
                .iter()
                .zip(&probabilities)
                .map(|(&count, &probability)| count as f64 * probability.log2())
                .sum::<f64>()
                / observations as f64
        };
        heldout_bits_per_umi.push(uniform_bits - model_bits);
        kernels.push(probabilities);
    }
    (pooled, per_sample, kernels, heldout_bits_per_umi)
}

fn kernel_likelihood(endpoint: u32, site: u32, rev: bool, kernel: &[f64]) -> Option<f64> {
    let distance = upstream_distance(endpoint, site, rev);
    if distance < KERNEL_MIN_BP as i64 || distance > KERNEL_MAX_BP as i64 {
        return None;
    }
    let bin = ((distance as i32 - KERNEL_MIN_BP) / KERNEL_BIN_BP) as usize;
    kernel.get(bin).copied()
}

fn compatible_site_range(endpoint: u32, rev: bool) -> (u32, u32) {
    let (lo, hi) = if rev {
        (
            endpoint as i64 - KERNEL_MAX_BP as i64,
            endpoint as i64 - KERNEL_MIN_BP as i64,
        )
    } else {
        (
            endpoint as i64 + KERNEL_MIN_BP as i64,
            endpoint as i64 + KERNEL_MAX_BP as i64,
        )
    };
    (
        lo.clamp(0, u32::MAX as i64) as u32,
        hi.clamp(0, u32::MAX as i64) as u32,
    )
}

fn em_site_counts(endpoints: &[(u32, u32)], sites: &[u32], rev: bool, kernel: &[f64]) -> Vec<f64> {
    if sites.is_empty() {
        return Vec::new();
    }
    // Candidate compatibility and the fixed observation kernel do not change between EM
    // iterations.  Materializing only non-zero likelihoods avoids two full site scans per
    // endpoint per iteration, which is the dominant cohort-analysis hot path.
    let likelihoods: Vec<Vec<(usize, f64)>> = endpoints
        .iter()
        .map(|&(endpoint, _)| {
            let (lo, hi) = compatible_site_range(endpoint, rev);
            let start = sites.partition_point(|&site| site < lo);
            let end = sites.partition_point(|&site| site <= hi);
            (start..end)
                .filter_map(|index| {
                    kernel_likelihood(endpoint, sites[index], rev, kernel)
                        .filter(|&value| value > 0.0)
                        .map(|value| (index, value))
                })
                .collect()
        })
        .collect();
    let mut theta = vec![1.0 / sites.len() as f64; sites.len()];
    let mut expected = vec![0.0; sites.len()];
    for _ in 0..50 {
        expected.fill(0.0);
        for ((_, count), compatible) in endpoints.iter().zip(&likelihoods) {
            let denominator = compatible
                .iter()
                .map(|&(site, likelihood)| likelihood * theta[site])
                .sum::<f64>();
            if denominator == 0.0 {
                continue;
            }
            for &(site, likelihood) in compatible {
                expected[site] += *count as f64 * likelihood * theta[site] / denominator;
            }
        }
        let total = expected.iter().sum::<f64>();
        if total == 0.0 {
            return vec![0.0; sites.len()];
        }
        let mut change: f64 = 0.0;
        for (prior, &count) in theta.iter_mut().zip(&expected) {
            let updated = count / total;
            change = change.max((updated - *prior).abs());
            *prior = updated;
        }
        if change < 1e-8 {
            break;
        }
    }
    expected
}

struct MixtureSite {
    gene: u32,
    point: u32,
    recurrent_samples: usize,
    counts: Vec<[f64; 2]>,
}

struct MixtureAudit {
    candidate_sites: usize,
    recurrent_sites: usize,
    assigned_umis: f64,
}

fn nearest(sorted: &[u32], point: u32) -> Option<u32> {
    let index = sorted.partition_point(|&value| value < point);
    [
        index.checked_sub(1),
        (index < sorted.len()).then_some(index),
    ]
    .into_iter()
    .flatten()
    .map(|i| sorted[i].abs_diff(point))
    .min()
}

fn transcript_upstream(seq: &[u8], point: u32, rev: bool) -> Vec<u8> {
    if !rev {
        let end = point.saturating_sub(5) as usize;
        let start = point.saturating_sub(60) as usize;
        seq.get(start.min(seq.len())..end.min(seq.len()))
            .unwrap_or_default()
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect()
    } else {
        let start = point.saturating_add(5) as usize;
        let end = point.saturating_add(60) as usize;
        seq.get(start.min(seq.len())..end.min(seq.len()))
            .unwrap_or_default()
            .iter()
            .rev()
            .map(|base| match base.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }
}

fn pas_motif(seq: &[u8], point: u32, rev: bool) -> Option<String> {
    let upstream = transcript_upstream(seq, point, rev);
    PAS_MOTIFS
        .iter()
        .find(|motif| upstream.windows(6).any(|window| window == motif.as_slice()))
        .map(|motif| String::from_utf8_lossy(motif.as_slice()).into_owned())
}

fn median(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let middle = values.len() / 2;
    Some(if values.len() % 2 == 1 {
        values[middle]
    } else {
        (values[middle - 1] + values[middle]) / 2.0
    })
}

fn site_catalogue(
    samples: &[SampleCounts],
    gene_count: usize,
    gene_rev: &[bool],
    site_gap: u32,
    min_site_umis: u64,
    min_site_samples: usize,
    max_sites: usize,
) -> Result<Vec<Site>> {
    let mut coordinates: Vec<Vec<u32>> = (0..gene_count).map(|_| Vec::new()).collect();
    for sample in samples {
        for &(key, _) in &sample.counts {
            coordinates[(key >> 32) as usize].push(key as u32);
        }
    }
    let mut sites = Vec::new();
    for (gene, points) in coordinates.iter_mut().enumerate() {
        if points.is_empty() {
            continue;
        }
        points.sort_unstable();
        points.dedup();
        let mut start = 0usize;
        while start < points.len() {
            let mut end = start + 1;
            while end < points.len() && points[end] - points[end - 1] <= site_gap {
                end += 1;
            }
            let cluster = &points[start..end];
            let mut counts = vec![[0u64; 2]; samples.len()];
            let mut coordinate_weights: BTreeMap<u32, u64> = BTreeMap::new();
            for (sample_index, sample) in samples.iter().enumerate() {
                for &point in cluster {
                    let key = ((gene as u64) << 32) | point as u64;
                    if let Ok(index) = sample.counts.binary_search_by_key(&key, |row| row.0) {
                        let group_counts = &sample.counts[index].1;
                        counts[sample_index][0] += group_counts[0] as u64;
                        counts[sample_index][1] += group_counts[1] as u64;
                        *coordinate_weights.entry(point).or_default() +=
                            (group_counts[0] + group_counts[1]) as u64;
                    }
                }
            }
            let recurrent_samples = counts
                .iter()
                .filter(|row| row[0] + row[1] >= min_site_umis)
                .count();
            if recurrent_samples >= min_site_samples {
                let representative = coordinate_weights
                    .iter()
                    .map(|(&point, &weight)| {
                        (
                            weight,
                            if gene_rev[gene] {
                                u32::MAX - point
                            } else {
                                point
                            },
                            point,
                        )
                    })
                    .max()
                    .map(|(_, _, point)| point)
                    .unwrap_or(cluster[0]);
                if sites.len() == max_sites {
                    bail!(
                        "transcript-end catalogue exceeds --max-sites {max_sites}; output was not truncated"
                    );
                }
                sites.push(Site {
                    gene: gene as u32,
                    lo: cluster[0],
                    hi: cluster[cluster.len() - 1],
                    representative,
                    recurrent_samples,
                    counts,
                    ip: false,
                    a20: 0,
                    arun: 0,
                    motif: None,
                    polyasite_distance: None,
                    high_confidence: false,
                });
            }
            start = end;
        }
    }
    Ok(sites)
}

fn score_sites(
    sites: &mut [Site],
    annotation: &anno::Annotation,
    gene_rev: &[bool],
    genome: &Path,
    polyasite: &BTreeMap<(String, bool), Vec<u32>>,
    reference_signature: Option<&evidence_io::genome::GenomeSig>,
    motif_min_samples: usize,
    capture_identity: bool,
) -> Result<Option<String>> {
    let mut chrom_names = vec![String::new(); annotation.chrom_ids.len()];
    for (name, &id) in &annotation.chrom_ids {
        chrom_names[id as usize] = name.clone();
    }
    let gene_chrom: Vec<u32> = {
        let mut value = vec![u32::MAX; annotation.gene_ids.len()];
        for transcript in &annotation.transcripts {
            value[transcript.gene as usize] = transcript.chrom;
        }
        value
    };
    let mut by_chrom: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (site_id, site) in sites.iter().enumerate() {
        by_chrom
            .entry(chrom_names[gene_chrom[site.gene as usize] as usize].clone())
            .or_default()
            .push(site_id);
    }
    let mut seen = BTreeSet::new();
    let mut verification_error = None;
    let mut observed_contigs = Vec::new();
    evidence_io::genome::for_each_contig(genome, |name, sequence| {
        if capture_identity {
            observed_contigs.push(evidence_io::genome::ContigSig {
                name: name.to_owned(),
                len: sequence.len() as u64,
                blake3: blake3::hash(sequence).to_hex().to_string(),
            });
        }
        let Some(ids) = by_chrom.get(name) else {
            return true;
        };
        if let Some(signature) = reference_signature {
            if let Err(error) = evidence_io::genome::verify_contig(signature, name, sequence) {
                verification_error = Some(error);
                return false;
            }
        }
        seen.insert(name.to_owned());
        for &site_id in ids {
            let site = &mut sites[site_id];
            let rev = gene_rev[site.gene as usize];
            let (a20, arun) = apastats::ip_stats(sequence, site.representative, rev);
            site.a20 = a20;
            site.arun = arun;
            site.ip = apastats::is_internal_priming(a20, arun);
            site.motif = pas_motif(sequence, site.representative, rev);
            site.polyasite_distance = polyasite
                .get(&(name.to_owned(), rev))
                .and_then(|positions| nearest(positions, site.representative));
            site.high_confidence = !site.ip
                && (site
                    .polyasite_distance
                    .is_some_and(|distance| distance <= 50)
                    || (site.motif.is_some() && site.recurrent_samples >= motif_min_samples));
        }
        true
    })?;
    if let Some(error) = verification_error {
        return Err(error);
    }
    if let Some(missing) = by_chrom.keys().find(|name| !seen.contains(*name)) {
        bail!("annotation/site chromosome {missing} is absent from --genome");
    }
    if capture_identity {
        if observed_contigs.is_empty() {
            bail!("{} contains no FASTA records", genome.display());
        }
        Ok(Some(format!(
            "{}:{}",
            evidence_io::genome::GENOME_SIG_ALGO,
            evidence_io::genome::GenomeSig::combined_digest(&observed_contigs)
        )))
    } else {
        Ok(None)
    }
}

fn fit_genes(
    sites: &[Site],
    gene_count: usize,
    gene_rev: &[bool],
    min_group_gene_umis: u64,
    min_samples: usize,
    min_distal_umis: u64,
    site_eligible: fn(&Site) -> bool,
) -> Vec<GeneTest> {
    let sample_count = sites.first().map_or(0, |site| site.counts.len());
    let mut by_gene: Vec<Vec<usize>> = (0..gene_count).map(|_| Vec::new()).collect();
    for (site_id, site) in sites.iter().enumerate() {
        if site_eligible(site) {
            by_gene[site.gene as usize].push(site_id);
        }
    }
    let mut tests = Vec::new();
    for (gene, site_ids) in by_gene.iter_mut().enumerate() {
        if site_ids.len() < 2 {
            continue;
        }
        site_ids.sort_unstable_by_key(|&site_id| {
            let point = sites[site_id].representative;
            if gene_rev[gene] {
                u32::MAX - point
            } else {
                point
            }
        });
        let denominator = (site_ids.len() - 1) as f64;
        let mut usages = vec![None; sample_count];
        let mut eligible_samples = Vec::new();
        let mut differences = Vec::new();
        for sample in 0..sample_count {
            let mut totals = [0u64; 2];
            let mut weighted = [0.0f64; 2];
            for (rank, &site_id) in site_ids.iter().enumerate() {
                for group in 0..2 {
                    let count = sites[site_id].counts[sample][group];
                    totals[group] += count;
                    weighted[group] += count as f64 * rank as f64 / denominator;
                }
            }
            if totals.iter().all(|&total| total >= min_group_gene_umis) {
                let value = [
                    weighted[0] / totals[0] as f64,
                    weighted[1] / totals[1] as f64,
                ];
                usages[sample] = Some(value);
                eligible_samples.push(sample);
                differences.push(value[1] - value[0]);
            }
        }
        if differences.len() < min_samples {
            continue;
        }
        let distal = *site_ids.last().unwrap();
        let distal_umis = eligible_samples
            .iter()
            .map(|&sample| sites[distal].counts[sample][0] + sites[distal].counts[sample][1])
            .sum();
        if distal_umis < min_distal_umis {
            continue;
        }
        let effect = differences.iter().sum::<f64>() / differences.len() as f64;
        let Some((t, p, p_flip)) = apastats::paired_test(&differences) else {
            continue;
        };
        let concordant = differences
            .iter()
            .filter(|&&difference| difference != 0.0 && difference.signum() == effect.signum())
            .count();
        let lodo_stable = (0..differences.len()).all(|left_out| {
            let mean = differences
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != left_out)
                .map(|(_, value)| value)
                .sum::<f64>()
                / (differences.len() - 1) as f64;
            mean != 0.0 && mean.signum() == effect.signum()
        });
        tests.push(GeneTest {
            gene: gene as u32,
            site_ids: site_ids.clone(),
            eligible_samples,
            usages,
            effect,
            t,
            p,
            p_flip,
            q: 1.0,
            concordant,
            lodo_stable,
            distal_umis,
            reported: false,
        });
    }
    let q = apastats::bh_fdr(&tests.iter().map(|test| test.p).collect::<Vec<_>>());
    for (test, q_value) in tests.iter_mut().zip(q) {
        test.q = q_value;
        test.reported =
            test.q <= 0.05 && test.effect.abs() >= 0.10 && test.concordant >= 6 && test.lodo_stable;
    }
    tests.sort_by(|a, b| a.p.partial_cmp(&b.p).unwrap_or(Ordering::Equal));
    tests
}

fn registered_site(site: &Site) -> bool {
    site.high_confidence
}

fn polyasite_site(site: &Site) -> bool {
    !site.ip
        && site
            .polyasite_distance
            .is_some_and(|distance| distance <= 50)
}

fn fit_polyasite_mixture(
    samples: &[SampleCounts],
    candidates: &[Vec<u32>],
    gene_rev: &[bool],
    sample_kernels: &[Vec<f64>],
    min_site_umis: u64,
    min_site_samples: usize,
    min_group_gene_umis: u64,
    min_samples: usize,
    min_distal_umis: u64,
) -> (Vec<MixtureSite>, Vec<GeneTest>, MixtureAudit) {
    debug_assert_eq!(samples.len(), sample_kernels.len());
    let gene_count = candidates.len();
    let mut mixture_sites = Vec::new();
    let mut site_ids_by_gene: Vec<Vec<usize>> = (0..gene_count).map(|_| Vec::new()).collect();
    let mut candidate_sites = 0usize;
    let mut assigned_umis = 0.0;
    for gene in 0..gene_count {
        if candidates[gene].is_empty() {
            continue;
        }
        let mut pooled_endpoints: Vec<u32> = samples
            .iter()
            .flat_map(|sample| gene_count_rows(sample, gene).iter().map(|row| row.0 as u32))
            .collect();
        pooled_endpoints.sort_unstable();
        pooled_endpoints.dedup();
        let active: Vec<u32> = candidates[gene]
            .iter()
            .copied()
            .filter(|&site| {
                let (lo, hi) = if gene_rev[gene] {
                    (
                        site as i64 + KERNEL_MIN_BP as i64,
                        site as i64 + KERNEL_MAX_BP as i64,
                    )
                } else {
                    (
                        site as i64 - KERNEL_MAX_BP as i64,
                        site as i64 - KERNEL_MIN_BP as i64,
                    )
                };
                let lo = lo.max(0) as u32;
                let hi = hi.clamp(0, u32::MAX as i64) as u32;
                let index = pooled_endpoints.partition_point(|&endpoint| endpoint < lo);
                pooled_endpoints
                    .get(index)
                    .is_some_and(|&endpoint| endpoint <= hi)
            })
            .collect();
        if active.is_empty() {
            continue;
        }
        candidate_sites += active.len();
        let estimate = |sites: &[u32]| -> Vec<Vec<[f64; 2]>> {
            let mut counts = vec![vec![[0.0; 2]; samples.len()]; sites.len()];
            for (sample, sample_counts) in samples.iter().enumerate() {
                let gene_rows = gene_count_rows(sample_counts, gene);
                for group in 0..2 {
                    let endpoints: Vec<(u32, u32)> = gene_rows
                        .iter()
                        .filter_map(|&(key, values)| {
                            (values[group] > 0).then_some((key as u32, values[group]))
                        })
                        .collect();
                    for (site, value) in
                        em_site_counts(&endpoints, sites, gene_rev[gene], &sample_kernels[sample])
                            .into_iter()
                            .enumerate()
                    {
                        counts[site][sample][group] = value;
                    }
                }
            }
            counts
        };
        let initial = estimate(&active);
        let recurrent: Vec<u32> = active
            .iter()
            .enumerate()
            .filter(|(site, _)| {
                initial[*site]
                    .iter()
                    .filter(|row| row[0] + row[1] >= min_site_umis as f64)
                    .count()
                    >= min_site_samples
            })
            .map(|(_, &site)| site)
            .collect();
        if recurrent.is_empty() {
            continue;
        }
        let final_counts = estimate(&recurrent);
        for (site, (&point, counts)) in recurrent.iter().zip(final_counts).enumerate() {
            assigned_umis += counts.iter().map(|row| row[0] + row[1]).sum::<f64>();
            let recurrent_samples = counts
                .iter()
                .filter(|row| row[0] + row[1] >= min_site_umis as f64)
                .count();
            let id = mixture_sites.len();
            mixture_sites.push(MixtureSite {
                gene: gene as u32,
                point,
                recurrent_samples,
                counts,
            });
            debug_assert_eq!(site_ids_by_gene[gene].len(), site);
            site_ids_by_gene[gene].push(id);
        }
    }
    let mut tests = Vec::new();
    for (gene, site_ids) in site_ids_by_gene.iter_mut().enumerate() {
        if site_ids.len() < 2 {
            continue;
        }
        site_ids.sort_unstable_by_key(|&site| {
            let point = mixture_sites[site].point;
            if gene_rev[gene] {
                u32::MAX - point
            } else {
                point
            }
        });
        let denominator = (site_ids.len() - 1) as f64;
        let mut usages = vec![None; samples.len()];
        let mut eligible_samples = Vec::new();
        let mut differences = Vec::new();
        for sample in 0..samples.len() {
            let mut totals = [0.0; 2];
            let mut weighted = [0.0; 2];
            for (rank, &site) in site_ids.iter().enumerate() {
                for group in 0..2 {
                    let count = mixture_sites[site].counts[sample][group];
                    totals[group] += count;
                    weighted[group] += count * rank as f64 / denominator;
                }
            }
            if totals
                .iter()
                .all(|&total| total >= min_group_gene_umis as f64)
            {
                let value = [weighted[0] / totals[0], weighted[1] / totals[1]];
                usages[sample] = Some(value);
                eligible_samples.push(sample);
                differences.push(value[1] - value[0]);
            }
        }
        if differences.len() < min_samples {
            continue;
        }
        let distal = *site_ids.last().unwrap();
        let distal_umis = eligible_samples
            .iter()
            .map(|&sample| {
                mixture_sites[distal].counts[sample][0] + mixture_sites[distal].counts[sample][1]
            })
            .sum::<f64>();
        if distal_umis < min_distal_umis as f64 {
            continue;
        }
        let effect = differences.iter().sum::<f64>() / differences.len() as f64;
        let Some((t, p, p_flip)) = apastats::paired_test(&differences) else {
            continue;
        };
        let concordant = differences
            .iter()
            .filter(|&&difference| difference != 0.0 && difference.signum() == effect.signum())
            .count();
        let lodo_stable = (0..differences.len()).all(|left_out| {
            let mean = differences
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != left_out)
                .map(|(_, value)| value)
                .sum::<f64>()
                / (differences.len() - 1) as f64;
            mean != 0.0 && mean.signum() == effect.signum()
        });
        tests.push(GeneTest {
            gene: gene as u32,
            site_ids: site_ids.clone(),
            eligible_samples,
            usages,
            effect,
            t,
            p,
            p_flip,
            q: 1.0,
            concordant,
            lodo_stable,
            distal_umis: distal_umis.round() as u64,
            reported: false,
        });
    }
    let q = apastats::bh_fdr(&tests.iter().map(|test| test.p).collect::<Vec<_>>());
    for (test, q_value) in tests.iter_mut().zip(q) {
        test.q = q_value;
        test.reported =
            test.q <= 0.05 && test.effect.abs() >= 0.10 && test.concordant >= 6 && test.lodo_stable;
    }
    tests.sort_by(|a, b| a.p.partial_cmp(&b.p).unwrap_or(Ordering::Equal));
    let audit = MixtureAudit {
        candidate_sites,
        recurrent_sites: mixture_sites.len(),
        assigned_umis,
    };
    (mixture_sites, tests, audit)
}

fn write_gene_table(
    path: &Path,
    annotation: &anno::Annotation,
    samples: &[SampleCounts],
    tests: &[GeneTest],
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    write!(
        writer,
        "gene_id\tgene_name\tn_sites\teligible_samples\teffect_B_minus_A\tt\tp\tp_signflip\tq\tconcordant\tlodo_stable\tdistal_umis\treported"
    )?;
    for sample in samples {
        write!(
            writer,
            "\t{}:A\t{}:B\t{}:delta",
            sample.row.sample, sample.row.sample, sample.row.sample
        )?;
    }
    writeln!(writer)?;
    for test in tests {
        let gene = test.gene as usize;
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{:.8}\t{:.8}\t{:.8e}\t{:.8e}\t{:.8e}\t{}\t{}\t{}\t{}",
            annotation.gene_ids[gene],
            annotation.gene_names[gene],
            test.site_ids.len(),
            test.eligible_samples.len(),
            test.effect,
            test.t,
            test.p,
            test.p_flip,
            test.q,
            test.concordant,
            test.lodo_stable as u8,
            test.distal_umis,
            test.reported as u8
        )?;
        for usage in &test.usages {
            if let Some([a, b]) = usage {
                write!(writer, "\t{a:.8}\t{b:.8}\t{:.8}", b - a)?;
            } else {
                write!(writer, "\tNA\tNA\tNA")?;
            }
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_kernel_table(path: &Path, kernel: &KernelAudit) -> Result<u64> {
    let total: u64 = kernel.bins.iter().sum();
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "distance_start\tdistance_end\tumis\tfraction")?;
    for (bin, &count) in kernel.bins.iter().enumerate() {
        let start = KERNEL_MIN_BP + bin as i32 * KERNEL_BIN_BP;
        let fraction = if total == 0 {
            0.0
        } else {
            count as f64 / total as f64
        };
        writeln!(
            writer,
            "{start}\t{}\t{count}\t{fraction:.10}",
            start + KERNEL_BIN_BP
        )?;
    }
    writer.flush()?;
    Ok(total)
}

fn donor_direction(
    tests: &[GeneTest],
    sample_count: usize,
) -> (Vec<Option<f64>>, usize, Option<f64>, bool) {
    let mut donor_medians = Vec::with_capacity(sample_count);
    for sample in 0..sample_count {
        let mut values: Vec<f64> = tests
            .iter()
            .filter_map(|test| test.usages[sample].map(|usage| usage[1] - usage[0]))
            .collect();
        donor_medians.push(median(&mut values));
    }
    let positive = donor_medians
        .iter()
        .filter(|value| value.is_some_and(|value| value > 0.0))
        .count();
    let mut observed: Vec<f64> = donor_medians.iter().flatten().copied().collect();
    let across = median(&mut observed);
    let pass = positive >= 7 && across.is_some_and(|value| value >= 0.02);
    (donor_medians, positive, across, pass)
}

fn write_mixture_site_table(
    path: &Path,
    annotation: &anno::Annotation,
    gene_chrom: &[String],
    gene_rev: &[bool],
    contrast: &[String; 2],
    samples: &[SampleCounts],
    sites: &[MixtureSite],
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    write!(
        writer,
        "site_id\tgene_id\tgene_name\tchrom\tposition\tstrand\trecurrent_samples"
    )?;
    for sample in samples {
        write!(
            writer,
            "\t{}:{}\t{}:{}",
            sample.row.sample, contrast[0], sample.row.sample, contrast[1]
        )?;
    }
    writeln!(writer)?;
    for (site_id, site) in sites.iter().enumerate() {
        let gene = site.gene as usize;
        write!(
            writer,
            "{site_id}\t{}\t{}\t{}\t{}\t{}\t{}",
            annotation.gene_ids[gene],
            annotation.gene_names[gene],
            gene_chrom[gene],
            site.point,
            if gene_rev[gene] { '-' } else { '+' },
            site.recurrent_samples,
        )?;
        for row in &site.counts {
            write!(writer, "\t{:.6}\t{:.6}", row[0], row[1])?;
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

const TRANSCRIPT_END_REPORT_SCHEMA: &str = "gravlax.cohort.transcript-ends.result.v1";
const POLYASITE_MIXTURE_REPORT_SCHEMA: &str = "gravlax.cohort.polyasite-mixture.result.v1";

#[derive(Serialize)]
struct ArtifactPublicationSummary<'a> {
    output_directory: &'a str,
    directory_transactional: bool,
    completion_marker: &'static str,
    report_file_atomic_no_clobber: bool,
}

#[derive(Serialize)]
struct TranscriptEndReportSummary<'a> {
    artifact: ArtifactPublicationSummary<'a>,
    samples: u64,
    recurrent_endpoint_sites: u64,
    high_confidence_endpoint_sites: u64,
    internal_priming_flagged_endpoint_sites: u64,
    tested_endpoint_genes: u64,
    reported_endpoint_genes: u64,
    tested_polyasite_only_genes: u64,
    reported_polyasite_only_genes: u64,
    mixture_candidate_sites: u64,
    mixture_recurrent_sites: u64,
    mixture_assigned_umis: f64,
    tested_mixture_genes: u64,
    reported_mixture_genes: u64,
    fragment_kernel_umis: u64,
}

#[derive(Serialize)]
struct PolyasiteMixtureReportSummary<'a> {
    artifact: ArtifactPublicationSummary<'a>,
    samples: u64,
    external_candidate_sites: u64,
    internal_priming_candidates_dropped: u64,
    compatible_candidate_sites: u64,
    recurrent_sites: u64,
    assigned_umis: f64,
    tested_genes: u64,
    reported_genes: u64,
    fragment_kernel_umis: u64,
}

struct ArtifactRow {
    kind: &'static str,
    path: String,
    bytes: u64,
    records: u64,
}

fn prospective_path(path: &Path, label: &str) -> Result<PathBuf> {
    let file_name = path
        .file_name()
        .with_context(|| format!("{label} must name a file or directory"))?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let parent = std::fs::canonicalize(parent)
        .with_context(|| format!("resolving {label} parent {}", parent.display()))?;
    if !parent.is_dir() {
        bail!("{label} parent is not a directory: {}", parent.display());
    }
    Ok(parent.join(file_name))
}

fn require_uniform_path<'a>(path: &'a Path, label: &str) -> Result<&'a str> {
    path.to_str().with_context(|| {
        format!(
            "uniform reports require UTF-8 {label} paths; {} is not representable",
            path.display()
        )
    })
}

fn preflight_uniform_report(args: &Args) -> Result<bool> {
    match (args.report_format, args.report_output.as_deref()) {
        (None, None) => Ok(false),
        (None, Some(_)) => bail!("--report-output requires --report-format"),
        (Some(_), output) => {
            for (label, path) in [
                ("design", args.design.as_path()),
                ("annotation", args.gtf.as_path()),
                ("genome", args.genome.as_path()),
                ("PolyASite catalogue", args.polyasite.as_path()),
                ("output-directory", args.out_dir.as_path()),
            ] {
                require_uniform_path(path, label)?;
            }
            let out_dir = prospective_path(&args.out_dir, "--out-dir")?;
            if let Some(output) = output {
                require_uniform_path(output, "report-output")?;
                match std::fs::symlink_metadata(output) {
                    Ok(_) => bail!(
                        "refusing to replace existing report output {}",
                        output.display()
                    ),
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                    Err(error) => {
                        return Err(error).with_context(|| {
                            format!("checking report output {}", output.display())
                        });
                    }
                }
                let report = prospective_path(output, "--report-output")?;
                if report == out_dir {
                    bail!("--report-output must differ from --out-dir");
                }
            }
            Ok(true)
        }
    }
}

fn validate_uniform_design_paths(design: &[DesignRow]) -> Result<()> {
    for row in design {
        require_uniform_path(&row.archive, "archive")?;
        require_uniform_path(&row.groups, "group-source")?;
    }
    Ok(())
}

fn load_annotation(
    path: &Path,
    capture_identity: bool,
) -> Result<(anno::Annotation, Option<String>)> {
    if !capture_identity {
        return Ok((anno::Annotation::from_path(path)?, None));
    }
    let annotation_label = require_uniform_path(path, "annotation")?;
    let identity =
        anno::intent::AnnotationIdentity::new("archive-stamped-reference", annotation_label)?;
    let bound = anno::intent::BoundAnnotation::from_path(path, identity)?;
    let (annotation, identity) = bound.into_parts();
    Ok((annotation, identity.digest))
}

fn uniform_context(
    args: &Args,
    samples: &[SampleCounts],
    design_identity: &str,
    annotation_identity: &str,
    genome_identity: &str,
    polyasite_identity: &str,
    include_motif_threshold: bool,
) -> Result<ResultContext> {
    let mut parameters = BTreeMap::new();
    parameters.insert(
        "design_path".into(),
        json!(require_uniform_path(&args.design, "design")?),
    );
    parameters.insert("design_identity".into(), json!(design_identity));
    parameters.insert(
        "genome_path".into(),
        json!(require_uniform_path(&args.genome, "genome")?),
    );
    parameters.insert("genome_identity".into(), json!(genome_identity));
    parameters.insert(
        "polyasite_path".into(),
        json!(require_uniform_path(
            &args.polyasite,
            "PolyASite catalogue"
        )?),
    );
    parameters.insert(
        "polyasite_catalogue_identity".into(),
        json!(polyasite_identity),
    );
    parameters.insert("group_contrast".into(), json!(args.group_contrast));
    parameters.insert("site_gap".into(), json!(args.site_gap));
    parameters.insert("tail_extend".into(), json!(args.tail_extend));
    parameters.insert("min_site_umis".into(), json!(args.min_site_umis));
    parameters.insert("min_site_samples".into(), json!(args.min_site_samples));
    if include_motif_threshold {
        parameters.insert("motif_min_samples".into(), json!(args.motif_min_samples));
    }
    parameters.insert(
        "min_group_gene_umis".into(),
        json!(args.min_group_gene_umis),
    );
    parameters.insert("min_samples".into(), json!(args.min_samples));
    parameters.insert("min_distal_umis".into(), json!(args.min_distal_umis));
    parameters.insert("max_sites".into(), json!(args.max_sites));
    if let Some(seed) = args.shuffle_seed {
        parameters.insert("shuffle_seed".into(), json!(seed));
    }
    let group_sources = samples
        .iter()
        .map(|sample| {
            Ok(json!({
                "sample": sample.row.sample,
                "path": require_uniform_path(&sample.row.groups, "group-source")?,
                "identity": sample.groups_identity.as_deref().context(
                    "uniform report did not capture a group-source identity"
                )?,
            }))
        })
        .collect::<Result<Vec<_>>>()?;
    parameters.insert("group_sources".into(), json!(group_sources));
    if let Some(signature) = samples
        .first()
        .and_then(|sample| sample.reference_signature.as_ref())
    {
        parameters.insert(
            "archive_genome_signature".into(),
            json!(format!("{}:{}", signature.algo, signature.digest)),
        );
    }
    let mut seen_archive_identities = BTreeSet::new();
    let mut archives = Vec::new();
    for sample in samples {
        let identity = sample
            .archive_identity
            .clone()
            .context("uniform report did not capture an archive identity")?;
        if seen_archive_identities.insert(identity.clone()) {
            archives.push(identity);
        }
    }
    let warnings = if samples
        .first()
        .is_some_and(|sample| sample.reference_signature.is_none())
    {
        vec![
            "archives are not stamped with a reference signature; sequence identity was not authenticated"
                .into(),
        ]
    } else {
        Vec::new()
    };
    Ok(ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            archives,
            annotation: Some(require_uniform_path(&args.gtf, "annotation")?.into()),
            annotation_digest: Some(annotation_identity.into()),
            parameters,
            ..Default::default()
        },
        warnings,
    })
}

fn schema(
    id: &str,
    fields: Vec<Field>,
    row_semantics: RowSemantics,
    key: Option<&[&str]>,
) -> std::result::Result<TableSchema, OutputError> {
    let semantics = match key {
        Some(key) => TableSemantics::new(row_semantics).with_key(key.iter().copied()),
        None => TableSemantics::new(row_semantics),
    };
    TableSchema::new(id, fields)?.with_semantics(semantics)
}

fn sample_schema(command: &str) -> std::result::Result<TableSchema, OutputError> {
    schema(
        &format!("gravlax.cohort.{command}.samples.v1"),
        vec![
            Field::new("sample", DataType::String),
            Field::new("condition", DataType::String),
            Field::new("archive_path", DataType::String),
            Field::new("archive_identity", DataType::String),
            Field::new("archive_format_version", DataType::UInt64),
            Field::new("groups_path", DataType::String),
            Field::new("groups_identity", DataType::String),
            Field::new("group_a", DataType::String),
            Field::new("group_b", DataType::String),
            Field::new("selected_group_a_cells", DataType::UInt64),
            Field::new("selected_group_b_cells", DataType::UInt64),
            Field::new("absent_group_a_design_cells", DataType::UInt64),
            Field::new("absent_group_b_design_cells", DataType::UInt64),
            Field::new("endpoint_records", DataType::UInt64),
            Field::new("deduplicated_gene_umis", DataType::UInt64),
            Field::new("outside_terminal_region_records", DataType::UInt64),
            Field::new("ambiguous_terminal_region_records", DataType::UInt64),
            Field::new("archive_chunks", DataType::UInt64),
        ],
        RowSemantics::Sequence,
        Some(&["sample"]),
    )
}

fn artifact_schema(command: &str) -> std::result::Result<TableSchema, OutputError> {
    schema(
        &format!("gravlax.cohort.{command}.artifacts.v1"),
        vec![
            Field::new("artifact_kind", DataType::String),
            Field::new("path", DataType::String),
            Field::new("bytes", DataType::UInt64),
            Field::new("records", DataType::UInt64),
        ],
        RowSemantics::Set,
        Some(&["artifact_kind", "path"]),
    )
}

fn endpoint_site_schema() -> std::result::Result<TableSchema, OutputError> {
    schema(
        "gravlax.cohort.transcript-ends.sites.v1",
        vec![
            Field::new("site_id", DataType::UInt64),
            Field::new("gene_id", DataType::String),
            Field::new("gene_name", DataType::String),
            Field::new("chrom", DataType::String),
            Field::new("start", DataType::UInt64),
            Field::new("end", DataType::UInt64),
            Field::new("endpoint", DataType::UInt64),
            Field::new("strand", DataType::String),
            Field::new("recurrent_samples", DataType::UInt64),
            Field::new("high_confidence", DataType::Boolean),
            Field::new("polyasite_distance", DataType::UInt64).nullable(),
            Field::new("motif", DataType::String).nullable(),
            Field::new("internal_priming", DataType::Boolean),
            Field::new("a20", DataType::UInt64),
            Field::new("arun", DataType::UInt64),
        ],
        RowSemantics::Set,
        Some(&["site_id"]),
    )
}

fn endpoint_count_schema() -> std::result::Result<TableSchema, OutputError> {
    schema(
        "gravlax.cohort.transcript-ends.site-counts.v1",
        vec![
            Field::new("site_id", DataType::UInt64),
            Field::new("sample", DataType::String),
            Field::new("group", DataType::String),
            Field::new("umis", DataType::UInt64),
        ],
        RowSemantics::Set,
        Some(&["site_id", "sample", "group"]),
    )
}

fn mixture_site_schema(command: &str) -> std::result::Result<TableSchema, OutputError> {
    schema(
        &format!("gravlax.cohort.{command}.sites.v1"),
        vec![
            Field::new("site_id", DataType::UInt64),
            Field::new("gene_id", DataType::String),
            Field::new("gene_name", DataType::String),
            Field::new("chrom", DataType::String),
            Field::new("position", DataType::UInt64),
            Field::new("strand", DataType::String),
            Field::new("recurrent_samples", DataType::UInt64),
        ],
        RowSemantics::Set,
        Some(&["site_id"]),
    )
}

fn mixture_count_schema(command: &str) -> std::result::Result<TableSchema, OutputError> {
    schema(
        &format!("gravlax.cohort.{command}.site-counts.v1"),
        vec![
            Field::new("site_id", DataType::UInt64),
            Field::new("sample", DataType::String),
            Field::new("group", DataType::String),
            Field::new("expected_umis", DataType::Float64),
        ],
        RowSemantics::Set,
        Some(&["site_id", "sample", "group"]),
    )
}

fn gene_schema(command: &str) -> std::result::Result<TableSchema, OutputError> {
    schema(
        &format!("gravlax.cohort.{command}.genes.v1"),
        vec![
            Field::new("analysis", DataType::String),
            Field::new("gene_id", DataType::String),
            Field::new("gene_name", DataType::String),
            Field::new("n_sites", DataType::UInt64),
            Field::new("eligible_samples", DataType::UInt64),
            Field::new("effect_b_minus_a", DataType::Float64),
            Field::new("t", DataType::Float64),
            Field::new("p", DataType::Float64),
            Field::new("p_signflip", DataType::Float64),
            Field::new("q", DataType::Float64),
            Field::new("concordant", DataType::UInt64),
            Field::new("lodo_stable", DataType::Boolean),
            Field::new("distal_umis", DataType::UInt64),
            Field::new("reported", DataType::Boolean),
        ],
        RowSemantics::Set,
        Some(&["analysis", "gene_id"]),
    )
}

fn gene_usage_schema(command: &str) -> std::result::Result<TableSchema, OutputError> {
    schema(
        &format!("gravlax.cohort.{command}.gene-usages.v1"),
        vec![
            Field::new("analysis", DataType::String),
            Field::new("gene_id", DataType::String),
            Field::new("sample", DataType::String),
            Field::new("group", DataType::String),
            Field::new("usage", DataType::Float64).nullable(),
        ],
        RowSemantics::Set,
        Some(&["analysis", "gene_id", "sample", "group"]),
    )
}

fn kernel_schema(command: &str) -> std::result::Result<TableSchema, OutputError> {
    let semantics = TableSemantics::new(RowSemantics::Sequence).ordered_by([OrderKey {
        field: "distance_start".into(),
        direction: SortDirection::Ascending,
    }]);
    TableSchema::new(
        format!("gravlax.cohort.{command}.fragment-kernel.v1"),
        vec![
            Field::new("distance_start", DataType::Int64),
            Field::new("distance_end", DataType::Int64),
            Field::new("umis", DataType::UInt64),
            Field::new("fraction", DataType::Float64),
        ],
    )?
    .with_semantics(semantics)
}

fn heldout_schema() -> std::result::Result<TableSchema, OutputError> {
    schema(
        "gravlax.cohort.polyasite-mixture.heldout-kernel.v1",
        vec![
            Field::new("sample", DataType::String),
            Field::new("heldout_unique_candidate_umis", DataType::UInt64),
            Field::new(
                "predictive_gain_over_uniform_bits_per_umi",
                DataType::Float64,
            )
            .nullable(),
        ],
        RowSemantics::Sequence,
        Some(&["sample"]),
    )
}

fn artifact_rows(
    out_dir: &Path,
    reported_out_dir: &str,
    definitions: &[(&'static str, &'static str, u64)],
) -> Result<Vec<ArtifactRow>> {
    definitions
        .iter()
        .map(|&(kind, name, records)| {
            let path = out_dir.join(name);
            let reported_path = Path::new(reported_out_dir).join(name);
            Ok(ArtifactRow {
                kind,
                path: require_uniform_path(&reported_path, "artifact")?.to_owned(),
                bytes: std::fs::metadata(&path)
                    .with_context(|| format!("reading artifact size {}", path.display()))?
                    .len(),
                records,
            })
        })
        .collect()
}

fn write_samples<W: Write>(
    bundle: &mut StreamingBundleWriter<W>,
    command: &str,
    samples: &[SampleCounts],
    contrast: &[String; 2],
) -> std::result::Result<(), OutputError> {
    let schema = sample_schema(command)?;
    bundle.write_table("samples", &schema, None, |rows| {
        for sample in samples {
            rows.write_row_with(|row| {
                row.string(&sample.row.sample)?;
                row.string(&sample.row.condition)?;
                row.string(sample.row.archive.to_str().ok_or_else(|| {
                    OutputError::Sink("archive path ceased to be valid UTF-8".into())
                })?)?;
                row.string(
                    sample
                        .archive_identity
                        .as_deref()
                        .ok_or_else(|| OutputError::Sink("archive identity is absent".into()))?,
                )?;
                row.uint64(
                    sample.archive_format_version.ok_or_else(|| {
                        OutputError::Sink("archive format version is absent".into())
                    })? as u64,
                )?;
                row.string(sample.row.groups.to_str().ok_or_else(|| {
                    OutputError::Sink("groups path ceased to be valid UTF-8".into())
                })?)?;
                row.string(
                    sample
                        .groups_identity
                        .as_deref()
                        .ok_or_else(|| OutputError::Sink("groups identity is absent".into()))?,
                )?;
                row.string(&contrast[0])?;
                row.string(&contrast[1])?;
                row.uint64(sample.selected_cells[0] as u64)?;
                row.uint64(sample.selected_cells[1] as u64)?;
                row.uint64(sample.absent_design_cells[0] as u64)?;
                row.uint64(sample.absent_design_cells[1] as u64)?;
                row.uint64(sample.endpoint_records as u64)?;
                row.uint64(sample.deduplicated_gene_umis as u64)?;
                row.uint64(sample.outside_terminal_region_records as u64)?;
                row.uint64(sample.ambiguous_terminal_region_records as u64)?;
                row.uint64(sample.chunks as u64)?;
                Ok(())
            })?;
        }
        Ok(())
    })
}

fn write_artifacts<W: Write>(
    bundle: &mut StreamingBundleWriter<W>,
    command: &str,
    artifacts: &[ArtifactRow],
) -> std::result::Result<(), OutputError> {
    let schema = artifact_schema(command)?;
    bundle.write_table("artifacts", &schema, None, |rows| {
        for artifact in artifacts {
            rows.write_row_with(|row| {
                row.string(artifact.kind)?;
                row.string(&artifact.path)?;
                row.uint64(artifact.bytes)?;
                row.uint64(artifact.records)?;
                Ok(())
            })?;
        }
        Ok(())
    })
}

fn write_mixture_sites<W: Write>(
    bundle: &mut StreamingBundleWriter<W>,
    table_name: &str,
    command: &str,
    annotation: &anno::Annotation,
    gene_chrom: &[String],
    gene_rev: &[bool],
    contrast: &[String; 2],
    samples: &[SampleCounts],
    sites: &[MixtureSite],
) -> std::result::Result<(), OutputError> {
    let site_schema = mixture_site_schema(command)?;
    bundle.write_table(table_name, &site_schema, None, |rows| {
        for (site_id, site) in sites.iter().enumerate() {
            let gene = site.gene as usize;
            rows.write_row_with(|row| {
                row.uint64(site_id as u64)?;
                row.string(&annotation.gene_ids[gene])?;
                row.string(&annotation.gene_names[gene])?;
                row.string(&gene_chrom[gene])?;
                row.uint64(site.point as u64)?;
                row.string(if gene_rev[gene] { "-" } else { "+" })?;
                row.uint64(site.recurrent_samples as u64)?;
                Ok(())
            })?;
        }
        Ok(())
    })?;
    let count_schema = mixture_count_schema(command)?;
    let count_table_name = format!(
        "{}_counts",
        table_name.strip_suffix('s').unwrap_or(table_name)
    );
    bundle.write_table(&count_table_name, &count_schema, None, |rows| {
        for (site_id, site) in sites.iter().enumerate() {
            for (sample_index, sample) in samples.iter().enumerate() {
                for (group_index, group) in contrast.iter().enumerate() {
                    rows.write_row_with(|row| {
                        row.uint64(site_id as u64)?;
                        row.string(&sample.row.sample)?;
                        row.string(group)?;
                        row.float64(site.counts[sample_index][group_index])?;
                        Ok(())
                    })?;
                }
            }
        }
        Ok(())
    })
}

fn write_genes<W: Write>(
    bundle: &mut StreamingBundleWriter<W>,
    command: &str,
    annotation: &anno::Annotation,
    contrast: &[String; 2],
    samples: &[SampleCounts],
    analyses: &[(&str, &[GeneTest])],
) -> std::result::Result<(), OutputError> {
    let schema = gene_schema(command)?;
    bundle.write_table("genes", &schema, None, |rows| {
        for &(analysis, tests) in analyses {
            for test in tests {
                let gene = test.gene as usize;
                rows.write_row_with(|row| {
                    row.string(analysis)?;
                    row.string(&annotation.gene_ids[gene])?;
                    row.string(&annotation.gene_names[gene])?;
                    row.uint64(test.site_ids.len() as u64)?;
                    row.uint64(test.eligible_samples.len() as u64)?;
                    row.float64(test.effect)?;
                    row.float64(test.t)?;
                    row.float64(test.p)?;
                    row.float64(test.p_flip)?;
                    row.float64(test.q)?;
                    row.uint64(test.concordant as u64)?;
                    row.boolean(test.lodo_stable)?;
                    row.uint64(test.distal_umis)?;
                    row.boolean(test.reported)?;
                    Ok(())
                })?;
            }
        }
        Ok(())
    })?;
    let usage_schema = gene_usage_schema(command)?;
    bundle.write_table("gene_usages", &usage_schema, None, |rows| {
        for &(analysis, tests) in analyses {
            for test in tests {
                let gene_id = &annotation.gene_ids[test.gene as usize];
                for (sample_index, sample) in samples.iter().enumerate() {
                    for (group_index, group) in contrast.iter().enumerate() {
                        rows.write_row_with(|row| {
                            row.string(analysis)?;
                            row.string(gene_id)?;
                            row.string(&sample.row.sample)?;
                            row.string(group)?;
                            if let Some(usage) = test.usages[sample_index] {
                                row.float64(usage[group_index])?;
                            } else {
                                row.null()?;
                            }
                            Ok(())
                        })?;
                    }
                }
            }
        }
        Ok(())
    })
}

fn write_kernel<W: Write>(
    bundle: &mut StreamingBundleWriter<W>,
    command: &str,
    kernel: &KernelAudit,
) -> std::result::Result<(), OutputError> {
    let schema = kernel_schema(command)?;
    let total: u64 = kernel.bins.iter().sum();
    bundle.write_table("fragment_kernel", &schema, None, |rows| {
        for (bin, &count) in kernel.bins.iter().enumerate() {
            let start = KERNEL_MIN_BP + bin as i32 * KERNEL_BIN_BP;
            let fraction = if total == 0 {
                0.0
            } else {
                count as f64 / total as f64
            };
            rows.write_row_with(|row| {
                row.int64(start as i64)?;
                row.int64((start + KERNEL_BIN_BP) as i64)?;
                row.uint64(count)?;
                row.float64(fraction)?;
                Ok(())
            })?;
        }
        Ok(())
    })
}

#[allow(clippy::too_many_arguments)]
fn write_transcript_end_report<W: Write>(
    writer: W,
    format: OutputFormat,
    context: &ResultContext,
    summary: &TranscriptEndReportSummary<'_>,
    annotation: &anno::Annotation,
    contrast: &[String; 2],
    gene_rev: &[bool],
    gene_chrom: &[String],
    samples: &[SampleCounts],
    sites: &[Site],
    tests: &[GeneTest],
    polyasite_tests: &[GeneTest],
    mixture_sites: &[MixtureSite],
    mixture_tests: &[GeneTest],
    kernel: &KernelAudit,
    artifacts: &[ArtifactRow],
) -> std::result::Result<(), OutputError> {
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        TRANSCRIPT_END_REPORT_SCHEMA,
        format,
        context,
        summary,
    )?;
    write_samples(&mut bundle, "transcript-ends", samples, contrast)?;
    let site_schema = endpoint_site_schema()?;
    bundle.write_table("sites", &site_schema, None, |rows| {
        for (site_id, site) in sites.iter().enumerate() {
            let gene = site.gene as usize;
            rows.write_row_with(|row| {
                row.uint64(site_id as u64)?;
                row.string(&annotation.gene_ids[gene])?;
                row.string(&annotation.gene_names[gene])?;
                row.string(&gene_chrom[gene])?;
                row.uint64(site.lo as u64)?;
                row.uint64(site.hi.saturating_add(1) as u64)?;
                row.uint64(site.representative as u64)?;
                row.string(if gene_rev[gene] { "-" } else { "+" })?;
                row.uint64(site.recurrent_samples as u64)?;
                row.boolean(site.high_confidence)?;
                if let Some(distance) = site.polyasite_distance {
                    row.uint64(distance as u64)?;
                } else {
                    row.null()?;
                }
                if let Some(motif) = site.motif.as_deref() {
                    row.string(motif)?;
                } else {
                    row.null()?;
                }
                row.boolean(site.ip)?;
                row.uint64(site.a20 as u64)?;
                row.uint64(site.arun as u64)?;
                Ok(())
            })?;
        }
        Ok(())
    })?;
    let count_schema = endpoint_count_schema()?;
    bundle.write_table("site_counts", &count_schema, None, |rows| {
        for (site_id, site) in sites.iter().enumerate() {
            for (sample_index, sample) in samples.iter().enumerate() {
                for (group_index, group) in contrast.iter().enumerate() {
                    rows.write_row_with(|row| {
                        row.uint64(site_id as u64)?;
                        row.string(&sample.row.sample)?;
                        row.string(group)?;
                        row.uint64(site.counts[sample_index][group_index])?;
                        Ok(())
                    })?;
                }
            }
        }
        Ok(())
    })?;
    write_mixture_sites(
        &mut bundle,
        "mixture_sites",
        "transcript-ends.polyasite-mixture",
        annotation,
        gene_chrom,
        gene_rev,
        contrast,
        samples,
        mixture_sites,
    )?;
    write_genes(
        &mut bundle,
        "transcript-ends",
        annotation,
        contrast,
        samples,
        &[
            ("endpoint", tests),
            ("polyasite_only", polyasite_tests),
            ("polyasite_mixture", mixture_tests),
        ],
    )?;
    write_kernel(&mut bundle, "transcript-ends", kernel)?;
    write_artifacts(&mut bundle, "transcript-ends", artifacts)?;
    bundle.finish()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_polyasite_mixture_report<W: Write>(
    writer: W,
    format: OutputFormat,
    context: &ResultContext,
    summary: &PolyasiteMixtureReportSummary<'_>,
    annotation: &anno::Annotation,
    contrast: &[String; 2],
    gene_rev: &[bool],
    gene_chrom: &[String],
    samples: &[SampleCounts],
    sites: &[MixtureSite],
    tests: &[GeneTest],
    kernel: &KernelAudit,
    per_sample_kernels: &[KernelAudit],
    heldout_bits: &[f64],
    artifacts: &[ArtifactRow],
) -> std::result::Result<(), OutputError> {
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        POLYASITE_MIXTURE_REPORT_SCHEMA,
        format,
        context,
        summary,
    )?;
    write_samples(&mut bundle, "polyasite-mixture", samples, contrast)?;
    write_mixture_sites(
        &mut bundle,
        "sites",
        "polyasite-mixture",
        annotation,
        gene_chrom,
        gene_rev,
        contrast,
        samples,
        sites,
    )?;
    write_genes(
        &mut bundle,
        "polyasite-mixture",
        annotation,
        contrast,
        samples,
        &[("polyasite_mixture", tests)],
    )?;
    write_kernel(&mut bundle, "polyasite-mixture", kernel)?;
    let heldout_schema = heldout_schema()?;
    bundle.write_table("heldout_kernel", &heldout_schema, None, |rows| {
        for ((sample, audit), &bits) in samples.iter().zip(per_sample_kernels).zip(heldout_bits) {
            rows.write_row_with(|row| {
                row.string(&sample.row.sample)?;
                row.uint64(audit.unique_candidate_umis)?;
                if bits.is_finite() {
                    row.float64(bits)?;
                } else {
                    row.null()?;
                }
                Ok(())
            })?;
        }
        Ok(())
    })?;
    write_artifacts(&mut bundle, "polyasite-mixture", artifacts)?;
    bundle.finish()?;
    Ok(())
}

fn publish_uniform_report<F>(
    format: UniformEndReportFormat,
    output: Option<&Path>,
    render: F,
) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> std::result::Result<(), OutputError>,
{
    if let Some(output) = output {
        let mut render = Some(render);
        let outcome = publish_file_no_clobber(output, Durability::Flush, |writer| {
            render.take().expect("uniform report renderer called once")(&mut *writer)
        })?;
        for warning in outcome.warnings {
            eprintln!("warning: {warning}");
        }
    } else {
        let stdout = std::io::stdout();
        let mut lock = stdout.lock();
        render(&mut lock)?;
    }
    let _ = format;
    Ok(())
}

fn write_outputs(
    args: &Args,
    annotation: &anno::Annotation,
    contrast: &[String; 2],
    gene_rev: &[bool],
    gene_chrom: &[String],
    samples: &[SampleCounts],
    sites: &[Site],
    tests: &[GeneTest],
    polyasite_tests: &[GeneTest],
    mixture_sites: &[MixtureSite],
    mixture_tests: &[GeneTest],
    mixture_audit: &MixtureAudit,
    external_candidate_count: usize,
    external_ip_dropped: usize,
    kernel: &KernelAudit,
    elapsed: f64,
) -> Result<serde_json::Value> {
    let sites_path = args.out_dir.join("sites.tsv");
    let mut site_writer = BufWriter::new(File::create(&sites_path)?);
    write!(
        site_writer,
        "site_id\tgene_id\tgene_name\tchrom\tstart\tend\tendpoint\tstrand\trecurrent_samples\thigh_confidence\tpolyasite_distance\tmotif\tip\ta20\tarun"
    )?;
    for sample in samples {
        write!(
            site_writer,
            "\t{}:{}\t{}:{}",
            sample.row.sample, contrast[0], sample.row.sample, contrast[1]
        )?;
    }
    writeln!(site_writer)?;
    for (site_id, site) in sites.iter().enumerate() {
        let gene = site.gene as usize;
        write!(
            site_writer,
            "{site_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            annotation.gene_ids[gene],
            annotation.gene_names[gene],
            gene_chrom[gene],
            site.lo,
            site.hi.saturating_add(1),
            site.representative,
            if gene_rev[gene] { '-' } else { '+' },
            site.recurrent_samples,
            site.high_confidence as u8,
            site.polyasite_distance
                .map(|value| value.to_string())
                .unwrap_or_else(|| "NA".into()),
            site.motif.as_deref().unwrap_or("NA"),
            site.ip as u8,
            site.a20,
            site.arun
        )?;
        for row in &site.counts {
            write!(site_writer, "\t{}\t{}", row[0], row[1])?;
        }
        writeln!(site_writer)?;
    }
    site_writer.flush()?;

    write_gene_table(&args.out_dir.join("genes.tsv"), annotation, samples, tests)?;
    write_gene_table(
        &args.out_dir.join("genes.polyasite.tsv"),
        annotation,
        samples,
        polyasite_tests,
    )?;
    write_mixture_site_table(
        &args.out_dir.join("polyasite-mixture-sites.tsv"),
        annotation,
        gene_chrom,
        gene_rev,
        contrast,
        samples,
        mixture_sites,
    )?;
    write_gene_table(
        &args.out_dir.join("polyasite-mixture-genes.tsv"),
        annotation,
        samples,
        mixture_tests,
    )?;
    let kernel_total = write_kernel_table(&args.out_dir.join("fragment-kernel.tsv"), kernel)?;

    let non_ip: Vec<&Site> = sites.iter().filter(|site| !site.ip).collect();
    let flagged: Vec<&Site> = sites.iter().filter(|site| site.ip).collect();
    let match_rate = |values: &[&Site]| -> Option<f64> {
        (!values.is_empty()).then(|| {
            values
                .iter()
                .filter(|site| {
                    site.polyasite_distance
                        .is_some_and(|distance| distance <= 50)
                })
                .count() as f64
                / values.len() as f64
        })
    };
    let (donor_medians, positive_donors, across_donor_median, global_pass) =
        donor_direction(tests, samples.len());
    let (polyasite_medians, polyasite_positive, polyasite_across, polyasite_global_pass) =
        donor_direction(polyasite_tests, samples.len());
    let (mixture_medians, mixture_positive, mixture_across, mixture_global_pass) =
        donor_direction(mixture_tests, samples.len());
    let non_ip_rate = match_rate(&non_ip);
    let flagged_rate = match_rate(&flagged);
    let accuracy_separation = non_ip_rate
        .zip(flagged_rate)
        .is_some_and(|(pass, ip)| pass > ip);
    let summary = json!({
        "schema": "gravlax.cohort.transcript-ends.v1",
        "semantics": {
            "replicate": "one unique donor/archive design row",
            "cells_or_molecules_are_replicates": false,
            "endpoint": "transcript-oriented terminal aligned base, not an asserted cleavage site",
            "gene_assignment": "unique same-strand union of annotated terminal exons, each with transcript-oriented tail extension",
            "molecule_count": "deduplicated once per gene and archive UMI class",
            "site_construction": "one common single-link endpoint catalogue across all donors",
            "effect": format!("{} minus {} distal usage index", contrast[1], contrast[0]),
            "shuffle": args.shuffle_seed.map(|seed| json!({"within_donor": true, "seed": seed})),
        },
        "design": {
            "path": args.design,
            "reference_digest": samples.first().and_then(|sample| sample.reference_signature.as_ref().map(|signature| signature.digest.clone())),
            "samples": samples.iter().map(|sample| json!({
                "sample": sample.row.sample, "condition": sample.row.condition,
                "archive": sample.row.archive_label, "groups": sample.row.groups_label,
                "selected_cells": {contrast[0].clone(): sample.selected_cells[0], contrast[1].clone(): sample.selected_cells[1]},
                "absent_design_cells": {contrast[0].clone(): sample.absent_design_cells[0], contrast[1].clone(): sample.absent_design_cells[1]},
                "endpoint_records": sample.endpoint_records,
                "deduplicated_gene_umis": sample.deduplicated_gene_umis,
                "outside_terminal_region_records": sample.outside_terminal_region_records,
                "ambiguous_terminal_region_records": sample.ambiguous_terminal_region_records,
                "chunks": sample.chunks,
            })).collect::<Vec<_>>(),
        },
        "thresholds": {
            "site_gap": args.site_gap, "tail_extend": args.tail_extend,
            "min_site_umis": args.min_site_umis, "min_site_samples": args.min_site_samples,
            "motif_min_samples": args.motif_min_samples,
            "min_group_gene_umis": args.min_group_gene_umis,
            "min_samples": args.min_samples, "min_distal_umis": args.min_distal_umis,
            "max_sites": args.max_sites,
        },
        "catalogue": {
            "recurrent_sites": sites.len(),
            "high_confidence_sites": sites.iter().filter(|site| site.high_confidence).count(),
            "internal_priming_flagged_sites": flagged.len(),
            "non_ip_polyasite_50bp_rate": non_ip_rate,
            "ip_flagged_polyasite_50bp_rate": flagged_rate,
            "accuracy_separation_pass": accuracy_separation,
        },
        "protocol_kernel": {
            "role": "label-blind audit for catalogue-constrained 10x 3-prime fragment deconvolution",
            "external_candidates_in_unique_terminal_regions": external_candidate_count,
            "internal_priming_candidates_dropped": external_ip_dropped,
            "retained_external_candidates": external_candidate_count - external_ip_dropped,
            "distance_definition": "transcript-upstream distance from aligned fragment endpoint to PolyASite candidate",
            "distance_range_bp": [KERNEL_MIN_BP, KERNEL_MAX_BP],
            "bin_bp": KERNEL_BIN_BP,
            "no_candidate_umis": kernel.no_candidate_umis,
            "unique_candidate_umis": kernel.unique_candidate_umis,
            "multiple_candidate_umis": kernel.multiple_candidate_umis,
            "kernel_umis": kernel_total,
        },
        "inference": {
            "tested_genes": tests.len(),
            "reported_genes": tests.iter().filter(|test| test.reported).count(),
            "test": "two-sided paired donor Student t test",
            "calibration": "exact paired sign-flip p-value",
            "multiplicity": "Benjamini-Hochberg over tested genes",
            "donor_median_effects": samples.iter().enumerate().map(|(index, sample)| json!({
                "sample": sample.row.sample, "condition": sample.row.condition,
                "median_effect": donor_medians[index],
            })).collect::<Vec<_>>(),
            "positive_donor_medians": positive_donors,
            "across_donor_median": across_donor_median,
            "registered_global_directional_pass": global_pass,
            "top_reported": tests.iter().filter(|test| test.reported).take(100).map(|test| json!({
                "gene_id": annotation.gene_ids[test.gene as usize],
                "gene_name": annotation.gene_names[test.gene as usize],
                "effect": test.effect, "q": test.q, "p_signflip": test.p_flip,
                "eligible_samples": test.eligible_samples.len(), "concordant": test.concordant,
            })).collect::<Vec<_>>(),
        },
        "polyasite_only_inference": {
            "role": "endpoint-cluster/external-catalogue audit retained for traceability; fragmented-library mixture inference is preferred",
            "tested_genes": polyasite_tests.len(),
            "reported_genes": polyasite_tests.iter().filter(|test| test.reported).count(),
            "donor_median_effects": samples.iter().enumerate().map(|(index, sample)| json!({
                "sample": sample.row.sample, "condition": sample.row.condition,
                "median_effect": polyasite_medians[index],
            })).collect::<Vec<_>>(),
            "positive_donor_medians": polyasite_positive,
            "across_donor_median": polyasite_across,
            "registered_global_directional_threshold_pass": polyasite_global_pass,
            "top_reported": polyasite_tests.iter().filter(|test| test.reported).take(100).map(|test| json!({
                "gene_id": annotation.gene_ids[test.gene as usize],
                "gene_name": annotation.gene_names[test.gene as usize],
                "effect": test.effect, "q": test.q, "p_signflip": test.p_flip,
                "eligible_samples": test.eligible_samples.len(), "concordant": test.concordant,
            })).collect::<Vec<_>>(),
        },
        "polyasite_mixture_inference": {
            "role": "protocol-aware preferred analysis: empirical fragment kernel plus PolyASite-constrained EM",
            "candidate_sites_with_compatible_evidence": mixture_audit.candidate_sites,
            "recurrent_sites": mixture_audit.recurrent_sites,
            "assigned_umis": mixture_audit.assigned_umis,
            "tested_genes": mixture_tests.len(),
            "reported_genes": mixture_tests.iter().filter(|test| test.reported).count(),
            "donor_median_effects": samples.iter().enumerate().map(|(index, sample)| json!({
                "sample": sample.row.sample, "condition": sample.row.condition,
                "median_effect": mixture_medians[index],
            })).collect::<Vec<_>>(),
            "positive_donor_medians": mixture_positive,
            "across_donor_median": mixture_across,
            "registered_global_directional_threshold_pass": mixture_global_pass,
            "top_reported": mixture_tests.iter().filter(|test| test.reported).take(100).map(|test| json!({
                "gene_id": annotation.gene_ids[test.gene as usize],
                "gene_name": annotation.gene_names[test.gene as usize],
                "effect": test.effect, "q": test.q, "p_signflip": test.p_flip,
                "eligible_samples": test.eligible_samples.len(), "concordant": test.concordant,
            })).collect::<Vec<_>>(),
        },
        "performance": {"wall_seconds": elapsed},
        "outputs": {"sites": "sites.tsv", "genes": "genes.tsv", "polyasite_genes": "genes.polyasite.tsv", "fragment_kernel": "fragment-kernel.tsv", "polyasite_mixture_sites": "polyasite-mixture-sites.tsv", "polyasite_mixture_genes": "polyasite-mixture-genes.tsv"},
    });
    std::fs::write(
        args.out_dir.join("summary.json"),
        serde_json::to_string_pretty(&summary)? + "\n",
    )?;
    Ok(summary)
}

pub fn run(args: Args) -> Result<()> {
    let started = std::time::Instant::now();
    for (name, value) in [
        ("--site-gap", args.site_gap as usize),
        ("--min-site-umis", args.min_site_umis as usize),
        ("--min-site-samples", args.min_site_samples),
        ("--motif-min-samples", args.motif_min_samples),
        ("--min-group-gene-umis", args.min_group_gene_umis as usize),
        ("--min-samples", args.min_samples),
        ("--min-distal-umis", args.min_distal_umis as usize),
        ("--max-sites", args.max_sites),
    ] {
        if value == 0 {
            bail!("{name} must be at least 1");
        }
    }
    if args.out_dir.exists() {
        bail!(
            "refusing to overwrite existing --out-dir {}",
            args.out_dir.display()
        );
    }
    let uniform_report = preflight_uniform_report(&args)?;
    let (design, design_identity) = parse_design(&args.design, uniform_report)?;
    if uniform_report {
        validate_uniform_design_paths(&design)?;
    }
    if args.min_site_samples > design.len() || args.min_samples > design.len() {
        bail!(
            "sample thresholds cannot exceed the {} design rows",
            design.len()
        );
    }
    let contrast = parse_contrast(&args.group_contrast)?;
    let (annotation, annotation_identity) = load_annotation(&args.gtf, uniform_report)?;
    let gene_index = Arc::new(build_gene_index(&annotation, args.tail_extend)?);
    let samples = reduce_design_bounded(
        &design,
        &gene_index,
        &contrast,
        args.shuffle_seed,
        uniform_report,
    )?;
    let signatures: BTreeSet<&str> = samples
        .iter()
        .filter_map(|sample| {
            sample
                .reference_signature
                .as_ref()
                .map(|signature| signature.digest.as_str())
        })
        .collect();
    if signatures.len() > 1
        || (signatures.len() == 1
            && samples
                .iter()
                .any(|sample| sample.reference_signature.is_none()))
    {
        bail!("transcript-end samples have incompatible or mixed stamped reference identities");
    }
    let reference_signature = samples
        .first()
        .and_then(|sample| sample.reference_signature.as_ref());
    let mut sites = site_catalogue(
        &samples,
        annotation.gene_ids.len(),
        &gene_index.rev,
        args.site_gap,
        args.min_site_umis,
        args.min_site_samples,
        args.max_sites,
    )?;
    let (polyasite, polyasite_identity) = load_polyasite(&args.polyasite, uniform_report)?;
    let mut external_candidates = external_catalogue_by_gene(
        &gene_index,
        &polyasite,
        annotation.gene_ids.len(),
        args.site_gap,
    );
    let external_candidate_count = external_candidates.iter().map(Vec::len).sum();
    let (external_ip_dropped, genome_identity) = remove_internal_priming_candidates(
        &mut external_candidates,
        &gene_index,
        &args.genome,
        reference_signature,
        uniform_report,
    )?;
    let (kernel, _per_sample_kernels, sample_kernels, _heldout_bits) =
        cross_fitted_kernels(&samples, &external_candidates, &gene_index.rev);
    let (mixture_sites, mixture_tests, mixture_audit) = fit_polyasite_mixture(
        &samples,
        &external_candidates,
        &gene_index.rev,
        &sample_kernels,
        args.min_site_umis,
        args.min_site_samples,
        args.min_group_gene_umis,
        args.min_samples,
        args.min_distal_umis,
    );
    let scoring_genome_identity = score_sites(
        &mut sites,
        &annotation,
        &gene_index.rev,
        &args.genome,
        &polyasite,
        reference_signature,
        args.motif_min_samples,
        uniform_report,
    )?;
    if genome_identity != scoring_genome_identity {
        bail!("--genome changed between transcript-end sequence passes");
    }
    let tests = fit_genes(
        &sites,
        annotation.gene_ids.len(),
        &gene_index.rev,
        args.min_group_gene_umis,
        args.min_samples,
        args.min_distal_umis,
        registered_site,
    );
    let polyasite_tests = fit_genes(
        &sites,
        annotation.gene_ids.len(),
        &gene_index.rev,
        args.min_group_gene_umis,
        args.min_samples,
        args.min_distal_umis,
        polyasite_site,
    );
    std::fs::create_dir(&args.out_dir)?;
    let summary = write_outputs(
        &args,
        &annotation,
        &contrast,
        &gene_index.rev,
        &gene_index.chrom,
        &samples,
        &sites,
        &tests,
        &polyasite_tests,
        &mixture_sites,
        &mixture_tests,
        &mixture_audit,
        external_candidate_count,
        external_ip_dropped,
        &kernel,
        started.elapsed().as_secs_f64(),
    )?;
    if let Some(format) = args.report_format {
        let design_identity = design_identity
            .as_deref()
            .context("uniform report did not capture a design identity")?;
        let annotation_identity = annotation_identity
            .as_deref()
            .context("uniform report did not capture an annotation identity")?;
        let genome_identity = genome_identity
            .as_deref()
            .context("uniform report did not capture a genome identity")?;
        let polyasite_identity = polyasite_identity
            .as_deref()
            .context("uniform report did not capture a PolyASite catalogue identity")?;
        let context = uniform_context(
            &args,
            &samples,
            design_identity,
            annotation_identity,
            genome_identity,
            polyasite_identity,
            true,
        )?;
        let reported_out_dir = reported_output_path(&args.out_dir)?;
        let artifacts = artifact_rows(
            &args.out_dir,
            &reported_out_dir,
            &[
                ("endpoint-sites", "sites.tsv", sites.len() as u64),
                ("endpoint-genes", "genes.tsv", tests.len() as u64),
                (
                    "polyasite-only-genes",
                    "genes.polyasite.tsv",
                    polyasite_tests.len() as u64,
                ),
                (
                    "polyasite-mixture-sites",
                    "polyasite-mixture-sites.tsv",
                    mixture_sites.len() as u64,
                ),
                (
                    "polyasite-mixture-genes",
                    "polyasite-mixture-genes.tsv",
                    mixture_tests.len() as u64,
                ),
                (
                    "fragment-kernel",
                    "fragment-kernel.tsv",
                    kernel.bins.len() as u64,
                ),
                ("completion-summary", "summary.json", 1),
            ],
        )?;
        let report_summary = TranscriptEndReportSummary {
            artifact: ArtifactPublicationSummary {
                output_directory: &reported_out_dir,
                directory_transactional: false,
                completion_marker: "summary.json",
                report_file_atomic_no_clobber: args.report_output.is_some(),
            },
            samples: samples.len() as u64,
            recurrent_endpoint_sites: sites.len() as u64,
            high_confidence_endpoint_sites: sites.iter().filter(|site| site.high_confidence).count()
                as u64,
            internal_priming_flagged_endpoint_sites: sites.iter().filter(|site| site.ip).count()
                as u64,
            tested_endpoint_genes: tests.len() as u64,
            reported_endpoint_genes: tests.iter().filter(|test| test.reported).count() as u64,
            tested_polyasite_only_genes: polyasite_tests.len() as u64,
            reported_polyasite_only_genes: polyasite_tests
                .iter()
                .filter(|test| test.reported)
                .count() as u64,
            mixture_candidate_sites: mixture_audit.candidate_sites as u64,
            mixture_recurrent_sites: mixture_audit.recurrent_sites as u64,
            mixture_assigned_umis: mixture_audit.assigned_umis,
            tested_mixture_genes: mixture_tests.len() as u64,
            reported_mixture_genes: mixture_tests.iter().filter(|test| test.reported).count()
                as u64,
            fragment_kernel_umis: kernel.bins.iter().sum(),
        };
        publish_uniform_report(format, args.report_output.as_deref(), |writer| {
            write_transcript_end_report(
                writer,
                format.into(),
                &context,
                &report_summary,
                &annotation,
                &contrast,
                &gene_index.rev,
                &gene_index.chrom,
                &samples,
                &sites,
                &tests,
                &polyasite_tests,
                &mixture_sites,
                &mixture_tests,
                &kernel,
                &artifacts,
            )
        })?;
    } else {
        println!("{}", serde_json::to_string_pretty(&summary)?);
    }
    eprintln!(
        "cohort transcript-ends: {} recurrent sites, {} high-confidence, {} tested genes in {:.2}s",
        sites.len(),
        sites.iter().filter(|site| site.high_confidence).count(),
        tests.len(),
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}

/// Production analysis for fragmented 3'-tag libraries. Unlike `transcript-ends`, this command
/// never interprets a fragment boundary as a cleavage coordinate. It learns a label-blind
/// fragment-distance kernel and deconvolves only orthogonally catalogued, non-internal-priming
/// PolyASite candidates.
pub fn run_polyasite_mixture(args: MixtureArgs) -> Result<()> {
    let args = Args::from(args);
    let started = std::time::Instant::now();
    for (name, value) in [
        ("--site-gap", args.site_gap as usize),
        ("--min-site-umis", args.min_site_umis as usize),
        ("--min-site-samples", args.min_site_samples),
        ("--min-group-gene-umis", args.min_group_gene_umis as usize),
        ("--min-samples", args.min_samples),
        ("--min-distal-umis", args.min_distal_umis as usize),
        ("--max-sites", args.max_sites),
    ] {
        if value == 0 {
            bail!("{name} must be at least 1");
        }
    }
    if args.out_dir.exists() {
        bail!(
            "refusing to overwrite existing --out-dir {}",
            args.out_dir.display()
        );
    }
    let uniform_report = preflight_uniform_report(&args)?;
    let (design, design_identity) = parse_design(&args.design, uniform_report)?;
    if uniform_report {
        validate_uniform_design_paths(&design)?;
    }
    if args.min_site_samples > design.len() || args.min_samples > design.len() {
        bail!(
            "sample thresholds cannot exceed the {} design rows",
            design.len()
        );
    }
    let contrast = parse_contrast(&args.group_contrast)?;
    let (annotation, annotation_identity) = load_annotation(&args.gtf, uniform_report)?;
    let gene_index = Arc::new(build_gene_index(&annotation, args.tail_extend)?);
    let samples = reduce_design_bounded(
        &design,
        &gene_index,
        &contrast,
        args.shuffle_seed,
        uniform_report,
    )?;
    let signatures: BTreeSet<&str> = samples
        .iter()
        .filter_map(|sample| {
            sample
                .reference_signature
                .as_ref()
                .map(|signature| signature.digest.as_str())
        })
        .collect();
    if signatures.len() > 1
        || (signatures.len() == 1
            && samples
                .iter()
                .any(|sample| sample.reference_signature.is_none()))
    {
        bail!("polyasite-mixture samples have incompatible or mixed stamped reference identities");
    }
    let reference_signature = samples
        .first()
        .and_then(|sample| sample.reference_signature.as_ref());
    let (polyasite, polyasite_identity) = load_polyasite(&args.polyasite, uniform_report)?;
    let mut candidates = external_catalogue_by_gene(
        &gene_index,
        &polyasite,
        annotation.gene_ids.len(),
        args.site_gap,
    );
    let external_candidate_count: usize = candidates.iter().map(Vec::len).sum();
    if external_candidate_count > args.max_sites {
        bail!(
            "external candidate catalogue has {external_candidate_count} sites, exceeding --max-sites {}; output was not truncated",
            args.max_sites
        );
    }
    let (external_ip_dropped, genome_identity) = remove_internal_priming_candidates(
        &mut candidates,
        &gene_index,
        &args.genome,
        reference_signature,
        uniform_report,
    )?;
    let (kernel, per_sample_kernels, sample_kernels, heldout_bits) =
        cross_fitted_kernels(&samples, &candidates, &gene_index.rev);
    if kernel.unique_candidate_umis == 0 {
        bail!("no unambiguous molecules are available to learn the fragment kernel");
    }
    let (sites, tests, audit) = fit_polyasite_mixture(
        &samples,
        &candidates,
        &gene_index.rev,
        &sample_kernels,
        args.min_site_umis,
        args.min_site_samples,
        args.min_group_gene_umis,
        args.min_samples,
        args.min_distal_umis,
    );
    if audit.recurrent_sites > args.max_sites {
        bail!(
            "recurrent mixture catalogue has {} sites, exceeding --max-sites {}; output was not truncated",
            audit.recurrent_sites,
            args.max_sites
        );
    }
    std::fs::create_dir(&args.out_dir)?;
    write_mixture_site_table(
        &args.out_dir.join("sites.tsv"),
        &annotation,
        &gene_index.chrom,
        &gene_index.rev,
        &contrast,
        &samples,
        &sites,
    )?;
    write_gene_table(
        &args.out_dir.join("genes.tsv"),
        &annotation,
        &samples,
        &tests,
    )?;
    let kernel_total = write_kernel_table(&args.out_dir.join("fragment-kernel.tsv"), &kernel)?;
    let (donor_medians, positive_donors, across_donor_median, global_pass) =
        donor_direction(&tests, samples.len());
    let peak_bin = kernel
        .bins
        .iter()
        .enumerate()
        .max_by_key(|(_, count)| *count)
        .map(|(index, _)| KERNEL_MIN_BP + index as i32 * KERNEL_BIN_BP);
    let summary = json!({
        "schema": "gravlax.cohort.polyasite-mixture.v1",
        "semantics": {
            "replicate": "one unique donor/archive design row",
            "cells_or_molecules_are_replicates": false,
            "site_identity": "PolyASite candidate in a unique same-strand annotated terminal region; never inferred from a raw fragment endpoint",
            "gene_assignment": "unique same-strand union of annotated terminal exons, each with transcript-oriented tail extension",
            "molecule_count": "deduplicated once per gene and archive UMI class",
            "kernel": "group-label-blind empirical transcript-upstream fragment distance learned from exactly-one-candidate molecules; each donor is deconvolved with a leave-one-donor-out kernel",
            "deconvolution": "50-iteration maximum EM, tolerance 1e-8, fixed smoothed empirical kernel",
            "effect": format!("{} minus {} distal PolyASite usage index", contrast[1], contrast[0]),
            "shuffle": args.shuffle_seed.map(|seed| json!({"within_donor": true, "seed": seed})),
        },
        "design": {
            "path": args.design,
            "reference_digest": samples.first().and_then(|sample| sample.reference_signature.as_ref().map(|signature| signature.digest.clone())),
            "samples": samples.iter().map(|sample| json!({
                "sample": sample.row.sample, "condition": sample.row.condition,
                "archive": sample.row.archive_label, "groups": sample.row.groups_label,
                "selected_cells": {contrast[0].clone(): sample.selected_cells[0], contrast[1].clone(): sample.selected_cells[1]},
                "absent_design_cells": {contrast[0].clone(): sample.absent_design_cells[0], contrast[1].clone(): sample.absent_design_cells[1]},
                "endpoint_records": sample.endpoint_records,
                "deduplicated_gene_umis": sample.deduplicated_gene_umis,
                "outside_terminal_region_records": sample.outside_terminal_region_records,
                "ambiguous_terminal_region_records": sample.ambiguous_terminal_region_records,
                "chunks": sample.chunks,
            })).collect::<Vec<_>>(),
        },
        "thresholds": {
            "candidate_merge_gap": args.site_gap, "tail_extend": args.tail_extend,
            "min_site_umis": args.min_site_umis, "min_site_samples": args.min_site_samples,
            "min_group_gene_umis": args.min_group_gene_umis,
            "min_samples": args.min_samples, "min_distal_umis": args.min_distal_umis,
            "max_sites": args.max_sites,
        },
        "accuracy": {
            "external_candidates_in_unique_terminal_regions": external_candidate_count,
            "internal_priming_candidates_dropped": external_ip_dropped,
            "retained_external_candidates": external_candidate_count - external_ip_dropped,
            "kernel_distance_range_bp": [KERNEL_MIN_BP, KERNEL_MAX_BP],
            "kernel_bin_bp": KERNEL_BIN_BP,
            "kernel_peak_bin_start_bp": peak_bin,
            "kernel_umis": kernel_total,
            "cross_fitted": true,
            "heldout_predictive_gain_over_uniform_bits_per_umi": samples.iter().enumerate().map(|(index, sample)| json!({
                "sample": sample.row.sample,
                "heldout_unique_candidate_umis": per_sample_kernels[index].unique_candidate_umis,
                "bits_per_umi": heldout_bits[index],
            })).collect::<Vec<_>>(),
            "no_candidate_umis": kernel.no_candidate_umis,
            "unique_candidate_umis": kernel.unique_candidate_umis,
            "multiple_candidate_umis": kernel.multiple_candidate_umis,
        },
        "mixture": {
            "candidate_sites_with_compatible_evidence": audit.candidate_sites,
            "recurrent_sites": audit.recurrent_sites,
            "assigned_umis": audit.assigned_umis,
        },
        "inference": {
            "tested_genes": tests.len(),
            "reported_genes": tests.iter().filter(|test| test.reported).count(),
            "test": "two-sided paired donor Student t test",
            "calibration": "exact paired sign-flip p-value",
            "multiplicity": "Benjamini-Hochberg over tested genes",
            "donor_median_effects": samples.iter().enumerate().map(|(index, sample)| json!({
                "sample": sample.row.sample, "condition": sample.row.condition,
                "median_effect": donor_medians[index],
            })).collect::<Vec<_>>(),
            "positive_donor_medians": positive_donors,
            "across_donor_median": across_donor_median,
            "registered_global_directional_threshold_pass": global_pass,
            "top_reported": tests.iter().filter(|test| test.reported).take(100).map(|test| json!({
                "gene_id": annotation.gene_ids[test.gene as usize],
                "gene_name": annotation.gene_names[test.gene as usize],
                "effect": test.effect, "q": test.q, "p_signflip": test.p_flip,
                "eligible_samples": test.eligible_samples.len(), "concordant": test.concordant,
            })).collect::<Vec<_>>(),
        },
        "performance": {"wall_seconds": started.elapsed().as_secs_f64()},
        "outputs": {"sites": "sites.tsv", "genes": "genes.tsv", "fragment_kernel": "fragment-kernel.tsv"},
    });
    std::fs::write(
        args.out_dir.join("summary.json"),
        serde_json::to_string_pretty(&summary)? + "\n",
    )?;
    if let Some(format) = args.report_format {
        let design_identity = design_identity
            .as_deref()
            .context("uniform report did not capture a design identity")?;
        let annotation_identity = annotation_identity
            .as_deref()
            .context("uniform report did not capture an annotation identity")?;
        let genome_identity = genome_identity
            .as_deref()
            .context("uniform report did not capture a genome identity")?;
        let polyasite_identity = polyasite_identity
            .as_deref()
            .context("uniform report did not capture a PolyASite catalogue identity")?;
        let context = uniform_context(
            &args,
            &samples,
            design_identity,
            annotation_identity,
            genome_identity,
            polyasite_identity,
            false,
        )?;
        let reported_out_dir = reported_output_path(&args.out_dir)?;
        let artifacts = artifact_rows(
            &args.out_dir,
            &reported_out_dir,
            &[
                ("polyasite-mixture-sites", "sites.tsv", sites.len() as u64),
                ("polyasite-mixture-genes", "genes.tsv", tests.len() as u64),
                (
                    "fragment-kernel",
                    "fragment-kernel.tsv",
                    kernel.bins.len() as u64,
                ),
                ("completion-summary", "summary.json", 1),
            ],
        )?;
        let report_summary = PolyasiteMixtureReportSummary {
            artifact: ArtifactPublicationSummary {
                output_directory: &reported_out_dir,
                directory_transactional: false,
                completion_marker: "summary.json",
                report_file_atomic_no_clobber: args.report_output.is_some(),
            },
            samples: samples.len() as u64,
            external_candidate_sites: external_candidate_count as u64,
            internal_priming_candidates_dropped: external_ip_dropped as u64,
            compatible_candidate_sites: audit.candidate_sites as u64,
            recurrent_sites: audit.recurrent_sites as u64,
            assigned_umis: audit.assigned_umis,
            tested_genes: tests.len() as u64,
            reported_genes: tests.iter().filter(|test| test.reported).count() as u64,
            fragment_kernel_umis: kernel_total,
        };
        publish_uniform_report(format, args.report_output.as_deref(), |writer| {
            write_polyasite_mixture_report(
                writer,
                format.into(),
                &context,
                &report_summary,
                &annotation,
                &contrast,
                &gene_index.rev,
                &gene_index.chrom,
                &samples,
                &sites,
                &tests,
                &kernel,
                &per_sample_kernels,
                &heldout_bits,
                &artifacts,
            )
        })?;
    } else {
        println!("{}", serde_json::to_string_pretty(&summary)?);
    }
    eprintln!(
        "cohort polyasite-mixture: {} recurrent sites, {} tested genes in {:.2}s",
        sites.len(),
        tests.len(),
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contrast_and_design_ids_are_strict() {
        assert_eq!(
            parse_contrast("astro_nsc:mature_neuron").unwrap()[1],
            "mature_neuron"
        );
        assert!(parse_contrast("x:x").is_err());
        assert!(parse_contrast("x:y:z").is_err());
    }

    #[test]
    fn paired_test_orders_null_and_shift() {
        let null = apastats::paired_test(&[0.1, -0.1, 0.05, -0.05, 0.02, -0.02]).unwrap();
        let shifted = apastats::paired_test(&[0.20, 0.15, 0.24, 0.18, 0.30, 0.12]).unwrap();
        assert!(shifted.1 < null.1);
        assert!(shifted.2 < null.2);
    }

    #[test]
    fn reverse_pas_motif_is_oriented() {
        let mut sequence = vec![b'C'; 200];
        sequence[120..126].copy_from_slice(b"TTTATT"); // reverse complement of AATAAA
        assert_eq!(pas_motif(&sequence, 100, true).as_deref(), Some("AATAAA"));
        assert!(pas_motif(&sequence, 100, false).is_none());
    }

    #[test]
    fn terminal_window_union_allows_isoforms_but_rejects_gene_overlap() {
        let windows = vec![
            GeneWindow {
                gene: 3,
                rev: false,
                lo: 100,
                hi: 220,
            },
            GeneWindow {
                gene: 3,
                rev: false,
                lo: 180,
                hi: 300,
            },
            GeneWindow {
                gene: 4,
                rev: false,
                lo: 260,
                hi: 340,
            },
        ];
        let index = ChromWindows {
            windows,
            bins: vec![vec![0, 1, 2]],
        };
        assert!(matches!(
            unique_gene(&index, 150, false),
            GeneMatch::Unique(3)
        ));
        assert!(matches!(
            unique_gene(&index, 200, false),
            GeneMatch::Unique(3)
        ));
        assert!(matches!(
            unique_gene(&index, 280, false),
            GeneMatch::Ambiguous
        ));
        assert!(matches!(unique_gene(&index, 80, false), GeneMatch::None));
    }

    #[test]
    fn catalogue_constrained_em_recovers_separated_site_mixture() {
        let mut kernel = vec![1e-9; ((KERNEL_MAX_BP - KERNEL_MIN_BP) / KERNEL_BIN_BP + 1) as usize];
        let mode = ((100 - KERNEL_MIN_BP) / KERNEL_BIN_BP) as usize;
        kernel[mode] = 1.0;
        let counts = em_site_counts(&[(900, 80), (1_400, 20)], &[1_000, 1_500], false, &kernel);
        assert!((counts.iter().sum::<f64>() - 100.0).abs() < 1e-6);
        assert!(counts[0] > 79.9);
        assert!(counts[1] > 19.9);
    }

    #[test]
    fn sparse_kernel_window_preserves_transcript_orientation_and_bounds() {
        assert_eq!(compatible_site_range(1_000, false), (950, 3_000));
        assert_eq!(compatible_site_range(1_000, true), (0, 1_050));
        assert_eq!(compatible_site_range(10, false), (0, 2_010));
        assert_eq!(
            compatible_site_range(u32::MAX - 10, true),
            (u32::MAX - 2_010, u32::MAX)
        );
    }

    #[test]
    fn endpoint_collapse_is_gene_class_deduplicated_and_strand_oriented() {
        let mut records = vec![
            EndpointRecord { gene: 0, class: 7, endpoint: 100, group: 0 },
            EndpointRecord { gene: 0, class: 7, endpoint: 140, group: 0 },
            EndpointRecord { gene: 0, class: 8, endpoint: 140, group: 1 },
            EndpointRecord { gene: 1, class: 9, endpoint: 500, group: 1 },
            EndpointRecord { gene: 1, class: 9, endpoint: 450, group: 1 },
            EndpointRecord { gene: 1, class: 10, endpoint: 450, group: 0 },
        ];
        let (counts, molecules) = collapse_endpoint_records(&mut records, &[false, true]).unwrap();
        assert_eq!(molecules, 4);
        assert_eq!(
            counts,
            vec![(140, [1, 1]), ((1u64 << 32) | 450, [1, 1])]
        );
    }

    #[test]
    fn flat_gene_ranges_preserve_sorted_coordinate_rows() {
        let sample = SampleCounts {
            row: DesignRow {
                sample: "A".into(),
                condition: "test".into(),
                archive: PathBuf::from("unused.aie"),
                archive_label: "unused.aie".into(),
                groups: PathBuf::from("unused.tsv"),
                groups_label: "unused.tsv".into(),
            },
            selected_cells: [1, 1],
            absent_design_cells: [0, 0],
            endpoint_records: 3,
            deduplicated_gene_umis: 3,
            outside_terminal_region_records: 0,
            ambiguous_terminal_region_records: 0,
            counts: vec![
                (10, [1, 0]),
                ((2u64 << 32) | 20, [0, 1]),
                ((2u64 << 32) | 30, [1, 1]),
            ],
            reference_signature: None,
            chunks: 1,
            archive_format_version: None,
            archive_identity: None,
            groups_identity: None,
        };
        assert!(gene_count_rows(&sample, 1).is_empty());
        assert_eq!(gene_count_rows(&sample, 2).len(), 2);
        assert_eq!(gene_count_rows(&sample, 2)[1].0 as u32, 30);
    }

    fn sample_with_counts(sample: &str, values: &[(u32, [u32; 2])]) -> SampleCounts {
        SampleCounts {
            row: DesignRow {
                sample: sample.to_owned(),
                condition: "test".into(),
                archive: PathBuf::from("unused.aie"),
                archive_label: "unused.aie".into(),
                groups: PathBuf::from("unused.tsv"),
                groups_label: "unused.tsv".into(),
            },
            selected_cells: [10, 10],
            absent_design_cells: [0, 0],
            endpoint_records: values
                .iter()
                .map(|(_, counts)| counts[0] as usize + counts[1] as usize)
                .sum(),
            deduplicated_gene_umis: values
                .iter()
                .map(|(_, counts)| counts[0] as usize + counts[1] as usize)
                .sum(),
            outside_terminal_region_records: 0,
            ambiguous_terminal_region_records: 0,
            counts: values
                .iter()
                .map(|(endpoint, counts)| (*endpoint as u64, *counts))
                .collect(),
            reference_signature: None,
            chunks: 1,
            archive_format_version: None,
            archive_identity: None,
            groups_identity: None,
        }
    }

    #[test]
    fn common_catalogue_preserves_explicit_sample_group_counts() {
        let samples = vec![
            sample_with_counts("A", &[(100, [12, 0]), (110, [0, 8]), (200, [5, 6])]),
            sample_with_counts("B", &[(100, [4, 6]), (205, [0, 12])]),
            sample_with_counts("C", &[(110, [5, 5]), (200, [9, 0])]),
        ];
        let sites = site_catalogue(&samples, 1, &[false], 24, 10, 3, 10).unwrap();
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].lo, 100);
        assert_eq!(sites[0].hi, 110);
        assert_eq!(sites[0].counts, vec![[12, 8], [4, 6], [5, 5]]);
    }

    #[test]
    fn fragment_kernels_leave_each_sample_out() {
        let samples = vec![
            sample_with_counts("A", &[(50, [12, 0])]),
            sample_with_counts("B", &[(80, [0, 8])]),
        ];
        let (pooled, per_sample, kernels, gains) =
            cross_fitted_kernels(&samples, &[vec![100]], &[false]);
        assert_eq!(pooled.unique_candidate_umis, 20);
        assert_eq!(per_sample[0].unique_candidate_umis, 12);
        assert_eq!(per_sample[1].unique_candidate_umis, 8);
        assert_eq!(kernels.len(), 2);
        assert_eq!(gains.len(), 2);
        assert!(
            kernels
                .iter()
                .all(|kernel| (kernel.iter().sum::<f64>() - 1.0).abs() < 1e-12)
        );
    }

    #[test]
    fn paired_gene_fit_uses_donors_and_preserves_effect_direction() {
        let mut proximal = Site {
            gene: 0,
            lo: 100,
            hi: 100,
            representative: 100,
            recurrent_samples: 8,
            counts: vec![[80, 20]; 8],
            ip: false,
            a20: 0,
            arun: 0,
            motif: Some("AATAAA".into()),
            polyasite_distance: Some(0),
            high_confidence: true,
        };
        let mut distal = Site {
            gene: 0,
            lo: 200,
            hi: 200,
            representative: 200,
            recurrent_samples: 8,
            counts: vec![[20, 80]; 8],
            ip: false,
            a20: 0,
            arun: 0,
            motif: Some("AATAAA".into()),
            polyasite_distance: Some(0),
            high_confidence: true,
        };
        // Break zero-variance symmetry while retaining a strong positive paired effect.
        for sample in 0..8 {
            proximal.counts[sample][1] += sample as u64;
            distal.counts[sample][0] += sample as u64;
        }
        let tests = fit_genes(&[proximal, distal], 1, &[false], 20, 6, 20, registered_site);
        assert_eq!(tests.len(), 1);
        assert!(tests[0].effect > 0.5);
        assert_eq!(tests[0].eligible_samples.len(), 8);
        assert_eq!(tests[0].concordant, 8);
        assert!(tests[0].lodo_stable);
    }
}
