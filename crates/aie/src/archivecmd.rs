//! `aie ingest-archive` — BAM → `.aie`, and `aie replay-rows` — quantify any GTF from either a
//! BAM (via row extraction) or an `.aie` archive. Both paths share `rows::{extract_rows,
//! replay_rows}` and must produce byte-identical matrices from equivalent evidence.

pub(crate) mod annotationcompare;
pub(crate) mod transcriptec;

use crate::rows::{
    extract_rows, extract_rows_for_archive, extract_rows_with_identity, identity_of_consumed_file,
    replay_rows_stranded, ConsumedFileIdentity, Extracted, ExtractedTerminalTails, MolChain,
    MolRec, PatAlt, ReplayRowsAccumulator, SAME_SHAPE,
};
use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use evidence_io::alignment_provenance::{
    AlignmentDeclaration, AlignmentInputs, AlignmentProvenanceManifest, DeclaredAlignmentFile,
    GenomeBindingAction, GenomeReferenceBinding, IngestProvenance, JunctionCatalogue,
    JunctionCatalogueRole, JunctionDiscoveryMode, ProvenanceStatus, VerifiedFileIdentity,
    ALIGNMENT_PROVENANCE_SCHEMA, ALIGNMENT_PROVENANCE_SECTION, JUNCTION_CATALOGUE_SECTION,
    MOLECULAR_EVIDENCE_SCHEMA,
};
use evidence_io::archive::{put_svarint, put_varint, Shape};
use evidence_io::format::{Cursor, SectionReader, SectionWriter};
use evidence_io::terminal_tail::{
    self, EncodedTerminalTailEvent, EncodedTerminalTailMolecule, TerminalTailMetadata,
    TerminalTailRoute, TerminalTailSignal,
};
use evidence_io::umi;
use gravlax_output::{
    canonical_destination_key, install_open_file_no_clobber, publish_file_no_clobber,
    reported_output_path, DataType, Durability, Field, MexManifest, OutputError, OutputFormat,
    Provenance, ResultContext, ResultEnvelope, RowSemantics, StreamingBundleWriter, TableSchema,
    TableSemantics,
};
use hashbrown::HashMap as HbMap;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use smallvec::SmallVec;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum UniformArchiveFormat {
    Text,
    Tsv,
    Json,
}

impl From<UniformArchiveFormat> for OutputFormat {
    fn from(value: UniformArchiveFormat) -> Self {
        match value {
            UniformArchiveFormat::Text => Self::Text,
            UniformArchiveFormat::Tsv => Self::Tsv,
            UniformArchiveFormat::Json => Self::Json,
        }
    }
}

/// Opt-in operation report for commands whose primary result is an archive or MEX artifact.
/// Omitting these flags preserves the established stdout and artifact behavior byte-for-byte.
#[derive(Parser, Clone, Debug, Default)]
pub struct UniformArchiveReportArgs {
    /// Emit a versioned uniform operation report instead of the legacy stdout summary.
    #[arg(long, value_enum)]
    pub report_format: Option<UniformArchiveFormat>,
    /// Atomically publish the operation report here without replacing an existing file.
    #[arg(long, requires = "report_format")]
    pub report_output: Option<PathBuf>,
}

#[derive(Parser)]
pub struct IngestArgs {
    /// Tagged, coordinate-sorted ingest BAM (CR/UR/CY, secondaries included).
    pub bam: PathBuf,
    #[arg(long)]
    pub whitelist: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long, default_value_t = 2_000)]
    pub locus_gap: u32,
    #[arg(long, default_value_t = 19)]
    pub zstd_level: i32,
    /// Genomic chunk size in megabases. Smaller chunks trade archive overhead for more selective
    /// regional reads; the ingest report prints the resulting archive size.
    #[arg(long, default_value_t = 4)]
    pub chunk_mb: u32,
    /// Reference FASTA (plain or gzipped) to bind for sequence-consulting queries. Gravlax hashes
    /// its exact bytes and normalized contigs; its relationship to the alignment is explicitly
    /// caller-declared rather than inferred. Hashing runs concurrently with extraction.
    #[arg(long)]
    pub genome: Option<PathBuf>,
    /// Retain the frozen, sequence-free terminal-tail observable from records with an explicit
    /// integer `NH=1` tag in forward-stranded 10x 3′ cDNA: trailing A on forward alignments and
    /// leading T on reverse alignments. This is not a chemistry-generic tail detector. Inspected
    /// clipped sequence is never stored.
    #[arg(long)]
    pub terminal_tails: bool,
    /// Declare how splice junctions were supplied or learned for this alignment. Gravlax records
    /// this declaration verbatim and does not infer it from an aligner command line.
    #[arg(long, value_enum, default_value_t = JunctionDiscoveryArg::Unspecified)]
    pub junction_discovery: JunctionDiscoveryArg,
    /// Exact STAR-style junction table supplied to pass 2. For per-library two-pass, this is the
    /// pass-1 output (for STAR Basic, `_STARpass1/SJ.out.tab`), not pass 2's final `SJ.out.tab`.
    /// Required for per-library-two-pass and frozen-catalogue declarations; the file's role is
    /// caller-declared, while its exact bytes and parsed data-row count are archived.
    #[arg(long)]
    pub junction_catalogue: Option<PathBuf>,
    /// Annotation supplied while building the alignment index or injecting splice junctions.
    /// Its exact-byte identity and locator are recorded; omission means no such file was declared.
    #[arg(long)]
    pub alignment_annotation: Option<PathBuf>,
    /// Caller-declared content identity or reproducible locator for the aligner index.
    #[arg(long)]
    pub alignment_index_identity: Option<String>,
    /// Ordered source-read or other files supplied to the aligner. Repeat in original argument
    /// order. Gravlax verifies each file's current bytes; their aligner role is caller-declared.
    #[arg(long = "alignment-input")]
    pub alignment_inputs: Vec<PathBuf>,
    /// Aligner log containing resolved defaults, if available.
    #[arg(long)]
    pub alignment_log: Option<PathBuf>,
    /// Caller-declared library chemistry used for alignment and strand interpretation.
    #[arg(long)]
    pub alignment_chemistry: Option<String>,
    #[command(flatten)]
    pub report: UniformArchiveReportArgs,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum JunctionDiscoveryArg {
    Unspecified,
    OnePass,
    PerLibraryTwoPass,
    FrozenCatalogue,
}

impl From<JunctionDiscoveryArg> for JunctionDiscoveryMode {
    fn from(value: JunctionDiscoveryArg) -> Self {
        match value {
            JunctionDiscoveryArg::Unspecified => Self::Unspecified,
            JunctionDiscoveryArg::OnePass => Self::OnePass,
            JunctionDiscoveryArg::PerLibraryTwoPass => Self::PerLibraryTwoPass,
            JunctionDiscoveryArg::FrozenCatalogue => Self::FrozenCatalogue,
        }
    }
}

fn junction_discovery_name(value: JunctionDiscoveryArg) -> &'static str {
    match value {
        JunctionDiscoveryArg::Unspecified => "unspecified",
        JunctionDiscoveryArg::OnePass => "one-pass",
        JunctionDiscoveryArg::PerLibraryTwoPass => "per-library-two-pass",
        JunctionDiscoveryArg::FrozenCatalogue => "frozen-catalogue",
    }
}

#[derive(Parser)]
pub struct ReplayRowsArgs {
    /// `.aie` archive, raw ingest BAM with --from-bam, or post-correction BAM with
    /// --from-molecule-bam.
    pub input: PathBuf,
    #[arg(long)]
    pub gtf: PathBuf,
    /// Interpret `input` as a BAM and extract rows on the fly (the regression reference path).
    #[arg(long, conflicts_with = "from_molecule_bam")]
    pub from_bam: bool,
    /// Interpret input as `export-molecule-bam` output. This path consumes opaque UMI-class ids
    /// and explicit 1MM edges; no whitelist or unavailable nucleotide UMI value is invented.
    #[arg(long, conflicts_with = "from_bam")]
    pub from_molecule_bam: bool,
    /// Whitelist; required with --from-bam.
    #[arg(long)]
    pub whitelist: Option<PathBuf>,
    /// Barcode list defining output column order (for example, STARsolo raw barcodes.tsv).
    #[arg(long)]
    pub barcodes: PathBuf,
    #[arg(long)]
    pub out_dir: PathBuf,
    #[arg(long, default_value_t = 2_000)]
    pub locus_gap: u32,
    /// Emit STARsolo Velocyto semantics (spliced/unspliced/ambiguous matrices) instead of Gene.
    #[arg(long)]
    pub velocity: bool,
    /// STARsolo cDNA-alignment/transcript strand relationship. 10x 3' uses forward; R2-only 10x
    /// 5' uses reverse.
    #[arg(long, value_enum, default_value_t = SoloStrandArg::Forward)]
    pub solo_strand: SoloStrandArg,
    /// Print the multi-gene ambiguity audit (the EM upside bound) instead of emitting a matrix.
    #[arg(long)]
    pub audit_multigene: bool,
    /// Diagnostic reference path: decode the entire archive before Gene replay. Normal archive
    /// replay streams bounded chunk batches and produces byte-identical output at lower peak RSS.
    #[arg(long)]
    pub eager: bool,
    #[command(flatten)]
    pub report: UniformArchiveReportArgs,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum SoloStrandArg {
    Forward,
    Reverse,
    Unstranded,
}

impl From<SoloStrandArg> for anno::assign::SoloStrand {
    fn from(value: SoloStrandArg) -> Self {
        match value {
            SoloStrandArg::Forward => Self::Forward,
            SoloStrandArg::Reverse => Self::Reverse,
            SoloStrandArg::Unstranded => Self::Unstranded,
        }
    }
}

/// Legacy/test archive writer. Production ingest uses [`write_archive_sections`] with a logical
/// v2 provenance manifest; this helper intentionally preserves the prior core-only byte contract
/// for compatibility fixtures.
#[allow(dead_code)]
pub fn write_archive(
    x: &Extracted,
    out: &Path,
    level: i32,
    chunk_bp: u32,
    genome_sig: Option<&evidence_io::genome::GenomeSig>,
) -> Result<Vec<(String, u64, u64)>> {
    write_archive_sections(
        x,
        SectionWriter::create(out, level)?,
        level,
        chunk_bp,
        genome_sig,
        ArchiveExtensions::default(),
    )?
    .finish()
}

struct PreparedTerminalTails {
    metadata: TerminalTailMetadata,
    routes: Vec<TerminalTailRoute>,
    sections: Vec<(String, Vec<u8>)>,
}

fn prepare_terminal_tails(
    x: &Extracted,
    tails: &ExtractedTerminalTails,
    chunk_bp: u32,
) -> Result<PreparedTerminalTails> {
    let mut sections = Vec::new();
    let mut routes = Vec::new();
    let mut tail_cursor = 0usize;
    let mut molecule_start = 0usize;
    let mut chunk = 0u32;
    while molecule_start < x.mols.len() {
        let chrom = x.mols[molecule_start].chrom;
        let bin_start = (x.mols[molecule_start].anchor() / chunk_bp) * chunk_bp;
        let bin_end = bin_start
            .checked_add(chunk_bp)
            .context("terminal-tail molecule chunk boundary overflow")?;
        let mut molecule_end = molecule_start;
        while molecule_end < x.mols.len()
            && x.mols[molecule_end].chrom == chrom
            && x.mols[molecule_end].anchor() < bin_end
        {
            molecule_end += 1;
        }
        if tail_cursor < tails.molecules.len()
            && tails.molecules[tail_cursor].molecule_ordinal < molecule_start as u64
        {
            bail!("terminal-tail attachments are not in molecule order");
        }
        let selected_start = tail_cursor;
        while tail_cursor < tails.molecules.len()
            && tails.molecules[tail_cursor].molecule_ordinal < molecule_end as u64
        {
            tail_cursor += 1;
        }
        if selected_start != tail_cursor {
            let selected = &tails.molecules[selected_start..tail_cursor];
            let mut encoded = Vec::with_capacity(selected.len());
            let mut min_anchor = u32::MAX;
            let mut max_anchor = 0u32;
            let mut event_count = 0usize;
            for molecule_events in selected {
                let ordinal = usize::try_from(molecule_events.molecule_ordinal)
                    .context("terminal-tail molecule ordinal exceeds usize")?;
                let molecule = x
                    .mols
                    .get(ordinal)
                    .context("terminal-tail attachment references an absent molecule")?;
                if molecule.chrom != chrom {
                    bail!("terminal-tail attachment chromosome differs from its molecule");
                }
                let mut events = Vec::with_capacity(molecule_events.events.len());
                for &(anchor, signal) in &molecule_events.events {
                    let anchor_delta = i64::from(anchor) - i64::from(molecule.anchor());
                    events.push(EncodedTerminalTailEvent {
                        anchor_delta,
                        reverse: molecule.strand_rev,
                        signal,
                    });
                    min_anchor = min_anchor.min(anchor);
                    max_anchor = max_anchor.max(anchor);
                    event_count += 1;
                }
                encoded.push(EncodedTerminalTailMolecule {
                    local_ordinal: u32::try_from(ordinal - molecule_start)
                        .context("terminal-tail local molecule ordinal exceeds u32")?,
                    events,
                });
            }
            let molecule_count = u32::try_from(molecule_end - molecule_start)
                .context("terminal-tail ordinary chunk size exceeds u32")?;
            sections.push((
                format!("tail.c{chunk}"),
                terminal_tail::encode_chunk(molecule_count, &encoded)?,
            ));
            routes.push(TerminalTailRoute {
                chunk,
                chrom,
                min_anchor,
                max_anchor,
                selected_molecules: u32::try_from(selected.len())
                    .context("terminal-tail selected count exceeds u32")?,
                events: u32::try_from(event_count)
                    .context("terminal-tail event count exceeds u32")?,
            });
        }
        molecule_start = molecule_end;
        chunk = chunk
            .checked_add(1)
            .context("terminal-tail chunk id overflow")?;
    }
    if tail_cursor != tails.molecules.len() {
        bail!("terminal-tail attachment references an absent molecule");
    }
    let selected_molecules = tails.molecules.len() as u64;
    let events = tails.molecules.iter().try_fold(0u64, |sum, molecule| {
        sum.checked_add(molecule.events.len() as u64)
            .context("terminal-tail event total overflow")
    })?;
    let metadata = TerminalTailMetadata::new(
        selected_molecules,
        events,
        u32::try_from(routes.len()).context("terminal-tail route count exceeds u32")?,
    );
    metadata.validate()?;
    Ok(PreparedTerminalTails {
        metadata,
        routes,
        sections,
    })
}

#[derive(Clone, Copy, Default)]
struct ArchiveExtensions<'a> {
    alignment_provenance: Option<&'a AlignmentProvenanceManifest>,
    junction_catalogue_bytes: Option<&'a [u8]>,
    terminal_tails: Option<&'a ExtractedTerminalTails>,
    genome_reference_binding: Option<&'a GenomeReferenceBinding>,
}

fn write_archive_sections(
    x: &Extracted,
    mut w: SectionWriter,
    level: i32,
    chunk_bp: u32,
    genome_sig: Option<&evidence_io::genome::GenomeSig>,
    extensions: ArchiveExtensions<'_>,
) -> Result<SectionWriter> {
    // Sections are built first and compressed in one rayon pass at the end — zstd-19 over the
    // chunk payloads is the ingest wall-clock, and it parallelizes embarrassingly.
    let mut sections: Vec<(String, Vec<u8>)> = Vec::new();

    if chunk_bp == 0 {
        bail!("archive chunk size must be positive");
    }
    if extensions.terminal_tails.is_some() && extensions.alignment_provenance.is_none() {
        bail!("terminal-tail evidence requires a logical v2 alignment-provenance manifest");
    }
    if let Some(manifest) = extensions.alignment_provenance {
        let rule_declared = manifest.ingest.terminal_tail_rule.as_deref()
            == Some(terminal_tail::TERMINAL_TAIL_RULE);
        if rule_declared != extensions.terminal_tails.is_some() {
            bail!("terminal-tail capability and provenance extraction rule disagree");
        }
    }
    match (
        extensions
            .alignment_provenance
            .and_then(|manifest| manifest.alignment.junction_catalogue.as_ref()),
        extensions.junction_catalogue_bytes,
    ) {
        (Some(catalogue), Some(bytes)) => {
            if catalogue.identity.bytes != bytes.len() as u64
                || catalogue.identity.blake3 != blake3::hash(bytes).to_hex().as_str()
            {
                bail!("junction catalogue bytes disagree with their provenance identity");
            }
        }
        (None, None) => {}
        _ => bail!("junction catalogue manifest and exact-byte section must occur together"),
    }
    if extensions.alignment_provenance.is_some() {
        match (genome_sig, extensions.genome_reference_binding) {
            (Some(signature), Some(binding)) if &binding.signature == signature => {
                binding.validate()?;
            }
            (None, None) => {}
            _ => bail!("logical-v2 genome signature and reference binding must occur together"),
        }
    }
    let prepared_tails = extensions
        .terminal_tails
        .map(|tails| prepare_terminal_tails(x, tails, chunk_bp))
        .transpose()?;

    let mut meta = serde_json::Map::new();
    meta.insert("mols".into(), serde_json::json!(x.mols.len()));
    meta.insert("edges".into(), serde_json::json!(x.edges.len()));
    meta.insert("cells".into(), serde_json::json!(x.cells.len()));
    meta.insert("shapes".into(), serde_json::json!(x.shapes.len()));
    meta.insert("patterns".into(), serde_json::json!(x.patterns.len()));
    meta.insert("classes".into(), serde_json::json!(x.n_classes));
    meta.insert("chunk_bp".into(), serde_json::json!(chunk_bp));
    meta.insert("chunk_streams".into(), serde_json::json!(10u32));
    meta.insert("coc_block".into(), serde_json::json!(COC_BLOCK));
    meta.insert("codec".into(), serde_json::json!("rans2"));
    if let Some(sig) = genome_sig {
        meta.insert("genome_sig".into(), serde_json::to_value(sig)?);
    }
    if let Some(binding) = extensions.genome_reference_binding {
        meta.insert(
            "genome_reference_binding".into(),
            serde_json::to_value(binding)?,
        );
    }
    if extensions.alignment_provenance.is_some() {
        meta.insert(
            "evidence_schema".into(),
            serde_json::json!(MOLECULAR_EVIDENCE_SCHEMA),
        );
        meta.insert(
            "alignment_provenance".into(),
            serde_json::json!(ALIGNMENT_PROVENANCE_SCHEMA),
        );
    }
    if let Some(prepared) = &prepared_tails {
        meta.insert(
            "terminal_tail".into(),
            serde_json::to_value(&prepared.metadata)?,
        );
    }
    sections.push(("meta".into(), serde_json::to_string(&meta)?.into_bytes()));
    if let Some(manifest) = extensions.alignment_provenance {
        sections.push((
            ALIGNMENT_PROVENANCE_SECTION.into(),
            manifest.to_canonical_json()?,
        ));
    }
    if let Some(bytes) = extensions.junction_catalogue_bytes {
        sections.push((JUNCTION_CATALOGUE_SECTION.into(), bytes.to_vec()));
    }
    sections.push(("chroms".into(), x.chrom_names.join("\n").into_bytes()));

    let mut cells = Vec::with_capacity(x.cells.len() * 4);
    for c in &x.cells {
        cells.extend_from_slice(&c.to_le_bytes());
    }
    sections.push(("cells".into(), cells));

    let mut shp = Vec::new();
    for sdef in &x.shapes {
        put_varint(&mut shp, sdef.blocks.len() as u64);
        for (off, len) in &sdef.blocks {
            put_varint(&mut shp, *off as u64);
            put_varint(&mut shp, *len as u64);
        }
    }
    sections.push(("shapes".into(), shp));

    // Patterns with the measured same-shape flag: 78.7% of alternatives share the applying row's
    // anchor shape, so those store one flag bit instead of a shape id.
    let mut pat = Vec::new();
    for pdef in &x.patterns {
        put_varint(&mut pat, pdef.len() as u64);
        for a in pdef {
            let same = a.shape == SAME_SHAPE;
            put_varint(&mut pat, ((a.chrom as u64) << 2) | ((same as u64) << 1) | (a.strand_flip as u64));
            put_svarint(&mut pat, a.offset);
            if !same {
                put_varint(&mut pat, a.shape as u64);
            }
        }
    }
    sections.push(("patterns".into(), pat));

    // ---- Chunks: contiguous runs of molecules with one (chrom, bin) each; per-chunk streams are
    // length-prefixed and concatenated into a single zstd frame (the chunk is the access unit). ----
    struct ChunkMeta {
        chrom: u32,
        bin_start: u32,
        n_mols: u32,
        class_base: u32,
        max_anchor: u32,
        n_cells: u32,
    }
    let mut chunk_meta: Vec<ChunkMeta> = Vec::new();
    let mut junctions: FxHashMap<(u32, u32, u32), u32> = FxHashMap::default();
    let mut j_postings: Vec<Vec<u32>> = Vec::new();
    let mut cell_postings: FxHashMap<u32, Vec<u32>> = FxHashMap::default();

    let mut next_class = 0u32;
    let mut i = 0usize;
    let mut chunk_idx = 0u32;
    // Cell is a pure function of the UMI class because classes are global (cell, value) sets. A
    // 22.16-million-molecule benchmark archive had zero violations of this invariant, so cell is
    // stored once per class in fresh-introduction order. COC_BLOCK-sized sections ("coc.<b>",
    // first value absolute, rest deltas) let queries decompress only the blocks containing classes
    // they touch, while full loads decode blocks in parallel.
    // Absolute values collected; each block is encoded BOTH as delta-varints and as rANS over
    // absolute ids at the end, and the smaller wins (a one-byte codec tag per block). Benchmarks
    // showed delta coding runs 15% above the independent-symbol bound when cell locality vanishes.
    let mut coc_vals: Vec<u32> = Vec::new();
    struct PendingChunk {
        idx: u32,
        anchor_s: Vec<u8>,
        layout_s: Vec<u8>,
        rep_shape_s: Vec<u8>,
        mm_shape_s: Vec<u8>,
        mm_pat_s: Vec<u8>,
        class_v: Vec<u64>,
        weight_v: Vec<u64>,
        rep_pos_v: Vec<u64>,
        mm_pos_v: Vec<u64>,
        mm_w_v: Vec<u64>,
    }
    let mut pending: Vec<PendingChunk> = Vec::new();
    let mut rans_counts = [[0u64; evidence_io::rans::NSYM]; 5];
    while i < x.mols.len() {
        let chrom = x.mols[i].chrom;
        let bin_start = (x.mols[i].anchor() / chunk_bp) * chunk_bp;
        let bin_end = bin_start + chunk_bp;
        let mut j = i;
        while j < x.mols.len() && x.mols[j].chrom == chrom && x.mols[j].anchor() < bin_end {
            j += 1;
        }
        let mols = &x.mols[i..j];

        let (mut anchor_s, mut layout_s) = (Vec::new(), Vec::new());
        let (mut rep_shape_s, mut mm_shape_s, mut mm_pat_s) = (Vec::new(), Vec::new(), Vec::new());
        // Value-collected streams use rANS with global static tables because zstd sits above the
        // order-0 value bound on these streams. Context-structured streams remain varint+zstd,
        // where zstd beats memoryless coding.
        let (mut class_v, mut weight_v, mut rep_pos_v, mut mm_pos_v, mut mm_w_v) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let class_base = next_class;
        let mut lanchor = bin_start;
        for m in mols {
            let anchor = m.anchor();
            put_varint(&mut anchor_s, (anchor - lanchor) as u64);
            lanchor = anchor;
            if m.umi_class == next_class {
                class_v.push(0u64);
                coc_vals.push(m.cell);
                next_class += 1;
            } else {
                class_v.push((next_class - m.umi_class) as u64);
            }
            // Chains are position-sorted at extraction, so for any molecule with chains the first
            // stored rep IS the anchor — its offset is elided with no flag (implied by
            // n_chains > 0). The invariant is enforced here, not assumed.
            let elide_first_rep = !m.chains.is_empty();
            if elide_first_rep && m.chains[0].reps[0].0 != anchor {
                anyhow::bail!(
                    "chain order invariant broken: first rep {} != anchor {anchor}",
                    m.chains[0].reps[0].0
                );
            }
            layout_s.push(m.strand_rev as u8);
            put_varint(&mut layout_s, m.chains.len() as u64);
            put_varint(&mut layout_s, m.mms.len() as u64);
            cell_postings.entry(m.cell).or_default().push(chunk_idx);
            let mut note_junctions = |pos: u32, shape: u32| {
                let sdef = &x.shapes[shape as usize];
                for wpair in sdef.blocks.windows(2) {
                    let key = (chrom, pos + wpair[0].0 + wpair[0].1, pos + wpair[1].0);
                    let next = junctions.len() as u32;
                    let jid = *junctions.entry(key).or_insert_with(|| {
                        j_postings.push(Vec::new());
                        next
                    });
                    j_postings[jid as usize].push(chunk_idx);
                }
            };
            let mut first_rep = true;
            for ch in &m.chains {
                weight_v.push(((ch.weight as u64) << 1) | ((ch.reps.len() == 2) as u64));
                for (pos, shape) in &ch.reps {
                    if !(first_rep && elide_first_rep) {
                        rep_pos_v.push((*pos - anchor) as u64);
                    }
                    first_rep = false;
                    put_varint(&mut rep_shape_s, *shape as u64);
                    note_junctions(*pos, *shape);
                }
            }
            for (pos, shape, pattern, weight) in &m.mms {
                mm_pos_v.push((*pos - anchor) as u64);
                put_varint(&mut mm_shape_s, *shape as u64);
                put_varint(&mut mm_pat_s, *pattern as u64);
                mm_w_v.push(*weight as u64);
                note_junctions(*pos, *shape);
            }
        }

        evidence_io::rans::count(&class_v, &mut rans_counts[0]);
        evidence_io::rans::count(&weight_v, &mut rans_counts[1]);
        evidence_io::rans::count(&rep_pos_v, &mut rans_counts[2]);
        evidence_io::rans::count(&mm_pos_v, &mut rans_counts[3]);
        evidence_io::rans::count(&mm_w_v, &mut rans_counts[4]);
        pending.push(PendingChunk {
            idx: chunk_idx,
            anchor_s, layout_s, rep_shape_s, mm_shape_s, mm_pat_s,
            class_v, weight_v, rep_pos_v, mm_pos_v, mm_w_v,
        });
        let max_anchor = mols.last().map(|m| m.anchor()).unwrap_or(bin_start);
        let n_cells_chunk = {
            let set: FxHashSet<u32> = mols.iter().map(|m| m.cell).collect();
            set.len() as u32
        };
        chunk_meta.push(ChunkMeta {
            chrom, bin_start, n_mols: mols.len() as u32, class_base,
            max_anchor, n_cells: n_cells_chunk,
        });
        chunk_idx += 1;
        i = j;
    }

    // Global static tables from whole-dataset counts, then parallel chunk assembly. Table 6
    // models the cell-of-class values (absolute ids) for the per-block codec choice.
    let mut coc_counts = [0u64; evidence_io::rans::NSYM];
    let coc_u64: Vec<u64> = coc_vals.iter().map(|v| *v as u64).collect();
    evidence_io::rans::count(&coc_u64, &mut coc_counts);
    let tables: Vec<evidence_io::rans::Table> = rans_counts
        .iter()
        .chain(std::iter::once(&coc_counts))
        .map(evidence_io::rans::Table::from_counts)
        .collect::<Result<_>>()?;
    let mut tbl_raw = Vec::new();
    for t in &tables {
        t.serialize(&mut tbl_raw);
    }
    sections.push(("rans.tables".into(), tbl_raw));
    let coc_blocks: Vec<(String, Vec<u8>)> = coc_u64
        .par_chunks(COC_BLOCK as usize)
        .enumerate()
        .map(|(b, vals)| {
            // Candidate A: absolute first, zigzag deltas after (the v1.3 coding).
            let mut a = vec![0u8];
            put_varint(&mut a, vals[0]);
            let mut last = vals[0] as i64;
            for v in &vals[1..] {
                put_svarint(&mut a, *v as i64 - last);
                last = *v as i64;
            }
            // Candidate B: rANS over absolute ids with the global coc table.
            let mut bpay = vec![1u8];
            evidence_io::rans::encode(vals, &tables[5], &mut bpay);
            (format!("coc.{b}"), if bpay.len() < a.len() { bpay } else { a })
        })
        .collect();
    let assembled: Vec<(String, Vec<u8>)> = pending
        .par_iter()
        .map(|p| {
            let mut chunk = Vec::new();
            let rans_seg = |vals: &[u64], t: &evidence_io::rans::Table, chunk: &mut Vec<u8>| {
                let mut seg = Vec::new();
                evidence_io::rans::encode(vals, t, &mut seg);
                put_varint(chunk, seg.len() as u64);
                chunk.extend_from_slice(&seg);
            };
            let byte_seg = |bytes: &[u8], chunk: &mut Vec<u8>| {
                put_varint(chunk, bytes.len() as u64);
                chunk.extend_from_slice(bytes);
            };
            byte_seg(&p.anchor_s, &mut chunk);
            rans_seg(&p.class_v, &tables[0], &mut chunk);
            byte_seg(&p.layout_s, &mut chunk);
            rans_seg(&p.weight_v, &tables[1], &mut chunk);
            rans_seg(&p.rep_pos_v, &tables[2], &mut chunk);
            byte_seg(&p.rep_shape_s, &mut chunk);
            rans_seg(&p.mm_pos_v, &tables[3], &mut chunk);
            byte_seg(&p.mm_shape_s, &mut chunk);
            byte_seg(&p.mm_pat_s, &mut chunk);
            rans_seg(&p.mm_w_v, &tables[4], &mut chunk);
            (format!("c{}", p.idx), chunk)
        })
        .collect();
    sections.extend(assembled);
    drop(pending);

    let mut cm = Vec::new();
    for c in &chunk_meta {
        put_varint(&mut cm, c.chrom as u64);
        put_varint(&mut cm, c.bin_start as u64);
        put_varint(&mut cm, c.n_mols as u64);
        put_varint(&mut cm, c.class_base as u64);
        put_varint(&mut cm, (c.max_anchor - c.bin_start) as u64);
        put_varint(&mut cm, c.n_cells as u64);
    }
    sections.push(("index.chunks".into(), cm));
    sections.extend(coc_blocks);
    if let Some(prepared) = prepared_tails {
        sections.push((
            terminal_tail::TERMINAL_TAIL_INDEX_SECTION.into(),
            terminal_tail::encode_index(&prepared.routes)?,
        ));
        sections.extend(prepared.sections);
    }

    // Junction catalogue + postings, sorted by coordinate for range-friendly lookup.
    {
        let mut cat: Vec<((u32, u32, u32), u32)> = junctions.into_iter().collect();
        cat.sort_unstable_by_key(|(k, _)| *k);
        let mut jc = Vec::new();
        let (mut lch, mut ld) = (u32::MAX, 0u32);
        for ((ch, d, a), _) in &cat {
            if *ch != lch {
                lch = *ch;
                ld = 0;
            }
            put_varint(&mut jc, *ch as u64);
            put_varint(&mut jc, (*d - ld) as u64);
            ld = *d;
            put_varint(&mut jc, (*a - *d) as u64);
        }
        sections.push(("index.junctions".into(), jc));
        let mut jp = Vec::new();
        for (_, old_id) in &cat {
            let posts = &j_postings[*old_id as usize];
            // Total supporting children (a cheap genome-wide support count served from the index
            // alone); chunk list deduplicated for the postings walk.
            put_varint(&mut jp, posts.len() as u64);
            let mut uniq: Vec<u32> = posts.clone();
            uniq.dedup();
            put_varint(&mut jp, uniq.len() as u64);
            let mut last = 0u32;
            for c in uniq {
                put_varint(&mut jp, (c - last) as u64);
                last = c;
            }
        }
        sections.push(("index.jpost".into(), jp));
    }

    // Cell postings: for each cell id in order, the chunks it appears in.
    {
        let mut cp = Vec::new();
        for cell in 0..x.cells.len() as u32 {
            match cell_postings.get(&cell) {
                Some(posts) => {
                    let mut uniq = posts.clone();
                    uniq.dedup();
                    put_varint(&mut cp, uniq.len() as u64);
                    let mut last = 0u32;
                    for c in uniq {
                        put_varint(&mut cp, (c - last) as u64);
                        last = c;
                    }
                }
                None => put_varint(&mut cp, 0),
            }
        }
        sections.push(("index.cellpost".into(), cp));
    }

    let mut edg = Vec::new();
    let mut la = 0i64;
    for (a, b) in &x.edges {
        put_svarint(&mut edg, *a as i64 - la);
        la = *a as i64;
        put_varint(&mut edg, (*b - *a) as u64);
    }
    sections.push(("edges".into(), edg));

    // One parallel compression pass, then a serial ordered write.
    let compressed: Vec<(String, usize, Vec<u8>)> = sections
        .into_par_iter()
        .map(|(name, raw)| {
            let comp = evidence_io::format::compress(&raw, level)?;
            Ok((name, raw.len(), comp))
        })
        .collect::<Result<_>>()?;
    for (name, raw_len, comp) in &compressed {
        w.section_precompressed(name, *raw_len as u64, comp)?;
    }

    Ok(w)
}

/// Report pattern-dictionary measurements at ingest without changing the archive format.
pub fn pattern_stats(x: &Extracted) {
    use rustc_hash::FxHashSet;
    let mut alts = 0u64;
    let mut same_shape = 0u64;
    let mut hops: FxHashSet<(u32, i64)> = FxHashSet::default();
    let mut anchor_shape_hits = 0u64;
    for (pi, pdef) in x.patterns.iter().enumerate() {
        let _ = pi;
        // The anchor is the offset-0 entry; its shape is the reference for the same-shape flag.
        let anchor_shape = pdef.iter().find(|a| a.offset == 0).map(|a| a.shape);
        for a in pdef {
            alts += 1;
            hops.insert((a.chrom, a.offset));
            if Some(a.shape) == anchor_shape {
                same_shape += 1;
                if a.offset != 0 {
                    anchor_shape_hits += 1;
                }
            }
        }
    }
    let _ = anchor_shape_hits;
    // Projected dictionary sizes: current entry = (chrom+flip varint ~1-4B, offset svarint ~1-5B,
    // shape varint ~1-3B) vs factored = hop-id varint into a hop vocabulary + same-shape bit.
    println!(
        "pattern #3 stats: {} patterns, {} alt entries, {} distinct (chrom,offset) hops ({:.2}x hop sharing), {:.1}% alts share the anchor's shape",
        x.patterns.len(),
        alts,
        hops.len(),
        alts as f64 / hops.len().max(1) as f64,
        100.0 * same_shape as f64 / alts.max(1) as f64
    );
}

/// Classes per cell-of-class block; each block is its own section ("coc.<b>") so it can be
/// decompressed independently.
pub const COC_BLOCK: u32 = 1 << 16;

/// Layout guard: decoding the wrong per-chunk stream count or cell-of-class layout produces
/// garbage rather than an error, so reject every unexpected value.
fn check_layout(meta: &serde_json::Value) -> Result<()> {
    let streams = meta["chunk_streams"].as_u64().unwrap_or(11);
    if streams != 10 {
        bail!("archive chunk layout has {streams} streams; this reader expects 10 — re-ingest with the current binary");
    }
    let coc = meta["coc_block"].as_u64().unwrap_or(0);
    if coc != COC_BLOCK as u64 {
        bail!("archive coc_block is {coc}; this reader expects {COC_BLOCK} — re-ingest with the current binary");
    }
    let codec = meta["codec"].as_str().unwrap_or("varint");
    if codec != "rans2" {
        bail!("archive stream codec is {codec}; this reader expects rans2 — re-ingest with the current binary");
    }
    Ok(())
}

/// Fields repeated in `meta` and the root-bound ingest manifest describe the same physical
/// layout. Validate both copies so a re-rooted archive cannot retain decodable payloads while
/// making a false construction claim in its provenance manifest.
fn check_provenance_layout(meta: &serde_json::Value, ingest: &IngestProvenance) -> Result<()> {
    let chunk_bp = required_meta_u32(meta, "chunk_bp")?;
    if chunk_bp != ingest.chunk_bp {
        bail!(
            "alignment provenance ingest layout disagrees with archive metadata: ingest.chunk_bp={} but meta.chunk_bp={chunk_bp}",
            ingest.chunk_bp
        );
    }
    let chunk_streams = required_meta_u32(meta, "chunk_streams")?;
    if chunk_streams != ingest.molecule_chunk_streams {
        bail!(
            "alignment provenance ingest layout disagrees with archive metadata: ingest.molecule_chunk_streams={} but meta.chunk_streams={chunk_streams}",
            ingest.molecule_chunk_streams
        );
    }
    let codec = meta
        .get("codec")
        .and_then(serde_json::Value::as_str)
        .context("archive meta.codec is missing or is not a string")?;
    if codec != ingest.molecule_codec {
        bail!(
            "alignment provenance ingest layout disagrees with archive metadata: ingest.molecule_codec={} but meta.codec={codec}",
            ingest.molecule_codec
        );
    }
    Ok(())
}

/// The six rANS tables: five chunk streams (class, weight, rep.pos, mm.pos, mm.weight) plus the
/// cell-of-class value table used by rANS-coded coc blocks.
pub fn read_rans_tables(r: &mut SectionReader) -> Result<Vec<evidence_io::rans::Table>> {
    let raw = r.read("rans.tables")?;
    let mut c = Cursor::new(&raw);
    let tables: Vec<_> =
        (0..6).map(|_| evidence_io::rans::Table::deserialize(&mut c)).collect::<Result<_>>()?;
    if !c.is_empty() {
        bail!("rans.tables has {} trailing bytes", raw.len() - c.position());
    }
    Ok(tables)
}

/// Decode one coc block. First byte tags the codec the writer chose for this block:
/// 0 = absolute-then-deltas, 1 = rANS over absolute ids (global coc table).
fn decode_coc_block(comp: &[u8], raw_len: usize, coc_table: &evidence_io::rans::Table) -> Result<Vec<u32>> {
    let raw = evidence_io::format::decompress(comp, raw_len)?;
    if raw.is_empty() {
        return Ok(Vec::new());
    }
    match raw[0] {
        0 => {
            let mut c = Cursor::new(&raw[1..]);
            let mut out = Vec::with_capacity(COC_BLOCK as usize);
            if c.is_empty() {
                return Ok(out);
            }
            let mut last = i64::try_from(c.varint()?).context("cell id exceeds i64")?;
            out.push(u32::try_from(last).context("cell id exceeds u32")?);
            while !c.is_empty() {
                last = last
                    .checked_add(c.svarint()?)
                    .context("cell-id delta overflow")?;
                out.push(u32::try_from(last).context("cell id is negative or exceeds u32")?);
            }
            Ok(out)
        }
        1 => evidence_io::rans::decode_limited(&raw[1..], coc_table, COC_BLOCK as usize)?
            .into_iter()
            .map(|v| u32::try_from(v).context("cell id exceeds u32"))
            .collect(),
        b => bail!("unknown coc block codec {b}"),
    }
}

fn decode_shapes(raw: &[u8]) -> Result<Vec<Shape>> {
    let mut c = Cursor::new(raw);
    let mut shapes = Vec::new();
    while !c.is_empty() {
        let n = usize::try_from(c.varint()?).context("shape block count exceeds usize")?;
        let remaining = raw.len() - c.position();
        if n > remaining / 2 {
            bail!("shape declares {n} blocks in only {remaining} remaining bytes");
        }
        let mut blocks = Vec::with_capacity(n);
        for _ in 0..n {
            let offset = u32::try_from(c.varint()?).context("shape block offset exceeds u32")?;
            let length = u32::try_from(c.varint()?).context("shape block length exceeds u32")?;
            offset
                .checked_add(length)
                .context("shape block coordinate overflow")?;
            blocks.push((offset, length));
        }
        shapes.push(Shape { blocks });
    }
    Ok(shapes)
}

fn decode_patterns(raw: &[u8]) -> Result<Vec<Vec<PatAlt>>> {
    let mut c = Cursor::new(raw);
    let mut patterns = Vec::new();
    while !c.is_empty() {
        let n = usize::try_from(c.varint()?).context("pattern length exceeds usize")?;
        let remaining = raw.len() - c.position();
        if n > remaining / 2 {
            bail!("pattern declares {n} alternatives in only {remaining} remaining bytes");
        }
        let mut alts = Vec::with_capacity(n);
        for _ in 0..n {
            let cf = c.varint()?;
            let same = (cf >> 1) & 1 == 1;
            alts.push(PatAlt {
                chrom: u32::try_from(cf >> 2).context("pattern chromosome exceeds u32")?,
                strand_flip: cf & 1 == 1,
                offset: c.svarint()?,
                shape: if same {
                    SAME_SHAPE
                } else {
                    u32::try_from(c.varint()?).context("pattern shape exceeds u32")?
                },
            });
        }
        patterns.push(alts);
    }
    Ok(patterns)
}

fn required_meta_u32(meta: &serde_json::Value, key: &str) -> Result<u32> {
    let value = meta
        .get(key)
        .and_then(serde_json::Value::as_u64)
        .with_context(|| format!("archive meta.{key} is missing or is not an unsigned integer"))?;
    u32::try_from(value).with_context(|| format!("archive meta.{key} exceeds u32"))
}

fn required_meta_usize(meta: &serde_json::Value, key: &str) -> Result<usize> {
    let value = meta
        .get(key)
        .and_then(serde_json::Value::as_u64)
        .with_context(|| format!("archive meta.{key} is missing or is not an unsigned integer"))?;
    usize::try_from(value).with_context(|| format!("archive meta.{key} exceeds usize"))
}

/// Shared decode of the dictionary sections.
pub struct Dicts {
    pub cells: Vec<u32>,
    pub shapes: Vec<Shape>,
    pub patterns: Vec<Vec<PatAlt>>,
    pub chrom_names: Vec<String>,
    /// Cell id per UMI class, in class order. Cell is a pure function of class, so it is stored
    /// once per class instead of once per molecule in each chunk.
    pub cell_of_class: Vec<u32>,
    pub rans_tables: Vec<evidence_io::rans::Table>,
    pub n_classes: u32,
    pub n_mols: usize,
}

pub fn read_dicts(r: &mut SectionReader) -> Result<Dicts> {
    let meta: serde_json::Value = serde_json::from_slice(&r.read("meta")?)?;
    // Layout guard: a reader decoding the wrong number of per-chunk streams produces garbage, not
    // an error — learned the hard way when a rebuilt binary read an older archive mid-analysis.
    check_layout(&meta)?;
    let chrom_names: Vec<String> =
        String::from_utf8_lossy(&r.read("chroms")?).lines().map(|s| s.to_string()).collect();
    let cells_raw = r.read("cells")?;
    if cells_raw.len() % 4 != 0 {
        bail!("cells section has {} trailing byte(s)", cells_raw.len() % 4);
    }
    let cells: Vec<u32> = cells_raw
        .as_chunks::<4>()
        .0
        .iter()
        .map(|&bytes| u32::from_le_bytes(bytes))
        .collect();
    let shapes = decode_shapes(&r.read("shapes")?)?;
    let patterns = decode_patterns(&r.read("patterns")?)?;
    let rans_tables = read_rans_tables(r)?;
    let n_classes = required_meta_u32(&meta, "classes")?;
    let n_blocks = n_classes.div_ceil(COC_BLOCK) as usize;
    let blocks: Vec<Vec<u32>> = (0..n_blocks)
        .into_par_iter()
        .map(|b| {
            let (comp, raw_len) = r.read_compressed_at(&format!("coc.{b}"))?;
            decode_coc_block(&comp, raw_len, &rans_tables[5])
        })
        .collect::<Result<_>>()?;
    let cell_of_class: Vec<u32> = blocks.concat();
    if cell_of_class.len() != n_classes as usize {
        anyhow::bail!(
            "coc blocks hold {} entries for {} classes",
            cell_of_class.len(),
            n_classes
        );
    }
    Ok(Dicts {
        cells,
        shapes,
        patterns,
        chrom_names,
        cell_of_class,
        rans_tables,
        n_classes,
        n_mols: required_meta_usize(&meta, "mols")?,
    })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChunkInfo {
    pub chrom: u32,
    pub bin_start: u32,
    pub n_mols: u32,
    pub class_base: u32,
    pub max_anchor: u32,
    pub n_cells: u32,
}

/// One decoded terminal-tail event attached to a stable serialized molecule record.
/// `molecule_ordinal` is global within the archive; `(chunk, local_molecule_ordinal)` is the
/// equivalent routing key for bounded readers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TerminalTailRecord {
    pub molecule_ordinal: u64,
    pub chunk: u32,
    pub local_molecule_ordinal: u32,
    pub cell: u32,
    pub umi_class: u32,
    pub chrom: u32,
    pub strand_rev: bool,
    pub anchor: u32,
    pub signal: TerminalTailSignal,
}

pub fn read_chunk_index(r: &mut SectionReader) -> Result<Vec<ChunkInfo>> {
    let raw = r.read("index.chunks")?;
    let mut c = Cursor::new(&raw);
    let mut out = Vec::new();
    while !c.is_empty() {
        let chrom = u32::try_from(c.varint()?).context("chunk chromosome id exceeds u32")?;
        let bin_start = u32::try_from(c.varint()?).context("chunk start exceeds u32")?;
        let n_mols = u32::try_from(c.varint()?).context("chunk molecule count exceeds u32")?;
        let class_base = u32::try_from(c.varint()?).context("chunk class base exceeds u32")?;
        let max_delta = u32::try_from(c.varint()?).context("chunk max-anchor delta exceeds u32")?;
        let max_anchor = bin_start
            .checked_add(max_delta)
            .context("chunk max-anchor overflow")?;
        let n_cells = u32::try_from(c.varint()?).context("chunk cell count exceeds u32")?;
        out.push(ChunkInfo {
            chrom,
            bin_start,
            n_mols,
            class_base,
            max_anchor,
            n_cells,
        });
    }
    Ok(out)
}

/// Decode one chunk section into molecule records. `cell_of_class` is the global per-class cell
/// table from the dictionaries; backref molecules routinely reference classes introduced in other
/// chunks (same cell+UMI at another locus), so cell resolution is global by construction. Pass
/// `None` to defer cell resolution (lazy queries): `cell` is then `u32::MAX` and callers resolve
/// via `LazyArchive::cell_of` for the molecules they actually touch.
pub fn decode_chunk(
    raw: &[u8],
    info: &ChunkInfo,
    cell_of_class: Option<&[u32]>,
    tables: &[evidence_io::rans::Table],
) -> Result<Vec<MolRec>> {
    if tables.len() < 5 {
        bail!("chunk decoder needs five rANS tables, got {}", tables.len());
    }
    let mut top = Cursor::new(raw);
    let mut streams: Vec<&[u8]> = Vec::with_capacity(10);
    for _ in 0..10 {
        let len = usize::try_from(top.varint()?).context("chunk stream length is too large")?;
        streams.push(top.take(len)?);
    }
    if !top.is_empty() {
        bail!("chunk has {} trailing bytes", raw.len() - top.position());
    }
    let (mut anchor_c, mut layout_c) = (Cursor::new(streams[0]), Cursor::new(streams[2]));
    let (mut rep_shape_c, mut mm_shape_c, mut mm_pat_c) =
        (Cursor::new(streams[5]), Cursor::new(streams[7]), Cursor::new(streams[8]));
    // The byte-layout stream gives tight bounds for every rANS stream. Probe it before allocating
    // decoded vectors so a corrupt count cannot amplify a small chunk into unbounded work.
    let mut layout_probe = Cursor::new(streams[2]);
    if info.n_mols as usize > streams[2].len() / 3 {
        bail!("layout stream is too short for {} molecules", info.n_mols);
    }
    let (mut n_chains_total, mut n_mms_total, mut mols_with_chains) = (0usize, 0usize, 0usize);
    for _ in 0..info.n_mols {
        layout_probe.byte()?;
        let n_chains = usize::try_from(layout_probe.varint()?).context("chain count is too large")?;
        let n_mms = usize::try_from(layout_probe.varint()?).context("multimapper count is too large")?;
        n_chains_total = n_chains_total.checked_add(n_chains).context("chain count overflow")?;
        n_mms_total = n_mms_total.checked_add(n_mms).context("multimapper count overflow")?;
        mols_with_chains += usize::from(n_chains > 0);
    }
    if !layout_probe.is_empty() {
        bail!("layout stream has trailing bytes");
    }
    // rANS streams pre-decode fully; the chunk is the access unit so this is the working set.
    let class_v = evidence_io::rans::decode_limited(streams[1], &tables[0], info.n_mols as usize)?;
    let weight_v = evidence_io::rans::decode_limited(streams[3], &tables[1], n_chains_total)?;
    if weight_v.len() != n_chains_total {
        bail!("weight stream has {} values for {n_chains_total} chains", weight_v.len());
    }
    let n_reps_total = weight_v.iter().try_fold(0usize, |sum, w| {
        sum.checked_add(if w & 1 == 1 { 2 } else { 1 })
            .context("representative count overflow")
    })?;
    let n_rep_positions = n_reps_total
        .checked_sub(mols_with_chains)
        .context("representative elision count underflow")?;
    let rep_pos_v =
        evidence_io::rans::decode_limited(streams[4], &tables[2], n_rep_positions)?;
    let mm_pos_v = evidence_io::rans::decode_limited(streams[6], &tables[3], n_mms_total)?;
    let mm_w_v = evidence_io::rans::decode_limited(streams[9], &tables[4], n_mms_total)?;
    if class_v.len() != info.n_mols as usize {
        bail!(
            "class stream has {} values for {} molecules",
            class_v.len(),
            info.n_mols
        );
    }
    let (mut ci, mut wi, mut rpi, mut mpi, mut mwi) = (0usize, 0usize, 0usize, 0usize, 0usize);

    let mut mols = Vec::with_capacity(info.n_mols as usize);
    let (mut lanchor, mut next_class) = (info.bin_start, info.class_base);
    for _ in 0..info.n_mols {
        let anchor_delta = u32::try_from(anchor_c.varint()?).context("anchor delta exceeds u32")?;
        let anchor = lanchor.checked_add(anchor_delta).context("anchor coordinate overflow")?;
        lanchor = anchor;
        let ctok = *class_v.get(ci).context("class stream underrun")?;
        ci += 1;
        let umi_class = if ctok == 0 {
            let id = next_class;
            next_class = next_class.checked_add(1).context("UMI class id overflow")?;
            id
        } else {
            let backref = u32::try_from(ctok).context("UMI class back-reference exceeds u32")?;
            next_class.checked_sub(backref).context("UMI class back-reference underflow")?
        };
        let cell = match cell_of_class {
            Some(t) => *t
                .get(umi_class as usize)
                .ok_or_else(|| anyhow::anyhow!("class {umi_class} beyond cellofclass table"))?,
            None => u32::MAX,
        };
        let strand_rev = layout_c.byte()? != 0;
        let n_chains = usize::try_from(layout_c.varint()?).context("chain count is too large")?;
        let elide_first_rep = n_chains > 0;
        let n_mms = usize::try_from(layout_c.varint()?).context("multimapper count is too large")?;
        if n_chains > weight_v.len().saturating_sub(wi) {
            bail!("chunk declares more chains than the weight stream contains");
        }
        if n_mms > mm_pos_v.len().saturating_sub(mpi)
            || n_mms > mm_w_v.len().saturating_sub(mwi)
        {
            bail!("chunk declares more multimappers than its value streams contain");
        }
        let mut chains: SmallVec<[MolChain; 1]> = SmallVec::with_capacity(n_chains);
        let mut first_rep = true;
        for _ in 0..n_chains {
            let wv = *weight_v.get(wi).context("chain-weight stream underrun")?;
            wi += 1;
            let weight = u32::try_from(wv >> 1).context("chain weight exceeds u32")?;
            let two = wv & 1 == 1;
            let mut reps: SmallVec<[(u32, u32); 2]> = SmallVec::new();
            for _ in 0..if two { 2 } else { 1 } {
                let pos = if first_rep && elide_first_rep {
                    anchor
                } else {
                    let delta = u32::try_from(
                        *rep_pos_v.get(rpi).context("representative-position stream underrun")?,
                    )
                    .context("representative-position delta exceeds u32")?;
                    rpi += 1;
                    anchor.checked_add(delta).context("representative position overflow")?
                };
                first_rep = false;
                let shape = u32::try_from(rep_shape_c.varint()?)
                    .context("representative shape id exceeds u32")?;
                reps.push((pos, shape));
            }
            chains.push(MolChain { weight, reps });
        }
        let mut mms: SmallVec<[(u32, u32, u32, u32); 1]> = SmallVec::with_capacity(n_mms);
        for _ in 0..n_mms {
            let delta = u32::try_from(
                *mm_pos_v.get(mpi).context("multimapper-position stream underrun")?,
            )
            .context("multimapper-position delta exceeds u32")?;
            let mp = anchor.checked_add(delta).context("multimapper position overflow")?;
            mpi += 1;
            let mw = u32::try_from(
                *mm_w_v.get(mwi).context("multimapper-weight stream underrun")?,
            )
            .context("multimapper weight exceeds u32")?;
            mwi += 1;
            let shape = u32::try_from(mm_shape_c.varint()?)
                .context("multimapper shape id exceeds u32")?;
            let pattern = u32::try_from(mm_pat_c.varint()?)
                .context("multimapper pattern id exceeds u32")?;
            mms.push((mp, shape, pattern, mw));
        }
        mols.push(MolRec { cell, umi_class, chrom: info.chrom, strand_rev, chains, mms });
    }
    if ci != class_v.len()
        || wi != weight_v.len()
        || rpi != rep_pos_v.len()
        || mpi != mm_pos_v.len()
        || mwi != mm_w_v.len()
    {
        bail!(
            "chunk value-stream length mismatch: class {ci}/{}, weight {wi}/{}, rep.pos {rpi}/{}, mm.pos {mpi}/{}, mm.weight {mwi}/{}",
            class_v.len(),
            weight_v.len(),
            rep_pos_v.len(),
            mm_pos_v.len(),
            mm_w_v.len()
        );
    }
    for (name, cursor) in [
        ("anchor", &anchor_c),
        ("layout", &layout_c),
        ("rep.shape", &rep_shape_c),
        ("mm.shape", &mm_shape_c),
        ("mm.pattern", &mm_pat_c),
    ] {
        if !cursor.is_empty() {
            bail!("chunk {name} stream has trailing bytes");
        }
    }
    Ok(mols)
}

/// Lazy archive access for queries: opens with meta + chrom names only; every dictionary loads on
/// first use, and cell-of-class blocks decompress individually as classes are touched. This keeps
/// query-open cost flat as archives grow; eager loading took 0.69 s on an approximately
/// 100-million-molecule benchmark archive.
pub struct LazyArchive {
    r: SectionReader,
    pub chrom_names: Vec<String>,
    pub chrom_digest: String,
    /// Present iff the archive was stamped with the reference genome's signature.
    pub genome_sig: Option<evidence_io::genome::GenomeSig>,
    /// Loaded eagerly: every chunk decode needs them and they are ~4 KB.
    pub rans_tables: Vec<evidence_io::rans::Table>,
    n_classes: u32,
    coc_cache: HbMap<u32, Vec<u32>>,
    cells: Option<Vec<u32>>,
    /// Arc so query code can hold the dictionary while the reader is borrowed mutably —
    /// previously callers deep-cloned the whole shape table per query.
    shapes: Option<std::sync::Arc<Vec<Shape>>>,
    /// Loaded only for queries that expand multimapper alternatives.
    patterns: Option<std::sync::Arc<Vec<Vec<PatAlt>>>>,
    /// `None` means the extraction rule was not run and callers must report unavailable, never
    /// biological zero. An empty, present capability is represented by zero metadata counts.
    terminal_tail: Option<TerminalTailMetadata>,
    terminal_tail_routes: Option<Vec<TerminalTailRoute>>,
}

fn archive_capabilities(
    reader: &mut SectionReader,
    meta: &serde_json::Value,
    verify_declared_payloads: bool,
) -> Result<(
    Option<AlignmentProvenanceManifest>,
    Option<TerminalTailMetadata>,
    Option<GenomeReferenceBinding>,
)> {
    let has_provenance_section = reader.has(ALIGNMENT_PROVENANCE_SECTION);
    let has_catalogue_section = reader.has(JUNCTION_CATALOGUE_SECTION);
    let has_tail_index = reader.has(terminal_tail::TERMINAL_TAIL_INDEX_SECTION);
    let tail_section_count = reader
        .names()
        .filter(|name| name.starts_with("tail.c"))
        .count();
    let evidence_schema = meta.get("evidence_schema");
    let provenance_schema = meta.get("alignment_provenance");
    let tail_value = meta.get("terminal_tail");

    let alignment_provenance = match evidence_schema {
        None => {
            if provenance_schema.is_some()
                || has_provenance_section
                || has_catalogue_section
                || tail_value.is_some()
                || meta.get("genome_reference_binding").is_some()
                || has_tail_index
                || tail_section_count != 0
            {
                bail!("archive has partial logical-v2 capability sections without an evidence_schema declaration");
            }
            None
        }
        Some(value) if value.as_str() == Some(MOLECULAR_EVIDENCE_SCHEMA) => {
            if reader.content_commitment().is_none() {
                bail!("logical molecular-evidence v2 requires a root-committed archive container");
            }
            check_layout(meta)?;
            if provenance_schema.and_then(serde_json::Value::as_str)
                != Some(ALIGNMENT_PROVENANCE_SCHEMA)
                || !has_provenance_section
            {
                bail!("logical molecular-evidence v2 lacks its declared alignment provenance");
            }
            let raw = reader.read(ALIGNMENT_PROVENANCE_SECTION)?;
            let manifest = AlignmentProvenanceManifest::from_json(&raw)?;
            check_provenance_layout(meta, &manifest.ingest)?;
            match manifest.alignment.junction_catalogue.as_ref() {
                Some(catalogue) => {
                    if !has_catalogue_section {
                        bail!("alignment provenance declares a junction catalogue whose exact-byte section is absent");
                    }
                    let raw_len = reader
                        .section_metadata()
                        .find(|section| section.name == JUNCTION_CATALOGUE_SECTION)
                        .map(|section| section.raw_len)
                        .expect("catalogue section presence was checked");
                    if raw_len != catalogue.identity.bytes {
                        bail!("junction catalogue section length disagrees with its provenance identity");
                    }
                    if verify_declared_payloads {
                        let catalogue_bytes = reader.read(JUNCTION_CATALOGUE_SECTION)?;
                        if blake3::hash(&catalogue_bytes).to_hex().as_str()
                            != catalogue.identity.blake3
                        {
                            bail!("junction catalogue section digest disagrees with its provenance identity");
                        }
                        let rows = junction_catalogue_data_rows(
                            &catalogue_bytes,
                            JUNCTION_CATALOGUE_SECTION,
                        )?;
                        if rows != catalogue.data_rows {
                            bail!("junction catalogue parsed row count disagrees with its provenance manifest");
                        }
                    }
                }
                None if has_catalogue_section => {
                    bail!("archive contains an undeclared alignment junction catalogue section")
                }
                None => {}
            }
            Some(manifest)
        }
        Some(value) => bail!(
            "unsupported molecular-evidence schema {}",
            value.as_str().unwrap_or("<non-string>")
        ),
    };

    let terminal_tail = match tail_value {
        Some(value) => {
            if alignment_provenance.is_none() || !has_tail_index {
                bail!("terminal-tail capability lacks its logical-v2 manifest or sparse index");
            }
            let metadata: TerminalTailMetadata = serde_json::from_value(value.clone())?;
            metadata.validate()?;
            if tail_section_count != metadata.chunks as usize {
                bail!(
                    "terminal-tail metadata declares {} routed chunks but archive contains {tail_section_count}",
                    metadata.chunks
                );
            }
            Some(metadata)
        }
        None => {
            if has_tail_index || tail_section_count != 0 {
                bail!("archive has terminal-tail sections without a capability declaration");
            }
            None
        }
    };
    if let Some(manifest) = &alignment_provenance {
        let rule_declared = manifest.ingest.terminal_tail_rule.as_deref()
            == Some(terminal_tail::TERMINAL_TAIL_RULE);
        if rule_declared != terminal_tail.is_some() {
            bail!("terminal-tail capability and provenance extraction rule disagree");
        }
    }
    if verify_declared_payloads {
        if let Some(metadata) = &terminal_tail {
            let chunks = read_chunk_index(reader)?;
            let routes = terminal_tail::decode_index(
                &reader.read(terminal_tail::TERMINAL_TAIL_INDEX_SECTION)?,
                chunks.len(),
            )?;
            let selected_molecules = routes.iter().try_fold(0u64, |sum, route| {
                sum.checked_add(u64::from(route.selected_molecules))
                    .context("terminal-tail selected-molecule total overflow")
            })?;
            let events = routes.iter().try_fold(0u64, |sum, route| {
                sum.checked_add(u64::from(route.events))
                    .context("terminal-tail event total overflow")
            })?;
            if routes.len() != metadata.chunks as usize
                || selected_molecules != metadata.selected_molecules
                || events != metadata.events
            {
                bail!("terminal-tail index cardinalities disagree with capability metadata");
            }
            let routed_sections: BTreeSet<String> = routes
                .iter()
                .map(|route| format!("tail.c{}", route.chunk))
                .collect();
            let actual_sections: BTreeSet<String> = reader
                .names()
                .filter(|name| name.starts_with("tail.c"))
                .map(str::to_owned)
                .collect();
            if routed_sections != actual_sections {
                bail!("terminal-tail sparse index and routed sections disagree");
            }
            for route in routes {
                let info = chunks
                    .get(route.chunk as usize)
                    .context("terminal-tail route references an absent molecule chunk")?;
                if route.chrom != info.chrom {
                    bail!("terminal-tail route chromosome disagrees with its molecule chunk");
                }
                let decoded = terminal_tail::decode_chunk(
                    &reader.read(&format!("tail.c{}", route.chunk))?,
                    info.n_mols,
                )?;
                let decoded_events = decoded.iter().try_fold(0u32, |sum, molecule| {
                    sum.checked_add(u32::try_from(molecule.events.len()).map_err(|_| {
                        anyhow::anyhow!("terminal-tail per-molecule event count exceeds u32")
                    })?)
                    .context("terminal-tail routed event count overflow")
                })?;
                if decoded.len() != route.selected_molecules as usize
                    || decoded_events != route.events
                {
                    bail!("terminal-tail routed section cardinalities disagree with its index");
                }
            }
        }
    }
    let meta_genome_signature = meta
        .get("genome_sig")
        .map(|value| serde_json::from_value::<evidence_io::genome::GenomeSig>(value.clone()))
        .transpose()?;
    let genome_reference_binding = match meta.get("genome_reference_binding") {
        Some(value) => {
            let binding: GenomeReferenceBinding = serde_json::from_value(value.clone())?;
            binding.validate()?;
            if alignment_provenance.is_none()
                || meta_genome_signature.as_ref() != Some(&binding.signature)
            {
                bail!("genome reference binding lacks matching logical-v2 metadata");
            }
            if binding.bound_by == GenomeBindingAction::IngestArchive {
                let manifest = alignment_provenance
                    .as_ref()
                    .expect("logical-v2 binding checked above");
                if manifest.inputs.genome_fasta.as_ref() != Some(&binding.identity)
                    || manifest.inputs.genome_signature.as_ref() != Some(&binding.signature)
                {
                    bail!("ingest-time genome binding disagrees with alignment provenance");
                }
            }
            Some(binding)
        }
        None => {
            if alignment_provenance.is_some() && meta_genome_signature.is_some() {
                bail!("logical-v2 genome signature lacks its reference-binding declaration");
            }
            if alignment_provenance
                .as_ref()
                .is_some_and(|manifest| manifest.inputs.genome_signature.is_some())
            {
                bail!("alignment provenance contains an ingest genome but no current reference binding");
            }
            None
        }
    };
    Ok((
        alignment_provenance,
        terminal_tail,
        genome_reference_binding,
    ))
}

impl LazyArchive {
    pub fn open(path: &Path) -> Result<LazyArchive> {
        let mut r = SectionReader::open(path)?;
        let meta: serde_json::Value = serde_json::from_slice(&r.read("meta")?)?;
        check_layout(&meta)?;
        let (_, terminal_tail, _) = archive_capabilities(&mut r, &meta, false)?;
        let chrom_bytes = r.read("chroms")?;
        let chrom_text = std::str::from_utf8(&chrom_bytes)
            .context("archive chromosome dictionary is not UTF-8")?;
        let chrom_names = chrom_text.lines().map(|s| s.to_string()).collect();
        let chrom_digest = blake3::hash(&chrom_bytes).to_hex().to_string();
        let rans_tables = read_rans_tables(&mut r)?;
        let genome_sig = meta
            .get("genome_sig")
            .and_then(|v| serde_json::from_value(v.clone()).ok());
        Ok(LazyArchive {
            r,
            chrom_names,
            chrom_digest,
            genome_sig,
            rans_tables,
            n_classes: required_meta_u32(&meta, "classes")?,
            coc_cache: HbMap::new(),
            cells: None,
            shapes: None,
            patterns: None,
            terminal_tail,
            terminal_tail_routes: None,
        })
    }

    /// Whether the terminal-tail extraction rule was evaluated. `Some` with zero events is a
    /// measured zero; `None` is unavailable evidence and must fail a requested tail predicate.
    pub fn terminal_tail_capability(&self) -> Option<&TerminalTailMetadata> {
        self.terminal_tail.as_ref()
    }

    /// Sparse terminal-tail routes, decoded only when requested.
    pub fn terminal_tail_routes(&mut self) -> Result<Option<&[TerminalTailRoute]>> {
        let Some(metadata) = self.terminal_tail.as_ref() else {
            return Ok(None);
        };
        if self.terminal_tail_routes.is_none() {
            let chunks = read_chunk_index(&mut self.r)?;
            let raw = self.r.read(terminal_tail::TERMINAL_TAIL_INDEX_SECTION)?;
            let routes = terminal_tail::decode_index(&raw, chunks.len())?;
            let selected = routes.iter().try_fold(0u64, |sum, route| {
                sum.checked_add(u64::from(route.selected_molecules))
                    .context("terminal-tail selected-molecule total overflow")
            })?;
            let events = routes.iter().try_fold(0u64, |sum, route| {
                sum.checked_add(u64::from(route.events))
                    .context("terminal-tail event total overflow")
            })?;
            if routes.len() != metadata.chunks as usize
                || selected != metadata.selected_molecules
                || events != metadata.events
            {
                bail!("terminal-tail index cardinalities disagree with capability metadata");
            }
            let routed_sections: BTreeSet<String> = routes
                .iter()
                .map(|route| format!("tail.c{}", route.chunk))
                .collect();
            let actual_sections: BTreeSet<String> = self
                .r
                .names()
                .filter(|name| name.starts_with("tail.c"))
                .map(str::to_owned)
                .collect();
            if routed_sections != actual_sections {
                bail!("terminal-tail sparse index and routed sections disagree");
            }
            self.terminal_tail_routes = Some(routes);
        }
        Ok(self.terminal_tail_routes.as_deref())
    }

    /// Decode one routed tail section and attach it to its already-decoded ordinary molecule
    /// chunk. This validates chromosome, strand, route envelope, counts, and local identity before
    /// returning query-friendly records; corrupt or mismatched side evidence therefore fails
    /// closed rather than being silently detached from its molecule.
    pub fn terminal_tail_records(
        &mut self,
        route: TerminalTailRoute,
        info: &ChunkInfo,
        molecule_base: u64,
        molecules: &[MolRec],
    ) -> Result<Vec<TerminalTailRecord>> {
        let known_route = self
            .terminal_tail_routes()?
            .context("terminal-tail capability is unavailable")?
            .contains(&route);
        if !known_route {
            bail!("terminal-tail route is not present in the archive index");
        }
        if route.chrom != info.chrom || molecules.len() != info.n_mols as usize {
            bail!("terminal-tail route does not match its ordinary molecule chunk");
        }
        let raw = self.r.read(&format!("tail.c{}", route.chunk))?;
        let decoded = terminal_tail::decode_chunk(&raw, info.n_mols)?;
        if decoded.len() != route.selected_molecules as usize {
            bail!("terminal-tail routed selected-molecule count mismatch");
        }
        let mut records = Vec::with_capacity(route.events as usize);
        let mut min_anchor = u32::MAX;
        let mut max_anchor = 0u32;
        for selected in decoded {
            let molecule = &molecules[selected.local_ordinal as usize];
            if molecule.chrom != info.chrom {
                bail!("terminal-tail event is attached to a molecule on another chromosome");
            }
            let cell = if molecule.cell == u32::MAX {
                self.cell_of(molecule.umi_class)?
            } else {
                molecule.cell
            };
            for event in selected.events {
                if event.reverse != molecule.strand_rev {
                    bail!("terminal-tail event strand disagrees with its attached molecule");
                }
                let anchor = i64::from(molecule.anchor())
                    .checked_add(event.anchor_delta)
                    .and_then(|value| u32::try_from(value).ok())
                    .context("terminal-tail anchor is outside the u32 coordinate range")?;
                min_anchor = min_anchor.min(anchor);
                max_anchor = max_anchor.max(anchor);
                records.push(TerminalTailRecord {
                    molecule_ordinal: molecule_base
                        .checked_add(u64::from(selected.local_ordinal))
                        .context("terminal-tail global molecule ordinal overflow")?,
                    chunk: route.chunk,
                    local_molecule_ordinal: selected.local_ordinal,
                    cell,
                    umi_class: molecule.umi_class,
                    chrom: molecule.chrom,
                    strand_rev: molecule.strand_rev,
                    anchor,
                    signal: event.signal,
                });
            }
        }
        if records.len() != route.events as usize
            || min_anchor != route.min_anchor
            || max_anchor != route.max_anchor
        {
            bail!("terminal-tail route envelope or event count disagrees with its section");
        }
        Ok(records)
    }

    pub fn reader(&mut self) -> &mut SectionReader {
        &mut self.r
    }

    /// Split borrow for parallel chunk decode: the reader (mutable, for compressed reads) and the
    /// rANS tables (shared, for the decode workers) at the same time.
    pub fn reader_and_tables(&mut self) -> (&mut SectionReader, &[evidence_io::rans::Table]) {
        (&mut self.r, &self.rans_tables)
    }

    /// Packed barcode per cell id (loads the cells dict on first call).
    pub fn cells(&mut self) -> Result<&[u32]> {
        if self.cells.is_none() {
            let raw = self.r.read("cells")?;
            if raw.len() % 4 != 0 {
                bail!("cells section has {} trailing byte(s)", raw.len() % 4);
            }
            self.cells = Some(
                raw.as_chunks::<4>()
                    .0
                    .iter()
                    .map(|&bytes| u32::from_le_bytes(bytes))
                    .collect(),
            );
        }
        Ok(self.cells.as_deref().unwrap())
    }

    pub fn shapes(&mut self) -> Result<std::sync::Arc<Vec<Shape>>> {
        if self.shapes.is_none() {
            let raw = self.r.read("shapes")?;
            self.shapes = Some(std::sync::Arc::new(decode_shapes(&raw)?));
        }
        Ok(self.shapes.clone().unwrap())
    }

    pub fn patterns(&mut self) -> Result<std::sync::Arc<Vec<Vec<PatAlt>>>> {
        if self.patterns.is_none() {
            let raw = self.r.read("patterns")?;
            self.patterns = Some(std::sync::Arc::new(decode_patterns(&raw)?));
        }
        Ok(self.patterns.clone().unwrap())
    }

    /// Decode (in parallel) every cell-of-class block the given classes touch, into the cache.
    /// `cell_of` then serves the callers' per-class lookups without serial block decodes — the
    /// serial path was 60%+ of a junction query's wall (rANS over 64Ki ids per touched block).
    pub fn prefetch_coc(&mut self, classes: impl Iterator<Item = u32>) -> Result<()> {
        let mut need: Vec<u32> = classes.map(|c| c / COC_BLOCK).collect();
        need.sort_unstable();
        need.dedup();
        need.retain(|b| !self.coc_cache.contains_key(b));
        let (reader, tables) = (&self.r, &self.rans_tables);
        let decoded: Vec<(u32, Vec<u32>)> = need
            .par_iter()
            .map(|b| {
                let (comp, raw_len) = reader.read_compressed_at(&format!("coc.{b}"))?;
                Ok((*b, decode_coc_block(&comp, raw_len, &tables[5])?))
            })
            .collect::<Result<_>>()?;
        for (b, block) in decoded {
            self.coc_cache.insert(b, block);
        }
        Ok(())
    }

    /// Release decoded cell-of-class blocks between bounded full-scan batches. Targeted queries
    /// benefit from retaining this cache, while a cohort reducer can safely re-prefetch a block
    /// if a later chromosome touches it again.
    pub fn clear_coc_cache(&mut self) {
        self.coc_cache.clear();
    }

    /// Cell id for a class whose COC block has already been prefetched. This read-only variant
    /// lets full-scan reducers classify independent decoded chunks in parallel without putting a
    /// lock around the lazy cache.
    pub fn cell_of_cached(&self, class: u32) -> Result<u32> {
        if class >= self.n_classes {
            bail!("class {class} beyond {} classes", self.n_classes);
        }
        let block = class / COC_BLOCK;
        self.coc_cache
            .get(&block)
            .and_then(|values| values.get((class % COC_BLOCK) as usize))
            .copied()
            .with_context(|| {
                format!("cell-of-class block {block} was not prefetched for class {class}")
            })
    }

    /// Cell id for one class, decoding only that class's coc block.
    pub fn cell_of(&mut self, class: u32) -> Result<u32> {
        if class >= self.n_classes {
            bail!("class {class} beyond {} classes", self.n_classes);
        }
        let b = class / COC_BLOCK;
        if !self.coc_cache.contains_key(&b) {
            let (comp, raw_len) = self.r.read_compressed(&format!("coc.{b}"))?;
            let block = decode_coc_block(&comp, raw_len, &self.rans_tables[5])?;
            self.coc_cache.insert(b, block);
        }
        self.coc_cache
            .get(&b)
            .and_then(|block| block.get((class % COC_BLOCK) as usize))
            .copied()
            .with_context(|| format!("cell-of-class block {b} is truncated for class {class}"))
    }
}

fn read_edges(r: &mut SectionReader) -> Result<Vec<(u32, u32)>> {
    let raw = r.read("edges")?;
    let mut c = Cursor::new(&raw);
    let mut edges = Vec::new();
    let mut la = 0i64;
    while !c.is_empty() {
        let a = (la + c.svarint()?) as u32;
        la = a as i64;
        let b = a + c.varint()? as u32;
        edges.push((a, b));
    }
    Ok(edges)
}

/// Dictionaries and directory for bounded Gene replay. `x.mols` stays empty: decoded records are
/// classified and released batch by batch.
struct StreamingReplayArchive {
    reader: SectionReader,
    chunks: Vec<ChunkInfo>,
    x: Extracted,
    cell_of_class: Vec<u32>,
    rans_tables: Vec<evidence_io::rans::Table>,
    n_mols: usize,
}

type StreamingReplayResult = (FxHashMap<(u32, u32), u32>, u64, u64);

impl StreamingReplayArchive {
    fn open(path: &Path) -> Result<Self> {
        let mut reader = SectionReader::open(path)?;
        let d = read_dicts(&mut reader)?;
        let chunks = read_chunk_index(&mut reader)?;
        let edges = read_edges(&mut reader)?;
        let Dicts { cells, shapes, patterns, chrom_names, cell_of_class, rans_tables,
            n_classes, n_mols } = d;
        let indexed_mols: usize = chunks.iter().map(|c| c.n_mols as usize).sum();
        if indexed_mols != n_mols {
            bail!("molecule count mismatch: {indexed_mols} indexed vs {n_mols} in meta");
        }
        Ok(Self {
            reader, chunks,
            x: Extracted { mols: Vec::new(), edges, cells, shapes, patterns, n_classes,
                chrom_names },
            cell_of_class, rans_tables, n_mols,
        })
    }

    fn replay(
        &self,
        anno: &anno::Annotation,
        solo_strand: anno::assign::SoloStrand,
    ) -> Result<StreamingReplayResult> {
        let mut replay = ReplayRowsAccumulator::with_strand(&self.x, anno, solo_strand);
        replay.reserve_assignments(self.n_mols);
        // Two access units per worker smooth decode density and amortize reducer barriers while
        // retaining bounded memory use in complete-command peak-RSS measurements.
        let batch_size = rayon::current_num_threads().max(1) * 2;
        let mut decoded_mols = 0usize;
        for (batch_no, batch) in self.chunks.chunks(batch_size).enumerate() {
            let first = batch_no * batch_size;
            let decoded: Vec<Vec<MolRec>> = batch.par_iter().enumerate().map(|(j, info)| {
                let i = first + j;
                let (comp, raw_len) = self.reader.read_compressed_at(&format!("c{i}"))?;
                let raw = evidence_io::format::decompress(&comp, raw_len)?;
                decode_chunk(&raw, info, Some(&self.cell_of_class), &self.rans_tables)
            }).collect::<Result<_>>()?;
            decoded_mols += decoded.iter().map(Vec::len).sum::<usize>();
            replay.add_archive_chunks(&decoded);
        }
        if decoded_mols != self.n_mols {
            bail!("molecule count mismatch: {decoded_mols} decoded vs {} in meta", self.n_mols);
        }
        Ok(replay.finish())
    }
    /// Decode the archive in the same bounded batches as streaming Gene replay, but reduce each
    /// batch immediately to packed (class,gene,evidence-kind) support words.  The final EM owns
    /// only CSR labels and shard-local numeric state; `x.mols` remains empty throughout.
    fn packed_em_accumulator(
        &self,
        anno: &anno::Annotation,
    ) -> Result<crate::rows::PackedEmAccumulator> {
        let mut em = crate::rows::PackedEmAccumulator::new(&self.x, anno, self.n_mols)?;
        let threads = rayon::current_num_threads().max(1);
        // The disk-backed path exists to bound large-archive memory. Eight access units keep its
        // decoded MolRec window small; classification still subdivides them into 65k tasks and
        // can use the full Rayon pool. Small in-memory archives retain the faster, benchmarked
        // two-per-thread decode window.
        let batch_size = if em.spills_supports() {
            threads.min(8)
        } else {
            threads * 2
        };
        let mut decoded_mols = 0usize;
        for (batch_no, batch) in self.chunks.chunks(batch_size).enumerate() {
            let first = batch_no * batch_size;
            let decoded: Vec<Vec<MolRec>> = batch
                .par_iter()
                .enumerate()
                .map(|(j, info)| {
                    let i = first + j;
                    let (comp, raw_len) = self.reader.read_compressed_at(&format!("c{i}"))?;
                    let raw = evidence_io::format::decompress(&comp, raw_len)?;
                    decode_chunk(&raw, info, Some(&self.cell_of_class), &self.rans_tables)
                })
                .collect::<Result<_>>()?;
            decoded_mols += decoded.iter().map(Vec::len).sum::<usize>();
            em.add_archive_chunks(&decoded, &self.x, anno)?;
        }
        if decoded_mols != self.n_mols {
            bail!(
                "molecule count mismatch: {decoded_mols} decoded vs {} in meta",
                self.n_mols
            );
        }
        Ok(em)
    }
}

pub fn read_archive(path: &Path) -> Result<Extracted> {
    let mut r = SectionReader::open(path)?;
    read_archive_from_reader(&mut r)
}

/// Decode the complete archive and bind provenance to that exact open reader. Rooted v2 identity
/// is directory-only; the additional legacy full-file scan occurs only on this opt-in path.
pub(crate) fn read_archive_with_identity(
    path: &Path,
) -> Result<(Extracted, ArchiveContentIdentity)> {
    let mut reader = SectionReader::open(path)?;
    let extracted = read_archive_from_reader(&mut reader)?;
    let identity = archive_identity(&reader)?;
    Ok((extracted, identity))
}

fn read_archive_from_reader(r: &mut SectionReader) -> Result<Extracted> {
    let d = read_dicts(r)?;
    let chunks = read_chunk_index(r)?;

    // Parallel pread + decompress + decode: one fused worker per chunk.
    let per_chunk: Vec<Vec<MolRec>> = (0..chunks.len())
        .into_par_iter()
        .map(|i| {
            let (c, raw_len) = r.read_compressed_at(&format!("c{i}"))?;
            let raw = evidence_io::format::decompress(&c, raw_len)?;
            decode_chunk(&raw, &chunks[i], Some(&d.cell_of_class), &d.rans_tables)
        })
        .collect::<Result<_>>()?;
    let total: usize = per_chunk.iter().map(|v| v.len()).sum();
    if total != d.n_mols {
        bail!("molecule count mismatch: {total} decoded vs {} in meta", d.n_mols);
    }
    // Parallel move-concat into one Vec: on an approximately 100-million-molecule benchmark
    // archive, serial `extend` walked about 9 GB of records.
    // Disjoint destination ranges per chunk; sources are forgotten without dropping the moved-out
    // elements (their allocations are still freed).
    let mut offsets = Vec::with_capacity(per_chunk.len());
    let mut acc = 0usize;
    for v in &per_chunk {
        offsets.push(acc);
        acc += v.len();
    }
    let mut mols: Vec<MolRec> = Vec::with_capacity(total);
    {
        let base_addr = mols.as_mut_ptr() as usize;
        per_chunk.into_par_iter().zip(offsets).for_each(|(v, off)| {
            let mut v = std::mem::ManuallyDrop::new(v);
            unsafe {
                std::ptr::copy_nonoverlapping(v.as_ptr(), (base_addr as *mut MolRec).add(off), v.len());
                // Free the source allocation without dropping the moved-out records.
                let (ptr, cap) = (v.as_mut_ptr(), v.capacity());
                drop(Vec::from_raw_parts(ptr as *mut std::mem::MaybeUninit<MolRec>, 0, cap));
            }
        });
        unsafe { mols.set_len(total) };
    }

    let edges = read_edges(r)?;

    Ok(Extracted {
        mols,
        edges,
        cells: d.cells,
        shapes: d.shapes,
        patterns: d.patterns,
        n_classes: d.n_classes,
        chrom_names: d.chrom_names,
    })
}

#[derive(Parser)]
pub struct EmArgs {
    pub archive: PathBuf,
    #[arg(long)]
    pub gtf: PathBuf,
    /// Fraction of mixed classes to mask for the labeled evaluation.
    #[arg(long, default_value_t = 0.2)]
    pub mask: f64,
    #[arg(long, default_value_t = 7)]
    pub seed: u64,
    /// Global-prior weight for blend mode.
    #[arg(long, default_value_t = 20.0)]
    pub alpha: f64,
    /// Two-column BARCODE/GROUP map enabling group-only and hierarchical evaluation modes.  Cells
    /// absent from the map retain the global fallback and are excluded from all mode scores.
    #[arg(long)]
    pub groups: Option<PathBuf>,
    /// Hierarchical pseudo-count mass borrowed from the supplied group distribution.
    #[arg(long, default_value_t = 20.0)]
    pub group_alpha: f64,
    /// Hierarchical pseudo-count mass borrowed from the whole-sample distribution.
    #[arg(long, default_value_t = 5.0)]
    pub global_alpha: f64,
    /// Candidate-normalized convex weight on the target cell distribution.
    #[arg(long, default_value_t = 0.10)]
    pub convex_cell_weight: f64,
    /// Candidate-normalized convex weight on the leave-one-cell-out group distribution. The
    /// whole-sample weight is the simplex remainder after cell and group weights.
    #[arg(long, default_value_t = 0.45)]
    pub convex_group_weight: f64,
    /// Candidate-level sample pseudo-count mass used to shrink the leave-one-cell-out group
    /// distribution. Unlike --group-alpha, this mass is applied after candidate normalization.
    #[arg(long, default_value_t = 20.0)]
    pub convex_group_prior: f64,
    /// Evaluate only the candidate-normalized convex mode. This avoids recomputing the eight
    /// reference modes during weight/prior grids and cannot be used for production emission.
    #[arg(long, requires = "groups", conflicts_with = "emit")]
    pub convex_only: bool,
    /// Candidate-level group-posterior mass borrowed by each cell in the posterior-mean
    /// hierarchical Dirichlet screening approximation.
    #[arg(long, default_value_t = 16.0)]
    pub dirichlet_cell_prior: f64,
    /// Candidate-level sample-posterior mass borrowed by each leave-one-cell-out group in the
    /// posterior-mean hierarchical Dirichlet screening approximation.
    #[arg(long, default_value_t = 20.0)]
    pub dirichlet_group_prior: f64,
    /// Evaluate only the posterior-mean Dirichlet proxy. This screens evidence-adaptive pooling
    /// before any full variational model is implemented and cannot emit production counts.
    #[arg(long, requires = "groups", conflicts_with_all = ["emit", "convex_only"])]
    pub dirichlet_only: bool,
    /// Half-transition evidence depth for hybrid weighting. The weight is
    /// depth^power / (depth^power + scale^power), using fitted unique candidate mass.
    #[arg(long, default_value_t = 8.0)]
    pub hybrid_depth_scale: f64,
    /// Positive Hill power controlling the sharpness of the monotone hybrid transition.
    #[arg(long, default_value_t = 8.0)]
    pub hybrid_depth_power: f64,
    /// Evaluate only the preferred monotone evidence-depth hybrid. This masked-recovery mode
    /// cannot emit production counts.
    #[arg(
        long,
        requires = "groups",
        conflicts_with_all = ["emit", "convex_only", "dirichlet_only"]
    )]
    pub hybrid_only: bool,
    /// Assign every archive cell to one group, but retain the --groups barcodes as the evaluation
    /// population. This is the registered collapse-to-pooled control.
    #[arg(long, requires = "groups")]
    pub collapse_groups: bool,
    /// Write proper-score and calibration metrics as JSON.
    #[arg(long)]
    pub metrics_json: Option<PathBuf>,
    /// Write the sorted union of genes occurring in masked evaluation candidate sets.
    #[arg(long)]
    pub candidate_genes_out: Option<PathBuf>,
    /// Stop after writing --candidate-genes-out, without running EM iterations.
    #[arg(long, requires = "candidate_genes_out")]
    pub candidate_genes_only: bool,
    /// With --mask 0: write the additive fractional recovered-counts layer (em.mtx, real-valued)
    /// into this directory. Requires --barcodes for column order.
    #[arg(long)]
    pub emit: Option<PathBuf>,
    /// Barcode list defining emitted column order (for example, STARsolo raw barcodes.tsv).
    #[arg(long)]
    pub barcodes: Option<PathBuf>,
    /// EM-0: replicate STARsolo --soloMultiMappers EM exactly (per-cell, intersection candidate
    /// sets, STAR's init/zeroing/convergence) and emit UniqueAndMult-EM.mtx into --emit.
    #[arg(long)]
    pub star: bool,
    /// Use the v1 eager row/vector/hash-table implementation as a semantic and performance
    /// reference.  The default recovery EM streams archive batches into packed cell shards.
    #[arg(long)]
    pub eager: bool,
    /// Write a reliability diagram (per-mode empirical accuracy by responsibility decile) from
    /// the masked run's calibration counts.
    #[arg(long)]
    pub plot: Option<PathBuf>,
}

fn read_em_groups(path: &PathBuf, cells: &[u32], collapse: bool) -> Result<crate::rows::EmGroups> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("read EM group map {}", path.display()))?;
    let mut cell_id: HbMap<u32, u32> = HbMap::with_capacity(cells.len());
    for (id, &barcode) in cells.iter().enumerate() {
        cell_id.insert(barcode, id as u32);
    }
    let mut rows: Vec<(u32, String)> = Vec::new();
    let mut seen: HbMap<u32, usize> = HbMap::new();
    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 2 {
            bail!(
                "{}:{}: expected exactly two columns: BARCODE GROUP",
                path.display(),
                line_no + 1
            );
        }
        let barcode = fields[0].strip_suffix("-1").unwrap_or(fields[0]);
        if barcode.len() != 16 {
            bail!(
                "{}:{}: barcode must contain exactly 16 A/C/G/T bases",
                path.display(),
                line_no + 1
            );
        }
        let packed = umi::pack(barcode.as_bytes()).with_context(|| {
            format!(
                "{}:{}: barcode must contain exactly 16 A/C/G/T bases",
                path.display(),
                line_no + 1
            )
        })?;
        let Some(&id) = cell_id.get(&packed) else {
            bail!(
                "{}:{}: barcode {} is not in the archive cell dictionary",
                path.display(),
                line_no + 1,
                fields[0]
            );
        };
        if let Some(first_line) = seen.insert(id, line_no + 1) {
            bail!(
                "{}:{}: duplicate barcode {} (first seen on line {})",
                path.display(),
                line_no + 1,
                fields[0],
                first_line
            );
        }
        rows.push((id, fields[1].to_string()));
    }
    if rows.is_empty() {
        bail!("EM group map {} contains no data rows", path.display());
    }

    let mut eval_cell = vec![false; cells.len()];
    for (cell, _) in &rows {
        eval_cell[*cell as usize] = true;
    }
    if collapse {
        return Ok(crate::rows::EmGroups {
            cell_group: vec![0; cells.len()],
            eval_cell,
            names: vec!["all-archive-cells".to_string()],
        });
    }

    let names: Vec<String> = rows
        .iter()
        .map(|(_, name)| name.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let name_to_id: HbMap<&str, u32> = names
        .iter()
        .enumerate()
        .map(|(id, name)| (name.as_str(), id as u32))
        .collect();
    let mut cell_group = vec![crate::rows::EM_NO_GROUP; cells.len()];
    for (cell, name) in rows {
        cell_group[cell as usize] = name_to_id[&name.as_str()];
    }
    drop(name_to_id);
    Ok(crate::rows::EmGroups {
        cell_group,
        eval_cell,
        names,
    })
}

/// Exit the process now, skipping destructors: the command is finished, and serially freeing a
/// 100-million-record molecule table plus row and assignment vectors costs several seconds.
/// Both stdio streams are flushed first. Callers MUST have flushed (or dropped) every BufWriter
/// they hold before calling this — destructors do not run past this point.
fn exit_without_teardown() -> ! {
    let _ = std::io::stdout().flush();
    let _ = std::io::stderr().flush();
    std::process::exit(0);
}

/// The packed EM allocates from many Rayon workers while building bounded decode windows. glibc's
/// default arena count retains hundreds of MiB of otherwise free pages on large archives. Apply
/// the same command-local bound as `MALLOC_ARENA_MAX=2` before EM performs its first large
/// allocation; other subcommands and non-glibc targets are unaffected.
#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn configure_em_allocator() {
    const M_ARENA_MAX: libc::c_int = -8;
    unsafe extern "C" {
        fn mallopt(param: libc::c_int, value: libc::c_int) -> libc::c_int;
    }
    unsafe {
        mallopt(M_ARENA_MAX, 2);
    }
}

#[cfg(not(all(target_os = "linux", target_env = "gnu")))]
fn configure_em_allocator() {}

pub fn run_em(args: EmArgs) -> Result<()> {
    configure_em_allocator();
    if !args.mask.is_finite() || !(0.0..=1.0).contains(&args.mask) {
        bail!("--mask must be finite and between 0 and 1");
    }
    if !args.alpha.is_finite() || args.alpha < 0.0 {
        bail!("--alpha must be finite and non-negative");
    }
    if !args.group_alpha.is_finite() || args.group_alpha < 0.0 {
        bail!("--group-alpha must be finite and non-negative");
    }
    if !args.global_alpha.is_finite() || args.global_alpha < 0.0 {
        bail!("--global-alpha must be finite and non-negative");
    }
    if !args.convex_cell_weight.is_finite()
        || !(0.0..=1.0).contains(&args.convex_cell_weight)
    {
        bail!("--convex-cell-weight must be finite and between 0 and 1");
    }
    if !args.convex_group_weight.is_finite()
        || !(0.0..=1.0).contains(&args.convex_group_weight)
    {
        bail!("--convex-group-weight must be finite and between 0 and 1");
    }
    if args.convex_cell_weight + args.convex_group_weight > 1.0 {
        bail!("--convex-cell-weight plus --convex-group-weight must not exceed 1");
    }
    if !args.convex_group_prior.is_finite() || args.convex_group_prior < 0.0 {
        bail!("--convex-group-prior must be finite and non-negative");
    }
    if !args.dirichlet_cell_prior.is_finite() || args.dirichlet_cell_prior < 0.0 {
        bail!("--dirichlet-cell-prior must be finite and non-negative");
    }
    if !args.dirichlet_group_prior.is_finite() || args.dirichlet_group_prior < 0.0 {
        bail!("--dirichlet-group-prior must be finite and non-negative");
    }
    if !args.hybrid_depth_scale.is_finite() || args.hybrid_depth_scale <= 0.0 {
        bail!("--hybrid-depth-scale must be finite and greater than zero");
    }
    if !args.hybrid_depth_power.is_finite() || args.hybrid_depth_power <= 0.0 {
        bail!("--hybrid-depth-power must be finite and greater than zero");
    }
    if args.eager
        && (args.groups.is_some()
            || args.metrics_json.is_some()
            || args.candidate_genes_out.is_some())
    {
        bail!("--eager does not support hierarchical groups or machine-readable EM outputs");
    }
    if args.star
        && (args.groups.is_some()
            || args.metrics_json.is_some()
            || args.candidate_genes_out.is_some())
    {
        bail!("--star does not support hierarchical groups or machine-readable EM outputs");
    }
    let anno = anno::Annotation::from_path(&args.gtf)?;
    if args.star {
        let x = read_archive(&args.archive)?;
        let m = crate::rows::em_star_matrix(&x, &anno);
        let out_dir = args.emit.as_ref().context("--star requires --emit")?;
        let barcodes = args.barcodes.as_ref().context("--star requires --barcodes")?;
        let bc_text = std::fs::read_to_string(barcodes)?;
        let out_barcodes: Vec<&str> =
            bc_text.lines().map(|l| l.trim()).filter(|l| !l.is_empty()).collect();
        let mut bc_col: HbMap<u32, u32> = HbMap::new(); // get-only after build: order never observed
        for (i, b) in out_barcodes.iter().enumerate() {
            if let Some(p) = umi::pack(b.as_bytes()) {
                bc_col.insert(p, i as u32);
            }
        }
        std::fs::create_dir_all(out_dir)?;
        let mut feat = std::io::BufWriter::new(std::fs::File::create(out_dir.join("features.tsv"))?);
        for (id, name) in anno.gene_ids.iter().zip(&anno.gene_names) {
            writeln!(feat, "{id}\t{name}\tGene Expression")?;
        }
        std::fs::write(out_dir.join("barcodes.tsv"), out_barcodes.join("\n") + "\n")?;
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(m.len());
        for ((cell, gene), v) in m {
            let packed = x.cells[cell as usize];
            let Some(col) = bc_col.get(&packed) else { bail!("barcode not in output list") };
            triplets.push((gene, *col, v));
        }
        triplets.par_sort_unstable_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        let mut mtx =
            std::io::BufWriter::new(std::fs::File::create(out_dir.join("UniqueAndMult-EM.mtx"))?);
        writeln!(mtx, "%%MatrixMarket matrix coordinate real general")?;
        writeln!(mtx, "%")?;
        writeln!(mtx, "{} {} {}", anno.gene_ids.len(), out_barcodes.len(), triplets.len())?;
        // Parallel line formatting; the float keeps std's {:.6} formatter for identical bytes.
        let blocks: Vec<String> = triplets
            .par_chunks(1 << 20)
            .map(|ch| {
                let mut s = String::with_capacity(ch.len() * 24);
                for (g, c, v) in ch {
                    push_u32(&mut s, g + 1);
                    s.push(' ');
                    push_u32(&mut s, c + 1);
                    s.push_str(&format!(" {v:.6}\n"));
                }
                s
            })
            .collect();
        for b in &blocks {
            mtx.write_all(b.as_bytes())?;
        }
        mtx.flush()?;
        feat.flush()?;
        eprintln!("wrote UniqueAndMult-EM ({} entries) to {}", triplets.len(), out_dir.display());
        exit_without_teardown();
    }
    let mut cals: Vec<(String, [[u64; 2]; 10], f64, u64)> = Vec::new();
    let mut metrics: Vec<crate::rows::PackedModeMetrics> = Vec::new();
    let mut paired_metrics: Vec<crate::rows::PackedPairedMetrics> = Vec::new();
    let mut group_names: Vec<String> = Vec::new();
    let mut group_eval_cells = 0usize;
    let (recovered, cells) = if args.eager {
        let x = read_archive(&args.archive)?;
        let recovered = crate::rows::em_experiment(
            &x, &anno, args.mask, args.seed, args.alpha,
            args.plot.is_some().then_some(&mut cals),
        );
        (recovered, x.cells)
    } else {
        let mut archive = StreamingReplayArchive::open(&args.archive)?;
        let em = archive.packed_em_accumulator(&anno)?;
        let cells = std::mem::take(&mut archive.x.cells);
        let cell_of_class = std::mem::take(&mut archive.cell_of_class);
        drop(archive);
        let packed = em.finish(&cell_of_class, args.mask, args.seed)?;
        drop(cell_of_class);
        if let Some(path) = &args.candidate_genes_out {
            let genes = packed.candidate_gene_ids(&anno);
            let text = if genes.is_empty() {
                String::new()
            } else {
                genes.join("\n") + "\n"
            };
            std::fs::write(path, text)
                .with_context(|| format!("write candidate genes {}", path.display()))?;
            if args.candidate_genes_only {
                exit_without_teardown();
            }
        }
        let groups = args
            .groups
            .as_ref()
            .map(|path| read_em_groups(path, &cells, args.collapse_groups))
            .transpose()?;
        if let Some(groups) = &groups {
            group_names = groups.names.clone();
            group_eval_cells = groups.eval_cell.iter().filter(|v| **v).count();
            eprintln!(
                "hierarchical EM: {} groups, {} scored cells{}",
                groups.names.len(),
                group_eval_cells,
                if args.collapse_groups { " (collapse control)" } else { "" }
            );
        } else {
            group_eval_cells = cells.len();
        }
        let recovered = packed.run(
            &anno,
            args.alpha,
            groups.as_ref(),
            args.group_alpha,
            args.global_alpha,
            crate::rows::ConvexEmParams {
                cell_weight: args.convex_cell_weight,
                group_weight: args.convex_group_weight,
                group_prior: args.convex_group_prior,
            },
            args.convex_only,
            crate::rows::DirichletProxyParams {
                cell_prior: args.dirichlet_cell_prior,
                group_prior: args.dirichlet_group_prior,
            },
            args.dirichlet_only,
            crate::rows::DepthHybridParams {
                depth_scale: args.hybrid_depth_scale,
                depth_power: args.hybrid_depth_power,
            },
            args.hybrid_only,
            args.plot.is_some().then_some(&mut cals),
            args.metrics_json.is_some().then_some(&mut metrics),
            args.metrics_json.is_some().then_some(&mut paired_metrics),
        );
        (recovered, cells)
    };
    if let Some(path) = &args.metrics_json {
        let modes: Vec<serde_json::Value> = metrics
            .iter()
            .map(|m| {
                let denom = m.n.max(1) as f64;
                let evidence_depth_strata: Vec<serde_json::Value> = m
                    .evidence_depth
                    .iter()
                    .map(|depth| {
                        let depth_denom = depth.n.max(1) as f64;
                        serde_json::json!({
                            "name": depth.name,
                            "n": depth.n,
                            "top1_percent": 100.0 * depth.top1 as f64 / depth_denom,
                            "expected_accuracy_percent": 100.0 * depth.expected / depth_denom,
                            "negative_log_loss": depth.negative_log_loss / depth_denom,
                            "multiclass_brier": depth.brier / depth_denom,
                        })
                    })
                    .collect();
                serde_json::json!({
                    "name": m.name,
                    "n": m.n,
                    "top1_percent": 100.0 * m.top1 as f64 / denom,
                    "expected_accuracy_percent": 100.0 * m.expected / denom,
                    "negative_log_loss": m.negative_log_loss / denom,
                    "multiclass_brier": m.brier / denom,
                    "calibration_max_probability_deciles": m.calibration,
                    "evidence_depth_strata": evidence_depth_strata,
                })
            })
            .collect();
        let paired_comparisons: Vec<serde_json::Value> = paired_metrics
            .iter()
            .map(|pair| {
                serde_json::json!({
                    "candidate": pair.candidate,
                    "reference": pair.reference,
                    "n": pair.n,
                    "cells": pair.cells,
                    "mean_negative_log_loss_difference": pair.mean_negative_log_loss_difference,
                    "negative_log_loss_clustered_se": pair.negative_log_loss_clustered_se,
                    "negative_log_loss_ci95": pair.negative_log_loss_ci95,
                    "negative_log_loss_wins": pair.negative_log_loss_wins,
                    "mean_brier_difference": pair.mean_brier_difference,
                    "brier_clustered_se": pair.brier_clustered_se,
                    "brier_ci95": pair.brier_ci95,
                    "brier_wins": pair.brier_wins,
                    "top1_percentage_point_difference": pair.top1_percentage_point_difference,
                })
            })
            .collect();
        let document = serde_json::json!({
            "schema_version": 2,
            "archive": args.archive,
            "gtf": args.gtf,
            "mask_fraction": args.mask,
            "seed": args.seed,
            "blend_alpha": args.alpha,
            "group_alpha": args.group_alpha,
            "global_alpha": args.global_alpha,
            "convex_cell_weight": args.convex_cell_weight,
            "convex_group_weight": args.convex_group_weight,
            "convex_global_weight": 1.0 - args.convex_cell_weight - args.convex_group_weight,
            "convex_group_prior": args.convex_group_prior,
            "convex_only": args.convex_only,
            "convex_candidate_normalized": true,
            "convex_leave_one_cell_out_group": true,
            "dirichlet_proxy_cell_prior": args.dirichlet_cell_prior,
            "dirichlet_proxy_group_prior": args.dirichlet_group_prior,
            "dirichlet_proxy_only": args.dirichlet_only,
            "dirichlet_proxy_candidate_normalized": true,
            "dirichlet_proxy_leave_one_cell_out_group": true,
            "dirichlet_proxy_inference": "posterior_mean_screen_not_full_variational",
            "depth_hybrid_scale": args.hybrid_depth_scale,
            "depth_hybrid_power": args.hybrid_depth_power,
            "depth_hybrid_only": args.hybrid_only,
            "depth_hybrid_gate": "d^power/(d^power+scale^power) using mode-independent fitted unique candidate mass",
            "depth_hybrid_fixed_component": "candidate-normalized convex",
            "depth_hybrid_fixed_cell_weight": 0.2,
            "depth_hybrid_fixed_group_weight": 0.6,
            "depth_hybrid_fixed_sample_weight": 0.2,
            "depth_hybrid_fixed_group_prior": 80.0,
            "depth_hybrid_adaptive_component": "posterior-mean Dirichlet proxy",
            "depth_hybrid_adaptive_cell_prior": 64.0,
            "depth_hybrid_adaptive_group_prior": 80.0,
            "evidence_depth_definition": "sum of mode-independent fitted unique counts over the target candidate genes: 0-1, (1,4], (4,16], 16+",
            "groups_file": args.groups,
            "collapse_groups": args.collapse_groups,
            "group_names": group_names,
            "scored_cells": group_eval_cells,
            "modes": modes,
            "paired_comparisons": paired_comparisons,
        });
        std::fs::write(path, serde_json::to_vec_pretty(&document)?)
            .with_context(|| format!("write EM metrics {}", path.display()))?;
    }
    if let Some(plot_path) = &args.plot {
        if cals.is_empty() {
            bail!("--plot needs a masked run (--mask > 0) to have calibration counts");
        }
        crate::plots::em_plot(plot_path, &cals, "EM masked recovery: reliability")?;
        eprintln!("wrote {}", plot_path.display());
    }
    if let Some(out_dir) = &args.emit {
        let Some(rec) = recovered else {
            bail!("--emit requires --mask 0 (the impact run produces the layer)");
        };
        let barcodes = args.barcodes.as_ref().context("--emit requires --barcodes")?;
        let bc_text = std::fs::read_to_string(barcodes)?;
        let out_barcodes: Vec<&str> =
            bc_text.lines().map(|l| l.trim()).filter(|l| !l.is_empty()).collect();
        let mut bc_col: HbMap<u32, u32> = HbMap::new(); // get-only after build: order never observed
        for (i, b) in out_barcodes.iter().enumerate() {
            if let Some(p) = umi::pack(b.as_bytes()) {
                bc_col.insert(p, i as u32);
            }
        }
        std::fs::create_dir_all(out_dir)?;
        let mut feat = std::io::BufWriter::new(std::fs::File::create(out_dir.join("features.tsv"))?);
        for (id, name) in anno.gene_ids.iter().zip(&anno.gene_names) {
            writeln!(feat, "{id}\t{name}\tGene Expression")?;
        }
        std::fs::write(out_dir.join("barcodes.tsv"), out_barcodes.join("\n") + "\n")?;
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(rec.len());
        for ((cell, gene), v) in rec {
            let packed = cells[cell as usize];
            let Some(col) = bc_col.get(&packed) else { bail!("barcode not in output list") };
            triplets.push((gene, *col, v));
        }
        triplets.par_sort_unstable_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        let mut mtx = std::io::BufWriter::new(std::fs::File::create(out_dir.join("em.mtx"))?);
        writeln!(mtx, "%%MatrixMarket matrix coordinate real general")?;
        writeln!(mtx, "% additive EM-recovered multimapper layer; exact-replay matrices unchanged")?;
        writeln!(mtx, "{} {} {}", anno.gene_ids.len(), out_barcodes.len(), triplets.len())?;
        // Parallel line formatting; the float keeps std's {:.4} formatter for identical bytes.
        let blocks: Vec<String> = triplets
            .par_chunks(1 << 20)
            .map(|ch| {
                let mut s = String::with_capacity(ch.len() * 20);
                for (g, c, v) in ch {
                    push_u32(&mut s, g + 1);
                    s.push(' ');
                    push_u32(&mut s, c + 1);
                    s.push_str(&format!(" {v:.4}\n"));
                }
                s
            })
            .collect();
        for b in &blocks {
            mtx.write_all(b.as_bytes())?;
        }
        mtx.flush()?;
        feat.flush()?;
        eprintln!("wrote {} recovered entries to {}", triplets.len(), out_dir.display());
    }
    exit_without_teardown();
}

struct DigestingWriter<W> {
    inner: W,
    hasher: blake3::Hasher,
    bytes: u64,
}

impl<W> DigestingWriter<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            hasher: blake3::Hasher::new(),
            bytes: 0,
        }
    }
}

impl<W: Write> Write for DigestingWriter<W> {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buffer)?;
        self.hasher.update(&buffer[..written]);
        self.bytes = self.bytes.saturating_add(written as u64);
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

fn write_file_with_identity<F>(path: &Path, name: &str, produce: F) -> Result<ArtifactFileIdentity>
where
    F: FnOnce(&mut DigestingWriter<std::io::BufWriter<std::fs::File>>) -> Result<()>,
{
    let file = std::fs::File::create(path)
        .with_context(|| format!("creating replay artifact {}", path.display()))?;
    let mut writer = DigestingWriter::new(std::io::BufWriter::new(file));
    produce(&mut writer)?;
    writer.flush()?;
    Ok(ArtifactFileIdentity {
        name: name.to_owned(),
        bytes: writer.bytes,
        blake3: writer.hasher.finalize().to_hex().to_string(),
    })
}

fn identity_for_bytes(name: &str, bytes: &[u8]) -> ArtifactFileIdentity {
    ArtifactFileIdentity {
        name: name.to_owned(),
        bytes: bytes.len() as u64,
        blake3: blake3::hash(bytes).to_hex().to_string(),
    }
}

fn write_features<W: Write>(writer: &mut W, anno: &anno::Annotation) -> Result<()> {
    for (id, name) in anno.gene_ids.iter().zip(&anno.gene_names) {
        writeln!(writer, "{id}\t{name}\tGene Expression")?;
    }
    Ok(())
}

fn write_velocity_matrix<W: Write>(
    writer: &mut W,
    anno: &anno::Annotation,
    barcode_count: usize,
    triplets: &[(u32, u32, [u32; 3], bool)],
    component: usize,
) -> Result<()> {
    writeln!(writer, "%%MatrixMarket matrix coordinate integer general")?;
    writeln!(writer, "%")?;
    writeln!(
        writer,
        "{} {} {}",
        anno.gene_ids.len(),
        barcode_count,
        triplets.len()
    )?;
    for (gene, cell, values, _) in triplets {
        writeln!(writer, "{} {} {}", gene + 1, cell + 1, values[component])?;
    }
    Ok(())
}

fn write_velocity_completeness<W: Write>(
    writer: &mut W,
    triplets: &[(u32, u32, [u32; 3], bool)],
) -> Result<()> {
    for (gene, cell, _, complete) in triplets {
        writeln!(writer, "{} {} {}", gene + 1, cell + 1, *complete as u8)?;
    }
    Ok(())
}

fn write_gene_matrix<W: Write>(
    writer: &mut W,
    anno: &anno::Annotation,
    barcode_count: usize,
    triplet_count: usize,
    blocks: &[String],
) -> Result<()> {
    writeln!(writer, "%%MatrixMarket matrix coordinate integer general")?;
    writeln!(writer, "%")?;
    writeln!(
        writer,
        "{} {} {}",
        anno.gene_ids.len(),
        barcode_count,
        triplet_count
    )?;
    for block in blocks {
        writer.write_all(block.as_bytes())?;
    }
    Ok(())
}

/// STARsolo Velocyto layout: features/barcodes plus spliced/unspliced/ambiguous.mtx sharing one
/// (gene, cell) entry set. Uniform reporting hashes bytes during the existing write pass.
fn emit_velocity(
    counts: &FxHashMap<(u32, u32), ([u32; 3], bool)>,
    x: &Extracted,
    anno: &anno::Annotation,
    barcodes: &PathBuf,
    barcodes_text: Option<&str>,
    out_dir: &PathBuf,
    report: bool,
) -> Result<Option<ReplayArtifactManifest>> {
    let owned_barcodes;
    let bc_text = if let Some(text) = barcodes_text {
        text
    } else {
        owned_barcodes = std::fs::read_to_string(barcodes)?;
        &owned_barcodes
    };
    let out_barcodes: Vec<&str> = bc_text
        .lines()
        .map(|line| line.trim())
        .filter(|line| !line.is_empty())
        .collect();
    let mut bc_col: HbMap<u32, u32> = HbMap::new();
    for (index, barcode) in out_barcodes.iter().enumerate() {
        if let Some(packed) = umi::pack(barcode.as_bytes()) {
            bc_col.insert(packed, index as u32);
        }
    }
    std::fs::create_dir_all(out_dir)?;
    let mut files = Vec::new();
    if report {
        files.push(write_file_with_identity(
            &out_dir.join("features.tsv"),
            "features.tsv",
            |writer| write_features(writer, anno),
        )?);
    } else {
        let mut writer =
            std::io::BufWriter::new(std::fs::File::create(out_dir.join("features.tsv"))?);
        write_features(&mut writer, anno)?;
        writer.flush()?;
    }
    let barcode_bytes = out_barcodes.join("\n") + "\n";
    std::fs::write(out_dir.join("barcodes.tsv"), barcode_bytes.as_bytes())?;
    if report {
        files.push(identity_for_bytes("barcodes.tsv", barcode_bytes.as_bytes()));
    }

    let mut triplets: Vec<(u32, u32, [u32; 3], bool)> = Vec::with_capacity(counts.len());
    for ((cell, gene), (counts, complete)) in counts {
        let packed = x.cells[*cell as usize];
        let Some(column) = bc_col.get(&packed) else {
            bail!("barcode not in output list")
        };
        triplets.push((*gene, *column, *counts, *complete));
    }
    triplets.par_sort_unstable();
    let matrix_names = ["spliced.mtx", "unspliced.mtx", "ambiguous.mtx"];
    for (component, name) in matrix_names.iter().enumerate() {
        if report {
            files.push(write_file_with_identity(
                &out_dir.join(name),
                name,
                |writer| {
                    write_velocity_matrix(writer, anno, out_barcodes.len(), &triplets, component)
                },
            )?);
        } else {
            let mut writer = std::io::BufWriter::new(std::fs::File::create(out_dir.join(name))?);
            write_velocity_matrix(&mut writer, anno, out_barcodes.len(), &triplets, component)?;
        }
    }
    if report {
        files.push(write_file_with_identity(
            &out_dir.join("entries.complete.tsv"),
            "entries.complete.tsv",
            |writer| write_velocity_completeness(writer, &triplets),
        )?);
    } else {
        let mut writer =
            std::io::BufWriter::new(std::fs::File::create(out_dir.join("entries.complete.tsv"))?);
        write_velocity_completeness(&mut writer, &triplets)?;
    }
    Ok(report.then(|| {
        let matrices: Vec<MexManifest> = matrix_names
            .iter()
            .map(|name| MexManifest {
                format: "matrix_market_coordinate".into(),
                matrix: (*name).to_owned(),
                features: "features.tsv".into(),
                barcodes: "barcodes.tsv".into(),
                feature_count: anno.gene_ids.len(),
                barcode_count: out_barcodes.len(),
                nonzero_count: triplets.len(),
                index_base: 1,
                value_type: "integer".into(),
            })
            .collect();
        ReplayArtifactManifest {
            // Keep the ordinary MEX manifest at metadata.data's top level so existing
            // readers can open the spliced component directly. The complete component list is
            // additive metadata for velocity-aware clients.
            primary: matrices[0].clone(),
            transactional_publication: false,
            completion_marker: "metadata.json",
            matrices,
            files,
        }
    }))
}

fn emit_matrix(
    counts: &FxHashMap<(u32, u32), u32>,
    x: &Extracted,
    anno: &anno::Annotation,
    barcodes: &PathBuf,
    barcodes_text: Option<&str>,
    out_dir: &PathBuf,
    report: bool,
) -> Result<Option<ReplayArtifactManifest>> {
    let owned_barcodes;
    let bc_text = if let Some(text) = barcodes_text {
        text
    } else {
        owned_barcodes = std::fs::read_to_string(barcodes)?;
        &owned_barcodes
    };
    let out_barcodes: Vec<&str> = bc_text
        .lines()
        .map(|line| line.trim())
        .filter(|line| !line.is_empty())
        .collect();
    let mut bc_col: HbMap<u32, u32> = HbMap::new();
    for (index, barcode) in out_barcodes.iter().enumerate() {
        if let Some(packed) = umi::pack(barcode.as_bytes()) {
            bc_col.insert(packed, index as u32);
        }
    }
    std::fs::create_dir_all(out_dir)?;
    let mut files = Vec::new();
    if report {
        files.push(write_file_with_identity(
            &out_dir.join("features.tsv"),
            "features.tsv",
            |writer| write_features(writer, anno),
        )?);
    } else {
        let mut writer =
            std::io::BufWriter::new(std::fs::File::create(out_dir.join("features.tsv"))?);
        write_features(&mut writer, anno)?;
        writer.flush()?;
    }
    let barcode_bytes = out_barcodes.join("\n") + "\n";
    std::fs::write(out_dir.join("barcodes.tsv"), barcode_bytes.as_bytes())?;
    if report {
        files.push(identity_for_bytes("barcodes.tsv", barcode_bytes.as_bytes()));
    }

    let mut triplets: Vec<(u32, u32, u32)> = Vec::with_capacity(counts.len());
    for ((cell, gene), count) in counts {
        let packed = x.cells[*cell as usize];
        let Some(column) = bc_col.get(&packed) else {
            bail!("barcode not in output list")
        };
        triplets.push((*gene, *column, *count));
    }
    triplets.par_sort_unstable();
    let blocks: Vec<String> = triplets
        .par_chunks(1 << 20)
        .map(|chunk| {
            let mut text = String::with_capacity(chunk.len() * 16);
            for (gene, cell, count) in chunk {
                push_u32(&mut text, gene + 1);
                text.push(' ');
                push_u32(&mut text, cell + 1);
                text.push(' ');
                push_u32(&mut text, *count);
                text.push('\n');
            }
            text
        })
        .collect();
    if report {
        files.push(write_file_with_identity(
            &out_dir.join("matrix.mtx"),
            "matrix.mtx",
            |writer| write_gene_matrix(writer, anno, out_barcodes.len(), triplets.len(), &blocks),
        )?);
    } else {
        let mut writer =
            std::io::BufWriter::new(std::fs::File::create(out_dir.join("matrix.mtx"))?);
        write_gene_matrix(
            &mut writer,
            anno,
            out_barcodes.len(),
            triplets.len(),
            &blocks,
        )?;
        writer.flush()?;
    }
    Ok(report.then(|| {
        let primary = MexManifest {
            format: "matrix_market_coordinate".into(),
            matrix: "matrix.mtx".into(),
            features: "features.tsv".into(),
            barcodes: "barcodes.tsv".into(),
            feature_count: anno.gene_ids.len(),
            barcode_count: out_barcodes.len(),
            nonzero_count: triplets.len(),
            index_base: 1,
            value_type: "integer".into(),
        };
        ReplayArtifactManifest {
            // This flattened manifest is the established gravlax-output MEX metadata shape.
            primary: primary.clone(),
            transactional_publication: false,
            completion_marker: "metadata.json",
            matrices: vec![primary],
            files,
        }
    }))
}

/// Append a u32's decimal digits — byte-identical to `format!("{v}")`, without the formatter.
#[inline]
fn push_u32(s: &mut String, mut v: u32) {
    let mut buf = [0u8; 10];
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    s.push_str(std::str::from_utf8(&buf[i..]).unwrap());
}

pub fn run_ingest(args: IngestArgs) -> Result<()> {
    validate_alignment_args(&args)?;
    let chunk_bp = args
        .chunk_mb
        .checked_mul(1_000_000)
        .context("--chunk-mb exceeds the supported coordinate range")?;
    let reporting = preflight_artifact_report(&args.report, &[&args.out])?;
    if reporting {
        path_parameter(&args.bam)?;
        path_parameter(&args.whitelist)?;
        if let Some(genome) = &args.genome {
            path_parameter(genome)?;
        }
    }
    let t0 = std::time::Instant::now();
    // Every newly ingested archive is logical molecular-evidence v2. Snapshot the whitelist and
    // derive the BAM identity from the same held file that supplies both extraction passes so the
    // root-bound manifest never describes a path reopened after the evidence was produced.
    let whitelist_snapshot = stable_utf8_file(&args.whitelist)?;
    let genome_input_thread = args.genome.as_deref().map(start_genome_input).transpose()?;
    // Genome hashing overlaps BAM extraction; joined (and validated) before the meta is written.
    let archive_extraction = extract_rows_for_archive(
        &args.bam,
        &whitelist_snapshot.0,
        args.locus_gap,
        args.terminal_tails,
    )?;
    let x = archive_extraction.evidence;
    let terminal_tails = archive_extraction.terminal_tails;
    let bam_programs = archive_extraction.bam_programs;
    let bam_identity = FileContentIdentity {
        scheme: "full-file-blake3-v1",
        blake3: archive_extraction.bam_identity.blake3,
        bytes: archive_extraction.bam_identity.bytes,
    };
    eprintln!(
        "extracted: {} molecules, {} edges, {} cells, {} shapes, {} patterns, {} classes ({:.1}s)",
        x.mols.len(), x.edges.len(), x.cells.len(), x.shapes.len(), x.patterns.len(), x.n_classes,
        t0.elapsed().as_secs_f32()
    );
    if !reporting {
        pattern_stats(&x);
    }
    let mut genome_file_identity = None;
    let genome_sig = match genome_input_thread {
        Some(thread) => {
            let (sig, identity) = thread
                .join()
                .map_err(|_| anyhow::anyhow!("genome input thread panicked"))??;
            genome_file_identity = Some(identity);
            Some(sig)
        }
        None => None,
    };
    let genome_sig = match genome_sig {
        Some(sig) => {
            let missing: Vec<&String> = x
                .chrom_names
                .iter()
                .filter(|n| sig.contig(n).is_none())
                .collect();
            if !missing.is_empty() {
                bail!(
                    "--genome does not look like the alignment reference: {} BAM contig(s) absent \
                     from the FASTA (first: {})",
                    missing.len(),
                    missing[0]
                );
            }
            eprintln!("genome signature: {} contigs, digest {}", sig.contigs.len(), &sig.digest[..16]);
            Some(sig)
        }
        None => None,
    };

    let catalogue_role = match args.junction_discovery {
        JunctionDiscoveryArg::PerLibraryTwoPass => Some(JunctionCatalogueRole::PerLibraryPass1),
        JunctionDiscoveryArg::FrozenCatalogue => Some(JunctionCatalogueRole::FrozenExternal),
        JunctionDiscoveryArg::OnePass | JunctionDiscoveryArg::Unspecified => None,
    };
    let (catalogue, junction_catalogue_bytes) = match (&args.junction_catalogue, catalogue_role) {
        (Some(path), Some(role)) => {
            let (catalogue, bytes) = junction_catalogue(path, role)?;
            (Some(catalogue), Some(bytes))
        }
        (None, None) => (None, None),
        _ => unreachable!("alignment argument validation enforces catalogue mode"),
    };
    let alignment_annotation = args
        .alignment_annotation
        .as_deref()
        .map(declared_alignment_file)
        .transpose()?;
    let ordered_inputs = args
        .alignment_inputs
        .iter()
        .map(|path| declared_alignment_file(path))
        .collect::<Result<Vec<_>>>()?;
    let alignment_log = args
        .alignment_log
        .as_deref()
        .map(declared_alignment_file)
        .transpose()?;
    let junction_discovery = args.junction_discovery.into();
    let genome_reference_binding = match (&genome_file_identity, &genome_sig) {
        (Some(identity), Some(signature)) => Some(GenomeReferenceBinding::new(
            GenomeBindingAction::IngestArchive,
            verified_identity(identity),
            signature.clone(),
        )),
        (None, None) => None,
        _ => unreachable!("genome signature and exact identity are produced together"),
    };
    let manifest = AlignmentProvenanceManifest {
        schema: ALIGNMENT_PROVENANCE_SCHEMA.into(),
        molecular_evidence_schema: MOLECULAR_EVIDENCE_SCHEMA.into(),
        alignment: AlignmentDeclaration {
            status: if junction_discovery == JunctionDiscoveryMode::Unspecified {
                ProvenanceStatus::Unspecified
            } else {
                ProvenanceStatus::DeclaredByCaller
            },
            junction_discovery,
            programs: bam_programs,
            junction_catalogue: catalogue,
            alignment_annotation,
            ordered_inputs,
            alignment_log,
            chemistry: args.alignment_chemistry.clone(),
            chemistry_status: if args.alignment_chemistry.is_some() {
                ProvenanceStatus::DeclaredByCaller
            } else {
                ProvenanceStatus::Unspecified
            },
            index_identity: args.alignment_index_identity.clone(),
            index_identity_status: if args.alignment_index_identity.is_some() {
                ProvenanceStatus::DeclaredByCaller
            } else {
                ProvenanceStatus::Unspecified
            },
        },
        inputs: AlignmentInputs {
            bam: verified_identity(&bam_identity),
            whitelist: verified_identity(&whitelist_snapshot.1),
            genome_fasta: genome_file_identity.as_ref().map(verified_identity),
            genome_signature: genome_sig.clone(),
            genome_relationship_status: if genome_sig.is_some() {
                ProvenanceStatus::DeclaredByCaller
            } else {
                ProvenanceStatus::Unspecified
            },
        },
        ingest: IngestProvenance {
            program: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
            locus_gap: args.locus_gap,
            chunk_bp,
            zstd_level: args.zstd_level,
            molecule_chunk_streams: 10,
            molecule_codec: "rans2".into(),
            barcode_correction: "unique-hamming1-quality-pseudocount-v1".into(),
            umi_classes: "global-cell-umi-equivalence-with-1mm-edges-v1".into(),
            unique_chain_reduction: "junction-chain-span-extremes-v1".into(),
            multimapper_reduction: "primary-relative-placement-pattern-v1".into(),
            terminal_tail_rule: args
                .terminal_tails
                .then(|| terminal_tail::TERMINAL_TAIL_RULE.into()),
        },
    };
    manifest.validate()?;
    let terminal_tail_summary = if reporting {
        terminal_tails
            .as_ref()
            .map(|tails| {
                prepare_terminal_tails(&x, tails, chunk_bp).map(|prepared| prepared.metadata)
            })
            .transpose()?
    } else {
        None
    };
    let (acct, output_identity) = if reporting {
        let (temporary, writer) = temporary_ingest_writer(&args.out, args.zstd_level)?;
        let staging_guard = writer.try_clone_file()?;
        let result = (|| -> Result<_> {
            let writer = write_archive_sections(
                &x,
                writer,
                args.zstd_level,
                chunk_bp,
                genome_sig.as_ref(),
                ArchiveExtensions {
                    alignment_provenance: Some(&manifest),
                    junction_catalogue_bytes: junction_catalogue_bytes.as_deref(),
                    terminal_tails: terminal_tails.as_ref(),
                    genome_reference_binding: genome_reference_binding.as_ref(),
                },
            )?;
            let (accounting, output_file, commitment) = writer.finish_with_file()?;
            let output_reader = SectionReader::from_file(output_file.try_clone()?)?;
            if output_reader.content_commitment() != Some(commitment) {
                bail!("finished ingest archive differs from its computed root commitment");
            }
            let identity = archive_identity(&output_reader)?;
            let outcome = install_open_file_no_clobber(
                &output_file,
                &temporary,
                &args.out,
                Durability::FileAndDirectory,
            )?;
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
            Ok((accounting, identity))
        })();
        if result.is_err() {
            warn_on_failed_staging_cleanup(&temporary, &staging_guard);
        }
        let (accounting, identity) = result?;
        (accounting, Some(identity))
    } else {
        let writer = write_archive_sections(
            &x,
            SectionWriter::create(&args.out, args.zstd_level)?,
            args.zstd_level,
            chunk_bp,
            genome_sig.as_ref(),
            ArchiveExtensions {
                alignment_provenance: Some(&manifest),
                junction_catalogue_bytes: junction_catalogue_bytes.as_deref(),
                terminal_tails: terminal_tails.as_ref(),
                genome_reference_binding: genome_reference_binding.as_ref(),
            },
        )?;
        (writer.finish()?, None)
    };
    let total: u64 = acct.iter().map(|(_, _, c)| *c).sum();
    if !reporting {
        println!("{:<14} {:>12} {:>12}", "section", "raw", "zstd");
        for (n, r, c) in &acct {
            println!("{n:<14} {r:>12} {c:>12}");
        }
        println!(
            "{:<14} {:>12} {:>12}  ({:.2} bits/molecule over {} molecules)",
            "TOTAL",
            "",
            total,
            8.0 * total as f64 / x.mols.len().max(1) as f64,
            x.mols.len()
        );
        println!(
            "wrote {} ({} bytes)",
            args.out.display(),
            std::fs::metadata(&args.out)?.len()
        );
        return Ok(());
    }

    let output_identity = output_identity.expect("reporting ingest captured held output identity");
    let mut inputs = BTreeMap::new();
    inputs.insert("bam", bam_identity);
    inputs.insert("whitelist", whitelist_snapshot.1);
    if let Some(identity) = genome_file_identity {
        inputs.insert("genome", identity);
    }
    let raw_section_bytes = acct.iter().try_fold(0u64, |sum, (_, raw, _)| {
        sum.checked_add(*raw)
            .context("raw section byte count overflow")
    })?;
    #[derive(Serialize)]
    struct IngestSummary<'a> {
        inputs: &'a BTreeMap<&'static str, FileContentIdentity>,
        output_archive: &'a ArchiveContentIdentity,
        molecules: u64,
        edges: u64,
        cells: u64,
        shapes: u64,
        patterns: u64,
        umi_classes: u64,
        sections: u64,
        raw_section_bytes: u64,
        compressed_payload_bytes: u64,
        bits_per_molecule: f64,
        genome_signature: Option<&'a evidence_io::genome::GenomeSig>,
        molecular_evidence_schema: &'static str,
        alignment_provenance: &'a AlignmentProvenanceManifest,
        terminal_tail_status: &'static str,
        terminal_tail: Option<&'a TerminalTailMetadata>,
        genome_reference_binding: Option<&'a GenomeReferenceBinding>,
    }
    let summary = IngestSummary {
        inputs: &inputs,
        output_archive: &output_identity,
        molecules: x.mols.len() as u64,
        edges: x.edges.len() as u64,
        cells: x.cells.len() as u64,
        shapes: x.shapes.len() as u64,
        patterns: x.patterns.len() as u64,
        umi_classes: x.n_classes as u64,
        sections: acct.len() as u64,
        raw_section_bytes,
        compressed_payload_bytes: total,
        bits_per_molecule: 8.0 * total as f64 / x.mols.len().max(1) as f64,
        genome_signature: genome_sig.as_ref(),
        molecular_evidence_schema: MOLECULAR_EVIDENCE_SCHEMA,
        alignment_provenance: &manifest,
        terminal_tail_status: if terminal_tail_summary.is_some() {
            "available"
        } else {
            "unavailable"
        },
        terminal_tail: terminal_tail_summary.as_ref(),
        genome_reference_binding: genome_reference_binding.as_ref(),
    };
    let mut parameters = BTreeMap::new();
    parameters.insert("bam".into(), path_parameter(&args.bam)?);
    parameters.insert("whitelist".into(), path_parameter(&args.whitelist)?);
    parameters.insert("output".into(), output_path_parameter(&args.out)?);
    parameters.insert("locus_gap".into(), serde_json::json!(args.locus_gap));
    parameters.insert("zstd_level".into(), serde_json::json!(args.zstd_level));
    parameters.insert("chunk_mb".into(), serde_json::json!(args.chunk_mb));
    parameters.insert(
        "terminal_tails".into(),
        serde_json::json!(args.terminal_tails),
    );
    parameters.insert(
        "junction_discovery".into(),
        serde_json::json!(junction_discovery_name(args.junction_discovery)),
    );
    if let Some(genome) = &args.genome {
        parameters.insert("genome".into(), path_parameter(genome)?);
    }
    if let Some(catalogue) = &args.junction_catalogue {
        parameters.insert("junction_catalogue".into(), path_parameter(catalogue)?);
    }
    if let Some(annotation) = &args.alignment_annotation {
        parameters.insert("alignment_annotation".into(), path_parameter(annotation)?);
    }
    if let Some(identity) = &args.alignment_index_identity {
        parameters.insert(
            "alignment_index_identity".into(),
            serde_json::json!(identity),
        );
    }
    if !args.alignment_inputs.is_empty() {
        let inputs = args
            .alignment_inputs
            .iter()
            .map(|path| path_parameter(path))
            .collect::<Result<Vec<_>>>()?;
        parameters.insert("alignment_inputs".into(), serde_json::Value::Array(inputs));
    }
    if let Some(log) = &args.alignment_log {
        parameters.insert("alignment_log".into(), path_parameter(log)?);
    }
    if let Some(chemistry) = &args.alignment_chemistry {
        parameters.insert("alignment_chemistry".into(), serde_json::json!(chemistry));
    }
    let context = report_context(
        vec![output_identity.provenance_identity()],
        parameters,
        Vec::new(),
    );
    let schema = section_table_schema()?;
    let format = args
        .report
        .report_format
        .expect("report format preflighted");
    write_artifact_report(&args.report, |writer| {
        let mut bundle = StreamingBundleWriter::new_with_summary(
            writer,
            INGEST_REPORT_SCHEMA,
            format.into(),
            &context,
            &summary,
        )?;
        bundle.write_table("sections", &schema, None, |table| {
            for (name, raw, compressed) in &acct {
                table.write_row_with(|row| {
                    row.string(name)?;
                    row.uint64(*raw)?;
                    row.uint64(*compressed)
                })?;
            }
            Ok(())
        })?;
        bundle.finish()?;
        Ok(())
    })?;
    Ok(())
}

pub fn run_replay_rows(args: ReplayRowsArgs) -> Result<()> {
    let reporting = preflight_artifact_report(&args.report, &[&args.out_dir])?;
    if reporting {
        path_parameter(&args.input)?;
        path_parameter(&args.gtf)?;
        path_parameter(&args.barcodes)?;
        if let Some(whitelist) = &args.whitelist {
            path_parameter(whitelist)?;
        }
    }
    let metadata_path = if reporting && !args.audit_multigene {
        Some(preflight_replay_metadata(
            &args.out_dir,
            args.report.report_output.as_deref(),
            args.velocity,
        )?)
    } else {
        None
    };
    let barcode_snapshot = if reporting && !args.audit_multigene {
        Some(stable_utf8_file(&args.barcodes)?)
    } else {
        None
    };
    let whitelist_snapshot = if reporting && args.from_bam {
        let whitelist = args
            .whitelist
            .as_ref()
            .context("--whitelist required with --from-bam")?;
        Some(stable_utf8_file(whitelist)?)
    } else {
        None
    };
    let t0 = std::time::Instant::now();
    let solo_strand = args.solo_strand.into();
    let streaming_gene = !args.from_bam && !args.from_molecule_bam && !args.eager
        && !args.velocity && !args.audit_multigene;
    if streaming_gene {
        let archive = StreamingReplayArchive::open(&args.input)?;
        let t_open = t0.elapsed().as_secs_f32();
        let (anno, annotation_identity) = load_replay_annotation(&args.gtf, reporting)?;
        let t_anno = t0.elapsed().as_secs_f32();
        let (counts, n_assigned, total) = archive.replay(&anno, solo_strand)?;
        let t_replay = t0.elapsed().as_secs_f32();
        let artifact = emit_matrix(
            &counts,
            &archive.x,
            &anno,
            &args.barcodes,
            barcode_snapshot.as_ref().map(|(text, _)| text.as_str()),
            &args.out_dir,
            reporting,
        )?;
        if reporting {
            let source_identity = archive_identity(&archive.reader)?;
            let annotation_identity = annotation_identity
                .as_ref()
                .expect("reporting replay bound annotation identity");
            let warnings = vec![
                "MEX component files retain the legacy in-place directory behavior; metadata.json is written last as a completion marker, but the directory is not atomically published."
                    .into(),
            ];
            let context = replay_report_context(
                &args,
                Some(&source_identity),
                annotation_identity,
                warnings,
            )?;
            let summary = ReplayReportSummary {
                source_kind: "archive",
                input_identity: serde_json::to_value(&source_identity)?,
                annotation_identity: annotation_identity.clone(),
                barcode_identity: barcode_snapshot
                    .as_ref()
                    .map(|(_, identity)| identity.clone()),
                whitelist_identity: None,
                strand_policy: replay_strand_name(args.solo_strand),
                velocity: false,
                multigene_audit: false,
                molecules: archive.n_mols as u64,
                assigned_molecules: Some(n_assigned),
                counted_umis: Some(total),
                matrix_entries: Some(counts.len() as u64),
                audit: None,
                artifact,
            };
            emit_replay_report(&args, &summary, &context, metadata_path.as_deref())?;
        }
        eprintln!(
            "molecules {} -> assigned {} -> {} UMIs in {} entries | streaming | open+dict {t_open:.1}s, +anno {:.1}s, +decode/replay {:.1}s, total {:.1}s",
            archive.n_mols, n_assigned, total, counts.len(), t_anno - t_open,
            t_replay - t_anno, t0.elapsed().as_secs_f32()
        );
        exit_without_teardown();
    }

    let mut decoded_archive_identity = None;
    let mut raw_input_identity = None;
    let x = if args.from_bam {
        let wl = args
            .whitelist
            .as_ref()
            .context("--whitelist required with --from-bam")?;
        if let Some((whitelist_text, _)) = &whitelist_snapshot {
            let (extracted, identity) =
                extract_rows_with_identity(&args.input, whitelist_text, args.locus_gap)?;
            raw_input_identity = Some(report_file_identity(identity));
            extracted
        } else {
            extract_rows(&args.input, wl, args.locus_gap)?
        }
    } else if args.from_molecule_bam {
        if args.whitelist.is_some() {
            bail!("--whitelist is not used with --from-molecule-bam");
        }
        if reporting {
            let (extracted, identity) =
                crate::moleculebam::read_molecule_bam_with_identity(&args.input)?;
            raw_input_identity = Some(report_file_identity(identity));
            extracted
        } else {
            crate::moleculebam::read_molecule_bam(&args.input)?
        }
    } else {
        if args.whitelist.is_some() {
            bail!("--whitelist requires --from-bam");
        }
        if reporting {
            let (extracted, identity) = read_archive_with_identity(&args.input)?;
            decoded_archive_identity = Some(identity);
            extracted
        } else {
            read_archive(&args.input)?
        }
    };
    let t_load = t0.elapsed().as_secs_f32();
    let (anno, annotation_identity) = load_replay_annotation(&args.gtf, reporting)?;
    let t_anno = t0.elapsed().as_secs_f32();
    if args.audit_multigene {
        if !reporting {
            crate::rows::multigene_audit_stranded(&x, &anno, solo_strand);
            return Ok(());
        }
        let audit = crate::rows::multigene_audit_summary_stranded(&x, &anno, solo_strand);
        let source_archive_identity = decoded_archive_identity.as_ref();
        let input_identity = if let Some(identity) = source_archive_identity {
            serde_json::to_value(identity)?
        } else {
            serde_json::to_value(
                raw_input_identity
                    .as_ref()
                    .context("uniform replay audit lacks an input identity")?,
            )?
        };
        let annotation_identity = annotation_identity
            .as_ref()
            .expect("reporting replay bound annotation identity");
        let context = replay_report_context(
            &args,
            source_archive_identity,
            annotation_identity,
            Vec::new(),
        )?;
        let summary = ReplayReportSummary {
            source_kind: if args.from_bam {
                "bam"
            } else if args.from_molecule_bam {
                "molecule_bam"
            } else {
                "archive"
            },
            input_identity,
            annotation_identity: annotation_identity.clone(),
            barcode_identity: None,
            whitelist_identity: whitelist_snapshot
                .as_ref()
                .map(|(_, identity)| identity.clone()),
            strand_policy: replay_strand_name(args.solo_strand),
            velocity: false,
            multigene_audit: true,
            molecules: x.mols.len() as u64,
            assigned_molecules: None,
            counted_umis: None,
            matrix_entries: None,
            audit: Some(audit),
            artifact: None,
        };
        emit_replay_report(&args, &summary, &context, None)?;
        return Ok(());
    }
    if args.velocity {
        let (vc, n_counted) = crate::rows::velocity_rows_stranded(&x, &anno, solo_strand);
        let artifact = emit_velocity(
            &vc,
            &x,
            &anno,
            &args.barcodes,
            barcode_snapshot.as_ref().map(|(text, _)| text.as_str()),
            &args.out_dir,
            reporting,
        )?;
        if reporting {
            let source_archive_identity = decoded_archive_identity.as_ref();
            let input_identity = if let Some(identity) = source_archive_identity {
                serde_json::to_value(identity)?
            } else {
                serde_json::to_value(
                    raw_input_identity
                        .as_ref()
                        .context("uniform velocity replay lacks an input identity")?,
                )?
            };
            let annotation_identity = annotation_identity
                .as_ref()
                .expect("reporting replay bound annotation identity");
            let context = replay_report_context(
                &args,
                source_archive_identity,
                annotation_identity,
                vec![
                    "MEX component files retain the legacy in-place directory behavior; metadata.json is written last as a completion marker, but the directory is not atomically published."
                        .into(),
                ],
            )?;
            let summary = ReplayReportSummary {
                source_kind: if args.from_bam {
                    "bam"
                } else if args.from_molecule_bam {
                    "molecule_bam"
                } else {
                    "archive"
                },
                input_identity,
                annotation_identity: annotation_identity.clone(),
                barcode_identity: barcode_snapshot
                    .as_ref()
                    .map(|(_, identity)| identity.clone()),
                whitelist_identity: whitelist_snapshot
                    .as_ref()
                    .map(|(_, identity)| identity.clone()),
                strand_policy: replay_strand_name(args.solo_strand),
                velocity: true,
                multigene_audit: false,
                molecules: x.mols.len() as u64,
                assigned_molecules: None,
                counted_umis: Some(n_counted),
                matrix_entries: Some(vc.len() as u64),
                audit: None,
                artifact,
            };
            emit_replay_report(&args, &summary, &context, metadata_path.as_deref())?;
        }
        eprintln!(
            "velocity: {} molecules -> {} UMIs in {} (cell,gene) entries | load {t_load:.1}s, total {:.1}s",
            x.mols.len(), n_counted, vc.len(), t0.elapsed().as_secs_f32()
        );
        return Ok(());
    }
    let (counts, n_assigned, total) = replay_rows_stranded(&x, &anno, solo_strand);
    let t_replay = t0.elapsed().as_secs_f32();
    let artifact = emit_matrix(
        &counts,
        &x,
        &anno,
        &args.barcodes,
        barcode_snapshot.as_ref().map(|(text, _)| text.as_str()),
        &args.out_dir,
        reporting,
    )?;
    if reporting {
        let source_archive_identity = decoded_archive_identity.as_ref();
        let input_identity = if let Some(identity) = source_archive_identity {
            serde_json::to_value(identity)?
        } else {
            serde_json::to_value(
                raw_input_identity
                    .as_ref()
                    .context("uniform gene replay lacks an input identity")?,
            )?
        };
        let annotation_identity = annotation_identity
            .as_ref()
            .expect("reporting replay bound annotation identity");
        let context = replay_report_context(
            &args,
            source_archive_identity,
            annotation_identity,
            vec![
                "MEX component files retain the legacy in-place directory behavior; metadata.json is written last as a completion marker, but the directory is not atomically published."
                    .into(),
            ],
        )?;
        let summary = ReplayReportSummary {
            source_kind: if args.from_bam {
                "bam"
            } else if args.from_molecule_bam {
                "molecule_bam"
            } else {
                "archive"
            },
            input_identity,
            annotation_identity: annotation_identity.clone(),
            barcode_identity: barcode_snapshot
                .as_ref()
                .map(|(_, identity)| identity.clone()),
            whitelist_identity: whitelist_snapshot
                .as_ref()
                .map(|(_, identity)| identity.clone()),
            strand_policy: replay_strand_name(args.solo_strand),
            velocity: false,
            multigene_audit: false,
            molecules: x.mols.len() as u64,
            assigned_molecules: Some(n_assigned),
            counted_umis: Some(total),
            matrix_entries: Some(counts.len() as u64),
            audit: None,
            artifact,
        };
        emit_replay_report(&args, &summary, &context, metadata_path.as_deref())?;
    }
    eprintln!(
        "molecules {} -> assigned {} -> {} UMIs in {} entries | eager | load {t_load:.1}s, +anno {:.1}s, +replay {:.1}s, total {:.1}s",
        x.mols.len(), n_assigned, total, counts.len(),
        t_anno - t_load, t_replay - t_anno, t0.elapsed().as_secs_f32()
    );
    exit_without_teardown();
}

#[derive(Parser)]
pub struct StampGenomeArgs {
    /// Archive to stamp.
    pub archive: PathBuf,
    /// Reference genome FASTA (plain or gzipped) the archive's reads were aligned to.
    #[arg(long)]
    pub genome: PathBuf,
    /// Write the stamped archive here instead of replacing the input in place.
    #[arg(long)]
    pub out: Option<PathBuf>,
    #[command(flatten)]
    pub report: UniformArchiveReportArgs,
}

#[derive(Parser)]
pub struct SealArchiveArgs {
    /// Legacy seekable v1 archive whose compressed sections will be copied exactly.
    pub archive: PathBuf,
    /// New root-committed v2 archive. Existing files are never overwritten.
    #[arg(long)]
    pub out: PathBuf,
    /// Emit one versioned JSON record instead of a human-readable summary.
    #[arg(long, conflicts_with = "report_format")]
    pub json: bool,
    #[command(flatten)]
    pub report: UniformArchiveReportArgs,
}

#[derive(Parser)]
pub struct InspectArchiveArgs {
    pub archive: PathBuf,
    /// Read and verify every compressed payload. Normal inspection authenticates only the v2
    /// directory/root; ordinary archive operations verify each selected payload lazily.
    #[arg(long)]
    pub verify_content: bool,
    #[arg(long, conflicts_with = "format")]
    pub json: bool,
    /// Use the versioned uniform result contract instead of the legacy presentation.
    #[arg(long, value_enum)]
    pub format: Option<UniformArchiveFormat>,
    /// Atomically publish uniform output here without replacing an existing file.
    #[arg(long, requires = "format")]
    pub output: Option<PathBuf>,
}

const INGEST_REPORT_SCHEMA: &str = "gravlax.archive.ingest-report.v1";
const REPLAY_REPORT_SCHEMA: &str = "gravlax.archive.replay-report.v1";
const SEAL_REPORT_SCHEMA: &str = "gravlax.archive.seal-report.v1";
const INSPECT_REPORT_SCHEMA: &str = "gravlax.archive.inspect-report.v1";
const STAMP_REPORT_SCHEMA: &str = "gravlax.archive.stamp-genome-report.v1";
const SECTION_TABLE_SCHEMA: &str = "gravlax.archive.section-accounting.v1";
const REPLAY_FILE_TABLE_SCHEMA: &str = "gravlax.archive.replay-artifact-files.v1";
const REPLAY_MEX_SCHEMA: &str = "gravlax.replay.mex-artifact.v1";

#[derive(Clone, Debug, Serialize)]
pub(crate) struct DigestIdentity {
    pub(crate) scheme: &'static str,
    pub(crate) blake3: String,
}

#[derive(Clone, Debug, Serialize)]
struct FileContentIdentity {
    scheme: &'static str,
    blake3: String,
    bytes: u64,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct ArchiveContentIdentity {
    pub(crate) format_version: u32,
    pub(crate) file_bytes: u64,
    pub(crate) native_identity: DigestIdentity,
    pub(crate) encoded_sections_identity: DigestIdentity,
}

impl ArchiveContentIdentity {
    pub(crate) fn provenance_identity(&self) -> String {
        format!(
            "{}:{}",
            self.native_identity.scheme, self.native_identity.blake3
        )
    }
}

#[derive(Clone, Debug, Serialize)]
struct ArtifactFileIdentity {
    name: String,
    bytes: u64,
    blake3: String,
}

#[derive(Clone, Debug, Serialize)]
struct ReplayArtifactManifest {
    /// The canonical single-matrix view consumed by the shared MEX reader. For velocity this is
    /// the spliced component; `matrices` retains the complete three-component manifest.
    #[serde(flatten)]
    primary: MexManifest,
    transactional_publication: bool,
    completion_marker: &'static str,
    matrices: Vec<MexManifest>,
    files: Vec<ArtifactFileIdentity>,
}

#[derive(Clone, Debug, Serialize)]
struct ReplayReportSummary {
    source_kind: &'static str,
    input_identity: serde_json::Value,
    annotation_identity: DigestIdentity,
    barcode_identity: Option<FileContentIdentity>,
    whitelist_identity: Option<FileContentIdentity>,
    strand_policy: &'static str,
    velocity: bool,
    multigene_audit: bool,
    molecules: u64,
    assigned_molecules: Option<u64>,
    counted_umis: Option<u64>,
    matrix_entries: Option<u64>,
    audit: Option<crate::rows::MultigeneAuditSummary>,
    artifact: Option<ReplayArtifactManifest>,
}

fn path_parameter(path: &Path) -> Result<serde_json::Value> {
    let value = path.to_str().with_context(|| {
        format!(
            "uniform reports require UTF-8 paths; {} is not representable",
            path.display()
        )
    })?;
    Ok(serde_json::Value::String(value.to_owned()))
}

fn output_path_parameter(path: &Path) -> Result<serde_json::Value> {
    Ok(serde_json::json!(reported_output_path(path)?))
}

fn preflight_new_report_path(path: &Path) -> Result<()> {
    if path.file_name().is_none() {
        bail!("uniform report output must name a file: {}", path.display());
    }
    match std::fs::symlink_metadata(path) {
        Ok(_) => bail!(
            "refusing to replace existing report output {}",
            path.display()
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error).with_context(|| format!("checking report output {}", path.display()))
        }
    }
    let parent = path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let metadata = std::fs::metadata(parent)
        .with_context(|| format!("checking report output directory {}", parent.display()))?;
    if !metadata.is_dir() {
        bail!(
            "report output parent is not a directory: {}",
            parent.display()
        );
    }
    Ok(())
}

fn preflight_artifact_report(
    report: &UniformArchiveReportArgs,
    primary_outputs: &[&Path],
) -> Result<bool> {
    match (report.report_format, report.report_output.as_deref()) {
        (None, None) => Ok(false),
        (None, Some(_)) => bail!("--report-output requires --report-format"),
        (Some(_), output) => {
            for path in primary_outputs {
                path_parameter(path)?;
            }
            if let Some(output) = output {
                let output_key = canonical_destination_key(output)?;
                if primary_outputs.iter().try_fold(false, |collision, primary| {
                    Ok::<_, OutputError>(
                        collision || canonical_destination_key(primary)? == output_key,
                    )
                })? {
                    bail!(
                        "operation report path must differ from the primary artifact path {}",
                        output.display()
                    );
                }
                path_parameter(output)?;
                preflight_new_report_path(output)?;
            }
            Ok(true)
        }
    }
}

fn preflight_inspect_output(args: &InspectArchiveArgs) -> Result<bool> {
    match (args.format, args.output.as_deref()) {
        (None, None) => Ok(false),
        (None, Some(_)) => bail!("--output requires --format"),
        (Some(_), output) => {
            path_parameter(&args.archive)?;
            if let Some(output) = output {
                path_parameter(output)?;
                preflight_new_report_path(output)?;
            }
            Ok(true)
        }
    }
}

fn write_uniform_output<F>(
    format: UniformArchiveFormat,
    output: Option<&Path>,
    render: F,
) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> std::result::Result<(), OutputError>,
{
    if let Some(path) = output {
        let mut render = Some(render);
        let outcome = publish_file_no_clobber(path, Durability::Flush, |writer| {
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

fn write_artifact_report<F>(report: &UniformArchiveReportArgs, render: F) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> std::result::Result<(), OutputError>,
{
    let format = report
        .report_format
        .context("uniform operation report requires --report-format")?;
    write_uniform_output(format, report.report_output.as_deref(), render)
}

fn report_file_identity(identity: ConsumedFileIdentity) -> FileContentIdentity {
    FileContentIdentity {
        scheme: "full-file-blake3-v1",
        blake3: identity.blake3,
        bytes: identity.bytes,
    }
}

type GenomeInputHandle = std::thread::JoinHandle<
    Result<(evidence_io::genome::GenomeSig, FileContentIdentity)>,
>;

fn start_genome_input(path: &Path) -> Result<GenomeInputHandle> {
    // Open both streams before the consuming work starts and prove they name the same immutable
    // regular file. The FASTA signature and exact encoded-byte identity can then run concurrently
    // without either reopening a replaceable pathname.
    let signature_file = std::fs::File::open(path)
        .with_context(|| format!("opening {} for genome signature", path.display()))?;
    let signature_before = signature_file.metadata()?;
    if !signature_before.is_file() {
        bail!("genome input is not a regular file: {}", path.display());
    }
    let signature_guard = signature_file.try_clone()?;
    let identity_file = std::fs::File::open(path)
        .with_context(|| format!("opening {} for genome identity", path.display()))?;
    let identity_before = identity_file.metadata()?;
    if !archive_metadata_matches(&signature_before, &identity_before)? {
        bail!("{} changed while its input snapshot was opened", path.display());
    }
    let path_before = std::fs::metadata(path)
        .with_context(|| format!("checking genome input {}", path.display()))?;
    if !archive_metadata_matches(&signature_before, &path_before)? {
        bail!("{} was replaced while its input snapshot was opened", path.display());
    }
    let path = path.to_owned();
    Ok(std::thread::spawn(move || {
        let identity_path = path.clone();
        let identity_thread = std::thread::spawn(move || {
            identity_of_consumed_file(identity_file, identity_before, &identity_path)
        });
        let signature = evidence_io::genome::sig_from_fasta_file(signature_file);
        let signature_after = signature_guard.metadata();
        let path_after = std::fs::metadata(&path)
            .with_context(|| format!("rechecking genome input {}", path.display()));
        let identity = identity_thread
            .join()
            .map_err(|_| anyhow::anyhow!("genome identity thread panicked"))?;
        let signature = signature?;
        let signature_after = signature_after?;
        if !archive_metadata_matches(&signature_before, &signature_after)? {
            bail!("{} changed while its genome signature was computed", path.display());
        }
        if !archive_metadata_matches(&signature_after, &path_after?)? {
            bail!("{} was replaced while its genome signature was computed", path.display());
        }
        Ok((signature, report_file_identity(identity?)))
    }))
}

fn stable_utf8_file(path: &Path) -> Result<(String, FileContentIdentity)> {
    let mut file = std::fs::File::open(path)
        .with_context(|| format!("opening {} for a stable text snapshot", path.display()))?;
    let before = file.metadata()?;
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)?;
    let after = file.metadata()?;
    if !archive_metadata_matches(&before, &after)? || bytes.len() as u64 != after.len() {
        bail!(
            "{} changed while its text snapshot was loaded",
            path.display()
        );
    }
    let path_metadata = std::fs::metadata(path)
        .with_context(|| format!("rechecking text input {}", path.display()))?;
    if !archive_inode_matches(&after, &path_metadata) {
        bail!(
            "{} was replaced while its text snapshot was loaded",
            path.display()
        );
    }
    let identity = FileContentIdentity {
        scheme: "full-file-blake3-v1",
        blake3: blake3::hash(&bytes).to_hex().to_string(),
        bytes: bytes.len() as u64,
    };
    let text = String::from_utf8(bytes)
        .with_context(|| format!("{} is not valid UTF-8", path.display()))?;
    Ok((text, identity))
}

fn verified_identity(identity: &FileContentIdentity) -> VerifiedFileIdentity {
    VerifiedFileIdentity::full_file_blake3(identity.blake3.clone(), identity.bytes)
}

fn stable_verified_file(path: &Path) -> Result<VerifiedFileIdentity> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("opening {} for content identity", path.display()))?;
    let before = file.metadata()?;
    if !before.is_file() {
        bail!(
            "alignment provenance input is not a regular file: {}",
            path.display()
        );
    }
    let identity = identity_of_consumed_file(file, before, path)?;
    Ok(VerifiedFileIdentity::full_file_blake3(
        identity.blake3,
        identity.bytes,
    ))
}

fn declared_alignment_file(path: &Path) -> Result<DeclaredAlignmentFile> {
    let locator = path
        .to_str()
        .with_context(|| format!("alignment provenance path is not UTF-8: {}", path.display()))?
        .to_owned();
    Ok(DeclaredAlignmentFile {
        relationship_status: ProvenanceStatus::DeclaredByCaller,
        locator,
        identity: stable_verified_file(path)?,
    })
}

fn junction_catalogue_data_rows(bytes: &[u8], label: &str) -> Result<u64> {
    let text = std::str::from_utf8(bytes)
        .with_context(|| format!("junction catalogue {label} is not valid UTF-8"))?;
    let mut data_rows = 0u64;
    for line in text.lines() {
        let row = line.trim_end_matches('\r');
        if row.trim().is_empty() || row.trim_start().starts_with('#') {
            continue;
        }
        if row.split('\t').count() < 4 {
            bail!("junction catalogue {label} contains a non-tabular data row");
        }
        data_rows = data_rows
            .checked_add(1)
            .context("junction catalogue row count overflow")?;
    }
    Ok(data_rows)
}

fn junction_catalogue(
    path: &Path,
    role: JunctionCatalogueRole,
) -> Result<(JunctionCatalogue, Vec<u8>)> {
    let mut file = std::fs::File::open(path)
        .with_context(|| format!("opening junction catalogue {}", path.display()))?;
    let before = file.metadata()?;
    if !before.is_file() {
        bail!(
            "junction catalogue is not a regular file: {}",
            path.display()
        );
    }
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)?;
    let data_rows = junction_catalogue_data_rows(&bytes, &path.display().to_string())?;
    let identity = identity_of_consumed_file(file, before, path)?;
    Ok((
        JunctionCatalogue {
            relationship_status: ProvenanceStatus::DeclaredByCaller,
            role,
            section: JUNCTION_CATALOGUE_SECTION.into(),
            identity: VerifiedFileIdentity::full_file_blake3(identity.blake3, identity.bytes),
            data_rows,
        },
        bytes,
    ))
}

fn validate_alignment_args(args: &IngestArgs) -> Result<()> {
    match (args.junction_discovery, args.junction_catalogue.as_ref()) {
        (JunctionDiscoveryArg::PerLibraryTwoPass | JunctionDiscoveryArg::FrozenCatalogue, None) => {
            bail!(
                "--junction-discovery {} requires --junction-catalogue",
                match args.junction_discovery {
                    JunctionDiscoveryArg::PerLibraryTwoPass => "per-library-two-pass",
                    JunctionDiscoveryArg::FrozenCatalogue => "frozen-catalogue",
                    _ => unreachable!(),
                }
            )
        }
        (JunctionDiscoveryArg::OnePass | JunctionDiscoveryArg::Unspecified, Some(_)) => {
            bail!(
                "--junction-catalogue requires per-library-two-pass or frozen-catalogue discovery"
            )
        }
        _ => {}
    }
    if args
        .alignment_index_identity
        .as_deref()
        .is_some_and(|value| value.trim().is_empty())
    {
        bail!("--alignment-index-identity must not be empty");
    }
    if args
        .alignment_chemistry
        .as_deref()
        .is_some_and(|value| value.trim().is_empty())
    {
        bail!("--alignment-chemistry must not be empty");
    }
    Ok(())
}

fn preflight_replay_metadata(
    out_dir: &Path,
    report_output: Option<&Path>,
    velocity: bool,
) -> Result<PathBuf> {
    let metadata = out_dir.join("metadata.json");
    if let Some(output) = report_output {
        let mut component_names = vec!["metadata.json", "features.tsv", "barcodes.tsv"];
        if velocity {
            component_names.extend([
                "spliced.mtx",
                "unspliced.mtx",
                "ambiguous.mtx",
                "entries.complete.tsv",
            ]);
        } else {
            component_names.push("matrix.mtx");
        }
        let exact_collision = component_names
            .iter()
            .map(|name| out_dir.join(name))
            .any(|component| output == component);
        let physical_collision = if out_dir.exists() {
            let root = std::fs::canonicalize(out_dir)?;
            let output_name = output.file_name();
            let output_parent = output
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            let output_key = output_name
                .map(|name| std::fs::canonicalize(output_parent).map(|parent| parent.join(name)))
                .transpose()?;
            output_key.is_some_and(|key| {
                component_names
                    .iter()
                    .any(|name| key == root.join(name))
            })
        } else {
            false
        };
        if exact_collision || physical_collision {
            bail!(
                "--report-output must differ from every replay artifact component under {}",
                out_dir.display()
            );
        }
    }
    match std::fs::symlink_metadata(&metadata) {
        Ok(_) => bail!(
            "refusing to replace existing replay metadata {}",
            metadata.display()
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => {
            return Err(error)
                .with_context(|| format!("checking replay metadata {}", metadata.display()))
        }
    }
    if let Ok(existing) = std::fs::metadata(out_dir) {
        if !existing.is_dir() {
            bail!("replay output is not a directory: {}", out_dir.display());
        }
    }
    Ok(metadata)
}

fn publish_replay_metadata(
    path: &Path,
    context: ResultContext,
    manifest: &ReplayArtifactManifest,
) -> Result<()> {
    let envelope = ResultEnvelope::new(REPLAY_MEX_SCHEMA, context, manifest)?;
    let outcome = publish_file_no_clobber(path, Durability::Flush, |writer| {
        serde_json::to_writer_pretty(&mut *writer, &envelope)?;
        writer.write_all(b"\n")?;
        Ok(())
    })?;
    for warning in outcome.warnings {
        eprintln!("warning: {warning}");
    }
    Ok(())
}

fn archive_identity(reader: &SectionReader) -> Result<ArchiveContentIdentity> {
    let before = reader.file_metadata()?;
    let version = reader.archive_version();
    let (native_scheme, native_digest, encoded_digest) =
        if let Some(root) = reader.content_commitment() {
            let encoded = reader
                .encoded_content_identity()?
                .context("rooted archive lacks encoded-section identity")?;
            ("aie-directory-root-v2", root.digest, encoded)
        } else {
            let scan = reader.scan_legacy_identities()?;
            (
                "full-file-blake3-v1",
                scan.full_file_blake3,
                scan.encoded_sections_blake3,
            )
        };
    let after = reader.file_metadata()?;
    if !archive_metadata_matches(&before, &after)? {
        bail!("archive changed while its exact identity was computed");
    }
    Ok(ArchiveContentIdentity {
        format_version: version,
        file_bytes: after.len(),
        native_identity: DigestIdentity {
            scheme: native_scheme,
            blake3: digest_hex(native_digest),
        },
        encoded_sections_identity: DigestIdentity {
            scheme: "aie-encoded-sections-v1",
            blake3: digest_hex(encoded_digest),
        },
    })
}

fn validate_archive_input(
    reader: &SectionReader,
    before: &std::fs::Metadata,
    path: &Path,
) -> Result<()> {
    let after = reader.file_metadata()?;
    if !archive_metadata_matches(before, &after)? {
        bail!("{} changed while its archive content was consumed", path.display());
    }
    let path_after = std::fs::metadata(path)
        .with_context(|| format!("rechecking consumed archive {}", path.display()))?;
    if !archive_metadata_matches(&after, &path_after)? {
        bail!("{} was replaced while its archive content was consumed", path.display());
    }
    Ok(())
}

fn report_context(
    archives: Vec<String>,
    parameters: BTreeMap<String, serde_json::Value>,
    warnings: Vec<String>,
) -> ResultContext {
    ResultContext {
        provenance: Provenance {
            archives,
            parameters,
            ..Provenance::default()
        },
        warnings,
        ..ResultContext::default()
    }
}

fn section_table_schema() -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(
        SECTION_TABLE_SCHEMA,
        vec![
            Field::new("section", DataType::String),
            Field::new("raw_bytes", DataType::UInt64),
            Field::new("compressed_bytes", DataType::UInt64),
        ],
    )?
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["section"]))
}

fn replay_file_table_schema() -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(
        REPLAY_FILE_TABLE_SCHEMA,
        vec![
            Field::new("file", DataType::String),
            Field::new("bytes", DataType::UInt64),
            Field::new("blake3", DataType::String),
        ],
    )?
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["file"]))
}

fn replay_strand_name(strand: SoloStrandArg) -> &'static str {
    match strand {
        SoloStrandArg::Forward => "forward",
        SoloStrandArg::Reverse => "reverse",
        SoloStrandArg::Unstranded => "unstranded",
    }
}

fn load_replay_annotation(
    path: &Path,
    reporting: bool,
) -> Result<(anno::Annotation, Option<DigestIdentity>)> {
    if !reporting {
        return Ok((anno::Annotation::from_path(path)?, None));
    }
    let label = path
        .to_str()
        .context("uniform replay reports require a UTF-8 annotation path")?;
    let identity = anno::intent::AnnotationIdentity::new("unspecified", label)?;
    let bound = anno::intent::BoundAnnotation::from_path(path, identity)?;
    let (annotation, identity) = bound.into_parts();
    let digest = identity
        .digest
        .context("bound replay annotation did not report a content digest")?;
    let digest = digest
        .strip_prefix("blake3:")
        .context("bound replay annotation returned a non-BLAKE3 digest")?;
    Ok((
        annotation,
        Some(DigestIdentity {
            scheme: "full-file-blake3-v1",
            blake3: digest.to_owned(),
        }),
    ))
}

fn replay_report_context(
    args: &ReplayRowsArgs,
    archive_identity: Option<&ArchiveContentIdentity>,
    annotation_identity: &DigestIdentity,
    warnings: Vec<String>,
) -> Result<ResultContext> {
    let mut parameters = BTreeMap::new();
    parameters.insert("input".into(), path_parameter(&args.input)?);
    parameters.insert("annotation".into(), path_parameter(&args.gtf)?);
    parameters.insert("barcodes".into(), path_parameter(&args.barcodes)?);
    parameters.insert(
        "output_directory".into(),
        output_path_parameter(&args.out_dir)?,
    );
    parameters.insert("from_bam".into(), serde_json::json!(args.from_bam));
    parameters.insert(
        "from_molecule_bam".into(),
        serde_json::json!(args.from_molecule_bam),
    );
    parameters.insert("velocity".into(), serde_json::json!(args.velocity));
    parameters.insert(
        "audit_multigene".into(),
        serde_json::json!(args.audit_multigene),
    );
    parameters.insert("eager".into(), serde_json::json!(args.eager));
    parameters.insert("locus_gap".into(), serde_json::json!(args.locus_gap));
    parameters.insert(
        "solo_strand".into(),
        serde_json::json!(replay_strand_name(args.solo_strand)),
    );
    if let Some(whitelist) = &args.whitelist {
        parameters.insert("whitelist".into(), path_parameter(whitelist)?);
    }
    let archives = archive_identity
        .map(|identity| vec![identity.provenance_identity()])
        .unwrap_or_default();
    let mut context = report_context(archives, parameters, warnings);
    context.provenance.assembly = Some("unspecified".into());
    context.provenance.annotation = Some(
        args.gtf
            .to_str()
            .context("uniform replay reports require a UTF-8 annotation path")?
            .to_owned(),
    );
    context.provenance.annotation_digest = Some(format!("blake3:{}", annotation_identity.blake3));
    Ok(context)
}

fn emit_replay_report(
    args: &ReplayRowsArgs,
    summary: &ReplayReportSummary,
    context: &ResultContext,
    metadata_path: Option<&Path>,
) -> Result<()> {
    if let (Some(path), Some(artifact)) = (metadata_path, summary.artifact.as_ref()) {
        publish_replay_metadata(path, context.clone(), artifact)?;
    }
    let schema = replay_file_table_schema()?;
    let format = args
        .report
        .report_format
        .expect("report format preflighted");
    write_artifact_report(&args.report, |writer| {
        let mut bundle = StreamingBundleWriter::new_with_summary(
            writer,
            REPLAY_REPORT_SCHEMA,
            format.into(),
            context,
            summary,
        )?;
        if let Some(artifact) = &summary.artifact {
            bundle.write_table("artifact_files", &schema, None, |table| {
                for file in &artifact.files {
                    table.write_row_with(|row| {
                        row.string(&file.name)?;
                        row.uint64(file.bytes)?;
                        row.string(&file.blake3)
                    })?;
                }
                Ok(())
            })?;
        }
        bundle.finish()?;
        Ok(())
    })
}

fn emit_stamp_report(
    args: &StampGenomeArgs,
    changed: bool,
    source: &ArchiveContentIdentity,
    output: &ArchiveContentIdentity,
    genome_file: &FileContentIdentity,
    genome_signature: &evidence_io::genome::GenomeSig,
    sections: &[(String, u64, u64, u64)],
) -> Result<()> {
    #[derive(Serialize)]
    struct StampSummary<'a> {
        changed: bool,
        source_archive: &'a ArchiveContentIdentity,
        output_archive: &'a ArchiveContentIdentity,
        genome_file: &'a FileContentIdentity,
        genome_signature: &'a evidence_io::genome::GenomeSig,
        sections: u64,
        compressed_sections_copied_without_reencoding: u64,
    }
    let summary = StampSummary {
        changed,
        source_archive: source,
        output_archive: output,
        genome_file,
        genome_signature,
        sections: sections.len() as u64,
        compressed_sections_copied_without_reencoding: if changed {
            sections
                .iter()
                .filter(|(name, _, _, _)| name != "meta")
                .count() as u64
        } else {
            sections.len() as u64
        },
    };
    let out = args.out.as_deref().unwrap_or(&args.archive);
    let mut parameters = BTreeMap::new();
    parameters.insert("archive".into(), path_parameter(&args.archive)?);
    parameters.insert("genome".into(), path_parameter(&args.genome)?);
    parameters.insert("output".into(), output_path_parameter(out)?);
    parameters.insert("in_place".into(), serde_json::json!(args.out.is_none()));
    let mut archives = vec![source.provenance_identity()];
    let output_identity = output.provenance_identity();
    if !archives.contains(&output_identity) {
        archives.push(output_identity);
    }
    let context = report_context(archives, parameters, Vec::new());
    let schema = section_table_schema()?;
    let format = args
        .report
        .report_format
        .expect("report format preflighted");
    write_artifact_report(&args.report, |writer| {
        let mut bundle = StreamingBundleWriter::new_with_summary(
            writer,
            STAMP_REPORT_SCHEMA,
            format.into(),
            &context,
            &summary,
        )?;
        bundle.write_table("sections", &schema, None, |table| {
            for (name, _, raw, compressed) in sections {
                table.write_row_with(|row| {
                    row.string(name)?;
                    row.uint64(*raw)?;
                    row.uint64(*compressed)
                })?;
            }
            Ok(())
        })?;
        bundle.finish()?;
        Ok(())
    })
}

fn digest_hex(digest: [u8; 32]) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn temporary_section_writer(
    out: &Path,
    tag: &str,
    level: i32,
) -> Result<(PathBuf, SectionWriter)> {
    let file_name = out
        .file_name()
        .context("archive output has no file name")?
        .to_string_lossy();
    for attempt in 0u32..1_000 {
        let temporary = out.with_file_name(format!(
            ".{file_name}.{tag}-{}-{attempt}",
            std::process::id()
        ));
        match SectionWriter::create_new(&temporary, level) {
            Ok(writer) => return Ok((temporary, writer)),
            Err(error)
                if error
                    .downcast_ref::<std::io::Error>()
                    .is_some_and(|io| io.kind() == std::io::ErrorKind::AlreadyExists) => {}
            Err(error) => return Err(error),
        }
    }
    bail!(
        "could not allocate a unique temporary archive beside {}",
        out.display()
    )
}

fn temporary_ingest_writer(out: &Path, level: i32) -> Result<(PathBuf, SectionWriter)> {
    temporary_section_writer(out, "ingest-tmp", level)
}

fn temporary_archive_writer(out: &Path) -> Result<(PathBuf, SectionWriter)> {
    temporary_section_writer(out, "seal-tmp", 0)
}

fn temporary_stamp_writer(out: &Path) -> Result<(PathBuf, SectionWriter)> {
    temporary_section_writer(out, "stamp-tmp", 19)
}

fn temporary_stamp_commit_path(out: &Path) -> Result<PathBuf> {
    let file_name = out
        .file_name()
        .context("stamp output has no file name")?
        .to_string_lossy();
    for attempt in 0u32..1_000 {
        let path = out.with_file_name(format!(
            ".{file_name}.stamp-commit-{}-{attempt}",
            std::process::id()
        ));
        match std::fs::symlink_metadata(&path) {
            Ok(_) => continue,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(path),
            Err(error) => return Err(error.into()),
        }
    }
    bail!(
        "could not allocate a unique stamp commit path beside {}",
        out.display()
    )
}

fn preflight_stamp_output(args: &StampGenomeArgs) -> Result<()> {
    let Some(out) = args.out.as_deref() else {
        return Ok(());
    };
    if canonical_destination_key(out)? == canonical_destination_key(&args.archive)? {
        bail!("--out must differ from the source archive; omit it for in-place stamping");
    }
    match std::fs::symlink_metadata(out) {
        Ok(_) => bail!("refusing to overwrite stamp output {}", out.display()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error)
            .with_context(|| format!("checking stamp output destination {}", out.display())),
    }
}

fn archive_metadata_matches(before: &std::fs::Metadata, after: &std::fs::Metadata) -> Result<bool> {
    if before.len() != after.len() || before.modified()? != after.modified()? {
        return Ok(false);
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        Ok(before.dev() == after.dev()
            && before.ino() == after.ino()
            && before.mtime() == after.mtime()
            && before.mtime_nsec() == after.mtime_nsec()
            && before.ctime() == after.ctime()
            && before.ctime_nsec() == after.ctime_nsec())
    }
    #[cfg(not(unix))]
    Ok(true)
}

/// Whether two metadata records name the same open regular file. Link creation updates ctime, so
/// this deliberately checks inode identity rather than the full immutability guard above.
fn archive_inode_matches(expected: &std::fs::Metadata, observed: &std::fs::Metadata) -> bool {
    if !expected.is_file() || !observed.is_file() || expected.len() != observed.len() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        expected.dev() == observed.dev() && expected.ino() == observed.ino()
    }
    #[cfg(not(unix))]
    {
        expected.modified().ok() == observed.modified().ok()
    }
}

/// Remove a staging name only when it still refers to the held produced file. A concurrent
/// replacement is foreign state and must be preserved, including on an error path.
pub(crate) fn remove_staging_if_owned(path: &Path, file: &std::fs::File) -> Result<bool> {
    let expected = file.metadata()?;
    let observed = match std::fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(true),
        Err(error) => return Err(error.into()),
    };
    if !archive_inode_matches(&expected, &observed) {
        return Ok(false);
    }
    std::fs::remove_file(path)?;
    Ok(true)
}

fn warn_on_failed_staging_cleanup(path: &Path, file: &std::fs::File) {
    match remove_staging_if_owned(path, file) {
        Ok(true) => {}
        Ok(false) => eprintln!(
            "warning: temporary path {} now names a different file and was preserved",
            path.display()
        ),
        Err(error) => eprintln!(
            "warning: could not remove owned temporary path {}: {error:#}",
            path.display()
        ),
    }
}

/// Byte-preserving v1 -> v2 migration. The one-time pass reads each compressed frame and writes
/// it unchanged; only the authenticated directory/footer are new.
pub fn run_seal_archive(args: SealArchiveArgs) -> Result<()> {
    let reporting = preflight_artifact_report(&args.report, &[&args.out])?;
    if reporting {
        path_parameter(&args.archive)?;
    }
    let started = std::time::Instant::now();
    if args.out.exists() {
        bail!("refusing to overwrite {}", args.out.display());
    }
    let source_file = std::fs::File::open(&args.archive)
        .with_context(|| format!("opening {}", args.archive.display()))?;
    let source_before = source_file.metadata()?;
    let mut reader = SectionReader::from_file(source_file)?;
    if reader.archive_version() != evidence_io::format::SEEKABLE_VERSION {
        bail!(
            "seal-archive requires a legacy v1 source; {} is archive version {}",
            args.archive.display(),
            reader.archive_version()
        );
    }
    let input_bytes = source_before.len();
    let entries = reader.entries().to_vec();
    for (name, _, raw_len, comp_len) in &entries {
        if *raw_len > evidence_io::format::MAX_V2_RAW_SECTION_BYTES {
            bail!("source section {name} raw length exceeds the v2 safety limit");
        }
        if *comp_len > evidence_io::format::MAX_V2_COMPRESSED_SECTION_BYTES {
            bail!("source section {name} compressed length exceeds the v2 safety limit");
        }
        let permitted_raw = comp_len
            .saturating_mul(4_096)
            .saturating_add(64 * 1024 * 1024);
        if *raw_len > permitted_raw {
            bail!("source section {name} declares an unsafe v2 compression ratio");
        }
    }
    let source_identity = reader.scan_legacy_identities()?;
    let (temporary, mut writer) = temporary_archive_writer(&args.out)?;
    let staging_guard = writer.try_clone_file()?;
    let result = (|| -> Result<_> {
        let mut compressed_bytes = 0u64;
        for (name, _, raw_len, expected_comp_len) in &entries {
            let (compressed, observed_raw_len) = reader.read_compressed(name)?;
            if observed_raw_len as u64 != *raw_len || compressed.len() as u64 != *expected_comp_len
            {
                bail!("source section {name} changed while sealing");
            }
            compressed_bytes = compressed_bytes
                .checked_add(compressed.len() as u64)
                .context("sealed payload byte count overflow")?;
            evidence_io::format::decompress(&compressed, observed_raw_len)
                .with_context(|| format!("validating legacy source section {name}"))?;
            writer.section_precompressed(name, *raw_len, &compressed)?;
        }
        let source_after = reader.file_metadata()?;
        if !archive_metadata_matches(&source_before, &source_after)? {
            bail!("source archive changed while it was being sealed");
        }
        let (finished_entries, output_file, root) = writer.finish_with_file()?;
        if finished_entries.len() != entries.len() {
            bail!("sealed archive section accounting differs from its v1 source");
        }
        let candidate = SectionReader::from_file(output_file.try_clone()?)?;
        let candidate_before = candidate.file_metadata()?;
        let parsed_root = candidate
            .content_commitment()
            .context("sealed archive lacks a v2 root")?;
        if parsed_root != root {
            bail!("sealed archive differs from its computed root commitment");
        }
        let content = candidate
            .encoded_content_identity()?
            .context("sealed archive lacks encoded-section identity")?;
        if content != source_identity.encoded_sections_blake3 {
            bail!("sealed archive encoded-section identity differs from its v1 source");
        }
        let verified_payload_bytes = candidate.verify_all_payloads()?;
        if verified_payload_bytes != compressed_bytes {
            bail!("sealed archive payload accounting differs after verification");
        }
        let candidate_after = candidate.file_metadata()?;
        if !archive_metadata_matches(&candidate_before, &candidate_after)? {
            bail!("temporary archive changed while it was being verified");
        }
        let output_bytes = candidate_after.len();
        let outcome = install_open_file_no_clobber(
            &output_file,
            &temporary,
            &args.out,
            Durability::FileAndDirectory,
        )?;
        for warning in outcome.warnings {
            eprintln!("warning: {warning}");
        }
        Ok((compressed_bytes, entries.len(), root, content, output_bytes))
    })();
    if result.is_err() {
        warn_on_failed_staging_cleanup(&temporary, &staging_guard);
    }
    let (compressed_bytes, sections, root, content, output_bytes) = result?;
    let value = serde_json::json!({
        "schema": "gravlax.archive.seal.v1",
        "input": args.archive,
        "output": args.out,
        "source_format_version": evidence_io::format::SEEKABLE_VERSION,
        "output_format_version": root.version,
        "sections": sections,
        "compressed_payload_bytes_copied": compressed_bytes,
        "input_bytes": input_bytes,
        "output_bytes": output_bytes,
        "source_full_file_blake3": digest_hex(source_identity.full_file_blake3),
        "source_identity_content_bytes_read": source_identity.bytes_read,
        "archive_root": {"scheme": "aie-directory-root-v2", "blake3": root.to_hex()},
        "encoded_sections_identity": {
            "scheme": "aie-encoded-sections-v1",
            "blake3": digest_hex(content),
        },
        "payload_verification_bytes_read": compressed_bytes,
        "elapsed_seconds": started.elapsed().as_secs_f64(),
    });
    if reporting {
        let source_archive = ArchiveContentIdentity {
            format_version: evidence_io::format::SEEKABLE_VERSION,
            file_bytes: input_bytes,
            native_identity: DigestIdentity {
                scheme: "full-file-blake3-v1",
                blake3: digest_hex(source_identity.full_file_blake3),
            },
            encoded_sections_identity: DigestIdentity {
                scheme: "aie-encoded-sections-v1",
                blake3: digest_hex(source_identity.encoded_sections_blake3),
            },
        };
        let output_archive = ArchiveContentIdentity {
            format_version: root.version,
            file_bytes: output_bytes,
            native_identity: DigestIdentity {
                scheme: "aie-directory-root-v2",
                blake3: root.to_hex(),
            },
            encoded_sections_identity: DigestIdentity {
                scheme: "aie-encoded-sections-v1",
                blake3: digest_hex(content),
            },
        };
        #[derive(Serialize)]
        struct SealSummary<'a> {
            source_archive: &'a ArchiveContentIdentity,
            output_archive: &'a ArchiveContentIdentity,
            sections: u64,
            compressed_payload_bytes_copied: u64,
            payload_verification_bytes_read: u64,
            source_identity_content_bytes_read: u64,
        }
        let summary = SealSummary {
            source_archive: &source_archive,
            output_archive: &output_archive,
            sections: sections as u64,
            compressed_payload_bytes_copied: compressed_bytes,
            payload_verification_bytes_read: compressed_bytes,
            source_identity_content_bytes_read: source_identity.bytes_read,
        };
        let mut parameters = BTreeMap::new();
        parameters.insert("source".into(), path_parameter(&args.archive)?);
        parameters.insert("output".into(), output_path_parameter(&args.out)?);
        let context = report_context(
            vec![
                source_archive.provenance_identity(),
                output_archive.provenance_identity(),
            ],
            parameters,
            Vec::new(),
        );
        let schema = section_table_schema()?;
        let format = args
            .report
            .report_format
            .expect("report format preflighted");
        write_artifact_report(&args.report, |writer| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                SEAL_REPORT_SCHEMA,
                format.into(),
                &context,
                &summary,
            )?;
            bundle.write_table("sections", &schema, None, |table| {
                for (name, _, raw, compressed) in &entries {
                    table.write_row_with(|row| {
                        row.string(name)?;
                        row.uint64(*raw)?;
                        row.uint64(*compressed)
                    })?;
                }
                Ok(())
            })?;
            bundle.finish()?;
            Ok(())
        })?;
    } else if args.json {
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        println!(
            "sealed {} v1 sections into {} ({} -> {} bytes, root {}) in {:.2}s",
            sections,
            args.out.display(),
            input_bytes,
            output_bytes,
            root.to_hex(),
            started.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

pub fn run_inspect_archive(args: InspectArchiveArgs) -> Result<()> {
    let uniform = preflight_inspect_output(&args)?;
    let started = std::time::Instant::now();
    let file = std::fs::File::open(&args.archive)
        .with_context(|| format!("opening {}", args.archive.display()))?;
    // Capture this before parsing: legacy inspection always performs a complete identity scan,
    // even without --verify-content, while explicit verification scans all payloads for either
    // format. Both long reads must describe one stable open-file snapshot.
    let inspection_before = file.metadata()?;
    let mut reader = SectionReader::from_file(file)?;
    let version = reader.archive_version();
    let file_bytes = reader.file_metadata()?.len();
    let section_count = reader.entries().len();
    let compressed_payload_bytes = reader.entries().iter().try_fold(0u64, |sum, entry| {
        sum.checked_add(entry.3)
            .context("compressed payload byte count overflow")
    })?;
    let (native_scheme, native_digest, encoded_digest, identity_content_bytes_read) =
        if let Some(root) = reader.content_commitment() {
            let encoded = reader
                .encoded_content_identity()?
                .context("rooted archive lacks encoded-section identity")?;
            let verified = if args.verify_content {
                reader.verify_all_payloads()?
            } else {
                0
            };
            ("aie-directory-root-v2", root.digest, encoded, verified)
        } else {
            let scan = reader.scan_legacy_identities()?;
            if args.verify_content {
                // Hash identity alone cannot validate the declared decompressed lengths. Explicit
                // full verification also exercises every legacy zstd frame.
                let names: Vec<String> = reader.names().map(str::to_owned).collect();
                for name in names {
                    reader.read(&name)?;
                }
            }
            (
                "full-file-blake3-v1",
                scan.full_file_blake3,
                scan.encoded_sections_blake3,
                scan.bytes_read,
            )
        };
    // `inspect-archive` also serves as the generic container inspector used for sealed legacy
    // section archives. Those archives need not carry Gravlax's application-level `meta`
    // section. Treat its absence as having no declared molecular-evidence capabilities; the
    // capability validator below still rejects any partial logical-v2 sections or declarations.
    let meta: serde_json::Value = if reader.has("meta") {
        serde_json::from_slice(&reader.read("meta")?)?
    } else {
        serde_json::json!({})
    };
    let (alignment_provenance, terminal_tail, genome_reference_binding) =
        archive_capabilities(&mut reader, &meta, true)?;
    let molecular_evidence_schema = meta
        .get("evidence_schema")
        .and_then(serde_json::Value::as_str);
    let alignment_provenance_status = if alignment_provenance.is_some() {
        "available"
    } else {
        "unavailable"
    };
    let terminal_tail_status = if terminal_tail.is_some() {
        "available"
    } else {
        "unavailable"
    };
    let genome_reference_binding_status = if genome_reference_binding.is_some() {
        "available"
    } else if molecular_evidence_schema.is_some() {
        "unavailable"
    } else {
        "legacy_unattributed"
    };
    if args.verify_content || version == evidence_io::format::SEEKABLE_VERSION {
        let after = reader.file_metadata()?;
        if !archive_metadata_matches(&inspection_before, &after)? {
            bail!("archive changed while its complete content was being inspected");
        }
    }
    let value = serde_json::json!({
        "schema": "gravlax.archive.identity.v1",
        "archive": args.archive,
        "format_version": version,
        "file_bytes": file_bytes,
        "sections": section_count,
        "native_identity": {"scheme": native_scheme, "blake3": digest_hex(native_digest)},
        "encoded_sections_identity": {
            "scheme": "aie-encoded-sections-v1",
            "blake3": digest_hex(encoded_digest),
        },
        "verification": {
            "directory_and_root": version == evidence_io::format::VERSION,
            "all_payloads": args.verify_content,
            "identity_content_bytes_read": identity_content_bytes_read,
            "ordinary_reads_verify_selected_payloads_only": version == evidence_io::format::VERSION,
        },
        "molecular_evidence": {
            "schema": molecular_evidence_schema,
            "alignment_provenance_status": alignment_provenance_status,
            "alignment_provenance": alignment_provenance,
            "terminal_tail_status": terminal_tail_status,
            "terminal_tail": terminal_tail,
            "genome_reference_binding_status": genome_reference_binding_status,
            "genome_reference_binding": genome_reference_binding,
        },
        "elapsed_seconds": started.elapsed().as_secs_f64(),
    });
    if uniform {
        let archive_identity = ArchiveContentIdentity {
            format_version: version,
            file_bytes,
            native_identity: DigestIdentity {
                scheme: native_scheme,
                blake3: digest_hex(native_digest),
            },
            encoded_sections_identity: DigestIdentity {
                scheme: "aie-encoded-sections-v1",
                blake3: digest_hex(encoded_digest),
            },
        };
        #[derive(Serialize)]
        struct InspectSummary<'a> {
            archive: &'a ArchiveContentIdentity,
            sections: u64,
            compressed_payload_bytes: u64,
            directory_and_root_verified: bool,
            all_payloads_verified: bool,
            identity_content_bytes_read: u64,
            ordinary_reads_verify_selected_payloads_only: bool,
            molecular_evidence_schema: Option<&'a str>,
            alignment_provenance_status: &'static str,
            alignment_provenance: Option<&'a AlignmentProvenanceManifest>,
            terminal_tail_status: &'static str,
            terminal_tail: Option<&'a TerminalTailMetadata>,
            genome_reference_binding_status: &'static str,
            genome_reference_binding: Option<&'a GenomeReferenceBinding>,
        }
        let summary = InspectSummary {
            archive: &archive_identity,
            sections: section_count as u64,
            compressed_payload_bytes,
            directory_and_root_verified: version == evidence_io::format::VERSION,
            all_payloads_verified: args.verify_content,
            identity_content_bytes_read,
            ordinary_reads_verify_selected_payloads_only: version == evidence_io::format::VERSION,
            molecular_evidence_schema,
            alignment_provenance_status,
            alignment_provenance: alignment_provenance.as_ref(),
            terminal_tail_status,
            terminal_tail: terminal_tail.as_ref(),
            genome_reference_binding_status,
            genome_reference_binding: genome_reference_binding.as_ref(),
        };
        let mut parameters = BTreeMap::new();
        parameters.insert("archive".into(), path_parameter(&args.archive)?);
        parameters.insert(
            "verify_content".into(),
            serde_json::json!(args.verify_content),
        );
        let context = report_context(
            vec![archive_identity.provenance_identity()],
            parameters,
            Vec::new(),
        );
        let schema = section_table_schema()?;
        let format = args.format.expect("uniform inspect format preflighted");
        write_uniform_output(format, args.output.as_deref(), |writer| {
            let mut bundle = StreamingBundleWriter::new_with_summary(
                writer,
                INSPECT_REPORT_SCHEMA,
                format.into(),
                &context,
                &summary,
            )?;
            bundle.write_table("sections", &schema, None, |table| {
                for (name, _, raw, compressed) in reader.entries() {
                    table.write_row_with(|row| {
                        row.string(name)?;
                        row.uint64(*raw)?;
                        row.uint64(*compressed)
                    })?;
                }
                Ok(())
            })?;
            bundle.finish()?;
            Ok(())
        })?;
    } else if args.json {
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        println!(
            "{}: archive v{} with {} sections and {} bytes",
            args.archive.display(),
            version,
            section_count,
            file_bytes
        );
        println!(
            "native identity: {native_scheme}:{}",
            digest_hex(native_digest)
        );
        println!(
            "encoded sections: aie-encoded-sections-v1:{}",
            digest_hex(encoded_digest)
        );
        match (&alignment_provenance, molecular_evidence_schema) {
            (Some(provenance), Some(schema)) => {
                println!("molecular evidence schema: {schema}");
                println!(
                    "alignment provenance: {} ({:?}, {:?})",
                    ALIGNMENT_PROVENANCE_SCHEMA,
                    provenance.alignment.junction_discovery,
                    provenance.alignment.status
                );
            }
            _ => {
                println!("molecular evidence schema: unavailable (legacy archive)");
                println!("alignment provenance: unavailable; junction discovery is unknown");
            }
        }
        if let Some(tails) = &terminal_tail {
            println!(
                "terminal tails: available ({} events on {} molecules in {} routed chunks; {})",
                tails.events, tails.selected_molecules, tails.chunks, tails.extraction_rule
            );
        } else {
            println!("terminal tails: unavailable (extraction rule was not recorded as evaluated)");
        }
        if let Some(binding) = &genome_reference_binding {
            println!(
                "genome reference binding: available ({:?}; caller-declared relationship)",
                binding.bound_by
            );
        } else if molecular_evidence_schema.is_some() {
            println!("genome reference binding: unavailable");
        } else {
            println!("genome reference binding: legacy/unattributed");
        }
        if args.verify_content {
            if version == evidence_io::format::VERSION {
                println!(
                    "verified all committed compressed payloads ({identity_content_bytes_read} bytes)"
                );
            } else {
                println!(
                    "completed a {identity_content_bytes_read}-byte full-file identity scan and decoded all {compressed_payload_bytes} compressed payload bytes"
                );
            }
        } else if version == evidence_io::format::VERSION {
            println!("verified directory/root; payloads will be verified when selected");
        } else {
            println!("legacy identity required a complete {identity_content_bytes_read}-byte scan");
        }
    }
    Ok(())
}

/// Retrofit a genome signature into an existing archive. Evidence sections are copied compressed
/// byte-for-byte. Logical-v2 archives record a separate current-reference binding in `meta`;
/// original alignment provenance is never rewritten by a later stamp.
pub fn run_stamp_genome(args: StampGenomeArgs) -> Result<()> {
    preflight_stamp_output(&args)?;
    let destination = args.out.as_deref().unwrap_or(&args.archive);
    let reporting = preflight_artifact_report(&args.report, &[destination])?;
    if reporting {
        path_parameter(&args.archive)?;
        path_parameter(&args.genome)?;
    }
    let t0 = std::time::Instant::now();
    let genome_thread = start_genome_input(&args.genome)?;
    let (sig, genome_identity) = genome_thread
        .join()
        .map_err(|_| anyhow::anyhow!("genome input thread panicked"))??;
    let mut r = SectionReader::open(&args.archive)?;
    let source_before = r.file_metadata()?;
    let source_identity = reporting.then(|| archive_identity(&r)).transpose()?;
    let mut meta: serde_json::Map<String, serde_json::Value> =
        serde_json::from_slice(&r.read("meta")?)?;
    let logical_v2 = meta
        .get("evidence_schema")
        .and_then(serde_json::Value::as_str)
        == Some(MOLECULAR_EVIDENCE_SCHEMA);
    let (_, _, previous_binding) =
        archive_capabilities(&mut r, &serde_json::Value::Object(meta.clone()), false)?;
    let target_binding = logical_v2.then(|| {
        GenomeReferenceBinding::new(
            GenomeBindingAction::StampGenome,
            verified_identity(&genome_identity),
            sig.clone(),
        )
    });
    let binding_changed = match (&previous_binding, &target_binding) {
        (Some(previous), Some(target)) => {
            previous.identity != target.identity || previous.signature != target.signature
        }
        (None, None) => false,
        _ => true,
    };
    let chroms: Vec<String> = String::from_utf8_lossy(&r.read("chroms")?)
        .lines()
        .map(|s| s.to_string())
        .collect();
    let missing: Vec<&String> = chroms.iter().filter(|n| sig.contig(n).is_none()).collect();
    if !missing.is_empty() {
        bail!(
            "--genome does not cover this archive: {} archive contig(s) absent from the FASTA \
             (first: {})",
            missing.len(),
            missing[0]
        );
    }
    if let Some(prev) = meta.get("genome_sig") {
        let prev: evidence_io::genome::GenomeSig = serde_json::from_value(prev.clone())?;
        if prev.digest == sig.digest && !binding_changed {
            let sections = r.entries().to_vec();
            if let Some(out) = args.out.as_deref() {
                let mut source = r.try_clone_file()?;
                source.seek(SeekFrom::Start(0))?;
                let outcome = publish_file_no_clobber(
                    out,
                    Durability::FileAndDirectory,
                    |writer| {
                        std::io::copy(&mut source, writer)?;
                        validate_archive_input(&r, &source_before, &args.archive).map_err(
                            |error| {
                                OutputError::Sink(format!(
                                    "source archive changed while it was copied: {error:#}"
                                ))
                            },
                        )?;
                        Ok(())
                    },
                )?;
                for warning in outcome.warnings {
                    eprintln!("warning: {warning}");
                }
            } else {
                validate_archive_input(&r, &source_before, &args.archive)?;
            }
            if reporting {
                let source_identity = source_identity
                    .as_ref()
                    .expect("reporting stamp captured source identity");
                emit_stamp_report(
                    &args,
                    false,
                    source_identity,
                    source_identity,
                    &genome_identity,
                    &sig,
                    &sections,
                )?;
            } else if let Some(out) = args.out.as_deref() {
                println!(
                    "archive already stamped with this genome (digest {}); copied unchanged to {}",
                    &sig.digest[..16],
                    out.display()
                );
            } else {
                println!(
                    "archive already stamped with this genome (digest {}); nothing to do",
                    &sig.digest[..16]
                );
            }
            return Ok(());
        }
        eprintln!("replacing existing genome signature (digest {} -> {})", &prev.digest[..16], &sig.digest[..16]);
    }
    meta.insert("genome_sig".into(), serde_json::to_value(&sig)?);
    if let Some(binding) = &target_binding {
        meta.insert(
            "genome_reference_binding".into(),
            serde_json::to_value(binding)?,
        );
    }
    let out = args.out.clone().unwrap_or_else(|| args.archive.clone());
    let (tmp, mut w) = temporary_stamp_writer(&out)?;
    let staging_guard = w.try_clone_file()?;
    let stamp_result = (|| -> Result<_> {
        let entries: Vec<(String, u64, u64, u64)> = r.entries().to_vec();
        for (name, _, raw_len, _) in &entries {
            if name == "meta" {
                w.section("meta", serde_json::to_string(&meta)?.as_bytes())?;
            } else {
                let (comp, rl) = r.read_compressed(name)?;
                debug_assert_eq!(rl as u64, *raw_len);
                w.section_precompressed(name, *raw_len, &comp)?;
            }
        }
        let (_accounting, output_file, commitment) = w.finish_with_file()?;
        validate_archive_input(&r, &source_before, &args.archive)?;
        let output_reader = SectionReader::from_file(output_file.try_clone()?)?;
        if output_reader.content_commitment() != Some(commitment) {
            bail!("stamped archive differs from its computed root commitment");
        }
        let output_identity = reporting.then(|| archive_identity(&output_reader)).transpose()?;
        let sections = reporting.then(|| output_reader.entries().to_vec());
        if args.out.is_some() {
            let outcome = install_open_file_no_clobber(
                &output_file,
                &tmp,
                &out,
                Durability::FileAndDirectory,
            )?;
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
        } else {
            // `rename` cannot name an open descriptor directly. First create a fresh link from
            // the held descriptor, so replacing the original predictable staging name can never
            // choose the inode installed over the source archive.
            let commit_path = temporary_stamp_commit_path(&out)?;
            let outcome = install_open_file_no_clobber(
                &output_file,
                &tmp,
                &commit_path,
                Durability::FileAndDirectory,
            )?;
            for warning in outcome.warnings {
                eprintln!("warning: {warning}");
            }
            let commit_metadata = std::fs::symlink_metadata(&commit_path).with_context(|| {
                format!("checking held stamp commit link {}", commit_path.display())
            })?;
            if !archive_inode_matches(&output_file.metadata()?, &commit_metadata) {
                bail!(
                    "stamp commit path {} was replaced before installation; preserving it",
                    commit_path.display()
                );
            }
            if let Err(error) = std::fs::rename(&commit_path, &out) {
                warn_on_failed_staging_cleanup(&commit_path, &output_file);
                return Err(error).with_context(|| {
                    format!("atomically replacing stamped archive {}", out.display())
                });
            }
            let installed = std::fs::symlink_metadata(&out)
                .with_context(|| format!("checking installed stamp output {}", out.display()))?;
            if !archive_inode_matches(&output_file.metadata()?, &installed) {
                bail!(
                    "installed stamp output is not the produced inode; preserving {} because a \
                     concurrent process may have replaced that path",
                    out.display()
                );
            }
        }
        Ok((output_identity, sections))
    })();
    if stamp_result.is_err() {
        warn_on_failed_staging_cleanup(&tmp, &staging_guard);
    }
    let (output_identity, sections) = stamp_result?;
    if reporting {
        emit_stamp_report(
            &args,
            true,
            source_identity
                .as_ref()
                .expect("reporting stamp captured source identity"),
            output_identity
                .as_ref()
                .expect("reporting stamp captured output identity"),
            &genome_identity,
            &sig,
            sections
                .as_ref()
                .expect("reporting stamp captured output sections"),
        )?;
    } else {
        println!(
            "stamped {} with genome signature ({} contigs, digest {}) in {:.1}s",
            out.display(),
            sig.contigs.len(),
            &sig.digest[..16],
            t0.elapsed().as_secs_f32()
        );
    }
    Ok(())
}

#[cfg(all(test, unix))]
mod archive_input_stability_tests {
    use super::{read_em_groups, validate_archive_input, SectionReader, SectionWriter};

    #[test]
    fn strict_em_group_barcodes_cannot_alias_a_short_packed_value() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gravlax-em-group-short-barcode-{}-{nonce}.tsv",
            std::process::id()
        ));
        std::fs::write(&path, "A\tgroup\n").unwrap();
        // `umi::pack("A")` and `umi::pack("AAAAAAAAAAAAAAAA")` have the same integer value;
        // external cell scopes must validate the barcode width before packing.
        let cells = [evidence_io::umi::pack(b"AAAAAAAAAAAAAAAA").unwrap()];
        let error = match read_em_groups(&path, &cells, false) {
            Ok(_) => panic!("short barcode unexpectedly entered the EM group scope"),
            Err(error) => error,
        };
        assert!(error
            .to_string()
            .contains("barcode must contain exactly 16 A/C/G/T bases"));
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn held_archive_reader_rejects_path_replacement() {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "gravlax-archive-input-stability-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&directory).unwrap();
        let path = directory.join("source.aie");
        let held = directory.join("held.aie");
        let mut writer = SectionWriter::create(&path, 1).unwrap();
        writer.section("meta", b"{}").unwrap();
        writer.section("chroms", b"chr1\n").unwrap();
        writer.finish().unwrap();

        let reader = SectionReader::open(&path).unwrap();
        let before = reader.file_metadata().unwrap();
        std::fs::rename(&path, &held).unwrap();
        std::fs::copy(&held, &path).unwrap();
        assert!(validate_archive_input(&reader, &before, &path).is_err());

        std::fs::remove_dir_all(directory).ok();
    }
}
