//! Standards-compliant BAM export of the post-correction archive abstraction.
//!
//! The archive no longer contains nucleotide UMI values: it contains opaque UMI-class ids plus
//! the 1MM adjacency graph consumed by replay. Pretending those ids are `UB` strings would be a
//! lossy and sometimes impossible embedding: an evaluated PBMC archive contains more global
//! classes than the 4^12 possible 12-base UMI strings. The export therefore uses local `X?` SAM
//! tags and explicit edge records while retaining ordinary BAM alignment fields for every
//! representative placement. This is a neutral-container baseline; generic BAM tools do not
//! interpret Gravlax's molecule semantics.

use crate::archivecmd::{read_archive, read_archive_with_identity, ArchiveContentIdentity};
use crate::build::to_ops;
use crate::rows::{
    identity_of_consumed_file, ConsumedFileIdentity, Extracted, MolChain, MolRec, PatAlt,
    SAME_SHAPE,
};
use anyhow::{bail, Context, Result};
use clap::ValueEnum;
use evidence_io::archive::Shape;
use evidence_io::umi;
use gravlax_output::{
    canonical_destination_key, publish_file_no_clobber, reported_output_path, DataType,
    Durability, Field, OutputError, OutputFormat, Producer, Provenance, ResultContext,
    RowSemantics, SelectionSummary, StreamingBundleWriter, TableSchema, TableSemantics,
};
use ingest::placement_from_alignment;
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam as sam;
use rustc_hash::FxHashMap;
use sam::alignment::{
    io::Write,
    record::{
        cigar::{op::Kind, Op},
        data::field::{Tag, Value as RecordValue},
        Flags, MappingQuality,
    },
    record_buf::{data::field::Value, Cigar, Data, RecordBuf},
};
use sam::header::record::value::{map::ReferenceSequence, Map};
use serde::Serialize;
use serde_json::json;
use smallvec::SmallVec;
use std::collections::BTreeMap;
use std::cell::RefCell;
use std::io::Write as IoWrite;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::rc::Rc;

const TAG_CELL: Tag = Tag::new(b'X', b'C');
const TAG_CLASS: Tag = Tag::new(b'X', b'I');
const TAG_MOLECULE: Tag = Tag::new(b'X', b'M');
const TAG_WEIGHT: Tag = Tag::new(b'X', b'W');
const TAG_KIND: Tag = Tag::new(b'X', b'K');
const TAG_GROUP: Tag = Tag::new(b'X', b'G');
const TAG_ALT: Tag = Tag::new(b'X', b'A');
const TAG_ANCHOR: Tag = Tag::new(b'X', b'P');
const TAG_EDGE_B: Tag = Tag::new(b'X', b'J');
const TAG_CB: Tag = Tag::new(b'C', b'B');

#[derive(clap::Parser)]
pub struct ExportArgs {
    /// Input `.aie` archive.
    pub archive: PathBuf,
    /// FASTA index supplying standards-compliant @SQ lengths.
    #[arg(long)]
    pub fai: PathBuf,
    /// Sequence-free BAM carrying the exact post-correction molecule abstraction.
    #[arg(long)]
    pub out: PathBuf,
    /// Emit a versioned export report while preserving the BAM artifact at --out.
    #[arg(long, value_enum)]
    pub report_format: Option<ExportReportFormat>,
    /// Publish the uniform report atomically without replacing an existing file.
    #[arg(long, requires = "report_format")]
    pub report_output: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum ExportReportFormat {
    Text,
    Tsv,
    Json,
}

impl From<ExportReportFormat> for OutputFormat {
    fn from(value: ExportReportFormat) -> Self {
        match value {
            ExportReportFormat::Text => Self::Text,
            ExportReportFormat::Tsv => Self::Tsv,
            ExportReportFormat::Json => Self::Json,
        }
    }
}

#[derive(Serialize)]
struct ExportSummary {
    molecules: u64,
    umi_edges: u64,
    bam_records: u64,
    output_bytes: u64,
    archive_format_version: u32,
    archive_bytes: u64,
}

fn artifact_schema() -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(
        "gravlax.molecule-bam.export.artifacts.v1",
        vec![
            Field::new("artifact_kind", DataType::String),
            Field::new("path", DataType::String),
            Field::new("bytes", DataType::UInt64),
            Field::new("records", DataType::UInt64),
            Field::new("identity", DataType::String),
        ],
    )?
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["artifact_kind", "path"]))
}

fn write_export_report<W: IoWrite>(
    writer: W,
    format: OutputFormat,
    archive_path: &str,
    fai_path: &str,
    fai_digest: &str,
    output_path: &str,
    output_digest: &str,
    archive_identity: &ArchiveContentIdentity,
    summary: &ExportSummary,
) -> std::result::Result<(), OutputError> {
    let mut parameters = BTreeMap::new();
    parameters.insert("archive_path".into(), json!(archive_path));
    parameters.insert("fasta_index".into(), json!(fai_path));
    parameters.insert("fasta_index_digest".into(), json!(fai_digest));
    parameters.insert("output".into(), json!(output_path));
    let context = ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            archives: vec![archive_identity.provenance_identity()],
            parameters,
            ..Default::default()
        },
        warnings: Vec::new(),
    };
    let schema = artifact_schema()?;
    let selection = SelectionSummary::complete(1);
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        "gravlax.molecule-bam.export.result.v1",
        format,
        &context,
        summary,
    )?;
    bundle.write_table("artifacts", &schema, Some(&selection), |rows| {
        rows.write_row_with(|row| {
            row.string("sequence-free-bam")?;
            row.string(output_path)?;
            row.uint64(summary.output_bytes)?;
            row.uint64(summary.bam_records)?;
            row.string(output_digest)?;
            Ok(())
        })
    })?;
    bundle.finish()?;
    Ok(())
}

fn parse_fai(path: &Path, text: &str) -> Result<FxHashMap<String, usize>> {
    let mut out = FxHashMap::default();
    for (i, line) in text.lines().enumerate() {
        let mut fields = line.split('\t');
        let name = fields
            .next()
            .filter(|s| !s.is_empty())
            .with_context(|| format!("{}:{}: missing contig name", path.display(), i + 1))?;
        let len: usize = fields
            .next()
            .with_context(|| format!("{}:{}: missing contig length", path.display(), i + 1))?
            .parse()
            .with_context(|| format!("{}:{}: invalid contig length", path.display(), i + 1))?;
        if len == 0 {
            bail!("{}:{}: zero-length contig", path.display(), i + 1);
        }
        if out.insert(name.to_string(), len).is_some() {
            bail!("{}:{}: duplicate contig {name}", path.display(), i + 1);
        }
    }
    Ok(out)
}

fn read_fai(path: &Path) -> Result<FxHashMap<String, usize>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading FASTA index {}", path.display()))?;
    parse_fai(path, &text)
}

fn read_fai_with_digest(path: &Path) -> Result<(FxHashMap<String, usize>, String)> {
    let bytes = std::fs::read(path)
        .with_context(|| format!("reading FASTA index {}", path.display()))?;
    let digest = format!("blake3:{}", blake3::hash(&bytes).to_hex());
    let text = std::str::from_utf8(&bytes)
        .with_context(|| format!("FASTA index {} is not UTF-8", path.display()))?;
    Ok((parse_fai(path, text)?, digest))
}

struct DigestingWriter<W> {
    inner: W,
    hasher: Rc<RefCell<blake3::Hasher>>,
}

impl<W> DigestingWriter<W> {
    fn new(inner: W) -> (Self, Rc<RefCell<blake3::Hasher>>) {
        let hasher = Rc::new(RefCell::new(blake3::Hasher::new()));
        (
            Self {
                inner,
                hasher: Rc::clone(&hasher),
            },
            hasher,
        )
    }
}

impl<W: IoWrite> IoWrite for DigestingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buf)?;
        self.hasher.borrow_mut().update(&buf[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

fn shape_cigar(shape: &Shape) -> Result<Cigar> {
    if shape.blocks.is_empty() {
        bail!("cannot export an empty placement shape");
    }
    let mut ops = Vec::with_capacity(shape.blocks.len() * 2 - 1);
    let mut prev_end = 0u32;
    for (i, &(off, len)) in shape.blocks.iter().enumerate() {
        if len == 0 {
            bail!("cannot export a zero-length alignment block");
        }
        if i == 0 {
            if off != 0 {
                bail!("first shape block starts at nonzero offset {off}");
            }
        } else {
            let gap = off
                .checked_sub(prev_end)
                .context("overlapping or reversed shape blocks")?;
            if gap == 0 {
                bail!("adjacent shape blocks must be coalesced");
            }
            ops.push(Op::new(Kind::Skip, gap as usize));
        }
        ops.push(Op::new(Kind::Match, len as usize));
        prev_end = off.checked_add(len).context("shape end overflow")?;
    }
    Ok(ops.into_iter().collect())
}

#[derive(Clone, Copy)]
struct PlacementTags {
    molecule: u32,
    cell: u32,
    class: u32,
    weight: u32,
    kind: u8,
    group: u32,
    alt: u32,
    anchor: bool,
    nh: u32,
}

fn placement_record(
    x: &Extracted,
    chrom: u32,
    pos: u32,
    strand_rev: bool,
    shape: u32,
    tags: PlacementTags,
) -> Result<RecordBuf> {
    let shape_def = x
        .shapes
        .get(shape as usize)
        .context("shape id out of bounds")?;
    let cigar = shape_cigar(shape_def)?;
    let start = usize::try_from(pos)
        .context("alignment start exceeds usize")?
        .checked_add(1)
        .context("one-based alignment start overflow")?;
    let barcode = x
        .cells
        .get(tags.cell as usize)
        .context("cell id out of bounds")?;
    let barcode = String::from_utf8(umi::unpack(*barcode, 16)).unwrap();
    let data: Data = [
        (TAG_CB, Value::from(barcode)),
        (TAG_CELL, Value::UInt32(tags.cell)),
        (TAG_CLASS, Value::UInt32(tags.class)),
        (TAG_MOLECULE, Value::UInt32(tags.molecule)),
        (TAG_WEIGHT, Value::UInt32(tags.weight)),
        (TAG_KIND, Value::Character(tags.kind)),
        (TAG_GROUP, Value::UInt32(tags.group)),
        (TAG_ALT, Value::UInt32(tags.alt)),
        (TAG_ANCHOR, Value::UInt8(tags.anchor as u8)),
        (Tag::ALIGNMENT_HIT_COUNT, Value::UInt32(tags.nh)),
    ]
    .into_iter()
    .collect();
    let mut flags = if strand_rev {
        Flags::REVERSE_COMPLEMENTED
    } else {
        Flags::empty()
    };
    if tags.kind == b'M' && !tags.anchor {
        flags |= Flags::SECONDARY;
    }
    let name = format!("gx{:x}", tags.molecule);
    Ok(RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(chrom as usize)
        .set_alignment_start(Position::try_from(start).context("invalid alignment start")?)
        .set_cigar(cigar)
        .set_data(data)
        .build())
}

fn edge_record(a: u32, b: u32) -> RecordBuf {
    let data: Data = [
        (TAG_KIND, Value::Character(b'E')),
        (TAG_CLASS, Value::UInt32(a)),
        (TAG_EDGE_B, Value::UInt32(b)),
    ]
    .into_iter()
    .collect();
    let name = format!("gxe{a:x}_{b:x}");
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::UNMAPPED)
        .set_mapping_quality(MappingQuality::MIN)
        .set_data(data)
        .build()
}

fn write_molecule_bam<W: IoWrite>(
    x: &Extracted,
    fai: &FxHashMap<String, usize>,
    sink: W,
) -> Result<u64> {
    let mut hb = sam::Header::builder();
    for name in &x.chrom_names {
        let len = *fai
            .get(name)
            .with_context(|| format!("contig {name} is absent from the FASTA index"))?;
        hb = hb.add_reference_sequence(
            name.as_str(),
            Map::<ReferenceSequence>::new(NonZero::new(len).unwrap()),
        );
    }
    let header = hb.build();
    let mut writer = bam::io::Writer::new(sink);
    writer.write_header(&header)?;
    let mut n_records = 0u64;
    for (mi, mol) in x.mols.iter().enumerate() {
        let molecule = u32::try_from(mi).context("more than u32::MAX molecules")?;
        for (gi, chain) in mol.chains.iter().enumerate() {
            for (ai, &(pos, shape)) in chain.reps.iter().enumerate() {
                let record = placement_record(
                    x,
                    mol.chrom,
                    pos,
                    mol.strand_rev,
                    shape,
                    PlacementTags {
                        molecule,
                        cell: mol.cell,
                        class: mol.umi_class,
                        weight: chain.weight,
                        kind: b'C',
                        group: gi as u32,
                        alt: ai as u32,
                        anchor: ai == 0,
                        nh: 1,
                    },
                )?;
                writer.write_alignment_record(&header, &record)?;
                n_records += 1;
            }
        }
        for (gi, &(anchor_pos, anchor_shape, pattern, weight)) in mol.mms.iter().enumerate() {
            let pdef = x
                .patterns
                .get(pattern as usize)
                .context("pattern id out of bounds")?;
            let anchor_i = pdef
                .iter()
                .enumerate()
                .find_map(|(i, alt)| {
                    let shape = if alt.shape == SAME_SHAPE {
                        anchor_shape
                    } else {
                        alt.shape
                    };
                    (alt.chrom == mol.chrom
                        && alt.offset == 0
                        && !alt.strand_flip
                        && shape == anchor_shape)
                        .then_some(i)
                })
                .with_context(|| {
                    format!("molecule {mi} multimapper group {gi} has no anchor alternative")
                })?;
            for (ai, alt) in pdef.iter().enumerate() {
                let apos = i64::from(anchor_pos)
                    .checked_add(alt.offset)
                    .context("multimapper position overflow")?;
                let pos = u32::try_from(apos).context("multimapper position outside u32")?;
                let shape = if alt.shape == SAME_SHAPE {
                    anchor_shape
                } else {
                    alt.shape
                };
                let record = placement_record(
                    x,
                    alt.chrom,
                    pos,
                    mol.strand_rev != alt.strand_flip,
                    shape,
                    PlacementTags {
                        molecule,
                        cell: mol.cell,
                        class: mol.umi_class,
                        weight,
                        kind: b'M',
                        group: gi as u32,
                        alt: ai as u32,
                        anchor: ai == anchor_i,
                        nh: pdef.len() as u32,
                    },
                )?;
                writer.write_alignment_record(&header, &record)?;
                n_records += 1;
            }
        }
    }
    for &(a, b) in &x.edges {
        writer.write_alignment_record(&header, &edge_record(a, b))?;
        n_records += 1;
    }
    writer.try_finish()?;
    // Preserve the historical writer lifecycle (and therefore artifact bytes) exactly.
    drop(writer);
    Ok(n_records)
}

pub fn run_export(args: ExportArgs) -> Result<()> {
    if let Some(path) = &args.report_output {
        if std::fs::symlink_metadata(path).is_ok() {
            bail!("refusing to overwrite {}", path.display());
        }
        if canonical_destination_key(path)? == canonical_destination_key(&args.out)? {
            bail!("BAM output and report output must differ");
        }
    }
    let uniform_paths = args
        .report_format
        .map(|_| {
            Ok::<_, anyhow::Error>((
                args.archive
                    .to_str()
                    .context("uniform report requires a UTF-8 archive path")?
                    .to_owned(),
                args.fai
                    .to_str()
                    .context("uniform report requires a UTF-8 FASTA-index path")?
                    .to_owned(),
                reported_output_path(&args.out)?,
            ))
        })
        .transpose()?;
    let (x, archive_identity) = if args.report_format.is_some() {
        let (x, identity) = read_archive_with_identity(&args.archive)?;
        (x, Some(identity))
    } else {
        (read_archive(&args.archive)?, None)
    };
    let (fai, fai_digest) = if args.report_format.is_some() {
        let (fai, digest) = read_fai_with_digest(&args.fai)?;
        (fai, Some(digest))
    } else {
        (read_fai(&args.fai)?, None)
    };
    let file = std::fs::File::create(&args.out)
        .with_context(|| format!("creating {}", args.out.display()))?;
    let (n_records, output_digest) = if args.report_format.is_some() {
        let (sink, hasher) = DigestingWriter::new(file);
        let n_records = write_molecule_bam(&x, &fai, sink)?;
        let output_digest = format!("blake3:{}", hasher.borrow().finalize().to_hex());
        (n_records, Some(output_digest))
    } else {
        (write_molecule_bam(&x, &fai, file)?, None)
    };
    let output_bytes = std::fs::metadata(&args.out)?.len();
    if let (Some(format), Some((archive_path, fai_path, output_path)), Some(identity)) =
        (args.report_format, uniform_paths, archive_identity.as_ref())
    {
        let summary = ExportSummary {
            molecules: x.mols.len() as u64,
            umi_edges: x.edges.len() as u64,
            bam_records: n_records,
            output_bytes,
            archive_format_version: identity.format_version,
            archive_bytes: identity.file_bytes,
        };
        match &args.report_output {
            Some(path) => {
                let outcome = publish_file_no_clobber(path, Durability::Flush, |writer| {
                    write_export_report(
                        writer,
                        format.into(),
                        &archive_path,
                        &fai_path,
                        fai_digest.as_deref().expect("uniform FASTA-index identity"),
                        &output_path,
                        output_digest.as_deref().expect("uniform BAM identity"),
                        identity,
                        &summary,
                    )
                })?;
                for warning in outcome.warnings {
                    eprintln!("warning: {warning}");
                }
            }
            None => write_export_report(
                std::io::stdout().lock(),
                format.into(),
                &archive_path,
                &fai_path,
                fai_digest.as_deref().expect("uniform FASTA-index identity"),
                &output_path,
                output_digest.as_deref().expect("uniform BAM identity"),
                identity,
                &summary,
            )?,
        }
        eprintln!(
            "exported {} molecules and {} 1MM edges as {n_records} BAM records",
            x.mols.len(),
            x.edges.len(),
        );
    } else {
        println!(
            "exported {} molecules and {} 1MM edges as {n_records} BAM records -> {} ({} bytes)",
            x.mols.len(),
            x.edges.len(),
            args.out.display(),
            output_bytes
        );
    }
    Ok(())
}

fn uint_tag(value: RecordValue<'_>, name: &str) -> Result<u32> {
    let n = value
        .as_int()
        .with_context(|| format!("{name} must be an integer"))?;
    u32::try_from(n).with_context(|| format!("{name} is outside u32"))
}

#[derive(Debug)]
struct RawPlacement {
    kind: u8,
    group: u32,
    alt: u32,
    anchor: bool,
    chrom: u32,
    pos: u32,
    strand_rev: bool,
    shape: u32,
    weight: u32,
    nh: u32,
    secondary: bool,
}

#[derive(Debug)]
struct RawMol {
    id: u32,
    cell: u32,
    class: u32,
    entries: Vec<RawPlacement>,
}

fn set_parent(parent: &mut Option<(u32, bool)>, value: (u32, bool)) -> Result<()> {
    match parent {
        Some(old) if *old != value => bail!("molecule groups disagree on anchor chromosome/strand"),
        Some(_) => {}
        None => *parent = Some(value),
    }
    Ok(())
}

fn finalize_molecule(
    raw: RawMol,
    patterns: &mut Vec<Vec<PatAlt>>,
    pattern_intern: &mut FxHashMap<Vec<PatAlt>, u32>,
) -> Result<MolRec> {
    let mut chains: BTreeMap<u32, Vec<RawPlacement>> = BTreeMap::new();
    let mut mms: BTreeMap<u32, Vec<RawPlacement>> = BTreeMap::new();
    for entry in raw.entries {
        match entry.kind {
            b'C' => chains.entry(entry.group).or_default().push(entry),
            b'M' => mms.entry(entry.group).or_default().push(entry),
            k => bail!(
                "molecule {} has invalid placement kind {}",
                raw.id,
                char::from(k)
            ),
        }
    }
    let mut parent = None;
    let mut out_chains: SmallVec<[MolChain; 1]> = SmallVec::new();
    for (expected_group, (group, mut entries)) in chains.into_iter().enumerate() {
        if group != expected_group as u32 {
            bail!("molecule {} has noncontiguous chain group {group}", raw.id);
        }
        entries.sort_unstable_by_key(|e| e.alt);
        let weight = entries.first().context("empty chain group")?.weight;
        let mut reps = SmallVec::new();
        for (i, entry) in entries.into_iter().enumerate() {
            if entry.alt != i as u32 || entry.weight != weight {
                bail!(
                    "molecule {} chain {group} has inconsistent indices or weights",
                    raw.id
                );
            }
            if entry.nh != 1 || entry.secondary {
                bail!(
                    "molecule {} chain {group} has invalid NH/secondary flags",
                    raw.id
                );
            }
            set_parent(&mut parent, (entry.chrom, entry.strand_rev))?;
            reps.push((entry.pos, entry.shape));
        }
        if reps.is_empty() || reps.len() > 2 {
            bail!(
                "molecule {} chain {group} has {} representatives",
                raw.id,
                reps.len()
            );
        }
        out_chains.push(MolChain { weight, reps });
    }
    let mut out_mms: SmallVec<[(u32, u32, u32, u32); 1]> = SmallVec::new();
    for (expected_group, (group, mut entries)) in mms.into_iter().enumerate() {
        if group != expected_group as u32 {
            bail!(
                "molecule {} has noncontiguous multimapper group {group}",
                raw.id
            );
        }
        entries.sort_unstable_by_key(|e| e.alt);
        let weight = entries.first().context("empty multimapper group")?.weight;
        let nh = u32::try_from(entries.len()).context("multimapper group exceeds u32")?;
        let mut anchor = None;
        for (i, entry) in entries.iter().enumerate() {
            if entry.alt != i as u32 || entry.weight != weight {
                bail!(
                    "molecule {} multimapper {group} has inconsistent indices or weights",
                    raw.id
                );
            }
            if entry.nh != nh || entry.secondary == entry.anchor {
                bail!(
                    "molecule {} multimapper {group} has invalid NH/anchor/secondary flags",
                    raw.id
                );
            }
            if entry.anchor {
                if anchor
                    .replace((entry.chrom, entry.pos, entry.strand_rev, entry.shape))
                    .is_some()
                {
                    bail!(
                        "molecule {} multimapper {group} has multiple anchors",
                        raw.id
                    );
                }
            }
        }
        let (achrom, apos, arev, ashape) = anchor
            .with_context(|| format!("molecule {} multimapper {group} has no anchor", raw.id))?;
        set_parent(&mut parent, (achrom, arev))?;
        let pdef: Vec<PatAlt> = entries
            .iter()
            .map(|entry| PatAlt {
                chrom: entry.chrom,
                offset: i64::from(entry.pos) - i64::from(apos),
                strand_flip: entry.strand_rev != arev,
                shape: if entry.shape == ashape {
                    SAME_SHAPE
                } else {
                    entry.shape
                },
            })
            .collect();
        let pattern = match pattern_intern.get(&pdef) {
            Some(&id) => id,
            None => {
                let id = u32::try_from(patterns.len()).context("more than u32::MAX patterns")?;
                pattern_intern.insert(pdef.clone(), id);
                patterns.push(pdef);
                id
            }
        };
        out_mms.push((apos, ashape, pattern, weight));
    }
    let (chrom, strand_rev) = parent.context("molecule has no placement records")?;
    Ok(MolRec {
        cell: raw.cell,
        umi_class: raw.class,
        chrom,
        strand_rev,
        chains: out_chains,
        mms: out_mms,
    })
}

pub fn read_molecule_bam(path: &PathBuf) -> Result<Extracted> {
    Ok(read_molecule_bam_inner(path, false)?.0)
}

pub(crate) fn read_molecule_bam_with_identity(
    path: &PathBuf,
) -> Result<(Extracted, ConsumedFileIdentity)> {
    let (extracted, identity) = read_molecule_bam_inner(path, true)?;
    Ok((
        extracted,
        identity.expect("reporting molecule-BAM read requested an identity"),
    ))
}

fn read_molecule_bam_inner(
    path: &PathBuf,
    capture_identity: bool,
) -> Result<(Extracted, Option<ConsumedFileIdentity>)> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let source_before = capture_identity.then(|| file.metadata()).transpose()?;
    if source_before.as_ref().is_some_and(|metadata| !metadata.is_file()) {
        bail!("molecule-BAM input is not a regular file: {}", path.display());
    }
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header()?;
    let chrom_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .collect();
    let mut rec = bam::Record::default();
    let mut cells: Vec<Option<u32>> = Vec::new();
    let mut class_cells: Vec<Option<u32>> = Vec::new();
    let mut shape_intern: FxHashMap<Shape, u32> = FxHashMap::default();
    let mut shapes = Vec::new();
    let mut pattern_intern: FxHashMap<Vec<PatAlt>, u32> = FxHashMap::default();
    let mut patterns = Vec::new();
    let mut mols = Vec::new();
    let mut edges = Vec::new();
    let mut current: Option<RawMol> = None;
    let mut edge_phase = false;

    while reader.read_record(&mut rec)? != 0 {
        let mut kind = None;
        let mut cell = None;
        let mut class = None;
        let mut molecule = None;
        let mut weight = None;
        let mut group = None;
        let mut alt = None;
        let mut anchor = false;
        let mut edge_b = None;
        let mut packed_cb = None;
        let mut nh = None;
        for field in rec.data().iter() {
            let (tag, value) = field?;
            match tag {
                TAG_KIND => {
                    kind = match value {
                        RecordValue::Character(v) => Some(v),
                        _ => bail!("XK must be a character"),
                    }
                }
                TAG_CELL => cell = Some(uint_tag(value, "XC")?),
                TAG_CLASS => class = Some(uint_tag(value, "XI")?),
                TAG_MOLECULE => molecule = Some(uint_tag(value, "XM")?),
                TAG_WEIGHT => weight = Some(uint_tag(value, "XW")?),
                TAG_GROUP => group = Some(uint_tag(value, "XG")?),
                TAG_ALT => alt = Some(uint_tag(value, "XA")?),
                TAG_ANCHOR => anchor = uint_tag(value, "XP")? != 0,
                TAG_EDGE_B => edge_b = Some(uint_tag(value, "XJ")?),
                TAG_CB => {
                    packed_cb = match value {
                        RecordValue::String(v) if v.len() == 16 => umi::pack(v.as_ref()),
                        _ => bail!("CB must be a 16-base string"),
                    }
                }
                Tag::ALIGNMENT_HIT_COUNT => nh = Some(uint_tag(value, "NH")?),
                _ => {}
            }
        }
        let kind = kind.context("record is missing XK")?;
        if kind == b'E' {
            edge_phase = true;
            if let Some(raw) = current.take() {
                mols.push(finalize_molecule(raw, &mut patterns, &mut pattern_intern)?);
            }
            if !Flags::from(rec.flags().bits()).is_unmapped() {
                bail!("edge record must be unmapped");
            }
            let a = class.context("edge record is missing XI")?;
            let b = edge_b.context("edge record is missing XJ")?;
            if a >= b {
                bail!("edge endpoints must satisfy XI < XJ");
            }
            edges.push((a, b));
            continue;
        }
        if edge_phase {
            bail!("placement record appears after edge records");
        }
        if kind != b'C' && kind != b'M' {
            bail!("invalid XK placement kind {}", char::from(kind));
        }
        let cell = cell.context("placement record is missing XC")?;
        let class = class.context("placement record is missing XI")?;
        let molecule = molecule.context("placement record is missing XM")?;
        let weight = weight.context("placement record is missing XW")?;
        if weight == 0 {
            bail!("placement record has zero XW");
        }
        let packed_cb = packed_cb.context("placement record has invalid or missing CB")?;
        if cells.len() <= cell as usize {
            cells.resize(cell as usize + 1, None);
        }
        match cells[cell as usize] {
            Some(old) if old != packed_cb => bail!("cell {cell} maps to multiple CB values"),
            Some(_) => {}
            None => cells[cell as usize] = Some(packed_cb),
        }
        if class_cells.len() <= class as usize {
            class_cells.resize(class as usize + 1, None);
        }
        match class_cells[class as usize] {
            Some(old) if old != cell => bail!("class {class} maps to multiple cells"),
            Some(_) => {}
            None => class_cells[class as usize] = Some(cell),
        }
        let flags = Flags::from(rec.flags().bits());
        if flags.is_unmapped() || flags.is_supplementary() {
            bail!("placement record is unmapped or supplementary");
        }
        let chrom = u32::try_from(
            rec.reference_sequence_id()
                .transpose()?
                .context("placement record has no reference")?,
        )
        .context("reference id exceeds u32")?;
        if chrom as usize >= chrom_names.len() {
            bail!("placement reference id is outside the header");
        }
        let pos = rec
            .alignment_start()
            .transpose()?
            .context("placement has no alignment start")?;
        let pos = u32::try_from(usize::from(pos) - 1).context("alignment start exceeds u32")?;
        if rec.sequence().len() != 0 || rec.quality_scores().len() != 0 {
            bail!("post-correction placement records must be sequence-free");
        }
        for op in rec.cigar().iter() {
            if !matches!(op?.kind(), Kind::Match | Kind::Skip) {
                bail!("post-correction placement CIGAR may contain only M and N operations");
            }
        }
        let pl = placement_from_alignment(
            chrom,
            pos,
            flags.is_reverse_complemented(),
            &to_ops(rec.cigar().iter())?,
            0,
            0,
            1,
        );
        let shape_def = Shape::of(&pl);
        let shape = match shape_intern.get(&shape_def) {
            Some(&id) => id,
            None => {
                let id = u32::try_from(shapes.len()).context("more than u32::MAX shapes")?;
                shape_intern.insert(shape_def.clone(), id);
                shapes.push(shape_def);
                id
            }
        };
        if current.as_ref().is_some_and(|m| m.id != molecule) {
            let raw = current.take().unwrap();
            if molecule != raw.id.checked_add(1).context("molecule id overflow")? {
                bail!(
                    "molecule ids are not contiguous: {} followed by {molecule}",
                    raw.id
                );
            }
            mols.push(finalize_molecule(raw, &mut patterns, &mut pattern_intern)?);
        }
        if current.is_none() && mols.is_empty() && molecule != 0 {
            bail!("first molecule id must be 0, found {molecule}");
        }
        let raw = current.get_or_insert_with(|| RawMol {
            id: molecule,
            cell,
            class,
            entries: Vec::new(),
        });
        if raw.cell != cell || raw.class != class {
            bail!("molecule {molecule} has inconsistent cell/class tags");
        }
        raw.entries.push(RawPlacement {
            kind,
            group: group.context("placement record is missing XG")?,
            alt: alt.context("placement record is missing XA")?,
            anchor,
            chrom,
            pos,
            strand_rev: flags.is_reverse_complemented(),
            shape,
            weight,
            nh: nh.context("placement record is missing NH")?,
            secondary: flags.is_secondary(),
        });
    }
    let identity_thread = source_before.map(|before| {
        let file = reader.into_inner().into_inner();
        let path = path.clone();
        std::thread::spawn(move || identity_of_consumed_file(file, before, &path))
    });
    if let Some(raw) = current.take() {
        mols.push(finalize_molecule(raw, &mut patterns, &mut pattern_intern)?);
    }
    let cells: Vec<u32> = cells
        .into_iter()
        .enumerate()
        .map(|(i, v)| v.with_context(|| format!("missing dense cell id {i}")))
        .collect::<Result<_>>()?;
    for (i, cell) in class_cells.iter().enumerate() {
        if cell.is_none() {
            bail!("missing dense class id {i}");
        }
    }
    for (i, &(a, b)) in edges.iter().enumerate() {
        if i > 0 && edges[i - 1] >= (a, b) {
            bail!("edge records must be strictly sorted and unique");
        }
        if b as usize >= class_cells.len() {
            bail!("edge ({a},{b}) exceeds class count {}", class_cells.len());
        }
        if class_cells[a as usize] != class_cells[b as usize] {
            bail!("edge ({a},{b}) crosses cells");
        }
    }
    let extracted = Extracted {
        mols,
        edges,
        cells,
        shapes,
        patterns,
        n_classes: u32::try_from(class_cells.len()).context("more than u32::MAX classes")?,
        chrom_names,
    };
    let identity = identity_thread
        .map(|thread| {
            thread
                .join()
                .map_err(|_| anyhow::anyhow!("molecule-BAM identity thread panicked"))?
        })
        .transpose()?;
    Ok((extracted, identity))
}
