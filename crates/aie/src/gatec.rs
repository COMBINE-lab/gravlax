//! Compare assignment-relevant evidence from annotation-free and annotation-aware alignments.
//!
//! Both input BAMs are produced from the same reads in the same order, so this streams them in
//! parallel and joins on read name rather than building a hash of the whole file. If the order
//! ever diverges the join detects it and fails loudly rather than silently reporting nonsense.

use anyhow::{bail, Context, Result};
use clap::Parser;
use ingest::{classify, cigar::Op, placement_from_alignment, Divergence};
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::alignment::record::Flags;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::path::PathBuf;

/// The concrete reader type produced by noodles' BAM builder.
type BamReader = bam::io::Reader<noodles_bgzf::io::Reader<File>>;

#[derive(Parser)]
pub struct Args {
    /// Baseline BAM (config i: annotation-free ingest).
    pub a: PathBuf,
    /// Comparison BAM (config ii or iii: annotation-aware).
    pub b: PathBuf,
    /// Restrict the comparison to reads assigned to a gene in this STARsolo BAM, using its `GX`
    /// tag. Reads that no annotation claims cannot affect a gene count.
    ///
    /// Read names are stored as 64-bit hashes (~270 MB for 34M reads rather than several GB); the
    /// collision probability over this many names is negligible.
    #[arg(long)]
    pub gene_assigned_from: Option<PathBuf>,
    /// Write the per-category counts here as JSON.
    #[arg(long)]
    pub json_out: Option<PathBuf>,
}

/// The primary alignment of one read, or `None` if the read is unmapped in this configuration.
struct Primary {
    name: String,
    placement: Option<evidence_io::Placement>,
}

fn to_ops(cigar: impl Iterator<Item = std::io::Result<noodles_sam::alignment::record::cigar::Op>>) -> Result<Vec<Op>> {
    let mut ops = Vec::new();
    for op in cigar {
        let op = op?;
        let n = op.len() as u32;
        ops.push(match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => Op::Match(n),
            Kind::Insertion => Op::Ins(n),
            Kind::Deletion => Op::Del(n),
            Kind::Skip => Op::Skip(n),
            Kind::SoftClip => Op::SoftClip(n),
            Kind::HardClip => Op::HardClip(n),
            Kind::Pad => Op::Pad(n),
        });
    }
    Ok(ops)
}

fn int_tag(rec: &bam::Record, tag: Tag) -> Option<i64> {
    match rec.data().get(&tag)? {
        Ok(Value::Int8(v)) => Some(v as i64),
        Ok(Value::UInt8(v)) => Some(v as i64),
        Ok(Value::Int16(v)) => Some(v as i64),
        Ok(Value::UInt16(v)) => Some(v as i64),
        Ok(Value::Int32(v)) => Some(v as i64),
        Ok(Value::UInt32(v)) => Some(v as i64),
        _ => None,
    }
}

/// Read the next primary (non-secondary, non-supplementary) record.
fn next_primary(
    reader: &mut BamReader,
    rec: &mut bam::Record,
) -> Result<Option<Primary>> {
    loop {
        if reader.read_record(rec)? == 0 {
            return Ok(None);
        }
        let flags = Flags::from(rec.flags().bits());
        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let name = String::from_utf8_lossy(rec.name().map(|n| n.as_ref()).unwrap_or(b"")).into_owned();

        if flags.is_unmapped() {
            return Ok(Some(Primary { name, placement: None }));
        }

        let chrom = rec.reference_sequence_id().transpose()?.unwrap_or(0) as u32;
        let pos = rec
            .alignment_start()
            .transpose()?
            .map(|p| usize::from(p) as u32 - 1)
            .unwrap_or(0);
        let ops = to_ops(rec.cigar().iter())?;
        let nm = int_tag(rec, Tag::MISMATCHED_POSITIONS).unwrap_or(0) as u16;
        let score = int_tag(rec, Tag::ALIGNMENT_SCORE).unwrap_or(0) as i32;
        let nh = int_tag(rec, Tag::ALIGNMENT_HIT_COUNT).unwrap_or(1) as u16;

        return Ok(Some(Primary {
            name,
            placement: Some(placement_from_alignment(
                chrom,
                pos,
                flags.is_reverse_complemented(),
                &ops,
                nm,
                score,
                nh,
            )),
        }));
    }
}

fn open(path: &PathBuf) -> Result<BamReader> {
    let mut r = bam::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;
    r.read_header()?;
    Ok(r)
}

/// Hash a read name to 64 bits.
fn name_hash(name: &str) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = rustc_hash::FxHasher::default();
    name.hash(&mut h);
    h.finish()
}

/// Collect name hashes of reads the oracle assigned to exactly one gene.
fn gene_assigned_set(path: &PathBuf) -> Result<rustc_hash::FxHashSet<u64>> {
    let mut r = open(path)?;
    let mut rec = bam::Record::default();
    let mut set = rustc_hash::FxHashSet::default();
    while r.read_record(&mut rec)? != 0 {
        let flags = Flags::from(rec.flags().bits());
        if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped() {
            continue;
        }
        let mut assigned = false;
        for field in rec.data().iter() {
            let (tag, value) = field?;
            if <[u8; 2]>::from(tag) != *b"GX" {
                continue;
            }
            if let Value::String(s) = value {
                let g: &[u8] = s.as_ref();
                assigned = !g.is_empty() && g != b"-" && !g.contains(&b',');
            }
            break;
        }
        if assigned {
            let name = String::from_utf8_lossy(rec.name().map(|n| n.as_ref()).unwrap_or(b""));
            set.insert(name_hash(&name));
        }
    }
    Ok(set)
}

pub fn run(args: Args) -> Result<()> {
    let restrict: Option<rustc_hash::FxHashSet<u64>> = match &args.gene_assigned_from {
        Some(p) => {
            eprintln!("collecting gene-assigned reads from {}", p.display());
            let s = gene_assigned_set(p)?;
            eprintln!("  {} gene-assigned reads", s.len());
            Some(s)
        }
        None => None,
    };

    let mut ra = open(&args.a)?;
    let mut rb = open(&args.b)?;
    let (mut reca, mut recb) = (bam::Record::default(), bam::Record::default());

    let mut counts: FxHashMap<&'static str, u64> = FxHashMap::default();
    let mut total = 0u64;

    loop {
        let (pa, pb) = (next_primary(&mut ra, &mut reca)?, next_primary(&mut rb, &mut recb)?);
        let (pa, pb) = match (pa, pb) {
            (Some(a), Some(b)) => (a, b),
            (None, None) => break,
            _ => bail!("BAMs contain different numbers of primary records; they must come from the same reads"),
        };
        if pa.name != pb.name {
            bail!("read-name mismatch at record {total}: {:?} vs {:?} — inputs must be in the same order", pa.name, pb.name);
        }
        if let Some(set) = &restrict {
            if !set.contains(&name_hash(&pa.name)) {
                continue;
            }
        }
        total += 1;
        let d = classify(pa.placement.as_ref(), pb.placement.as_ref());
        *counts.entry(label(d)).or_insert(0) += 1;
    }

    let identical = *counts.get("identical").unwrap_or(&0);
    let pct = if total > 0 { 100.0 * identical as f64 / total as f64 } else { 0.0 };

    println!("reads compared:      {total}");
    println!("identical evidence:  {identical}  ({pct:.4}%)");
    let mut rest: Vec<_> = counts.iter().filter(|(k, _)| **k != "identical").collect();
    rest.sort_by_key(|(k, _)| *k);
    for (k, v) in rest {
        let p = if total > 0 { 100.0 * *v as f64 / total as f64 } else { 0.0 };
        println!("  {k:<14} {v:>12}  ({p:.4}%)");
    }

    if let Some(p) = &args.json_out {
        let mut obj = serde_json::Map::new();
        obj.insert("total".into(), serde_json::json!(total));
        obj.insert("identical_pct".into(), serde_json::json!(pct));
        for (k, v) in &counts {
            obj.insert((*k).into(), serde_json::json!(v));
        }
        std::fs::write(p, serde_json::to_string_pretty(&obj)?)?;
    }
    Ok(())
}

fn label(d: Divergence) -> &'static str {
    match d {
        Divergence::Identical => "identical",
        Divergence::Locus => "locus",
        Divergence::BlockBoundary => "block_boundary",
        Divergence::Junction => "junction",
        Divergence::Multiplicity => "multiplicity",
        Divergence::Presence => "presence",
    }
}
