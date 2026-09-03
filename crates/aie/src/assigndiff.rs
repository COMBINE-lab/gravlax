//! `aie assign-diff` — read-level confusion between STARsolo's gene assignment (the `GX` tag) and
//! our `anno::assign` port, on IDENTICAL alignments.
//!
//! In an evaluated PBMC comparison, the assignment port accounted for a 1.95% difference, compared
//! with 0.66 and 0.13 percentage points from barcode processing and alignment, respectively.
//! Matrix-level comparisons cannot identify the differing rule, so this command counts each
//! disagreement class and prints a bounded sample with full geometry for inspection against STAR's
//! source.

use anyhow::{Context, Result};
use clap::Parser;
use ingest::{cigar::Op, placement_from_alignment};
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Value;
use noodles_sam::alignment::record::Flags;
use rustc_hash::FxHashMap;
use std::path::PathBuf;

#[derive(Parser)]
pub struct Args {
    /// STARsolo reference BAM with GX tags and alignments identical to those being classified.
    pub bam: PathBuf,
    /// The same GTF used for the STARsolo reference run.
    #[arg(long)]
    pub gtf: PathBuf,
    /// Print up to this many examples per disagreement class.
    #[arg(long, default_value_t = 5)]
    pub examples: usize,
    #[arg(long)]
    pub json_out: Option<PathBuf>,
}

pub fn run(args: Args) -> Result<()> {
    let anno = anno::Annotation::from_path(&args.gtf)?;
    let gene_lookup: FxHashMap<&str, u32> = anno
        .gene_ids
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i as u32))
        .collect();

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&args.bam)
        .with_context(|| format!("opening {}", args.bam.display()))?;
    let header = reader.read_header()?;
    let chrom_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .collect();
    let bam2anno: Vec<Option<u32>> = chrom_names.iter().map(|n| anno.chrom_ids.get(n).copied()).collect();

    let mut counts: FxHashMap<&'static str, u64> = FxHashMap::default();
    let mut shown: FxHashMap<&'static str, usize> = FxHashMap::default();
    let mut rec = bam::Record::default();
    let mut n = 0u64;

    while reader.read_record(&mut rec)? != 0 {
        let flags = Flags::from(rec.flags().bits());
        if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped() {
            continue;
        }

        let (mut gx, mut nh) = (None::<String>, 1u16);
        for field in rec.data().iter() {
            let (tag, value) = field?;
            let key = <[u8; 2]>::from(tag);
            match (&key, value) {
                (b"GX", Value::String(s)) => gx = Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
                (b"NH", v) => {
                    nh = match v {
                        Value::Int8(x) => x as u16,
                        Value::UInt8(x) => x as u16,
                        Value::Int16(x) => x as u16,
                        Value::UInt16(x) => x as u16,
                        Value::Int32(x) => x as u16,
                        Value::UInt32(x) => x as u16,
                        _ => 1,
                    }
                }
                _ => {}
            }
        }
        if nh != 1 {
            continue; // STARsolo assigns genes only to unique alignments; so do we.
        }
        n += 1;

        let chrom_bam = rec.reference_sequence_id().transpose()?.unwrap_or(0) as usize;
        let pos = rec.alignment_start().transpose()?.map(|p| usize::from(p) as u32 - 1).unwrap_or(0);
        let mut ops = Vec::new();
        for op in rec.cigar().iter() {
            let op = op?;
            let l = op.len() as u32;
            ops.push(match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => Op::Match(l),
                Kind::Insertion => Op::Ins(l),
                Kind::Deletion => Op::Del(l),
                Kind::Skip => Op::Skip(l),
                Kind::SoftClip => Op::SoftClip(l),
                Kind::HardClip => Op::HardClip(l),
                Kind::Pad => Op::Pad(l),
            });
        }
        let p = placement_from_alignment(
            chrom_bam as u32, pos, flags.is_reverse_complemented(), &ops, 0, 0, nh,
        );

        let ours: Vec<u32> = match bam2anno.get(chrom_bam).copied().flatten() {
            Some(ac) => anno::assign::concordant_genes(&p, &anno, ac),
            None => Vec::new(),
        };

        // The oracle's verdict: a specific gene, "-" (none), or a comma list (multi → uncounted).
        let oracle: Option<Option<u32>> = match gx.as_deref() {
            Some("-") | None => Some(None),
            Some(g) if g.contains(',') => None, // multi-gene: STARsolo does not count it either
            Some(g) => Some(Some(*gene_lookup.get(g).unwrap_or(&u32::MAX))),
        };

        let class: &'static str = match (oracle, ours.as_slice()) {
            (Some(Some(og)), [g]) if *g == og => "agree_assigned",
            (Some(None), []) => "agree_none",
            (Some(None), [_]) => "oracle_none_ours_one",
            (Some(Some(_)), []) => "oracle_one_ours_none",
            (Some(Some(_)), [_]) => "both_one_different_gene",
            (Some(Some(_)), _) => "oracle_one_ours_multi",
            (Some(None), _) => "oracle_none_ours_multi",
            (None, _) => "oracle_multi",
        };
        *counts.entry(class).or_insert(0) += 1;

        let show = shown.entry(class).or_insert(0);
        if *show < args.examples && !class.starts_with("agree") && class != "oracle_multi" {
            *show += 1;
            let ours_names: Vec<&str> = ours.iter().map(|&g| anno.gene_ids[g as usize].as_str()).collect();
            eprintln!(
                "[{class}] {}:{} blocks={:?} junctions={:?} strand={:?} GX={:?} ours={:?}",
                chrom_names[chrom_bam], pos, p.blocks, p.junctions, p.strand, gx, ours_names
            );
        }
    }

    println!("=== assign-diff: {} unique-alignment reads ===", n);
    let mut rows: Vec<_> = counts.iter().collect();
    rows.sort_by_key(|(_, v)| std::cmp::Reverse(**v));
    for (k, v) in &rows {
        println!("  {k:<26} {v:>12}  ({:.4}%)", 100.0 * **v as f64 / n as f64);
    }

    if let Some(p) = &args.json_out {
        if let Some(d) = p.parent() {
            std::fs::create_dir_all(d)?;
        }
        let obj: serde_json::Map<String, serde_json::Value> =
            counts.iter().map(|(k, v)| ((*k).to_string(), serde_json::json!(v))).collect();
        std::fs::write(p, serde_json::to_string_pretty(&serde_json::Value::Object(obj))?)?;
    }
    Ok(())
}
