//! `aie sig-stats` — measure the evidence structures that determine archive size.
//!
//! The format's molecule payload is defined so fidelity is inherited BY CONSTRUCTION rather than
//! approximated: a molecule stores its distinct *read signatures* with counts, where a signature
//! is the complete placement set of a read (one placement for a unique read, the whole alternative
//! set for a multimapper). Replay over signatures reproduces read-level semantics exactly — the
//! union rule and containment behavior. What remains is a storage question, which this command
//! measures on real data:
//!
//!  1. Signature multiplicity — distinct signatures per molecule (how much diversity the payload
//!     actually carries).
//!  2. Paralog patterns — a multimapper's alternatives sit at fixed genomic offsets from its
//!     anchor, conserved across every read of a repeat family. If the distinct (chrom, offset)
//!     PATTERNS are few, an interned pattern dictionary makes E4 multimapper evidence nearly free.
//!  3. Position structure — molecule starts pile at 3' sites. Two-level (site anchor + offset)
//!     coding is measured against the flat delta coding the current encoder uses.
//!  4. Cell stream orderings — (pos, cell) vs (site, cell) secondary sort, since the cell stream
//!     is the largest non-UMI stream (11.3 bits/molecule).
//!
//! Everything is reported as real zstd-compressed stream sizes, comparable to `archive::encode`.

use crate::build::{correct_barcode, load_whitelist, to_ops};
use anyhow::{Context, Result};
use clap::Parser;
use evidence_io::archive::{put_svarint, put_varint};
use evidence_io::{umi, Placement};
use ingest::placement_from_alignment;
use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Value;
use noodles_sam::alignment::record::Flags;
use rustc_hash::FxHashMap;
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser)]
pub struct Args {
    /// Annotation-free ingest BAM with CR/UR, secondaries included.
    pub bam: PathBuf,
    #[arg(long)]
    pub whitelist: PathBuf,
    #[arg(long, default_value_t = 2_000)]
    pub locus_gap: u32,
    /// Gap between molecule starts that closes a 3' site (positions within a site are offsets).
    #[arg(long, default_value_t = 64)]
    pub site_gap: u32,
    #[arg(long, default_value_t = 19)]
    pub zstd_level: i32,
    #[arg(long)]
    pub json_out: Option<PathBuf>,
}

fn zbytes(buf: &[u8], level: i32) -> Result<u64> {
    let mut e = zstd::Encoder::new(Vec::new(), level)?;
    e.write_all(buf)?;
    Ok(e.finish()?.len() as u64)
}

fn name_hash(name: &[u8]) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = rustc_hash::FxHasher::default();
    name.hash(&mut h);
    h.finish()
}

fn placement_key(p: &Placement) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = rustc_hash::FxHasher::default();
    p.assignment_key().hash(&mut h);
    h.finish()
}

/// One read reduced to signature-building parts.
struct SRead {
    cell: u32,
    umi: u32,
    /// Primary placement (molecule anchor and locus key).
    prim: Placement,
    /// For multimappers: the full alternative set including the primary, sorted.
    alts: Option<Vec<Placement>>,
}

pub fn run(args: Args) -> Result<()> {
    let wl = load_whitelist(&args.whitelist)?;
    let mut reader = bam::io::reader::Builder
        .build_from_path(&args.bam)
        .with_context(|| format!("opening {}", args.bam.display()))?;
    reader.read_header()?;

    // Multimapper records arrive scattered (coordinate order), so they are collected by read name
    // and assembled after the scan. Primaries of unique reads stream straight into `reads`.
    type MultimapRead = (Option<(u32, u32)>, Option<usize>, Vec<Placement>);
    let mut mm: FxHashMap<u64, MultimapRead> = FxHashMap::default();
    let mut reads: Vec<SRead> = Vec::new();
    let mut bc_cache: FxHashMap<u32, Option<u32>> = FxHashMap::default();
    let mut rec = bam::Record::default();
    let (mut n_prim, mut n_sec) = (0u64, 0u64);

    while reader.read_record(&mut rec)? != 0 {
        let flags = Flags::from(rec.flags().bits());
        if flags.is_supplementary() || flags.is_unmapped() {
            continue;
        }
        let (mut cr, mut ur, mut nm, mut score, mut nh) = (None, None, 0u16, 0i32, 1u16);
        for field in rec.data().iter() {
            let (tag, value) = field?;
            let key = <[u8; 2]>::from(tag);
            match value {
                Value::String(s) => match &key {
                    b"CR" => cr = umi::pack(s.as_ref()),
                    b"UR" => ur = umi::pack(s.as_ref()),
                    _ => {}
                },
                v => {
                    let iv = match v {
                        Value::Int8(x) => x as i64,
                        Value::UInt8(x) => x as i64,
                        Value::Int16(x) => x as i64,
                        Value::UInt16(x) => x as i64,
                        Value::Int32(x) => x as i64,
                        Value::UInt32(x) => x as i64,
                        _ => continue,
                    };
                    match &key {
                        b"nM" => nm = iv as u16,
                        b"AS" => score = iv as i32,
                        b"NH" => nh = iv as u16,
                        _ => {}
                    }
                }
            }
        }
        let chrom = rec.reference_sequence_id().transpose()?.unwrap_or(0) as u32;
        let pos = rec.alignment_start().transpose()?.map(|p| usize::from(p) as u32 - 1).unwrap_or(0);
        let p = placement_from_alignment(
            chrom, pos, flags.is_reverse_complemented(), &to_ops(rec.cigar().iter())?, nm, score, nh,
        );

        if nh > 1 {
            let key = name_hash(rec.name().map(|n| n.as_ref()).unwrap_or(b""));
            let e = mm.entry(key).or_insert((None, None, Vec::new()));
            if !flags.is_secondary() {
                n_prim += 1;
                let cb = match (cr, ur) {
                    (Some(c), Some(u)) => bc_cache
                        .entry(c)
                        .or_insert_with(|| correct_barcode(c, &wl))
                        .map(|cc| (cc, u)),
                    _ => None,
                };
                e.0 = cb;
                e.1 = Some(e.2.len());
            } else {
                n_sec += 1;
            }
            e.2.push(p);
            continue;
        }

        if flags.is_secondary() {
            continue;
        }
        n_prim += 1;
        let (Some(c), Some(u)) = (cr, ur) else { continue };
        let Some(cb) = *bc_cache.entry(c).or_insert_with(|| correct_barcode(c, &wl)) else { continue };
        reads.push(SRead { cell: cb, umi: u, prim: p, alts: None });
    }
    eprintln!("primaries {n_prim}, secondaries {n_sec}, mm reads {}", mm.len());

    for (_k, (cb, prim_idx, mut alts)) in mm.drain() {
        let (Some((cell, u)), Some(pi)) = (cb, prim_idx) else { continue };
        let prim = alts[pi].clone();
        alts.sort_by_key(|p| (p.chrom, p.start()));
        reads.push(SRead { cell, umi: u, prim, alts: Some(alts) });
    }

    // Dense cell interning by first occurrence — the archive's real coding. Using raw packed
    // barcodes here inflated the first run's cell streams (32 random bits vs ~log2(#cells)).
    {
        let mut intern: FxHashMap<u32, u32> = FxHashMap::default();
        for r in reads.iter_mut() {
            let next = intern.len() as u32;
            r.cell = *intern.entry(r.cell).or_insert(next);
        }
        eprintln!("distinct cells: {}", intern.len());
    }

    // ---- Molecules: (cell, chrom, strand) → locus clusters → exact UMI → signature multiset ----
    reads.sort_by_key(|r| (r.cell, r.prim.chrom, matches!(r.prim.strand, evidence_io::Strand::Reverse), r.prim.start()));

    struct MolSig {
        anchor_chrom: u32,
        anchor_pos: u32,
        cell: u32,
        sigs: Vec<(u64, u32, bool)>, // (signature hash, read count, is_multimapper)
        n_chains: u32,               // distinct junction chains among unique-read signatures
    }
    let mut molecules: Vec<MolSig> = Vec::new();
    // Pattern dictionary for multimapper signatures: offsets from the read's own anchor.
    let mut patterns: FxHashMap<Vec<(u32, i64, bool)>, u32> = FxHashMap::default();
    let mut pattern_uses: Vec<(u32, u64)> = Vec::new(); // (pattern id, signature hash) per mm sig
    let sig_of_read = |r: &SRead, patterns: &mut FxHashMap<Vec<(u32, i64, bool)>, u32>| -> (u64, bool, Option<u32>) {
        match &r.alts {
            None => (placement_key(&r.prim), false, None),
            Some(alts) => {
                use std::hash::{Hash, Hasher};
                let mut h = rustc_hash::FxHasher::default();
                let a0 = &alts[0];
                let pat: Vec<(u32, i64, bool)> = alts
                    .iter()
                    .map(|p| (
                        p.chrom,
                        p.start() as i64 - a0.start() as i64,
                        matches!(p.strand, evidence_io::Strand::Reverse) != matches!(a0.strand, evidence_io::Strand::Reverse),
                    ))
                    .collect();
                for p in alts {
                    p.assignment_key().hash(&mut h);
                }
                let next = patterns.len() as u32;
                let pid = *patterns.entry(pat).or_insert(next);
                (h.finish(), true, Some(pid))
            }
        }
    };

    let mut i = 0usize;
    while i < reads.len() {
        let (cell, chrom, strand) = (
            reads[i].cell,
            reads[i].prim.chrom,
            matches!(reads[i].prim.strand, evidence_io::Strand::Reverse),
        );
        let mut j = i;
        while j < reads.len()
            && reads[j].cell == cell
            && reads[j].prim.chrom == chrom
            && matches!(reads[j].prim.strand, evidence_io::Strand::Reverse) == strand
        {
            j += 1;
        }
        let group = &reads[i..j];
        let mut ls = 0usize;
        for k in 1..=group.len() {
            let split = k == group.len() || group[k].prim.start() - group[k - 1].prim.start() > args.locus_gap;
            if !split {
                continue;
            }
            let locus = &group[ls..k];
            ls = k;
            let mut by_umi: FxHashMap<u32, FxHashMap<u64, (u32, bool)>> = FxHashMap::default();
            let mut anchor: FxHashMap<u32, u32> = FxHashMap::default();
            let mut chains: FxHashMap<u32, rustc_hash::FxHashSet<u64>> = FxHashMap::default();
            for r in locus {
                let (sig, is_mm, pid) = sig_of_read(r, &mut patterns);
                if let Some(pid) = pid {
                    pattern_uses.push((pid, sig));
                }
                let e = by_umi.entry(r.umi).or_default().entry(sig).or_insert((0, is_mm));
                e.0 += 1;
                anchor.entry(r.umi).or_insert(r.prim.start());
                if !is_mm {
                    use std::hash::{Hash, Hasher};
                    let mut h = rustc_hash::FxHasher::default();
                    r.prim.junctions.hash(&mut h);
                    chains.entry(r.umi).or_default().insert(h.finish());
                }
            }
            for (u, sigs) in by_umi {
                molecules.push(MolSig {
                    anchor_chrom: chrom,
                    anchor_pos: anchor[&u],
                    cell,
                    sigs: sigs.into_iter().map(|(s, (c, m))| (s, c, m)).collect(),
                    n_chains: chains.get(&u).map(|c| c.len() as u32).unwrap_or(0),
                });
            }
        }
        i = j;
    }

    // ---- 1. Signature multiplicity ----
    let n_mol = molecules.len() as u64;
    let mut hist = [0u64; 5];
    let mut mm_mols = 0u64;
    for m in &molecules {
        hist[m.sigs.len().min(4)] += 1;
        if m.sigs.iter().any(|(_, _, mm)| *mm) {
            mm_mols += 1;
        }
    }
    // Chain diversity is tracked at signature-build time via chain hashes carried alongside.
    let total_sigs: u64 = molecules.iter().map(|m| m.sigs.len() as u64).sum();
    let mut chain_hist = [0u64; 4];
    for m in &molecules {
        chain_hist[(m.n_chains as usize).min(3)] += 1;
    }
    let total_chains: u64 = molecules.iter().map(|m| m.n_chains as u64).sum();

    // ---- 2. Paralog patterns ----
    let mm_sigs = pattern_uses.len() as u64;
    let distinct_patterns = patterns.len() as u64;
    let mut use_count: FxHashMap<u32, u64> = FxHashMap::default();
    for (pid, _) in &pattern_uses {
        *use_count.entry(*pid).or_insert(0) += 1;
    }
    let mut uses: Vec<u64> = use_count.values().copied().collect();
    uses.sort_unstable_by(|a, b| b.cmp(a));
    let top100: u64 = uses.iter().take(100).sum();
    let pattern_dict_entries: u64 = patterns.keys().map(|p| p.len() as u64).sum();
    // Encoded pattern dictionary: per entry (chrom varint, offset svarint, flip bit in chrom lsb).
    let mut pd = Vec::new();
    for pat in patterns.keys() {
        put_varint(&mut pd, pat.len() as u64);
        for (c, off, flip) in pat {
            put_varint(&mut pd, ((*c as u64) << 1) | (*flip as u64));
            put_svarint(&mut pd, *off);
        }
    }
    let pattern_dict_bytes = zbytes(&pd, args.zstd_level)?;
    // Per-mm-signature pattern-id stream.
    let mut pid_stream = Vec::new();
    for (pid, _) in &pattern_uses {
        put_varint(&mut pid_stream, *pid as u64);
    }
    let pid_bytes = zbytes(&pid_stream, args.zstd_level)?;

    // ---- 3 & 4. Position/site and cell-ordering experiments (over molecule anchors) ----
    molecules.sort_by_key(|m| (m.anchor_chrom, m.anchor_pos, m.cell));
    let enc_flat = {
        let (mut pos, mut cellb) = (Vec::new(), Vec::new());
        let (mut lc, mut lp, mut lcell) = (u32::MAX, 0u32, 0u32);
        for m in &molecules {
            if m.anchor_chrom != lc {
                lc = m.anchor_chrom;
                lp = 0;
            }
            put_svarint(&mut pos, m.anchor_pos as i64 - lp as i64);
            lp = m.anchor_pos;
            put_svarint(&mut cellb, m.cell as i64 - lcell as i64);
            lcell = m.cell;
        }
        (zbytes(&pos, args.zstd_level)?, zbytes(&cellb, args.zstd_level)?)
    };

    // Two-level: site anchors + per-molecule offsets; then cell-sorted within site.
    let (enc_site, n_sites) = {
        let (mut site_anchor, mut offset, mut cellb) = (Vec::new(), Vec::new(), Vec::new());
        let mut n_sites = 0u64;
        let (mut i, mut last_anchor, mut lc, mut lcell) = (0usize, 0u32, u32::MAX, 0u32);
        while i < molecules.len() {
            let (chrom, start) = (molecules[i].anchor_chrom, molecules[i].anchor_pos);
            let mut j = i;
            while j + 1 < molecules.len()
                && molecules[j + 1].anchor_chrom == chrom
                && molecules[j + 1].anchor_pos - molecules[j].anchor_pos <= args.site_gap
            {
                j += 1;
            }
            n_sites += 1;
            if chrom != lc {
                lc = chrom;
                last_anchor = 0;
            }
            put_varint(&mut site_anchor, (start - last_anchor) as u64);
            last_anchor = start;
            // Within the site, order by cell: offsets lose monotonicity (svarint), cells gain runs.
            let mut site: Vec<(u32, u32)> = molecules[i..=j].iter().map(|m| (m.cell, m.anchor_pos - start)).collect();
            site.sort_unstable();
            put_varint(&mut offset, site.len() as u64);
            let mut lo = 0i64;
            for (c, o) in site {
                put_svarint(&mut offset, o as i64 - lo);
                lo = o as i64;
                put_svarint(&mut cellb, c as i64 - lcell as i64);
                lcell = c;
            }
            i = j + 1;
        }
        ((zbytes(&site_anchor, args.zstd_level)?, zbytes(&offset, args.zstd_level)?, zbytes(&cellb, args.zstd_level)?), n_sites)
    };

    let bits = |b: u64| 8.0 * b as f64 / n_mol.max(1) as f64;
    println!("=== signature statistics ({n_mol} molecules, {total_sigs} signatures) ===");
    println!("signatures per molecule: 1: {} ({:.3}%)  2: {}  3: {}  >=4: {}",
        hist[1], 100.0 * hist[1] as f64 / n_mol as f64, hist[2], hist[3], hist[4]);
    println!("molecules containing any multimapper signature: {mm_mols} ({:.3}%)", 100.0 * mm_mols as f64 / n_mol as f64);
    println!("junction chains per molecule: 0: {}  1: {} ({:.3}%)  2: {}  >=3: {}   total chains {}",
        chain_hist[0], chain_hist[1], 100.0 * chain_hist[1] as f64 / n_mol as f64, chain_hist[2], chain_hist[3], total_chains);
    println!();
    println!("=== paralog patterns ===");
    println!("multimapper signatures: {mm_sigs}   distinct offset-patterns: {distinct_patterns}  ({:.2}x sharing)",
        mm_sigs as f64 / distinct_patterns.max(1) as f64);
    println!("top-100 patterns cover: {:.2}% of mm signatures", 100.0 * top100 as f64 / mm_sigs.max(1) as f64);
    println!("pattern dictionary: {pattern_dict_entries} entries -> {pattern_dict_bytes} B zstd; pattern-id stream {pid_bytes} B ({:.3} bits/mol)",
        bits(pattern_dict_bytes + pid_bytes));
    println!();
    println!("=== position & cell stream experiments ===");
    println!("flat  (pos,cell) sort:  pos {} B ({:.3} b/mol)   cell {} B ({:.3} b/mol)", enc_flat.0, bits(enc_flat.0), enc_flat.1, bits(enc_flat.1));
    println!("sites ({} sites, gap {}): anchors+offsets {}+{} B ({:.3} b/mol)   cell {} B ({:.3} b/mol)",
        n_sites, args.site_gap, enc_site.0, enc_site.1, bits(enc_site.0 + enc_site.1), enc_site.2, bits(enc_site.2));
    println!("total pos+cell: flat {:.3} vs site-ordered {:.3} bits/mol",
        bits(enc_flat.0 + enc_flat.1), bits(enc_site.0 + enc_site.1 + enc_site.2));

    if let Some(p) = &args.json_out {
        if let Some(d) = p.parent() {
            std::fs::create_dir_all(d)?;
        }
        let obj = serde_json::json!({
            "molecules": n_mol, "signatures": total_sigs, "sig_hist": hist.to_vec(),
            "mm_molecules": mm_mols, "mm_signatures": mm_sigs,
            "chain_hist": chain_hist.to_vec(), "total_chains": total_chains,
            "distinct_patterns": distinct_patterns, "top100_cover": top100 as f64 / mm_sigs.max(1) as f64,
            "pattern_dict_bytes": pattern_dict_bytes, "pattern_id_bytes": pid_bytes,
            "flat_pos_bytes": enc_flat.0, "flat_cell_bytes": enc_flat.1,
            "n_sites": n_sites, "site_anchor_bytes": enc_site.0, "site_offset_bytes": enc_site.1, "site_cell_bytes": enc_site.2,
        });
        std::fs::write(p, serde_json::to_string_pretty(&obj)?)?;
        println!("wrote {}", p.display());
    }
    Ok(())
}
