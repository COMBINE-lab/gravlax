//! Compare annotation-independent UMI grouping with STARsolo's per-gene grouping.
//!
//! STARsolo collapses UMIs **per (cell, gene)** — `SoloFeature_collapseUMIall.cpp` loops
//! `for (iG=0; iG<nGenes; iG++)`. That grouping is annotation-dependent, whereas a reusable
//! molecule archive must group reads before a gene annotation is supplied.
//!
//! Two groupings are compared over the *same* alignments:
//!
//! * **STARsolo grouping** — a molecule is a distinct corrected `UB` within a `(cell, gene)` group.
//! * **Annotation-independent grouping** — a molecule is keyed by cell, UMI, genomic locus, and
//!   orientation. Loci come from single-linkage clustering of read positions rather than
//!   fixed-width bucketing, which can split molecules that straddle a boundary.
//!
//! The annotation-independent grouping reads `UR` (raw UMI), never `UB`: STARsolo corrects UMIs
//! per gene, so grouping on `UB` would reintroduce annotation-derived information.
//!
//! The reported difference quantifies the effect of grouping without a gene model and can be
//! compared with changes caused by choosing a different annotation version.

use anyhow::{Context, Result};
use clap::Parser;
use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Value;
use noodles_sam::alignment::record::Flags;
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::PathBuf;

#[derive(Parser)]
pub struct Args {
    /// Coordinate-sorted STARsolo reference BAM carrying CB/UB/UR/GX.
    pub bam: PathBuf,
    /// Maximum gap, in bp, between consecutive read starts of one molecule. Reads of a single 10x
    /// molecule pile at a shared 3' end, so this only has to absorb fragment-length jitter.
    #[arg(long, default_value_t = 500)]
    pub locus_gap: u32,
    /// STARsolo `Solo.out/Gene/filtered/barcodes.tsv`. Per-cell statistics are meaningless without
    /// it: the unfiltered whitelist is ~65k empty droplets carrying a handful of counts each, on
    /// which G1 and G2 agree trivially, so the median over all barcodes reads 0% however bad the
    /// grouping is.
    #[arg(long)]
    pub cells: Option<PathBuf>,
    /// Merge UMIs within one substitution of a more abundant UMI at the same locus. STARsolo
    /// applies `1MM_CR` per (cell, gene); without an equivalent, G2 over-splits and the comparison
    /// measures missing error-correction rather than missing annotation. The correction runs inside
    /// a genomic locus and consults no annotation.
    #[arg(long)]
    pub correct_umi: bool,
    #[arg(long)]
    pub json_out: Option<PathBuf>,
}

#[derive(Default)]
struct Tags {
    gx: Option<String>,
    cb: Option<String>,
    ub: Option<String>,
    ur: Option<String>,
}

/// Pull the string tags we need in one pass. `data().get()` does not match STAR's tags reliably
/// here, so fields are scanned directly.
fn scan_tags(rec: &bam::Record) -> Result<Tags> {
    let mut t = Tags::default();
    for field in rec.data().iter() {
        let (tag, value) = field?;
        let key = <[u8; 2]>::from(tag);
        let Value::String(s) = value else { continue };
        let s = String::from_utf8_lossy(s.as_ref()).into_owned();
        match &key {
            b"GX" => t.gx = Some(s),
            b"CB" => t.cb = Some(s),
            b"UB" => t.ub = Some(s),
            b"UR" => t.ur = Some(s),
            _ => {}
        }
    }
    Ok(t)
}

/// True when two equal-length sequences differ in exactly one position.
fn hamming1(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut diff = 0;
    for (x, y) in a.iter().zip(b) {
        if x != y {
            diff += 1;
            if diff > 1 {
                return false;
            }
        }
    }
    diff == 1
}

/// Map each UMI to the canonical UMI it collapses into, most-abundant-first, within one
/// substitution. Mirrors the shape of STARsolo's `1MM_CR`, but scoped to a genomic locus.
fn correct_umis_1mm(counts: &FxHashMap<String, u32>) -> FxHashMap<String, String> {
    let mut order: Vec<(&String, &u32)> = counts.iter().collect();
    // Ties broken lexicographically so the result never depends on hash iteration order.
    order.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));

    let mut map: FxHashMap<String, String> = FxHashMap::default();
    for (umi, n) in &order {
        let mut target: Option<String> = None;
        for (cand, cn) in &order {
            if cand == umi {
                break; // nothing more abundant remains
            }
            if *cn > *n && hamming1(cand.as_bytes(), umi.as_bytes()) {
                target = Some(map.get(*cand).cloned().unwrap_or_else(|| (*cand).clone()));
                break;
            }
        }
        map.insert((*umi).clone(), target.unwrap_or_else(|| (*umi).clone()));
    }
    map
}

/// One aligned read's grouping evidence, plus the gene the oracle assigned it.
struct Read {
    pos: u32,
    umi: String,
    gene: String,
}

pub fn run(args: Args) -> Result<()> {
    let cells: Option<FxHashSet<String>> = match &args.cells {
        Some(p) => Some(
            std::fs::read_to_string(p)?
                .lines()
                .map(|l| l.trim().to_string())
                .filter(|l| !l.is_empty())
                .collect(),
        ),
        None => None,
    };

    let mut reader = bam::io::reader::Builder
        .build_from_path(&args.bam)
        .with_context(|| format!("opening {}", args.bam.display()))?;
    reader.read_header()?;

    // Reference grouping: (cell, gene) -> distinct corrected UMIs. Reproduces the reference
    // matrix entry.
    let mut g1: FxHashMap<(String, String), FxHashSet<String>> = FxHashMap::default();
    // Annotation-independent input: (cell, chrom, strand) -> reads.
    let mut g2_in: FxHashMap<(String, i32, bool), Vec<Read>> = FxHashMap::default();

    let mut rec = bam::Record::default();
    let (mut n_reads, mut n_used) = (0u64, 0u64);

    while reader.read_record(&mut rec)? != 0 {
        let flags = Flags::from(rec.flags().bits());
        if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped() {
            continue;
        }
        n_reads += 1;

        // GX is "-" when unassigned and comma-joined when multi-gene; STARsolo excludes both from
        // the Gene matrix, so skip them and let the two groupings see one population.
        let t = scan_tags(&rec)?;
        let gene = match t.gx {
            Some(g) if !g.is_empty() && g != "-" && !g.contains(',') => g,
            _ => continue,
        };
        let (cb, ub, ur) = match (t.cb, t.ub, t.ur) {
            (Some(c), Some(u), Some(r)) if c != "-" && u != "-" => (c, u, r),
            _ => continue,
        };
        n_used += 1;

        g1.entry((cb.clone(), gene.clone())).or_default().insert(ub);

        let chrom = rec.reference_sequence_id().transpose()?.unwrap_or(0) as i32;
        let pos = rec
            .alignment_start()
            .transpose()?
            .map(|p| usize::from(p) as u32 - 1)
            .unwrap_or(0);
        g2_in
            .entry((cb, chrom, flags.is_reverse_complemented()))
            .or_default()
            .push(Read { pos, umi: ur, gene });
    }

    // Resolve G2: cluster reads into loci by single linkage, then into molecules by UMI.
    let mut g2_matrix: FxHashMap<(String, String), u64> = FxHashMap::default();
    let (mut n_molecules, mut multi_gene, mut merged) = (0u64, 0u64, 0u64);

    for ((cb, _chrom, _strand), mut reads) in g2_in {
        reads.sort_by_key(|r| r.pos);
        let mut start = 0usize;
        for i in 1..=reads.len() {
            let split = i == reads.len() || reads[i].pos - reads[i - 1].pos > args.locus_gap;
            if !split {
                continue;
            }
            let locus = &reads[start..i];
            start = i;

            let mut by_umi: FxHashMap<&str, (u32, FxHashSet<&str>)> = FxHashMap::default();
            for r in locus {
                let e = by_umi.entry(&r.umi).or_insert_with(|| (0, FxHashSet::default()));
                e.0 += 1;
                e.1.insert(&r.gene);
            }

            let mut molecules: FxHashMap<String, FxHashSet<&str>> = FxHashMap::default();
            if args.correct_umi {
                let counts: FxHashMap<String, u32> =
                    by_umi.iter().map(|(u, (n, _))| ((*u).to_string(), *n)).collect();
                let map = correct_umis_1mm(&counts);
                for (u, (_, genes)) in &by_umi {
                    let canon = map.get(*u).cloned().unwrap_or_else(|| (*u).to_string());
                    if canon != **u {
                        merged += 1;
                    }
                    molecules.entry(canon).or_default().extend(genes.iter().copied());
                }
            } else {
                for (u, (_, genes)) in &by_umi {
                    molecules.insert((*u).to_string(), genes.clone());
                }
            }

            for genes in molecules.values() {
                n_molecules += 1;
                if genes.len() > 1 {
                    multi_gene += 1;
                }
                for g in genes {
                    *g2_matrix.entry((cb.clone(), (*g).to_string())).or_insert(0) += 1;
                }
            }
        }
    }

    // Compare the two matrices.
    let mut g1_total = 0u64;
    let mut per_cell_g1: FxHashMap<String, u64> = FxHashMap::default();
    let mut per_cell_err: FxHashMap<String, u64> = FxHashMap::default();
    let mut l1 = 0u64;

    for (key, umis) in &g1 {
        let a = umis.len() as u64;
        let b = *g2_matrix.get(key).unwrap_or(&0);
        g1_total += a;
        *per_cell_g1.entry(key.0.clone()).or_insert(0) += a;
        *per_cell_err.entry(key.0.clone()).or_insert(0) += a.abs_diff(b);
        l1 += a.abs_diff(b);
    }
    let mut g2_only = 0u64;
    for (key, b) in &g2_matrix {
        if !g1.contains_key(key) {
            g2_only += b;
            l1 += b;
            *per_cell_err.entry(key.0.clone()).or_insert(0) += *b;
        }
    }

    let mut rels: Vec<f64> = per_cell_g1
        .iter()
        .filter(|(c, &t)| t > 0 && cells.as_ref().is_none_or(|s| s.contains(*c)))
        .map(|(c, &t)| *per_cell_err.get(c).unwrap_or(&0) as f64 / t as f64)
        .collect();
    rels.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let pick = |q: f64| {
        if rels.is_empty() { 0.0 } else { rels[((rels.len() - 1) as f64 * q) as usize] }
    };
    let (median, p90) = (pick(0.5), pick(0.9));

    let g2_total: u64 = g2_matrix.values().sum();
    let overall = if g1_total > 0 { l1 as f64 / g1_total as f64 } else { 0.0 };

    println!("=== Annotation-independent UMI grouping comparison ===");
    println!("locus gap: {} bp   UMI 1MM correction: {}", args.locus_gap, args.correct_umi);
    println!("mapped primary reads:          {n_reads}");
    println!("unique-gene reads used:        {n_used}");
    println!("G1 molecules (STARsolo):       {g1_total}");
    println!("G2 molecules (annot-free):     {n_molecules}  (matrix mass {g2_total})");
    println!("G2 molecules spanning >1 gene: {multi_gene}");
    println!("UMIs merged by 1MM:            {merged}");
    println!("count mass only in G2:         {g2_only}");
    println!("overall L1 / G1 total:         {:.4}%", 100.0 * overall);
    println!(
        "cells in per-cell stats:       {}  ({})",
        rels.len(),
        if cells.is_some() { "filtered cell list" } else { "ALL barcodes — not meaningful" }
    );
    println!("median per-cell L1 rel err:    {:.4}%   [PASS < 1%]", 100.0 * median);
    println!("p90 per-cell L1 rel err:       {:.4}%", 100.0 * p90);

    if let Some(p) = &args.json_out {
        if let Some(d) = p.parent() {
            std::fs::create_dir_all(d)?;
        }
        let obj = serde_json::json!({
            "locus_gap": args.locus_gap,
            "correct_umi": args.correct_umi,
            "mapped_primary_reads": n_reads,
            "unique_gene_reads": n_used,
            "g1_molecules": g1_total,
            "g2_molecules": n_molecules,
            "g2_matrix_mass": g2_total,
            "g2_multi_gene_molecules": multi_gene,
            "umis_merged_1mm": merged,
            "mass_only_in_g2": g2_only,
            "overall_l1_rel": overall,
            "cells_in_stats": rels.len(),
            "median_per_cell_l1_rel": median,
            "p90_per_cell_l1_rel": p90,
        });
        std::fs::write(p, serde_json::to_string_pretty(&obj)?)?;
        println!("wrote {}", p.display());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hamming1_is_exactly_one_substitution() {
        assert!(hamming1(b"ACGT", b"ACGA"));
        assert!(!hamming1(b"ACGT", b"ACGT"), "identical is zero differences, not one");
        assert!(!hamming1(b"ACGT", b"ACAA"));
        assert!(!hamming1(b"ACGT", b"ACG"), "different lengths never match");
    }

    #[test]
    fn correction_absorbs_rare_umi_into_abundant_neighbour() {
        let mut c = FxHashMap::default();
        c.insert("AAAAAAAAAAAA".to_string(), 100u32);
        c.insert("AAAAAAAAAAAC".to_string(), 1u32);
        let m = correct_umis_1mm(&c);
        assert_eq!(m["AAAAAAAAAAAC"], "AAAAAAAAAAAA");
        assert_eq!(m["AAAAAAAAAAAA"], "AAAAAAAAAAAA");
    }

    #[test]
    fn correction_leaves_equally_abundant_umis_alone() {
        // Without a strict abundance ordering, merging would be arbitrary and order-dependent.
        let mut c = FxHashMap::default();
        c.insert("AAAAAAAAAAAA".to_string(), 5u32);
        c.insert("AAAAAAAAAAAC".to_string(), 5u32);
        let m = correct_umis_1mm(&c);
        assert_eq!(m["AAAAAAAAAAAA"], "AAAAAAAAAAAA");
        assert_eq!(m["AAAAAAAAAAAC"], "AAAAAAAAAAAC");
    }

    #[test]
    fn correction_is_transitive_to_the_most_abundant() {
        // C(1) -> B(5) -> A(100): C must land on A, not stop at B.
        let mut c = FxHashMap::default();
        c.insert("AAAAAAAAAAAA".to_string(), 100u32);
        c.insert("AAAAAAAAAAAC".to_string(), 5u32);
        c.insert("AAAAAAAAAACC".to_string(), 1u32);
        let m = correct_umis_1mm(&c);
        assert_eq!(m["AAAAAAAAAAAC"], "AAAAAAAAAAAA");
        assert_eq!(m["AAAAAAAAAACC"], "AAAAAAAAAAAA");
    }

    #[test]
    fn two_substitutions_are_not_merged() {
        let mut c = FxHashMap::default();
        c.insert("AAAAAAAAAAAA".to_string(), 100u32);
        c.insert("AAAAAAAAAACC".to_string(), 1u32);
        let m = correct_umis_1mm(&c);
        assert_eq!(m["AAAAAAAAAACC"], "AAAAAAAAAACC");
    }
}
