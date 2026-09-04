//! Measure whether the *collision structure* of UMIs, rather than their nucleotide values,
//! preserves per-gene UMI collapse across annotations.
//!
//! The UMI stream occupies 56% of an evaluated archive at 27 bits per molecule, but replay never
//! reads a UMI *value*. Values are consulted for exactly two decisions: (a) molecules with the
//! same value in one `(cell, gene)` group are one molecule; (b) a value within one substitution of
//! a more abundant value in that group is absorbed into it. Both decisions depend only on the
//! sparse equivalence and adjacency relations between molecules. Per cell, the UMI space
//! (4^12 ≈ 16.7M) is much larger than the molecule count (~10^4), so collisions need not be a
//! per-molecule payload.
//!
//! This module measures three things on real data:
//!  1. Graph statistics — same-value links and 1-substitution edges per cell, as a function of the
//!     genomic window they span (the window bounds which future gene spans the graph can serve).
//!  2. Preservation — replay per-(cell, gene) collapse using only molecules + read counts + the
//!     graph, and compare it with STARsolo's own collapse (the `UB` tag). Two 1MM policies are
//!     compared because a strictly-more-abundant policy left a 0.9% difference, while
//!     CellRanger-style 1MM_CR also merges ties.
//!  3. Bytes — a real encoding of the graph, to set against the 67.2 MB UMI stream it replaces.
//!
//! On a BAM without `GX` (the annotation-free ingest), only (1) and (3) run.

use anyhow::{Context, Result};
use clap::Parser;
use evidence_io::archive::put_varint;
use evidence_io::umi;
use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Value;
use noodles_sam::alignment::record::Flags;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser)]
pub struct Args {
    /// BAM with CR/UR and, when collapse comparison is requested, reference CB/UB/GX tags.
    pub bam: PathBuf,
    /// Filtered cell list for per-cell statistics.
    #[arg(long)]
    pub cells: Option<PathBuf>,
    /// Locus gap for the base (ingest-time) molecule definition. Reads of one physical molecule
    /// pile within a fragment length, so this is deliberately small; everything beyond it is the
    /// graph's job at replay time, not the ingest's.
    #[arg(long, default_value_t = 2_000)]
    pub locus_gap: u32,
    /// Window for edges retained in the STORED graph. Bounds the gene spans replay can serve;
    /// statistics are reported across a range of windows regardless.
    #[arg(long, default_value_t = 3_000_000)]
    pub store_window: u32,
    #[arg(long)]
    pub json_out: Option<PathBuf>,
}

#[derive(Default)]
struct Tags {
    gx: Option<String>,
    cb: Option<String>,
    ub: Option<u32>,
    ur: Option<u32>,
    cr: Option<u32>,
}

fn scan_tags(rec: &bam::Record) -> Result<Tags> {
    let mut t = Tags::default();
    for field in rec.data().iter() {
        let (tag, value) = field?;
        let key = <[u8; 2]>::from(tag);
        let Value::String(s) = value else { continue };
        match &key {
            b"GX" => t.gx = Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
            b"CB" => t.cb = Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
            b"UB" => t.ub = umi::pack(s.as_ref()),
            b"UR" => t.ur = umi::pack(s.as_ref()),
            b"CR" => t.cr = umi::pack(s.as_ref()),
            _ => {}
        }
    }
    Ok(t)
}

/// One read, reduced to grouping evidence. `gene == u32::MAX` means unassigned.
struct Read {
    cell: u32,
    chrom: i32,
    strand_rev: bool,
    pos: u32,
    umi: u32,
    gene: u32,
}

/// A base molecule: one exact UMI value at one locus. No 1MM correction is applied at ingest in
/// graph mode — every merge is deferred to replay so no future annotation's grouping is foreclosed.
struct Mol {
    umi: u32,
    pos: u32,
    /// (gene, reads) pairs; genes are available only from the reference BAM.
    genes: Vec<(u32, u32)>,
    n_reads: u32,
}

/// The 36 packed values within one substitution over the 12 positions used by 10x v3. Fixed 10 bp
/// UMIs are represented in the low 20 bits: their 30 real neighbours are included, while the six
/// mutations of the two zero-padded high positions cannot match another fixed 10 bp UMI.
pub(crate) fn neighbors_1mm(u: u32) -> impl Iterator<Item = u32> {
    (0..12u32).flat_map(move |pos| {
        let shift = 2 * pos;
        let cur = (u >> shift) & 0b11;
        (0..4u32).filter_map(move |b| {
            if b == cur {
                None
            } else {
                Some((u & !(0b11 << shift)) | (b << shift))
            }
        })
    })
}

/// Greedy most-abundant-first 1MM collapse over (node, count) pairs; returns surviving node count.
///
/// Each node maps to the canonical survivor of its best (earliest-ranked) 1MM neighbour that
/// outranks it, transitively. These are the same semantics as `correct_umis_1mm`, making the two
/// measurements directly comparable. `merge_ties` additionally lets an equally abundant,
/// earlier-sorting neighbour absorb, matching CellRanger-style 1MM_CR and addressing the 0.9%
/// difference left by the strictly-more-abundant policy.
/// Public alias for cross-module use.
pub(crate) fn neighbors_1mm_pub(u: u32) -> impl Iterator<Item = u32> {
    neighbors_1mm(u)
}

pub(crate) fn collapse_nodes(counts: &FxHashMap<u32, u32>, merge_ties: bool) -> u32 {
    let mut order: Vec<(u32, u32)> = counts.iter().map(|(u, n)| (*u, *n)).collect();
    order.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let rank: FxHashMap<u32, usize> = order.iter().enumerate().map(|(i, (u, _))| (*u, i)).collect();

    let mut canon: FxHashMap<u32, u32> = FxHashMap::default();
    let mut roots = 0u32;
    for (i, (u, n)) in order.iter().enumerate() {
        let mut best: Option<usize> = None;
        for v in neighbors_1mm(*u) {
            if let Some(&r) = rank.get(&v) {
                if r < i {
                    let cn = order[r].1;
                    if cn > *n || (merge_ties && cn == *n) {
                        best = Some(best.map_or(r, |b| b.min(r)));
                    }
                }
            }
        }
        match best {
            Some(r) => {
                let target = order[r].0;
                let c = *canon.get(&target).expect("earlier node always has a canonical");
                canon.insert(*u, c);
            }
            None => {
                canon.insert(*u, *u);
                roots += 1;
            }
        }
    }
    roots
}

pub fn run(args: Args) -> Result<()> {
    let cells_filter: Option<FxHashSet<String>> = match &args.cells {
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

    let mut cell_ids: FxHashMap<String, u32> = FxHashMap::default();
    let mut gene_ids: FxHashMap<String, u32> = FxHashMap::default();
    let mut cell_in_filter: Vec<bool> = Vec::new();

    // G1 oracle: distinct corrected UB per (cell, gene). Only populated when the BAM carries GX/UB.
    let mut g1: FxHashMap<(u32, u32), FxHashSet<u32>> = FxHashMap::default();

    let mut reads: Vec<Read> = Vec::new();
    let mut rec = bam::Record::default();
    let (mut n_total, mut n_used) = (0u64, 0u64);
    let mut have_oracle = false;

    while reader.read_record(&mut rec)? != 0 {
        let flags = Flags::from(rec.flags().bits());
        if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped() {
            continue;
        }
        n_total += 1;
        let t = scan_tags(&rec)?;

        // Cell identity: corrected CB on the oracle BAM; raw CR otherwise (the ingest path corrects
        // barcodes itself in `build`, but for graph statistics the raw barcode is an acceptable
        // cell key — an uncorrected barcode only fragments a cell, never merges two).
        let cell_key = match (&t.cb, &t.cr) {
            (Some(c), _) if c != "-" => c.clone(),
            (_, Some(r)) => format!("{r:08x}"),
            _ => continue,
        };
        let Some(ur) = t.ur else { continue };

        let gene = match &t.gx {
            Some(g) if !g.is_empty() && g != "-" && !g.contains(',') => {
                have_oracle = true;
                let next = gene_ids.len() as u32;
                *gene_ids.entry(g.clone()).or_insert(next)
            }
            _ => u32::MAX,
        };

        let next = cell_ids.len() as u32;
        let cell = *cell_ids.entry(cell_key.clone()).or_insert_with(|| {
            cell_in_filter.push(
                cells_filter.as_ref().is_none_or(|s| s.contains(&cell_key)),
            );
            next
        });

        if gene != u32::MAX {
            if let Some(ub) = t.ub {
                g1.entry((cell, gene)).or_default().insert(ub);
            }
        }

        let chrom = rec.reference_sequence_id().transpose()?.unwrap_or(0) as i32;
        let pos = rec
            .alignment_start()
            .transpose()?
            .map(|p| usize::from(p) as u32 - 1)
            .unwrap_or(0);
        n_used += 1;
        reads.push(Read { cell, chrom, strand_rev: flags.is_reverse_complemented(), pos, umi: ur, gene });
    }
    eprintln!(
        "reads used: {n_used} / {n_total}   cells: {}   reference tags: {have_oracle}",
        cell_ids.len()
    );

    // ---- Base molecules: exact UMI within a locus, no 1MM at ingest. ----
    reads.sort_unstable_by_key(|r| (r.cell, r.chrom, r.strand_rev, r.pos));

    // Per (cell, chrom, strand) group: molecules and the graph over them.
    let windows: [u32; 5] = [10_000, 100_000, 1_000_000, 3_000_000, u32::MAX];
    let mut same_links_by_w = [0u64; 5];
    let mut edges_1mm_by_w = [0u64; 5];
    let mut n_molecules = 0u64;
    let mut n_classes = 0u64;
    let mut cells_with_any_edge: FxHashSet<u32> = FxHashSet::default();

    // Encoded stored-graph stream (window = args.store_window): per group, links and edges as
    // varint-coded molecule-index pairs.
    let mut graph_stream: Vec<u8> = Vec::new();
    let mut stored_links = 0u64;
    let mut stored_edges = 0u64;

    // Replay accumulator: (cell, gene) -> umi -> reads. Aggregating by value IS the exact-collision
    // part of the graph: within (cell, chrom, strand) an exact class is one value, and a gene lives
    // on one (chrom, strand), so this consults nothing the stored graph would not provide.
    let mut replay: FxHashMap<(u32, u32), FxHashMap<u32, u32>> = FxHashMap::default();
    // (cell, umi) -> (gene -> reads). STARsolo's --soloUMIfiltering MultiGeneUMI_CR keeps a
    // multi-gene UMI only for its best-supported gene; replaying that needs this cross-gene view,
    // which is a property of the molecule set, not of the UMI *values*.
    let mut umi_genes: FxHashMap<(u32, u32), FxHashMap<u32, u32>> = FxHashMap::default();

    let mut i = 0usize;
    while i < reads.len() {
        let (cell, chrom, strand) = (reads[i].cell, reads[i].chrom, reads[i].strand_rev);
        let mut j = i;
        while j < reads.len()
            && reads[j].cell == cell
            && reads[j].chrom == chrom
            && reads[j].strand_rev == strand
        {
            j += 1;
        }
        let group = &reads[i..j];
        i = j;

        // Loci by single linkage, then molecules by exact UMI within each locus.
        let mut mols: Vec<Mol> = Vec::new();
        let mut ls = 0usize;
        for k in 1..=group.len() {
            let split = k == group.len() || group[k].pos - group[k - 1].pos > args.locus_gap;
            if !split {
                continue;
            }
            let locus = &group[ls..k];
            ls = k;
            let mut by_umi: FxHashMap<u32, usize> = FxHashMap::default();
            for r in locus {
                let idx = *by_umi.entry(r.umi).or_insert_with(|| {
                    mols.push(Mol { umi: r.umi, pos: r.pos, genes: Vec::new(), n_reads: 0 });
                    mols.len() - 1
                });
                let m = &mut mols[idx];
                m.n_reads += 1;
                if r.gene != u32::MAX {
                    match m.genes.iter_mut().find(|(g, _)| *g == r.gene) {
                        Some((_, c)) => *c += 1,
                        None => m.genes.push((r.gene, 1)),
                    }
                    *replay
                        .entry((cell, r.gene))
                        .or_default()
                        .entry(r.umi)
                        .or_insert(0) += 1;
                    *umi_genes
                        .entry((cell, r.umi))
                        .or_default()
                        .entry(r.gene)
                        .or_insert(0) += 1;
                }
            }
        }
        n_molecules += mols.len() as u64;

        // Same-value links: chain each value's molecules in position order.
        let mut by_value: FxHashMap<u32, Vec<usize>> = FxHashMap::default();
        for (idx, m) in mols.iter().enumerate() {
            by_value.entry(m.umi).or_default().push(idx);
        }
        n_classes += by_value.len() as u64;
        for idxs in by_value.values() {
            for w in idxs.windows(2) {
                let d = mols[w[1]].pos - mols[w[0]].pos;
                for (wi, win) in windows.iter().enumerate() {
                    if d <= *win {
                        same_links_by_w[wi] += 1;
                    }
                }
                cells_with_any_edge.insert(cell);
                if d <= args.store_window {
                    stored_links += 1;
                    put_varint(&mut graph_stream, w[0] as u64);
                    put_varint(&mut graph_stream, (w[1] - w[0]) as u64);
                }
            }
        }

        // 1MM edges between value classes: one edge per value pair, at the minimum inter-molecule
        // distance (which decides the window it falls in).
        for (v, idxs) in &by_value {
            for nb in neighbors_1mm(*v) {
                if nb <= *v {
                    continue; // count each unordered pair once
                }
                let Some(nidxs) = by_value.get(&nb) else { continue };
                let mut dmin = u32::MAX;
                for a in idxs {
                    for b in nidxs {
                        dmin = dmin.min(mols[*a].pos.abs_diff(mols[*b].pos));
                    }
                }
                for (wi, win) in windows.iter().enumerate() {
                    if dmin <= *win {
                        edges_1mm_by_w[wi] += 1;
                    }
                }
                cells_with_any_edge.insert(cell);
                if dmin <= args.store_window {
                    stored_edges += 1;
                    put_varint(&mut graph_stream, idxs[0] as u64);
                    put_varint(&mut graph_stream, nidxs[0] as u64);
                }
            }
        }
    }

    let mut enc = zstd::Encoder::new(Vec::new(), 19)?;
    enc.write_all(&graph_stream)?;
    let graph_bytes = enc.finish()?.len() as u64;

    println!("=== UMI adjacency graph ===");
    println!("base molecules (exact-UMI, {} bp loci): {n_molecules}", args.locus_gap);
    println!("distinct (cell,chrom,strand,value) classes: {n_classes}");
    println!("cells touched by any edge: {}", cells_with_any_edge.len());
    println!();
    println!("{:>12} {:>16} {:>16}", "window", "same-value links", "1MM edges");
    for (wi, w) in windows.iter().enumerate() {
        let wname = if *w == u32::MAX { "chrom-wide".into() } else { format!("{w}") };
        println!("{wname:>12} {:>16} {:>16}", same_links_by_w[wi], edges_1mm_by_w[wi]);
    }
    println!();
    println!(
        "stored graph @ {} bp: {} links + {} edges = {} B zstd  ({:.4} bits/molecule; UMI stream was 27.00)",
        args.store_window,
        stored_links,
        stored_edges,
        graph_bytes,
        8.0 * graph_bytes as f64 / n_molecules.max(1) as f64
    );

    let mut out = serde_json::Map::new();
    out.insert("base_molecules".into(), serde_json::json!(n_molecules));
    out.insert("classes".into(), serde_json::json!(n_classes));
    out.insert("graph_bytes_zstd".into(), serde_json::json!(graph_bytes));
    out.insert("stored_links".into(), serde_json::json!(stored_links));
    out.insert("stored_edges".into(), serde_json::json!(stored_edges));
    out.insert("store_window".into(), serde_json::json!(args.store_window));
    for (wi, w) in windows.iter().enumerate() {
        out.insert(format!("same_links_w{w}"), serde_json::json!(same_links_by_w[wi]));
        out.insert(format!("edges_1mm_w{w}"), serde_json::json!(edges_1mm_by_w[wi]));
    }

    // ---- Sufficiency: replay per-(cell, gene) collapse from the graph, diff against G1. ----
    if have_oracle {
        // Multi-gene UMI statistics, and the filtered replay view.
        let mut multi_gene_umis = 0u64;
        let mut best_gene: FxHashMap<(u32, u32), u32> = FxHashMap::default();
        for ((cell, u), genes) in &umi_genes {
            if genes.len() > 1 {
                multi_gene_umis += 1;
            }
            let best = genes.iter().max_by_key(|(g, c)| (**c, std::cmp::Reverse(**g))).map(|(g, _)| *g).unwrap();
            best_gene.insert((*cell, *u), best);
        }
        println!();
        println!("(cell,UMI) pairs spanning >1 gene: {multi_gene_umis} of {}", umi_genes.len());
        out.insert("multi_gene_umi_pairs".into(), serde_json::json!(multi_gene_umis));
        out.insert("total_cell_umi_pairs".into(), serde_json::json!(umi_genes.len()));

        let mut replay_filtered: FxHashMap<(u32, u32), FxHashMap<u32, u32>> = FxHashMap::default();
        for ((cell, gene), counts) in &replay {
            for (u, n) in counts {
                if best_gene.get(&(*cell, *u)) == Some(gene) {
                    replay_filtered.entry((*cell, *gene)).or_default().insert(*u, *n);
                }
            }
        }

        for (multigene_filter, replay) in [(false, &replay), (true, &replay_filtered)] {
        for merge_ties in [false, true] {
            let mut l1 = 0u64;
            let mut g1_total = 0u64;
            let mut replay_total = 0u64;
            let mut per_cell_g1: FxHashMap<u32, u64> = FxHashMap::default();
            let mut per_cell_err: FxHashMap<u32, u64> = FxHashMap::default();

            for ((cell, gene), ubs) in &g1 {
                let a = ubs.len() as u64;
                let b = match replay.get(&(*cell, *gene)) {
                    Some(counts) => collapse_nodes(counts, merge_ties) as u64,
                    None => 0,
                };
                g1_total += a;
                replay_total += b;
                l1 += a.abs_diff(b);
                *per_cell_g1.entry(*cell).or_insert(0) += a;
                *per_cell_err.entry(*cell).or_insert(0) += a.abs_diff(b);
            }
            for ((cell, gene), counts) in replay.iter() {
                if !g1.contains_key(&(*cell, *gene)) {
                    let b = collapse_nodes(counts, merge_ties) as u64;
                    replay_total += b;
                    l1 += b;
                    *per_cell_err.entry(*cell).or_insert(0) += b;
                }
            }

            let mut rels: Vec<f64> = per_cell_g1
                .iter()
                .filter(|(c, &t)| t > 0 && cell_in_filter[**c as usize])
                .map(|(c, &t)| *per_cell_err.get(c).unwrap_or(&0) as f64 / t as f64)
                .collect();
            rels.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let pick = |q: f64| if rels.is_empty() { 0.0 } else { rels[((rels.len() - 1) as f64 * q) as usize] };
            let policy = if merge_ties { "merge_ties (CR-like)" } else { "strict (G-D greedy)" };
            let filt = if multigene_filter { " + MultiGeneUMI_CR" } else { "" };

            println!();
            println!("--- graph replay vs STARsolo collapse [{policy}{filt}] ---");
            println!("G1 molecules:     {g1_total}");
            println!("replay molecules: {replay_total}");
            println!("overall L1 / G1:  {:.4}%", 100.0 * l1 as f64 / g1_total.max(1) as f64);
            println!("median per-cell:  {:.4}%", 100.0 * pick(0.5));
            println!("p90 per-cell:     {:.4}%", 100.0 * pick(0.9));

            let key = format!("{}{}", if merge_ties { "ties" } else { "strict" },
                              if multigene_filter { "_mgfilt" } else { "" });
            out.insert(format!("replay_{key}_g1"), serde_json::json!(g1_total));
            out.insert(format!("replay_{key}_ours"), serde_json::json!(replay_total));
            out.insert(format!("replay_{key}_overall_l1"), serde_json::json!(l1 as f64 / g1_total.max(1) as f64));
            out.insert(format!("replay_{key}_median"), serde_json::json!(pick(0.5)));
            out.insert(format!("replay_{key}_p90"), serde_json::json!(pick(0.9)));
        }
        }
    } else {
        println!("(no GX tags — sufficiency test skipped; statistics and bytes only)");
    }

    if let Some(p) = &args.json_out {
        if let Some(d) = p.parent() {
            std::fs::create_dir_all(d)?;
        }
        std::fs::write(p, serde_json::to_string_pretty(&out)?)?;
        println!("wrote {}", p.display());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neighbors_are_exactly_36_distinct_values() {
        let u = umi::pack(b"ACGTACGTACGT").unwrap();
        let n: FxHashSet<u32> = neighbors_1mm(u).collect();
        assert_eq!(n.len(), 36);
        assert!(!n.contains(&u));
    }

    #[test]
    fn all_10bp_neighbours_are_enumerated() {
        let u = umi::pack(b"ACGTACGTAC").unwrap();
        let n: FxHashSet<u32> = neighbors_1mm(u).collect();
        for pos in 0..10 {
            let mut seq = *b"ACGTACGTAC";
            let original = seq[pos];
            for &base in b"ACGT" {
                if base == original {
                    continue;
                }
                seq[pos] = base;
                assert!(n.contains(&umi::pack(&seq).unwrap()));
            }
        }
    }

    #[test]
    fn strict_policy_keeps_equal_count_neighbours_separate() {
        let a = umi::pack(b"AAAAAAAAAAAA").unwrap();
        let b = umi::pack(b"AAAAAAAAAAAC").unwrap();
        let mut c = FxHashMap::default();
        c.insert(a, 5u32);
        c.insert(b, 5u32);
        assert_eq!(collapse_nodes(&c, false), 2);
        assert_eq!(collapse_nodes(&c, true), 1, "tie policy must merge them");
    }

    #[test]
    fn abundant_absorbs_rare_under_both_policies() {
        let a = umi::pack(b"AAAAAAAAAAAA").unwrap();
        let b = umi::pack(b"AAAAAAAAAAAC").unwrap();
        let mut c = FxHashMap::default();
        c.insert(a, 100u32);
        c.insert(b, 1u32);
        assert_eq!(collapse_nodes(&c, false), 1);
        assert_eq!(collapse_nodes(&c, true), 1);
    }

    #[test]
    fn two_substitutions_never_merge() {
        let a = umi::pack(b"AAAAAAAAAAAA").unwrap();
        let b = umi::pack(b"AAAAAAAAAACC").unwrap();
        let mut c = FxHashMap::default();
        c.insert(a, 100u32);
        c.insert(b, 1u32);
        assert_eq!(collapse_nodes(&c, true), 2);
    }

    #[test]
    fn chain_collapses_transitively_toward_most_abundant() {
        // A(100) - B(5) - C(1) in a 1MM chain: B into A, C into B's survivor chain => 1 node.
        let a = umi::pack(b"AAAAAAAAAAAA").unwrap();
        let b = umi::pack(b"AAAAAAAAAAAC").unwrap();
        let cc = umi::pack(b"AAAAAAAAAACC").unwrap();
        let mut c = FxHashMap::default();
        c.insert(a, 100u32);
        c.insert(b, 5u32);
        c.insert(cc, 1u32);
        assert_eq!(collapse_nodes(&c, false), 1);
    }
}
