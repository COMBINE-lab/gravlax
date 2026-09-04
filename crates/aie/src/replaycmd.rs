//! `aie replay` — quantify an annotation-free BAM against a supplied annotation to produce a
//! cell-by-gene matrix.
//!
//! Everything upstream of this command is annotation-free (alignment, barcode correction, molecule
//! grouping). The annotation enters here and only here, through the `anno` compiler:
//!
//!   molecule placement → concordant genes (`anno::assign`, the alignToTranscript port)
//!   → uniqueness rules (one gene, NH == 1)
//!   → per-(cell, gene) UMI collapse using tie-merging 1MM and MultiGeneUMI_CR
//!   → MatrixMarket output shaped like STARsolo's raw matrix for entry-by-entry comparison.
//!
//! This BAM-backed command builds molecules in memory. `aie replay-rows` provides the corresponding
//! path for compact `.aie` archives.

use crate::build::{correct_barcode, load_whitelist, molecules_for_chrom, placement_from_parts, to_ops, BcCorrector, Read};
use crate::graph::collapse_nodes;
use anyhow::{bail, Context, Result};
use clap::Parser;
use evidence_io::{umi, Molecule};
use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Value;
use noodles_sam::alignment::record::Flags;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser)]
pub struct Args {
    /// Annotation-free ingest BAM (coordinate-sorted, CR/UR tags).
    pub bam: PathBuf,
    /// GTF to replay — the annotation the archive has never seen.
    #[arg(long)]
    pub gtf: PathBuf,
    /// 10x barcode whitelist.
    #[arg(long)]
    pub whitelist: PathBuf,
    /// Barcode list defining output column order. Pass the reference run's raw `barcodes.tsv` to
    /// compare matrices entry-for-entry.
    #[arg(long)]
    pub barcodes: PathBuf,
    /// Output directory; receives matrix.mtx, features.tsv, barcodes.tsv.
    #[arg(long)]
    pub out_dir: PathBuf,
    #[arg(long, default_value_t = 2_000)]
    pub locus_gap: u32,
    /// Take cell identity from the BAM's corrected `CB` tag instead of correcting `CR` ourselves.
    /// Diagnostic only: reuse corrected `CB` tags from the reference BAM instead of exercising
    /// Gravlax barcode correction. This separates barcode differences from assignment differences.
    #[arg(long)]
    pub use_cb_tag: bool,
    /// Diagnostic: skip 1MM UMI collapse (count distinct raw UMIs).
    #[arg(long)]
    pub no_collapse: bool,
    /// Diagnostic: skip the MultiGeneUMI_CR filter.
    #[arg(long)]
    pub no_multigene_filter: bool,
    /// Classify every READ against the annotation instead of one representative per molecule.
    /// This is STARsolo's own semantics (a UMI counts for a gene if any of its reads is
    /// concordant there) and the upper bound a molecule-level archive must approximate.
    #[arg(long)]
    pub read_level: bool,
    /// Use the simple exact-or-unique-1MM barcode rule instead of the STARsolo-style
    /// pseudocount-posterior corrector. Diagnostic comparison only.
    #[arg(long)]
    pub simple_bc: bool,
    /// With --read-level: keep only ONE unique read per (cell, UMI, locus, junction-chain) — the
    /// most-contained (max start, tie min end). Emulates replay from an archive whose molecule
    /// payload is one representative per chain plus multimapper patterns, pricing the fidelity
    /// cost of not storing every read's exact span.
    #[arg(long)]
    pub chain_representative: bool,
    /// With --chain-representative: also keep the most-EXTENDED read of each chain when distinct
    /// from the most-contained one. The two span extremes bracket the containment behaviour of
    /// every read in between, at the cost of one extra placement for span-diverse chains.
    #[arg(long)]
    pub two_reps: bool,
    /// Skip multi-locus (NH>1) reads entirely. Default is to count a multimapper whose alignments
    /// all resolve to ONE gene, which is what STARsolo's `Gene` feature actually does ("unique
    /// gene", not "unique alignment") — verified on EEF1A1, where a single cell shows 106 distinct
    /// UMIs across 123 mostly-secondary alignments, a matrix count of 102, and primary-GX tags of
    /// just 4. On the evaluated PBMC sample this behavior accounts for about 2.2% of total matrix
    /// mass.
    #[arg(long)]
    pub no_multimappers: bool,
}

pub fn run(args: Args) -> Result<()> {
    let t0 = std::time::Instant::now();
    let anno = anno::Annotation::from_path(&args.gtf)?;
    eprintln!(
        "annotation compiled: {} genes, {} transcripts ({:.1}s)",
        anno.gene_ids.len(),
        anno.transcripts.len(),
        t0.elapsed().as_secs_f32()
    );

    let wl = load_whitelist(&args.whitelist)?;

    // Output barcode order, and packed-barcode → column lookup.
    let bc_text = std::fs::read_to_string(&args.barcodes)?;
    let out_barcodes: Vec<&str> = bc_text.lines().map(|l| l.trim()).filter(|l| !l.is_empty()).collect();
    let mut bc_col: FxHashMap<u32, u32> = FxHashMap::default();
    for (i, b) in out_barcodes.iter().enumerate() {
        if let Some(p) = umi::pack(b.as_bytes()) {
            bc_col.insert(p, i as u32);
        }
    }

    // Pass 1 for the pseudocount corrector: exact-whitelist frequencies of raw barcodes across the
    // run. Annotation-free by construction; skipped when the oracle's CB tags are used directly.
    let corrector: Option<BcCorrector> = if !args.use_cb_tag && !args.simple_bc {
        let t = std::time::Instant::now();
        let mut counts: FxHashMap<u32, u32> = FxHashMap::default();
        let mut r0 = bam::io::reader::Builder
            .build_from_path(&args.bam)
            .with_context(|| format!("opening {}", args.bam.display()))?;
        r0.read_header()?;
        let mut rec0 = bam::Record::default();
        while r0.read_record(&mut rec0)? != 0 {
            let f = Flags::from(rec0.flags().bits());
            if f.is_secondary() || f.is_supplementary() {
                continue;
            }
            for field in rec0.data().iter() {
                let (tag, value) = field?;
                if <[u8; 2]>::from(tag) == *b"CR" {
                    if let Value::String(sv) = value {
                        if let Some(pk) = umi::pack(sv.as_ref()) {
                            if wl.contains(&pk) {
                                *counts.entry(pk).or_insert(0) += 1;
                            }
                        }
                    }
                    break;
                }
            }
        }
        eprintln!("barcode frequency pass: {} exact whitelist barcodes ({:.1}s)", counts.len(), t.elapsed().as_secs_f32());
        Some(BcCorrector::new(wl.clone(), counts))
    } else {
        None
    };

    let mut reader = bam::io::reader::Builder
        .build_from_path(&args.bam)
        .with_context(|| format!("opening {}", args.bam.display()))?;
    let header = reader.read_header()?;
    // BAM reference id → annotation chromosome id, joined by name. Never assume the orders agree.
    let bam2anno: Vec<Option<u32>> = header
        .reference_sequences()
        .keys()
        .map(|name| anno.chrom_ids.get(std::str::from_utf8(name.as_ref()).unwrap_or("")).copied())
        .collect();

    // Accumulators, filled chromosome by chromosome so peak memory stays bounded.
    let mut per_cg: FxHashMap<(u32, u32), FxHashMap<u32, u32>> = FxHashMap::default();
    let mut umi_genes: FxHashMap<(u32, u32), FxHashMap<u32, u32>> = FxHashMap::default();
    let (mut n_reads_used, mut n_mol, mut n_assigned, mut n_multigene, mut n_multimap) =
        (0u64, 0u64, 0u64, 0u64, 0u64);

    let mut rec = bam::Record::default();
    let mut cur_chrom: i32 = -2;
    let mut buf: Vec<Read> = Vec::new();
    let mut bc_cache: FxHashMap<u32, Option<u32>> = FxHashMap::default();
    // NH>1 reads, keyed by read-name hash across ALL their alignment records (secondaries appear
    // on other chromosomes, so this cannot live inside the per-chromosome flush). gene == u32::MAX
    // marks "some alignment was multi-gene", which poisons the read.
    let mut mm: FxHashMap<u64, (u32, u32, u32, bool)> = FxHashMap::default();
    let t_align = std::time::Instant::now();

    fn name_hash(name: &[u8]) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = rustc_hash::FxHasher::default();
        name.hash(&mut h);
        h.finish()
    }

    let flush = |chrom: i32,
                     buf: &mut Vec<Read>,
                     per_cg: &mut FxHashMap<(u32, u32), FxHashMap<u32, u32>>,
                     umi_genes: &mut FxHashMap<(u32, u32), FxHashMap<u32, u32>>,
                     n_mol: &mut u64,
                     n_assigned: &mut u64,
                     n_multigene: &mut u64,
                     n_multimap: &mut u64|
     -> Result<()> {
        if buf.is_empty() {
            return Ok(());
        }
        let Some(anno_chrom) = bam2anno.get(chrom as usize).copied().flatten() else {
            // Chromosome absent from this annotation (e.g. a scaffold): molecules exist but no
            // gene can claim them under this GTF.
            buf.clear();
            return Ok(());
        };
        // The annotation-dependent step, parallelised: placement → concordant gene set.
        let assigned: Vec<(u32, u32, u32, u32)> = if args.read_level {
            // STARsolo's own per-read semantics — the ceiling a molecule archive approximates.
            let mut reads = std::mem::take(buf);
            if args.chain_representative {
                // Reduce to one representative per (cell, strand, locus, UMI, junction-chain).
                use std::hash::{Hash, Hasher};
                reads.sort_by_key(|r| (r.cb, r.strand_rev, r.pos));
                let mut kept: Vec<Read> = Vec::with_capacity(reads.len() / 2);
                let mut i = 0usize;
                while i < reads.len() {
                    let (cb, st) = (reads[i].cb, reads[i].strand_rev);
                    let mut j = i;
                    while j < reads.len() && reads[j].cb == cb && reads[j].strand_rev == st {
                        j += 1;
                    }
                    let grp = &reads[i..j];
                    let mut ls = 0usize;
                    for k in 1..=grp.len() {
                        let split = k == grp.len() || grp[k].pos - grp[k - 1].pos > args.locus_gap;
                        if !split {
                            continue;
                        }
                        let locus = &grp[ls..k];
                        ls = k;
                        // (umi, chain) -> (most-contained idx, most-extended idx, chain read count).
                        let mut best: FxHashMap<(u32, u64), (usize, usize, u32)> = FxHashMap::default();
                        for (li, r) in locus.iter().enumerate() {
                            let mut h = rustc_hash::FxHasher::default();
                            let pl = crate::build::placement_of(0, r);
                            pl.junctions.hash(&mut h);
                            let key = (r.umi, h.finish());
                            let end = pl.end();
                            let e = best.entry(key).or_insert((li, li, 0));
                            e.2 += 1;
                            let cin = &locus[e.0];
                            let cin_end = crate::build::placement_of(0, cin).end();
                            if (r.pos, std::cmp::Reverse(end)) > (cin.pos, std::cmp::Reverse(cin_end)) {
                                e.0 = li;
                            }
                            let cex = &locus[e.1];
                            let cex_end = crate::build::placement_of(0, cex).end();
                            if (r.pos, std::cmp::Reverse(end)) < (cex.pos, std::cmp::Reverse(cex_end)) {
                                e.1 = li;
                            }
                        }
                        let mut idxs: Vec<(usize, usize, u32)> = best.values().copied().collect();
                        idxs.sort_unstable();
                        for (ci, ei, w) in idxs {
                            for (n, li) in [(0usize, ci), (1, ei)] {
                                if n == 1 && (!args.two_reps || li == ci) {
                                    continue;
                                }
                                let r = &locus[li];
                                kept.push(Read {
                                    cb: r.cb, umi: r.umi, pos: r.pos, strand_rev: r.strand_rev,
                                    ops: r.ops.clone(), nm: r.nm,
                                    // The chain's read count rides in `score`, which gene
                                    // assignment never reads; the archive stores it explicitly.
                                    score: w as i32, nh: r.nh,
                                });
                            }
                        }
                    }
                    i = j;
                }
                reads = kept;
            }
            let reads = reads;
            reads
                .par_iter()
                .filter_map(|r| {
                    if r.nh != 1 {
                        return None;
                    }
                    let w = if args.chain_representative { (r.score.max(1)) as u32 } else { 1 };
                    let p = crate::build::placement_of(chrom as u32, r);
                    let genes = anno::assign::concordant_genes(&p, &anno, anno_chrom);
                    match genes.as_slice() {
                        [g] => Some((r.cb, *g, r.umi, w)),
                        [] => None,
                        _ => Some((r.cb, u32::MAX, r.umi, w)),
                    }
                })
                .collect()
        } else {
            let mut mols: Vec<Molecule> = Vec::new();
            molecules_for_chrom(chrom as u32, std::mem::take(buf), args.locus_gap, false, &mut mols);
            *n_mol += mols.len() as u64;
            mols.par_iter()
                .filter_map(|m| {
                    let p = m.placements.first()?;
                    if p.nh != 1 {
                        return None; // STARsolo Gene counts unique alignments only
                    }
                    let genes = anno::assign::concordant_genes(p, &anno, anno_chrom);
                    match genes.as_slice() {
                        [g] => Some((m.cb, *g, m.umi, m.n_reads)),
                        [] => None,
                        _ => Some((m.cb, u32::MAX, m.umi, m.n_reads)), // multi-gene marker
                    }
                })
                .collect()
        };

        for (cb, gene, u, reads) in assigned {
            if gene == u32::MAX {
                *n_multigene += 1;
                continue;
            }
            *n_assigned += 1;
            per_cg
                .entry((cb, gene))
                .or_default()
                .entry(u)
                .and_modify(|r| *r += reads)
                .or_insert(reads);
            umi_genes
                .entry((cb, u))
                .or_default()
                .entry(gene)
                .and_modify(|r| *r += reads)
                .or_insert(reads);
        }
        let _ = n_multimap;
        Ok(())
    };

    loop {
        let more = reader.read_record(&mut rec)? != 0;
        let this_chrom = if more {
            let f = Flags::from(rec.flags().bits());
            if f.is_supplementary() || f.is_unmapped() {
                continue; // secondaries flow through: multimapper evidence lives on them
            }
            rec.reference_sequence_id().transpose()?.unwrap_or(0) as i32
        } else {
            -1
        };
        if this_chrom != cur_chrom && cur_chrom >= 0 {
            flush(cur_chrom, &mut buf, &mut per_cg, &mut umi_genes, &mut n_mol, &mut n_assigned, &mut n_multigene, &mut n_multimap)?;
        }
        if !more {
            break;
        }
        cur_chrom = this_chrom;

        let (mut cr, mut cr_bytes, mut cy_bytes, mut ur, mut cbtag, mut nm, mut score, mut nh) =
            (None, None::<Vec<u8>>, None::<Vec<u8>>, None, None, 0u16, 0i32, 1u16);
        for field in rec.data().iter() {
            let (tag, value) = field?;
            let key = <[u8; 2]>::from(tag);
            match value {
                Value::String(s) => match &key {
                    b"CR" => {
                        cr = umi::pack(s.as_ref());
                        cr_bytes = Some(<[u8] as ToOwned>::to_owned(s.as_ref()));
                    }
                    b"CY" => cy_bytes = Some(<[u8] as ToOwned>::to_owned(s.as_ref())),
                    b"UR" => ur = umi::pack(s.as_ref()),
                    b"CB" => cbtag = umi::pack(s.as_ref()),
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
        let Some(ur) = ur else { continue };
        let cb = if args.use_cb_tag {
            match cbtag {
                Some(c) => c,
                None => continue,
            }
        } else if let Some(corr) = &corrector {
            // Exact hits shortcut the posterior machinery; everything else is scored per read
            // with its own base qualities.
            match (cr, &cr_bytes) {
                (Some(pk), _) if corr.wl_contains(pk) => pk,
                (_, Some(bytes)) => match corr.correct(bytes, cy_bytes.as_deref()) {
                    Some(c) => c,
                    None => continue,
                },
                _ => continue,
            }
        } else {
            let Some(cr) = cr else { continue };
            match *bc_cache.entry(cr).or_insert_with(|| correct_barcode(cr, &wl)) {
                Some(c) => c,
                None => continue,
            }
        };
        let flags = Flags::from(rec.flags().bits());
        let pos = rec.alignment_start().transpose()?.map(|p| usize::from(p) as u32 - 1).unwrap_or(0);

        if nh > 1 {
            if args.no_multimappers {
                continue;
            }
            // Classify this alignment now and fold it into the read's running gene union.
            let anno_chrom = bam2anno.get(this_chrom as usize).copied().flatten();
            let genes: Vec<u32> = match anno_chrom {
                Some(ac) => {
                    let p = placement_from_parts(
                        this_chrom as u32, pos, flags.is_reverse_complemented(),
                        &to_ops(rec.cigar().iter())?, nm, score, nh,
                    );
                    anno::assign::concordant_genes(&p, &anno, ac)
                }
                None => Vec::new(),
            };
            let key = name_hash(rec.name().map(|n| n.as_ref()).unwrap_or(b""));
            let e = mm.entry(key).or_insert((cb, ur, u32::MAX - 1, false)); // MAX-1 = none yet
            for g in genes {
                if e.2 == u32::MAX - 1 || e.2 == g {
                    e.2 = g;
                } else {
                    e.3 = true; // a second distinct gene: poisoned
                }
            }
            continue;
        }

        if flags.is_secondary() {
            continue; // NH==1 has no secondaries; anything else here is malformed
        }
        n_reads_used += 1;
        buf.push(Read {
            cb,
            umi: ur,
            pos,
            strand_rev: flags.is_reverse_complemented(),
            ops: to_ops(rec.cigar().iter())?,
            nm,
            score,
            nh,
        });
    }

    // Resolve multimappers: a read whose alignments all reduced to exactly one gene counts there.
    let mut n_mm_counted = 0u64;
    for (_k, (cb, ur, gene, poisoned)) in mm {
        if poisoned || gene >= u32::MAX - 1 {
            continue;
        }
        n_mm_counted += 1;
        per_cg
            .entry((cb, gene))
            .or_default()
            .entry(ur)
            .and_modify(|r| *r += 1)
            .or_insert(1);
        umi_genes
            .entry((cb, ur))
            .or_default()
            .entry(gene)
            .and_modify(|r| *r += 1)
            .or_insert(1);
    }
    eprintln!("multimapper reads counted (single-gene union): {n_mm_counted}");
    eprintln!(
        "molecules assigned: {n_assigned} of {n_mol} ({n_multigene} multi-gene) from {n_reads_used} reads ({:.1}s)",
        t_align.elapsed().as_secs_f32()
    );

    // MultiGeneUMI_CR: a (cell, UMI) seen in several genes counts only for its best-supported one.
    let t_collapse = std::time::Instant::now();
    let mut best_gene: FxHashMap<(u32, u32), u32> = FxHashMap::default();
    for ((cell, u), genes) in &umi_genes {
        let best = genes
            .iter()
            .max_by_key(|(g, c)| (**c, std::cmp::Reverse(**g)))
            .map(|(g, _)| *g)
            .unwrap();
        best_gene.insert((*cell, *u), best);
    }

    // Per-(cell, gene) collapse, in parallel; policies fixed by the graph experiment.
    let entries: Vec<((u32, u32), FxHashMap<u32, u32>)> = per_cg
        .into_iter()
        .map(|((cell, gene), counts)| {
            let kept: FxHashMap<u32, u32> = if args.no_multigene_filter {
                counts
            } else {
                counts
                    .into_iter()
                    .filter(|(u, _)| best_gene.get(&(cell, *u)) == Some(&gene))
                    .collect()
            };
            ((cell, gene), kept)
        })
        .collect();
    let counted: Vec<((u32, u32), u32)> = entries
        .par_iter()
        .filter_map(|((cell, gene), counts)| {
            if counts.is_empty() {
                return None;
            }
            let n = if args.no_collapse { counts.len() as u32 } else { collapse_nodes(counts, true) };
            Some(((*cell, *gene), n))
        })
        .collect();
    eprintln!("collapse done ({:.1}s)", t_collapse.elapsed().as_secs_f32());

    // Emit the matrix in STARsolo raw shape: rows = features (GTF order), cols = barcodes.
    std::fs::create_dir_all(&args.out_dir)?;
    let mut feat = std::io::BufWriter::new(std::fs::File::create(args.out_dir.join("features.tsv"))?);
    for (id, name) in anno.gene_ids.iter().zip(&anno.gene_names) {
        writeln!(feat, "{id}\t{name}\tGene Expression")?;
    }
    std::fs::write(args.out_dir.join("barcodes.tsv"), out_barcodes.join("\n") + "\n")?;

    let mut cells_seen = 0u64;
    let mut triplets: Vec<(u32, u32, u32)> = Vec::with_capacity(counted.len());
    for ((cell, gene), n) in counted {
        // `cell` is a packed corrected barcode; translate to the output column.
        let Some(col) = bc_col.get(&cell) else { cells_seen += 1; continue };
        triplets.push((gene, *col, n));
    }
    triplets.sort_unstable();
    let mut mtx = std::io::BufWriter::new(std::fs::File::create(args.out_dir.join("matrix.mtx"))?);
    writeln!(mtx, "%%MatrixMarket matrix coordinate integer general")?;
    writeln!(mtx, "%")?;
    writeln!(mtx, "{} {} {}", anno.gene_ids.len(), out_barcodes.len(), triplets.len())?;
    for (g, c, n) in &triplets {
        writeln!(mtx, "{} {} {}", g + 1, c + 1, n)?;
    }
    if cells_seen > 0 {
        bail!("{cells_seen} corrected barcodes were not in the output barcode list — whitelist mismatch");
    }
    let total: u64 = triplets.iter().map(|(_, _, n)| *n as u64).sum();
    eprintln!(
        "replay complete: {} nonzero entries, {} total UMIs -> {} ({:.1}s wall total)",
        triplets.len(),
        total,
        args.out_dir.display(),
        t0.elapsed().as_secs_f32()
    );
    Ok(())
}
