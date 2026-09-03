//! The shared row abstraction both replay paths run on.
//!
//! `extract_rows` turns the ingest BAM into the evidence stored by the archive: each junction chain
//! has two span-extreme representatives and a read count; each multimapper signature has an anchor
//! and an interned paralog pattern; UMIs are represented as **classes and 1MM edges**, never
//! values. `replay_rows` quantifies any annotation from those rows alone.
//!
//! Both the BAM-backed and archive-backed replay commands call these same functions, so the
//! archive regression reduces to: serialize → deserialize → identical rows → identical matrix.
//!
//! Value-free collapse: 1MM collapse needs an ordering among equal-count neighbours. The
//! value-based reference broke ties lexicographically on the UMI string; rows carry no values, so
//! ties break on **class id** (global first-occurrence order in the ingest scan), which both paths
//! share by construction. The substitution has been benchmarked against value-based collapse.

use crate::build::{load_whitelist, to_ops, BcCorrector};
use anyhow::{bail, Context, Result};
use evidence_io::archive::Shape;
use evidence_io::{umi, Block, Junction, Placement, Strand};
use ingest::placement_from_alignment;
use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Value;
use noodles_sam::alignment::record::Flags;
use rayon::prelude::*;
// Hash-map policy: the hot replay pipeline is hash-free (sorted runs + CSR adjacency); the maps
// that remain here are rustc-hash FxHashMap — whose iteration order is deterministic run-to-run —
// because their iteration order leaks into observables: EM float accumulation order (targets are
// built by map iteration and f64 addition is not associative) and the ingest-side interning that
// fixes archive bytes. hashbrown's default foldhash seeds per-process, so putting it on these
// paths would make outputs nondeterministic, not just different.
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::fs::File;
use std::io::{BufWriter, Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

pub const EDGE_WINDOW: u32 = 3_000_000;

/// Full-file identity captured from the same held file that supplied a raw replay/ingest parse.
/// Keeping this dependency-neutral lets the archive command wrap it in its report schema without
/// making the hot row representation depend on the output crate.
#[derive(Debug)]
pub(crate) struct ConsumedFileIdentity {
    pub(crate) blake3: String,
    pub(crate) bytes: u64,
}

fn stable_file_metadata(before: &std::fs::Metadata, after: &std::fs::Metadata) -> Result<bool> {
    if !before.is_file()
        || !after.is_file()
        || before.len() != after.len()
        || before.modified()? != after.modified()?
    {
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

/// Hash a held input from byte zero and reject mutation or pathname replacement across the
/// computation that consumed it. The caller captures `before` before parsing the same `File`.
pub(crate) fn identity_of_consumed_file(
    mut file: File,
    before: std::fs::Metadata,
    path: &Path,
) -> Result<ConsumedFileIdentity> {
    file.rewind()
        .with_context(|| format!("rewinding {} for content identity", path.display()))?;
    let mut hasher = blake3::Hasher::new();
    let mut buffer = vec![0u8; 1 << 20];
    let mut bytes = 0u64;
    loop {
        let read = file.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
        bytes = bytes
            .checked_add(read as u64)
            .context("content identity byte count overflow")?;
    }
    let after = file.metadata()?;
    if !stable_file_metadata(&before, &after)? || bytes != after.len() {
        bail!(
            "{} changed while it was consumed and identified",
            path.display()
        );
    }
    let path_after = std::fs::metadata(path)
        .with_context(|| format!("rechecking consumed input {}", path.display()))?;
    if !stable_file_metadata(&after, &path_after)? {
        bail!("{} was replaced while it was consumed", path.display());
    }
    Ok(ConsumedFileIdentity {
        blake3: hasher.finalize().to_hex().to_string(),
        bytes,
    })
}

#[cfg(all(test, unix))]
mod consumed_file_identity_tests {
    use super::identity_of_consumed_file;
    use std::io::Read;

    fn scratch_file(name: &str) -> std::path::PathBuf {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "gravlax-consumed-identity-{}-{nonce}-{name}",
            std::process::id()
        ))
    }

    #[test]
    fn identity_uses_consumed_inode_and_rejects_path_replacement() {
        let path = scratch_file("input");
        let held_path = scratch_file("held");
        std::fs::write(&path, b"consumed bytes").unwrap();
        let mut file = std::fs::File::open(&path).unwrap();
        let before = file.metadata().unwrap();
        let mut consumed = Vec::new();
        file.read_to_end(&mut consumed).unwrap();
        assert_eq!(consumed, b"consumed bytes");

        std::fs::rename(&path, &held_path).unwrap();
        std::fs::write(&path, b"replacement!!").unwrap();
        let error = identity_of_consumed_file(file, before, &path).unwrap_err();
        assert!(
            error.to_string().contains("replaced") || error.to_string().contains("changed"),
            "unexpected stability error: {error:#}"
        );

        std::fs::remove_file(path).ok();
        std::fs::remove_file(held_path).ok();
    }
}

/// One archive row. `pattern == u32::MAX` marks a unique-read chain representative (weight = the
/// chain's read count, carried by both span extremes); otherwise the row is an aggregated
/// multimapper signature whose alternatives are `patterns[pattern]` applied to the anchor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub cell: u32,
    pub umi_class: u32,
    pub chrom: u32,
    pub pos: u32,
    pub strand_rev: bool,
    pub shape: u32,
    pub weight: u32,
    pub pattern: u32,
}

/// One alternative placement of a paralog pattern, relative to the anchor.
/// `shape == u32::MAX` means "same shape as the applying row's anchor" — measured at 78.7% of
/// alternatives, and making it relative also merges patterns that differ only in absolute shape,
/// increasing dictionary sharing.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct PatAlt {
    pub chrom: u32,
    pub offset: i64,
    pub strand_flip: bool,
    pub shape: u32,
}

pub const SAME_SHAPE: u32 = u32::MAX;

/// One junction chain of a molecule: its read count and one or two span-extreme representatives.
/// `reps` is inline (never more than 2 entries by construction), so a chain costs no heap.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolChain {
    pub weight: u32,
    pub reps: SmallVec<[(u32, u32); 2]>, // (pos, shape), 1 or 2 entries
}

/// A molecule-major record: cell and UMI class paid once; children carry only local geometry.
/// The child vectors are SmallVecs sized for the dominant case (one chain, zero-or-one mm
/// signature). On an approximately 100-million-molecule benchmark archive, the previous
/// per-molecule `Vec` pair created about 200 million small heap allocations, dominating archive
/// load time and teardown.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolRec {
    pub cell: u32,
    pub umi_class: u32,
    pub chrom: u32,
    pub strand_rev: bool,
    pub chains: SmallVec<[MolChain; 1]>,
    /// Aggregated multimapper signatures: (anchor pos, anchor shape, pattern id, weight).
    pub mms: SmallVec<[(u32, u32, u32, u32); 1]>,
}

impl MolRec {
    pub fn anchor(&self) -> u32 {
        let c = self.chains.iter().flat_map(|c| c.reps.iter().map(|r| r.0)).min();
        let m = self.mms.iter().map(|m| m.0).min();
        match (c, m) {
            (Some(a), Some(b)) => a.min(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => 0,
        }
    }
}

pub struct Extracted {
    pub mols: Vec<MolRec>,
    /// Undirected 1MM edges between UMI classes (`a < b`), scoped to one cell.
    pub edges: Vec<(u32, u32)>,
    /// Dense cell id → packed 16 bp barcode.
    pub cells: Vec<u32>,
    pub shapes: Vec<Shape>,
    pub patterns: Vec<Vec<PatAlt>>,
    pub n_classes: u32,
    pub chrom_names: Vec<String>,
}

/// Canonical flat rows from molecule records. Both replay paths call this, so any serialization
/// that reproduces the molecule records byte-for-byte reproduces the matrices byte-for-byte.
pub fn flatten(mols: &[MolRec]) -> Vec<Row> {
    let mut rows = flatten_unsorted(mols);
    rows.par_sort_unstable_by_key(|r| {
        (
            r.chrom,
            r.pos,
            r.cell,
            r.umi_class,
            r.pattern,
            r.shape,
            r.weight,
        )
    });
    rows
}

/// Rows in molecule order, without the canonical sort. `replay_rows` uses this directly: its
/// per-row classification is independent and its aggregation sorts full tuples itself, so row
/// order cannot reach the matrix and the big sort is pure cost there. The EM paths keep the
/// sorted [`flatten`]: their candidate-list construction order feeds f64 accumulation, where
/// order is part of the stable, byte-checked output.
pub fn flatten_unsorted(mols: &[MolRec]) -> Vec<Row> {
    mols.par_iter()
        .flat_map_iter(|m| {
            let chains = m.chains.iter().flat_map(move |ch| {
                ch.reps.iter().map(move |(pos, shape)| Row {
                    cell: m.cell, umi_class: m.umi_class, chrom: m.chrom, pos: *pos,
                    strand_rev: m.strand_rev, shape: *shape, weight: ch.weight, pattern: u32::MAX,
                })
            });
            let mms = m.mms.iter().map(move |(pos, shape, pattern, weight)| Row {
                cell: m.cell, umi_class: m.umi_class, chrom: m.chrom, pos: *pos,
                strand_rev: m.strand_rev, shape: *shape, weight: *weight, pattern: *pattern,
            });
            chains.chain(mms)
        })
        .collect()
}

/// STARsolo Velocyto semantics from archive rows (port of `SoloFeature_countVelocyto.cpp` and
/// `alignToTranscriptMinOverlap`). Unique evidence only (chains; multimapper
/// molecules contribute nothing, as STAR streams only nAlignG==1 reads). Per (cell, class):
/// transcript sets intersect across representatives with type-bit OR — and because classes are
/// global (cell, value), a same-cell UMI at two loci intersects to empty and drops, exactly as
/// STAR's cell-scoped `cuTrTypes[iCB][umi]` does. Known evidence substitution, measured not
/// assumed: the two span extremes stand in for "all reads of the UMI".
/// Returns ((cell, gene) → [spliced, unspliced, ambiguous], classes counted).
/// Strand-parameterized velocity replay for chemistries such as 10x 5' R2-only whose cDNA
/// alignments are antisense to the counted transcript.
pub fn velocity_rows_stranded(
    x: &Extracted,
    anno: &anno::Annotation,
    solo_strand: anno::assign::SoloStrand,
) -> (FxHashMap<(u32, u32), ([u32; 3], bool)>, u64) {
    let bam2anno: Vec<Option<u32>> =
        x.chrom_names.iter().map(|n| anno.chrom_ids.get(n).copied()).collect();

    // Per molecule: the rep-intersected (transcript, bits) list, sorted by transcript, plus a
    // completeness flag — a chain storing 2 extremes for >2 reads has missing middles, everything
    // else carries every read (1 rep = identical placements at any weight). VELO-2's stratifier.
    let per_mol: Vec<Option<(u32, u32, Vec<(u32, u8)>, bool)>> = x
        .mols
        .par_iter()
        .map(|m| {
            if m.chains.is_empty() {
                return None;
            }
            let complete = m.chains.iter().all(|c| !(c.reps.len() == 2 && c.weight > 2));
            let ac = (*bam2anno.get(m.chrom as usize)?)?;
            let mut acc: Option<Vec<(u32, u8)>> = None;
            for ch in &m.chains {
                for (pos, shape) in &ch.reps {
                    let p = placement_from_parts(m.chrom, *pos, m.strand_rev, &x.shapes[*shape as usize], 1);
                    let mut cur: Vec<(u32, u8)> = Vec::new();
                    for txi in anno.overlapping(ac, p.start(), p.end()) {
                        let t = &anno.transcripts[txi as usize];
                        if !solo_strand.accepts(m.strand_rev, t.strand_rev) {
                            continue;
                        }
                        if let Some(v) = anno::assign::align_vs_transcript_min_overlap(&p, t) {
                            cur.push((txi, anno::assign::velocyto_bits(v)));
                        }
                    }
                    if cur.is_empty() {
                        continue; // no compatible transcript: the read is absent from STAR's
                                  // velocity stream and must NOT blank the UMI's intersection
                    }
                    cur.sort_unstable_by_key(|(t, _)| *t);
                    acc = Some(match acc {
                        None => cur,
                        Some(old) => intersect_or(&old, &cur),
                    });
                    if acc.as_ref().unwrap().is_empty() {
                        break; // a REAL conflict between non-empty sets kills the UMI, as in STAR
                    }
                }
            }
            // None = no velocity evidence at all (contributes nothing);
            // Some(empty) = conflict-dead (kills the class).
            acc.map(|trs| (m.cell, m.umi_class, trs, complete))
        })
        .collect();

    // Class-level intersection across molecules (cell-scoped by construction of classes).
    let mut per_class: FxHashMap<(u32, u32), (Option<Vec<(u32, u8)>>, bool)> = FxHashMap::default();
    for e in per_mol.into_iter().flatten() {
        let (cell, cls, trs, complete) = e;
        per_class
            .entry((cell, cls))
            .and_modify(|(old, cflag)| {
                let merged = match old.take() {
                    Some(o) => intersect_or(&o, &trs),
                    None => Vec::new(),
                };
                *old = Some(merged);
                *cflag &= complete;
            })
            .or_insert((Some(trs), complete));
    }

    // 1MM UMI correction: STAR velocity consumes the GENE feature's corrected UMIs, so apply
    // exactly the merges Gene collapse performs (gene_canon_map) and nothing else — intronic-only
    // UMIs were never in Gene's collapse and stay uncorrected (the first cut merged them too and
    // lost 10% of unspliced mass; the second cut skipped merging and over-counted spliced).
    let canon = gene_canon_map_stranded(x, anno, solo_strand);
    let mut merged: FxHashMap<(u32, u32), (Option<Vec<(u32, u8)>>, bool)> = FxHashMap::default();
    for ((cell, cls), (trs, complete)) in per_class {
        let Some(trs) = trs else { continue };
        let cls_eff = canon.get(&(cell, cls)).copied().unwrap_or(cls);
        merged
            .entry((cell, cls_eff))
            .and_modify(|(old, cflag)| {
                let m2 = match old.take() {
                    Some(o) => intersect_or(&o, &trs),
                    None => Vec::new(),
                };
                *old = Some(m2);
                *cflag &= complete;
            })
            .or_insert((Some(trs), complete));
    }

    let mut counts: FxHashMap<(u32, u32), ([u32; 3], bool)> = FxHashMap::default();
    let mut n_counted = 0u64;
    for ((cell, _cls), (trs, complete)) in merged {
        let trs = trs.unwrap_or_default();
        if trs.is_empty() {
            continue;
        }
        let gene = anno.transcripts[trs[0].0 as usize].gene;
        if trs.iter().any(|(t, _)| anno.transcripts[*t as usize].gene != gene) {
            continue; // multigene
        }
        use anno::assign::{VB_CONCORDANT, VB_EXONINTRON, VB_INTRON, VB_SPAN};
        let (mut exon_m, mut intron_m, mut mixed_m) = (false, false, false);
        let mut span_m = true;
        for (_, b) in &trs {
            let (i, e, x2, s) = (
                b & VB_INTRON != 0,
                b & VB_CONCORDANT != 0,
                b & VB_EXONINTRON != 0,
                b & VB_SPAN != 0,
            );
            mixed_m |= ((i && e) || x2) && !s;
            span_m &= s;
            exon_m |= e && !i && !x2;
            intron_m |= i && !x2 && !e;
        }
        let comp = if exon_m && !intron_m && !mixed_m {
            0 // spliced
        } else if span_m || ((intron_m || mixed_m) && !exon_m) {
            1 // unspliced
        } else {
            2 // ambiguous
        };
        let e = counts.entry((cell, gene)).or_insert(([0; 3], true));
        e.0[comp] += 1;
        e.1 &= complete;
        n_counted += 1;
    }
    (counts, n_counted)
}

/// EM-0: STARsolo `--soloMultiMappers EM` semantics, ported from `SoloFeature_collapseUMIall.cpp`
/// (lines 270–500, read on this host). Per cell: multi-gene UMIs are classes whose every row is
/// multi-gene (a UMI seen among unique-mapped reads is skipped — our mixed classes — and UMIs are
/// raw values, uncorrected — our classes); a UMI's candidate genes are the INTERSECTION of its
/// reads' gene sets; EM initializes at unique+uniform, zeroes counts < 0.01 each iteration, and
/// stops at max-abs-change < 0.01 or 100 iterations. Returns unique + EM-multimapper totals per
/// (cell, gene) — the UniqueAndMult-EM matrix.
pub fn em_star_matrix(x: &Extracted, anno: &anno::Annotation) -> FxHashMap<(u32, u32), f64> {
    // Unique side = the byte-exact Gene replay.
    let (unique_counts, _, _) = replay_rows(x, anno);

    let bam2anno: Vec<Option<u32>> =
        x.chrom_names.iter().map(|n| anno.chrom_ids.get(n).copied()).collect();
    let rows = flatten(&x.mols);
    let per_row: Vec<Option<(u32, u32, Vec<u32>)>> = rows
        .par_iter()
        .map_init(RowScratch::default, |s, r| {
            row_genes(r, x, anno, &bam2anno, MmMissing::SkipAlt, s)?;
            if s.genes.is_empty() { None } else { Some((r.cell, r.umi_class, s.genes.clone())) }
        })
        .collect();

    // Per (cell, class): has_single kills the UMI for the mult side; else intersect gene sets.
    let mut per_class: FxHashMap<(u32, u32), (bool, Option<Vec<u32>>)> = FxHashMap::default();
    for e in per_row.into_iter().flatten() {
        let (cell, cls, mut genes) = e;
        genes.sort_unstable();
        let ent = per_class.entry((cell, cls)).or_insert((false, None));
        if genes.len() == 1 {
            ent.0 = true;
        } else {
            ent.1 = Some(match ent.1.take() {
                None => genes,
                Some(old) => old.into_iter().filter(|g| genes.binary_search(g).is_ok()).collect(),
            });
        }
    }
    // Cell -> the multi-gene UMIs' candidate sets.
    let mut cell_umis: FxHashMap<u32, Vec<Vec<u32>>> = FxHashMap::default();
    for ((cell, _), (has_single, inter)) in per_class {
        if has_single {
            continue; // skipped: this UMI was among uniquely-mapped
        }
        if let Some(v) = inter {
            // STAR keeps ANY non-empty intersection: a size-1 set assigns its UMI wholly to that
            // gene inside the EM (SoloFeature_collapseUMIall pushes vg regardless of size).
            if !v.is_empty() {
                cell_umis.entry(cell).or_default().push(v);
            }
        }
    }

    let cells_em: Vec<(u32, Vec<Vec<u32>>)> = cell_umis.into_iter().collect();
    let results: Vec<Vec<((u32, u32), f64)>> = cells_em
        .par_iter()
        .map(|(cell, umis)| {
            // genesM: union of candidates, indexed.
            let mut gidx: FxHashMap<u32, usize> = FxHashMap::default();
            let mut glist: Vec<u32> = Vec::new();
            for ug in umis {
                for g in ug {
                    gidx.entry(*g).or_insert_with(|| {
                        glist.push(*g);
                        glist.len() - 1
                    });
                }
            }
            let n = glist.len();
            let mut g_uniform = vec![0.0f64; n];
            for ug in umis {
                for g in ug {
                    g_uniform[gidx[g]] += 1.0 / ug.len() as f64;
                }
            }
            let mut g_u = vec![0.0f64; n];
            for (i, g) in glist.iter().enumerate() {
                g_u[i] = unique_counts.get(&(*cell, *g)).copied().unwrap_or(0) as f64;
            }
            let mut em_old: Vec<f64> = g_uniform.iter().zip(&g_u).map(|(a, b)| a + b).collect();
            let mut em_new = vec![0.0f64; n];
            let idx_sets: Vec<Vec<usize>> =
                umis.iter().map(|ug| ug.iter().map(|g| gidx[g]).collect()).collect();
            for _iter in 0..100 {
                em_new.copy_from_slice(&g_u);
                for v in em_old.iter_mut() {
                    if *v < 0.01 {
                        *v = 0.0;
                    }
                }
                for ug in &idx_sets {
                    let norm: f64 = ug.iter().map(|&i| em_old[i]).sum();
                    if norm == 0.0 {
                        continue;
                    }
                    for &i in ug {
                        em_new[i] += em_old[i] / norm;
                    }
                }
                let max_change = em_new
                    .iter()
                    .zip(&em_old)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0f64, f64::max);
                std::mem::swap(&mut em_old, &mut em_new);
                if max_change < 0.01 {
                    break;
                }
            }
            // em_old now holds the final state: unique + mult. Emit the mult part per gene.
            glist
                .iter()
                .enumerate()
                .map(|(i, g)| ((*cell, *g), (em_old[i] - g_u[i]).max(0.0)))
                .filter(|(_, v)| *v > 1e-9)
                .collect()
        })
        .collect();

    // UniqueAndMult = unique + mult.
    let mut out: FxHashMap<(u32, u32), f64> = FxHashMap::default();
    for ((cell, gene), c) in unique_counts {
        out.insert((cell, gene), c as f64);
    }
    for chunk in results {
        for (k, v) in chunk {
            *out.entry(k).or_insert(0.0) += v;
        }
    }
    out
}

/// Measure how much evidence never reaches the Gene matrix because of multi-gene ambiguity.
/// Buckets every row by its concordant-gene count and every `(cell, class)` by its fate; the
/// "recoverable pool" is what an equivalence-class EM could fractionally assign.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize)]
pub struct MultigeneAuditSummary {
    pub no_gene_read_weight: u64,
    pub single_gene_read_weight: u64,
    pub multi_gene_read_weight: u64,
    pub no_gene_rows: u64,
    pub single_gene_rows: u64,
    pub multi_gene_rows: u64,
    pub gene_evidence_classes: u64,
    pub single_gene_counted_classes: u64,
    pub multi_only_classes: u64,
    pub mixed_classes: u64,
    pub no_gene_classes: u64,
}

impl MultigeneAuditSummary {
    pub fn total_read_weight(&self) -> u64 {
        self.no_gene_read_weight + self.single_gene_read_weight + self.multi_gene_read_weight
    }
}

pub fn multigene_audit_summary_stranded(
    x: &Extracted,
    anno: &anno::Annotation,
    solo_strand: anno::assign::SoloStrand,
) -> MultigeneAuditSummary {
    let bam2anno: Vec<Option<u32>> = x
        .chrom_names
        .iter()
        .map(|n| anno.chrom_ids.get(n).copied())
        .collect();
    let rows = flatten(&x.mols);
    // Per row: number of concordant genes (unique rows) or union genes (mm rows), plus weight.
    let per_row: Vec<(u32, u32, usize, u32)> = rows
        .par_iter()
        .map_init(RowScratch::default, |s, r| {
            let ngenes = match row_genes_stranded(
                r, x, anno, &bam2anno, MmMissing::SkipAlt, solo_strand, s,
            ) {
                Some(()) => s.genes.len(),
                None => 0, // unique row on an unannotated chromosome counted as 0 genes
            };
            (r.cell, r.umi_class, ngenes, r.weight)
        })
        .collect();

    let (mut w0, mut w1, mut wm) = (0u64, 0u64, 0u64);
    let (mut r0, mut r1, mut rm) = (0u64, 0u64, 0u64);
    // Class fate: has_single (>=1 single-gene row), has_multi (>=1 multi-gene row), has_any.
    let mut class_state: FxHashMap<(u32, u32), (bool, bool)> = FxHashMap::default();
    for (cell, cls, ng, w) in &per_row {
        match ng {
            0 => {
                w0 += *w as u64;
                r0 += 1;
            }
            1 => {
                w1 += *w as u64;
                r1 += 1;
            }
            _ => {
                wm += *w as u64;
                rm += 1;
            }
        }
        let e = class_state.entry((*cell, *cls)).or_insert((false, false));
        e.0 |= *ng == 1;
        e.1 |= *ng >= 2;
    }
    let n_classes = class_state.len() as u64;
    let mut cls_single = 0u64; // counted today
    let mut cls_multi_only = 0u64; // gene-informative, fully discarded — the EM pool
    let mut cls_mixed = 0u64; // counted, but multi-gene evidence ignored
    for (_, (s, m)) in &class_state {
        match (s, m) {
            (true, false) => cls_single += 1,
            (false, true) => cls_multi_only += 1,
            (true, true) => cls_mixed += 1,
            (false, false) => {}
        }
    }
    MultigeneAuditSummary {
        no_gene_read_weight: w0,
        single_gene_read_weight: w1,
        multi_gene_read_weight: wm,
        no_gene_rows: r0,
        single_gene_rows: r1,
        multi_gene_rows: rm,
        gene_evidence_classes: n_classes,
        single_gene_counted_classes: cls_single,
        multi_only_classes: cls_multi_only,
        mixed_classes: cls_mixed,
        no_gene_classes: n_classes - cls_single - cls_multi_only - cls_mixed,
    }
}

fn write_multigene_audit<W: Write>(
    writer: &mut W,
    summary: &MultigeneAuditSummary,
) -> std::io::Result<()> {
    let weight = summary.total_read_weight();
    writeln!(writer, "== multigene audit (EM upside bound) ==")?;
    writeln!(
        writer,
        "rows by concordant-gene count (read-weighted): 0 genes {:.2}% | 1 gene {:.2}% | >=2 genes {:.2}%  ({} / {} / {} rows)",
        100.0 * summary.no_gene_read_weight as f64 / weight as f64,
        100.0 * summary.single_gene_read_weight as f64 / weight as f64,
        100.0 * summary.multi_gene_read_weight as f64 / weight as f64,
        summary.no_gene_rows,
        summary.single_gene_rows,
        summary.multi_gene_rows,
    )?;
    writeln!(
        writer,
        "classes with gene evidence: single-gene-counted {} | MULTI-ONLY (fully discarded, the EM pool) {} ({:.2}% of counted) | mixed {}; no-gene classes {}",
        summary.single_gene_counted_classes,
        summary.multi_only_classes,
        100.0 * summary.multi_only_classes as f64
            / (summary.single_gene_counted_classes + summary.mixed_classes).max(1) as f64,
        summary.mixed_classes,
        summary.no_gene_classes,
    )
}

pub fn multigene_audit_stranded(
    x: &Extracted,
    anno: &anno::Annotation,
    solo_strand: anno::assign::SoloStrand,
) {
    let summary = multigene_audit_summary_stranded(x, anno, solo_strand);
    let stdout = std::io::stdout();
    write_multigene_audit(&mut stdout.lock(), &summary).expect("writing multigene audit");
}

// The construction layout below follows the same principle as piscem-infer's PackedEqMap and
// Salmon's PackedEqClasses: equivalence-class labels live in one flat u32 array and u64 CSR
// offsets delimit classes.  The archive is decoded in bounded batches and only two u32 words per
// observed (class, gene) support item survive classification.  Cell shards make the sort/reduce
// working set bounded and give the EM a deterministic reduction order.
const EM_SHARDS: usize = 64;
const EM_IN_MEMORY_SUPPORT_BUDGET: usize = 512 << 20;
const EM_GENE_MASK: u32 = (1 << 30) - 1;
const EM_MULTI_FLAG: u32 = 1 << 30;
const EM_UNIQUE_FLAG: u32 = 1 << 31;
const EM_NO_TRUTH: u32 = u32::MAX;
pub(crate) const EM_NO_GROUP: u32 = u32::MAX;

#[derive(Clone, Copy)]
struct EmSupport {
    class: u32,
    gene_flags: u32,
}

#[derive(Default, Clone, Copy)]
struct EmCounters {
    single: u64,
    mixed: u64,
    multi_only: u64,
    masked: u64,
    truth_lost: u64,
}

impl std::ops::AddAssign for EmCounters {
    fn add_assign(&mut self, rhs: Self) {
        self.single += rhs.single;
        self.mixed += rhs.mixed;
        self.multi_only += rhs.multi_only;
        self.masked += rhs.masked;
        self.truth_lost += rhs.truth_lost;
    }
}

struct PackedEmShard {
    // Sorted packed (cell,gene) universe needed by either the base counts or a target.
    cell_gene_keys: Vec<u64>,
    unique_base: Vec<f64>,
    // Target i has candidates in target_genes/target_local[starts[i]..starts[i+1]].
    starts: Vec<u64>,
    target_cells: Vec<u32>,
    target_genes: Vec<u32>,
    target_local: Vec<u32>,
    target_truth: Vec<u32>,
}

impl PackedEmShard {
    #[inline]
    fn target_range(&self, i: usize) -> std::ops::Range<usize> {
        self.starts[i] as usize..self.starts[i + 1] as usize
    }
}

/// Packed, cell-sharded input to the recovery EM.  It contains no molecule/row table, no hash
/// table, and no per-class or per-target heap allocation.
pub(crate) struct PackedEm {
    shards: Vec<PackedEmShard>,
    global_base: Vec<f64>,
    counters: EmCounters,
    mask_frac: f64,
}

/// Optional query-time hierarchy.  `cell_group` covers the archive cell dictionary, while
/// `eval_cell` independently fixes the population on which every mode is scored.  Keeping those
/// roles separate permits an exact all-archive-cells/one-group collapse control without changing
/// the called-cell denominator.
pub(crate) struct EmGroups {
    pub(crate) cell_group: Vec<u32>,
    pub(crate) eval_cell: Vec<bool>,
    pub(crate) names: Vec<String>,
}

/// Candidate-normalized convex partial-pooling parameters. `global_weight` is derived as the
/// simplex remainder so CLI inputs cannot silently be renormalized. `group_prior` is an effective
/// candidate-set pseudo-count mass, unlike the whole-transcriptome alphas retained by the
/// registered additive evaluator.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ConvexEmParams {
    pub(crate) cell_weight: f64,
    pub(crate) group_weight: f64,
    pub(crate) group_prior: f64,
}

impl ConvexEmParams {
    #[inline]
    fn global_weight(self) -> f64 {
        1.0 - self.cell_weight - self.group_weight
    }
}

/// Posterior-mean screening approximation to the reserved hierarchical Dirichlet model. The
/// group posterior borrows `group_prior` candidate-level counts from the sample, and the cell
/// posterior borrows `cell_prior` counts from its leave-one-cell-out group posterior. This makes
/// pooling strength adapt to the evidence depth of each candidate set without fitting latent
/// concentration parameters.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DirichletProxyParams {
    pub(crate) cell_prior: f64,
    pub(crate) group_prior: f64,
}

/// Monotone weight blending the fixed-convex and posterior-mean Dirichlet-proxy probabilities.
/// The weight uses mode-independent fitted unique evidence over the target candidate set. A
/// squared Hill curve keeps shallow targets close to fixed convex while approaching the proxy for
/// deep targets without adding a discontinuous threshold.
#[derive(Clone, Copy, Debug)]
pub(crate) struct DepthHybridParams {
    pub(crate) depth_scale: f64,
    pub(crate) depth_power: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct PackedDepthMetrics {
    pub(crate) name: &'static str,
    pub(crate) top1: u64,
    pub(crate) expected: f64,
    pub(crate) negative_log_loss: f64,
    pub(crate) brier: f64,
    pub(crate) n: u64,
}

#[derive(Clone, Debug)]
pub(crate) struct PackedModeMetrics {
    pub(crate) name: String,
    pub(crate) top1: u64,
    pub(crate) expected: f64,
    pub(crate) negative_log_loss: f64,
    pub(crate) brier: f64,
    pub(crate) n: u64,
    pub(crate) calibration: [[u64; 2]; 10],
    pub(crate) evidence_depth: Vec<PackedDepthMetrics>,
}

#[derive(Clone, Debug)]
struct PackedTargetObservation {
    cell: u32,
    negative_log_loss: f64,
    brier: f64,
    correct: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct PackedPairedMetrics {
    pub(crate) candidate: String,
    pub(crate) reference: String,
    pub(crate) n: u64,
    pub(crate) cells: u64,
    pub(crate) mean_negative_log_loss_difference: f64,
    pub(crate) negative_log_loss_clustered_se: f64,
    pub(crate) negative_log_loss_ci95: [f64; 2],
    pub(crate) negative_log_loss_wins: u64,
    pub(crate) mean_brier_difference: f64,
    pub(crate) brier_clustered_se: f64,
    pub(crate) brier_ci95: [f64; 2],
    pub(crate) brier_wins: u64,
    pub(crate) top1_percentage_point_difference: f64,
}

/// Streaming builder for [`PackedEm`].  Archive batches may be discarded immediately after
/// `add_archive_chunks`; only compact support words remain.
pub(crate) struct PackedEmAccumulator {
    n_classes: u32,
    gene_count: usize,
    bam2anno: Vec<Option<u32>>,
    shards: Vec<Vec<EmSupport>>,
    spill: Option<EmSupportSpill>,
}

struct EmSupportSpill {
    dir: PathBuf,
    writers: Vec<Option<BufWriter<File>>>,
    counts: Vec<usize>,
}

impl EmSupportSpill {
    fn new() -> Result<Self> {
        let base = std::env::temp_dir();
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .context("system clock precedes Unix epoch")?
            .as_nanos();
        let mut dir = None;
        for attempt in 0..100u32 {
            let candidate = base.join(format!(
                "gravlax-em-{}-{nonce}-{attempt}",
                std::process::id()
            ));
            match std::fs::create_dir(&candidate) {
                Ok(()) => {
                    dir = Some(candidate);
                    break;
                }
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
                Err(error) => {
                    return Err(error).with_context(|| {
                        format!("create EM spill directory below {}", base.display())
                    });
                }
            }
        }
        let dir = dir.context("could not allocate a unique EM spill directory")?;
        let writers = match (0..EM_SHARDS)
            .map(|shard| {
                File::create(dir.join(format!("shard-{shard:02}.bin")))
                    .map(|file| Some(BufWriter::with_capacity(1 << 20, file)))
            })
            .collect::<std::io::Result<Vec<_>>>()
        {
            Ok(writers) => writers,
            Err(error) => {
                let _ = std::fs::remove_dir_all(&dir);
                return Err(error)
                    .with_context(|| format!("create EM support shards in {}", dir.display()));
            }
        };
        Ok(Self {
            dir,
            writers,
            counts: vec![0; EM_SHARDS],
        })
    }

    fn append(&mut self, shard: usize, support: &[EmSupport]) -> Result<()> {
        let writer = self.writers[shard]
            .as_mut()
            .context("EM support spill writer is already closed")?;
        let mut bytes = Vec::with_capacity(support.len() * 8);
        for row in support {
            bytes.extend_from_slice(&row.class.to_le_bytes());
            bytes.extend_from_slice(&row.gene_flags.to_le_bytes());
        }
        writer
            .write_all(&bytes)
            .with_context(|| format!("write EM support shard {shard}"))?;
        self.counts[shard] += support.len();
        Ok(())
    }

    fn close_writers(&mut self) -> Result<()> {
        for (shard, writer) in self.writers.iter_mut().enumerate() {
            if let Some(mut writer) = writer.take() {
                writer
                    .flush()
                    .with_context(|| format!("flush EM support shard {shard}"))?;
            }
        }
        Ok(())
    }

    fn read_shard(&self, shard: usize) -> Result<Vec<EmSupport>> {
        let path = self.dir.join(format!("shard-{shard:02}.bin"));
        let bytes = std::fs::read(&path)
            .with_context(|| format!("read EM support shard {}", path.display()))?;
        if bytes.len() != self.counts[shard] * 8 {
            bail!(
                "EM support shard {shard} has {} bytes for {} records",
                bytes.len(),
                self.counts[shard]
            );
        }
        Ok(bytes
            .chunks_exact(8)
            .map(|record| EmSupport {
                class: u32::from_le_bytes(record[..4].try_into().unwrap()),
                gene_flags: u32::from_le_bytes(record[4..].try_into().unwrap()),
            })
            .collect())
    }
}

impl Drop for EmSupportSpill {
    fn drop(&mut self) {
        self.writers.clear();
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

impl PackedEmAccumulator {
    pub(crate) fn new(
        x: &Extracted,
        anno: &anno::Annotation,
        expected_molecules: usize,
    ) -> Result<Self> {
        if anno.gene_ids.len() > EM_GENE_MASK as usize + 1 {
            bail!("annotation has too many genes for packed EM labels");
        }
        let bam2anno = x
            .chrom_names
            .iter()
            .map(|n| anno.chrom_ids.get(n).copied())
            .collect();
        Ok(Self {
            n_classes: x.n_classes,
            gene_count: anno.gene_ids.len(),
            bam2anno,
            shards: (0..EM_SHARDS).map(|_| Vec::new()).collect(),
            spill: (expected_molecules.saturating_mul(std::mem::size_of::<EmSupport>())
                > EM_IN_MEMORY_SUPPORT_BUDGET)
                .then(EmSupportSpill::new)
                .transpose()?,
        })
    }

    pub(crate) fn spills_supports(&self) -> bool {
        self.spill.is_some()
    }

    fn classify(
        &self,
        mols: &[MolRec],
        x: &Extracted,
        anno: &anno::Annotation,
        s: &mut RowScratch,
    ) -> Vec<Vec<EmSupport>> {
        let mut shards: Vec<Vec<EmSupport>> = (0..EM_SHARDS).map(|_| Vec::new()).collect();
        for m in mols {
            let shard = m.cell as usize % EM_SHARDS;
            let mut handle = |r: Row| {
                if row_genes(&r, x, anno, &self.bam2anno, MmMissing::DropRow, s).is_none()
                    || s.genes.is_empty()
                {
                    return;
                }
                let flag = if r.weight == 0 {
                    0
                } else if s.genes.len() == 1 {
                    EM_UNIQUE_FLAG
                } else {
                    EM_MULTI_FLAG
                };
                for &gene in &s.genes {
                    shards[shard].push(EmSupport {
                        class: r.umi_class,
                        gene_flags: gene | flag,
                    });
                }
            };
            for ch in &m.chains {
                for (pos, shape) in &ch.reps {
                    handle(Row {
                        cell: m.cell,
                        umi_class: m.umi_class,
                        chrom: m.chrom,
                        pos: *pos,
                        strand_rev: m.strand_rev,
                        shape: *shape,
                        weight: ch.weight,
                        pattern: u32::MAX,
                    });
                }
            }
            for (pos, shape, pattern, weight) in &m.mms {
                handle(Row {
                    cell: m.cell,
                    umi_class: m.umi_class,
                    chrom: m.chrom,
                    pos: *pos,
                    strand_rev: m.strand_rev,
                    shape: *shape,
                    weight: *weight,
                    pattern: *pattern,
                });
            }
        }
        shards
    }

    pub(crate) fn add_archive_chunks(
        &mut self,
        chunks: &[Vec<MolRec>],
        x: &Extracted,
        anno: &anno::Annotation,
    ) -> Result<()> {
        let tasks: Vec<&[MolRec]> = chunks
            .iter()
            .flat_map(|chunk| chunk.chunks(1 << 16))
            .collect();
        let mut parts: Vec<Vec<Vec<EmSupport>>> = tasks
            .into_par_iter()
            .map_init(RowScratch::default, |s, part| {
                self.classify(part, x, anno, s)
            })
            .collect();
        parts.par_iter_mut().for_each(|worker_shards| {
            for support in worker_shards {
                compact_em_supports(support);
            }
        });
        if let Some(spill) = &mut self.spill {
            for worker_shards in &parts {
                for (shard, support) in worker_shards.iter().enumerate() {
                    spill.append(shard, support)?;
                }
            }
        } else {
            for k in 0..EM_SHARDS {
                self.shards[k].reserve_exact(parts.iter().map(|p| p[k].len()).sum());
            }
            for worker_shards in &mut parts {
                for (dst, src) in self.shards.iter_mut().zip(worker_shards) {
                    dst.append(src);
                }
            }
        }
        Ok(())
    }

    pub(crate) fn finish(
        mut self,
        cell_of_class: &[u32],
        mask_frac: f64,
        seed: u64,
    ) -> Result<PackedEm> {
        if cell_of_class.len() != self.n_classes as usize {
            bail!(
                "cell-of-class table has {} entries for {} classes",
                cell_of_class.len(),
                self.n_classes
            );
        }
        let mut shards = Vec::with_capacity(EM_SHARDS);
        let mut counters = EmCounters::default();
        if let Some(mut spill) = self.spill.take() {
            spill.close_writers()?;
            for shard_id in 0..EM_SHARDS {
                let mut support = spill.read_shard(shard_id)?;
                support.sort_unstable_by_key(|r| (r.class, r.gene_flags & EM_GENE_MASK));
                let (shard, c) =
                    build_em_shard(shard_id, support, cell_of_class, mask_frac, seed)?;
                shards.push(shard);
                counters += c;
            }
        } else {
            let built: Vec<Result<(PackedEmShard, EmCounters)>> = self
                .shards
                .into_par_iter()
                .enumerate()
                .map(|(shard_id, mut support)| {
                    support.sort_unstable_by_key(|r| (r.class, r.gene_flags & EM_GENE_MASK));
                    build_em_shard(shard_id, support, cell_of_class, mask_frac, seed)
                })
                .collect();
            for item in built {
                let (shard, c) = item?;
                shards.push(shard);
                counters += c;
            }
        }
        let mut global_base = vec![0.0f64; self.gene_count];
        for shard in &shards {
            for (&key, &count) in shard.cell_gene_keys.iter().zip(&shard.unique_base) {
                if count != 0.0 {
                    global_base[(key & u32::MAX as u64) as usize] += count;
                }
            }
        }
        Ok(PackedEm {
            shards,
            global_base,
            counters,
            mask_frac,
        })
    }
}

/// Within one bounded decode task, classification may see the same class/gene through several
/// representatives. Finalization consumes only the OR of their evidence-kind flags, so collapse
/// those duplicates before they enter the resident shard buffers. The final sort/reduction still
/// combines identical keys that span decode tasks or batches.
fn compact_em_supports(support: &mut Vec<EmSupport>) {
    support.sort_unstable_by_key(|row| (row.class, row.gene_flags & EM_GENE_MASK));
    let mut write = 0usize;
    for read in 0..support.len() {
        let row = support[read];
        if write > 0
            && support[write - 1].class == row.class
            && support[write - 1].gene_flags & EM_GENE_MASK == row.gene_flags & EM_GENE_MASK
        {
            support[write - 1].gene_flags |= row.gene_flags & !EM_GENE_MASK;
        } else {
            support[write] = row;
            write += 1;
        }
    }
    support.truncate(write);
}

#[inline]
fn cell_gene_key(cell: u32, gene: u32) -> u64 {
    ((cell as u64) << 32) | gene as u64
}

/// SplitMix64 keyed holdout: selection depends only on the declared seed and biological identity,
/// never on hash-table capacity, allocation order, worker count, or shard scheduling.
#[inline]
fn em_masked(seed: u64, cell: u32, class: u32, frac: f64) -> bool {
    if frac <= 0.0 {
        return false;
    }
    if frac >= 1.0 {
        return true;
    }
    let mut z = seed ^ cell_gene_key(cell, class);
    z = z.wrapping_add(0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^= z >> 31;
    let u = (z >> 11) as f64 * (1.0 / ((1u64 << 53) as f64));
    u < frac
}

fn push_packed_target(
    cell: u32,
    genes: impl Iterator<Item = u32>,
    truth: u32,
    starts: &mut Vec<u64>,
    target_cells: &mut Vec<u32>,
    target_genes: &mut Vec<u32>,
    target_truth: &mut Vec<u32>,
    key_events: &mut Vec<u64>,
) -> Result<()> {
    target_cells.push(cell);
    target_truth.push(truth);
    for gene in genes {
        target_genes.push(gene);
        key_events.push(cell_gene_key(cell, gene));
    }
    starts.push(u64::try_from(target_genes.len()).context("packed EM incidence count overflow")?);
    Ok(())
}

fn build_em_shard(
    shard_id: usize,
    support: Vec<EmSupport>,
    cell_of_class: &[u32],
    mask_frac: f64,
    seed: u64,
) -> Result<(PackedEmShard, EmCounters)> {
    let mut starts = vec![0u64];
    let mut target_cells = Vec::new();
    let mut target_genes = Vec::new();
    let mut target_truth = Vec::new();
    let mut unique_events: Vec<u64> = Vec::new();
    let mut key_events: Vec<u64> = Vec::new();
    let mut counters = EmCounters::default();
    let mut class_genes: SmallVec<[(u32, u32); 4]> = SmallVec::new();

    let mut i = 0usize;
    while i < support.len() {
        let class = support[i].class;
        let cell = *cell_of_class
            .get(class as usize)
            .with_context(|| format!("EM support references class {class} beyond cell table"))?;
        if cell as usize % EM_SHARDS != shard_id {
            bail!("class {class} is in EM shard {shard_id} but belongs to cell {cell}");
        }
        class_genes.clear();
        while i < support.len() && support[i].class == class {
            let gene = support[i].gene_flags & EM_GENE_MASK;
            let mut flags = 0u32;
            while i < support.len()
                && support[i].class == class
                && support[i].gene_flags & EM_GENE_MASK == gene
            {
                flags |= support[i].gene_flags & !EM_GENE_MASK;
                i += 1;
            }
            class_genes.push((gene, flags));
        }

        let unique_genes: SmallVec<[u32; 2]> = class_genes
            .iter()
            .filter(|(_, flags)| flags & EM_UNIQUE_FLAG != 0)
            .map(|(gene, _)| *gene)
            .collect();
        if unique_genes.len() == 1 && class_genes.len() == 1 {
            counters.single += 1;
            let key = cell_gene_key(cell, class_genes[0].0);
            unique_events.push(key);
            key_events.push(key);
        } else if unique_genes.len() == 1 {
            counters.mixed += 1;
            let truth = unique_genes[0];
            if em_masked(seed, cell, class, mask_frac) {
                counters.masked += 1;
                let cands: SmallVec<[u32; 4]> = class_genes
                    .iter()
                    .filter(|(_, flags)| flags & EM_MULTI_FLAG != 0)
                    .map(|(gene, _)| *gene)
                    .collect();
                if cands.len() >= 2 && cands.contains(&truth) {
                    push_packed_target(
                        cell,
                        cands.into_iter(),
                        truth,
                        &mut starts,
                        &mut target_cells,
                        &mut target_genes,
                        &mut target_truth,
                        &mut key_events,
                    )?;
                } else {
                    counters.truth_lost += 1;
                }
            } else {
                let key = cell_gene_key(cell, truth);
                unique_events.push(key);
                key_events.push(key);
            }
        } else if unique_genes.is_empty() && class_genes.len() >= 2 {
            counters.multi_only += 1;
            push_packed_target(
                cell,
                class_genes.iter().map(|(gene, _)| *gene),
                EM_NO_TRUTH,
                &mut starts,
                &mut target_cells,
                &mut target_genes,
                &mut target_truth,
                &mut key_events,
            )?;
        }
        // unique_genes >= 2 is contradictory and remains excluded, matching the v1 experiment.
    }
    drop(support);

    unique_events.sort_unstable();
    key_events.sort_unstable();
    key_events.dedup();
    let cell_gene_keys = key_events;
    let mut unique_base = vec![0.0f64; cell_gene_keys.len()];
    let mut j = 0usize;
    while j < unique_events.len() {
        let key = unique_events[j];
        let mut e = j + 1;
        while e < unique_events.len() && unique_events[e] == key {
            e += 1;
        }
        let at = cell_gene_keys
            .binary_search(&key)
            .expect("unique key inserted into universe");
        unique_base[at] = (e - j) as f64;
        j = e;
    }

    let mut target_local = Vec::with_capacity(target_genes.len());
    for (&cell, range) in target_cells.iter().zip(starts.windows(2)) {
        for &gene in &target_genes[range[0] as usize..range[1] as usize] {
            let at = cell_gene_keys
                .binary_search(&cell_gene_key(cell, gene))
                .expect("target key inserted into universe");
            target_local.push(u32::try_from(at).context("EM shard has more than u32::MAX keys")?);
        }
    }
    Ok((
        PackedEmShard {
            cell_gene_keys,
            unique_base,
            starts,
            target_cells,
            target_genes,
            target_local,
            target_truth,
        },
        counters,
    ))
}

#[derive(Default)]
struct PackedModeStats {
    top1: u64,
    expected: f64,
    negative_log_loss: f64,
    brier: f64,
    n: u64,
    cal: [[u64; 2]; 10],
}

impl PackedModeStats {
    fn add(&mut self, other: &Self) {
        self.top1 += other.top1;
        self.expected += other.expected;
        self.negative_log_loss += other.negative_log_loss;
        self.brier += other.brier;
        self.n += other.n;
        for (a, b) in self.cal.iter_mut().zip(&other.cal) {
            a[0] += b[0];
            a[1] += b[1];
        }
    }

    fn observe(&mut self, correct: bool, truth_r: f64, sum_r_squared: f64, max_r: f64) {
        self.n += 1;
        self.top1 += correct as u64;
        self.expected += truth_r;
        self.negative_log_loss -= truth_r.max(f64::MIN_POSITIVE).ln();
        self.brier += sum_r_squared + 1.0 - 2.0 * truth_r;
        let b = ((max_r * 10.0) as usize).min(9);
        self.cal[b][0] += 1;
        self.cal[b][1] += correct as u64;
    }
}

const EM_EVIDENCE_DEPTH_NAMES: [&str; 4] = ["0-1", "1-4", "4-16", "16+"];

#[inline]
fn evidence_depth_stratum(total: f64) -> usize {
    if total <= 1.0 {
        0
    } else if total <= 4.0 {
        1
    } else if total <= 16.0 {
        2
    } else {
        3
    }
}

#[inline]
fn packed_weight(mode: u8, alpha: f64, total_unique_global: f64, local: f64, global: f64) -> f64 {
    match mode {
        1 => local + 1e-9,
        2 => global + 1e-9,
        _ => local + alpha * (global / total_unique_global) + 1e-9,
    }
}

#[inline]
fn packed_hierarchical_weight(
    mode: u8,
    group_alpha: f64,
    global_alpha: f64,
    total_group: f64,
    total_global: f64,
    local: f64,
    group: f64,
    global: f64,
) -> f64 {
    match mode {
        // Raw group and global counts differ only by a candidate-independent scale when all cells
        // occupy one group, so this makes the collapse control numerically identical to pooled.
        4 => group + 1e-9,
        _ => {
            local
                + group_alpha * (group / total_group.max(f64::MIN_POSITIVE))
                + global_alpha * (global / total_global.max(f64::MIN_POSITIVE))
                + 1e-9
        }
    }
}

/// One candidate probability under convex partial pooling. Component distributions are
/// normalized over the target's candidate set before mixing. The group component removes the
/// target cell's current abundance, then borrows `group_prior` candidate-level pseudo-counts from
/// the sample distribution. An empty cell or leave-one-cell-out group falls back to the sample
/// distribution. For an ungrouped cell, the group weight is also transferred to the sample.
#[allow(clippy::too_many_arguments)]
#[inline]
fn packed_convex_probability(
    params: ConvexEmParams,
    has_group: bool,
    local: f64,
    group: f64,
    global: f64,
    local_candidate_total: f64,
    loo_group_candidate_total: f64,
    global_candidate_total: f64,
) -> f64 {
    let p_global = (global + 1e-9) / global_candidate_total;
    let p_cell = if local_candidate_total > 0.0 {
        local / local_candidate_total
    } else {
        p_global
    };
    let p_group = if has_group {
        let denom = loo_group_candidate_total + params.group_prior;
        if denom > 0.0 {
            ((group - local).max(0.0) + params.group_prior * p_global) / denom
        } else {
            p_global
        }
    } else {
        p_global
    };
    let group_weight = if has_group { params.group_weight } else { 0.0 };
    let global_weight = params.global_weight() + if has_group { 0.0 } else { params.group_weight };
    params.cell_weight * p_cell + group_weight * p_group + global_weight * p_global
}

/// One posterior-mean probability under the hierarchical Dirichlet screening approximation.
/// Both posterior updates operate on the target candidate set, and the target cell is removed
/// from the group counts. Empty cell or leave-one-cell-out group posteriors fall back toward the
/// next level. This is a cheap test of evidence-adaptive pooling, not full variational inference.
#[allow(clippy::too_many_arguments)]
#[inline]
fn packed_dirichlet_proxy_probability(
    params: DirichletProxyParams,
    has_group: bool,
    local: f64,
    group: f64,
    global: f64,
    local_candidate_total: f64,
    loo_group_candidate_total: f64,
    global_candidate_total: f64,
) -> f64 {
    let p_global = (global + 1e-9) / global_candidate_total;
    let p_group = if has_group {
        let denom = loo_group_candidate_total + params.group_prior;
        if denom > 0.0 {
            ((group - local).max(0.0) + params.group_prior * p_global) / denom
        } else {
            p_global
        }
    } else {
        p_global
    };
    let denom = local_candidate_total + params.cell_prior;
    if denom > 0.0 {
        (local + params.cell_prior * p_group) / denom
    } else {
        p_group
    }
}

#[inline]
fn monotone_depth_gate(base_candidate_total: f64, depth_scale: f64, depth_power: f64) -> f64 {
    if base_candidate_total <= 0.0 {
        return 0.0;
    }
    // Stable logistic form of d^p / (d^p + D^p), avoiding pow overflow at large depth/power.
    let logit = depth_power * (base_candidate_total.ln() - depth_scale.ln());
    if logit >= 0.0 {
        1.0 / (1.0 + (-logit).exp())
    } else {
        let exp = logit.exp();
        exp / (1.0 + exp)
    }
}

/// Blend two already candidate-normalized models with a weight fixed by pre-mode evidence depth.
/// Because the weight is candidate-independent within a target, the result remains normalized.
#[allow(clippy::too_many_arguments)]
#[inline]
fn packed_depth_hybrid_probability(
    hybrid: DepthHybridParams,
    base_candidate_total: f64,
    has_group: bool,
    local: f64,
    group: f64,
    global: f64,
    local_candidate_total: f64,
    loo_group_candidate_total: f64,
    global_candidate_total: f64,
) -> f64 {
    let convex = ConvexEmParams {
        cell_weight: 0.2,
        group_weight: 0.6,
        group_prior: 80.0,
    };
    let dirichlet = DirichletProxyParams {
        cell_prior: 64.0,
        group_prior: 80.0,
    };
    let fixed = packed_convex_probability(
        convex,
        has_group,
        local,
        group,
        global,
        local_candidate_total,
        loo_group_candidate_total,
        global_candidate_total,
    );
    let adaptive = packed_dirichlet_proxy_probability(
        dirichlet,
        has_group,
        local,
        group,
        global,
        local_candidate_total,
        loo_group_candidate_total,
        global_candidate_total,
    );
    let gate = monotone_depth_gate(
        base_candidate_total,
        hybrid.depth_scale,
        hybrid.depth_power,
    );
    fixed + gate * (adaptive - fixed)
}

fn clustered_standard_error(differences: &[f64], cells: &[u32], mean: f64) -> (u64, f64) {
    assert_eq!(differences.len(), cells.len());
    if differences.is_empty() {
        return (0, 0.0);
    }
    let mut cluster_sums = vec![0.0f64; cells.iter().copied().max().unwrap_or(0) as usize + 1];
    let mut used = vec![false; cluster_sums.len()];
    for (&difference, &cell) in differences.iter().zip(cells) {
        cluster_sums[cell as usize] += difference - mean;
        used[cell as usize] = true;
    }
    let clusters = used.iter().filter(|value| **value).count() as u64;
    if clusters <= 1 {
        return (clusters, 0.0);
    }
    let sum_squares: f64 = cluster_sums.iter().map(|value| value * value).sum();
    let correction = clusters as f64 / (clusters - 1) as f64;
    (
        clusters,
        (correction * sum_squares).sqrt() / differences.len() as f64,
    )
}

fn paired_metrics(
    candidate_name: &str,
    candidate: &[PackedTargetObservation],
    reference_name: &str,
    reference: &[PackedTargetObservation],
) -> PackedPairedMetrics {
    assert_eq!(candidate.len(), reference.len());
    let mut cells = Vec::with_capacity(candidate.len());
    let mut nll = Vec::with_capacity(candidate.len());
    let mut brier = Vec::with_capacity(candidate.len());
    let mut nll_wins = 0u64;
    let mut brier_wins = 0u64;
    let mut top1_difference = 0i64;
    for (candidate, reference) in candidate.iter().zip(reference) {
        assert_eq!(candidate.cell, reference.cell);
        cells.push(candidate.cell);
        let nll_difference = candidate.negative_log_loss - reference.negative_log_loss;
        let brier_difference = candidate.brier - reference.brier;
        nll.push(nll_difference);
        brier.push(brier_difference);
        nll_wins += (nll_difference < 0.0) as u64;
        brier_wins += (brier_difference < 0.0) as u64;
        top1_difference += candidate.correct as i64 - reference.correct as i64;
    }
    let denominator = candidate.len().max(1) as f64;
    let mean_nll = nll.iter().sum::<f64>() / denominator;
    let mean_brier = brier.iter().sum::<f64>() / denominator;
    let (clusters, nll_se) = clustered_standard_error(&nll, &cells, mean_nll);
    let (_, brier_se) = clustered_standard_error(&brier, &cells, mean_brier);
    PackedPairedMetrics {
        candidate: candidate_name.to_string(),
        reference: reference_name.to_string(),
        n: candidate.len() as u64,
        cells: clusters,
        mean_negative_log_loss_difference: mean_nll,
        negative_log_loss_clustered_se: nll_se,
        negative_log_loss_ci95: [mean_nll - 1.96 * nll_se, mean_nll + 1.96 * nll_se],
        negative_log_loss_wins: nll_wins,
        mean_brier_difference: mean_brier,
        brier_clustered_se: brier_se,
        brier_ci95: [mean_brier - 1.96 * brier_se, mean_brier + 1.96 * brier_se],
        brier_wins,
        top1_percentage_point_difference: 100.0 * top1_difference as f64 / denominator,
    }
}

impl PackedEm {
    /// Run the four fixed sharing modes, plus optional group, additive-hierarchical,
    /// candidate-normalized convex, posterior-mean Dirichlet-proxy, and monotone depth-hybrid
    /// modes, over packed CSR labels. Each iteration walks shards in a fixed order and classes in
    /// sorted class order, so
    /// f64 accumulation is reproducible across worker counts. Responsibilities are recomputed from
    /// the previous iteration's fixed state and never stored as one allocation per target.
    pub(crate) fn run(
        &self,
        anno: &anno::Annotation,
        alpha: f64,
        groups: Option<&EmGroups>,
        group_alpha: f64,
        global_alpha: f64,
        convex: ConvexEmParams,
        convex_only: bool,
        dirichlet: DirichletProxyParams,
        dirichlet_only: bool,
        hybrid: DepthHybridParams,
        hybrid_only: bool,
        cal_out: Option<&mut Vec<(String, [[u64; 2]; 10], f64, u64)>>,
        metrics_out: Option<&mut Vec<PackedModeMetrics>>,
        paired_metrics_out: Option<&mut Vec<PackedPairedMetrics>>,
    ) -> Option<FxHashMap<(u32, u32), f64>> {
        if let Some(groups) = groups {
            assert_eq!(groups.cell_group.len(), groups.eval_cell.len());
        }
        let c = self.counters;
        let eval_targets: usize = self
            .shards
            .iter()
            .map(|s| {
                s.target_truth
                    .iter()
                    .zip(&s.target_cells)
                    .filter(|(truth, cell)| {
                        **truth != EM_NO_TRUTH
                            && groups.is_none_or(|g| g.eval_cell[**cell as usize])
                    })
                    .count()
            })
            .sum();
        println!(
            "classes: single {}, mixed {} (masked {}, truth-lost {}), multi-only {}, eval targets {}",
            c.single, c.mixed, c.masked, c.truth_lost, c.multi_only, eval_targets
        );

        let total_unique_global: f64 = self.global_base.iter().sum();
        assert!(
            [convex_only, dirichlet_only, hybrid_only]
                .into_iter()
                .filter(|value| *value)
                .count()
                <= 1
        );
        let mut modes = if convex_only {
            vec![("convex", 6u8)]
        } else if dirichlet_only {
            vec![("dirichlet-proxy", 7u8)]
        } else if hybrid_only {
            vec![("depth-hybrid", 8u8)]
        } else {
            vec![("uniform", 0u8), ("cell", 1), ("pooled", 2), ("blend", 3)]
        };
        if groups.is_some() && !convex_only && !dirichlet_only && !hybrid_only {
            modes.extend([
                ("group", 4),
                ("hierarchical", 5),
                ("convex", 6),
                ("dirichlet-proxy", 7),
                ("depth-hybrid", 8),
            ]);
        }
        let mut cal_rows = Vec::new();
        let mut metric_rows = Vec::new();
        let collect_paired = paired_metrics_out.is_some()
            && !convex_only
            && !dirichlet_only
            && !hybrid_only
            && groups.is_some();
        let mut target_observations: Vec<(String, Vec<PackedTargetObservation>)> = Vec::new();
        let mut impact: Option<(Vec<Vec<f64>>, Vec<f64>, u64)> = None;

        let group_base = groups.map(|groups| {
            let mut base = vec![0.0f64; groups.names.len() * self.global_base.len()];
            for shard in &self.shards {
                for (&key, &count) in shard.cell_gene_keys.iter().zip(&shard.unique_base) {
                    let group = groups.cell_group[(key >> 32) as usize];
                    if group != EM_NO_GROUP && count != 0.0 {
                        let gene = key as u32 as usize;
                        base[group as usize * self.global_base.len() + gene] += count;
                    }
                }
            }
            base
        });

        for (name, mode) in modes {
            let mut stats = PackedModeStats::default();
            let mut depth_stats: [PackedModeStats; 4] = std::array::from_fn(|_| Default::default());
            let mut observations = (collect_paired && matches!(mode, 2 | 6 | 7 | 8))
                .then(|| Vec::with_capacity(eval_targets));
            if mode == 0 {
                for shard in &self.shards {
                    for i in 0..shard.target_cells.len() {
                        let truth = shard.target_truth[i];
                        if truth == EM_NO_TRUTH
                            || groups.is_some_and(|g| !g.eval_cell[shard.target_cells[i] as usize])
                        {
                            continue;
                        }
                        let range = shard.target_range(i);
                        let n = range.len();
                        let ti = shard.target_genes[range.clone()]
                            .iter()
                            .position(|g| *g == truth)
                            .unwrap();
                        // Iterator::max_by keeps the last equal item; uniform therefore chooses
                        // the final packed candidate, matching Rust's v1 tie rule.
                        let argmax = n - 1;
                        let r = 1.0 / n as f64;
                        let correct = argmax == ti;
                        let base_candidate_total: f64 = range
                            .clone()
                            .map(|k| shard.unique_base[shard.target_local[k] as usize])
                            .sum();
                        let depth = evidence_depth_stratum(base_candidate_total);
                        stats.observe(correct, r, 1.0 / n as f64, r);
                        depth_stats[depth].observe(correct, r, 1.0 / n as f64, r);
                    }
                }
            } else {
                let mut pi_cell: Vec<Vec<f64>> =
                    self.shards.iter().map(|s| s.unique_base.clone()).collect();
                let mut pi_global = self.global_base.clone();
                let mut pi_group = (mode >= 4).then(|| {
                    group_base
                        .as_ref()
                        .expect("group mode requires groups")
                        .clone()
                });
                let capture_impact = mode == 2 && self.mask_frac == 0.0;
                let mut recovered = capture_impact.then(|| {
                    self.shards
                        .iter()
                        .map(|s| vec![0.0; s.cell_gene_keys.len()])
                        .collect::<Vec<_>>()
                });
                let mut gene_gain = capture_impact.then(|| vec![0.0f64; self.global_base.len()]);
                let mut confident = 0u64;

                for iter in 0..10 {
                    let final_iter = iter == 9;
                    let mut next_global = self.global_base.clone();
                    let mut next_group = pi_group.as_ref().map(|_| {
                        group_base.as_ref().expect("group mode requires groups").clone()
                    });
                    let iter_total_global: f64 = pi_global.iter().sum();
                    let group_totals = pi_group.as_ref().map(|values| {
                        values
                            .chunks_exact(self.global_base.len())
                            .map(|slice| slice.iter().sum::<f64>())
                            .collect::<Vec<_>>()
                    });
                    for (shard_id, shard) in self.shards.iter().enumerate() {
                        let old_cell = &pi_cell[shard_id];
                        let mut next_cell = shard.unique_base.clone();
                        let mut shard_stats = PackedModeStats::default();
                        let mut shard_depth_stats: [PackedModeStats; 4] =
                            std::array::from_fn(|_| Default::default());
                        for ti in 0..shard.target_cells.len() {
                            let range = shard.target_range(ti);
                            let cell = shard.target_cells[ti] as usize;
                            let group = groups.map_or(EM_NO_GROUP, |g| g.cell_group[cell]);
                            let (denom, candidate_totals) = if mode >= 6 {
                                let mut local_total = 0.0f64;
                                let mut loo_group_total = 0.0f64;
                                let mut global_total = 0.0f64;
                                let mut base_total = 0.0f64;
                                for k in range.clone() {
                                    let gene = shard.target_genes[k] as usize;
                                    let local_idx = shard.target_local[k] as usize;
                                    let local = old_cell[local_idx];
                                    local_total += local;
                                    base_total += shard.unique_base[local_idx];
                                    global_total += pi_global[gene] + 1e-9;
                                    if group != EM_NO_GROUP {
                                        let at = group as usize * self.global_base.len() + gene;
                                        loo_group_total +=
                                            (pi_group.as_ref().unwrap()[at] - local).max(0.0);
                                    }
                                }
                                (
                                    1.0,
                                    Some((local_total, loo_group_total, global_total, base_total)),
                                )
                            } else {
                                let mut denom = 0.0f64;
                                for k in range.clone() {
                                    let gene = shard.target_genes[k] as usize;
                                    let local = old_cell[shard.target_local[k] as usize];
                                    denom += if mode < 4 || group == EM_NO_GROUP {
                                        packed_weight(
                                            if mode == 4 { 2 } else if mode == 5 { 3 } else { mode },
                                            if mode == 5 { global_alpha } else { alpha },
                                            total_unique_global,
                                            local,
                                            pi_global[gene],
                                        )
                                    } else {
                                        let at = group as usize * self.global_base.len() + gene;
                                        packed_hierarchical_weight(
                                            mode,
                                            group_alpha,
                                            global_alpha,
                                            group_totals.as_ref().unwrap()[group as usize],
                                            iter_total_global,
                                            local,
                                            pi_group.as_ref().unwrap()[at],
                                            pi_global[gene],
                                        )
                                    };
                                }
                                (denom, None)
                            };
                            let truth = shard.target_truth[ti];
                            let mut argmax = 0usize;
                            let mut max_r = f64::NEG_INFINITY;
                            let mut truth_r = 0.0f64;
                            let mut sum_r_squared = 0.0f64;
                            for (within, k) in range.clone().enumerate() {
                                let gene = shard.target_genes[k] as usize;
                                let local_idx = shard.target_local[k] as usize;
                                let w = if mode == 6 {
                                    let (local_total, loo_group_total, global_total, _) =
                                        candidate_totals.unwrap();
                                    let group_value = if group == EM_NO_GROUP {
                                        0.0
                                    } else {
                                        let at =
                                            group as usize * self.global_base.len() + gene;
                                        pi_group.as_ref().unwrap()[at]
                                    };
                                    packed_convex_probability(
                                        convex,
                                        group != EM_NO_GROUP,
                                        old_cell[local_idx],
                                        group_value,
                                        pi_global[gene],
                                        local_total,
                                        loo_group_total,
                                        global_total,
                                    )
                                } else if mode == 7 {
                                    let (local_total, loo_group_total, global_total, _) =
                                        candidate_totals.unwrap();
                                    let group_value = if group == EM_NO_GROUP {
                                        0.0
                                    } else {
                                        let at =
                                            group as usize * self.global_base.len() + gene;
                                        pi_group.as_ref().unwrap()[at]
                                    };
                                    packed_dirichlet_proxy_probability(
                                        dirichlet,
                                        group != EM_NO_GROUP,
                                        old_cell[local_idx],
                                        group_value,
                                        pi_global[gene],
                                        local_total,
                                        loo_group_total,
                                        global_total,
                                    )
                                } else if mode == 8 {
                                    let (local_total, loo_group_total, global_total, base_total) =
                                        candidate_totals.unwrap();
                                    let group_value = if group == EM_NO_GROUP {
                                        0.0
                                    } else {
                                        let at =
                                            group as usize * self.global_base.len() + gene;
                                        pi_group.as_ref().unwrap()[at]
                                    };
                                    packed_depth_hybrid_probability(
                                        hybrid,
                                        base_total,
                                        group != EM_NO_GROUP,
                                        old_cell[local_idx],
                                        group_value,
                                        pi_global[gene],
                                        local_total,
                                        loo_group_total,
                                        global_total,
                                    )
                                } else if mode < 4 || group == EM_NO_GROUP {
                                    packed_weight(
                                        if mode == 4 { 2 } else if mode == 5 { 3 } else { mode },
                                        if mode == 5 { global_alpha } else { alpha },
                                        total_unique_global,
                                        old_cell[local_idx],
                                        pi_global[gene],
                                    )
                                } else {
                                    let at = group as usize * self.global_base.len() + gene;
                                    packed_hierarchical_weight(
                                        mode,
                                        group_alpha,
                                        global_alpha,
                                        group_totals.as_ref().unwrap()[group as usize],
                                        iter_total_global,
                                        old_cell[local_idx],
                                        pi_group.as_ref().unwrap()[at],
                                        pi_global[gene],
                                    )
                                };
                                let r = w / denom;
                                next_cell[local_idx] += r;
                                next_global[gene] += r;
                                if let Some(values) = &mut next_group {
                                    if group != EM_NO_GROUP {
                                        values[group as usize * self.global_base.len() + gene] += r;
                                    }
                                }
                                sum_r_squared += r * r;
                                if r >= max_r {
                                    max_r = r;
                                    argmax = within;
                                }
                                if shard.target_genes[k] == truth {
                                    truth_r = r;
                                }
                                if final_iter && capture_impact {
                                    recovered.as_mut().unwrap()[shard_id][local_idx] += r;
                                    gene_gain.as_mut().unwrap()[gene] += r;
                                }
                            }
                            if final_iter
                                && truth != EM_NO_TRUTH
                                && groups.is_none_or(|g| g.eval_cell[cell])
                            {
                                let truth_i = shard.target_genes[range.clone()]
                                    .iter()
                                    .position(|g| *g == truth)
                                    .unwrap();
                                let correct = argmax == truth_i;
                                let base_candidate_total: f64 = range
                                    .clone()
                                    .map(|k| shard.unique_base[shard.target_local[k] as usize])
                                    .sum();
                                let depth = evidence_depth_stratum(base_candidate_total);
                                shard_stats.observe(correct, truth_r, sum_r_squared, max_r);
                                shard_depth_stats[depth]
                                    .observe(correct, truth_r, sum_r_squared, max_r);
                                if let Some(observations) = &mut observations {
                                    observations.push(PackedTargetObservation {
                                        cell: cell as u32,
                                        negative_log_loss: -truth_r.max(f64::MIN_POSITIVE).ln(),
                                        brier: sum_r_squared + 1.0 - 2.0 * truth_r,
                                        correct,
                                    });
                                }
                            }
                            if final_iter && capture_impact {
                                confident += (max_r > 0.8) as u64;
                            }
                        }
                        pi_cell[shard_id] = next_cell;
                        if final_iter {
                            stats.add(&shard_stats);
                            for (total, shard_total) in
                                depth_stats.iter_mut().zip(&shard_depth_stats)
                            {
                                total.add(shard_total);
                            }
                        }
                    }
                    pi_global = next_global;
                    pi_group = next_group;
                }
                if capture_impact {
                    impact = Some((recovered.unwrap(), gene_gain.unwrap(), confident));
                }
            }

            println!(
                "mode {name:8}: top-1 {:.1}%  expected-accuracy {:.1}%  (n={})",
                100.0 * stats.top1 as f64 / stats.n.max(1) as f64,
                100.0 * stats.expected / stats.n.max(1) as f64,
                stats.n
            );
            metric_rows.push(PackedModeMetrics {
                name: name.to_string(),
                top1: stats.top1,
                expected: stats.expected,
                negative_log_loss: stats.negative_log_loss,
                brier: stats.brier,
                n: stats.n,
                calibration: stats.cal,
                evidence_depth: EM_EVIDENCE_DEPTH_NAMES
                    .into_iter()
                    .zip(depth_stats)
                    .map(|(name, depth)| PackedDepthMetrics {
                        name,
                        top1: depth.top1,
                        expected: depth.expected,
                        negative_log_loss: depth.negative_log_loss,
                        brier: depth.brier,
                        n: depth.n,
                    })
                    .collect(),
            });
            if let Some(observations) = observations {
                target_observations.push((name.to_string(), observations));
            }
            if stats.n > 0 {
                cal_rows.push((
                    name.to_string(),
                    stats.cal,
                    stats.top1 as f64 / stats.n as f64,
                    stats.n,
                ));
            }
            if mode == 2 {
                print!("  calibration (max-r decile: empirical accuracy | n): ");
                for (i, [cn, cc]) in stats.cal.iter().enumerate() {
                    if *cn > 0 {
                        print!(
                            "{:.1}-: {:.0}%|{} ",
                            i as f64 / 10.0,
                            100.0 * *cc as f64 / *cn as f64,
                            cn
                        );
                    }
                }
                println!();
            }
        }
        if let Some(out) = cal_out {
            *out = cal_rows;
        }
        if let Some(out) = metrics_out {
            *out = metric_rows;
        }
        if let Some(out) = paired_metrics_out {
            if collect_paired {
                let find = |name: &str| {
                    target_observations
                        .iter()
                        .find(|(candidate, _)| candidate == name)
                        .map(|(_, observations)| observations.as_slice())
                        .expect("paired EM mode was not collected")
                };
                *out = [
                    ("depth-hybrid", "pooled"),
                    ("depth-hybrid", "convex"),
                    ("depth-hybrid", "dirichlet-proxy"),
                    ("dirichlet-proxy", "pooled"),
                    ("convex", "pooled"),
                ]
                .into_iter()
                .map(|(candidate, reference)| {
                    paired_metrics(candidate, find(candidate), reference, find(reference))
                })
                .collect();
            }
        }

        if let Some((recovered, gene_gain, confident)) = impact {
            let counted_today = c.single + c.mixed;
            let n_targets: usize = self.shards.iter().map(|s| s.target_cells.len()).sum();
            println!(
                "EM-2: recovered {} multi-only classes (+{:.2}% of {} counted); {:.1}% assigned with r > 0.8",
                n_targets,
                100.0 * n_targets as f64 / counted_today.max(1) as f64,
                counted_today,
                100.0 * confident as f64 / n_targets.max(1) as f64
            );
            let mut top: Vec<(u32, f64)> = gene_gain
                .iter()
                .enumerate()
                .filter(|(_, mass)| **mass > 0.0)
                .map(|(gene, mass)| (gene as u32, *mass))
                .collect();
            top.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            for (gene, mass) in top.iter().take(15) {
                println!("  +{:8.0}  {}", mass, anno.gene_names[*gene as usize]);
            }
            let mut out = FxHashMap::default();
            for (shard, values) in self.shards.iter().zip(recovered) {
                for (&key, value) in shard.cell_gene_keys.iter().zip(values) {
                    if value > 1e-6 {
                        out.insert(((key >> 32) as u32, key as u32), value);
                    }
                }
            }
            return Some(out);
        }
        None
    }

    pub(crate) fn candidate_gene_ids(&self, anno: &anno::Annotation) -> Vec<String> {
        let mut used = vec![false; self.global_base.len()];
        for shard in &self.shards {
            for i in 0..shard.target_cells.len() {
                if shard.target_truth[i] == EM_NO_TRUTH {
                    continue;
                }
                for &gene in &shard.target_genes[shard.target_range(i)] {
                    used[gene as usize] = true;
                }
            }
        }
        let mut genes: Vec<String> = used.into_iter()
            .zip(&anno.gene_ids)
            .filter_map(|(used, gene)| used.then(|| gene.clone()))
            .collect();
        genes.sort_unstable();
        genes
    }
}

#[cfg(test)]
mod hierarchical_em_tests {
    use super::{
        compact_em_supports, monotone_depth_gate, packed_convex_probability,
        packed_depth_hybrid_probability, packed_dirichlet_proxy_probability,
        packed_hierarchical_weight, packed_weight, paired_metrics, write_multigene_audit,
        ConvexEmParams, DepthHybridParams, DirichletProxyParams, EmSupport, EmSupportSpill,
        MultigeneAuditSummary, PackedTargetObservation, EM_MULTI_FLAG, EM_UNIQUE_FLAG,
    };

    #[test]
    fn multigene_audit_legacy_rendering_is_frozen() {
        let summary = MultigeneAuditSummary {
            no_gene_read_weight: 2,
            single_gene_read_weight: 6,
            multi_gene_read_weight: 2,
            no_gene_rows: 1,
            single_gene_rows: 3,
            multi_gene_rows: 1,
            gene_evidence_classes: 9,
            single_gene_counted_classes: 4,
            multi_only_classes: 2,
            mixed_classes: 1,
            no_gene_classes: 2,
        };
        let mut bytes = Vec::new();
        write_multigene_audit(&mut bytes, &summary).unwrap();
        assert_eq!(
            String::from_utf8(bytes).unwrap(),
            "== multigene audit (EM upside bound) ==\n\
rows by concordant-gene count (read-weighted): 0 genes 20.00% | 1 gene 60.00% | >=2 genes 20.00%  (1 / 3 / 1 rows)\n\
classes with gene evidence: single-gene-counted 4 | MULTI-ONLY (fully discarded, the EM pool) 2 (40.00% of counted) | mixed 1; no-gene classes 2\n"
        );
    }

    #[test]
    fn support_spill_roundtrips_and_cleans_up() {
        let rows = [
            EmSupport {
                class: 1,
                gene_flags: 7 | EM_UNIQUE_FLAG,
            },
            EmSupport {
                class: 2,
                gene_flags: 9 | EM_MULTI_FLAG,
            },
        ];
        let mut spill = EmSupportSpill::new().unwrap();
        let dir = spill.dir.clone();
        spill.append(3, &rows).unwrap();
        spill.close_writers().unwrap();
        let restored = spill.read_shard(3).unwrap();
        assert_eq!(restored.len(), rows.len());
        for (actual, expected) in restored.iter().zip(rows) {
            assert_eq!((actual.class, actual.gene_flags), (expected.class, expected.gene_flags));
        }
        drop(spill);
        assert!(!dir.exists());
    }

    #[test]
    fn batch_support_compaction_ors_flags_without_merging_genes_or_classes() {
        let mut rows = vec![
            EmSupport {
                class: 2,
                gene_flags: 7 | EM_MULTI_FLAG,
            },
            EmSupport {
                class: 1,
                gene_flags: 7 | EM_UNIQUE_FLAG,
            },
            EmSupport {
                class: 1,
                gene_flags: 9 | EM_MULTI_FLAG,
            },
            EmSupport {
                class: 1,
                gene_flags: 7 | EM_MULTI_FLAG,
            },
        ];
        compact_em_supports(&mut rows);
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].class, 1);
        assert_eq!(rows[0].gene_flags, 7 | EM_UNIQUE_FLAG | EM_MULTI_FLAG);
        assert_eq!(
            (rows[1].class, rows[1].gene_flags),
            (1, 9 | EM_MULTI_FLAG)
        );
        assert_eq!(
            (rows[2].class, rows[2].gene_flags),
            (2, 7 | EM_MULTI_FLAG)
        );
    }

    #[test]
    fn one_group_weight_is_exactly_pooled_weight() {
        for count in [0.0, 1.0, 17.5, 1_000_000.0] {
            assert_eq!(
                packed_hierarchical_weight(4, 20.0, 5.0, 10.0, 10.0, 3.0, count, count),
                packed_weight(2, 20.0, 10.0, 3.0, count)
            );
        }
    }

    #[test]
    fn convex_probabilities_are_candidate_normalized() {
        let params = ConvexEmParams {
            cell_weight: 0.1,
            group_weight: 0.45,
            group_prior: 20.0,
        };
        let local = [1.0, 3.0, 0.0];
        let group = [6.0, 5.0, 1.0];
        let global = [20.0, 10.0, 5.0];
        let local_total = local.iter().sum();
        let loo_group_total: f64 = group
            .iter()
            .zip(local)
            .map(|(group, local)| f64::max(*group - local, 0.0))
            .sum();
        let global_total: f64 = global.iter().map(|value| value + 1e-9).sum();
        let sum: f64 = local
            .into_iter()
            .zip(group)
            .zip(global)
            .map(|((local, group), global)| {
                packed_convex_probability(
                    params,
                    true,
                    local,
                    group,
                    global,
                    local_total,
                    loo_group_total,
                    global_total,
                )
            })
            .sum();
        assert!((sum - 1.0).abs() < 1e-12, "candidate probabilities sum to {sum}");
    }

    #[test]
    fn convex_global_vertex_matches_pooled_probability() {
        let params = ConvexEmParams {
            cell_weight: 0.0,
            group_weight: 0.0,
            group_prior: 20.0,
        };
        let global = [7.0, 2.0, 0.0];
        let denom: f64 = global.iter().map(|value| value + 1e-9).sum();
        for value in global {
            let convex = packed_convex_probability(
                params, true, 99.0, 123.0, value, 99.0, 24.0, denom,
            );
            assert_eq!(convex, packed_weight(2, 0.0, 1.0, 99.0, value) / denom);
        }
    }

    #[test]
    fn convex_group_component_leaves_target_cell_out() {
        let params = ConvexEmParams {
            cell_weight: 0.0,
            group_weight: 1.0,
            group_prior: 0.0,
        };
        // With no other candidate evidence in the group, group-only convex sharing falls back to
        // the sample distribution rather than recycling the target cell's 9:1 local evidence.
        let global = [2.0, 8.0];
        let global_total: f64 = global.iter().map(|value| value + 1e-9).sum();
        for (local, value) in [9.0, 1.0].into_iter().zip(global) {
            let probability = packed_convex_probability(
                params,
                true,
                local,
                local,
                value,
                10.0,
                0.0,
                global_total,
            );
            assert_eq!(probability, (value + 1e-9) / global_total);
        }
    }

    #[test]
    fn dirichlet_proxy_probabilities_are_candidate_normalized() {
        let params = DirichletProxyParams {
            cell_prior: 16.0,
            group_prior: 20.0,
        };
        let local = [1.0, 3.0, 0.0];
        let group = [6.0, 5.0, 1.0];
        let global = [20.0, 10.0, 5.0];
        let local_total = local.iter().sum();
        let loo_group_total: f64 = group
            .iter()
            .zip(local)
            .map(|(group, local)| f64::max(*group - local, 0.0))
            .sum();
        let global_total: f64 = global.iter().map(|value| value + 1e-9).sum();
        let sum: f64 = local
            .into_iter()
            .zip(group)
            .zip(global)
            .map(|((local, group), global)| {
                packed_dirichlet_proxy_probability(
                    params,
                    true,
                    local,
                    group,
                    global,
                    local_total,
                    loo_group_total,
                    global_total,
                )
            })
            .sum();
        assert!((sum - 1.0).abs() < 1e-12, "candidate probabilities sum to {sum}");
    }

    #[test]
    fn dirichlet_proxy_adapts_pooling_to_cell_evidence_depth() {
        let params = DirichletProxyParams {
            cell_prior: 10.0,
            group_prior: 20.0,
        };
        let global_total = 2.0 + 2e-9;
        let shallow = packed_dirichlet_proxy_probability(
            params,
            false,
            9.0,
            0.0,
            1.0,
            10.0,
            0.0,
            global_total,
        );
        let deep = packed_dirichlet_proxy_probability(
            params,
            false,
            90.0,
            0.0,
            1.0,
            100.0,
            0.0,
            global_total,
        );
        assert!(deep > shallow, "deeper evidence should receive more cell weight");
        assert!((shallow - 0.7).abs() < 1e-9);
        assert!((deep - 95.0 / 110.0).abs() < 1e-9);
    }

    #[test]
    fn depth_hybrid_gate_is_monotone_and_has_locked_midpoint() {
        let scale = 16.0;
        let values = [0.0, 1.0, 4.0, 16.0, 64.0, 1024.0]
            .map(|depth| monotone_depth_gate(depth, scale, 2.0));
        assert_eq!(values[0], 0.0);
        assert_eq!(values[3], 0.5);
        assert!(values.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(values[5] < 1.0 && values[5] > 0.999);
    }

    #[test]
    fn depth_hybrid_probabilities_are_candidate_normalized() {
        let hybrid = DepthHybridParams {
            depth_scale: 16.0,
            depth_power: 2.0,
        };
        let local = [1.0, 3.0, 0.0];
        let group = [6.0, 5.0, 1.0];
        let global = [20.0, 10.0, 5.0];
        let local_total = local.iter().sum();
        let loo_group_total: f64 = group
            .iter()
            .zip(local)
            .map(|(group, local)| f64::max(*group - local, 0.0))
            .sum();
        let global_total: f64 = global.iter().map(|value| value + 1e-9).sum();
        let sum: f64 = local
            .into_iter()
            .zip(group)
            .zip(global)
            .map(|((local, group), global)| {
                packed_depth_hybrid_probability(
                    hybrid,
                    16.0,
                    true,
                    local,
                    group,
                    global,
                    local_total,
                    loo_group_total,
                    global_total,
                )
            })
            .sum();
        assert!((sum - 1.0).abs() < 1e-12, "candidate probabilities sum to {sum}");
    }

    #[test]
    fn depth_hybrid_moves_from_fixed_toward_proxy_with_depth() {
        let convex = ConvexEmParams {
            cell_weight: 0.2,
            group_weight: 0.6,
            group_prior: 80.0,
        };
        let dirichlet = DirichletProxyParams {
            cell_prior: 64.0,
            group_prior: 80.0,
        };
        let hybrid = DepthHybridParams {
            depth_scale: 16.0,
            depth_power: 2.0,
        };
        let args = (true, 8.0, 30.0, 2.0, 10.0, 20.0, 10.0 + 2e-9);
        let fixed = packed_convex_probability(
            convex, args.0, args.1, args.2, args.3, args.4, args.5, args.6,
        );
        let proxy = packed_dirichlet_proxy_probability(
            dirichlet, args.0, args.1, args.2, args.3, args.4, args.5, args.6,
        );
        let shallow = packed_depth_hybrid_probability(
            hybrid, 1.0, args.0, args.1, args.2, args.3, args.4, args.5, args.6,
        );
        let deep = packed_depth_hybrid_probability(
            hybrid, 256.0, args.0, args.1, args.2, args.3, args.4, args.5, args.6,
        );
        assert!((shallow - fixed).abs() < (deep - fixed).abs());
        assert!((deep - proxy).abs() < (shallow - proxy).abs());
    }

    #[test]
    fn paired_metrics_preserve_target_pairing_and_cluster_by_cell() {
        let reference = [
            PackedTargetObservation {
                cell: 0,
                negative_log_loss: 0.3,
                brier: 0.2,
                correct: false,
            },
            PackedTargetObservation {
                cell: 0,
                negative_log_loss: 0.2,
                brier: 0.1,
                correct: true,
            },
            PackedTargetObservation {
                cell: 1,
                negative_log_loss: 0.4,
                brier: 0.3,
                correct: false,
            },
        ];
        let candidate = [
            PackedTargetObservation {
                cell: 0,
                negative_log_loss: 0.2,
                brier: 0.1,
                correct: true,
            },
            PackedTargetObservation {
                cell: 0,
                negative_log_loss: 0.3,
                brier: 0.2,
                correct: true,
            },
            PackedTargetObservation {
                cell: 1,
                negative_log_loss: 0.1,
                brier: 0.2,
                correct: true,
            },
        ];
        let metrics = paired_metrics("candidate", &candidate, "reference", &reference);
        assert_eq!(metrics.n, 3);
        assert_eq!(metrics.cells, 2);
        assert!((metrics.mean_negative_log_loss_difference + 0.1).abs() < 1e-12);
        assert!((metrics.mean_brier_difference + 1.0 / 30.0).abs() < 1e-12);
        assert_eq!(metrics.negative_log_loss_wins, 2);
        assert_eq!(metrics.brier_wins, 2);
        assert!((metrics.top1_percentage_point_difference - 200.0 / 3.0).abs() < 1e-12);
        assert!(metrics.negative_log_loss_clustered_se > 0.0);
        assert!(metrics.brier_clustered_se > 0.0);
    }
}

/// Evaluate EM multimapper recovery with cross-cell sharing. Builds per-`(cell, class)` gene-support
/// vectors, masks a seeded 20% of mixed classes down to their multi-gene evidence (truth is the
/// gene supported by their unique rows), and scores per-cell, pooled, and blended EM against a
/// uniform baseline.
pub fn em_experiment(
    x: &Extracted,
    anno: &anno::Annotation,
    mask_frac: f64,
    seed: u64,
    alpha: f64,
    cal_out: Option<&mut Vec<(String, [[u64; 2]; 10], f64, u64)>>,
) -> Option<FxHashMap<(u32, u32), f64>> {
    let bam2anno: Vec<Option<u32>> =
        x.chrom_names.iter().map(|n| anno.chrom_ids.get(n).copied()).collect();
    let rows = flatten(&x.mols);
    // Per row: candidate genes + weight, tagged unique(1 gene) / multi(>1). DropRow preserves the
    // historical semantics: an mm alternative on an unannotated chromosome discards the row.
    // SmallVec: candidate lists are 1–2 genes for the overwhelming majority of rows — inline.
    let per_row: Vec<Option<(u32, u32, SmallVec<[u32; 4]>, u32)>> = rows
        .par_iter()
        .map_init(RowScratch::default, |s, r| {
            row_genes(r, x, anno, &bam2anno, MmMissing::DropRow, s)?;
            if s.genes.is_empty() {
                None
            } else {
                Some((r.cell, r.umi_class, SmallVec::from_slice(&s.genes), r.weight))
            }
        })
        .collect();

    // Class-level support: gene -> (unique weight, multi weight).
    struct Cls {
        support: Vec<(u32, u32, u32)>, // gene, u_w, m_w
    }
    let mut classes: FxHashMap<(u32, u32), Cls> = FxHashMap::default();
    for e in per_row.into_iter().flatten() {
        let (cell, cls, genes, w) = e;
        let uniq = genes.len() == 1;
        let c = classes.entry((cell, cls)).or_insert(Cls { support: Vec::new() });
        for g in genes {
            match c.support.iter_mut().find(|(gg, _, _)| *gg == g) {
                Some(s) => {
                    if uniq {
                        s.1 += w;
                    } else {
                        s.2 += w;
                    }
                }
                None => c.support.push((g, if uniq { w } else { 0 }, if !uniq { w } else { 0 })),
            }
        }
    }

    // Partition: single-gene (π estimators), mixed (labeled truth pool), multi-only (targets).
    // Truth for a mixed class = the gene its unique rows support (unique across genes in
    // practice; classes with unique support for >1 gene are excluded from the eval).
    let mut lcg = seed | 1;
    let mut rng = move || {
        lcg = lcg.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (lcg >> 33) as f64 / (1u64 << 31) as f64
    };
    let mut unique_counts_cell: FxHashMap<(u32, u32), f64> = FxHashMap::default(); // (cell,gene)
    let mut unique_counts_global: FxHashMap<u32, f64> = FxHashMap::default();
    struct Target {
        cell: u32,
        cands: Vec<u32>,
        truth: Option<u32>, // Some = masked eval item
    }
    let mut targets: Vec<Target> = Vec::new();
    let (mut n_single, mut n_mixed, mut n_multi_only, mut n_masked, mut n_truth_lost) =
        (0u64, 0u64, 0u64, 0u64, 0u64);
    for ((cell, _cls), c) in &classes {
        let with_u: Vec<&(u32, u32, u32)> = c.support.iter().filter(|(_, u, _)| *u > 0).collect();
        let all_genes: Vec<u32> = c.support.iter().map(|(g, _, _)| *g).collect();
        if with_u.len() == 1 && all_genes.len() == 1 {
            n_single += 1;
            *unique_counts_cell.entry((*cell, all_genes[0])).or_insert(0.0) += 1.0;
            *unique_counts_global.entry(all_genes[0]).or_insert(0.0) += 1.0;
        } else if with_u.len() == 1 {
            n_mixed += 1;
            let truth = with_u[0].0;
            if mask_frac > 0.0 && rng() < mask_frac {
                // Masked: candidate set = genes with MULTI evidence only (the unique row is
                // hidden). Truth must survive in the candidate set to be recoverable.
                let cands: Vec<u32> =
                    c.support.iter().filter(|(_, _, m)| *m > 0).map(|(g, _, _)| *g).collect();
                n_masked += 1;
                if cands.len() >= 2 && cands.contains(&truth) {
                    targets.push(Target { cell: *cell, cands, truth: Some(truth) });
                } else {
                    n_truth_lost += 1;
                }
            } else {
                // Unmasked mixed classes count as today (best gene = truth) and feed π.
                *unique_counts_cell.entry((*cell, truth)).or_insert(0.0) += 1.0;
                *unique_counts_global.entry(truth).or_insert(0.0) += 1.0;
            }
        } else if with_u.is_empty() && all_genes.len() >= 2 {
            n_multi_only += 1;
            targets.push(Target { cell: *cell, cands: all_genes, truth: None });
        }
        // with_u >= 2: unique evidence for two genes — contradictory, excluded (rare).
    }
    println!(
        "classes: single {n_single}, mixed {n_mixed} (masked {n_masked}, truth-lost {n_truth_lost}), multi-only {n_multi_only}, eval targets {}",
        targets.iter().filter(|t| t.truth.is_some()).count()
    );

    // EM per mode. Responsibilities over candidates from π scope; 10 iterations updating the
    // scope's π with assigned mass.
    //
    // Blend-mode normalizer, hoisted: the old code recomputed
    // `unique_counts_global.values().sum()` inside the per-target per-gene closure —
    // O(targets × candidates × genes) map traversals per iteration, which dominated EM wall time
    // on an approximately 100-million-molecule benchmark archive.
    // The map is never mutated after construction and FxHashMap iteration is deterministic, so
    // this single sum is bit-identical to every one of the sums it replaces.
    let total_unique_global: f64 = unique_counts_global.values().sum();
    let modes: Vec<(&str, u8)> = vec![("uniform", 0), ("cell", 1), ("pooled", 2), ("blend", 3)];
    let mut cal_rows: Vec<(String, [[u64; 2]; 10], f64, u64)> = Vec::new();
    for (name, mode) in modes {
        let mut pi_cell = unique_counts_cell.clone();
        let mut pi_global = unique_counts_global.clone();
        let mut resp: Vec<Vec<f64>> = targets
            .iter()
            .map(|t| vec![1.0 / t.cands.len() as f64; t.cands.len()])
            .collect();
        let iters = if mode == 0 { 1 } else { 10 };
        for _ in 0..iters {
            if mode != 0 {
                // E-step: each target's responsibilities depend only on the fixed π maps, so
                // targets parallelize with bit-identical results (per-target sums keep their
                // serial order; nothing accumulates across targets here). Written in place —
                // the collect-into-fresh-Vec version was ~35% of EM wall in malloc/free.
                targets.par_iter().zip(resp.par_iter_mut()).for_each(|(t, r)| {
                    for (v, g) in r.iter_mut().zip(t.cands.iter()) {
                        *v = match mode {
                            1 => pi_cell.get(&(t.cell, *g)).copied().unwrap_or(0.0) + 1e-9,
                            2 => pi_global.get(g).copied().unwrap_or(0.0) + 1e-9,
                            _ => {
                                pi_cell.get(&(t.cell, *g)).copied().unwrap_or(0.0)
                                    + alpha
                                        * (pi_global.get(g).copied().unwrap_or(0.0)
                                            / total_unique_global)
                                    + 1e-9
                            }
                        };
                    }
                    let s: f64 = r.iter().sum();
                    for v in r.iter_mut() {
                        *v /= s;
                    }
                });
                // M-step: rebuild π = unique base + assigned mass (clone_from reuses the tables).
                pi_cell.clone_from(&unique_counts_cell);
                pi_global.clone_from(&unique_counts_global);
                for (t, r) in targets.iter().zip(resp.iter()) {
                    for (g, v) in t.cands.iter().zip(r.iter()) {
                        *pi_cell.entry((t.cell, *g)).or_insert(0.0) += v;
                        *pi_global.entry(*g).or_insert(0.0) += v;
                    }
                }
            }
        }
        // Score the masked targets, with a calibration curve: does max-responsibility r mean
        // "correct with probability r"?
        let (mut top1, mut expacc, mut n) = (0u64, 0.0f64, 0u64);
        let mut cal = [[0u64; 2]; 10]; // decile of max-r -> [n, n_correct]
        for (t, r) in targets.iter().zip(resp.iter()) {
            let Some(truth) = t.truth else { continue };
            n += 1;
            let ti = t.cands.iter().position(|g| *g == truth).unwrap();
            expacc += r[ti];
            let argmax = r
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .map(|(i, _)| i)
                .unwrap();
            let correct = argmax == ti;
            top1 += correct as u64;
            let mx = r[argmax];
            let b2 = ((mx * 10.0) as usize).min(9);
            cal[b2][0] += 1;
            cal[b2][1] += correct as u64;
        }
        println!(
            "mode {name:8}: top-1 {:.1}%  expected-accuracy {:.1}%  (n={n})",
            100.0 * top1 as f64 / n.max(1) as f64,
            100.0 * expacc / n.max(1) as f64
        );
        if n > 0 {
            cal_rows.push((name.to_string(), cal, top1 as f64 / n as f64, n));
        }
        if mode == 2 {
            print!("  calibration (max-r decile: empirical accuracy | n): ");
            for (i, [cn, cc]) in cal.iter().enumerate() {
                if *cn > 0 {
                    print!("{:.1}-: {:.0}%|{} ", i as f64 / 10.0, 100.0 * *cc as f64 / *cn as f64, cn);
                }
            }
            println!();
        }
    }
    if let Some(out) = cal_out {
        *out = cal_rows;
    }
    // EM-2 (mask = 0): matrix impact with the winning (pooled) mode over the real multi-only
    // pool — recovered mass, its confidence profile, and where it lands.
    if mask_frac == 0.0 {
        let mut pi_global = unique_counts_global.clone();
        let mut resp: Vec<Vec<f64>> = targets
            .iter()
            .map(|t| vec![1.0 / t.cands.len() as f64; t.cands.len()])
            .collect();
        for _ in 0..10 {
            // Parallel in-place E-step, same bit-identical argument as above.
            targets.par_iter().zip(resp.par_iter_mut()).for_each(|(t, r)| {
                for (v, g) in r.iter_mut().zip(t.cands.iter()) {
                    *v = pi_global.get(g).copied().unwrap_or(0.0) + 1e-9;
                }
                let s: f64 = r.iter().sum();
                for v in r.iter_mut() {
                    *v /= s;
                }
            });
            pi_global.clone_from(&unique_counts_global);
            for (t, r) in targets.iter().zip(resp.iter()) {
                for (g, v) in t.cands.iter().zip(r.iter()) {
                    *pi_global.entry(*g).or_insert(0.0) += v;
                }
            }
        }
        let counted_today = n_single + n_mixed;
        let mut confident = 0u64;
        let mut gene_gain: FxHashMap<u32, f64> = FxHashMap::default();
        for (t, r) in targets.iter().zip(resp.iter()) {
            let mx = r.iter().cloned().fold(0.0f64, f64::max);
            confident += (mx > 0.8) as u64;
            for (g, v) in t.cands.iter().zip(r.iter()) {
                *gene_gain.entry(*g).or_insert(0.0) += v;
            }
        }
        let mut top: Vec<(u32, f64)> = gene_gain.into_iter().collect();
        top.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        println!(
            "EM-2: recovered {} multi-only classes (+{:.2}% of {} counted); {:.1}% assigned with r > 0.8",
            targets.len(),
            100.0 * targets.len() as f64 / counted_today as f64,
            counted_today,
            100.0 * confident as f64 / targets.len().max(1) as f64
        );
        for (g, m) in top.iter().take(15) {
            println!("  +{:8.0}  {}", m, anno.gene_names[*g as usize]);
        }
        // Per-(cell, gene) fractional recovered mass, for the additive matrix layer.
        let mut recovered: FxHashMap<(u32, u32), f64> = FxHashMap::default();
        for (t, r) in targets.iter().zip(resp.iter()) {
            for (g, v) in t.cands.iter().zip(r.iter()) {
                if *v > 1e-6 {
                    *recovered.entry((t.cell, *g)).or_insert(0.0) += v;
                }
            }
        }
        return Some(recovered);
    }
    None
}

/// The (cell, class) → canonical-class merges that the Gene collapse performs: replay_rows'
/// assignment → aggregation → best-gene filter → tie-merge, replicated to return the merge map
/// instead of counts (velocity needs the corrected-UMI identities, not the numbers). Only merged
/// (non-root) classes appear in the map.
/// Strand-parameterized Gene correction map, kept in lockstep with the Gene feature that supplies
/// STAR velocity's corrected UMI identities.
pub fn gene_canon_map_stranded(
    x: &Extracted,
    anno: &anno::Annotation,
    solo_strand: anno::assign::SoloStrand,
) -> FxHashMap<(u32, u32), u32> {
    let bam2anno: Vec<Option<u32>> =
        x.chrom_names.iter().map(|n| anno.chrom_ids.get(n).copied()).collect();
    let rows = flatten(&x.mols);
    let assigned: Vec<Option<(u32, u32, u32, u32)>> = rows
        .par_iter()
        .map_init(RowScratch::default, |s, r| {
            row_genes_stranded(
                r, x, anno, &bam2anno, MmMissing::SkipAlt, solo_strand, s,
            )?;
            match s.genes.as_slice() {
                [g] => Some((r.cell, r.umi_class, *g, r.weight)),
                _ => None,
            }
        })
        .collect();

    const NSH: usize = 64;
    let mut shards: Vec<Vec<(u32, u32, u32, u32)>> = vec![Vec::new(); NSH];
    for a in assigned.into_iter().flatten() {
        shards[a.0 as usize % NSH].push(a);
    }
    let mut adj: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for (a, b) in &x.edges {
        adj.entry(*a).or_default().push(*b);
        adj.entry(*b).or_default().push(*a);
    }

    shards
        .into_par_iter()
        .flat_map_iter(|shard| {
            let mut per_cg: FxHashMap<(u32, u32), FxHashMap<u32, u32>> = FxHashMap::default();
            let mut class_genes: FxHashMap<(u32, u32), FxHashMap<u32, u32>> = FxHashMap::default();
            for (cell, cls, gene, w) in shard {
                *per_cg.entry((cell, gene)).or_default().entry(cls).or_insert(0) += w;
                *class_genes.entry((cell, cls)).or_default().entry(gene).or_insert(0) += w;
            }
            let mut best: FxHashMap<(u32, u32), u32> = FxHashMap::default();
            for ((cell, cls), genes) in &class_genes {
                let b = genes.iter().max_by_key(|(g, c)| (**c, std::cmp::Reverse(**g))).map(|(g, _)| *g).unwrap();
                best.insert((*cell, *cls), b);
            }
            let mut out: Vec<((u32, u32), u32)> = Vec::new();
            for ((cell, gene), cnts) in per_cg {
                let kept: Vec<(u32, u32)> = cnts
                    .into_iter()
                    .filter(|(cls, _)| best.get(&(cell, *cls)) == Some(&gene))
                    .collect();
                if kept.is_empty() {
                    continue;
                }
                let mut order = kept;
                order.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
                let rank: FxHashMap<u32, usize> =
                    order.iter().enumerate().map(|(i, (c, _))| (*c, i)).collect();
                let mut canon: FxHashMap<u32, u32> = FxHashMap::default();
                for (idx, (c, _)) in order.iter().enumerate() {
                    let mut tgt: Option<usize> = None;
                    if let Some(nbs) = adj.get(c) {
                        for nb in nbs {
                            if let Some(&r2) = rank.get(nb) {
                                if r2 < idx {
                                    tgt = Some(tgt.map_or(r2, |t| t.min(r2)));
                                }
                            }
                        }
                    }
                    match tgt {
                        Some(r2) => {
                            let cc = *canon.get(&order[r2].0).unwrap();
                            canon.insert(*c, cc);
                            out.push(((cell, *c), cc));
                        }
                        None => {
                            canon.insert(*c, *c);
                        }
                    }
                }
            }
            out
        })
        .collect()
}

/// Intersect two transcript-sorted (tr, bits) lists, OR-ing bits — `countVelocyto` lines 66–77.
fn intersect_or(a: &[(u32, u8)], b: &[(u32, u8)]) -> Vec<(u32, u8)> {
    let mut out = Vec::with_capacity(a.len().min(b.len()));
    let mut j = 0usize;
    for (t, ba) in a {
        while j < b.len() && b[j].0 < *t {
            j += 1;
        }
        if j == b.len() {
            break;
        }
        if b[j].0 == *t {
            out.push((*t, ba | b[j].1));
        }
    }
    out
}

/// Reconstruct a placement from (chrom, pos, strand, shape). Junctions are the gaps between the
/// shape's blocks — they were never stored separately, which is why shapes suffice.
pub fn placement_from_parts(chrom: u32, pos: u32, strand_rev: bool, shape: &Shape, nh: u16) -> Placement {
    let blocks: Vec<Block> = shape
        .blocks
        .iter()
        .map(|&(off, len)| Block { start: pos + off, end: pos + off + len })
        .collect();
    let junctions: Vec<Junction> = blocks
        .windows(2)
        .map(|w| Junction { donor: w[0].end, acceptor: w[1].start })
        .collect();
    Placement {
        chrom,
        strand: if strand_rev { Strand::Reverse } else { Strand::Forward },
        blocks,
        junctions,
        nm: 0,
        score: 0,
        nh,
        clip: (0, 0),
    }
}

/// In-place variant of [`placement_from_parts`]: reuses the caller's block/junction buffers, so
/// the per-row classifiers pay zero allocations per placement.
pub fn placement_from_parts_into(
    p: &mut Placement,
    chrom: u32,
    pos: u32,
    strand_rev: bool,
    shape: &Shape,
    nh: u16,
) {
    p.chrom = chrom;
    p.strand = if strand_rev { Strand::Reverse } else { Strand::Forward };
    p.nh = nh;
    p.blocks.clear();
    p.junctions.clear();
    for &(off, len) in &shape.blocks {
        p.blocks.push(Block { start: pos + off, end: pos + off + len });
    }
    for w in shape.blocks.windows(2) {
        p.junctions.push(Junction { donor: pos + w[0].0 + w[0].1, acceptor: pos + w[1].0 });
    }
}

/// Per-worker scratch for the row classifiers (one per rayon split via `map_init`): placement
/// geometry, the transcript-overlap buffer and the gene accumulators are reused across rows
/// instead of allocated per row — this churn was the assignment stage's dominant cost.
pub struct RowScratch {
    place: Placement,
    txbuf: Vec<u32>,
    genes: Vec<u32>,
    alt_genes: Vec<u32>,
    /// (chrom, pos, shape, pattern, strand_rev) of the previous row: those five fields fully
    /// determine the gene set, and rows arrive position-clustered (UMIs pile up at hot loci), so
    /// a one-entry memo skips the whole classification for exact repeats.
    last_key: Option<(u32, u32, u32, u32, bool)>,
    last_none: bool,
}

impl Default for RowScratch {
    fn default() -> RowScratch {
        RowScratch {
            place: Placement {
                chrom: 0,
                strand: Strand::Forward,
                blocks: Vec::new(),
                junctions: Vec::new(),
                nm: 0,
                score: 0,
                nh: 1,
                clip: (0, 0),
            },
            txbuf: Vec::new(),
            genes: Vec::new(),
            alt_genes: Vec::new(),
            last_key: None,
            last_none: false,
        }
    }
}

/// How a multimapper alternative on a chromosome absent from the annotation is treated: `SkipAlt`
/// drops just that alternative (replay/canon/audit/EM-0 semantics); `DropRow` discards the whole
/// row (`em_experiment`'s historical semantics, preserved exactly).
#[derive(Clone, Copy, PartialEq, Eq)]
enum MmMissing {
    SkipAlt,
    DropRow,
}

/// Candidate genes for one row, into `s.genes` — same gene order as the old allocating code.
/// Returns `None` exactly when the old code early-returned before inspecting genes: a unique row
/// on an unannotated chromosome, or (with `DropRow`) any mm alternative on one.
fn row_genes(
    r: &Row,
    x: &Extracted,
    anno: &anno::Annotation,
    bam2anno: &[Option<u32>],
    mm_missing: MmMissing,
    s: &mut RowScratch,
) -> Option<()> {
    row_genes_stranded(
        r, x, anno, bam2anno, mm_missing, anno::assign::SoloStrand::Forward, s,
    )
}

fn row_genes_stranded(
    r: &Row,
    x: &Extracted,
    anno: &anno::Annotation,
    bam2anno: &[Option<u32>],
    mm_missing: MmMissing,
    solo_strand: anno::assign::SoloStrand,
    s: &mut RowScratch,
) -> Option<()> {
    // One-entry memo: same placement key → same genes (s.genes still holds them). The mm-missing
    // mode is constant across one pass, so it cannot alias between cached entries.
    let key = (r.chrom, r.pos, r.shape, r.pattern, r.strand_rev);
    if s.last_key == Some(key) {
        return if s.last_none { None } else { Some(()) };
    }
    s.last_key = Some(key);
    s.last_none = true; // corrected on the success path below
    s.genes.clear();
    if r.pattern == u32::MAX {
        let ac = (*bam2anno.get(r.chrom as usize)?)?;
        placement_from_parts_into(&mut s.place, r.chrom, r.pos, r.strand_rev, &x.shapes[r.shape as usize], 1);
        anno::assign::concordant_genes_stranded_into(
            &s.place, anno, ac, solo_strand, &mut s.txbuf, &mut s.genes,
        );
    } else {
        for alt in &x.patterns[r.pattern as usize] {
            let ac = match bam2anno.get(alt.chrom as usize).copied().flatten() {
                Some(ac) => ac,
                None => match mm_missing {
                    MmMissing::SkipAlt => continue,
                    MmMissing::DropRow => return None,
                },
            };
            let apos = (r.pos as i64 + alt.offset) as u32;
            let arev = r.strand_rev != alt.strand_flip;
            let shape_id = if alt.shape == SAME_SHAPE { r.shape } else { alt.shape };
            placement_from_parts_into(&mut s.place, alt.chrom, apos, arev, &x.shapes[shape_id as usize], 2);
            anno::assign::concordant_genes_stranded_into(
                &s.place, anno, ac, solo_strand, &mut s.txbuf, &mut s.alt_genes,
            );
            for &g in &s.alt_genes {
                if !s.genes.contains(&g) {
                    s.genes.push(g);
                }
            }
        }
    }
    s.last_none = false;
    Some(())
}

fn chain_hash(pos: u32, shape: &Shape) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = rustc_hash::FxHasher::default();
    for w in shape.blocks.windows(2) {
        (pos + w[0].0 + w[0].1, pos + w[1].0).hash(&mut h);
    }
    h.finish()
}

fn end_of(pos: u32, shape: &Shape) -> u32 {
    shape.blocks.last().map(|&(o, l)| pos + o + l).unwrap_or(pos)
}

struct URead {
    cell: u32,
    umi: u32,
    chrom: u32,
    pos: u32,
    strand_rev: bool,
    shape: u32,
}

/// BAM reader over a multithreaded BGZF decoder: block inflation runs on a worker pool while the
/// caller parses records. Worker count is deliberately below the rayon cap — beyond ~8 the serial
/// record-parsing consumer is the bottleneck, not decompression.
fn open_bam_mt(
    bam_path: &PathBuf,
) -> Result<bam::io::Reader<noodles_bgzf::io::MultithreadedReader<std::fs::File>>> {
    let file =
        std::fs::File::open(bam_path).with_context(|| format!("opening {}", bam_path.display()))?;
    Ok(bam_reader_mt(file))
}

fn bam_reader_mt(
    file: File,
) -> bam::io::Reader<noodles_bgzf::io::MultithreadedReader<std::fs::File>> {
    let workers = std::num::NonZeroUsize::new(rayon::current_num_threads().clamp(1, 8)).unwrap();
    bam::io::Reader::from(
        noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, file),
    )
}

pub fn extract_rows(bam_path: &PathBuf, whitelist: &PathBuf, locus_gap: u32) -> Result<Extracted> {
    let wl = load_whitelist(whitelist)?;
    Ok(extract_rows_inner(bam_path, wl, locus_gap, false)?.0)
}

/// Reporting path: parse the exact whitelist snapshot supplied by the caller and derive the BAM
/// identity from the same held file used by both extraction passes.
pub(crate) fn extract_rows_with_identity(
    bam_path: &PathBuf,
    whitelist_text: &str,
    locus_gap: u32,
) -> Result<(Extracted, ConsumedFileIdentity)> {
    let mut wl = FxHashSet::default();
    for line in whitelist_text.lines() {
        let bases = line.trim().as_bytes();
        if bases.len() == 16 {
            if let Some(packed) = umi::pack(bases) {
                wl.insert(packed);
            }
        }
    }
    let (extracted, identity) = extract_rows_inner(bam_path, wl, locus_gap, true)?;
    Ok((
        extracted,
        identity.expect("reporting extraction requested a BAM identity"),
    ))
}

fn extract_rows_inner(
    bam_path: &PathBuf,
    wl: FxHashSet<u32>,
    locus_gap: u32,
    capture_identity: bool,
) -> Result<(Extracted, Option<ConsumedFileIdentity>)> {
    let source = if capture_identity {
        let file = File::open(bam_path)
            .with_context(|| format!("opening {}", bam_path.display()))?;
        let before = file.metadata()?;
        if !before.is_file() {
            bail!("raw alignment input is not a regular file: {}", bam_path.display());
        }
        Some((file, before))
    } else {
        None
    };

    // Pass 1: exact-whitelist frequencies for the pseudocount corrector (annotation-free).
    let mut freq: FxHashMap<u32, u32> = FxHashMap::default();
    let source = {
        let (mut r0, before) = match source {
            Some((file, before)) => (bam_reader_mt(file), Some(before)),
            None => (open_bam_mt(bam_path)?, None),
        };
        r0.read_header()?;
        let mut rec = bam::Record::default();
        while r0.read_record(&mut rec)? != 0 {
            let f = Flags::from(rec.flags().bits());
            if f.is_secondary() || f.is_supplementary() {
                continue;
            }
            for field in rec.data().iter() {
                let (tag, value) = field?;
                if <[u8; 2]>::from(tag) == *b"CR" {
                    if let Value::String(s) = value {
                        if let Some(pk) = umi::pack(s.as_ref()) {
                            if wl.contains(&pk) {
                                *freq.entry(pk).or_insert(0) += 1;
                            }
                        }
                    }
                    break;
                }
            }
        }
        match before {
            Some(before) => {
                let mut bgzf = r0.into_inner();
                let mut file = bgzf.finish()?;
                file.rewind().with_context(|| {
                    format!("rewinding {} for extraction pass 2", bam_path.display())
                })?;
                Some((file, before))
            }
            None => None,
        }
    };
    let corrector = BcCorrector::new(wl, freq);

    // Pass 2: all alignment records.
    let (mut reader, source_before) = match source {
        Some((file, before)) => (bam_reader_mt(file), Some(before)),
        None => (open_bam_mt(bam_path)?, None),
    };
    let header = reader.read_header()?;
    let chrom_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .collect();

    let mut cell_intern: FxHashMap<u32, u32> = FxHashMap::default();
    let mut cells: Vec<u32> = Vec::new();
    let mut shape_intern: FxHashMap<Shape, u32> = FxHashMap::default();
    let mut shapes: Vec<Shape> = Vec::new();
    let mut bc_cache: FxHashMap<u64, Option<u32>> = FxHashMap::default();

    let mut ureads: Vec<URead> = Vec::new();
    type MmEntry = (Option<(u32, u32)>, Option<usize>, Vec<(u32, u32, bool, u32)>); // (cb,umi), prim idx, alts (chrom,pos,rev,shape)
    let mut mm: FxHashMap<u64, MmEntry> = FxHashMap::default();

    let mut rec = bam::Record::default();
    while reader.read_record(&mut rec)? != 0 {
        let flags = Flags::from(rec.flags().bits());
        if flags.is_supplementary() || flags.is_unmapped() {
            continue;
        }
        let (mut cr_b, mut cy_b, mut ur, mut nh) = (None::<Vec<u8>>, None::<Vec<u8>>, None, 1u16);
        for field in rec.data().iter() {
            let (tag, value) = field?;
            let key = <[u8; 2]>::from(tag);
            match value {
                Value::String(s) => match &key {
                    b"CR" => cr_b = Some(<[u8] as ToOwned>::to_owned(s.as_ref())),
                    b"CY" => cy_b = Some(<[u8] as ToOwned>::to_owned(s.as_ref())),
                    b"UR" => ur = umi::pack(s.as_ref()),
                    _ => {}
                },
                v => {
                    if key == *b"NH" {
                        nh = match v {
                            Value::Int8(x) => x as u16,
                            Value::UInt8(x) => x as u16,
                            Value::Int16(x) => x as u16,
                            Value::UInt16(x) => x as u16,
                            Value::Int32(x) => x as u16,
                            Value::UInt32(x) => x as u16,
                            _ => 1,
                        };
                    }
                }
            }
        }
        let chrom = rec.reference_sequence_id().transpose()?.unwrap_or(0) as u32;
        let pos = rec
            .alignment_start()
            .transpose()?
            .map(|p| usize::from(p) as u32 - 1)
            .unwrap_or(0);
        let pl = placement_from_alignment(
            chrom, pos, flags.is_reverse_complemented(), &to_ops(rec.cigar().iter())?, 0, 0, nh,
        );
        let shape_id = {
            let s = Shape::of(&pl);
            match shape_intern.get(&s) {
                Some(&id) => id,
                None => {
                    let id = shapes.len() as u32;
                    shape_intern.insert(s.clone(), id);
                    shapes.push(s);
                    id
                }
            }
        };

        let mut correct = |crb: &Vec<u8>, cyb: &Option<Vec<u8>>| -> Option<u32> {
            use std::hash::{Hash, Hasher};
            let mut h = rustc_hash::FxHasher::default();
            crb.hash(&mut h);
            cyb.hash(&mut h);
            *bc_cache
                .entry(h.finish())
                .or_insert_with(|| corrector.correct(crb, cyb.as_deref()))
        };

        if nh > 1 {
            let key = {
                use std::hash::{Hash, Hasher};
                let mut h = rustc_hash::FxHasher::default();
                match rec.name() {
                    Some(n) => {
                        let nb: &[u8] = n.as_ref();
                        nb.hash(&mut h);
                    }
                    None => (b"" as &[u8]).hash(&mut h),
                }
                h.finish()
            };
            let e = mm.entry(key).or_insert((None, None, Vec::new()));
            if !flags.is_secondary() {
                e.0 = match (&cr_b, ur) {
                    (Some(crb), Some(u)) => correct(crb, &cy_b).map(|c| (c, u)),
                    _ => None,
                };
                e.1 = Some(e.2.len());
            }
            e.2.push((chrom, pl.start(), flags.is_reverse_complemented(), shape_id));
            continue;
        }
        if flags.is_secondary() {
            continue;
        }
        let (Some(crb), Some(u)) = (&cr_b, ur) else {
            continue;
        };
        let Some(packed) = correct(crb, &cy_b) else {
            continue;
        };
        let next = cell_intern.len() as u32;
        let cell = *cell_intern.entry(packed).or_insert_with(|| {
            cells.push(packed);
            next
        });
        ureads.push(URead {
            cell, umi: u, chrom, pos: pl.start(), strand_rev: flags.is_reverse_complemented(), shape: shape_id,
        });
    }

    // All BAM records have now been reduced to owned row inputs. Hash the same held file while
    // the CPU-heavy deterministic grouping below runs, and validate its pre-parse metadata at the
    // end. This closes the report-only pathname race without changing the legacy extraction path.
    let identity_thread = match source_before {
        Some(before) => {
            let mut bgzf = reader.into_inner();
            let file = bgzf.finish()?;
            let path = bam_path.clone();
            Some(std::thread::spawn(move || {
                identity_of_consumed_file(file, before, &path)
            }))
        }
        None => None,
    };

    // Multimapper reads become (cell, umi, anchor, pattern) tuples; the pattern is interned.
    let mut pattern_intern: FxHashMap<Vec<PatAlt>, u32> = FxHashMap::default();
    let mut patterns: Vec<Vec<PatAlt>> = Vec::new();
    struct MRead {
        cell: u32,
        umi: u32,
        chrom: u32,
        pos: u32,
        strand_rev: bool,
        shape: u32,
        pattern: u32,
    }
    let mut mreads: Vec<MRead> = Vec::new();
    for (_k, (cb, prim_idx, mut alts)) in mm.drain() {
        let (Some((packed_cell, u)), Some(pi)) = (cb, prim_idx) else {
            continue;
        };
        let next = cell_intern.len() as u32;
        let cell = *cell_intern.entry(packed_cell).or_insert_with(|| {
            cells.push(packed_cell);
            next
        });
        let (a_chrom, a_pos, a_rev, a_shape) = alts[pi];
        alts.sort_unstable_by_key(|&(c, p, _, _)| (c, p));
        let pat: Vec<PatAlt> = alts
            .iter()
            .map(|&(c, p, r, s)| PatAlt {
                chrom: c,
                offset: p as i64 - a_pos as i64,
                strand_flip: r != a_rev,
                shape: if s == a_shape { SAME_SHAPE } else { s },
            })
            .collect();
        let pid = match pattern_intern.get(&pat) {
            Some(&id) => id,
            None => {
                let id = patterns.len() as u32;
                pattern_intern.insert(pat.clone(), id);
                patterns.push(pat);
                id
            }
        };
        mreads.push(MRead { cell, umi: u, chrom: a_chrom, pos: a_pos, strand_rev: a_rev, shape: a_shape, pattern: pid });
    }

    // ---- Class assignment + chain reduction, per (cell, chrom, strand) in sorted order. ----
    // Class ids are global, assigned in first-occurrence order of this deterministic scan; the
    // same scan order is what the collapse's tie-break relies on.
    ureads.sort_unstable_by_key(|r| (r.cell, r.chrom, r.strand_rev, r.pos));
    mreads.sort_unstable_by_key(|r| (r.cell, r.chrom, r.strand_rev, r.pos));

    let mut mols: Vec<MolRec> = Vec::new();
    let mut edges: Vec<(u32, u32)> = Vec::new();
    let mut n_classes = 0u32;
    // Global class map: (cell, umi value) -> class id. Value semantics, not locus semantics.
    let mut global_classes: FxHashMap<(u32, u32), u32> = FxHashMap::default();
    // Per cell: the set of values seen, for cell-scoped 1MM edge enumeration at the end.
    let mut cell_values: FxHashMap<u32, Vec<u32>> = FxHashMap::default();

    let mut mread_cursor = 0usize;
    let mut i = 0usize;
    while i < ureads.len() || mread_cursor < mreads.len() {
        // Advance over the group key present in whichever stream comes first.
        let key = {
            let ku = ureads.get(i).map(|r| (r.cell, r.chrom, r.strand_rev));
            let km = mreads.get(mread_cursor).map(|r| (r.cell, r.chrom, r.strand_rev));
            match (ku, km) {
                (Some(a), Some(b)) => a.min(b),
                (Some(a), None) => a,
                (None, Some(b)) => b,
                (None, None) => break,
            }
        };
        let (cell, chrom, strand_rev) = key;

        let ustart = i;
        while i < ureads.len() && (ureads[i].cell, ureads[i].chrom, ureads[i].strand_rev) == key {
            i += 1;
        }
        let mstart = mread_cursor;
        while mread_cursor < mreads.len()
            && (mreads[mread_cursor].cell, mreads[mread_cursor].chrom, mreads[mread_cursor].strand_rev) == key
        {
            mread_cursor += 1;
        }
        let ugroup = &ureads[ustart..i];
        let mgroup = &mreads[mstart..mread_cursor];

        let mut class_of = |u: u32,
                            _pos: u32,
                            global_classes: &mut FxHashMap<(u32, u32), u32>,
                            cell_values: &mut FxHashMap<u32, Vec<u32>>,
                            n_classes: &mut u32|
         -> u32 {
            match global_classes.get(&(cell, u)) {
                Some(&id) => id,
                None => {
                    let id = *n_classes;
                    *n_classes += 1;
                    global_classes.insert((cell, u), id);
                    cell_values.entry(cell).or_default().push(u);
                    id
                }
            }
        };

        // Unique reads: loci by single linkage, then per (umi, chain) the two span extremes.
        let mut ls = 0usize;
        for k in 1..=ugroup.len() {
            let split = k == ugroup.len() || ugroup[k].pos - ugroup[k - 1].pos > locus_gap;
            if !split {
                continue;
            }
            let locus = &ugroup[ls..k];
            ls = k;
            // (umi, chain) -> (contained idx, extended idx, count)
            let mut chains: FxHashMap<(u32, u64), (usize, usize, u32)> = FxHashMap::default();
            for (li, r) in locus.iter().enumerate() {
                let sh = &shapes[r.shape as usize];
                let ch = chain_hash(r.pos, sh);
                let end = end_of(r.pos, sh);
                let e = chains.entry((r.umi, ch)).or_insert((li, li, 0));
                e.2 += 1;
                let (cin, cex) = (&locus[e.0], &locus[e.1]);
                let cin_end = end_of(cin.pos, &shapes[cin.shape as usize]);
                let cex_end = end_of(cex.pos, &shapes[cex.shape as usize]);
                if (r.pos, std::cmp::Reverse(end)) > (cin.pos, std::cmp::Reverse(cin_end)) {
                    e.0 = li;
                }
                if (r.pos, std::cmp::Reverse(end)) < (cex.pos, std::cmp::Reverse(cex_end)) {
                    e.1 = li;
                }
            }
            // One MolRec per UMI in this locus, gathering all its chains.
            let mut per_umi: FxHashMap<u32, SmallVec<[MolChain; 1]>> = FxHashMap::default();
            let mut items: Vec<((u32, u64), (usize, usize, u32))> = chains.into_iter().collect();
            items.sort_unstable_by_key(|((u, ch), _)| (*u, *ch));
            for ((u, _ch), (ci, ei, w)) in items {
                // Span-minimum (extended) rep first: with chains position-sorted below, the
                // molecule's first stored rep is then its anchor and the serialization elides it.
                let mut reps = SmallVec::new();
                reps.push((locus[ei].pos, locus[ei].shape));
                if ei != ci {
                    reps.push((locus[ci].pos, locus[ci].shape));
                }
                per_umi.entry(u).or_default().push(MolChain { weight: w, reps });
            }
            let mut umis: Vec<(u32, SmallVec<[MolChain; 1]>)> = per_umi.into_iter().collect();
            umis.sort_unstable_by_key(|(u, _)| *u);
            for (u, mut chains) in umis {
                // Position-canonical chain order. Within a chain reps[0] is the span minimum, so
                // after this sort chains[0].reps[0].0 IS the molecule anchor — the serialization
                // elides that first representative's offset with no flag.
                chains.sort_unstable_by_key(|c| (c.reps[0], c.reps.last().copied(), c.weight));
                let anchor = chains.iter().flat_map(|c| c.reps.iter().map(|r| r.0)).min().unwrap_or(0);
                let cls = class_of(u, anchor, &mut global_classes, &mut cell_values, &mut n_classes);
                mols.push(MolRec { cell, umi_class: cls, chrom, strand_rev, chains, mms: SmallVec::new() });
            }
        }

        // Multimapper signatures: aggregate identical (umi, anchor, pattern) tuples.
        let mut agg: FxHashMap<(u32, u32, u32, u32), u32> = FxHashMap::default();
        for r in mgroup {
            *agg.entry((r.umi, r.pos, r.shape, r.pattern)).or_insert(0) += 1;
        }
        let mut magg: Vec<((u32, u32, u32, u32), u32)> = agg.into_iter().collect();
        magg.sort_unstable();
        let mut per_umi_mm: FxHashMap<u32, Vec<(u32, u32, u32, u32)>> = FxHashMap::default();
        for ((u, pos, shape, pattern), w) in magg {
            per_umi_mm.entry(u).or_default().push((pos, shape, pattern, w));
        }
        let mut mm_umis: Vec<(u32, Vec<(u32, u32, u32, u32)>)> = per_umi_mm.into_iter().collect();
        mm_umis.sort_unstable_by_key(|(u, _)| *u);
        for (u, mms) in mm_umis {
            let anchor = mms.iter().map(|m| m.0).min().unwrap_or(0);
            let cls = class_of(u, anchor, &mut global_classes, &mut cell_values, &mut n_classes);
            mols.push(MolRec { cell, umi_class: cls, chrom, strand_rev, chains: SmallVec::new(), mms: SmallVec::from_vec(mms) });
        }

    }

    // 1MM edges, cell-scoped over values — mirroring the value-based collapse, which enumerates a
    // value's neighbours within the (cell, gene) group with no positional cutoff. A cell-wide edge
    // that never co-occurs in one gene is simply never consulted at replay.
    for (cell, values) in &cell_values {
        let lookup: FxHashMap<u32, u32> = values
            .iter()
            .map(|v| (*v, *global_classes.get(&(*cell, *v)).unwrap()))
            .collect();
        for v in values {
            let ca = lookup[v];
            for nb in crate::graph::neighbors_1mm_pub(*v) {
                if nb <= *v {
                    continue;
                }
                if let Some(&cb2) = lookup.get(&nb) {
                    edges.push((ca.min(cb2), ca.max(cb2)));
                }
            }
        }
    }

    // Cell ids renumbered by molecule frequency (#2): common cells get small ids, shrinking the
    // per-molecule cell delta stream. The dictionary records the new order.
    {
        let mut freq: Vec<u64> = vec![0; cells.len()];
        for m in &mols {
            freq[m.cell as usize] += 1;
        }
        let mut order: Vec<u32> = (0..cells.len() as u32).collect();
        order.sort_unstable_by_key(|&c| (std::cmp::Reverse(freq[c as usize]), c));
        let mut remap: Vec<u32> = vec![0; cells.len()];
        let mut new_cells: Vec<u32> = Vec::with_capacity(cells.len());
        for (new_id, &old) in order.iter().enumerate() {
            remap[old as usize] = new_id as u32;
            new_cells.push(cells[old as usize]);
        }
        for m in mols.iter_mut() {
            m.cell = remap[m.cell as usize];
        }
        cells = new_cells;
    }

    // Serialized molecule order: genome-major by anchor. Class ids are then renumbered by first
    // occurrence in THIS order. Cell-major ids delta-code poorly against genome-major storage, so
    // renumbering is applied at molecule granularity.
    mols.sort_unstable_by_key(|m| (m.chrom, m.anchor(), m.cell, m.umi_class, m.strand_rev));
    {
        let mut remap: FxHashMap<u32, u32> = FxHashMap::default();
        for m in mols.iter_mut() {
            let next = remap.len() as u32;
            m.umi_class = *remap.entry(m.umi_class).or_insert(next);
        }
        for (a, b) in edges.iter_mut() {
            let na = *remap.get(a).expect("edge references class with no molecule");
            let nb = *remap.get(b).expect("edge references class with no molecule");
            (*a, *b) = (na.min(nb), na.max(nb));
        }
        n_classes = remap.len() as u32;
    }
    edges.sort_unstable();
    edges.dedup();

    let extracted = Extracted {
        mols,
        edges,
        cells,
        shapes,
        patterns,
        n_classes,
        chrom_names,
    };
    let identity = identity_thread
        .map(|thread| {
            thread
                .join()
                .map_err(|_| anyhow::anyhow!("BAM identity thread panicked"))?
        })
        .transpose()?;
    Ok((extracted, identity))
}

type ReplayTuple = (u32, u32, u32, u32); // (cell, class, gene, weight)
const REPLAY_SHARDS: usize = 64;

/// Bounded input-side reducer. Only assignment tuples survive each decoded archive batch; final
/// aggregation remains global, so chunk boundaries cannot affect UMI classes or counts.
pub(crate) struct ReplayRowsAccumulator<'a> {
    x: &'a Extracted,
    anno: &'a anno::Annotation,
    bam2anno: Vec<Option<u32>>,
    solo_strand: anno::assign::SoloStrand,
    shards: Vec<Vec<ReplayTuple>>,
    n_assigned: u64,
}

impl<'a> ReplayRowsAccumulator<'a> {
    pub(crate) fn with_strand(
        x: &'a Extracted,
        anno: &'a anno::Annotation,
        solo_strand: anno::assign::SoloStrand,
    ) -> Self {
        let bam2anno: Vec<Option<u32>> =
            x.chrom_names.iter().map(|n| anno.chrom_ids.get(n).copied()).collect();
        Self { x, anno, bam2anno, solo_strand,
            shards: (0..REPLAY_SHARDS).map(|_| Vec::new()).collect(), n_assigned: 0 }
    }

    pub(crate) fn reserve_assignments(&mut self, expected: usize) {
        let per_shard = expected.div_ceil(REPLAY_SHARDS);
        for shard in &mut self.shards {
            shard.reserve(per_shard);
        }
    }

    fn classify(&self, mols: &[MolRec], s: &mut RowScratch)
        -> (u64, Vec<Vec<ReplayTuple>>)
    {
        let mut shards: Vec<Vec<ReplayTuple>> =
            (0..REPLAY_SHARDS).map(|_| Vec::new()).collect();
        let mut n = 0u64;
        for m in mols {
            let mut handle = |r: Row| {
                if row_genes_stranded(
                    &r, self.x, self.anno, &self.bam2anno, MmMissing::SkipAlt,
                    self.solo_strand, s,
                ).is_some() {
                    if let [g] = s.genes.as_slice() {
                        n += 1;
                        shards[r.cell as usize % REPLAY_SHARDS]
                            .push((r.cell, r.umi_class, *g, r.weight));
                    }
                }
            };
            for ch in &m.chains {
                for (pos, shape) in &ch.reps {
                    handle(Row { cell: m.cell, umi_class: m.umi_class, chrom: m.chrom, pos: *pos,
                        strand_rev: m.strand_rev, shape: *shape, weight: ch.weight,
                        pattern: u32::MAX });
                }
            }
            for (pos, shape, pattern, weight) in &m.mms {
                handle(Row { cell: m.cell, umi_class: m.umi_class, chrom: m.chrom, pos: *pos,
                    strand_rev: m.strand_rev, shape: *shape, weight: *weight, pattern: *pattern });
            }
        }
        (n, shards)
    }

    fn append_parts(&mut self, mut parts: Vec<(u64, Vec<Vec<ReplayTuple>>)>) {
        self.n_assigned += parts.iter().map(|(n, _)| *n).sum::<u64>();
        for k in 0..REPLAY_SHARDS {
            self.shards[k].reserve(parts.iter().map(|(_, p)| p[k].len()).sum());
        }
        for (_, worker_shards) in &mut parts {
            for (dst, src) in self.shards.iter_mut().zip(worker_shards) {
                dst.append(src);
            }
        }
    }

    /// Fine reducer tasks let work stealing smooth gene-density differences within the fixed
    /// decode window without concatenating molecule records.
    pub(crate) fn add_archive_chunks(&mut self, chunks: &[Vec<MolRec>]) {
        let tasks: Vec<&[MolRec]> = chunks.iter()
            .flat_map(|chunk| chunk.chunks(1 << 16)).collect();
        let parts: Vec<_> = tasks.into_par_iter()
            .map_init(RowScratch::default, |s, part| self.classify(part, s)).collect();
        self.append_parts(parts);
    }

    /// Preserve the original parallel shard assembly for a performance-faithful `--eager` arm.
    fn add_molecules_eager(&mut self, mols: &[MolRec]) {
        let parts: Vec<_> = mols.par_chunks(1 << 19)
            .map_init(RowScratch::default, |s, part| self.classify(part, s)).collect();
        self.n_assigned += parts.iter().map(|(n, _)| *n).sum::<u64>();
        let mut shards: Vec<Vec<ReplayTuple>> = (0..REPLAY_SHARDS)
            .map(|k| Vec::with_capacity(parts.iter().map(|(_, p)| p[k].len()).sum())).collect();
        shards.par_iter_mut().enumerate().for_each(|(k, dst)| {
            for (_, worker_shards) in &parts {
                dst.extend_from_slice(&worker_shards[k]);
            }
        });
        self.shards = shards;
    }

    pub(crate) fn finish(self) -> (FxHashMap<(u32, u32), u32>, u64, u64) {
        let Self {
            x,
            shards,
            n_assigned,
            ..
        } = self;

        // Adjacency over classes as CSR arrays: deterministic, index-addressed, zero hashing. The
        // hash-map version of this stage (nested maps keyed by (cell, gene)/(cell, class)) was 50%+
        // of the replay profile — sorting the shard makes every group a contiguous run instead.
        let nc = x.n_classes as usize;
        let mut off = vec![0u32; nc + 1];
        for (a, b) in &x.edges {
            off[*a as usize + 1] += 1;
            off[*b as usize + 1] += 1;
        }
        for i in 0..nc {
            off[i + 1] += off[i];
        }
        let mut cur: Vec<u32> = off[..nc].to_vec();
        let mut nbrs = vec![0u32; *off.last().unwrap() as usize];
        for (a, b) in &x.edges {
            // Same per-node neighbour order as the old push-per-edge map build; only the min rank of
            // the neighbours is consumed, so order is immaterial anyway.
            nbrs[cur[*a as usize] as usize] = *b;
            cur[*a as usize] += 1;
            nbrs[cur[*b as usize] as usize] = *a;
            cur[*b as usize] += 1;
        }
        drop(cur);

        let counted: Vec<((u32, u32), u32)> = shards
            .into_par_iter()
            .flat_map_iter(|mut shard| {
                // Sorting groups (cell, class) runs with genes as sub-runs; weight sums, the
                // MultiGeneUMI_CR best-gene choice and the per-(cell, gene) collapse all read off
                // contiguous slices. Semantics identical to the old nested-map pipeline.
                shard.sort_unstable();
                // Per (cell, class): best gene by (weight total, Reverse(gene)); that class then
                // counts only toward its best gene, carrying that gene's weight total.
                let mut kept: Vec<(u32, u32, u32, u32)> = Vec::with_capacity(shard.len()); // (cell, gene, cls, n)
                let mut i = 0usize;
                while i < shard.len() {
                    let (cell, cls) = (shard[i].0, shard[i].1);
                    let (mut best_g, mut best_w, mut have) = (0u32, 0u32, false);
                    while i < shard.len() && shard[i].0 == cell && shard[i].1 == cls {
                        let gene = shard[i].2;
                        let mut w = 0u32;
                        while i < shard.len()
                            && shard[i].0 == cell
                            && shard[i].1 == cls
                            && shard[i].2 == gene
                        {
                            w += shard[i].3;
                            i += 1;
                        }
                        if !have || w > best_w || (w == best_w && gene < best_g) {
                            (best_g, best_w, have) = (gene, w, true);
                        }
                    }
                    kept.push((cell, best_g, cls, best_w));
                }
                // Group by (cell, gene); tie-merging greedy over (class, count), most-abundant first,
                // class-id tie order, transitive to the earliest survivor — collapse_nodes' semantics
                // minus UMI values. Rank lookup is a binary search over the (small) group.
                kept.sort_unstable();
                let mut out: Vec<((u32, u32), u32)> = Vec::new();
                let mut order: Vec<(u32, u32)> = Vec::new(); // (cls, n), reused
                let mut byclass: Vec<(u32, u32)> = Vec::new(); // (cls, rank), reused
                let mut canon: Vec<u32> = Vec::new(); // canonical rank per rank, reused
                let mut g0 = 0usize;
                while g0 < kept.len() {
                    let (cell, gene) = (kept[g0].0, kept[g0].1);
                    let mut g1 = g0;
                    while g1 < kept.len() && kept[g1].0 == cell && kept[g1].1 == gene {
                        g1 += 1;
                    }
                    order.clear();
                    order.extend(kept[g0..g1].iter().map(|k| (k.2, k.3)));
                    order.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
                    byclass.clear();
                    byclass.extend(order.iter().enumerate().map(|(r, (c, _))| (*c, r as u32)));
                    byclass.sort_unstable();
                    canon.clear();
                    canon.resize(order.len(), 0);
                    let mut roots = 0u32;
                    for (idx, (c, _n)) in order.iter().enumerate() {
                        let mut tgt: Option<u32> = None;
                        for nb in &nbrs[off[*c as usize] as usize..off[*c as usize + 1] as usize] {
                            if let Ok(p) = byclass.binary_search_by_key(nb, |(c2, _)| *c2) {
                                let r = byclass[p].1;
                                if (r as usize) < idx {
                                    tgt = Some(tgt.map_or(r, |t| t.min(r)));
                                }
                            }
                        }
                        match tgt {
                            Some(r) => canon[idx] = canon[r as usize],
                            None => {
                                canon[idx] = idx as u32;
                                roots += 1;
                            }
                        }
                    }
                    out.push(((cell, gene), roots));
                    g0 = g1;
                }
                out
            })
            .collect();

    let total: u64 = counted.iter().map(|(_, n)| *n as u64).sum();
    (counted.into_iter().collect(), n_assigned, total)
    }
}

/// Quantify one annotation from rows alone. The eager entry point remains the reference used by
/// BAM, molecule-BAM, EM, velocity helpers, and `replay-rows --eager`.
pub fn replay_rows(
    x: &Extracted,
    anno: &anno::Annotation,
) -> (FxHashMap<(u32, u32), u32>, u64, u64) {
    replay_rows_stranded(x, anno, anno::assign::SoloStrand::Forward)
}

/// Strand-parameterized [`replay_rows`]. Existing callers retain `Forward` through the wrapper;
/// the CLI selects this explicitly for non-3' chemistries.
pub fn replay_rows_stranded(
    x: &Extracted,
    anno: &anno::Annotation,
    solo_strand: anno::assign::SoloStrand,
) -> (FxHashMap<(u32, u32), u32>, u64, u64) {
    let mut replay = ReplayRowsAccumulator::with_strand(x, anno, solo_strand);
    replay.add_molecules_eager(&x.mols);
    replay.finish()
}
