//! Build the annotation-independent molecule archive from an annotation-free alignment, and
//! report where its bytes go.
//!
//! Molecules are keyed by `(cell, UMI, locus, orientation)` and carry only evidence explained by
//! the reference genome. No gene model is consulted while they are built.
//!
//! Read starts are grouped by single-linkage clustering within `(cell, chrom, strand)`, using a
//! 50 kb maximum gap, followed by one-mismatch UMI correction inside each locus. Fixed-width
//! bucketing is deliberately avoided because it splits molecules that straddle a boundary.
//!
//! Barcodes are corrected here rather than taken from STARsolo, which refuses to emit `CB` without
//! a gene model. Whitelist matching is exact-or-one-mismatch and consults no annotation, so ingest
//! remains annotation-free end to end.
//!
//! The archive is encoded twice — with and without the UMI stream — because the UMI is over half
//! the record and is needed only if a future annotation must redo the per-gene collapse.

use anyhow::{Context, Result};
use clap::Parser;
use evidence_io::{archive, umi, Ev, Molecule};
use ingest::{cigar::Op, placement_from_alignment};
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Value;
use noodles_sam::alignment::record::Flags;
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::PathBuf;

#[derive(Parser)]
pub struct Args {
    /// Annotation-free ingest BAM, coordinate-sorted, carrying raw `CR`/`UR`.
    pub bam: PathBuf,
    /// 10x barcode whitelist (one 16 bp barcode per line).
    #[arg(long)]
    pub whitelist: PathBuf,
    /// Locus gap in bp for single-linkage molecule grouping.
    #[arg(long, default_value_t = 50_000)]
    pub locus_gap: u32,
    /// zstd level for the size report.
    #[arg(long, default_value_t = 19)]
    pub zstd_level: i32,
    /// Skip 1-mismatch UMI collapse at ingest, keeping one molecule per exact UMI value per locus.
    /// This is the graph-variant base: every merge decision is deferred to replay (served by the
    /// UMI-adjacency graph), so none is baked in. Slightly more molecules, strictly more replayable.
    #[arg(long)]
    pub no_umi_collapse: bool,
    #[arg(long)]
    pub json_out: Option<PathBuf>,
}

/// A read reduced to what the archive may keep.
pub(crate) struct Read {
    pub(crate) cb: u32,
    pub(crate) umi: u32,
    pub(crate) pos: u32,
    pub(crate) strand_rev: bool,
    pub(crate) ops: Vec<Op>,
    pub(crate) nm: u16,
    pub(crate) score: i32,
    pub(crate) nh: u16,
}

pub(crate) fn to_ops(
    cigar: impl Iterator<Item = std::io::Result<noodles_sam::alignment::record::cigar::Op>>,
) -> Result<Vec<Op>> {
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

/// Load the whitelist as 2-bit packed 16-mers.
pub(crate) fn load_whitelist(path: &PathBuf) -> Result<FxHashSet<u32>> {
    let txt = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut set = FxHashSet::default();
    for line in txt.lines() {
        let b = line.trim().as_bytes();
        if b.len() == 16 {
            if let Some(p) = umi::pack(b) {
                set.insert(p);
            }
        }
    }
    Ok(set)
}

/// STARsolo-style `1MM_multi_Nbase_pseudocounts` barcode correction.
///
/// The simple exact-or-unique-1MM rule (below) DROPS any read whose barcode sits within one
/// substitution of two whitelist entries. STARsolo instead resolves the ambiguity by posterior:
/// candidate weight = (exact-match read count of that whitelist barcode + 1 pseudocount) ×
/// P(sequencing error at the mismatched base, from its quality). The winner is accepted when it
/// holds > 0.975 of the total posterior — CellRanger's threshold. A barcode with a single `N` is
/// treated as a free mismatch at that position. The whole procedure consults only the whitelist
/// and the run's own barcode frequencies: annotation-free.
///
/// On the evaluated PBMC data, the simple rule accounted for an additional 0.68 percentage points
/// of changed UMI mass. This posterior-based procedure replaces that rule.
pub(crate) struct BcCorrector {
    wl: FxHashSet<u32>,
    /// Exact-whitelist-hit read counts from a first pass over the run.
    counts: FxHashMap<u32, u32>,
}

impl BcCorrector {
    pub(crate) fn new(wl: FxHashSet<u32>, counts: FxHashMap<u32, u32>) -> Self {
        BcCorrector { wl, counts }
    }

    pub(crate) fn wl_contains(&self, packed: u32) -> bool {
        self.wl.contains(&packed)
    }

    fn weight(&self, cand: u32, perr: f64) -> f64 {
        (*self.counts.get(&cand).unwrap_or(&0) as f64 + 1.0) * perr
    }

    /// `raw` is the 16 bp barcode string; `qual` its Phred+33 qualities when available.
    pub(crate) fn correct(&self, raw: &[u8], qual: Option<&[u8]>) -> Option<u32> {
        if raw.len() != 16 {
            return None;
        }
        let n_positions: Vec<usize> = raw.iter().enumerate().filter(|(_, b)| !matches!(b, b'A'|b'C'|b'G'|b'T'|b'a'|b'c'|b'g'|b't')).map(|(i, _)| i).collect();
        if n_positions.len() > 1 {
            return None;
        }
        let perr_at = |i: usize| -> f64 {
            match qual.and_then(|q| q.get(i)) {
                // Phred+33; cap as STARsolo does so a Q2 base cannot dominate the posterior.
                Some(&q) => 10f64.powf(-(((q.saturating_sub(33)).min(33)) as f64) / 10.0),
                None => 0.01,
            }
        };

        if let [npos] = n_positions.as_slice() {
            // One N: try the four bases there; whitelist membership decides the candidate set.
            let mut cands: Vec<(u32, f64)> = Vec::new();
            let mut fixed = raw.to_vec();
            for b in [b'A', b'C', b'G', b'T'] {
                fixed[*npos] = b;
                if let Some(p) = umi::pack(&fixed) {
                    if self.wl.contains(&p) {
                        cands.push((p, self.weight(p, 1.0)));
                    }
                }
            }
            return pick(&cands, 0.975);
        }

        let packed = umi::pack(raw)?;
        if self.wl.contains(&packed) {
            return Some(packed);
        }
        // One substitution at any position; posterior over whitelist neighbours.
        let mut cands: Vec<(u32, f64)> = Vec::new();
        for str_pos in 0..16usize {
            let shift = 2 * (15 - str_pos);
            let cur = (packed >> shift) & 0b11;
            for base in 0..4u32 {
                if base == cur {
                    continue;
                }
                let cand = (packed & !(0b11 << shift)) | (base << shift);
                if self.wl.contains(&cand) {
                    cands.push((cand, self.weight(cand, perr_at(str_pos))));
                }
            }
        }
        pick(&cands, 0.975)
    }
}

fn pick(cands: &[(u32, f64)], threshold: f64) -> Option<u32> {
    let total: f64 = cands.iter().map(|(_, w)| w).sum();
    if total <= 0.0 {
        return None;
    }
    let best = cands.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())?;
    if best.1 / total > threshold {
        Some(best.0)
    } else {
        None
    }
}

/// Exact-or-one-substitution whitelist match. Returns the corrected barcode, or `None` when the
/// read is ambiguous or unmatched — those reads are dropped rather than guessed at, since a
/// mis-corrected barcode silently moves counts between cells.
pub(crate) fn correct_barcode(raw: u32, wl: &FxHashSet<u32>) -> Option<u32> {
    if wl.contains(&raw) {
        return Some(raw);
    }
    let mut hit = None;
    for pos in 0..16 {
        let shift = 2 * pos;
        let cur = (raw >> shift) & 0b11;
        for base in 0..4u32 {
            if base == cur {
                continue;
            }
            let cand = (raw & !(0b11 << shift)) | (base << shift);
            if wl.contains(&cand) {
                if hit.is_some() {
                    return None; // ambiguous within one substitution
                }
                hit = Some(cand);
            }
        }
    }
    hit
}

fn hamming1_u32(a: u32, b: u32, len: usize) -> bool {
    let mut diff = 0;
    for pos in 0..len {
        let shift = 2 * pos;
        if (a >> shift) & 0b11 != (b >> shift) & 0b11 {
            diff += 1;
            if diff > 1 {
                return false;
            }
        }
    }
    diff == 1
}

/// Collapse a locus's UMIs, absorbing each into a strictly more abundant neighbour within one
/// substitution. Ties break on the UMI value so the result is order-independent.
fn collapse_umis(counts: &FxHashMap<u32, u32>) -> FxHashMap<u32, u32> {
    let mut order: Vec<(u32, u32)> = counts.iter().map(|(u, n)| (*u, *n)).collect();
    order.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let mut map: FxHashMap<u32, u32> = FxHashMap::default();
    for (u, n) in &order {
        let mut target = None;
        for (cand, cn) in &order {
            if cand == u {
                break;
            }
            if cn > n && hamming1_u32(*cand, *u, 12) {
                target = Some(*map.get(cand).unwrap_or(cand));
                break;
            }
        }
        map.insert(*u, target.unwrap_or(*u));
    }
    map
}

/// Turn one chromosome's reads into molecules.
/// Placement of a single read, for read-level (STARsolo-semantics) replay.
pub(crate) fn placement_of(chrom: u32, r: &Read) -> evidence_io::Placement {
    placement_from_alignment(chrom, r.pos, r.strand_rev, &r.ops, r.nm, r.score, r.nh)
}

/// Placement from loose parts, for records that never become a `Read` (multimapper secondaries).
#[allow(clippy::too_many_arguments)]
pub(crate) fn placement_from_parts(
    chrom: u32, pos: u32, rev: bool, ops: &[Op], nm: u16, score: i32, nh: u16,
) -> evidence_io::Placement {
    placement_from_alignment(chrom, pos, rev, ops, nm, score, nh)
}

pub(crate) fn molecules_for_chrom(chrom: u32, reads: Vec<Read>, gap: u32, collapse: bool, out: &mut Vec<Molecule>) {
    // Bucket by (cell, strand); positions within a bucket are then clustered by single linkage.
    let mut by_cell: FxHashMap<(u32, bool), Vec<Read>> = FxHashMap::default();
    for r in reads {
        by_cell.entry((r.cb, r.strand_rev)).or_default().push(r);
    }

    for ((cb, strand_rev), mut rs) in by_cell {
        rs.sort_by_key(|r| r.pos);
        let mut start = 0usize;
        for i in 1..=rs.len() {
            let split = i == rs.len() || rs[i].pos - rs[i - 1].pos > gap;
            if !split {
                continue;
            }
            let locus = &rs[start..i];
            start = i;

            let canon = if collapse {
                let mut counts: FxHashMap<u32, u32> = FxHashMap::default();
                for r in locus {
                    *counts.entry(r.umi).or_insert(0) += 1;
                }
                collapse_umis(&counts)
            } else {
                FxHashMap::default() // identity: every exact value is its own molecule
            };

            // One molecule per corrected UMI. Its placement is taken from the read with the
            // longest alignment, and the count of reads collapsed into it is retained: that ratio
            // is the entire storage thesis.
            let mut groups: FxHashMap<u32, (&Read, u32)> = FxHashMap::default();
            for r in locus {
                let key = *canon.get(&r.umi).unwrap_or(&r.umi);
                let e = groups.entry(key).or_insert((r, 0));
                e.1 += 1;
                let cur_len: u32 = e.0.ops.iter().map(|o| match o {
                    Op::Match(n) | Op::Del(n) => *n,
                    _ => 0,
                }).sum();
                let new_len: u32 = r.ops.iter().map(|o| match o {
                    Op::Match(n) | Op::Del(n) => *n,
                    _ => 0,
                }).sum();
                // Representative selection is surprisingly load-bearing: the LONGEST read is the
                // one most likely to protrude past an annotated transcript end and fail STAR's
                // containment check, silently dropping the whole molecule at replay. Prefer the
                // SHORTEST — most likely to be contained — pending a proper consensus scheme.
                if new_len < cur_len {
                    e.0 = r;
                }
            }

            for (u, (rep, n_reads)) in groups {
                let p = placement_from_alignment(
                    chrom, rep.pos, strand_rev, &rep.ops, rep.nm, rep.score, rep.nh,
                );
                out.push(Molecule {
                    cb,
                    umi: u,
                    placements: vec![p],
                    n_reads,
                    residual: None,
                });
            }
        }
    }
}

pub fn run(args: Args) -> Result<()> {
    let wl = load_whitelist(&args.whitelist)?;
    eprintln!("whitelist: {} barcodes", wl.len());

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&args.bam)
        .with_context(|| format!("opening {}", args.bam.display()))?;
    reader.read_header()?;

    let mut molecules: Vec<Molecule> = Vec::new();
    // Barcodes are interned to dense indices: the 2-bit packed barcode is 32 essentially random
    // bits per molecule, whereas an index into the ~1M barcodes actually observed is ~20 bits and
    // delta-codes far better. The dictionary is a one-off cost of 4 bytes per distinct barcode.
    let mut cb_intern: FxHashMap<u32, u32> = FxHashMap::default();
    let mut rec = bam::Record::default();
    let (mut n_reads, mut n_kept, mut n_bc_fail) = (0u64, 0u64, 0u64);

    // The BAM is coordinate-sorted, so reads can be accumulated one chromosome at a time and
    // flushed, which bounds peak memory to a single chromosome's reads rather than the whole run.
    let mut cur_chrom: i32 = -2;
    let mut buf: Vec<Read> = Vec::new();
    let mut bc_cache: FxHashMap<u32, Option<u32>> = FxHashMap::default();

    loop {
        let more = reader.read_record(&mut rec)? != 0;
        let this_chrom = if more {
            let f = Flags::from(rec.flags().bits());
            if f.is_secondary() || f.is_supplementary() || f.is_unmapped() {
                continue;
            }
            rec.reference_sequence_id().transpose()?.unwrap_or(0) as i32
        } else {
            -1
        };

        if this_chrom != cur_chrom && cur_chrom >= 0 {
            molecules_for_chrom(cur_chrom as u32, std::mem::take(&mut buf), args.locus_gap, !args.no_umi_collapse, &mut molecules);
        }
        if !more {
            break;
        }
        cur_chrom = this_chrom;
        n_reads += 1;

        let (mut cr, mut ur) = (None, None);
        for field in rec.data().iter() {
            let (tag, value) = field?;
            let key = <[u8; 2]>::from(tag);
            if let Value::String(s) = value {
                match &key {
                    b"CR" => cr = umi::pack(s.as_ref()),
                    b"UR" => ur = umi::pack(s.as_ref()),
                    _ => {}
                }
            }
        }
        let (Some(cr), Some(ur)) = (cr, ur) else { continue };

        let cb_packed = match *bc_cache.entry(cr).or_insert_with(|| correct_barcode(cr, &wl)) {
            Some(c) => c,
            None => {
                n_bc_fail += 1;
                continue;
            }
        };
        let next_id = cb_intern.len() as u32;
        let cb = *cb_intern.entry(cb_packed).or_insert(next_id);

        let flags = Flags::from(rec.flags().bits());
        let pos = rec.alignment_start().transpose()?.map(|p| usize::from(p) as u32 - 1).unwrap_or(0);
        let mut nm = 0u16;
        let mut score = 0i32;
        let mut nh = 1u16;
        for field in rec.data().iter() {
            let (tag, value) = field?;
            let k = <[u8; 2]>::from(tag);
            let v = match value {
                Value::Int8(v) => v as i64,
                Value::UInt8(v) => v as i64,
                Value::Int16(v) => v as i64,
                Value::UInt16(v) => v as i64,
                Value::Int32(v) => v as i64,
                Value::UInt32(v) => v as i64,
                _ => continue,
            };
            match &k {
                b"nM" => nm = v as u16,
                b"AS" => score = v as i32,
                b"NH" => nh = v as u16,
                _ => {}
            }
        }

        n_kept += 1;
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

    let n_cells = cb_intern.len() as u64;
    // Dictionary cost of interning, charged honestly against the archive.
    let cb_dict_bytes = 4 * n_cells;
    let n_mol = molecules.len() as u64;
    let total_reads_in_mol: u64 = molecules.iter().map(|m| m.n_reads as u64).sum();

    println!("=== archive build ===");
    println!("mapped primary reads:     {n_reads}");
    println!("reads with usable barcode:{n_kept}");
    println!("reads dropped (barcode):  {n_bc_fail}");
    println!("distinct barcodes:        {n_cells}  (dictionary {cb_dict_bytes} B)");
    println!("molecules:                {n_mol}");
    println!("reads per molecule:       {:.3}", total_reads_in_mol as f64 / n_mol.max(1) as f64);
    println!();

    // E1 retains the assignment-relevant geometry available from the alignment; E0 is reported
    // alongside it to show the encoded cost of splice structure.
    let mut report = serde_json::Map::new();
    report.insert("molecules".into(), serde_json::json!(n_mol));
    report.insert("distinct_barcodes".into(), serde_json::json!(n_cells));
    report.insert("cb_dictionary_bytes".into(), serde_json::json!(cb_dict_bytes));
    report.insert("mapped_primary_reads".into(), serde_json::json!(n_reads));
    report.insert("reads_per_molecule".into(),
                  serde_json::json!(total_reads_in_mol as f64 / n_mol.max(1) as f64));

    for (level, name) in [(Ev::E0, "E0"), (Ev::E1, "E1")] {
        let mut projected: Vec<Molecule> = molecules.iter().map(|m| m.project(level)).collect();
        for (with_umi, tag) in [(true, "with_umi"), (false, "no_umi")] {
            let r = archive::encode(&mut projected, with_umi, args.zstd_level)?;
            let total_with_dict = r.total + cb_dict_bytes;
            println!(
                "{name:<3} {tag:<9} total {:>12} B  (+dict {:>12} B = {:.2} bits/molecule, {} shapes)",
                r.total, total_with_dict,
                8.0 * total_with_dict as f64 / n_mol.max(1) as f64,
                r.distinct_shapes
            );
            println!(
                "      streams: pos {} shape {} flags {} cb {} umi {} nreads {} chrom {} dict {}",
                r.pos_delta, r.shape_id, r.flags, r.cb_delta, r.umi, r.n_reads, r.chrom, r.shape_dict
            );
            report.insert(format!("{name}_{tag}"), serde_json::to_value(&r)?);
        }
    }

    if let Some(p) = &args.json_out {
        if let Some(d) = p.parent() {
            std::fs::create_dir_all(d)?;
        }
        std::fs::write(p, serde_json::to_string_pretty(&report)?)?;
        println!("wrote {}", p.display());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_whitelist_match_is_returned_unchanged() {
        let mut wl = FxHashSet::default();
        let bc = umi::pack(b"ACGTACGTACGTACGT").unwrap();
        wl.insert(bc);
        assert_eq!(correct_barcode(bc, &wl), Some(bc));
    }

    #[test]
    fn one_substitution_is_corrected() {
        let mut wl = FxHashSet::default();
        let good = umi::pack(b"ACGTACGTACGTACGT").unwrap();
        wl.insert(good);
        let mutated = umi::pack(b"ACGTACGTACGTACGA").unwrap();
        assert_eq!(correct_barcode(mutated, &wl), Some(good));
    }

    #[test]
    fn ambiguous_barcode_is_dropped_not_guessed() {
        // Two whitelist entries one substitution away: correcting either way would silently move
        // counts between cells, so the read must be discarded.
        let mut wl = FxHashSet::default();
        wl.insert(umi::pack(b"ACGTACGTACGTACGA").unwrap());
        wl.insert(umi::pack(b"ACGTACGTACGTACGC").unwrap());
        let mutated = umi::pack(b"ACGTACGTACGTACGT").unwrap();
        assert_eq!(correct_barcode(mutated, &wl), None);
    }

    #[test]
    fn two_substitutions_are_not_corrected() {
        let mut wl = FxHashSet::default();
        wl.insert(umi::pack(b"ACGTACGTACGTACGT").unwrap());
        let mutated = umi::pack(b"ACGTACGTACGTACAA").unwrap();
        assert_eq!(correct_barcode(mutated, &wl), None);
    }

    #[test]
    fn umi_collapse_absorbs_rare_into_abundant() {
        let a = umi::pack(b"AAAAAAAAAAAA").unwrap();
        let b = umi::pack(b"AAAAAAAAAAAC").unwrap();
        let mut c = FxHashMap::default();
        c.insert(a, 100);
        c.insert(b, 1);
        let m = collapse_umis(&c);
        assert_eq!(m[&b], a);
        assert_eq!(m[&a], a);
    }

    #[test]
    fn single_linkage_splits_only_on_a_real_gap() {
        // Two reads 60 kb apart with a 50 kb gap are distinct molecules even with the same UMI;
        // two reads 100 bp apart are one. This is the behaviour fixed-width bucketing gets wrong.
        let mk = |pos| Read {
            cb: 1, umi: 7, pos, strand_rev: false,
            ops: vec![Op::Match(91)], nm: 0, score: 0, nh: 1,
        };
        let mut out = Vec::new();
        molecules_for_chrom(0, vec![mk(1000), mk(1100)], 50_000, true, &mut out);
        assert_eq!(out.len(), 1, "nearby reads of one UMI are one molecule");
        assert_eq!(out[0].n_reads, 2);

        out.clear();
        molecules_for_chrom(0, vec![mk(1000), mk(61_000)], 50_000, true, &mut out);
        assert_eq!(out.len(), 2, "reads beyond the gap are separate molecules");
    }
}
