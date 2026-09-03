//! `aie debug`: measurement instrument for the coding stack. Read-only; touches no writer code.
//!
//! Three questions it answers, all against a real archive:
//! 1. Where do the bytes go, section by section (raw vs zstd), and stream by stream inside chunks?
//! 2. What would an ideal memoryless *value* coder achieve per stream (order-0 entropy over the
//!    decoded varint values)? zstd below that bound means it exploits cross-value structure; zstd
//!    far above it means the entropy layer is leaving money; both are model headroom signals.
//! 3. What sharing ratios do the candidate factorizations have — junction-chain equivalence
//!    classes over shapes, hop (intron-length) vocabulary, pattern-alt vocabulary — and what
//!    would each save in the dictionaries?

use anyhow::Result;
use clap::Args as ClapArgs;
use evidence_io::format::{Cursor, SectionReader};
use rustc_hash::FxHashMap;
use std::io::Write;
use std::path::PathBuf;

use crate::archivecmd::read_chunk_index;
use crate::rows::SAME_SHAPE;

#[derive(ClapArgs)]
pub struct Args {
    archive: PathBuf,
    /// Write each chunk stream (concatenated across chunks) and each dictionary/index section as
    /// a raw .bin into this directory, for external recompression experiments.
    #[arg(long)]
    dump_dir: Option<PathBuf>,
}

const STREAM_NAMES: [&str; 10] = [
    "anchor", "class", "layout", "weight",
    "rep.pos", "rep.shape", "mm.pos", "mm.shape", "mm.pattern", "mm.weight",
];

fn entropy_report(name: &str, hist: &FxHashMap<u64, u64>, stream_bytes: usize) {
    let total: u64 = hist.values().sum();
    if total == 0 {
        return;
    }
    let mut bits = 0f64;
    for &c in hist.values() {
        let p = c as f64 / total as f64;
        bits -= (c as f64) * p.log2();
    }
    println!(
        "  value-entropy {name}: {total} values, {} distinct, order-0 bound {:.2} MB ({:.2} bits/value); varint stream {:.2} MB",
        hist.len(),
        bits / 8.0 / 1e6,
        bits / total as f64,
        stream_bytes as f64 / 1e6
    );
}

pub fn run(args: Args) -> Result<()> {
    let mut r = SectionReader::open(&args.archive)?;

    // ---- 1. Section accounting ----
    let (mut chunk_raw, mut chunk_comp, mut n_chunks) = (0u64, 0u64, 0u64);
    let (mut coc_raw, mut coc_comp, mut n_coc) = (0u64, 0u64, 0u64);
    println!("== sections (raw -> zstd) ==");
    let entries: Vec<_> = r.entries().to_vec();
    for (name, _, raw, comp) in &entries {
        if name.starts_with('c') && name[1..].chars().all(|c| c.is_ascii_digit()) {
            chunk_raw += raw;
            chunk_comp += comp;
            n_chunks += 1;
        } else if name.starts_with("coc.") {
            coc_raw += raw;
            coc_comp += comp;
            n_coc += 1;
        } else {
            println!("  {name}: {:.3} MB -> {:.3} MB", *raw as f64 / 1e6, *comp as f64 / 1e6);
        }
    }
    println!("  coc x{n_coc}: {:.3} MB -> {:.3} MB", coc_raw as f64 / 1e6, coc_comp as f64 / 1e6);
    println!("  chunks x{n_chunks}: {:.3} MB -> {:.3} MB", chunk_raw as f64 / 1e6, chunk_comp as f64 / 1e6);

    // ---- 2. Split chunk streams ----
    let chunks = read_chunk_index(&mut r)?;
    let tables = crate::archivecmd::read_rans_tables(&mut r)?;
    // rANS-coded streams (class, weight, rep.pos, mm.pos, mm.weight) are decoded back to varints
    // per chunk so every downstream analysis keeps its value semantics.
    const RANS_AT: [(usize, usize); 5] = [(1, 0), (3, 1), (4, 2), (6, 3), (9, 4)];
    let mut streams: Vec<Vec<u8>> = vec![Vec::new(); STREAM_NAMES.len()];
    for i in 0..chunks.len() {
        let raw = r.read(&format!("c{i}"))?;
        let mut c = Cursor::new(&raw);
        for (si, s) in streams.iter_mut().enumerate() {
            let len = c.varint()? as usize;
            let start = c.position();
            let seg = &raw[start..start + len];
            if let Some((_, ti)) = RANS_AT.iter().find(|(a, _)| *a == si) {
                for v in evidence_io::rans::decode(seg, &tables[*ti])? {
                    evidence_io::archive::put_varint(s, v);
                }
            } else {
                s.extend_from_slice(seg);
            }
            c.set_position(start + len);
        }
    }
    println!("== chunk streams (raw, concatenated) ==");
    for (name, s) in STREAM_NAMES.iter().zip(&streams) {
        println!("  {name}: {:.3} MB", s.len() as f64 / 1e6);
    }

    if let Some(dir) = &args.dump_dir {
        std::fs::create_dir_all(dir)?;
        for (name, s) in STREAM_NAMES.iter().zip(&streams) {
            std::fs::File::create(dir.join(format!("stream.{name}.bin")))?.write_all(s)?;
        }
        let mut coc_all = Vec::new();
        for (name, _, _, _) in &entries {
            if name.starts_with('c') && name[1..].chars().all(|c| c.is_ascii_digit()) {
                continue;
            }
            let raw = r.read(name)?;
            if name.starts_with("coc.") {
                coc_all.extend_from_slice(&raw);
                continue;
            }
            std::fs::File::create(dir.join(format!("section.{name}.bin")))?.write_all(&raw)?;
        }
        std::fs::File::create(dir.join("section.cellofclass.bin"))?.write_all(&coc_all)?;
        println!("dumped to {}", dir.display());
    }

    // ---- 3. Value-level order-0 entropy per stream ----
    println!("== value entropy (memoryless coder bound per stream) ==");
    for (i, name) in STREAM_NAMES.iter().enumerate() {
        if *name == "layout" {
            // Repeating triple: strand byte, n_chains varint, n_mms varint. Entropy of the triple.
            let mut hist: FxHashMap<u64, u64> = FxHashMap::default();
            let mut c = Cursor::new(&streams[i]);
            while !c.is_empty() {
                let strand = c.byte()? as u64;
                let nch = c.varint()?;
                let nmm = c.varint()?;
                *hist.entry(strand | (nch << 1) | (nmm << 20)).or_default() += 1;
            }
            entropy_report(name, &hist, streams[i].len());
        } else if *name == "cell" {
            let mut hist: FxHashMap<u64, u64> = FxHashMap::default();
            let mut c = Cursor::new(&streams[i]);
            while !c.is_empty() {
                *hist.entry(c.svarint()? as u64).or_default() += 1;
            }
            entropy_report(name, &hist, streams[i].len());
        } else {
            let mut hist: FxHashMap<u64, u64> = FxHashMap::default();
            let mut c = Cursor::new(&streams[i]);
            while !c.is_empty() {
                *hist.entry(c.varint()?).or_default() += 1;
            }
            entropy_report(name, &hist, streams[i].len());
            if *name == "class" {
                let fresh = {
                    let mut c = Cursor::new(&streams[i]);
                    let mut f = 0u64;
                    let mut t = 0u64;
                    while !c.is_empty() {
                        t += 1;
                        f += (c.varint()? == 0) as u64;
                    }
                    (f, t)
                };
                println!("    class: {} of {} tokens fresh ({:.1}%)", fresh.0, fresh.1, 100.0 * fresh.0 as f64 / fresh.1 as f64);
            }
            if *name == "rep.pos" {
                let mut c = Cursor::new(&streams[i]);
                let (mut z, mut t) = (0u64, 0u64);
                while !c.is_empty() {
                    t += 1;
                    z += (c.varint()? == 0) as u64;
                }
                println!("    rep.pos: {} of {} offsets are 0 ({:.1}%) — anchor-equal representatives", z, t, 100.0 * z as f64 / t as f64);
            }
        }
    }

    // ---- 4. Factorization candidates ----
    let x = crate::archivecmd::read_archive(&args.archive)?;

    // 4a. Shapes -> junction-chain equivalence classes (terminal block lengths factored out).
    // Chain key: every block boundary in donor0-relative coordinates, except the shape's outer
    // ends. Free parameters per shape: first and last block length.
    let mut chains: FxHashMap<Vec<i64>, u64> = FxHashMap::default();
    let mut single_block = 0u64;
    let mut factored_shape_bytes = 0u64; // per-shape residual: chain id + left len + right len
    let mut current_shape_bytes = 0u64;
    for s in &x.shapes {
        let mut cur = Vec::new();
        for (off, len) in &s.blocks {
            evidence_io::archive::put_varint(&mut cur, *off as u64);
            evidence_io::archive::put_varint(&mut cur, *len as u64);
        }
        current_shape_bytes += cur.len() as u64 + 1;
        if s.blocks.len() == 1 {
            single_block += 1;
            factored_shape_bytes += 3; // chain id (unspliced class) + one length
            continue;
        }
        let donor0 = (s.blocks[0].0 + s.blocks[0].1) as i64;
        let mut key = Vec::new();
        for w in s.blocks.windows(2) {
            key.push((w[0].0 + w[0].1) as i64 - donor0); // donor
            key.push(w[1].0 as i64 - donor0); // acceptor
        }
        *chains.entry(key).or_default() += 1;
        factored_shape_bytes += 3 + 2; // ~varint chain id + left len + right len
    }
    let mut chain_dict_bytes = 0u64;
    for key in chains.keys() {
        let mut b = Vec::new();
        for v in key {
            evidence_io::archive::put_svarint(&mut b, *v);
        }
        chain_dict_bytes += b.len() as u64 + 1;
    }
    let multi = x.shapes.len() as u64 - single_block;
    println!("== factorization: shapes -> junction chains ==");
    println!(
        "  {} shapes ({} single-block, {} spliced) -> {} distinct chains ({:.2}x sharing among spliced)",
        x.shapes.len(), single_block, multi, chains.len(),
        multi as f64 / chains.len().max(1) as f64
    );
    println!(
        "  shape dict now ~{:.3} MB raw; factored ~{:.3} MB raw (chain dict {:.3} + residuals {:.3})",
        current_shape_bytes as f64 / 1e6,
        (chain_dict_bytes + factored_shape_bytes) as f64 / 1e6,
        chain_dict_bytes as f64 / 1e6,
        factored_shape_bytes as f64 / 1e6
    );

    // 4b. Hop vocabulary: intron lengths across the shape dictionary.
    let mut hops: FxHashMap<u32, u64> = FxHashMap::default();
    let mut hop_occ = 0u64;
    for s in &x.shapes {
        for w in s.blocks.windows(2) {
            *hops.entry(w[1].0 - (w[0].0 + w[0].1)).or_default() += 1;
            hop_occ += 1;
        }
    }
    println!("== factorization: hop (intron length) vocabulary ==");
    println!("  {} intron occurrences in shape dict -> {} distinct lengths ({:.2}x sharing)", hop_occ, hops.len(), hop_occ as f64 / hops.len().max(1) as f64);

    // 4c. Pattern-alt vocabulary: distinct (chrom, offset, flip, same_shape?) across all patterns.
    let mut alt_vocab: FxHashMap<(u32, i64, bool, bool), u64> = FxHashMap::default();
    let mut alt_occ = 0u64;
    let mut same_shape_alts = 0u64;
    for p in &x.patterns {
        for a in p {
            let same = a.shape == SAME_SHAPE;
            same_shape_alts += same as u64;
            *alt_vocab.entry((a.chrom, a.offset, a.strand_flip, same)).or_default() += 1;
            alt_occ += 1;
        }
    }
    let mut offsets: FxHashMap<i64, u64> = FxHashMap::default();
    for (_, off, _, _) in alt_vocab.keys() {
        *offsets.entry(*off).or_default() += 1;
    }
    println!("== factorization: pattern alternatives ==");
    println!(
        "  {} patterns, {} alt entries ({:.1}% same-shape) -> {} distinct (chrom,offset,flip) alts ({:.2}x sharing), {} distinct offsets",
        x.patterns.len(), alt_occ, 100.0 * same_shape_alts as f64 / alt_occ.max(1) as f64,
        alt_vocab.len(), alt_occ as f64 / alt_vocab.len().max(1) as f64, offsets.len()
    );

    Ok(())
}
