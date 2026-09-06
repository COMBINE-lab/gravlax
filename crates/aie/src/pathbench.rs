//! Experimental native geometry decoders. Outputs are NOT production archives.
//! Preserve molecule order, shape IDs, MM evidence and all auxiliary sections.
use crate::{
    archivecmd::{self, ChunkInfo, Dicts},
    rows::{MolChain, MolRec},
};
use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use evidence_io::{
    archive::{put_varint, Shape},
    format::{compress, decompress, Cursor, SectionReader, SectionWriter},
    rans,
};
use std::{collections::BTreeMap, path::PathBuf, time::Instant};

#[derive(ClapArgs)]
pub struct Args {
    archive: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value_t = 19)]
    level: i32,
    #[arg(long, default_value_t = 7, value_parser = clap::value_parser!(u32).range(1..=100))]
    repeats: u32,
}
type Columns = Vec<Vec<u8>>;
// Lookup only: hash iteration never influences the encoded bytes.
type ShapeMap = rustc_hash::FxHashMap<Vec<(u32, u32)>, u32>;
type Path = Vec<u32>;
type Groups = BTreeMap<Path, Vec<(u32, u32, u32)>>;
const MAX_ITEMS: usize = 16_000_000;
fn put(c: &mut [Vec<u8>], i: usize, v: u64) {
    put_varint(&mut c[i], v);
}
fn u32v(c: &mut Cursor<'_>) -> Result<u32> {
    Ok(u32::try_from(c.varint()?)?)
}
fn count(c: &mut Cursor<'_>) -> Result<usize> {
    let n = usize::try_from(c.varint()?)?;
    if n > MAX_ITEMS {
        bail!("experiment item limit exceeded");
    }
    Ok(n)
}
fn pack(cols: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for col in cols {
        put_varint(&mut out, col.len() as u64);
        out.extend(col);
    }
    out
}
fn unpack(raw: &[u8], n: usize) -> Result<Columns> {
    let mut c = Cursor::new(raw);
    let mut out = Vec::new();
    for _ in 0..n {
        let len = usize::try_from(c.varint()?)?;
        out.push(c.take(len)?.to_vec());
    }
    if !c.is_empty() {
        bail!("trailing column bytes");
    }
    Ok(out)
}
fn path_length(pos: u32, shape: &Shape) -> Result<(Path, u32)> {
    if shape.blocks.first().map(|b| b.0) != Some(0) {
        bail!("nonzero initial shape offset");
    }
    let mut path = Vec::new();
    let mut total = 0u32;
    for (i, &(off, len)) in shape.blocks.iter().enumerate() {
        if len == 0 {
            bail!("zero block");
        }
        total = total.checked_add(len).context("aligned length overflow")?;
        let end = pos
            .checked_add(off)
            .and_then(|p| p.checked_add(len))
            .context("coordinate overflow")?;
        if let Some(&(next, _)) = shape.blocks.get(i + 1) {
            let acceptor = pos.checked_add(next).context("coordinate overflow")?;
            if acceptor <= end {
                bail!("invalid junction order");
            }
            path.extend([end, acceptor]);
        }
    }
    Ok((path, total))
}
fn shape_from_path(
    start: u32,
    length: u32,
    path: &[u32],
) -> Result<smallvec::SmallVec<[(u32, u32); 4]>> {
    if !path.len().is_multiple_of(2) || length == 0 {
        bail!("invalid path or length");
    }
    let mut blocks = smallvec::SmallVec::new();
    let mut here = start;
    let mut remaining = length;
    for pair in path.chunks_exact(2) {
        let n = pair[0].checked_sub(here).context("junction before start")?;
        if n == 0 || pair[1] <= pair[0] {
            bail!("invalid junction");
        }
        remaining = remaining
            .checked_sub(n)
            .context("path exceeds aligned length")?;
        blocks.push((here - start, n));
        here = pair[1];
    }
    if remaining == 0 {
        bail!("empty terminal block");
    }
    here.checked_add(remaining)
        .context("terminal coordinate overflow")?;
    blocks.push((here - start, remaining));
    Ok(blocks)
}
fn grouped(mols: &[MolRec], shapes: &[Shape]) -> Result<Vec<Groups>> {
    mols.iter()
        .map(|m| {
            let mut groups: Groups = BTreeMap::new();
            for chain in &m.chains {
                if chain.reps.len() != 1 || chain.weight == 0 {
                    bail!("requires positive full-fidelity singleton chains");
                }
                let (pos, id) = chain.reps[0];
                let (path, length) =
                    path_length(pos, shapes.get(id as usize).context("missing shape")?)?;
                groups
                    .entry(path)
                    .or_default()
                    .push((pos, length, chain.weight));
            }
            for entries in groups.values_mut() {
                entries.sort_unstable();
            }
            Ok(groups)
        })
        .collect()
}
// Absolute observed junctions, not transcript annotations. Length differences are
// relative to a stored mode; no library length is assumed by the decoder.
fn encode_paths(mols: &[MolRec], shapes: &[Shape], bin: u32, constrained: bool) -> Result<Columns> {
    let groups = grouped(mols, shapes)?;
    let mut paths = BTreeMap::new();
    let mut lengths = BTreeMap::new();
    for g in &groups {
        for (path, entries) in g {
            paths.insert(path.clone(), 0u64);
            for &(_, len, _) in entries {
                *lengths.entry(len).or_insert(0usize) += 1;
            }
        }
    }
    for (i, id) in paths.values_mut().enumerate() {
        *id = i as u64;
    }
    let mode = lengths
        .into_iter()
        .max_by_key(|&(len, n)| (n, std::cmp::Reverse(len)))
        .map_or(1, |(l, _)| l);
    let mut c = vec![Vec::new(); 8];
    put(&mut c, 0, u64::from(mode));
    put(&mut c, 0, paths.len() as u64);
    for path in paths.keys() {
        put(&mut c, 0, (path.len() / 2) as u64);
        let mut last = bin;
        for &p in path {
            put(
                &mut c,
                0,
                u64::from(p.checked_sub(last).context("path precedes chunk")?),
            );
            last = p;
        }
    }
    for (m, g) in mols.iter().zip(&groups) {
        put(&mut c, 1, g.len() as u64);
        let mut last_id = 0;
        for (path, entries) in g {
            let id = paths[path];
            put(&mut c, 2, id - last_id);
            last_id = id;
            put(&mut c, 3, entries.len() as u64);
            let extra: u64 = entries.iter().map(|e| u64::from(e.2 - 1)).sum();
            if constrained {
                put(&mut c, 6, extra);
            }
            let mut last = m.anchor();
            for (i, &(pos, len, w)) in entries.iter().enumerate() {
                put(
                    &mut c,
                    4,
                    u64::from(pos.checked_sub(last).context("unsorted start")?),
                );
                last = pos;
                evidence_io::archive::put_svarint(&mut c[5], i64::from(len) - i64::from(mode));
                if !constrained || (extra > 0 && i + 1 < entries.len()) {
                    put(&mut c, 7, u64::from(w - 1));
                }
            }
        }
    }
    Ok(c)
}
fn decode_paths(
    cols: &[Vec<u8>],
    anchors: &[u32],
    bin: u32,
    shapes: &ShapeMap,
    constrained: bool,
) -> Result<Vec<Vec<MolChain>>> {
    if cols.len() != 8 {
        bail!("wrong path column count");
    }
    let mut c: Vec<_> = cols.iter().map(|x| Cursor::new(x)).collect();
    let mode = u32v(&mut c[0])?;
    let np = count(&mut c[0])?;
    let mut paths = Vec::new();
    for _ in 0..np {
        let nj = count(&mut c[0])?;
        let mut path = Vec::new();
        let mut last = bin;
        for _ in 0..nj.checked_mul(2).context("path length overflow")? {
            last = last
                .checked_add(u32v(&mut c[0])?)
                .context("path overflow")?;
            path.push(last);
        }
        paths.push(path);
    }
    let mut out = Vec::new();
    let mut seen = 0usize;
    for &anchor in anchors {
        let ng = count(&mut c[1])?;
        let mut id = 0usize;
        let mut chains = Vec::new();
        for _ in 0..ng {
            id = id
                .checked_add(count(&mut c[2])?)
                .context("path id overflow")?;
            let path = paths.get(id).context("missing path")?;
            let n = count(&mut c[3])?;
            if n == 0 {
                bail!("empty path group");
            }
            seen = seen.checked_add(n).context("geometry count overflow")?;
            if seen > MAX_ITEMS {
                bail!("geometry budget exceeded");
            }
            let mut remaining = if constrained { c[6].varint()? } else { 0 };
            let has_extra = remaining > 0;
            let mut pos = anchor;
            for i in 0..n {
                pos = pos
                    .checked_add(u32v(&mut c[4])?)
                    .context("start overflow")?;
                let len = u32::try_from(
                    i64::from(mode)
                        .checked_add(c[5].svarint()?)
                        .context("length overflow")?,
                )?;
                let blocks = shape_from_path(pos, len, path)?;
                let shape = *shapes
                    .get(blocks.as_slice())
                    .context("reconstructed shape missing from source dictionary")?;
                let extra = if !constrained {
                    c[7].varint()?
                } else if !has_extra {
                    0
                } else if i + 1 == n {
                    remaining
                } else {
                    let e = c[7].varint()?;
                    remaining = remaining
                        .checked_sub(e)
                        .context("multiplicity exceeds total")?;
                    e
                };
                let weight = u32::try_from(extra.checked_add(1).context("weight overflow")?)?;
                chains.push(MolChain {
                    weight,
                    reps: smallvec::smallvec![(pos, shape)],
                });
            }
        }
        chains.sort_unstable_by_key(|c| (c.reps[0], c.reps.last().copied(), c.weight));
        out.push(chains);
    }
    if c.iter().any(|c| !c.is_empty()) {
        bail!("trailing path values");
    }
    Ok(out)
}
// Share individual geometries, not whole molecule templates. Inline one-off entries.
fn encode_shared(mols: &[MolRec], bin: u32, deltas: bool, anchored: bool) -> Result<Columns> {
    let mut freq = BTreeMap::new();
    for m in mols {
        for ch in &m.chains {
            if ch.reps.len() != 1 || ch.weight == 0 {
                bail!("not full fidelity");
            }
            *freq.entry(ch.reps[0]).or_insert(0usize) += 1;
        }
    }
    let ids: BTreeMap<_, _> = freq
        .iter()
        .filter(|(_, n)| **n > 1)
        .enumerate()
        .map(|(i, (&g, _))| (g, i as u64 + 1))
        .collect();
    let mut c = vec![Vec::new(); 6];
    put(&mut c, 0, ids.len() as u64);
    let mut last = bin;
    for &(p, s) in ids.keys() {
        put(&mut c, 0, u64::from(p - last));
        put(&mut c, 0, u64::from(s));
        last = p;
    }
    let positions: Vec<_> = ids.keys().map(|g| g.0).collect();
    let mut last_id = 0i64;
    for m in mols {
        if anchored {
            last_id = positions.partition_point(|&p| p < m.anchor()) as i64;
        }
        put(&mut c, 1, m.chains.len() as u64);
        for ch in &m.chains {
            let g = ch.reps[0];
            let id = ids.get(&g).copied().unwrap_or(0);
            let token = if anchored && id > 0 {
                let delta = (id as i64)
                    .checked_sub(last_id)
                    .context("shared rank underflow")?;
                if delta <= 0 {
                    bail!("anchored references require distinct canonical geometries");
                }
                last_id = id as i64;
                delta as u64
            } else if deltas && id > 0 {
                let delta = i64::try_from(id)? - last_id;
                last_id = id as i64;
                ((delta << 1) ^ (delta >> 63)) as u64 + 1
            } else {
                id
            };
            put(&mut c, 2, token);
            if id == 0 {
                put(&mut c, 3, u64::from(g.0 - m.anchor()));
                put(&mut c, 4, u64::from(g.1));
            }
            put(&mut c, 5, u64::from(ch.weight - 1));
        }
    }
    Ok(c)
}
fn decode_shared(
    cols: &[Vec<u8>],
    anchors: &[u32],
    bin: u32,
    nshapes: usize,
    deltas: bool,
    anchored: bool,
) -> Result<Vec<Vec<MolChain>>> {
    if cols.len() != 6 {
        bail!("wrong shared column count");
    }
    let mut c: Vec<_> = cols.iter().map(|x| Cursor::new(x)).collect();
    let n = count(&mut c[0])?;
    let mut dict = Vec::new();
    let mut p = bin;
    for _ in 0..n {
        p = p
            .checked_add(u32v(&mut c[0])?)
            .context("dictionary position overflow")?;
        let s = u32v(&mut c[0])?;
        if s as usize >= nshapes {
            bail!("invalid shape");
        }
        dict.push((p, s));
    }
    let mut out = Vec::new();
    let mut seen = 0usize;
    let mut last_id = 0i64;
    for &anchor in anchors {
        if anchored {
            last_id = dict.partition_point(|g| g.0 < anchor) as i64;
        }
        let n = count(&mut c[1])?;
        seen = seen.checked_add(n).context("count overflow")?;
        if seen > MAX_ITEMS {
            bail!("geometry budget exceeded");
        }
        let mut chains = Vec::new();
        for _ in 0..n {
            let token = c[2].varint()?;
            let id = if anchored && token > 0 {
                last_id = last_id
                    .checked_add(i64::try_from(token)?)
                    .context("shared rank overflow")?;
                usize::try_from(last_id)?
            } else if deltas && token > 0 {
                let z = token - 1;
                let delta = ((z >> 1) as i64) ^ -((z & 1) as i64);
                last_id = last_id
                    .checked_add(delta)
                    .context("geometry reference overflow")?;
                let id = usize::try_from(last_id)?;
                if id == 0 {
                    bail!("invalid zero shared reference");
                }
                id
            } else {
                usize::try_from(token)?
            };
            let g = if id == 0 {
                (
                    anchor
                        .checked_add(u32v(&mut c[3])?)
                        .context("inline position overflow")?,
                    u32v(&mut c[4])?,
                )
            } else {
                *dict.get(id - 1).context("invalid geometry reference")?
            };
            if g.1 as usize >= nshapes {
                bail!("invalid inline shape");
            }
            let weight = u32v(&mut c[5])?.checked_add(1).context("weight overflow")?;
            chains.push(MolChain {
                weight,
                reps: smallvec::smallvec![g],
            });
        }
        out.push(chains);
    }
    if c.iter().any(|c| !c.is_empty()) {
        bail!("trailing shared values");
    }
    Ok(out)
}
fn skeleton(raw: &[u8]) -> Result<Vec<u8>> {
    let mut c = unpack(raw, 10)?;
    for i in [3, 4, 5] {
        c[i].clear();
    }
    Ok(pack(&c))
}
// Per-column entropy ablation using the existing Rust backend. Charge the mode,
// count and full local table; retain varints when their final zstd frame is smaller.
fn entropy_encode(raw: &[u8], level: i32) -> Result<Vec<u8>> {
    let mut c = Cursor::new(raw);
    let mut values = Vec::new();
    while !c.is_empty() {
        values.push(c.varint()?);
    }
    let mut counts = [0; rans::NSYM];
    rans::count(&values, &mut counts);
    let table = rans::Table::from_counts(&counts)?;
    let mut coded = vec![1];
    put_varint(&mut coded, values.len() as u64);
    table.serialize(&mut coded);
    rans::encode(&values, &table, &mut coded);
    let plain = [&[0][..], raw].concat();
    Ok(
        if compress(&coded, level)?.len() < compress(&plain, level)?.len() {
            coded
        } else {
            plain
        },
    )
}
fn entropy_decode(raw: &[u8]) -> Result<Vec<u8>> {
    let mut c = Cursor::new(raw);
    match c.byte()? {
        0 => Ok(raw[1..].to_vec()),
        1 => {
            let n = count(&mut c)?;
            let table = rans::Table::deserialize(&mut c)?;
            let values = rans::decode_limited(&raw[c.position()..], &table, n)?;
            if values.len() != n {
                bail!("entropy count mismatch");
            }
            let mut out = Vec::new();
            for value in values {
                put_varint(&mut out, value);
            }
            Ok(out)
        }
        _ => bail!("unknown entropy column codec"),
    }
}
// Directly decode unchanged identity/MM columns; never rebuild old rANS geometry.
type Skeleton = (Vec<MolRec>, Vec<u32>, Vec<usize>);
fn decode_skeleton(raw: &[u8], info: &ChunkInfo, d: &Dicts) -> Result<Skeleton> {
    let s = unpack(raw, 10)?;
    let mut a = Cursor::new(&s[0]);
    let mut l = Cursor::new(&s[2]);
    let ids = rans::decode_limited(&s[1], &d.rans_tables[0], info.n_mols as usize)?;
    if ids.len() != info.n_mols as usize {
        bail!("identity count mismatch");
    }
    let mut probe = Cursor::new(&s[2]);
    let mut nm = 0usize;
    for _ in &ids {
        probe.byte()?;
        count(&mut probe)?;
        nm = nm
            .checked_add(count(&mut probe)?)
            .context("MM count overflow")?;
    }
    if nm > MAX_ITEMS {
        bail!("MM budget exceeded");
    }
    let mp = rans::decode_limited(&s[6], &d.rans_tables[3], nm)?;
    let mw = rans::decode_limited(&s[9], &d.rans_tables[4], nm)?;
    if mp.len() != nm || mw.len() != nm {
        bail!("MM cardinality mismatch");
    }
    let mut ms = Cursor::new(&s[7]);
    let mut mt = Cursor::new(&s[8]);
    let (mut anchor, mut next, mut mi) = (info.bin_start, info.class_base, 0usize);
    let (mut out, mut anchors, mut sizes) = (Vec::new(), Vec::new(), Vec::new());
    for token in ids {
        anchor = anchor
            .checked_add(u32v(&mut a)?)
            .context("anchor overflow")?;
        let class = if token == 0 {
            let id = next;
            next = next.checked_add(1).context("class overflow")?;
            id
        } else {
            next.checked_sub(u32::try_from(token)?)
                .context("class underflow")?
        };
        let strand = l.byte()?;
        if strand > 1 {
            bail!("invalid strand");
        }
        sizes.push(count(&mut l)?);
        let n = count(&mut l)?;
        let mut mms = smallvec::SmallVec::new();
        for _ in 0..n {
            let pos = anchor
                .checked_add(u32::try_from(*mp.get(mi).context("missing MM")?)?)
                .context("MM position overflow")?;
            mms.push((pos, u32v(&mut ms)?, u32v(&mut mt)?, u32::try_from(mw[mi])?));
            mi += 1;
        }
        out.push(MolRec {
            cell: *d
                .cell_of_class
                .get(class as usize)
                .context("invalid class")?,
            umi_class: class,
            chrom: info.chrom,
            strand_rev: strand == 1,
            chains: smallvec::SmallVec::new(),
            mms,
        });
        anchors.push(anchor);
    }
    if [a, l, ms, mt].iter().any(|c| !c.is_empty())
        || mi != nm
        || [3, 4, 5].iter().any(|&i| !s[i].is_empty())
    {
        bail!("trailing skeleton values");
    }
    Ok((out, anchors, sizes))
}
#[derive(Clone, Copy)]
struct Variant {
    name: &'static str,
    mode: u8,
    split: u8,
    constrained: bool,
}
const VARIANTS: &[Variant] = &[
    Variant {
        name: "original-joint",
        mode: 0,
        split: 0,
        constrained: false,
    },
    Variant {
        name: "original-split",
        mode: 0,
        split: 2,
        constrained: false,
    },
    Variant {
        name: "path-joint",
        mode: 1,
        split: 0,
        constrained: false,
    },
    Variant {
        name: "path-grouped",
        mode: 1,
        split: 1,
        constrained: false,
    },
    Variant {
        name: "path-split",
        mode: 1,
        split: 2,
        constrained: false,
    },
    Variant {
        name: "path-entropy-split",
        mode: 1,
        split: 3,
        constrained: false,
    },
    Variant {
        name: "path-counts-grouped",
        mode: 1,
        split: 1,
        constrained: true,
    },
    Variant {
        name: "shared-grouped",
        mode: 2,
        split: 1,
        constrained: false,
    },
    Variant {
        name: "shared-deltas-grouped",
        mode: 3,
        split: 1,
        constrained: false,
    },
    Variant {
        name: "shared-anchor-grouped",
        mode: 4,
        split: 1,
        constrained: false,
    },
];
fn encode(
    raw: &[u8],
    mols: &[MolRec],
    d: &Dicts,
    info: &ChunkInfo,
    v: Variant,
    level: i32,
) -> Result<Columns> {
    if v.mode == 0 {
        return if v.split == 0 {
            Ok(vec![raw.to_vec()])
        } else {
            unpack(raw, 10)
        };
    }
    let sk = skeleton(raw)?;
    let mut c = if v.mode == 1 {
        encode_paths(mols, &d.shapes, info.bin_start, v.constrained)?
    } else {
        encode_shared(mols, info.bin_start, v.mode == 3, v.mode == 4)?
    };
    if v.split == 3 {
        c = c
            .iter()
            .map(|col| entropy_encode(col, level))
            .collect::<Result<_>>()?;
    }
    Ok(match v.split {
        0 => vec![pack(&[vec![sk], c].concat())],
        1 => vec![sk, pack(&c[..1]), pack(&c[1..4]), pack(&c[4..])],
        _ => [vec![sk], c].concat(),
    })
}
fn decode(
    parts: &[Vec<u8>],
    d: &Dicts,
    info: &ChunkInfo,
    shapes: &ShapeMap,
    v: Variant,
) -> Result<Vec<MolRec>> {
    if v.mode == 0 {
        let raw = if v.split == 0 {
            parts[0].clone()
        } else {
            pack(parts)
        };
        return archivecmd::decode_chunk(&raw, info, Some(&d.cell_of_class), &d.rans_tables);
    }
    let n = if v.mode == 1 { 8 } else { 6 };
    let mut all = match v.split {
        0 => unpack(&parts[0], n + 1)?,
        1 => [
            vec![parts[0].clone()],
            unpack(&parts[1], 1)?,
            unpack(&parts[2], 3)?,
            unpack(&parts[3], n - 4)?,
        ]
        .concat(),
        _ => parts.to_vec(),
    };
    if v.split == 3 {
        for col in &mut all[1..] {
            *col = entropy_decode(col)?;
        }
    }
    let (mut mols, anchors, sizes) = decode_skeleton(&all[0], info, d)?;
    let chains = if v.mode == 1 {
        decode_paths(&all[1..], &anchors, info.bin_start, shapes, v.constrained)?
    } else {
        decode_shared(
            &all[1..],
            &anchors,
            info.bin_start,
            d.shapes.len(),
            v.mode == 3,
            v.mode == 4,
        )?
    };
    for ((m, c), n) in mols.iter_mut().zip(chains).zip(sizes) {
        if c.len() != n {
            bail!("chain count mismatch");
        }
        m.chains = c.into();
    }
    Ok(mols)
}
fn median(v: &[f64]) -> f64 {
    let mut x = v.to_vec();
    x.sort_by(f64::total_cmp);
    x[x.len() / 2]
}
pub fn run(args: Args) -> Result<()> {
    if args.out.exists() {
        bail!("output directory already exists");
    }
    let mut reader = SectionReader::open(&args.archive)?;
    reader.verify_all_payloads()?;
    let d = archivecmd::read_dicts(&mut reader)?;
    let infos = archivecmd::read_chunk_index(&mut reader)?;
    let shapes: ShapeMap = d
        .shapes
        .iter()
        .enumerate()
        .map(|(i, s)| (s.blocks.clone(), i as u32))
        .collect();
    if shapes.len() != d.shapes.len() {
        bail!("ambiguous duplicate shape IDs");
    }
    let mut chunks = Vec::new();
    let mut originals = Vec::new();
    for (i, info) in infos.iter().enumerate() {
        let raw = reader.read(&format!("c{i}"))?;
        originals.push(archivecmd::decode_chunk(
            &raw,
            info,
            Some(&d.cell_of_class),
            &d.rans_tables,
        )?);
        chunks.push(raw);
    }
    let mut geo = BTreeMap::new();
    let mut lens = BTreeMap::new();
    let mut entries = 0usize;
    let mut groups = 0usize;
    for mols in &originals {
        groups += grouped(mols, &d.shapes)?
            .iter()
            .map(|g| g.len())
            .sum::<usize>();
        for m in mols {
            for ch in &m.chains {
                entries += 1;
                let (p, s) = ch.reps[0];
                *geo.entry((m.chrom, p, s)).or_insert(0usize) += 1;
                *lens
                    .entry(path_length(p, &d.shapes[s as usize])?.1)
                    .or_insert(0usize) += 1;
            }
        }
    }
    std::fs::create_dir(&args.out)?;
    let names: Vec<_> = reader.names().map(str::to_owned).collect();
    let mut results = Vec::new();
    for &v in VARIANTS {
        let path = args.out.join(format!("{}.gexp", v.name));
        let mut w = SectionWriter::create_new(&path, args.level)?;
        let manifest = serde_json::json!({"schema":"gravlax.path-experiment.v1","variant":v.name,"source_root":reader.content_commitment().context("source root required")?.to_hex(),"production_reader":false,"retained_dictionaries":true});
        w.section("experiment", &serde_json::to_vec(&manifest)?)?;
        for name in &names {
            if name
                .strip_prefix('c')
                .is_some_and(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()))
            {
                continue;
            }
            let (comp, len) = reader.read_compressed(name)?;
            w.section_precompressed(name, len as u64, &comp)?;
        }
        let mut cached = Vec::new();
        let mut encode_seconds = 0.0;
        let mut payload = 0usize;
        let mut group_bytes = Vec::<usize>::new();
        for (i, ((raw, mols), info)) in chunks.iter().zip(&originals).zip(&infos).enumerate() {
            let t = Instant::now();
            let parts = encode(raw, mols, &d, info, v, args.level)?;
            let compressed: Vec<_> = parts
                .iter()
                .map(|p| compress(p, args.level))
                .collect::<Result<_>>()?;
            encode_seconds += t.elapsed().as_secs_f64();
            if decode(&parts, &d, info, &shapes, v)? != *mols {
                bail!("{} reconstruction differs at chunk {i}", v.name);
            }
            if group_bytes.len() < parts.len() {
                group_bytes.resize(parts.len(), 0);
            }
            for (j, (p, c)) in parts.iter().zip(&compressed).enumerate() {
                w.section_precompressed(&format!("c{i}.g{j}"), p.len() as u64, c)?;
                payload += c.len();
                group_bytes[j] += c.len();
                if v.name == "path-split" && j > 0 {
                    let mut cur = Cursor::new(p);
                    let mut fixed = Vec::new();
                    while !cur.is_empty() {
                        fixed.extend(cur.varint()?.to_le_bytes());
                    }
                    std::fs::write(args.out.join(format!("typed-c{i}-s{}.u64le", j - 1)), fixed)?;
                }
            }
            cached.push(
                parts
                    .iter()
                    .zip(compressed)
                    .map(|(p, c)| (p.len(), c))
                    .collect::<Vec<_>>(),
            );
        }
        w.finish()?;
        let mut check = SectionReader::open(&path)?;
        check.verify_all_payloads()?;
        for (i, info) in infos.iter().enumerate() {
            let parts = (0..cached[i].len())
                .map(|j| check.read(&format!("c{i}.g{j}")))
                .collect::<Result<Vec<_>>>()?;
            if decode(&parts, &d, info, &shapes, v)? != originals[i] {
                bail!("persisted reconstruction differs");
            }
        }
        let mut times = Vec::new();
        for _ in 0..args.repeats {
            let t = Instant::now();
            for (i, info) in infos.iter().enumerate() {
                let parts = cached[i]
                    .iter()
                    .map(|(n, c)| decompress(c, *n))
                    .collect::<Result<Vec<_>>>()?;
                std::hint::black_box(decode(&parts, &d, info, &shapes, v)?);
            }
            times.push(t.elapsed().as_secs_f64());
        }
        let before = check.bytes_read();
        let t = Instant::now();
        for (i, info) in infos.iter().enumerate() {
            let identity_raw = if v.mode == 0 && v.split == 2 {
                let mut cols = vec![Vec::new(); 10];
                for (j, col) in cols.iter_mut().enumerate().take(3) {
                    *col = check.read(&format!("c{i}.g{j}"))?;
                }
                pack(&cols)
            } else {
                let raw = check.read(&format!("c{i}.g0"))?;
                if v.mode != 0 && v.split == 0 {
                    unpack(&raw, 9)?[0].clone()
                } else {
                    raw
                }
            };
            let identities =
                archivecmd::decode_chunk_identities(&identity_raw, info, &d.rans_tables)?;
            for (id, m) in identities.iter().zip(&originals[i]) {
                if (id.class, id.anchor, id.reverse) != (m.umi_class, m.anchor(), m.strand_rev) {
                    bail!("identity projection differs");
                }
            }
        }
        let identity_seconds = t.elapsed().as_secs_f64();
        let identity_read_bytes = check.bytes_read() - before;
        results.push(serde_json::json!({"name":v.name,"file_bytes":std::fs::metadata(&path)?.len(),"chunk_payload_bytes":payload,"group_payload_bytes":group_bytes,"encode_seconds":encode_seconds,"decode_median_seconds":median(&times),"decode_seconds":times,"identity_scan_read_bytes":identity_read_bytes,"identity_scan_seconds":identity_seconds,"exact_rows_verified":true}));
        eprintln!(
            "{}: {} bytes, {:.4}s decode",
            v.name,
            std::fs::metadata(&path)?.len(),
            median(&times)
        );
    }
    let report = serde_json::json!({"schema":"gravlax.path-benchmark.v1","source":args.archive,"source_bytes":std::fs::metadata(&args.archive)?.len(),"source_root":reader.content_commitment().context("root missing")?.to_hex(),"chunks":infos.len(),"level":args.level,"repeats":args.repeats,"geometry_entries":entries,"path_groups":groups,"distinct_absolute_geometries":geo.len(),"repeated_geometry_entries":geo.values().filter(|&&n|n>1).sum::<usize>(),"aligned_lengths":lens,"variants":results,"limitations":"Experimental containers, not production archives. All original dictionaries/indexes/auxiliary sections retained. Exact complete MolRec equality checked before and after persistence. Timings are warm in-memory decompression plus native row reconstruction, not production query or replay; dictionary setup excluded equally. Identity I/O counts authenticated section reads excluding already-loaded directory/dictionaries. No new biological evidence is captured."});
    let report = serde_json::to_string_pretty(&report)?;
    std::fs::write(args.out.join("report.json"), &report)?;
    println!("{report}");
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (Vec<MolRec>, Vec<Shape>) {
        let shapes = vec![
            Shape {
                blocks: vec![(0, 30), (130, 61)],
            },
            Shape {
                blocks: vec![(0, 20), (120, 70)],
            },
            Shape {
                blocks: vec![(0, 91)],
            },
            Shape {
                blocks: vec![(0, 10), (100, 20), (200, 62)],
            },
        ];
        let make = |cell, chains| MolRec {
            cell,
            umi_class: cell,
            chrom: 0,
            strand_rev: cell == 1,
            chains,
            mms: smallvec::smallvec![(100, 2, 0, 3)],
        };
        let ch = |pos, shape, weight| MolChain {
            weight,
            reps: smallvec::smallvec![(pos, shape)],
        };
        (
            vec![
                make(
                    0,
                    smallvec::smallvec![ch(100, 0, 4), ch(110, 1, 2), ch(115, 2, 1), ch(150, 3, 1)],
                ),
                make(1, smallvec::smallvec![ch(100, 0, 1), ch(110, 1, 1)]),
                make(2, smallvec::smallvec![]),
            ],
            shapes,
        )
    }
    #[test]
    fn exact_paths_counts_and_shared_geometries() {
        let (mols, shapes) = fixture();
        let map = shapes
            .iter()
            .enumerate()
            .map(|(i, s)| (s.blocks.clone(), i as u32))
            .collect();
        let anchors = mols.iter().map(MolRec::anchor).collect::<Vec<_>>();
        let expected = mols.iter().map(|m| m.chains.to_vec()).collect::<Vec<_>>();
        for constrained in [false, true] {
            let c = encode_paths(&mols, &shapes, 0, constrained).unwrap();
            assert_eq!(
                decode_paths(&c, &anchors, 0, &map, constrained).unwrap(),
                expected
            );
            for i in 0..c.len() {
                if !c[i].is_empty() {
                    let mut bad = c.clone();
                    bad[i].pop();
                    assert!(decode_paths(&bad, &anchors, 0, &map, constrained).is_err());
                }
                let mut bad = c.clone();
                bad[i].push(0);
                assert!(decode_paths(&bad, &anchors, 0, &map, constrained).is_err());
            }
        }
        for (deltas, anchored) in [(false, false), (true, false), (false, true)] {
            let c = encode_shared(&mols, 0, deltas, anchored).unwrap();
            assert_eq!(
                decode_shared(&c, &anchors, 0, shapes.len(), deltas, anchored).unwrap(),
                expected
            );
            for i in 0..c.len() {
                let mut bad = c.clone();
                bad[i].push(0);
                assert!(decode_shared(&bad, &anchors, 0, shapes.len(), deltas, anchored).is_err());
            }
        }
    }
    #[test]
    fn geometry_edge_cases() {
        assert_eq!(
            shape_from_path(100, 91, &[130, 230]).unwrap().as_slice(),
            &[(0, 30), (130, 61)]
        );
        for (s, l, p) in [
            (100, 0, vec![]),
            (100, 30, vec![130, 230]),
            (100, 91, vec![90, 230]),
            (100, 91, vec![130, 120]),
            (u32::MAX, 1, vec![]),
        ] {
            assert!(shape_from_path(s, l, &p).is_err());
        }
        let (mols, shapes) = fixture();
        let mut bad = mols;
        bad[0].chains[0].reps.push((100, 0));
        assert!(encode_paths(&bad, &shapes, 0, false).is_err());
    }
    #[test]
    fn variable_length_multijunction_randomized_roundtrips() {
        // Deliberately not a 91-base library: 24 deterministic seeds, mixed paths,
        // broad lengths, same-start geometries, MM-only rows and maximum counts.
        for seed in 1..=24u64 {
            let mut state = seed;
            let mut next = || {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                (state >> 32) as u32
            };
            let mut shapes = Vec::new();
            let mut ids = ShapeMap::default();
            let mut mols = Vec::new();
            for i in 0..80u32 {
                let anchor = 1000 + i * 10000;
                let mut chains = Vec::new();
                for j in 0..(next() % 12) {
                    let start = anchor + (next() % 5) * 4;
                    let nblocks = 1 + next() % 7;
                    let mut blocks = Vec::new();
                    let mut off = 0;
                    for _ in 0..nblocks {
                        let len = 1 + next() % 200;
                        blocks.push((off, len));
                        off += len + 1 + next() % 1000;
                    }
                    let id = *ids.entry(blocks.clone()).or_insert_with(|| {
                        let id = shapes.len() as u32;
                        shapes.push(Shape { blocks });
                        id
                    });
                    chains.push(MolChain {
                        weight: if j == 0 { u32::MAX } else { 1 + next() % 1000 },
                        reps: smallvec::smallvec![(start, id)],
                    });
                }
                chains.sort_unstable_by_key(|c| (c.reps[0], c.reps.last().copied(), c.weight));
                let m = MolRec {
                    cell: i % 3,
                    umi_class: i,
                    chrom: 2,
                    strand_rev: i % 2 == 1,
                    chains: chains.into(),
                    mms: smallvec::smallvec![(anchor, 0, 0, 1)],
                };
                mols.push(m.clone());
                if i % 3 == 0 {
                    mols.push(m);
                }
            }
            let anchors = mols.iter().map(MolRec::anchor).collect::<Vec<_>>();
            let expected = mols.iter().map(|m| m.chains.to_vec()).collect::<Vec<_>>();
            for flag in [false, true] {
                let c = encode_paths(&mols, &shapes, 1000, flag).unwrap();
                assert_eq!(c, encode_paths(&mols, &shapes, 1000, flag).unwrap());
                assert_eq!(
                    decode_paths(&c, &anchors, 1000, &ids, flag).unwrap(),
                    expected
                );
                let c = encode_shared(&mols, 1000, flag, false).unwrap();
                assert_eq!(
                    decode_shared(&c, &anchors, 1000, shapes.len(), flag, false).unwrap(),
                    expected
                );
            }
            let c = encode_shared(&mols, 1000, false, true).unwrap();
            assert_eq!(
                decode_shared(&c, &anchors, 1000, shapes.len(), false, true).unwrap(),
                expected
            );
        }
    }
    #[test]
    fn malformed_envelopes_and_count_conservation() {
        let (mols, shapes) = fixture();
        let anchors = mols.iter().map(MolRec::anchor).collect::<Vec<_>>();
        let ids = shapes
            .iter()
            .enumerate()
            .map(|(i, s)| (s.blocks.clone(), i as u32))
            .collect();
        let mut cols = encode_paths(&mols, &shapes, 0, true).unwrap();
        // First group's total extra is zero: a forged total larger than u32 must
        // fail when reconstructed into positive per-geometry counts.
        cols[6] = vec![255, 255, 255, 255, 255, 255, 255, 255, 255, 1];
        assert!(decode_paths(&cols, &anchors, 0, &ids, true).is_err());
        assert!(unpack(&[255; 10], 1).is_err());
        assert!(unpack(&[4, 1], 1).is_err());
        let mut c = Cursor::new(&[255, 255, 255, 127]);
        assert!(count(&mut c).is_err());
    }
    #[test]
    fn adaptive_entropy_columns_roundtrip() {
        for values in [
            vec![],
            vec![0; 20000],
            (0..20000).map(|i| i % 7).collect(),
            vec![0, 1, u64::MAX, u64::MAX - 1],
        ] {
            let mut raw = Vec::new();
            for v in values {
                put_varint(&mut raw, v);
            }
            let coded = entropy_encode(&raw, 1).unwrap();
            assert_eq!(entropy_decode(&coded).unwrap(), raw);
            if coded[0] == 1 {
                let mut bad = coded.clone();
                bad.pop();
                assert!(entropy_decode(&bad).is_err());
            }
        }
        assert!(entropy_decode(&[]).is_err());
        assert!(entropy_decode(&[2]).is_err());
    }
}
