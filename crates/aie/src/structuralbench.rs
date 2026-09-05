//! Read-only, reversible structural compression experiments. Never writes an archive.
use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use evidence_io::{
    archive::{put_svarint, put_varint},
    format::{compress, decompress, Cursor, SectionReader},
};
use std::{collections::BTreeMap, path::PathBuf, time::Instant};

#[derive(ClapArgs)]
pub struct Args {
    archive: PathBuf,
    #[arg(long, default_value_t = 19)]
    level: i32,
    #[arg(long, default_value_t = 7, value_parser = clap::value_parser!(u32).range(1..=100))]
    repeats: u32,
}

fn split(raw: &[u8], n: usize) -> Result<Vec<Vec<u8>>> {
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
fn join(streams: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for stream in streams {
        put_varint(&mut out, stream.len() as u64);
        out.extend_from_slice(stream);
    }
    out
}
fn values(raw: &[u8]) -> Result<Vec<u64>> {
    let mut c = Cursor::new(raw);
    let mut out = Vec::new();
    while !c.is_empty() {
        out.push(c.varint()?);
    }
    Ok(out)
}
fn varints(vals: &[u64]) -> Vec<u8> {
    let mut out = Vec::new();
    for &v in vals {
        put_varint(&mut out, v);
    }
    out
}

// Birth positions are a separate bitmap. Only repeats pay for an identity.
// Mode 2 retains backward distance; mode 3 differences successive repeated absolute IDs.
fn class_encode(vals: &[u64], base: u32, mode: u8) -> Result<Vec<u8>> {
    if mode == 1 {
        return Ok(varints(vals));
    }
    let mut bitmap = vec![0u8; vals.len().div_ceil(8)];
    let mut rest = Vec::new();
    let (mut next, mut last) = (u64::from(base), 0i64);
    for (i, &v) in vals.iter().enumerate() {
        if v == 0 {
            next += 1;
        } else {
            bitmap[i / 8] |= 1 << (i % 8);
            let id = i64::try_from(
                next.checked_sub(v)
                    .context("invalid class back-reference")?,
            )?;
            if mode == 2 {
                put_varint(&mut rest, v);
            } else {
                put_svarint(&mut rest, id - last);
                last = id;
            }
        }
    }
    bitmap.extend(rest);
    Ok(bitmap)
}
fn class_decode(raw: &[u8], n: usize, base: u32, mode: u8) -> Result<Vec<u64>> {
    if mode == 1 {
        let out = values(raw)?;
        if out.len() != n {
            bail!("class count mismatch");
        }
        return Ok(out);
    }
    let bytes = n.div_ceil(8);
    let bitmap = raw.get(..bytes).context("truncated class bitmap")?;
    let mut c = Cursor::new(&raw[bytes..]);
    let (mut next, mut last) = (u64::from(base), 0i64);
    let mut out = Vec::new();
    for i in 0..n {
        if bitmap[i / 8] & (1 << (i % 8)) == 0 {
            out.push(0);
            next += 1;
        } else {
            let v = if mode == 2 {
                c.varint()?
            } else {
                last = last
                    .checked_add(c.svarint()?)
                    .context("class delta overflow")?;
                next.checked_sub(u64::try_from(last)?)
                    .context("invalid repeated class")?
            };
            if v == 0 || v > next {
                bail!("invalid class back-reference");
            }
            out.push(v);
        }
    }
    if !c.is_empty() {
        bail!("trailing class bytes");
    }
    Ok(out)
}

// Share a frequency-ranked local shape vocabulary across unique and MM streams.
fn shape_encode(a: &[u8], b: &[u8]) -> Result<(Vec<u8>, Vec<u8>)> {
    let (a, b) = (values(a)?, values(b)?);
    let mut counts = BTreeMap::new();
    for &v in a.iter().chain(&b) {
        *counts.entry(v).or_insert(0u64) += 1;
    }
    let mut ranked: Vec<_> = counts.into_iter().collect();
    ranked.sort_unstable_by_key(|&(v, n)| (std::cmp::Reverse(n), v));
    let ids: BTreeMap<_, _> = ranked
        .iter()
        .enumerate()
        .map(|(i, &(v, _))| (v, i as u64))
        .collect();
    let mut first = Vec::new();
    put_varint(&mut first, ranked.len() as u64);
    for &(v, _) in &ranked {
        put_varint(&mut first, v);
    }
    for v in a {
        put_varint(&mut first, ids[&v]);
    }
    Ok((
        first,
        varints(&b.iter().map(|v| ids[v]).collect::<Vec<_>>()),
    ))
}
fn shape_decode(a: &[u8], b: &[u8]) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut c = Cursor::new(a);
    let n = usize::try_from(c.varint()?)?;
    if n > a.len() {
        bail!("invalid local shape dictionary count");
    }
    let mut dict = Vec::new();
    for _ in 0..n {
        dict.push(c.varint()?);
    }
    let expand = |raw: &[u8]| -> Result<Vec<u8>> {
        let mut out = Vec::new();
        for v in values(raw)? {
            put_varint(
                &mut out,
                *dict
                    .get(usize::try_from(v)?)
                    .context("invalid local shape id")?,
            );
        }
        Ok(out)
    };
    Ok((expand(&a[c.position()..])?, expand(b)?))
}

// Separate strand, chain counts and MM counts, instead of interleaving triples.
fn layout_encode(raw: &[u8]) -> Result<Vec<u8>> {
    let mut c = Cursor::new(raw);
    let mut cols = vec![Vec::new(); 3];
    while !c.is_empty() {
        cols[0].extend_from_slice(c.take(1)?);
        for col in &mut cols[1..] {
            put_varint(col, c.varint()?);
        }
    }
    Ok(join(&cols))
}
fn layout_decode(raw: &[u8]) -> Result<Vec<u8>> {
    let cols = split(raw, 3)?;
    let (a, b) = (values(&cols[1])?, values(&cols[2])?);
    if a.len() != cols[0].len() || b.len() != a.len() {
        bail!("layout count mismatch");
    }
    let mut out = Vec::new();
    for (i, &strand) in cols[0].iter().enumerate() {
        out.push(strand);
        put_varint(&mut out, a[i]);
        put_varint(&mut out, b[i]);
    }
    Ok(out)
}

#[derive(Clone, Copy)]
struct Variant {
    name: &'static str,
    class: u8,
    shapes: bool,
    layout: bool,
}
const VARIANTS: &[Variant] = &[
    Variant {
        name: "position-deltas-local-rans",
        class: 0,
        shapes: false,
        layout: true,
    },
    Variant {
        name: "position-deltas-chain-flanks-local-rans",
        class: 0,
        shapes: false,
        layout: true,
    },
    Variant {
        name: "chain-flanks",
        class: 0,
        shapes: false,
        layout: false,
    },
    Variant {
        name: "position-deltas-plus-chain-flanks",
        class: 0,
        shapes: false,
        layout: true,
    },
    Variant {
        name: "representative-position-deltas",
        class: 0,
        shapes: false,
        layout: false,
    },
    Variant {
        name: "position-deltas-plus-layout",
        class: 0,
        shapes: false,
        layout: true,
    },
    Variant {
        name: "unchanged",
        class: 0,
        shapes: false,
        layout: false,
    },
    Variant {
        name: "class-varints",
        class: 1,
        shapes: false,
        layout: false,
    },
    Variant {
        name: "class-bitmap-backrefs",
        class: 2,
        shapes: false,
        layout: false,
    },
    Variant {
        name: "class-bitmap-repeat-deltas",
        class: 3,
        shapes: false,
        layout: false,
    },
    Variant {
        name: "local-shape-ids",
        class: 0,
        shapes: true,
        layout: false,
    },
    Variant {
        name: "layout-columns",
        class: 0,
        shapes: false,
        layout: true,
    },
    Variant {
        name: "combined",
        class: 3,
        shapes: true,
        layout: true,
    },
    Variant {
        name: "class-plus-layout",
        class: 3,
        shapes: false,
        layout: true,
    },
    Variant {
        name: "geometry-templates",
        class: 0,
        shapes: false,
        layout: false,
    },
    Variant {
        name: "class-plus-geometry-templates",
        class: 3,
        shapes: false,
        layout: false,
    },
];

fn representative_counts(s: &[Vec<u8>], tables: &[evidence_io::rans::Table]) -> Result<Vec<usize>> {
    let mut layout = Cursor::new(&s[2]);
    let mut counts = Vec::new();
    let mut total = 0usize;
    while !layout.is_empty() {
        layout.take(1)?;
        let count = usize::try_from(layout.varint()?)?;
        layout.varint()?;
        total = total.checked_add(count).context("chain count overflow")?;
        counts.push(count);
    }
    let weights = evidence_io::rans::decode_limited(&s[3], &tables[1], total)?;
    if weights.len() != total {
        bail!("weight count mismatch");
    }
    let mut wi = 0;
    let mut reps = Vec::new();
    for count in counts {
        reps.push(
            weights[wi..wi + count]
                .iter()
                .map(|w| 1 + (w & 1) as usize)
                .sum(),
        );
        wi += count;
    }
    Ok(reps)
}

fn position_deltas(
    s: &mut [Vec<u8>],
    tables: &[evidence_io::rans::Table],
    inverse: bool,
) -> Result<()> {
    let counts = representative_counts(s, tables)?;
    let mut out = Vec::new();
    let input = if inverse {
        Vec::new()
    } else {
        evidence_io::rans::decode(&s[4], &tables[2])?
    };
    let mut input = input.into_iter();
    let mut c = Cursor::new(&s[4]);
    let mut restored = Vec::new();
    for reps in counts {
        let mut last = 0i64;
        for _ in 1..reps {
            if inverse {
                last = last
                    .checked_add(c.svarint()?)
                    .context("position delta overflow")?;
                restored.push(u64::try_from(last)?);
            } else {
                let pos = i64::try_from(input.next().context("position underrun")?)?;
                put_svarint(&mut out, pos - last);
                last = pos;
            }
        }
    }
    if inverse {
        if !c.is_empty() {
            bail!("trailing position deltas");
        }
        evidence_io::rans::encode(&restored, &tables[2], &mut out);
    } else if input.next().is_some() {
        bail!("trailing positions");
    }
    s[4] = out;
    Ok(())
}

fn chain_flanks(
    s: &mut [Vec<u8>],
    tables: &[evidence_io::rans::Table],
    shapes: &[evidence_io::archive::Shape],
    inverse: bool,
) -> Result<()> {
    let counts = representative_counts(s, tables)?;
    let mut positions = evidence_io::rans::decode(&s[4], &tables[2])?.into_iter();
    let shape_ids: BTreeMap<_, _> = shapes
        .iter()
        .enumerate()
        .map(|(i, s)| (s.blocks.clone(), i as u64))
        .collect();
    let mut c = Cursor::new(&s[5]);
    let mut out = Vec::new();
    for count in counts {
        let mut previous: Vec<(u64, u64)> = Vec::new();
        let mut chains = BTreeMap::new();
        for i in 0..count {
            let pos = if i == 0 {
                0
            } else {
                positions.next().context("position underrun")?
            };
            let shape = if inverse {
                let back = usize::try_from(c.varint()?)?;
                if back == 0 {
                    c.varint()?
                } else {
                    let &(old_pos, old_shape) = previous
                        .get(i.checked_sub(back).context("invalid chain reference")?)
                        .context("missing chain reference")?;
                    let old = shapes
                        .get(usize::try_from(old_shape)?)
                        .context("missing shape")?;
                    if old.blocks.first().map(|b| b.0) != Some(0) {
                        bail!("unsupported nonzero initial block offset");
                    }
                    let delta = c.svarint()?;
                    let mut blocks = old.blocks.clone();
                    let displacement = i64::try_from(old_pos)? - i64::try_from(pos)?;
                    if blocks.len() > 1 {
                        blocks[0].1 = u32::try_from(i64::from(blocks[0].1) + displacement)?;
                    }
                    for block in &mut blocks[1..] {
                        block.0 = u32::try_from(i64::from(block.0) + displacement)?;
                    }
                    let last = blocks.last_mut().context("empty shape")?;
                    last.1 = u32::try_from(i64::from(last.1) + delta)?;
                    *shape_ids
                        .get(&blocks)
                        .context("reconstructed geometry missing from shape dictionary")?
                }
            } else {
                c.varint()?
            };
            let definition = shapes
                .get(usize::try_from(shape)?)
                .context("missing shape")?;
            if inverse {
                put_varint(&mut out, shape);
            } else {
                let mut chain = Vec::new();
                for pair in definition.blocks.windows(2) {
                    chain.push(pos + u64::from(pair[0].0) + u64::from(pair[0].1));
                    chain.push(pos + u64::from(pair[1].0));
                }
                let last = definition.blocks.last().context("empty shape")?;
                let end = i64::from(last.1);
                if definition.blocks[0].0 != 0 {
                    put_varint(&mut out, 0);
                    put_varint(&mut out, shape);
                } else {
                    if let Some(&(ordinal, old_end)) = chains.get(&chain) {
                        put_varint(&mut out, (i - ordinal) as u64);
                        put_svarint(&mut out, end - old_end);
                    } else {
                        put_varint(&mut out, 0);
                        put_varint(&mut out, shape);
                    }
                    chains.insert(chain, (i, end));
                }
            }
            previous.push((pos, shape));
        }
    }
    if !c.is_empty() || positions.next().is_some() {
        bail!("trailing chain geometry");
    }
    s[5] = out;
    Ok(())
}

fn geometry_encode(mols: &[crate::rows::MolRec]) -> Vec<u8> {
    let mut dictionary = BTreeMap::new();
    let mut templates = Vec::new();
    let mut ids = Vec::new();
    for m in mols {
        let anchor = m.anchor();
        let mut key = vec![m.strand_rev as u8];
        put_varint(&mut key, m.chains.len() as u64);
        put_varint(&mut key, m.mms.len() as u64);
        for ch in &m.chains {
            put_varint(&mut key, ch.reps.len() as u64);
            for &(pos, shape) in &ch.reps {
                put_varint(&mut key, (pos - anchor) as u64);
                put_varint(&mut key, shape as u64);
            }
        }
        for &(pos, shape, pattern, _) in &m.mms {
            for v in [pos - anchor, shape, pattern] {
                put_varint(&mut key, v as u64);
            }
        }
        let next = templates.len() as u64;
        let id = *dictionary.entry(key.clone()).or_insert_with(|| {
            templates.push(key);
            next
        });
        ids.push(id);
    }
    let mut out = Vec::new();
    put_varint(&mut out, templates.len() as u64);
    out.extend(join(&templates));
    out.extend(varints(&ids));
    out
}

fn geometry_decode(
    raw: &[u8],
    n: usize,
    tables: &[evidence_io::rans::Table],
    s: &mut [Vec<u8>],
) -> Result<()> {
    let mut c = Cursor::new(raw);
    let count = usize::try_from(c.varint()?)?;
    if count > raw.len() {
        bail!("invalid geometry template count");
    }
    let mut templates = Vec::new();
    for _ in 0..count {
        let len = usize::try_from(c.varint()?)?;
        templates.push(c.take(len)?);
    }
    for i in [2, 4, 5, 6, 7, 8] {
        s[i].clear();
    }
    let (mut rep_pos, mut mm_pos) = (Vec::new(), Vec::new());
    for _ in 0..n {
        let id = usize::try_from(c.varint()?)?;
        let raw = templates.get(id).context("invalid geometry template id")?;
        let mut t = Cursor::new(raw);
        s[2].extend_from_slice(t.take(1)?);
        let (chains, mms) = (t.varint()?, t.varint()?);
        if chains > raw.len() as u64 || mms > raw.len() as u64 {
            bail!("invalid geometry child count");
        }
        put_varint(&mut s[2], chains);
        put_varint(&mut s[2], mms);
        let mut first = true;
        for _ in 0..chains {
            let reps = t.varint()?;
            if !(1..=2).contains(&reps) {
                bail!("invalid template representative count");
            }
            for _ in 0..reps {
                let pos = t.varint()?;
                if !first {
                    rep_pos.push(pos);
                } else if pos != 0 {
                    bail!("first representative is not the anchor");
                }
                first = false;
                put_varint(&mut s[5], t.varint()?);
            }
        }
        for _ in 0..mms {
            mm_pos.push(t.varint()?);
            put_varint(&mut s[7], t.varint()?);
            put_varint(&mut s[8], t.varint()?);
        }
        if !t.is_empty() {
            bail!("trailing geometry template bytes");
        }
    }
    if !c.is_empty() {
        bail!("trailing geometry IDs");
    }
    evidence_io::rans::encode(&rep_pos, &tables[2], &mut s[4]);
    evidence_io::rans::encode(&mm_pos, &tables[3], &mut s[6]);
    Ok(())
}

fn transform(
    raw: &[u8],
    info: &crate::archivecmd::ChunkInfo,
    tables: &[evidence_io::rans::Table],
    shapes: &[evidence_io::archive::Shape],
    v: Variant,
    inverse: bool,
) -> Result<Vec<u8>> {
    let raw = if inverse && v.name != "unchanged" {
        if raw.first() != Some(&1) {
            bail!("unsupported experiment envelope");
        }
        &raw[1..]
    } else {
        raw
    };
    let mut s = split(raw, 10)?;
    if inverse && v.name.contains("local-rans") {
        let mut c = Cursor::new(&s[4]);
        let table = evidence_io::rans::Table::deserialize(&mut c)?;
        let decoded = evidence_io::rans::decode(&s[4][c.position()..], &table)?;
        s[4] = varints(&decoded);
    }
    if v.name.contains("chain-flanks") && !inverse {
        chain_flanks(&mut s, tables, shapes, false)?;
    }
    if v.name.contains("position-deltas") && !inverse {
        position_deltas(&mut s, tables, false)?;
    }
    if !inverse && v.name.contains("local-rans") {
        let decoded = values(&s[4])?;
        let mut counts = [0; evidence_io::rans::NSYM];
        evidence_io::rans::count(&decoded, &mut counts);
        let table = evidence_io::rans::Table::from_counts(&counts)?;
        s[4].clear();
        table.serialize(&mut s[4]);
        evidence_io::rans::encode(&decoded, &table, &mut s[4]);
    }
    if v.name.contains("geometry-templates") {
        if inverse {
            let geometry = s[2].clone();
            geometry_decode(&geometry, info.n_mols as usize, tables, &mut s)?;
        } else {
            let mols = crate::archivecmd::decode_chunk(raw, info, None, tables)?;
            s[2] = geometry_encode(&mols);
            for i in [4, 5, 6, 7, 8] {
                s[i].clear();
            }
        }
    }
    let table = &tables[0];
    if v.class != 0 {
        s[1] = if inverse {
            let vals = class_decode(&s[1], info.n_mols as usize, info.class_base, v.class)?;
            let mut out = Vec::new();
            evidence_io::rans::encode(&vals, table, &mut out);
            out
        } else {
            class_encode(
                &evidence_io::rans::decode_limited(&s[1], table, info.n_mols as usize)?,
                info.class_base,
                v.class,
            )?
        };
    }
    if v.shapes {
        let (a, b) = if inverse {
            shape_decode(&s[5], &s[7])?
        } else {
            shape_encode(&s[5], &s[7])?
        };
        s[5] = a;
        s[7] = b;
    }
    if v.layout {
        s[2] = if inverse {
            layout_decode(&s[2])?
        } else {
            layout_encode(&s[2])?
        };
    }
    if v.name.contains("position-deltas") && inverse {
        position_deltas(&mut s, tables, true)?;
    }
    if v.name.contains("chain-flanks") && inverse {
        chain_flanks(&mut s, tables, shapes, true)?;
    }
    let mut out = join(&s);
    if !inverse && v.name != "unchanged" {
        out.insert(0, 1);
    }
    Ok(out)
}

pub fn run(args: Args) -> Result<()> {
    let mut reader = SectionReader::open(&args.archive)?;
    reader.verify_all_payloads()?;
    let identity = reader.content_commitment().map(|v| format!("{v:?}"));
    let bytes = std::fs::metadata(&args.archive)?.len();
    let dicts = crate::archivecmd::read_dicts(&mut reader)?;
    let chunks = crate::archivecmd::read_chunk_index(&mut reader)?;
    let mut totals = vec![0u64; VARIANTS.len()];
    let mut encode_times = vec![vec![0f64; args.repeats as usize]; VARIANTS.len()];
    let mut decode_times = encode_times.clone();
    let mut original = 0u64;
    for (i, info) in chunks.iter().enumerate() {
        let name = format!("c{i}");
        original += reader
            .entries()
            .iter()
            .find(|e| e.0 == name)
            .context("missing chunk")?
            .3;
        let raw = reader.read(&name)?;
        let source = crate::archivecmd::decode_chunk(
            &raw,
            info,
            Some(&dicts.cell_of_class),
            &dicts.rans_tables,
        )?;
        for (vi, &variant) in VARIANTS.iter().enumerate() {
            for rep in 0..args.repeats as usize {
                let start = Instant::now();
                let encoded = transform(
                    &raw,
                    info,
                    &dicts.rans_tables,
                    &dicts.shapes,
                    variant,
                    false,
                )?;
                let compressed = compress(&encoded, args.level)?;
                encode_times[vi][rep] += start.elapsed().as_secs_f64();
                if rep == 0 {
                    totals[vi] += compressed.len() as u64;
                }
                let start = Instant::now();
                let inflated = decompress(&compressed, encoded.len())?;
                let restored = if variant.name == "unchanged" {
                    inflated
                } else {
                    transform(
                        &inflated,
                        info,
                        &dicts.rans_tables,
                        &dicts.shapes,
                        variant,
                        true,
                    )?
                };
                let decoded = crate::archivecmd::decode_chunk(
                    &restored,
                    info,
                    Some(&dicts.cell_of_class),
                    &dicts.rans_tables,
                )?;
                decode_times[vi][rep] += start.elapsed().as_secs_f64();
                if restored != raw || decoded != source {
                    bail!("{} failed exact restoration", variant.name);
                }
            }
        }
    }
    let mut shape_raw = Vec::new();
    for shape in &dicts.shapes {
        put_varint(&mut shape_raw, shape.blocks.len() as u64);
        for &(off, len) in &shape.blocks {
            put_varint(&mut shape_raw, off as u64);
            put_varint(&mut shape_raw, len as u64);
        }
    }
    let factored = crate::shapecodec::encode(&dicts.shapes)?;
    if crate::shapecodec::decode(&factored)? != dicts.shapes {
        bail!("shape factoring changed geometry");
    }
    let median = |v: &[f64]| {
        let mut v = v.to_vec();
        v.sort_by(f64::total_cmp);
        v[v.len() / 2]
    };
    let variants:Vec<_>=VARIANTS.iter().enumerate().map(|(i,v)|serde_json::json!({
        "name":v.name,"chunk_compressed_bytes":totals[i],"projected_archive_bytes":bytes-original+totals[i],
        "encode_median_seconds":median(&encode_times[i]),"restore_and_decode_median_seconds":median(&decode_times[i]),
        "encode_seconds":encode_times[i],"restore_and_decode_seconds":decode_times[i]
    })).collect();
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({
            "schema":"gravlax.structural-experiment.v1","archive":args.archive,"source_root":identity,
            "source_file_bytes":bytes,"source_chunk_compressed_bytes":original,"chunks":chunks.len(),"level":args.level,
            "repeats":args.repeats,"exact_restoration_verified":true,"variants":variants,
            "shapes":{"ordinary_frame_bytes":compress(&shape_raw,args.level)?.len(),"factored_frame_bytes":compress(&factored,args.level)?.len(),"additional_section_name_bytes":18},
            "limitations":"Read-only research transforms, not supported archive codecs. Projected file sizes retain existing dictionaries/tables/indexes and charge an experiment-envelope byte per changed chunk. Restore timing includes canonical stream reconstruction (including rANS re-encoding), not a native candidate reader or query/replay benchmark."
        }))?
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn flanks_and_position_deltas_preserve_nonoverlapping_and_spliced_reads() {
        use evidence_io::archive::Shape;
        let shapes = vec![
            Shape {
                blocks: vec![(0, 20)],
            },
            Shape {
                blocks: vec![(0, 25)],
            },
            Shape {
                blocks: vec![(0, 10), (100, 10)],
            },
            Shape {
                blocks: vec![(0, 7), (97, 15)],
            },
        ];
        let mut counts = [0; evidence_io::rans::NSYM];
        evidence_io::rans::count(&[2, 4, 25, 100, 103], &mut counts);
        let tables: Vec<_> = (0..4)
            .map(|_| evidence_io::rans::Table::from_counts(&counts).unwrap())
            .collect();
        let mut s = vec![Vec::new(); 10];
        s[2] = vec![0, 4, 0];
        evidence_io::rans::encode(&[2, 4, 2, 2], &tables[1], &mut s[3]);
        evidence_io::rans::encode(&[25, 100, 103], &tables[2], &mut s[4]);
        s[5] = varints(&[0, 1, 2, 3]);
        let expected = s.clone();
        chain_flanks(&mut s, &tables, &shapes, false).unwrap();
        position_deltas(&mut s, &tables, false).unwrap();
        position_deltas(&mut s, &tables, true).unwrap();
        chain_flanks(&mut s, &tables, &shapes, true).unwrap();
        assert_eq!(s, expected);
        s[5] = vec![1];
        assert!(chain_flanks(&mut s, &tables, &shapes, true).is_err());
    }
    #[test]
    fn geometry_templates_restore_unique_and_multimapper_columns() {
        use crate::rows::{MolChain, MolRec};
        use smallvec::smallvec;
        let molecule = MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: true,
            chains: smallvec![MolChain {
                weight: 7,
                reps: smallvec![(100, 1), (102, 1)]
            }],
            mms: smallvec![(103, 2, 0, 4)],
        };
        let raw = geometry_encode(&[molecule.clone(), molecule]);
        let mut counts = [0; evidence_io::rans::NSYM];
        evidence_io::rans::count(&[2, 3], &mut counts);
        let tables: Vec<_> = (0..4)
            .map(|_| evidence_io::rans::Table::from_counts(&counts).unwrap())
            .collect();
        let mut s = vec![Vec::new(); 10];
        geometry_decode(&raw, 2, &tables, &mut s).unwrap();
        assert_eq!(s[2], vec![1, 1, 1, 1, 1, 1]);
        assert_eq!(values(&s[5]).unwrap(), vec![1, 1, 1, 1]);
        assert_eq!(values(&s[7]).unwrap(), vec![2, 2]);
        assert_eq!(values(&s[8]).unwrap(), vec![0, 0]);
        assert_eq!(
            evidence_io::rans::decode(&s[4], &tables[2]).unwrap(),
            vec![2, 2]
        );
        assert_eq!(
            evidence_io::rans::decode(&s[6], &tables[3]).unwrap(),
            vec![3, 3]
        );
        for end in 0..raw.len() {
            assert!(geometry_decode(&raw[..end], 2, &tables, &mut s).is_err());
        }
    }
    #[test]
    fn structures_restore_exact_values() {
        for vals in [vec![], vec![0, 0, 1, 0, 3, 2, 0, 4], vec![2, 1, 0, 3]] {
            for mode in 1..=3 {
                assert_eq!(
                    class_decode(&class_encode(&vals, 2, mode).unwrap(), vals.len(), 2, mode)
                        .unwrap(),
                    vals
                );
            }
        }
        let (a, b) = (varints(&[900, 900, 7, 1, 900]), varints(&[7, 900]));
        let (x, y) = shape_encode(&a, &b).unwrap();
        assert_eq!(shape_decode(&x, &y).unwrap(), (a, b));
        let layout = vec![0, 1, 0, 1, 2, 1, 0, 0, 1];
        assert_eq!(
            layout_decode(&layout_encode(&layout).unwrap()).unwrap(),
            layout
        );
        assert!(class_decode(&[], 1, 0, 2).is_err());
        assert!(shape_decode(&[1], &[]).is_err());
        assert!(layout_decode(&[1, 0],).is_err());
    }
}
