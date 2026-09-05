//! Experimental compact-base + sparse corrections capsule, not a production archive codec.
use crate::{
    archivecmd::{decode_chunk, read_chunk_index, read_dicts},
    rows::{MolChain, MolRec},
};
use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use evidence_io::{
    archive::{put_svarint, put_varint, Shape},
    format::{Cursor, SectionReader, SectionWriter},
};
use std::{collections::BTreeMap, path::PathBuf};

#[derive(ClapArgs)]
pub struct Args {
    compact: PathBuf,
    fidelity: PathBuf,
    /// New directory for experimental correction capsules and their JSON report.
    #[arg(long)]
    out: PathBuf,
    #[arg(long, default_value_t = 19)]
    level: i32,
}
type Geometry = (u32, u32);
fn key(pos: u32, shape: &Shape) -> Vec<u64> {
    shape
        .blocks
        .windows(2)
        .flat_map(|p| {
            [
                u64::from(pos) + u64::from(p[0].0) + u64::from(p[0].1),
                u64::from(pos) + u64::from(p[1].0),
            ]
        })
        .collect()
}
fn retained(chain: &MolChain) -> Vec<Geometry> {
    let mut out = chain.reps.to_vec();
    out.sort_unstable();
    out.dedup();
    out
}
fn pack(columns: &[Vec<u8>]) -> Vec<u8> {
    let mut out = b"SPARSE01".to_vec();
    for c in columns {
        put_varint(&mut out, c.len() as u64);
        out.extend(c);
    }
    out
}
fn entropy_pack(raw: &[u8], level: i32) -> Result<Vec<u8>> {
    use evidence_io::{format::compress, rans};
    let mut c = Cursor::new(raw);
    if c.take(8)? != b"SPARSE01" {
        bail!("invalid sparse input");
    }
    let mut out = b"SPARSE02".to_vec();
    for _ in 0..7 {
        let n = usize::try_from(c.varint()?)?;
        let column = c.take(n)?;
        let mut v = Cursor::new(column);
        let mut values = Vec::new();
        while !v.is_empty() {
            values.push(v.varint()?);
        }
        let mut counts = [0; rans::NSYM];
        rans::count(&values, &mut counts);
        let table = rans::Table::from_counts(&counts)?;
        let mut encoded = Vec::new();
        table.serialize(&mut encoded);
        rans::encode(&values, &table, &mut encoded);
        let use_rans = compress(&encoded, level)?.len() < compress(column, level)?.len();
        out.push(u8::from(use_rans));
        put_varint(&mut out, n as u64);
        let bytes = if use_rans { encoded.as_slice() } else { column };
        put_varint(&mut out, bytes.len() as u64);
        out.extend(bytes);
    }
    if !c.is_empty() {
        bail!("trailing sparse input");
    }
    Ok(out)
}

fn entropy_unpack(raw: &[u8]) -> Result<Vec<u8>> {
    use evidence_io::rans;
    let mut c = Cursor::new(raw);
    if c.take(8)? != b"SPARSE02" {
        bail!("invalid entropy envelope");
    }
    let mut columns = Vec::new();
    for _ in 0..7 {
        let tag = c.take(1)?[0];
        let n = usize::try_from(c.varint()?)?;
        if n > 256 * 1024 * 1024 {
            bail!("correction column exceeds experiment budget");
        }
        let len = usize::try_from(c.varint()?)?;
        let bytes = c.take(len)?;
        let column = match tag {
            0 => bytes.to_vec(),
            1 => {
                let mut tc = Cursor::new(bytes);
                let table = rans::Table::deserialize(&mut tc)?;
                let values = rans::decode_limited(&bytes[tc.position()..], &table, n)?;
                let mut out = Vec::new();
                for value in values {
                    put_varint(&mut out, value);
                }
                out
            }
            _ => bail!("unsupported correction column codec"),
        };
        if column.len() != n {
            bail!("correction column length mismatch");
        }
        columns.push(column);
    }
    if !c.is_empty() {
        bail!("trailing entropy envelope");
    }
    Ok(pack(&columns))
}

fn encode(
    base: &[MolRec],
    full: &[MolRec],
    shapes: &[Shape],
    delta_shapes: bool,
) -> Result<(Vec<u8>, usize, usize)> {
    if base.len() != full.len() {
        bail!("record counts differ");
    }
    let mut cols = vec![Vec::new(); 7];
    let (mut ordinal, mut last_route, mut changed, mut omitted) = (0usize, 0usize, 0usize, 0usize);
    for (b, f) in base.iter().zip(full) {
        if (
            b.cell,
            b.umi_class,
            b.chrom,
            b.strand_rev,
            b.anchor(),
            &b.mms,
        ) != (
            f.cell,
            f.umi_class,
            f.chrom,
            f.strand_rev,
            f.anchor(),
            &f.mms,
        ) {
            bail!("base/fidelity record identity differs");
        }
        let mut groups: BTreeMap<Vec<u64>, BTreeMap<Geometry, u32>> = BTreeMap::new();
        for ch in &f.chains {
            if ch.reps.len() != 1 || ch.weight == 0 {
                bail!("full input must have positive singleton geometry weights");
            }
            let (pos, shape) = ch.reps[0];
            if groups
                .entry(key(
                    pos,
                    shapes.get(shape as usize).context("missing shape")?,
                ))
                .or_default()
                .insert((pos, shape), ch.weight)
                .is_some()
            {
                bail!("duplicate fidelity geometry");
            }
        }
        for ch in &b.chains {
            ordinal += 1;
            let reps = retained(ch);
            let first = *reps.first().context("empty base chain")?;
            let group = groups
                .remove(&key(first.0, &shapes[first.1 as usize]))
                .context("base chain absent from fidelity")?;
            if group.values().map(|&w| u64::from(w)).sum::<u64>() != u64::from(ch.weight) {
                bail!("chain multiplicities differ");
            }
            for rep in &reps {
                if !group.contains_key(rep) {
                    bail!("base representative absent from full geometry");
                }
            }
            let extras: Vec<_> = group
                .keys()
                .filter(|g| reps.binary_search(g).is_err())
                .copied()
                .collect();
            // Counts default to one; the final geometry receives the conserved remainder.
            let patches: Vec<_> = group
                .values()
                .take(group.len().saturating_sub(1))
                .enumerate()
                .filter(|(_, w)| **w != 1)
                .map(|(i, &w)| (i, w))
                .collect();
            if extras.is_empty() && patches.is_empty() {
                continue;
            }
            changed += 1;
            omitted += extras.len();
            put_varint(&mut cols[0], (ordinal - last_route) as u64);
            last_route = ordinal;
            put_varint(&mut cols[1], extras.len() as u64);
            put_varint(&mut cols[2], patches.len() as u64);
            let (mut pos, mut shape) = (first.0, first.1);
            for (p, s) in extras {
                put_varint(
                    &mut cols[3],
                    u64::from(
                        p.checked_sub(pos)
                            .context("extra precedes retained anchor")?,
                    ),
                );
                pos = p;
                if delta_shapes {
                    put_svarint(&mut cols[4], i64::from(s) - i64::from(shape));
                } else {
                    put_varint(&mut cols[4], s as u64);
                }
                shape = s;
            }
            let mut last = 0;
            for (i, w) in patches {
                put_varint(&mut cols[5], (i - last) as u64);
                last = i;
                put_varint(&mut cols[6], w as u64);
            }
        }
        if !groups.is_empty() {
            bail!("full geometry has chains absent from compact base");
        }
    }
    Ok((pack(&cols), changed, omitted))
}

fn decode(
    base: &[MolRec],
    raw: &[u8],
    shapes: &[Shape],
    delta_shapes: bool,
) -> Result<Vec<MolRec>> {
    if raw.starts_with(b"SPARSE02") {
        return decode(base, &entropy_unpack(raw)?, shapes, delta_shapes);
    }
    let mut c = Cursor::new(raw);
    if c.take(8)? != b"SPARSE01" {
        bail!("unsupported corrections");
    }
    let mut cols = Vec::new();
    for _ in 0..7 {
        let n = usize::try_from(c.varint()?)?;
        cols.push(Cursor::new(c.take(n)?));
    }
    if !c.is_empty() {
        bail!("trailing correction columns");
    }
    let mut patches = BTreeMap::new();
    let total: usize = base.iter().map(|m| m.chains.len()).sum();
    let mut ordinal = 0usize;
    while !cols[0].is_empty() {
        let delta = usize::try_from(cols[0].varint()?)?;
        ordinal = ordinal.checked_add(delta).context("route overflow")?;
        if delta == 0 || ordinal > total {
            bail!("invalid correction route");
        }
        let extra = usize::try_from(cols[1].varint()?)?;
        let counts = usize::try_from(cols[2].varint()?)?;
        if extra > raw.len() || counts > raw.len() || extra + counts == 0 {
            bail!("invalid correction count");
        }
        let mut geos = Vec::new();
        for _ in 0..extra {
            geos.push((
                cols[3].varint()?,
                if delta_shapes {
                    cols[4].svarint()?
                } else {
                    i64::try_from(cols[4].varint()?)?
                },
            ));
        }
        let mut weights = Vec::new();
        let mut index = 0usize;
        for i in 0..counts {
            let delta = usize::try_from(cols[5].varint()?)?;
            if i > 0 && delta == 0 {
                bail!("duplicate count correction");
            }
            index = index.checked_add(delta).context("count index overflow")?;
            let weight = u32::try_from(cols[6].varint()?)?;
            if weight <= 1 {
                bail!("noncanonical count correction");
            }
            weights.push((index, weight));
        }
        patches.insert(ordinal, (geos, weights));
    }
    if cols.iter().any(|c| !c.is_empty()) {
        bail!("unused correction bytes");
    }
    let mut out = Vec::new();
    let mut ordinal = 0;
    for b in base {
        let mut m = b.clone();
        m.chains.clear();
        for chain in &b.chains {
            ordinal += 1;
            let mut reps = retained(chain);
            let first = *reps.first().context("empty base chain")?;
            let chain_key = key(first.0, &shapes[first.1 as usize]);
            let mut corrected = Vec::new();
            if let Some((geos, weights)) = patches.remove(&ordinal) {
                let (mut pos, mut shape) = (first.0, first.1);
                for (delta, s) in geos {
                    pos = pos
                        .checked_add(u32::try_from(delta)?)
                        .context("geometry position overflow")?;
                    shape = if delta_shapes {
                        u32::try_from(
                            i64::from(shape)
                                .checked_add(s)
                                .context("shape delta overflow")?,
                        )?
                    } else {
                        u32::try_from(s)?
                    };
                    if key(
                        pos,
                        shapes
                            .get(shape as usize)
                            .context("missing correction shape")?,
                    ) != chain_key
                    {
                        bail!("correction changes junction chain");
                    }
                    reps.push((pos, shape));
                }
                corrected = weights;
            }
            reps.sort_unstable();
            if reps.windows(2).any(|p| p[0] == p[1]) {
                bail!("duplicate correction geometry");
            }
            let mut weights = vec![1u32; reps.len()];
            for (i, w) in corrected {
                if i >= weights.len() - 1 {
                    bail!("invalid count correction index");
                }
                weights[i] = w;
            }
            let sum: u64 = weights[..weights.len() - 1]
                .iter()
                .map(|&w| u64::from(w))
                .sum();
            let last = u64::from(chain.weight)
                .checked_sub(sum)
                .context("multiplicities exceed chain total")?;
            if last == 0 {
                bail!("zero residual multiplicity");
            }
            *weights.last_mut().unwrap() = u32::try_from(last)?;
            for (rep, weight) in reps.into_iter().zip(weights) {
                m.chains.push(MolChain {
                    weight,
                    reps: smallvec::smallvec![rep],
                });
            }
        }
        m.chains
            .sort_unstable_by_key(|ch| (ch.reps[0], ch.reps.last().copied(), ch.weight));
        out.push(m);
    }
    if !patches.is_empty() {
        bail!("unapplied corrections");
    }
    Ok(out)
}

pub fn run(args: Args) -> Result<()> {
    let mut base = SectionReader::open(&args.compact)?;
    let mut full = SectionReader::open(&args.fidelity)?;
    base.verify_all_payloads()?;
    full.verify_all_payloads()?;
    let bd = read_dicts(&mut base)?;
    let fd = read_dicts(&mut full)?;
    if bd.shapes != fd.shapes
        || bd.cells != fd.cells
        || bd.patterns != fd.patterns
        || bd.cell_of_class != fd.cell_of_class
        || bd.chrom_names != fd.chrom_names
    {
        bail!("source dictionaries differ");
    }
    let bc = read_chunk_index(&mut base)?;
    let fc = read_chunk_index(&mut full)?;
    if bc.len() != fc.len() {
        bail!("chunk counts differ");
    }
    let base_index = if base.has(crate::accessindex::SECTION) {
        Some(crate::accessindex::Index::decode(
            &base.read(crate::accessindex::SECTION)?,
            bc.len(),
            bd.n_classes,
            bd.chrom_names.len(),
        )?)
    } else {
        None
    };
    let full_index = if full.has(crate::accessindex::SECTION) {
        Some(crate::accessindex::Index::decode(
            &full.read(crate::accessindex::SECTION)?,
            fc.len(),
            fd.n_classes,
            fd.chrom_names.len(),
        )?)
    } else {
        None
    };
    let supplement = match (&base_index, &full_index) {
        (Some(b), Some(f)) => f.supplement_to(b)?,
        (None, None) => None,
        _ => bail!("both inputs must have the same access-index capability"),
    };
    let manifest = serde_json::to_vec(
        &serde_json::json!({"schema":"gravlax.experimental-sparse-fidelity.v1","base_root":base.content_commitment().context("compact archive needs a root")?.to_hex(),"fidelity_root":full.content_commitment().context("fidelity archive needs a root")?.to_hex(),"rule":"one per retained geometry, last gets conserved remainder; explicit missing geometry and nondefault nonfinal counts","routing":"include every correction chunk in full-fidelity geometry queries; base geometry postings alone are insufficient","production_reader":false}),
    )?;
    std::fs::create_dir(&args.out)?;
    let base_bytes = std::fs::metadata(&args.compact)?.len();
    let mut results = Vec::new();
    for (delta, entropy, name) in [
        (false, false, "absolute-shapes"),
        (true, false, "delta-shapes"),
        (true, true, "delta-shapes-entropy"),
    ] {
        let path = args.out.join(format!("{name}.capsule"));
        let mut writer = SectionWriter::create_new(&path, args.level)?;
        let mut declaration: serde_json::Value = serde_json::from_slice(&manifest)?;
        declaration["delta_shapes"] = serde_json::json!(delta);
        declaration["entropy_columns"] = serde_json::json!(entropy);
        declaration["routing"] = serde_json::json!(if base_index.is_some() {
            "union base routes with index.access.supplement; apply corrections to selected chunks before predicates"
        } else {
            "full scan required for fidelity geometry queries"
        });
        writer.section("manifest", &serde_json::to_vec(&declaration)?)?;
        if let Some(raw) = &supplement {
            writer.section("index.access.supplement", raw)?;
        }
        let (mut changed, mut omitted, mut correction_chunks) = (0usize, 0usize, 0usize);
        let (mut base_posts, mut full_posts) = (0usize, 0usize);
        for (i, (b, f)) in bc.iter().zip(&fc).enumerate() {
            if (b.chrom, b.bin_start, b.class_base, b.n_mols)
                != (f.chrom, f.bin_start, f.class_base, f.n_mols)
            {
                bail!("chunk identity mismatch");
            }
            let bm = decode_chunk(
                &base.read(&format!("c{i}"))?,
                b,
                Some(&bd.cell_of_class),
                &bd.rans_tables,
            )?;
            let fm = decode_chunk(
                &full.read(&format!("c{i}"))?,
                f,
                Some(&fd.cell_of_class),
                &fd.rans_tables,
            )?;
            if let Some(index) = &base_index {
                base_posts +=
                    index.verify_chunk(i as u32, b.class_base, &bm, &bd.shapes, &bd.patterns)?;
            }
            if let Some(index) = &full_index {
                full_posts +=
                    index.verify_chunk(i as u32, f.class_base, &fm, &fd.shapes, &fd.patterns)?;
            }
            let (raw, n, extra) = encode(&bm, &fm, &bd.shapes, delta)?;
            let raw = if entropy {
                entropy_pack(&raw, args.level)?
            } else {
                raw
            };
            if decode(&bm, &raw, &bd.shapes, delta)? != fm {
                bail!("sparse reconstruction differs from full fidelity");
            }
            if n > 0 {
                writer.section(&format!("correction.c{i}"), &raw)?;
                correction_chunks += 1;
            }
            changed += n;
            omitted += extra;
        }
        if base_index
            .as_ref()
            .is_some_and(|i| i.posting_count() != base_posts)
            || full_index
                .as_ref()
                .is_some_and(|i| i.posting_count() != full_posts)
        {
            bail!("access index has extra routes");
        }
        writer.finish()?;
        let capsule_bytes = std::fs::metadata(&path)?.len();
        let mut check = SectionReader::open(&path)?;
        check.verify_all_payloads()?;
        if serde_json::from_slice::<serde_json::Value>(&check.read("manifest")?)? != declaration {
            bail!("capsule manifest differs");
        }
        for (i, (b, f)) in bc.iter().zip(&fc).enumerate() {
            let bm = decode_chunk(
                &base.read(&format!("c{i}"))?,
                b,
                Some(&bd.cell_of_class),
                &bd.rans_tables,
            )?;
            let fm = decode_chunk(
                &full.read(&format!("c{i}"))?,
                f,
                Some(&fd.cell_of_class),
                &fd.rans_tables,
            )?;
            let name = format!("correction.c{i}");
            let raw = if check.has(&name) {
                check.read(&name)?
            } else {
                pack(&vec![Vec::new(); 7])
            };
            if decode(&bm, &raw, &bd.shapes, delta)? != fm {
                bail!("published correction reconstruction mismatch");
            }
        }
        results.push(serde_json::json!({"delta_shapes":delta,"entropy_columns":entropy,"capsule_bytes":capsule_bytes,"base_plus_capsule_bytes":base_bytes+capsule_bytes,"premium_percent":100.0*capsule_bytes as f64/base_bytes as f64,"changed_chains":changed,"omitted_geometries":omitted,"correction_chunks":correction_chunks,"reconstruction_verified":true}));
    }
    let report = serde_json::to_string_pretty(
        &serde_json::json!({"schema":"gravlax.sparse-fidelity-benchmark.v1","base_bytes":base_bytes,"full_fidelity_bytes":std::fs::metadata(&args.fidelity)?.len(),"level":args.level,"variants":results,"supplemental_routing_required":supplement.is_some(),"base_access_index":base_index.is_some(),"limitation":"Experimental capsules, no production query adapter. Route by union of base and included supplemental postings; without an index use full scans."}),
    )?;
    std::fs::write(args.out.join("report.json"), &report)?;
    println!("{report}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    #[test]
    fn entropy_columns_roundtrip_and_reject_truncation() {
        let mut columns = vec![Vec::new(); 7];
        let mut state = 1u64;
        for _ in 0..2000 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            put_varint(&mut columns[0], state >> 54);
        }
        let raw = pack(&columns);
        let encoded = entropy_pack(&raw, 1).unwrap();
        assert_eq!(entropy_unpack(&encoded).unwrap(), raw);
        for end in [0, 7, 8, encoded.len() - 1] {
            assert!(entropy_unpack(&encoded[..end]).is_err());
        }
        let mut trailing = encoded;
        trailing.push(0);
        assert!(entropy_unpack(&trailing).is_err());
    }
    #[test]
    fn omitted_middle_and_multiplicities_reconstruct() {
        let shapes = vec![Shape {
            blocks: vec![(0, 20)],
        }];
        let base = MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: false,
            chains: smallvec![MolChain {
                weight: 7,
                reps: smallvec![(100, 0), (140, 0)]
            }],
            mms: smallvec![],
        };
        let mut full = base.clone();
        full.chains = smallvec![
            MolChain {
                weight: 2,
                reps: smallvec![(100, 0)]
            },
            MolChain {
                weight: 4,
                reps: smallvec![(120, 0)]
            },
            MolChain {
                weight: 1,
                reps: smallvec![(140, 0)]
            }
        ];
        for delta in [false, true] {
            let (raw, n, extra) = encode(
                std::slice::from_ref(&base),
                std::slice::from_ref(&full),
                &shapes,
                delta,
            )
            .unwrap();
            assert_eq!((n, extra), (1, 1));
            assert_eq!(
                decode(std::slice::from_ref(&base), &raw, &shapes, delta).unwrap(),
                vec![full.clone()]
            );
            for end in 0..raw.len() {
                assert!(decode(std::slice::from_ref(&base), &raw[..end], &shapes, delta).is_err());
            }
        }
    }

    #[test]
    fn counts_defaults_sparse_routes_and_duplicate_extremes() {
        let shapes = vec![Shape {
            blocks: vec![(0, 20)],
        }];
        for a in 1..=4 {
            for b in 1..=4 {
                for c in 1..=4 {
                    let unchanged = MolRec {
                        cell: 0,
                        umi_class: 0,
                        chrom: 0,
                        strand_rev: false,
                        chains: smallvec![MolChain {
                            weight: 5,
                            reps: smallvec![(10, 0), (10, 0)]
                        }],
                        mms: smallvec![],
                    };
                    let base = MolRec {
                        cell: 0,
                        umi_class: 0,
                        chrom: 0,
                        strand_rev: false,
                        chains: smallvec![MolChain {
                            weight: a + b + c,
                            reps: smallvec![(5000, 0), (5040, 0)]
                        }],
                        mms: smallvec![],
                    };
                    let mut full = base.clone();
                    full.chains = [(5000, a), (5020, b), (5040, c)]
                        .into_iter()
                        .map(|(p, w)| MolChain {
                            weight: w,
                            reps: smallvec![(p, 0)],
                        })
                        .collect();
                    let mut first = unchanged.clone();
                    first.chains[0].reps.pop();
                    let bases = [unchanged, base];
                    let expected = [first, full];
                    for delta in [false, true] {
                        let (raw, n, extra) = encode(&bases, &expected, &shapes, delta).unwrap();
                        assert_eq!((n, extra), (1, 1));
                        assert_eq!(decode(&bases, &raw, &shapes, delta).unwrap(), expected);
                    }
                }
            }
        }
    }
}
