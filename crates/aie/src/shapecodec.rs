//! Optional shape dictionary factoring. IDs and decoded block geometry are unchanged.
use anyhow::{bail, Context, Result};
use evidence_io::{
    archive::{put_varint, Shape},
    format::Cursor,
};
use std::collections::BTreeMap;

pub const SECTION: &str = "shapes.factored";
const MAGIC: &[u8] = b"SHPFACT1";
const MAX_BLOCKS: usize = 10_000_000;

pub fn encode(shapes: &[Shape]) -> Result<Vec<u8>> {
    let mut ids = BTreeMap::new();
    let mut skeletons = Vec::new();
    let mut records = Vec::new();
    for shape in shapes {
        let first = shape.blocks.first().context("empty shape")?;
        let mut skeleton = Vec::new();
        for (i, pair) in shape.blocks.windows(2).enumerate() {
            let end = pair[0].0.checked_add(pair[0].1).context("shape overflow")?;
            skeleton.push(
                pair[1]
                    .0
                    .checked_sub(end)
                    .context("overlapping shape blocks")?,
            );
            if i + 2 < shape.blocks.len() {
                skeleton.push(pair[1].1);
            }
        }
        let next = ids.len();
        let id = *ids.entry(skeleton.clone()).or_insert_with(|| {
            skeletons.push(skeleton);
            next
        });
        records.push((id, first.0, first.1, shape.blocks.last().unwrap().1));
    }
    let mut out = MAGIC.to_vec();
    put_varint(&mut out, skeletons.len() as u64);
    for skeleton in &skeletons {
        put_varint(&mut out, skeleton.len() as u64);
        for v in skeleton {
            put_varint(&mut out, *v as u64);
        }
    }
    put_varint(&mut out, records.len() as u64);
    for (id, off, first, last) in records {
        for v in [id as u64, off as u64, first as u64] {
            put_varint(&mut out, v);
        }
        if !skeletons[id].is_empty() {
            put_varint(&mut out, last as u64);
        }
    }
    Ok(out)
}

pub fn decode(raw: &[u8]) -> Result<Vec<Shape>> {
    if !raw.starts_with(MAGIC) {
        bail!("unsupported factored shape codec");
    }
    let mut c = Cursor::new(&raw[MAGIC.len()..]);
    let n = usize::try_from(c.varint()?)?;
    if n > raw.len() {
        bail!("invalid shape skeleton count");
    }
    let mut skeletons = Vec::new();
    for _ in 0..n {
        let len = usize::try_from(c.varint()?)?;
        if len > raw.len() || (len != 0 && len % 2 == 0) {
            bail!("invalid shape skeleton length");
        }
        let mut skeleton = Vec::new();
        for _ in 0..len {
            skeleton.push(u32::try_from(c.varint()?)?);
        }
        skeletons.push(skeleton);
    }
    let n = usize::try_from(c.varint()?)?;
    if n > raw.len() / 3 {
        bail!("invalid factored shape count");
    }
    let mut shapes = Vec::new();
    let mut total = 0usize;
    for _ in 0..n {
        let id = usize::try_from(c.varint()?)?;
        let skeleton = skeletons.get(id).context("invalid shape skeleton id")?;
        total = total
            .checked_add(1 + skeleton.len().div_ceil(2))
            .context("shape count overflow")?;
        if total > MAX_BLOCKS {
            bail!("factored shapes exceed decoded block budget");
        }
        let mut off = u32::try_from(c.varint()?)?;
        let mut len = u32::try_from(c.varint()?)?;
        if len == 0 {
            bail!("zero-length shape block");
        }
        off.checked_add(len).context("shape coordinate overflow")?;
        let mut blocks = vec![(off, len)];
        let last = if skeleton.is_empty() {
            0
        } else {
            u32::try_from(c.varint()?)?
        };
        for i in (0..skeleton.len()).step_by(2) {
            off = off
                .checked_add(len)
                .and_then(|v| v.checked_add(skeleton[i]))
                .context("shape coordinate overflow")?;
            len = skeleton.get(i + 1).copied().unwrap_or(last);
            if len == 0 {
                bail!("zero-length shape block");
            }
            off.checked_add(len).context("shape coordinate overflow")?;
            blocks.push((off, len));
        }
        shapes.push(Shape { blocks });
    }
    if !c.is_empty() {
        bail!("trailing factored shape bytes");
    }
    Ok(shapes)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_roundtrip_and_truncation() {
        let shapes = vec![
            Shape {
                blocks: vec![(0, 20)],
            },
            Shape {
                blocks: vec![(0, 20), (100, 30), (200, 40)],
            },
            Shape {
                blocks: vec![(3, 25), (108, 30), (208, 50)],
            },
        ];
        let raw = encode(&shapes).unwrap();
        let decoded = decode(&raw).unwrap();
        for (a, b) in shapes.iter().zip(decoded) {
            assert_eq!(a.blocks, b.blocks);
        }
        for n in 0..raw.len() {
            assert!(decode(&raw[..n]).is_err());
        }
        let mut trailing = raw;
        trailing.push(0);
        assert!(decode(&trailing).is_err());
        assert!(decode(b"SHPFACT1\x00\x01\x00\x00\x01").is_err());
    }
}
