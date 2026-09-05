//! Optional root-bound routes over retained evidence. Routes only prune chunks;
//! ordinary source records still determine every predicate and count.
use crate::rows::{Extracted, MolRec, PatAlt, SAME_SHAPE};
use anyhow::{bail, Context, Result};
use evidence_io::archive::Shape;
use evidence_io::{archive::put_varint, format::Cursor};
use std::collections::{BTreeMap, BTreeSet};

pub const SECTION: &str = "index.access";
pub const SCHEMA: &str = "gravlax.access-index.v1";
const MAGIC: &[u8] = b"ACCESS01";
const TILE: u32 = 65536;
const MAX_POSTS: usize = 10_000_000;
// key: kind (0 repeated class, 1 interval tile, 2 junction), placement mode,
// chromosome/class, start/tile, end. mode: unique=0, direct=1, all=2.
type Key = [u32; 5];

pub struct Index {
    bases: Vec<u32>,
    chroms: u32,
    posts: BTreeMap<Key, Vec<u32>>,
}

impl Index {
    pub fn build(x: &Extracted, counts: impl Iterator<Item = u32>) -> Result<Self> {
        let mut index = Self {
            bases: Vec::new(),
            chroms: x.chrom_names.len() as u32,
            posts: BTreeMap::new(),
        };
        let (mut offset, mut next, mut total) = (0usize, 0u32, 0usize);
        for (chunk, count) in counts.enumerate() {
            let chunk = u32::try_from(chunk)?;
            index.bases.push(next);
            let end = offset
                .checked_add(count as usize)
                .context("access-index record count overflow")?;
            for molecule in x
                .mols
                .get(offset..end)
                .context("access-index chunk exceeds source records")?
            {
                if molecule.umi_class == next {
                    next += 1;
                } else if molecule.umi_class < index.bases[chunk as usize] {
                    index.add([0, 0, molecule.umi_class, 0, 0], chunk, &mut total)?;
                }
                for chain in &molecule.chains {
                    for &(pos, shape) in &chain.reps {
                        for mode in 0..3 {
                            index.placement(
                                &x.shapes,
                                mode,
                                molecule.chrom,
                                pos,
                                shape,
                                chunk,
                                &mut total,
                            )?;
                        }
                    }
                }
                for &(pos, shape, pattern, _) in &molecule.mms {
                    index.placement(&x.shapes, 1, molecule.chrom, pos, shape, chunk, &mut total)?;
                    for alt in x
                        .patterns
                        .get(pattern as usize)
                        .context("access-index missing pattern")?
                    {
                        let pos = i64::from(pos)
                            .checked_add(alt.offset)
                            .and_then(|v| u32::try_from(v).ok())
                            .context("access-index alternative coordinate overflow")?;
                        index.placement(
                            &x.shapes,
                            2,
                            alt.chrom,
                            pos,
                            if alt.shape == SAME_SHAPE {
                                shape
                            } else {
                                alt.shape
                            },
                            chunk,
                            &mut total,
                        )?;
                    }
                }
            }
            offset = end;
        }
        if offset != x.mols.len() || next != x.n_classes {
            bail!("access-index class introduction or chunk count mismatch");
        }
        index.bases.push(next);
        Ok(index)
    }

    fn add(&mut self, key: Key, chunk: u32, total: &mut usize) -> Result<()> {
        if self.posts.get(&key).and_then(|p| p.last()) == Some(&chunk) {
            return Ok(());
        }
        if *total == MAX_POSTS {
            bail!("access index exceeds {MAX_POSTS} distinct postings; omit --access-index");
        }
        self.posts.entry(key).or_default().push(chunk);
        *total += 1;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn placement(
        &mut self,
        shapes: &[Shape],
        mode: u32,
        chrom: u32,
        pos: u32,
        shape: u32,
        chunk: u32,
        total: &mut usize,
    ) -> Result<()> {
        let shape = shapes
            .get(shape as usize)
            .context("access-index missing shape")?;
        let mut previous = None;
        for &(offset, len) in &shape.blocks {
            let start = pos
                .checked_add(offset)
                .context("access-index coordinate overflow")?;
            let end = start
                .checked_add(len)
                .context("access-index coordinate overflow")?;
            if len == 0 || previous.is_some_and(|p| p > start) {
                bail!("access-index malformed aligned blocks");
            }
            for tile in start / TILE..=(end - 1) / TILE {
                self.add([1, mode, chrom, tile, 0], chunk, total)?;
            }
            if let Some(donor) = previous.filter(|&p| p < start) {
                self.add([2, mode, chrom, donor, start], chunk, total)?;
            }
            previous = Some(end);
        }
        Ok(())
    }

    /// Reconstruct one chunk's routes during full semantic inspection. Only this
    /// chunk's routes are materialized; the caller compares the total at the end.
    pub fn verify_chunk(
        &self,
        chunk: u32,
        base: u32,
        molecules: &[MolRec],
        shapes: &[Shape],
        patterns: &[Vec<PatAlt>],
    ) -> Result<usize> {
        if self.bases.get(chunk as usize) != Some(&base) {
            bail!("access-index class base differs from source chunk");
        }
        let mut expected = Self {
            bases: vec![],
            chroms: self.chroms,
            posts: BTreeMap::new(),
        };
        let mut total = 0;
        for m in molecules {
            if m.umi_class < base {
                expected.add([0, 0, m.umi_class, 0, 0], chunk, &mut total)?;
            }
            for chain in &m.chains {
                for &(pos, shape) in &chain.reps {
                    for mode in 0..3 {
                        expected.placement(shapes, mode, m.chrom, pos, shape, chunk, &mut total)?;
                    }
                }
            }
            for &(pos, shape, pattern, _) in &m.mms {
                expected.placement(shapes, 1, m.chrom, pos, shape, chunk, &mut total)?;
                for alt in patterns
                    .get(pattern as usize)
                    .context("missing access-index pattern")?
                {
                    let apos = i64::from(pos)
                        .checked_add(alt.offset)
                        .and_then(|p| u32::try_from(p).ok())
                        .context("alternative coordinate overflow")?;
                    expected.placement(
                        shapes,
                        2,
                        alt.chrom,
                        apos,
                        if alt.shape == SAME_SHAPE {
                            shape
                        } else {
                            alt.shape
                        },
                        chunk,
                        &mut total,
                    )?;
                }
            }
        }
        for key in expected.posts.keys() {
            if self
                .posts
                .get(key)
                .is_none_or(|p| p.binary_search(&chunk).is_err())
            {
                bail!("access index omits a source-derived route");
            }
        }
        Ok(total)
    }

    pub fn posting_count(&self) -> usize {
        self.posts.values().map(Vec::len).sum()
    }

    pub fn validate_bases(&self, bases: impl Iterator<Item = u32>) -> Result<()> {
        if !self.bases[..self.bases.len()-1].iter().copied().eq(bases) {
            bail!("access-index class bases differ from the chunk directory");
        }
        Ok(())
    }

    pub fn encode(&self) -> Vec<u8> {
        let mut out = MAGIC.to_vec();
        put_varint(&mut out, self.chroms.into());
        put_varint(&mut out, self.bases.len() as u64);
        for &base in &self.bases {
            put_varint(&mut out, base.into());
        }
        put_varint(&mut out, self.posts.len() as u64);
        for (key, posts) in &self.posts {
            for &v in key {
                put_varint(&mut out, v.into());
            }
            put_varint(&mut out, posts.len() as u64);
            let mut last = 0;
            for &p in posts {
                put_varint(&mut out, (p - last).into());
                last = p;
            }
        }
        out
    }

    pub fn decode(raw: &[u8], chunks: usize, classes: u32, chroms: usize) -> Result<Self> {
        let mut c = Cursor::new(raw);
        if c.take(MAGIC.len())? != MAGIC {
            bail!("unsupported access index");
        }
        let read = |c: &mut Cursor<'_>| -> Result<u32> { Ok(u32::try_from(c.varint()?)?) };
        if read(&mut c)? as usize != chroms {
            bail!("access-index chromosome count mismatch");
        }
        let n = read(&mut c)? as usize;
        if n != chunks + 1 || n > raw.len() {
            bail!("access-index chunk count mismatch");
        }
        let mut bases = Vec::with_capacity(n);
        for _ in 0..n {
            bases.push(read(&mut c)?);
        }
        if bases.first() != Some(&0)
            || bases.last() != Some(&classes)
            || bases.windows(2).any(|w| w[0] > w[1])
        {
            bail!("invalid access-index class bases");
        }
        let rows = read(&mut c)? as usize;
        if rows > MAX_POSTS || rows > raw.len() / 7 {
            bail!("access-index rows exceed bounds");
        }
        let (mut posts, mut previous, mut total) = (BTreeMap::new(), None, 0usize);
        for _ in 0..rows {
            let mut key = [0; 5];
            for v in &mut key {
                *v = read(&mut c)?;
            }
            if previous.is_some_and(|p| p >= key) || key[0] > 2 || key[1] > 2 {
                bail!("invalid access-index key order");
            }
            if key[0] == 0 {
                if key[1] != 0 || key[2] >= classes || key[3] != 0 || key[4] != 0 {
                    bail!("invalid class route");
                }
            } else if key[2] as usize >= chroms
                || (key[0] == 1 && (key[3] > u32::MAX / TILE || key[4] != 0))
                || (key[0] == 2 && key[3] >= key[4])
            {
                bail!("invalid geometry route");
            }
            previous = Some(key);
            let count = read(&mut c)? as usize;
            total = total
                .checked_add(count)
                .context("access-index posting overflow")?;
            if count == 0 || count > chunks || count > raw.len() || total > MAX_POSTS {
                bail!("access-index postings exceed bounds");
            }
            let mut row = Vec::with_capacity(count);
            let mut last = 0u32;
            for i in 0..count {
                let delta = read(&mut c)?;
                last = last
                    .checked_add(delta)
                    .context("access-index chunk overflow")?;
                if last as usize >= chunks || (i > 0 && delta == 0) {
                    bail!("invalid access-index chunk route");
                }
                row.push(last);
            }
            posts.insert(key, row);
        }
        if !c.is_empty() {
            bail!("trailing access-index bytes");
        }
        Ok(Self {
            bases,
            chroms: chroms as u32,
            posts,
        })
    }

    pub fn class_chunks(&self, classes: impl Iterator<Item = u32>) -> Result<Vec<usize>> {
        let mut out = BTreeSet::new();
        for class in classes {
            if class >= *self.bases.last().unwrap() {
                bail!("access-index class out of range");
            }
            out.insert(self.bases.partition_point(|&b| b <= class) - 1);
            if let Some(posts) = self.posts.get(&[0, 0, class, 0, 0]) {
                out.extend(posts.iter().map(|&p| p as usize));
            }
        }
        Ok(out.into_iter().collect())
    }

    pub fn geometry_chunks(
        &self,
        junction: bool,
        mode: u32,
        chrom: u32,
        start: u32,
        end: u32,
    ) -> Vec<usize> {
        if start >= end {
            return Vec::new();
        }
        let mut out = BTreeSet::new();
        if junction {
            if let Some(p) = self.posts.get(&[2, mode, chrom, start, end]) {
                out.extend(p.iter().map(|&v| v as usize));
            }
        } else {
            for (_, p) in self
                .posts
                .range([1, mode, chrom, start / TILE, 0]..=[1, mode, chrom, (end - 1) / TILE, 0])
            {
                out.extend(p.iter().map(|&v| v as usize));
            }
        }
        out.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rows::{MolChain, MolRec, PatAlt};
    use evidence_io::archive::Shape;
    use smallvec::smallvec;

    fn fixture() -> Extracted {
        let mol = |class, pos| MolRec {
            cell: 0,
            umi_class: class,
            chrom: 0,
            strand_rev: false,
            chains: smallvec![MolChain {
                weight: 1,
                reps: smallvec![(pos, 0)]
            }],
            mms: smallvec![],
        };
        let mut mm = mol(2, 900);
        mm.chains.clear();
        mm.mms.push((900, 0, 0, 1));
        Extracted {
            mols: vec![mol(0, 100), mol(0, 70000), mol(1, 140000), mm],
            cells: vec![0],
            edges: vec![],
            n_classes: 3,
            chrom_names: vec!["chr1".into(), "chr2".into()],
            shapes: vec![Shape {
                blocks: vec![(0, 10), (100, 10)],
            }],
            patterns: vec![vec![PatAlt {
                chrom: 1,
                offset: 1000000,
                strand_flip: false,
                shape: SAME_SHAPE,
            }]],
        }
    }

    #[test]
    fn sparse_repeats_and_births_route_exact_classes() {
        let x = fixture();
        let index = Index::build(&x, [1, 1, 1, 1].into_iter()).unwrap();
        let raw = index.encode();
        let index = Index::decode(&raw, 4, 3, 2).unwrap();
        assert_eq!(index.class_chunks([0].into_iter()).unwrap(), vec![0, 1]);
        assert_eq!(index.class_chunks([1, 2].into_iter()).unwrap(), vec![2, 3]);
        assert!(index.class_chunks([3].into_iter()).is_err());
        assert_eq!(index.encode(), raw);
    }

    #[test]
    fn interval_and_alternative_junction_routes_preserve_scope() {
        let index = Index::build(&fixture(), [1, 1, 1, 1].into_iter()).unwrap();
        assert_eq!(index.geometry_chunks(false, 0, 0, 70001, 70002), vec![1]);
        assert_eq!(index.geometry_chunks(true, 2, 1, 1000910, 1001000), vec![3]);
        assert!(index
            .geometry_chunks(true, 0, 1, 1000910, 1001000)
            .is_empty());
        assert_eq!(index.geometry_chunks(true, 1, 0, 910, 1000), vec![3]);
        assert!(index.geometry_chunks(true, 2, 0, 910, 1000).is_empty());
    }

    #[test]
    fn malformed_and_truncated_indexes_fail_closed() {
        let raw = Index::build(&fixture(), [1, 1, 1, 1].into_iter())
            .unwrap()
            .encode();
        for end in 0..raw.len() {
            assert!(Index::decode(&raw[..end], 4, 3, 2).is_err());
        }
        assert!(Index::decode(&raw, 5, 3, 2).is_err());
        assert!(Index::decode(&raw, 4, 4, 2).is_err());
        let mut extra = raw;
        extra.push(0);
        assert!(Index::decode(&extra, 4, 3, 2).is_err());
    }

    #[test]
    fn semantic_verification_detects_omissions_and_extras() {
        let x = fixture();
        let mut index = Index::build(&x, [1,1,1,1].into_iter()).unwrap();
        index.validate_bases([0,1,1,2].into_iter()).unwrap();
        assert!(index.validate_bases([0,0,1,2].into_iter()).is_err());
        let verify = |index: &Index| -> Result<usize> {
            [0,1,1,2].into_iter().enumerate().map(|(i,base)|
                index.verify_chunk(i as u32,base,&x.mols[i..i+1],&x.shapes,&x.patterns)).sum()
        };
        assert_eq!(verify(&index).unwrap(), index.posting_count());
        let key = *index.posts.keys().next().unwrap();
        let saved = index.posts.remove(&key).unwrap();
        assert!(verify(&index).is_err());
        index.posts.insert(key, saved);
        index.posts.insert([1,0,0,60000,0],vec![0]);
        assert_ne!(verify(&index).unwrap(),index.posting_count());
    }
}
