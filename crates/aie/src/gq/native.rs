//! Bounded memoization of pure placement predicates, using the built-in shape kernel
//! where its semantics match GQ. Never cache extent-sensitive omitted geometry.
use super::{eval::*, model::*};
use crate::querycmd::ShapePredicate;
use anyhow::{bail, Context, Result};
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
};

const MAX_ENTRIES: usize = 1_024;
type Entry = Option<((usize, Alignment), Truth)>;
#[derive(Default)]
pub(super) struct Cache {
    next: usize,
    entries: Vec<Entry>,
}
pub(super) type SharedCache = Rc<RefCell<Cache>>;

pub(super) struct Leaf<'a> {
    node: &'a Node,
    matcher: Option<ShapePredicate>,
    cache: SharedCache,
    id: usize,
    validated_shape: Cell<Option<(u32, u32)>>,
}
impl<'a> Leaf<'a> {
    pub fn compile(node: &'a Node, cache: &SharedCache) -> Self {
        let matcher = match &node.op {
            Op::Geometry {
                kind,
                region,
                minimum,
                total,
            } if region.intervals.len() == 1 && !total => {
                let (start, end) = region.intervals[0];
                Some(ShapePredicate::new(kind, start, end, *minimum, vec![]))
            }
            Op::Match(p) if p.junctions.len() == 1 && p.left == p.right => {
                let (start, end) = p.junctions[0];
                Some(ShapePredicate::new("junction", start, end, p.left, vec![]))
            }
            // The native path matcher has different restart behavior for overlapping
            // patterns; only single-junction patterns are lowered for now.
            _ => None,
        };
        let id = cache.borrow().next;
        cache.borrow_mut().next += 1;
        Self {
            node,
            matcher,
            cache: Rc::clone(cache),
            id,
            validated_shape: Cell::new(None),
        }
    }

    pub fn eval(&self, ev: &Evaluator<'_>, scope: Scope<'_>) -> Result<Truth> {
        let alignment = scope.alignment.context("predicate requires alignment")?;
        let (chrom, reverse) = match &self.node.op {
            Op::Geometry {
                kind,
                region,
                minimum,
                total,
            } => {
                if scope.omitted.is_some() {
                    return ev.geometry(scope, kind, region, *minimum, *total);
                }
                (&region.chrom, region.reverse)
            }
            Op::Match(p) => (&p.chrom, p.reverse),
            _ => unreachable!(),
        };
        if ev.data.chroms.get(alignment.chrom as usize) != Some(chrom)
            || reverse.is_some_and(|r| r != alignment.reverse)
        {
            return Ok(Truth::False);
        }
        let key = (self.id, alignment);
        let slot = (alignment.position as usize).wrapping_mul(0x9e3779b1)
            ^ (alignment.shape as usize).wrapping_mul(0x85ebca6b)
            ^ self.id.wrapping_mul(0xc2b2ae35);
        let slot = (slot ^ (slot >> 16)) & (MAX_ENTRIES - 1);
        if let Some(Some((stored, value))) = self.cache.borrow().entries.get(slot) {
            if *stored == key {
                return Ok(*value);
            }
        }
        let value = if let Some(matcher) = &self.matcher {
            let shape = ev
                .data
                .shapes
                .get(alignment.shape as usize)
                .context("shape identifier out of bounds")?;
            // Validate every block without materializing absolute coordinates. The
            // maximum relative end proves every addition safe at this position.
            let end = if let Some((_, end)) = self
                .validated_shape
                .get()
                .filter(|&(id, _)| id == alignment.shape)
            {
                end
            } else {
                let mut end = 0;
                for &(offset, length) in &shape.blocks {
                    if length == 0 {
                        bail!("zero-length aligned block");
                    }
                    end = end.max(offset.checked_add(length).context("block end overflow")?);
                }
                self.validated_shape.set(Some((alignment.shape, end)));
                end
            };
            alignment
                .position
                .checked_add(end)
                .context("block end overflow")?;
            Truth::from_bool(matcher.matches(alignment.position, shape)?)
        } else {
            match &self.node.op {
                Op::Geometry {
                    kind,
                    region,
                    minimum,
                    total,
                } => ev.geometry(scope, kind, region, *minimum, *total)?,
                Op::Match(p) => Truth::from_bool(ev.matches(alignment, p)?),
                _ => unreachable!(),
            }
        };
        let mut cache = self.cache.borrow_mut();
        if cache.entries.is_empty() {
            cache.entries.resize(MAX_ENTRIES, None);
        }
        // Direct-mapped replacement: fixed allocation, no growth or hash-table probes.
        cache.entries[slot] = Some((key, value));
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rows::MolChain;
    use evidence_io::archive::Shape;
    use smallvec::smallvec;

    #[test]
    fn shared_predicates_match_interpreter_and_do_not_cache_omitted_extents() {
        let shapes = vec![Shape {
            blocks: vec![(0, 10), (30, 20)],
        }];
        let chroms = vec!["1".into()];
        let ev = Evaluator {
            data: Data {
                records: &[],
                shapes: &shapes,
                patterns: &[],
                chroms: &chroms,
            },
            steps: 0,
            max_steps: u64::MAX,
            unknown_omitted: 0,
        };
        let cache = SharedCache::default();
        let region = Region {
            chrom: "1".into(),
            reverse: None,
            intervals: vec![(103, 116)],
        };
        for kind in ["overlaps", "start_in", "end_in"] {
            for total in [false, true] {
                let n = Node {
                    at: 0,
                    ty: Type::Truth,
                    op: Op::Geometry {
                        kind: kind.into(),
                        region: region.clone(),
                        minimum: 4,
                        total,
                    },
                };
                let leaf = Leaf::compile(&n, &cache);
                for position in 90..125 {
                    let a = Alignment {
                        chrom: 0,
                        reverse: false,
                        position,
                        shape: 0,
                    };
                    let scope = Scope {
                        records: &[],
                        record: None,
                        alignment: Some(a),
                        signature: None,
                        omitted: None,
                    };
                    for _ in 0..2 {
                        assert_eq!(
                            leaf.eval(&ev, scope).unwrap(),
                            ev.geometry(scope, kind, &region, 4, total).unwrap()
                        );
                    }
                    let chain = MolChain {
                        weight: 8,
                        reps: smallvec![(position, 0), (position + 10, 0)],
                    };
                    let omitted = Scope {
                        omitted: Some(&chain),
                        ..scope
                    };
                    assert_eq!(
                        leaf.eval(&ev, omitted).unwrap(),
                        ev.geometry(omitted, kind, &region, 4, total).unwrap()
                    );
                }
            }
        }
        for (left, right) in [(0, 0), (3, 3), (2, 7)] {
            let p = Pattern {
                chrom: "1".into(),
                reverse: Some(false),
                junctions: vec![(110, 130)],
                subsequence: false,
                left,
                right,
            };
            let n = Node {
                at: 0,
                ty: Type::Truth,
                op: Op::Match(p.clone()),
            };
            let leaf = Leaf::compile(&n, &cache);
            for position in 90..125 {
                for reverse in [false, true] {
                    let a = Alignment {
                        chrom: 0,
                        reverse,
                        position,
                        shape: 0,
                    };
                    let scope = Scope {
                        records: &[],
                        record: None,
                        alignment: Some(a),
                        signature: None,
                        omitted: None,
                    };
                    for _ in 0..2 {
                        assert_eq!(
                            leaf.eval(&ev, scope).unwrap(),
                            Truth::from_bool(ev.matches(a, &p).unwrap())
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn cache_is_bounded_across_leaves_and_does_not_hide_invalid_blocks() {
        let shapes = vec![
            Shape {
                blocks: vec![(0, 10)],
            },
            Shape {
                blocks: vec![(0, 10), (20, 0)],
            },
        ];
        let chroms = vec!["1".into()];
        let ev = Evaluator {
            data: Data {
                records: &[],
                shapes: &shapes,
                patterns: &[],
                chroms: &chroms,
            },
            steps: 0,
            max_steps: u64::MAX,
            unknown_omitted: 0,
        };
        let n = Node {
            at: 0,
            ty: Type::Truth,
            op: Op::Geometry {
                kind: "overlaps".into(),
                region: Region {
                    chrom: "1".into(),
                    reverse: None,
                    intervals: vec![(0, 100_000)],
                },
                minimum: 1,
                total: false,
            },
        };
        let cache = SharedCache::default();
        let leaves = [Leaf::compile(&n, &cache), Leaf::compile(&n, &cache)];
        let scope = Scope {
            records: &[],
            record: None,
            alignment: Some(Alignment {
                chrom: 0,
                reverse: false,
                position: 0,
                shape: 0,
            }),
            signature: None,
            omitted: None,
        };
        for position in 0..40_000 {
            let s = Scope {
                alignment: Some(Alignment {
                    position,
                    ..scope.alignment.unwrap()
                }),
                ..scope
            };
            for leaf in &leaves {
                assert_eq!(leaf.eval(&ev, s).unwrap(), Truth::True);
            }
        }
        assert_eq!(cache.borrow().entries.len(), MAX_ENTRIES);
        for a in [
            Alignment {
                shape: 1,
                ..scope.alignment.unwrap()
            },
            Alignment {
                position: u32::MAX,
                ..scope.alignment.unwrap()
            },
            Alignment {
                shape: 9,
                ..scope.alignment.unwrap()
            },
        ] {
            assert!(leaves[0]
                .eval(
                    &ev,
                    Scope {
                        alignment: Some(a),
                        ..scope
                    }
                )
                .is_err());
        }
    }
}
