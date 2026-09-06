//! Native interpreter and sound abstract interpretation of omitted unique geometry.
use super::model::*;
use crate::rows::{MolRec, PatAlt, SAME_SHAPE};
use anyhow::{bail, Context, Result};
use evidence_io::archive::Shape;
use std::collections::BTreeMap;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Alignment {
    pub chrom: u32,
    pub reverse: bool,
    pub position: u32,
    pub shape: u32,
}
pub struct Evidence {
    pub molecule: MolRec,
    pub tails: Vec<u32>,
    pub ordinal: u64,
}
pub struct Data<'a> {
    pub records: &'a [Evidence],
    pub shapes: &'a [Shape],
    pub patterns: &'a [Vec<PatAlt>],
    pub chroms: &'a [String],
}
#[derive(Clone, Copy)]
pub struct Scope<'a> {
    pub records: &'a [usize],
    pub record: Option<usize>,
    pub alignment: Option<Alignment>,
    pub signature: Option<usize>,
    pub omitted: Option<&'a crate::rows::MolChain>,
}
pub struct Evaluator<'a> {
    pub data: Data<'a>,
    pub steps: u64,
    pub max_steps: u64,
    pub unknown_omitted: u64,
}
impl Evaluator<'_> {
    pub fn eval(
        &mut self,
        n: &Node,
        s: Scope<'_>,
        fields: &BTreeMap<String, Value>,
    ) -> Result<Value> {
        self.steps += 1;
        if self.steps > self.max_steps {
            bail!("expression evaluation budget exceeded");
        }
        Ok(match &n.op {
            Op::Constant(v) => v.clone(),
            Op::Field(name) => {
                if name.starts_with("reads.") {
                    let mut unique = 0u64;
                    let mut mm = 0u64;
                    let one;
                    let indices = if let Some(i) = s.record {
                        one = [i];
                        &one[..]
                    } else {
                        s.records
                    };
                    for &i in indices {
                        let m = &self.data.records[i].molecule;
                        for c in &m.chains {
                            unique = unique
                                .checked_add(u64::from(c.weight))
                                .context("read total overflow")?;
                        }
                        for &(_, _, _, w) in &m.mms {
                            mm = mm
                                .checked_add(u64::from(w))
                                .context("read total overflow")?;
                        }
                    }
                    Value::Count(match name.as_str() {
                        "reads.unique" => unique,
                        "reads.multimapping" => mm,
                        _ => unique.checked_add(mm).context("read total overflow")?,
                    })
                } else {
                    fields
                        .get(name)
                        .cloned()
                        .with_context(|| format!("field {name} unavailable in execution scope"))?
                }
            }
            Op::Unary(op, a) => {
                let a = self.eval(a, s, fields)?;
                if op == "!" {
                    Value::Truth(a.truth()?.not())
                } else {
                    match a {
                        Value::Number(v) => Value::Number(-v),
                        Value::Null => Value::Null,
                        _ => bail!("invalid unary numeric input"),
                    }
                }
            }
            Op::Binary(op, a, b) => {
                let a = self.eval(a, s, fields)?;
                if (op == "&" && a.truth().ok() == Some(Truth::False))
                    || (op == "|" && a.truth().ok() == Some(Truth::True))
                {
                    return Ok(a);
                }
                let b = self.eval(b, s, fields)?;
                scalar_binary(op, a, b)?
            }
            Op::Geometry {
                kind,
                region,
                minimum,
                total,
            } => Value::Truth(self.geometry(s, kind, region, *minimum, *total)?),
            Op::Match(pattern) => Value::Truth(Truth::from_bool(
                self.matches(s.alignment.context("pattern requires alignment")?, pattern)?,
            )),
            Op::TxStrand(reverse) => Value::Truth(Truth::from_bool(
                s.alignment.context("strand requires alignment")?.reverse == *reverse,
            )),
            Op::Terminal(region) => {
                let record = &self.data.records[s.record.context("terminal requires record")?];
                let m = &record.molecule;
                Value::Truth(Truth::from_bool(
                    self.same(region, m.chrom, m.strand_rev)
                        && record.tails.iter().any(|&p| region.contains(p)),
                ))
            }
            Op::Nonempty(domain) => {
                let yes = if domain == "record" {
                    !s.records.is_empty()
                } else {
                    let m =
                        &self.data.records[s.record.context("nonempty requires record")?].molecule;
                    match domain.as_str() {
                        "unique" | "stored" => !m.chains.is_empty(),
                        "multimap" => !m.mms.is_empty(),
                        "alternative" => !self
                            .alternatives(
                                m,
                                s.signature.context("alternative requires signature")?,
                            )?
                            .is_empty(),
                        _ => bail!("invalid domain"),
                    }
                };
                Value::Truth(Truth::from_bool(yes))
            }
            Op::Quant { all, domain, body } => {
                let identity = if *all { Truth::True } else { Truth::False };
                let combine = |a: Truth, b: Truth| if *all { a.and(b) } else { a.or(b) };
                let mut state = identity;
                if domain == "record" {
                    for &i in s.records {
                        state = combine(
                            state,
                            self.eval(
                                body,
                                Scope {
                                    record: Some(i),
                                    alignment: None,
                                    signature: None,
                                    omitted: None,
                                    ..s
                                },
                                fields,
                            )?
                            .truth()?,
                        );
                    }
                } else {
                    let records = self.data.records;
                    let m = &records[s.record.context("quantifier requires record")?].molecule;
                    match domain.as_str() {
                        "unique" | "stored" => {
                            for chain in &m.chains {
                                validate_chain(chain)?;
                                for &(position, shape) in &chain.reps {
                                    let alignment = Some(Alignment {
                                        chrom: m.chrom,
                                        reverse: m.strand_rev,
                                        position,
                                        shape,
                                    });
                                    state = combine(
                                        state,
                                        self.eval(
                                            body,
                                            Scope {
                                                alignment,
                                                omitted: None,
                                                ..s
                                            },
                                            fields,
                                        )?
                                        .truth()?,
                                    );
                                }
                                if domain == "unique"
                                    && chain.reps.len() > 1
                                    && chain.weight as usize > chain.reps.len()
                                {
                                    let (position, shape) = chain.reps[0];
                                    let alignment = Some(Alignment {
                                        chrom: m.chrom,
                                        reverse: m.strand_rev,
                                        position,
                                        shape,
                                    });
                                    let value = self
                                        .eval(
                                            body,
                                            Scope {
                                                alignment,
                                                omitted: Some(chain),
                                                ..s
                                            },
                                            fields,
                                        )?
                                        .truth()?;
                                    state = combine(state, value);
                                    if value == Truth::Unknown {
                                        self.unknown_omitted += 1;
                                    }
                                }
                            }
                        }
                        "multimap" => {
                            for i in 0..m.mms.len() {
                                state = combine(
                                    state,
                                    self.eval(
                                        body,
                                        Scope {
                                            signature: Some(i),
                                            ..s
                                        },
                                        fields,
                                    )?
                                    .truth()?,
                                );
                            }
                        }
                        "alternative" => {
                            for alignment in self.alternatives(
                                m,
                                s.signature.context("alternative requires signature")?,
                            )? {
                                state = combine(
                                    state,
                                    self.eval(
                                        body,
                                        Scope {
                                            alignment: Some(alignment),
                                            omitted: None,
                                            ..s
                                        },
                                        fields,
                                    )?
                                    .truth()?,
                                );
                            }
                        }
                        _ => bail!("invalid quantifier domain"),
                    }
                }
                Value::Truth(state)
            }
            Op::Support(body) => {
                let records = self.data.records;
                let m = &records[s.record.context("read support requires record")?].molecule;
                let mut lower = 0u64;
                let mut upper = 0u64;
                let mut add = |v: Truth, w: u64| -> Result<()> {
                    if v == Truth::True {
                        lower = lower.checked_add(w).context("support overflow")?;
                    }
                    if v != Truth::False {
                        upper = upper.checked_add(w).context("support overflow")?;
                    }
                    Ok(())
                };
                for chain in &m.chains {
                    validate_chain(chain)?;
                    for &(position, shape) in &chain.reps {
                        let alignment = Some(Alignment {
                            chrom: m.chrom,
                            reverse: m.strand_rev,
                            position,
                            shape,
                        });
                        let weight = if chain.reps.len() == 1 {
                            u64::from(chain.weight)
                        } else {
                            1
                        };
                        add(
                            self.eval(
                                body,
                                Scope {
                                    alignment,
                                    omitted: None,
                                    ..s
                                },
                                fields,
                            )?
                            .truth()?,
                            weight,
                        )?;
                    }
                    if chain.reps.len() > 1 && chain.weight as usize > chain.reps.len() {
                        let (position, shape) = chain.reps[0];
                        let alignment = Some(Alignment {
                            chrom: m.chrom,
                            reverse: m.strand_rev,
                            position,
                            shape,
                        });
                        add(
                            self.eval(
                                body,
                                Scope {
                                    alignment,
                                    omitted: Some(chain),
                                    ..s
                                },
                                fields,
                            )?
                            .truth()?,
                            u64::from(chain.weight) - chain.reps.len() as u64,
                        )?;
                    }
                }
                for (i, &(_, _, _, weight)) in m.mms.iter().enumerate() {
                    let alternatives = self.alternatives(m, i)?;
                    let mut any = Truth::False;
                    let mut all = Truth::True;
                    for alignment in alternatives {
                        let v = self
                            .eval(
                                body,
                                Scope {
                                    alignment: Some(alignment),
                                    omitted: None,
                                    ..s
                                },
                                fields,
                            )?
                            .truth()?;
                        any = any.or(v);
                        all = all.and(v);
                    }
                    add(
                        if all == Truth::True {
                            Truth::True
                        } else if any == Truth::False {
                            Truth::False
                        } else {
                            Truth::Unknown
                        },
                        u64::from(weight),
                    )?;
                }
                Value::Bounds { lower, upper }
            }
            Op::Call(name, args) => match name.as_str() {
                "fraction" => {
                    let a = self.eval(&args[0], s, fields)?;
                    let b = self.eval(&args[1], s, fields)?;
                    let numeric = |v: Value| match v {
                        Value::Number(n) => Some(n),
                        Value::Count(n) => Some(n as f64),
                        _ => None,
                    };
                    match (numeric(a), numeric(b)) {
                        (Some(a), Some(b)) if b != 0.0 => {
                            let ratio = a / b;
                            if !ratio.is_finite() {
                                bail!("non-finite fraction");
                            }
                            Value::Number(ratio)
                        }
                        _ => Value::Null,
                    }
                }
                "bounds_lower" | "bounds_upper" | "bounds_status" => {
                    let Value::Bounds { lower, upper } = self.eval(&args[0], s, fields)? else {
                        bail!("expected read support bounds");
                    };
                    match name.as_str() {
                        "bounds_lower" => Value::Count(lower),
                        "bounds_upper" => Value::Count(upper),
                        _ => Value::String(if lower == upper { "exact" } else { "bounded" }.into()),
                    }
                }
                "is_unknown" => Value::Truth(Truth::from_bool(
                    self.eval(&args[0], s, fields)?.truth()? == Truth::Unknown,
                )),
                "is_null" => Value::Truth(Truth::from_bool(
                    self.eval(&args[0], s, fields)? == Value::Null,
                )),
                _ => bail!("aggregate {name} is only valid in summarize"),
            },
        })
    }
    pub fn same(&self, r: &Region, chrom: u32, reverse: bool) -> bool {
        self.data.chroms.get(chrom as usize) == Some(&r.chrom)
            && r.reverse.is_none_or(|v| v == reverse)
    }
    pub fn blocks(&self, a: Alignment) -> Result<smallvec::SmallVec<[(u32, u32); 8]>> {
        self.data
            .shapes
            .get(a.shape as usize)
            .context("shape identifier out of bounds")?
            .blocks
            .iter()
            .map(|&(offset, len)| {
                let start = a
                    .position
                    .checked_add(offset)
                    .context("block start overflow")?;
                let end = start.checked_add(len).context("block end overflow")?;
                if len == 0 {
                    bail!("zero-length aligned block");
                }
                Ok((start, end))
            })
            .collect()
    }
    pub fn alternatives(&self, m: &MolRec, index: usize) -> Result<Vec<Alignment>> {
        let &(pos, shape, pattern, _) = m.mms.get(index).context("invalid signature")?;
        let alternatives = self
            .data
            .patterns
            .get(pattern as usize)
            .context("invalid alternative pattern")?;
        if alternatives.is_empty() {
            bail!("multimapping signature has no alternatives");
        }
        alternatives
            .iter()
            .map(|alt| {
                Ok(Alignment {
                    chrom: alt.chrom,
                    reverse: m.strand_rev ^ alt.strand_flip,
                    position: i64::from(pos)
                        .checked_add(alt.offset)
                        .and_then(|p| u32::try_from(p).ok())
                        .context("alternative position overflow")?,
                    shape: if alt.shape == SAME_SHAPE {
                        shape
                    } else {
                        alt.shape
                    },
                })
            })
            .collect()
    }
    pub fn matches(&self, a: Alignment, p: &Pattern) -> Result<bool> {
        if self.data.chroms.get(a.chrom as usize) != Some(&p.chrom)
            || p.reverse.is_some_and(|v| v != a.reverse)
        {
            return Ok(false);
        }
        let blocks = self.blocks(a)?;
        let matches = |a: (u32, u32), b: (u32, u32)| {
            a.0.abs_diff(b.0) <= p.left && a.1.abs_diff(b.1) <= p.right
        };
        // Single-junction predicates need no intermediate junction vector.
        if p.junctions.len() == 1 {
            return Ok(blocks
                .windows(2)
                .any(|w| matches((w[0].1, w[1].0), p.junctions[0])));
        }
        let observed: smallvec::SmallVec<[(u32, u32); 8]> =
            blocks.windows(2).map(|w| (w[0].1, w[1].0)).collect();
        if p.subsequence {
            let mut next = 0;
            for j in observed {
                if matches(j, p.junctions[next]) {
                    next += 1;
                    if next == p.junctions.len() {
                        return Ok(true);
                    }
                }
            }
            Ok(false)
        } else {
            Ok(observed
                .windows(p.junctions.len())
                .any(|w| w.iter().zip(&p.junctions).all(|(&a, &b)| matches(a, b))))
        }
    }
    pub(super) fn geometry(
        &self,
        s: Scope<'_>,
        kind: &str,
        r: &Region,
        minimum: u32,
        total: bool,
    ) -> Result<Truth> {
        let a = s.alignment.context("geometry requires alignment")?;
        if !self.same(r, a.chrom, a.reverse) {
            return Ok(Truth::False);
        }
        if let Some(chain) = s.omitted {
            let starts = chain
                .reps
                .iter()
                .map(|&(position, shape)| {
                    let first = self
                        .data
                        .shapes
                        .get(shape as usize)
                        .and_then(|s| s.blocks.first())
                        .context("empty chain shape")?;
                    if first.0 != 0 {
                        bail!("start interval proof requires normalized shape offset zero");
                    }
                    Ok(position)
                })
                .collect::<Result<smallvec::SmallVec<[_; 2]>>>()?;
            let lo = *starts.iter().min().unwrap();
            let hi = *starts.iter().max().unwrap();
            if kind == "start_in" {
                if r.intervals.iter().any(|&(l, h)| l <= lo && hi < h) {
                    return Ok(Truth::True);
                }
                if !r.intervals.iter().any(|&(l, h)| l <= hi && lo < h) {
                    return Ok(Truth::False);
                }
            } else {
                // Junctions are invariant, so intronic gaps cannot gain aligned overlap.
                // Only starts are globally bracketed. Ends are bracketed *only when all
                // starts agree*, because end is the extractor's lexicographic tie-breaker.
                let blocks = self.blocks(a)?;
                let junctions: smallvec::SmallVec<[_; 8]> =
                    blocks.windows(2).map(|w| (w[0].1, w[1].0)).collect();
                let (min_end, max_end) = if lo == hi {
                    let ends = chain
                        .reps
                        .iter()
                        .map(|&(position, shape)| {
                            self.blocks(Alignment {
                                position,
                                shape,
                                ..a
                            })?
                            .last()
                            .map(|b| b.1)
                            .context("empty shape")
                        })
                        .collect::<Result<smallvec::SmallVec<[_; 2]>>>()?;
                    (*ends.iter().min().unwrap(), *ends.iter().max().unwrap())
                } else {
                    (
                        junctions
                            .last()
                            .map_or(lo, |j| j.1)
                            .checked_add(1)
                            .context("end bound overflow")?,
                        u32::MAX,
                    )
                };
                if kind == "end_in" {
                    if r.intervals
                        .iter()
                        .any(|&(l, h)| l <= min_end && max_end < h)
                    {
                        return Ok(Truth::True);
                    }
                    if !r
                        .intervals
                        .iter()
                        .any(|&(l, h)| l <= max_end && min_end < h)
                    {
                        return Ok(Truth::False);
                    }
                } else if kind == "overlaps" {
                    let mut upper = smallvec::SmallVec::<[_; 8]>::new();
                    let mut lower = smallvec::SmallVec::<[_; 8]>::new();
                    let mut left = lo;
                    let mut inner_left = hi;
                    for &(donor, acceptor) in &junctions {
                        upper.push((left, donor));
                        lower.push((inner_left, donor));
                        left = acceptor;
                        inner_left = acceptor;
                    }
                    upper.push((left, max_end));
                    if inner_left < min_end {
                        lower.push((inner_left, min_end));
                    }
                    let overlap = |blocks: &[(u32, u32)]| {
                        let counts = blocks.iter().map(|&(a, b)| {
                            r.intervals
                                .iter()
                                .map(|&(c, d)| u64::from(b.min(d).saturating_sub(a.max(c))))
                                .sum::<u64>()
                        });
                        if total {
                            counts.sum::<u64>()
                        } else {
                            counts.max().unwrap_or(0)
                        }
                    };
                    if overlap(&upper) < u64::from(minimum) {
                        return Ok(Truth::False);
                    }
                    if overlap(&lower) >= u64::from(minimum) {
                        return Ok(Truth::True);
                    }
                }
            }
            return Ok(Truth::Unknown);
        }
        let blocks = self.blocks(a)?;
        let answer = match kind {
            "start_in" => blocks.first().is_some_and(|&(p, _)| r.contains(p)),
            "end_in" => blocks.last().is_some_and(|&(_, p)| r.contains(p)),
            "overlaps" => {
                let mut counts = smallvec::SmallVec::<[u64; 8]>::new();
                for &(a, b) in &blocks {
                    let mut n = 0u64;
                    for &(c, d) in &r.intervals {
                        n += u64::from(b.min(d).saturating_sub(a.max(c)));
                    }
                    counts.push(n);
                }
                if total {
                    counts.iter().sum::<u64>() >= u64::from(minimum)
                } else {
                    counts.into_iter().max().unwrap_or(0) >= u64::from(minimum)
                }
            }
            _ => bail!("unknown geometry operation"),
        };
        Ok(Truth::from_bool(answer))
    }
}
pub(super) fn validate_chain(chain: &crate::rows::MolChain) -> Result<()> {
    if chain.reps.is_empty() || chain.weight < (chain.reps.len() as u32) {
        bail!("invalid chain multiplicity or representatives");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rows::MolChain;
    use smallvec::smallvec;
    fn node(q: &str) -> Node {
        let doc=super::super::parser::parse(&format!("header {{gq=1,assembly=\"test\"}} from @x.records |> within all |> derive {{a={q}}} |> select {{a}}")).unwrap();
        let p = super::super::check::Compiler::new(&doc, BTreeMap::new())
            .unwrap()
            .compile()
            .unwrap();
        let super::super::model::Stage::Derive(fields) = &p.stages[0] else {
            panic!()
        };
        fields[0].1.clone()
    }
    fn evaluate(q: &str, weight: u32) -> Value {
        let shapes = vec![
            Shape {
                blocks: vec![(0, 20)],
            },
            Shape {
                blocks: vec![(0, 20)],
            },
        ];
        let records = vec![Evidence {
            molecule: MolRec {
                cell: 0,
                umi_class: 0,
                chrom: 0,
                strand_rev: false,
                chains: smallvec![MolChain {
                    weight,
                    reps: smallvec![(100, 0), (110, 1)]
                }],
                mms: smallvec![],
            },
            tails: vec![],
            ordinal: 0,
        }];
        let mut evaluator = Evaluator {
            data: Data {
                records: &records,
                shapes: &shapes,
                patterns: &[],
                chroms: &["1".into()],
            },
            steps: 0,
            max_steps: 1000,
            unknown_omitted: 0,
        };
        evaluator
            .eval(
                &node(q),
                Scope {
                    records: &[0],
                    record: Some(0),
                    alignment: None,
                    signature: None,
                    omitted: None,
                },
                &BTreeMap::new(),
            )
            .unwrap()
    }
    #[test]
    fn omitted_end_counterexample() {
        assert_eq!(
            evaluate("any unique {overlaps(g[1:180..190])}", 3),
            Value::Truth(Truth::Unknown)
        );
        assert_eq!(
            evaluate("any stored {overlaps(g[1:180..190])}", 3),
            Value::Truth(Truth::False)
        );
        assert_eq!(
            evaluate("any unique {overlaps(g[1:180..190])}", 2),
            Value::Truth(Truth::False)
        );
        assert_eq!(
            evaluate("support_reads(overlaps(g[1:180..190]))", 3),
            Value::Bounds { lower: 0, upper: 1 }
        );
    }
    #[test]
    fn start_proofs_and_read_witnesses() {
        assert_eq!(
            evaluate("all unique {start_in(g[1:90..120])}", 10),
            Value::Truth(Truth::True)
        );
        assert_eq!(
            evaluate("any unique {start_in(g[1:104..106])}", 10),
            Value::Truth(Truth::Unknown)
        );
        assert_eq!(
            evaluate("support_reads(start_in(g[1:90..120]))", 10),
            Value::Bounds {
                lower: 10,
                upper: 10
            }
        );
        assert_eq!(
            evaluate("support_reads(start_in(g[1:100..105]))", 10),
            Value::Bounds { lower: 1, upper: 9 }
        );
    }
    #[test]
    fn kleene_laws() {
        for p in [Truth::True, Truth::False, Truth::Unknown] {
            for q in [Truth::True, Truth::False, Truth::Unknown] {
                assert_eq!(p.and(q).not(), p.not().or(q.not()));
            }
        }
        assert_eq!(Truth::Unknown.or(Truth::Unknown.not()), Truth::Unknown);
    }
    #[test]
    fn vacuous_all() {
        assert_eq!(
            evaluate("all multimap {all alternative {true}}", 2),
            Value::Truth(Truth::True)
        );
        assert_eq!(
            evaluate(
                "nonempty(multimap) & all multimap {all alternative {true}}",
                2
            ),
            Value::Truth(Truth::False)
        );
    }
    #[test]
    fn exhaustive_compact_bounds_are_sound() {
        let geometries: Vec<_> = (1..5)
            .flat_map(|start| (start + 1..8).map(move |end| (start, end)))
            .collect();
        for i in 0..geometries.len() {
            for j in i..geometries.len() {
                for k in j..geometries.len() {
                    let reads = [geometries[i], geometries[j], geometries[k]];
                    let shapes: Vec<_> = reads
                        .iter()
                        .map(|&(a, b)| Shape {
                            blocks: vec![(0, b - a)],
                        })
                        .collect();
                    let lo = (0..3)
                        .min_by_key(|&i| (reads[i].0, std::cmp::Reverse(reads[i].1)))
                        .unwrap();
                    let hi = (0..3)
                        .max_by_key(|&i| (reads[i].0, std::cmp::Reverse(reads[i].1)))
                        .unwrap();
                    let mut reps = smallvec![(reads[lo].0, lo as u32)];
                    if reads[hi] != reads[lo] {
                        reps.push((reads[hi].0, hi as u32));
                    }
                    let records = vec![Evidence {
                        molecule: MolRec {
                            cell: 0,
                            umi_class: 0,
                            chrom: 0,
                            strand_rev: false,
                            chains: smallvec![MolChain { weight: 3, reps }],
                            mms: smallvec![],
                        },
                        tails: vec![],
                        ordinal: 0,
                    }];
                    let mut ev = Evaluator {
                        data: Data {
                            records: &records,
                            shapes: &shapes,
                            patterns: &[],
                            chroms: &["1".into()],
                        },
                        steps: 0,
                        max_steps: 100_000,
                        unknown_omitted: 0,
                    };
                    let scope = Scope {
                        records: &[0],
                        record: Some(0),
                        alignment: None,
                        signature: None,
                        omitted: None,
                    };
                    for start in 0..8 {
                        for kind in ["overlaps", "start_in", "end_in"] {
                            let predicate = Node {
                                at: 0,
                                ty: Type::Truth,
                                op: Op::Geometry {
                                    kind: kind.into(),
                                    region: Region {
                                        chrom: "1".into(),
                                        intervals: vec![(start, start + 2)],
                                        reverse: None,
                                    },
                                    minimum: 1,
                                    total: false,
                                },
                            };
                            let actual = reads
                                .iter()
                                .enumerate()
                                .map(|(i, &(position, _))| {
                                    ev.eval(
                                        &predicate,
                                        Scope {
                                            alignment: Some(Alignment {
                                                chrom: 0,
                                                reverse: false,
                                                position,
                                                shape: i as u32,
                                            }),
                                            ..scope
                                        },
                                        &BTreeMap::new(),
                                    )
                                    .unwrap()
                                    .truth()
                                    .unwrap()
                                })
                                .filter(|v| *v == Truth::True)
                                .count() as u64;
                            let Value::Bounds { lower, upper } = ev
                                .eval(
                                    &Node {
                                        at: 0,
                                        ty: Type::Bounds,
                                        op: Op::Support(Box::new(predicate)),
                                    },
                                    scope,
                                    &BTreeMap::new(),
                                )
                                .unwrap()
                            else {
                                panic!()
                            };
                            assert!(
                                lower <= actual && actual <= upper,
                                "{reads:?}, {kind}, {start}: {lower} <= {actual} <= {upper}"
                            );
                        }
                    }
                }
            }
        }
    }
}
