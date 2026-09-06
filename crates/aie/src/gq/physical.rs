//! Conservative physical lowering. Logical scopes and evaluation order are unchanged.
//! Truth slots avoid materializing a named row between filter, derive and aggregate.
use super::{
    eval::*,
    model::*,
    native::{Leaf, SharedCache},
};
use anyhow::{bail, Context, Result};
use smallvec::SmallVec;
use std::collections::BTreeMap;

type Row = BTreeMap<String, Value>;
pub(super) type States = SmallVec<[Truth; 8]>;
pub(super) type Values = SmallVec<[Value; 8]>;

enum Expr<'a> {
    Native(&'a Node),
    Constant(Truth),
    Slot(usize),
    Not(Box<Self>),
    And(Box<Self>, Box<Self>),
    Or(Box<Self>, Box<Self>),
    Quant {
        all: bool,
        domain: &'a str,
        body: Box<Self>,
    },
    Placement(Leaf<'a>),
    Strand(bool),
}

fn references_slot(n: &Node, slots: &BTreeMap<String, usize>) -> bool {
    match &n.op {
        Op::Field(name) => slots.contains_key(name),
        Op::Unary(_, a) | Op::Quant { body: a, .. } | Op::Support(a) => references_slot(a, slots),
        Op::Binary(_, a, b) => references_slot(a, slots) || references_slot(b, slots),
        Op::Call(_, args) => args.iter().any(|a| references_slot(a, slots)),
        _ => false,
    }
}
impl<'a> Expr<'a> {
    fn compile(n: &'a Node, slots: &BTreeMap<String, usize>, cache: &SharedCache) -> Option<Self> {
        if n.ty != Type::Truth {
            return None;
        }
        Some(match &n.op {
            Op::Constant(Value::Truth(t)) => Self::Constant(*t),
            Op::Field(name) if slots.contains_key(name) => Self::Slot(slots[name]),
            Op::Unary(op, a) if op == "!" => Self::Not(Box::new(Self::compile(a, slots, cache)?)),
            Op::Binary(op, a, b) if op == "&" || op == "|" => {
                let (a, b) = (
                    Box::new(Self::compile(a, slots, cache)?),
                    Box::new(Self::compile(b, slots, cache)?),
                );
                if op == "&" {
                    Self::And(a, b)
                } else {
                    Self::Or(a, b)
                }
            }
            Op::Quant { all, domain, body } => Self::Quant {
                all: *all,
                domain,
                body: Box::new(Self::compile(body, slots, cache)?),
            },
            Op::Geometry { .. } | Op::Match(_) => Self::Placement(Leaf::compile(n, cache)),
            Op::TxStrand(r) => Self::Strand(*r),
            _ if !references_slot(n, slots) => Self::Native(n),
            _ => return None,
        })
    }

    fn eval(
        &self,
        ev: &mut Evaluator<'_>,
        scope: Scope<'_>,
        row: &Row,
        slots: &[Truth],
    ) -> Result<Truth> {
        // Native leaves account for their own logical work. All other nodes charge exactly
        // one step as the reference interpreter does, including slot reads and quantifiers.
        if let Self::Native(n) = self {
            return ev.eval(n, scope, row)?.truth();
        }
        ev.steps += 1;
        if ev.steps > ev.max_steps {
            bail!("expression evaluation budget exceeded");
        }
        Ok(match self {
            Self::Native(_) => unreachable!(),
            Self::Constant(t) => *t,
            Self::Slot(i) => slots[*i],
            Self::Not(a) => a.eval(ev, scope, row, slots)?.not(),
            Self::And(a, b) => {
                let a = a.eval(ev, scope, row, slots)?;
                if a == Truth::False {
                    a
                } else {
                    a.and(b.eval(ev, scope, row, slots)?)
                }
            }
            Self::Or(a, b) => {
                let a = a.eval(ev, scope, row, slots)?;
                if a == Truth::True {
                    a
                } else {
                    a.or(b.eval(ev, scope, row, slots)?)
                }
            }
            Self::Placement(leaf) => leaf.eval(ev, scope)?,
            Self::Strand(r) => Truth::from_bool(
                scope
                    .alignment
                    .context("strand requires alignment")?
                    .reverse
                    == *r,
            ),
            Self::Quant { all, domain, body } => {
                let mut state = if *all { Truth::True } else { Truth::False };
                let combine = |a: Truth, b| if *all { a.and(b) } else { a.or(b) };
                if *domain == "record" {
                    for &i in scope.records {
                        state = combine(
                            state,
                            body.eval(
                                ev,
                                Scope {
                                    record: Some(i),
                                    alignment: None,
                                    signature: None,
                                    omitted: None,
                                    ..scope
                                },
                                row,
                                slots,
                            )?,
                        );
                    }
                } else {
                    let records = ev.data.records;
                    let m = &records[scope.record.context("quantifier requires record")?].molecule;
                    match *domain {
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
                                        body.eval(
                                            ev,
                                            Scope {
                                                alignment,
                                                omitted: None,
                                                ..scope
                                            },
                                            row,
                                            slots,
                                        )?,
                                    );
                                }
                                if *domain == "unique"
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
                                    let v = body.eval(
                                        ev,
                                        Scope {
                                            alignment,
                                            omitted: Some(chain),
                                            ..scope
                                        },
                                        row,
                                        slots,
                                    )?;
                                    state = combine(state, v);
                                    if v == Truth::Unknown {
                                        ev.unknown_omitted += 1;
                                    }
                                }
                            }
                        }
                        "multimap" => {
                            for i in 0..m.mms.len() {
                                state = combine(
                                    state,
                                    body.eval(
                                        ev,
                                        Scope {
                                            signature: Some(i),
                                            ..scope
                                        },
                                        row,
                                        slots,
                                    )?,
                                );
                            }
                        }
                        "alternative" => {
                            for alignment in ev.alternatives(
                                m,
                                scope.signature.context("alternative requires signature")?,
                            )? {
                                state = combine(
                                    state,
                                    body.eval(
                                        ev,
                                        Scope {
                                            alignment: Some(alignment),
                                            omitted: None,
                                            ..scope
                                        },
                                        row,
                                        slots,
                                    )?,
                                );
                            }
                        }
                        _ => bail!("invalid quantifier domain"),
                    }
                }
                // Do not short-circuit domain traversal: omitted-cause counters and logical
                // budgets are part of the reference contract, even after a decisive witness.
                state
            }
        })
    }
}
enum Step<'a> {
    Filter(Expr<'a>),
    Derive(usize, Expr<'a>),
}
enum Aggregate<'a> {
    Count,
    Truth(Expr<'a>, Truth),
    Sum(&'a Node),
    Distinct(&'a Node),
}
enum Sink<'a> {
    Tally(Vec<Expr<'a>>),
    Summary(Vec<Aggregate<'a>>),
}
pub(super) enum Outcome {
    Dropped(Truth),
    Tally(States),
    Summary,
}
pub(super) struct Kernel<'a> {
    steps: Vec<Step<'a>>,
    sink: Sink<'a>,
    slots: Vec<Truth>,
    pub values: Values,
}
impl<'a> Kernel<'a> {
    pub fn compile(plan: &'a Plan) -> Option<Self> {
        if plan.enumeration {
            return None;
        }
        let mut names = BTreeMap::new();
        let cache = SharedCache::default();
        let mut steps = Vec::new();
        let mut sink = None;
        let mut size = 0;
        for stage in &plan.stages {
            if sink.is_some() {
                if matches!(stage, Stage::Select(_) | Stage::Sort(_) | Stage::Take(_)) {
                    continue;
                }
                return None;
            }
            match stage {
                Stage::Where(n) => steps.push(Step::Filter(Expr::compile(n, &names, &cache)?)),
                Stage::Derive(fields) => {
                    for (name, n) in fields {
                        steps.push(Step::Derive(size, Expr::compile(n, &names, &cache)?));
                        names.insert(name.clone(), size);
                        size += 1;
                    }
                }
                Stage::Tally(fields, _) => {
                    sink = Some(Sink::Tally(
                        fields
                            .iter()
                            .map(|(_, n)| Expr::compile(n, &names, &cache))
                            .collect::<Option<_>>()?,
                    ))
                }
                Stage::Summarize(fields, _) => {
                    let mut aggregates = Vec::new();
                    for (_, n) in fields {
                        let Op::Call(name, args) = &n.op else {
                            return None;
                        };
                        aggregates.push(match name.as_str() {
                            "count" => Aggregate::Count,
                            "count_true" | "count_false" | "count_unknown" => Aggregate::Truth(
                                Expr::compile(&args[0], &names, &cache)?,
                                match name.as_str() {
                                    "count_true" => Truth::True,
                                    "count_false" => Truth::False,
                                    _ => Truth::Unknown,
                                },
                            ),
                            "sum" if !references_slot(&args[0], &names) => Aggregate::Sum(&args[0]),
                            "count_distinct" if !references_slot(&args[0], &names) => {
                                Aggregate::Distinct(&args[0])
                            }
                            _ => return None,
                        });
                    }
                    sink = Some(Sink::Summary(aggregates));
                }
                _ => return None,
            }
        }
        Some(Self {
            steps,
            sink: sink?,
            slots: vec![Truth::False; size],
            values: Values::new(),
        })
    }
    pub fn run(&mut self, ev: &mut Evaluator<'_>, scope: Scope<'_>, row: &Row) -> Result<Outcome> {
        for step in &self.steps {
            match step {
                Step::Filter(expr) => {
                    let state = expr.eval(ev, scope, row, &self.slots)?;
                    if state != Truth::True {
                        return Ok(Outcome::Dropped(state));
                    }
                }
                Step::Derive(slot, expr) => {
                    self.slots[*slot] = expr.eval(ev, scope, row, &self.slots)?
                }
            }
        }
        Ok(match &self.sink {
            Sink::Tally(fields) => Outcome::Tally(
                fields
                    .iter()
                    .map(|e| e.eval(ev, scope, row, &self.slots))
                    .collect::<Result<_>>()?,
            ),
            Sink::Summary(fields) => {
                self.values.clear();
                for e in fields {
                    self.values.push(match e {
                        Aggregate::Count => Value::Count(1),
                        Aggregate::Truth(expr, want) => Value::Count(u64::from(
                            expr.eval(ev, scope, row, &self.slots)? == *want,
                        )),
                        Aggregate::Sum(n) | Aggregate::Distinct(n) => ev.eval(n, scope, row)?,
                    });
                }
                Outcome::Summary
            }
        })
    }
}

pub(super) fn description(plan: &Plan, reference: bool, parallel_decode: bool) -> String {
    if plan.enumeration {
        return "dedicated exact junction enumeration (both engines)".into();
    }
    let pipeline = if reference {
        "reference interpreter (requested)"
    } else if Kernel::compile(plan).is_some() {
        "fused filter/derive/aggregate; shared native shape predicates; bounded placement cache; Truth slots; sparse typed aggregation"
    } else {
        "reference interpreter (unsupported fused shape)"
    };
    let top_k = !reference
        && plan
            .stages
            .windows(2)
            .any(|s| matches!((&s[0], &s[1]), (Stage::Sort(_), Stage::Take(_))));
    format!(
        "{pipeline}; ordered parallel chunk decode={}; member-constant grouping={}; direct record iteration={}; unit.id materialized={}; stable top-k={top_k}",
        !reference && parallel_decode,
        !reference && plan.stages.iter().any(|stage| match stage {
            Stage::Tally(_, fields) | Stage::Summarize(_, fields) => fields.iter().all(|(_, n)| matches!(&n.op, Op::Field(name) if name == "sample")),
            _ => false,
        }),
        !reference && plan.unit == Unit::Record && !plan.enumeration,
        reference || uses_field(plan, "unit.id")
    )
}

pub(super) fn uses_field(plan: &Plan, field: &str) -> bool {
    let names = BTreeMap::from([(field.to_owned(), 0)]);
    plan.stages.iter().any(|stage| match stage {
        Stage::Where(n) => references_slot(n, &names),
        Stage::Derive(fields) | Stage::Select(fields) | Stage::Sort(fields) => {
            fields.iter().any(|(_, n)| references_slot(n, &names))
        }
        Stage::Tally(fields, by) | Stage::Summarize(fields, by) => fields
            .iter()
            .chain(by)
            .any(|(_, n)| references_slot(n, &names)),
        Stage::Take(_) => false,
    })
}
