use super::{
    ast::{self, Expr, Kind},
    model::*,
};
use anyhow::{bail, Context, Result};
use std::collections::{BTreeMap, BTreeSet};

pub struct Compiler<'a> {
    doc: &'a ast::Document,
    pub fields: BTreeMap<String, Type>,
    definitions: BTreeMap<String, Expr>,
    active: BTreeSet<String>,
    nodes: usize,
    pub terminal: bool,
    pub stored: bool,
    library: Option<String>,
    resources: Option<&'a super::resources::Resources>,
    aggregate: bool,
    assembly: Option<String>,
    depth: usize,
}
impl<'a> Compiler<'a> {
    pub fn new(doc: &'a ast::Document, fields: BTreeMap<String, Type>) -> Result<Self> {
        let mut definitions = BTreeMap::new();
        for (n, e) in &doc.bindings {
            if definitions.insert(n.clone(), e.clone()).is_some() {
                bail!("duplicate let {n}");
            }
        }
        let mut names: BTreeSet<_> = definitions.keys().cloned().collect();
        for f in &doc.functions {
            if !names.insert(f.name.clone()) {
                bail!("duplicate definition {}", f.name);
            }
        }
        validate_definitions(doc)?;
        Ok(Self {
            doc,
            fields,
            definitions,
            active: BTreeSet::new(),
            nodes: 0,
            terminal: false,
            stored: false,
            library: None,
            resources: None,
            aggregate: false,
            assembly: None,
            depth: 0,
        })
    }
    pub fn with_resources(mut self, resources: &'a super::resources::Resources) -> Self {
        self.resources = Some(resources);
        self
    }
    fn constant(&mut self, e: &Expr) -> Result<Value> {
        let node = self.node(e, Unit::Alignment, &BTreeMap::new())?;
        if let Op::Constant(v) = node.op {
            Ok(v)
        } else {
            bail!("byte {}: expected compile-time value", e.at)
        }
    }
    fn local_constant(
        &mut self,
        e: &Expr,
        u: Unit,
        locals: &BTreeMap<String, Node>,
    ) -> Result<Value> {
        let node = self.node(e, u, locals)?;
        match node.op {
            Op::Constant(v) => Ok(v),
            _ => bail!("byte {}: option requires a compile-time value", e.at),
        }
    }
    fn predicate(&mut self, e: &Expr, u: Unit, locals: &BTreeMap<String, Node>) -> Result<Node> {
        let mut n = self.node(e, u, locals)?;
        if let Op::Constant(Value::Pattern(p)) = &n.op {
            require_unit(e.at, u, Unit::Alignment)?;
            n.op = Op::Match(p.clone());
            n.ty = Type::Truth;
        }
        if let Op::Constant(Value::PatternSet(patterns)) = &n.op {
            require_unit(e.at, u, Unit::Alignment)?;
            let mut tree = Node {
                at: e.at,
                ty: Type::Truth,
                op: Op::Constant(Value::Truth(Truth::False)),
            };
            if patterns.len() > 64 {
                bail!("pattern set exceeds 64 predicates");
            }
            for p in patterns {
                tree = Node {
                    at: e.at,
                    ty: Type::Truth,
                    op: Op::Binary(
                        "|".into(),
                        Box::new(tree),
                        Box::new(Node {
                            at: e.at,
                            ty: Type::Truth,
                            op: Op::Match(p.clone()),
                        }),
                    ),
                };
            }
            n = tree;
        }
        if n.ty != Type::Truth {
            bail!("byte {}: expected Truth predicate, found {:?}", e.at, n.ty);
        }
        Ok(n)
    }
    pub fn node(&mut self, e: &Expr, u: Unit, locals: &BTreeMap<String, Node>) -> Result<Node> {
        if self.depth >= 128 {
            bail!("byte {}: expanded expression depth exceeds 128", e.at);
        }
        self.depth += 1;
        let result = self.node_inner(e, u, locals);
        self.depth -= 1;
        result
    }
    fn node_inner(&mut self, e: &Expr, u: Unit, locals: &BTreeMap<String, Node>) -> Result<Node> {
        self.nodes += 1;
        if self.nodes > 8192 || self.active.len() > 32 {
            bail!("expanded expression budget exceeded");
        }
        let at = e.at;
        let (ty, op) = match &e.kind {
            Kind::Bool(v) => (
                Type::Truth,
                Op::Constant(Value::Truth(Truth::from_bool(*v))),
            ),
            Kind::Unknown => (Type::Truth, Op::Constant(Value::Truth(Truth::Unknown))),
            Kind::String(v) => (Type::String, Op::Constant(Value::String(v.clone()))),
            Kind::Number(v) => (Type::Number, Op::Constant(Value::Number(*v))),
            Kind::Count(v) => (Type::Count, Op::Constant(Value::Count(*v))),
            Kind::Bases(v) => (Type::Bases, Op::Constant(Value::Bases(*v))),
            Kind::Null => (Type::Null, Op::Constant(Value::Null)),
            Kind::Locus {
                junction,
                chrom,
                start,
                end,
                reverse,
            } => {
                if *junction {
                    (
                        Type::Pattern,
                        Op::Constant(Value::Pattern(Pattern {
                            chrom: chrom.clone(),
                            reverse: *reverse,
                            junctions: vec![(*start, *end)],
                            subsequence: false,
                            left: 0,
                            right: 0,
                        })),
                    )
                } else {
                    (
                        Type::Region,
                        Op::Constant(Value::Region(Region {
                            chrom: chrom.clone(),
                            intervals: vec![(*start, *end)],
                            reverse: *reverse,
                        })),
                    )
                }
            }
            Kind::Name(n) => {
                if let Some(value) = locals.get(n) {
                    return Ok(value.clone());
                }
                if let Some(expr) = self.definitions.get(n).cloned() {
                    if !self.active.insert(n.clone()) {
                        bail!("byte {at}: recursive definition {n}");
                    }
                    let result = self.node(&expr, u, locals);
                    self.active.remove(n);
                    return result;
                }
                (
                    self.fields.get(n).cloned().with_context(|| {
                        format!("byte {at}: unknown name or undeclared metadata column {n}")
                    })?,
                    Op::Field(n.clone()),
                )
            }
            Kind::Field(base, field) => {
                let name = super::parser::field_name(e).unwrap_or_default();
                if !name.starts_with("reads.") && !self.fields.contains_key(&name) {
                    let base = self.node(base, u, locals)?;
                    if let Op::Constant(Value::Feature(feature)) = base.op {
                        let Feature {
                            span,
                            exons,
                            junctions,
                            junction_path,
                            ..
                        } = *feature;
                        let (ty, value) = match field.as_str() {
                            "span" => (Type::Region, Value::Region(span)),
                            "exons" => (Type::Region, Value::Region(exons)),
                            "junctions" => (Type::PatternSet, Value::PatternSet(junctions)),
                            "junction_path" => (
                                Type::Pattern,
                                Value::Pattern(junction_path.context(
                                    "junction_path requires a spliced transcript, not a gene union",
                                )?),
                            ),
                            _ => bail!("unknown feature projection {field}"),
                        };
                        return Ok(Node {
                            at,
                            ty,
                            op: Op::Constant(value),
                        });
                    }
                    if base.ty == Type::Bounds
                        && ["lower", "upper", "status"].contains(&field.as_str())
                    {
                        return Ok(Node {
                            at,
                            ty: if field == "status" {
                                Type::String
                            } else {
                                Type::Count
                            },
                            op: Op::Call(format!("bounds_{field}"), vec![base]),
                        });
                    }
                    bail!("byte {at}: undeclared metadata column or invalid projection {name}");
                }
                let ty = if name.starts_with("reads.") {
                    if !["reads.unique", "reads.multimapping", "reads.total"]
                        .contains(&name.as_str())
                    {
                        bail!("unknown multiplicity {name}");
                    }
                    if ![Unit::Record, Unit::Class, Unit::Cell].contains(&u) {
                        bail!("reads totals require record, class, or cell scope");
                    }
                    Type::Count
                } else {
                    self.fields.get(&name).cloned().with_context(||format!("byte {at}: undeclared metadata column {name}; declare optional columns explicitly"))?
                };
                (ty, Op::Field(name))
            }
            Kind::List(items) => {
                let mut values = Vec::new();
                let mut ty = None;
                for item in items {
                    let n = self.node(item, u, locals)?;
                    if ty.as_ref().is_some_and(|t| t != &n.ty)
                        || matches!(
                            n.ty,
                            Type::Null
                                | Type::List(_)
                                | Type::Region
                                | Type::Pattern
                                | Type::Bounds
                        )
                    {
                        bail!("byte {at}: list elements must have one non-null scalar type");
                    }
                    ty = Some(n.ty);
                    let Op::Constant(v) = n.op else {
                        bail!("v1 list elements must be literals");
                    };
                    values.push(v);
                }
                (
                    Type::List(Box::new(ty.unwrap_or(Type::Null))),
                    Op::Constant(Value::List(values)),
                )
            }
            Kind::Unary(op, a) => {
                let a = if op == "!" {
                    self.predicate(a, u, locals)?
                } else {
                    self.node(a, u, locals)?
                };
                if op == "-" && a.ty != Type::Number {
                    bail!("byte {at}: unary minus requires Number");
                }
                (a.ty.clone(), Op::Unary(op.clone(), Box::new(a)))
            }
            Kind::Binary(op, _, _) if op == ">>" || op == "~>" => {
                return self.path(e, None, locals)
            }
            Kind::Binary(op, a, b) => {
                let boolean = op == "&" || op == "|";
                let a = if boolean {
                    self.predicate(a, u, locals)?
                } else {
                    self.node(a, u, locals)?
                };
                let b = if boolean {
                    self.predicate(b, u, locals)?
                } else {
                    self.node(b, u, locals)?
                };
                let arithmetic = ["+", "-", "*", "/"].contains(&op.as_str());
                if op == "in" || op == "not in" {
                    if !matches!(&b.ty,Type::List(t) if **t==a.ty || **t==Type::Null) {
                        bail!("byte {at}: membership requires a homogeneous list matching its operand");
                    }
                } else if arithmetic {
                    if a.ty != b.ty || !matches!(a.ty, Type::Number | Type::Count) {
                        bail!("byte {at}: arithmetic needs matching numeric types; no implicit conversions");
                    }
                } else if !boolean
                    && (a.ty != b.ty && a.ty != Type::Null && b.ty != Type::Null
                        || matches!(
                            a.ty,
                            Type::Region | Type::Pattern | Type::List(_) | Type::Bounds
                        ))
                {
                    bail!("byte {at}: comparison requires matching scalar types");
                }
                (
                    if arithmetic {
                        a.ty.clone()
                    } else {
                        Type::Truth
                    },
                    Op::Binary(op.clone(), Box::new(a), Box::new(b)),
                )
            }
            Kind::Quant { all, domain, body } => {
                let inner = domain_unit(at, domain, u)?;
                self.stored |= domain == "stored";
                (
                    Type::Truth,
                    Op::Quant {
                        all: *all,
                        domain: domain.clone(),
                        body: Box::new(self.predicate(body, inner, locals)?),
                    },
                )
            }
            Kind::Path { order, body } => return self.path(body, Some(order), locals),
            Kind::Call(name, args) => {
                if let Some(f) = self.doc.functions.iter().find(|f| &f.name == name).cloned() {
                    if args.len() != f.parameters.len() || args.iter().any(|(n, _)| n.is_some()) {
                        bail!(
                            "function {name} requires {} positional arguments",
                            f.parameters.len()
                        );
                    }
                    if !self.active.insert(name.clone()) {
                        bail!("recursive function {name}");
                    }
                    let mut bound = BTreeMap::new();
                    for ((param, ty), (_, arg)) in f.parameters.iter().zip(args) {
                        let value = self.node(arg, u, locals)?;
                        if let Some(ty) = ty {
                            self.check_type(ty, &value, u, f.exported)?;
                        }
                        bound.insert(param.clone(), value);
                    }
                    let result = if f
                        .returns
                        .as_ref()
                        .is_some_and(|t| t.starts_with("Predicate<"))
                    {
                        self.predicate(&f.body, u, &bound)?
                    } else {
                        self.node(&f.body, u, &bound)?
                    };
                    if let Some(ty) = &f.returns {
                        self.check_type(ty, &result, u, f.exported)?;
                    }
                    self.active.remove(name);
                    return Ok(result);
                }
                return self.builtin(at, name, args, u, locals);
            }
            Kind::Resource(_) => bail!("byte {at}: resource alias is not an evidence expression"),
        };
        Ok(Node { at, ty, op })
    }
    fn check_type(&self, declared: &str, node: &Node, unit: Unit, exported: bool) -> Result<()> {
        if exported && node.ty == Type::Truth && !declared.starts_with("Predicate<") {
            bail!("exported Truth predicates must declare Predicate<Unit,Assembly,StrandFrame>");
        }
        let (base, qualifiers) = if let Some((base, args)) = declared.split_once('<') {
            (
                base,
                args.strip_suffix('>')
                    .context("malformed type")?
                    .split(',')
                    .collect::<Vec<_>>(),
            )
        } else {
            (declared, Vec::new())
        };
        if base == "Predicate" {
            if node.ty != Type::Truth || qualifiers.len() != 3 {
                bail!("Predicate requires <Unit,Assembly,alignment|transcript|unstranded>");
            }
            if qualifiers[0] != format!("{unit:?}") {
                bail!(
                    "exported predicate requires {}, found {unit:?}",
                    qualifiers[0]
                );
            }
        } else if base != format!("{:?}", node.ty) {
            bail!("declared type {declared} does not match {:?}", node.ty);
        }
        let coordinate = matches!(
            node.ty,
            Type::Region
                | Type::Pattern
                | Type::PatternSet
                | Type::GeneFeature
                | Type::TranscriptFeature
        ) || base == "Predicate";
        if coordinate {
            let args = if base == "Predicate" {
                &qualifiers[1..]
            } else {
                &qualifiers[..]
            };
            if args.len() != 2 {
                if exported {
                    bail!("exported coordinate types require explicit <Assembly,StrandFrame>");
                }
                return Ok(());
            }
            if self.assembly.as_deref() != Some(args[0]) {
                bail!("type assembly {} disagrees with query assembly", args[0]);
            }
            if !["alignment", "transcript", "unstranded"].contains(&args[1]) {
                bail!("unknown strand-frame requirement {}", args[1]);
            }
            if args[1] == "transcript" && self.library.is_none() {
                bail!("transcript-frame function requires bound library rule");
            }
            if args[1] == "unstranded"
                && matches!(&node.op,Op::Constant(Value::Region(r)) if r.reverse.is_some())
            {
                bail!("function requires unstranded region");
            }
        } else if !qualifiers.is_empty() {
            bail!("scalar types do not accept coordinate qualifiers");
        }
        Ok(())
    }
    fn path(
        &mut self,
        e: &Expr,
        order: Option<&str>,
        locals: &BTreeMap<String, Node>,
    ) -> Result<Node> {
        fn collect<'b>(
            e: &'b Expr,
            out: &mut Vec<&'b Expr>,
            mode: &mut Option<String>,
        ) -> Result<()> {
            if let Kind::Binary(op, a, b) = &e.kind {
                if op == ">>" || op == "~>" {
                    if mode.as_ref().is_some_and(|m| m != op) {
                        bail!(
                            "mixing consecutive and subsequence operators needs separate patterns"
                        );
                    }
                    *mode = Some(op.clone());
                    collect(a, out, mode)?;
                    collect(b, out, mode)?;
                    return Ok(());
                }
            }
            out.push(e);
            Ok(())
        }
        let mut items = Vec::new();
        let mut mode = None;
        collect(e, &mut items, &mut mode)?;
        if items.len() > 64 {
            bail!("path exceeds 64 elements");
        }
        let mut patterns = Vec::new();
        for item in items {
            let n = self.node(item, Unit::Alignment, locals)?;
            let Op::Constant(Value::Pattern(p)) = n.op else {
                bail!("path composition requires constant junction patterns");
            };
            patterns.push(p);
        }
        let mut p = patterns.first().context("empty path")?.clone();
        let descending = match order {
            Some("genomic") => false,
            Some("-") => true,
            Some("+") => false,
            None => p.reverse.context(
                "unstranded composition requires path(order: genomic) or an explicit strand",
            )?,
            _ => bail!("invalid path order"),
        };
        for next in &patterns[1..] {
            if next.chrom != p.chrom || next.reverse != p.reverse {
                bail!("path junctions must share chromosome and strand");
            }
            p.junctions.extend_from_slice(&next.junctions);
        }
        if p.junctions.windows(2).any(|w| {
            if descending {
                w[1].1 >= w[0].0
            } else {
                w[0].1 >= w[1].0
            }
        }) {
            bail!("byte {}: path is not in declared {} order; written junctions are never silently reordered",e.at,if descending{"descending"}else{"genomic"});
        }
        if order == Some("-") || order == Some("+") {
            let opposite = self.library.as_deref().context(
                "transcriptional path requires header library = \"same\" or \"opposite\"",
            )? == "opposite";
            let expected = descending ^ opposite;
            if p.reverse.is_some_and(|r| r != expected) {
                bail!("path alignment strand disagrees with transcript frame/library rule");
            }
            p.reverse = Some(expected);
        }
        if descending {
            p.junctions.reverse();
        }
        p.subsequence = mode.as_deref() == Some("~>");
        Ok(Node {
            at: e.at,
            ty: Type::Pattern,
            op: Op::Constant(Value::Pattern(p)),
        })
    }
    fn builtin(
        &mut self,
        at: usize,
        name: &str,
        args: &[(Option<String>, Expr)],
        u: Unit,
        locals: &BTreeMap<String, Node>,
    ) -> Result<Node> {
        if name == "gene" || name == "transcript" {
            if args.len() != 2 || args.iter().any(|(n, _)| n.is_some()) {
                bail!("{name} requires annotation alias and stable name/id");
            }
            let Kind::Resource(alias) = &args[0].1.kind else {
                bail!("{name} first argument must be @annotation");
            };
            let Op::Constant(Value::String(label)) = self.node(&args[1].1, u, locals)?.op else {
                bail!("feature label must be constant string");
            };
            let value = self
                .resources
                .context("feature resolution needs bound project annotation")?
                .feature(alias, &label, name == "transcript", self.library.as_deref())?;
            return Ok(Node {
                at,
                ty: if name == "gene" {
                    Type::GeneFeature
                } else {
                    Type::TranscriptFeature
                },
                op: Op::Constant(value),
            });
        }
        let allowed: &[&str] = match name {
            "fraction" => &["denominator"],
            "overlaps" => &["min", "blocks"],
            "near" => &["left", "right"],
            _ => &[],
        };
        let mut seen = BTreeSet::new();
        for (n, _) in args {
            if let Some(n) = n {
                if !allowed.contains(&n.as_str()) || !seen.insert(n) {
                    bail!("byte {at}: unexpected or repeated {name} argument {n}");
                }
            }
        }
        let positional: Vec<_> = args
            .iter()
            .filter(|(n, _)| n.is_none())
            .map(|(_, e)| e)
            .collect();
        let arity = if name == "count" { 0 } else { 1 };
        let is_aggregate = [
            "count",
            "sum",
            "count_distinct",
            "count_true",
            "count_false",
            "count_unknown",
        ]
        .contains(&name);
        if is_aggregate && !self.aggregate {
            bail!("byte {at}: aggregate {name} is only valid as a summarize field");
        }
        if positional.len() != arity {
            bail!("byte {at}: {name} requires {arity} positional arguments");
        }
        let arg = positional.first().copied();
        let option = |key: &str| {
            args.iter()
                .find(|(n, _)| n.as_deref() == Some(key))
                .map(|(_, e)| e)
        };
        let (ty, op) = match name {
            "fraction" => {
                let numerator = self.node(arg.unwrap(), u, locals)?;
                let denominator = self.node(
                    option("denominator").context("fraction requires explicit denominator:")?,
                    u,
                    locals,
                )?;
                if numerator.ty != denominator.ty
                    || !matches!(numerator.ty, Type::Count | Type::Number)
                {
                    bail!("fraction requires matching numeric numerator and denominator types");
                }
                (
                    Type::Number,
                    Op::Call(name.into(), vec![numerator, denominator]),
                )
            }
            "overlaps" | "start_in" | "end_in" | "terminal" => {
                require_unit(
                    at,
                    u,
                    if name == "terminal" {
                        Unit::Record
                    } else {
                        Unit::Alignment
                    },
                )?;
                let n = self.node(arg.unwrap(), u, locals)?;
                let Op::Constant(Value::Region(mut region)) = n.op else {
                    bail!("{name} requires Region or RegionSet");
                };
                region.union();
                let minimum = if let Some(e) = option("min") {
                    bases(&self.local_constant(e, u, locals)?)?
                } else {
                    1
                };
                if minimum == 0 {
                    bail!("overlap minimum must be positive");
                }
                let total = if let Some(e) = option("blocks") {
                    match &e.kind {
                        Kind::Name(n) if n == "one" => false,
                        Kind::Name(n) if n == "total" => true,
                        _ => bail!("blocks must be one or total"),
                    }
                } else {
                    false
                };
                if name == "terminal" {
                    self.terminal = true;
                    (Type::Truth, Op::Terminal(region))
                } else {
                    (
                        Type::Truth,
                        Op::Geometry {
                            kind: name.into(),
                            region,
                            minimum,
                            total,
                        },
                    )
                }
            }
            "near" => {
                let n = self.node(arg.unwrap(), u, locals)?;
                let Op::Constant(Value::Pattern(mut p)) = n.op else {
                    bail!("near requires JunctionPattern");
                };
                if p.junctions.len() != 1 {
                    bail!("near requires a single junction");
                }
                p.left = if let Some(e) = option("left") {
                    bases(&self.local_constant(e, u, locals)?)?
                } else {
                    0
                };
                p.right = if let Some(e) = option("right") {
                    bases(&self.local_constant(e, u, locals)?)?
                } else {
                    0
                };
                (Type::Pattern, Op::Constant(Value::Pattern(p)))
            }
            "tx_strand" => {
                require_unit(at, u, Unit::Alignment)?;
                let opposite = self
                    .library
                    .as_deref()
                    .context("tx_strand requires header library = \"same\" or \"opposite\"")?
                    == "opposite";
                let v = self.local_constant(arg.unwrap(), u, locals)?;
                let reverse = match v {
                    Value::String(s) if s == "+" => false,
                    Value::String(s) if s == "-" => true,
                    _ => bail!("tx_strand requires \"+\" or \"-\""),
                };
                (Type::Truth, Op::TxStrand(reverse ^ opposite))
            }
            "nonempty" => {
                let Kind::Name(domain) = &arg.unwrap().kind else {
                    bail!("nonempty requires a domain name");
                };
                domain_unit(at, domain, u)?;
                (Type::Truth, Op::Nonempty(domain.clone()))
            }
            "support_reads" => {
                require_unit(at, u, Unit::Record)?;
                (
                    Type::Bounds,
                    Op::Support(Box::new(self.predicate(
                        arg.unwrap(),
                        Unit::Alignment,
                        locals,
                    )?)),
                )
            }
            "is_unknown" | "is_null" => {
                let n = if name == "is_unknown" {
                    self.predicate(arg.unwrap(), u, locals)?
                } else {
                    self.node(arg.unwrap(), u, locals)?
                };
                (Type::Truth, Op::Call(name.into(), vec![n]))
            }
            "count" => (Type::Count, Op::Call(name.into(), vec![])),
            "sum" | "count_distinct" | "count_true" | "count_false" | "count_unknown" => {
                self.aggregate = false;
                let n = if name.starts_with("count_") && name != "count_distinct" {
                    self.predicate(arg.unwrap(), u, locals)?
                } else {
                    self.node(arg.unwrap(), u, locals)?
                };
                self.aggregate = true;
                if name == "sum" && !matches!(n.ty, Type::Number | Type::Count) {
                    bail!("sum requires numeric input");
                }
                (
                    if name == "sum" {
                        n.ty.clone()
                    } else {
                        Type::Count
                    },
                    Op::Call(name.into(), vec![n]),
                )
            }
            _ => bail!("byte {at}: unknown function {name}"),
        };
        Ok(Node { at, ty, op })
    }
    pub fn compile(mut self) -> Result<Plan> {
        if let Kind::Call(name, args) = &self.doc.source.kind {
            if name != "junctions"
                || args.len() != 2
                || args[0].0.is_some()
                || args[1].0.as_deref() != Some("within")
            {
                bail!("enumeration source requires junctions(@source, within: Region)");
            }
            let Kind::Resource(alias) = &args[0].1.kind else {
                bail!("junctions requires source alias");
            };
            if self.doc.stages.len() != 1 {
                bail!("v1 junction enumeration requires one support(unit:, by:) stage");
            }
            let ast::Stage::Support(unit, by) = &self.doc.stages[0] else {
                bail!("junction enumeration requires support stage");
            };
            let plural = match unit.as_str() {
                "class" => "classes",
                "record" => "records",
                "cell" => "cells",
                _ => bail!("support unit must be record, class, or cell"),
            };
            let mut doc = self.doc.clone();
            doc.source = Expr {
                at: 0,
                kind: Kind::Field(
                    Box::new(Expr {
                        at: 0,
                        kind: Kind::Resource(alias.clone()),
                    }),
                    plural.into(),
                ),
            };
            doc.stages = vec![
                ast::Stage::Within(args[1].1.clone()),
                ast::Stage::Summarize(
                    vec![(
                        "support".into(),
                        Expr {
                            at: 0,
                            kind: Kind::Call("count".into(), vec![]),
                        },
                    )],
                    by.clone(),
                ),
            ];
            let mut compiler = Compiler::new(&doc, self.fields.clone())?;
            if let Some(r) = self.resources {
                compiler = compiler.with_resources(r);
            }
            let mut plan = compiler.compile()?;
            plan.enumeration = true;
            if !matches!(plan.within, Some(Value::Region(_))) {
                bail!("junction enumeration within requires Region or RegionSet");
            }
            return Ok(plan);
        }
        let mut assembly = None;
        let mut version = false;
        for (name, e) in &self.doc.header {
            match (name.as_str(), self.constant(e)?) {
                ("gq", Value::Number(1.0)) => version = true,
                ("assembly", Value::String(s)) => assembly = Some(s),
                ("library", Value::String(s)) if s == "same" || s == "opposite" => {
                    self.library = Some(s)
                }
                _ => bail!("unsupported header field/value {name}"),
            }
        }
        if !version {
            bail!("header must declare gq = 1");
        }
        self.assembly = assembly.clone();
        let Kind::Field(base, unit) = &self.doc.source.kind else {
            bail!("source must be @name.records, .classes, or .cells");
        };
        let Kind::Resource(source) = &base.kind else {
            bail!("source requires resource alias");
        };
        let unit = match unit.as_str() {
            "records" => Unit::Record,
            "classes" => Unit::Class,
            "cells" => Unit::Cell,
            _ => bail!("unknown source unit {unit}"),
        };
        let mut within = None;
        let mut has_within = false;
        let mut stages = Vec::new();
        let mut reduced = false;
        let mut selected = false;
        let mut sorted = false;
        for (i, stage) in self.doc.stages.iter().enumerate() {
            if matches!(stage, ast::Stage::Take(_)) && i + 1 != self.doc.stages.len() {
                bail!("take must be the final pipeline stage");
            }
            let empty = BTreeMap::new();
            stages.push(match stage {
                ast::Stage::Support(..)=>bail!("support requires junctions source"),
                ast::Stage::Within(e)=>{
                    if i!=0||has_within {bail!("within must be the first pipeline stage, exactly once");} has_within=true;
                    if !matches!(&e.kind,Kind::Name(n) if n=="all") {
                        let v=self.constant(e)?; if !matches!(v,Value::Region(_)|Value::Pattern(_)){bail!("within requires region, junction, or all");} within=Some(v);
                    } continue;
                }
                ast::Stage::Where(e)=> {if reduced||selected{bail!("where must precede aggregation/projection");} Stage::Where(self.predicate(e,unit,&empty)?)},
                ast::Stage::Derive(fields)=> {
                    if reduced||selected{bail!("derive must precede aggregation/projection");}
                    let mut out=Vec::new(); for (name,e) in fields { if self.fields.contains_key(name)||self.definitions.contains_key(name){bail!("derive cannot shadow {name}");} let n=self.node(e,unit,&empty)?; self.fields.insert(name.clone(),n.ty.clone()); out.push((name.clone(),n)); } Stage::Derive(out)
                }
                ast::Stage::Tally(fields,by)|ast::Stage::Summarize(fields,by)=> {
                    if reduced||selected{bail!("only one aggregation is allowed, before projection");} reduced=true;
                    let tally=matches!(stage,ast::Stage::Tally(..));
                    self.aggregate = !tally;
                    let fields=fields.iter().map(|(k,e)|Ok((k.clone(),if tally{self.predicate(e,unit,&empty)?}else{self.node(e,unit,&empty)?}))).collect::<Result<Vec<_>>>()?;
                    self.aggregate=false;
                    if !tally&&fields.iter().any(|(_,n)|!matches!(&n.op,Op::Call(name,_) if ["count","sum","count_distinct","count_true","count_false","count_unknown"].contains(&name.as_str()))){bail!("summarize fields must be aggregate calls");}
                    let by=by.iter().map(|(k,e)|Ok((k.clone(),self.node(e,unit,&empty)?))).collect::<Result<Vec<_>>>()?;
                    if by.iter().any(|(_,n)|!matches!(&n.op,Op::Field(name) if name=="sample"||name=="cell"||name=="cell.id"||self.resources.is_some_and(|r|r.metadata.columns.contains_key(name)))){bail!("v1 grouping uses base sample/cell metadata so pre-filter denominators remain defined");}
                    if by.iter().any(|(_,n)|!matches!(n.ty,Type::String|Type::Number|Type::Count|Type::Truth)){bail!("group keys must be scalars");}
                    let mut schema:BTreeMap<String,Type>=by.iter().map(|(k,n)|(k.clone(),n.ty.clone())).collect();
                    for (name,n) in &fields{let key=if tally{format!("{name}_state")}else{name.clone()};if schema.insert(key.clone(),n.ty.clone()).is_some(){bail!("duplicate output column {key}");}}
                    if tally&&schema.insert("count".into(),Type::Count).is_some(){bail!("group key conflicts with tally count");}
                    self.fields=schema;
                    if tally{Stage::Tally(fields,by)}else{Stage::Summarize(fields,by)}
                }
                ast::Stage::Select(fields)=>{
                    if selected{bail!("only one select stage is allowed");}selected=true;
                    let fields=fields.iter().map(|(k,e)|Ok((k.clone(),self.node(e,unit,&empty)?))).collect::<Result<Vec<_>>>()?;
                    if fields.iter().any(|(_,n)|matches!(n.ty,Type::Region|Type::Pattern|Type::GeneFeature|Type::TranscriptFeature|Type::PatternSet)){bail!("select emits scalar values and support bounds, not unevaluated patterns/features");}
                    self.fields=fields.iter().map(|(k,n)|(k.clone(),n.ty.clone())).collect();Stage::Select(fields)
                }
                ast::Stage::Sort(fields)=>{if !reduced&&!selected{bail!("sort follows aggregation or select");}sorted=true;let fields=fields.iter().map(|(k,e)|Ok((k.clone(),self.node(e,unit,&empty)?))).collect::<Result<Vec<_>>>()?;if fields.iter().any(|(name,n)|!matches!(&n.op,Op::Field(field) if field==name)){bail!("sort requires output field names, without aliases");}Stage::Sort(fields)},
                ast::Stage::Take(n)=>{if !sorted{bail!("take requires explicit deterministic sort");}Stage::Take(*n)},
            });
        }
        if !has_within {
            bail!("explicit within region/junction/all is required");
        }
        if !reduced && !stages.iter().any(|s| matches!(s, Stage::Select(_))) {
            bail!("query must end in tally, summarize, or select");
        }
        Ok(Plan {
            enumeration: false,
            source: source.clone(),
            unit,
            assembly: assembly.context("header assembly = \"...\" is required")?,
            library: self.library,
            within,
            stages,
            terminal: self.terminal,
            diagnostic_stored: self.stored,
            fields: self.fields,
        })
    }
}

fn validate_definitions(doc: &ast::Document) -> Result<()> {
    fn references(e: &Expr, out: &mut BTreeSet<String>) {
        match &e.kind {
            Kind::Name(name) => {
                out.insert(name.clone());
            }
            Kind::Call(name, args) => {
                out.insert(name.clone());
                for (_, e) in args {
                    references(e, out);
                }
            }
            Kind::Binary(_, a, b) => {
                references(a, out);
                references(b, out);
            }
            Kind::Unary(_, a)
            | Kind::Field(a, _)
            | Kind::Quant { body: a, .. }
            | Kind::Path { body: a, .. } => references(a, out),
            Kind::List(values) => {
                for e in values {
                    references(e, out);
                }
            }
            _ => {}
        }
    }
    let mut graph = BTreeMap::new();
    for (name, e) in &doc.bindings {
        let mut refs = BTreeSet::new();
        references(e, &mut refs);
        graph.insert(name.clone(), refs);
    }
    for f in &doc.functions {
        let mut params = BTreeSet::new();
        for (name, _) in &f.parameters {
            if !params.insert(name) {
                bail!("duplicate parameter {name} in {}", f.name);
            }
        }
        let mut refs = BTreeSet::new();
        references(&f.body, &mut refs);
        for name in params {
            refs.remove(name);
        }
        graph.insert(f.name.clone(), refs);
        if f.exported {
            for declared in f
                .parameters
                .iter()
                .filter_map(|(_, t)| t.as_deref())
                .chain(f.returns.as_deref())
            {
                if [
                    "Region",
                    "Pattern",
                    "PatternSet",
                    "GeneFeature",
                    "TranscriptFeature",
                    "Truth",
                    "Predicate",
                ]
                .contains(&declared)
                {
                    bail!("exported scientific type {declared} requires explicit unit/assembly/strand qualifiers");
                }
            }
        }
    }
    fn visit(
        name: &str,
        graph: &BTreeMap<String, BTreeSet<String>>,
        active: &mut BTreeSet<String>,
        done: &mut BTreeSet<String>,
    ) -> Result<()> {
        if done.contains(name) || !graph.contains_key(name) {
            return Ok(());
        }
        if !active.insert(name.into()) {
            bail!("recursive definition {name}");
        }
        if active.len() > 32 {
            bail!("definition dependency depth exceeds 32");
        }
        for next in &graph[name] {
            visit(next, graph, active, done)?;
        }
        active.remove(name);
        done.insert(name.into());
        Ok(())
    }
    let mut done = BTreeSet::new();
    for name in graph.keys() {
        visit(name, &graph, &mut BTreeSet::new(), &mut done)?;
    }
    Ok(())
}
fn require_unit(at: usize, actual: Unit, expected: Unit) -> Result<()> {
    if actual != expected {
        bail!("byte {at}: predicate requires {expected:?}; found {actual:?}; use an explicit quantifier (terminal linkage is Record only)");
    }
    Ok(())
}
fn domain_unit(at: usize, domain: &str, u: Unit) -> Result<Unit> {
    match domain {
        "record" if u == Unit::Class || u == Unit::Cell => Ok(Unit::Record),
        "unique" | "stored" if u == Unit::Record => Ok(Unit::Alignment),
        "multimap" if u == Unit::Record => Ok(Unit::Signature),
        "alternative" if u == Unit::Signature => Ok(Unit::Alignment),
        _ => bail!("byte {at}: {domain} is not a domain of {u:?}"),
    }
}
fn bases(v: &Value) -> Result<u32> {
    match v {
        Value::Bases(v) => Ok(*v),
        Value::Number(v) if *v >= 0.0 && *v <= u32::MAX as f64 && v.fract() == 0.0 => Ok(*v as u32),
        _ => bail!("expected nonnegative integral bp quantity"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn compile(q: &str) -> Result<Plan> {
        let d = super::super::parser::parse(q)?;
        Compiler::new(&d, BTreeMap::from([("sample".into(), Type::String)]))?.compile()
    }
    #[test]
    fn scope_and_legacy_diagnostics() {
        for expression in [
            "terminal(g[1:10..20])",
            "any unique { terminal(g[1:10..20]) }",
            "j[1:10..20] & true",
        ] {
            let unit = if expression.starts_with("terminal") {
                "cells"
            } else {
                "records"
            };
            assert!(compile(&format!("header {{gq=1,assembly=\"test\"}} from @x.{unit} |> within all |> derive {{a={expression}}} |> tally {{a}}")).is_err());
        }
    }
    #[test]
    fn explicit_domains_compile() {
        assert!(compile("header {gq=1,assembly=\"test\"} from @x.records |> within all |> derive {a=any unique{j[1:10..20]},b=any multimap{all alternative{j[1:10..20]}}} |> tally {a,b}").is_ok());
    }
    #[test]
    fn path_order_and_recursion() {
        assert!(compile("header {gq=1,assembly=\"test\"} let A=A from @x.records |> within all |> derive {a=any unique{A}} |> tally {a}").is_err());
        assert!(compile("header {gq=1,assembly=\"test\"} let A=j[1:10..20] >> j[1:30..40] from @x.records |> within all |> derive {a=any unique{A}} |> tally {a}").is_err());
    }
}
