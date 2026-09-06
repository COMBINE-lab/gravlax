use super::{
    eval::*,
    model::*,
    physical::{self, Kernel, Outcome, States},
    resources::{Member, Resources},
    Engine, Input,
};
use crate::archivecmd::{decode_chunk, read_chunk_index, ChunkInfo, LazyArchive};
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::time::Instant;

type Row = BTreeMap<String, Value>;
pub struct Table {
    pub columns: Vec<(String, Type)>,
    pub rows: Vec<Row>,
    pub summary: serde_json::Value,
}
#[derive(Default)]
struct Accumulator {
    key: Row,
    values: Vec<Value>,
    distinct: Vec<HashSet<String>>,
}
struct PartialAggregate {
    key: String,
    values: Vec<Value>,
    distinct: Vec<HashSet<String>>,
}
#[derive(Default)]
struct Totals {
    chunks: usize,
    records: usize,
    steps: u64,
    unknown: u64,
    where_unknown: u64,
    population: u64,
    source_cells: u64,
    source_units: u64,
    timings: BTreeMap<&'static str, f64>,
}
impl Totals {
    fn timed(&mut self, name: &'static str, start: Option<Instant>) {
        if let Some(start) = start {
            *self.timings.entry(name).or_default() += start.elapsed().as_secs_f64() * 1000.0;
        }
    }
}

fn routes(
    plan: &Plan,
    la: &mut LazyArchive,
    chunks: &[ChunkInfo],
    allow: bool,
) -> Result<(Vec<usize>, String)> {
    if plan.enumeration {
        if !allow {
            bail!("v1 exact junction enumeration requires --allow-full-scan; aligned-block postings do not route intron interiors");
        }
        return Ok((
            (0..chunks.len()).collect(),
            "explicit full scan for complete junction enumeration and unit support".into(),
        ));
    }
    if let Some(value) = &plan.within {
        let (chrom, intervals, junction) = match value {
            Value::Region(r) => (&r.chrom, r.intervals.clone(), false),
            Value::Pattern(p) if p.junctions.len() == 1 && p.left == 0 && p.right == 0 => {
                (&p.chrom, p.junctions.clone(), true)
            }
            _ => bail!("within requires a RegionSet or one exact junction"),
        };
        let chrom = la
            .chrom_names
            .iter()
            .position(|s| s == chrom)
            .with_context(|| {
                format!("source lacks chromosome {chrom}; no implicit contig aliases")
            })? as u32;
        if let Some(index) = la.access_index()? {
            let mut routes = BTreeSet::new();
            for (start, end) in intervals {
                routes.extend(index.geometry_chunks(junction, 2, chrom, start, end));
            }
            return Ok((
                routes.into_iter().collect(),
                "access-index retained-geometry routes".into(),
            ));
        }
    }
    if !allow {
        bail!("initial witnessed-universe routing requires --allow-full-scan: no complete access index (or within all)");
    }
    Ok(((0..chunks.len()).collect(), "explicit full scan".into()))
}
pub fn explain(plan: &Plan, resources: &Resources, input: &Input) -> Result<serde_json::Value> {
    let mut sources = Vec::new();
    let mut roots = BTreeSet::new();
    for member in resources.members(&plan.source)? {
        let mut la = LazyArchive::open(&member.path)?;
        let root = la
            .reader()
            .content_commitment()
            .context("GQ requires authenticated v2 archives")?
            .to_hex();
        if !roots.insert(root.clone()) {
            bail!("duplicate archive content under different samples");
        }
        if plan.terminal && la.terminal_tail_capability().is_none() {
            bail!(
                "sample {} lacks required terminal-tail capability",
                member.sample
            );
        }
        let chunks = read_chunk_index(la.reader())?;
        check_chromosomes(plan, &la.chrom_names)?;
        let (initial, path) = routes(plan, &mut la, &chunks, input.allow_full_scan)?;
        if initial.len() > input.max_chunks {
            bail!("initial routes exceed chunk budget");
        }
        let closure = match plan.unit {
            Unit::Record => "none: complete selected records",
            Unit::Class if la.access_index()?.is_some() => "data-dependent class index closure",
            _ if input.allow_full_scan => "explicit full scan for complete unit children",
            _ => bail!(
                "{:?} closure needs --allow-full-scan on this archive",
                plan.unit
            ),
        };
        let denominator = if needs_denominator_scan(plan) {
            if !input.allow_full_scan {
                bail!("exact source-unit denominators by cell metadata require --allow-full-scan");
            }
            "explicit full scan for per-group source-unit totals"
        } else {
            "archive metadata plus cell dictionary"
        };
        sources.push(json!({"sample":member.sample,"archive":member.path,"content_root":root,"initial_chunks":initial.len(),"initial_routing":path,"closure":closure,"denominators":denominator,"total_chunks":chunks.len()}));
    }
    Ok(
        json!({"schema":"gravlax.gq.explain.v1","logical_plan":plan,"physical_strategy":physical::description(plan, input.engine == Engine::Reference, input.parallel_decode),"sources":sources,"execution_policy":{"allow_full_scan":input.allow_full_scan,"max_chunks":input.max_chunks,"max_records":input.max_records,"max_steps":input.max_steps,"max_rows":input.max_rows},"semantics":{"universe":"positive retained witness across unique representatives and grouped multimapper alternatives; children are not cropped","all":"classical; true on empty domains","unknown":"Kleene logic; omitted ends are not bounded by retained ends","predicate_effects":{"junctions_and_paths":"chain-invariant","start_in":"extent-sensitive; sound start-interval proof","overlaps_and_end_in":"extent-sensitive; conservative when omitted geometry is unresolved"},"stored":"diagnostic literal representatives","support_reads":"accepted-observation bounds; unique weights plus each multimapping signature once","sum_reads_after_where":"sum selected units' total accepted observations, not reads satisfying the filter"}}),
    )
}
fn needs_denominator_scan(plan: &Plan) -> bool {
    !plan.enumeration
        && plan.unit != Unit::Cell
        && plan.stages.iter().any(|s| match s {
            Stage::Tally(_, by) | Stage::Summarize(_, by) => by
                .iter()
                .any(|(_, n)| !matches!(&n.op,Op::Field(name) if name=="sample")),
            _ => false,
        })
}
fn check_chromosomes(plan: &Plan, chroms: &[String]) -> Result<()> {
    let mut stack = Vec::new();
    for stage in &plan.stages {
        match stage {
            Stage::Where(n) => stack.push(n),
            Stage::Derive(fields) | Stage::Select(fields) | Stage::Sort(fields) => {
                stack.extend(fields.iter().map(|(_, n)| n))
            }
            Stage::Tally(fields, by) | Stage::Summarize(fields, by) => {
                stack.extend(fields.iter().chain(by).map(|(_, n)| n))
            }
            _ => {}
        }
    }
    while let Some(n) = stack.pop() {
        let chrom = match &n.op {
            Op::Geometry { region, .. } | Op::Terminal(region) => Some(&region.chrom),
            Op::Match(p) => Some(&p.chrom),
            Op::Unary(_, a) | Op::Quant { body: a, .. } | Op::Support(a) => {
                stack.push(a);
                None
            }
            Op::Binary(_, a, b) => {
                stack.push(a);
                stack.push(b);
                None
            }
            Op::Call(_, args) => {
                stack.extend(args);
                None
            }
            _ => None,
        };
        if chrom.is_some_and(|c| !chroms.contains(c)) {
            bail!(
                "byte {}: predicate chromosome {:?} is absent from archive",
                n.at,
                chrom.unwrap()
            );
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn load_chunks(
    la: &mut LazyArchive,
    chunks: &[ChunkInfo],
    selected: &[usize],
    seen: &mut BTreeSet<usize>,
    records: &mut Vec<Evidence>,
    totals: &mut Totals,
    input: &Input,
    terminal: bool,
) -> Result<()> {
    let start = input.profile.then(Instant::now);
    let mut bases = Vec::with_capacity(chunks.len());
    let mut base = 0u64;
    for c in chunks {
        bases.push(base);
        base += u64::from(c.n_mols);
    }
    let tails = if terminal {
        la.terminal_tail_routes()?
            .context("missing terminal stream")?
            .to_vec()
    } else {
        Vec::new()
    };
    if tails.iter().map(|t| u64::from(t.events)).sum::<u64>() > input.max_terminal_events {
        bail!("terminal-event budget exceeded");
    }
    let mut pending = Vec::new();
    for &i in selected {
        if !seen.insert(i) {
            continue;
        }
        totals.chunks += 1;
        totals.records = totals
            .records
            .checked_add(chunks[i].n_mols as usize)
            .context("record count overflow")?;
        if totals.chunks > input.max_chunks || totals.records > input.max_records {
            bail!("decoded chunk/record budget exceeded; no partial result published");
        }
        pending.push(i);
    }
    // Validate the entire pending route against budgets before speculative I/O.
    // Indexed parallel collection retains route order; cell/terminal attachment and
    // aggregation stay ordered, including floating point sums across federations.
    let window = if input.engine == Engine::Auto && input.parallel_decode {
        rayon::current_num_threads().saturating_mul(2).max(1)
    } else {
        1
    };
    for selected_window in pending.chunks(window) {
        let decoded: Vec<_> = {
            let (reader, tables) = la.reader_and_tables();
            let reader = &*reader;
            let decode = |&i: &usize| -> Result<_> {
                let (compressed, raw_len) = reader.read_compressed_at(&format!("c{i}"))?;
                let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                Ok((i, decode_chunk(&raw, &chunks[i], None, tables)?))
            };
            if input.engine == Engine::Auto && selected_window.len() > 1 {
                selected_window
                    .par_iter()
                    .map(decode)
                    .collect::<Result<_>>()?
            } else {
                selected_window.iter().map(decode).collect::<Result<_>>()?
            }
        };
        for (i, mut molecules) in decoded {
            la.prefetch_coc(molecules.iter().map(|m| m.umi_class))?;
            for m in &mut molecules {
                m.cell = la.cell_of_cached(m.umi_class)?;
            }
            let mut events: BTreeMap<usize, Vec<u32>> = BTreeMap::new();
            if let Some(route) = tails.iter().find(|r| r.chunk as usize == i) {
                for e in la.terminal_tail_records(*route, &chunks[i], bases[i], &molecules)? {
                    events
                        .entry(e.local_molecule_ordinal as usize)
                        .or_default()
                        .push(e.anchor);
                }
            }
            for (j, molecule) in molecules.into_iter().enumerate() {
                records.push(Evidence {
                    molecule,
                    tails: events.remove(&j).unwrap_or_default(),
                    ordinal: bases[i] + j as u64,
                });
            }
        }
    }
    totals.timed("chunk_read_decode_attach", start);
    Ok(())
}
fn witness(evaluator: &mut Evaluator<'_>, i: usize, within: Option<&Node>) -> Result<bool> {
    let Some(node) = within else { return Ok(true) };
    let m = &evaluator.data.records[i].molecule;
    let mut alignments = smallvec::SmallVec::<[Alignment; 8]>::new();
    for c in &m.chains {
        for &(position, shape) in &c.reps {
            alignments.push(Alignment {
                chrom: m.chrom,
                reverse: m.strand_rev,
                position,
                shape,
            });
        }
    }
    for j in 0..m.mms.len() {
        alignments.extend(evaluator.alternatives(m, j)?);
    }
    for alignment in alignments {
        if evaluator
            .eval(
                node,
                Scope {
                    records: &[],
                    record: Some(i),
                    alignment: Some(alignment),
                    signature: None,
                    omitted: None,
                },
                &BTreeMap::new(),
            )?
            .truth()?
            == Truth::True
        {
            return Ok(true);
        }
    }
    Ok(false)
}
fn key_value(v: &Value) -> Result<String> {
    Ok(serde_json::to_string(&v.json())?)
}
fn key(row: &Row) -> Result<String> {
    use std::fmt::Write;
    if row.is_empty() {
        return Ok("{}".into());
    }
    let mut key = String::new();
    for (name, value) in row {
        write!(key, "{}:{name}", name.len())?;
        match value {
            Value::Truth(t) => key.push(match t {
                Truth::True => 'T',
                Truth::False => 'F',
                Truth::Unknown => 'U',
            }),
            Value::Null => key.push('N'),
            Value::String(s) => write!(key, "S{}:{s}", s.len())?,
            Value::Number(n) => write!(key, "D{:016x}", if *n == 0.0 { 0 } else { n.to_bits() })?,
            Value::Count(n) => write!(key, "C{n};")?,
            _ => {
                let json = serde_json::to_string(&value.json())?;
                write!(key, "J{}:{json}", json.len())?;
            }
        }
    }
    Ok(key)
}
fn evaluate_fields(
    fields: &[(String, Node)],
    ev: &mut Evaluator<'_>,
    scope: Scope<'_>,
    row: &Row,
) -> Result<Row> {
    fields
        .iter()
        .map(|(k, n)| Ok((k.clone(), ev.eval(n, scope, row)?)))
        .collect()
}

pub fn execute(plan: &Plan, resources: &Resources, input: &Input) -> Result<Table> {
    if plan.enumeration {
        return enumerate(plan, resources, input);
    }
    let start = input.profile.then(Instant::now);
    let explanation = explain(plan, resources, input)?;
    let mut totals = Totals::default();
    totals.timed("preflight_planning", start);
    let mut output = Vec::new();
    let mut accumulators: BTreeMap<String, Accumulator> = BTreeMap::new();
    let mut populations: BTreeMap<String, (Row, u64, u64)> = BTreeMap::new();
    let mut access = Vec::new();
    let mut columns = Vec::new();
    let mut state_fields = Vec::new();
    for stage in &plan.stages {
        match stage {
            Stage::Tally(fields, by) => {
                columns = by.iter().map(|(k, n)| (k.clone(), n.ty.clone())).collect();
                state_fields = fields.iter().map(|(k, _)| format!("{k}_state")).collect();
                columns.extend(state_fields.iter().map(|k| (k.clone(), Type::Truth)));
                columns.push(("count".into(), Type::Count));
            }
            Stage::Summarize(fields, by) => {
                columns = by
                    .iter()
                    .chain(fields)
                    .map(|(k, n)| (k.clone(), n.ty.clone()))
                    .collect();
                if by.is_empty() {
                    accumulators.insert(
                        "{}".into(),
                        Accumulator {
                            key: BTreeMap::new(),
                            values: fields
                                .iter()
                                .map(|(_, n)| {
                                    if n.ty == Type::Number {
                                        Value::Number(0.0)
                                    } else {
                                        Value::Count(0)
                                    }
                                })
                                .collect(),
                            distinct: vec![HashSet::new(); fields.len()],
                        },
                    );
                }
            }
            Stage::Select(fields) => {
                columns = fields
                    .iter()
                    .map(|(k, n)| (k.clone(), n.ty.clone()))
                    .collect()
            }
            _ => {}
        }
    }
    for member in resources.members(&plan.source)? {
        execute_member(
            plan,
            resources,
            input,
            member,
            &mut totals,
            &mut output,
            &mut accumulators,
            &mut populations,
            &mut access,
            &mut columns,
            &explanation,
        )?;
    }
    let start = input.profile.then(Instant::now);
    for accumulator in accumulators.into_values() {
        let mut row = accumulator.key;
        let values = accumulator.values;
        let aggregations = plan
            .stages
            .iter()
            .find_map(|s| match s {
                Stage::Tally(fields, _) | Stage::Summarize(fields, _) => {
                    Some((matches!(s, Stage::Tally(..)), fields))
                }
                _ => None,
            })
            .context("missing aggregation")?;
        if aggregations.0 {
            row.insert("count".into(), values[0].clone());
        } else {
            for (i, (name, n)) in aggregations.1.iter().enumerate() {
                row.insert(
                    name.clone(),
                    if matches!(&n.op,Op::Call(name,_) if name=="count_distinct") {
                        Value::Count(accumulator.distinct[i].len() as u64)
                    } else {
                        values[i].clone()
                    },
                );
            }
        }
        output.push(row);
    }
    let mut available = output.len();
    let mut sorted = false;
    let mut aggregated = false;
    for (stage_index, stage) in plan.stages.iter().enumerate() {
        match stage {
            Stage::Tally(..) | Stage::Summarize(..) => aggregated = true,
            Stage::Select(fields) if aggregated => {
                let mut ev = Evaluator {
                    data: Data {
                        records: &[],
                        shapes: &[],
                        patterns: &[],
                        chroms: &[],
                    },
                    steps: totals.steps,
                    max_steps: input.max_steps,
                    unknown_omitted: 0,
                };
                let scope = Scope {
                    records: &[],
                    record: None,
                    alignment: None,
                    signature: None,
                    omitted: None,
                };
                output = output
                    .iter()
                    .map(|row| evaluate_fields(fields, &mut ev, scope, row))
                    .collect::<Result<_>>()?;
                columns = fields
                    .iter()
                    .map(|(k, n)| (k.clone(), n.ty.clone()))
                    .collect();
                totals.steps = ev.steps;
            }
            Stage::Sort(fields) => {
                // Stable scalar ordering, with null before present values. Complete key tie-break
                // makes pagination deterministic even when the requested key is non-unique.
                let compare = |a: &Row, b: &Row| {
                    for (name, _) in fields {
                        let order = compare_value(a.get(name), b.get(name));
                        if !order.is_eq() {
                            return order;
                        }
                    }
                    key(a).unwrap_or_default().cmp(&key(b).unwrap_or_default())
                };
                if let Some(Stage::Take(n)) = plan
                    .stages
                    .get(stage_index + 1)
                    .filter(|_| input.engine == Engine::Auto)
                {
                    if *n < output.len() {
                        // The original ordinal is the last tie-break, preserving stable-sort
                        // behavior even for semantically equal rows (e.g. signed zero).
                        let mut ranked: Vec<_> = output.into_iter().enumerate().collect();
                        let order = |a: &(usize, Row), b: &(usize, Row)| {
                            compare(&a.1, &b.1).then_with(|| a.0.cmp(&b.0))
                        };
                        if *n > 0 {
                            ranked.select_nth_unstable_by(*n, order);
                            ranked[..*n].sort_unstable_by(order);
                        }
                        output = ranked.into_iter().map(|(_, row)| row).collect();
                    } else {
                        output.sort_by(compare);
                    }
                } else {
                    output.sort_by(compare);
                }
                sorted = true;
            }
            Stage::Take(n) => {
                if !sorted {
                    bail!("take requires sort");
                }
                available = output.len();
                output.truncate(*n);
            }
            _ => {}
        }
    }
    if output.len() > input.max_rows {
        bail!("result-row budget exceeded");
    }
    let populations:Vec<_>=populations.into_values().map(|(k,n,source)|json!({"group":k.iter().map(|(k,v)|(k,v.json())).collect::<BTreeMap<_,_>>(),"population_units":n,"source_scope_units":source})).collect();
    let truncated = available > output.len();
    totals.timed("finalize_sort_project", start);
    Ok(Table {
        columns,
        rows: output,
        summary: json!({"language_version":1,"unit":plan.unit,"state_fields":state_fields,"population_units":totals.population,"source_scope_units":totals.source_units,"source_scope_cells":totals.source_cells,"population_by_group":populations,"where_dropped_unknown":totals.where_unknown,"unknown_evaluations_by_cause":{"omitted_unique_geometry":totals.unknown},"decoded_chunks":totals.chunks,"decoded_records":totals.records,"expression_steps":totals.steps,"archive_access":access,"explain":explanation,"available_rows":available,"truncated":truncated,"timings_ms":totals.timings}),
    })
}

fn enumerate(plan: &Plan, resources: &Resources, input: &Input) -> Result<Table> {
    let start = input.profile.then(Instant::now);
    let explanation = explain(plan, resources, input)?;
    let Some(Value::Region(region)) = &plan.within else {
        bail!("junction enumeration requires region");
    };
    let Some(Stage::Summarize(_, by)) = plan.stages.first() else {
        bail!("invalid support plan");
    };
    let mut counts: BTreeMap<String, (Row, u64)> = BTreeMap::new();
    let mut totals = Totals::default();
    totals.timed("preflight_planning", start);
    for member in resources.members(&plan.source)? {
        let start = input.profile.then(Instant::now);
        let mut la = LazyArchive::open(&member.path)?;
        verify_root(&mut la, member, &explanation)?;
        let chunks = read_chunk_index(la.reader())?;
        let shapes = la.shapes()?;
        let chroms = la.chrom_names.clone();
        let cells = la.cells()?.to_vec();
        let mut records = Vec::new();
        let mut seen = BTreeSet::new();
        totals.timed("archive_setup", start);
        load_chunks(
            &mut la,
            &chunks,
            &(0..chunks.len()).collect::<Vec<_>>(),
            &mut seen,
            &mut records,
            &mut totals,
            input,
            false,
        )?;
        let start = input.profile.then(Instant::now);
        let mut units: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        for (i, r) in records.iter().enumerate() {
            units
                .entry(match plan.unit {
                    Unit::Class => u64::from(r.molecule.umi_class),
                    Unit::Cell => u64::from(r.molecule.cell),
                    _ => r.ordinal,
                })
                .or_default()
                .push(i);
        }
        let mut ev = Evaluator {
            data: Data {
                records: &records,
                shapes: &shapes,
                patterns: &[],
                chroms: &chroms,
            },
            steps: totals.steps,
            max_steps: input.max_steps,
            unknown_omitted: 0,
        };
        for indices in units.into_values() {
            let mut junctions = BTreeSet::new();
            for &i in &indices {
                let m = &records[i].molecule;
                if !ev.same(region, m.chrom, m.strand_rev) {
                    continue;
                }
                for chain in &m.chains {
                    let &(position, shape) = chain.reps.first().context("empty chain")?;
                    for pair in ev
                        .blocks(Alignment {
                            chrom: m.chrom,
                            reverse: m.strand_rev,
                            position,
                            shape,
                        })?
                        .windows(2)
                    {
                        ev.steps += 1;
                        if ev.steps > input.max_steps {
                            bail!("junction enumeration work budget exceeded");
                        }
                        let (donor, acceptor) = (pair[0].1, pair[1].0);
                        if region
                            .intervals
                            .iter()
                            .any(|&(a, b)| a <= donor && acceptor <= b)
                        {
                            junctions.insert((donor, acceptor, m.strand_rev));
                        }
                    }
                }
            }
            if junctions.is_empty() {
                continue;
            }
            let cell = records[indices[0]].molecule.cell;
            let barcode = String::from_utf8(
                crate::querycmd::unpack_cell_bytes(cells[cell as usize]).to_vec(),
            )?;
            let fields = resources.row(member, &barcode)?;
            let scope = Scope {
                records: &indices,
                record: None,
                alignment: None,
                signature: None,
                omitted: None,
            };
            let group = evaluate_fields(by, &mut ev, scope, &fields)?;
            for (donor, acceptor, reverse) in junctions {
                let mut row = group.clone();
                row.insert("chrom".into(), Value::String(region.chrom.clone()));
                row.insert("donor".into(), Value::Count(u64::from(donor)));
                row.insert("acceptor".into(), Value::Count(u64::from(acceptor)));
                row.insert(
                    "alignment_strand".into(),
                    Value::String(if reverse { "-" } else { "+" }.into()),
                );
                let key = key(&row)?;
                if !counts.contains_key(&key) && counts.len() >= input.max_rows {
                    bail!("junction candidate/group budget exceeded");
                }
                counts.entry(key).or_insert((row, 0)).1 += 1;
            }
        }
        totals.steps = ev.steps;
        totals.timed("enumerate_and_aggregate", start);
    }
    let start = input.profile.then(Instant::now);
    let rows = counts
        .into_values()
        .map(|(mut row, count)| {
            row.insert("support".into(), Value::Count(count));
            row
        })
        .collect();
    let mut columns: Vec<_> = by.iter().map(|(k, n)| (k.clone(), n.ty.clone())).collect();
    columns.extend([
        ("chrom".into(), Type::String),
        ("donor".into(), Type::Count),
        ("acceptor".into(), Type::Count),
        ("alignment_strand".into(), Type::String),
        ("support".into(), Type::Count),
    ]);
    totals.timed("finalize_sort_project", start);
    Ok(Table {
        columns,
        rows,
        summary: json!({"unit":plan.unit,"domain":"unique","support_status":"exact","candidate_scope":"both intron boundaries contained in one target interval","catalogue_counts":"not emitted; no unverified upper-bound reinterpretation","decoded_chunks":totals.chunks,"decoded_records":totals.records,"expression_steps":totals.steps,"explain":explanation,"timings_ms":totals.timings}),
    })
}

#[allow(clippy::too_many_arguments)]
fn execute_member(
    plan: &Plan,
    resources: &Resources,
    input: &Input,
    member: &Member,
    totals: &mut Totals,
    output: &mut Vec<Row>,
    accumulators: &mut BTreeMap<String, Accumulator>,
    populations: &mut BTreeMap<String, (Row, u64, u64)>,
    access: &mut Vec<serde_json::Value>,
    columns: &mut Vec<(String, Type)>,
    explanation: &serde_json::Value,
) -> Result<()> {
    let start = input.profile.then(Instant::now);
    let mut kernel = if input.engine == Engine::Auto {
        Kernel::compile(plan)
    } else {
        None
    };
    let mut la = LazyArchive::open(&member.path)?;
    verify_root(&mut la, member, explanation)?;
    let chunks = read_chunk_index(la.reader())?;
    let (initial, routing) = routes(plan, &mut la, &chunks, input.allow_full_scan)?;
    let shapes = la.shapes()?;
    let patterns = la.patterns()?;
    let chroms = la.chrom_names.clone();
    let cells = la.cells()?.to_vec();
    totals.source_cells += cells.len() as u64;
    let source_units = match plan.unit {
        Unit::Record => chunks.iter().map(|c| u64::from(c.n_mols)).sum(),
        Unit::Class => u64::from(la.class_count()),
        Unit::Cell => cells.len() as u64,
        _ => unreachable!(),
    };
    totals.source_units += source_units;
    let mut seen = BTreeSet::new();
    let mut records = Vec::new();
    totals.timed("archive_setup", start);
    load_chunks(
        &mut la,
        &chunks,
        &initial,
        &mut seen,
        &mut records,
        totals,
        input,
        plan.terminal,
    )?;
    let start = input.profile.then(Instant::now);
    let mut ev = Evaluator {
        data: Data {
            records: &records,
            shapes: &shapes,
            patterns: &patterns,
            chroms: &chroms,
        },
        steps: totals.steps,
        max_steps: input.max_steps,
        unknown_omitted: 0,
    };
    let mut selected_records = BTreeSet::new();
    let direct_records = input.engine == Engine::Auto && plan.unit == Unit::Record;
    let mut selected_indices = Vec::new();
    let mut selected_classes = BTreeSet::new();
    let mut selected_cells = BTreeSet::new();
    let within = plan
        .within
        .as_ref()
        .map(|value| {
            Ok(Node {
                at: 0,
                ty: Type::Truth,
                op: match value {
                    Value::Region(r) => Op::Geometry {
                        kind: "overlaps".into(),
                        region: r.clone(),
                        minimum: 1,
                        total: false,
                    },
                    Value::Pattern(p) => Op::Match(p.clone()),
                    _ => bail!("invalid universe"),
                },
            })
        })
        .transpose()?;
    for (i, r) in records.iter().enumerate() {
        if witness(&mut ev, i, within.as_ref())? {
            if direct_records {
                selected_indices.push((r.ordinal, smallvec::smallvec![i]));
            } else {
                selected_records.insert(r.ordinal);
                selected_classes.insert(r.molecule.umi_class);
                selected_cells.insert(r.molecule.cell);
            }
        }
    }
    totals.steps = ev.steps;
    totals.timed("universe_witness", start);
    if plan.within.is_none() && plan.unit == Unit::Cell {
        selected_cells.extend(0..cells.len() as u32);
    }
    let mut closure = "none";
    if plan.unit != Unit::Record && !selected_cells.is_empty() {
        let closure_chunks = if plan.unit == Unit::Class {
            if let Some(index) = la.access_index()? {
                closure = "class index";
                index.class_chunks(selected_classes.iter().copied())?
            } else {
                if !input.allow_full_scan {
                    bail!("class closure requires full scan allowance");
                }
                closure = "explicit full scan";
                (0..chunks.len()).collect()
            }
        } else {
            if !input.allow_full_scan {
                bail!("cell closure requires full scan allowance");
            }
            closure = "explicit full scan";
            (0..chunks.len()).collect()
        };
        load_chunks(
            &mut la,
            &chunks,
            &closure_chunks,
            &mut seen,
            &mut records,
            totals,
            input,
            plan.terminal,
        )?;
    }
    if needs_denominator_scan(plan) {
        load_chunks(
            &mut la,
            &chunks,
            &(0..chunks.len()).collect::<Vec<_>>(),
            &mut seen,
            &mut records,
            totals,
            input,
            false,
        )?;
    }
    access.push(json!({"sample":member.sample,"initial_routing":routing,"initial_chunks":initial.len(),"closure":closure,"actual_chunks":seen.len(),"allow_full_scan":input.allow_full_scan}));
    let start = input.profile.then(Instant::now);
    let mut units: BTreeMap<u64, smallvec::SmallVec<[usize; 1]>> = BTreeMap::new();
    for (i, r) in records.iter().enumerate() {
        if direct_records {
            break;
        }
        let m = &r.molecule;
        let selected = match plan.unit {
            Unit::Record => selected_records.contains(&r.ordinal),
            Unit::Class => selected_classes.contains(&m.umi_class),
            Unit::Cell => selected_cells.contains(&m.cell),
            _ => false,
        };
        if selected {
            units
                .entry(match plan.unit {
                    Unit::Record => r.ordinal,
                    Unit::Class => u64::from(m.umi_class),
                    _ => u64::from(m.cell),
                })
                .or_default()
                .push(i);
        }
    }
    if plan.unit == Unit::Cell {
        for &cell in &selected_cells {
            units.entry(u64::from(cell)).or_default();
        }
    }
    let mut ev = Evaluator {
        data: Data {
            records: &records,
            shapes: &shapes,
            patterns: &patterns,
            chroms: &chroms,
        },
        steps: totals.steps,
        max_steps: input.max_steps,
        unknown_omitted: 0,
    };
    if let Some(by) = plan.stages.iter().find_map(|s| match s {
        Stage::Tally(_, by) | Stage::Summarize(_, by) => Some(by),
        _ => None,
    }) {
        let mut source_groups: BTreeMap<u32, u64> = BTreeMap::new();
        if plan.unit == Unit::Cell {
            for cell in 0..cells.len() as u32 {
                source_groups.insert(cell, 1);
            }
        } else if needs_denominator_scan(plan) {
            let mut classes = BTreeSet::new();
            for r in &records {
                if plan.unit == Unit::Record || classes.insert(r.molecule.umi_class) {
                    *source_groups.entry(r.molecule.cell).or_default() += 1;
                }
            }
        } else {
            source_groups.insert(u32::MAX, source_units);
        }
        for (cell, count) in source_groups {
            let row = if cell == u32::MAX {
                BTreeMap::from([("sample".into(), Value::String(member.sample.clone()))])
            } else {
                let barcode = String::from_utf8(
                    crate::querycmd::unpack_cell_bytes(cells[cell as usize]).to_vec(),
                )?;
                resources.row(member, &barcode)?
            };
            let scope = Scope {
                records: &[],
                record: None,
                alignment: None,
                signature: None,
                omitted: None,
            };
            let group = evaluate_fields(by, &mut ev, scope, &row)?;
            let k = key(&group)?;
            if !populations.contains_key(&k) && populations.len() >= input.max_rows {
                bail!("source-scope group budget exceeded");
            }
            populations.entry(k).or_insert((group, 0, 0)).2 += count;
        }
    }
    let mut cell_rows: BTreeMap<u32, Row> = BTreeMap::new();
    let mut cell_groups: BTreeMap<u32, (Row, String)> = BTreeMap::new();
    let by = plan.stages.iter().find_map(|s| match s {
        Stage::Tally(_, by) | Stage::Summarize(_, by) => Some(by),
        _ => None,
    });
    let aggregate_fields = plan.stages.iter().find_map(|s| match s {
        Stage::Tally(fields, _) | Stage::Summarize(fields, _) => Some(fields),
        _ => None,
    });
    // An empty or sample-only grouping has one key per member, not one key per
    // cell. Still charge each cell's original grouping evaluation exactly once.
    let member_group = input.engine == Engine::Auto
        && by.is_some_and(|fields| {
            fields
                .iter()
                .all(|(_, n)| matches!(&n.op, Op::Field(name) if name == "sample"))
        });
    let mut grouped_cells = if member_group {
        vec![false; cells.len()]
    } else {
        Vec::new()
    };
    let mut group_ids = BTreeMap::<String, usize>::new();
    let mut cell_group_ids = BTreeMap::<u32, usize>::new();
    let mut fast_keys = HashMap::<(usize, States), usize>::new();
    let mut fast_values: Vec<PartialAggregate> = Vec::new();
    let distinct_fields: Vec<_> = aggregate_fields
        .map(|fields| {
            fields
                .iter()
                .map(|(_, n)| matches!(&n.op, Op::Call(name, _) if name == "count_distinct"))
                .collect()
        })
        .unwrap_or_default();
    totals.timed("unit_and_denominator_setup", start);
    let start = input.profile.then(Instant::now);
    let needs_unit_id = input.engine == Engine::Reference || physical::uses_field(plan, "unit.id");
    // No metadata validation is skipped: this applies only when there are no declared
    // metadata columns and neither cell identity field is referenced anywhere.
    let sample_only = input.engine == Engine::Auto
        && resources.metadata.columns.is_empty()
        && !physical::uses_field(plan, "cell")
        && !physical::uses_field(plan, "cell.id");
    if sample_only {
        cell_rows.insert(
            0,
            BTreeMap::from([("sample".into(), Value::String(member.sample.clone()))]),
        );
    }
    let units: Box<dyn Iterator<Item = (u64, smallvec::SmallVec<[usize; 1]>)>> = if direct_records {
        Box::new(selected_indices.into_iter())
    } else {
        Box::new(units.into_iter())
    };
    for (id, indices) in units {
        totals.population += 1;
        let cell = if plan.unit == Unit::Cell {
            id as u32
        } else {
            records[indices[0]].molecule.cell
        };
        let packed = *cells.get(cell as usize).context("invalid cell id")?;
        let row_cell = if sample_only { 0 } else { cell };
        if let std::collections::btree_map::Entry::Vacant(entry) = cell_rows.entry(row_cell) {
            let barcode = String::from_utf8(crate::querycmd::unpack_cell_bytes(packed).to_vec())?;
            entry.insert(resources.row(member, &barcode)?);
        }
        let row = cell_rows.get_mut(&row_cell).expect("initialized cell row");
        if needs_unit_id {
            row.insert(
                "unit.id".into(),
                Value::String(format!("{}:{:?}:{id}", member.sample, plan.unit)),
            );
        }
        let scope = Scope {
            records: &indices,
            record: if plan.unit == Unit::Record {
                Some(indices[0])
            } else {
                None
            },
            alignment: None,
            signature: None,
            omitted: None,
        };
        // Denominators are computed before evidence where, using base metadata group keys.
        if let Some(by) = by {
            let group_cell = if member_group { 0 } else { cell };
            if member_group && !grouped_cells[cell as usize] {
                grouped_cells[cell as usize] = true;
                // The first key's normal evaluation below already charges its steps.
                if cell_groups.contains_key(&group_cell) {
                    for _ in by {
                        ev.steps += 1;
                        if ev.steps > ev.max_steps {
                            bail!("expression evaluation budget exceeded");
                        }
                    }
                }
            }
            if let std::collections::btree_map::Entry::Vacant(entry) = cell_groups.entry(group_cell)
            {
                let group = evaluate_fields(by, &mut ev, scope, row)?;
                let key = key(&group)?;
                let next = group_ids.len();
                let id = *group_ids.entry(key.clone()).or_insert(next);
                cell_group_ids.insert(group_cell, id);
                entry.insert((group, key));
            }
            let (group, key) = &cell_groups[&group_cell];
            if !populations.contains_key(key) && populations.len() >= input.max_rows {
                bail!("population group budget exceeded");
            }
            if let Some(population) = populations.get_mut(key) {
                population.1 += 1;
            } else {
                populations.insert(key.clone(), (group.clone(), 1, 0));
            }
        }
        if let Some(kernel) = &mut kernel {
            let outcome = kernel.run(&mut ev, scope, row)?;
            let (states, tally) = match outcome {
                Outcome::Dropped(state) => {
                    totals.where_unknown += u64::from(state == Truth::Unknown);
                    continue;
                }
                Outcome::Tally(states) => (states, true),
                Outcome::Summary => (States::new(), false),
            };
            let group_cell = if member_group { 0 } else { cell };
            let group_id = cell_group_ids[&group_cell];
            let lookup = (group_id, states);
            let index = if let Some(&i) = fast_keys.get(&lookup) {
                i
            } else {
                let (base_group, base_key) = &cell_groups[&group_cell];
                let mut k = base_key.clone();
                for state in &lookup.1 {
                    k.push(match state {
                        Truth::True => 'T',
                        Truth::False => 'F',
                        Truth::Unknown => 'U',
                    });
                }
                if !accumulators.contains_key(&k) && accumulators.len() >= input.max_rows {
                    bail!("aggregation group budget exceeded");
                }
                let fields = aggregate_fields.context("missing aggregation")?;
                let initial: Vec<_> = if tally {
                    vec![Value::Count(0)]
                } else {
                    fields
                        .iter()
                        .map(|(_, n)| {
                            if n.ty == Type::Number {
                                Value::Number(0.0)
                            } else {
                                Value::Count(0)
                            }
                        })
                        .collect()
                };
                let accumulator = accumulators
                    .entry(k.clone())
                    .or_insert_with(|| Accumulator {
                        key: base_group
                            .iter()
                            .map(|(k, v)| (k.clone(), v.clone()))
                            .chain(fields.iter().zip(&lookup.1).map(|((name, _), state)| {
                                (format!("{name}_state"), Value::Truth(*state))
                            }))
                            .collect(),
                        values: initial,
                        distinct: vec![HashSet::new(); fields.len()],
                    });
                let index = fast_values.len();
                // Continue the existing accumulator in source/unit order. In particular,
                // do not reassociate floating-point sums across federation members.
                fast_values.push(PartialAggregate {
                    key: k,
                    values: std::mem::take(&mut accumulator.values),
                    distinct: std::mem::take(&mut accumulator.distinct),
                });
                fast_keys.insert(lookup, index);
                index
            };
            let a = &mut fast_values[index];
            if tally {
                add(&mut a.values[0], Value::Count(1))?;
            }
            for (i, b) in kernel.values.drain(..).enumerate() {
                if b == Value::Null {
                    continue;
                }
                if distinct_fields[i] {
                    insert_distinct(&mut a.distinct[i], b, input.max_records)?;
                } else {
                    add(&mut a.values[i], b)?;
                }
            }
            continue;
        }
        let mut aggregated = false;
        for stage in &plan.stages {
            match stage {
                Stage::Where(n) => {
                    let state = ev.eval(n, scope, row)?.truth()?;
                    if state != Truth::True {
                        if state == Truth::Unknown {
                            totals.where_unknown += 1;
                        }
                        break;
                    }
                }
                Stage::Derive(fields) => {
                    for (name, n) in fields {
                        let value = ev.eval(n, scope, row)?;
                        row.insert(name.clone(), value);
                    }
                }
                Stage::Select(fields) => {
                    if aggregated {
                        break;
                    }
                    let selected = evaluate_fields(fields, &mut ev, scope, row)?;
                    if output.len() >= input.max_rows {
                        bail!("result row budget exceeded");
                    }
                    *columns = fields
                        .iter()
                        .map(|(k, n)| (k.clone(), n.ty.clone()))
                        .collect();
                    output.push(selected);
                }
                Stage::Tally(fields, _) | Stage::Summarize(fields, _) => {
                    aggregated = true;
                    let tally = matches!(stage, Stage::Tally(..));
                    let group_cell = if member_group { 0 } else { cell };
                    let (base_group, base_key) = cell_groups
                        .get(&group_cell)
                        .context("missing cached group")?;
                    let mut k = base_key.clone();
                    let mut states = smallvec::SmallVec::<[Value; 8]>::new();
                    if tally {
                        for (_, n) in fields {
                            let value = Value::Truth(ev.eval(n, scope, row)?.truth()?);
                            k.push(match value.truth()? {
                                Truth::True => 'T',
                                Truth::False => 'F',
                                Truth::Unknown => 'U',
                            });
                            states.push(value);
                        }
                    }
                    if !accumulators.contains_key(&k) && accumulators.len() >= input.max_rows {
                        bail!("aggregation group budget exceeded");
                    }
                    let a = accumulators.entry(k).or_insert_with(|| Accumulator {
                        key: base_group
                            .iter()
                            .map(|(k, v)| (k.clone(), v.clone()))
                            .chain(
                                fields
                                    .iter()
                                    .zip(states)
                                    .map(|((name, _), value)| (format!("{name}_state"), value)),
                            )
                            .collect(),
                        values: if tally {
                            vec![Value::Count(0)]
                        } else {
                            fields
                                .iter()
                                .map(|(_, n)| {
                                    if n.ty == Type::Number {
                                        Value::Number(0.0)
                                    } else {
                                        Value::Count(0)
                                    }
                                })
                                .collect()
                        },
                        distinct: vec![HashSet::new(); fields.len()],
                    });
                    if tally {
                        add(&mut a.values[0], Value::Count(1))?;
                    } else {
                        for (i, (_, n)) in fields.iter().enumerate() {
                            let Op::Call(name, args) = &n.op else {
                                bail!("invalid aggregate")
                            };
                            match name.as_str() {
                                "count" => add(&mut a.values[i], Value::Count(1))?,
                                "sum" => {
                                    let value = ev.eval(&args[0], scope, row)?;
                                    if value != Value::Null {
                                        add(&mut a.values[i], value)?;
                                    }
                                }
                                "count_distinct" => {
                                    let value = ev.eval(&args[0], scope, row)?;
                                    if value != Value::Null {
                                        insert_distinct(
                                            &mut a.distinct[i],
                                            value,
                                            input.max_records,
                                        )?;
                                    }
                                }
                                _ => {
                                    let v = ev.eval(&args[0], scope, row)?.truth()?;
                                    let want = match name.as_str() {
                                        "count_true" => Truth::True,
                                        "count_false" => Truth::False,
                                        _ => Truth::Unknown,
                                    };
                                    add(&mut a.values[i], Value::Count(u64::from(v == want)))?;
                                }
                            }
                        }
                    }
                }
                Stage::Sort(_) | Stage::Take(_) => {}
            }
        }
    }
    for PartialAggregate {
        key,
        values,
        distinct,
    } in fast_values
    {
        let a = accumulators
            .get_mut(&key)
            .context("missing fused accumulator")?;
        a.values = values;
        a.distinct = distinct;
    }
    totals.timed("evaluate_and_aggregate", start);
    totals.steps = ev.steps;
    totals.unknown += ev.unknown_omitted;
    Ok(())
}
fn add(a: &mut Value, b: Value) -> Result<()> {
    if let (Value::Count(a), Value::Count(b)) = (&mut *a, &b) {
        *a = a.checked_add(*b).context("Count overflow")?;
        return Ok(());
    }
    *a = scalar_binary("+", a.clone(), b)?;
    Ok(())
}
fn insert_distinct(set: &mut HashSet<String>, value: Value, limit: usize) -> Result<()> {
    let key = key_value(&value)?;
    if set.len() >= limit && !set.contains(&key) {
        bail!("distinct-value budget exceeded");
    }
    set.insert(key);
    Ok(())
}
fn verify_root(
    la: &mut LazyArchive,
    member: &Member,
    explanation: &serde_json::Value,
) -> Result<()> {
    let root = la
        .reader()
        .content_commitment()
        .context("archive lost authenticated root")?
        .to_hex();
    let expected = explanation["sources"]
        .as_array()
        .context("missing resolved sources")?
        .iter()
        .find(|s| s["sample"].as_str() == Some(&member.sample))
        .and_then(|s| s["content_root"].as_str())
        .context("missing expected archive root")?;
    if root != expected {
        bail!(
            "archive changed between resource validation and execution: {}",
            member.path.display()
        );
    }
    Ok(())
}
fn compare_value(a: Option<&Value>, b: Option<&Value>) -> std::cmp::Ordering {
    match (a, b) {
        (Some(Value::Null), Some(Value::Null)) => std::cmp::Ordering::Equal,
        (Some(Value::Null), _) => std::cmp::Ordering::Less,
        (_, Some(Value::Null)) => std::cmp::Ordering::Greater,
        (Some(Value::String(a)), Some(Value::String(b))) => a.cmp(b),
        (Some(Value::Number(a)), Some(Value::Number(b))) => a.total_cmp(b),
        (Some(Value::Count(a)), Some(Value::Count(b))) => a.cmp(b),
        (Some(a), Some(b)) => a.json().to_string().cmp(&b.json().to_string()),
        _ => a.is_some().cmp(&b.is_some()),
    }
}
