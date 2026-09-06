use super::*;
use super::{model::*, resources::Resources};
use crate::rows::{Extracted, MolChain, MolRec, PatAlt, SAME_SHAPE};
use evidence_io::archive::Shape;
use smallvec::smallvec;

struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "gravlax-gq-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}
fn fixture(path: &std::path::Path, chunk_bp: u32) {
    let record = |cell, class, pos, shape| MolRec {
        cell,
        umi_class: class,
        chrom: 0,
        strand_rev: false,
        chains: smallvec![MolChain {
            weight: 1,
            reps: smallvec![(pos, shape)]
        }],
        mms: smallvec![],
    };
    let mut compact = record(0, 0, 100, 0);
    compact.chains[0] = MolChain {
        weight: 10,
        reps: smallvec![(100, 0), (110, 0)],
    };
    let mm = MolRec {
        cell: 0,
        umi_class: 2,
        chrom: 0,
        strand_rev: false,
        chains: smallvec![],
        mms: smallvec![(700, 1, 0, 5)],
    };
    let x = Extracted {
        mols: vec![compact, record(1, 1, 150, 1), record(0, 0, 500, 1), mm],
        edges: vec![],
        cells: vec![0, 1, 2],
        shapes: vec![
            Shape {
                blocks: vec![(0, 20)],
            },
            Shape {
                blocks: vec![(0, 10), (30, 10)],
            },
        ],
        patterns: vec![vec![
            PatAlt {
                chrom: 0,
                offset: 0,
                strand_flip: false,
                shape: SAME_SHAPE,
            },
            PatAlt {
                chrom: 0,
                offset: 100,
                strand_flip: false,
                shape: 0,
            },
        ]],
        n_classes: 3,
        chrom_names: vec!["1".into()],
    };
    crate::archivecmd::write_archive(&x, path, 3, chunk_bp, None).unwrap();
}

#[test]
#[ignore = "export a deterministic release smoke fixture to a new caller-supplied file"]
fn export_release_smoke_fixture() {
    use std::io::Write;
    let scratch = Scratch::new();
    let archive = scratch.0.join("smoke.aie");
    let chunk_bp = std::env::var("GRAVLAX_GQ_SMOKE_CHUNK_BP")
        .map(|s| s.parse().unwrap())
        .unwrap_or(100);
    fixture(&archive, chunk_bp);
    let path = std::env::var_os("GRAVLAX_GQ_SMOKE_OUTPUT").expect("set GRAVLAX_GQ_SMOKE_OUTPUT");
    let mut output = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .unwrap();
    output.write_all(&std::fs::read(archive).unwrap()).unwrap();
}

#[test]
#[ignore = "export a clearly synthetic 600k-record three-chromosome stress archive"]
fn export_release_scale_fixture() {
    use std::io::Write;
    let path = std::env::var_os("GRAVLAX_GQ_SCALE_OUTPUT").expect("set GRAVLAX_GQ_SCALE_OUTPUT");
    let mut output = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .unwrap();
    let scratch = Scratch::new();
    let archive = scratch.0.join("scale.aie");
    let mols = (0..600_000_u32)
        .map(|i| {
            let chrom = i / 200_000;
            let position = (i % 200_000) * 10;
            MolRec {
                cell: i % 60_000,
                umi_class: i,
                chrom,
                strand_rev: i % 2 == 0,
                chains: smallvec![MolChain {
                    weight: 3,
                    reps: smallvec![(position, 1), (position + 2, 2)]
                }],
                mms: if i % 5 == 0 {
                    smallvec![(position, 1, chrom, 2)]
                } else {
                    smallvec![]
                },
            }
        })
        .collect();
    let data = Extracted {
        mols,
        edges: vec![],
        cells: (0..60_000).collect(),
        shapes: vec![
            Shape {
                blocks: vec![(0, 20)],
            },
            Shape {
                blocks: vec![(0, 10), (30, 10)],
            },
            Shape {
                blocks: vec![(0, 8), (28, 10)],
            },
        ],
        patterns: (0..3)
            .map(|chrom| {
                vec![
                    PatAlt {
                        chrom,
                        offset: 0,
                        strand_flip: false,
                        shape: SAME_SHAPE,
                    },
                    PatAlt {
                        chrom: (chrom + 1) % 3,
                        offset: 100,
                        strand_flip: true,
                        shape: SAME_SHAPE,
                    },
                ]
            })
            .collect(),
        n_classes: 600_000,
        chrom_names: vec!["1".into(), "2".into(), "3".into()],
    };
    crate::archivecmd::write_archive(&data, &archive, 3, 100_000, None).unwrap();
    output.write_all(&std::fs::read(archive).unwrap()).unwrap();
}
fn query(text: &str, archive: &std::path::Path) -> execute::Table {
    let source = format!("header {{gq=1,assembly=\"test\"}} {text}");
    let doc = parser::parse(&source).unwrap();
    let mut input = Input {
        engine: Engine::Auto,
        parallel_decode: false,
        profile: false,
        query: PathBuf::new(),
        bind: vec![format!("x={}", archive.display())],
        project: None,
        metadata: None,
        allow_full_scan: true,
        max_chunks: 100,
        max_records: 100,
        max_steps: 1_000_000,
        max_rows: 100,
        max_terminal_events: 100,
        output: None,
    };
    let resources = Resources::open(&input.bind, None, None, "test").unwrap();
    let plan = check::Compiler::new(&doc, resources.fields().unwrap())
        .unwrap()
        .with_resources(&resources)
        .compile()
        .unwrap();
    let optimized = execute::execute(&plan, &resources, &input).unwrap();
    input.engine = Engine::Reference;
    let reference = execute::execute(&plan, &resources, &input).unwrap();
    assert_equivalent(&optimized, &reference);
    input.engine = Engine::Auto;
    input.parallel_decode = true;
    assert_equivalent(
        &optimized,
        &execute::execute(&plan, &resources, &input).unwrap(),
    );
    input.engine = Engine::Reference;
    if let Some(steps) = reference.summary["expression_steps"]
        .as_u64()
        .filter(|n| *n > 1)
    {
        input.max_steps = steps - 1;
        assert!(execute::execute(&plan, &resources, &input).is_err());
        input.engine = Engine::Auto;
        assert!(execute::execute(&plan, &resources, &input).is_err());
    }
    optimized
}
fn assert_equivalent(a: &execute::Table, b: &execute::Table) {
    assert_eq!(a.columns, b.columns);
    assert_eq!(a.rows, b.rows);
    let clean = |t: &execute::Table| {
        let mut s = t.summary.clone();
        s.as_object_mut().unwrap().remove("timings_ms");
        s["explain"]
            .as_object_mut()
            .unwrap()
            .remove("physical_strategy");
        s
    };
    assert_eq!(clean(a), clean(b));
}

#[test]
fn fused_summaries_filters_scopes_and_fallbacks() {
    let scratch = Scratch::new();
    let archive = scratch.0.join("fixture.aie");
    fixture(&archive, 100);
    for source in [
        "from @x.records |> within all |> derive {a=any unique{start_in(g[1:105..115])}, b=!a} |> summarize {n=count(),t=count_true(a),f=count_false(a),u=count_unknown(a),reads=sum(reads.total)}",
        "from @x.records |> within all |> derive {a=any unique{overlaps(g[1:100..120])}} |> where a |> derive {b=all unique{overlaps(g[1:100..120])}} |> tally {a,b}",
        "from @x.records |> within all |> derive {a=any multimap{all alternative{overlaps(g[1:700..820])}}} |> tally {a}",
        "from @x.classes |> within all |> derive {a=any record{any unique{overlaps(g[1:100..120])}},b=any record{any unique{j[1:510..530]}},joint=any record{any unique{overlaps(g[1:100..120]) & j[1:510..530]}}} |> tally {a,b,joint}",
        "from @x.cells |> within all |> derive {a=all record{nonempty(unique)}} |> summarize {t=count_true(a),f=count_false(a),reads=sum(reads.total)}",
        "from @x.records |> within all |> summarize {n=count_distinct(cell.id)}",
        "from @x.records |> within all |> derive {a=any unique{start_in(g[1:105..115])},b=is_unknown(a)} |> tally {a,b}",
        "from @x.records |> within all |> where any unique{overlaps(g[1:900..1000])} |> summarize {n=count(),r=sum(reads.total)}",
    ] { query(source, &archive); }
}

#[test]
fn parallel_decode_preserves_order_and_aggregate_semantics() {
    let scratch = Scratch::new();
    let archive = scratch.0.join("fixture.aie");
    fixture(&archive, 1);
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let four = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    for source in [
        "from @x.records |> within all |> derive {a=any unique{start_in(g[1:105..115])},b=all unique{end_in(g[1:120..131])}} |> tally {a,b} by {cell.id}",
        "from @x.classes |> within all |> derive {a=any record{any unique{j[1:160..180]}}} |> summarize {n=count(),t=count_true(a),r=sum(reads.total)} by {sample}",
        "from @x.cells |> within all |> summarize {n=count(),r=sum(reads.total)}",
        "from @x.records |> within all |> select {unit.id} |> sort {unit.id} |> take 3",
    ] {
        let a = one.install(|| query(source, &archive));
        let b = four.install(|| query(source, &archive));
        assert_equivalent(&a, &b);
    }
}

#[test]
fn physical_top_k_and_wide_aggregates_match_reference() {
    let scratch = Scratch::new();
    let archive = scratch.0.join("fixture.aie");
    fixture(&archive, 100);
    for n in [0, 1, 2, 4, 10] {
        query(&format!("from @x.records |> within all |> select {{cell.id,unit.id,reads.total}} |> sort {{cell.id}} |> take {n}"), &archive);
        query(&format!("from @x.records |> within all |> tally {{a=any unique{{start_in(g[1:105..115])}}}} |> sort {{count}} |> take {n}"), &archive);
    }
    let fields = (0..12)
        .map(|i| format!("q{i}=any unique{{overlaps(g[1:100..120])}}"))
        .collect::<Vec<_>>()
        .join(",");
    query(
        &format!("from @x.records |> within all |> tally {{{fields}}}"),
        &archive,
    );
    let fields = (0..12)
        .map(|i| format!("n{i}=count_distinct(cell.id)"))
        .collect::<Vec<_>>()
        .join(",");
    query(
        &format!("from @x.records |> within all |> summarize {{{fields}}}"),
        &archive,
    );
}
#[test]
fn complete_class_children_and_chunk_invariance() {
    let scratch = Scratch::new();
    let fine = scratch.0.join("fine.aie");
    let coarse = scratch.0.join("coarse.aie");
    fixture(&fine, 100);
    fixture(&coarse, 1000);
    let text="from @x.classes |> within g[1:100..130] |> derive {a=any record{any unique{j[1:510..530]}}} |> tally {a}";
    let a = query(text, &fine);
    let b = query(text, &coarse);
    assert_eq!(a.rows, b.rows);
    assert_eq!(a.rows.len(), 1);
    assert_eq!(a.rows[0]["a_state"], Value::Truth(Truth::True));
    let records=query("from @x.records |> within g[1:100..130] |> derive {a=any unique{j[1:510..530]}} |> tally {a}",&fine);
    assert_eq!(records.rows[0]["a_state"], Value::Truth(Truth::False));
}
#[test]
fn empty_denominators_multiplicity_and_grouped_alternatives() {
    let scratch = Scratch::new();
    let path = scratch.0.join("archive.aie");
    fixture(&path, 100);
    let empty = query(
        "from @x.records |> within g[1:900..950] |> summarize {n=count()}",
        &path,
    );
    assert_eq!(empty.rows[0]["n"], Value::Count(0));
    let tally = query(
        "from @x.records |> within g[1:900..950] |> derive {a=true} |> tally {a}",
        &path,
    );
    assert!(tally.rows.is_empty());
    assert_eq!(tally.columns.len(), 2);
    let sum = query(
        "from @x.records |> within all |> summarize {n=sum(reads.total)}",
        &path,
    );
    assert_eq!(sum.rows[0]["n"], Value::Count(17));
    let mm=query("from @x.records |> within g[1:700..750] |> derive {possible=any multimap{any alternative{j[1:710..730]}}, certain=any multimap{all alternative{j[1:710..730]}}} |> tally {possible,certain}",&path);
    assert_eq!(mm.rows[0]["possible_state"], Value::Truth(Truth::True));
    assert_eq!(mm.rows[0]["certain_state"], Value::Truth(Truth::False));
    let cells = query(
        "from @x.cells |> within all |> derive {a=all record{nonempty(unique)}} |> tally {a}",
        &path,
    );
    assert_eq!(cells.summary["population_units"], 3);
}
#[test]
fn support_enumeration_and_projection() {
    let scratch = Scratch::new();
    let path = scratch.0.join("archive.aie");
    fixture(&path, 100);
    let table = query(
        "from junctions(@x,within:g[1:1..900]) |> support(unit:class,by:{sample})",
        &path,
    );
    assert_eq!(table.rows.len(), 2);
    assert!(table.rows.iter().all(|r| r["support"] == Value::Count(1)));
    let projected=query("from @x.records |> within all |> summarize {n=count()} |> select {n} |> sort {n} |> take 1",&path);
    assert_eq!(projected.rows[0]["n"], Value::Count(4));
}

#[test]
fn functions_require_scientific_types_and_preserve_option_parameters() {
    let source="header {gq=1,assembly=\"test\",library=\"opposite\"} export fn supported(exon:Region<test,alignment>,minimum:Bases) -> Predicate<Alignment,test,alignment> = overlaps(exon,min:minimum) from @x.records |> within all |> derive {a=any unique{supported(g[1:10..30],12bp)&tx_strand(-)}} |> tally {a}";
    let doc = parser::parse(source).unwrap();
    assert!(check::Compiler::new(&doc, BTreeMap::new())
        .unwrap()
        .compile()
        .is_ok());
    let wrong = source.replace(
        "Predicate<Alignment,test,alignment>",
        "Predicate<Record,test,alignment>",
    );
    let doc = parser::parse(&wrong).unwrap();
    assert!(check::Compiler::new(&doc, BTreeMap::new())
        .unwrap()
        .compile()
        .is_err());
    let recursion="header {gq=1,assembly=\"test\"} fn unused() = unused() from @x.records |> within all |> summarize {n=count()}";
    let doc = parser::parse(recursion).unwrap();
    assert!(check::Compiler::new(&doc, BTreeMap::new()).is_err());
}

#[test]
fn pinned_features_and_union_overlap() {
    let scratch = Scratch::new();
    let archive = scratch.0.join("archive.aie");
    fixture(&archive, 100);
    let annotation = scratch.0.join("anno.gtf");
    std::fs::write(&annotation,"1\tx\texon\t101\t110\t.\t+\t.\tgene_id \"g\"; gene_name \"G\"; transcript_id \"t1\";\n1\tx\texon\t131\t140\t.\t+\t.\tgene_id \"g\"; gene_name \"G\"; transcript_id \"t1\";\n1\tx\texon\t101\t110\t.\t+\t.\tgene_id \"g\"; gene_name \"G\"; transcript_id \"t2\";\n").unwrap();
    let project = scratch.0.join("aie-project.yaml");
    std::fs::write(&project,"schema_version: 1\nname: test\nresources:\n  anno:\n    kind: annotation\n    path: anno.gtf\n    annotation_identity:\n      assembly: test\n      annotation: release1\n").unwrap();
    let resources = Resources::open(&[], None, Some(&project), "test").unwrap();
    let feature = resources.feature("anno", "G", false, Some("same")).unwrap();
    let Value::Feature(feature) = feature else {
        panic!()
    };
    assert_eq!(feature.exons.intervals, vec![(100, 110), (130, 140)]);
    assert_eq!(feature.junctions.len(), 1);
    let source="header {gq=1,assembly=\"test\",library=\"same\"} let G=gene(@anno,\"G\") from @x.records |> within G.span |> derive {a=any unique{overlaps(G.exons,blocks:total,min:12bp)}} |> tally {a}";
    let doc = parser::parse(source).unwrap();
    assert!(check::Compiler::new(&doc, resources.fields().unwrap())
        .unwrap()
        .with_resources(&resources)
        .compile()
        .is_ok());
    assert!(resources
        .feature("anno", "missing", false, Some("same"))
        .is_err());
    assert!(resources.feature("anno", "G", false, None).is_err());
}

#[test]
fn federation_partitioning_metadata_nulls_and_budgets() {
    let scratch = Scratch::new();
    let a = scratch.0.join("a.aie");
    let b = scratch.0.join("b.aie");
    fixture(&a, 100);
    fixture(&b, 1000);
    let federation = scratch.0.join("cohort.json");
    std::fs::write(&federation,serde_json::to_vec(&serde_json::json!({"schema_version":1,"assembly":"test","archives":[{"sample":"A","path":"a.aie","metadata":{"score":1e16}},{"sample":"B","path":"b.aie","metadata":{"score":2.0}}]})).unwrap()).unwrap();
    let metadata = scratch.0.join("metadata.json");
    std::fs::write(&metadata,r#"{"columns":{"score":{"type":"number"},"cell.type":{"type":"string","optional":true},"cell.good":{"type":"truth","optional":true}},"cells":{"A":{"AAAAAAAAAAAAAAAA":{"cell.type":"Astro","cell.good":true}}}}"#).unwrap();
    let input = Input {
        engine: Engine::Auto,
        parallel_decode: false,
        profile: false,
        query: PathBuf::new(),
        bind: vec![format!("x={}", federation.display())],
        project: None,
        metadata: Some(metadata.clone()),
        allow_full_scan: true,
        max_chunks: 100,
        max_records: 100,
        max_steps: 1_000_000,
        max_rows: 100,
        max_terminal_events: 100,
        output: None,
    };
    let resources = Resources::open(&input.bind, Some(&metadata), None, "test").unwrap();
    let doc=parser::parse("header {gq=1,assembly=\"test\"} from @x.records |> within all |> tally {good=cell.good} by {sample,cell.type}").unwrap();
    let plan = check::Compiler::new(&doc, resources.fields().unwrap())
        .unwrap()
        .with_resources(&resources)
        .compile()
        .unwrap();
    let table = execute::execute(&plan, &resources, &input).unwrap();
    let mut reference_input = input.clone();
    reference_input.engine = Engine::Reference;
    assert_equivalent(
        &table,
        &execute::execute(&plan, &resources, &reference_input).unwrap(),
    );
    // A partial-sum-per-member implementation would add 8.0 at the end and change
    // this result. Each 2.0 must instead be added in the original source/unit order.
    let summary_doc = parser::parse("header {gq=1,assembly=\"test\"} from @x.records |> within all |> summarize {n=count(),s=sum(score),d=count_distinct(cell.id)}").unwrap();
    let summary_plan = check::Compiler::new(&summary_doc, resources.fields().unwrap())
        .unwrap()
        .with_resources(&resources)
        .compile()
        .unwrap();
    let optimized = execute::execute(&summary_plan, &resources, &input).unwrap();
    assert_equivalent(
        &optimized,
        &execute::execute(&summary_plan, &resources, &reference_input).unwrap(),
    );
    assert_eq!(optimized.rows[0]["s"], Value::Number(4e16));
    assert_eq!(optimized.rows[0]["n"], Value::Count(8));
    assert_eq!(optimized.rows[0]["d"], Value::Count(2));
    assert_eq!(
        table
            .rows
            .iter()
            .map(|r| match r["count"] {
                Value::Count(n) => n,
                _ => panic!(),
            })
            .sum::<u64>(),
        8
    );
    assert!(table.rows.iter().all(|r| r["good_state"] != Value::Null));
    assert_eq!(table.summary["source_scope_units"], 8);
    let mut low = input;
    low.max_records = 1;
    assert!(execute::execute(&plan, &resources, &low).is_err());
    let duplicate = Resources::open(
        &[format!("x={}", federation.display())],
        None,
        None,
        "wrong",
    );
    assert!(duplicate.is_err());
}

#[test]
fn logical_digest_ignores_source_locations() {
    let a="header {gq=1,assembly=\"test\"} from @x.records |> within all |> derive {a=any unique{j[1:10..20]}} |> tally {a}";
    let b = a.replace("from", "\n\nfrom");
    let compile = |source: &str| {
        let doc = parser::parse(source).unwrap();
        check::Compiler::new(&doc, BTreeMap::new())
            .unwrap()
            .compile()
            .unwrap()
    };
    assert_eq!(
        logical_digest(&compile(a)).unwrap(),
        logical_digest(&compile(&b)).unwrap()
    );
}

#[test]
fn count_literals_ratios_and_null_sorting() {
    let scratch = Scratch::new();
    let archive = scratch.0.join("archive.aie");
    fixture(&archive, 100);
    let table=query("from @x.records |> within all |> where reads.total >= 2u64 |> summarize {n=count(),reads=sum(reads.total)} |> select {n,ratio=fraction(n,denominator:reads)}",&archive);
    assert_eq!(table.rows[0]["n"], Value::Count(2));
    assert_eq!(table.rows[0]["ratio"], Value::Number(2.0 / 15.0));
    let table=query("from @x.records |> within all |> summarize {n=count()} |> select {ratio=fraction(n,denominator:0u64)}",&archive);
    assert_eq!(table.rows[0]["ratio"], Value::Null);
}
