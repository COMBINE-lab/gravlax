//! Isolate the former mixed-strand placement lookup against direct block/strand lookup.
//! cargo run --release -p gravlax-anno --example genefull_index_benchmark -- annotation.aic
//! Deterministic gene-local queries; does not measure archive decoding or EM iterations.
use anno::assign::{GeneFullIndex, SoloStrand};
use evidence_io::{Block, Junction, Placement, Strand};
use std::{hint::black_box, time::Instant};

struct Span {
    start: u32,
    end: u32,
    max_end: u32,
    gene: u32,
    reverse: bool,
}

fn reference_index(annotation: &anno::Annotation) -> hashbrown::HashMap<u32, Vec<Span>> {
    let mut bounds = hashbrown::HashMap::<(u32, u32, bool), (u32, u32)>::new();
    for t in &annotation.transcripts {
        let (start, end) = t.span();
        if start >= end {
            continue;
        }
        let entry = bounds
            .entry((t.chrom, t.gene, t.strand_rev))
            .or_insert((start, end));
        entry.0 = entry.0.min(start);
        entry.1 = entry.1.max(end);
    }
    let mut chroms: hashbrown::HashMap<u32, Vec<Span>> = hashbrown::HashMap::new();
    for ((chrom, gene, reverse), (start, end)) in bounds {
        chroms.entry(chrom).or_default().push(Span {
            start,
            end,
            max_end: 0,
            gene,
            reverse,
        });
    }
    for spans in chroms.values_mut() {
        spans.sort_unstable_by_key(|s| (s.start, s.end, s.gene, s.reverse));
        let mut max_end = 0;
        for span in spans {
            max_end = max_end.max(span.end);
            span.max_end = max_end;
        }
    }
    chroms
}

fn reference_query(
    index: &hashbrown::HashMap<u32, Vec<Span>>,
    p: &Placement,
    strand: SoloStrand,
    out: &mut Vec<u32>,
) {
    out.clear();
    if let Some(spans) = index.get(&p.chrom) {
        for b in &p.blocks {
            if b.start >= b.end {
                continue;
            }
            let hi = spans.partition_point(|s| s.start < b.end);
            for span in spans[..hi].iter().rev() {
                if span.max_end <= b.start {
                    break;
                }
                if span.end > b.start
                    && strand.accepts(matches!(p.strand, Strand::Reverse), span.reverse)
                {
                    out.push(span.gene);
                }
            }
        }
    }
    out.sort_unstable();
    out.dedup();
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args_os().nth(1).expect("annotation path");
    let annotation = anno::Annotation::from_path(std::path::Path::new(&path))?;
    anyhow::ensure!(
        !annotation.transcripts.is_empty(),
        "annotation has no transcripts"
    );
    let old = reference_index(&annotation);
    let new = GeneFullIndex::new(&annotation);
    let mut seed = 2718u64;
    let mut random = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 32) as u32
    };
    let queries: Vec<_> = (0..200_000)
        .map(|_| {
            let t = &annotation.transcripts[random() as usize % annotation.transcripts.len()];
            let (start, end) = t.span();
            let pos = start + random() % (end - start).max(1);
            (t.chrom, random() % 2 == 0, pos, random() % 2 == 0)
        })
        .collect();
    let mut placement = Placement {
        chrom: 0,
        strand: Strand::Forward,
        blocks: Vec::new(),
        junctions: Vec::new(),
        nm: 0,
        score: 0,
        nh: 1,
        clip: (0, 0),
    };
    let mut observed = Vec::new();
    let mut expected = Vec::new();
    println!("strand\trepetition\timplementation\tqueries\tseconds\tgene_hits");
    for strand in [
        SoloStrand::Forward,
        SoloStrand::Reverse,
        SoloStrand::Unstranded,
    ] {
        // Check every query before timing; use direct reconstruction as the reference.
        for &(chrom, reverse, pos, spliced) in &queries {
            let blocks = [(pos, pos + 50), (pos + 150, pos + 200)];
            let blocks = &blocks[..if spliced { 2 } else { 1 }];
            placement.chrom = chrom;
            placement.strand = if reverse {
                Strand::Reverse
            } else {
                Strand::Forward
            };
            placement.blocks.clear();
            placement
                .blocks
                .extend(blocks.iter().map(|&(start, end)| Block { start, end }));
            reference_query(&old, &placement, strand, &mut expected);
            new.genes_for_blocks_into(
                chrom,
                reverse,
                strand,
                blocks.iter().copied(),
                &mut observed,
            );
            assert_eq!(observed, expected);
        }
        for rep in 0..3 {
            for optimized in if rep % 2 == 0 {
                [false, true]
            } else {
                [true, false]
            } {
                let start = Instant::now();
                let mut hits = 0usize;
                for _ in 0..5 {
                    for &(chrom, reverse, pos, spliced) in &queries {
                        let blocks = [(pos, pos + 50), (pos + 150, pos + 200)];
                        let blocks = &blocks[..if spliced { 2 } else { 1 }];
                        if optimized {
                            new.genes_for_blocks_into(
                                chrom,
                                reverse,
                                strand,
                                blocks.iter().copied(),
                                &mut observed,
                            );
                        } else {
                            placement.chrom = chrom;
                            placement.strand = if reverse {
                                Strand::Reverse
                            } else {
                                Strand::Forward
                            };
                            placement.blocks.clear();
                            placement
                                .blocks
                                .extend(blocks.iter().map(|&(start, end)| Block { start, end }));
                            placement.junctions.clear();
                            placement
                                .junctions
                                .extend(blocks.windows(2).map(|w| Junction {
                                    donor: w[0].1,
                                    acceptor: w[1].0,
                                }));
                            reference_query(&old, &placement, strand, &mut observed);
                        }
                        hits += black_box(&observed).len();
                    }
                }
                println!(
                    "{strand:?}\t{rep}\t{}\t{}\t{:.6}\t{hits}",
                    if optimized {
                        "direct-strand-index"
                    } else {
                        "placement-mixed-index"
                    },
                    queries.len() * 5,
                    start.elapsed().as_secs_f64()
                );
            }
        }
    }
    Ok(())
}
