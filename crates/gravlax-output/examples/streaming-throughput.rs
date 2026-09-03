//! Small dependency-free comparison of the legacy allocating row adapter and the borrowed writer.
//! Run with, for example:
//! `cargo run --release -p gravlax-output --example streaming-throughput -- 500000`.

use gravlax_output::{
    write_streaming_table, write_table, DataType, Field, OutputFormat, ResultContext, ScalarValue,
    TableSchema,
};
use std::hint::black_box;
use std::io::{sink, BufWriter, Write};
use std::time::{Duration, Instant};

struct SourceRow {
    gene: String,
    class: String,
    count: u64,
    delta: i64,
    score: f64,
    retained: bool,
}

fn schema() -> TableSchema {
    TableSchema::new(
        "gravlax.benchmark.rows.v1",
        vec![
            Field::new("gene", DataType::String),
            Field::new("class", DataType::String),
            Field::new("count", DataType::UInt64),
            Field::new("delta", DataType::Int64),
            Field::new("score", DataType::Float64),
            Field::new("retained", DataType::Boolean),
        ],
    )
    .unwrap()
}

fn legacy(rows: &[SourceRow], format: OutputFormat) -> Duration {
    let started = Instant::now();
    let mut output = BufWriter::with_capacity(64 * 1024, sink());
    write_table(
        &mut output,
        &schema(),
        rows.iter().map(|row| {
            vec![
                ScalarValue::String(row.gene.clone()),
                ScalarValue::String(row.class.clone()),
                ScalarValue::UInt64(row.count),
                ScalarValue::Int64(row.delta),
                ScalarValue::Float64(row.score),
                ScalarValue::Boolean(row.retained),
            ]
        }),
        format,
        &ResultContext::default(),
    )
    .unwrap();
    output.flush().unwrap();
    started.elapsed()
}

fn streaming(rows: &[SourceRow], format: OutputFormat) -> Duration {
    let schema = schema();
    let started = Instant::now();
    let mut output = BufWriter::with_capacity(64 * 1024, sink());
    write_streaming_table(
        &mut output,
        &schema,
        format,
        &ResultContext::default(),
        None,
        |output| {
            for source in rows {
                output.write_row_with(|row| {
                    row.string(&source.gene)?;
                    row.string(&source.class)?;
                    row.uint64(source.count)?;
                    row.int64(source.delta)?;
                    row.float64(source.score)?;
                    row.boolean(source.retained)
                })?;
            }
            Ok(())
        },
    )
    .unwrap();
    started.elapsed()
}

fn best_of_three(mut run: impl FnMut() -> Duration) -> Duration {
    (0..3).map(|_| run()).min().unwrap()
}

fn main() {
    let count = std::env::args()
        .nth(1)
        .map(|value| {
            value
                .parse::<usize>()
                .expect("row count must be an integer")
        })
        .unwrap_or(250_000);
    let rows: Vec<_> = (0..count)
        .map(|index| SourceRow {
            gene: format!("ENSG{:011}", index % 10_000),
            class: format!("evidence-class-{}", index % 7),
            count: index as u64 * 3,
            delta: index as i64 % 19 - 9,
            score: (index % 1_000) as f64 / 37.0,
            retained: index % 3 != 0,
        })
        .collect();
    let rows = black_box(rows);
    println!("rows\tformat\tlegacy_ms\tstreaming_ms\tstreaming_over_legacy");
    for format in [OutputFormat::Tsv, OutputFormat::Json] {
        // Warm both paths before recording the minimum of three runs.
        black_box(legacy(&rows[..rows.len().min(1_000)], format));
        black_box(streaming(&rows[..rows.len().min(1_000)], format));
        let legacy = best_of_three(|| legacy(&rows, format));
        let streaming = best_of_three(|| streaming(&rows, format));
        println!(
            "{}\t{:?}\t{:.3}\t{:.3}\t{:.4}",
            rows.len(),
            format,
            legacy.as_secs_f64() * 1_000.0,
            streaming.as_secs_f64() * 1_000.0,
            streaming.as_secs_f64() / legacy.as_secs_f64()
        );
    }
}
