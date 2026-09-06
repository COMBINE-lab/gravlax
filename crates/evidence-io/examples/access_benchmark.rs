//! Deterministic warm-cache microbenchmark; run release builds before and after changes.
use anyhow::Result;
use evidence_io::{
    format::{SectionReader, SectionWriter},
    rans,
};
use std::{hint::black_box, time::Instant};

fn main() -> Result<()> {
    let dir = std::env::temp_dir().join(format!("gravlax-access-bench-{}", std::process::id()));
    std::fs::create_dir(&dir)?;
    let path = dir.join("sections.aie");
    let mut writer = SectionWriter::create_new(&path, 1)?;
    for i in 0..4096 {
        writer.section(&format!("c{i}"), &(i as u32).to_le_bytes())?;
    }
    writer.finish()?;
    let mut times = Vec::new();
    for _ in 0..7 {
        let start = Instant::now();
        for _ in 0..5 {
            black_box(SectionReader::open(&path)?);
        }
        let open = start.elapsed().as_secs_f64() / 5.0;
        let reader = SectionReader::open(&path)?;
        let start = Instant::now();
        for i in 0..20_000 {
            black_box(reader.read_compressed_at(&format!("c{}", i % 4096))?);
        }
        times.push(
            serde_json::json!({"open_seconds":open,"reads_seconds":start.elapsed().as_secs_f64()}),
        );
    }
    let values: Vec<u64> = (0..1_000_000)
        .map(|i| if i % 3 == 0 { i * 193 } else { i % 11 })
        .collect();
    let mut counts = [0; rans::NSYM];
    rans::count(&values, &mut counts);
    let table = rans::Table::from_counts(&counts)?;
    let mut encoded = Vec::new();
    rans::encode(&values, &table, &mut encoded);
    let mut decode = Vec::new();
    for _ in 0..7 {
        let start = Instant::now();
        let result = rans::decode(&encoded, &table)?;
        decode.push(start.elapsed().as_secs_f64());
        assert_eq!(result, values);
    }
    println!(
        "{}",
        serde_json::json!({"sections":4096,"reads":20000,"samples":times,"rans_seconds":decode,"rans_bytes":encoded.len()})
    );
    std::fs::remove_file(path)?;
    std::fs::remove_dir(dir)?;
    Ok(())
}
