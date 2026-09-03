use noodles_bam as bam;
fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).unwrap();
    let mut r = bam::io::reader::Builder::default().build_from_path(&path)?;
    r.read_header()?;
    let mut rec = bam::Record::default();
    let mut n = 0;
    while r.read_record(&mut rec)? != 0 && n < 400000 {
        n += 1;
        let mut keys = Vec::new();
        for f in rec.data().iter() {
            let (t, v) = f?;
            keys.push(format!("{}={:?}", String::from_utf8_lossy(&<[u8;2]>::from(t)), v));
        }
        if keys.iter().any(|k| k.starts_with("GX")) {
            println!("{}", keys.join("  "));
            if n > 0 { break; }
        }
    }
    println!("scanned {n}");
    Ok(())
}
