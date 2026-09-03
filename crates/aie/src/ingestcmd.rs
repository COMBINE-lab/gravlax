//! Preflight and recipe UX for annotation-independent archive ingestion.

use anyhow::{bail, Context, Result};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::record::data::field::Value;
use sam::alignment::record::Flags;
use serde::Serialize;
use std::collections::{BTreeMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Scan a BAM and whitelist before spending time building an archive.
    Check(CheckArgs),
    /// Print an explicit annotation-free STAR recipe for a supported chemistry.
    Recipe(RecipeArgs),
}

#[derive(Clone, Copy, Debug, ValueEnum, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Chemistry {
    /// 10x Chromium Single Cell 3' v2: 16 bp cell barcode and 10 bp UMI.
    #[value(name = "10x-3p-v2", alias = "chromium-3p-v2")]
    #[serde(rename = "10x-3p-v2")]
    TenX3pV2,
    /// 10x Chromium Single Cell 3' v3/v3.1: 16 bp cell barcode and 12 bp UMI.
    #[value(name = "10x-3p-v3", alias = "10x-3p-v3.1", alias = "chromium-3p-v3")]
    #[serde(rename = "10x-3p-v3")]
    TenX3pV3,
}

impl Chemistry {
    fn label(self) -> &'static str {
        match self {
            Self::TenX3pV2 => "10x-3p-v2",
            Self::TenX3pV3 => "10x-3p-v3",
        }
    }

    fn umi_len(self) -> usize {
        match self {
            Self::TenX3pV2 => 10,
            Self::TenX3pV3 => 12,
        }
    }

    fn default_whitelist(self) -> &'static str {
        match self {
            Self::TenX3pV2 => "737K-august-2016.txt",
            Self::TenX3pV3 => "3M-february-2018.txt",
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum ReportFormat {
    Text,
    Json,
}

#[derive(ClapArgs)]
struct CheckArgs {
    /// Annotation-free, coordinate-sorted STAR BAM intended for `aie ingest-archive`.
    #[arg(value_name = "BAM")]
    bam: PathBuf,

    /// Barcode whitelist intended for ingest (one 16 bp A/C/G/T barcode per line).
    #[arg(long, value_name = "FILE")]
    whitelist: PathBuf,

    /// Validate the observed UMI length against a supported library chemistry.
    #[arg(long, value_enum)]
    chemistry: Option<Chemistry>,

    /// Treat warnings, including an unexercised secondary-alignment check, as failure.
    #[arg(long)]
    strict: bool,

    /// Report representation.
    #[arg(long, value_enum, default_value_t = ReportFormat::Text)]
    format: ReportFormat,
}

#[derive(ClapArgs)]
struct RecipeArgs {
    /// Supported library chemistry and read structure.
    #[arg(long, value_enum)]
    chemistry: Chemistry,

    /// STAR genome directory built without a GTF or annotation-derived splice junctions.
    #[arg(long, default_value = "star-index-nogtf", value_name = "DIR")]
    genome_dir: PathBuf,

    /// Barcode/UMI read (R1 for supported 10x 3' chemistries).
    #[arg(long, default_value = "sample_R1.fastq.gz", value_name = "FASTQ")]
    read1: PathBuf,

    /// cDNA read (R2 for supported 10x 3' chemistries).
    #[arg(long, default_value = "sample_R2.fastq.gz", value_name = "FASTQ")]
    read2: PathBuf,

    /// Chemistry-matched barcode whitelist; defaults to the conventional 10x filename.
    #[arg(long, value_name = "FILE")]
    whitelist: Option<PathBuf>,

    /// STAR output prefix.
    #[arg(long, default_value = "align/", value_name = "PREFIX")]
    out_prefix: PathBuf,

    /// STAR worker threads.
    #[arg(long, default_value_t = 24)]
    threads: usize,

    /// Inputs are plain FASTQ; omit STAR's `--readFilesCommand zcat`.
    #[arg(long)]
    plain_fastq: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
enum Status {
    Pass,
    Warn,
    Fail,
}

impl Status {
    fn marker(self) -> &'static str {
        match self {
            Self::Pass => "ok",
            Self::Warn => "warning",
            Self::Fail => "failed",
        }
    }
}

#[derive(Debug, Serialize)]
struct Check {
    id: &'static str,
    status: Status,
    summary: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    detail: Option<String>,
}

impl Check {
    fn new(id: &'static str, status: Status, summary: impl Into<String>) -> Self {
        Self {
            id,
            status,
            summary: summary.into(),
            detail: None,
        }
    }

    fn detail(mut self, detail: impl Into<String>) -> Self {
        self.detail = Some(detail.into());
        self
    }
}

#[derive(Debug, Serialize)]
struct ScanCounts {
    records: u64,
    mapped_primary_records: u64,
    secondary_records: u64,
    supplementary_records: u64,
    unmapped_records: u64,
}

#[derive(Debug, Serialize)]
struct Report {
    schema: &'static str,
    ok: bool,
    strict: bool,
    bam: String,
    whitelist: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    chemistry: Option<Chemistry>,
    counts: ScanCounts,
    checks: Vec<Check>,
}

#[derive(Default)]
struct TagStats {
    primary: u64,
    cr_missing: u64,
    ur_missing: u64,
    cy_missing: u64,
    nh_missing: u64,
    nh_invalid: u64,
    cy_wrong_length: u64,
    invalid_cr: u64,
    invalid_ur: u64,
    cr_lengths: BTreeMap<usize, u64>,
    ur_lengths: BTreeMap<usize, u64>,
}

#[derive(Default)]
struct MultimapStats {
    expected_secondaries: u64,
    observed_secondaries: u64,
    invalid_secondary_nh: u64,
    missing_names: u64,
    expected_fingerprint: [u64; 4],
    observed_fingerprint: [u64; 4],
}

fn sequence_is_acgtn(value: &[u8], allow_n: bool) -> bool {
    let mut n_bases = 0usize;
    value.iter().all(|base| match base.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => true,
        b'N' if allow_n => {
            n_bases += 1;
            n_bases <= 1
        }
        _ => false,
    })
}

fn integer(value: Value<'_>) -> Option<u64> {
    match value {
        Value::Int8(value) => u64::try_from(value).ok(),
        Value::UInt8(value) => Some(u64::from(value)),
        Value::Int16(value) => u64::try_from(value).ok(),
        Value::UInt16(value) => Some(u64::from(value)),
        Value::Int32(value) => u64::try_from(value).ok(),
        Value::UInt32(value) => Some(u64::from(value)),
        _ => None,
    }
}

fn inspect_tags(record: &bam::Record, stats: &mut TagStats) -> Result<Option<u64>> {
    let mut cr: Option<&[u8]> = None;
    let mut ur: Option<&[u8]> = None;
    let mut cy: Option<&[u8]> = None;
    let mut nh: Option<u64> = None;
    for field in record.data().iter() {
        let (tag, value) = field?;
        match (<[u8; 2]>::from(tag), value) {
            (tag, Value::String(value)) if tag == *b"CR" => cr = Some(value.as_ref()),
            (tag, Value::String(value)) if tag == *b"UR" => ur = Some(value.as_ref()),
            (tag, Value::String(value)) if tag == *b"CY" => cy = Some(value.as_ref()),
            (tag, value) if tag == *b"NH" => nh = integer(value),
            _ => {}
        }
    }
    stats.primary += 1;
    match cr {
        Some(value) => {
            *stats.cr_lengths.entry(value.len()).or_default() += 1;
            if !sequence_is_acgtn(value, true) {
                stats.invalid_cr += 1;
            }
        }
        None => stats.cr_missing += 1,
    }
    match ur {
        Some(value) => {
            *stats.ur_lengths.entry(value.len()).or_default() += 1;
            if !sequence_is_acgtn(value, false) {
                stats.invalid_ur += 1;
            }
        }
        None => stats.ur_missing += 1,
    }
    match (cr, cy) {
        (_, None) => stats.cy_missing += 1,
        (Some(cr), Some(cy)) if cr.len() != cy.len() => stats.cy_wrong_length += 1,
        _ => {}
    }
    match nh {
        Some(value) => {
            if value == 0 {
                stats.nh_invalid += 1;
            }
        }
        None => stats.nh_missing += 1,
    }
    Ok(nh)
}

fn record_nh(record: &bam::Record) -> Result<Option<u64>> {
    for field in record.data().iter() {
        let (tag, value) = field?;
        if <[u8; 2]>::from(tag) == *b"NH" {
            return Ok(integer(value));
        }
    }
    Ok(None)
}

fn alignment_fingerprint(record: &bam::Record, nh: u64) -> Option<[u64; 4]> {
    let name = record.name()?;
    let mut hasher = blake3::Hasher::new();
    hasher.update(name.as_ref());
    hasher.update(&[0]);
    hasher.update(&nh.to_le_bytes());
    let digest = hasher.finalize();
    let mut fingerprint = [0u64; 4];
    for (index, chunk) in digest.as_bytes().chunks_exact(8).enumerate() {
        fingerprint[index] = u64::from_le_bytes(chunk.try_into().unwrap());
    }
    Some(fingerprint)
}

fn add_fingerprint(total: &mut [u64; 4], value: [u64; 4], weight: u64) {
    for (total, value) in total.iter_mut().zip(value) {
        *total = total.wrapping_add(value.wrapping_mul(weight));
    }
}

fn whitelist_checks(path: &Path) -> Result<Vec<Check>> {
    let source = std::fs::read_to_string(path)
        .with_context(|| format!("reading whitelist {}", path.display()))?;
    let mut valid = 0usize;
    let mut invalid = Vec::new();
    let mut duplicates = 0usize;
    let mut seen = HashSet::new();
    for (index, line) in source.lines().enumerate() {
        let barcode = line.trim();
        if barcode.len() != 16
            || !barcode
                .bytes()
                .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
        {
            if invalid.len() < 5 {
                invalid.push(index + 1);
            }
            continue;
        }
        valid += 1;
        if !seen.insert(barcode.to_ascii_uppercase()) {
            duplicates += 1;
        }
    }
    let mut checks = Vec::new();
    if valid == 0 {
        checks.push(Check::new(
            "whitelist",
            Status::Fail,
            "whitelist contains no valid 16 bp A/C/G/T barcodes",
        ));
    } else if invalid.is_empty() {
        checks.push(Check::new(
            "whitelist",
            Status::Pass,
            format!("{valid} valid 16 bp barcodes"),
        ));
    } else {
        checks.push(
            Check::new(
                "whitelist",
                Status::Fail,
                format!("{valid} valid barcode(s), with malformed lines present"),
            )
            .detail(format!("first malformed line(s): {invalid:?}")),
        );
    }
    if duplicates > 0 {
        checks.push(Check::new(
            "whitelist_duplicates",
            Status::Warn,
            format!("{duplicates} duplicate whitelist line(s) have no effect at ingest"),
        ));
    }
    Ok(checks)
}

fn length_summary(lengths: &BTreeMap<usize, u64>) -> String {
    lengths
        .iter()
        .map(|(length, count)| format!("{length} bp: {count}"))
        .collect::<Vec<_>>()
        .join(", ")
}

fn tag_checks(stats: &TagStats, chemistry: Option<Chemistry>) -> Vec<Check> {
    if stats.primary == 0 {
        return vec![Check::new(
            "raw_tags",
            Status::Fail,
            "BAM contains no mapped primary records to ingest",
        )];
    }
    let mut checks = Vec::new();
    let required_missing = stats.cr_missing + stats.ur_missing + stats.nh_missing;
    if required_missing == 0 && stats.nh_invalid == 0 {
        checks.push(Check::new(
            "raw_tags",
            Status::Pass,
            format!(
                "CR, UR, and NH are present on all {} mapped primary records",
                stats.primary
            ),
        ));
    } else {
        checks.push(
            Check::new(
                "raw_tags",
                Status::Fail,
                "required raw barcode, UMI, or multiplicity tags are missing",
            )
            .detail(format!(
                "mapped primary records missing CR: {}; UR: {}; NH: {}; with invalid NH=0: {}",
                stats.cr_missing, stats.ur_missing, stats.nh_missing, stats.nh_invalid
            )),
        );
    }
    if stats.cy_wrong_length > 0 {
        checks.push(
            Check::new(
                "barcode_qualities",
                Status::Fail,
                format!(
                    "{} mapped primary record(s) have CY length different from CR",
                    stats.cy_wrong_length
                ),
            )
            .detail("CY must carry one barcode quality per CR base"),
        );
    } else if stats.cy_missing == 0 {
        checks.push(Check::new(
            "barcode_qualities",
            Status::Pass,
            "CY barcode qualities are present for posterior barcode correction",
        ));
    } else {
        checks.push(
            Check::new(
                "barcode_qualities",
                Status::Warn,
                format!(
                    "{} mapped primary record(s) lack CY; ingest will use its quality fallback",
                    stats.cy_missing
                ),
            )
            .detail(
                "rerun STAR with `--outSAMattributes ... CR CY UR UY` for full correction fidelity",
            ),
        );
    }
    let cr_lengths_ok =
        !stats.cr_lengths.is_empty() && stats.cr_lengths.keys().all(|length| *length == 16);
    if cr_lengths_ok {
        checks.push(Check::new(
            "barcode_length",
            Status::Pass,
            "observed CR values are 16 bp",
        ));
    } else {
        checks.push(
            Check::new(
                "barcode_length",
                Status::Fail,
                "observed CR length is incompatible with the ingest corrector",
            )
            .detail(length_summary(&stats.cr_lengths)),
        );
    }
    let expected = chemistry.map(Chemistry::umi_len);
    let umi_lengths_ok = match expected {
        Some(expected) => {
            !stats.ur_lengths.is_empty()
                && stats.ur_lengths.keys().all(|length| *length == expected)
        }
        None => {
            stats.ur_lengths.len() == 1
                && stats
                    .ur_lengths
                    .keys()
                    .next()
                    .is_some_and(|length| matches!(length, 10 | 12))
        }
    };
    if umi_lengths_ok {
        checks.push(Check::new(
            "umi_length",
            Status::Pass,
            match chemistry {
                Some(chemistry) => format!(
                    "observed UR values match {} ({} bp)",
                    chemistry.label(),
                    chemistry.umi_len()
                ),
                None => format!(
                    "one supported UR length observed ({})",
                    length_summary(&stats.ur_lengths)
                ),
            },
        ));
    } else {
        checks.push(
            Check::new(
                "umi_length",
                Status::Fail,
                "observed UR lengths do not match the selected/supported chemistry",
            )
            .detail(length_summary(&stats.ur_lengths)),
        );
    }
    if stats.invalid_cr > 0 || stats.invalid_ur > 0 {
        checks.push(
            Check::new(
                "raw_sequences",
                Status::Warn,
                "some raw barcode/UMI values contain unsupported characters and will not yield molecules",
            )
            .detail(format!(
                "CR values outside A/C/G/T plus at most one N: {}; non-ACGT UR values: {}",
                stats.invalid_cr, stats.invalid_ur
            )),
        );
    } else {
        checks.push(Check::new(
            "raw_sequences",
            Status::Pass,
            "observed raw barcode and UMI alphabets are supported",
        ));
    }
    checks
}

fn scan_bam(path: &Path, chemistry: Option<Chemistry>) -> Result<(ScanCounts, Vec<Check>)> {
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("opening BAM {}", path.display()))?;
    let header = reader
        .read_header()
        .with_context(|| format!("reading BAM header {}", path.display()))?;
    let mut counts = ScanCounts {
        records: 0,
        mapped_primary_records: 0,
        secondary_records: 0,
        supplementary_records: 0,
        unmapped_records: 0,
    };
    let mut stats = TagStats::default();
    let mut record = bam::Record::default();
    let mut last_coordinate: Option<(usize, usize)> = None;
    let mut saw_unmapped = false;
    let mut out_of_order = 0u64;
    let mut multimaps = MultimapStats::default();
    while reader.read_record(&mut record)? != 0 {
        counts.records += 1;
        let flags = Flags::from(record.flags().bits());
        if flags.is_unmapped() {
            counts.unmapped_records += 1;
            saw_unmapped = true;
            continue;
        }
        let coordinate = record
            .reference_sequence_id()
            .transpose()?
            .zip(record.alignment_start().transpose()?)
            .map(|(reference, position)| (reference, usize::from(position)));
        if let Some(coordinate) = coordinate {
            if saw_unmapped || last_coordinate.is_some_and(|last| coordinate < last) {
                out_of_order += 1;
            }
            last_coordinate = Some(coordinate);
        } else {
            out_of_order += 1;
        }
        if flags.is_supplementary() {
            counts.supplementary_records += 1;
            continue;
        }
        if flags.is_secondary() {
            counts.secondary_records += 1;
            match record_nh(&record)? {
                Some(nh) if nh > 1 => {
                    multimaps.observed_secondaries += 1;
                    if let Some(fingerprint) = alignment_fingerprint(&record, nh) {
                        add_fingerprint(&mut multimaps.observed_fingerprint, fingerprint, 1);
                    } else {
                        multimaps.missing_names += 1;
                    }
                }
                _ => multimaps.invalid_secondary_nh += 1,
            }
            continue;
        }
        counts.mapped_primary_records += 1;
        if let Some(nh) = inspect_tags(&record, &mut stats)? {
            if nh > 1 {
                let expected = nh - 1;
                multimaps.expected_secondaries =
                    multimaps.expected_secondaries.saturating_add(expected);
                if let Some(fingerprint) = alignment_fingerprint(&record, nh) {
                    add_fingerprint(&mut multimaps.expected_fingerprint, fingerprint, expected);
                } else {
                    multimaps.missing_names += 1;
                }
            }
        }
    }

    let mut checks = Vec::new();
    if header.reference_sequences().is_empty() {
        checks.push(Check::new(
            "reference_dictionary",
            Status::Fail,
            "BAM header has no reference sequences",
        ));
    } else {
        checks.push(Check::new(
            "reference_dictionary",
            Status::Pass,
            format!(
                "BAM header declares {} reference sequence(s)",
                header.reference_sequences().len()
            ),
        ));
    }
    if out_of_order == 0 {
        checks.push(Check::new(
            "coordinate_order",
            Status::Pass,
            format!(
                "all {} BAM records were scanned in coordinate order",
                counts.records
            ),
        ));
    } else {
        checks.push(Check::new(
            "coordinate_order",
            Status::Fail,
            format!("{out_of_order} record(s) violate coordinate order"),
        ));
    }
    checks.extend(tag_checks(&stats, chemistry));
    if multimaps.invalid_secondary_nh > 0 || multimaps.missing_names > 0 {
        checks.push(
            Check::new(
                "secondary_alignments",
                Status::Fail,
                "secondary alignment identity or multiplicity metadata is incomplete",
            )
            .detail(format!(
                "secondary records without numeric NH>1: {}; multimapping records without read names: {}",
                multimaps.invalid_secondary_nh, multimaps.missing_names
            )),
        );
    } else if multimaps.expected_secondaries > 0
        && multimaps.expected_secondaries == multimaps.observed_secondaries
        && multimaps.expected_fingerprint == multimaps.observed_fingerprint
    {
        checks.push(Check::new(
            "secondary_alignments",
            Status::Pass,
            format!(
                "all {} NH-declared secondary alignment(s) are retained; read-name/multiplicity fingerprints match",
                multimaps.observed_secondaries
            ),
        ));
    } else if multimaps.expected_secondaries > 0 || multimaps.observed_secondaries > 0 {
        checks.push(
            Check::new(
                "secondary_alignments",
                Status::Fail,
                "secondary alignments do not match the NH-declared primary alignments",
            )
            .detail(format!(
                "expected {} secondary record(s) from primaries; observed {} with numeric NH>1",
                multimaps.expected_secondaries, multimaps.observed_secondaries
            )),
        );
    } else {
        checks.push(
            Check::new(
                "secondary_alignments",
                Status::Warn,
                "no multimapping reads were observed, so secondary retention could not be exercised",
            )
            .detail("this is valid for a uniquely mapping fixture, but unusual for a complete scRNA-seq run"),
        );
    }
    Ok((counts, checks))
}

fn emit_report(report: &Report, format: ReportFormat) -> Result<()> {
    let stdout = std::io::stdout();
    let mut writer = stdout.lock();
    match format {
        ReportFormat::Json => {
            serde_json::to_writer_pretty(&mut writer, report)?;
            writeln!(writer)?;
        }
        ReportFormat::Text => {
            writeln!(
                writer,
                "ingest preflight: {}",
                if report.ok { "ready" } else { "not ready" }
            )?;
            writeln!(writer, "BAM: {}", report.bam)?;
            writeln!(writer, "whitelist: {}", report.whitelist)?;
            if let Some(chemistry) = report.chemistry {
                writeln!(writer, "chemistry: {}", chemistry.label())?;
            }
            for check in &report.checks {
                writeln!(
                    writer,
                    "[{}] {}: {}",
                    check.status.marker(),
                    check.id,
                    check.summary
                )?;
                if let Some(detail) = &check.detail {
                    writeln!(writer, "    {detail}")?;
                }
            }
        }
    }
    Ok(())
}

fn run_check(args: CheckArgs) -> Result<()> {
    let mut checks = whitelist_checks(&args.whitelist)?;
    let (counts, bam_checks) = scan_bam(&args.bam, args.chemistry)?;
    checks.extend(bam_checks);
    let failed = checks.iter().any(|check| check.status == Status::Fail);
    let warned = checks.iter().any(|check| check.status == Status::Warn);
    let report = Report {
        schema: "gravlax.ingest.preflight.v1",
        ok: !failed && (!args.strict || !warned),
        strict: args.strict,
        bam: args.bam.to_string_lossy().into_owned(),
        whitelist: args.whitelist.to_string_lossy().into_owned(),
        chemistry: args.chemistry,
        counts,
        checks,
    };
    emit_report(&report, args.format)?;
    if !report.ok {
        bail!(
            "ingest preflight failed; fix the reported inputs before running `aie ingest-archive`"
        );
    }
    Ok(())
}

fn shell_quote(path: &Path) -> String {
    let value = path.to_string_lossy();
    if !value.is_empty()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'/' | b'.' | b'_' | b'-'))
    {
        return value.into_owned();
    }
    format!("'{}'", value.replace('\'', "'\\''"))
}

fn run_recipe(args: RecipeArgs) -> Result<()> {
    if args.threads == 0 {
        bail!("--threads must be at least 1");
    }
    let whitelist = args
        .whitelist
        .unwrap_or_else(|| PathBuf::from(args.chemistry.default_whitelist()));
    println!(
        "# Gravlax annotation-free STAR recipe for {}",
        args.chemistry.label()
    );
    println!(
        "# Read order is cDNA R2, then barcode/UMI R1. The genome directory must have been built without a GTF or annotation-derived splice junctions."
    );
    println!("STAR \\");
    println!("  --runThreadN {} \\", args.threads);
    println!("  --genomeDir {} \\", shell_quote(&args.genome_dir));
    println!(
        "  --readFilesIn {} {} \\",
        shell_quote(&args.read2),
        shell_quote(&args.read1)
    );
    if !args.plain_fastq {
        println!("  --readFilesCommand zcat \\");
    }
    println!("  --twopassMode Basic \\");
    println!("  --soloType CB_UMI_Simple \\");
    println!("  --soloCBstart 1 \\");
    println!("  --soloCBlen 16 \\");
    println!("  --soloUMIstart 17 \\");
    println!("  --soloUMIlen {} \\", args.chemistry.umi_len());
    println!("  --soloFeatures SJ \\");
    println!("  --soloCBwhitelist {} \\", shell_quote(&whitelist));
    println!("  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\");
    println!("  --outSAMattributes NH HI AS nM CR CY UR UY \\");
    println!("  --outSAMtype BAM SortedByCoordinate \\");
    println!("  --outSAMunmapped None \\");
    println!("  --outFilterMultimapNmax 50 \\");
    println!("  --outSAMmultNmax 50 \\");
    println!("  --outFileNamePrefix {}", shell_quote(&args.out_prefix));
    println!();
    println!("# Then verify and ingest:");
    println!(
        "aie ingest check {}Aligned.sortedByCoord.out.bam --whitelist {} --chemistry {}",
        shell_quote(&args.out_prefix),
        shell_quote(&whitelist),
        args.chemistry.label()
    );
    println!(
        "aie ingest-archive {}Aligned.sortedByCoord.out.bam --whitelist {} --out sample.aie",
        shell_quote(&args.out_prefix),
        shell_quote(&whitelist)
    );
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Check(args) => run_check(args),
        Command::Recipe(args) => run_recipe(args),
    }
}
