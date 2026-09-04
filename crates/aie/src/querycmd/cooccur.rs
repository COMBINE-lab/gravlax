//! Exact Boolean predicates over one retained molecule record or one archive UMI class.
//!
//! A query always names a positive universe predicate.  That makes the candidate population and
//! its I/O route explicit, and gives `!` the narrow, reproducible meaning "not observed in this
//! archive evidence unit" rather than biological absence.

use super::{
    load_query_scope, read_junction_metadata, region_selects_chunk, uniform_query_context,
    unpack_cell_bytes, write_uniform_bundle_output, QueryAggregation, QueryScopeArgs,
    UniformQueryOutputArgs,
};
use crate::archivecmd::{decode_chunk, ChunkInfo, LazyArchive};
use crate::rows::{MolRec, PatAlt, SAME_SHAPE};
use anyhow::{bail, Context, Result};
use clap::ValueEnum;
use evidence_io::archive::Shape;
use evidence_io::terminal_tail::TerminalTailRoute;
use gravlax_output::{
    DataType, Field, OutputError, OutputFormat, RowSemantics, SelectionSummary,
    StreamingBundleWriter, TableSchema, TableSemantics,
};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use serde_json::json;
use std::collections::BTreeMap;
use std::path::Path;

const MAX_PREDICATES: usize = 64;
const NEGATIVE_TERM_SEMANTICS: &str =
    "not observed in the selected retained archive evidence unit; not biological absence";
const UMI_CLASS_IDENTITY: &str = "the default molecule_record evaluates only uniquely mapped, locus-resolved archive records for a barcode-corrected cell and exact raw UMI value; archive-wide exact raw-UMI-value-class union is explicit because same-value collisions can combine distinct physical molecules; one-mismatch UMI edges are not collapsed";

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum EvidenceUnitArg {
    /// Evaluate the expression against one serialized molecule record.
    MoleculeRecord,
    /// Merge every record carrying one barcode-corrected cell and exact raw UMI value;
    /// requires a complete archive pass and does not collapse one-mismatch UMI edges.
    UmiClass,
}

impl EvidenceUnitArg {
    fn name(self) -> &'static str {
        match self {
            Self::MoleculeRecord => "molecule_record",
            Self::UmiClass => "umi_class",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum RegionMatchArg {
    /// Match the archive anchor used by the range index.
    Anchor,
    /// Match any aligned block among the retained direct placement representatives.
    AlignedBlock,
}

impl RegionMatchArg {
    fn name(self) -> &'static str {
        match self {
            Self::Anchor => "archive_anchor",
            Self::AlignedBlock => "retained_aligned_block_overlap",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum PlacementScopeArg {
    /// Inspect only uniquely mapped chain evidence; this is the conservative default for
    /// same-molecule claims.
    Unique,
    /// Inspect unique chains plus the stored direct placement of a multimapper for junction and
    /// aligned-block predicates; archive-anchor regions remain record-anchor tests.
    Direct,
    /// Inspect every retained multimapper alternative for junction and aligned-block predicates;
    /// archive-anchor regions remain record-anchor tests. An indexed junction universe then
    /// requires --allow-full-scan because archive v2 has no alternative-placement postings.
    All,
}

impl PlacementScopeArg {
    fn name(self) -> &'static str {
        match self {
            Self::Unique => "unique_chains",
            Self::Direct => "direct",
            Self::All => "all_retained_alternatives",
        }
    }
}

#[derive(clap::Args)]
pub(crate) struct Args {
    /// Named predicate NAME=KIND:chrom:start-end[:+|-]; KIND is region, junction, or terminal.
    #[arg(long = "predicate", required = true)]
    predicates: Vec<String>,
    /// Boolean expression over predicate names using !, &, |, and parentheses.
    #[arg(long = "where")]
    expression: String,
    /// Named positive predicate defining the population in which the expression is evaluated.
    #[arg(long)]
    universe: String,
    /// Evidence identity on which predicates must co-occur.
    #[arg(long, value_enum, default_value_t = EvidenceUnitArg::MoleculeRecord)]
    unit: EvidenceUnitArg,
    /// How region predicates match retained evidence.
    #[arg(long, value_enum, default_value_t = RegionMatchArg::Anchor)]
    region_match: RegionMatchArg,
    /// Use uniquely mapped chains, or explicitly inspect a multimapper's stored direct placement
    /// or every retained alternative. Multimapper modes are not physical-molecule proof.
    #[arg(long, value_enum, default_value_t = PlacementScopeArg::Unique)]
    placements: PlacementScopeArg,
    /// Permit the exact complete pass required by UMI-class, aligned-block-universe, or
    /// all-alternative junction-universe evaluation.
    #[arg(long)]
    allow_full_scan: bool,
    /// Hard bound on routed archive chunks; exceeding it is an error, never truncation.
    #[arg(long, default_value_t = 100_000)]
    max_chunks: usize,
    /// Hard bound on decoded evidence records retained for exact evaluation.
    #[arg(long, default_value_t = 10_000_000)]
    max_evidence_records: usize,
    /// Hard bound checked against declared terminal-tail events before sparse tail decoding.
    #[arg(long, default_value_t = 10_000_000)]
    max_terminal_events: u64,
    /// Include one row per evidence unit in the named `memberships` table.
    #[arg(long)]
    emit_membership: bool,
    /// Hard membership-row bound; exceeding it is an error, never truncation.
    #[arg(long, default_value_t = 1_000_000, requires = "emit_membership")]
    max_memberships: usize,
    /// Hard bound on aggregate pattern rows; exceeding it is an error, never truncation.
    #[arg(long, default_value_t = 1_000_000)]
    max_pattern_rows: usize,
    #[command(flatten)]
    scope: QueryScopeArgs,
    #[command(flatten)]
    output: UniformQueryOutputArgs,
}

pub(super) fn validate(args: &Args) -> Result<()> {
    if args.output.format.is_none() {
        bail!("cooccur requires --format text, --format tsv, or --format json");
    }
    super::validate_uniform_output_flags(&args.output, false, false)?;
    if args.max_chunks == 0 {
        bail!("--max-chunks must be at least 1");
    }
    if args.max_memberships == 0 {
        bail!("--max-memberships must be at least 1");
    }
    if args.max_evidence_records == 0 {
        bail!("--max-evidence-records must be at least 1");
    }
    if args.max_pattern_rows == 0 {
        bail!("--max-pattern-rows must be at least 1");
    }
    if args.unit == EvidenceUnitArg::UmiClass && !args.allow_full_scan {
        bail!("--unit umi-class requires --allow-full-scan for exact cross-record evaluation");
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PredicateKind {
    Region,
    Junction,
    Terminal,
}

impl PredicateKind {
    fn name(self) -> &'static str {
        match self {
            Self::Region => "region",
            Self::Junction => "junction",
            Self::Terminal => "terminal",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Predicate {
    name: String,
    kind: PredicateKind,
    locus: String,
    chrom: String,
    chrom_id: u32,
    start: u32,
    end: u32,
    strand_rev: Option<bool>,
}

fn valid_name(name: &str) -> bool {
    let mut chars = name.chars();
    chars
        .next()
        .is_some_and(|first| first == '_' || first.is_ascii_alphabetic())
        && chars.all(|character| {
            character == '_'
                || character == '-'
                || character == '.'
                || character.is_ascii_alphanumeric()
        })
}

fn parse_stranded_locus(value: &str) -> Result<(String, u32, u32, Option<bool>)> {
    let (locus, strand_rev) = if let Some(locus) = value.strip_suffix(":+") {
        (locus, Some(false))
    } else if let Some(locus) = value.strip_suffix(":-") {
        (locus, Some(true))
    } else {
        (value, None)
    };
    let (chrom, start, end) = super::parse_locus(locus)?;
    Ok((chrom, start, end, strand_rev))
}

fn parse_predicates(values: &[String], chrom_names: &[String]) -> Result<Vec<Predicate>> {
    if values.len() > MAX_PREDICATES {
        bail!("cooccur supports at most {MAX_PREDICATES} predicates");
    }
    let mut seen = FxHashSet::default();
    values
        .iter()
        .map(|value| -> Result<Predicate> {
            let (name, descriptor) = value
                .split_once('=')
                .context("predicate must be NAME=KIND:chrom:start-end[:+|-]")?;
            if !valid_name(name) {
                bail!(
                    "predicate name {name:?} must start with an ASCII letter or '_' and contain only letters, digits, '_', '-', or '.'"
                );
            }
            if !seen.insert(name.to_owned()) {
                bail!("duplicate predicate name {name}");
            }
            let (kind, locus) = descriptor
                .split_once(':')
                .context("predicate must include region:, junction:, or terminal:")?;
            let kind = match kind {
                "region" => PredicateKind::Region,
                "junction" => PredicateKind::Junction,
                "terminal" | "terminal-tail" => PredicateKind::Terminal,
                _ => bail!("predicate {name} has unknown kind {kind:?}"),
            };
            let (chrom, start, end, strand_rev) = parse_stranded_locus(locus)
                .with_context(|| format!("invalid predicate {name}"))?;
            let chrom_id = chrom_names
                .iter()
                .position(|candidate| candidate == &chrom)
                .map(|index| index as u32)
                .with_context(|| format!("predicate {name} names unknown chromosome {chrom}"))?;
            Ok(Predicate {
                name: name.to_owned(),
                kind,
                locus: locus.to_owned(),
                chrom,
                chrom_id,
                start,
                end,
                strand_rev,
            })
        })
        .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Expression {
    Predicate(usize),
    Not(Box<Expression>),
    And(Box<Expression>, Box<Expression>),
    Or(Box<Expression>, Box<Expression>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum TruthValue {
    False,
    Unknown,
    True,
}

impl TruthValue {
    fn name(self) -> &'static str {
        match self {
            Self::False => "false",
            Self::Unknown => "unknown",
            Self::True => "true",
        }
    }

    fn selected(self) -> Option<bool> {
        match self {
            Self::False => Some(false),
            Self::Unknown => None,
            Self::True => Some(true),
        }
    }

    fn not(self) -> Self {
        match self {
            Self::False => Self::True,
            Self::Unknown => Self::Unknown,
            Self::True => Self::False,
        }
    }

    fn and(self, other: Self) -> Self {
        match (self, other) {
            (Self::False, _) | (_, Self::False) => Self::False,
            (Self::True, Self::True) => Self::True,
            _ => Self::Unknown,
        }
    }

    fn or(self, other: Self) -> Self {
        match (self, other) {
            (Self::True, _) | (_, Self::True) => Self::True,
            (Self::False, Self::False) => Self::False,
            _ => Self::Unknown,
        }
    }
}

impl Expression {
    fn evaluate(&self, observed_mask: u64, completeness_mask: u64) -> TruthValue {
        match self {
            Self::Predicate(index) => {
                let bit = 1_u64 << index;
                if observed_mask & bit != 0 {
                    TruthValue::True
                } else if completeness_mask & bit != 0 {
                    TruthValue::False
                } else {
                    TruthValue::Unknown
                }
            }
            Self::Not(value) => value.evaluate(observed_mask, completeness_mask).not(),
            Self::And(left, right) => left
                .evaluate(observed_mask, completeness_mask)
                .and(right.evaluate(observed_mask, completeness_mask)),
            Self::Or(left, right) => left
                .evaluate(observed_mask, completeness_mask)
                .or(right.evaluate(observed_mask, completeness_mask)),
        }
    }
}

struct ExpressionParser<'a> {
    source: &'a str,
    position: usize,
    names: &'a FxHashMap<String, usize>,
}

impl<'a> ExpressionParser<'a> {
    fn new(source: &'a str, names: &'a FxHashMap<String, usize>) -> Self {
        Self {
            source,
            position: 0,
            names,
        }
    }

    fn parse(mut self) -> Result<Expression> {
        let expression = self.parse_or()?;
        self.skip_space();
        if self.position != self.source.len() {
            bail!("unexpected token in --where at byte {}", self.position);
        }
        Ok(expression)
    }

    fn parse_or(&mut self) -> Result<Expression> {
        let mut expression = self.parse_and()?;
        loop {
            self.skip_space();
            if !self.consume('|') {
                break;
            }
            expression = Expression::Or(Box::new(expression), Box::new(self.parse_and()?));
        }
        Ok(expression)
    }

    fn parse_and(&mut self) -> Result<Expression> {
        let mut expression = self.parse_unary()?;
        loop {
            self.skip_space();
            if !self.consume('&') {
                break;
            }
            expression = Expression::And(Box::new(expression), Box::new(self.parse_unary()?));
        }
        Ok(expression)
    }

    fn parse_unary(&mut self) -> Result<Expression> {
        self.skip_space();
        if self.consume('!') {
            return Ok(Expression::Not(Box::new(self.parse_unary()?)));
        }
        if self.consume('(') {
            let expression = self.parse_or()?;
            self.skip_space();
            if !self.consume(')') {
                bail!("unclosed '(' in --where");
            }
            return Ok(expression);
        }
        let start = self.position;
        while self.peek().is_some_and(valid_expression_name_character) {
            self.position += 1;
        }
        if start == self.position {
            bail!("expected predicate name in --where at byte {start}");
        }
        let name = &self.source[start..self.position];
        let index = self
            .names
            .get(name)
            .copied()
            .with_context(|| format!("--where references unknown predicate {name}"))?;
        Ok(Expression::Predicate(index))
    }

    fn skip_space(&mut self) {
        while self.peek().is_some_and(char::is_whitespace) {
            self.position += 1;
        }
    }

    fn peek(&self) -> Option<char> {
        self.source
            .as_bytes()
            .get(self.position)
            .copied()
            .map(char::from)
    }

    fn consume(&mut self, expected: char) -> bool {
        if self.peek() == Some(expected) {
            self.position += 1;
            true
        } else {
            false
        }
    }
}

fn valid_expression_name_character(character: char) -> bool {
    character == '_' || character == '-' || character == '.' || character.is_ascii_alphanumeric()
}

#[derive(Clone, Copy, Debug)]
struct TailObservation {
    local_record: u32,
    chrom_id: u32,
    cleavage_anchor: u32,
    strand_rev: bool,
}

struct TerminalEvidence {
    available: bool,
    routes: Vec<TerminalTailRoute>,
}

impl TerminalEvidence {
    fn open(la: &mut LazyArchive, required: bool, max_events: u64) -> Result<Self> {
        let capability = la.terminal_tail_capability();
        let available = capability.is_some();
        if !available {
            if required {
                bail!(
                    "terminal predicates require the lossless terminal-tail capability; this archive does not declare meta.terminal_tail"
                );
            }
            return Ok(Self {
                available: false,
                routes: Vec::new(),
            });
        }
        if !required {
            return Ok(Self {
                available: true,
                routes: Vec::new(),
            });
        }
        let declared_events = capability
            .expect("terminal-tail availability was checked above")
            .events;
        if declared_events > max_events {
            bail!(
                "terminal predicates require decoding {declared_events} declared events, exceeding --max-terminal-events {max_events}; raise the explicit bound"
            );
        }
        let routes = la
            .terminal_tail_routes()?
            .context("terminal-tail capability disappeared while loading its sparse routes")?
            .to_vec();
        Ok(Self {
            available: true,
            routes,
        })
    }

    fn routed_chunks(&self, predicate: &Predicate) -> Result<Vec<usize>> {
        Ok(self
            .routes
            .iter()
            .filter(|route| {
                route.chrom == predicate.chrom_id
                    && route.min_anchor < predicate.end
                    && route.max_anchor >= predicate.start
            })
            .map(|route| route.chunk as usize)
            .collect())
    }

    fn read_chunk(
        &self,
        la: &mut LazyArchive,
        chunk: usize,
        info: &ChunkInfo,
        molecule_base: u64,
        molecules: &[MolRec],
    ) -> Result<Vec<TailObservation>> {
        let Some(route) = self
            .routes
            .iter()
            .find(|route| route.chunk as usize == chunk)
        else {
            return Ok(Vec::new());
        };
        la.terminal_tail_records(*route, info, molecule_base, molecules)?
            .into_iter()
            .map(|record| {
                if record.chunk as usize != chunk {
                    bail!("terminal-tail API returned an event for the wrong chunk");
                }
                Ok(TailObservation {
                    local_record: record.local_molecule_ordinal,
                    chrom_id: record.chrom,
                    cleavage_anchor: record.anchor,
                    strand_rev: record.strand_rev,
                })
            })
            .collect()
    }
}

fn strand_matches(requested: Option<bool>, observed: bool) -> bool {
    requested.is_none_or(|strand| strand == observed)
}

fn interval_contains(predicate: &Predicate, coordinate: u32) -> bool {
    coordinate >= predicate.start && coordinate < predicate.end
}

struct MatchContext<'a> {
    predicates: &'a [Predicate],
    region_match: RegionMatchArg,
    placements: PlacementScopeArg,
    shapes: &'a [Shape],
    patterns: Option<&'a [Vec<PatAlt>]>,
}

impl MatchContext<'_> {
    fn masks(&self, molecule: &MolRec, tails: &[TailObservation]) -> Result<(u64, u64)> {
        let mut observed = 0_u64;
        let mut complete = 0_u64;
        for (index, predicate) in self.predicates.iter().enumerate() {
            let bit = 1_u64 << index;
            if self.matches(predicate, molecule, tails)? {
                observed |= bit;
            }
            if self.absence_is_complete(predicate, molecule) {
                complete |= bit;
            }
        }
        Ok((observed, complete))
    }

    fn absence_is_complete(&self, predicate: &Predicate, molecule: &MolRec) -> bool {
        if predicate.kind != PredicateKind::Region
            || self.region_match != RegionMatchArg::AlignedBlock
            || predicate.chrom_id != molecule.chrom
            || !strand_matches(predicate.strand_rev, molecule.strand_rev)
        {
            return true;
        }
        molecule
            .chains
            .iter()
            .all(|chain| !(chain.reps.len() == 2 && chain.weight > 2))
    }

    fn matches(
        &self,
        predicate: &Predicate,
        molecule: &MolRec,
        tails: &[TailObservation],
    ) -> Result<bool> {
        if self.placements == PlacementScopeArg::Unique && molecule.chains.is_empty() {
            return Ok(false);
        }
        match predicate.kind {
            PredicateKind::Region if self.region_match == RegionMatchArg::Anchor => {
                let anchor = if self.placements == PlacementScopeArg::Unique {
                    molecule
                        .chains
                        .iter()
                        .flat_map(|chain| chain.reps.iter().map(|&(position, _)| position))
                        .min()
                        .context("unique-chain molecule record has no retained representative")?
                } else {
                    molecule.anchor()
                };
                Ok(molecule.chrom == predicate.chrom_id
                    && strand_matches(predicate.strand_rev, molecule.strand_rev)
                    && interval_contains(predicate, anchor))
            }
            PredicateKind::Terminal => Ok(tails.iter().any(|tail| {
                tail.chrom_id == predicate.chrom_id
                    && strand_matches(predicate.strand_rev, tail.strand_rev)
                    && interval_contains(predicate, tail.cleavage_anchor)
            })),
            PredicateKind::Region | PredicateKind::Junction => {
                let mut matched = false;
                let mut inspect = |chrom: u32,
                                   position: u32,
                                   strand_rev: bool,
                                   shape_id: u32|
                 -> Result<()> {
                    if matched
                        || chrom != predicate.chrom_id
                        || !strand_matches(predicate.strand_rev, strand_rev)
                    {
                        return Ok(());
                    }
                    let shape = self
                        .shapes
                        .get(shape_id as usize)
                        .with_context(|| format!("molecule references missing shape {shape_id}"))?;
                    matched = match predicate.kind {
                        PredicateKind::Region => {
                            let mut overlaps = false;
                            for &(offset, length) in &shape.blocks {
                                let start = position.checked_add(offset).with_context(|| {
                                        format!(
                                            "shape {shape_id} block start overflows genomic coordinates"
                                        )
                                    })?;
                                let end = start.checked_add(length).with_context(|| {
                                    format!(
                                        "shape {shape_id} block end overflows genomic coordinates"
                                    )
                                })?;
                                overlaps |= start < predicate.end && end > predicate.start;
                            }
                            overlaps
                        }
                        PredicateKind::Junction => {
                            let mut contains = false;
                            for blocks in shape.blocks.windows(2) {
                                let donor = position
                                        .checked_add(blocks[0].0)
                                        .and_then(|start| start.checked_add(blocks[0].1))
                                        .with_context(|| {
                                            format!(
                                                "shape {shape_id} junction donor overflows genomic coordinates"
                                            )
                                        })?;
                                let acceptor = position.checked_add(blocks[1].0).with_context(
                                        || {
                                            format!(
                                                "shape {shape_id} junction acceptor overflows genomic coordinates"
                                            )
                                        },
                                    )?;
                                contains |= donor == predicate.start && acceptor == predicate.end;
                            }
                            contains
                        }
                        PredicateKind::Terminal => unreachable!(),
                    };
                    Ok(())
                };
                for chain in &molecule.chains {
                    for &(position, shape) in &chain.reps {
                        inspect(molecule.chrom, position, molecule.strand_rev, shape)?;
                    }
                }
                if self.placements == PlacementScopeArg::Unique {
                    return Ok(matched);
                }
                for &(position, shape, pattern, _) in &molecule.mms {
                    if self.placements == PlacementScopeArg::Direct {
                        inspect(molecule.chrom, position, molecule.strand_rev, shape)?;
                        continue;
                    }
                    let alternatives = self
                        .patterns
                        .context("all-alternative matching requires the pattern dictionary")?
                        .get(pattern as usize)
                        .with_context(|| {
                            format!("molecule references missing pattern {pattern}")
                        })?;
                    for alternative in alternatives {
                        let alternative_position = u32::try_from(
                            i64::from(position)
                                .checked_add(alternative.offset)
                                .context("multimapper alternative position overflow")?,
                        )
                        .context("multimapper alternative position is negative")?;
                        let alternative_shape = if alternative.shape == SAME_SHAPE {
                            shape
                        } else {
                            alternative.shape
                        };
                        inspect(
                            alternative.chrom,
                            alternative_position,
                            molecule.strand_rev != alternative.strand_flip,
                            alternative_shape,
                        )?;
                    }
                }
                Ok(matched)
            }
        }
    }
}

#[derive(Clone, Debug)]
struct RecordHit {
    chunk: u32,
    local_record: u32,
    global_record: u64,
    class: u32,
    mask: u64,
    completeness_mask: u64,
}

struct ScanRequest<'a> {
    chunks: &'a [ChunkInfo],
    selected_chunks: &'a [usize],
    global_offsets: &'a [u64],
    predicates: &'a [Predicate],
    candidate_classes: Option<&'a FxHashSet<u32>>,
    region_match: RegionMatchArg,
    placements: PlacementScopeArg,
    shapes: &'a [Shape],
    patterns: Option<&'a [Vec<PatAlt>]>,
    terminal: &'a TerminalEvidence,
    max_records: usize,
}

fn scan_records(la: &mut LazyArchive, request: &ScanRequest<'_>) -> Result<Vec<RecordHit>> {
    let decode_window = rayon::current_num_threads().saturating_mul(2).max(1);
    let mut hits = Vec::new();
    for selected_window in request.selected_chunks.chunks(decode_window) {
        let decoded_chunks: Vec<(usize, Vec<MolRec>)> = {
            let (reader, tables) = la.reader_and_tables();
            let reader = &*reader;
            selected_window
                .par_iter()
                .map(|&chunk_index| -> Result<(usize, Vec<MolRec>)> {
                    let (compressed, raw_len) =
                        reader.read_compressed_at(&format!("c{chunk_index}"))?;
                    let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                    let molecules = decode_chunk(&raw, &request.chunks[chunk_index], None, tables)?;
                    Ok((chunk_index, molecules))
                })
                .collect::<Result<_>>()?
        };
        for (chunk_index, molecules) in decoded_chunks {
            let tails = request.terminal.read_chunk(
                la,
                chunk_index,
                &request.chunks[chunk_index],
                request.global_offsets[chunk_index],
                &molecules,
            )?;
            let mut tails_by_record: FxHashMap<u32, Vec<TailObservation>> = FxHashMap::default();
            for tail in tails {
                tails_by_record
                    .entry(tail.local_record)
                    .or_default()
                    .push(tail);
            }
            let matcher = MatchContext {
                predicates: request.predicates,
                region_match: request.region_match,
                placements: request.placements,
                shapes: request.shapes,
                patterns: request.patterns,
            };
            let mut chunk_hits: Vec<RecordHit> = molecules
                .iter()
                .enumerate()
                .filter(|(_, molecule)| {
                    (request.placements != PlacementScopeArg::Unique || !molecule.chains.is_empty())
                        && request
                            .candidate_classes
                            .is_none_or(|classes| classes.contains(&molecule.umi_class))
                })
                .map(|(local_record, molecule)| -> Result<RecordHit> {
                    let local_record = u32::try_from(local_record)
                        .context("chunk molecule ordinal exceeds u32")?;
                    let (mask, completeness_mask) = matcher.masks(
                        molecule,
                        tails_by_record
                            .get(&local_record)
                            .map(Vec::as_slice)
                            .unwrap_or_default(),
                    )?;
                    Ok(RecordHit {
                        chunk: chunk_index as u32,
                        local_record,
                        global_record: request.global_offsets[chunk_index]
                            .checked_add(u64::from(local_record))
                            .context("global molecule-record ordinal overflow")?,
                        class: molecule.umi_class,
                        mask,
                        completeness_mask,
                    })
                })
                .collect::<Result<_>>()?;
            hits.append(&mut chunk_hits);
            if hits.len() > request.max_records {
                bail!(
                    "cooccur retained more than {} evidence records, exceeding --max-evidence-records {}; narrow --universe or raise the explicit bound",
                    hits.len(),
                    request.max_records
                );
            }
        }
    }
    Ok(hits)
}

fn global_record_offsets(chunks: &[ChunkInfo]) -> Result<Vec<u64>> {
    let mut next = 0_u64;
    chunks
        .iter()
        .map(|chunk| {
            let offset = next;
            next = next
                .checked_add(u64::from(chunk.n_mols))
                .context("archive molecule-record count overflow")?;
            Ok(offset)
        })
        .collect()
}

fn full_scan_reason(args: &Args, universe: &Predicate) -> Option<&'static str> {
    if args.unit == EvidenceUnitArg::UmiClass {
        Some("exact raw-UMI-value-class union across all molecule records")
    } else if universe.kind == PredicateKind::Region
        && args.region_match == RegionMatchArg::AlignedBlock
    {
        Some("aligned-block universe has no complete archive posting index")
    } else if universe.kind == PredicateKind::Junction && args.placements == PlacementScopeArg::All
    {
        Some("multimapper alternative junctions have no complete archive posting index")
    } else {
        None
    }
}

fn universe_route_full_scan_reason(args: &Args, universe: &Predicate) -> Option<&'static str> {
    if universe.kind == PredicateKind::Region && args.region_match == RegionMatchArg::AlignedBlock {
        Some("aligned-block universe has no complete archive posting index")
    } else if universe.kind == PredicateKind::Junction && args.placements == PlacementScopeArg::All
    {
        Some("multimapper alternative junctions have no complete archive posting index")
    } else {
        None
    }
}

fn route_universe(
    args: &Args,
    universe: &Predicate,
    la: &mut LazyArchive,
    chunks: &[ChunkInfo],
    terminal: &TerminalEvidence,
) -> Result<(Vec<usize>, Option<&'static str>)> {
    let reason = universe_route_full_scan_reason(args, universe);
    if !args.allow_full_scan {
        if let Some(message) = reason {
            bail!("{message}; rerun with --allow-full-scan");
        }
    }
    let mut selected = if reason.is_some() {
        (0..chunks.len()).collect::<Vec<_>>()
    } else {
        match universe.kind {
            PredicateKind::Region => chunks
                .iter()
                .enumerate()
                .filter(|(_, chunk)| {
                    region_selects_chunk(chunk, universe.chrom_id, universe.start, universe.end)
                })
                .map(|(index, _)| index)
                .collect(),
            PredicateKind::Junction => {
                let metadata = read_junction_metadata(la)?;
                metadata
                    .iter()
                    .find(|row| {
                        row.chrom == universe.chrom_id
                            && row.donor == universe.start
                            && row.acceptor == universe.end
                    })
                    .map(|row| row.posts.iter().map(|post| *post as usize).collect())
                    .unwrap_or_default()
            }
            PredicateKind::Terminal => terminal.routed_chunks(universe)?,
        }
    };
    selected.sort_unstable();
    selected.dedup();
    if let Some(&bad) = selected.iter().find(|&&index| index >= chunks.len()) {
        bail!("predicate route references missing archive chunk {bad}");
    }
    if selected.len() > args.max_chunks {
        bail!(
            "cooccur route selects {} chunks, exceeding --max-chunks {}; narrow --universe or raise the explicit bound",
            selected.len(),
            args.max_chunks
        );
    }
    Ok((selected, reason))
}

#[derive(Clone, Debug)]
struct EvaluatedUnit {
    cell: u32,
    class: u32,
    chunk: Option<u32>,
    local_record: Option<u32>,
    global_record: Option<u64>,
    contributing_records: u64,
    mask: u64,
    completeness_mask: u64,
    selection: TruthValue,
}

fn evaluate_record_units(
    records: Vec<RecordHit>,
    universe_mask: u64,
    expression: &Expression,
    la: &mut LazyArchive,
) -> Result<Vec<EvaluatedUnit>> {
    let records: Vec<RecordHit> = records
        .into_iter()
        .filter(|record| record.mask & universe_mask != 0)
        .collect();
    la.prefetch_coc(records.iter().map(|record| record.class))?;
    records
        .into_iter()
        .map(|record| {
            Ok(EvaluatedUnit {
                cell: la.cell_of_cached(record.class)?,
                class: record.class,
                chunk: Some(record.chunk),
                local_record: Some(record.local_record),
                global_record: Some(record.global_record),
                contributing_records: 1,
                mask: record.mask,
                completeness_mask: record.completeness_mask,
                selection: expression.evaluate(record.mask, record.completeness_mask),
            })
        })
        .collect()
}

fn evaluate_class_units(
    candidates: Vec<RecordHit>,
    all_records: Vec<RecordHit>,
    universe_mask: u64,
    expression: &Expression,
    la: &mut LazyArchive,
) -> Result<Vec<EvaluatedUnit>> {
    let merged = merge_class_evidence(&candidates, all_records, universe_mask);
    la.prefetch_coc(merged.keys().copied())?;
    merged
        .into_iter()
        .map(|(class, evidence)| {
            Ok(EvaluatedUnit {
                cell: la.cell_of_cached(class)?,
                class,
                chunk: None,
                local_record: None,
                global_record: None,
                contributing_records: evidence.records,
                mask: evidence.mask,
                completeness_mask: evidence.completeness_mask,
                selection: expression.evaluate(evidence.mask, evidence.completeness_mask),
            })
        })
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ClassEvidence {
    mask: u64,
    completeness_mask: u64,
    records: u64,
}

fn merge_class_evidence(
    candidates: &[RecordHit],
    all_records: Vec<RecordHit>,
    universe_mask: u64,
) -> FxHashMap<u32, ClassEvidence> {
    let candidate_classes: FxHashSet<u32> = candidates
        .iter()
        .filter(|record| record.mask & universe_mask != 0)
        .map(|record| record.class)
        .collect();
    let mut merged: FxHashMap<u32, ClassEvidence> = FxHashMap::default();
    for record in all_records {
        if candidate_classes.contains(&record.class) {
            let aggregate = merged.entry(record.class).or_insert(ClassEvidence {
                mask: 0,
                completeness_mask: !0_u64,
                records: 0,
            });
            aggregate.mask |= record.mask;
            aggregate.completeness_mask &= record.completeness_mask;
            aggregate.records += 1;
        }
    }
    merged
}

#[derive(Clone, Debug, Default)]
struct PatternAggregate {
    units: u64,
    cells: FxHashSet<u32>,
}

#[derive(Clone, Debug)]
struct PatternRow {
    entity: String,
    mask: u64,
    completeness_mask: u64,
    selection: TruthValue,
    units: u64,
    cells: Option<u64>,
    scope_cells: Option<u64>,
}

fn aggregate_patterns(
    units: &[EvaluatedUnit],
    scope: &super::QueryScope,
    cell_dictionary: &[u32],
    max_rows: usize,
) -> Result<Vec<PatternRow>> {
    let mut aggregates: BTreeMap<(String, u64, u64, TruthValue), PatternAggregate> =
        BTreeMap::new();
    for unit in units.iter().filter(|unit| scope.includes(unit.cell)) {
        let entity = match scope.aggregation {
            QueryAggregation::Cell => {
                let barcode = unpack_cell_bytes(cell_dictionary[unit.cell as usize]);
                std::str::from_utf8(&barcode)
                    .expect("packed archive barcode decodes to ASCII")
                    .to_owned()
            }
            QueryAggregation::Group => match scope.group_of.get(&unit.cell) {
                Some(group) => scope.group_names[*group as usize].clone(),
                None => continue,
            },
            QueryAggregation::Bulk => "all".to_owned(),
        };
        let key = (entity, unit.mask, unit.completeness_mask, unit.selection);
        if !aggregates.contains_key(&key) && aggregates.len() == max_rows {
            bail!(
                "cooccur would produce more than --max-pattern-rows {max_rows}; narrow the scope or raise the explicit bound"
            );
        }
        let aggregate = aggregates.entry(key).or_default();
        aggregate.units += 1;
        aggregate.cells.insert(unit.cell);
    }
    Ok(aggregates
        .into_iter()
        .map(
            |((entity, mask, completeness_mask, selection), aggregate)| {
                let scope_cells = match scope.aggregation {
                    QueryAggregation::Group => scope
                        .group_names
                        .iter()
                        .position(|name| name == &entity)
                        .map(|group| scope.selected_per_group[group] as u64),
                    QueryAggregation::Bulk => Some(scope.selected_cells as u64),
                    QueryAggregation::Cell => None,
                };
                PatternRow {
                    entity,
                    mask,
                    completeness_mask,
                    selection,
                    units: aggregate.units,
                    cells: (scope.aggregation != QueryAggregation::Cell)
                        .then_some(aggregate.cells.len() as u64),
                    scope_cells,
                }
            },
        )
        .collect())
}

fn mask_names(mask: u64, predicates: &[Predicate]) -> Vec<&str> {
    predicates
        .iter()
        .enumerate()
        .filter(|(index, _)| mask & (1_u64 << index) != 0)
        .map(|(_, predicate)| predicate.name.as_str())
        .collect()
}

fn mask_text(mask: u64) -> String {
    format!("0x{mask:016x}")
}

#[derive(Serialize)]
struct Summary<'a> {
    coordinates: &'static str,
    evidence_unit: &'a str,
    expression: &'a str,
    universe: &'a str,
    predicates: u64,
    candidate_units: u64,
    selected_units: u64,
    indeterminate_units: u64,
    candidate_cells: u64,
    selected_cells: u64,
    pattern_rows: u64,
    chunks_read: u64,
    archive_chunks: u64,
    full_scan: bool,
    terminal_tail_available: bool,
    negative_term_semantics: &'static str,
    negative_terms_complete_within_declared_evidence_scope: bool,
    positive_matches_are_witnessed: bool,
    biological_absence_claimed: bool,
    omitted_middle_placement_semantics: &'static str,
    multimapper_alternative_semantics: &'static str,
    jointly_realizable_multimapper_placement_claimed: bool,
    umi_class_identity: &'static str,
}

pub(super) struct RunContext<'a> {
    pub(super) archive: &'a Path,
    pub(super) la: &'a mut LazyArchive,
    pub(super) chunks: &'a [ChunkInfo],
    pub(super) chrom_names: &'a [String],
    pub(super) t0: std::time::Instant,
    pub(super) t_open: f32,
}

pub(super) fn run(context: RunContext<'_>, args: Args) -> Result<()> {
    let RunContext {
        archive,
        la,
        chunks,
        chrom_names,
        t0,
        t_open,
    } = context;
    let predicates = parse_predicates(&args.predicates, chrom_names)?;
    let names: FxHashMap<String, usize> = predicates
        .iter()
        .enumerate()
        .map(|(index, predicate)| (predicate.name.clone(), index))
        .collect();
    let expression = ExpressionParser::new(&args.expression, &names).parse()?;
    let universe_index = names
        .get(&args.universe)
        .copied()
        .with_context(|| format!("--universe names unknown predicate {}", args.universe))?;
    let universe = &predicates[universe_index];
    let universe_mask = 1_u64 << universe_index;
    let terminal_required = predicates
        .iter()
        .any(|predicate| predicate.kind == PredicateKind::Terminal);
    let terminal = TerminalEvidence::open(la, terminal_required, args.max_terminal_events)?;
    let (selected_chunks, _) = route_universe(&args, universe, la, chunks, &terminal)?;
    let full_scan_reason = full_scan_reason(&args, universe);
    if full_scan_reason.is_some() && chunks.len() > args.max_chunks {
        bail!(
            "cooccur exact full scan needs {} chunks, exceeding --max-chunks {}; raise the explicit bound or choose a routed evidence unit",
            chunks.len(),
            args.max_chunks
        );
    }
    let global_offsets = global_record_offsets(chunks)?;
    let shapes = la.shapes()?;
    let patterns = (args.placements == PlacementScopeArg::All)
        .then(|| la.patterns())
        .transpose()?;
    let candidate_scan = ScanRequest {
        chunks,
        selected_chunks: &selected_chunks,
        global_offsets: &global_offsets,
        predicates: &predicates,
        candidate_classes: None,
        region_match: args.region_match,
        placements: args.placements,
        shapes: &shapes,
        patterns: patterns.as_deref().map(Vec::as_slice),
        terminal: &terminal,
        max_records: args.max_evidence_records,
    };
    let candidate_records = scan_records(la, &candidate_scan)?;
    let (mut units, chunks_read) = match args.unit {
        EvidenceUnitArg::MoleculeRecord => (
            evaluate_record_units(candidate_records, universe_mask, &expression, la)?,
            selected_chunks.len(),
        ),
        EvidenceUnitArg::UmiClass => {
            let classes: FxHashSet<u32> = candidate_records
                .iter()
                .filter(|record| record.mask & universe_mask != 0)
                .map(|record| record.class)
                .collect();
            let routed: FxHashSet<usize> = selected_chunks.iter().copied().collect();
            let remaining_chunks: Vec<usize> = (0..chunks.len())
                .filter(|chunk| !routed.contains(chunk))
                .collect();
            let mut all_records = candidate_records.clone();
            let chunks_read = if classes.is_empty() {
                selected_chunks.len()
            } else {
                let all_scan = ScanRequest {
                    chunks,
                    selected_chunks: &remaining_chunks,
                    global_offsets: &global_offsets,
                    predicates: &predicates,
                    candidate_classes: Some(&classes),
                    region_match: args.region_match,
                    placements: args.placements,
                    shapes: &shapes,
                    patterns: patterns.as_deref().map(Vec::as_slice),
                    terminal: &terminal,
                    max_records: args
                        .max_evidence_records
                        .saturating_sub(candidate_records.len()),
                };
                all_records.extend(scan_records(la, &all_scan)?);
                selected_chunks.len() + remaining_chunks.len()
            };
            (
                evaluate_class_units(
                    candidate_records,
                    all_records,
                    universe_mask,
                    &expression,
                    la,
                )?,
                chunks_read,
            )
        }
    };
    let mut scope = load_query_scope(la, &args.scope)?;
    scope.ensure_resolved_mapping_digest()?;
    units.retain(|unit| scope.includes(unit.cell));
    let cell_dictionary = la.cells()?.to_vec();
    units.sort_unstable_by_key(|unit| {
        (
            cell_dictionary[unit.cell as usize],
            unit.class,
            unit.global_record.unwrap_or_default(),
        )
    });
    let pattern_rows = aggregate_patterns(&units, &scope, &cell_dictionary, args.max_pattern_rows)?;
    if args.emit_membership && units.len() > args.max_memberships {
        bail!(
            "cooccur produced {} membership rows, exceeding --max-memberships {}; narrow the scope or raise the explicit bound",
            units.len(),
            args.max_memberships
        );
    }
    let candidate_cells: FxHashSet<u32> = units.iter().map(|unit| unit.cell).collect();
    let selected_cells: FxHashSet<u32> = units
        .iter()
        .filter(|unit| unit.selection == TruthValue::True)
        .map(|unit| unit.cell)
        .collect();
    let all_predicates_mask = if predicates.len() == 64 {
        u64::MAX
    } else {
        (1_u64 << predicates.len()) - 1
    };
    let summary = Summary {
        coordinates: "0-based half-open",
        evidence_unit: args.unit.name(),
        expression: &args.expression,
        universe: &args.universe,
        predicates: predicates.len() as u64,
        candidate_units: units.len() as u64,
        selected_units: units
            .iter()
            .filter(|unit| unit.selection == TruthValue::True)
            .count() as u64,
        indeterminate_units: units
            .iter()
            .filter(|unit| unit.selection == TruthValue::Unknown)
            .count() as u64,
        candidate_cells: candidate_cells.len() as u64,
        selected_cells: selected_cells.len() as u64,
        pattern_rows: pattern_rows.len() as u64,
        chunks_read: chunks_read as u64,
        archive_chunks: chunks.len() as u64,
        full_scan: chunks_read == chunks.len(),
        terminal_tail_available: terminal.available,
        negative_term_semantics: NEGATIVE_TERM_SEMANTICS,
        negative_terms_complete_within_declared_evidence_scope: units
            .iter()
            .all(|unit| unit.completeness_mask & all_predicates_mask == all_predicates_mask),
        positive_matches_are_witnessed: true,
        biological_absence_claimed: false,
        omitted_middle_placement_semantics: "when a junction chain represents more than two distinct reads with two span-extreme placements, an unobserved aligned-block predicate is unknown",
        multimapper_alternative_semantics: match args.placements {
            PlacementScopeArg::Unique => {
                "only uniquely mapped chain evidence is inspected; multimapper records and alternatives are excluded"
            }
            PlacementScopeArg::Direct => {
                "opt-in diagnostic mode: junction/aligned-block predicates inspect retained direct-placement children, while archive-anchor regions remain record-anchor tests; a mask asserts one evidence unit, not one physical molecule, read, or placement"
            }
            PlacementScopeArg::All => {
                "opt-in diagnostic mode: junction/aligned-block predicates are existential across retained alternatives, while archive-anchor regions remain record-anchor tests; a multi-predicate mask does not assert one physical molecule or one jointly realizable placement"
            }
        },
        jointly_realizable_multimapper_placement_claimed: false,
        umi_class_identity: UMI_CLASS_IDENTITY,
    };
    let predicate_schema = TableSchema::new(
        "gravlax.query.cooccur.predicates.v1",
        vec![
            Field::new("predicate_index", DataType::UInt64),
            Field::new("name", DataType::String),
            Field::new("predicate_kind", DataType::String),
            Field::new("locus", DataType::String),
            Field::new("chrom", DataType::String),
            Field::new("start", DataType::UInt64),
            Field::new("end", DataType::UInt64),
            Field::new("strand", DataType::String).nullable(),
            Field::new("match_semantics", DataType::String),
        ],
    )?
    .with_semantics(
        TableSemantics::new(RowSemantics::Sequence)
            .ordered_by([gravlax_output::OrderKey::ascending("predicate_index")]),
    )?;
    let pattern_schema = TableSchema::new(
        "gravlax.query.cooccur.patterns.v1",
        vec![
            Field::new("aggregation", DataType::String),
            Field::new("entity", DataType::String),
            Field::new("pattern_mask", DataType::String),
            Field::new("completeness_mask", DataType::String),
            Field::new("matched_predicates", DataType::Json),
            Field::new("selection_state", DataType::String),
            Field::new("selected", DataType::Boolean).nullable(),
            Field::new("evidence_units", DataType::UInt64),
            Field::new("cells", DataType::UInt64).nullable(),
            Field::new("scope_cells", DataType::UInt64).nullable(),
        ],
    )?
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key([
        "aggregation",
        "entity",
        "pattern_mask",
        "completeness_mask",
    ]))?;
    let membership_schema = TableSchema::new(
        "gravlax.query.cooccur.memberships.v1",
        vec![
            Field::new("cell_id", DataType::UInt64),
            Field::new("unit_id", DataType::String),
            Field::new("barcode", DataType::String),
            Field::new("umi_class", DataType::UInt64),
            Field::new("chunk", DataType::UInt64).nullable(),
            Field::new("local_record", DataType::UInt64).nullable(),
            Field::new("global_record", DataType::UInt64).nullable(),
            Field::new("contributing_records", DataType::UInt64),
            Field::new("pattern_mask", DataType::String),
            Field::new("completeness_mask", DataType::String),
            Field::new("matched_predicates", DataType::Json),
            Field::new("selection_state", DataType::String),
            Field::new("selected", DataType::Boolean).nullable(),
        ],
    )?
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["unit_id"]))?;
    let mut parameters = BTreeMap::new();
    parameters.insert("expression".into(), json!(args.expression));
    parameters.insert("universe".into(), json!(args.universe));
    parameters.insert("evidence_unit".into(), json!(args.unit.name()));
    parameters.insert("region_match".into(), json!(args.region_match.name()));
    parameters.insert("placement_scope".into(), json!(args.placements.name()));
    parameters.insert("cell_scope".into(), scope.provenance_json());
    parameters.insert("aggregation".into(), json!(scope.aggregation_name()));
    parameters.insert("max_chunks".into(), json!(args.max_chunks));
    parameters.insert(
        "max_evidence_records".into(),
        json!(args.max_evidence_records),
    );
    parameters.insert(
        "max_terminal_events".into(),
        json!(args.max_terminal_events),
    );
    parameters.insert("max_pattern_rows".into(), json!(args.max_pattern_rows));
    parameters.insert("emit_membership".into(), json!(args.emit_membership));
    parameters.insert("max_memberships".into(), json!(args.max_memberships));
    parameters.insert(
        "archive_access".into(),
        json!(full_scan_reason.unwrap_or("universe-routed archive chunks")),
    );
    let context = uniform_query_context(
        archive,
        la.reader().archive_version(),
        la.reader()
            .content_commitment()
            .map(|commitment| commitment.to_hex()),
        parameters,
    );
    let predicate_selection = SelectionSummary::complete(predicates.len() as u64);
    let pattern_selection = SelectionSummary::complete(pattern_rows.len() as u64);
    let membership_selection = SelectionSummary::complete(units.len() as u64);
    write_uniform_bundle_output(&args.output, |writer, format| {
        let mut bundle = StreamingBundleWriter::new_with_summary(
            writer,
            "gravlax.query.cooccur.result.v1",
            OutputFormat::from(format),
            &context,
            &summary,
        )?;
        bundle.write_table(
            "predicates",
            &predicate_schema,
            Some(&predicate_selection),
            |rows| {
                for (index, predicate) in predicates.iter().enumerate() {
                    rows.write_row_with(|row| {
                        row.uint64(index as u64)?;
                        row.string(&predicate.name)?;
                        row.string(predicate.kind.name())?;
                        row.string(&predicate.locus)?;
                        row.string(&predicate.chrom)?;
                        row.uint64(predicate.start as u64)?;
                        row.uint64(predicate.end as u64)?;
                        match predicate.strand_rev {
                            Some(false) => row.string("+")?,
                            Some(true) => row.string("-")?,
                            None => row.null()?,
                        }
                        row.string(match predicate.kind {
                            PredicateKind::Region => args.region_match.name(),
                            PredicateKind::Junction => "exact retained junction boundaries",
                            PredicateKind::Terminal => {
                                "lossless terminal-tail cleavage anchor in interval"
                            }
                        })?;
                        Ok(())
                    })?;
                }
                Ok(())
            },
        )?;
        bundle.write_table(
            "patterns",
            &pattern_schema,
            Some(&pattern_selection),
            |rows| {
                for pattern in &pattern_rows {
                    rows.write_row_with(|row| {
                        row.string(scope.aggregation_name())?;
                        row.string(&pattern.entity)?;
                        row.string(&mask_text(pattern.mask))?;
                        row.string(&mask_text(pattern.completeness_mask))?;
                        row.json(&json!(mask_names(pattern.mask, &predicates)))?;
                        row.string(pattern.selection.name())?;
                        match pattern.selection.selected() {
                            Some(selected) => row.boolean(selected)?,
                            None => row.null()?,
                        }
                        row.uint64(pattern.units)?;
                        match pattern.cells {
                            Some(cells) => row.uint64(cells)?,
                            None => row.null()?,
                        }
                        match pattern.scope_cells {
                            Some(cells) => row.uint64(cells)?,
                            None => row.null()?,
                        }
                        Ok(())
                    })?;
                }
                Ok(())
            },
        )?;
        if args.emit_membership {
            bundle.write_table(
                "memberships",
                &membership_schema,
                Some(&membership_selection),
                |rows| {
                    for unit in &units {
                        let barcode = unpack_cell_bytes(cell_dictionary[unit.cell as usize]);
                        let barcode = std::str::from_utf8(&barcode)
                            .expect("packed archive barcode decodes to ASCII");
                        rows.write_row_with(|row| {
                            row.uint64(unit.cell as u64)?;
                            row.string(&match unit.global_record {
                                Some(record) => format!("record:{record}"),
                                None => format!("umi-class:{}:{}", unit.cell, unit.class),
                            })?;
                            row.string(barcode)?;
                            row.uint64(unit.class as u64)?;
                            match unit.chunk {
                                Some(value) => row.uint64(value as u64)?,
                                None => row.null()?,
                            }
                            match unit.local_record {
                                Some(value) => row.uint64(value as u64)?,
                                None => row.null()?,
                            }
                            match unit.global_record {
                                Some(value) => row.uint64(value)?,
                                None => row.null()?,
                            }
                            row.uint64(unit.contributing_records)?;
                            row.string(&mask_text(unit.mask))?;
                            row.string(&mask_text(unit.completeness_mask))?;
                            row.json(&json!(mask_names(unit.mask, &predicates)))?;
                            row.string(unit.selection.name())?;
                            match unit.selection.selected() {
                                Some(selected) => row.boolean(selected)?,
                                None => row.null()?,
                            }
                            Ok(())
                        })?;
                    }
                    Ok(())
                },
            )?;
        }
        bundle.finish()?;
        Ok::<(), OutputError>(())
    })?;
    eprintln!(
        "cooccur: {} selected / {} candidate {}s; {} / {} chunks read (open {:.2}s, total {:.2}s)",
        summary.selected_units,
        summary.candidate_units,
        args.unit.name(),
        summary.chunks_read,
        summary.archive_chunks,
        t_open,
        t0.elapsed().as_secs_f32()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    fn predicates() -> Vec<Predicate> {
        parse_predicates(
            &[
                "j= junction:chr1:100-200".replace("= ", "="),
                "r=region:chr1:50-250:+".to_owned(),
                "t=terminal:chr1:300-301:-".to_owned(),
            ],
            &["chr1".to_owned()],
        )
        .unwrap()
    }

    #[test]
    fn parses_named_stranded_predicates() {
        let parsed = predicates();
        assert_eq!(parsed[0].kind, PredicateKind::Junction);
        assert_eq!(parsed[1].strand_rev, Some(false));
        assert_eq!(parsed[2].strand_rev, Some(true));
        assert_eq!((parsed[2].start, parsed[2].end), (300, 301));
    }

    #[test]
    fn boolean_precedence_and_parentheses_are_explicit() {
        let names = [
            ("a".to_owned(), 0),
            ("b".to_owned(), 1),
            ("c".to_owned(), 2),
        ]
        .into_iter()
        .collect();
        let expression = ExpressionParser::new("a | b & !c", &names).parse().unwrap();
        assert_eq!(expression.evaluate(0b001, !0), TruthValue::True);
        assert_eq!(expression.evaluate(0b010, !0), TruthValue::True);
        assert_eq!(expression.evaluate(0b110, !0), TruthValue::False);
        let grouped = ExpressionParser::new("(a | b) & !c", &names)
            .parse()
            .unwrap();
        assert_eq!(grouped.evaluate(0b101, !0), TruthValue::False);
        assert_eq!(grouped.evaluate(0b001, 0b001), TruthValue::Unknown);
        assert!(ExpressionParser::new("missing", &names).parse().is_err());
        assert!(NEGATIVE_TERM_SEMANTICS.contains("not biological absence"));
    }

    #[test]
    fn predicate_names_are_strict_and_duplicates_fail() {
        assert!(
            parse_predicates(&["1bad=region:chr1:1-2".to_owned()], &["chr1".to_owned()]).is_err()
        );
        assert!(parse_predicates(
            &[
                "x=region:chr1:1-2".to_owned(),
                "x=junction:chr1:1-2".to_owned()
            ],
            &["chr1".to_owned()]
        )
        .is_err());
    }

    #[test]
    fn interval_and_strand_matching_are_half_open() {
        let predicate = &predicates()[2];
        assert!(interval_contains(predicate, 300));
        assert!(!interval_contains(predicate, 301));
        assert!(strand_matches(predicate.strand_rev, true));
        assert!(!strand_matches(predicate.strand_rev, false));
    }

    #[test]
    fn terminal_match_requires_the_requested_chromosome() {
        let predicate = &predicates()[2];
        let molecule = MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: true,
            chains: smallvec::smallvec![crate::rows::MolChain {
                weight: 1,
                reps: smallvec::smallvec![(300, 0)],
            }],
            mms: Default::default(),
        };
        let matcher = MatchContext {
            predicates: std::slice::from_ref(predicate),
            region_match: RegionMatchArg::Anchor,
            placements: PlacementScopeArg::Unique,
            shapes: &[],
            patterns: None,
        };
        assert!(!matcher
            .matches(
                predicate,
                &molecule,
                &[TailObservation {
                    local_record: 0,
                    chrom_id: 1,
                    cleavage_anchor: 300,
                    strand_rev: true,
                }],
            )
            .unwrap());
        assert!(matcher
            .matches(
                predicate,
                &molecule,
                &[TailObservation {
                    local_record: 0,
                    chrom_id: 0,
                    cleavage_anchor: 300,
                    strand_rev: true,
                }],
            )
            .unwrap());
    }

    #[test]
    fn omitted_middle_read_makes_aligned_block_absence_unknown() {
        let predicate = parse_predicates(
            &["r=region:chr1:500-501:+".to_owned()],
            &["chr1".to_owned()],
        )
        .unwrap()
        .remove(0);
        let mut molecule = MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: false,
            chains: Default::default(),
            mms: Default::default(),
        };
        let mut reps = smallvec::SmallVec::new();
        reps.push((100, 0));
        reps.push((900, 0));
        molecule
            .chains
            .push(crate::rows::MolChain { weight: 3, reps });
        let matcher = MatchContext {
            predicates: std::slice::from_ref(&predicate),
            region_match: RegionMatchArg::AlignedBlock,
            placements: PlacementScopeArg::Unique,
            shapes: &[Shape {
                blocks: vec![(0, 10)],
            }],
            patterns: None,
        };
        let (observed, complete) = matcher.masks(&molecule, &[]).unwrap();
        assert_eq!(observed, 0);
        assert_eq!(complete, 0);
        let expression = Expression::Not(Box::new(Expression::Predicate(0)));
        assert_eq!(expression.evaluate(observed, complete), TruthValue::Unknown);
    }

    #[test]
    fn malformed_shape_coordinate_overflow_fails_closed() {
        let predicate = parse_predicates(&["r=region:chr1:1-2:+".to_owned()], &["chr1".to_owned()])
            .unwrap()
            .remove(0);
        let molecule = MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: false,
            chains: smallvec::smallvec![crate::rows::MolChain {
                weight: 1,
                reps: smallvec::smallvec![(1, 0)],
            }],
            mms: Default::default(),
        };
        let matcher = MatchContext {
            predicates: std::slice::from_ref(&predicate),
            region_match: RegionMatchArg::AlignedBlock,
            placements: PlacementScopeArg::Unique,
            shapes: &[Shape {
                blocks: vec![(u32::MAX, 1)],
            }],
            patterns: None,
        };
        let error = matcher.matches(&predicate, &molecule, &[]).unwrap_err();
        assert!(error.to_string().contains("overflows genomic coordinates"));
    }

    #[test]
    fn molecule_record_is_the_conservative_cli_default() {
        let cli = crate::Cli::try_parse_from([
            "aie",
            "query",
            "sample.aie",
            "cooccur",
            "--predicate",
            "u=region:chr1:0-1",
            "--where",
            "u",
            "--universe",
            "u",
            "--format",
            "json",
        ])
        .unwrap();
        let crate::Cmd::Query(query) = cli.cmd else {
            panic!("parsed the wrong top-level command")
        };
        let super::super::What::Cooccur(args) = query.what else {
            panic!("parsed the wrong query command")
        };
        assert_eq!(args.unit, EvidenceUnitArg::MoleculeRecord);
        assert_eq!(args.placements, PlacementScopeArg::Unique);
    }

    #[test]
    fn unique_placement_scope_excludes_multimapper_only_records() {
        let predicate = parse_predicates(
            &["r=region:chr1:100-101:+".to_owned()],
            &["chr1".to_owned()],
        )
        .unwrap()
        .remove(0);
        let molecule = MolRec {
            cell: 0,
            umi_class: 7,
            chrom: 0,
            strand_rev: false,
            chains: Default::default(),
            mms: smallvec::smallvec![(100, 0, 0, 1), (1_000_000, 0, 1, 1)],
        };
        let shapes = [Shape {
            blocks: vec![(0, 1)],
        }];
        let unique = MatchContext {
            predicates: std::slice::from_ref(&predicate),
            region_match: RegionMatchArg::Anchor,
            placements: PlacementScopeArg::Unique,
            shapes: &shapes,
            patterns: None,
        };
        assert!(!unique.matches(&predicate, &molecule, &[]).unwrap());

        let direct = MatchContext {
            placements: PlacementScopeArg::Direct,
            ..unique
        };
        assert!(direct.matches(&predicate, &molecule, &[]).unwrap());
    }

    #[test]
    fn pattern_row_limit_is_enforced_during_aggregation() {
        let units = vec![
            EvaluatedUnit {
                cell: 0,
                class: 0,
                chunk: Some(0),
                local_record: Some(0),
                global_record: Some(0),
                contributing_records: 1,
                mask: 1,
                completeness_mask: 3,
                selection: TruthValue::True,
            },
            EvaluatedUnit {
                cell: 1,
                class: 1,
                chunk: Some(0),
                local_record: Some(1),
                global_record: Some(1),
                contributing_records: 1,
                mask: 2,
                completeness_mask: 3,
                selection: TruthValue::False,
            },
        ];
        let scope = super::super::QueryScope {
            selected: None,
            group_names: Vec::new(),
            group_of: FxHashMap::default(),
            selected_per_group: Vec::new(),
            selected_cells: 2,
            aggregation: QueryAggregation::Bulk,
            source: "all_archive_cells",
            source_path: None,
            source_content_blake3: None,
            resolved_mapping_blake3: None,
            archive_cells: 2,
            active: false,
        };
        let error = aggregate_patterns(&units, &scope, &[0, 1], 1).unwrap_err();
        assert!(error.to_string().contains("--max-pattern-rows 1"));
    }

    #[test]
    fn distinct_locus_records_do_not_cooccur_unless_class_union_is_explicit() {
        let expression = Expression::And(
            Box::new(Expression::Predicate(0)),
            Box::new(Expression::Predicate(1)),
        );
        let first_locus = RecordHit {
            chunk: 0,
            local_record: 0,
            global_record: 0,
            class: 7,
            mask: 0b01,
            completeness_mask: 0b11,
        };
        let second_locus = RecordHit {
            chunk: 1,
            local_record: 0,
            global_record: 1,
            class: 7,
            mask: 0b10,
            completeness_mask: 0b11,
        };
        assert_eq!(
            expression.evaluate(first_locus.mask, first_locus.completeness_mask),
            TruthValue::False
        );
        let merged = merge_class_evidence(
            std::slice::from_ref(&first_locus),
            vec![first_locus.clone(), second_locus],
            0b01,
        );
        let evidence = merged.get(&7).unwrap();
        assert_eq!(evidence.mask, 0b11);
        assert_eq!(evidence.records, 2);
        assert_eq!(
            expression.evaluate(evidence.mask, evidence.completeness_mask),
            TruthValue::True
        );
        assert!(UMI_CLASS_IDENTITY.contains("collisions"));
        assert!(UMI_CLASS_IDENTITY.contains("explicit"));
    }
}
