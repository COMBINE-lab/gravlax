//! Whole-collection reverse search over rooted molecular evidence.
//!
//! Planning starts from the collection's compact junction catalogue. Candidate predicates are
//! applied before any source molecule chunk is decoded; surviving coordinates are then routed
//! back to the authoritative molecule chunks and reduced exactly. The collection never becomes a
//! count store.

use super::*;
use crate::archivecmd::TerminalTailRecord;
use crate::rows::MolRec;
use clap::ValueEnum;
use smallvec::SmallVec;
use std::collections::BTreeSet;

const RESULT_SCHEMA: &str = "gravlax.collection.find-events.result.v1";
const ENTITY_SCHEMA: &str = "gravlax.collection.find-events.entities.v1";
const COMPONENT_SCHEMA: &str = "gravlax.collection.find-events.components.v1";
const COUNT_SCHEMA: &str = "gravlax.collection.find-events.counts.v1";
const CAPABILITY_SCHEMA: &str = "gravlax.collection.find-events.capabilities.v1";
const TERMINAL_ANCHOR_SCHEMA: &str = "gravlax.collection.find-events.terminal-anchors.v1";
const TERMINAL_COUNT_SCHEMA: &str = "gravlax.collection.find-events.terminal-counts.v1";

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum SearchKind {
    Junction,
    AltAcceptor,
    AltDonor,
    Cassette,
    TerminalTail,
}

impl SearchKind {
    fn name(self) -> &'static str {
        match self {
            Self::Junction => "junction",
            Self::AltAcceptor => "alt_acceptor",
            Self::AltDonor => "alt_donor",
            Self::Cassette => "cassette",
            Self::TerminalTail => "terminal_tail",
        }
    }

    fn is_junction(self) -> bool {
        self == Self::Junction
    }

    fn is_alternative_site(self) -> bool {
        matches!(self, Self::AltAcceptor | Self::AltDonor)
    }
}

#[derive(clap::Args)]
pub(super) struct Args {
    /// Authenticated collection whose indexes supply reverse-search candidates and routes.
    pub(super) collection: PathBuf,
    /// Entity kind to search; repeat to select several. The default selects every supported kind.
    #[arg(long = "kind", value_enum)]
    kinds: Vec<SearchKind>,
    /// Strict sample-to-donor TSV with header `sample<TAB>donor`. Without it, each sample is its
    /// own biological donor; this explicit default is recorded in result provenance.
    #[arg(long)]
    design: Option<PathBuf>,
    /// Strict cell-group TSV with header `sample<TAB>barcode<TAB>group`. When supplied, only listed
    /// cells are in scope. Without it, every archive cell contributes to the `bulk` group.
    #[arg(long)]
    groups: Option<PathBuf>,
    /// Require exact support in each named group; repeat for multiple groups.
    #[arg(long = "require-group", requires = "groups")]
    required_groups: Vec<String>,
    /// Minimum exact archive UMI-class count in each --require-group.
    #[arg(
        long = "min-group-umi-classes",
        visible_alias = "min-group-umis",
        default_value_t = 1,
        requires = "required_groups"
    )]
    min_group_umis: usize,
    /// Minimum distinct donors with exact unique-chain UMI-class support.
    #[arg(long, default_value_t = 1)]
    min_donors: usize,
    /// Minimum distinct source archives with exact unique-chain UMI-class support.
    #[arg(long, default_value_t = 1)]
    min_samples: usize,
    /// Minimum exact archive UMI classes pooled across the selected cells.
    #[arg(
        long = "min-umi-classes",
        visible_alias = "min-umis",
        default_value_t = 1
    )]
    min_umis: usize,
    /// For alternative-splicing entities, require at least this many exact include-only and
    /// exclude-only archive UMI classes. This prevents one-sided candidates from surviving.
    #[arg(
        long = "min-side-umi-classes",
        visible_alias = "min-side-umis",
        default_value_t = 1
    )]
    min_side_umis: usize,
    /// Minimum catalogue route/support upper bound for every splice component.
    #[arg(long, default_value_t = 2)]
    min_support: u64,
    /// Maximum gap between consecutive exact terminal-tail anchors in one strand-aware cluster.
    /// Zero reports one entity per exact anchor; exact anchor rows are always retained in output.
    #[arg(long, default_value_t = 25)]
    terminal_cluster_bp: u32,
    /// Fail before decoding tail-bearing molecule chunks when capable archives declare more
    /// terminal events than this bound. Results are never truncated.
    #[arg(long, default_value_t = 10_000_000)]
    max_terminal_events: u64,
    /// Uncompressed GTF or compiled AIC used to classify annotation gaps.
    #[arg(long, requires_all = ["assembly", "annotation_label"])]
    annotation: Option<PathBuf>,
    /// Caller-declared reference assembly for --annotation; recorded but not inferred as
    /// compatible with the collection's stamped genome.
    #[arg(long, requires = "annotation")]
    assembly: Option<String>,
    /// Immutable release or descriptive label for --annotation.
    #[arg(long = "annotation-label", requires = "annotation")]
    annotation_label: Option<String>,
    /// Optional expected annotation snapshot digest in blake3:<64 lowercase hex> form.
    #[arg(long = "annotation-digest", requires = "annotation")]
    annotation_digest: Option<String>,
    /// Keep only entities incompatible with every eligible transcript; requires --annotation.
    #[arg(long, requires = "annotation")]
    novel_only: bool,
    /// STARsolo alignment/transcript strand relationship used for annotation compatibility.
    #[arg(long, value_enum, default_value_t = crate::archivecmd::SoloStrandArg::Forward)]
    solo_strand: crate::archivecmd::SoloStrandArg,
    /// Hard limit on retained routed splice candidates plus terminal clusters. Never truncates.
    #[arg(long, default_value_t = 100_000)]
    max_candidates: usize,
    /// Hard limit on attempted splice-candidate definitions before recurrence filtering.
    #[arg(long, default_value_t = 1_000_000)]
    max_candidates_considered: usize,
    /// Hard limit on materialized candidate-to-archive target associations plus routed chunk
    /// postings. This bounds exact-routing memory; results are never truncated.
    #[arg(long, default_value_t = 10_000_000)]
    max_routed_entries: usize,
    /// Hard limit on exact candidate-target checks against molecule junctions. Checked before
    /// each target list is expanded; results are never truncated.
    #[arg(long, default_value_t = 25_000_000)]
    max_exact_match_attempts: u64,
    /// Hard limit on indexed annotation transcript-comparison work. Checked before each local
    /// candidate set is classified; results are never truncated.
    #[arg(long, default_value_t = 10_000_000)]
    max_annotation_comparisons: u64,
    /// Re-hash every source archive before exact routing.
    #[arg(long)]
    verify_content: bool,
    #[command(flatten)]
    pub(super) uniform_output: CollectionIoArgs,
}

#[derive(Clone, Debug)]
struct Design {
    donor_of_sample: Vec<usize>,
    donor_names: Vec<String>,
    source: Option<PathBuf>,
    content_blake3: Option<String>,
}

#[derive(Clone, Debug)]
struct Groups {
    names: Vec<String>,
    /// Packed barcode -> group, one map per collection archive.
    by_sample: Vec<FxHashMap<u32, usize>>,
    source: Option<PathBuf>,
    content_blake3: Option<String>,
    explicit: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Coordinate {
    chrom: u32,
    donor: u32,
    acceptor: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum ComponentSide {
    Support,
    Include,
    Exclude,
}

impl ComponentSide {
    fn name(self) -> &'static str {
        match self {
            Self::Support => "support",
            Self::Include => "include",
            Self::Exclude => "exclude",
        }
    }

    fn mask(self) -> u8 {
        match self {
            Self::Support | Self::Include => 1,
            Self::Exclude => 2,
        }
    }

    fn output_name(self, kind: SearchKind) -> &'static str {
        if kind.is_alternative_site() {
            match self {
                Self::Include => "side_a",
                Self::Exclude => "side_b",
                Self::Support => "support",
            }
        } else {
            self.name()
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Component {
    coordinate: Coordinate,
    side: ComponentSide,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct EntityKey {
    kind: SearchKind,
    chrom: u32,
    /// Strand stored on source alignments and used while exact-routing molecules. `None` is the
    /// explicit unstranded wildcard and coalesces both alignment orientations into one entity.
    strand_rev: Option<bool>,
    /// Transcript/evidence strand after applying STARsolo semantics; `None` is reported as `.`.
    reported_strand_rev: Option<bool>,
    /// At most three: a junction has one, an alternative site two, a cassette three. Keeping them
    /// inline avoids one heap allocation per attempted definition, of which a cohort-wide search
    /// makes millions.
    components: SmallVec<[Component; 3]>,
}

impl EntityKey {
    fn id(&self, chroms: &[String]) -> String {
        let coordinates = self
            .components
            .iter()
            .map(|component| {
                format!(
                    "{}-{}",
                    component.coordinate.donor, component.coordinate.acceptor
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "{}:{}:{}:{coordinates}",
            self.kind.name(),
            chroms[self.chrom as usize],
            match self.reported_strand_rev {
                Some(true) => '-',
                Some(false) => '+',
                None => '.',
            },
        )
    }
}

#[derive(Clone, Debug)]
struct Candidate {
    key: EntityKey,
    catalogue_support_upper_bound: u64,
    catalogue_samples: usize,
    catalogue_donors: usize,
}

/// One candidate component's exact-routing target. Targets live in a single shared arena, so a
/// coordinate's target list is never cloned per archive or per chunk.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Target {
    entity: u32,
    component_bit: u8,
    side: u8,
    /// Requested alignment strand: 0 accepts either orientation, 1 requires forward, 2 requires
    /// reverse. The encoding keeps the sort order of the `Option<bool>` it replaces.
    strand: u8,
}

impl Target {
    fn strand_code(strand_rev: Option<bool>) -> u8 {
        match strand_rev {
            None => 0,
            Some(false) => 1,
            Some(true) => 2,
        }
    }

    fn accepts(self, strand_rev: bool) -> bool {
        self.strand == 0 || (self.strand == 2) == strand_rev
    }
}

/// Half-open range of one coordinate's targets inside the routing arena.
type TargetRange = (u32, u32);

#[derive(Clone, Debug)]
struct RoutedTarget {
    coordinate: Coordinate,
    targets: TargetRange,
    posts: SmallVec<[u32; 4]>,
}

struct RoutedPlan {
    /// Shared target arena; every `RoutedTarget::targets` names a range of it.
    targets: Vec<Target>,
    per_archive: Vec<Vec<RoutedTarget>>,
    target_associations: usize,
    chunk_postings: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct MaskCounts {
    support: usize,
    include_only: usize,
    exclude_only: usize,
    both: usize,
}

impl MaskCounts {
    fn add_mask(&mut self, kind: SearchKind, mask: u8) {
        if kind.is_junction() {
            if mask & 1 != 0 {
                self.support += 1;
            }
        } else {
            match mask & 3 {
                1 => self.include_only += 1,
                2 => self.exclude_only += 1,
                3 => self.both += 1,
                _ => {}
            }
        }
    }

    fn metric(self, kind: SearchKind) -> usize {
        if kind.is_junction() {
            self.support
        } else {
            self.include_only + self.exclude_only
        }
    }

    fn total(self, kind: SearchKind) -> usize {
        if kind.is_junction() {
            self.support
        } else {
            self.include_only + self.exclude_only + self.both
        }
    }

    fn add(&mut self, other: Self) {
        self.support += other.support;
        self.include_only += other.include_only;
        self.exclude_only += other.exclude_only;
        self.both += other.both;
    }
}

fn mask_is_informative(kind: SearchKind, mask: u8) -> bool {
    kind.is_junction() || matches!(mask & 3, 1 | 2)
}

/// One archive's exact counts for a candidate entity in one cell group. The archive and its donor
/// come from the owning `ArchiveExact`, and distinct cells are already reduced to a count, so a
/// cohort run holds tens of millions of these in tens of bytes each.
#[derive(Clone, Copy, Debug)]
struct SampleCount {
    entity: u32,
    group: u32,
    counts: SampleMaskCounts,
    cells: u32,
}

/// `MaskCounts` narrowed for storage. Each field counts distinct archive UMI classes within one
/// archive, which is itself bounded by a `u32` class dictionary, so nothing can overflow.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct SampleMaskCounts {
    support: u32,
    include_only: u32,
    exclude_only: u32,
    both: u32,
}

impl SampleMaskCounts {
    fn narrow(counts: MaskCounts) -> Result<Self> {
        Ok(Self {
            support: u32::try_from(counts.support).context("exact class count exceeds u32")?,
            include_only: u32::try_from(counts.include_only)
                .context("exact class count exceeds u32")?,
            exclude_only: u32::try_from(counts.exclude_only)
                .context("exact class count exceeds u32")?,
            both: u32::try_from(counts.both).context("exact class count exceeds u32")?,
        })
    }

    fn metric(self, kind: SearchKind) -> u32 {
        if kind.is_junction() {
            self.support
        } else {
            self.include_only + self.exclude_only
        }
    }
}

/// Cohort-wide exact evidence for one candidate entity. Recurrence is counted rather than
/// collected: archives contribute once each, and reducing them in donor order makes the distinct
/// donor count a counter plus the previous donor. Per-group counts live in `ExactAggregate`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct EntityExact {
    counts: MaskCounts,
    component_umi_classes: [u32; 3],
    strand_umis: [u32; 2],
    cells: u32,
    samples: u32,
    donors: u32,
    last_donor_plus_one: u32,
}

#[derive(Clone, Debug, Default, Serialize)]
struct GapClassification {
    incompatible_with_every_transcript: bool,
    compatible_transcripts: usize,
    missing_junction: bool,
    boundary: bool,
    strand: bool,
    overlap: bool,
    primary_class: Option<&'static str>,
    overlapping_gene_ids: Vec<String>,
}

fn valid_identifier(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"_.-".contains(&byte))
}

fn read_bound_text(path: &Path, label: &str) -> Result<(String, String)> {
    let bytes =
        std::fs::read(path).with_context(|| format!("reading {label} {}", path.display()))?;
    let digest = format!("blake3:{}", blake3::hash(&bytes).to_hex());
    let text = String::from_utf8(bytes)
        .with_context(|| format!("{label} {} is not UTF-8", path.display()))?;
    Ok((text, digest))
}

fn load_design(collection: &Collection, path: Option<&Path>) -> Result<Design> {
    let Some(path) = path else {
        return Ok(Design {
            donor_of_sample: (0..collection.archives.len()).collect(),
            donor_names: collection
                .archives
                .iter()
                .map(|row| row.id.clone())
                .collect(),
            source: None,
            content_blake3: None,
        });
    };
    let (text, digest) = read_bound_text(path, "donor design")?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .map(|line| line.trim_end_matches('\r'))
        .context("donor design is empty")?;
    if header != "sample\tdonor" {
        bail!("donor design header must be exactly: sample<TAB>donor");
    }
    let sample_index: BTreeMap<&str, usize> = collection
        .archives
        .iter()
        .enumerate()
        .map(|(index, archive)| (archive.id.as_str(), index))
        .collect();
    let mut donor_names = Vec::new();
    let mut donor_of_sample = vec![None; collection.archives.len()];
    for (line_index, raw) in lines.enumerate() {
        let line_no = line_index + 2;
        let line = raw.trim_end_matches('\r');
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 || fields.iter().any(|field| !valid_identifier(field)) {
            bail!("donor design line {line_no} must contain valid sample<TAB>donor identifiers");
        }
        let &sample = sample_index.get(fields[0]).with_context(|| {
            format!(
                "donor design line {line_no} names unknown sample {}",
                fields[0]
            )
        })?;
        if donor_of_sample[sample].is_some() {
            bail!("donor design repeats sample {}", fields[0]);
        }
        let donor = donor_names
            .iter()
            .position(|name| name == fields[1])
            .unwrap_or_else(|| {
                donor_names.push(fields[1].to_owned());
                donor_names.len() - 1
            });
        donor_of_sample[sample] = Some(donor);
    }
    let missing: Vec<&str> = donor_of_sample
        .iter()
        .enumerate()
        .filter_map(|(index, donor)| {
            donor
                .is_none()
                .then_some(collection.archives[index].id.as_str())
        })
        .collect();
    if !missing.is_empty() {
        bail!(
            "donor design omits collection sample(s): {}",
            missing.join(", ")
        );
    }
    Ok(Design {
        donor_of_sample: donor_of_sample.into_iter().map(Option::unwrap).collect(),
        donor_names,
        source: Some(path.to_path_buf()),
        content_blake3: Some(digest),
    })
}

fn load_groups(collection: &Collection, path: Option<&Path>) -> Result<Groups> {
    let Some(path) = path else {
        return Ok(Groups {
            names: vec!["bulk".to_owned()],
            by_sample: (0..collection.archives.len())
                .map(|_| FxHashMap::default())
                .collect(),
            source: None,
            content_blake3: None,
            explicit: false,
        });
    };
    let (text, digest) = read_bound_text(path, "cell-group map")?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .map(|line| line.trim_end_matches('\r'))
        .context("cell-group map is empty")?;
    if header != "sample\tbarcode\tgroup" {
        bail!("cell-group map header must be exactly: sample<TAB>barcode<TAB>group");
    }
    let sample_index: BTreeMap<&str, usize> = collection
        .archives
        .iter()
        .enumerate()
        .map(|(index, archive)| (archive.id.as_str(), index))
        .collect();
    let mut names = Vec::new();
    let mut by_sample: Vec<FxHashMap<u32, usize>> = (0..collection.archives.len())
        .map(|_| FxHashMap::default())
        .collect();
    for (line_index, raw) in lines.enumerate() {
        let line_no = line_index + 2;
        let line = raw.trim_end_matches('\r');
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 3 || !valid_identifier(fields[0]) || !valid_identifier(fields[2]) {
            bail!("cell-group map line {line_no} must contain sample<TAB>barcode<TAB>group");
        }
        let &sample = sample_index.get(fields[0]).with_context(|| {
            format!(
                "cell-group map line {line_no} names unknown sample {}",
                fields[0]
            )
        })?;
        let packed = crate::querycmd::pack_cell_barcode_16(fields[1])
            .with_context(|| format!("cell-group map line {line_no} has an invalid barcode"))?;
        let group = names
            .iter()
            .position(|name| name == fields[2])
            .unwrap_or_else(|| {
                names.push(fields[2].to_owned());
                names.len() - 1
            });
        if by_sample[sample].insert(packed, group).is_some() {
            bail!(
                "cell-group map repeats sample {} barcode {}",
                fields[0],
                fields[1]
            );
        }
    }
    if names.is_empty() {
        bail!("cell-group map contains no assignments");
    }
    Ok(Groups {
        names,
        by_sample,
        source: Some(path.to_path_buf()),
        content_blake3: Some(digest),
        explicit: true,
    })
}

fn selected_kinds(requested: &[SearchKind]) -> BTreeSet<SearchKind> {
    if requested.is_empty() {
        [
            SearchKind::Junction,
            SearchKind::AltAcceptor,
            SearchKind::AltDonor,
            SearchKind::Cassette,
            SearchKind::TerminalTail,
        ]
        .into_iter()
        .collect()
    } else {
        requested.iter().copied().collect()
    }
}

fn validate_kind_strand(
    kinds: &BTreeSet<SearchKind>,
    solo_strand: crate::archivecmd::SoloStrandArg,
) -> Result<()> {
    if kinds.contains(&SearchKind::TerminalTail)
        && !matches!(solo_strand, crate::archivecmd::SoloStrandArg::Forward)
    {
        bail!(
            "--kind terminal-tail requires --solo-strand forward because the archived extraction capability has fixed forward-cDNA strand semantics"
        );
    }
    if matches!(solo_strand, crate::archivecmd::SoloStrandArg::Unstranded)
        && (kinds.contains(&SearchKind::AltAcceptor) || kinds.contains(&SearchKind::AltDonor))
    {
        bail!(
            "--solo-strand unstranded cannot orient --kind alt-acceptor or alt-donor; request junction/cassette or supply the library's STARsolo strand relationship"
        );
    }
    Ok(())
}

fn scan_chain_junctions(chain: &CollectionChain) -> Result<Vec<GlobalJunction>> {
    let mut merged: BTreeMap<Coordinate, GlobalJunction> = BTreeMap::new();
    for layer in &chain.layers {
        for name in layer.file.names().filter(|name| name.starts_with("j.")) {
            let rows = decode_junction_rows(
                &layer.file.read(name)?,
                &layer.manifest.chroms,
                &layer.manifest.archives,
            )?;
            for mut row in rows {
                for route in &mut row.routes {
                    route.archive = u32::try_from(layer.local_to_global[route.archive as usize])
                        .context("global archive index exceeds u32")?;
                }
                let coordinate = Coordinate {
                    chrom: row.chrom,
                    donor: row.donor,
                    acceptor: row.acceptor,
                };
                match merged.entry(coordinate) {
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        entry.insert(row);
                    }
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        let current = entry.get_mut();
                        current.support_upper_bound = current
                            .support_upper_bound
                            .checked_add(row.support_upper_bound)
                            .context("junction support upper bound overflow across layers")?;
                        current.routes.extend(row.routes);
                    }
                }
            }
        }
    }
    let mut rows = Vec::with_capacity(merged.len());
    for (_, mut row) in merged {
        row.routes.sort_unstable_by_key(|route| route.archive);
        if row
            .routes
            .windows(2)
            .any(|pair| pair[0].archive == pair[1].archive)
        {
            bail!("collection chain contains duplicate archive routes for one junction");
        }
        rows.push(row);
    }
    Ok(rows)
}

fn coordinate_of_row(row: &GlobalJunction) -> Coordinate {
    Coordinate {
        chrom: row.chrom,
        donor: row.donor,
        acceptor: row.acceptor,
    }
}

/// The catalogue merges collection layers through a coordinate-keyed map, so its rows are unique
/// and ascending. Planning relies on that to look a coordinate up by binary search and to keep
/// each archive's routes an append-only sorted vector, so the invariant is checked once.
fn require_sorted_catalogue(rows: &[GlobalJunction]) -> Result<()> {
    if rows
        .windows(2)
        .any(|pair| coordinate_of_row(&pair[1]) <= coordinate_of_row(&pair[0]))
    {
        bail!("junction catalogue rows must be unique and ascending by coordinate");
    }
    Ok(())
}

fn catalogue_row(rows: &[GlobalJunction], coordinate: Coordinate) -> Option<&GlobalJunction> {
    rows.binary_search_by(|row| coordinate_of_row(row).cmp(&coordinate))
        .ok()
        .map(|index| &rows[index])
}

/// Catalogue support bound and route recurrence for one candidate definition. `samples` and
/// `donors` are caller-owned scratch buffers: a cohort-wide search evaluates millions of
/// definitions, and a fresh set per definition dominated both time and memory.
fn candidate_catalogue_stats(
    components: &[Component],
    rows: &[GlobalJunction],
    design: &Design,
    samples: &mut Vec<u32>,
    donors: &mut Vec<u32>,
) -> Result<(u64, usize, usize)> {
    if components.is_empty() {
        bail!("candidate has no components");
    }
    samples.clear();
    donors.clear();
    let mut minimum_support = u64::MAX;
    for component in components {
        let row = catalogue_row(rows, component.coordinate)
            .context("candidate component is absent from the junction catalogue")?;
        minimum_support = minimum_support.min(row.support_upper_bound);
        samples.extend(row.routes.iter().map(|route| route.archive));
    }
    // Exact informative support is an OR over a candidate's component junctions (the molecule
    // mask later determines its side). Their route union is therefore the safe recurrence upper
    // bound, including when alternative evidence is distributed across samples or donors.
    samples.sort_unstable();
    samples.dedup();
    for &sample in samples.iter() {
        let donor = design
            .donor_of_sample
            .get(sample as usize)
            .context("junction route references missing collection archive")?;
        donors.push(u32::try_from(*donor).context("biological donor index exceeds u32")?);
    }
    donors.sort_unstable();
    donors.dedup();
    Ok((minimum_support, samples.len(), donors.len()))
}

/// Collects every attempted splice-event definition. Deduplication is a sort over the collected
/// keys rather than an ordered set of heap-allocated keys, which halved the stage's peak memory:
/// a cohort-wide cassette search attempts millions of definitions.
struct CandidateBuilder {
    max_candidates_considered: usize,
    attempted: usize,
    keys: Vec<EntityKey>,
}

#[derive(Clone, Copy)]
struct CandidateThresholds {
    min_support: u64,
    min_samples: usize,
    min_donors: usize,
    max_candidates: usize,
    max_candidates_considered: usize,
}

struct CandidateDiscovery {
    candidates: Vec<Candidate>,
    attempted: usize,
    distinct: usize,
}

impl CandidateBuilder {
    fn push(&mut self, key: EntityKey) -> Result<()> {
        if self.attempted == self.max_candidates_considered {
            bail!(
                "reverse search attempted more than --max-candidates-considered {} splice-event definitions; narrow --kind, increase --min-support, or raise the explicit limit",
                self.max_candidates_considered
            );
        }
        self.attempted += 1;
        self.keys.push(key);
        Ok(())
    }
}

fn component(coordinate: Coordinate, side: ComponentSide) -> Component {
    Component { coordinate, side }
}

fn reported_strand(
    solo_strand: crate::archivecmd::SoloStrandArg,
    alignment_strand_rev: bool,
) -> bool {
    match solo_strand {
        crate::archivecmd::SoloStrandArg::Forward
        | crate::archivecmd::SoloStrandArg::Unstranded => alignment_strand_rev,
        crate::archivecmd::SoloStrandArg::Reverse => !alignment_strand_rev,
    }
}

fn candidate_strands(
    solo_strand: crate::archivecmd::SoloStrandArg,
) -> Vec<(Option<bool>, Option<bool>)> {
    match solo_strand {
        crate::archivecmd::SoloStrandArg::Unstranded => vec![(None, None)],
        crate::archivecmd::SoloStrandArg::Forward => {
            vec![(Some(false), Some(false)), (Some(true), Some(true))]
        }
        crate::archivecmd::SoloStrandArg::Reverse => {
            vec![(Some(false), Some(true)), (Some(true), Some(false))]
        }
    }
}

fn discover_candidates(
    rows: &[GlobalJunction],
    kinds: &BTreeSet<SearchKind>,
    design: &Design,
    solo_strand: crate::archivecmd::SoloStrandArg,
    thresholds: CandidateThresholds,
) -> Result<CandidateDiscovery> {
    if thresholds.max_candidates == 0 {
        bail!("--max-candidates must be at least 1");
    }
    if thresholds.max_candidates_considered == 0 {
        bail!("--max-candidates-considered must be at least 1");
    }
    require_sorted_catalogue(rows)?;
    let eligible: Vec<Coordinate> = rows
        .iter()
        .filter(|row| row.support_upper_bound >= thresholds.min_support)
        .map(coordinate_of_row)
        .collect();
    let mut builder = CandidateBuilder {
        max_candidates_considered: thresholds.max_candidates_considered,
        attempted: 0,
        keys: Vec::new(),
    };

    if kinds.contains(&SearchKind::Junction) {
        for coordinate in &eligible {
            for (strand_rev, reported_strand_rev) in candidate_strands(solo_strand) {
                builder.push(EntityKey {
                    kind: SearchKind::Junction,
                    chrom: coordinate.chrom,
                    strand_rev,
                    reported_strand_rev,
                    components: smallvec::smallvec![component(*coordinate, ComponentSide::Support)],
                })?;
            }
        }
    }

    let mut by_donor: BTreeMap<(u32, u32), Vec<u32>> = BTreeMap::new();
    let mut by_acceptor: BTreeMap<(u32, u32), Vec<u32>> = BTreeMap::new();
    for coordinate in &eligible {
        by_donor
            .entry((coordinate.chrom, coordinate.donor))
            .or_default()
            .push(coordinate.acceptor);
        by_acceptor
            .entry((coordinate.chrom, coordinate.acceptor))
            .or_default()
            .push(coordinate.donor);
    }
    for acceptors in by_donor.values_mut() {
        acceptors.sort_unstable();
        acceptors.dedup();
    }
    for donors in by_acceptor.values_mut() {
        donors.sort_unstable();
        donors.dedup();
    }

    if kinds.contains(&SearchKind::AltAcceptor) || kinds.contains(&SearchKind::AltDonor) {
        for (&(chrom, donor), acceptors) in &by_donor {
            for left in 0..acceptors.len() {
                for right in left + 1..acceptors.len() {
                    let components: SmallVec<[Component; 3]> = smallvec::smallvec![
                        component(
                            Coordinate {
                                chrom,
                                donor,
                                acceptor: acceptors[left],
                            },
                            ComponentSide::Include,
                        ),
                        component(
                            Coordinate {
                                chrom,
                                donor,
                                acceptor: acceptors[right],
                            },
                            ComponentSide::Exclude,
                        ),
                    ];
                    for alignment_strand_rev in [false, true] {
                        let transcript_rev = reported_strand(solo_strand, alignment_strand_rev);
                        let kind = if transcript_rev {
                            SearchKind::AltDonor
                        } else {
                            SearchKind::AltAcceptor
                        };
                        if kinds.contains(&kind) {
                            builder.push(EntityKey {
                                kind,
                                chrom,
                                strand_rev: Some(alignment_strand_rev),
                                reported_strand_rev: Some(transcript_rev),
                                components: components.clone(),
                            })?;
                        }
                    }
                }
            }
        }
    }

    if kinds.contains(&SearchKind::AltDonor) || kinds.contains(&SearchKind::AltAcceptor) {
        for (&(chrom, acceptor), donors) in &by_acceptor {
            for left in 0..donors.len() {
                for right in left + 1..donors.len() {
                    let components: SmallVec<[Component; 3]> = smallvec::smallvec![
                        component(
                            Coordinate {
                                chrom,
                                donor: donors[left],
                                acceptor,
                            },
                            ComponentSide::Include,
                        ),
                        component(
                            Coordinate {
                                chrom,
                                donor: donors[right],
                                acceptor,
                            },
                            ComponentSide::Exclude,
                        ),
                    ];
                    for alignment_strand_rev in [false, true] {
                        let transcript_rev = reported_strand(solo_strand, alignment_strand_rev);
                        let kind = if transcript_rev {
                            SearchKind::AltAcceptor
                        } else {
                            SearchKind::AltDonor
                        };
                        if kinds.contains(&kind) {
                            builder.push(EntityKey {
                                kind,
                                chrom,
                                strand_rev: Some(alignment_strand_rev),
                                reported_strand_rev: Some(transcript_rev),
                                components: components.clone(),
                            })?;
                        }
                    }
                }
            }
        }
    }

    if kinds.contains(&SearchKind::Cassette) {
        for skip in &eligible {
            let Some(left_acceptors) = by_donor.get(&(skip.chrom, skip.donor)) else {
                continue;
            };
            let Some(right_donors) = by_acceptor.get(&(skip.chrom, skip.acceptor)) else {
                continue;
            };
            let left_begin = left_acceptors.partition_point(|&value| value <= skip.donor);
            let left_end = left_acceptors.partition_point(|&value| value < skip.acceptor);
            for &inner_acceptor in &left_acceptors[left_begin..left_end] {
                let right_begin = right_donors.partition_point(|&value| value <= inner_acceptor);
                let right_end = right_donors.partition_point(|&value| value < skip.acceptor);
                for &inner_donor in &right_donors[right_begin..right_end] {
                    let left = Coordinate {
                        chrom: skip.chrom,
                        donor: skip.donor,
                        acceptor: inner_acceptor,
                    };
                    let right = Coordinate {
                        chrom: skip.chrom,
                        donor: inner_donor,
                        acceptor: skip.acceptor,
                    };
                    // Both flanks come directly from the two sorted catalogue maps, so every
                    // tuple reaching this loop is an actual candidate definition and consumes
                    // the same hard attempted-candidate budget as the quadratic event kinds.
                    for (strand_rev, reported_strand_rev) in candidate_strands(solo_strand) {
                        builder.push(EntityKey {
                            kind: SearchKind::Cassette,
                            chrom: skip.chrom,
                            strand_rev,
                            reported_strand_rev,
                            components: smallvec::smallvec![
                                component(left, ComponentSide::Include),
                                component(right, ComponentSide::Include),
                                component(*skip, ComponentSide::Exclude),
                            ],
                        })?;
                    }
                }
            }
        }
    }

    // Deduplicate by sorting the attempted keys, then apply the catalogue predicates in key
    // order. The retained list is therefore already in its reported order, and `--max-candidates`
    // still fails on exactly the same searches the incremental check failed on.
    let mut keys = builder.keys;
    keys.par_sort_unstable();
    keys.dedup();
    let distinct = keys.len();
    let mut candidates = Vec::new();
    let mut samples = Vec::new();
    let mut donors = Vec::new();
    for key in keys {
        let (support, catalogue_samples, catalogue_donors) =
            candidate_catalogue_stats(&key.components, rows, design, &mut samples, &mut donors)?;
        if support < thresholds.min_support
            || catalogue_samples < thresholds.min_samples
            || catalogue_donors < thresholds.min_donors
        {
            continue;
        }
        if candidates.len() == thresholds.max_candidates {
            bail!(
                "reverse-search catalogue exceeds --max-candidates {}; increase --min-support/--min-samples/--min-donors or raise the explicit limit",
                thresholds.max_candidates
            );
        }
        candidates.push(Candidate {
            key,
            catalogue_support_upper_bound: support,
            catalogue_samples,
            catalogue_donors,
        });
    }
    Ok(CandidateDiscovery {
        candidates,
        attempted: builder.attempted,
        distinct,
    })
}

fn route_candidates(
    collection: &Collection,
    rows: &[GlobalJunction],
    candidates: &[Candidate],
    max_routed_entries: usize,
) -> Result<RoutedPlan> {
    if max_routed_entries == 0 {
        bail!("--max-routed-entries must be at least 1");
    }
    if candidates.len() > MAX_PACKED_ENTITIES {
        bail!(
            "reverse search retained {} candidates, beyond the {MAX_PACKED_ENTITIES} the packed exact reducer addresses; strengthen catalogue predicates or lower --max-candidates",
            candidates.len()
        );
    }
    // Targets are built as one flat `(coordinate, target)` list, sorted and deduplicated once,
    // then flattened into a shared arena. Every later consumer names a half-open range of that
    // arena, so no target list is copied per archive, per chunk, or per junction posting.
    let mut pairs: Vec<(Coordinate, Target)> = Vec::new();
    for (entity, candidate) in candidates.iter().enumerate() {
        for (ordinal, component) in candidate.key.components.iter().enumerate() {
            if ordinal >= EXACT_HIT_COMPONENT_BITS as usize {
                bail!("event has too many components for exact-support tracking");
            }
            pairs.push((
                component.coordinate,
                Target {
                    entity: entity as u32,
                    component_bit: 1u8 << ordinal,
                    side: component.side.mask(),
                    strand: Target::strand_code(candidate.key.strand_rev),
                },
            ));
        }
    }
    pairs.par_sort_unstable();
    pairs.dedup();
    let mut arena: Vec<Target> = Vec::with_capacity(pairs.len());
    let mut coordinates: Vec<Coordinate> = Vec::new();
    let mut ranges: Vec<TargetRange> = Vec::new();
    for (coordinate, target) in pairs {
        let start = u32::try_from(arena.len()).context("routed target arena exceeds u32")?;
        if coordinates.last() != Some(&coordinate) {
            coordinates.push(coordinate);
            ranges.push((start, start));
        }
        arena.push(target);
        let end = u32::try_from(arena.len()).context("routed target arena exceeds u32")?;
        ranges
            .last_mut()
            .context("routed target arena lost its coordinate range")?
            .1 = end;
    }
    require_sorted_catalogue(rows)?;
    let mut per_archive: Vec<Vec<RoutedTarget>> =
        (0..collection.archives.len()).map(|_| Vec::new()).collect();
    let mut target_associations = 0usize;
    let mut chunk_postings = 0usize;
    for row in rows {
        let coordinate = coordinate_of_row(row);
        let Ok(index) = coordinates.binary_search(&coordinate) else {
            continue;
        };
        let range = ranges[index];
        let targets = (range.1 - range.0) as usize;
        for route in &row.routes {
            let added = targets
                .checked_add(route.posts.len())
                .context("routed-entry count overflow")?;
            let routed_entries = target_associations
                .checked_add(chunk_postings)
                .and_then(|value| value.checked_add(added))
                .context("routed-entry count overflow")?;
            if routed_entries > max_routed_entries {
                bail!(
                    "reverse-search exact plan requires more than --max-routed-entries {max_routed_entries} candidate/archive associations and chunk postings; strengthen catalogue predicates or raise the explicit limit"
                );
            }
            target_associations = target_associations
                .checked_add(targets)
                .context("routed target-association count overflow")?;
            chunk_postings = chunk_postings
                .checked_add(route.posts.len())
                .context("routed chunk-posting count overflow")?;
            let archive = per_archive
                .get_mut(route.archive as usize)
                .context("junction route references missing collection archive")?;
            match archive.last_mut() {
                Some(last) if last.coordinate == coordinate => {
                    last.posts.extend_from_slice(&route.posts);
                }
                _ => archive.push(RoutedTarget {
                    coordinate,
                    targets: range,
                    posts: SmallVec::from_slice(&route.posts),
                }),
            }
        }
    }
    for archive in &mut per_archive {
        for target in archive.iter_mut() {
            target.posts.sort_unstable();
            target.posts.dedup();
        }
    }
    Ok(RoutedPlan {
        targets: arena,
        per_archive,
        target_associations,
        chunk_postings,
    })
}

#[derive(Debug)]
struct ArchiveExact {
    sample: usize,
    rows: Vec<SampleCount>,
    unique_chunks: usize,
    independent_chunk_decodes: usize,
    planned_bytes: u64,
    actual_bytes: u64,
}

#[derive(Clone, Copy, Debug, Default)]
struct ArchiveEntityExact {
    entity: u32,
    counts: MaskCounts,
    component_umi_classes: [u32; 3],
    cells: u32,
    strands: [u32; 2],
}

/// Routed `(donor, acceptor)` junctions a chunk must match, each naming a range of the shared
/// target arena rather than its own copy of the target list.
type ChunkTargets = FxHashMap<(u32, u32), TargetRange>;

struct ExactMatchBudget {
    limit: u64,
    attempted: std::sync::atomic::AtomicU64,
}

impl ExactMatchBudget {
    fn new(limit: u64) -> Result<Self> {
        if limit == 0 {
            bail!("--max-exact-match-attempts must be at least 1");
        }
        Ok(Self {
            limit,
            attempted: std::sync::atomic::AtomicU64::new(0),
        })
    }

    fn claim(&self, count: usize) -> Result<()> {
        let count = u64::try_from(count).context("exact match-attempt count exceeds u64")?;
        self.attempted
            .fetch_update(
                std::sync::atomic::Ordering::Relaxed,
                std::sync::atomic::Ordering::Relaxed,
                |attempted| {
                    attempted
                        .checked_add(count)
                        .filter(|&next| next <= self.limit)
                },
            )
            .map(|_| ())
            .map_err(|attempted| {
                anyhow::anyhow!(
                    "reverse-search exact reduction would exceed --max-exact-match-attempts {} after {} candidate-target checks; strengthen predicates or raise the explicit limit",
                    self.limit,
                    attempted,
                )
            })
    }

    fn attempted(&self) -> u64 {
        self.attempted.load(std::sync::atomic::Ordering::Relaxed)
    }
}

struct AnnotationComparisonBudget {
    limit: u64,
    attempted: std::sync::atomic::AtomicU64,
}

impl AnnotationComparisonBudget {
    fn new(limit: u64) -> Result<Self> {
        if limit == 0 {
            bail!("--max-annotation-comparisons must be at least 1");
        }
        Ok(Self {
            limit,
            attempted: std::sync::atomic::AtomicU64::new(0),
        })
    }

    fn claim(&self, count: usize) -> Result<()> {
        let count = u64::try_from(count).context("annotation comparison count exceeds u64")?;
        self.attempted
            .fetch_update(
                std::sync::atomic::Ordering::Relaxed,
                std::sync::atomic::Ordering::Relaxed,
                |attempted| {
                    attempted
                        .checked_add(count)
                        .filter(|&next| next <= self.limit)
                },
            )
            .map(|_| ())
            .map_err(|attempted| {
                anyhow::anyhow!(
                    "annotation-gap classification would exceed --max-annotation-comparisons {} after {} indexed transcript-comparison units; strengthen predicates or raise the explicit limit",
                    self.limit,
                    attempted,
                )
            })
    }

    fn attempted(&self) -> u64 {
        self.attempted.load(std::sync::atomic::Ordering::Relaxed)
    }
}

/// Packed exact hit for one `(entity, archive UMI class)` pair. Layout from the low bit: side
/// mask (2 bits), component mask (3 bits), observed alignment-strand mask (2 bits), UMI class
/// (32 bits), entity ordinal (25 bits). Sorting the `u64` orders by entity then class, so a run of
/// equal high bits reduces with one bitwise OR of the seven payload bits. This replaces a hash map
/// keyed by `(entity, class)`: one cohort archive reduces to tens of millions of hits, and eight
/// bytes each is an order of magnitude below what the map cost.
const EXACT_HIT_PAYLOAD_BITS: u32 = 7;
const EXACT_HIT_PAYLOAD_MASK: u64 = (1 << EXACT_HIT_PAYLOAD_BITS) - 1;
const EXACT_HIT_COMPONENT_SHIFT: u32 = 2;
const EXACT_HIT_COMPONENT_BITS: u32 = 3;
const EXACT_HIT_STRAND_SHIFT: u32 = 5;
const EXACT_HIT_CLASS_SHIFT: u32 = EXACT_HIT_PAYLOAD_BITS;
const EXACT_HIT_ENTITY_SHIFT: u32 = EXACT_HIT_CLASS_SHIFT + 32;
/// Entities addressable by the packed reducer. `--max-candidates` is checked against this before
/// any molecule chunk is decoded, so an oversized search fails while planning.
const MAX_PACKED_ENTITIES: usize = 1 << (64 - EXACT_HIT_ENTITY_SHIFT);

fn pack_exact_hit(entity: u32, class: u32, payload: u8) -> u64 {
    debug_assert!(u64::from(payload) <= EXACT_HIT_PAYLOAD_MASK);
    debug_assert!((entity as usize) < MAX_PACKED_ENTITIES);
    (u64::from(entity) << EXACT_HIT_ENTITY_SHIFT)
        | (u64::from(class) << EXACT_HIT_CLASS_SHIFT)
        | u64::from(payload)
}

/// Entity, archive UMI class, side mask, component mask, observed-strand mask.
fn unpack_exact_hit(hit: u64) -> (u32, u32, u8, u8, u8) {
    let payload = (hit & EXACT_HIT_PAYLOAD_MASK) as u8;
    (
        (hit >> EXACT_HIT_ENTITY_SHIFT) as u32,
        ((hit >> EXACT_HIT_CLASS_SHIFT) & u64::from(u32::MAX)) as u32,
        payload & 3,
        (payload >> EXACT_HIT_COMPONENT_SHIFT) & 7,
        payload >> EXACT_HIT_STRAND_SHIFT,
    )
}

/// OR-reduce runs of equal `(entity, class)` in a sorted packed-hit list, in place.
fn reduce_sorted_exact_hits(hits: &mut Vec<u64>) {
    hits.dedup_by(|later, earlier| {
        if *later >> EXACT_HIT_PAYLOAD_BITS == *earlier >> EXACT_HIT_PAYLOAD_BITS {
            *earlier |= *later & EXACT_HIT_PAYLOAD_MASK;
            true
        } else {
            false
        }
    });
}

/// Split sorted packed hits into roughly equal parts that never cut an entity's run, so the
/// per-entity aggregation below can run in parallel and still emit entities in order.
fn entity_hit_slices(hits: &[u64], parts: usize) -> Vec<&[u64]> {
    if hits.is_empty() {
        return Vec::new();
    }
    let target = hits.len().div_ceil(parts.max(1));
    let mut slices = Vec::new();
    let mut start = 0usize;
    while start < hits.len() {
        let mut end = start.saturating_add(target).min(hits.len());
        while end < hits.len()
            && hits[end] >> EXACT_HIT_ENTITY_SHIFT == hits[end - 1] >> EXACT_HIT_ENTITY_SHIFT
        {
            end += 1;
        }
        slices.push(&hits[start..end]);
        start = end;
    }
    slices
}

/// Everything one decoded chunk's molecules are matched against.
#[derive(Clone, Copy)]
struct ChunkMatch<'a> {
    shapes: &'a [evidence_io::archive::Shape],
    wanted: &'a ChunkTargets,
    arena: &'a [Target],
    budget: &'a ExactMatchBudget,
}

/// Accumulate one placement's `(entity, side | component)` hits into `scratch`.
fn collect_placement_targets(
    position: u32,
    shape_id: u32,
    strand_rev: bool,
    chunk: ChunkMatch<'_>,
    scratch: &mut Vec<(u32, u8)>,
) -> Result<()> {
    let shape = chunk
        .shapes
        .get(shape_id as usize)
        .with_context(|| format!("molecule references missing shape {shape_id}"))?;
    for blocks in shape.blocks.windows(2) {
        let donor = position
            .checked_add(blocks[0].0)
            .and_then(|value| value.checked_add(blocks[0].1))
            .context("junction donor coordinate overflow")?;
        let acceptor = position
            .checked_add(blocks[1].0)
            .context("junction acceptor coordinate overflow")?;
        if let Some(&(start, end)) = chunk.wanted.get(&(donor, acceptor)) {
            let targets = chunk
                .arena
                .get(start as usize..end as usize)
                .context("routed target range is outside the shared target arena")?;
            chunk.budget.claim(targets.len())?;
            for target in targets {
                if target.accepts(strand_rev) {
                    scratch.push((
                        target.entity,
                        target.side | (target.component_bit << EXACT_HIT_COMPONENT_SHIFT),
                    ));
                }
            }
        }
    }
    Ok(())
}

/// Append one molecule's packed per-class hits, one per distinct entity it supports.
fn molecule_packed_hits(
    molecule: &MolRec,
    chunk: ChunkMatch<'_>,
    scratch: &mut Vec<(u32, u8)>,
    hits: &mut Vec<u64>,
) -> Result<()> {
    scratch.clear();
    for chain in &molecule.chains {
        for &(position, shape) in &chain.reps {
            collect_placement_targets(position, shape, molecule.strand_rev, chunk, scratch)?;
        }
    }
    // Coordinate/event search deliberately uses only unique-read chain representatives. An MM
    // tuple's stored anchor is the BAM-designated primary placement, while its complete placement
    // set lives in the pattern dictionary. Counting only that primary would make support depend on
    // an arbitrary representative; expanding every alternative would require a different,
    // explicitly declared ambiguity policy and a complete candidate universe.
    if scratch.is_empty() {
        return Ok(());
    }
    scratch.sort_unstable();
    let mut write = 0usize;
    for index in 0..scratch.len() {
        if write > 0 && scratch[write - 1].0 == scratch[index].0 {
            scratch[write - 1].1 |= scratch[index].1;
        } else {
            scratch[write] = scratch[index];
            write += 1;
        }
    }
    scratch.truncate(write);
    let strand = if molecule.strand_rev { 2u8 } else { 1u8 } << EXACT_HIT_STRAND_SHIFT;
    hits.extend(
        scratch
            .iter()
            .map(|&(entity, payload)| pack_exact_hit(entity, molecule.umi_class, payload | strand)),
    );
    Ok(())
}

/// One archive's exact hits, already OR-reduced per `(entity, class)` and sorted.
/// Entity-range buckets used to merge one archive's chunk hits. Each decoded chunk sorts its own
/// hits once and records where each bucket starts; the merge then gathers one bucket at a time, so
/// the transient copy is a fraction of the archive rather than a second copy of all of it.
const EXACT_MERGE_BUCKETS: usize = 32;

/// One worker's reduced hits for a group of chunks, partitioned into `EXACT_MERGE_BUCKETS`
/// entity ranges.
struct ChunkHits {
    hits: Vec<u64>,
    /// `EXACT_MERGE_BUCKETS + 1` ascending offsets into `hits`.
    bucket_offsets: Vec<u32>,
}

/// Entity bits dropped to name a merge bucket, so that every entity falls in one of the buckets
/// and buckets stay in entity order.
fn merge_bucket_shift(entities: usize) -> u32 {
    let bits = usize::BITS - entities.saturating_sub(1).leading_zeros();
    bits.saturating_sub(EXACT_MERGE_BUCKETS.trailing_zeros())
}

fn bucket_offsets(hits: &[u64], shift: u32) -> Result<Vec<u32>> {
    let mut offsets = Vec::with_capacity(EXACT_MERGE_BUCKETS + 1);
    offsets.push(0u32);
    for bucket in 1..EXACT_MERGE_BUCKETS {
        let limit = (bucket as u64)
            .checked_shl(shift + EXACT_HIT_ENTITY_SHIFT)
            .unwrap_or(u64::MAX);
        let offset = hits.partition_point(|&hit| hit < limit);
        offsets.push(u32::try_from(offset).context("chunk exact hit count exceeds u32")?);
    }
    offsets.push(u32::try_from(hits.len()).context("chunk exact hit count exceeds u32")?);
    Ok(offsets)
}

/// Decode every selected chunk in parallel and reduce it to sorted, bucketed packed hits.
fn archive_packed_hits(
    archive: &mut LazyArchive,
    chunks: &[ChunkInfo],
    selected: &[usize],
    chunk_wanted: &[ChunkTargets],
    arena: &[Target],
    budget: &ExactMatchBudget,
    shift: u32,
) -> Result<Vec<ChunkHits>> {
    let shapes = archive.shapes()?;
    let (reader, tables) = archive.reader_and_tables();
    let reader = &*reader;
    // Chunks are decoded in groups that share one growing hit buffer. An archive's hits live
    // until it has been merged, so accumulating them in a few dozen buffers instead of one per
    // chunk keeps the allocator's growth slack to a fraction of the archive.
    let group = selected
        .len()
        .div_ceil(rayon::current_num_threads().saturating_mul(4).max(1))
        .max(1);
    selected
        .par_chunks(group)
        .map(|group| -> Result<ChunkHits> {
            let mut hits = Vec::new();
            let mut scratch = Vec::new();
            for &chunk_index in group {
                // The compressed and decompressed chunk buffers are released before the molecules
                // are matched, so a worker holds one chunk's worth rather than three.
                let molecules = {
                    let (compressed, raw_len) =
                        reader.read_compressed_at(&format!("c{chunk_index}"))?;
                    let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                    decode_chunk(&raw, &chunks[chunk_index], None, tables)?
                };
                let matcher = ChunkMatch {
                    shapes: &shapes,
                    wanted: &chunk_wanted[chunk_index],
                    arena,
                    budget,
                };
                for molecule in &molecules {
                    molecule_packed_hits(molecule, matcher, &mut scratch, &mut hits)?;
                }
            }
            hits.sort_unstable();
            reduce_sorted_exact_hits(&mut hits);
            // These live until the whole archive is merged, so give back the growth slack.
            hits.shrink_to_fit();
            let bucket_offsets = bucket_offsets(&hits, shift)?;
            Ok(ChunkHits {
                hits,
                bucket_offsets,
            })
        })
        .collect()
}

/// Gather one entity-range bucket from every chunk, then sort and OR-reduce it.
fn merge_bucket(chunk_hits: &[ChunkHits], bucket: usize) -> Result<Vec<u64>> {
    let total = chunk_hits
        .iter()
        .try_fold(0usize, |total, chunk| {
            let span = (chunk.bucket_offsets[bucket + 1] - chunk.bucket_offsets[bucket]) as usize;
            total.checked_add(span)
        })
        .context("exact hit count overflow")?;
    let mut merged = Vec::with_capacity(total);
    for chunk in chunk_hits {
        let start = chunk.bucket_offsets[bucket] as usize;
        let end = chunk.bucket_offsets[bucket + 1] as usize;
        merged.extend_from_slice(&chunk.hits[start..end]);
    }
    merged.par_sort_unstable();
    reduce_sorted_exact_hits(&mut merged);
    Ok(merged)
}

/// One archive's per-entity exact aggregate plus its sparse per-group count rows.
#[derive(Debug, Default)]
struct ArchiveAggregate {
    entities: Vec<ArchiveEntityExact>,
    rows: Vec<SampleCount>,
}

/// Reduce one entity-aligned slice of packed hits into per-entity and per-group aggregates.
/// Distinct cells are counted by sorting `(group, cell)` pairs, which needs no per-entity or
/// per-group hash set: a cell belongs to exactly one group, so the deduplicated pair count is the
/// entity's distinct informative cell count and each group's share is its own.
fn aggregate_exact_slice(
    hits: &[u64],
    archive: &LazyArchive,
    group_of_cell: &[u32],
    group_count: usize,
    candidates: &[Candidate],
) -> Result<ArchiveAggregate> {
    let mut out = ArchiveAggregate::default();
    let mut group_counts = vec![MaskCounts::default(); group_count];
    let mut group_touched = vec![false; group_count];
    let mut touched: Vec<u32> = Vec::new();
    let mut cells: Vec<(u32, u32)> = Vec::new();
    let mut index = 0usize;
    while index < hits.len() {
        let entity = (hits[index] >> EXACT_HIT_ENTITY_SHIFT) as u32;
        let mut end = index + 1;
        while end < hits.len() && (hits[end] >> EXACT_HIT_ENTITY_SHIFT) as u32 == entity {
            end += 1;
        }
        let kind = candidates
            .get(entity as usize)
            .context("packed exact hit references an unknown candidate entity")?
            .key
            .kind;
        let mut exact = ArchiveEntityExact {
            entity,
            ..ArchiveEntityExact::default()
        };
        let mut supported = false;
        cells.clear();
        touched.clear();
        for &hit in &hits[index..end] {
            let (_, class, side, components, strands) = unpack_exact_hit(hit);
            let cell = archive.cell_of_cached(class)?;
            let group = *group_of_cell
                .get(cell as usize)
                .with_context(|| format!("cell id {cell} is outside the archive dictionary"))?;
            if group == u32::MAX {
                continue;
            }
            supported = true;
            exact.counts.add_mask(kind, side);
            for (ordinal, count) in exact.component_umi_classes.iter_mut().enumerate() {
                if components & (1u8 << ordinal) != 0 {
                    *count += 1;
                }
            }
            if !group_touched[group as usize] {
                group_touched[group as usize] = true;
                touched.push(group);
            }
            group_counts[group as usize].add_mask(kind, side);
            if mask_is_informative(kind, side) {
                if strands & 1 != 0 {
                    exact.strands[0] += 1;
                }
                if strands & 2 != 0 {
                    exact.strands[1] += 1;
                }
                cells.push((group, cell));
            }
        }
        if supported {
            cells.sort_unstable();
            cells.dedup();
            exact.cells = u32::try_from(cells.len()).context("exact cell count exceeds u32")?;
            out.entities.push(exact);
            touched.sort_unstable();
            for &group in &touched {
                group_touched[group as usize] = false;
                let counts = std::mem::take(&mut group_counts[group as usize]);
                if counts.total(kind) == 0 {
                    continue;
                }
                let group_cells = cells.iter().filter(|(name, _)| *name == group).count();
                out.rows.push(SampleCount {
                    entity,
                    group,
                    counts: SampleMaskCounts::narrow(counts)?,
                    cells: u32::try_from(group_cells).context("group cell count exceeds u32")?,
                });
            }
        }
        index = end;
    }
    Ok(out)
}

fn exact_archive(
    collection: &Collection,
    sample: usize,
    routes: &[RoutedTarget],
    arena: &[Target],
    groups: &Groups,
    candidates: &[Candidate],
    budget: &ExactMatchBudget,
) -> Result<(ArchiveExact, Vec<ArchiveEntityExact>)> {
    let entry = &collection.archives[sample];
    let chunks = chunk_infos(entry);
    let mut chunk_wanted: Vec<ChunkTargets> =
        (0..chunks.len()).map(|_| FxHashMap::default()).collect();
    let mut independent_chunk_decodes = 0usize;
    for route in routes {
        for &post in &route.posts {
            let chunk = chunks
                .get(post as usize)
                .with_context(|| format!("reverse-search route references missing chunk {post}"))?;
            if chunk.chrom != route.coordinate.chrom {
                bail!("reverse-search route references a chunk on the wrong chromosome");
            }
            // One archive chunk covers one chromosome, so a `(donor, acceptor)` pair names exactly
            // one routed coordinate and can carry that coordinate's arena range directly. Nothing
            // copies, sorts, or deduplicates a target list per chunk any more.
            chunk_wanted[post as usize].insert(
                (route.coordinate.donor, route.coordinate.acceptor),
                route.targets,
            );
            independent_chunk_decodes += 1;
        }
    }
    let selected: Vec<usize> = chunk_wanted
        .iter()
        .enumerate()
        .filter_map(|(index, wanted)| (!wanted.is_empty()).then_some(index))
        .collect();
    let planned_bytes = planned_bytes(
        entry,
        &selected
            .iter()
            .map(|&index| index as u32)
            .collect::<Vec<_>>(),
    )?;
    let mut archive = open_source(collection, sample, None)?;
    let shift = merge_bucket_shift(candidates.len());
    let chunk_hits = archive_packed_hits(
        &mut archive,
        &chunks,
        &selected,
        &chunk_wanted,
        arena,
        budget,
        shift,
    )?;
    drop(chunk_wanted);
    archive.prefetch_coc(
        chunk_hits
            .iter()
            .flat_map(|chunk| chunk.hits.iter().map(|&hit| unpack_exact_hit(hit).1)),
    )?;
    // Resolving each cell's group through a dense per-cell table keeps the aggregation free of
    // barcode hashing and lets out-of-scope cells be skipped with one comparison.
    let cell_dictionary = archive.cells()?;
    let group_of_cell: Vec<u32> = if groups.explicit {
        cell_dictionary
            .iter()
            .map(|packed| {
                groups.by_sample[sample]
                    .get(packed)
                    .map_or(u32::MAX, |&group| group as u32)
            })
            .collect()
    } else {
        vec![0u32; cell_dictionary.len()]
    };
    let group_count = groups.names.len();
    let mut entities = Vec::new();
    let mut rows = Vec::new();
    for bucket in 0..EXACT_MERGE_BUCKETS {
        let merged = merge_bucket(&chunk_hits, bucket)?;
        if merged.is_empty() {
            continue;
        }
        let parts = entity_hit_slices(&merged, rayon::current_num_threads().saturating_mul(4));
        let aggregates: Vec<ArchiveAggregate> = parts
            .par_iter()
            .map(|slice| {
                aggregate_exact_slice(slice, &archive, &group_of_cell, group_count, candidates)
            })
            .collect::<Result<_>>()?;
        for mut part in aggregates {
            entities.append(&mut part.entities);
            rows.append(&mut part.rows);
        }
    }
    drop(chunk_hits);
    release_free_heap();
    let actual_bytes = archive.reader().bytes_read();
    Ok((
        ArchiveExact {
            sample,
            rows,
            unique_chunks: selected.len(),
            independent_chunk_decodes,
            planned_bytes,
            actual_bytes,
        },
        entities,
    ))
}

/// Cohort-wide exact aggregates, indexed by candidate entity.
struct ExactAggregate {
    entities: Vec<EntityExact>,
    /// Informative-metric total per `(entity, --require-group ordinal)`, flattened. A dense row of
    /// small counters replaces a per-entity map from group to counts.
    required_group_metrics: Vec<u32>,
    required_groups: usize,
}

impl ExactAggregate {
    fn new(entities: usize, required_groups: usize) -> Self {
        Self {
            entities: vec![EntityExact::default(); entities],
            required_group_metrics: vec![0u32; entities * required_groups],
            required_groups,
        }
    }

    fn required_metrics(&self, entity: usize) -> &[u32] {
        let start = entity * self.required_groups;
        &self.required_group_metrics[start..start + self.required_groups]
    }
}

/// Fold one archive's aggregate for a candidate into its cohort-wide total. Archives contribute
/// at most once each, so distinct samples is a counter; archives are reduced in donor order, so
/// distinct donors is a counter plus the donor that last contributed.
fn merge_archive_entity(
    slot: &mut EntityExact,
    archive_exact: &ArchiveEntityExact,
    kind: SearchKind,
    donor: u32,
) -> Result<()> {
    for (total, count) in slot
        .component_umi_classes
        .iter_mut()
        .zip(archive_exact.component_umi_classes)
    {
        *total = total
            .checked_add(count)
            .context("exact component UMI-class count overflow")?;
    }
    if archive_exact.counts.metric(kind) == 0 {
        return Ok(());
    }
    slot.counts.add(archive_exact.counts);
    slot.samples += 1;
    if slot.last_donor_plus_one != donor + 1 {
        slot.donors += 1;
        slot.last_donor_plus_one = donor + 1;
    }
    slot.cells += archive_exact.cells;
    slot.strand_umis[0] += archive_exact.strands[0];
    slot.strand_umis[1] += archive_exact.strands[1];
    Ok(())
}

fn reduce_exact(
    collection: &Collection,
    plan: &RoutedPlan,
    groups: &Groups,
    design: &Design,
    candidates: &[Candidate],
    required_groups: &[usize],
    budget: &ExactMatchBudget,
) -> Result<(Vec<ArchiveExact>, ExactAggregate)> {
    let mut group_ordinal = vec![usize::MAX; groups.names.len()];
    for (ordinal, &group) in required_groups.iter().enumerate() {
        group_ordinal[group] = ordinal;
    }
    let mut exact = ExactAggregate::new(candidates.len(), required_groups.len());
    // Archives are reduced one at a time, with their chunks decoded in parallel, so only one
    // archive's packed hits are ever resident. Samples of one donor stay contiguous, which lets
    // distinct-donor recurrence be a counter plus the previous donor instead of a per-entity set.
    // Donor groups then go in decreasing routed-posting order, so the archive with the most exact
    // evidence is reduced while the fewest per-group count rows have accumulated.
    let mut donor_cost = vec![0usize; design.donor_names.len()];
    for (sample, routes) in plan.per_archive.iter().enumerate() {
        let postings: usize = routes.iter().map(|route| route.posts.len()).sum();
        donor_cost[design.donor_of_sample[sample]] += postings;
    }
    let mut order: Vec<usize> = (0..collection.archives.len())
        .filter(|&sample| !plan.per_archive[sample].is_empty())
        .collect();
    order.sort_unstable_by_key(|&sample| {
        let donor = design.donor_of_sample[sample];
        (std::cmp::Reverse(donor_cost[donor]), donor, sample)
    });
    let mut archives = Vec::with_capacity(order.len());
    for sample in order {
        let (archive, entities) = exact_archive(
            collection,
            sample,
            &plan.per_archive[sample],
            &plan.targets,
            groups,
            candidates,
            budget,
        )?;
        let donor = u32::try_from(design.donor_of_sample[sample])
            .context("biological donor index exceeds u32")?;
        for archive_exact in &entities {
            let entity = archive_exact.entity as usize;
            let kind = candidates[entity].key.kind;
            merge_archive_entity(&mut exact.entities[entity], archive_exact, kind, donor)?;
        }
        drop(entities);
        if !required_groups.is_empty() {
            for row in &archive.rows {
                let ordinal = group_ordinal[row.group as usize];
                if ordinal == usize::MAX {
                    continue;
                }
                let kind = candidates[row.entity as usize].key.kind;
                let metric = row.counts.metric(kind);
                let slot = &mut exact.required_group_metrics
                    [row.entity as usize * required_groups.len() + ordinal];
                *slot = slot
                    .checked_add(metric)
                    .context("exact group UMI-class count overflow")?;
            }
        }
        release_free_heap();
        archives.push(archive);
    }
    archives.sort_unstable_by_key(|archive| archive.sample);
    Ok((archives, exact))
}

#[derive(Clone, Copy, Debug)]
struct TailHit {
    sample: usize,
    donor: usize,
    group: usize,
    cell: u32,
    class: u32,
    chrom: u32,
    strand_rev: bool,
    anchor: u32,
    signal: evidence_io::terminal_tail::TerminalTailSignal,
}

#[derive(Debug)]
struct TailArchivePlan {
    sample: usize,
    available: bool,
    declared_selected_molecules: u64,
    declared_events: u64,
    planning_bytes: u64,
}

#[derive(Debug)]
struct TailArchiveScan {
    hits: Vec<TailHit>,
    routed_chunks: usize,
    actual_bytes: u64,
    annotation_excluded_routes: usize,
    annotation_excluded_events: u64,
    annotation_unmatched_chroms: BTreeSet<u32>,
}

#[derive(Clone, Debug)]
struct TailCount {
    sample: usize,
    donor: usize,
    group: usize,
    umis: usize,
    cells: usize,
}

#[derive(Clone, Debug)]
struct TailAnchorCount {
    anchor: u32,
    umis: usize,
    cells: usize,
    samples: usize,
    donors: usize,
    max_clip_len: u8,
    max_tail_bases: u8,
    max_terminal_run: u8,
    counts: Vec<TailCount>,
    gap: Option<GapClassification>,
}

#[derive(Clone, Debug)]
struct TailEntity {
    id: String,
    chrom: u32,
    strand_rev: bool,
    start: u32,
    end: u32,
    summit: u32,
    umis: usize,
    cells: usize,
    samples: usize,
    donors: usize,
    group_umis: FxHashMap<usize, usize>,
    counts: Vec<TailCount>,
    anchors: Vec<TailAnchorCount>,
    gap: Option<GapClassification>,
}

#[derive(Clone, Debug, Default)]
struct TailCapabilitySummary {
    requested: bool,
    archive_available: Vec<bool>,
    available_archives: usize,
    unavailable_archives: usize,
    available_donors: usize,
    declared_selected_molecules: u64,
    declared_events: u64,
    routed_chunks: usize,
    actual_bytes: u64,
    candidate_clusters: usize,
    annotation_excluded_routes: usize,
    annotation_excluded_events: u64,
    annotation_unmatched_chroms: BTreeSet<u32>,
}

struct TailCountSummary {
    counts: Vec<TailCount>,
    umis: usize,
    cells: usize,
    samples: usize,
    donors: usize,
    group_umis: FxHashMap<usize, usize>,
}

type TailGroupKey = (usize, usize, usize);
type TailDistinctEvidence = (FxHashSet<u32>, FxHashSet<u32>);
type TailGroupMap = BTreeMap<TailGroupKey, TailDistinctEvidence>;

fn decode_tail_chunk(
    archive: &mut LazyArchive,
    chunks: &[ChunkInfo],
    route: evidence_io::terminal_tail::TerminalTailRoute,
    molecule_base: u64,
) -> Result<Vec<TerminalTailRecord>> {
    let info = chunks.get(route.chunk as usize).with_context(|| {
        format!(
            "terminal-tail route references missing chunk {}",
            route.chunk
        )
    })?;
    let molecules = {
        let (reader, tables) = archive.reader_and_tables();
        let (compressed, raw_len) = reader.read_compressed_at(&format!("c{}", route.chunk))?;
        let raw = evidence_io::format::decompress(&compressed, raw_len)?;
        crate::archivecmd::decode_chunk_identities(&raw, info, tables)?
    };
    archive.terminal_tail_records_projected(route, info, molecule_base, &molecules)
}

fn plan_tail_archive(collection: &Collection, sample: usize) -> Result<TailArchivePlan> {
    let mut archive = open_source(collection, sample, None)?;
    let Some(metadata) = archive.terminal_tail_capability().cloned() else {
        return Ok(TailArchivePlan {
            sample,
            available: false,
            declared_selected_molecules: 0,
            declared_events: 0,
            planning_bytes: archive.reader().bytes_read(),
        });
    };
    Ok(TailArchivePlan {
        sample,
        available: true,
        declared_selected_molecules: metadata.selected_molecules,
        declared_events: metadata.events,
        planning_bytes: archive.reader().bytes_read(),
    })
}

fn scan_tail_archive(
    collection: &Collection,
    plan: &TailArchivePlan,
    groups: &Groups,
    design: &Design,
    annotation: Option<&AnnotationIndex>,
) -> Result<TailArchiveScan> {
    debug_assert!(plan.available);
    let sample = plan.sample;
    let mut archive = open_source(collection, sample, None)?;
    let routes = archive
        .terminal_tail_routes()?
        .context("terminal-tail capability has no route index")?
        .to_vec();
    let mut retained_routes = Vec::with_capacity(routes.len());
    let mut annotation_excluded_routes = 0usize;
    let mut annotation_excluded_events = 0u64;
    let mut annotation_unmatched_chroms = BTreeSet::new();
    for route in routes {
        if let Some(annotation) = annotation {
            if !annotation.chrom_is_matched(collection, route.chrom, "terminal-tail")? {
                annotation_excluded_routes += 1;
                annotation_excluded_events = annotation_excluded_events
                    .checked_add(u64::from(route.events))
                    .context("annotation-excluded terminal-event count overflow")?;
                annotation_unmatched_chroms.insert(route.chrom);
                continue;
            }
        }
        retained_routes.push(route);
    }
    let routes = retained_routes;
    let chunks = chunk_infos(&collection.archives[plan.sample]);
    let mut molecule_bases = Vec::with_capacity(chunks.len());
    let mut molecule_base = 0u64;
    for chunk in &chunks {
        molecule_bases.push(molecule_base);
        molecule_base = molecule_base
            .checked_add(u64::from(chunk.n_mols))
            .context("terminal-tail molecule ordinal overflow")?;
    }
    let mut records = Vec::new();
    for route in &routes {
        records.extend(decode_tail_chunk(
            &mut archive,
            &chunks,
            *route,
            molecule_bases[route.chunk as usize],
        )?);
    }
    let cell_dictionary = archive.cells()?.to_vec();
    let mut hits = Vec::with_capacity(records.len());
    for record in records {
        let packed = *cell_dictionary.get(record.cell as usize).with_context(|| {
            format!(
                "terminal-tail cell {} is outside the dictionary",
                record.cell
            )
        })?;
        let group = if groups.explicit {
            let Some(&group) = groups.by_sample[sample].get(&packed) else {
                continue;
            };
            group
        } else {
            0
        };
        hits.push(TailHit {
            sample,
            donor: design.donor_of_sample[sample],
            group,
            cell: record.cell,
            class: record.umi_class,
            chrom: record.chrom,
            strand_rev: record.strand_rev,
            anchor: record.anchor,
            signal: record.signal,
        });
    }
    hits.sort_unstable_by_key(|hit| (hit.chrom, hit.strand_rev, hit.anchor, hit.sample, hit.class));
    let mut unique: Vec<TailHit> = Vec::with_capacity(hits.len());
    for hit in hits {
        match unique.last_mut() {
            Some(previous)
                if (
                    previous.chrom,
                    previous.strand_rev,
                    previous.anchor,
                    previous.sample,
                    previous.class,
                ) == (hit.chrom, hit.strand_rev, hit.anchor, hit.sample, hit.class) =>
            {
                if hit.signal.stronger_than(previous.signal) {
                    previous.signal = hit.signal;
                }
            }
            _ => unique.push(hit),
        }
    }
    Ok(TailArchiveScan {
        hits: unique,
        routed_chunks: routes.len(),
        actual_bytes: archive.reader().bytes_read(),
        annotation_excluded_routes,
        annotation_excluded_events,
        annotation_unmatched_chroms,
    })
}

fn tail_counts(hits: &[TailHit]) -> TailCountSummary {
    let mut grouped = TailGroupMap::new();
    for hit in hits {
        let (classes, cells) = grouped
            .entry((hit.sample, hit.donor, hit.group))
            .or_insert_with(|| (FxHashSet::default(), FxHashSet::default()));
        classes.insert(hit.class);
        cells.insert(hit.cell);
    }
    let counts: Vec<TailCount> = grouped
        .into_iter()
        .map(|((sample, donor, group), (classes, cells))| TailCount {
            sample,
            donor,
            group,
            umis: classes.len(),
            cells: cells.len(),
        })
        .collect();
    let umis = counts.iter().map(|count| count.umis).sum();
    let cells = counts.iter().map(|count| count.cells).sum();
    let samples: FxHashSet<usize> = counts.iter().map(|count| count.sample).collect();
    let donors: FxHashSet<usize> = counts.iter().map(|count| count.donor).collect();
    let mut group_umis = FxHashMap::default();
    for count in &counts {
        *group_umis.entry(count.group).or_insert(0) += count.umis;
    }
    TailCountSummary {
        counts,
        umis,
        cells,
        samples: samples.len(),
        donors: donors.len(),
        group_umis,
    }
}

fn build_tail_entity(
    chroms: &[String],
    chrom: u32,
    strand_rev: bool,
    hits: &[TailHit],
    annotation: Option<(&AnnotationIndex, &AnnotationComparisonBudget)>,
) -> Result<TailEntity> {
    let start = hits
        .first()
        .context("terminal-tail cluster is empty")?
        .anchor;
    let last = hits
        .last()
        .context("terminal-tail cluster is empty")?
        .anchor;
    let end = last
        .checked_add(1)
        .context("terminal-tail cluster end coordinate overflow")?;
    let mut by_anchor: BTreeMap<u32, Vec<TailHit>> = BTreeMap::new();
    for hit in hits {
        by_anchor.entry(hit.anchor).or_default().push(*hit);
    }
    let mut anchors = Vec::with_capacity(by_anchor.len());
    for (anchor, anchor_hits) in by_anchor {
        let summary = tail_counts(&anchor_hits);
        anchors.push(TailAnchorCount {
            anchor,
            umis: summary.umis,
            cells: summary.cells,
            samples: summary.samples,
            donors: summary.donors,
            max_clip_len: anchor_hits
                .iter()
                .map(|hit| hit.signal.clip_len)
                .max()
                .unwrap_or(0),
            max_tail_bases: anchor_hits
                .iter()
                .map(|hit| hit.signal.tail_bases)
                .max()
                .unwrap_or(0),
            max_terminal_run: anchor_hits
                .iter()
                .map(|hit| hit.signal.terminal_run)
                .max()
                .unwrap_or(0),
            counts: summary.counts,
            gap: match annotation {
                Some((annotation, budget)) => Some(classify_terminal_anchor(
                    annotation, chrom, strand_rev, anchor, budget,
                )?),
                None => None,
            },
        });
    }
    let summit = anchors
        .iter()
        .max_by_key(|anchor| (anchor.umis, std::cmp::Reverse(anchor.anchor)))
        .context("terminal-tail cluster has no exact anchors")?
        .anchor;
    let summary = tail_counts(hits);
    let group_umis = summary.group_umis;
    let id = format!(
        "terminal_tail:{}:{}:{start}-{end}:summit={summit}",
        chroms[chrom as usize],
        if strand_rev { '-' } else { '+' }
    );
    let mut entity = TailEntity {
        id,
        chrom,
        strand_rev,
        start,
        end,
        summit,
        umis: summary.umis,
        cells: summary.cells,
        samples: summary.samples,
        donors: summary.donors,
        group_umis,
        counts: summary.counts,
        anchors,
        gap: None,
    };
    entity.gap = match annotation {
        Some((annotation, budget)) => Some(classify_terminal_gap(annotation, &entity, budget)?),
        None => None,
    };
    Ok(entity)
}

fn retain_incompatible_tail_anchors(
    hits: &mut Vec<TailHit>,
    annotation: &AnnotationIndex,
    budget: &AnnotationComparisonBudget,
) -> Result<()> {
    let mut compatibility = BTreeMap::new();
    for hit in hits.iter() {
        let key = (hit.chrom, hit.strand_rev, hit.anchor);
        if let std::collections::btree_map::Entry::Vacant(entry) = compatibility.entry(key) {
            entry.insert(
                classify_terminal_anchor(
                    annotation,
                    hit.chrom,
                    hit.strand_rev,
                    hit.anchor,
                    budget,
                )?
                .incompatible_with_every_transcript,
            );
        }
    }
    hits.retain(|hit| compatibility[&(hit.chrom, hit.strand_rev, hit.anchor)]);
    Ok(())
}

fn reserve_terminal_candidate(
    splice_candidate_count: usize,
    terminal_candidate_count: &mut usize,
    max_candidates: usize,
) -> Result<()> {
    let combined_candidates = splice_candidate_count
        .checked_add(*terminal_candidate_count)
        .context("reverse-search candidate count overflow")?;
    if combined_candidates >= max_candidates {
        bail!(
            "reverse-search catalogue contains at least {} splice/junction candidates and {} terminal-tail clusters, exceeding --max-candidates {}; strengthen exact/catalogue thresholds or raise the explicit limit",
            splice_candidate_count,
            terminal_candidate_count
                .checked_add(1)
                .context("terminal-tail candidate count overflow")?,
            max_candidates,
        );
    }
    *terminal_candidate_count = terminal_candidate_count
        .checked_add(1)
        .context("terminal-tail candidate count overflow")?;
    Ok(())
}

fn scan_terminal_tails(
    collection: &Collection,
    groups: &Groups,
    design: &Design,
    splice_candidate_count: usize,
    args: &Args,
    required_groups: &[usize],
    annotation: Option<(&AnnotationIndex, &AnnotationComparisonBudget)>,
) -> Result<(Vec<TailEntity>, TailCapabilitySummary)> {
    let mut plans: Vec<TailArchivePlan> = (0..collection.archives.len())
        .into_par_iter()
        .map(|sample| plan_tail_archive(collection, sample))
        .collect::<Result<_>>()?;
    plans.sort_unstable_by_key(|plan| plan.sample);
    let declared_selected_molecules = plans.iter().try_fold(0u64, |sum, plan| {
        sum.checked_add(plan.declared_selected_molecules)
            .context("terminal-tail declared selected-molecule count overflow")
    })?;
    let declared_events = plans.iter().try_fold(0u64, |sum, plan| {
        sum.checked_add(plan.declared_events)
            .context("terminal-tail declared event count overflow")
    })?;
    let planning_bytes = plans.iter().try_fold(0u64, |sum, plan| {
        sum.checked_add(plan.planning_bytes)
            .context("terminal-tail planning byte count overflow")
    })?;
    let mut capability = TailCapabilitySummary {
        requested: true,
        archive_available: plans.iter().map(|plan| plan.available).collect(),
        available_archives: plans.iter().filter(|plan| plan.available).count(),
        unavailable_archives: plans.iter().filter(|plan| !plan.available).count(),
        declared_selected_molecules,
        declared_events,
        routed_chunks: 0,
        actual_bytes: planning_bytes,
        ..TailCapabilitySummary::default()
    };
    capability.available_donors = plans
        .iter()
        .filter(|plan| plan.available)
        .map(|plan| design.donor_of_sample[plan.sample])
        .collect::<FxHashSet<_>>()
        .len();
    if capability.declared_events > args.max_terminal_events {
        bail!(
            "terminal-tail capability indexes declare {} events, exceeding --max-terminal-events {}; raise the explicit limit to authorize exact decoding",
            capability.declared_events,
            args.max_terminal_events,
        );
    }
    if capability.available_archives < args.min_samples
        || capability.available_donors < args.min_donors
        || capability.declared_selected_molecules < args.min_umis as u64
    {
        return Ok((Vec::new(), capability));
    }
    let scans: Vec<TailArchiveScan> = plans
        .par_iter()
        .filter(|plan| plan.available)
        .map(|plan| scan_tail_archive(collection, plan, groups, design, annotation.map(|v| v.0)))
        .collect::<Result<_>>()?;
    capability.routed_chunks = scans.iter().try_fold(0usize, |sum, scan| {
        sum.checked_add(scan.routed_chunks)
            .context("terminal-tail routed-chunk count overflow")
    })?;
    capability.annotation_excluded_routes = scans.iter().try_fold(0usize, |sum, scan| {
        sum.checked_add(scan.annotation_excluded_routes)
            .context("annotation-excluded terminal-route count overflow")
    })?;
    capability.annotation_excluded_events = scans.iter().try_fold(0u64, |sum, scan| {
        sum.checked_add(scan.annotation_excluded_events)
            .context("annotation-excluded terminal-event count overflow")
    })?;
    for scan in &scans {
        capability
            .annotation_unmatched_chroms
            .extend(scan.annotation_unmatched_chroms.iter().copied());
    }
    capability.actual_bytes = scans
        .iter()
        .try_fold(capability.actual_bytes, |sum, scan| {
            sum.checked_add(scan.actual_bytes)
                .context("terminal-tail archive byte count overflow")
        })?;
    let mut hits: Vec<TailHit> = scans.into_iter().flat_map(|scan| scan.hits).collect();
    hits.sort_unstable_by_key(|hit| (hit.chrom, hit.strand_rev, hit.anchor, hit.sample, hit.class));
    if args.novel_only {
        let (annotation, budget) = annotation.context("--novel-only requires --annotation")?;
        retain_incompatible_tail_anchors(&mut hits, annotation, budget)?;
    }
    let mut entities = Vec::new();
    let mut begin = 0usize;
    while begin < hits.len() {
        let chrom = hits[begin].chrom;
        let strand_rev = hits[begin].strand_rev;
        let mut end = begin + 1;
        let mut previous_anchor = hits[begin].anchor;
        while end < hits.len()
            && hits[end].chrom == chrom
            && hits[end].strand_rev == strand_rev
            && hits[end].anchor <= previous_anchor.saturating_add(args.terminal_cluster_bp)
        {
            previous_anchor = hits[end].anchor;
            end += 1;
        }
        reserve_terminal_candidate(
            splice_candidate_count,
            &mut capability.candidate_clusters,
            args.max_candidates,
        )?;
        let entity = build_tail_entity(
            &collection.chroms,
            chrom,
            strand_rev,
            &hits[begin..end],
            annotation,
        )?;
        if entity.umis >= args.min_umis
            && entity.samples >= args.min_samples
            && entity.donors >= args.min_donors
            && required_groups.iter().all(|group| {
                entity.group_umis.get(group).copied().unwrap_or(0) >= args.min_group_umis
            })
        {
            entities.push(entity);
        }
        begin = end;
    }
    Ok((entities, capability))
}

#[derive(Clone, Debug)]
struct AnnotationTranscript {
    gene_id: String,
    strand_rev: bool,
    start: u32,
    end: u32,
    junctions: FxHashSet<(u32, u32)>,
    donors: FxHashSet<u32>,
    acceptors: FxHashSet<u32>,
}

#[derive(Debug, Default)]
struct AnnotationIntervalIndex {
    transcript_by_start: Vec<usize>,
    prefix_max_end: Vec<u32>,
}

impl AnnotationIntervalIndex {
    fn new(transcripts: &[AnnotationTranscript]) -> Self {
        let mut transcript_by_start: Vec<usize> = (0..transcripts.len()).collect();
        transcript_by_start.sort_unstable_by_key(|&index| {
            (transcripts[index].start, transcripts[index].end, index)
        });
        let mut prefix_max_end = Vec::with_capacity(transcript_by_start.len());
        let mut max_end = 0u32;
        for &index in &transcript_by_start {
            max_end = max_end.max(transcripts[index].end);
            prefix_max_end.push(max_end);
        }
        Self {
            transcript_by_start,
            prefix_max_end,
        }
    }

    fn overlapping_interval(
        &self,
        transcripts: &[AnnotationTranscript],
        start: u32,
        end: u32,
        budget: &AnnotationComparisonBudget,
    ) -> Result<Vec<usize>> {
        let upper = self
            .transcript_by_start
            .partition_point(|&index| transcripts[index].start < end);
        let mut overlapping = Vec::new();
        for position in (0..upper).rev() {
            budget.claim(1)?;
            if self.prefix_max_end[position] <= start {
                break;
            }
            let index = self.transcript_by_start[position];
            if transcripts[index].end > start {
                overlapping.push(index);
            }
        }
        Ok(overlapping)
    }

    fn overlapping_point(
        &self,
        transcripts: &[AnnotationTranscript],
        point: u32,
        budget: &AnnotationComparisonBudget,
    ) -> Result<Vec<usize>> {
        let upper = self
            .transcript_by_start
            .partition_point(|&index| transcripts[index].start <= point);
        let mut overlapping = Vec::new();
        for position in (0..upper).rev() {
            budget.claim(1)?;
            if self.prefix_max_end[position] < point {
                break;
            }
            let index = self.transcript_by_start[position];
            if transcripts[index].end >= point {
                overlapping.push(index);
            }
        }
        Ok(overlapping)
    }
}

#[derive(Debug)]
struct AnnotationIndex {
    by_chrom: Vec<Vec<AnnotationTranscript>>,
    intervals: Vec<AnnotationIntervalIndex>,
    matched_chroms: Vec<bool>,
    identity: anno::intent::AnnotationIdentity,
    content_blake3: String,
}

impl AnnotationIndex {
    fn chrom_is_matched(&self, collection: &Collection, chrom: u32, kind: &str) -> Result<bool> {
        let index = chrom as usize;
        collection
            .chroms
            .get(index)
            .with_context(|| format!("{kind} evidence references chromosome id {chrom} outside the collection dictionary"))?;
        self.matched_chroms.get(index).copied().with_context(|| {
            format!("annotation contig-match index is missing collection chromosome id {chrom}")
        })
    }

    fn interval_transcripts(
        &self,
        chrom: u32,
        start: u32,
        end: u32,
        budget: &AnnotationComparisonBudget,
    ) -> Result<Vec<&AnnotationTranscript>> {
        let chrom = chrom as usize;
        self.intervals[chrom]
            .overlapping_interval(&self.by_chrom[chrom], start, end, budget)
            .map(|indices| {
                indices
                    .into_iter()
                    .map(|index| &self.by_chrom[chrom][index])
                    .collect()
            })
    }

    fn point_transcripts(
        &self,
        chrom: u32,
        point: u32,
        budget: &AnnotationComparisonBudget,
    ) -> Result<Vec<&AnnotationTranscript>> {
        let chrom = chrom as usize;
        self.intervals[chrom]
            .overlapping_point(&self.by_chrom[chrom], point, budget)
            .map(|indices| {
                indices
                    .into_iter()
                    .map(|index| &self.by_chrom[chrom][index])
                    .collect()
            })
    }
}

#[derive(Debug, Default)]
struct AnnotationCandidateExclusions {
    candidates: usize,
    unmatched_chroms: BTreeSet<u32>,
}

fn exclude_annotation_unmatched_candidates(
    collection: &Collection,
    annotation: &AnnotationIndex,
    candidates: &mut Vec<Candidate>,
) -> Result<AnnotationCandidateExclusions> {
    let mut retained = Vec::with_capacity(candidates.len());
    let mut excluded = AnnotationCandidateExclusions::default();
    for candidate in candidates.drain(..) {
        if annotation.chrom_is_matched(collection, candidate.key.chrom, "splice-event")? {
            retained.push(candidate);
        } else {
            excluded.candidates += 1;
            excluded.unmatched_chroms.insert(candidate.key.chrom);
        }
    }
    *candidates = retained;
    Ok(excluded)
}

fn load_annotation(
    path: &Path,
    collection: &Collection,
    identity: anno::intent::AnnotationIdentity,
) -> Result<AnnotationIndex> {
    let bound = anno::intent::BoundAnnotation::from_path(path, identity)
        .with_context(|| format!("loading bound annotation {}", path.display()))?;
    let identity = bound.identity().clone();
    let content_blake3 = bound
        .identity()
        .digest
        .clone()
        .context("bound annotation lacks its observed content digest")?;
    let annotation = bound.annotation();
    let annotation_chrom_of: Vec<Option<u32>> = collection
        .chroms
        .iter()
        .map(|name| annotation.chrom_ids.get(name).copied())
        .collect();
    let matched_chroms: Vec<bool> = annotation_chrom_of.iter().map(Option::is_some).collect();
    if !matched_chroms.iter().any(|matched| *matched) {
        bail!(
            "annotation {} shares no exact contig names with the collection (for example, `1` and `chr1` are not treated as equivalent)",
            path.display()
        );
    }
    let mut collection_chrom_of = FxHashMap::default();
    for (collection_chrom, annotation_chrom) in annotation_chrom_of.iter().enumerate() {
        if let Some(annotation_chrom) = annotation_chrom {
            collection_chrom_of.insert(*annotation_chrom, collection_chrom);
        }
    }
    let mut by_chrom: Vec<Vec<AnnotationTranscript>> =
        (0..collection.chroms.len()).map(|_| Vec::new()).collect();
    for transcript in &annotation.transcripts {
        let Some(&collection_chrom) = collection_chrom_of.get(&transcript.chrom) else {
            continue;
        };
        let (start, end) = transcript.span();
        let mut junctions = FxHashSet::default();
        let mut donors = FxHashSet::default();
        let mut acceptors = FxHashSet::default();
        for pair in transcript.exons.windows(2) {
            junctions.insert((pair[0].end, pair[1].start));
            donors.insert(pair[0].end);
            acceptors.insert(pair[1].start);
        }
        by_chrom[collection_chrom].push(AnnotationTranscript {
            gene_id: annotation.gene_ids[transcript.gene as usize].clone(),
            strand_rev: transcript.strand_rev,
            start,
            end,
            junctions,
            donors,
            acceptors,
        });
    }
    let intervals = by_chrom
        .iter()
        .map(|transcripts| AnnotationIntervalIndex::new(transcripts))
        .collect();
    Ok(AnnotationIndex {
        by_chrom,
        intervals,
        matched_chroms,
        identity,
        content_blake3,
    })
}

fn transcript_accepted(
    solo_strand: anno::assign::SoloStrand,
    observed_strands: [u32; 2],
    transcript_rev: bool,
) -> bool {
    (observed_strands[0] > 0 && solo_strand.accepts(false, transcript_rev))
        || (observed_strands[1] > 0 && solo_strand.accepts(true, transcript_rev))
}

fn side_components(key: &EntityKey) -> Vec<Vec<Coordinate>> {
    if key.kind.is_junction() {
        return vec![key
            .components
            .iter()
            .map(|component| component.coordinate)
            .collect()];
    }
    let include: Vec<Coordinate> = key
        .components
        .iter()
        .filter(|component| component.side == ComponentSide::Include)
        .map(|component| component.coordinate)
        .collect();
    let exclude: Vec<Coordinate> = key
        .components
        .iter()
        .filter(|component| component.side == ComponentSide::Exclude)
        .map(|component| component.coordinate)
        .collect();
    vec![include, exclude]
}

fn assign_gap_primary(gap: &mut GapClassification) {
    gap.primary_class = if gap.strand {
        Some("strand")
    } else if gap.boundary {
        Some("boundary")
    } else if gap.overlap {
        Some("overlap")
    } else if gap.missing_junction {
        Some("missing_junction")
    } else {
        None
    };
}

fn classify_gap(
    annotation: &AnnotationIndex,
    candidate: &Candidate,
    exact: &EntityExact,
    solo_strand: anno::assign::SoloStrand,
    budget: &AnnotationComparisonBudget,
) -> Result<GapClassification> {
    let start = candidate
        .key
        .components
        .iter()
        .map(|component| component.coordinate.donor)
        .min()
        .context("annotation classification candidate has no components")?;
    let end = candidate
        .key
        .components
        .iter()
        .map(|component| component.coordinate.acceptor)
        .max()
        .context("annotation classification candidate has no components")?;
    let transcripts = annotation.interval_transcripts(candidate.key.chrom, start, end, budget)?;
    let comparison_passes = candidate
        .key
        .components
        .len()
        .checked_add(side_components(&candidate.key).len())
        .and_then(|value| value.checked_add(3))
        .context("annotation comparison-pass count overflow")?;
    budget.claim(
        transcripts
            .len()
            .checked_mul(comparison_passes)
            .context("annotation comparison count overflow")?,
    )?;
    let accepted = |transcript: &AnnotationTranscript| {
        transcript_accepted(solo_strand, exact.strand_umis, transcript.strand_rev)
    };
    let mut incompatible_path_overlap = false;
    let mut side_match_sets = Vec::new();
    for side in side_components(&candidate.key) {
        let side_matches: Vec<usize> = transcripts
            .iter()
            .enumerate()
            .filter(|(_, transcript)| {
                accepted(transcript)
                    && side.iter().all(|coordinate| {
                        transcript
                            .junctions
                            .contains(&(coordinate.donor, coordinate.acceptor))
                    })
            })
            .map(|(index, _)| index)
            .collect();
        if side_matches.is_empty()
            && side.iter().all(|coordinate| {
                transcripts.iter().any(|transcript| {
                    accepted(transcript)
                        && transcript
                            .junctions
                            .contains(&(coordinate.donor, coordinate.acceptor))
                })
            })
        {
            incompatible_path_overlap = true;
        }
        side_match_sets.push(side_matches);
    }

    let every_side_has_matches = side_match_sets.iter().all(|matches| !matches.is_empty());
    let mut common_genes: Option<BTreeSet<&str>> = None;
    for matches in &side_match_sets {
        let side_genes: BTreeSet<&str> = matches
            .iter()
            .map(|&index| transcripts[index].gene_id.as_str())
            .collect();
        if let Some(common) = common_genes.as_mut() {
            common.retain(|gene| side_genes.contains(gene));
        } else {
            common_genes = Some(side_genes);
        }
    }
    let common_genes = common_genes.unwrap_or_default();
    let event_compatible = every_side_has_matches && !common_genes.is_empty();
    let cross_gene_overlap = every_side_has_matches && common_genes.is_empty();
    let compatible: FxHashSet<usize> = if event_compatible {
        side_match_sets
            .iter()
            .flatten()
            .copied()
            .filter(|&index| common_genes.contains(transcripts[index].gene_id.as_str()))
            .collect()
    } else {
        FxHashSet::default()
    };

    let overlapping_gene_ids: Vec<String> = transcripts
        .iter()
        .filter(|transcript| {
            candidate.key.components.iter().any(|component| {
                transcript.start < component.coordinate.acceptor
                    && transcript.end > component.coordinate.donor
            })
        })
        .map(|transcript| transcript.gene_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    let mut gap = GapClassification {
        incompatible_with_every_transcript: !event_compatible,
        compatible_transcripts: compatible.len(),
        overlap: incompatible_path_overlap || cross_gene_overlap,
        overlapping_gene_ids,
        ..GapClassification::default()
    };
    if event_compatible {
        return Ok(gap);
    }

    for component in &candidate.key.components {
        let coordinate = component.coordinate;
        let mut accepted_exact = false;
        let mut accepted_donor = false;
        let mut accepted_acceptor = false;
        let mut accepted_overlap = false;
        let mut opposite_exact = false;
        let mut opposite_overlap = false;
        for transcript in &transcripts {
            let overlaps =
                transcript.start < coordinate.acceptor && transcript.end > coordinate.donor;
            if accepted(transcript) {
                accepted_exact |= transcript
                    .junctions
                    .contains(&(coordinate.donor, coordinate.acceptor));
                accepted_donor |= transcript.donors.contains(&coordinate.donor);
                accepted_acceptor |= transcript.acceptors.contains(&coordinate.acceptor);
                accepted_overlap |= overlaps;
            } else {
                opposite_exact |= transcript
                    .junctions
                    .contains(&(coordinate.donor, coordinate.acceptor));
                opposite_overlap |= overlaps;
            }
        }
        if !accepted_exact {
            if accepted_donor && accepted_acceptor {
                gap.missing_junction = true;
            } else {
                gap.boundary = true;
            }
            gap.strand |= opposite_exact || (opposite_overlap && !accepted_overlap);
        }
    }
    assign_gap_primary(&mut gap);
    Ok(gap)
}

fn classify_terminal_anchor(
    annotation: &AnnotationIndex,
    chrom: u32,
    strand_rev: bool,
    anchor: u32,
    budget: &AnnotationComparisonBudget,
) -> Result<GapClassification> {
    let transcripts = annotation.point_transcripts(chrom, anchor, budget)?;
    budget.claim(
        transcripts
            .len()
            .checked_mul(2)
            .context("terminal annotation comparison count overflow")?,
    )?;
    let mut compatible = FxHashSet::default();
    let mut accepted_overlap = false;
    let mut opposite_exact = false;
    let mut opposite_overlap = false;
    let mut genes = BTreeSet::new();
    for (index, transcript) in transcripts.iter().enumerate() {
        let terminal = if transcript.strand_rev {
            transcript.start
        } else {
            transcript.end
        };
        let overlaps = transcript.start <= anchor && anchor <= transcript.end;
        let exact = anchor == terminal;
        if overlaps || exact {
            genes.insert(transcript.gene_id.clone());
        }
        if transcript.strand_rev == strand_rev {
            accepted_overlap |= overlaps;
            if exact {
                compatible.insert(index);
            }
        } else {
            opposite_exact |= exact;
            opposite_overlap |= overlaps;
        }
    }
    let incompatible = compatible.is_empty();
    let mut gap = GapClassification {
        incompatible_with_every_transcript: incompatible,
        compatible_transcripts: compatible.len(),
        boundary: incompatible,
        strand: incompatible && (opposite_exact || (opposite_overlap && !accepted_overlap)),
        overlap: incompatible && accepted_overlap,
        overlapping_gene_ids: genes.into_iter().collect(),
        ..GapClassification::default()
    };
    assign_gap_primary(&mut gap);
    Ok(gap)
}

fn classify_terminal_gap(
    annotation: &AnnotationIndex,
    entity: &TailEntity,
    budget: &AnnotationComparisonBudget,
) -> Result<GapClassification> {
    let transcripts = annotation.interval_transcripts(
        entity.chrom,
        entity.start.saturating_sub(1),
        entity.end.saturating_add(1),
        budget,
    )?;
    budget.claim(
        transcripts
            .len()
            .checked_mul(entity.anchors.len().saturating_add(2))
            .context("terminal annotation comparison count overflow")?,
    )?;
    let exact_anchors: FxHashSet<u32> = entity.anchors.iter().map(|anchor| anchor.anchor).collect();
    let mut compatible = FxHashSet::default();
    let mut accepted_overlap = false;
    let mut opposite_exact = false;
    let mut opposite_overlap = false;
    let mut genes = BTreeSet::new();
    for (index, transcript) in transcripts.iter().enumerate() {
        let terminal = if transcript.strand_rev {
            transcript.start
        } else {
            transcript.end
        };
        let overlaps = exact_anchors
            .iter()
            .any(|&anchor| transcript.start <= anchor && anchor <= transcript.end);
        let exact = exact_anchors.contains(&terminal);
        if overlaps || exact {
            genes.insert(transcript.gene_id.clone());
        }
        if transcript.strand_rev == entity.strand_rev {
            accepted_overlap |= overlaps;
            if exact {
                compatible.insert(index);
            }
        } else {
            opposite_exact |= exact;
            opposite_overlap |= overlaps;
        }
    }
    let incompatible = compatible.is_empty();
    let mut gap = GapClassification {
        incompatible_with_every_transcript: incompatible,
        compatible_transcripts: compatible.len(),
        boundary: incompatible,
        strand: incompatible && (opposite_exact || (opposite_overlap && !accepted_overlap)),
        overlap: incompatible && accepted_overlap,
        overlapping_gene_ids: genes.into_iter().collect(),
        ..GapClassification::default()
    };
    assign_gap_primary(&mut gap);
    Ok(gap)
}

/// `required_metrics` holds this entity's informative metric in each `--require-group`, in the
/// order the groups were requested.
fn evidence_passes_predicates(
    candidate: &Candidate,
    exact: &EntityExact,
    required_metrics: &[u32],
    args: &Args,
) -> bool {
    let side_support = candidate.key.kind.is_junction()
        || (exact.counts.include_only >= args.min_side_umis
            && exact.counts.exclude_only >= args.min_side_umis
            && exact
                .component_umi_classes
                .iter()
                .take(candidate.key.components.len())
                .all(|&count| count as usize >= args.min_side_umis));
    exact.counts.metric(candidate.key.kind) >= args.min_umis
        && side_support
        && exact.samples as usize >= args.min_samples
        && exact.donors as usize >= args.min_donors
        && required_metrics
            .iter()
            .all(|&metric| metric as usize >= args.min_group_umis)
}

fn keep_entity(
    candidate: &Candidate,
    exact: &EntityExact,
    required_metrics: &[u32],
    gap: &GapClassification,
    args: &Args,
) -> bool {
    evidence_passes_predicates(candidate, exact, required_metrics, args)
        && (!args.novel_only || gap.incompatible_with_every_transcript)
}

/// Per-stage wall clock, recorded in execution order. The window matches `total_seconds`: both
/// are sampled before the result tables are streamed, so `output` covers assembling the result
/// summary aggregates rather than serializing the tables.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize)]
struct StageSeconds {
    load_catalogue: f64,
    discover_candidates: f64,
    route_candidates: f64,
    exact_counting: f64,
    terminal_tails: f64,
    annotation_classification: f64,
    output: f64,
}

/// Process peak resident set size observed at the end of each stage. The kernel reports a high
/// water mark, so these values are monotonically non-decreasing and attribute a peak to the last
/// stage that could have caused it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize)]
struct StagePeakRssBytes {
    load_catalogue: u64,
    discover_candidates: u64,
    route_candidates: u64,
    exact_counting: u64,
    terminal_tails: u64,
    annotation_classification: u64,
    output: u64,
}

/// Return freed heap pages to the operating system. Reducing one cohort archive allocates and
/// frees hundreds of megabytes of packed hits; without this the allocator keeps those pages, and
/// the next archive's peak is measured on top of them. This is advisory and a no-op elsewhere.
fn release_free_heap() {
    #[cfg(target_env = "gnu")]
    // SAFETY: `malloc_trim` takes no pointers and only releases already-free heap pages.
    unsafe {
        libc::malloc_trim(0);
    }
}

/// Current process peak RSS in bytes, or zero where the platform does not publish one. Reading
/// `VmHWM` costs one small file read per stage boundary and never fails the search.
fn peak_rss_bytes() -> u64 {
    let Ok(status) = std::fs::read_to_string("/proc/self/status") else {
        return 0;
    };
    status
        .lines()
        .find_map(|line| {
            let rest = line.strip_prefix("VmHWM:")?;
            let kib: u64 = rest.split_whitespace().next()?.parse().ok()?;
            kib.checked_mul(1024)
        })
        .unwrap_or(0)
}

/// Stopwatch that attributes each elapsed interval to one named stage.
struct StageClock {
    last: std::time::Instant,
    seconds: StageSeconds,
    peak_rss: StagePeakRssBytes,
}

impl StageClock {
    fn new(started: std::time::Instant) -> Self {
        Self {
            last: started,
            seconds: StageSeconds::default(),
            peak_rss: StagePeakRssBytes::default(),
        }
    }

    fn lap(
        &mut self,
        seconds: fn(&mut StageSeconds) -> &mut f64,
        rss: fn(&mut StagePeakRssBytes) -> &mut u64,
    ) {
        let now = std::time::Instant::now();
        *seconds(&mut self.seconds) = now.duration_since(self.last).as_secs_f64();
        *rss(&mut self.peak_rss) = peak_rss_bytes();
        self.last = now;
    }
}

#[derive(Serialize)]
struct SearchSummary<'a> {
    coordinates: &'static str,
    candidate_source: &'static str,
    exact_counting: &'static str,
    missing_count_rows: &'static str,
    event_usage: &'static str,
    alternative_site_side_semantics: &'static str,
    cassette_component_semantics: &'static str,
    donor_semantics: &'static str,
    cell_scope: &'static str,
    evidence_placement_policy: &'static str,
    multimapper_placements_included: bool,
    multimapper_alternatives_available_to_search: bool,
    requested_kinds: Vec<&'static str>,
    catalogue_junctions: u64,
    splice_candidate_definitions_attempted: u64,
    splice_candidate_definitions_distinct: u64,
    candidate_entities: u64,
    retained_entities: u64,
    archives_total: u64,
    archives_opened: u64,
    donors_total: u64,
    groups_total: u64,
    min_support: u64,
    min_samples: u64,
    min_donors: u64,
    min_umi_classes: u64,
    min_side_umi_classes: u64,
    min_group_umi_classes: u64,
    max_candidates: u64,
    max_candidates_considered: u64,
    max_routed_entries: u64,
    routed_target_associations: u64,
    routed_chunk_postings: u64,
    max_exact_match_attempts: u64,
    exact_match_attempts: u64,
    max_annotation_comparisons: u64,
    annotation_comparisons: u64,
    terminal_cluster_bp: u64,
    max_terminal_events: u64,
    required_groups: &'a [String],
    annotation_gap_flags_nonexclusive: bool,
    annotation_gap_primary_precedence: [&'static str; 4],
    annotation_identity: Option<&'a anno::intent::AnnotationIdentity>,
    annotation_collection_compatibility: Option<&'static str>,
    annotation_exact_contig_names_matched: Option<u64>,
    annotation_collection_contigs: Option<u64>,
    annotation_unmatched_evidence_policy: Option<&'static str>,
    annotation_excluded_splice_candidates: u64,
    annotation_excluded_terminal_routes: u64,
    annotation_excluded_terminal_events: u64,
    annotation_unmatched_evidence_contigs: Vec<String>,
    collection_genome_algo: Option<&'a str>,
    collection_genome_digest: Option<&'a str>,
    unique_chunks_decoded: u64,
    independent_chunk_decodes: u64,
    planned_compressed_bytes: u64,
    actual_archive_bytes_read: u64,
    source_archive_identity_bytes_read: u64,
    collection_sidecar_bytes_read: u64,
    total_seconds: f64,
    stage_seconds: StageSeconds,
    stage_peak_rss_bytes: StagePeakRssBytes,
    terminal_tail_available_archives: u64,
    terminal_tail_unavailable_archives: u64,
    terminal_tail_available_donors: u64,
    terminal_tail_declared_selected_molecules: u64,
    terminal_tail_declared_events: u64,
    terminal_tail_routed_chunks: u64,
}

fn capability_schema() -> std::result::Result<TableSchema, OutputError> {
    set_table_schema(
        CAPABILITY_SCHEMA,
        vec![
            Field::new("kind", DataType::String),
            Field::new("scope", DataType::String),
            Field::new("sample", DataType::String),
            Field::new("donor", DataType::String),
            Field::new("status", DataType::String),
            Field::new("archives_available", DataType::UInt64).nullable(),
            Field::new("archives_unavailable", DataType::UInt64).nullable(),
            Field::new("included_in_denominator", DataType::Boolean).nullable(),
            Field::new("denominator", DataType::String),
        ],
        &["kind", "scope", "sample"],
    )
}

fn entity_schema() -> std::result::Result<TableSchema, OutputError> {
    set_table_schema(
        ENTITY_SCHEMA,
        vec![
            Field::new("entity_id", DataType::String),
            Field::new("kind", DataType::String),
            Field::new("chrom", DataType::String),
            Field::new("strand", DataType::String),
            Field::new("start", DataType::UInt64),
            Field::new("end", DataType::UInt64),
            Field::new("summit", DataType::UInt64).nullable(),
            Field::new(
                "catalogue_min_component_route_upper_bound",
                DataType::UInt64,
            )
            .nullable(),
            Field::new("catalogue_sample_route_upper_bound", DataType::UInt64).nullable(),
            Field::new("catalogue_donor_route_upper_bound", DataType::UInt64).nullable(),
            Field::new("exact_umi_classes", DataType::UInt64),
            Field::new("exact_cells", DataType::UInt64),
            Field::new("exact_samples", DataType::UInt64),
            Field::new("exact_donors", DataType::UInt64),
            Field::new("forward_alignment_umi_classes", DataType::UInt64),
            Field::new("reverse_alignment_umi_classes", DataType::UInt64),
            Field::new("annotation_incompatible", DataType::Boolean).nullable(),
            Field::new("compatible_transcripts", DataType::UInt64).nullable(),
            Field::new("gap_primary_class", DataType::String).nullable(),
            Field::new("gap_missing_junction", DataType::Boolean).nullable(),
            Field::new("gap_boundary", DataType::Boolean).nullable(),
            Field::new("gap_strand", DataType::Boolean).nullable(),
            Field::new("gap_overlap", DataType::Boolean).nullable(),
            Field::new("overlapping_gene_ids", DataType::Json).nullable(),
        ],
        &["entity_id"],
    )
}

fn component_schema() -> std::result::Result<TableSchema, OutputError> {
    set_table_schema(
        COMPONENT_SCHEMA,
        vec![
            Field::new("entity_id", DataType::String),
            Field::new("role", DataType::String),
            Field::new("ordinal", DataType::UInt64),
            Field::new("chrom", DataType::String),
            Field::new("coordinate_kind", DataType::String),
            Field::new("position", DataType::UInt64).nullable(),
            Field::new("donor", DataType::UInt64).nullable(),
            Field::new("acceptor", DataType::UInt64).nullable(),
            Field::new("exact_umi_classes", DataType::UInt64),
        ],
        &["entity_id", "role", "ordinal"],
    )
}

fn terminal_anchor_schema() -> std::result::Result<TableSchema, OutputError> {
    set_table_schema(
        TERMINAL_ANCHOR_SCHEMA,
        vec![
            Field::new("entity_id", DataType::String),
            Field::new("chrom", DataType::String),
            Field::new("strand", DataType::String),
            Field::new("anchor", DataType::UInt64),
            Field::new("exact_umi_classes", DataType::UInt64),
            Field::new("exact_cells", DataType::UInt64),
            Field::new("exact_samples", DataType::UInt64),
            Field::new("exact_donors", DataType::UInt64),
            Field::new("max_clip_len", DataType::UInt64),
            Field::new("max_tail_bases", DataType::UInt64),
            Field::new("max_terminal_run", DataType::UInt64),
            Field::new("annotation_incompatible", DataType::Boolean).nullable(),
            Field::new("compatible_transcripts", DataType::UInt64).nullable(),
            Field::new("gap_primary_class", DataType::String).nullable(),
            Field::new("gap_missing_junction", DataType::Boolean).nullable(),
            Field::new("gap_boundary", DataType::Boolean).nullable(),
            Field::new("gap_strand", DataType::Boolean).nullable(),
            Field::new("gap_overlap", DataType::Boolean).nullable(),
            Field::new("overlapping_gene_ids", DataType::Json).nullable(),
        ],
        &["entity_id", "anchor"],
    )
}

fn terminal_count_schema() -> std::result::Result<TableSchema, OutputError> {
    set_table_schema(
        TERMINAL_COUNT_SCHEMA,
        vec![
            Field::new("entity_id", DataType::String),
            Field::new("anchor", DataType::UInt64),
            Field::new("sample", DataType::String),
            Field::new("donor", DataType::String),
            Field::new("group", DataType::String),
            Field::new("umi_classes", DataType::UInt64),
            Field::new("cells", DataType::UInt64),
        ],
        &["entity_id", "anchor", "sample", "group"],
    )
}

fn count_schema() -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(
        COUNT_SCHEMA,
        vec![
            Field::new("entity_id", DataType::String),
            Field::new("sample", DataType::String),
            Field::new("donor", DataType::String),
            Field::new("group", DataType::String),
            Field::new("support_umi_classes", DataType::UInt64).nullable(),
            Field::new("side_a_only_umi_classes", DataType::UInt64).nullable(),
            Field::new("side_b_only_umi_classes", DataType::UInt64).nullable(),
            Field::new("include_only_umi_classes", DataType::UInt64).nullable(),
            Field::new("exclude_only_umi_classes", DataType::UInt64).nullable(),
            Field::new("both_umi_classes", DataType::UInt64).nullable(),
            Field::new("informative_umi_classes", DataType::UInt64),
            Field::new("cells", DataType::UInt64),
        ],
    )?
    .with_semantics(TableSemantics::new(RowSemantics::Set).with_key([
        "entity_id",
        "sample",
        "group",
    ]))
}

struct OutputData<'a> {
    collection: &'a Collection,
    candidates: &'a [Candidate],
    exact: &'a [EntityExact],
    gap: &'a [Option<GapClassification>],
    retained: &'a [usize],
    archives: &'a [ArchiveExact],
    design: &'a Design,
    groups: &'a Groups,
    tails: &'a [TailEntity],
    tail_capability: &'a TailCapabilitySummary,
}

fn stream_result<W: Write>(
    writer: W,
    format: CollectionOutputFormat,
    context: &ResultContext,
    summary: &SearchSummary<'_>,
    data: &OutputData<'_>,
) -> std::result::Result<W, OutputError> {
    let capability_schema = capability_schema()?;
    let entity_schema = entity_schema()?;
    let component_schema = component_schema()?;
    let count_schema = count_schema()?;
    let terminal_anchor_schema = terminal_anchor_schema()?;
    let terminal_count_schema = terminal_count_schema()?;
    let retained_set: FxHashSet<usize> = data.retained.iter().copied().collect();
    let component_count: usize = data
        .retained
        .iter()
        .map(|&entity| data.candidates[entity].key.components.len())
        .sum::<usize>()
        + data
            .tails
            .iter()
            .map(|entity| entity.anchors.len())
            .sum::<usize>();
    let count_rows: usize = data
        .archives
        .iter()
        .flat_map(|archive| &archive.rows)
        .filter(|row| retained_set.contains(&(row.entity as usize)))
        .count()
        + data
            .tails
            .iter()
            .map(|entity| entity.counts.len())
            .sum::<usize>();
    let terminal_anchor_rows: usize = data.tails.iter().map(|entity| entity.anchors.len()).sum();
    let terminal_count_rows: usize = data
        .tails
        .iter()
        .flat_map(|entity| &entity.anchors)
        .map(|anchor| anchor.counts.len())
        .sum();
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        RESULT_SCHEMA,
        OutputFormat::from(format),
        context,
        summary,
    )?;
    bundle.write_table(
        "capabilities",
        &capability_schema,
        Some(&SelectionSummary::complete(
            (2 + data.tail_capability.archive_available.len()) as u64,
        )),
        |rows| {
            rows.write_row_with(|row| {
                row.string("junction")?;
                row.string("aggregate")?;
                row.string("*")?;
                row.string("*")?;
                row.string("available")?;
                row.uint64(data.collection.archives.len() as u64)?;
                row.uint64(0)?;
                row.null()?;
                row.string("all collection archives")?;
                Ok(())
            })?;
            rows.write_row_with(|row| {
                row.string("terminal_tail")?;
                row.string("aggregate")?;
                row.string("*")?;
                row.string("*")?;
                row.string(if !data.tail_capability.requested {
                    "not_requested"
                } else if data.tail_capability.available_archives == 0 {
                    "unavailable"
                } else if data.tail_capability.unavailable_archives == 0 {
                    "available"
                } else {
                    "partially_available"
                })?;
                row.uint64(data.tail_capability.available_archives as u64)?;
                row.uint64(data.tail_capability.unavailable_archives as u64)?;
                row.null()?;
                row.string(if data.tail_capability.requested {
                    "only archives declaring the typed terminal-tail capability; unsupported archives are unavailable, never zero"
                } else {
                    "not evaluated because terminal-tail was not requested"
                })?;
                Ok(())
            })?;
            for (sample, &available) in data
                .tail_capability
                .archive_available
                .iter()
                .enumerate()
            {
                rows.write_row_with(|row| {
                    row.string("terminal_tail")?;
                    row.string("archive")?;
                    row.string(&data.collection.archives[sample].id)?;
                    row.string(
                        &data.design.donor_names[data.design.donor_of_sample[sample]],
                    )?;
                    row.string(if available { "available" } else { "unavailable" })?;
                    row.null()?;
                    row.null()?;
                    row.boolean(available)?;
                    row.string(if available {
                        "included in terminal-tail recurrence and zero-count denominator"
                    } else {
                        "excluded from terminal-tail recurrence and zero-count denominator"
                    })?;
                    Ok(())
                })?;
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "entities",
        &entity_schema,
        Some(&SelectionSummary::complete(
            (data.retained.len() + data.tails.len()) as u64,
        )),
        |rows| {
            for &entity in data.retained {
                let candidate = &data.candidates[entity];
                let exact = &data.exact[entity];
                let gap = data.gap[entity].as_ref();
                rows.write_row_with(|row| {
                    row.string(&candidate.key.id(&data.collection.chroms))?;
                    row.string(candidate.key.kind.name())?;
                    row.string(&data.collection.chroms[candidate.key.chrom as usize])?;
                    row.string(match candidate.key.reported_strand_rev {
                        Some(true) => "-",
                        Some(false) => "+",
                        None => ".",
                    })?;
                    let start = candidate
                        .key
                        .components
                        .iter()
                        .map(|component| component.coordinate.donor)
                        .min()
                        .unwrap_or(0);
                    let end = candidate
                        .key
                        .components
                        .iter()
                        .map(|component| component.coordinate.acceptor)
                        .max()
                        .unwrap_or(start);
                    row.uint64(start as u64)?;
                    row.uint64(end as u64)?;
                    row.null()?;
                    row.uint64(candidate.catalogue_support_upper_bound)?;
                    row.uint64(candidate.catalogue_samples as u64)?;
                    row.uint64(candidate.catalogue_donors as u64)?;
                    row.uint64(exact.counts.metric(candidate.key.kind) as u64)?;
                    row.uint64(u64::from(exact.cells))?;
                    row.uint64(u64::from(exact.samples))?;
                    row.uint64(u64::from(exact.donors))?;
                    row.uint64(u64::from(exact.strand_umis[0]))?;
                    row.uint64(u64::from(exact.strand_umis[1]))?;
                    if let Some(gap) = gap {
                        row.boolean(gap.incompatible_with_every_transcript)?;
                        row.uint64(gap.compatible_transcripts as u64)?;
                        if let Some(primary) = gap.primary_class {
                            row.string(primary)?;
                        } else {
                            row.null()?;
                        }
                        row.boolean(gap.missing_junction)?;
                        row.boolean(gap.boundary)?;
                        row.boolean(gap.strand)?;
                        row.boolean(gap.overlap)?;
                        row.json(&json!(gap.overlapping_gene_ids))?;
                    } else {
                        for _ in 0..8 {
                            row.null()?;
                        }
                    }
                    Ok(())
                })?;
            }
            for entity in data.tails {
                let gap = entity.gap.as_ref();
                rows.write_row_with(|row| {
                    row.string(&entity.id)?;
                    row.string("terminal_tail")?;
                    row.string(&data.collection.chroms[entity.chrom as usize])?;
                    row.string(if entity.strand_rev { "-" } else { "+" })?;
                    row.uint64(entity.start as u64)?;
                    row.uint64(entity.end as u64)?;
                    row.uint64(entity.summit as u64)?;
                    row.null()?;
                    row.null()?;
                    row.null()?;
                    row.uint64(entity.umis as u64)?;
                    row.uint64(entity.cells as u64)?;
                    row.uint64(entity.samples as u64)?;
                    row.uint64(entity.donors as u64)?;
                    row.uint64(if entity.strand_rev {
                        0
                    } else {
                        entity.umis as u64
                    })?;
                    row.uint64(if entity.strand_rev {
                        entity.umis as u64
                    } else {
                        0
                    })?;
                    if let Some(gap) = gap {
                        row.boolean(gap.incompatible_with_every_transcript)?;
                        row.uint64(gap.compatible_transcripts as u64)?;
                        if let Some(primary) = gap.primary_class {
                            row.string(primary)?;
                        } else {
                            row.null()?;
                        }
                        row.boolean(gap.missing_junction)?;
                        row.boolean(gap.boundary)?;
                        row.boolean(gap.strand)?;
                        row.boolean(gap.overlap)?;
                        row.json(&json!(gap.overlapping_gene_ids))?;
                    } else {
                        for _ in 0..8 {
                            row.null()?;
                        }
                    }
                    Ok(())
                })?;
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "components",
        &component_schema,
        Some(&SelectionSummary::complete(component_count as u64)),
        |rows| {
            for &entity in data.retained {
                let candidate = &data.candidates[entity];
                let id = candidate.key.id(&data.collection.chroms);
                let mut ordinal_of: BTreeMap<ComponentSide, usize> = BTreeMap::new();
                for (component_index, component) in candidate.key.components.iter().enumerate() {
                    let ordinal = ordinal_of.entry(component.side).or_insert(0);
                    rows.write_row_with(|row| {
                        row.string(&id)?;
                        row.string(component.side.output_name(candidate.key.kind))?;
                        row.uint64(*ordinal as u64)?;
                        row.string(&data.collection.chroms[component.coordinate.chrom as usize])?;
                        row.string("junction")?;
                        row.null()?;
                        row.uint64(component.coordinate.donor as u64)?;
                        row.uint64(component.coordinate.acceptor as u64)?;
                        row.uint64(u64::from(
                            data.exact[entity].component_umi_classes[component_index],
                        ))?;
                        Ok(())
                    })?;
                    *ordinal += 1;
                }
            }
            for entity in data.tails {
                for (ordinal, anchor) in entity.anchors.iter().enumerate() {
                    rows.write_row_with(|row| {
                        row.string(&entity.id)?;
                        row.string("support")?;
                        row.uint64(ordinal as u64)?;
                        row.string(&data.collection.chroms[entity.chrom as usize])?;
                        row.string("terminal_anchor")?;
                        row.uint64(anchor.anchor as u64)?;
                        row.null()?;
                        row.null()?;
                        row.uint64(anchor.umis as u64)?;
                        Ok(())
                    })?;
                }
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "counts",
        &count_schema,
        Some(&SelectionSummary::complete(count_rows as u64)),
        |rows| {
            for archive in data.archives {
                let donor = data.design.donor_of_sample[archive.sample];
                for count in &archive.rows {
                    if !retained_set.contains(&(count.entity as usize)) {
                        continue;
                    }
                    let candidate = &data.candidates[count.entity as usize];
                    rows.write_row_with(|row| {
                        row.string(&candidate.key.id(&data.collection.chroms))?;
                        row.string(&data.collection.archives[archive.sample].id)?;
                        row.string(&data.design.donor_names[donor])?;
                        row.string(&data.groups.names[count.group as usize])?;
                        if candidate.key.kind.is_junction() {
                            row.uint64(u64::from(count.counts.support))?;
                            for _ in 0..5 {
                                row.null()?;
                            }
                        } else if candidate.key.kind.is_alternative_site() {
                            row.null()?;
                            row.uint64(u64::from(count.counts.include_only))?;
                            row.uint64(u64::from(count.counts.exclude_only))?;
                            row.null()?;
                            row.null()?;
                            row.uint64(u64::from(count.counts.both))?;
                        } else {
                            row.null()?;
                            row.null()?;
                            row.null()?;
                            row.uint64(u64::from(count.counts.include_only))?;
                            row.uint64(u64::from(count.counts.exclude_only))?;
                            row.uint64(u64::from(count.counts.both))?;
                        }
                        row.uint64(u64::from(count.counts.metric(candidate.key.kind)))?;
                        row.uint64(u64::from(count.cells))?;
                        Ok(())
                    })?;
                }
            }
            for entity in data.tails {
                for count in &entity.counts {
                    rows.write_row_with(|row| {
                        row.string(&entity.id)?;
                        row.string(&data.collection.archives[count.sample].id)?;
                        row.string(&data.design.donor_names[count.donor])?;
                        row.string(&data.groups.names[count.group])?;
                        row.uint64(count.umis as u64)?;
                        for _ in 0..5 {
                            row.null()?;
                        }
                        row.uint64(count.umis as u64)?;
                        row.uint64(count.cells as u64)?;
                        Ok(())
                    })?;
                }
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "terminal_anchors",
        &terminal_anchor_schema,
        Some(&SelectionSummary::complete(terminal_anchor_rows as u64)),
        |rows| {
            for entity in data.tails {
                for anchor in &entity.anchors {
                    rows.write_row_with(|row| {
                        row.string(&entity.id)?;
                        row.string(&data.collection.chroms[entity.chrom as usize])?;
                        row.string(if entity.strand_rev { "-" } else { "+" })?;
                        row.uint64(anchor.anchor as u64)?;
                        row.uint64(anchor.umis as u64)?;
                        row.uint64(anchor.cells as u64)?;
                        row.uint64(anchor.samples as u64)?;
                        row.uint64(anchor.donors as u64)?;
                        row.uint64(anchor.max_clip_len as u64)?;
                        row.uint64(anchor.max_tail_bases as u64)?;
                        row.uint64(anchor.max_terminal_run as u64)?;
                        if let Some(gap) = anchor.gap.as_ref() {
                            row.boolean(gap.incompatible_with_every_transcript)?;
                            row.uint64(gap.compatible_transcripts as u64)?;
                            if let Some(primary) = gap.primary_class {
                                row.string(primary)?;
                            } else {
                                row.null()?;
                            }
                            row.boolean(gap.missing_junction)?;
                            row.boolean(gap.boundary)?;
                            row.boolean(gap.strand)?;
                            row.boolean(gap.overlap)?;
                            row.json(&json!(gap.overlapping_gene_ids))?;
                        } else {
                            for _ in 0..8 {
                                row.null()?;
                            }
                        }
                        Ok(())
                    })?;
                }
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "terminal_counts",
        &terminal_count_schema,
        Some(&SelectionSummary::complete(terminal_count_rows as u64)),
        |rows| {
            for entity in data.tails {
                for anchor in &entity.anchors {
                    for count in &anchor.counts {
                        rows.write_row_with(|row| {
                            row.string(&entity.id)?;
                            row.uint64(anchor.anchor as u64)?;
                            row.string(&data.collection.archives[count.sample].id)?;
                            row.string(&data.design.donor_names[count.donor])?;
                            row.string(&data.groups.names[count.group])?;
                            row.uint64(count.umis as u64)?;
                            row.uint64(count.cells as u64)?;
                            Ok(())
                        })?;
                    }
                }
            }
            Ok(())
        },
    )?;
    bundle.finish()
}

fn validate_group_scope(collection: &Collection, groups: &Groups) -> Result<u64> {
    if !groups.explicit {
        return Ok(0);
    }
    (0..collection.archives.len())
        .into_par_iter()
        .map(|sample| -> Result<u64> {
            let mut archive = open_source(collection, sample, None)?;
            let dictionary: FxHashSet<u32> = archive.cells()?.iter().copied().collect();
            if let Some(&unknown) = groups.by_sample[sample]
                .keys()
                .find(|packed| !dictionary.contains(packed))
            {
                bail!(
                    "cell-group map barcode {} is not in collection sample {}",
                    unpack_barcode(unknown),
                    collection.archives[sample].id
                );
            }
            Ok(archive.reader().bytes_read())
        })
        .try_reduce(
            || 0u64,
            |left, right| {
                left.checked_add(right)
                    .context("group-scope archive byte count overflow")
            },
        )
}

/// Report the stage profile on stderr so the human output carries the same accounting the JSON
/// summary records in `stage_seconds`/`stage_peak_rss_bytes`.
fn eprint_stage_profile(seconds: &StageSeconds, peak_rss: &StagePeakRssBytes) {
    let stages: [(&str, f64, u64); 7] = [
        (
            "load_catalogue",
            seconds.load_catalogue,
            peak_rss.load_catalogue,
        ),
        (
            "discover_candidates",
            seconds.discover_candidates,
            peak_rss.discover_candidates,
        ),
        (
            "route_candidates",
            seconds.route_candidates,
            peak_rss.route_candidates,
        ),
        (
            "exact_counting",
            seconds.exact_counting,
            peak_rss.exact_counting,
        ),
        (
            "terminal_tails",
            seconds.terminal_tails,
            peak_rss.terminal_tails,
        ),
        (
            "annotation_classification",
            seconds.annotation_classification,
            peak_rss.annotation_classification,
        ),
        ("output", seconds.output, peak_rss.output),
    ];
    for (name, elapsed, rss) in stages {
        eprintln!(
            "collection find-events stage {name}: {elapsed:.3}s, peak rss {:.3} GiB",
            rss as f64 / (1024.0 * 1024.0 * 1024.0)
        );
    }
}

fn human_output(data: &OutputData<'_>) {
    if data.tail_capability.requested {
        eprintln!(
            "terminal-tail capability: {} available, {} unavailable; only available archives are in the recurrence/zero denominator",
            data.tail_capability.available_archives, data.tail_capability.unavailable_archives
        );
        for (sample, &available) in data.tail_capability.archive_available.iter().enumerate() {
            eprintln!(
                "terminal-tail archive {}: {}",
                data.collection.archives[sample].id,
                if available {
                    "available (included)"
                } else {
                    "unavailable (excluded)"
                }
            );
        }
    }
    println!(
        "entity_id\tkind\tchrom\tstrand\texact_umi_classes\tcells\tsamples\tdonors\tgap_primary_class"
    );
    for &entity in data.retained {
        let candidate = &data.candidates[entity];
        let exact = &data.exact[entity];
        let primary = data.gap[entity]
            .as_ref()
            .and_then(|gap| gap.primary_class)
            .unwrap_or("NA");
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            candidate.key.id(&data.collection.chroms),
            candidate.key.kind.name(),
            data.collection.chroms[candidate.key.chrom as usize],
            match candidate.key.reported_strand_rev {
                Some(true) => "-",
                Some(false) => "+",
                None => ".",
            },
            exact.counts.metric(candidate.key.kind),
            exact.cells,
            exact.samples,
            exact.donors,
            primary,
        );
    }
    for entity in data.tails {
        let primary = entity
            .gap
            .as_ref()
            .and_then(|gap| gap.primary_class)
            .unwrap_or("NA");
        println!(
            "{}\tterminal_tail\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            entity.id,
            data.collection.chroms[entity.chrom as usize],
            if entity.strand_rev { "-" } else { "+" },
            entity.umis,
            entity.cells,
            entity.samples,
            entity.donors,
            primary,
        );
    }
}

pub(super) fn run(args: Args) -> Result<()> {
    let started = std::time::Instant::now();
    let mut clock = StageClock::new(started);
    let chain = open_collection_chain_with_locations(
        &args.collection,
        args.uniform_output.locations.as_deref(),
    )?;
    let collection = &chain.collection;
    if collection.archives.is_empty() {
        bail!("collection contains no source archives");
    }
    if args.min_samples == 0 || args.min_samples > collection.archives.len() {
        bail!(
            "--min-samples must be between 1 and the {} collection archives",
            collection.archives.len()
        );
    }
    if args.min_umis == 0 {
        bail!("--min-umi-classes must be at least 1");
    }
    if args.min_side_umis == 0 {
        bail!("--min-side-umi-classes must be at least 1");
    }
    if args.min_group_umis == 0 {
        bail!("--min-group-umi-classes must be at least 1");
    }
    let design = load_design(collection, args.design.as_deref())?;
    if args.min_donors == 0 || args.min_donors > design.donor_names.len() {
        bail!(
            "--min-donors must be between 1 and the {} biological donors",
            design.donor_names.len()
        );
    }
    let groups = load_groups(collection, args.groups.as_deref())?;
    let mut seen_required = FxHashSet::default();
    let required_groups: Vec<usize> = args
        .required_groups
        .iter()
        .map(|name| {
            if !seen_required.insert(name.as_str()) {
                bail!("duplicate --require-group {name}");
            }
            groups
                .names
                .iter()
                .position(|known| known == name)
                .with_context(|| {
                    format!("--require-group {name} is absent from the cell-group map")
                })
        })
        .collect::<Result<_>>()?;
    let kinds = selected_kinds(&args.kinds);
    validate_kind_strand(&kinds, args.solo_strand)?;

    let source_identity = validate_collection_sources(collection, args.verify_content)?;
    let identity_bytes = source_archive_bytes_read(&source_identity)?;
    let scope_bytes = validate_group_scope(collection, &groups)?;
    let catalogue = scan_chain_junctions(&chain)?;
    let annotation = if let Some(path) = args.annotation.as_deref() {
        let assembly = args
            .assembly
            .as_deref()
            .context("--annotation requires --assembly")?;
        let label = args
            .annotation_label
            .as_deref()
            .context("--annotation requires --annotation-label")?;
        let identity = anno::intent::AnnotationIdentity::new(assembly, label)
            .context("constructing annotation identity")?;
        let identity = match args.annotation_digest.as_deref() {
            Some(digest) => identity
                .with_digest(digest)
                .context("validating --annotation-digest")?,
            None => identity,
        };
        Some(load_annotation(path, collection, identity)?)
    } else {
        None
    };
    let annotation_comparison_budget = annotation
        .as_ref()
        .map(|_| AnnotationComparisonBudget::new(args.max_annotation_comparisons))
        .transpose()?;
    let annotation_context = annotation
        .as_ref()
        .zip(annotation_comparison_budget.as_ref());
    clock.lap(|s| &mut s.load_catalogue, |r| &mut r.load_catalogue);
    let discovery = discover_candidates(
        &catalogue,
        &kinds,
        &design,
        args.solo_strand,
        CandidateThresholds {
            min_support: args.min_support,
            min_samples: args.min_samples,
            min_donors: args.min_donors,
            max_candidates: args.max_candidates,
            max_candidates_considered: args.max_candidates_considered,
        },
    )?;
    let CandidateDiscovery {
        mut candidates,
        attempted: splice_candidates_attempted,
        distinct: splice_candidates_distinct,
    } = discovery;
    let annotation_candidate_exclusions = if let Some(annotation) = annotation.as_ref() {
        exclude_annotation_unmatched_candidates(collection, annotation, &mut candidates)?
    } else {
        AnnotationCandidateExclusions::default()
    };
    clock.lap(
        |s| &mut s.discover_candidates,
        |r| &mut r.discover_candidates,
    );
    let routes = route_candidates(collection, &catalogue, &candidates, args.max_routed_entries)?;
    // The catalogue's rows and their route postings are not needed once the exact plan exists.
    let catalogue_junctions = catalogue.len() as u64;
    drop(catalogue);
    release_free_heap();
    clock.lap(|s| &mut s.route_candidates, |r| &mut r.route_candidates);
    let exact_match_budget = ExactMatchBudget::new(args.max_exact_match_attempts)?;
    let (archives, exact) = reduce_exact(
        collection,
        &routes,
        &groups,
        &design,
        &candidates,
        &required_groups,
        &exact_match_budget,
    )?;
    clock.lap(|s| &mut s.exact_counting, |r| &mut r.exact_counting);
    let (tails, tail_capability) = if kinds.contains(&SearchKind::TerminalTail) {
        scan_terminal_tails(
            collection,
            &groups,
            &design,
            candidates.len(),
            &args,
            &required_groups,
            annotation_context,
        )?
    } else {
        (Vec::new(), TailCapabilitySummary::default())
    };
    if candidates
        .len()
        .checked_add(tail_capability.candidate_clusters)
        .context("reverse-search candidate count overflow")?
        > args.max_candidates
    {
        bail!(
            "reverse-search catalogue contains {} splice/junction candidates and {} terminal-tail clusters, exceeding --max-candidates {}; strengthen exact/catalogue thresholds or raise the explicit limit",
            candidates.len(),
            tail_capability.candidate_clusters,
            args.max_candidates,
        );
    }
    clock.lap(|s| &mut s.terminal_tails, |r| &mut r.terminal_tails);

    let solo_strand: anno::assign::SoloStrand = args.solo_strand.into();
    let gap: Vec<Option<GapClassification>> = candidates
        .iter()
        .enumerate()
        .map(|(entity, candidate)| -> Result<Option<GapClassification>> {
            let Some((annotation, budget)) = annotation_context else {
                return Ok(None);
            };
            let entity_exact = &exact.entities[entity];
            if !evidence_passes_predicates(
                candidate,
                entity_exact,
                exact.required_metrics(entity),
                &args,
            ) {
                return Ok(None);
            }
            classify_gap(annotation, candidate, entity_exact, solo_strand, budget).map(Some)
        })
        .collect::<Result<_>>()?;
    let retained: Vec<usize> = candidates
        .iter()
        .enumerate()
        .filter_map(|(entity, candidate)| {
            let absent = GapClassification::default();
            keep_entity(
                candidate,
                &exact.entities[entity],
                exact.required_metrics(entity),
                gap[entity].as_ref().unwrap_or(&absent),
                &args,
            )
            .then_some(entity)
        })
        .collect();
    clock.lap(
        |s| &mut s.annotation_classification,
        |r| &mut r.annotation_classification,
    );

    let unique_chunks: usize = archives.iter().map(|archive| archive.unique_chunks).sum();
    let independent_chunks: usize = archives
        .iter()
        .map(|archive| archive.independent_chunk_decodes)
        .sum();
    let planned_bytes: u64 = archives.iter().map(|archive| archive.planned_bytes).sum();
    let exact_bytes: u64 = archives.iter().map(|archive| archive.actual_bytes).sum();
    let actual_archive_bytes = identity_bytes
        .checked_add(scope_bytes)
        .and_then(|value| value.checked_add(exact_bytes))
        .and_then(|value| value.checked_add(tail_capability.actual_bytes))
        .context("reverse-search archive byte count overflow")?;
    let sidecar_bytes = collection_sidecar_bytes_read(&chain)?;
    let requested_kinds: Vec<&str> = kinds.iter().map(|kind| kind.name()).collect();
    let mut annotation_unmatched_chroms = annotation_candidate_exclusions.unmatched_chroms.clone();
    annotation_unmatched_chroms.extend(tail_capability.annotation_unmatched_chroms.iter().copied());
    let annotation_unmatched_evidence_contigs: Vec<String> = annotation_unmatched_chroms
        .into_iter()
        .map(|chrom| collection.chroms[chrom as usize].clone())
        .collect();
    clock.lap(|s| &mut s.output, |r| &mut r.output);
    let summary = SearchSummary {
        coordinates:
            "0-based junction boundaries and strand-aware terminal cleavage anchors; entity intervals are half-open",
        candidate_source:
            "authenticated collection junction route superset (unique chains plus direct multimapper anchors) followed by unique-chain-only exact filtering, plus rooted terminal-tail capability/index routes",
        exact_counting:
            "distinct archive classes (corrected cell barcode, retained raw UMI value) supported by unique-read chains; stored 1MM edges are not collapsed",
        missing_count_rows:
            "logical zero over retained entity x capable sample x requested cell-group dimensions; consult capabilities",
        event_usage:
            "alternative-site informative_umi_classes = side_a_only_umi_classes + side_b_only_umi_classes; cassette informative_umi_classes = include_only_umi_classes + exclude_only_umi_classes; both_umi_classes is reported but excluded",
        alternative_site_side_semantics:
            "side_a is the lexicographically lower (donor,acceptor) junction and side_b the higher; these are neutral coordinate labels, not exon inclusion or PSI",
        cassette_component_semantics:
            "include-only is the distinct-class union witnessing either inclusion flank and not the skip; marginal component counts must each pass, possibly in disjoint classes; no same-class full inclusion path is claimed",
        donor_semantics: if design.source.is_some() {
            "explicit sample-to-donor design"
        } else {
            "each collection sample is one donor"
        },
        cell_scope: if groups.explicit {
            "only sample/barcode rows listed in --groups"
        } else {
            "all archive cells in group bulk"
        },
        evidence_placement_policy: "unique_chain_representatives_only",
        multimapper_placements_included: false,
        multimapper_alternatives_available_to_search: false,
        requested_kinds,
        catalogue_junctions,
        splice_candidate_definitions_attempted: splice_candidates_attempted as u64,
        splice_candidate_definitions_distinct: splice_candidates_distinct as u64,
        candidate_entities: (candidates.len() + tail_capability.candidate_clusters) as u64,
        retained_entities: (retained.len() + tails.len()) as u64,
        archives_total: collection.archives.len() as u64,
        archives_opened: if tail_capability.requested {
            collection.archives.len() as u64
        } else {
            archives.len() as u64
        },
        donors_total: design.donor_names.len() as u64,
        groups_total: groups.names.len() as u64,
        min_support: args.min_support,
        min_samples: args.min_samples as u64,
        min_donors: args.min_donors as u64,
        min_umi_classes: args.min_umis as u64,
        min_side_umi_classes: args.min_side_umis as u64,
        min_group_umi_classes: args.min_group_umis as u64,
        max_candidates: args.max_candidates as u64,
        max_candidates_considered: args.max_candidates_considered as u64,
        max_routed_entries: args.max_routed_entries as u64,
        routed_target_associations: routes.target_associations as u64,
        routed_chunk_postings: routes.chunk_postings as u64,
        max_exact_match_attempts: args.max_exact_match_attempts,
        exact_match_attempts: exact_match_budget.attempted(),
        max_annotation_comparisons: args.max_annotation_comparisons,
        annotation_comparisons: annotation_comparison_budget
            .as_ref()
            .map_or(0, AnnotationComparisonBudget::attempted),
        terminal_cluster_bp: args.terminal_cluster_bp as u64,
        max_terminal_events: args.max_terminal_events,
        required_groups: &args.required_groups,
        annotation_gap_flags_nonexclusive: true,
        annotation_gap_primary_precedence: ["strand", "boundary", "overlap", "missing_junction"],
        annotation_identity: annotation.as_ref().map(|value| &value.identity),
        annotation_collection_compatibility: annotation
            .as_ref()
            .map(|_| "caller_declared_unverified"),
        annotation_exact_contig_names_matched: annotation.as_ref().map(|value| {
            value.matched_chroms.iter().filter(|matched| **matched).count() as u64
        }),
        annotation_collection_contigs: annotation
            .as_ref()
            .map(|_| collection.chroms.len() as u64),
        annotation_unmatched_evidence_policy: annotation
            .as_ref()
            .map(|_| "omit_unmatched_contig_evidence; never classify it as novel or zero"),
        annotation_excluded_splice_candidates: annotation_candidate_exclusions.candidates as u64,
        annotation_excluded_terminal_routes: tail_capability.annotation_excluded_routes as u64,
        annotation_excluded_terminal_events: tail_capability.annotation_excluded_events,
        annotation_unmatched_evidence_contigs,
        collection_genome_algo: collection.genome_algo.as_deref(),
        collection_genome_digest: collection.genome_digest.as_deref(),
        unique_chunks_decoded: unique_chunks as u64,
        independent_chunk_decodes: independent_chunks as u64,
        planned_compressed_bytes: planned_bytes,
        actual_archive_bytes_read: actual_archive_bytes,
        source_archive_identity_bytes_read: identity_bytes,
        collection_sidecar_bytes_read: sidecar_bytes,
        total_seconds: started.elapsed().as_secs_f64(),
        stage_seconds: clock.seconds,
        stage_peak_rss_bytes: clock.peak_rss,
        terminal_tail_available_archives: tail_capability.available_archives as u64,
        terminal_tail_unavailable_archives: tail_capability.unavailable_archives as u64,
        terminal_tail_available_donors: tail_capability.available_donors as u64,
        terminal_tail_declared_selected_molecules: tail_capability.declared_selected_molecules,
        terminal_tail_declared_events: tail_capability.declared_events,
        terminal_tail_routed_chunks: tail_capability.routed_chunks as u64,
    };
    let data = OutputData {
        collection,
        candidates: &candidates,
        exact: &exact.entities,
        gap: &gap,
        retained: &retained,
        archives: &archives,
        design: &design,
        groups: &groups,
        tails: &tails,
        tail_capability: &tail_capability,
    };
    if let Some(format) = args.uniform_output.format {
        let mut parameters = BTreeMap::new();
        parameters.insert("kinds".into(), json!(summary.requested_kinds));
        parameters.insert("min_support".into(), json!(args.min_support));
        parameters.insert("min_samples".into(), json!(args.min_samples));
        parameters.insert("min_donors".into(), json!(args.min_donors));
        parameters.insert("min_umi_classes".into(), json!(args.min_umis));
        parameters.insert("min_side_umi_classes".into(), json!(args.min_side_umis));
        parameters.insert(
            "terminal_cluster_bp".into(),
            json!(args.terminal_cluster_bp),
        );
        parameters.insert(
            "max_terminal_events".into(),
            json!(args.max_terminal_events),
        );
        parameters.insert("required_groups".into(), json!(&args.required_groups));
        parameters.insert("min_group_umi_classes".into(), json!(args.min_group_umis));
        parameters.insert("max_candidates".into(), json!(args.max_candidates));
        parameters.insert(
            "max_candidates_considered".into(),
            json!(args.max_candidates_considered),
        );
        parameters.insert("max_routed_entries".into(), json!(args.max_routed_entries));
        parameters.insert(
            "routed_target_associations".into(),
            json!(routes.target_associations),
        );
        parameters.insert("routed_chunk_postings".into(), json!(routes.chunk_postings));
        parameters.insert(
            "max_exact_match_attempts".into(),
            json!(args.max_exact_match_attempts),
        );
        parameters.insert(
            "exact_match_attempts".into(),
            json!(exact_match_budget.attempted()),
        );
        parameters.insert(
            "max_annotation_comparisons".into(),
            json!(args.max_annotation_comparisons),
        );
        parameters.insert(
            "annotation_comparisons".into(),
            json!(annotation_comparison_budget
                .as_ref()
                .map_or(0, AnnotationComparisonBudget::attempted)),
        );
        parameters.insert(
            "evidence_placement_policy".into(),
            json!("unique_chain_representatives_only"),
        );
        parameters.insert("multimapper_placements_included".into(), json!(false));
        parameters.insert(
            "multimapper_alternatives_available_to_search".into(),
            json!(false),
        );
        parameters.insert("verify_content".into(), json!(args.verify_content));
        parameters.insert("novel_only".into(), json!(args.novel_only));
        parameters.insert(
            "solo_strand".into(),
            json!(format!("{:?}", args.solo_strand).to_ascii_lowercase()),
        );
        if let Some(path) = design.source.as_deref() {
            parameters.insert(
                "design_path".into(),
                json!(uniform_path(path, "donor design")?),
            );
            parameters.insert("design_content_blake3".into(), json!(design.content_blake3));
        }
        if let Some(path) = groups.source.as_deref() {
            parameters.insert(
                "groups_path".into(),
                json!(uniform_path(path, "cell-group map")?),
            );
            parameters.insert("groups_content_blake3".into(), json!(groups.content_blake3));
        }
        if let (Some(path), Some(annotation)) = (args.annotation.as_deref(), annotation.as_ref()) {
            parameters.insert(
                "annotation_path".into(),
                json!(uniform_path(path, "annotation")?),
            );
            parameters.insert(
                "annotation_content_blake3".into(),
                json!(&annotation.content_blake3),
            );
            parameters.insert("annotation_identity".into(), json!(&annotation.identity));
            parameters.insert(
                "annotation_expected_digest".into(),
                json!(&args.annotation_digest),
            );
            parameters.insert(
                "annotation_collection_compatibility".into(),
                json!({
                    "status": "caller_declared_unverified",
                    "annotation_assembly": annotation.identity.assembly,
                    "collection_genome_algo": collection.genome_algo,
                    "collection_genome_digest": collection.genome_digest,
                    "exact_contig_names_matched": annotation
                        .matched_chroms
                        .iter()
                        .filter(|matched| **matched)
                        .count(),
                    "collection_contigs": collection.chroms.len(),
                    "unmatched_evidence_policy": "omit_unmatched_contig_evidence; never classify it as novel or zero",
                    "excluded_splice_candidates": summary.annotation_excluded_splice_candidates,
                    "excluded_terminal_routes": summary.annotation_excluded_terminal_routes,
                    "excluded_terminal_events": summary.annotation_excluded_terminal_events,
                    "unmatched_evidence_contigs": &summary.annotation_unmatched_evidence_contigs,
                    "note": "exact contig-name overlap is required for classified evidence but does not prove reference-sequence identity"
                }),
            );
        }
        parameters.insert(
            "archive_access".into(),
            json!("collection catalogue/capability prefilter followed by exact routed union chunk decoding"),
        );
        let context = collection_uniform_context(&args.collection, &chain, parameters)?;
        write_uniform_collection_result(&args.uniform_output, |writer| {
            stream_result(&mut *writer, format, &context, &summary, &data).map(|_| ())
        })?;
    } else {
        human_output(&data);
        eprint_stage_profile(&summary.stage_seconds, &summary.stage_peak_rss_bytes);
        if !summary.annotation_unmatched_evidence_contigs.is_empty() {
            eprintln!(
                "collection find-events: omitted {} splice candidates and {} terminal routes ({} declared events) on annotation-unmatched contigs: {}",
                summary.annotation_excluded_splice_candidates,
                summary.annotation_excluded_terminal_routes,
                summary.annotation_excluded_terminal_events,
                summary.annotation_unmatched_evidence_contigs.join(", "),
            );
        }
        eprintln!(
            "collection find-events: {} retained / {} candidates; {}/{} archives, {unique_chunks}/{independent_chunks} chunks ({:.3}s)",
            retained.len() + tails.len(),
            candidates.len() + tail_capability.candidate_clusters,
            if tail_capability.requested { collection.archives.len() } else { archives.len() },
            collection.archives.len(),
            started.elapsed().as_secs_f64(),
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {

    /// One routed coordinate whose targets occupy the whole arena, as the exact reducer sees it.
    fn one_chunk_target(
        donor: u32,
        acceptor: u32,
        targets: Vec<Target>,
    ) -> (ChunkTargets, Vec<Target>) {
        let range = (0u32, targets.len() as u32);
        ([((donor, acceptor), range)].into_iter().collect(), targets)
    }

    fn target(entity: u32, component_bit: u8, side: u8, strand_rev: Option<bool>) -> Target {
        Target {
            entity,
            component_bit,
            side,
            strand: Target::strand_code(strand_rev),
        }
    }

    /// Every packed hit one molecule contributes, in entity order.
    fn packed_hits(
        molecule: &MolRec,
        shapes: &[evidence_io::archive::Shape],
        wanted: &ChunkTargets,
        arena: &[Target],
        budget: &ExactMatchBudget,
    ) -> Result<Vec<u64>> {
        let mut hits = Vec::new();
        let mut scratch = Vec::new();
        let matcher = ChunkMatch {
            shapes,
            wanted,
            arena,
            budget,
        };
        molecule_packed_hits(molecule, matcher, &mut scratch, &mut hits)?;
        Ok(hits)
    }

    #[test]
    fn stage_clock_attributes_each_lap_and_records_a_peak() {
        let mut clock = StageClock::new(std::time::Instant::now());
        clock.lap(|s| &mut s.load_catalogue, |r| &mut r.load_catalogue);
        std::thread::sleep(std::time::Duration::from_millis(20));
        clock.lap(|s| &mut s.exact_counting, |r| &mut r.exact_counting);
        assert!(clock.seconds.exact_counting >= 0.02);
        assert!(clock.seconds.load_catalogue < clock.seconds.exact_counting);
        assert_eq!(clock.seconds.output, 0.0);
        // Linux publishes VmHWM; elsewhere the profile degrades to zero rather than failing.
        if std::path::Path::new("/proc/self/status").exists() {
            assert!(clock.peak_rss.load_catalogue > 0);
            assert!(clock.peak_rss.exact_counting >= clock.peak_rss.load_catalogue);
        }
    }
    use super::*;

    fn archive(id: &str) -> ArchiveEntry {
        ArchiveEntry {
            id: id.to_owned(),
            path: PathBuf::from(format!("{id}.aie")),
            identity: FileIdentity {
                len: 1,
                modified_secs: 0,
                modified_nanos: 0,
                changed_secs: 0,
                changed_nanos: 0,
                dev: 0,
                inode: 0,
                archive_format_version: evidence_io::format::VERSION,
                native_scheme: ROOTED_DIRECTORY_SCHEME.to_owned(),
                native_digest: "0".repeat(64),
                encoded_sections_digest: Some("1".repeat(64)),
            },
            chunks: Vec::new(),
            shape_routes: None,
        }
    }

    fn junction(donor: u32, acceptor: u32) -> GlobalJunction {
        GlobalJunction {
            chrom: 0,
            donor,
            acceptor,
            support_upper_bound: 10,
            routes: vec![ArchiveRoute {
                archive: 0,
                supporting_children: 10,
                posts: vec![0],
            }],
        }
    }

    fn design() -> Design {
        Design {
            donor_of_sample: vec![0],
            donor_names: vec!["D1".to_owned()],
            source: None,
            content_blake3: None,
        }
    }

    fn tail_hit(
        sample: usize,
        donor: usize,
        group: usize,
        cell: u32,
        class: u32,
        anchor: u32,
    ) -> TailHit {
        TailHit {
            sample,
            donor,
            group,
            cell,
            class,
            chrom: 0,
            strand_rev: false,
            anchor,
            signal: evidence_io::terminal_tail::TerminalTailSignal {
                clip_len: 10,
                tail_bases: 9,
                terminal_run: 7,
            },
        }
    }

    fn annotation_index(transcripts: Vec<AnnotationTranscript>) -> AnnotationIndex {
        let intervals = vec![AnnotationIntervalIndex::new(&transcripts)];
        AnnotationIndex {
            by_chrom: vec![transcripts],
            intervals,
            matched_chroms: vec![true],
            identity: anno::intent::AnnotationIdentity::new("GRCh38", "fixture-v1").unwrap(),
            content_blake3: format!("blake3:{}", "0".repeat(64)),
        }
    }

    fn match_budget() -> ExactMatchBudget {
        ExactMatchBudget::new(1_000_000).unwrap()
    }

    fn annotation_budget() -> AnnotationComparisonBudget {
        AnnotationComparisonBudget::new(1_000_000).unwrap()
    }

    #[test]
    fn biological_alt_site_names_follow_exact_strand() {
        assert_eq!(
            ComponentSide::Include.output_name(SearchKind::AltAcceptor),
            "side_a"
        );
        assert_eq!(
            ComponentSide::Exclude.output_name(SearchKind::AltDonor),
            "side_b"
        );
        let rows = vec![junction(100, 200), junction(100, 300), junction(150, 300)];
        let candidates = discover_candidates(
            &rows,
            &[SearchKind::AltAcceptor, SearchKind::AltDonor]
                .into_iter()
                .collect(),
            &design(),
            crate::archivecmd::SoloStrandArg::Forward,
            CandidateThresholds {
                min_support: 1,
                min_samples: 1,
                min_donors: 1,
                max_candidates: 100,
                max_candidates_considered: 1_000,
            },
        )
        .unwrap()
        .candidates;
        assert!(candidates.iter().any(|candidate| {
            candidate.key.kind == SearchKind::AltAcceptor
                && candidate.key.strand_rev == Some(false)
                && candidate.key.components[0].coordinate.donor == 100
                && candidate.key.components[1].coordinate.donor == 100
        }));
        assert!(candidates.iter().any(|candidate| {
            candidate.key.kind == SearchKind::AltDonor
                && candidate.key.strand_rev == Some(true)
                && candidate.key.components[0].coordinate.donor == 100
                && candidate.key.components[1].coordinate.donor == 100
        }));
        assert!(candidates.iter().any(|candidate| {
            candidate.key.kind == SearchKind::AltDonor
                && candidate.key.strand_rev == Some(false)
                && candidate.key.components[0].coordinate.acceptor == 300
                && candidate.key.components[1].coordinate.acceptor == 300
        }));
        assert!(candidates.iter().any(|candidate| {
            candidate.key.kind == SearchKind::AltAcceptor
                && candidate.key.strand_rev == Some(true)
                && candidate.key.components[0].coordinate.acceptor == 300
                && candidate.key.components[1].coordinate.acceptor == 300
        }));
    }

    #[test]
    fn catalogue_recurrence_unions_mutually_alternative_sides() {
        let rows = vec![
            GlobalJunction {
                chrom: 0,
                donor: 100,
                acceptor: 200,
                support_upper_bound: 10,
                routes: vec![ArchiveRoute {
                    archive: 0,
                    supporting_children: 10,
                    posts: vec![0],
                }],
            },
            GlobalJunction {
                chrom: 0,
                donor: 100,
                acceptor: 300,
                support_upper_bound: 10,
                routes: vec![ArchiveRoute {
                    archive: 1,
                    supporting_children: 10,
                    posts: vec![0],
                }],
            },
        ];
        let design = Design {
            donor_of_sample: vec![0, 1],
            donor_names: vec!["D1".to_owned(), "D2".to_owned()],
            source: None,
            content_blake3: None,
        };
        let candidates = discover_candidates(
            &rows,
            &[SearchKind::AltAcceptor].into_iter().collect(),
            &design,
            crate::archivecmd::SoloStrandArg::Forward,
            CandidateThresholds {
                min_support: 1,
                min_samples: 2,
                min_donors: 2,
                max_candidates: 100,
                max_candidates_considered: 1_000,
            },
        )
        .unwrap()
        .candidates;
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].catalogue_samples, 2);
        assert_eq!(candidates[0].catalogue_donors, 2);
    }

    #[test]
    fn unstranded_junction_coalesces_six_plus_six_before_thresholding() {
        let candidates = discover_candidates(
            &[junction(100, 200)],
            &[SearchKind::Junction].into_iter().collect(),
            &design(),
            crate::archivecmd::SoloStrandArg::Unstranded,
            CandidateThresholds {
                min_support: 1,
                min_samples: 1,
                min_donors: 1,
                max_candidates: 100,
                max_candidates_considered: 1_000,
            },
        )
        .unwrap()
        .candidates;
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].key.strand_rev, None);
        assert_eq!(candidates[0].key.reported_strand_rev, None);
        assert_eq!(
            candidates[0].key.id(&["chr1".to_owned()]),
            "junction:chr1:.:100-200"
        );

        let shapes = vec![evidence_io::archive::Shape {
            blocks: vec![(0, 10), (110, 10)],
        }];
        let (wanted, arena) = one_chunk_target(100, 200, vec![target(0, 1, 1, None)]);
        let mut counts = MaskCounts::default();
        let mut strand_umis = [0u32; 2];
        for class in 0..12u32 {
            let strand_rev = class >= 6;
            let molecule = MolRec {
                cell: class,
                umi_class: class,
                chrom: 0,
                strand_rev,
                chains: smallvec::smallvec![crate::rows::MolChain {
                    weight: 1,
                    reps: smallvec::smallvec![(90, 0)],
                }],
                mms: smallvec::smallvec![],
            };
            let hits = packed_hits(&molecule, &shapes, &wanted, &arena, &match_budget()).unwrap();
            assert_eq!(hits.len(), 1);
            counts.add_mask(SearchKind::Junction, unpack_exact_hit(hits[0]).2);
            strand_umis[usize::from(strand_rev)] += 1;
        }
        assert_eq!(counts.support, 12);
        assert_eq!(strand_umis, [6, 6]);

        let exact = EntityExact {
            counts,
            samples: 1,
            donors: 1,
            cells: 12,
            component_umi_classes: [12, 0, 0],
            strand_umis,
            last_donor_plus_one: 1,
        };
        let args = Args {
            collection: PathBuf::from("test.aicollection"),
            kinds: vec![SearchKind::Junction],
            design: None,
            groups: None,
            required_groups: Vec::new(),
            min_group_umis: 1,
            min_donors: 1,
            min_samples: 1,
            min_umis: 10,
            min_side_umis: 1,
            min_support: 1,
            terminal_cluster_bp: 25,
            max_terminal_events: 10_000_000,
            annotation: None,
            assembly: None,
            annotation_label: None,
            annotation_digest: None,
            novel_only: false,
            solo_strand: crate::archivecmd::SoloStrandArg::Unstranded,
            max_candidates: 100,
            max_candidates_considered: 1_000,
            max_routed_entries: 10_000_000,
            max_exact_match_attempts: 10_000_000,
            max_annotation_comparisons: 10_000_000,
            verify_content: false,
            uniform_output: CollectionIoArgs::default(),
        };
        assert!(keep_entity(
            &candidates[0],
            &exact,
            &[],
            &GapClassification::default(),
            &args,
        ));
    }

    #[test]
    fn annotation_gap_flags_are_nonexclusive_with_stable_primary_class() {
        let transcript = AnnotationTranscript {
            gene_id: "G1".to_owned(),
            strand_rev: false,
            start: 50,
            end: 350,
            junctions: [(100, 200)].into_iter().collect(),
            donors: [100].into_iter().collect(),
            acceptors: [200].into_iter().collect(),
        };
        let annotation = annotation_index(vec![transcript]);
        let candidate = Candidate {
            key: EntityKey {
                kind: SearchKind::Junction,
                chrom: 0,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![component(
                    Coordinate {
                        chrom: 0,
                        donor: 100,
                        acceptor: 250,
                    },
                    ComponentSide::Support,
                )],
            },
            catalogue_support_upper_bound: 2,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let exact = EntityExact {
            strand_umis: [2, 0],
            ..EntityExact::default()
        };
        let gap = classify_gap(
            &annotation,
            &candidate,
            &exact,
            anno::assign::SoloStrand::Forward,
            &annotation_budget(),
        )
        .unwrap();
        assert!(gap.incompatible_with_every_transcript);
        assert!(!gap.missing_junction);
        assert!(gap.boundary);
        assert!(!gap.overlap);
        assert!(!gap.strand);
        assert_eq!(gap.primary_class, Some("boundary"));
        assert_eq!(gap.overlapping_gene_ids, vec!["G1"]);
    }

    #[test]
    fn gap_primary_overlap_and_missing_junction_are_structurally_distinct() {
        let endpoint_transcript = AnnotationTranscript {
            gene_id: "G-overlap".to_owned(),
            strand_rev: false,
            start: 50,
            end: 300,
            junctions: [(100, 150), (180, 250)].into_iter().collect(),
            donors: [100, 180].into_iter().collect(),
            acceptors: [150, 250].into_iter().collect(),
        };
        let annotation = annotation_index(vec![endpoint_transcript.clone()]);
        let junction_candidate = |donor, acceptor| Candidate {
            key: EntityKey {
                kind: SearchKind::Junction,
                chrom: 0,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![component(
                    Coordinate {
                        chrom: 0,
                        donor,
                        acceptor,
                    },
                    ComponentSide::Support,
                )],
            },
            catalogue_support_upper_bound: 2,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let exact = EntityExact {
            strand_umis: [1, 0],
            ..EntityExact::default()
        };

        let missing = classify_gap(
            &annotation,
            &junction_candidate(100, 250),
            &exact,
            anno::assign::SoloStrand::Forward,
            &annotation_budget(),
        )
        .unwrap();
        assert!(missing.missing_junction);
        assert!(!missing.overlap);
        assert!(!missing.boundary);
        assert_eq!(missing.primary_class, Some("missing_junction"));

        let split_path_annotation = annotation_index(vec![
            AnnotationTranscript {
                gene_id: "G-overlap".to_owned(),
                strand_rev: false,
                start: 50,
                end: 170,
                junctions: [(100, 150)].into_iter().collect(),
                donors: [100].into_iter().collect(),
                acceptors: [150].into_iter().collect(),
            },
            AnnotationTranscript {
                gene_id: "G-overlap".to_owned(),
                strand_rev: false,
                start: 170,
                end: 300,
                junctions: [(180, 250)].into_iter().collect(),
                donors: [180].into_iter().collect(),
                acceptors: [250].into_iter().collect(),
            },
            AnnotationTranscript {
                gene_id: "G-overlap".to_owned(),
                strand_rev: false,
                start: 50,
                end: 300,
                junctions: [(100, 250)].into_iter().collect(),
                donors: [100].into_iter().collect(),
                acceptors: [250].into_iter().collect(),
            },
        ]);
        let path_candidate = Candidate {
            key: EntityKey {
                kind: SearchKind::Cassette,
                chrom: 0,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 100,
                            acceptor: 150,
                        },
                        ComponentSide::Include,
                    ),
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 180,
                            acceptor: 250,
                        },
                        ComponentSide::Include,
                    ),
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 100,
                            acceptor: 250,
                        },
                        ComponentSide::Exclude,
                    ),
                ],
            },
            catalogue_support_upper_bound: 2,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let overlap = classify_gap(
            &split_path_annotation,
            &path_candidate,
            &exact,
            anno::assign::SoloStrand::Forward,
            &annotation_budget(),
        )
        .unwrap();
        assert!(overlap.overlap);
        assert!(!overlap.missing_junction);
        assert!(!overlap.boundary);
        assert_eq!(overlap.primary_class, Some("overlap"));
    }

    #[test]
    fn opposite_exact_junction_is_a_strand_gap() {
        let transcript = AnnotationTranscript {
            gene_id: "G2".to_owned(),
            strand_rev: true,
            start: 50,
            end: 250,
            junctions: [(100, 200)].into_iter().collect(),
            donors: [100].into_iter().collect(),
            acceptors: [200].into_iter().collect(),
        };
        let annotation = annotation_index(vec![transcript]);
        let candidate = Candidate {
            key: EntityKey {
                kind: SearchKind::Junction,
                chrom: 0,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![component(
                    Coordinate {
                        chrom: 0,
                        donor: 100,
                        acceptor: 200,
                    },
                    ComponentSide::Support,
                )],
            },
            catalogue_support_upper_bound: 2,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let exact = EntityExact {
            strand_umis: [1, 0],
            ..EntityExact::default()
        };
        let gap = classify_gap(
            &annotation,
            &candidate,
            &exact,
            anno::assign::SoloStrand::Forward,
            &annotation_budget(),
        )
        .unwrap();
        assert!(gap.strand);
        assert_eq!(gap.primary_class, Some("strand"));
    }

    #[test]
    fn alternative_sides_require_a_common_annotated_gene() {
        let transcript = |gene: &str, donor: u32, acceptor: u32| AnnotationTranscript {
            gene_id: gene.to_owned(),
            strand_rev: false,
            start: 50,
            end: 350,
            junctions: [(donor, acceptor)].into_iter().collect(),
            donors: [donor].into_iter().collect(),
            acceptors: [acceptor].into_iter().collect(),
        };
        let candidate = Candidate {
            key: EntityKey {
                kind: SearchKind::AltAcceptor,
                chrom: 0,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 100,
                            acceptor: 200,
                        },
                        ComponentSide::Include,
                    ),
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 100,
                            acceptor: 300,
                        },
                        ComponentSide::Exclude,
                    ),
                ],
            },
            catalogue_support_upper_bound: 2,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let exact = EntityExact {
            strand_umis: [2, 0],
            ..EntityExact::default()
        };

        let same_gene =
            annotation_index(vec![transcript("G1", 100, 200), transcript("G1", 100, 300)]);
        let compatible = classify_gap(
            &same_gene,
            &candidate,
            &exact,
            anno::assign::SoloStrand::Forward,
            &annotation_budget(),
        )
        .unwrap();
        assert!(!compatible.incompatible_with_every_transcript);
        assert_eq!(compatible.compatible_transcripts, 2);
        assert_eq!(compatible.overlapping_gene_ids, vec!["G1"]);

        let disjoint_genes =
            annotation_index(vec![transcript("G1", 100, 200), transcript("G2", 100, 300)]);
        let incompatible = classify_gap(
            &disjoint_genes,
            &candidate,
            &exact,
            anno::assign::SoloStrand::Forward,
            &annotation_budget(),
        )
        .unwrap();
        assert!(incompatible.incompatible_with_every_transcript);
        assert_eq!(incompatible.compatible_transcripts, 0);
        assert!(incompatible.overlap);
        assert_eq!(incompatible.primary_class, Some("overlap"));
    }

    #[test]
    fn alternative_event_requires_exact_support_on_both_sides() {
        let candidate = Candidate {
            key: EntityKey {
                kind: SearchKind::AltAcceptor,
                chrom: 0,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 100,
                            acceptor: 200,
                        },
                        ComponentSide::Include,
                    ),
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 100,
                            acceptor: 300,
                        },
                        ComponentSide::Exclude,
                    ),
                ],
            },
            catalogue_support_upper_bound: 10,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let mut exact = EntityExact {
            counts: MaskCounts {
                include_only: 10,
                exclude_only: 0,
                ..MaskCounts::default()
            },
            component_umi_classes: [10, 10, 0],
            samples: 1,
            donors: 1,
            ..EntityExact::default()
        };
        let args = Args {
            collection: PathBuf::from("test.aicollection"),
            kinds: Vec::new(),
            design: None,
            groups: None,
            required_groups: Vec::new(),
            min_group_umis: 1,
            min_donors: 1,
            min_samples: 1,
            min_umis: 1,
            min_side_umis: 1,
            min_support: 1,
            terminal_cluster_bp: 25,
            max_terminal_events: 10_000_000,
            annotation: None,
            assembly: None,
            annotation_label: None,
            annotation_digest: None,
            novel_only: false,
            solo_strand: crate::archivecmd::SoloStrandArg::Forward,
            max_candidates: 100,
            max_candidates_considered: 1_000,
            max_routed_entries: 10_000_000,
            max_exact_match_attempts: 10_000_000,
            max_annotation_comparisons: 10_000_000,
            verify_content: false,
            uniform_output: CollectionIoArgs::default(),
        };
        assert!(!keep_entity(
            &candidate,
            &exact,
            &[],
            &GapClassification::default(),
            &args,
        ));
        exact.counts.exclude_only = 1;
        assert!(keep_entity(
            &candidate,
            &exact,
            &[],
            &GapClassification::default(),
            &args,
        ));

        let cassette = Candidate {
            key: EntityKey {
                kind: SearchKind::Cassette,
                chrom: 0,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 100,
                            acceptor: 150,
                        },
                        ComponentSide::Include,
                    ),
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 180,
                            acceptor: 250,
                        },
                        ComponentSide::Include,
                    ),
                    component(
                        Coordinate {
                            chrom: 0,
                            donor: 100,
                            acceptor: 250,
                        },
                        ComponentSide::Exclude,
                    ),
                ],
            },
            catalogue_support_upper_bound: 10,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let mut cassette_exact = EntityExact {
            counts: MaskCounts {
                include_only: 10,
                exclude_only: 10,
                ..MaskCounts::default()
            },
            samples: 1,
            donors: 1,
            component_umi_classes: [10, 0, 10],
            ..EntityExact::default()
        };
        assert!(!keep_entity(
            &cassette,
            &cassette_exact,
            &[],
            &GapClassification::default(),
            &args,
        ));
        cassette_exact.component_umi_classes[1] = 1;
        assert!(keep_entity(
            &cassette,
            &cassette_exact,
            &[],
            &GapClassification::default(),
            &args,
        ));
    }

    #[test]
    fn terminal_cluster_retains_exact_anchors_and_deduplicates_classes() {
        let hits = vec![
            tail_hit(0, 0, 0, 1, 10, 100),
            tail_hit(0, 0, 0, 1, 10, 101),
            tail_hit(1, 1, 1, 2, 10, 101),
        ];
        let entity = build_tail_entity(&["chr1".to_owned()], 0, false, &hits, None).unwrap();
        assert_eq!((entity.start, entity.end, entity.summit), (100, 102, 101));
        assert_eq!((entity.umis, entity.samples, entity.donors), (2, 2, 2));
        assert_eq!(entity.group_umis.get(&0), Some(&1));
        assert_eq!(entity.group_umis.get(&1), Some(&1));
        assert_eq!(entity.anchors.len(), 2);
        assert_eq!(entity.anchors[0].umis, 1);
        assert_eq!(entity.anchors[1].umis, 2);
    }

    #[test]
    fn terminal_candidate_cap_fails_before_building_a_third_separated_cluster() {
        let separated_anchors = [100, 200, 300];
        let mut terminal_candidates = 0;
        for anchor in separated_anchors.into_iter().take(2) {
            let _ = anchor;
            reserve_terminal_candidate(0, &mut terminal_candidates, 2).unwrap();
        }
        assert_eq!(terminal_candidates, 2);
        let error = reserve_terminal_candidate(0, &mut terminal_candidates, 2).unwrap_err();
        assert!(error.to_string().contains("--max-candidates 2"));
        assert_eq!(terminal_candidates, 2);
    }

    #[test]
    fn terminal_novel_only_filters_exact_anchors_before_reclustering() {
        let transcript = AnnotationTranscript {
            gene_id: "G-terminal".to_owned(),
            strand_rev: false,
            start: 50,
            end: 100,
            junctions: FxHashSet::default(),
            donors: FxHashSet::default(),
            acceptors: FxHashSet::default(),
        };
        let annotation = annotation_index(vec![transcript]);
        let budget = annotation_budget();
        let mut hits = vec![tail_hit(0, 0, 0, 1, 10, 100), tail_hit(0, 0, 0, 2, 11, 110)];
        retain_incompatible_tail_anchors(&mut hits, &annotation, &budget).unwrap();
        assert_eq!(
            hits.iter().map(|hit| hit.anchor).collect::<Vec<_>>(),
            vec![110]
        );
        let entity = build_tail_entity(
            &["chr1".to_owned()],
            0,
            false,
            &hits,
            Some((&annotation, &budget)),
        )
        .unwrap();
        assert_eq!((entity.start, entity.end, entity.umis), (110, 111, 1));
        assert_eq!(entity.anchors.len(), 1);
        assert!(
            entity.anchors[0]
                .gap
                .as_ref()
                .unwrap()
                .incompatible_with_every_transcript
        );
    }

    #[test]
    fn candidate_attempt_cap_applies_even_when_recurrence_rejects_everything() {
        let rows: Vec<GlobalJunction> = (0..8)
            .map(|index| junction(100, 200 + index * 10))
            .collect();
        let error = discover_candidates(
            &rows,
            &[SearchKind::AltAcceptor].into_iter().collect(),
            &design(),
            crate::archivecmd::SoloStrandArg::Forward,
            CandidateThresholds {
                min_support: 1,
                min_samples: 2,
                min_donors: 1,
                max_candidates: 100,
                max_candidates_considered: 5,
            },
        )
        .err()
        .expect("the attempted-candidate cap must fail before unbounded enumeration");
        assert!(error.to_string().contains("--max-candidates-considered 5"));

        let mut cassette_rows = vec![junction(100, 500)];
        for acceptor in [150, 200, 250] {
            cassette_rows.push(junction(100, acceptor));
        }
        for donor in [300, 350, 400] {
            cassette_rows.push(junction(donor, 500));
        }
        // The catalogue is always unique and ascending by coordinate; planning relies on it.
        cassette_rows.sort_unstable_by_key(|row| (row.chrom, row.donor, row.acceptor));
        let cassette_error = discover_candidates(
            &cassette_rows,
            &[SearchKind::Cassette].into_iter().collect(),
            &design(),
            crate::archivecmd::SoloStrandArg::Forward,
            CandidateThresholds {
                min_support: 1,
                min_samples: 2,
                min_donors: 1,
                max_candidates: 100,
                max_candidates_considered: 5,
            },
        )
        .err()
        .expect("the attempted-candidate cap must bound cubic cassette enumeration");
        assert!(cassette_error
            .to_string()
            .contains("--max-candidates-considered 5"));
    }

    #[test]
    fn direct_primary_multimapper_tuple_is_not_exact_search_support() {
        let shapes = vec![evidence_io::archive::Shape {
            blocks: vec![(0, 10), (110, 10)],
        }];
        let (wanted, arena) = one_chunk_target(100, 200, vec![target(0, 1, 1, Some(false))]);
        let multimapper = MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: false,
            chains: smallvec::smallvec![],
            // The BAM-designated primary placement directly matches the target.
            mms: smallvec::smallvec![(90, 0, 0, 1)],
        };
        assert!(
            packed_hits(&multimapper, &shapes, &wanted, &arena, &match_budget())
                .unwrap()
                .is_empty()
        );

        let unique = MolRec {
            chains: smallvec::smallvec![crate::rows::MolChain {
                weight: 1,
                reps: smallvec::smallvec![(90, 0)],
            }],
            mms: smallvec::smallvec![],
            ..multimapper
        };
        assert_eq!(
            packed_hits(&unique, &shapes, &wanted, &arena, &match_budget())
                .unwrap()
                .len(),
            1
        );
    }

    #[test]
    fn exact_match_budget_fails_before_expanding_a_hot_target_list() {
        let shapes = vec![evidence_io::archive::Shape {
            blocks: vec![(0, 10), (110, 10)],
        }];
        let (wanted, arena) = one_chunk_target(
            100,
            200,
            vec![
                target(0, 1, 1, Some(false)),
                target(1, 1, 1, Some(false)),
                target(2, 1, 1, Some(false)),
            ],
        );
        let molecule = MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: false,
            chains: smallvec::smallvec![crate::rows::MolChain {
                weight: 1,
                reps: smallvec::smallvec![(90, 0)],
            }],
            mms: smallvec::smallvec![],
        };
        let budget = ExactMatchBudget::new(2).unwrap();
        let error = packed_hits(&molecule, &shapes, &wanted, &arena, &budget).unwrap_err();
        assert!(error.to_string().contains("--max-exact-match-attempts 2"));
        assert_eq!(budget.attempted(), 0);
    }

    #[test]
    fn alternative_only_multimapper_pattern_is_not_searchable_support() {
        let shapes = vec![evidence_io::archive::Shape {
            blocks: vec![(0, 10), (110, 10)],
        }];
        let (wanted, arena) = one_chunk_target(100, 200, vec![target(0, 1, 1, Some(false))]);
        let alternative_only = MolRec {
            cell: 0,
            umi_class: 0,
            chrom: 0,
            strand_rev: false,
            chains: smallvec::smallvec![],
            // Its direct anchor is elsewhere. Pattern id 1 is deliberately not expanded by the
            // unique-chain-only search, even if that dictionary entry contains the target.
            mms: smallvec::smallvec![(1_000, 0, 1, 1)],
        };
        assert!(
            packed_hits(&alternative_only, &shapes, &wanted, &arena, &match_budget())
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn terminal_gap_uses_exact_strand_aware_three_prime_boundary() {
        let transcript = AnnotationTranscript {
            gene_id: "G3".to_owned(),
            strand_rev: false,
            start: 50,
            end: 200,
            junctions: FxHashSet::default(),
            donors: FxHashSet::default(),
            acceptors: FxHashSet::default(),
        };
        let annotation = annotation_index(vec![transcript]);
        let budget = annotation_budget();
        let exact = build_tail_entity(
            &["chr1".to_owned()],
            0,
            false,
            &[tail_hit(0, 0, 0, 1, 10, 200)],
            Some((&annotation, &budget)),
        )
        .unwrap();
        let exact_gap = classify_terminal_gap(&annotation, &exact, &budget).unwrap();
        assert!(!exact_gap.incompatible_with_every_transcript);
        assert_eq!(exact_gap.compatible_transcripts, 1);

        let internal = build_tail_entity(
            &["chr1".to_owned()],
            0,
            false,
            &[tail_hit(0, 0, 0, 1, 11, 190)],
            Some((&annotation, &budget)),
        )
        .unwrap();
        let internal_gap = classify_terminal_gap(&annotation, &internal, &budget).unwrap();
        assert!(internal_gap.incompatible_with_every_transcript);
        assert!(internal_gap.boundary);
        assert!(internal_gap.overlap);
        assert_eq!(internal_gap.primary_class, Some("boundary"));
        assert_eq!(internal_gap.overlapping_gene_ids, vec!["G3"]);
    }

    #[test]
    fn terminal_search_rejects_nonforward_strand_reinterpretation() {
        let kinds = [SearchKind::TerminalTail].into_iter().collect();
        validate_kind_strand(&kinds, crate::archivecmd::SoloStrandArg::Forward).unwrap();
        let error =
            validate_kind_strand(&kinds, crate::archivecmd::SoloStrandArg::Reverse).unwrap_err();
        assert!(error.to_string().contains("fixed forward-cDNA"));
        assert!(
            validate_kind_strand(&kinds, crate::archivecmd::SoloStrandArg::Unstranded).is_err()
        );
    }

    #[test]
    fn design_defaults_each_collection_sample_to_its_own_donor() {
        let collection = Collection {
            base: None,
            genome_algo: None,
            genome_digest: None,
            chroms: vec!["chr1".to_owned()],
            chroms_digest: "0".repeat(64),
            archives: vec![archive("A"), archive("B")],
            junctions: Vec::new(),
            junction_count: 0,
            route_count: 0,
            posting_count: 0,
            shape_route_blocks: Vec::new(),
            encoded_shape_route_blocks: Vec::new(),
        };
        let design = load_design(&collection, None).unwrap();
        assert_eq!(design.donor_of_sample, vec![0, 1]);
        assert_eq!(design.donor_names, vec!["A", "B"]);
    }

    #[test]
    fn annotation_contig_guard_omits_unmatched_evidence_contigs() {
        let collection = Collection {
            base: None,
            genome_algo: None,
            genome_digest: None,
            chroms: vec!["chr1".to_owned(), "chr2".to_owned()],
            chroms_digest: "0".repeat(64),
            archives: vec![archive("A")],
            junctions: Vec::new(),
            junction_count: 0,
            route_count: 0,
            posting_count: 0,
            shape_route_blocks: Vec::new(),
            encoded_shape_route_blocks: Vec::new(),
        };
        let annotation = AnnotationIndex {
            by_chrom: vec![Vec::new(), Vec::new()],
            intervals: vec![
                AnnotationIntervalIndex::default(),
                AnnotationIntervalIndex::default(),
            ],
            matched_chroms: vec![true, false],
            identity: anno::intent::AnnotationIdentity::new("GRCh38", "fixture-v1").unwrap(),
            content_blake3: format!("blake3:{}", "0".repeat(64)),
        };
        assert!(annotation
            .chrom_is_matched(&collection, 0, "splice-event")
            .unwrap());
        assert!(!annotation
            .chrom_is_matched(&collection, 1, "splice-event")
            .unwrap());

        let candidate = |chrom| Candidate {
            key: EntityKey {
                kind: SearchKind::Junction,
                chrom,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![component(
                    Coordinate {
                        chrom,
                        donor: 100,
                        acceptor: 200,
                    },
                    ComponentSide::Support,
                )],
            },
            catalogue_support_upper_bound: 10,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let mut candidates = vec![candidate(0), candidate(1)];
        let excluded =
            exclude_annotation_unmatched_candidates(&collection, &annotation, &mut candidates)
                .unwrap();
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].key.chrom, 0);
        assert_eq!(excluded.candidates, 1);
        assert_eq!(excluded.unmatched_chroms, [1].into_iter().collect());
    }

    #[test]
    fn annotation_interval_index_localizes_work_and_budget_fails_closed() {
        let transcripts: Vec<AnnotationTranscript> = (0..1_000)
            .map(|index| AnnotationTranscript {
                gene_id: format!("G{index}"),
                strand_rev: false,
                start: index * 100,
                end: index * 100 + 50,
                junctions: FxHashSet::default(),
                donors: FxHashSet::default(),
                acceptors: FxHashSet::default(),
            })
            .collect();
        let index = AnnotationIntervalIndex::new(&transcripts);
        let budget = AnnotationComparisonBudget::new(10).unwrap();
        let overlap = index
            .overlapping_interval(&transcripts, 50_010, 50_020, &budget)
            .unwrap();
        assert_eq!(overlap, vec![500]);
        assert!(budget.attempted() <= 3);

        let too_small = AnnotationComparisonBudget::new(1).unwrap();
        let error = index
            .overlapping_interval(&transcripts, 50_010, 50_020, &too_small)
            .unwrap_err();
        assert!(error.to_string().contains("--max-annotation-comparisons 1"));
    }

    #[test]
    fn exact_route_plan_is_bounded_before_materialization() {
        let collection = Collection {
            base: None,
            genome_algo: None,
            genome_digest: None,
            chroms: vec!["chr1".to_owned()],
            chroms_digest: "0".repeat(64),
            archives: vec![archive("A")],
            junctions: Vec::new(),
            junction_count: 0,
            route_count: 0,
            posting_count: 0,
            shape_route_blocks: Vec::new(),
            encoded_shape_route_blocks: Vec::new(),
        };
        let row = GlobalJunction {
            chrom: 0,
            donor: 100,
            acceptor: 200,
            support_upper_bound: 10,
            routes: vec![ArchiveRoute {
                archive: 0,
                supporting_children: 10,
                posts: vec![0, 1],
            }],
        };
        let candidate = Candidate {
            key: EntityKey {
                kind: SearchKind::Junction,
                chrom: 0,
                strand_rev: Some(false),
                reported_strand_rev: Some(false),
                components: smallvec::smallvec![component(
                    Coordinate {
                        chrom: 0,
                        donor: 100,
                        acceptor: 200,
                    },
                    ComponentSide::Support,
                )],
            },
            catalogue_support_upper_bound: 10,
            catalogue_samples: 1,
            catalogue_donors: 1,
        };
        let error = route_candidates(
            &collection,
            std::slice::from_ref(&row),
            std::slice::from_ref(&candidate),
            2,
        )
        .err()
        .expect("one association plus two postings must exceed a bound of two");
        assert!(error.to_string().contains("--max-routed-entries 2"));

        let plan = route_candidates(&collection, &[row], &[candidate], 3).unwrap();
        assert_eq!(plan.target_associations, 1);
        assert_eq!(plan.chunk_postings, 2);
    }

    #[test]
    fn group_scope_barcodes_are_exactly_sixteen_bases() {
        assert!(crate::querycmd::pack_cell_barcode_16("A").is_err());
        assert!(crate::querycmd::pack_cell_barcode_16("AAAAAAAAAAAAAAA").is_err());
        assert!(crate::querycmd::pack_cell_barcode_16("AAAAAAAAAAAAAAAA").is_ok());
        assert!(crate::querycmd::pack_cell_barcode_16("AAAAAAAAAAAAAAAN").is_err());
    }

    /// Deterministic xorshift, so the reducer equivalence test covers many shapes without a
    /// random-number dependency.
    fn next_random(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn packed_exact_reduction_matches_a_naive_hash_map() {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut packed = Vec::new();
        let mut naive: FxHashMap<(u32, u32), u8> = FxHashMap::default();
        for _ in 0..50_000 {
            let value = next_random(&mut state);
            let entity = (value % 97) as u32;
            let class = ((value >> 17) % 211) as u32;
            let payload = ((value >> 40) & EXACT_HIT_PAYLOAD_MASK) as u8;
            packed.push(pack_exact_hit(entity, class, payload));
            *naive.entry((entity, class)).or_insert(0) |= payload;
        }
        packed.sort_unstable();
        reduce_sorted_exact_hits(&mut packed);

        let observed: Vec<(u32, u32, u8)> = packed
            .iter()
            .map(|&hit| {
                let (entity, class, side, components, strands) = unpack_exact_hit(hit);
                (
                    entity,
                    class,
                    side | (components << EXACT_HIT_COMPONENT_SHIFT)
                        | (strands << EXACT_HIT_STRAND_SHIFT),
                )
            })
            .collect();
        let mut expected: Vec<(u32, u32, u8)> = naive
            .into_iter()
            .map(|((entity, class), payload)| (entity, class, payload))
            .collect();
        expected.sort_unstable();
        assert_eq!(observed, expected);
        assert!(observed.windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn packed_hit_round_trips_every_field_at_its_bounds() {
        let entity = (MAX_PACKED_ENTITIES - 1) as u32;
        let payload = EXACT_HIT_PAYLOAD_MASK as u8;
        let hit = pack_exact_hit(entity, u32::MAX, payload);
        assert_eq!(unpack_exact_hit(hit), (entity, u32::MAX, 3, 7, 3));
        assert_eq!(unpack_exact_hit(pack_exact_hit(0, 0, 0)), (0, 0, 0, 0, 0));
    }

    #[test]
    fn entity_hit_slices_cover_every_hit_without_splitting_an_entity() {
        let mut hits: Vec<u64> = Vec::new();
        for entity in 0..40u32 {
            for class in 0..(entity % 7 + 1) {
                hits.push(pack_exact_hit(entity, class, 1));
            }
        }
        for parts in [1usize, 2, 3, 8, 64] {
            let slices = entity_hit_slices(&hits, parts);
            assert_eq!(
                slices.iter().map(|slice| slice.len()).sum::<usize>(),
                hits.len()
            );
            let mut previous: Option<u32> = None;
            for slice in &slices {
                assert!(!slice.is_empty());
                let first = unpack_exact_hit(slice[0]).0;
                assert!(previous.is_none_or(|last| last < first));
                previous = Some(unpack_exact_hit(slice[slice.len() - 1]).0);
            }
        }
        assert!(entity_hit_slices(&[], 4).is_empty());
    }

    #[test]
    fn compact_recurrence_counts_samples_once_and_donors_in_donor_order() {
        let archive_exact = |metric: usize, cells: u32| ArchiveEntityExact {
            entity: 0,
            counts: MaskCounts {
                support: metric,
                ..MaskCounts::default()
            },
            component_umi_classes: [1, 0, 0],
            cells,
            strands: [metric as u32, 0],
        };
        let mut exact = EntityExact::default();
        // Two archives of donor 0, then one of donor 1: three samples, two donors.
        merge_archive_entity(&mut exact, &archive_exact(3, 2), SearchKind::Junction, 0).unwrap();
        merge_archive_entity(&mut exact, &archive_exact(4, 5), SearchKind::Junction, 0).unwrap();
        merge_archive_entity(&mut exact, &archive_exact(1, 1), SearchKind::Junction, 1).unwrap();
        assert_eq!(exact.samples, 3);
        assert_eq!(exact.donors, 2);
        assert_eq!(exact.cells, 8);
        assert_eq!(exact.counts.support, 8);
        assert_eq!(exact.strand_umis, [8, 0]);
        assert_eq!(exact.component_umi_classes, [3, 0, 0]);

        // An archive with component evidence but no informative metric is not recurrence.
        merge_archive_entity(&mut exact, &archive_exact(0, 0), SearchKind::Junction, 2).unwrap();
        assert_eq!(exact.samples, 3);
        assert_eq!(exact.donors, 2);
        assert_eq!(exact.component_umi_classes, [4, 0, 0]);
    }

    #[test]
    fn exact_group_metrics_are_dense_per_required_group_and_cells_are_informative() {
        // Per-group counts are kept only for the requested groups, as one dense row of counters
        // per entity, rather than a map from group to counts on every candidate.
        let mut aggregate = ExactAggregate::new(100, 2);
        assert!(aggregate
            .entities
            .iter()
            .all(|exact| *exact == EntityExact::default()));
        assert_eq!(aggregate.required_group_metrics.len(), 200);
        aggregate.required_group_metrics[7 * 2 + 1] = 5;
        assert_eq!(aggregate.required_metrics(7), &[0, 5]);
        assert_eq!(aggregate.required_metrics(8), &[0, 0]);
        assert_eq!(
            aggregate.required_group_metrics.iter().sum::<u32>(),
            5,
            "no other entity gained a group counter"
        );
        assert_eq!(ExactAggregate::new(4, 0).required_metrics(3), &[] as &[u32]);

        assert!(mask_is_informative(SearchKind::Cassette, 1));
        assert!(mask_is_informative(SearchKind::Cassette, 2));
        assert!(!mask_is_informative(SearchKind::Cassette, 3));
        assert!(mask_is_informative(SearchKind::Junction, 1));
    }
}
