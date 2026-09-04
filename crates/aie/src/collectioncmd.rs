//! Optional collection index and exact federated query planner.
//!
//! `.aicollection` is a derived sidecar. Its base index decodes only archive metadata, the chunk
//! directory, and junction catalogue/postings; optional local-shape routes additionally decode
//! the source shape dictionary, never molecule chunks. A rooted v2 archive supplies both content
//! identities from its authenticated directory, whereas a legacy v1 source requires one complete
//! identity scan. Query answers remain grounded in the source archives: the sidecar chooses local
//! shapes, archives, and chunks, then the ordinary archive decoder recomputes exact molecule-class
//! counts.

use crate::archivecmd::{
    decode_chunk, read_chunk_index, remove_staging_if_owned, ChunkInfo, LazyArchive,
};
use crate::querycmd::{junction_counts_routed_with_shape_route, parse_locus, region_selects_chunk};
use crate::shaperoute::{
    self, EncodedRouteBlock, EncodedShapeRoutes, RouteBlock, ShapeRouteBinding, SpanRoute,
};
use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use evidence_io::archive::put_varint;
use evidence_io::format::Cursor;
use gravlax_output::{
    install_open_file_no_clobber, publish_file_no_clobber, reported_output_path, DataType,
    Durability, Field, OutputError, OutputFormat, Producer, Provenance, ResultContext,
    RowSemantics, SelectionSummary, StreamingBundleWriter, TableSchema, TableSemantics,
};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use serde_json::json;
use std::collections::BTreeMap;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const MAGIC: &[u8; 8] = b"GRVLXCOL";
const VERSION: u32 = 4;
const PREVIOUS_VERSION: u32 = 3;
const LEGACY_VERSION: u32 = 2;
const SCHEMA: &str = "gravlax.collection.v4";
const LEGACY_FULL_FILE_SCHEME: &str = "full-file-blake3-v1";
const ROOTED_DIRECTORY_SCHEME: &str = "aie-directory-root-v2";
const MAX_ITEMS: usize = 10_000_000;
const MAX_SECTIONS: usize = 1_000_000;
const MAX_CHROMS: usize = 65_536;
const MAX_ARCHIVES_PER_LAYER: usize = 65_536;
const MAX_CHUNKS_PER_ARCHIVE: usize = 2_000_000;
const MAX_JUNCTION_ROWS_PER_SEGMENT: usize = 250_000;
const MAX_SECTION_RAW_BYTES: u64 = 256 * 1024 * 1024;
const MAX_COMPRESSION_RATIO: u64 = 4_096;
const JUNCTION_BIN_BP: u32 = 16_000_000;

#[derive(Parser)]
pub struct Args {
    #[command(subcommand)]
    what: What,
}

/// The opt-in collection result contract. There is intentionally no default: omitting
/// `--format` preserves each collection command's historical stdout byte contract.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum CollectionOutputFormat {
    Text,
    Tsv,
    Json,
}

impl From<CollectionOutputFormat> for OutputFormat {
    fn from(value: CollectionOutputFormat) -> Self {
        match value {
            CollectionOutputFormat::Text => Self::Text,
            CollectionOutputFormat::Tsv => Self::Tsv,
            CollectionOutputFormat::Json => Self::Json,
        }
    }
}

#[derive(clap::Args, Clone, Debug, Default)]
struct CollectionOutputArgs {
    /// Use the versioned uniform result contract instead of the legacy presentation.
    #[arg(long, value_enum)]
    format: Option<CollectionOutputFormat>,
    /// Publish the uniform result atomically without replacing an existing file; requires --format.
    #[arg(short = 'o', long, requires = "format")]
    output: Option<PathBuf>,
}

#[derive(clap::Subcommand)]
enum What {
    /// Build a deterministic sidecar from decoded indexes and an exact source identity.
    Build {
        /// Named archive as ID=PATH. IDs are sorted, so input order does not affect the file.
        #[arg(long = "sample", required_unless_present = "base")]
        samples: Vec<String>,
        /// Expected native BLAKE3 as ID=HEX (full-file digest for v1, directory root for v2).
        /// The builder always derives and verifies the identity from the source itself.
        #[arg(long = "source-digest", value_name = "ID=BLAKE3")]
        source_digests: Vec<String>,
        /// Extend an existing sidecar without rescanning its source archive indexes.
        #[arg(long)]
        base: Option<PathBuf>,
        /// Destination `.aicollection`; existing files are never overwritten.
        #[arg(long)]
        out: PathBuf,
        /// Permit archives lacking a stamped genome digest. Chromosome dictionaries must still
        /// match exactly, but reference-sequence identity then cannot be proven.
        #[arg(long)]
        allow_unstamped: bool,
        /// Derive source-root-bound exact intron-span routes from each archive's shape dictionary.
        /// Molecule chunks remain authoritative; collections built without this flag retain the
        /// exact full-shape fallback.
        #[arg(long)]
        shape_routes: bool,
        /// Emit one versioned JSON object with exact source-I/O accounting.
        #[arg(long, conflicts_with = "format")]
        json: bool,
        #[command(flatten)]
        uniform_output: CollectionOutputArgs,
    },
    /// Print the versioned collection manifest and index cardinalities.
    Inspect {
        collection: PathBuf,
        /// Re-hash every source archive in addition to the default exact filesystem guard.
        /// Routed collections also reconstruct every shape route from its bound source shapes.
        #[arg(long)]
        verify_content: bool,
        /// Decode every route block and reconstruct it exactly from the root-bound source shapes.
        #[arg(long)]
        verify_routes: bool,
        #[command(flatten)]
        uniform_output: CollectionOutputArgs,
    },
    /// Exact point-junction counts across routed archives.
    Junction {
        collection: PathBuf,
        /// Junction as chrom:donor-acceptor (0-based boundaries).
        locus: String,
        /// Skip all source archives when the global support upper bound is below this threshold.
        #[arg(long, default_value_t = 0)]
        min_support: u64,
        #[arg(long, default_value_t = 5)]
        top: usize,
        /// Include per-sample routing decisions and planned compressed bytes.
        #[arg(long)]
        explain: bool,
        /// Re-hash every source archive before querying.
        #[arg(long)]
        verify_content: bool,
        /// Emit one versioned JSON object.
        #[arg(long, conflicts_with = "format")]
        json: bool,
        #[command(flatten)]
        uniform_output: CollectionOutputArgs,
    },
    /// Exact anchor-region counts across routed archives.
    Region {
        collection: PathBuf,
        /// Window as chrom:start-end (0-based, half-open).
        locus: String,
        #[arg(long, default_value_t = 5)]
        top: usize,
        #[arg(long)]
        explain: bool,
        #[arg(long)]
        verify_content: bool,
        #[arg(long, conflicts_with = "format")]
        json: bool,
        #[command(flatten)]
        uniform_output: CollectionOutputArgs,
    },
    /// Exact inclusion/exclusion junction-set usage across routed archives.
    Jset {
        collection: PathBuf,
        #[arg(long = "include", required = true, num_args = 1..)]
        include: Vec<String>,
        #[arg(long = "exclude", required = true, num_args = 1..)]
        exclude: Vec<String>,
        /// Require every requested component's global support upper bound to reach this threshold.
        #[arg(long, default_value_t = 0)]
        min_support: u64,
        #[arg(long, default_value_t = 5)]
        top: usize,
        #[arg(long)]
        explain: bool,
        #[arg(long)]
        verify_content: bool,
        #[arg(long, conflicts_with = "format")]
        json: bool,
        #[command(flatten)]
        uniform_output: CollectionOutputArgs,
    },
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct FileIdentity {
    len: u64,
    modified_secs: u64,
    modified_nanos: u32,
    changed_secs: u64,
    changed_nanos: u32,
    dev: u64,
    inode: u64,
    archive_format_version: u32,
    native_scheme: String,
    native_digest: String,
    /// Exact ordered encoded-section identity. Legacy v2 collection manifests predate this field;
    /// collection v3 and later always carry it.
    encoded_sections_digest: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SourceIo {
    id: String,
    format_version: u32,
    identity_scheme: String,
    identity_content_bytes_read: u64,
    total_bytes_read: u64,
    shape_route_payload_bytes_read: u64,
    shape_route_source_bytes_read: u64,
    sections_read: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct IndexedChunk {
    info: ChunkInfo,
    compressed_bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ArchiveEntry {
    id: String,
    path: PathBuf,
    identity: FileIdentity,
    chunks: Vec<IndexedChunk>,
    /// Optional collection-local route binding. The archive ordinal is local to the immutable
    /// layer that owns this entry, never to the flattened collection chain.
    shape_routes: Option<ShapeRouteBinding>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ArchiveRoute {
    archive: u32,
    supporting_children: u64,
    posts: Vec<u32>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct GlobalJunction {
    chrom: u32,
    donor: u32,
    acceptor: u32,
    support_upper_bound: u64,
    routes: Vec<ArchiveRoute>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct BaseCollection {
    path: PathBuf,
    root_digest: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Collection {
    base: Option<BaseCollection>,
    genome_algo: Option<String>,
    genome_digest: Option<String>,
    chroms: Vec<String>,
    chroms_digest: String,
    archives: Vec<ArchiveEntry>,
    junctions: Vec<GlobalJunction>,
    junction_count: usize,
    route_count: usize,
    posting_count: usize,
    /// Derived payloads written only into this layer. They are not serialized in the manifest and
    /// are empty when a manifest is read or when a chain is flattened.
    shape_route_blocks: Vec<RouteBlock>,
    /// Canonically encoded build representation. Production construction uses this flat form to
    /// avoid retaining a separately allocated pair vector for every exact span.
    encoded_shape_route_blocks: Vec<EncodedRouteBlock>,
}

#[derive(Clone, Debug)]
struct LocalJunction {
    chrom: u32,
    donor: u32,
    acceptor: u32,
    supporting_children: u64,
    posts: Vec<u32>,
}

#[derive(Clone, Debug)]
struct LocalArchive {
    entry: ArchiveEntry,
    source_io: SourceIo,
    genome_algo: Option<String>,
    genome_digest: Option<String>,
    chroms: Vec<String>,
    chroms_digest: String,
    junctions: Vec<LocalJunction>,
    shape_routes: Option<LocalShapeRoutes>,
}

#[derive(Clone, Debug)]
struct LocalShapeRoutes {
    shapes_digest: [u8; 32],
    compressed_bytes_read: u64,
    derived: EncodedShapeRoutes,
}

#[derive(Clone, Debug)]
struct Counts {
    molecules: Option<u64>,
    per_cell: FxHashMap<u32, FxHashSet<u32>>,
}

fn checked_usize(value: u64, label: &str) -> Result<usize> {
    let value = usize::try_from(value).with_context(|| format!("{label} exceeds usize"))?;
    if value > MAX_ITEMS {
        bail!("{label} {value} exceeds safety limit {MAX_ITEMS}");
    }
    Ok(value)
}

fn enforce_count(value: usize, maximum: usize, label: &str) -> Result<()> {
    if value > maximum {
        bail!("{label} {value} exceeds safety limit {maximum}");
    }
    Ok(())
}

fn put_string(out: &mut Vec<u8>, value: &str) {
    put_varint(out, value.len() as u64);
    out.extend_from_slice(value.as_bytes());
}

fn get_string(cursor: &mut Cursor<'_>, label: &str) -> Result<String> {
    let len = checked_usize(cursor.varint()?, &format!("{label} length"))?;
    if len > 1_048_576 {
        bail!("{label} length {len} exceeds 1 MiB");
    }
    String::from_utf8(cursor.take(len)?.to_vec()).with_context(|| format!("{label} is not UTF-8"))
}

fn put_optional_string(out: &mut Vec<u8>, value: &Option<String>) {
    out.push(u8::from(value.is_some()));
    if let Some(value) = value {
        put_string(out, value);
    }
}

fn get_optional_string(cursor: &mut Cursor<'_>, label: &str) -> Result<Option<String>> {
    match cursor.byte()? {
        0 => Ok(None),
        1 => Ok(Some(get_string(cursor, label)?)),
        value => bail!("invalid {label} presence flag {value}"),
    }
}

fn digest_from_hex(value: &str, label: &str) -> Result<[u8; 32]> {
    if !valid_digest(value) {
        bail!("{label} is not a canonical lower-case BLAKE3 digest");
    }
    let mut digest = [0u8; 32];
    for (index, byte) in digest.iter_mut().enumerate() {
        *byte = u8::from_str_radix(&value[index * 2..index * 2 + 2], 16)
            .with_context(|| format!("{label} contains invalid hexadecimal"))?;
    }
    Ok(digest)
}

fn encode_shape_route_binding(out: &mut Vec<u8>, binding: &ShapeRouteBinding) -> Result<()> {
    put_varint(out, shaperoute::ROUTE_CODEC_VERSION as u64);
    out.extend_from_slice(&binding.source_root);
    out.extend_from_slice(&binding.shapes_digest);
    put_varint(out, binding.n_shapes as u64);
    put_varint(out, binding.blocks.len() as u64);
    for descriptor in &binding.blocks {
        put_varint(out, descriptor.first_span as u64);
        put_varint(out, descriptor.last_span as u64);
        put_string(out, &descriptor.section_name);
    }
    Ok(())
}

fn decode_shape_route_binding(
    cursor: &mut Cursor<'_>,
    archive_ordinal: u32,
    identity: &FileIdentity,
) -> Result<ShapeRouteBinding> {
    let version = decode_u32(cursor, "shape-route codec version")?;
    if version != shaperoute::ROUTE_CODEC_VERSION {
        bail!(
            "unsupported shape-route codec version {version}; expected {}",
            shaperoute::ROUTE_CODEC_VERSION
        );
    }
    let source_root: [u8; 32] = cursor.take(32)?.try_into().unwrap();
    let shapes_digest: [u8; 32] = cursor.take(32)?.try_into().unwrap();
    let n_shapes = decode_u32(cursor, "shape-route source shape count")?;
    let n_blocks = checked_usize(cursor.varint()?, "shape-route block count")?;
    enforce_count(n_blocks, MAX_SECTIONS, "shape-route block count")?;
    let mut blocks = Vec::with_capacity(n_blocks);
    for _ in 0..n_blocks {
        let first_span = decode_u32(cursor, "shape-route first span")?;
        let last_span = decode_u32(cursor, "shape-route last span")?;
        let section_name = get_string(cursor, "shape-route section name")?;
        blocks.push(shaperoute::RouteBlockDescriptor {
            first_span,
            last_span,
            section_name,
        });
    }
    let binding = ShapeRouteBinding {
        archive_ordinal,
        source_root,
        shapes_digest,
        n_shapes,
        blocks,
    };
    let expected_root = digest_from_hex(
        &identity.native_digest,
        &format!("archive {archive_ordinal} root"),
    )?;
    shaperoute::validate_binding(&binding, archive_ordinal, expected_root, shapes_digest)?;
    Ok(binding)
}

fn encode_collection(collection: &Collection) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    out.push(u8::from(collection.base.is_some()));
    if let Some(base) = &collection.base {
        put_string(
            &mut out,
            base.path
                .to_str()
                .with_context(|| format!("base path {} is not UTF-8", base.path.display()))?,
        );
        put_string(&mut out, &base.root_digest);
    }
    put_optional_string(&mut out, &collection.genome_algo);
    put_optional_string(&mut out, &collection.genome_digest);
    put_varint(&mut out, collection.chroms.len() as u64);
    for chrom in &collection.chroms {
        put_string(&mut out, chrom);
    }
    put_string(&mut out, &collection.chroms_digest);
    put_varint(&mut out, collection.archives.len() as u64);
    for archive in &collection.archives {
        put_string(&mut out, &archive.id);
        put_string(
            &mut out,
            archive
                .path
                .to_str()
                .with_context(|| format!("archive path {} is not UTF-8", archive.path.display()))?,
        );
        for value in [
            archive.identity.len,
            archive.identity.modified_secs,
            archive.identity.modified_nanos as u64,
            archive.identity.changed_secs,
            archive.identity.changed_nanos as u64,
            archive.identity.dev,
            archive.identity.inode,
        ] {
            put_varint(&mut out, value);
        }
        put_varint(&mut out, archive.identity.archive_format_version as u64);
        put_string(&mut out, &archive.identity.native_scheme);
        put_string(&mut out, &archive.identity.native_digest);
        put_string(
            &mut out,
            archive
                .identity
                .encoded_sections_digest
                .as_ref()
                .context("collection archive lacks an encoded-sections digest")?,
        );
        put_varint(&mut out, archive.chunks.len() as u64);
        for chunk in &archive.chunks {
            for value in [
                chunk.info.chrom as u64,
                chunk.info.bin_start as u64,
                chunk.info.n_mols as u64,
                chunk.info.class_base as u64,
                chunk.info.max_anchor as u64,
                chunk.info.n_cells as u64,
                chunk.compressed_bytes,
            ] {
                put_varint(&mut out, value);
            }
        }
        out.push(u8::from(archive.shape_routes.is_some()));
        if let Some(binding) = &archive.shape_routes {
            encode_shape_route_binding(&mut out, binding)?;
        }
    }
    put_varint(&mut out, collection.junction_count as u64);
    put_varint(&mut out, collection.route_count as u64);
    put_varint(&mut out, collection.posting_count as u64);
    put_varint(&mut out, collection.junctions.len() as u64);
    let (mut last_chrom, mut last_donor) = (u32::MAX, 0u32);
    for junction in &collection.junctions {
        if junction.chrom != last_chrom {
            last_chrom = junction.chrom;
            last_donor = 0;
        }
        put_varint(&mut out, junction.chrom as u64);
        put_varint(&mut out, (junction.donor - last_donor) as u64);
        last_donor = junction.donor;
        put_varint(&mut out, (junction.acceptor - junction.donor) as u64);
        put_varint(&mut out, junction.support_upper_bound);
        let presence = presence_bitmap(collection.archives.len(), &junction.routes);
        put_varint(&mut out, presence.len() as u64);
        for word in &presence {
            out.extend_from_slice(&word.to_le_bytes());
        }
        put_varint(&mut out, junction.routes.len() as u64);
        for route in &junction.routes {
            put_varint(&mut out, route.archive as u64);
            put_varint(&mut out, route.supporting_children);
            put_varint(&mut out, route.posts.len() as u64);
            let mut last = 0u32;
            for post in &route.posts {
                put_varint(&mut out, (*post - last) as u64);
                last = *post;
            }
        }
    }
    Ok(out)
}

fn decode_u32(cursor: &mut Cursor<'_>, label: &str) -> Result<u32> {
    u32::try_from(cursor.varint()?).with_context(|| format!("{label} exceeds u32"))
}

fn valid_digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
}

fn validate_identity(identity: &FileIdentity, id: &str, require_encoded: bool) -> Result<()> {
    let expected_scheme = match identity.archive_format_version {
        evidence_io::format::SEEKABLE_VERSION => LEGACY_FULL_FILE_SCHEME,
        evidence_io::format::VERSION => ROOTED_DIRECTORY_SCHEME,
        version => bail!("archive {id} records unsupported format version {version}"),
    };
    if identity.native_scheme != expected_scheme {
        bail!(
            "archive {id} format version {} requires native identity scheme {expected_scheme}, not {}",
            identity.archive_format_version,
            identity.native_scheme
        );
    }
    if !valid_digest(&identity.native_digest) {
        bail!("archive {id} has an invalid native BLAKE3 digest");
    }
    match &identity.encoded_sections_digest {
        Some(digest) if valid_digest(digest) => {}
        Some(_) => bail!("archive {id} has an invalid encoded-sections BLAKE3 digest"),
        None if require_encoded => {
            bail!("archive {id} lacks the encoded-sections digest required by collection v3+")
        }
        None => {}
    }
    Ok(())
}

fn decode_collection(raw: &[u8], format_version: u32) -> Result<Collection> {
    if !matches!(format_version, LEGACY_VERSION | PREVIOUS_VERSION | VERSION) {
        bail!("unsupported collection version {format_version}");
    }
    let mut cursor = Cursor::new(raw);
    let base = match cursor.byte()? {
        0 => None,
        1 => {
            let path = PathBuf::from(get_string(&mut cursor, "base collection path")?);
            let root_digest = get_string(&mut cursor, "base collection root digest")?;
            if root_digest.len() != 64 || !root_digest.bytes().all(|byte| byte.is_ascii_hexdigit())
            {
                bail!("base collection has invalid BLAKE3 root digest");
            }
            Some(BaseCollection { path, root_digest })
        }
        value => bail!("invalid base collection presence flag {value}"),
    };
    let genome_algo = get_optional_string(&mut cursor, "genome algorithm")?;
    let genome_digest = get_optional_string(&mut cursor, "genome digest")?;
    if genome_algo.is_some() != genome_digest.is_some() {
        bail!("collection has only one of genome algorithm and digest");
    }
    let n_chroms = checked_usize(cursor.varint()?, "chromosome count")?;
    enforce_count(n_chroms, MAX_CHROMS, "chromosome count")?;
    if n_chroms > raw.len().saturating_sub(cursor.position()) {
        bail!("chromosome count exceeds remaining manifest bytes");
    }
    let mut chroms = Vec::with_capacity(n_chroms);
    for index in 0..n_chroms {
        chroms.push(get_string(&mut cursor, &format!("chromosome {index}"))?);
    }
    let chroms_digest = get_string(&mut cursor, "chromosome dictionary digest")?;
    if chroms_digest.len() != 64 || !chroms_digest.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        bail!("collection has an invalid chromosome dictionary digest");
    }
    let n_archives = checked_usize(cursor.varint()?, "archive count")?;
    enforce_count(n_archives, MAX_ARCHIVES_PER_LAYER, "archive count")?;
    if n_archives > raw.len().saturating_sub(cursor.position()) {
        bail!("archive count exceeds remaining manifest bytes");
    }
    let mut archives = Vec::with_capacity(n_archives);
    for index in 0..n_archives {
        let id = get_string(&mut cursor, &format!("archive {index} id"))?;
        let path = PathBuf::from(get_string(&mut cursor, &format!("archive {index} path"))?);
        let len = cursor.varint()?;
        let modified_secs = cursor.varint()?;
        let modified_nanos = decode_u32(&mut cursor, "modified nanoseconds")?;
        if modified_nanos >= 1_000_000_000 {
            bail!("archive {id} has invalid modified nanoseconds {modified_nanos}");
        }
        let changed_secs = cursor.varint()?;
        let changed_nanos = decode_u32(&mut cursor, "changed nanoseconds")?;
        if changed_nanos >= 1_000_000_000 {
            bail!("archive {id} has invalid changed nanoseconds {changed_nanos}");
        }
        let dev = cursor.varint()?;
        let inode = cursor.varint()?;
        let (archive_format_version, native_scheme, native_digest, encoded_sections_digest) =
            if format_version == LEGACY_VERSION {
                (
                    evidence_io::format::SEEKABLE_VERSION,
                    LEGACY_FULL_FILE_SCHEME.to_owned(),
                    get_string(&mut cursor, &format!("archive {id} digest"))?.to_ascii_lowercase(),
                    None,
                )
            } else {
                (
                    decode_u32(&mut cursor, "archive format version")?,
                    get_string(&mut cursor, &format!("archive {id} identity scheme"))?,
                    get_string(&mut cursor, &format!("archive {id} native digest"))?,
                    Some(get_string(
                        &mut cursor,
                        &format!("archive {id} encoded-sections digest"),
                    )?),
                )
            };
        let n_chunks = checked_usize(cursor.varint()?, &format!("archive {id} chunk count"))?;
        enforce_count(
            n_chunks,
            MAX_CHUNKS_PER_ARCHIVE,
            &format!("archive {id} chunk count"),
        )?;
        if n_chunks > raw.len().saturating_sub(cursor.position()) / 7 {
            bail!("archive {id} chunk count exceeds remaining manifest bytes");
        }
        let mut chunks = Vec::with_capacity(n_chunks);
        for _ in 0..n_chunks {
            chunks.push(IndexedChunk {
                info: ChunkInfo {
                    chrom: decode_u32(&mut cursor, "chunk chromosome")?,
                    bin_start: decode_u32(&mut cursor, "chunk start")?,
                    n_mols: decode_u32(&mut cursor, "chunk molecule count")?,
                    class_base: decode_u32(&mut cursor, "chunk class base")?,
                    max_anchor: decode_u32(&mut cursor, "chunk max anchor")?,
                    n_cells: decode_u32(&mut cursor, "chunk cell count")?,
                },
                compressed_bytes: cursor.varint()?,
            });
        }
        let identity = FileIdentity {
            len,
            modified_secs,
            modified_nanos,
            changed_secs,
            changed_nanos,
            dev,
            inode,
            archive_format_version,
            native_scheme,
            native_digest,
            encoded_sections_digest,
        };
        validate_identity(&identity, &id, format_version != LEGACY_VERSION)?;
        let shape_routes = if format_version == VERSION {
            match cursor.byte()? {
                0 => None,
                1 => {
                    if identity.archive_format_version != evidence_io::format::VERSION {
                        bail!("shape routes require a root-committed v2 archive source");
                    }
                    Some(decode_shape_route_binding(
                        &mut cursor,
                        index as u32,
                        &identity,
                    )?)
                }
                value => bail!("invalid shape-route presence flag {value}"),
            }
        } else {
            None
        };
        archives.push(ArchiveEntry {
            id,
            path,
            identity,
            chunks,
            shape_routes,
        });
    }
    let junction_count = checked_usize(cursor.varint()?, "declared junction count")?;
    let route_count = checked_usize(cursor.varint()?, "declared route count")?;
    let posting_count = checked_usize(cursor.varint()?, "declared posting count")?;
    let n_junctions = checked_usize(cursor.varint()?, "junction payload row count")?;
    enforce_count(
        n_junctions,
        MAX_JUNCTION_ROWS_PER_SEGMENT,
        "junction payload row count",
    )?;
    let minimum_row_bytes = 6usize.saturating_add(n_archives.div_ceil(64).saturating_mul(8));
    if minimum_row_bytes == 0
        || n_junctions > raw.len().saturating_sub(cursor.position()) / minimum_row_bytes
    {
        bail!("junction payload row count exceeds remaining bytes");
    }
    let mut junctions = Vec::with_capacity(n_junctions);
    let (mut last_chrom, mut last_donor) = (u32::MAX, 0u32);
    for _ in 0..n_junctions {
        let chrom = decode_u32(&mut cursor, "junction chromosome")?;
        if chrom as usize >= chroms.len() {
            bail!("junction chromosome {chrom} is outside the chromosome dictionary");
        }
        if chrom != last_chrom {
            last_chrom = chrom;
            last_donor = 0;
        }
        let donor = last_donor
            .checked_add(decode_u32(&mut cursor, "junction donor delta")?)
            .context("junction donor overflow")?;
        last_donor = donor;
        let acceptor = donor
            .checked_add(decode_u32(&mut cursor, "junction span")?)
            .context("junction acceptor overflow")?;
        let support_upper_bound = cursor.varint()?;
        let n_words = checked_usize(cursor.varint()?, "presence bitmap words")?;
        let expected_words = n_archives.div_ceil(64);
        if n_words != expected_words {
            bail!("presence bitmap has {n_words} words; expected {expected_words}");
        }
        let mut presence = Vec::with_capacity(n_words);
        for _ in 0..n_words {
            presence.push(u64::from_le_bytes(cursor.take(8)?.try_into().unwrap()));
        }
        let n_routes = checked_usize(cursor.varint()?, "junction route count")?;
        if n_routes > n_archives || n_routes > raw.len().saturating_sub(cursor.position()) / 3 {
            bail!("junction route count exceeds archive or payload bounds");
        }
        let mut routes = Vec::with_capacity(n_routes);
        let mut previous_archive = None;
        for _ in 0..n_routes {
            let archive = decode_u32(&mut cursor, "route archive")?;
            if archive as usize >= n_archives {
                bail!("junction route references missing archive {archive}");
            }
            if previous_archive.is_some_and(|previous| archive <= previous) {
                bail!("junction routes are not strictly archive-sorted");
            }
            previous_archive = Some(archive);
            let supporting_children = cursor.varint()?;
            let n_posts = checked_usize(cursor.varint()?, "junction posting count")?;
            let chunk_count = archives[archive as usize].chunks.len();
            if n_posts == 0
                || n_posts > chunk_count
                || n_posts > raw.len().saturating_sub(cursor.position())
            {
                bail!("junction posting count exceeds chunk or payload bounds");
            }
            let mut posts = Vec::with_capacity(n_posts);
            let mut last = 0u32;
            for _ in 0..n_posts {
                let delta = decode_u32(&mut cursor, "junction posting delta")?;
                if !posts.is_empty() && delta == 0 {
                    bail!("junction postings are not strictly increasing");
                }
                last = last
                    .checked_add(delta)
                    .context("junction posting overflow")?;
                if last as usize >= archives[archive as usize].chunks.len() {
                    bail!("junction route references missing chunk {last}");
                }
                if archives[archive as usize].chunks[last as usize].info.chrom != chrom {
                    bail!("junction posting references a chunk on the wrong chromosome");
                }
                posts.push(last);
            }
            if presence[archive as usize / 64] & (1u64 << (archive as usize % 64)) == 0 {
                bail!("junction route is absent from its presence bitmap");
            }
            routes.push(ArchiveRoute {
                archive,
                supporting_children,
                posts,
            });
        }
        let support_sum = routes.iter().try_fold(0u64, |sum, route| {
            sum.checked_add(route.supporting_children)
                .context("junction support sum overflow")
        })?;
        if support_sum != support_upper_bound {
            bail!("junction support upper bound {support_upper_bound} differs from route sum {support_sum}");
        }
        if presence
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum::<usize>()
            != routes.len()
        {
            bail!("junction presence bitmap contains a route-free archive");
        }
        junctions.push(GlobalJunction {
            chrom,
            donor,
            acceptor,
            support_upper_bound,
            routes,
        });
    }
    if !cursor.is_empty() {
        bail!(
            "collection payload has {} trailing bytes",
            raw.len() - cursor.position()
        );
    }
    if archives.windows(2).any(|pair| pair[0].id >= pair[1].id) {
        bail!("collection archives are not strictly ID-sorted");
    }
    if junctions.windows(2).any(|pair| {
        (pair[0].chrom, pair[0].donor, pair[0].acceptor)
            >= (pair[1].chrom, pair[1].donor, pair[1].acceptor)
    }) {
        bail!("collection junctions are not strictly coordinate-sorted");
    }
    if !junctions.is_empty() {
        let observed_routes: usize = junctions.iter().map(|row| row.routes.len()).sum();
        let observed_postings: usize = junctions
            .iter()
            .flat_map(|row| &row.routes)
            .map(|route| route.posts.len())
            .sum();
        if junctions.len() != junction_count
            || observed_routes != route_count
            || observed_postings != posting_count
        {
            bail!("collection declared junction cardinalities do not match its payload");
        }
    }
    Ok(Collection {
        base,
        genome_algo,
        genome_digest,
        chroms,
        chroms_digest,
        archives,
        junction_count,
        junctions,
        route_count,
        posting_count,
        shape_route_blocks: Vec::new(),
        encoded_shape_route_blocks: Vec::new(),
    })
}

fn presence_bitmap(n_archives: usize, routes: &[ArchiveRoute]) -> Vec<u64> {
    let mut presence = vec![0u64; n_archives.div_ceil(64)];
    for route in routes {
        presence[route.archive as usize / 64] |= 1u64 << (route.archive as usize % 64);
    }
    presence
}

fn encode_junction_rows(rows: &[GlobalJunction], n_archives: usize) -> Vec<u8> {
    let mut out = Vec::new();
    put_varint(&mut out, rows.len() as u64);
    let (mut last_chrom, mut last_donor) = (u32::MAX, 0u32);
    for junction in rows {
        if junction.chrom != last_chrom {
            last_chrom = junction.chrom;
            last_donor = 0;
        }
        put_varint(&mut out, junction.chrom as u64);
        put_varint(&mut out, (junction.donor - last_donor) as u64);
        last_donor = junction.donor;
        put_varint(&mut out, (junction.acceptor - junction.donor) as u64);
        put_varint(&mut out, junction.support_upper_bound);
        let presence = presence_bitmap(n_archives, &junction.routes);
        put_varint(&mut out, presence.len() as u64);
        for word in &presence {
            out.extend_from_slice(&word.to_le_bytes());
        }
        put_varint(&mut out, junction.routes.len() as u64);
        for route in &junction.routes {
            put_varint(&mut out, route.archive as u64);
            put_varint(&mut out, route.supporting_children);
            put_varint(&mut out, route.posts.len() as u64);
            let mut last = 0u32;
            for post in &route.posts {
                put_varint(&mut out, (*post - last) as u64);
                last = *post;
            }
        }
    }
    out
}

fn decode_junction_rows(
    raw: &[u8],
    chroms: &[String],
    archives: &[ArchiveEntry],
) -> Result<Vec<GlobalJunction>> {
    let mut cursor = Cursor::new(raw);
    let n_rows = checked_usize(cursor.varint()?, "junction segment row count")?;
    enforce_count(
        n_rows,
        MAX_JUNCTION_ROWS_PER_SEGMENT,
        "junction segment row count",
    )?;
    let expected_words = archives.len().div_ceil(64);
    let minimum_row_bytes = 6usize.saturating_add(expected_words.saturating_mul(8));
    if minimum_row_bytes == 0
        || n_rows > raw.len().saturating_sub(cursor.position()) / minimum_row_bytes
    {
        bail!("junction segment row count exceeds remaining bytes");
    }
    let mut rows = Vec::with_capacity(n_rows);
    let (mut last_chrom, mut last_donor) = (u32::MAX, 0u32);
    for _ in 0..n_rows {
        let chrom = decode_u32(&mut cursor, "junction chromosome")?;
        if chrom as usize >= chroms.len() {
            bail!("junction chromosome {chrom} is outside the chromosome dictionary");
        }
        if chrom != last_chrom {
            last_chrom = chrom;
            last_donor = 0;
        }
        let donor = last_donor
            .checked_add(decode_u32(&mut cursor, "junction donor delta")?)
            .context("junction donor overflow")?;
        last_donor = donor;
        let acceptor = donor
            .checked_add(decode_u32(&mut cursor, "junction span")?)
            .context("junction acceptor overflow")?;
        let support_upper_bound = cursor.varint()?;
        let n_words = checked_usize(cursor.varint()?, "presence bitmap words")?;
        if n_words != expected_words {
            bail!("presence bitmap has {n_words} words; expected {expected_words}");
        }
        let mut presence = Vec::with_capacity(n_words);
        for _ in 0..n_words {
            presence.push(u64::from_le_bytes(cursor.take(8)?.try_into().unwrap()));
        }
        let n_routes = checked_usize(cursor.varint()?, "junction route count")?;
        if n_routes > archives.len() || n_routes > raw.len().saturating_sub(cursor.position()) / 3 {
            bail!("junction route count exceeds archive or payload bounds");
        }
        let mut routes = Vec::with_capacity(n_routes);
        let mut previous_archive = None;
        for _ in 0..n_routes {
            let archive = decode_u32(&mut cursor, "route archive")?;
            if archive as usize >= archives.len() {
                bail!("junction route references missing archive {archive}");
            }
            if previous_archive.is_some_and(|previous| archive <= previous) {
                bail!("junction routes are not strictly archive-sorted");
            }
            previous_archive = Some(archive);
            let supporting_children = cursor.varint()?;
            let n_posts = checked_usize(cursor.varint()?, "junction posting count")?;
            let chunk_count = archives[archive as usize].chunks.len();
            if n_posts == 0
                || n_posts > chunk_count
                || n_posts > raw.len().saturating_sub(cursor.position())
            {
                bail!("junction posting count exceeds chunk or payload bounds");
            }
            let mut posts = Vec::with_capacity(n_posts);
            let mut last = 0u32;
            for _ in 0..n_posts {
                let delta = decode_u32(&mut cursor, "junction posting delta")?;
                if !posts.is_empty() && delta == 0 {
                    bail!("junction postings are not strictly increasing");
                }
                last = last
                    .checked_add(delta)
                    .context("junction posting overflow")?;
                if last as usize >= archives[archive as usize].chunks.len() {
                    bail!("junction route references missing chunk {last}");
                }
                if archives[archive as usize].chunks[last as usize].info.chrom != chrom {
                    bail!("junction posting references a chunk on the wrong chromosome");
                }
                posts.push(last);
            }
            if presence[archive as usize / 64] & (1u64 << (archive as usize % 64)) == 0 {
                bail!("junction route is absent from its presence bitmap");
            }
            routes.push(ArchiveRoute {
                archive,
                supporting_children,
                posts,
            });
        }
        let support_sum = routes.iter().try_fold(0u64, |sum, route| {
            sum.checked_add(route.supporting_children)
                .context("junction support sum overflow")
        })?;
        if support_sum != support_upper_bound {
            bail!("junction support upper bound differs from route sum");
        }
        if presence
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum::<usize>()
            != routes.len()
        {
            bail!("junction presence bitmap contains a route-free archive");
        }
        rows.push(GlobalJunction {
            chrom,
            donor,
            acceptor,
            support_upper_bound,
            routes,
        });
    }
    if !cursor.is_empty() {
        bail!(
            "junction segment has {} trailing bytes",
            raw.len() - cursor.position()
        );
    }
    if rows.windows(2).any(|pair| {
        (pair[0].chrom, pair[0].donor, pair[0].acceptor)
            >= (pair[1].chrom, pair[1].donor, pair[1].acceptor)
    }) {
        bail!("junction segment is not strictly coordinate-sorted");
    }
    Ok(rows)
}

#[derive(Clone, Debug)]
struct CollectionSection {
    name: String,
    offset: u64,
    raw_len: u64,
    compressed_len: u64,
    digest: [u8; 32],
}

struct CollectionFile {
    file: std::fs::File,
    sections: Vec<CollectionSection>,
    root_digest: [u8; 32],
    version: u32,
    bytes_read: std::sync::atomic::AtomicU64,
}

fn parse_segment_name(name: &str) -> Result<(u32, u32)> {
    let fields: Vec<&str> = name.split('.').collect();
    if fields.len() != 3 || fields[0] != "j" {
        bail!("collection layer contains unknown section {name}");
    }
    let chrom: u32 = fields[1]
        .parse()
        .with_context(|| format!("invalid chromosome in collection section {name}"))?;
    let bin: u32 = fields[2]
        .parse()
        .with_context(|| format!("invalid genomic bin in collection section {name}"))?;
    let canonical_name = format!("j.{chrom}.{bin}");
    if name != canonical_name {
        bail!("collection section {name} is not canonically named {canonical_name}");
    }
    Ok((chrom, bin))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CollectionSectionKind {
    Junction(u32, u32),
    ShapeRoute(u32, u32),
}

fn parse_collection_section_name(name: &str) -> Result<CollectionSectionKind> {
    if name.starts_with("j.") {
        let (chrom, bin) = parse_segment_name(name)?;
        Ok(CollectionSectionKind::Junction(chrom, bin))
    } else if name.starts_with("s.") {
        let (archive, first_span) = shaperoute::parse_route_section_name(name)?;
        Ok(CollectionSectionKind::ShapeRoute(archive, first_span))
    } else {
        bail!("collection layer contains unknown section {name}")
    }
}

impl CollectionFile {
    fn open(path: &Path) -> Result<Self> {
        use std::io::{Seek, SeekFrom};
        let file =
            std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let file_len = file.metadata()?.len();
        if file_len < 48 {
            bail!("collection is truncated: {file_len} bytes");
        }
        // A route-bearing collection can have thousands of small directory records. Buffer both
        // authentication passes so opening it does not issue several tiny reads per record. The
        // authenticated bytes, validation order, and logical byte accounting remain unchanged.
        let mut file = BufReader::with_capacity(1024 * 1024, file);
        let mut fixed = [0u8; 48];
        file.read_exact(&mut fixed)?;
        if &fixed[..8] != MAGIC {
            bail!("not a .aicollection file (bad magic)");
        }
        let version = u32::from_le_bytes(fixed[8..12].try_into().unwrap());
        if !matches!(version, LEGACY_VERSION | PREVIOUS_VERSION | VERSION) {
            bail!(
                "unsupported collection version {version}; expected {LEGACY_VERSION}, {PREVIOUS_VERSION}, or {VERSION}"
            );
        }
        let n = checked_usize(
            u32::from_le_bytes(fixed[12..16].try_into().unwrap()) as u64,
            "collection section count",
        )?;
        if n > MAX_SECTIONS {
            bail!("collection section count {n} exceeds safety limit {MAX_SECTIONS}");
        }
        let mut root_digest = [0u8; 32];
        root_digest.copy_from_slice(&fixed[16..48]);
        let mut directory_hasher = blake3::Hasher::new();
        directory_hasher.update(&fixed[..16]);
        // Authenticate the complete variable-length directory before trusting its cardinalities,
        // names, lengths, or payload routes. The directory commits to each section's raw digest;
        // selected payloads are then verified lazily against that authenticated digest.
        for _ in 0..n {
            let mut name_len = [0u8; 1];
            file.read_exact(&mut name_len)?;
            directory_hasher.update(&name_len);
            let mut name = vec![0u8; name_len[0] as usize];
            file.read_exact(&mut name)?;
            directory_hasher.update(&name);
            let mut fields = [0u8; 48];
            file.read_exact(&mut fields)?;
            directory_hasher.update(&fields);
        }
        let payload_start = file.stream_position()?;
        if directory_hasher.finalize().as_bytes() != &root_digest {
            bail!("collection directory checksum mismatch");
        }
        file.seek(SeekFrom::Start(48))?;
        let mut pending = Vec::with_capacity(n);
        let mut names = FxHashSet::default();
        for _ in 0..n {
            let mut name_len = [0u8; 1];
            file.read_exact(&mut name_len)?;
            let mut name = vec![0u8; name_len[0] as usize];
            file.read_exact(&mut name)?;
            let name = String::from_utf8(name).context("collection section name is not UTF-8")?;
            if name.is_empty() || !names.insert(name.clone()) {
                bail!("collection has an empty or duplicate section name");
            }
            let mut value = [0u8; 8];
            file.read_exact(&mut value)?;
            let raw_len = u64::from_le_bytes(value);
            file.read_exact(&mut value)?;
            let compressed_len = u64::from_le_bytes(value);
            if raw_len > MAX_SECTION_RAW_BYTES {
                bail!("collection section {name} raw length exceeds safety limit");
            }
            let permitted_raw = compressed_len
                .saturating_mul(MAX_COMPRESSION_RATIO)
                .saturating_add(64 * 1024 * 1024);
            if raw_len > permitted_raw {
                bail!("collection section {name} declares an unsafe compression ratio");
            }
            let mut digest = [0u8; 32];
            file.read_exact(&mut digest)?;
            pending.push((name, raw_len, compressed_len, digest));
        }
        if pending.first().map(|row| row.0.as_str()) != Some("manifest") {
            bail!("collection manifest must be the first section");
        }
        let mut previous_junction = None;
        let mut previous_shape_route = None;
        let mut saw_shape_route = false;
        for (name, _, _, _) in pending.iter().skip(1) {
            match parse_collection_section_name(name)? {
                CollectionSectionKind::Junction(chrom, bin) => {
                    if saw_shape_route {
                        bail!("collection junction sections must precede all shape-route sections");
                    }
                    let segment = (chrom, bin);
                    if previous_junction.is_some_and(|previous| segment <= previous) {
                        bail!("collection junction sections are not in canonical genomic order");
                    }
                    previous_junction = Some(segment);
                }
                CollectionSectionKind::ShapeRoute(archive, first_span) => {
                    if version != VERSION {
                        bail!("shape-route sections require collection format v{VERSION}");
                    }
                    saw_shape_route = true;
                    let segment = (archive, first_span);
                    if previous_shape_route.is_some_and(|previous| segment <= previous) {
                        bail!("collection shape-route sections are not in canonical archive/span order");
                    }
                    previous_shape_route = Some(segment);
                }
            }
        }
        if file.stream_position()? != payload_start {
            bail!("collection directory length changed between validation passes");
        }
        let mut offset = payload_start;
        let mut sections = Vec::with_capacity(n);
        for (name, raw_len, compressed_len, digest) in pending {
            sections.push(CollectionSection {
                name,
                offset,
                raw_len,
                compressed_len,
                digest,
            });
            offset = offset
                .checked_add(compressed_len)
                .context("collection section extent overflow")?;
        }
        if offset != file_len {
            bail!("collection section lengths do not match file length");
        }
        file.seek(SeekFrom::Start(payload_start))?;
        Ok(Self {
            file: file.into_inner(),
            sections,
            root_digest,
            version,
            bytes_read: std::sync::atomic::AtomicU64::new(
                payload_start
                    .checked_mul(2)
                    .and_then(|value| value.checked_sub(48))
                    .context("collection directory byte count overflow")?,
            ),
        })
    }

    fn names(&self) -> impl Iterator<Item = &str> {
        self.sections.iter().map(|section| section.name.as_str())
    }

    fn root_digest_hex(&self) -> String {
        self.root_digest
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }

    fn has(&self, name: &str) -> bool {
        self.sections.iter().any(|section| section.name == name)
    }

    fn bytes_read(&self) -> u64 {
        self.bytes_read.load(std::sync::atomic::Ordering::Relaxed)
    }

    fn compressed_len(&self, name: &str) -> Result<u64> {
        self.sections
            .iter()
            .find(|section| section.name == name)
            .map(|section| section.compressed_len)
            .with_context(|| format!("collection missing section {name}"))
    }

    fn read(&self, name: &str) -> Result<Vec<u8>> {
        use std::io::{Seek, SeekFrom};
        let section = self
            .sections
            .iter()
            .find(|section| section.name == name)
            .with_context(|| format!("collection missing section {name}"))?;
        let compressed_len = usize::try_from(section.compressed_len)
            .context("collection section compressed length exceeds usize")?;
        let raw_len = usize::try_from(section.raw_len)
            .context("collection section raw length exceeds usize")?;
        let mut file = self.file.try_clone()?;
        file.seek(SeekFrom::Start(section.offset))?;
        let mut compressed = vec![0u8; compressed_len];
        file.read_exact(&mut compressed)?;
        self.bytes_read
            .fetch_add(section.compressed_len, std::sync::atomic::Ordering::Relaxed);
        let raw = evidence_io::format::decompress(&compressed, raw_len)?;
        if blake3::hash(&raw).as_bytes() != &section.digest {
            bail!("collection section {name} checksum mismatch");
        }
        Ok(raw)
    }
}

struct CollectionLayer {
    path: PathBuf,
    file: CollectionFile,
    manifest: Collection,
    local_to_global: Vec<usize>,
}

struct CollectionChain {
    layers: Vec<CollectionLayer>,
    collection: Collection,
    global_to_local: Vec<(usize, usize)>,
}

const COLLECTION_BUILD_RESULT_SCHEMA: &str = "gravlax.collection.build.result.v1";
const COLLECTION_BUILD_ARCHIVES_SCHEMA: &str = "gravlax.collection.build.archives.v1";
const COLLECTION_BUILD_SOURCE_IO_SCHEMA: &str = "gravlax.collection.build.source-io.v1";
const COLLECTION_BUILD_SOURCE_SECTIONS_SCHEMA: &str = "gravlax.collection.build.source-sections.v1";
const COLLECTION_INSPECT_RESULT_SCHEMA: &str = "gravlax.collection.inspect.result.v1";
const COLLECTION_INSPECT_LAYERS_SCHEMA: &str = "gravlax.collection.inspect.layers.v1";
const COLLECTION_INSPECT_CHROMS_SCHEMA: &str = "gravlax.collection.inspect.chromosomes.v1";
const COLLECTION_INSPECT_ARCHIVES_SCHEMA: &str = "gravlax.collection.inspect.archives.v1";
const COLLECTION_INSPECT_ROUTE_BLOCKS_SCHEMA: &str =
    "gravlax.collection.inspect.shape-route-blocks.v1";
const COLLECTION_JUNCTION_RESULT_SCHEMA: &str = "gravlax.collection.junction.result.v1";
const COLLECTION_JUNCTION_SAMPLES_SCHEMA: &str = "gravlax.collection.junction.samples.v1";
const COLLECTION_JUNCTION_CELLS_SCHEMA: &str = "gravlax.collection.junction.cells.v1";
const COLLECTION_REGION_RESULT_SCHEMA: &str = "gravlax.collection.region.result.v1";
const COLLECTION_REGION_SAMPLES_SCHEMA: &str = "gravlax.collection.region.samples.v1";
const COLLECTION_REGION_CELLS_SCHEMA: &str = "gravlax.collection.region.cells.v1";
const COLLECTION_JSET_RESULT_SCHEMA: &str = "gravlax.collection.jset.result.v1";
const COLLECTION_JSET_REQUESTS_SCHEMA: &str = "gravlax.collection.jset.requests.v1";
const COLLECTION_JSET_SAMPLES_SCHEMA: &str = "gravlax.collection.jset.samples.v1";
const COLLECTION_JSET_CELLS_SCHEMA: &str = "gravlax.collection.jset.cells.v1";

fn set_table_schema(
    id: &'static str,
    fields: Vec<Field>,
    key: &[&str],
) -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(id, fields)?
        .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(key.iter().copied()))
}

fn sequence_table_schema(
    id: &'static str,
    fields: Vec<Field>,
    key: &str,
) -> std::result::Result<TableSchema, OutputError> {
    TableSchema::new(id, fields)?.with_semantics(
        TableSemantics::new(RowSemantics::Sequence)
            .with_key([key])
            .ordered_by([gravlax_output::OrderKey::ascending(key)]),
    )
}

fn uniform_path<'a>(path: &'a Path, role: &str) -> Result<&'a str> {
    path.to_str()
        .with_context(|| format!("{role} path is not valid UTF-8: {}", path.display()))
}

/// Run one uniform serializer against either locked stdout or the shared atomic no-clobber file
/// publisher. The producer is invoked exactly once, and a failed producer never installs a file.
fn write_uniform_collection_result<F>(output: &CollectionOutputArgs, produce: F) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> std::result::Result<(), OutputError>,
{
    output
        .format
        .context("uniform collection output requires an explicit --format")?;
    if let Some(path) = output.output.as_deref() {
        let outcome =
            publish_file_no_clobber(path, Durability::Flush, |writer| produce(&mut *writer))?;
        for warning in outcome.warnings {
            eprintln!("warning: {warning}");
        }
    } else {
        let stdout = std::io::stdout();
        let mut writer = BufWriter::new(stdout.lock());
        produce(&mut writer)?;
    }
    Ok(())
}

/// Avoid an expensive collection scan when a requested result destination is already occupied.
/// The publisher remains the authoritative race-safe no-clobber check.
fn preflight_uniform_collection_output(output: &CollectionOutputArgs) -> Result<()> {
    let Some(path) = output.output.as_deref() else {
        return Ok(());
    };
    if path.file_name().is_none() {
        bail!("uniform output path must name a file: {}", path.display());
    }
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let parent_metadata = std::fs::metadata(parent).with_context(|| {
        format!(
            "checking parent directory {} for uniform output {}",
            parent.display(),
            path.display()
        )
    })?;
    if !parent_metadata.is_dir() {
        bail!(
            "uniform output parent is not a directory: {}",
            parent.display()
        );
    }
    match std::fs::symlink_metadata(path) {
        Ok(_) => bail!("refusing to replace existing output {}", path.display()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error)
            .with_context(|| format!("checking uniform output destination {}", path.display())),
    }
}

fn destination_key(path: &Path) -> Result<PathBuf> {
    let file_name = path
        .file_name()
        .with_context(|| format!("output path must name a file: {}", path.display()))?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    Ok(std::fs::canonicalize(parent)
        .with_context(|| format!("resolving output parent {}", parent.display()))?
        .join(file_name))
}

fn archive_provenance_identity(archive: &ArchiveEntry) -> String {
    format!(
        "{}:{}",
        archive.identity.native_scheme, archive.identity.native_digest
    )
}

fn uniform_source_context<'a, I>(
    archives: I,
    mut parameters: BTreeMap<String, serde_json::Value>,
) -> Result<ResultContext>
where
    I: IntoIterator<Item = &'a ArchiveEntry>,
{
    let mut by_id = archives
        .into_iter()
        .map(|archive| (archive.id.as_str(), archive_provenance_identity(archive)))
        .collect::<Vec<_>>();
    by_id.sort_unstable_by(|left, right| left.0.cmp(right.0));
    let mut sample_identities = serde_json::Map::new();
    for (sample, identity) in &by_id {
        if sample_identities
            .insert((*sample).to_owned(), json!(identity))
            .is_some()
        {
            bail!("duplicate sample id {sample} in uniform collection provenance");
        }
    }
    parameters.insert(
        "sample_identities".into(),
        serde_json::Value::Object(sample_identities),
    );
    let context = ResultContext {
        producer: Producer {
            name: "aie".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        },
        provenance: Provenance {
            archives: by_id.into_iter().map(|(_, identity)| identity).collect(),
            parameters,
            ..Default::default()
        },
        warnings: Vec::new(),
    };
    context.validate()?;
    Ok(context)
}

fn collection_uniform_context(
    requested_path: &Path,
    chain: &CollectionChain,
    mut parameters: BTreeMap<String, serde_json::Value>,
) -> Result<ResultContext> {
    parameters.insert(
        "collection_path".into(),
        json!(uniform_path(requested_path, "collection")?),
    );
    parameters.insert(
        "collection_layers".into(),
        serde_json::Value::Array(
            chain
                .layers
                .iter()
                .map(|layer| -> Result<_> {
                    Ok(json!({
                        "path": uniform_path(&layer.path, "collection layer")?,
                        "format_version": layer.file.version,
                        "root": format!("aicollection-directory-root-v1:{}", layer.file.root_digest_hex()),
                    }))
                })
                .collect::<Result<Vec<_>>>()?,
        ),
    );
    parameters.insert(
        "reference".into(),
        match (
            chain.collection.genome_algo.as_deref(),
            chain.collection.genome_digest.as_deref(),
        ) {
            (Some(algo), Some(digest)) => {
                json!({"stamped": true, "algo": algo, "digest": digest})
            }
            _ => json!({"stamped": false}),
        },
    );
    parameters.insert(
        "chromosome_dictionary_digest".into(),
        json!(&chain.collection.chroms_digest),
    );
    uniform_source_context(chain.collection.archives.iter(), parameters)
}

fn segment_name(chrom: u32, donor: u32) -> String {
    format!("j.{chrom}.{}", donor / JUNCTION_BIN_BP)
}

fn encode_manifest(collection: &Collection) -> Result<Vec<u8>> {
    let mut manifest = collection.clone();
    manifest.junctions.clear();
    manifest.shape_route_blocks.clear();
    manifest.encoded_shape_route_blocks.clear();
    encode_collection(&manifest)
}

fn read_manifest(file: &CollectionFile) -> Result<Collection> {
    let manifest = decode_collection(&file.read("manifest")?, file.version)?;
    if !manifest.junctions.is_empty() {
        bail!("collection manifest unexpectedly embeds junction rows");
    }
    let observed: Vec<&str> = file.names().filter(|name| name.starts_with("s.")).collect();
    let expected: Vec<&str> = manifest
        .archives
        .iter()
        .flat_map(|archive| {
            archive
                .shape_routes
                .iter()
                .flat_map(|binding| binding.blocks.iter())
                .map(|descriptor| descriptor.section_name.as_str())
        })
        .collect();
    if observed != expected {
        bail!("collection shape-route sections do not exactly match the manifest bindings");
    }
    Ok(manifest)
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct CollectionWriteStats {
    raw_bytes: u64,
    file_bytes: u64,
    root_digest: [u8; 32],
    shape_route_sections: usize,
    shape_route_raw_bytes: u64,
    shape_route_compressed_bytes: u64,
}

fn install_collection_output(
    output: &std::fs::File,
    temporary: &Path,
    destination: &Path,
) -> Result<u64> {
    let file_bytes = output.metadata()?.len();
    let outcome = install_open_file_no_clobber(
        output,
        temporary,
        destination,
        Durability::FileAndDirectory,
    )?;
    for warning in outcome.warnings {
        eprintln!("warning: {warning}");
    }
    Ok(file_bytes)
}

fn write_collection(path: &Path, collection: &Collection) -> Result<CollectionWriteStats> {
    if path.exists() {
        bail!("refusing to overwrite {}", path.display());
    }
    let mut junction_sections = Vec::new();
    let mut begin = 0usize;
    while begin < collection.junctions.len() {
        let row = &collection.junctions[begin];
        let key = (row.chrom, row.donor / JUNCTION_BIN_BP);
        let mut end = begin + 1;
        while end < collection.junctions.len()
            && (
                collection.junctions[end].chrom,
                collection.junctions[end].donor / JUNCTION_BIN_BP,
            ) == key
        {
            end += 1;
        }
        junction_sections.push((format!("j.{}.{}", key.0, key.1), begin, end));
        begin = end;
    }
    let expected_route_sections: Vec<(&ShapeRouteBinding, &shaperoute::RouteBlockDescriptor)> =
        collection
            .archives
            .iter()
            .filter_map(|archive| archive.shape_routes.as_ref())
            .flat_map(|binding| {
                binding
                    .blocks
                    .iter()
                    .map(move |descriptor| (binding, descriptor))
            })
            .collect();
    if !collection.shape_route_blocks.is_empty()
        && !collection.encoded_shape_route_blocks.is_empty()
    {
        bail!("collection mixes decoded and encoded shape-route build representations");
    }
    let route_block_count = collection
        .shape_route_blocks
        .len()
        .checked_add(collection.encoded_shape_route_blocks.len())
        .context("collection shape-route block count overflow")?;
    if expected_route_sections.len() != route_block_count {
        bail!(
            "collection has {} route descriptors but {} route blocks",
            expected_route_sections.len(),
            route_block_count
        );
    }
    for ((binding, descriptor), block) in expected_route_sections
        .iter()
        .zip(&collection.shape_route_blocks)
    {
        if block.archive_ordinal != binding.archive_ordinal
            || block.n_shapes != binding.n_shapes
            || block.spans.first().map(|row| row.span) != Some(descriptor.first_span)
            || block.spans.last().map(|row| row.span) != Some(descriptor.last_span)
            || descriptor.section_name
                != shaperoute::route_section_name(block.archive_ordinal, descriptor.first_span)
        {
            bail!("collection shape-route block differs from its manifest descriptor");
        }
    }
    for ((binding, descriptor), block) in expected_route_sections
        .iter()
        .zip(&collection.encoded_shape_route_blocks)
    {
        if block.archive_ordinal != binding.archive_ordinal
            || block.n_shapes != binding.n_shapes
            || block.first_span != descriptor.first_span
            || block.last_span != descriptor.last_span
            || block.descriptor() != **descriptor
        {
            bail!("encoded collection shape-route block differs from its manifest descriptor");
        }
    }
    let shape_route_sections = expected_route_sections.len();
    let mut section_names = Vec::with_capacity(
        1usize
            .checked_add(junction_sections.len())
            .and_then(|count| count.checked_add(shape_route_sections))
            .context("collection section count overflow")?,
    );
    section_names.push("manifest".to_owned());
    section_names.extend(junction_sections.iter().map(|(name, _, _)| name.clone()));
    section_names.extend(
        expected_route_sections
            .iter()
            .map(|(_, descriptor)| descriptor.section_name.clone()),
    );
    let section_count =
        u32::try_from(section_names.len()).context("too many collection sections")?;
    let mut header = Vec::with_capacity(16);
    header.extend_from_slice(MAGIC);
    header.extend_from_slice(&VERSION.to_le_bytes());
    header.extend_from_slice(&section_count.to_le_bytes());
    let directory_len = section_names.iter().try_fold(0usize, |sum, name| {
        u8::try_from(name.len()).context("collection section name is too long")?;
        sum.checked_add(1 + name.len() + 8 + 8 + 32)
            .context("collection directory length overflow")
    })?;
    let payload_start = header
        .len()
        .checked_add(32)
        .and_then(|value| value.checked_add(directory_len))
        .context("collection payload offset overflow")?;
    let file_name = path
        .file_name()
        .context("collection output has no file name")?
        .to_string_lossy();
    let (temporary, output) = (0u32..1_000)
        .find_map(|attempt| {
            let candidate =
                path.with_file_name(format!(".{file_name}.tmp-{}-{attempt}", std::process::id()));
            match std::fs::OpenOptions::new()
                .read(true)
                .write(true)
                .create_new(true)
                .open(&candidate)
            {
                Ok(file) => Some(Ok((candidate, file))),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => None,
                Err(error) => Some(Err(error)),
            }
        })
        .transpose()
        .with_context(|| format!("creating temporary output beside {}", path.display()))?
        .context("could not allocate a unique collection temporary file")?;
    let staging_guard = match output.try_clone() {
        Ok(file) => file,
        Err(error) => {
            let _ = remove_staging_if_owned(&temporary, &output);
            return Err(error).context("cloning collection temporary file for cleanup guard");
        }
    };
    let write_result = (|| -> Result<CollectionWriteStats> {
        use std::io::{Seek, SeekFrom};

        let mut writer = BufWriter::new(output);
        writer.seek(SeekFrom::Start(
            u64::try_from(payload_start).context("collection payload offset exceeds u64")?,
        ))?;
        let mut directory_rows = Vec::with_capacity(section_names.len());
        let mut raw_bytes = 0u64;
        let mut shape_route_raw_bytes = 0u64;
        let mut shape_route_compressed_bytes = 0u64;
        {
            let mut write_section = |name: &str, raw: &[u8]| -> Result<()> {
                let raw_len = u64::try_from(raw.len()).context("collection section exceeds u64")?;
                let digest = *blake3::hash(raw).as_bytes();
                let compressed = evidence_io::format::compress(raw, 9)?;
                let compressed_len =
                    u64::try_from(compressed.len()).context("compressed section exceeds u64")?;
                writer.write_all(&compressed)?;
                raw_bytes = raw_bytes
                    .checked_add(raw_len)
                    .context("collection raw section byte count overflow")?;
                if name.starts_with("s.") {
                    shape_route_raw_bytes = shape_route_raw_bytes
                        .checked_add(raw_len)
                        .context("shape-route raw byte count overflow")?;
                    shape_route_compressed_bytes = shape_route_compressed_bytes
                        .checked_add(compressed_len)
                        .context("shape-route compressed byte count overflow")?;
                }
                directory_rows.push((name.to_owned(), raw_len, compressed_len, digest));
                Ok(())
            };
            let manifest_raw = encode_manifest(collection)?;
            write_section("manifest", &manifest_raw)?;
            for (name, begin, end) in &junction_sections {
                let raw = encode_junction_rows(
                    &collection.junctions[*begin..*end],
                    collection.archives.len(),
                );
                write_section(name, &raw)?;
            }
            if collection.encoded_shape_route_blocks.is_empty() {
                for ((_, descriptor), block) in expected_route_sections
                    .iter()
                    .zip(&collection.shape_route_blocks)
                {
                    let raw = shaperoute::encode_block(block)?;
                    write_section(&descriptor.section_name, &raw)?;
                }
            } else {
                for ((_, descriptor), block) in expected_route_sections
                    .iter()
                    .zip(&collection.encoded_shape_route_blocks)
                {
                    write_section(&descriptor.section_name, &block.raw)?;
                }
            }
        }
        let observed_names: Vec<&str> = directory_rows
            .iter()
            .map(|(name, _, _, _)| name.as_str())
            .collect();
        if observed_names != section_names.iter().map(String::as_str).collect::<Vec<_>>() {
            bail!("collection sections were not written in their planned order");
        }
        let mut directory = Vec::with_capacity(directory_len);
        for (name, raw_len, compressed_len, digest) in &directory_rows {
            directory.push(u8::try_from(name.len()).context("collection section name is too long")?);
            directory.extend_from_slice(name.as_bytes());
            directory.extend_from_slice(&raw_len.to_le_bytes());
            directory.extend_from_slice(&compressed_len.to_le_bytes());
            directory.extend_from_slice(digest);
        }
        if directory.len() != directory_len {
            bail!("collection directory length differs from its planned length");
        }
        let mut directory_hasher = blake3::Hasher::new();
        directory_hasher.update(&header);
        directory_hasher.update(&directory);
        let root_digest = directory_hasher.finalize();
        writer.seek(SeekFrom::Start(0))?;
        writer.write_all(&header)?;
        writer.write_all(root_digest.as_bytes())?;
        writer.write_all(&directory)?;
        writer.flush()?;
        let output = writer.into_inner().map_err(|error| error.into_error())?;
        let file_bytes = install_collection_output(&output, &temporary, path)?;
        Ok(CollectionWriteStats {
            raw_bytes,
            file_bytes,
            root_digest: *root_digest.as_bytes(),
            shape_route_sections,
            shape_route_raw_bytes,
            shape_route_compressed_bytes,
        })
    })();
    if write_result.is_err() {
        match remove_staging_if_owned(&temporary, &staging_guard) {
            Ok(true) => {}
            Ok(false) => eprintln!(
                "warning: collection temporary path {} now names a different file and was preserved",
                temporary.display()
            ),
            Err(error) => eprintln!(
                "warning: could not remove owned collection temporary path {}: {error:#}",
                temporary.display()
            ),
        }
    }
    write_result
}

fn open_collection_manifest(path: &Path) -> Result<(CollectionFile, Collection)> {
    let file = CollectionFile::open(path)?;
    let collection = read_manifest(&file)?;
    Ok((file, collection))
}

/// Return the authenticated collection-directory identity without decoding section payloads.
/// This is the cheap provenance path used by resolved analysis plans.
pub fn native_collection_identity(path: &Path) -> Result<String> {
    Ok(CollectionFile::open(path)?.root_digest_hex())
}

fn load_collection_layers(
    path: &Path,
    depth: usize,
    seen: &mut FxHashSet<PathBuf>,
    layers: &mut Vec<CollectionLayer>,
) -> Result<String> {
    if depth >= 32 {
        bail!("collection base chain exceeds the 32-layer safety limit; compact it");
    }
    let canonical = std::fs::canonicalize(path)
        .with_context(|| format!("resolving collection layer {}", path.display()))?;
    if !seen.insert(canonical.clone()) {
        bail!("collection base chain contains a cycle or repeated layer");
    }
    let (file, manifest) = open_collection_manifest(&canonical)?;
    if let Some(base) = &manifest.base {
        let base_path = if base.path.is_absolute() {
            base.path.clone()
        } else {
            canonical
                .parent()
                .unwrap_or_else(|| Path::new("."))
                .join(&base.path)
        };
        let observed = load_collection_layers(&base_path, depth + 1, seen, layers)?;
        if observed != base.root_digest {
            bail!(
                "base collection digest mismatch for {}; expected {}, observed {}",
                base_path.display(),
                base.root_digest,
                observed
            );
        }
    }
    let root_digest = file.root_digest_hex();
    layers.push(CollectionLayer {
        path: canonical,
        file,
        local_to_global: vec![0; manifest.archives.len()],
        manifest,
    });
    Ok(root_digest)
}

fn open_collection_chain(path: &Path) -> Result<CollectionChain> {
    let mut layers = Vec::new();
    load_collection_layers(path, 0, &mut FxHashSet::default(), &mut layers)?;
    let has_v2 = layers
        .iter()
        .any(|layer| layer.file.version == LEGACY_VERSION);
    let has_encoded_identity = layers
        .iter()
        .any(|layer| layer.file.version != LEGACY_VERSION);
    if has_v2 && has_encoded_identity {
        bail!(
            "mixed collection v2 and v3/v4 chains cannot prove cross-container duplicate evidence; rebuild the v2 layers"
        );
    }
    let first = layers.first().context("collection chain is empty")?;
    let genome_algo = first.manifest.genome_algo.clone();
    let genome_digest = first.manifest.genome_digest.clone();
    let chroms = first.manifest.chroms.clone();
    let chroms_digest = first.manifest.chroms_digest.clone();
    for layer in &layers {
        if layer.manifest.genome_algo != genome_algo
            || layer.manifest.genome_digest != genome_digest
            || layer.manifest.chroms != chroms
            || layer.manifest.chroms_digest != chroms_digest
        {
            bail!(
                "collection layer {} has an incompatible reference identity",
                layer.path.display()
            );
        }
    }

    let mut by_id: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut by_path: FxHashMap<PathBuf, String> = FxHashMap::default();
    let mut by_inode: FxHashMap<(u64, u64), String> = FxHashMap::default();
    let mut by_digest: FxHashMap<String, String> = FxHashMap::default();
    for (layer_index, layer) in layers.iter().enumerate() {
        for (local_index, archive) in layer.manifest.archives.iter().enumerate() {
            if let Some((old_layer, _old_local)) =
                by_id.insert(archive.id.clone(), (layer_index, local_index))
            {
                bail!(
                    "duplicate sample id {} in collection layers {} and {}",
                    archive.id,
                    layers[old_layer].path.display(),
                    layer.path.display()
                );
            }
            if let Some(previous) = by_path.insert(archive.path.clone(), archive.id.clone()) {
                bail!(
                    "samples {previous} and {} reuse resolved archive {}",
                    archive.id,
                    archive.path.display()
                );
            }
            let identity = &archive.identity;
            if identity.dev != 0 || identity.inode != 0 {
                if let Some(previous) =
                    by_inode.insert((identity.dev, identity.inode), archive.id.clone())
                {
                    bail!(
                        "samples {previous} and {} reuse the same archive inode",
                        archive.id
                    );
                }
            }
            let duplicate_key = identity.encoded_sections_digest.clone().unwrap_or_else(|| {
                format!("{}:{}", identity.native_scheme, identity.native_digest)
            });
            if let Some(previous) = by_digest.insert(duplicate_key, archive.id.clone()) {
                bail!(
                    "samples {previous} and {} have byte-identical archive content",
                    archive.id
                );
            }
        }
    }
    let mut archives = Vec::with_capacity(by_id.len());
    let mut global_to_local = Vec::with_capacity(by_id.len());
    for (global, (_, (layer, local))) in by_id.iter().enumerate() {
        layers[*layer].local_to_global[*local] = global;
        archives.push(layers[*layer].manifest.archives[*local].clone());
        global_to_local.push((*layer, *local));
    }
    let junction_count = layers.iter().try_fold(0usize, |sum, layer| {
        sum.checked_add(layer.manifest.junction_count)
            .context("collection junction-row count overflow")
    })?;
    let route_count = layers.iter().try_fold(0usize, |sum, layer| {
        sum.checked_add(layer.manifest.route_count)
            .context("collection route count overflow")
    })?;
    let posting_count = layers.iter().try_fold(0usize, |sum, layer| {
        sum.checked_add(layer.manifest.posting_count)
            .context("collection posting count overflow")
    })?;
    Ok(CollectionChain {
        layers,
        collection: Collection {
            base: None,
            genome_algo,
            genome_digest,
            chroms,
            chroms_digest,
            archives,
            junctions: Vec::new(),
            junction_count,
            route_count,
            posting_count,
            shape_route_blocks: Vec::new(),
            encoded_shape_route_blocks: Vec::new(),
        },
        global_to_local,
    })
}

fn lookup_chain_junction(
    chain: &CollectionChain,
    chrom: u32,
    donor: u32,
    acceptor: u32,
) -> Result<Option<GlobalJunction>> {
    let mut routes = Vec::new();
    let mut support_upper_bound = 0u64;
    for layer in &chain.layers {
        if let Some(row) = lookup_junction(&layer.file, &layer.manifest, chrom, donor, acceptor)? {
            support_upper_bound = support_upper_bound
                .checked_add(row.support_upper_bound)
                .context("collection support upper bound overflow across layers")?;
            for mut route in row.routes {
                route.archive = u32::try_from(layer.local_to_global[route.archive as usize])
                    .context("global archive index exceeds u32")?;
                routes.push(route);
            }
        }
    }
    routes.sort_unstable_by_key(|route| route.archive);
    Ok((!routes.is_empty()).then_some(GlobalJunction {
        chrom,
        donor,
        acceptor,
        support_upper_bound,
        routes,
    }))
}

fn lookup_chain_junctions(
    chain: &CollectionChain,
    loci: &[(u32, u32, u32)],
) -> Result<Vec<Option<GlobalJunction>>> {
    let mut supports = vec![0u64; loci.len()];
    let mut routes: Vec<Vec<ArchiveRoute>> = (0..loci.len()).map(|_| Vec::new()).collect();
    for layer in &chain.layers {
        let mut bins: BTreeMap<(u32, u32), Vec<GlobalJunction>> = BTreeMap::new();
        for &(chrom, donor, _) in loci {
            let key = (chrom, donor / JUNCTION_BIN_BP);
            if let std::collections::btree_map::Entry::Vacant(entry) = bins.entry(key) {
                entry.insert(load_junction_segment(
                    &layer.file,
                    &layer.manifest,
                    chrom,
                    donor,
                )?);
            }
        }
        for (index, &(chrom, donor, acceptor)) in loci.iter().enumerate() {
            let rows = &bins[&(chrom, donor / JUNCTION_BIN_BP)];
            if let Ok(row_index) = rows.binary_search_by_key(&(chrom, donor, acceptor), |row| {
                (row.chrom, row.donor, row.acceptor)
            }) {
                let row = &rows[row_index];
                supports[index] = supports[index]
                    .checked_add(row.support_upper_bound)
                    .context("collection support upper bound overflow across layers")?;
                for route in &row.routes {
                    let mut route = route.clone();
                    route.archive = u32::try_from(layer.local_to_global[route.archive as usize])
                        .context("global archive index exceeds u32")?;
                    routes[index].push(route);
                }
            }
        }
    }
    Ok(loci
        .iter()
        .copied()
        .zip(supports)
        .zip(routes)
        .map(
            |(((chrom, donor, acceptor), support_upper_bound), mut routes)| {
                routes.sort_unstable_by_key(|route| route.archive);
                (!routes.is_empty()).then_some(GlobalJunction {
                    chrom,
                    donor,
                    acceptor,
                    support_upper_bound,
                    routes,
                })
            },
        )
        .collect())
}

#[derive(Clone, Debug)]
struct ArchiveShapeRoutePlan {
    binding: ShapeRouteBinding,
    spans: BTreeMap<u32, SpanRoute>,
    section_names: Vec<String>,
    sidecar_compressed_bytes: u64,
}

fn lookup_chain_shape_routes(
    chain: &CollectionChain,
    global_archive: usize,
    requested_spans: &[u32],
) -> Result<Option<ArchiveShapeRoutePlan>> {
    let &(layer_index, local_archive) = chain
        .global_to_local
        .get(global_archive)
        .with_context(|| format!("missing collection owner for archive {global_archive}"))?;
    let layer = &chain.layers[layer_index];
    let archive = layer
        .manifest
        .archives
        .get(local_archive)
        .context("collection owner refers to a missing local archive")?;
    let Some(binding) = &archive.shape_routes else {
        return Ok(None);
    };
    let mut descriptor_spans: BTreeMap<usize, Vec<u32>> = BTreeMap::new();
    for &exact_span in requested_spans {
        let descriptor_index = binding
            .blocks
            .partition_point(|descriptor| descriptor.first_span <= exact_span)
            .checked_sub(1);
        let descriptor_index = descriptor_index
            .filter(|&index| exact_span <= binding.blocks[index].last_span)
            .with_context(|| {
                format!(
                    "shape-route binding for sample {} lacks catalogued intron span {exact_span}",
                    archive.id
                )
            })?;
        descriptor_spans
            .entry(descriptor_index)
            .or_default()
            .push(exact_span);
    }
    let mut spans = BTreeMap::new();
    let mut section_names = Vec::new();
    let mut sidecar_compressed_bytes = 0u64;
    for (descriptor_index, mut exact_spans) in descriptor_spans {
        exact_spans.sort_unstable();
        exact_spans.dedup();
        let descriptor = &binding.blocks[descriptor_index];
        let raw = layer.file.read(&descriptor.section_name)?;
        let block = shaperoute::decode_block(
            &raw,
            local_archive as u32,
            binding.n_shapes,
            descriptor.first_span,
            descriptor.last_span,
        )?;
        for exact_span in exact_spans {
            let span = block.span(exact_span).cloned().with_context(|| {
                format!(
                    "shape-route section {} omits catalogued intron span {exact_span}",
                    descriptor.section_name
                )
            })?;
            spans.insert(exact_span, span);
        }
        section_names.push(descriptor.section_name.clone());
        sidecar_compressed_bytes = sidecar_compressed_bytes
            .checked_add(layer.file.compressed_len(&descriptor.section_name)?)
            .context("shape-route sidecar byte count overflow")?;
    }
    Ok(Some(ArchiveShapeRoutePlan {
        binding: binding.clone(),
        spans,
        section_names,
        sidecar_compressed_bytes,
    }))
}

fn collection_sidecar_bytes_read(chain: &CollectionChain) -> Result<u64> {
    chain.layers.iter().try_fold(0u64, |sum, layer| {
        sum.checked_add(layer.file.bytes_read())
            .context("collection-sidecar byte count overflow")
    })
}

fn load_junction_segment(
    file: &CollectionFile,
    collection: &Collection,
    chrom: u32,
    donor: u32,
) -> Result<Vec<GlobalJunction>> {
    let name = segment_name(chrom, donor);
    if !file.has(&name) {
        return Ok(Vec::new());
    }
    let rows = decode_junction_rows(&file.read(&name)?, &collection.chroms, &collection.archives)?;
    let bin = donor / JUNCTION_BIN_BP;
    if rows
        .iter()
        .any(|row| row.chrom != chrom || row.donor / JUNCTION_BIN_BP != bin)
    {
        bail!("collection section {name} contains a junction outside its routing bin");
    }
    Ok(rows)
}

fn lookup_junction(
    file: &CollectionFile,
    collection: &Collection,
    chrom: u32,
    donor: u32,
    acceptor: u32,
) -> Result<Option<GlobalJunction>> {
    let rows = load_junction_segment(file, collection, chrom, donor)?;
    Ok(rows
        .binary_search_by_key(&(chrom, donor, acceptor), |row| {
            (row.chrom, row.donor, row.acceptor)
        })
        .ok()
        .map(|index| rows[index].clone()))
}

fn digest_hex(digest: [u8; 32]) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn identity_from_metadata(
    metadata: std::fs::Metadata,
    archive_format_version: u32,
    native_scheme: String,
    native_digest: String,
    encoded_sections_digest: Option<String>,
) -> Result<FileIdentity> {
    if !metadata.is_file() {
        bail!("archive source is not a regular file");
    }
    let modified = metadata
        .modified()
        .context("archive modification time is unavailable")?
        .duration_since(std::time::UNIX_EPOCH)
        .context("archive modification time precedes Unix epoch")?;
    #[cfg(unix)]
    let (dev, inode, changed_secs, changed_nanos) = {
        use std::os::unix::fs::MetadataExt;
        let changed_secs =
            u64::try_from(metadata.ctime()).context("archive change time precedes Unix epoch")?;
        let changed_nanos = u32::try_from(metadata.ctime_nsec())
            .context("archive change-time nanoseconds are invalid")?;
        if changed_nanos >= 1_000_000_000 {
            bail!("archive change-time nanoseconds are invalid");
        }
        (metadata.dev(), metadata.ino(), changed_secs, changed_nanos)
    };
    #[cfg(not(unix))]
    let (dev, inode, changed_secs, changed_nanos) = (0, 0, 0, 0);
    Ok(FileIdentity {
        len: metadata.len(),
        modified_secs: modified.as_secs(),
        modified_nanos: modified.subsec_nanos(),
        changed_secs,
        changed_nanos,
        dev,
        inode,
        archive_format_version,
        native_scheme,
        native_digest,
        encoded_sections_digest,
    })
}

fn metadata_identity(
    path: &Path,
    archive_format_version: u32,
    native_scheme: String,
    native_digest: String,
    encoded_sections_digest: Option<String>,
) -> Result<FileIdentity> {
    let metadata = std::fs::metadata(path)
        .with_context(|| format!("reading metadata for {}", path.display()))?;
    identity_from_metadata(
        metadata,
        archive_format_version,
        native_scheme,
        native_digest,
        encoded_sections_digest,
    )
}

fn stat_identity(path: &Path) -> Result<FileIdentity> {
    metadata_identity(path, 0, String::new(), String::new(), None)
}

fn identity_stat_matches(expected: &FileIdentity, observed: &FileIdentity) -> bool {
    expected.len == observed.len
        && expected.modified_secs == observed.modified_secs
        && expected.modified_nanos == observed.modified_nanos
        && expected.changed_secs == observed.changed_secs
        && expected.changed_nanos == observed.changed_nanos
        && expected.dev == observed.dev
        && expected.inode == observed.inode
}

fn validate_archive_identity(archive: &ArchiveEntry, verify_content: bool) -> Result<SourceIo> {
    let file = std::fs::File::open(&archive.path)
        .with_context(|| format!("opening archive identity {}", archive.path.display()))?;
    validate_archive_identity_file(archive, file, verify_content)
}

fn validate_archive_identity_file(
    archive: &ArchiveEntry,
    file: std::fs::File,
    verify_content: bool,
) -> Result<SourceIo> {
    let observed = identity_from_metadata(file.metadata()?, 0, String::new(), String::new(), None)?;
    if !identity_stat_matches(&archive.identity, &observed) {
        bail!(
            "archive identity changed for sample {} at {}; rebuild the collection index",
            archive.id,
            archive.path.display()
        );
    }
    // A normal legacy-v1 guard remains stat-only. Its v3+ manifest nevertheless records the
    // scheme-independent identity established during construction. Re-reading that identity would
    // require the full-file scan that rooted archives are designed to eliminate.
    if archive.identity.archive_format_version == evidence_io::format::SEEKABLE_VERSION
        && !verify_content
    {
        return Ok(SourceIo {
            id: archive.id.clone(),
            format_version: archive.identity.archive_format_version,
            identity_scheme: archive.identity.native_scheme.clone(),
            identity_content_bytes_read: 0,
            total_bytes_read: 0,
            shape_route_payload_bytes_read: 0,
            shape_route_source_bytes_read: 0,
            sections_read: Vec::new(),
        });
    }

    // Parse and verify the same open file description whose stat tuple was checked above. This
    // prevents a replaceable archive path from splitting the stat and content guards across two
    // different inodes.
    let reader = evidence_io::format::SectionReader::from_file(file)?;
    if reader.archive_version() != archive.identity.archive_format_version {
        bail!(
            "archive format version changed for sample {} at {}; rebuild the collection index",
            archive.id,
            archive.path.display()
        );
    }
    let identity_content_bytes_read = if reader.archive_version() == evidence_io::format::VERSION {
        let root = reader
            .content_commitment()
            .context("rooted archive lacks its directory commitment")?;
        if archive.identity.native_scheme != ROOTED_DIRECTORY_SCHEME
            || root.to_hex() != archive.identity.native_digest
        {
            bail!(
                "archive root changed for sample {} at {}; rebuild the collection index",
                archive.id,
                archive.path.display()
            );
        }
        let encoded = reader
            .encoded_content_identity()?
            .context("rooted archive lacks its encoded-sections identity")?;
        let encoded = digest_hex(encoded);
        if archive.identity.encoded_sections_digest.as_deref() != Some(encoded.as_str()) {
            bail!(
                    "archive encoded-sections identity changed for sample {} at {}; rebuild the collection index",
                    archive.id,
                    archive.path.display()
                );
        }
        if let Some(binding) = &archive.shape_routes {
            let shapes_digest = reader
                .section_metadata()
                .find(|section| section.name == "shapes")
                .context("routed archive lacks its shapes section")?
                .compressed_blake3
                .context("routed archive shapes entry lacks its committed payload digest")?;
            shaperoute::validate_binding(
                binding,
                binding.archive_ordinal,
                root.digest,
                shapes_digest,
            )
            .with_context(|| {
                format!("validating shape-route source binding for sample {}", archive.id)
            })?;
        }
        if verify_content {
            reader.verify_all_payloads()?
        } else {
            0
        }
    } else {
        let scan = reader.scan_legacy_identities()?;
        let full = digest_hex(scan.full_file_blake3);
        let encoded = digest_hex(scan.encoded_sections_blake3);
        if archive.identity.native_scheme != LEGACY_FULL_FILE_SCHEME
            || full != archive.identity.native_digest
            || archive
                .identity
                .encoded_sections_digest
                .as_deref()
                .is_some_and(|expected| expected != encoded)
        {
            bail!(
                "archive content digest changed for sample {} at {}; rebuild the collection index",
                archive.id,
                archive.path.display()
            );
        }
        scan.bytes_read
    };
    let sections_read =
        if verify_content && reader.archive_version() == evidence_io::format::VERSION {
            reader.names().map(str::to_owned).collect()
        } else {
            Vec::new()
        };
    let after = identity_from_metadata(
        reader.file_metadata()?,
        0,
        String::new(),
        String::new(),
        None,
    )?;
    if !identity_stat_matches(&archive.identity, &after) {
        bail!(
            "archive identity changed while validating sample {} at {}; rebuild the collection index",
            archive.id,
            archive.path.display()
        );
    }
    Ok(SourceIo {
        id: archive.id.clone(),
        format_version: reader.archive_version(),
        identity_scheme: archive.identity.native_scheme.clone(),
        identity_content_bytes_read,
        total_bytes_read: reader.bytes_read(),
        shape_route_payload_bytes_read: 0,
        shape_route_source_bytes_read: 0,
        sections_read,
    })
}

fn parse_sample(value: &str) -> Result<(String, PathBuf)> {
    let (id, path) = value.split_once('=').context("--sample must be ID=PATH")?;
    if id.is_empty()
        || !id
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
    {
        bail!("sample id {id:?} must contain only ASCII letters, digits, '.', '_' or '-'");
    }
    if path.is_empty() {
        bail!("sample {id} has an empty archive path");
    }
    Ok((id.to_owned(), PathBuf::from(path)))
}

fn parse_source_digest(value: &str) -> Result<(String, String)> {
    let (id, digest) = value
        .split_once('=')
        .context("--source-digest must be ID=BLAKE3")?;
    if id.is_empty()
        || !id
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
    {
        bail!("source digest id {id:?} is invalid");
    }
    if digest.len() != 64 || !digest.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        bail!("source digest for {id} is not a 64-character BLAKE3 hex digest");
    }
    Ok((id.to_owned(), digest.to_ascii_lowercase()))
}

fn parse_local_junctions(
    archive: &mut LazyArchive,
    chunks: &[IndexedChunk],
) -> Result<Vec<LocalJunction>> {
    let catalogue = archive.reader().read("index.junctions")?;
    let postings = archive.reader().read("index.jpost")?;
    let mut coordinates = Cursor::new(&catalogue);
    let mut rows = Vec::new();
    let (mut last_chrom, mut last_donor) = (u32::MAX, 0u32);
    while !coordinates.is_empty() {
        let chrom = decode_u32(&mut coordinates, "junction chromosome")?;
        if chrom != last_chrom {
            last_chrom = chrom;
            last_donor = 0;
        }
        let donor = last_donor
            .checked_add(decode_u32(&mut coordinates, "junction donor delta")?)
            .context("junction donor overflow")?;
        last_donor = donor;
        let acceptor = donor
            .checked_add(decode_u32(&mut coordinates, "junction span")?)
            .context("junction acceptor overflow")?;
        rows.push(LocalJunction {
            chrom,
            donor,
            acceptor,
            supporting_children: 0,
            posts: Vec::new(),
        });
    }
    if rows.windows(2).any(|pair| {
        (pair[0].chrom, pair[0].donor, pair[0].acceptor)
            >= (pair[1].chrom, pair[1].donor, pair[1].acceptor)
    }) {
        bail!("archive junction catalogue is not strictly coordinate-sorted");
    }
    let mut posting_cursor = Cursor::new(&postings);
    for row in &mut rows {
        row.supporting_children = posting_cursor.varint()?;
        let n = checked_usize(posting_cursor.varint()?, "junction posting count")?;
        if n > chunks.len() || n > postings.len().saturating_sub(posting_cursor.position()) {
            bail!("junction posting count exceeds chunk or payload bounds");
        }
        let mut last = 0u32;
        let mut previous = None;
        row.posts.reserve(n);
        for _ in 0..n {
            last = last
                .checked_add(decode_u32(&mut posting_cursor, "junction posting delta")?)
                .context("junction posting overflow")?;
            if previous.is_some_and(|value| last <= value) {
                bail!("junction postings are not strictly chunk-sorted");
            }
            let chunk = chunks
                .get(last as usize)
                .with_context(|| format!("junction posting references missing chunk {last}"))?;
            if chunk.info.chrom != row.chrom {
                bail!(
                    "junction posting references chromosome {} chunk for chromosome {} junction",
                    chunk.info.chrom,
                    row.chrom
                );
            }
            row.posts.push(last);
            previous = Some(last);
        }
    }
    if !posting_cursor.is_empty() {
        bail!("junction postings contain more rows than the catalogue");
    }
    Ok(rows)
}

fn load_local_archive(
    id: String,
    path: PathBuf,
    expected_digest: Option<String>,
    build_shape_routes: bool,
    route_archive_ordinal: u32,
) -> Result<LocalArchive> {
    let canonical = std::fs::canonicalize(&path)
        .with_context(|| format!("resolving sample {id} archive {}", path.display()))?;
    let before = stat_identity(&canonical)?;
    let mut archive = LazyArchive::open(&canonical)
        .with_context(|| format!("opening sample {id} archive {}", canonical.display()))?;
    let chroms = archive.chrom_names.clone();
    let chroms_digest = archive.chrom_digest.clone();
    let genome_algo = archive.genome_sig.as_ref().map(|sig| sig.algo.clone());
    let genome_digest = archive.genome_sig.as_ref().map(|sig| sig.digest.clone());
    let chunk_info = read_chunk_index(archive.reader())?;
    let compressed: FxHashMap<usize, u64> = archive
        .reader()
        .entries()
        .iter()
        .filter_map(|(name, _, _, compressed)| {
            name.strip_prefix('c')
                .and_then(|suffix| suffix.parse::<usize>().ok())
                .map(|index| (index, *compressed))
        })
        .collect();
    let chunks: Vec<IndexedChunk> = chunk_info
        .into_iter()
        .enumerate()
        .map(|(index, info)| {
            Ok(IndexedChunk {
                info,
                compressed_bytes: *compressed
                    .get(&index)
                    .with_context(|| format!("archive missing chunk section c{index}"))?,
            })
        })
        .collect::<Result<_>>()?;
    if let Some(chunk) = chunks
        .iter()
        .find(|chunk| chunk.info.chrom as usize >= chroms.len())
    {
        bail!(
            "sample {id} chunk references chromosome id {} outside its dictionary",
            chunk.info.chrom
        );
    }
    let junctions = parse_local_junctions(&mut archive, &chunks)?;
    if let Some(row) = junctions
        .iter()
        .find(|row| row.chrom as usize >= chroms.len())
    {
        bail!(
            "sample {id} junction references chromosome id {} outside its dictionary",
            row.chrom
        );
    }
    if let Some(row) = junctions
        .iter()
        .find(|row| row.posts.iter().any(|&post| post as usize >= chunks.len()))
    {
        bail!(
            "sample {id} junction {}:{}-{} references a missing chunk",
            chroms
                .get(row.chrom as usize)
                .map(String::as_str)
                .unwrap_or("<bad-chrom>"),
            row.donor,
            row.acceptor
        );
    }
    let shape_routes = if build_shape_routes {
        if archive.reader().archive_version() != evidence_io::format::VERSION {
            bail!("sample {id} must be a root-committed v2 archive to derive shape routes");
        }
        let shapes_digest = archive
            .reader()
            .section_metadata()
            .find(|section| section.name == "shapes")
            .context("archive lacks its shapes section")?
            .compressed_blake3
            .context("rooted archive shapes entry lacks its committed payload digest")?;
        let compressed_bytes_read = archive
            .reader()
            .section_metadata()
            .find(|section| section.name == "shapes")
            .context("archive lacks its shapes section")?
            .compressed_len;
        let raw = archive.reader().read("shapes")?;
        let derived = shaperoute::derive_encoded_from_shapes(&raw, route_archive_ordinal)
            .with_context(|| format!("deriving local shape routes for sample {id}"))?;
        Some(LocalShapeRoutes {
            shapes_digest,
            compressed_bytes_read,
            derived,
        })
    } else {
        None
    };
    let reader = archive.reader();
    let format_version = reader.archive_version();
    let (identity_scheme, native_digest, encoded_sections_digest, identity_content_bytes_read) =
        if format_version == evidence_io::format::VERSION {
            let root = reader
                .content_commitment()
                .context("rooted archive lacks its directory commitment")?;
            let encoded = reader
                .encoded_content_identity()?
                .context("rooted archive lacks its encoded-sections identity")?;
            (
                ROOTED_DIRECTORY_SCHEME.to_owned(),
                root.to_hex(),
                digest_hex(encoded),
                0,
            )
        } else if format_version == evidence_io::format::SEEKABLE_VERSION {
            let scan = reader.scan_legacy_identities()?;
            (
                LEGACY_FULL_FILE_SCHEME.to_owned(),
                digest_hex(scan.full_file_blake3),
                digest_hex(scan.encoded_sections_blake3),
                scan.bytes_read,
            )
        } else {
            bail!("sample {id} archive version {format_version} is not seekable v1 or rooted v2");
        };
    if let Some(expected) = expected_digest {
        if native_digest != expected {
            bail!(
                "sample {id} source digest mismatch: expected {expected}, observed {native_digest}"
            );
        }
    }
    let total_bytes_read = reader.bytes_read();
    let compressed_by_name: FxHashMap<String, u64> = reader
        .entries()
        .iter()
        .map(|(name, _, _, compressed_len)| (name.clone(), *compressed_len))
        .collect();
    let after_metadata = reader.file_metadata()?;
    let identity = identity_from_metadata(
        after_metadata,
        format_version,
        identity_scheme.clone(),
        native_digest,
        Some(encoded_sections_digest),
    )?;
    if !identity_stat_matches(&before, &identity) {
        bail!("sample {id} archive changed while the collection index was being built");
    }
    validate_identity(&identity, &id, true)?;
    let mut sections_read: Vec<String> = [
        "chroms",
        "index.chunks",
        "index.jpost",
        "index.junctions",
        "meta",
        "rans.tables",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect();
    if build_shape_routes {
        sections_read.push("shapes".to_owned());
    }
    sections_read.sort_unstable();
    let selected_payload_bytes = sections_read.iter().try_fold(0u64, |sum, name| {
        sum.checked_add(
            *compressed_by_name
                .get(name)
                .with_context(|| format!("sample {id} lacks declared section {name}"))?,
        )
        .context("selected source-payload byte count overflow")
    })?;
    let structural_bytes = total_bytes_read
        .checked_sub(selected_payload_bytes)
        .context("source structural byte accounting underflow")?;
    let shape_route_payload_bytes_read = if build_shape_routes {
        shape_routes
            .as_ref()
            .context("shape-route source accounting lacks derived routes")?
            .compressed_bytes_read
    } else {
        0
    };
    let shape_route_source_bytes_read = if build_shape_routes {
        structural_bytes
            .checked_add(
                shape_route_payload_bytes_read,
            )
            .context("shape-route source byte count overflow")?
    } else {
        0
    };
    let source_io = SourceIo {
        id: id.clone(),
        format_version,
        identity_scheme,
        identity_content_bytes_read,
        total_bytes_read,
        shape_route_payload_bytes_read,
        shape_route_source_bytes_read,
        sections_read,
    };
    Ok(LocalArchive {
        entry: ArchiveEntry {
            id,
            path: canonical,
            identity,
            chunks,
            shape_routes: None,
        },
        source_io,
        genome_algo,
        genome_digest,
        chroms,
        chroms_digest,
        junctions,
        shape_routes,
    })
}

fn assemble_collection(
    mut locals: Vec<LocalArchive>,
    allow_unstamped: bool,
    minimum_archives: usize,
) -> Result<Collection> {
    if locals.len() < minimum_archives {
        bail!("this collection build requires at least {minimum_archives} archive(s)");
    }
    locals.sort_unstable_by(|left, right| left.entry.id.cmp(&right.entry.id));
    if let Some(pair) = locals
        .windows(2)
        .find(|pair| pair[0].entry.id == pair[1].entry.id)
    {
        bail!("duplicate sample id {}", pair[0].entry.id);
    }
    let mut source_ids: FxHashMap<&Path, &str> = FxHashMap::default();
    let mut inode_ids: FxHashMap<(u64, u64), &str> = FxHashMap::default();
    let mut digest_ids: FxHashMap<&str, &str> = FxHashMap::default();
    for local in &locals {
        if let Some(previous) =
            source_ids.insert(local.entry.path.as_path(), local.entry.id.as_str())
        {
            bail!(
                "samples {previous} and {} reuse resolved archive {}; each source archive may appear only once",
                local.entry.id,
                local.entry.path.display()
            );
        }
        let identity = &local.entry.identity;
        if identity.dev != 0 || identity.inode != 0 {
            if let Some(previous) =
                inode_ids.insert((identity.dev, identity.inode), local.entry.id.as_str())
            {
                bail!(
                    "samples {previous} and {} reuse the same archive inode; hard-link aliases are not independent samples",
                    local.entry.id
                );
            }
        }
        let encoded = identity
            .encoded_sections_digest
            .as_deref()
            .context("new collection source lacks an encoded-sections digest")?;
        if let Some(previous) = digest_ids.insert(encoded, local.entry.id.as_str()) {
            bail!(
                "samples {previous} and {} have identical archive content; duplicate evidence cannot be counted twice",
                local.entry.id
            );
        }
    }
    let chroms = locals[0].chroms.clone();
    let chroms_digest = locals[0].chroms_digest.clone();
    let genome_algo = locals[0].genome_algo.clone();
    let genome_digest = locals[0].genome_digest.clone();
    if genome_digest.is_none() && !allow_unstamped {
        bail!("archives lack a stamped genome identity; pass --allow-unstamped to accept chromosome-only guards");
    }
    for local in &locals {
        if local.chroms != chroms {
            bail!(
                "sample {} has an incompatible chromosome dictionary",
                local.entry.id
            );
        }
        if local.chroms_digest != chroms_digest {
            bail!(
                "sample {} has a non-identical chromosome dictionary encoding",
                local.entry.id
            );
        }
        if local.genome_algo != genome_algo || local.genome_digest != genome_digest {
            bail!(
                "sample {} has an incompatible stamped genome identity",
                local.entry.id
            );
        }
    }

    let mut encoded_shape_route_blocks = Vec::new();
    for (archive_ordinal, local) in locals.iter_mut().enumerate() {
        let Some(local_routes) = local.shape_routes.take() else {
            continue;
        };
        let archive_ordinal = u32::try_from(archive_ordinal)
            .context("collection-local archive ordinal exceeds u32")?;
        let n_shapes = local_routes.derived.n_shapes;
        if local_routes.derived.archive_ordinal != archive_ordinal {
            bail!("derived shape-route archive ordinal differs from sorted sample order");
        }
        let blocks = local_routes.derived.blocks;
        let descriptors = blocks
            .iter()
            .map(EncodedRouteBlock::descriptor)
            .collect::<Vec<_>>();
        let source_root = digest_from_hex(
            &local.entry.identity.native_digest,
            &format!("sample {} source root", local.entry.id),
        )?;
        let binding = ShapeRouteBinding {
            archive_ordinal,
            source_root,
            shapes_digest: local_routes.shapes_digest,
            n_shapes,
            blocks: descriptors,
        };
        shaperoute::validate_binding(
            &binding,
            archive_ordinal,
            source_root,
            local_routes.shapes_digest,
        )?;
        local.entry.shape_routes = Some(binding);
        encoded_shape_route_blocks.extend(blocks);
    }

    // Flatten, sort, and reduce rather than inserting 1M+ individually allocated BTreeMap nodes.
    // Posts move into the global routes, so the build never retains two copies of the catalogues.
    let local_rows = locals.iter().try_fold(0usize, |sum, local| {
        sum.checked_add(local.junctions.len())
            .context("collection local junction count overflow")
    })?;
    enforce_count(local_rows, MAX_ITEMS, "collection local junction count")?;
    let mut flat: Vec<((u32, u32, u32), ArchiveRoute)> = Vec::with_capacity(local_rows);
    for (archive, local) in locals.iter_mut().enumerate() {
        for junction in std::mem::take(&mut local.junctions) {
            flat.push((
                (junction.chrom, junction.donor, junction.acceptor),
                ArchiveRoute {
                    archive: archive as u32,
                    supporting_children: junction.supporting_children,
                    posts: junction.posts,
                },
            ));
        }
    }
    flat.sort_unstable_by_key(|(key, route)| (*key, route.archive));
    let mut junctions: Vec<GlobalJunction> = Vec::new();
    let mut rows = flat.into_iter().peekable();
    while let Some((key, route)) = rows.next() {
        let mut routes = vec![route];
        while rows.peek().is_some_and(|(next, _)| *next == key) {
            routes.push(rows.next().unwrap().1);
        }
        let mut support_upper_bound = 0u64;
        for route in &routes {
            support_upper_bound = support_upper_bound
                .checked_add(route.supporting_children)
                .context("global junction support upper bound overflow")?;
        }
        junctions.push(GlobalJunction {
            chrom: key.0,
            donor: key.1,
            acceptor: key.2,
            support_upper_bound,
            routes,
        });
    }
    let junction_count = junctions.len();
    let route_count = junctions.iter().try_fold(0usize, |sum, row| {
        sum.checked_add(row.routes.len())
            .context("collection route count overflow")
    })?;
    let posting_count = junctions.iter().try_fold(0usize, |sum, row| {
        row.routes.iter().try_fold(sum, |sum, route| {
            sum.checked_add(route.posts.len())
                .context("collection posting count overflow")
        })
    })?;
    enforce_count(route_count, MAX_ITEMS, "collection route count")?;
    enforce_count(posting_count, MAX_ITEMS, "collection posting count")?;
    Ok(Collection {
        base: None,
        genome_algo,
        genome_digest,
        chroms,
        chroms_digest,
        archives: locals.into_iter().map(|local| local.entry).collect(),
        junctions,
        junction_count,
        route_count,
        posting_count,
        shape_route_blocks: Vec::new(),
        encoded_shape_route_blocks,
    })
}

fn validate_extension(base: &Collection, extension: &Collection) -> Result<()> {
    if base.genome_algo != extension.genome_algo
        || base.genome_digest != extension.genome_digest
        || base.chroms != extension.chroms
        || base.chroms_digest != extension.chroms_digest
    {
        bail!("new collection segment has an incompatible reference identity");
    }
    let base_ids: FxHashSet<&str> = base.archives.iter().map(|row| row.id.as_str()).collect();
    let base_paths: FxHashSet<&Path> = base.archives.iter().map(|row| row.path.as_path()).collect();
    let base_inodes: FxHashSet<(u64, u64)> = base
        .archives
        .iter()
        .filter_map(|row| {
            let identity = &row.identity;
            (identity.dev != 0 || identity.inode != 0).then_some((identity.dev, identity.inode))
        })
        .collect();
    let base_digests: FxHashSet<&str> = base
        .archives
        .iter()
        .map(|row| {
            row.identity
                .encoded_sections_digest
                .as_deref()
                .context("base collection archive lacks an encoded-sections digest")
        })
        .collect::<Result<_>>()?;
    for archive in &extension.archives {
        if base_ids.contains(archive.id.as_str()) {
            bail!(
                "duplicate sample id {} across collection layers",
                archive.id
            );
        }
        if base_paths.contains(archive.path.as_path()) {
            bail!(
                "sample {} reuses resolved archive {} from the base collection",
                archive.id,
                archive.path.display()
            );
        }
        let identity = &archive.identity;
        if (identity.dev != 0 || identity.inode != 0)
            && base_inodes.contains(&(identity.dev, identity.inode))
        {
            bail!(
                "sample {} reuses an archive inode from the base collection",
                archive.id
            );
        }
        let encoded = identity
            .encoded_sections_digest
            .as_deref()
            .context("extension archive lacks an encoded-sections digest")?;
        if base_digests.contains(encoded) {
            bail!(
                "sample {} has byte-identical archive content to a base sample",
                archive.id
            );
        }
    }
    Ok(())
}

fn require_encoded_extension_base(chain: &CollectionChain) -> Result<()> {
    if chain
        .layers
        .iter()
        .any(|layer| layer.file.version == LEGACY_VERSION)
    {
        bail!(
            "collection v2 cannot be used as a v4 extension base because it lacks cross-container encoded identities; rebuild the base"
        );
    }
    Ok(())
}

#[derive(Serialize)]
struct BuildShapeRouteSummary {
    archives: u64,
    sections: u64,
    exact_spans: u64,
    pairs: u64,
    raw_bytes: u64,
    compressed_bytes: u64,
}

#[derive(Serialize)]
struct BuildSourceIoSummary {
    archives: u64,
    identity_content_bytes_read: u64,
    total_bytes_read: u64,
    shape_route_payload_bytes_read: u64,
    shape_route_source_bytes_read: u64,
}

#[derive(Serialize)]
struct BuildUniformSummary<'a> {
    collection_format_version: u32,
    output_path: &'a str,
    output_root: String,
    shape_routes_requested: bool,
    new_archives: u64,
    segment_junctions: u64,
    archive_routes: u64,
    chunk_postings: u64,
    raw_section_bytes: u64,
    file_bytes: u64,
    shape_routes: BuildShapeRouteSummary,
    source_io: BuildSourceIoSummary,
}

fn build_source_io_summary(rows: &[SourceIo]) -> Result<BuildSourceIoSummary> {
    Ok(BuildSourceIoSummary {
        archives: u64::try_from(rows.len()).context("build source count exceeds u64")?,
        identity_content_bytes_read: rows.iter().try_fold(0u64, |sum, row| {
            sum.checked_add(row.identity_content_bytes_read)
                .context("build identity-content byte count overflow")
        })?,
        total_bytes_read: rows.iter().try_fold(0u64, |sum, row| {
            sum.checked_add(row.total_bytes_read)
                .context("build source byte count overflow")
        })?,
        shape_route_payload_bytes_read: rows.iter().try_fold(0u64, |sum, row| {
            sum.checked_add(row.shape_route_payload_bytes_read)
                .context("build shape-route source byte count overflow")
        })?,
        shape_route_source_bytes_read: rows.iter().try_fold(0u64, |sum, row| {
            sum.checked_add(row.shape_route_source_bytes_read)
                .context("build shape-route attributable byte count overflow")
        })?,
    })
}

fn stream_uniform_build(
    writer: &mut dyn Write,
    format: CollectionOutputFormat,
    context: &ResultContext,
    summary: &BuildUniformSummary<'_>,
    collection: &Collection,
    source_io: &[SourceIo],
) -> std::result::Result<(), OutputError> {
    // Validate every fallible row invariant before writing an envelope prefix. Sink failures may
    // still interrupt stdout, but malformed result data cannot produce a partial logical result.
    for archive in &collection.archives {
        if archive.path.to_str().is_none() {
            return Err(OutputError::Sink(format!(
                "archive path is not valid UTF-8: {}",
                archive.path.display()
            )));
        }
    }
    let archive_schema = set_table_schema(
        COLLECTION_BUILD_ARCHIVES_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("path", DataType::String),
            Field::new("archive_format_version", DataType::UInt64),
            Field::new("identity_scheme", DataType::String),
            Field::new("identity_digest", DataType::String),
            Field::new("encoded_sections_digest", DataType::String).nullable(),
            Field::new("chunks", DataType::UInt64),
            Field::new("shape_routes", DataType::Boolean),
        ],
        &["sample"],
    )?;
    let source_io_schema = set_table_schema(
        COLLECTION_BUILD_SOURCE_IO_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("archive_format_version", DataType::UInt64),
            Field::new("identity_scheme", DataType::String),
            Field::new("identity_content_bytes_read", DataType::UInt64),
            Field::new("total_bytes_read", DataType::UInt64),
            Field::new("shape_route_payload_bytes_read", DataType::UInt64),
            Field::new("shape_route_source_bytes_read", DataType::UInt64),
        ],
        &["sample"],
    )?;
    let source_sections_schema = set_table_schema(
        COLLECTION_BUILD_SOURCE_SECTIONS_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("section", DataType::String),
        ],
        &["sample", "section"],
    )?;
    let archive_selection = SelectionSummary::complete(collection.archives.len() as u64);
    let source_selection = SelectionSummary::complete(source_io.len() as u64);
    let source_section_count = source_io.iter().try_fold(0u64, |sum, row| {
        sum.checked_add(row.sections_read.len() as u64)
            .ok_or(OutputError::InvalidSchema(
                "build source-section count exceeds u64".into(),
            ))
    })?;
    let source_section_selection = SelectionSummary::complete(source_section_count);
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        COLLECTION_BUILD_RESULT_SCHEMA,
        OutputFormat::from(format),
        context,
        summary,
    )?;
    bundle.write_table(
        "archives",
        &archive_schema,
        Some(&archive_selection),
        |rows| {
            for archive in &collection.archives {
                let path = archive.path.to_str().ok_or_else(|| {
                    OutputError::Sink(format!(
                        "archive path is not valid UTF-8: {}",
                        archive.path.display()
                    ))
                })?;
                rows.write_row_with(|row| {
                    row.string(&archive.id)?;
                    row.string(path)?;
                    row.uint64(archive.identity.archive_format_version as u64)?;
                    row.string(&archive.identity.native_scheme)?;
                    row.string(&archive.identity.native_digest)?;
                    match archive.identity.encoded_sections_digest.as_deref() {
                        Some(digest) => row.string(digest)?,
                        None => row.null()?,
                    }
                    row.uint64(archive.chunks.len() as u64)?;
                    row.boolean(archive.shape_routes.is_some())?;
                    Ok(())
                })?;
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "source_io",
        &source_io_schema,
        Some(&source_selection),
        |rows| {
            for source in source_io {
                rows.write_row_with(|row| {
                    row.string(&source.id)?;
                    row.uint64(source.format_version as u64)?;
                    row.string(&source.identity_scheme)?;
                    row.uint64(source.identity_content_bytes_read)?;
                    row.uint64(source.total_bytes_read)?;
                    row.uint64(source.shape_route_payload_bytes_read)?;
                    row.uint64(source.shape_route_source_bytes_read)?;
                    Ok(())
                })?;
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "source_sections",
        &source_sections_schema,
        Some(&source_section_selection),
        |rows| {
            for source in source_io {
                for section in &source.sections_read {
                    rows.write_row_with(|row| {
                        row.string(&source.id)?;
                        row.string(section)?;
                        Ok(())
                    })?;
                }
            }
            Ok(())
        },
    )?;
    bundle.finish()?;
    Ok(())
}

struct BuildRun {
    samples: Vec<String>,
    source_digests: Vec<String>,
    base: Option<PathBuf>,
    out: PathBuf,
    allow_unstamped: bool,
    build_shape_routes: bool,
    json_output: bool,
    uniform_output: CollectionOutputArgs,
}

fn run_build(args: BuildRun) -> Result<()> {
    let BuildRun {
        samples,
        source_digests,
        base,
        out,
        allow_unstamped,
        build_shape_routes,
        json_output,
        uniform_output,
    } = args;
    let started = std::time::Instant::now();
    if base.is_some() && samples.is_empty() {
        bail!("an incremental build requires at least one new --sample");
    }
    let base_chain = base.as_deref().map(open_collection_chain).transpose()?;
    let mut source_io = if let Some(chain) = &base_chain {
        require_encoded_extension_base(chain)?;
        validate_collection_sources(&chain.collection, false)?
    } else {
        Vec::new()
    };
    let mut parsed = samples
        .iter()
        .map(|sample| parse_sample(sample))
        .collect::<Result<Vec<_>>>()?;
    parsed.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    let mut digest_of = BTreeMap::new();
    for value in &source_digests {
        let (id, digest) = parse_source_digest(value)?;
        if digest_of.insert(id.clone(), digest).is_some() {
            bail!("duplicate --source-digest for sample {id}");
        }
    }
    let sample_ids: FxHashSet<&str> = parsed.iter().map(|(id, _)| id.as_str()).collect();
    if let Some(extra) = digest_of
        .keys()
        .find(|id| !sample_ids.contains(id.as_str()))
    {
        bail!("--source-digest supplied for unknown new sample {extra}");
    }
    // Derive archives in parallel. Each worker keeps a compact contiguous route table and its
    // canonical encoded blocks; it does not retain a per-span allocation hierarchy.
    let loaded = parsed
        .into_par_iter()
        .enumerate()
        .map(|(archive_ordinal, (id, path))| {
            let digest = digest_of.get(&id).cloned();
            let archive_ordinal = u32::try_from(archive_ordinal)
                .context("collection-local archive ordinal exceeds u32")?;
            load_local_archive(
                id,
                path,
                digest,
                build_shape_routes,
                archive_ordinal,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    source_io.extend(loaded.iter().map(|local| local.source_io.clone()));
    source_io.sort_unstable_by(|left, right| left.id.cmp(&right.id));
    let mut collection = assemble_collection(
        loaded,
        allow_unstamped,
        if base_chain.is_some() { 1 } else { 2 },
    )?;
    if let Some(chain) = &base_chain {
        validate_extension(&chain.collection, &collection)?;
        let root = chain
            .layers
            .last()
            .context("base collection chain is empty")?;
        collection.base = Some(BaseCollection {
            path: root.path.clone(),
            root_digest: root.file.root_digest_hex(),
        });
    }
    let archive_routes: usize = collection
        .junctions
        .iter()
        .map(|row| row.routes.len())
        .sum();
    let posting_routes: usize = collection
        .junctions
        .iter()
        .flat_map(|row| &row.routes)
        .map(|route| route.posts.len())
        .sum();
    let (shape_route_spans, shape_route_pairs): (usize, usize) = if collection
        .encoded_shape_route_blocks
        .is_empty()
    {
        (
            collection
                .shape_route_blocks
                .iter()
                .map(|block| block.spans.len())
                .sum(),
            collection
                .shape_route_blocks
                .iter()
                .flat_map(|block| &block.spans)
                .map(|span| span.pairs.len())
                .sum(),
        )
    } else {
        (
            collection
                .encoded_shape_route_blocks
                .iter()
                .map(|block| block.n_spans)
                .sum(),
            collection
                .encoded_shape_route_blocks
                .iter()
                .map(|block| block.n_pairs)
                .sum(),
        )
    };
    let write_stats = write_collection(&out, &collection)?;
    let elapsed_seconds = started.elapsed().as_secs_f64();
    if let Some(format) = uniform_output.format {
        let reported_out = reported_output_path(&out)?;
        let mut parameters = BTreeMap::new();
        parameters.insert("output_path".into(), json!(&reported_out));
        parameters.insert("allow_unstamped".into(), json!(allow_unstamped));
        parameters.insert("shape_routes".into(), json!(build_shape_routes));
        parameters.insert(
            "source_digest_assertions".into(),
            serde_json::to_value(&digest_of)?,
        );
        parameters.insert(
            "reference".into(),
            match (
                collection.genome_algo.as_deref(),
                collection.genome_digest.as_deref(),
            ) {
                (Some(algo), Some(digest)) => {
                    json!({"stamped": true, "algo": algo, "digest": digest})
                }
                _ => json!({"stamped": false}),
            },
        );
        parameters.insert(
            "chromosome_dictionary_digest".into(),
            json!(&collection.chroms_digest),
        );
        parameters.insert(
            "base_layers".into(),
            serde_json::Value::Array(
                base_chain
                    .as_ref()
                    .map(|chain| {
                        chain
                            .layers
                            .iter()
                            .map(|layer| -> Result<_> {
                                Ok(json!({
                                    "path": uniform_path(&layer.path, "base collection layer")?,
                                    "format_version": layer.file.version,
                                    "root": format!("aicollection-directory-root-v1:{}", layer.file.root_digest_hex()),
                                }))
                            })
                            .collect::<Result<Vec<_>>>()
                    })
                    .transpose()?
                    .unwrap_or_default(),
            ),
        );
        let mut source_archives = base_chain
            .as_ref()
            .map(|chain| chain.collection.archives.iter().collect::<Vec<_>>())
            .unwrap_or_default();
        source_archives.extend(collection.archives.iter());
        let context = uniform_source_context(source_archives, parameters)?;
        let summary = BuildUniformSummary {
            collection_format_version: VERSION,
            output_path: &reported_out,
            output_root: format!(
                "aicollection-directory-root-v1:{}",
                digest_hex(write_stats.root_digest)
            ),
            shape_routes_requested: build_shape_routes,
            new_archives: collection.archives.len() as u64,
            segment_junctions: collection.junctions.len() as u64,
            archive_routes: archive_routes as u64,
            chunk_postings: posting_routes as u64,
            raw_section_bytes: write_stats.raw_bytes,
            file_bytes: write_stats.file_bytes,
            shape_routes: BuildShapeRouteSummary {
                archives: collection
                    .archives
                    .iter()
                    .filter(|archive| archive.shape_routes.is_some())
                    .count() as u64,
                sections: write_stats.shape_route_sections as u64,
                exact_spans: shape_route_spans as u64,
                pairs: shape_route_pairs as u64,
                raw_bytes: write_stats.shape_route_raw_bytes,
                compressed_bytes: write_stats.shape_route_compressed_bytes,
            },
            source_io: build_source_io_summary(&source_io)?,
        };
        write_uniform_collection_result(&uniform_output, |writer| {
            stream_uniform_build(writer, format, &context, &summary, &collection, &source_io)
        })?;
        return Ok(());
    }
    if json_output {
        let identity_content_bytes_read = source_io.iter().try_fold(0u64, |sum, row| {
            sum.checked_add(row.identity_content_bytes_read)
                .context("build identity-content byte count overflow")
        })?;
        let total_bytes_read = source_io.iter().try_fold(0u64, |sum, row| {
            sum.checked_add(row.total_bytes_read)
                .context("build source byte count overflow")
        })?;
        let shape_route_payload_bytes_read = source_io.iter().try_fold(0u64, |sum, row| {
            sum.checked_add(row.shape_route_payload_bytes_read)
                .context("build shape-route source byte count overflow")
        })?;
        let shape_route_source_bytes_read = source_io.iter().try_fold(0u64, |sum, row| {
            sum.checked_add(row.shape_route_source_bytes_read)
                .context("build shape-route attributable byte count overflow")
        })?;
        let mut sections_read: Vec<String> = source_io
            .iter()
            .flat_map(|row| row.sections_read.iter().cloned())
            .collect();
        sections_read.sort_unstable();
        sections_read.dedup();
        let value = json!({
            "schema": "gravlax.collection.build.v2",
            "output": out,
            "collection_format_version": VERSION,
            "shape_routes_requested": build_shape_routes,
            "new_archives": collection.archives.len(),
            "segment_junctions": collection.junctions.len(),
            "archive_routes": archive_routes,
            "chunk_postings": posting_routes,
            "raw_section_bytes": write_stats.raw_bytes,
            "file_bytes": write_stats.file_bytes,
            "shape_routes": {
                "archives": collection.archives.iter().filter(|archive| archive.shape_routes.is_some()).count(),
                "sections": write_stats.shape_route_sections,
                "exact_spans": shape_route_spans,
                "pairs": shape_route_pairs,
                "raw_bytes": write_stats.shape_route_raw_bytes,
                "compressed_bytes": write_stats.shape_route_compressed_bytes,
            },
            "elapsed_seconds": elapsed_seconds,
            "source_io": {
                "identity_content_bytes_read": identity_content_bytes_read,
                "total_bytes_read": total_bytes_read,
                "shape_route_payload_bytes_read": shape_route_payload_bytes_read,
                "shape_route_source_bytes_read": shape_route_source_bytes_read,
                "sections_read": sections_read,
                "archives": source_io.iter().map(|row| json!({
                    "id": row.id,
                    "format_version": row.format_version,
                    "identity_scheme": row.identity_scheme,
                    "identity_content_bytes_read": row.identity_content_bytes_read,
                    "total_bytes_read": row.total_bytes_read,
                    "shape_route_payload_bytes_read": row.shape_route_payload_bytes_read,
                    "shape_route_source_bytes_read": row.shape_route_source_bytes_read,
                    "sections_read": row.sections_read,
                })).collect::<Vec<_>>(),
            },
        });
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        println!(
            "built {}: {} new archives, {} segment junctions, {} archive routes, {} chunk postings; {} raw / {} file bytes in {:.2}s",
            out.display(),
            collection.archives.len(),
            collection.junctions.len(),
            archive_routes,
            posting_routes,
            write_stats.raw_bytes,
            write_stats.file_bytes,
            elapsed_seconds
        );
    }
    Ok(())
}

fn validate_collection_sources(
    collection: &Collection,
    verify_content: bool,
) -> Result<Vec<SourceIo>> {
    collection
        .archives
        .par_iter()
        .map(|archive| validate_archive_identity(archive, verify_content))
        .collect()
}

fn identity_content_bytes(rows: &[SourceIo]) -> Result<u64> {
    rows.iter().try_fold(0u64, |sum, row| {
        sum.checked_add(row.identity_content_bytes_read)
            .context("identity-verification byte count overflow")
    })
}

fn source_archive_bytes_read(rows: &[SourceIo]) -> Result<u64> {
    rows.iter().try_fold(0u64, |sum, row| {
        sum.checked_add(row.total_bytes_read)
            .context("source-archive byte count overflow")
    })
}

fn chunk_infos(archive: &ArchiveEntry) -> Vec<ChunkInfo> {
    archive
        .chunks
        .iter()
        .map(|chunk| chunk.info.clone())
        .collect()
}

fn open_source(
    collection: &Collection,
    archive_index: usize,
    shape_route: Option<&ShapeRouteBinding>,
) -> Result<LazyArchive> {
    let entry = &collection.archives[archive_index];
    let mut archive = LazyArchive::open(&entry.path).with_context(|| {
        format!(
            "opening routed sample {} at {}",
            entry.id,
            entry.path.display()
        )
    })?;
    let observed = identity_from_metadata(
        archive.reader().file_metadata()?,
        0,
        String::new(),
        String::new(),
        None,
    )?;
    if !identity_stat_matches(&entry.identity, &observed) {
        bail!(
            "archive identity changed for sample {} at {}; rebuild the collection index",
            entry.id,
            entry.path.display()
        );
    }
    let reader = archive.reader();
    if reader.archive_version() != entry.identity.archive_format_version {
        bail!(
            "archive format version changed for sample {} at {}; rebuild the collection index",
            entry.id,
            entry.path.display()
        );
    }
    if reader.archive_version() == evidence_io::format::VERSION {
        let root = reader
            .content_commitment()
            .context("rooted archive lacks its directory commitment")?;
        let encoded = reader
            .encoded_content_identity()?
            .context("rooted archive lacks its encoded-sections identity")?;
        let encoded = digest_hex(encoded);
        if entry.identity.native_scheme != ROOTED_DIRECTORY_SCHEME
            || root.to_hex() != entry.identity.native_digest
            || entry.identity.encoded_sections_digest.as_deref() != Some(encoded.as_str())
        {
            bail!(
                "archive root identity changed for sample {} at {}; rebuild the collection index",
                entry.id,
                entry.path.display()
            );
        }
        if let Some(binding) = shape_route {
            let shapes_digest = reader
                .section_metadata()
                .find(|section| section.name == "shapes")
                .context("routed archive lacks its shapes section")?
                .compressed_blake3
                .context("routed archive shapes entry lacks its committed payload digest")?;
            shaperoute::validate_binding(
                binding,
                binding.archive_ordinal,
                root.digest,
                shapes_digest,
            )
            .with_context(|| format!("validating shape-route source for sample {}", entry.id))?;
        }
    } else if entry.identity.native_scheme != LEGACY_FULL_FILE_SCHEME {
        bail!(
            "legacy archive identity scheme changed for sample {}; rebuild the collection index",
            entry.id
        );
    }
    if shape_route.is_some() && reader.archive_version() != evidence_io::format::VERSION {
        bail!("shape routes cannot be used with legacy archive sample {}", entry.id);
    }
    if archive.chrom_names != collection.chroms || archive.chrom_digest != collection.chroms_digest
    {
        bail!(
            "sample {} chromosome dictionary bytes no longer match the collection",
            entry.id
        );
    }
    let observed = archive
        .genome_sig
        .as_ref()
        .map(|sig| (&sig.algo, &sig.digest));
    let expected = collection
        .genome_algo
        .as_ref()
        .zip(collection.genome_digest.as_ref());
    if observed != expected {
        bail!(
            "sample {} stamped genome identity no longer matches the collection",
            entry.id
        );
    }
    Ok(archive)
}

fn chrom_id(collection: &Collection, chrom: &str) -> Result<u32> {
    collection
        .chroms
        .iter()
        .position(|name| name == chrom)
        .map(|index| index as u32)
        .with_context(|| format!("unknown chromosome {chrom}"))
}

fn unpack_barcode(packed: u32) -> String {
    String::from_utf8(evidence_io::umi::unpack(packed, 16))
        .expect("packed barcode always decodes as ASCII")
}

type CountSummary = (usize, usize, Vec<(String, usize)>);

fn summarize_counts(
    archive: &mut LazyArchive,
    per_cell: &FxHashMap<u32, FxHashSet<u32>>,
    top: usize,
) -> Result<CountSummary> {
    let mut cells: Vec<(u32, usize)> = per_cell
        .iter()
        .map(|(cell, classes)| (*cell, classes.len()))
        .collect();
    cells.sort_unstable_by_key(|(cell, count)| (std::cmp::Reverse(*count), *cell));
    let umis = cells.iter().map(|(_, count)| count).sum();
    let dictionary = archive.cells()?;
    let limit = if top == 0 {
        cells.len()
    } else {
        top.min(cells.len())
    };
    let rows = cells
        .into_iter()
        .take(limit)
        .map(|(cell, count)| {
            let packed = dictionary
                .get(cell as usize)
                .copied()
                .with_context(|| format!("cell id {cell} is outside the barcode dictionary"))?;
            Ok((unpack_barcode(packed), count))
        })
        .collect::<Result<_>>()?;
    Ok((umis, per_cell.len(), rows))
}

#[derive(Debug)]
struct SampleCounts {
    archive: usize,
    present: bool,
    supporting_children: u64,
    chunks: usize,
    planned_bytes: u64,
    actual_bytes: u64,
    shape_route_used: bool,
    shape_route_section: Option<String>,
    shape_route_sidecar_bytes: u64,
    molecules: Option<u64>,
    umis: usize,
    cells: usize,
    top_cells: Vec<(String, usize)>,
}

fn zero_counts(archive: usize, present: bool) -> SampleCounts {
    SampleCounts {
        archive,
        present,
        supporting_children: 0,
        chunks: 0,
        planned_bytes: 0,
        actual_bytes: 0,
        shape_route_used: false,
        shape_route_section: None,
        shape_route_sidecar_bytes: 0,
        molecules: None,
        umis: 0,
        cells: 0,
        top_cells: Vec::new(),
    }
}

#[derive(Serialize)]
struct CollectionSelectionPolicy {
    requested_top_per_sample: usize,
    top_zero_means_all: bool,
    comparator: &'static str,
}

#[derive(Serialize)]
struct CollectionQueryPlanningSummary {
    collection_layers: u64,
    archives_total: u64,
    archives_opened: u64,
    archives_pruned: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    support_bound_pruned: Option<bool>,
    unique_chunks_decoded: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    independent_chunk_decodes: Option<u64>,
    planned_compressed_bytes: u64,
    actual_archive_bytes_read: u64,
    source_archive_identity_bytes_read: u64,
    source_archive_execution_bytes_read: u64,
    source_archive_bytes_read: u64,
    collection_sidecar_bytes_read: u64,
    shape_route_sidecar_payload_bytes_read: u64,
    total_logical_bytes_read: u64,
    route_blocks_loaded: u64,
    routed_archives: u64,
    fallback_archives: u64,
}

#[derive(Serialize)]
struct CollectionJunctionUniformSummary<'a> {
    coordinates: &'static str,
    chrom: &'a str,
    donor: u32,
    acceptor: u32,
    support_upper_bound: u64,
    min_support: u64,
    umis: u64,
    cells: u64,
    selection: CollectionSelectionPolicy,
    planning: CollectionQueryPlanningSummary,
}

#[derive(Serialize)]
struct CollectionRegionUniformSummary<'a> {
    coordinates: &'static str,
    anchor_semantics: bool,
    chrom: &'a str,
    start: u32,
    end: u32,
    molecules: u64,
    umis: u64,
    cells: u64,
    selection: CollectionSelectionPolicy,
    planning: CollectionQueryPlanningSummary,
}

fn collection_query_parameters(
    locus: &str,
    top: usize,
    explain: bool,
    verify_content: bool,
    archive_access: &'static str,
) -> BTreeMap<String, serde_json::Value> {
    let mut parameters = BTreeMap::new();
    parameters.insert("locus".into(), json!(locus));
    parameters.insert("top".into(), json!(top));
    parameters.insert("explain".into(), json!(explain));
    parameters.insert("verify_content".into(), json!(verify_content));
    parameters.insert("archive_access".into(), json!(archive_access));
    parameters
}

fn collection_count_cell_selection(
    results: &[SampleCounts],
) -> std::result::Result<SelectionSummary, OutputError> {
    let available = results.iter().try_fold(0u64, |sum, result| {
        sum.checked_add(result.cells as u64)
            .ok_or(OutputError::InvalidSchema(
                "collection cell availability exceeds u64".into(),
            ))
    })?;
    let emitted = results.iter().try_fold(0u64, |sum, result| {
        sum.checked_add(result.top_cells.len() as u64)
            .ok_or(OutputError::InvalidSchema(
                "collection emitted cell count exceeds u64".into(),
            ))
    })?;
    SelectionSummary::selected(available, emitted)
}

fn junction_sample_schema() -> std::result::Result<TableSchema, OutputError> {
    set_table_schema(
        COLLECTION_JUNCTION_SAMPLES_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("present", DataType::Boolean),
            Field::new("supporting_children", DataType::UInt64),
            Field::new("umis", DataType::UInt64),
            Field::new("cells", DataType::UInt64),
            Field::new("chunks_decoded", DataType::UInt64),
            Field::new("planned_compressed_bytes", DataType::UInt64),
            Field::new("actual_archive_bytes_read", DataType::UInt64),
            Field::new("shape_route_used", DataType::Boolean),
            Field::new("shape_route_section", DataType::String).nullable(),
            Field::new("shape_route_sidecar_bytes_read", DataType::UInt64),
            Field::new("decision", DataType::String).nullable(),
        ],
        &["sample"],
    )
}

fn region_sample_schema() -> std::result::Result<TableSchema, OutputError> {
    set_table_schema(
        COLLECTION_REGION_SAMPLES_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("present", DataType::Boolean),
            Field::new("molecules", DataType::UInt64),
            Field::new("umis", DataType::UInt64),
            Field::new("cells", DataType::UInt64),
            Field::new("chunks_decoded", DataType::UInt64),
            Field::new("planned_compressed_bytes", DataType::UInt64),
            Field::new("actual_archive_bytes_read", DataType::UInt64),
            Field::new("decision", DataType::String).nullable(),
        ],
        &["sample"],
    )
}

fn collection_cell_schema(id: &'static str) -> std::result::Result<TableSchema, OutputError> {
    set_table_schema(
        id,
        vec![
            Field::new("sample", DataType::String),
            Field::new("rank", DataType::UInt64),
            Field::new("barcode", DataType::String),
            Field::new("umis", DataType::UInt64),
        ],
        &["sample", "barcode"],
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "the streaming serializer keeps each schema-bound input explicit at the output boundary"
)]
fn stream_uniform_junction(
    writer: &mut dyn Write,
    format: CollectionOutputFormat,
    context: &ResultContext,
    summary: &CollectionJunctionUniformSummary<'_>,
    collection: &Collection,
    results: &[SampleCounts],
    support_bound_pruned: bool,
    explain: bool,
) -> std::result::Result<(), OutputError> {
    let sample_schema = junction_sample_schema()?;
    let cell_schema = collection_cell_schema(COLLECTION_JUNCTION_CELLS_SCHEMA)?;
    let sample_selection = SelectionSummary::complete(results.len() as u64);
    let cell_selection = collection_count_cell_selection(results)?;
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        COLLECTION_JUNCTION_RESULT_SCHEMA,
        OutputFormat::from(format),
        context,
        summary,
    )?;
    bundle.write_table("samples", &sample_schema, Some(&sample_selection), |rows| {
        for result in results {
            rows.write_row_with(|row| {
                row.string(&collection.archives[result.archive].id)?;
                row.boolean(result.present)?;
                row.uint64(result.supporting_children)?;
                row.uint64(result.umis as u64)?;
                row.uint64(result.cells as u64)?;
                row.uint64(result.chunks as u64)?;
                row.uint64(result.planned_bytes)?;
                row.uint64(result.actual_bytes)?;
                row.boolean(result.shape_route_used)?;
                match result.shape_route_section.as_deref() {
                    Some(section) => row.string(section)?,
                    None => row.null()?,
                }
                row.uint64(result.shape_route_sidecar_bytes)?;
                if explain {
                    row.string(if support_bound_pruned && result.present {
                        "prune_support_bound"
                    } else if result.present {
                        "open"
                    } else {
                        "prune_absent"
                    })?;
                } else {
                    row.null()?;
                }
                Ok(())
            })?;
        }
        Ok(())
    })?;
    bundle.write_table("cells", &cell_schema, Some(&cell_selection), |rows| {
        for result in results {
            let sample = &collection.archives[result.archive].id;
            for (rank, (barcode, umis)) in result.top_cells.iter().enumerate() {
                rows.write_row_with(|row| {
                    row.string(sample)?;
                    row.uint64(rank as u64 + 1)?;
                    row.string(barcode)?;
                    row.uint64(*umis as u64)?;
                    Ok(())
                })?;
            }
        }
        Ok(())
    })?;
    bundle.finish()?;
    Ok(())
}

fn stream_uniform_region(
    writer: &mut dyn Write,
    format: CollectionOutputFormat,
    context: &ResultContext,
    summary: &CollectionRegionUniformSummary<'_>,
    collection: &Collection,
    results: &[SampleCounts],
    explain: bool,
) -> std::result::Result<(), OutputError> {
    let sample_schema = region_sample_schema()?;
    let cell_schema = collection_cell_schema(COLLECTION_REGION_CELLS_SCHEMA)?;
    let sample_selection = SelectionSummary::complete(results.len() as u64);
    let cell_selection = collection_count_cell_selection(results)?;
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        COLLECTION_REGION_RESULT_SCHEMA,
        OutputFormat::from(format),
        context,
        summary,
    )?;
    bundle.write_table("samples", &sample_schema, Some(&sample_selection), |rows| {
        for result in results {
            rows.write_row_with(|row| {
                row.string(&collection.archives[result.archive].id)?;
                row.boolean(result.present)?;
                row.uint64(result.molecules.unwrap_or(0))?;
                row.uint64(result.umis as u64)?;
                row.uint64(result.cells as u64)?;
                row.uint64(result.chunks as u64)?;
                row.uint64(result.planned_bytes)?;
                row.uint64(result.actual_bytes)?;
                if explain {
                    row.string(if result.present {
                        "open"
                    } else {
                        "prune_no_overlapping_chunks"
                    })?;
                } else {
                    row.null()?;
                }
                Ok(())
            })?;
        }
        Ok(())
    })?;
    bundle.write_table("cells", &cell_schema, Some(&cell_selection), |rows| {
        for result in results {
            let sample = &collection.archives[result.archive].id;
            for (rank, (barcode, umis)) in result.top_cells.iter().enumerate() {
                rows.write_row_with(|row| {
                    row.string(sample)?;
                    row.uint64(rank as u64 + 1)?;
                    row.string(barcode)?;
                    row.uint64(*umis as u64)?;
                    Ok(())
                })?;
            }
        }
        Ok(())
    })?;
    bundle.finish()?;
    Ok(())
}

fn sample_counts_json(collection: &Collection, result: &SampleCounts) -> serde_json::Value {
    let archive = &collection.archives[result.archive];
    let mut value = json!({
        "sample": archive.id,
        "archive": archive.path,
        "present": result.present,
        "supporting_children": result.supporting_children,
        "chunks_decoded": result.chunks,
        "planned_compressed_bytes": result.planned_bytes,
        "actual_archive_bytes_read": result.actual_bytes,
        "shape_route_used": result.shape_route_used,
        "shape_route_section": result.shape_route_section,
        "shape_route_sidecar_bytes_read": result.shape_route_sidecar_bytes,
        "umis": result.umis,
        "cells": result.cells,
        "top_cells": result.top_cells.iter().map(|(barcode, umis)| json!({
            "barcode": barcode,
            "umis": umis,
        })).collect::<Vec<_>>(),
    });
    if let Some(molecules) = result.molecules {
        value
            .as_object_mut()
            .unwrap()
            .insert("molecules".into(), json!(molecules));
    }
    value
}

fn planned_bytes(entry: &ArchiveEntry, posts: &[u32]) -> Result<u64> {
    posts.iter().try_fold(0u64, |sum, &post| {
        let bytes = entry
            .chunks
            .get(post as usize)
            .with_context(|| format!("sample {} route references missing chunk {post}", entry.id))?
            .compressed_bytes;
        sum.checked_add(bytes)
            .context("planned compressed byte count overflow")
    })
}

fn region_counts_routed(
    archive: &mut LazyArchive,
    chunks: &[ChunkInfo],
    selected: &[u32],
    start: u32,
    end: u32,
) -> Result<Counts> {
    if let Some(&bad) = selected
        .iter()
        .find(|&&chunk| chunk as usize >= chunks.len())
    {
        bail!("region route references missing chunk {bad}");
    }
    let mut per_cell: FxHashMap<u32, FxHashSet<u32>> = FxHashMap::default();
    let mut molecules = 0u64;
    for chunk in selected {
        let decoded = {
            let (reader, tables) = archive.reader_and_tables();
            let (compressed, raw_len) = reader.read_compressed_at(&format!("c{chunk}"))?;
            let raw = evidence_io::format::decompress(&compressed, raw_len)?;
            decode_chunk(&raw, &chunks[*chunk as usize], None, tables)?
        };
        let classes: Vec<u32> = decoded
            .iter()
            .filter(|molecule| molecule.anchor() >= start && molecule.anchor() < end)
            .map(|molecule| molecule.umi_class)
            .collect();
        molecules = molecules
            .checked_add(classes.len() as u64)
            .context("region molecule count overflow")?;
        archive.prefetch_coc(classes.iter().copied())?;
        for class in classes {
            per_cell
                .entry(archive.cell_of(class)?)
                .or_default()
                .insert(class);
        }
    }
    Ok(Counts {
        molecules: Some(molecules),
        per_cell,
    })
}

struct JunctionQueryRun {
    path: PathBuf,
    locus: String,
    min_support: u64,
    top: usize,
    explain: bool,
    verify_content: bool,
    json_output: bool,
    uniform_output: CollectionOutputArgs,
}

fn run_junction_query(args: JunctionQueryRun) -> Result<()> {
    let JunctionQueryRun {
        path,
        locus,
        min_support,
        top,
        explain,
        verify_content,
        json_output,
        uniform_output,
    } = args;
    let started = std::time::Instant::now();
    let chain = open_collection_chain(&path)?;
    let collection = &chain.collection;
    let load_seconds = started.elapsed().as_secs_f64();
    let uniform_context = if uniform_output.format.is_some() {
        let mut parameters = collection_query_parameters(
            &locus,
            top,
            explain,
            verify_content,
            "collection catalogue/postings and optional source-root-bound local-shape routes",
        );
        parameters.insert("min_support".into(), json!(min_support));
        Some(collection_uniform_context(&path, &chain, parameters)?)
    } else {
        None
    };
    let identity_started = std::time::Instant::now();
    let identity_io = validate_collection_sources(collection, verify_content)?;
    let identity_seconds = identity_started.elapsed().as_secs_f64();
    let identity_bytes_read = identity_content_bytes(&identity_io)?;
    let identity_source_bytes = source_archive_bytes_read(&identity_io)?;
    let planning_started = std::time::Instant::now();
    let (chrom, donor, acceptor) = parse_locus(&locus)?;
    let chrom_id = chrom_id(collection, &chrom)?;
    let junction = lookup_chain_junction(&chain, chrom_id, donor, acceptor)?;
    let support_bound_pruned = junction
        .as_ref()
        .is_some_and(|row| row.support_upper_bound < min_support);
    let mut route_of: FxHashMap<usize, &ArchiveRoute> = FxHashMap::default();
    if !support_bound_pruned {
        if let Some(junction) = &junction {
            route_of.extend(
                junction
                    .routes
                    .iter()
                    .map(|route| (route.archive as usize, route)),
            );
        }
    }
    let mut routed: Vec<usize> = route_of.keys().copied().collect();
    routed.sort_unstable();
    let exact_span = acceptor
        .checked_sub(donor)
        .context("junction acceptor precedes donor")?;
    let mut shape_route_of = FxHashMap::default();
    for &archive_index in &routed {
        if let Some(route) = lookup_chain_shape_routes(&chain, archive_index, &[exact_span])? {
            shape_route_of.insert(archive_index, route);
        }
    }
    let planning_seconds = planning_started.elapsed().as_secs_f64();
    let execution_started = std::time::Instant::now();
    let mut computed: Vec<SampleCounts> = routed
        .par_iter()
        .map(|&archive_index| -> Result<SampleCounts> {
            let route = route_of[&archive_index];
            let entry = &collection.archives[archive_index];
            let bytes = planned_bytes(entry, &route.posts)?;
            let chunks = chunk_infos(entry);
            let shape_route = shape_route_of.get(&archive_index);
            let mut archive = open_source(
                collection,
                archive_index,
                shape_route.map(|route| &route.binding),
            )?;
            let (_, n_chunks, per_cell) = junction_counts_routed_with_shape_route(
                &mut archive,
                &chunks,
                chrom_id,
                donor,
                acceptor,
                route.supporting_children,
                &route.posts,
                shape_route.map(|route| {
                    (
                        route
                            .spans
                            .get(&exact_span)
                            .expect("planned exact span is present"),
                        route.binding.n_shapes,
                    )
                }),
            )?;
            let (umis, cells, top_cells) = summarize_counts(&mut archive, &per_cell, top)?;
            let actual_bytes = archive.reader().bytes_read();
            Ok(SampleCounts {
                archive: archive_index,
                present: true,
                supporting_children: route.supporting_children,
                chunks: n_chunks,
                planned_bytes: bytes,
                actual_bytes,
                shape_route_used: shape_route.is_some(),
                shape_route_section: shape_route
                    .and_then(|route| route.section_names.first().cloned()),
                shape_route_sidecar_bytes: shape_route
                    .map_or(0, |route| route.sidecar_compressed_bytes),
                molecules: None,
                umis,
                cells,
                top_cells,
            })
        })
        .collect::<Result<_>>()?;
    let execution_seconds = execution_started.elapsed().as_secs_f64();
    let mut results: Vec<SampleCounts> = (0..collection.archives.len())
        .map(|archive| {
            zero_counts(
                archive,
                junction.as_ref().is_some_and(|row| {
                    row.routes
                        .iter()
                        .any(|route| route.archive as usize == archive)
                }),
            )
        })
        .collect();
    for result in computed.drain(..) {
        let archive = result.archive;
        results[archive] = result;
    }
    let archives_opened = routed.len();
    let chunks_decoded: usize = results.iter().map(|result| result.chunks).sum();
    let bytes: u64 = results.iter().map(|result| result.planned_bytes).sum();
    let actual_bytes: u64 = results.iter().map(|result| result.actual_bytes).sum();
    let route_sidecar_bytes: u64 = results
        .iter()
        .map(|result| result.shape_route_sidecar_bytes)
        .sum();
    let collection_sidecar_bytes = collection_sidecar_bytes_read(&chain)?;
    let source_archive_bytes = identity_source_bytes
        .checked_add(actual_bytes)
        .context("junction source-archive byte count overflow")?;
    let total_logical_bytes = collection_sidecar_bytes
        .checked_add(source_archive_bytes)
        .context("junction total logical byte count overflow")?;
    let total_umis: usize = results.iter().map(|result| result.umis).sum();
    let total_cells: usize = results.iter().map(|result| result.cells).sum();
    let total_seconds = started.elapsed().as_secs_f64();
    if let Some(format) = uniform_output.format {
        let summary = CollectionJunctionUniformSummary {
            coordinates: "0-based junction boundaries",
            chrom: &chrom,
            donor,
            acceptor,
            support_upper_bound: junction.as_ref().map_or(0, |row| row.support_upper_bound),
            min_support,
            umis: total_umis as u64,
            cells: total_cells as u64,
            selection: CollectionSelectionPolicy {
                requested_top_per_sample: top,
                top_zero_means_all: true,
                comparator: "UMIs descending, archive cell id ascending",
            },
            planning: CollectionQueryPlanningSummary {
                collection_layers: chain.layers.len() as u64,
                archives_total: collection.archives.len() as u64,
                archives_opened: archives_opened as u64,
                archives_pruned: (collection.archives.len() - archives_opened) as u64,
                support_bound_pruned: Some(support_bound_pruned),
                unique_chunks_decoded: chunks_decoded as u64,
                independent_chunk_decodes: None,
                planned_compressed_bytes: bytes,
                actual_archive_bytes_read: actual_bytes,
                source_archive_identity_bytes_read: identity_source_bytes,
                source_archive_execution_bytes_read: actual_bytes,
                source_archive_bytes_read: source_archive_bytes,
                collection_sidecar_bytes_read: collection_sidecar_bytes,
                shape_route_sidecar_payload_bytes_read: route_sidecar_bytes,
                total_logical_bytes_read: total_logical_bytes,
                route_blocks_loaded: shape_route_of.len() as u64,
                routed_archives: shape_route_of.len() as u64,
                fallback_archives: (archives_opened - shape_route_of.len()) as u64,
            },
        };
        let context = uniform_context
            .as_ref()
            .expect("uniform collection junction context was prepared");
        write_uniform_collection_result(&uniform_output, |writer| {
            stream_uniform_junction(
                writer,
                format,
                context,
                &summary,
                collection,
                &results,
                support_bound_pruned,
                explain,
            )
        })?;
        return Ok(());
    }
    let mut value = json!({
        "schema": "gravlax.collection.junction.v2",
        "collection_schema": SCHEMA,
        "coordinates": "0-based junction boundaries",
        "chrom": chrom,
        "donor": donor,
        "acceptor": acceptor,
        "support_upper_bound": junction.as_ref().map_or(0, |row| row.support_upper_bound),
        "min_support": min_support,
        "totals": {"umis": total_umis, "cells": total_cells},
        "samples": results.iter().map(|result| sample_counts_json(collection, result)).collect::<Vec<_>>(),
        "planning": {
            "collection_layers": chain.layers.len(),
            "archives_total": collection.archives.len(),
            "archives_opened": archives_opened,
            "archives_pruned": collection.archives.len() - archives_opened,
            "support_bound_pruned": support_bound_pruned,
            "unique_chunks_decoded": chunks_decoded,
            "planned_compressed_bytes": bytes,
            "actual_archive_bytes_read": actual_bytes,
            "source_archive_identity_bytes_read": identity_source_bytes,
            "source_archive_execution_bytes_read": actual_bytes,
            "source_archive_bytes_read": source_archive_bytes,
            "collection_sidecar_bytes_read": collection_sidecar_bytes,
            "shape_route_sidecar_payload_bytes_read": route_sidecar_bytes,
            "total_logical_bytes_read": total_logical_bytes,
            "route_blocks_loaded": shape_route_of.len(),
            "routed_archives": shape_route_of.len(),
            "fallback_archives": archives_opened - shape_route_of.len(),
            "archive_catalogue_sections_read": 0,
            "archive_posting_sections_read": 0,
            "collection_load_seconds": load_seconds,
            "identity_guard_seconds": identity_seconds,
            "identity_content_bytes_read": identity_bytes_read,
            "route_planning_seconds": planning_seconds,
            "source_execution_seconds": execution_seconds,
            "total_seconds": total_seconds,
        },
    });
    if explain {
        value.as_object_mut().unwrap().insert(
            "explain".into(),
            json!(results
                .iter()
                .map(|result| json!({
                    "sample": collection.archives[result.archive].id,
                    "decision": if support_bound_pruned && result.present {
                        "prune_support_bound"
                    } else if result.present {
                        "open"
                    } else {
                        "prune_absent"
                    },
                    "posting_chunks": result.chunks,
                    "planned_compressed_bytes": result.planned_bytes,
                    "actual_archive_bytes_read": result.actual_bytes,
                }))
                .collect::<Vec<_>>()),
        );
    }
    if json_output {
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        for result in &results {
            println!(
                "{}\t{}\t{}\t{}",
                collection.archives[result.archive].id,
                result.umis,
                result.cells,
                if result.present { "routed" } else { "absent" }
            );
        }
        println!(
            "collection junction {locus}: {total_umis} UMIs / {total_cells} cells; {archives_opened}/{} archives, {chunks_decoded} chunks, {bytes} compressed bytes ({total_seconds:.3}s)",
            collection.archives.len()
        );
        if explain {
            eprintln!("{}", serde_json::to_string_pretty(&value["planning"])?);
        }
    }
    Ok(())
}

fn run_region_query(
    path: PathBuf,
    locus: String,
    top: usize,
    explain: bool,
    verify_content: bool,
    json_output: bool,
    uniform_output: CollectionOutputArgs,
) -> Result<()> {
    let started = std::time::Instant::now();
    let chain = open_collection_chain(&path)?;
    let collection = &chain.collection;
    let load_seconds = started.elapsed().as_secs_f64();
    let uniform_context = if uniform_output.format.is_some() {
        Some(collection_uniform_context(
            &path,
            &chain,
            collection_query_parameters(
                &locus,
                top,
                explain,
                verify_content,
                "collection interval routing and selected archive chunks",
            ),
        )?)
    } else {
        None
    };
    let identity_started = std::time::Instant::now();
    let identity_io = validate_collection_sources(collection, verify_content)?;
    let identity_seconds = identity_started.elapsed().as_secs_f64();
    let identity_bytes_read = identity_content_bytes(&identity_io)?;
    let identity_source_bytes = source_archive_bytes_read(&identity_io)?;
    let planning_started = std::time::Instant::now();
    let (chrom, start, end) = parse_locus(&locus)?;
    let chrom_id = chrom_id(collection, &chrom)?;
    let routes: Vec<Vec<u32>> = collection
        .archives
        .iter()
        .map(|archive| {
            archive
                .chunks
                .iter()
                .enumerate()
                .filter(|(_, chunk)| region_selects_chunk(&chunk.info, chrom_id, start, end))
                .map(|(index, _)| index as u32)
                .collect()
        })
        .collect();
    let routed: Vec<usize> = routes
        .iter()
        .enumerate()
        .filter_map(|(archive, chunks)| (!chunks.is_empty()).then_some(archive))
        .collect();
    let planning_seconds = planning_started.elapsed().as_secs_f64();
    let execution_started = std::time::Instant::now();
    let computed: Vec<SampleCounts> = routed
        .iter()
        .map(|&archive_index| -> Result<SampleCounts> {
            let entry = &collection.archives[archive_index];
            let bytes = planned_bytes(entry, &routes[archive_index])?;
            let chunks = chunk_infos(entry);
            let mut archive = open_source(collection, archive_index, None)?;
            let counts =
                region_counts_routed(&mut archive, &chunks, &routes[archive_index], start, end)?;
            let (umis, cells, top_cells) = summarize_counts(&mut archive, &counts.per_cell, top)?;
            let actual_bytes = archive.reader().bytes_read();
            Ok(SampleCounts {
                archive: archive_index,
                present: true,
                supporting_children: 0,
                chunks: routes[archive_index].len(),
                planned_bytes: bytes,
                actual_bytes,
                shape_route_used: false,
                shape_route_section: None,
                shape_route_sidecar_bytes: 0,
                molecules: counts.molecules,
                umis,
                cells,
                top_cells,
            })
        })
        .collect::<Result<_>>()?;
    let execution_seconds = execution_started.elapsed().as_secs_f64();
    let mut results: Vec<SampleCounts> = (0..collection.archives.len())
        .map(|archive| zero_counts(archive, false))
        .collect();
    for result in computed {
        let archive = result.archive;
        results[archive] = result;
    }
    let archives_opened = routed.len();
    let chunks_decoded: usize = results.iter().map(|result| result.chunks).sum();
    let bytes: u64 = results.iter().map(|result| result.planned_bytes).sum();
    let actual_bytes: u64 = results.iter().map(|result| result.actual_bytes).sum();
    let collection_sidecar_bytes = collection_sidecar_bytes_read(&chain)?;
    let source_archive_bytes = identity_source_bytes
        .checked_add(actual_bytes)
        .context("region source-archive byte count overflow")?;
    let total_logical_bytes = collection_sidecar_bytes
        .checked_add(source_archive_bytes)
        .context("region total logical byte count overflow")?;
    let total_molecules: u64 = results.iter().filter_map(|result| result.molecules).sum();
    let total_umis: usize = results.iter().map(|result| result.umis).sum();
    let total_cells: usize = results.iter().map(|result| result.cells).sum();
    let total_seconds = started.elapsed().as_secs_f64();
    if let Some(format) = uniform_output.format {
        let summary = CollectionRegionUniformSummary {
            coordinates: "0-based half-open",
            anchor_semantics: true,
            chrom: &chrom,
            start,
            end,
            molecules: total_molecules,
            umis: total_umis as u64,
            cells: total_cells as u64,
            selection: CollectionSelectionPolicy {
                requested_top_per_sample: top,
                top_zero_means_all: true,
                comparator: "UMIs descending, archive cell id ascending",
            },
            planning: CollectionQueryPlanningSummary {
                collection_layers: chain.layers.len() as u64,
                archives_total: collection.archives.len() as u64,
                archives_opened: archives_opened as u64,
                archives_pruned: (collection.archives.len() - archives_opened) as u64,
                support_bound_pruned: None,
                unique_chunks_decoded: chunks_decoded as u64,
                independent_chunk_decodes: None,
                planned_compressed_bytes: bytes,
                actual_archive_bytes_read: actual_bytes,
                source_archive_identity_bytes_read: identity_source_bytes,
                source_archive_execution_bytes_read: actual_bytes,
                source_archive_bytes_read: source_archive_bytes,
                collection_sidecar_bytes_read: collection_sidecar_bytes,
                shape_route_sidecar_payload_bytes_read: 0,
                total_logical_bytes_read: total_logical_bytes,
                route_blocks_loaded: 0,
                routed_archives: 0,
                fallback_archives: 0,
            },
        };
        let context = uniform_context
            .as_ref()
            .expect("uniform collection region context was prepared");
        write_uniform_collection_result(&uniform_output, |writer| {
            stream_uniform_region(
                writer, format, context, &summary, collection, &results, explain,
            )
        })?;
        return Ok(());
    }
    let mut value = json!({
        "schema": "gravlax.collection.region.v2",
        "collection_schema": SCHEMA,
        "coordinates": "0-based half-open",
        "anchor_semantics": true,
        "chrom": chrom,
        "start": start,
        "end": end,
        "totals": {"molecules": total_molecules, "umis": total_umis, "cells": total_cells},
        "samples": results.iter().map(|result| sample_counts_json(collection, result)).collect::<Vec<_>>(),
        "planning": {
            "collection_layers": chain.layers.len(),
            "archives_total": collection.archives.len(),
            "archives_opened": archives_opened,
            "archives_pruned": collection.archives.len() - archives_opened,
            "unique_chunks_decoded": chunks_decoded,
            "planned_compressed_bytes": bytes,
            "actual_archive_bytes_read": actual_bytes,
            "source_archive_identity_bytes_read": identity_source_bytes,
            "source_archive_execution_bytes_read": actual_bytes,
            "source_archive_bytes_read": source_archive_bytes,
            "collection_sidecar_bytes_read": collection_sidecar_bytes,
            "shape_route_sidecar_payload_bytes_read": 0,
            "total_logical_bytes_read": total_logical_bytes,
            "route_blocks_loaded": 0,
            "routed_archives": 0,
            "fallback_archives": 0,
            "archive_chunk_index_sections_read": 0,
            "collection_load_seconds": load_seconds,
            "identity_guard_seconds": identity_seconds,
            "identity_content_bytes_read": identity_bytes_read,
            "route_planning_seconds": planning_seconds,
            "source_execution_seconds": execution_seconds,
            "total_seconds": total_seconds,
        },
    });
    if explain {
        value.as_object_mut().unwrap().insert(
            "explain".into(),
            json!(results
                .iter()
                .map(|result| json!({
                    "sample": collection.archives[result.archive].id,
                    "decision": if result.present { "open" } else { "prune_no_overlapping_chunks" },
                    "chunks": result.chunks,
                    "planned_compressed_bytes": result.planned_bytes,
                    "actual_archive_bytes_read": result.actual_bytes,
                }))
                .collect::<Vec<_>>()),
        );
    }
    if json_output {
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        for result in &results {
            println!(
                "{}\t{}\t{}\t{}\t{}",
                collection.archives[result.archive].id,
                result.molecules.unwrap_or(0),
                result.umis,
                result.cells,
                if result.present { "routed" } else { "empty" }
            );
        }
        println!(
            "collection region {locus}: {total_molecules} molecules / {total_umis} UMIs / {total_cells} cells; {archives_opened}/{} archives, {chunks_decoded} chunks, {bytes} compressed bytes ({total_seconds:.3}s)",
            collection.archives.len()
        );
        if explain {
            eprintln!("{}", serde_json::to_string_pretty(&value["planning"])?);
        }
    }
    Ok(())
}

#[derive(Clone, Debug)]
struct RoutedComponent {
    chrom: u32,
    donor: u32,
    acceptor: u32,
    side: u8,
    posts: Vec<u32>,
    shape_route: Option<SpanRoute>,
}

#[derive(Clone, Copy, Debug, Default)]
struct CategoryCounts {
    include_only: usize,
    exclude_only: usize,
    both: usize,
}

impl CategoryCounts {
    fn add_mask(&mut self, mask: u8) {
        match mask {
            1 => self.include_only += 1,
            2 => self.exclude_only += 1,
            3 => self.both += 1,
            _ => {}
        }
    }

    fn add(&mut self, other: Self) {
        self.include_only += other.include_only;
        self.exclude_only += other.exclude_only;
        self.both += other.both;
    }

    fn informative(&self) -> usize {
        self.include_only + self.exclude_only
    }

    fn total(&self) -> usize {
        self.informative() + self.both
    }

    fn usage(&self) -> Option<f64> {
        (self.informative() > 0).then_some(self.include_only as f64 / self.informative() as f64)
    }

    fn json(&self) -> serde_json::Value {
        json!({
            "include_only": self.include_only,
            "exclude_only": self.exclude_only,
            "both": self.both,
            "informative_umis": self.informative(),
            "usage_fraction": self.usage(),
        })
    }
}

fn jset_class_masks_routed(
    archive: &mut LazyArchive,
    chunks: &[ChunkInfo],
    components: &[RoutedComponent],
    routed_n_shapes: Option<u32>,
) -> Result<(FxHashMap<u32, u8>, usize, usize)> {
    let use_shape_routes = components.iter().all(|component| component.shape_route.is_some());
    if use_shape_routes != routed_n_shapes.is_some() {
        bail!("junction-set shape-route plan is internally inconsistent");
    }
    let shapes = (!use_shape_routes).then(|| archive.shapes()).transpose()?;
    let mut chunk_wanted: Vec<FxHashMap<(u32, u32), u8>> =
        (0..chunks.len()).map(|_| FxHashMap::default()).collect();
    let mut chunk_components: Vec<Vec<usize>> = (0..chunks.len()).map(|_| Vec::new()).collect();
    let mut independent_chunk_decodes = 0usize;
    for (component_index, component) in components.iter().enumerate() {
        for &post in &component.posts {
            let chunk = chunks
                .get(post as usize)
                .with_context(|| format!("junction-set route references missing chunk {post}"))?;
            if chunk.chrom != component.chrom {
                bail!("junction-set route references a chunk on the wrong chromosome");
            }
            *chunk_wanted[post as usize]
                .entry((component.donor, component.acceptor))
                .or_insert(0) |= component.side;
            chunk_components[post as usize].push(component_index);
            independent_chunk_decodes += 1;
        }
    }
    let selected: Vec<usize> = chunk_wanted
        .iter()
        .enumerate()
        .filter_map(|(index, wanted)| (!wanted.is_empty()).then_some(index))
        .collect();
    let hits: Vec<Vec<(u32, u8)>> = {
        let (reader, tables) = archive.reader_and_tables();
        let reader = &*reader;
        selected
            .par_iter()
            .map(|&chunk_index| -> Result<Vec<(u32, u8)>> {
                let (compressed, raw_len) =
                    reader.read_compressed_at(&format!("c{chunk_index}"))?;
                let raw = evidence_io::format::decompress(&compressed, raw_len)?;
                let molecules = decode_chunk(&raw, &chunks[chunk_index], None, tables)?;
                let wanted = &chunk_wanted[chunk_index];
                let mut chunk_hits = Vec::new();
                for molecule in &molecules {
                    let mut mask = 0u8;
                    let mut inspect = |position: u32, shape: u32| -> Result<()> {
                        if let Some(n_shapes) = routed_n_shapes {
                            if shape >= n_shapes {
                                bail!(
                                    "molecule references shape {shape} outside the bound {n_shapes}-shape dictionary"
                                );
                            }
                            for &component_index in &chunk_components[chunk_index] {
                                let component = &components[component_index];
                                if component
                                    .shape_route
                                    .as_ref()
                                    .expect("route mode requires every component route")
                                    .matches(
                                        shape,
                                        position,
                                        component.donor,
                                        component.acceptor,
                                    )?
                                {
                                    mask |= component.side;
                                }
                            }
                        } else {
                            let shape = shapes
                                .as_ref()
                                .expect("fallback mode loads shapes")
                                .get(shape as usize)
                                .with_context(|| {
                                    format!("molecule references missing shape {shape}")
                                })?;
                            for blocks in shape.blocks.windows(2) {
                                let donor = position
                                    .checked_add(blocks[0].0)
                                    .and_then(|value| value.checked_add(blocks[0].1))
                                    .context("junction donor coordinate overflow")?;
                                let acceptor = position
                                    .checked_add(blocks[1].0)
                                    .context("junction acceptor coordinate overflow")?;
                                if let Some(side) = wanted.get(&(donor, acceptor)) {
                                    mask |= *side;
                                }
                            }
                        }
                        Ok(())
                    };
                    for chain in &molecule.chains {
                        for &(position, shape) in &chain.reps {
                            inspect(position, shape)?;
                        }
                    }
                    for &(position, shape, _, _) in &molecule.mms {
                        inspect(position, shape)?;
                    }
                    if mask != 0 {
                        chunk_hits.push((molecule.umi_class, mask));
                    }
                }
                chunk_hits.sort_unstable();
                let mut reduced: Vec<(u32, u8)> = Vec::with_capacity(chunk_hits.len());
                for (class, mask) in chunk_hits {
                    match reduced.last_mut() {
                        Some((previous, combined)) if *previous == class => *combined |= mask,
                        _ => reduced.push((class, mask)),
                    }
                }
                Ok(reduced)
            })
            .collect::<Result<_>>()?
    };
    archive.prefetch_coc(hits.iter().flatten().map(|(class, _)| *class))?;
    let mut masks = FxHashMap::default();
    for (class, mask) in hits.into_iter().flatten() {
        *masks.entry(class).or_insert(0) |= mask;
    }
    Ok((masks, selected.len(), independent_chunk_decodes))
}

#[derive(Debug)]
struct SampleJset {
    archive: usize,
    present_components: usize,
    unique_chunks: usize,
    independent_chunks: usize,
    planned_bytes: u64,
    actual_bytes: u64,
    shape_route_used: bool,
    shape_route_sections: Vec<String>,
    shape_route_sidecar_bytes: u64,
    totals: CategoryCounts,
    cells: usize,
    top_cells: Vec<(String, CategoryCounts)>,
}

struct JsetRequest {
    locus: String,
    chrom: u32,
    donor: u32,
    acceptor: u32,
    side: u8,
}

#[derive(Serialize)]
struct JsetCategorySummary {
    include_only: u64,
    exclude_only: u64,
    both: u64,
    informative_umis: u64,
    usage_fraction: Option<f64>,
}

impl From<CategoryCounts> for JsetCategorySummary {
    fn from(value: CategoryCounts) -> Self {
        Self {
            include_only: value.include_only as u64,
            exclude_only: value.exclude_only as u64,
            both: value.both as u64,
            informative_umis: value.informative() as u64,
            usage_fraction: value.usage(),
        }
    }
}

#[derive(Serialize)]
struct CollectionJsetUniformSummary {
    coordinates: &'static str,
    class_categories: [&'static str; 3],
    informative_umis: &'static str,
    usage_fraction: &'static str,
    both_in_usage_denominator: bool,
    min_support: u64,
    totals: JsetCategorySummary,
    selection: CollectionSelectionPolicy,
    planning: CollectionQueryPlanningSummary,
}

fn jset_cell_selection(
    results: &[Option<SampleJset>],
) -> std::result::Result<SelectionSummary, OutputError> {
    let available = results.iter().flatten().try_fold(0u64, |sum, result| {
        sum.checked_add(result.cells as u64)
            .ok_or(OutputError::InvalidSchema(
                "junction-set cell availability exceeds u64".into(),
            ))
    })?;
    let emitted = results.iter().flatten().try_fold(0u64, |sum, result| {
        sum.checked_add(result.top_cells.len() as u64)
            .ok_or(OutputError::InvalidSchema(
                "junction-set emitted cell count exceeds u64".into(),
            ))
    })?;
    SelectionSummary::selected(available, emitted)
}

#[expect(
    clippy::too_many_arguments,
    reason = "the streaming serializer keeps request, routing, selection, and schema inputs explicit at the output boundary"
)]
fn stream_uniform_jset(
    writer: &mut dyn Write,
    format: CollectionOutputFormat,
    context: &ResultContext,
    summary: &CollectionJsetUniformSummary,
    collection: &Collection,
    requests: &[JsetRequest],
    request_junctions: &[Option<GlobalJunction>],
    min_support: u64,
    results: &[Option<SampleJset>],
    support_bound_pruned: bool,
    explain: bool,
) -> std::result::Result<(), OutputError> {
    let request_schema = sequence_table_schema(
        COLLECTION_JSET_REQUESTS_SCHEMA,
        vec![
            Field::new("request_index", DataType::UInt64),
            Field::new("side", DataType::String),
            Field::new("locus", DataType::String),
            Field::new("present_samples", DataType::UInt64),
            Field::new("support_upper_bound", DataType::UInt64),
            Field::new("passes_min_support", DataType::Boolean),
        ],
        "request_index",
    )?;
    let sample_schema = set_table_schema(
        COLLECTION_JSET_SAMPLES_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("present_components", DataType::UInt64),
            Field::new("include_only", DataType::UInt64),
            Field::new("exclude_only", DataType::UInt64),
            Field::new("both", DataType::UInt64),
            Field::new("informative_umis", DataType::UInt64),
            Field::new("usage_fraction", DataType::Float64).nullable(),
            Field::new("cells", DataType::UInt64),
            Field::new("unique_chunks_decoded", DataType::UInt64),
            Field::new("independent_chunk_decodes", DataType::UInt64),
            Field::new("planned_compressed_bytes", DataType::UInt64),
            Field::new("actual_archive_bytes_read", DataType::UInt64),
            Field::new("shape_route_used", DataType::Boolean),
            Field::new("shape_route_sections", DataType::UInt64),
            Field::new("shape_route_sidecar_bytes_read", DataType::UInt64),
            Field::new("decision", DataType::String).nullable(),
        ],
        &["sample"],
    )?;
    let cell_schema = set_table_schema(
        COLLECTION_JSET_CELLS_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("rank", DataType::UInt64),
            Field::new("barcode", DataType::String),
            Field::new("include_only", DataType::UInt64),
            Field::new("exclude_only", DataType::UInt64),
            Field::new("both", DataType::UInt64),
            Field::new("informative_umis", DataType::UInt64),
            Field::new("usage_fraction", DataType::Float64).nullable(),
        ],
        &["sample", "barcode"],
    )?;
    let request_selection = SelectionSummary::complete(requests.len() as u64);
    let sample_selection = SelectionSummary::complete(collection.archives.len() as u64);
    let cell_selection = jset_cell_selection(results)?;
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        COLLECTION_JSET_RESULT_SCHEMA,
        OutputFormat::from(format),
        context,
        summary,
    )?;
    bundle.write_table(
        "requests",
        &request_schema,
        Some(&request_selection),
        |rows| {
            for (index, (request, junction)) in requests.iter().zip(request_junctions).enumerate() {
                let support = junction.as_ref().map_or(0, |row| row.support_upper_bound);
                rows.write_row_with(|row| {
                    row.uint64(index as u64)?;
                    row.string(if request.side == 1 {
                        "include"
                    } else {
                        "exclude"
                    })?;
                    row.string(&request.locus)?;
                    row.uint64(junction.as_ref().map_or(0, |row| row.routes.len()) as u64)?;
                    row.uint64(support)?;
                    row.boolean(support >= min_support)?;
                    Ok(())
                })?;
            }
            Ok(())
        },
    )?;
    bundle.write_table("samples", &sample_schema, Some(&sample_selection), |rows| {
        for (archive, result) in results.iter().enumerate() {
            let counts = result
                .as_ref()
                .map_or(CategoryCounts::default(), |row| row.totals);
            rows.write_row_with(|row| {
                row.string(&collection.archives[archive].id)?;
                row.uint64(result.as_ref().map_or(0, |row| row.present_components) as u64)?;
                row.uint64(counts.include_only as u64)?;
                row.uint64(counts.exclude_only as u64)?;
                row.uint64(counts.both as u64)?;
                row.uint64(counts.informative() as u64)?;
                match counts.usage() {
                    Some(usage) => row.float64(usage)?,
                    None => row.null()?,
                }
                row.uint64(result.as_ref().map_or(0, |row| row.cells) as u64)?;
                row.uint64(result.as_ref().map_or(0, |row| row.unique_chunks) as u64)?;
                row.uint64(result.as_ref().map_or(0, |row| row.independent_chunks) as u64)?;
                row.uint64(result.as_ref().map_or(0, |row| row.planned_bytes))?;
                row.uint64(result.as_ref().map_or(0, |row| row.actual_bytes))?;
                row.boolean(result.as_ref().is_some_and(|row| row.shape_route_used))?;
                row.uint64(
                    result
                        .as_ref()
                        .map_or(0, |row| row.shape_route_sections.len()) as u64,
                )?;
                row.uint64(
                    result
                        .as_ref()
                        .map_or(0, |row| row.shape_route_sidecar_bytes),
                )?;
                if explain {
                    row.string(if support_bound_pruned {
                        "prune_support_bound"
                    } else if result.is_some() {
                        "open"
                    } else {
                        "prune_all_components_absent"
                    })?;
                } else {
                    row.null()?;
                }
                Ok(())
            })?;
        }
        Ok(())
    })?;
    bundle.write_table("cells", &cell_schema, Some(&cell_selection), |rows| {
        for result in results.iter().flatten() {
            let sample = &collection.archives[result.archive].id;
            for (rank, (barcode, counts)) in result.top_cells.iter().enumerate() {
                rows.write_row_with(|row| {
                    row.string(sample)?;
                    row.uint64(rank as u64 + 1)?;
                    row.string(barcode)?;
                    row.uint64(counts.include_only as u64)?;
                    row.uint64(counts.exclude_only as u64)?;
                    row.uint64(counts.both as u64)?;
                    row.uint64(counts.informative() as u64)?;
                    match counts.usage() {
                        Some(usage) => row.float64(usage)?,
                        None => row.null()?,
                    }
                    Ok(())
                })?;
            }
        }
        Ok(())
    })?;
    bundle.finish()?;
    Ok(())
}

struct JsetQueryRun {
    path: PathBuf,
    includes: Vec<String>,
    excludes: Vec<String>,
    min_support: u64,
    top: usize,
    explain: bool,
    verify_content: bool,
    json_output: bool,
    uniform_output: CollectionOutputArgs,
}

fn run_jset_query(args: JsetQueryRun) -> Result<()> {
    let JsetQueryRun {
        path,
        includes,
        excludes,
        min_support,
        top,
        explain,
        verify_content,
        json_output,
        uniform_output,
    } = args;
    let started = std::time::Instant::now();
    let chain = open_collection_chain(&path)?;
    let collection = &chain.collection;
    let load_seconds = started.elapsed().as_secs_f64();
    let uniform_context = if uniform_output.format.is_some() {
        let mut parameters = BTreeMap::new();
        parameters.insert("include".into(), json!(&includes));
        parameters.insert("exclude".into(), json!(&excludes));
        parameters.insert("min_support".into(), json!(min_support));
        parameters.insert("top".into(), json!(top));
        parameters.insert("explain".into(), json!(explain));
        parameters.insert("verify_content".into(), json!(verify_content));
        parameters.insert(
            "archive_access".into(),
            json!("collection catalogue/postings with shared chunk decoding and optional source-root-bound local-shape routes"),
        );
        Some(collection_uniform_context(&path, &chain, parameters)?)
    } else {
        None
    };
    let identity_started = std::time::Instant::now();
    let identity_io = validate_collection_sources(collection, verify_content)?;
    let identity_seconds = identity_started.elapsed().as_secs_f64();
    let identity_bytes_read = identity_content_bytes(&identity_io)?;
    let identity_source_bytes = source_archive_bytes_read(&identity_io)?;
    let planning_started = std::time::Instant::now();

    let mut requests = Vec::new();
    let mut seen: BTreeMap<(u32, u32, u32), u8> = BTreeMap::new();
    for (label, side, loci) in [("include", 1u8, includes), ("exclude", 2u8, excludes)] {
        for locus in loci {
            let (chrom_name, donor, acceptor) =
                parse_locus(&locus).with_context(|| format!("invalid --{label} junction"))?;
            let chrom = chrom_id(collection, &chrom_name)?;
            let key = (chrom, donor, acceptor);
            if let Some(previous) = seen.insert(key, side) {
                if previous == side {
                    bail!("duplicate --{label} junction {locus}");
                }
                bail!("junction {locus} appears on both inclusion and exclusion sides");
            }
            requests.push(JsetRequest {
                locus,
                chrom,
                donor,
                acceptor,
                side,
            });
        }
    }
    let mut per_archive: Vec<Vec<RoutedComponent>> =
        (0..collection.archives.len()).map(|_| Vec::new()).collect();
    let request_loci: Vec<(u32, u32, u32)> = requests
        .iter()
        .map(|request| (request.chrom, request.donor, request.acceptor))
        .collect();
    let request_junctions = lookup_chain_junctions(&chain, &request_loci)?;
    let support_bound_pruned = request_junctions
        .iter()
        .any(|row| row.as_ref().map_or(0, |row| row.support_upper_bound) < min_support);
    if !support_bound_pruned {
        for (request, junction) in requests.iter().zip(&request_junctions) {
            if let Some(junction) = junction {
                for route in &junction.routes {
                    per_archive[route.archive as usize].push(RoutedComponent {
                        chrom: request.chrom,
                        donor: request.donor,
                        acceptor: request.acceptor,
                        side: request.side,
                        posts: route.posts.clone(),
                        shape_route: None,
                    });
                }
            }
        }
    }
    let routed: Vec<usize> = per_archive
        .iter()
        .enumerate()
        .filter_map(|(archive, requests)| (!requests.is_empty()).then_some(archive))
        .collect();
    let mut shape_route_plans: Vec<Option<ArchiveShapeRoutePlan>> =
        (0..collection.archives.len()).map(|_| None).collect();
    for &archive_index in &routed {
        let spans = per_archive[archive_index]
            .iter()
            .map(|component| {
                component
                    .acceptor
                    .checked_sub(component.donor)
                    .context("junction-set acceptor precedes donor")
            })
            .collect::<Result<Vec<_>>>()?;
        if let Some(plan) = lookup_chain_shape_routes(&chain, archive_index, &spans)? {
            for component in &mut per_archive[archive_index] {
                let span = component
                    .acceptor
                    .checked_sub(component.donor)
                    .context("junction-set acceptor precedes donor")?;
                component.shape_route = Some(
                    plan.spans
                        .get(&span)
                        .cloned()
                        .context("planned junction-set shape route is missing its exact span")?,
                );
            }
            shape_route_plans[archive_index] = Some(plan);
        }
    }
    let planning_seconds = planning_started.elapsed().as_secs_f64();
    let execution_started = std::time::Instant::now();
    let computed: Vec<SampleJset> = routed
        .par_iter()
        .map(|&archive_index| -> Result<SampleJset> {
            let entry = &collection.archives[archive_index];
            let chunks = chunk_infos(entry);
            let mut union_posts: Vec<u32> = per_archive[archive_index]
                .iter()
                .flat_map(|request| request.posts.iter().copied())
                .collect();
            union_posts.sort_unstable();
            union_posts.dedup();
            let bytes = planned_bytes(entry, &union_posts)?;
            let shape_plan = shape_route_plans[archive_index].as_ref();
            let mut archive = open_source(
                collection,
                archive_index,
                shape_plan.map(|plan| &plan.binding),
            )?;
            let (masks, unique_chunks, independent_chunks) =
                jset_class_masks_routed(
                    &mut archive,
                    &chunks,
                    &per_archive[archive_index],
                    shape_plan.map(|plan| plan.binding.n_shapes),
                )?;
            let mut per_cell: FxHashMap<u32, CategoryCounts> = FxHashMap::default();
            for (class, mask) in masks {
                per_cell
                    .entry(archive.cell_of(class)?)
                    .or_default()
                    .add_mask(mask);
            }
            let mut cells: Vec<(u32, CategoryCounts)> = per_cell.into_iter().collect();
            cells.sort_unstable_by_key(|(cell, counts)| (std::cmp::Reverse(counts.total()), *cell));
            let mut totals = CategoryCounts::default();
            for (_, counts) in &cells {
                totals.add(*counts);
            }
            let dictionary = archive.cells()?;
            let limit = if top == 0 {
                cells.len()
            } else {
                top.min(cells.len())
            };
            let top_cells = cells
                .iter()
                .take(limit)
                .map(|(cell, counts)| {
                    let packed = dictionary.get(*cell as usize).copied().with_context(|| {
                        format!("cell id {cell} is outside the barcode dictionary")
                    })?;
                    Ok((unpack_barcode(packed), *counts))
                })
                .collect::<Result<_>>()?;
            let actual_bytes = archive.reader().bytes_read();
            Ok(SampleJset {
                archive: archive_index,
                present_components: per_archive[archive_index].len(),
                unique_chunks,
                independent_chunks,
                planned_bytes: bytes,
                actual_bytes,
                shape_route_used: shape_plan.is_some(),
                shape_route_sections: shape_plan
                    .map_or_else(Vec::new, |plan| plan.section_names.clone()),
                shape_route_sidecar_bytes: shape_plan
                    .map_or(0, |plan| plan.sidecar_compressed_bytes),
                totals,
                cells: cells.len(),
                top_cells,
            })
        })
        .collect::<Result<_>>()?;
    let execution_seconds = execution_started.elapsed().as_secs_f64();
    let mut result_of: Vec<Option<SampleJset>> =
        (0..collection.archives.len()).map(|_| None).collect();
    for result in computed {
        let archive = result.archive;
        result_of[archive] = Some(result);
    }
    let mut grand = CategoryCounts::default();
    for result in result_of.iter().flatten() {
        grand.add(result.totals);
    }
    let unique_chunks: usize = result_of
        .iter()
        .flatten()
        .map(|result| result.unique_chunks)
        .sum();
    let independent_chunks: usize = result_of
        .iter()
        .flatten()
        .map(|result| result.independent_chunks)
        .sum();
    let bytes: u64 = result_of
        .iter()
        .flatten()
        .map(|result| result.planned_bytes)
        .sum();
    let actual_bytes: u64 = result_of
        .iter()
        .flatten()
        .map(|result| result.actual_bytes)
        .sum();
    let route_sidecar_bytes: u64 = result_of
        .iter()
        .flatten()
        .map(|result| result.shape_route_sidecar_bytes)
        .sum();
    let route_blocks_loaded: usize = result_of
        .iter()
        .flatten()
        .map(|result| result.shape_route_sections.len())
        .sum();
    let routed_archives = result_of
        .iter()
        .flatten()
        .filter(|result| result.shape_route_used)
        .count();
    let collection_sidecar_bytes = collection_sidecar_bytes_read(&chain)?;
    let source_archive_bytes = identity_source_bytes
        .checked_add(actual_bytes)
        .context("junction-set source-archive byte count overflow")?;
    let total_logical_bytes = collection_sidecar_bytes
        .checked_add(source_archive_bytes)
        .context("junction-set total logical byte count overflow")?;
    let total_seconds = started.elapsed().as_secs_f64();
    if let Some(format) = uniform_output.format {
        let summary = CollectionJsetUniformSummary {
            coordinates: "0-based junction boundaries",
            class_categories: ["include_only", "exclude_only", "both"],
            informative_umis: "include_only + exclude_only",
            usage_fraction: "include_only / informative_umis",
            both_in_usage_denominator: false,
            min_support,
            totals: grand.into(),
            selection: CollectionSelectionPolicy {
                requested_top_per_sample: top,
                top_zero_means_all: true,
                comparator: "total category UMIs descending, archive cell id ascending",
            },
            planning: CollectionQueryPlanningSummary {
                collection_layers: chain.layers.len() as u64,
                archives_total: collection.archives.len() as u64,
                archives_opened: routed.len() as u64,
                archives_pruned: (collection.archives.len() - routed.len()) as u64,
                support_bound_pruned: Some(support_bound_pruned),
                unique_chunks_decoded: unique_chunks as u64,
                independent_chunk_decodes: Some(independent_chunks as u64),
                planned_compressed_bytes: bytes,
                actual_archive_bytes_read: actual_bytes,
                source_archive_identity_bytes_read: identity_source_bytes,
                source_archive_execution_bytes_read: actual_bytes,
                source_archive_bytes_read: source_archive_bytes,
                collection_sidecar_bytes_read: collection_sidecar_bytes,
                shape_route_sidecar_payload_bytes_read: route_sidecar_bytes,
                total_logical_bytes_read: total_logical_bytes,
                route_blocks_loaded: route_blocks_loaded as u64,
                routed_archives: routed_archives as u64,
                fallback_archives: (routed.len() - routed_archives) as u64,
            },
        };
        let context = uniform_context
            .as_ref()
            .expect("uniform collection junction-set context was prepared");
        write_uniform_collection_result(&uniform_output, |writer| {
            stream_uniform_jset(
                writer,
                format,
                context,
                &summary,
                collection,
                &requests,
                &request_junctions,
                min_support,
                &result_of,
                support_bound_pruned,
                explain,
            )
        })?;
        return Ok(());
    }
    let request_rows = requests
        .iter()
        .zip(&request_junctions)
        .map(|(request, junction)| {
            json!({
                "side": if request.side == 1 { "include" } else { "exclude" },
                "locus": request.locus,
                "present_samples": junction.as_ref().map_or(0, |row| row.routes.len()),
                "support_upper_bound": junction.as_ref().map_or(0, |row| row.support_upper_bound),
                "passes_min_support": junction.as_ref().map_or(0, |row| row.support_upper_bound) >= min_support,
            })
        })
        .collect::<Vec<_>>();
    let sample_rows = result_of
        .iter()
        .enumerate()
        .map(|(archive, result)| match result {
            Some(result) => json!({
                "sample": collection.archives[archive].id,
                "archive": collection.archives[archive].path,
                "present_components": result.present_components,
                "totals": result.totals.json(),
                "cells": result.cells,
                "top_cells": result.top_cells.iter().map(|(barcode, counts)| {
                    let mut value = counts.json();
                    value.as_object_mut().unwrap().insert("barcode".into(), json!(barcode));
                    value
                }).collect::<Vec<_>>(),
                "unique_chunks_decoded": result.unique_chunks,
                "independent_chunk_decodes": result.independent_chunks,
                "planned_compressed_bytes": result.planned_bytes,
                "actual_archive_bytes_read": result.actual_bytes,
                "shape_route_used": result.shape_route_used,
                "shape_route_sections": result.shape_route_sections,
                "shape_route_sidecar_bytes_read": result.shape_route_sidecar_bytes,
            }),
            None => json!({
                "sample": collection.archives[archive].id,
                "archive": collection.archives[archive].path,
                "present_components": 0,
                "totals": CategoryCounts::default().json(),
                "cells": 0,
                "top_cells": [],
                "unique_chunks_decoded": 0,
                "independent_chunk_decodes": 0,
                "planned_compressed_bytes": 0,
                "actual_archive_bytes_read": 0,
                "shape_route_used": false,
                "shape_route_sections": [],
                "shape_route_sidecar_bytes_read": 0,
            }),
        })
        .collect::<Vec<_>>();
    let mut value = json!({
        "schema": "gravlax.collection.jset.v2",
        "collection_schema": SCHEMA,
        "coordinates": "0-based junction boundaries",
        "semantics": {
            "class_categories": ["include_only", "exclude_only", "both"],
            "informative_umis": "include_only + exclude_only",
            "usage_fraction": "include_only / informative_umis",
            "both_in_usage_denominator": false,
        },
        "junctions": request_rows,
        "min_support": min_support,
        "totals": grand.json(),
        "samples": sample_rows,
        "planning": {
            "collection_layers": chain.layers.len(),
            "archives_total": collection.archives.len(),
            "archives_opened": routed.len(),
            "archives_pruned": collection.archives.len() - routed.len(),
            "support_bound_pruned": support_bound_pruned,
            "unique_chunks_decoded": unique_chunks,
            "independent_chunk_decodes": independent_chunks,
            "chunk_decode_reduction_fraction": if independent_chunks == 0 { 0.0 } else {
                1.0 - unique_chunks as f64 / independent_chunks as f64
            },
            "planned_compressed_bytes": bytes,
            "actual_archive_bytes_read": actual_bytes,
            "source_archive_identity_bytes_read": identity_source_bytes,
            "source_archive_execution_bytes_read": actual_bytes,
            "source_archive_bytes_read": source_archive_bytes,
            "collection_sidecar_bytes_read": collection_sidecar_bytes,
            "shape_route_sidecar_payload_bytes_read": route_sidecar_bytes,
            "total_logical_bytes_read": total_logical_bytes,
            "route_blocks_loaded": route_blocks_loaded,
            "routed_archives": routed_archives,
            "fallback_archives": routed.len() - routed_archives,
            "archive_catalogue_sections_read": 0,
            "archive_posting_sections_read": 0,
            "collection_load_seconds": load_seconds,
            "identity_guard_seconds": identity_seconds,
            "identity_content_bytes_read": identity_bytes_read,
            "route_planning_seconds": planning_seconds,
            "source_execution_seconds": execution_seconds,
            "total_seconds": total_seconds,
        },
    });
    if explain {
        value.as_object_mut().unwrap().insert(
            "explain".into(),
            json!(result_of
                .iter()
                .enumerate()
                .map(|(archive, result)| json!({
                    "sample": collection.archives[archive].id,
                    "decision": if support_bound_pruned {
                        "prune_support_bound"
                    } else if result.is_some() {
                        "open"
                    } else {
                        "prune_all_components_absent"
                    },
                    "present_components": result.as_ref().map_or(0, |row| row.present_components),
                    "unique_chunks": result.as_ref().map_or(0, |row| row.unique_chunks),
                    "planned_compressed_bytes": result.as_ref().map_or(0, |row| row.planned_bytes),
                    "actual_archive_bytes_read": result.as_ref().map_or(0, |row| row.actual_bytes),
                }))
                .collect::<Vec<_>>()),
        );
    }
    if json_output {
        println!("{}", serde_json::to_string_pretty(&value)?);
    } else {
        for (archive, result) in result_of.iter().enumerate() {
            let counts = result
                .as_ref()
                .map_or(CategoryCounts::default(), |row| row.totals);
            println!(
                "{}\t{}\t{}\t{}\t{}",
                collection.archives[archive].id,
                counts.include_only,
                counts.exclude_only,
                counts.both,
                counts
                    .usage()
                    .map_or_else(|| "NA".to_owned(), |usage| format!("{usage:.6}"))
            );
        }
        println!(
            "collection jset: {} include-only / {} exclude-only / {} both; {}/{} archives, {unique_chunks}/{independent_chunks} chunks, {bytes} compressed bytes ({total_seconds:.3}s)",
            grand.include_only,
            grand.exclude_only,
            grand.both,
            routed.len(),
            collection.archives.len(),
        );
        if explain {
            eprintln!("{}", serde_json::to_string_pretty(&value["planning"])?);
        }
    }
    Ok(())
}

#[derive(Clone, Debug)]
struct InspectedArchiveRoutes {
    archive_ordinal: usize,
    blocks: Vec<RouteBlock>,
    compressed_bytes: Vec<u64>,
    verification_source_bytes: Option<u64>,
}

impl InspectedArchiveRoutes {
    fn section_count(&self) -> usize {
        self.blocks.len()
    }

    fn span_count(&self) -> usize {
        self.blocks.iter().map(|block| block.spans.len()).sum()
    }

    fn pair_count(&self) -> usize {
        self.blocks.iter().map(RouteBlock::pair_count).sum()
    }

    fn compressed_byte_count(&self) -> Result<u64> {
        self.compressed_bytes.iter().try_fold(0u64, |sum, bytes| {
            sum.checked_add(*bytes)
                .context("shape-route compressed byte count overflow")
        })
    }
}

#[derive(Clone, Debug, Default)]
struct CollectionLayerAudit {
    shape_routes: Vec<InspectedArchiveRoutes>,
}

impl CollectionLayerAudit {
    fn section_count(&self) -> usize {
        self.shape_routes
            .iter()
            .map(InspectedArchiveRoutes::section_count)
            .sum()
    }

    fn span_count(&self) -> usize {
        self.shape_routes
            .iter()
            .map(InspectedArchiveRoutes::span_count)
            .sum()
    }

    fn pair_count(&self) -> usize {
        self.shape_routes
            .iter()
            .map(InspectedArchiveRoutes::pair_count)
            .sum()
    }

    fn compressed_byte_count(&self) -> Result<u64> {
        self.shape_routes.iter().try_fold(0u64, |sum, routes| {
            sum.checked_add(routes.compressed_byte_count()?)
                .context("shape-route compressed byte count overflow")
        })
    }

    fn verification_source_byte_count(&self) -> Result<u64> {
        self.shape_routes.iter().try_fold(0u64, |sum, routes| {
            sum.checked_add(routes.verification_source_bytes.unwrap_or(0))
                .context("shape-route verification source byte count overflow")
        })
    }
}

fn audit_collection_layer(layer: &CollectionLayer) -> Result<CollectionLayerAudit> {
    let observed_shape_sections = layer
        .file
        .names()
        .filter(|name| name.starts_with("s."))
        .collect::<Vec<_>>();
    let expected_shape_sections = layer
        .manifest
        .archives
        .iter()
        .flat_map(|archive| {
            archive
                .shape_routes
                .iter()
                .flat_map(|binding| binding.blocks.iter())
                .map(|descriptor| descriptor.section_name.as_str())
        })
        .collect::<Vec<_>>();
    if observed_shape_sections != expected_shape_sections {
        bail!("collection shape-route sections do not exactly match the manifest bindings");
    }
    let mut rows = 0usize;
    let mut routes = 0usize;
    let mut postings = 0usize;
    let mut previous_section = None;
    let mut previous_coordinate = None;
    for name in layer.file.names() {
        if name == "manifest" {
            layer.file.read(name)?;
            continue;
        }
        if name.starts_with("s.") {
            continue;
        }
        let (chrom, bin) = parse_segment_name(name)?;
        let section_key = (chrom, bin);
        if previous_section.is_some_and(|previous| section_key <= previous) {
            bail!("collection junction sections are not in canonical genomic order");
        }
        previous_section = Some(section_key);
        let decoded = decode_junction_rows(
            &layer.file.read(name)?,
            &layer.manifest.chroms,
            &layer.manifest.archives,
        )?;
        if decoded.is_empty() {
            bail!("collection junction section {name} is empty");
        }
        for row in &decoded {
            if row.chrom != chrom || row.donor / JUNCTION_BIN_BP != bin {
                bail!("collection section {name} contains a row outside its routing bin");
            }
            let coordinate = (row.chrom, row.donor, row.acceptor);
            if previous_coordinate.is_some_and(|previous| coordinate <= previous) {
                bail!("collection junction rows are duplicated or globally unsorted");
            }
            previous_coordinate = Some(coordinate);
            routes = routes
                .checked_add(row.routes.len())
                .context("collection route audit count overflow")?;
            for route in &row.routes {
                postings = postings
                    .checked_add(route.posts.len())
                    .context("collection posting audit count overflow")?;
            }
        }
        rows = rows
            .checked_add(decoded.len())
            .context("collection junction audit count overflow")?;
    }
    if rows != layer.manifest.junction_count
        || routes != layer.manifest.route_count
        || postings != layer.manifest.posting_count
    {
        bail!(
            "collection layer cardinalities differ from its manifest: rows {rows}/{}, routes {routes}/{}, postings {postings}/{}",
            layer.manifest.junction_count,
            layer.manifest.route_count,
            layer.manifest.posting_count
        );
    }

    let mut shape_routes = Vec::new();
    for (archive_ordinal, archive) in layer.manifest.archives.iter().enumerate() {
        let Some(binding) = &archive.shape_routes else {
            continue;
        };
        let mut blocks = Vec::with_capacity(binding.blocks.len());
        let mut compressed_bytes = Vec::with_capacity(binding.blocks.len());
        for descriptor in &binding.blocks {
            let raw = layer.file.read(&descriptor.section_name)?;
            blocks.push(shaperoute::decode_block(
                &raw,
                archive_ordinal as u32,
                binding.n_shapes,
                descriptor.first_span,
                descriptor.last_span,
            )?);
            compressed_bytes.push(layer.file.compressed_len(&descriptor.section_name)?);
        }
        shaperoute::validate_blocks_against_binding(binding, &blocks).with_context(|| {
            format!(
                "validating shape-route directory for sample {}",
                archive.id
            )
        })?;
        shape_routes.push(InspectedArchiveRoutes {
            archive_ordinal,
            blocks,
            compressed_bytes,
            verification_source_bytes: None,
        });
    }
    Ok(CollectionLayerAudit { shape_routes })
}

fn verify_inspected_archive_routes(
    archive: &ArchiveEntry,
    binding: &ShapeRouteBinding,
    blocks: &[RouteBlock],
) -> Result<u64> {
    let file = std::fs::File::open(&archive.path)
        .with_context(|| format!("opening routed archive {}", archive.path.display()))?;
    let before = identity_from_metadata(file.metadata()?, 0, String::new(), String::new(), None)?;
    if !identity_stat_matches(&archive.identity, &before) {
        bail!(
            "archive identity changed for sample {} at {}; rebuild the collection index",
            archive.id,
            archive.path.display()
        );
    }
    let mut reader = evidence_io::format::SectionReader::from_file(file)?;
    if reader.archive_version() != evidence_io::format::VERSION {
        bail!(
            "shape-route verification requires a rooted v2 archive for sample {}",
            archive.id
        );
    }
    let root = reader
        .content_commitment()
        .context("rooted route source lacks its directory commitment")?;
    let encoded = reader
        .encoded_content_identity()?
        .context("rooted route source lacks its encoded-sections identity")?;
    if archive.identity.native_scheme != ROOTED_DIRECTORY_SCHEME
        || root.to_hex() != archive.identity.native_digest
        || archive.identity.encoded_sections_digest.as_deref()
            != Some(digest_hex(encoded).as_str())
    {
        bail!(
            "archive root identity changed for sample {} at {}; rebuild the collection index",
            archive.id,
            archive.path.display()
        );
    }
    let shapes_digest = reader
        .section_metadata()
        .find(|section| section.name == "shapes")
        .context("routed archive lacks its shapes section")?
        .compressed_blake3
        .context("routed archive shapes entry lacks its committed payload digest")?;
    shaperoute::validate_binding(
        binding,
        binding.archive_ordinal,
        root.digest,
        shapes_digest,
    )
    .with_context(|| format!("validating shape-route source for sample {}", archive.id))?;
    let shapes_raw = reader.read("shapes")?;
    shaperoute::verify_reconstruction(&shapes_raw, binding, blocks).with_context(|| {
        format!(
            "reconstructing shape routes for sample {} from its bound source",
            archive.id
        )
    })?;
    let after = identity_from_metadata(
        reader.file_metadata()?,
        0,
        String::new(),
        String::new(),
        None,
    )?;
    if !identity_stat_matches(&archive.identity, &after) {
        bail!(
            "archive identity changed while verifying shape routes for sample {} at {}; rebuild the collection index",
            archive.id,
            archive.path.display()
        );
    }
    Ok(reader.bytes_read())
}

fn route_binding_json(
    binding: &ShapeRouteBinding,
    inspected: &InspectedArchiveRoutes,
    reconstruction_verified: bool,
) -> Result<serde_json::Value> {
    if binding.blocks.len() != inspected.blocks.len()
        || binding.blocks.len() != inspected.compressed_bytes.len()
    {
        bail!("inspected shape-route block directory is internally inconsistent");
    }
    let blocks = binding
        .blocks
        .iter()
        .zip(&inspected.blocks)
        .zip(&inspected.compressed_bytes)
        .map(|((descriptor, block), compressed_bytes)| {
            json!({
                "first_span": descriptor.first_span,
                "last_span": descriptor.last_span,
                "section_name": descriptor.section_name,
                "exact_spans": block.spans.len(),
                "pairs": block.pair_count(),
                "compressed_bytes": compressed_bytes,
            })
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "codec_version": shaperoute::ROUTE_CODEC_VERSION,
        "archive_ordinal": binding.archive_ordinal,
        "source_root": digest_hex(binding.source_root),
        "shapes_digest": digest_hex(binding.shapes_digest),
        "n_shapes": binding.n_shapes,
        "sections": inspected.section_count(),
        "exact_spans": inspected.span_count(),
        "pairs": inspected.pair_count(),
        "compressed_bytes": inspected.compressed_byte_count()?,
        "reconstruction_verified": reconstruction_verified,
        "descriptors": blocks,
    }))
}

#[derive(Serialize)]
struct InspectReferenceSummary<'a> {
    stamped: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    algo: Option<&'a str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    digest: Option<&'a str>,
}

#[derive(Serialize)]
struct InspectIndexSummary {
    segment_junction_rows: u64,
    archive_routes: u64,
    chunk_postings: u64,
    interval_chunks: u64,
    junction_bin_bp: u32,
    junction_segments: u64,
    shape_route_archives: u64,
    shape_route_sections: u64,
    shape_route_exact_spans: u64,
    shape_route_pairs: u64,
    shape_route_compressed_bytes: u64,
}

#[derive(Serialize)]
struct InspectGuardSummary {
    filesystem_identity: &'static str,
    content_digest_recorded: bool,
    content_digest_verified: bool,
    shape_route_payloads_verified: bool,
    shape_route_reconstruction_verified: bool,
}

#[derive(Serialize)]
struct InspectIoSummary {
    collection_sidecar_bytes_read: u64,
    source_identity_content_bytes_read: u64,
    source_identity_total_bytes_read: u64,
    route_verification_source_bytes_read: u64,
}

#[derive(Serialize)]
struct InspectUniformSummary<'a> {
    collection_schema: &'static str,
    format_version: Option<u32>,
    path: &'a str,
    file_bytes: u64,
    reference: InspectReferenceSummary<'a>,
    index: InspectIndexSummary,
    guard: InspectGuardSummary,
    io: InspectIoSummary,
}

fn stream_uniform_inspect(
    writer: &mut dyn Write,
    format: CollectionOutputFormat,
    context: &ResultContext,
    summary: &InspectUniformSummary<'_>,
    chain: &CollectionChain,
    audits: &[CollectionLayerAudit],
    reconstruct_routes: bool,
) -> std::result::Result<(), OutputError> {
    for layer in &chain.layers {
        if layer.path.to_str().is_none() {
            return Err(OutputError::Sink(format!(
                "collection layer path is not valid UTF-8: {}",
                layer.path.display()
            )));
        }
    }
    for (global, archive) in chain.collection.archives.iter().enumerate() {
        if archive.path.to_str().is_none() {
            return Err(OutputError::Sink(format!(
                "archive path is not valid UTF-8: {}",
                archive.path.display()
            )));
        }
        let (layer_index, local_archive) = chain.global_to_local[global];
        let inspected = audits[layer_index]
            .shape_routes
            .iter()
            .find(|routes| routes.archive_ordinal == local_archive);
        match (&archive.shape_routes, inspected) {
            (Some(binding), Some(inspected)) => {
                if binding.blocks.len() != inspected.blocks.len()
                    || binding.blocks.len() != inspected.compressed_bytes.len()
                {
                    return Err(OutputError::Sink(
                        "inspected shape-route block directory is internally inconsistent".into(),
                    ));
                }
                inspected
                    .compressed_byte_count()
                    .map_err(|error| OutputError::Sink(error.to_string()))?;
            }
            (None, None) => {}
            _ => {
                return Err(OutputError::Sink(
                    "flattened collection route binding differs from its owning layer".into(),
                ));
            }
        }
    }
    let layer_schema = sequence_table_schema(
        COLLECTION_INSPECT_LAYERS_SCHEMA,
        vec![
            Field::new("layer_index", DataType::UInt64),
            Field::new("path", DataType::String),
            Field::new("format_version", DataType::UInt64),
            Field::new("root_digest", DataType::String),
            Field::new("archives", DataType::UInt64),
            Field::new("junction_rows", DataType::UInt64),
            Field::new("junction_segments", DataType::UInt64),
            Field::new("shape_route_archives", DataType::UInt64),
            Field::new("shape_route_sections", DataType::UInt64),
            Field::new("shape_route_exact_spans", DataType::UInt64),
            Field::new("shape_route_pairs", DataType::UInt64),
            Field::new("shape_route_compressed_bytes", DataType::UInt64),
            Field::new("shape_route_reconstruction_verified", DataType::Boolean),
        ],
        "layer_index",
    )?;
    let chrom_schema = sequence_table_schema(
        COLLECTION_INSPECT_CHROMS_SCHEMA,
        vec![
            Field::new("chrom_index", DataType::UInt64),
            Field::new("name", DataType::String),
        ],
        "chrom_index",
    )?;
    let archive_schema = set_table_schema(
        COLLECTION_INSPECT_ARCHIVES_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("path", DataType::String),
            Field::new("layer_index", DataType::UInt64),
            Field::new("layer_archive_index", DataType::UInt64),
            Field::new("bytes", DataType::UInt64),
            Field::new("archive_format_version", DataType::UInt64),
            Field::new("native_identity_scheme", DataType::String),
            Field::new("native_identity_digest", DataType::String),
            Field::new("encoded_sections_digest", DataType::String).nullable(),
            Field::new("chunks", DataType::UInt64),
            Field::new("shape_routes", DataType::Boolean),
            Field::new("shape_route_codec_version", DataType::UInt64).nullable(),
            Field::new("shape_route_archive_ordinal", DataType::UInt64).nullable(),
            Field::new("shape_route_source_root", DataType::String).nullable(),
            Field::new("shape_route_shapes_digest", DataType::String).nullable(),
            Field::new("shape_route_n_shapes", DataType::UInt64).nullable(),
            Field::new("shape_route_sections", DataType::UInt64),
            Field::new("shape_route_exact_spans", DataType::UInt64),
            Field::new("shape_route_pairs", DataType::UInt64),
            Field::new("shape_route_compressed_bytes", DataType::UInt64),
            Field::new("shape_route_reconstruction_verified", DataType::Boolean),
        ],
        &["sample"],
    )?;
    let route_block_schema = set_table_schema(
        COLLECTION_INSPECT_ROUTE_BLOCKS_SCHEMA,
        vec![
            Field::new("sample", DataType::String),
            Field::new("section_name", DataType::String),
            Field::new("first_span", DataType::UInt64),
            Field::new("last_span", DataType::UInt64),
            Field::new("exact_spans", DataType::UInt64),
            Field::new("pairs", DataType::UInt64),
            Field::new("compressed_bytes", DataType::UInt64),
        ],
        &["sample", "section_name"],
    )?;
    let layer_selection = SelectionSummary::complete(chain.layers.len() as u64);
    let chrom_selection = SelectionSummary::complete(chain.collection.chroms.len() as u64);
    let archive_selection = SelectionSummary::complete(chain.collection.archives.len() as u64);
    let route_block_count = audits.iter().try_fold(0u64, |sum, audit| {
        audit.shape_routes.iter().try_fold(sum, |sum, inspected| {
            sum.checked_add(inspected.blocks.len() as u64)
                .ok_or(OutputError::InvalidSchema(
                    "inspected shape-route block count exceeds u64".into(),
                ))
        })
    })?;
    let route_block_selection = SelectionSummary::complete(route_block_count);
    let mut bundle = StreamingBundleWriter::new_with_summary(
        writer,
        COLLECTION_INSPECT_RESULT_SCHEMA,
        OutputFormat::from(format),
        context,
        summary,
    )?;
    bundle.write_table("layers", &layer_schema, Some(&layer_selection), |rows| {
        for (index, (layer, audit)) in chain.layers.iter().zip(audits).enumerate() {
            let path = layer.path.to_str().ok_or_else(|| {
                OutputError::Sink(format!(
                    "collection layer path is not valid UTF-8: {}",
                    layer.path.display()
                ))
            })?;
            rows.write_row_with(|row| {
                row.uint64(index as u64)?;
                row.string(path)?;
                row.uint64(layer.file.version as u64)?;
                row.string(&layer.file.root_digest_hex())?;
                row.uint64(layer.manifest.archives.len() as u64)?;
                row.uint64(layer.manifest.junction_count as u64)?;
                row.uint64(
                    layer
                        .file
                        .names()
                        .filter(|name| name.starts_with("j."))
                        .count() as u64,
                )?;
                row.uint64(audit.shape_routes.len() as u64)?;
                row.uint64(audit.section_count() as u64)?;
                row.uint64(audit.span_count() as u64)?;
                row.uint64(audit.pair_count() as u64)?;
                row.uint64(
                    audit
                        .compressed_byte_count()
                        .map_err(|error| OutputError::Sink(error.to_string()))?,
                )?;
                row.boolean(reconstruct_routes)?;
                Ok(())
            })?;
        }
        Ok(())
    })?;
    bundle.write_table(
        "chromosomes",
        &chrom_schema,
        Some(&chrom_selection),
        |rows| {
            for (index, chrom) in chain.collection.chroms.iter().enumerate() {
                rows.write_row_with(|row| {
                    row.uint64(index as u64)?;
                    row.string(chrom)?;
                    Ok(())
                })?;
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "archives",
        &archive_schema,
        Some(&archive_selection),
        |rows| {
            for (global, archive) in chain.collection.archives.iter().enumerate() {
                let (layer_index, local_archive) = chain.global_to_local[global];
                let inspected = audits[layer_index]
                    .shape_routes
                    .iter()
                    .find(|routes| routes.archive_ordinal == local_archive);
                let path = archive.path.to_str().ok_or_else(|| {
                    OutputError::Sink(format!(
                        "archive path is not valid UTF-8: {}",
                        archive.path.display()
                    ))
                })?;
                rows.write_row_with(|row| {
                    row.string(&archive.id)?;
                    row.string(path)?;
                    row.uint64(layer_index as u64)?;
                    row.uint64(local_archive as u64)?;
                    row.uint64(archive.identity.len)?;
                    row.uint64(archive.identity.archive_format_version as u64)?;
                    row.string(&archive.identity.native_scheme)?;
                    row.string(&archive.identity.native_digest)?;
                    match archive.identity.encoded_sections_digest.as_deref() {
                        Some(digest) => row.string(digest)?,
                        None => row.null()?,
                    }
                    row.uint64(archive.chunks.len() as u64)?;
                    match (&archive.shape_routes, inspected) {
                        (Some(binding), Some(inspected)) => {
                            row.boolean(true)?;
                            row.uint64(shaperoute::ROUTE_CODEC_VERSION as u64)?;
                            row.uint64(binding.archive_ordinal as u64)?;
                            row.string(&digest_hex(binding.source_root))?;
                            row.string(&digest_hex(binding.shapes_digest))?;
                            row.uint64(binding.n_shapes as u64)?;
                            row.uint64(inspected.section_count() as u64)?;
                            row.uint64(inspected.span_count() as u64)?;
                            row.uint64(inspected.pair_count() as u64)?;
                            row.uint64(
                                inspected
                                    .compressed_byte_count()
                                    .map_err(|error| OutputError::Sink(error.to_string()))?,
                            )?;
                            row.boolean(reconstruct_routes)?;
                        }
                        (None, None) => {
                            row.boolean(false)?;
                            row.null()?;
                            row.null()?;
                            row.null()?;
                            row.null()?;
                            row.null()?;
                            row.uint64(0)?;
                            row.uint64(0)?;
                            row.uint64(0)?;
                            row.uint64(0)?;
                            row.boolean(false)?;
                        }
                        _ => {
                            return Err(OutputError::Sink(
                                "flattened collection route binding differs from its owning layer"
                                    .into(),
                            ));
                        }
                    }
                    Ok(())
                })?;
            }
            Ok(())
        },
    )?;
    bundle.write_table(
        "shape_route_blocks",
        &route_block_schema,
        Some(&route_block_selection),
        |rows| {
            for (layer_index, audit) in audits.iter().enumerate() {
                for inspected in &audit.shape_routes {
                    let archive =
                        &chain.layers[layer_index].manifest.archives[inspected.archive_ordinal];
                    let binding = archive.shape_routes.as_ref().ok_or_else(|| {
                        OutputError::Sink(
                            "inspected route archive lacks its manifest binding".into(),
                        )
                    })?;
                    for ((descriptor, block), compressed_bytes) in binding
                        .blocks
                        .iter()
                        .zip(&inspected.blocks)
                        .zip(&inspected.compressed_bytes)
                    {
                        rows.write_row_with(|row| {
                            row.string(&archive.id)?;
                            row.string(&descriptor.section_name)?;
                            row.uint64(descriptor.first_span as u64)?;
                            row.uint64(descriptor.last_span as u64)?;
                            row.uint64(block.spans.len() as u64)?;
                            row.uint64(block.pair_count() as u64)?;
                            row.uint64(*compressed_bytes)?;
                            Ok(())
                        })?;
                    }
                }
            }
            Ok(())
        },
    )?;
    bundle.finish()?;
    Ok(())
}

fn run_inspect(
    path: PathBuf,
    verify_content: bool,
    verify_routes: bool,
    uniform_output: CollectionOutputArgs,
) -> Result<()> {
    let started = std::time::Instant::now();
    let chain = open_collection_chain(&path)?;
    let collection = &chain.collection;
    let uniform_context = if uniform_output.format.is_some() {
        let mut parameters = BTreeMap::new();
        parameters.insert("verify_content".into(), json!(verify_content));
        parameters.insert("verify_routes".into(), json!(verify_routes));
        parameters.insert(
            "access".into(),
            json!("authenticated collection directories and all indexed payloads; source content only when explicitly verified"),
        );
        Some(collection_uniform_context(&path, &chain, parameters)?)
    } else {
        None
    };
    let mut layer_audits = Vec::with_capacity(chain.layers.len());
    for layer in &chain.layers {
        layer_audits.push(audit_collection_layer(layer)?);
    }
    let source_io = validate_collection_sources(collection, verify_content)?;
    let reconstruct_routes = verify_routes || verify_content;
    if reconstruct_routes {
        for (layer, audit) in chain.layers.iter().zip(&mut layer_audits) {
            for inspected in &mut audit.shape_routes {
                let archive = &layer.manifest.archives[inspected.archive_ordinal];
                let binding = archive
                    .shape_routes
                    .as_ref()
                    .context("inspected route archive lacks its manifest binding")?;
                inspected.verification_source_bytes = Some(verify_inspected_archive_routes(
                    archive,
                    binding,
                    &inspected.blocks,
                )?);
            }
        }
    }
    let junction_segments: usize = chain
        .layers
        .iter()
        .map(|layer| {
            layer
                .file
                .names()
                .filter(|name| name.starts_with("j."))
                .count()
        })
        .sum();
    let file_bytes: u64 = chain.layers.iter().try_fold(0u64, |sum, layer| {
        sum.checked_add(std::fs::metadata(&layer.path)?.len())
            .context("collection chain byte count overflow")
    })?;
    let shape_route_archives: usize = layer_audits
        .iter()
        .map(|audit| audit.shape_routes.len())
        .sum();
    let shape_route_sections: usize = layer_audits
        .iter()
        .map(CollectionLayerAudit::section_count)
        .sum();
    let shape_route_exact_spans: usize = layer_audits
        .iter()
        .map(CollectionLayerAudit::span_count)
        .sum();
    let shape_route_pairs: usize = layer_audits
        .iter()
        .map(CollectionLayerAudit::pair_count)
        .sum();
    let shape_route_compressed_bytes = layer_audits.iter().try_fold(0u64, |sum, audit| {
        sum.checked_add(audit.compressed_byte_count()?)
            .context("shape-route compressed byte count overflow")
    })?;
    let route_verification_source_bytes_read =
        layer_audits.iter().try_fold(0u64, |sum, audit| {
            sum.checked_add(audit.verification_source_byte_count()?)
                .context("shape-route verification source byte count overflow")
        })?;
    if let Some(format) = uniform_output.format {
        let (algo, digest) = collection
            .genome_algo
            .as_deref()
            .zip(collection.genome_digest.as_deref())
            .map_or((None, None), |(algo, digest)| (Some(algo), Some(digest)));
        let summary = InspectUniformSummary {
            collection_schema: SCHEMA,
            format_version: chain.layers.last().map(|layer| layer.file.version),
            path: uniform_path(&path, "collection")?,
            file_bytes,
            reference: InspectReferenceSummary {
                stamped: algo.is_some(),
                algo,
                digest,
            },
            index: InspectIndexSummary {
                segment_junction_rows: collection.junction_count as u64,
                archive_routes: collection.route_count as u64,
                chunk_postings: collection.posting_count as u64,
                interval_chunks: collection
                    .archives
                    .iter()
                    .map(|archive| archive.chunks.len() as u64)
                    .sum(),
                junction_bin_bp: JUNCTION_BIN_BP,
                junction_segments: junction_segments as u64,
                shape_route_archives: shape_route_archives as u64,
                shape_route_sections: shape_route_sections as u64,
                shape_route_exact_spans: shape_route_exact_spans as u64,
                shape_route_pairs: shape_route_pairs as u64,
                shape_route_compressed_bytes,
            },
            guard: InspectGuardSummary {
                filesystem_identity: "size + nanosecond mtime + nanosecond ctime + device + inode",
                content_digest_recorded: true,
                content_digest_verified: verify_content,
                shape_route_payloads_verified: true,
                shape_route_reconstruction_verified: reconstruct_routes,
            },
            io: InspectIoSummary {
                collection_sidecar_bytes_read: collection_sidecar_bytes_read(&chain)?,
                source_identity_content_bytes_read: identity_content_bytes(&source_io)?,
                source_identity_total_bytes_read: source_archive_bytes_read(&source_io)?,
                route_verification_source_bytes_read,
            },
        };
        let context = uniform_context
            .as_ref()
            .expect("uniform collection inspect context was prepared");
        write_uniform_collection_result(&uniform_output, |writer| {
            stream_uniform_inspect(
                writer,
                format,
                context,
                &summary,
                &chain,
                &layer_audits,
                reconstruct_routes,
            )
        })?;
        return Ok(());
    }
    let value = json!({
        "schema": SCHEMA,
        "format_version": chain.layers.last().map(|layer| layer.file.version),
        "path": path,
        "file_bytes": file_bytes,
        "layers": chain.layers.iter().zip(&layer_audits).map(|(layer, audit)| -> Result<_> {
            Ok(json!({
                "path": layer.path,
                "format_version": layer.file.version,
                "root_digest": layer.file.root_digest_hex(),
                "archives": layer.manifest.archives.len(),
                "junction_rows": layer.manifest.junction_count,
                "junction_segments": layer.file.names().filter(|name| name.starts_with("j.")).count(),
                "shape_routes": {
                    "archives": audit.shape_routes.len(),
                    "sections": audit.section_count(),
                    "exact_spans": audit.span_count(),
                    "pairs": audit.pair_count(),
                    "compressed_bytes": audit.compressed_byte_count()?,
                    "reconstruction_verified": reconstruct_routes,
                },
            }))
        }).collect::<Result<Vec<_>>>()?,
        "reference": match (&collection.genome_algo, &collection.genome_digest) {
            (Some(algo), Some(digest)) => json!({"stamped": true, "algo": algo, "digest": digest}),
            _ => json!({"stamped": false}),
        },
        "chromosomes": collection.chroms,
        "archives": collection.archives.iter().enumerate().map(|(global, archive)| -> Result<_> {
            let (layer_index, local_archive) = chain.global_to_local[global];
            let inspected = layer_audits[layer_index]
                .shape_routes
                .iter()
                .find(|routes| routes.archive_ordinal == local_archive);
            let shape_routes = match (&archive.shape_routes, inspected) {
                (Some(binding), Some(inspected)) => {
                    Some(route_binding_json(binding, inspected, reconstruct_routes)?)
                }
                (None, None) => None,
                _ => bail!("flattened collection route binding differs from its owning layer"),
            };
            Ok(json!({
                "id": archive.id,
                "path": archive.path,
                "bytes": archive.identity.len,
                "archive_format_version": archive.identity.archive_format_version,
                "native_identity": {
                    "scheme": archive.identity.native_scheme,
                    "blake3": archive.identity.native_digest,
                },
                "encoded_sections_identity": archive.identity.encoded_sections_digest.as_ref().map(|digest| json!({
                    "scheme": "aie-encoded-sections-v1",
                    "blake3": digest,
                })),
                "chunks": archive.chunks.len(),
                "shape_routes": shape_routes,
            }))
        }).collect::<Result<Vec<_>>>()?,
        "index": {
            "segment_junction_rows": collection.junction_count,
            "global_junctions": serde_json::Value::Null,
            "archive_routes": collection.route_count,
            "chunk_postings": collection.posting_count,
            "interval_chunks": collection.archives.iter().map(|archive| archive.chunks.len()).sum::<usize>(),
            "junction_bin_bp": JUNCTION_BIN_BP,
            "junction_segments": junction_segments,
            "shape_route_archives": shape_route_archives,
            "shape_route_sections": shape_route_sections,
            "shape_route_exact_spans": shape_route_exact_spans,
            "shape_route_pairs": shape_route_pairs,
            "shape_route_compressed_bytes": shape_route_compressed_bytes,
        },
        "guard": {
            "filesystem_identity": "size + nanosecond mtime + nanosecond ctime + device + inode",
            "content_digest_recorded": true,
            "content_digest_verified": verify_content,
            "shape_route_payloads_verified": true,
            "shape_route_reconstruction_verified": reconstruct_routes,
        },
        "io": {
            "collection_sidecar_bytes_read": collection_sidecar_bytes_read(&chain)?,
            "source_identity_content_bytes_read": identity_content_bytes(&source_io)?,
            "source_identity_total_bytes_read": source_archive_bytes_read(&source_io)?,
            "route_verification_source_bytes_read": route_verification_source_bytes_read,
        },
        "elapsed_seconds": started.elapsed().as_secs_f64(),
    });
    println!("{}", serde_json::to_string_pretty(&value)?);
    Ok(())
}

fn validate_collection_args(what: &What) -> Result<()> {
    let output = match what {
        What::Build { uniform_output, .. }
        | What::Inspect { uniform_output, .. }
        | What::Junction { uniform_output, .. }
        | What::Region { uniform_output, .. }
        | What::Jset { uniform_output, .. } => uniform_output,
    };
    preflight_uniform_collection_output(output)?;
    if output.format.is_some() {
        match what {
            What::Build {
                base,
                out,
                uniform_output,
                ..
            } => {
                uniform_path(out, "collection output")?;
                if let Some(base) = base {
                    uniform_path(base, "base collection")?;
                }
                if let Some(report) = uniform_output.output.as_deref() {
                    if destination_key(out)? == destination_key(report)? {
                        bail!(
                            "collection --out and uniform --output must be different destinations"
                        );
                    }
                }
            }
            What::Inspect { collection, .. }
            | What::Junction { collection, .. }
            | What::Region { collection, .. }
            | What::Jset { collection, .. } => {
                uniform_path(collection, "collection")?;
            }
        }
    }
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    validate_collection_args(&args.what)?;
    match args.what {
        What::Build {
            samples,
            source_digests,
            base,
            out,
            allow_unstamped,
            shape_routes,
            json,
            uniform_output,
        } => run_build(BuildRun {
            samples,
            source_digests,
            base,
            out,
            allow_unstamped,
            build_shape_routes: shape_routes,
            json_output: json,
            uniform_output,
        }),
        What::Inspect {
            collection,
            verify_content,
            verify_routes,
            uniform_output,
        } => run_inspect(collection, verify_content, verify_routes, uniform_output),
        What::Junction {
            collection,
            locus,
            min_support,
            top,
            explain,
            verify_content,
            json,
            uniform_output,
        } => run_junction_query(JunctionQueryRun {
            path: collection,
            locus,
            min_support,
            top,
            explain,
            verify_content,
            json_output: json,
            uniform_output,
        }),
        What::Region {
            collection,
            locus,
            top,
            explain,
            verify_content,
            json,
            uniform_output,
        } => run_region_query(
            collection,
            locus,
            top,
            explain,
            verify_content,
            json,
            uniform_output,
        ),
        What::Jset {
            collection,
            include,
            exclude,
            min_support,
            top,
            explain,
            verify_content,
            json,
            uniform_output,
        } => run_jset_query(JsetQueryRun {
            path: collection,
            includes: include,
            excludes: exclude,
            min_support,
            top,
            explain,
            verify_content,
            json_output: json,
            uniform_output,
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    fn test_path(tag: &str, extension: &str) -> PathBuf {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        std::env::temp_dir().join(format!(
            "gravlax-collection-{tag}-{}-{}.{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed),
            extension
        ))
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn held_collection_output_ignores_and_preserves_replaced_staging_path() {
        let temporary = test_path("held-stage", "tmp");
        let displaced = test_path("held-displaced", "tmp");
        let destination = test_path("held-destination", "aicollection");
        let mut output = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&temporary)
            .unwrap();
        output.write_all(b"produced collection bytes\n").unwrap();
        output.flush().unwrap();
        std::fs::rename(&temporary, &displaced).unwrap();
        std::fs::write(&temporary, b"concurrent replacement\n").unwrap();

        let bytes = install_collection_output(&output, &temporary, &destination).unwrap();
        assert_eq!(bytes, b"produced collection bytes\n".len() as u64);
        assert_eq!(
            std::fs::read(&destination).unwrap(),
            b"produced collection bytes\n"
        );
        assert_eq!(
            std::fs::read(&temporary).unwrap(),
            b"concurrent replacement\n"
        );
        assert_eq!(
            std::fs::read(&displaced).unwrap(),
            b"produced collection bytes\n"
        );

        for path in [temporary, displaced, destination] {
            std::fs::remove_file(path).ok();
        }
    }

    fn fixture_archive() -> ArchiveEntry {
        ArchiveEntry {
            id: "A".into(),
            path: PathBuf::from("A.aie"),
            identity: FileIdentity {
                len: 1,
                modified_secs: 0,
                modified_nanos: 0,
                changed_secs: 0,
                changed_nanos: 0,
                dev: 1,
                inode: 1,
                archive_format_version: evidence_io::format::VERSION,
                native_scheme: ROOTED_DIRECTORY_SCHEME.into(),
                native_digest: "0".repeat(64),
                encoded_sections_digest: Some("1".repeat(64)),
            },
            chunks: vec![
                IndexedChunk {
                    info: ChunkInfo {
                        chrom: 0,
                        bin_start: 0,
                        n_mols: 1,
                        class_base: 0,
                        max_anchor: 100,
                        n_cells: 1,
                    },
                    compressed_bytes: 1,
                },
                IndexedChunk {
                    info: ChunkInfo {
                        chrom: 1,
                        bin_start: 0,
                        n_mols: 1,
                        class_base: 1,
                        max_anchor: 100,
                        n_cells: 1,
                    },
                    compressed_bytes: 1,
                },
            ],
            shape_routes: None,
        }
    }

    fn fixture_row(posts: Vec<u32>) -> GlobalJunction {
        GlobalJunction {
            chrom: 0,
            donor: 100,
            acceptor: 200,
            support_upper_bound: 1,
            routes: vec![ArchiveRoute {
                archive: 0,
                supporting_children: 1,
                posts,
            }],
        }
    }

    fn fixture_collection() -> Collection {
        Collection {
            base: None,
            genome_algo: Some("blake3".into()),
            genome_digest: Some("2".repeat(64)),
            chroms: vec!["chr1".into(), "chr2".into()],
            chroms_digest: "3".repeat(64),
            archives: vec![fixture_archive()],
            junctions: Vec::new(),
            junction_count: 0,
            route_count: 0,
            posting_count: 0,
            shape_route_blocks: Vec::new(),
            encoded_shape_route_blocks: Vec::new(),
        }
    }

    fn scientific_archive_fixture() -> PathBuf {
        use crate::rows::{Extracted, MolChain, MolRec};
        use evidence_io::archive::Shape;
        use smallvec::smallvec;

        let chain = |weight, reps| MolChain { weight, reps };
        let molecule = |cell, umi_class, strand_rev, chains, mms| MolRec {
            cell,
            umi_class,
            chrom: 0,
            strand_rev,
            chains,
            mms,
        };
        // Shape 0 exposes the same 10-base intron at donor offsets 10 and 25. This lets the
        // fixture exercise both span-extreme representatives without making either one special.
        let extracted = Extracted {
            mols: vec![
                molecule(0, 0, false, smallvec![chain(1, smallvec![(85, 0)])], smallvec![]),
                molecule(
                    1,
                    1,
                    false,
                    smallvec![chain(2, smallvec![(85, 0), (100, 0)])],
                    smallvec![],
                ),
                molecule(2, 2, true, smallvec![chain(1, smallvec![(100, 0)])], smallvec![]),
                molecule(3, 3, false, smallvec![], smallvec![(100, 0, 0, 1)]),
                molecule(4, 4, false, smallvec![chain(1, smallvec![(100, 2)])], smallvec![]),
                // Class 0 reappears in a later chunk and contributes the exclusion junction.
                molecule(0, 0, false, smallvec![chain(1, smallvec![(205, 3)])], smallvec![]),
                molecule(5, 5, false, smallvec![chain(1, smallvec![(205, 3)])], smallvec![]),
            ],
            edges: Vec::new(),
            cells: vec![0, 1, 2, 3, 4, 5],
            shapes: vec![
                Shape {
                    blocks: vec![(0, 10), (20, 5), (35, 5)],
                },
                Shape {
                    blocks: vec![(0, 8), (18, 4)],
                },
                Shape {
                    blocks: vec![(0, 20)],
                },
                Shape {
                    blocks: vec![(0, 5), (25, 5)],
                },
            ],
            patterns: vec![Vec::new()],
            n_classes: 6,
            chrom_names: vec!["chr1".into()],
        };
        let path = test_path("scientific-routes", "aie");
        crate::archivecmd::write_archive(&extracted, &path, 1, 100, None).unwrap();
        path
    }

    fn routed_collection_fixture(source: &Path) -> (Collection, Vec<u8>) {
        let mut reader = evidence_io::format::SectionReader::open(source).unwrap();
        assert_eq!(reader.archive_version(), evidence_io::format::VERSION);
        let root = reader.content_commitment().unwrap();
        let encoded = digest_hex(reader.encoded_content_identity().unwrap().unwrap());
        let shapes_digest = reader
            .section_metadata()
            .find(|section| section.name == "shapes")
            .unwrap()
            .compressed_blake3
            .unwrap();
        let shapes_raw = reader.read("shapes").unwrap();
        let chunks = read_chunk_index(&mut reader).unwrap();
        let indexed_chunks = chunks
            .into_iter()
            .enumerate()
            .map(|(index, info)| IndexedChunk {
                info,
                compressed_bytes: reader
                    .section_metadata()
                    .find(|section| section.name == format!("c{index}"))
                    .unwrap()
                    .compressed_len,
            })
            .collect();
        let identity = identity_from_metadata(
            reader.file_metadata().unwrap(),
            evidence_io::format::VERSION,
            ROOTED_DIRECTORY_SCHEME.into(),
            root.to_hex(),
            Some(encoded),
        )
        .unwrap();
        let derived = shaperoute::derive_from_shapes(&shapes_raw, 0).unwrap();
        let binding = derived.binding(root.digest, shapes_digest).unwrap();
        let chrom_bytes = b"chr1";
        (
            Collection {
                base: None,
                genome_algo: None,
                genome_digest: None,
                chroms: vec!["chr1".into()],
                chroms_digest: blake3::hash(chrom_bytes).to_hex().to_string(),
                archives: vec![ArchiveEntry {
                    id: "A".into(),
                    path: source.to_path_buf(),
                    identity,
                    chunks: indexed_chunks,
                    shape_routes: Some(binding),
                }],
                junctions: Vec::new(),
                junction_count: 0,
                route_count: 0,
                posting_count: 0,
                shape_route_blocks: derived.blocks,
                encoded_shape_route_blocks: Vec::new(),
            },
            shapes_raw,
        )
    }

    fn minimal_local_route_fixture(
        id: &str,
        span: u32,
        archive_ordinal: u32,
    ) -> (PathBuf, LocalArchive) {
        let path = test_path(&format!("local-route-{id}"), "aie");
        let mut shapes_raw = Vec::new();
        put_varint(&mut shapes_raw, 2);
        put_varint(&mut shapes_raw, 0);
        put_varint(&mut shapes_raw, 5);
        put_varint(&mut shapes_raw, (5 + span) as u64);
        put_varint(&mut shapes_raw, 5);
        let mut writer = evidence_io::format::SectionWriter::create_new(&path, 1).unwrap();
        writer.section("shapes", &shapes_raw).unwrap();
        writer.section("fixture.id", id.as_bytes()).unwrap();
        writer.finish().unwrap();

        let reader = evidence_io::format::SectionReader::open(&path).unwrap();
        let root = reader.content_commitment().unwrap();
        let encoded = digest_hex(reader.encoded_content_identity().unwrap().unwrap());
        let shapes_meta = reader
            .section_metadata()
            .find(|section| section.name == "shapes")
            .unwrap();
        let identity = identity_from_metadata(
            reader.file_metadata().unwrap(),
            evidence_io::format::VERSION,
            ROOTED_DIRECTORY_SCHEME.into(),
            root.to_hex(),
            Some(encoded),
        )
        .unwrap();
        let shape_routes = LocalShapeRoutes {
            shapes_digest: shapes_meta.compressed_blake3.unwrap(),
            compressed_bytes_read: shapes_meta.compressed_len,
            derived: shaperoute::derive_encoded_from_shapes(&shapes_raw, archive_ordinal).unwrap(),
        };
        let local = LocalArchive {
            entry: ArchiveEntry {
                id: id.into(),
                path: path.clone(),
                identity: identity.clone(),
                chunks: Vec::new(),
                shape_routes: None,
            },
            source_io: SourceIo {
                id: id.into(),
                format_version: identity.archive_format_version,
                identity_scheme: identity.native_scheme.clone(),
                identity_content_bytes_read: 0,
                total_bytes_read: 0,
                shape_route_payload_bytes_read: shape_routes.compressed_bytes_read,
                shape_route_source_bytes_read: shape_routes.compressed_bytes_read,
                sections_read: vec!["shapes".into()],
            },
            genome_algo: None,
            genome_digest: None,
            chroms: vec!["chr1".into()],
            chroms_digest: blake3::hash(b"chr1").to_hex().to_string(),
            junctions: Vec::new(),
            shape_routes: Some(shape_routes),
        };
        (path, local)
    }

    fn encode_v2_manifest_fixture(collection: &Collection, digest: &str) -> Vec<u8> {
        let mut out = Vec::new();
        out.push(0); // no base
        out.push(u8::from(collection.genome_algo.is_some()));
        if let Some(value) = &collection.genome_algo {
            put_string(&mut out, value);
        }
        out.push(u8::from(collection.genome_digest.is_some()));
        if let Some(value) = &collection.genome_digest {
            put_string(&mut out, value);
        }
        put_varint(&mut out, collection.chroms.len() as u64);
        for chrom in &collection.chroms {
            put_string(&mut out, chrom);
        }
        put_string(&mut out, &collection.chroms_digest);
        put_varint(&mut out, collection.archives.len() as u64);
        for archive in &collection.archives {
            put_string(&mut out, &archive.id);
            put_string(&mut out, archive.path.to_str().unwrap());
            for value in [
                archive.identity.len,
                archive.identity.modified_secs,
                archive.identity.modified_nanos as u64,
                archive.identity.changed_secs,
                archive.identity.changed_nanos as u64,
                archive.identity.dev,
                archive.identity.inode,
            ] {
                put_varint(&mut out, value);
            }
            put_string(&mut out, digest);
            put_varint(&mut out, archive.chunks.len() as u64);
            for chunk in &archive.chunks {
                for value in [
                    chunk.info.chrom as u64,
                    chunk.info.bin_start as u64,
                    chunk.info.n_mols as u64,
                    chunk.info.class_base as u64,
                    chunk.info.max_anchor as u64,
                    chunk.info.n_cells as u64,
                    chunk.compressed_bytes,
                ] {
                    put_varint(&mut out, value);
                }
            }
        }
        for value in [
            collection.junction_count,
            collection.route_count,
            collection.posting_count,
            0,
        ] {
            put_varint(&mut out, value as u64);
        }
        out
    }

    fn encode_v3_manifest_fixture(collection: &Collection) -> Vec<u8> {
        assert!(collection.base.is_none());
        assert!(collection.junctions.is_empty());
        assert!(collection
            .archives
            .iter()
            .all(|archive| archive.shape_routes.is_none()));
        let mut out = Vec::new();
        out.push(0); // no base
        put_optional_string(&mut out, &collection.genome_algo);
        put_optional_string(&mut out, &collection.genome_digest);
        put_varint(&mut out, collection.chroms.len() as u64);
        for chrom in &collection.chroms {
            put_string(&mut out, chrom);
        }
        put_string(&mut out, &collection.chroms_digest);
        put_varint(&mut out, collection.archives.len() as u64);
        for archive in &collection.archives {
            put_string(&mut out, &archive.id);
            put_string(&mut out, archive.path.to_str().unwrap());
            for value in [
                archive.identity.len,
                archive.identity.modified_secs,
                archive.identity.modified_nanos as u64,
                archive.identity.changed_secs,
                archive.identity.changed_nanos as u64,
                archive.identity.dev,
                archive.identity.inode,
                archive.identity.archive_format_version as u64,
            ] {
                put_varint(&mut out, value);
            }
            put_string(&mut out, &archive.identity.native_scheme);
            put_string(&mut out, &archive.identity.native_digest);
            put_string(
                &mut out,
                archive.identity.encoded_sections_digest.as_ref().unwrap(),
            );
            put_varint(&mut out, archive.chunks.len() as u64);
            for chunk in &archive.chunks {
                for value in [
                    chunk.info.chrom as u64,
                    chunk.info.bin_start as u64,
                    chunk.info.n_mols as u64,
                    chunk.info.class_base as u64,
                    chunk.info.max_anchor as u64,
                    chunk.info.n_cells as u64,
                    chunk.compressed_bytes,
                ] {
                    put_varint(&mut out, value);
                }
            }
        }
        for value in [
            collection.junction_count,
            collection.route_count,
            collection.posting_count,
            0,
        ] {
            put_varint(&mut out, value as u64);
        }
        out
    }

    fn write_raw_collection_fixture(path: &Path, version: u32, sections: &[(&str, &[u8])]) {
        let compressed = sections
            .iter()
            .map(|(name, raw)| (*name, *raw, evidence_io::format::compress(raw, 1).unwrap()))
            .collect::<Vec<_>>();
        let mut header = Vec::new();
        header.extend_from_slice(MAGIC);
        header.extend_from_slice(&version.to_le_bytes());
        header.extend_from_slice(&(sections.len() as u32).to_le_bytes());
        let mut directory = Vec::new();
        for (name, raw, bytes) in &compressed {
            directory.push(name.len() as u8);
            directory.extend_from_slice(name.as_bytes());
            directory.extend_from_slice(&(raw.len() as u64).to_le_bytes());
            directory.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
            directory.extend_from_slice(blake3::hash(raw).as_bytes());
        }
        let mut hasher = blake3::Hasher::new();
        hasher.update(&header);
        hasher.update(&directory);
        let root = hasher.finalize();
        let mut bytes = header;
        bytes.extend_from_slice(root.as_bytes());
        bytes.extend_from_slice(&directory);
        for (_, _, compressed) in compressed {
            bytes.extend_from_slice(&compressed);
        }
        std::fs::write(path, bytes).unwrap();
    }

    fn write_manifest_collection_fixture(path: &Path, manifest: &[u8], version: u32) {
        write_raw_collection_fixture(path, version, &[("manifest", manifest)]);
    }

    #[test]
    fn v4_manifest_roundtrips_identities_and_v2_manifest_remains_readable() {
        let collection = fixture_collection();
        let raw = encode_collection(&collection).unwrap();
        assert_eq!(decode_collection(&raw, VERSION).unwrap(), collection);

        let v3 = encode_v3_manifest_fixture(&collection);
        assert_eq!(
            decode_collection(&v3, PREVIOUS_VERSION).unwrap(),
            collection
        );
        let v3_path = test_path("v3-sidecar", "aicollection");
        write_manifest_collection_fixture(&v3_path, &v3, PREVIOUS_VERSION);
        let (v3_file, v3_decoded) = open_collection_manifest(&v3_path).unwrap();
        assert_eq!(v3_file.version, PREVIOUS_VERSION);
        assert_eq!(v3_decoded, collection);
        drop(v3_file);
        std::fs::remove_file(v3_path).unwrap();

        // Collection v2 accepted upper- or lower-case hexadecimal. Decode to the canonical
        // lower-case representation used by current duplicate and content guards.
        let legacy_digest = "A".repeat(64);
        let raw = encode_v2_manifest_fixture(&collection, &legacy_digest);
        let decoded = decode_collection(&raw, LEGACY_VERSION).unwrap();
        let identity = &decoded.archives[0].identity;
        assert_eq!(
            identity.archive_format_version,
            evidence_io::format::SEEKABLE_VERSION
        );
        assert_eq!(identity.native_scheme, LEGACY_FULL_FILE_SCHEME);
        assert_eq!(identity.native_digest, legacy_digest.to_ascii_lowercase());
        assert_eq!(identity.encoded_sections_digest, None);

        let path = test_path("legacy-v2-sidecar", "aicollection");
        write_manifest_collection_fixture(&path, &raw, LEGACY_VERSION);
        let (file, decoded_file) = open_collection_manifest(&path).unwrap();
        assert_eq!(file.version, LEGACY_VERSION);
        assert_eq!(decoded_file, decoded);
        drop(file);
        let chain = open_collection_chain(&path).unwrap();
        assert!(require_encoded_extension_base(&chain)
            .unwrap_err()
            .to_string()
            .contains("cannot be used as a v4 extension base"));
        drop(chain);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn v3_identity_contract_rejects_wrong_scheme_digest_and_missing_encoded_identity() {
        let mut identity = fixture_archive().identity;
        assert!(validate_identity(&identity, "A", true).is_ok());
        identity.native_scheme = LEGACY_FULL_FILE_SCHEME.into();
        assert!(validate_identity(&identity, "A", true)
            .unwrap_err()
            .to_string()
            .contains("requires native identity scheme"));
        identity = fixture_archive().identity;
        identity.native_digest = "not-a-digest".into();
        assert!(validate_identity(&identity, "A", true).is_err());
        identity = fixture_archive().identity;
        identity.encoded_sections_digest = None;
        assert!(validate_identity(&identity, "A", true)
            .unwrap_err()
            .to_string()
            .contains("lacks the encoded-sections digest"));
    }

    #[test]
    fn rooted_identity_guard_uses_no_payload_bytes_and_rechecks_root_and_encoded_identity() {
        let path = test_path("rooted-source", "aie");
        let mut writer = evidence_io::format::SectionWriter::create_new(&path, 1).unwrap();
        writer.section("x", b"root-bound payload").unwrap();
        writer.finish().unwrap();
        let reader = evidence_io::format::SectionReader::open(&path).unwrap();
        let root = reader.content_commitment().unwrap().to_hex();
        let encoded = digest_hex(reader.encoded_content_identity().unwrap().unwrap());
        let identity = identity_from_metadata(
            reader.file_metadata().unwrap(),
            evidence_io::format::VERSION,
            ROOTED_DIRECTORY_SCHEME.into(),
            root,
            Some(encoded),
        )
        .unwrap();
        drop(reader);
        let valid_identity = identity.clone();
        let mut archive = ArchiveEntry {
            id: "rooted".into(),
            path: path.clone(),
            identity,
            chunks: Vec::new(),
            shape_routes: None,
        };
        let io = validate_archive_identity(&archive, false).unwrap();
        assert_eq!(io.identity_content_bytes_read, 0);
        assert!(io.total_bytes_read > 0);
        assert!(io.sections_read.is_empty());

        archive.identity.native_digest = "f".repeat(64);
        assert!(validate_archive_identity(&archive, false)
            .unwrap_err()
            .to_string()
            .contains("root changed"));
        archive.identity = valid_identity;
        archive.identity.encoded_sections_digest = Some("c".repeat(64));
        assert!(validate_archive_identity(&archive, false)
            .unwrap_err()
            .to_string()
            .contains("encoded-sections identity changed"));
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn identity_validation_keeps_stat_and_root_on_the_same_open_file() {
        let path = test_path("rooted-open-fd-path", "aie");
        let moved = test_path("rooted-open-fd-held", "aie");
        let mut writer = evidence_io::format::SectionWriter::create_new(&path, 1).unwrap();
        writer.section("x", b"original root-bound payload").unwrap();
        writer.finish().unwrap();

        let reader = evidence_io::format::SectionReader::open(&path).unwrap();
        let root = reader.content_commitment().unwrap().to_hex();
        let encoded = digest_hex(reader.encoded_content_identity().unwrap().unwrap());
        drop(reader);

        let held = std::fs::File::open(&path).unwrap();
        std::fs::rename(&path, &moved).unwrap();
        let mut replacement = evidence_io::format::SectionWriter::create_new(&path, 1).unwrap();
        replacement
            .section("x", b"replacement payload at the old path")
            .unwrap();
        replacement.finish().unwrap();

        // Capture the held inode's post-rename stat tuple. Validation must continue on this open
        // file description rather than following `archive.path` to the replacement.
        let identity = identity_from_metadata(
            held.metadata().unwrap(),
            evidence_io::format::VERSION,
            ROOTED_DIRECTORY_SCHEME.into(),
            root,
            Some(encoded),
        )
        .unwrap();
        let archive = ArchiveEntry {
            id: "held".into(),
            path: path.clone(),
            identity,
            chunks: Vec::new(),
            shape_routes: None,
        };
        let io = validate_archive_identity_file(&archive, held, false).unwrap();
        assert_eq!(io.identity_content_bytes_read, 0);
        assert!(io.total_bytes_read > 0);

        // A fresh validation follows the path and therefore rejects the replacement inode.
        assert!(validate_archive_identity(&archive, false).is_err());
        std::fs::remove_file(path).unwrap();
        std::fs::remove_file(moved).unwrap();
    }

    fn fixture_local(id: &str, path: &str, identity: FileIdentity) -> LocalArchive {
        LocalArchive {
            entry: ArchiveEntry {
                id: id.into(),
                path: PathBuf::from(path),
                identity: identity.clone(),
                chunks: Vec::new(),
                shape_routes: None,
            },
            source_io: SourceIo {
                id: id.into(),
                format_version: identity.archive_format_version,
                identity_scheme: identity.native_scheme.clone(),
                identity_content_bytes_read: 0,
                total_bytes_read: 0,
                shape_route_payload_bytes_read: 0,
                shape_route_source_bytes_read: 0,
                sections_read: Vec::new(),
            },
            genome_algo: None,
            genome_digest: None,
            chroms: vec!["chr1".into()],
            chroms_digest: "4".repeat(64),
            junctions: Vec::new(),
            shape_routes: None,
        }
    }

    #[test]
    fn duplicate_guard_uses_encoded_identity_across_legacy_and_rooted_containers() {
        let encoded = Some("e".repeat(64));
        let legacy = FileIdentity {
            len: 1,
            modified_secs: 0,
            modified_nanos: 0,
            changed_secs: 0,
            changed_nanos: 0,
            dev: 1,
            inode: 1,
            archive_format_version: evidence_io::format::SEEKABLE_VERSION,
            native_scheme: LEGACY_FULL_FILE_SCHEME.into(),
            native_digest: "a".repeat(64),
            encoded_sections_digest: encoded.clone(),
        };
        let rooted = FileIdentity {
            len: 2,
            modified_secs: 0,
            modified_nanos: 0,
            changed_secs: 0,
            changed_nanos: 0,
            dev: 2,
            inode: 2,
            archive_format_version: evidence_io::format::VERSION,
            native_scheme: ROOTED_DIRECTORY_SCHEME.into(),
            native_digest: "b".repeat(64),
            encoded_sections_digest: encoded,
        };
        let error = assemble_collection(
            vec![
                fixture_local("legacy", "legacy.aie", legacy),
                fixture_local("rooted", "rooted.aie", rooted),
            ],
            true,
            2,
        )
        .unwrap_err();
        assert!(error.to_string().contains("identical archive content"));
    }

    #[test]
    fn routed_v4_manifest_sections_are_deterministic_and_fully_reconstructable() {
        let source = scientific_archive_fixture();
        let (collection, shapes_raw) = routed_collection_fixture(&source);
        let binding = collection.archives[0].shape_routes.as_ref().unwrap();

        let manifest = encode_collection(&collection).unwrap();
        let decoded = decode_collection(&manifest, VERSION).unwrap();
        assert_eq!(decoded.archives[0].shape_routes.as_ref(), Some(binding));
        assert!(decoded.shape_route_blocks.is_empty());
        assert_eq!(
            decoded.archives[0].shape_routes.as_ref().unwrap().blocks[0].section_name,
            binding.blocks[0].section_name
        );

        let first = test_path("routed-v4-first", "aicollection");
        let second = test_path("routed-v4-second", "aicollection");
        let first_stats = write_collection(&first, &collection).unwrap();
        let second_stats = write_collection(&second, &collection).unwrap();
        assert_eq!(std::fs::read(&first).unwrap(), std::fs::read(&second).unwrap());
        assert_eq!(first_stats.shape_route_sections, binding.blocks.len());
        assert_eq!(first_stats, second_stats);

        let chain = open_collection_chain(&first).unwrap();
        assert_eq!(chain.layers.len(), 1);
        assert_eq!(chain.layers[0].file.version, VERSION);
        let audit = audit_collection_layer(&chain.layers[0]).unwrap();
        assert_eq!(audit.shape_routes.len(), 1);
        assert_eq!(audit.section_count(), binding.blocks.len());
        assert_eq!(audit.span_count(), 2);
        assert_eq!(audit.pair_count(), 4);
        assert_eq!(audit.shape_routes[0].blocks, collection.shape_route_blocks);
        let route_json = route_binding_json(binding, &audit.shape_routes[0], false).unwrap();
        assert_eq!(route_json["archive_ordinal"], 0);
        assert_eq!(route_json["sections"], binding.blocks.len());
        assert_eq!(route_json["exact_spans"], 2);
        assert_eq!(route_json["pairs"], 4);
        assert_eq!(
            route_json["descriptors"][0]["section_name"],
            binding.blocks[0].section_name
        );
        assert_eq!(route_json["reconstruction_verified"], false);
        assert!(verify_inspected_archive_routes(
            &collection.archives[0],
            binding,
            &audit.shape_routes[0].blocks,
        )
        .unwrap()
            > 0);
        shaperoute::verify_reconstruction(
            &shapes_raw,
            binding,
            &audit.shape_routes[0].blocks,
        )
        .unwrap();
        drop(chain);

        // A manifest may not silently omit a bound route section.
        let missing = test_path("routed-v4-missing-section", "aicollection");
        let manifest_only = encode_manifest(&collection).unwrap();
        write_raw_collection_fixture(&missing, VERSION, &[("manifest", &manifest_only)]);
        let missing_error = match open_collection_manifest(&missing) {
            Ok(_) => panic!("collection with an omitted route section was accepted"),
            Err(error) => error,
        };
        assert!(missing_error.to_string().contains("do not exactly match"));

        // A correctly named and checksummed, but malformed, route payload is rejected by the
        // full layer audit rather than surviving until its particular span is queried.
        let malformed = test_path("routed-v4-malformed-block", "aicollection");
        let mut bad_block = shaperoute::encode_block(&collection.shape_route_blocks[0]).unwrap();
        bad_block.push(0);
        let section_name = &binding.blocks[0].section_name;
        write_raw_collection_fixture(
            &malformed,
            VERSION,
            &[("manifest", &manifest_only), (section_name, &bad_block)],
        );
        let malformed_chain = open_collection_chain(&malformed).unwrap();
        assert!(audit_collection_layer(&malformed_chain.layers[0])
            .unwrap_err()
            .to_string()
            .contains("trailing"));
        drop(malformed_chain);

        let orphan = test_path("routed-v4-orphan-section", "aicollection");
        let good_block = shaperoute::encode_block(&collection.shape_route_blocks[0]).unwrap();
        write_raw_collection_fixture(
            &orphan,
            VERSION,
            &[
                ("manifest", &manifest_only),
                (section_name, &good_block),
                ("s.0.999", &good_block),
            ],
        );
        let orphan_error = match open_collection_manifest(&orphan) {
            Ok(_) => panic!("collection with an orphan route section was accepted"),
            Err(error) => error,
        };
        assert!(orphan_error.to_string().contains("do not exactly match"));

        for path in [source, first, second, missing, malformed, orphan] {
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn route_binding_codec_rejects_wrong_context_names_intervals_and_legacy_sources() {
        let source = scientific_archive_fixture();
        let (collection, _) = routed_collection_fixture(&source);
        let archive = &collection.archives[0];
        let binding = archive.shape_routes.as_ref().unwrap();
        let mut encoded = Vec::new();
        encode_shape_route_binding(&mut encoded, binding).unwrap();
        assert_eq!(
            decode_shape_route_binding(&mut Cursor::new(&encoded), 0, &archive.identity).unwrap(),
            *binding
        );

        let mut unknown_codec = encoded.clone();
        unknown_codec[0] = (shaperoute::ROUTE_CODEC_VERSION + 1) as u8;
        assert!(decode_shape_route_binding(
            &mut Cursor::new(&unknown_codec),
            0,
            &archive.identity,
        )
        .unwrap_err()
        .to_string()
        .contains("unsupported shape-route codec"));

        let mut wrong_root = encoded.clone();
        wrong_root[1] ^= 1;
        assert!(decode_shape_route_binding(
            &mut Cursor::new(&wrong_root),
            0,
            &archive.identity,
        )
        .unwrap_err()
        .to_string()
        .contains("source-root"));

        let mut wrong_name = binding.clone();
        wrong_name.blocks[0].section_name = "s.0.999".into();
        let mut wrong_name_raw = Vec::new();
        encode_shape_route_binding(&mut wrong_name_raw, &wrong_name).unwrap();
        assert!(decode_shape_route_binding(
            &mut Cursor::new(&wrong_name_raw),
            0,
            &archive.identity,
        )
        .unwrap_err()
        .to_string()
        .contains("disagrees with its archive/span binding"));

        let mut overlap = binding.clone();
        overlap.blocks.push(overlap.blocks[0].clone());
        assert!(shaperoute::validate_binding(
            &overlap,
            0,
            binding.source_root,
            binding.shapes_digest,
        )
        .unwrap_err()
        .to_string()
        .contains("overlaps"));

        let mut legacy_bound = collection.clone();
        legacy_bound.archives[0].identity.archive_format_version =
            evidence_io::format::SEEKABLE_VERSION;
        legacy_bound.archives[0].identity.native_scheme = LEGACY_FULL_FILE_SCHEME.into();
        let legacy_manifest = encode_collection(&legacy_bound).unwrap();
        assert!(decode_collection(&legacy_manifest, VERSION)
            .unwrap_err()
            .to_string()
            .contains("shape routes require a root-committed v2 archive"));

        let mut wrong_shapes_source = archive.clone();
        wrong_shapes_source
            .shape_routes
            .as_mut()
            .unwrap()
            .shapes_digest[0] ^= 1;
        let error = validate_archive_identity(&wrong_shapes_source, false).unwrap_err();
        assert!(format!("{error:#}").contains("compressed shapes-digest binding mismatch"));

        std::fs::remove_file(source).unwrap();
    }

    #[test]
    fn route_free_v4_is_deterministic_and_selects_the_exact_fallback() {
        let collection = fixture_collection();
        let first = test_path("route-free-v4-first", "aicollection");
        let second = test_path("route-free-v4-second", "aicollection");
        write_collection(&first, &collection).unwrap();
        write_collection(&second, &collection).unwrap();
        assert_eq!(std::fs::read(&first).unwrap(), std::fs::read(&second).unwrap());
        let chain = open_collection_chain(&first).unwrap();
        assert_eq!(chain.layers[0].file.version, VERSION);
        assert!(audit_collection_layer(&chain.layers[0])
            .unwrap()
            .shape_routes
            .is_empty());
        assert!(lookup_chain_shape_routes(&chain, 0, &[10]).unwrap().is_none());
        drop(chain);
        std::fs::remove_file(first).unwrap();
        std::fs::remove_file(second).unwrap();
    }

    #[test]
    fn routed_build_is_identical_for_repeated_and_reversed_sample_inputs() {
        let (source_a, local_a) = minimal_local_route_fixture("A", 10, 0);
        let (source_b, local_b) = minimal_local_route_fixture("B", 20, 1);
        let forward =
            assemble_collection(vec![local_a.clone(), local_b.clone()], true, 2).unwrap();
        let reversed = assemble_collection(vec![local_b, local_a], true, 2).unwrap();
        assert_eq!(forward, reversed);
        assert_eq!(
            forward
                .archives
                .iter()
                .map(|archive| archive.id.as_str())
                .collect::<Vec<_>>(),
            ["A", "B"]
        );
        assert_eq!(
            forward
                .archives
                .iter()
                .map(|archive| archive.shape_routes.as_ref().unwrap().archive_ordinal)
                .collect::<Vec<_>>(),
            [0, 1]
        );

        let first = test_path("reversed-build-first", "aicollection");
        let second = test_path("reversed-build-second", "aicollection");
        write_collection(&first, &forward).unwrap();
        write_collection(&second, &reversed).unwrap();
        assert_eq!(std::fs::read(&first).unwrap(), std::fs::read(&second).unwrap());
        for path in [source_a, source_b, first, second] {
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn routed_and_fallback_reducers_match_on_geometry_multimappers_and_cross_chunk_umis() {
        let source = scientific_archive_fixture();
        let mut reader = evidence_io::format::SectionReader::open(&source).unwrap();
        let shapes_raw = reader.read("shapes").unwrap();
        let chunks = read_chunk_index(&mut reader).unwrap();
        assert!(chunks.len() >= 3);
        let derived = shaperoute::derive_from_shapes(&shapes_raw, 0).unwrap();
        let span = |exact| {
            derived
                .blocks
                .iter()
                .find_map(|block| block.span(exact))
                .unwrap()
                .clone()
        };
        let include_route = span(10);
        let exclude_route = span(20);
        let posts = (0..chunks.len() as u32).collect::<Vec<_>>();

        let mut fallback = LazyArchive::open(&source).unwrap();
        let fallback_point = junction_counts_routed_with_shape_route(
            &mut fallback,
            &chunks,
            0,
            110,
            120,
            4,
            &posts,
            None,
        )
        .unwrap();
        let mut routed = LazyArchive::open(&source).unwrap();
        let routed_point = junction_counts_routed_with_shape_route(
            &mut routed,
            &chunks,
            0,
            110,
            120,
            4,
            &posts,
            Some((&include_route, derived.n_shapes)),
        )
        .unwrap();
        assert_eq!(fallback_point, routed_point);
        let hit_classes = routed_point
            .2
            .values()
            .flat_map(|classes| classes.iter().copied())
            .collect::<FxHashSet<_>>();
        // Class 1 has both containment representatives, class 2 is reverse-strand, class 3 is a
        // multimapper, and class 0 later reappears in another chunk.
        assert_eq!(hit_classes, FxHashSet::from_iter([0, 1, 2, 3]));

        let fallback_components = vec![
            RoutedComponent {
                chrom: 0,
                donor: 110,
                acceptor: 120,
                side: 1,
                posts: posts.clone(),
                shape_route: None,
            },
            RoutedComponent {
                chrom: 0,
                donor: 210,
                acceptor: 230,
                side: 2,
                posts: posts.clone(),
                shape_route: None,
            },
        ];
        let routed_components = vec![
            RoutedComponent {
                shape_route: Some(include_route.clone()),
                ..fallback_components[0].clone()
            },
            RoutedComponent {
                shape_route: Some(exclude_route),
                ..fallback_components[1].clone()
            },
        ];
        let mut fallback = LazyArchive::open(&source).unwrap();
        let fallback_masks =
            jset_class_masks_routed(&mut fallback, &chunks, &fallback_components, None).unwrap();
        let mut routed = LazyArchive::open(&source).unwrap();
        let routed_masks = jset_class_masks_routed(
            &mut routed,
            &chunks,
            &routed_components,
            Some(derived.n_shapes),
        )
        .unwrap();
        assert_eq!(fallback_masks, routed_masks);
        assert_eq!(routed_masks.0.get(&0), Some(&3));
        for class in [1, 2, 3] {
            assert_eq!(routed_masks.0.get(&class), Some(&1));
        }
        assert_eq!(routed_masks.0.get(&5), Some(&2));
        assert_eq!(routed_masks.1, chunks.len());
        assert_eq!(routed_masks.2, chunks.len() * 2);

        let absent = [RoutedComponent {
            chrom: 0,
            donor: 999,
            acceptor: 1009,
            side: 1,
            posts,
            shape_route: Some(include_route),
        }];
        let mut routed = LazyArchive::open(&source).unwrap();
        assert!(jset_class_masks_routed(
            &mut routed,
            &chunks,
            &absent,
            Some(derived.n_shapes),
        )
        .unwrap()
        .0
        .is_empty());

        std::fs::remove_file(source).unwrap();
    }

    #[test]
    fn junction_segment_rejects_duplicate_and_wrong_chromosome_postings() {
        let chroms = vec!["chr1".into(), "chr2".into()];
        let archives = vec![fixture_archive()];

        let valid = encode_junction_rows(&[fixture_row(vec![0])], 1);
        assert!(decode_junction_rows(&valid, &chroms, &archives).is_ok());

        let duplicate = encode_junction_rows(&[fixture_row(vec![0, 0])], 1);
        let error = decode_junction_rows(&duplicate, &chroms, &archives).unwrap_err();
        assert!(error.to_string().contains("strictly increasing"));

        let wrong_chromosome = encode_junction_rows(&[fixture_row(vec![1])], 1);
        let error = decode_junction_rows(&wrong_chromosome, &chroms, &archives).unwrap_err();
        assert!(error.to_string().contains("wrong chromosome"));
    }

    fn fixture_planning() -> CollectionQueryPlanningSummary {
        CollectionQueryPlanningSummary {
            collection_layers: 1,
            archives_total: 1,
            archives_opened: 1,
            archives_pruned: 0,
            support_bound_pruned: Some(false),
            unique_chunks_decoded: 2,
            independent_chunk_decodes: None,
            planned_compressed_bytes: 10,
            actual_archive_bytes_read: 11,
            source_archive_identity_bytes_read: 0,
            source_archive_execution_bytes_read: 11,
            source_archive_bytes_read: 11,
            collection_sidecar_bytes_read: 12,
            shape_route_sidecar_payload_bytes_read: 3,
            total_logical_bytes_read: 23,
            route_blocks_loaded: 1,
            routed_archives: 1,
            fallback_archives: 0,
        }
    }

    fn named_table<'a>(value: &'a serde_json::Value, name: &str) -> &'a serde_json::Value {
        value["data"]["tables"]
            .as_array()
            .unwrap()
            .iter()
            .find(|table| table["name"] == name)
            .unwrap()
    }

    #[test]
    fn uniform_junction_stream_has_typed_science_selection_and_set_semantics() {
        let collection = fixture_collection();
        let results = vec![SampleCounts {
            archive: 0,
            present: true,
            supporting_children: 4,
            chunks: 2,
            planned_bytes: 10,
            actual_bytes: 11,
            shape_route_used: true,
            shape_route_section: Some("s.0.100".into()),
            shape_route_sidecar_bytes: 3,
            molecules: None,
            umis: 3,
            cells: 2,
            top_cells: vec![("AAAAAAAAAAAAAAAA".into(), 2)],
        }];
        let context = uniform_source_context(collection.archives.iter(), BTreeMap::new()).unwrap();
        let summary = CollectionJunctionUniformSummary {
            coordinates: "0-based junction boundaries",
            chrom: "chr1",
            donor: 100,
            acceptor: 200,
            support_upper_bound: 4,
            min_support: 0,
            umis: 3,
            cells: 2,
            selection: CollectionSelectionPolicy {
                requested_top_per_sample: 1,
                top_zero_means_all: true,
                comparator: "UMIs descending, archive cell id ascending",
            },
            planning: fixture_planning(),
        };
        let mut bytes = Vec::new();
        stream_uniform_junction(
            &mut bytes,
            CollectionOutputFormat::Json,
            &context,
            &summary,
            &collection,
            &results,
            false,
            true,
        )
        .unwrap();
        let value: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(value["$schema"], gravlax_output::ENVELOPE_SCHEMA);
        assert_eq!(value["result_schema"], COLLECTION_JUNCTION_RESULT_SCHEMA);
        assert_eq!(value["data"]["summary"]["umis"], 3);
        assert_eq!(value["data"]["summary"]["cells"], 2);
        let samples = named_table(&value, "samples");
        assert_eq!(samples["schema"]["semantics"]["row_semantics"], "set");
        assert_eq!(samples["schema"]["semantics"]["key"], json!(["sample"]));
        assert_eq!(samples["rows"][0][3], 3);
        assert_eq!(samples["rows"][0][4], 2);
        assert_eq!(samples["rows"][0][11], "open");
        let cells = named_table(&value, "cells");
        assert_eq!(cells["selection"]["available_rows"], 2);
        assert_eq!(cells["selection"]["emitted_rows"], 1);
        assert_eq!(cells["selection"]["truncated"], true);
        assert_eq!(cells["rows"][0], json!(["A", 1, "AAAAAAAAAAAAAAAA", 2]));
    }

    #[test]
    fn collection_uniform_preflight_rejects_occupied_destination() {
        let output = test_path("uniform-occupied", "json");
        std::fs::write(&output, b"keep\n").unwrap();
        let args = CollectionOutputArgs {
            format: Some(CollectionOutputFormat::Json),
            output: Some(output.clone()),
        };
        let error = preflight_uniform_collection_output(&args).unwrap_err();
        assert!(error.to_string().contains("refusing to replace"));
        assert_eq!(std::fs::read(&output).unwrap(), b"keep\n");
        std::fs::remove_file(output).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn uniform_build_rejects_non_utf8_rows_before_writing_any_bytes() {
        use std::os::unix::ffi::OsStringExt;

        let mut collection = fixture_collection();
        collection.archives[0].path = PathBuf::from(std::ffi::OsString::from_vec(vec![0xff]));
        let summary = BuildUniformSummary {
            collection_format_version: VERSION,
            output_path: "atlas.aicollection",
            output_root: "0".repeat(64),
            shape_routes_requested: false,
            new_archives: 1,
            segment_junctions: 1,
            archive_routes: 1,
            chunk_postings: 1,
            raw_section_bytes: 1,
            file_bytes: 1,
            shape_routes: BuildShapeRouteSummary {
                archives: 0,
                sections: 0,
                exact_spans: 0,
                pairs: 0,
                raw_bytes: 0,
                compressed_bytes: 0,
            },
            source_io: BuildSourceIoSummary {
                archives: 0,
                identity_content_bytes_read: 0,
                total_bytes_read: 0,
                shape_route_payload_bytes_read: 0,
                shape_route_source_bytes_read: 0,
            },
        };
        let mut bytes = Vec::new();
        let error = stream_uniform_build(
            &mut bytes,
            CollectionOutputFormat::Json,
            &ResultContext::default(),
            &summary,
            &collection,
            &[],
        )
        .unwrap_err();
        assert!(error.to_string().contains("not valid UTF-8"));
        assert!(bytes.is_empty());
    }
}
