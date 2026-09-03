//! The `.aie` container: named, individually-zstd'd sections behind a tiny directory.
//!
//! v0 deliberately stops at whole-archive columnar sections (genome-ordered rows, one stream per
//! column). v1 adds a seekable tail directory. v2 commits that directory, including one digest per
//! compressed section payload, so a reader can authenticate the exact encoded sections it selects
//! without scanning unrelated molecule evidence. The commitment is not a publisher signature.

use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::io::{Read, Write};
use std::path::Path;

pub const MAGIC: &[u8; 4] = b"AIE0";
/// v1 adds a seekable directory at the end of the file. v2 extends each directory entry with a
/// digest of the exact compressed payload and carries a root commitment in its footer. Readers
/// retain the v1 path verbatim; writers emit only v2.
pub const VERSION: u32 = 2;
pub const SEEKABLE_VERSION: u32 = 1;
/// Domain separator for the rooted archive commitment. Changing these bytes changes every v2
/// archive identity.
pub const ROOT_DOMAIN: &[u8] = b"gravlax-aie-directory-root-v2\0";
/// Scheme-independent exact encoded-section identity used to reject the same evidence when one
/// copy is legacy v1 and another is its byte-preserving v2 seal.
pub const ENCODED_CONTENT_DOMAIN: &[u8] = b"gravlax-aie-encoded-sections-v1\0";
const V1_FOOTER_LEN: u64 = 12;
const V2_FOOTER_LEN: u64 = 44;
const MAX_V2_SECTIONS: usize = 1_000_000;
/// The largest raw section in the benchmark archive corpus was 100,768,566 bytes. This limit
/// leaves more than fivefold headroom while preventing a tiny hostile frame from expanding toward
/// host-scale memory.
pub const MAX_V2_RAW_SECTION_BYTES: u64 = 512 * 1024 * 1024;
/// Includes ample zstd framing/compress-bound headroom for an incompressible maximum-size raw
/// section while retaining a fixed allocation ceiling.
pub const MAX_V2_COMPRESSED_SECTION_BYTES: u64 = MAX_V2_RAW_SECTION_BYTES + 16 * 1024 * 1024;
const MAX_V2_COMPRESSION_RATIO: u64 = 4_096;
const V2_COMPRESSION_SLACK: u64 = 64 * 1024 * 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ArchiveCommitment {
    pub version: u32,
    pub digest: [u8; 32],
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LegacyIdentityScan {
    pub full_file_blake3: [u8; 32],
    pub encoded_sections_blake3: [u8; 32],
    pub bytes_read: u64,
}

/// Authenticated directory metadata for one section.
///
/// `compressed_blake3` is present only in a root-committed v2 archive.  Its absence for v1 is
/// deliberate: no digest read from an unauthenticated legacy directory could provide the same
/// guarantee.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SectionMetadata<'a> {
    pub name: &'a str,
    pub offset: u64,
    pub raw_len: u64,
    pub compressed_len: u64,
    pub compressed_blake3: Option<[u8; 32]>,
}

impl ArchiveCommitment {
    pub fn to_hex(self) -> String {
        self.digest.iter().map(|byte| format!("{byte:02x}")).collect()
    }
}

fn root_hasher(directory_offset: u64) -> blake3::Hasher {
    let mut hasher = blake3::Hasher::new();
    hasher.update(ROOT_DOMAIN);
    hasher.update(MAGIC);
    hasher.update(&VERSION.to_le_bytes());
    hasher.update(&directory_offset.to_le_bytes());
    hasher
}

fn encoded_content_identity(
    directory: &[(String, u64, u64, u64)],
    payload_digests: &[[u8; 32]],
) -> Result<[u8; 32]> {
    if directory.len() != payload_digests.len() {
        bail!("section directory and payload-digest cardinalities differ");
    }
    let count = u32::try_from(directory.len()).context("too many sections for content identity")?;
    let mut hasher = blake3::Hasher::new();
    hasher.update(ENCODED_CONTENT_DOMAIN);
    hasher.update(&count.to_le_bytes());
    for ((name, _, raw_len, comp_len), digest) in directory.iter().zip(payload_digests) {
        let name_len = u8::try_from(name.len()).context("section name exceeds u8")?;
        hasher.update(&[name_len]);
        hasher.update(name.as_bytes());
        hasher.update(&raw_len.to_le_bytes());
        hasher.update(&comp_len.to_le_bytes());
        hasher.update(digest);
    }
    Ok(*hasher.finalize().as_bytes())
}

fn validate_section_lengths(version: u32, name: &str, raw_len: u64, comp_len: u64) -> Result<()> {
    if raw_len > MAX_V2_RAW_SECTION_BYTES {
        bail!("archive v{version} section {name} raw length {raw_len} exceeds safety limit");
    }
    if comp_len > MAX_V2_COMPRESSED_SECTION_BYTES {
        bail!(
            "archive v{version} section {name} compressed length {comp_len} exceeds safety limit"
        );
    }
    let permitted_raw = comp_len
        .saturating_mul(MAX_V2_COMPRESSION_RATIO)
        .saturating_add(V2_COMPRESSION_SLACK);
    if raw_len > permitted_raw {
        bail!("archive v{version} section {name} declares an unsafe compression ratio");
    }
    Ok(())
}

pub struct SectionWriter {
    out: std::io::BufWriter<std::fs::File>,
    directory: Vec<(String, u64, u64, [u8; 32])>,
    offsets: Vec<u64>,
    level: i32,
}

impl SectionWriter {
    pub fn create(path: &Path, level: i32) -> Result<Self> {
        let file = std::fs::File::create(path)
            .with_context(|| format!("creating {}", path.display()))?;
        Self::from_file(file, level)
    }

    /// Create a writer without permitting an existing destination to be truncated.
    pub fn create_new(path: &Path, level: i32) -> Result<Self> {
        let file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(path)
            .with_context(|| format!("creating new archive {}", path.display()))?;
        Self::from_file(file, level)
    }

    fn from_file(file: std::fs::File, level: i32) -> Result<Self> {
        let mut out = std::io::BufWriter::new(file);
        out.write_all(MAGIC)?;
        out.write_all(&VERSION.to_le_bytes())?;
        Ok(SectionWriter { out, directory: Vec::new(), offsets: Vec::new(), level })
    }

    /// Clone the open file description used by this writer. This is intended for ownership
    /// guards around staging cleanup; section bytes remain buffered until `finish`.
    pub fn try_clone_file(&self) -> Result<std::fs::File> {
        Ok(self.out.get_ref().try_clone()?)
    }

    pub fn section(&mut self, name: &str, raw: &[u8]) -> Result<()> {
        if raw.len() as u64 > MAX_V2_RAW_SECTION_BYTES {
            bail!("section {name} raw length {} exceeds v2 safety limit", raw.len());
        }
        let comp = compress(raw, self.level)?;
        self.section_precompressed(name, raw.len() as u64, &comp)
    }

    /// Write a section whose zstd payload was produced elsewhere (e.g. compressed in parallel
    /// with rayon). `comp` must be exactly what `compress(raw, level)` yields.
    pub fn section_precompressed(&mut self, name: &str, raw_len: u64, comp: &[u8]) -> Result<()> {
        use std::io::Seek;
        if name.is_empty() {
            bail!("section name must not be empty");
        }
        if name.len() > u8::MAX as usize {
            bail!("section name is {} bytes; maximum is {}", name.len(), u8::MAX);
        }
        if self.directory.iter().any(|(existing, _, _, _)| existing == name) {
            bail!("duplicate section name {name}");
        }
        if self.directory.len() >= MAX_V2_SECTIONS {
            bail!("archive section count exceeds v2 safety limit {MAX_V2_SECTIONS}");
        }
        let comp_len = comp.len() as u64;
        validate_section_lengths(VERSION, name, raw_len, comp_len)?;
        self.offsets.push(self.out.stream_position()?);
        self.out.write_all(&(name.len() as u8).to_le_bytes())?;
        self.out.write_all(name.as_bytes())?;
        self.out.write_all(&raw_len.to_le_bytes())?;
        self.out.write_all(&comp_len.to_le_bytes())?;
        self.out.write_all(comp)?;
        self.directory.push((
            name.to_string(),
            raw_len,
            comp_len,
            *blake3::hash(comp).as_bytes(),
        ));
        Ok(())
    }

    /// Finish a canonical v2 file: terminator, authenticated physical-order directory, and fixed
    /// 44-byte footer `(directory_offset, root, "AIED")`. Returns the historical three-column
    /// section accounting used by ingest reports.
    pub fn finish(self) -> Result<Vec<(String, u64, u64)>> {
        let (accounting, _file, _commitment) = self.finish_with_file()?;
        Ok(accounting)
    }

    /// Finish while retaining the exact open file description and its computed root commitment.
    /// This supports publication and provenance without reopening a replaceable staging path.
    pub fn finish_with_file(
        mut self,
    ) -> Result<(Vec<(String, u64, u64)>, std::fs::File, ArchiveCommitment)> {
        use std::io::Seek;
        self.out.write_all(&0u8.to_le_bytes())?; // terminator for linear readers
        let dir_offset = self.out.stream_position()?;
        let section_count = u32::try_from(self.offsets.len()).context("too many archive sections")?;
        let mut directory = Vec::new();
        directory.extend_from_slice(&section_count.to_le_bytes());
        for ((name, raw, comp, digest), off) in self.directory.iter().zip(&self.offsets) {
            directory.push(name.len() as u8);
            directory.extend_from_slice(name.as_bytes());
            directory.extend_from_slice(&off.to_le_bytes());
            directory.extend_from_slice(&raw.to_le_bytes());
            directory.extend_from_slice(&comp.to_le_bytes());
            directory.extend_from_slice(digest);
        }
        let mut hasher = root_hasher(dir_offset);
        hasher.update(&directory);
        let root = hasher.finalize();
        self.out.write_all(&directory)?;
        self.out.write_all(&dir_offset.to_le_bytes())?;
        self.out.write_all(root.as_bytes())?;
        self.out.write_all(b"AIED")?;
        self.out.flush()?;
        let accounting = self
            .directory
            .into_iter()
            .map(|(name, raw, comp, _)| (name, raw, comp))
            .collect();
        let file = self.out.into_inner().map_err(|error| error.into_error())?;
        Ok((
            accounting,
            file,
            ArchiveCommitment {
                version: VERSION,
                digest: *root.as_bytes(),
            },
        ))
    }
}

/// One zstd frame, as sections store them.
pub fn compress(raw: &[u8], level: i32) -> Result<Vec<u8>> {
    let mut enc = zstd::Encoder::new(Vec::new(), level)?;
    enc.write_all(raw)?;
    Ok(enc.finish()?)
}

/// Inverse of `compress`, verifying the expected length.
pub fn decompress(comp: &[u8], raw_len: usize) -> Result<Vec<u8>> {
    let dec = zstd::Decoder::new(comp)?;
    let limit = raw_len.checked_add(1).context("declared raw section length overflow")?;
    let mut raw = Vec::with_capacity(raw_len.min(8 * 1024 * 1024));
    dec.take(limit as u64).read_to_end(&mut raw)?;
    if raw.len() != raw_len {
        bail!("decompressed {} bytes, expected {raw_len}", raw.len());
    }
    Ok(raw)
}

/// Read every section into memory. Archives are tens of MB by design; streaming readers can come
/// with the chunked layout.
pub fn read_sections(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let mut f = std::io::BufReader::new(
        std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?,
    );
    let mut magic = [0u8; 4];
    f.read_exact(&mut magic)?;
    if &magic != MAGIC {
        bail!("not an .aie file (bad magic)");
    }
    let mut v = [0u8; 4];
    f.read_exact(&mut v)?;
    let version = u32::from_le_bytes(v);
    if version > VERSION {
        bail!("archive version {version} is newer than this reader ({VERSION})");
    }
    if version >= SEEKABLE_VERSION {
        // No seekable path may bypass directory validation and the shared allocation bounds.  In
        // addition, v2 reads authenticate the directory and each selected compressed payload.
        let file = f.into_inner();
        let mut reader = SectionReader::from_file(file)?;
        let names: Vec<String> = reader.names().map(str::to_owned).collect();
        return names
            .into_iter()
            .map(|name| {
                let raw = reader.read(&name)?;
                Ok((name, raw))
            })
            .collect();
    }
    let mut out = Vec::new();
    loop {
        let mut nl = [0u8; 1];
        if f.read_exact(&mut nl).is_err() {
            break;
        }
        if nl[0] == 0 {
            break;
        }
        let mut name = vec![0u8; nl[0] as usize];
        f.read_exact(&mut name)?;
        let mut sz = [0u8; 8];
        f.read_exact(&mut sz)?;
        let raw_len = u64::from_le_bytes(sz) as usize;
        f.read_exact(&mut sz)?;
        let comp_len = u64::from_le_bytes(sz) as usize;
        let mut comp = vec![0u8; comp_len];
        f.read_exact(&mut comp)?;
        let raw = decompress(&comp, raw_len).with_context(|| {
            format!("section {:?}", String::from_utf8_lossy(&name))
        })?;
        out.push((String::from_utf8_lossy(&name).into_owned(), raw));
    }
    Ok(out)
}

/// Seekable section reader: loads only the directory, then serves sections on demand.
pub struct SectionReader {
    file: std::fs::File,
    directory: Vec<(String, u64, u64, u64)>, // name, offset, raw_len, comp_len
    payload_digests: Vec<Option<[u8; 32]>>,
    version: u32,
    root: Option<[u8; 32]>,
    bytes_read: std::sync::atomic::AtomicU64,
}

impl SectionReader {
    pub fn open(path: &Path) -> Result<SectionReader> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("opening {}", path.display()))?;
        Self::from_file(file)
    }

    /// Parse one already-open file description. Callers that first inspect a version byte can use
    /// this entry point without reopening a replaceable path.
    pub fn from_file(mut file: std::fs::File) -> Result<SectionReader> {
        use std::io::{Seek, SeekFrom};
        file.seek(SeekFrom::Start(0))?;
        let file_len = file.metadata()?.len();
        if file_len < 20 {
            bail!("archive is truncated: {file_len} bytes");
        }
        let mut head = [0u8; 8];
        file.read_exact(&mut head)?;
        if &head[..4] != MAGIC {
            bail!("not an .aie file");
        }
        let version = u32::from_le_bytes(head[4..8].try_into().unwrap());
        if version < 1 {
            bail!("archive version {version} predates the seekable directory; use read_sections");
        }
        if version > VERSION {
            bail!("archive version {version} is newer than this reader ({VERSION})");
        }
        match version {
            SEEKABLE_VERSION => Self::open_v1(file, file_len),
            VERSION => Self::open_v2(file, file_len, head),
            _ => unreachable!(),
        }
    }

    /// The v1 parser remains behavior-compatible with the original archive reader.
    fn open_v1(mut file: std::fs::File, file_len: u64) -> Result<Self> {
        use std::io::{Seek, SeekFrom};
        file.seek(SeekFrom::End(-(V1_FOOTER_LEN as i64)))?;
        let mut tail = [0u8; 12];
        file.read_exact(&mut tail)?;
        if &tail[8..] != b"AIED" {
            bail!("missing directory footer");
        }
        let dir_offset = u64::from_le_bytes(tail[..8].try_into().unwrap());
        let footer_offset = file_len - 12;
        if dir_offset < 9 || dir_offset > footer_offset.saturating_sub(4) {
            bail!("invalid directory offset {dir_offset} for {file_len}-byte archive");
        }
        file.seek(SeekFrom::Start(dir_offset))?;
        let mut n4 = [0u8; 4];
        file.read_exact(&mut n4)?;
        let n = u32::from_le_bytes(n4) as usize;
        if n > MAX_V2_SECTIONS {
            bail!("v1 directory declares {n} sections; limit is {MAX_V2_SECTIONS}");
        }
        let directory_bytes = footer_offset - dir_offset - 4;
        if n as u64 > directory_bytes / 25 {
            bail!("directory declares {n} entries in only {directory_bytes} bytes");
        }
        let mut directory = Vec::with_capacity(n);
        let mut names = HashSet::with_capacity(n);
        for _ in 0..n {
            let mut nl = [0u8; 1];
            file.read_exact(&mut nl)?;
            let mut name = vec![0u8; nl[0] as usize];
            file.read_exact(&mut name)?;
            let mut u = [0u8; 8];
            file.read_exact(&mut u)?;
            let offset = u64::from_le_bytes(u);
            file.read_exact(&mut u)?;
            let raw_len = u64::from_le_bytes(u);
            file.read_exact(&mut u)?;
            let comp_len = u64::from_le_bytes(u);
            let name = String::from_utf8(name).context("section name is not valid UTF-8")?;
            if name.is_empty() {
                bail!("directory contains an empty section name");
            }
            if !names.insert(name.clone()) {
                bail!("directory contains duplicate section {name}");
            }
            validate_section_lengths(SEEKABLE_VERSION, &name, raw_len, comp_len)?;
            let header_len = 1u64
                .checked_add(name.len() as u64)
                .and_then(|v| v.checked_add(16))
                .context("section header length overflow")?;
            let section_end = offset
                .checked_add(header_len)
                .and_then(|v| v.checked_add(comp_len))
                .context("section extent overflow")?;
            if offset < 8 || section_end > dir_offset - 1 {
                bail!(
                    "section {name} extent {offset}..{section_end} lies outside the section area"
                );
            }
            directory.push((name, offset, raw_len, comp_len));
        }
        if file.stream_position()? != footer_offset {
            bail!("directory length does not match footer offset");
        }
        // The directory is not trusted until it agrees with each inline section header.
        let mut bytes_read = 8u64
            .checked_add(12)
            .and_then(|value| value.checked_add(footer_offset - dir_offset))
            .context("archive open byte accounting overflow")?;
        for (name, offset, raw_len, comp_len) in &directory {
            file.seek(SeekFrom::Start(*offset))?;
            let mut nl = [0u8; 1];
            file.read_exact(&mut nl)?;
            if nl[0] as usize != name.len() {
                bail!("section {name} has inconsistent inline name length");
            }
            let mut inline_name = vec![0u8; nl[0] as usize];
            file.read_exact(&mut inline_name)?;
            if inline_name != name.as_bytes() {
                bail!("section {name} has inconsistent inline name");
            }
            let mut u = [0u8; 8];
            file.read_exact(&mut u)?;
            if u64::from_le_bytes(u) != *raw_len {
                bail!("section {name} has inconsistent raw length");
            }
            file.read_exact(&mut u)?;
            if u64::from_le_bytes(u) != *comp_len {
                bail!("section {name} has inconsistent compressed length");
            }
            bytes_read = bytes_read
                .checked_add(1 + name.len() as u64 + 16)
                .context("archive open byte accounting overflow")?;
        }
        Ok(SectionReader {
            file,
            payload_digests: vec![None; directory.len()],
            directory,
            version: SEEKABLE_VERSION,
            root: None,
            bytes_read: std::sync::atomic::AtomicU64::new(bytes_read),
        })
    }

    fn open_v2(mut file: std::fs::File, file_len: u64, head: [u8; 8]) -> Result<Self> {
        use std::io::{Seek, SeekFrom};
        let minimum = 8 + 1 + 4 + V2_FOOTER_LEN;
        if file_len < minimum {
            bail!("v2 archive is truncated: {file_len} bytes");
        }
        file.seek(SeekFrom::End(-(V2_FOOTER_LEN as i64)))?;
        let mut footer = [0u8; V2_FOOTER_LEN as usize];
        file.read_exact(&mut footer)?;
        if &footer[40..] != b"AIED" {
            bail!("missing v2 directory footer");
        }
        let dir_offset = u64::from_le_bytes(footer[..8].try_into().unwrap());
        let footer_offset = file_len - V2_FOOTER_LEN;
        if dir_offset < 9 || dir_offset > footer_offset.saturating_sub(4) {
            bail!("invalid v2 directory offset {dir_offset} for {file_len}-byte archive");
        }
        let mut expected_root = [0u8; 32];
        expected_root.copy_from_slice(&footer[8..40]);
        file.seek(SeekFrom::Start(dir_offset))?;
        let mut hasher = root_hasher(dir_offset);
        let mut n4 = [0u8; 4];
        file.read_exact(&mut n4)?;
        hasher.update(&n4);
        let n = u32::from_le_bytes(n4) as usize;
        if n > MAX_V2_SECTIONS {
            bail!("v2 directory declares {n} sections; limit is {MAX_V2_SECTIONS}");
        }
        let directory_bytes = footer_offset - dir_offset - 4;
        // One non-empty one-byte name plus offset/raw/compressed/digest fields.
        if n as u64 > directory_bytes / 58 {
            bail!("v2 directory declares {n} entries in only {directory_bytes} bytes");
        }
        let mut directory = Vec::with_capacity(n);
        let mut payload_digests = Vec::with_capacity(n);
        let mut names = HashSet::with_capacity(n);
        let mut expected_offset = 8u64;
        for _ in 0..n {
            let mut nl = [0u8; 1];
            file.read_exact(&mut nl)?;
            hasher.update(&nl);
            if nl[0] == 0 {
                bail!("v2 directory contains an empty section name");
            }
            let mut name_bytes = vec![0u8; nl[0] as usize];
            file.read_exact(&mut name_bytes)?;
            hasher.update(&name_bytes);
            let name = String::from_utf8(name_bytes)
                .context("v2 section name is not valid UTF-8")?;
            if !names.insert(name.clone()) {
                bail!("v2 directory contains duplicate section {name}");
            }
            let mut fields = [0u8; 56];
            file.read_exact(&mut fields)?;
            hasher.update(&fields);
            let offset = u64::from_le_bytes(fields[..8].try_into().unwrap());
            let raw_len = u64::from_le_bytes(fields[8..16].try_into().unwrap());
            let comp_len = u64::from_le_bytes(fields[16..24].try_into().unwrap());
            let mut digest = [0u8; 32];
            digest.copy_from_slice(&fields[24..]);
            validate_section_lengths(VERSION, &name, raw_len, comp_len)?;
            let header_len = 1u64
                .checked_add(name.len() as u64)
                .and_then(|value| value.checked_add(16))
                .context("v2 section header length overflow")?;
            let section_end = offset
                .checked_add(header_len)
                .and_then(|value| value.checked_add(comp_len))
                .context("v2 section extent overflow")?;
            if offset != expected_offset {
                bail!(
                    "v2 section {name} begins at {offset}, expected contiguous physical offset {expected_offset}"
                );
            }
            if section_end > dir_offset.saturating_sub(1) {
                bail!("v2 section {name} extends outside the canonical section area");
            }
            expected_offset = section_end;
            directory.push((name, offset, raw_len, comp_len));
            payload_digests.push(Some(digest));
        }
        if file.stream_position()? != footer_offset {
            bail!("v2 directory length does not match footer offset");
        }
        if expected_offset.checked_add(1) != Some(dir_offset) {
            bail!("v2 section area is not contiguous with its directory");
        }
        if hasher.finalize().as_bytes() != &expected_root {
            bail!("v2 archive directory root mismatch");
        }
        file.seek(SeekFrom::Start(dir_offset - 1))?;
        let mut terminator = [0u8; 1];
        file.read_exact(&mut terminator)?;
        if terminator[0] != 0 {
            bail!("v2 archive is missing its canonical section terminator");
        }
        // Root authentication precedes interpretation of inline section headers. Exact agreement
        // then proves that the committed physical route points at the intended payload bytes.
        let mut bytes_read = 8u64
            .checked_add(V2_FOOTER_LEN)
            .and_then(|value| value.checked_add(footer_offset - dir_offset))
            .and_then(|value| value.checked_add(1))
            .context("v2 archive open byte accounting overflow")?;
        for (name, offset, raw_len, comp_len) in &directory {
            file.seek(SeekFrom::Start(*offset))?;
            let mut nl = [0u8; 1];
            file.read_exact(&mut nl)?;
            if nl[0] as usize != name.len() {
                bail!("v2 section {name} has inconsistent inline name length");
            }
            let mut inline_name = vec![0u8; nl[0] as usize];
            file.read_exact(&mut inline_name)?;
            if inline_name != name.as_bytes() {
                bail!("v2 section {name} has inconsistent inline name");
            }
            let mut u = [0u8; 8];
            file.read_exact(&mut u)?;
            if u64::from_le_bytes(u) != *raw_len {
                bail!("v2 section {name} has inconsistent raw length");
            }
            file.read_exact(&mut u)?;
            if u64::from_le_bytes(u) != *comp_len {
                bail!("v2 section {name} has inconsistent compressed length");
            }
            bytes_read = bytes_read
                .checked_add(1 + name.len() as u64 + 16)
                .context("v2 archive open byte accounting overflow")?;
        }
        Ok(Self {
            file,
            directory,
            payload_digests,
            version: u32::from_le_bytes(head[4..8].try_into().unwrap()),
            root: Some(expected_root),
            bytes_read: std::sync::atomic::AtomicU64::new(bytes_read),
        })
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.directory.iter().map(|(n, _, _, _)| n.as_str())
    }

    /// (name, offset, raw_len, comp_len) for every section — the per-section byte accounting.
    pub fn entries(&self) -> &[(String, u64, u64, u64)] {
        &self.directory
    }

    /// Authenticated per-section directory records in physical order.  For v1, the structural
    /// fields are validated against inline headers but `compressed_blake3` is necessarily absent.
    pub fn section_metadata(&self) -> impl ExactSizeIterator<Item = SectionMetadata<'_>> + '_ {
        self.directory.iter().enumerate().map(|(index, (name, offset, raw, compressed))| {
            SectionMetadata {
                name,
                offset: *offset,
                raw_len: *raw,
                compressed_len: *compressed,
                compressed_blake3: self.payload_digests[index],
            }
        })
    }

    pub fn archive_version(&self) -> u32 {
        self.version
    }

    /// Exact encoded-content commitment for a rooted v2 archive. Legacy v1 archives return None
    /// and require an external full-file digest when a persistent identity is needed.
    pub fn content_commitment(&self) -> Option<ArchiveCommitment> {
        self.root.map(|digest| ArchiveCommitment { version: self.version, digest })
    }

    /// Identity of ordered exact encoded section payloads, independent of the v1/v2 container
    /// directory/footer. Rooted v2 derives it without reading payload bytes.
    pub fn encoded_content_identity(&self) -> Result<Option<[u8; 32]>> {
        let Some(digests) = self
            .payload_digests
            .iter()
            .copied()
            .collect::<Option<Vec<_>>>()
        else {
            return Ok(None);
        };
        Ok(Some(encoded_content_identity(&self.directory, &digests)?))
    }

    /// One sequential legacy-v1 pass computes both the historical full-file identity and the
    /// scheme-independent exact encoded-section identity. This is used only where v1 lacks a
    /// committed directory; v2 callers must use [`Self::encoded_content_identity`].
    pub fn scan_legacy_identities(&self) -> Result<LegacyIdentityScan> {
        if self.version != SEEKABLE_VERSION {
            bail!("legacy identity scan is only valid for archive v1");
        }
        let file_len = self.file.metadata()?.len();
        // v1 readers historically allow directory order to differ from physical order. Scan in
        // physical order but deposit each digest at its original directory ordinal so sealing a
        // legacy archive retains its exact ordered encoded-section identity.
        let mut spans = Vec::with_capacity(self.directory.len());
        for (ordinal, (name, offset, _, comp_len)) in self.directory.iter().enumerate() {
            let start = offset
                .checked_add(1 + name.len() as u64 + 16)
                .context("legacy payload offset overflow")?;
            let end = start.checked_add(*comp_len).context("legacy payload extent overflow")?;
            spans.push((start, end, ordinal, name.as_str()));
        }
        spans.sort_unstable_by_key(|&(start, end, ordinal, _)| (start, end, ordinal));
        for pair in spans.windows(2) {
            if pair[1].0 < pair[0].1 {
                bail!(
                    "legacy sections {} and {} have overlapping payload extents",
                    pair[0].3,
                    pair[1].3
                );
            }
        }
        let mut full_hasher = blake3::Hasher::new();
        let mut section_hasher = blake3::Hasher::new();
        let mut span_index = 0usize;
        let mut payload_digests = vec![[0u8; 32]; self.directory.len()];
        let mut absolute = 0u64;
        let mut buffer = vec![0u8; 4 * 1024 * 1024];
        while absolute < file_len {
            let n = usize::try_from((file_len - absolute).min(buffer.len() as u64)).unwrap();
            self.read_exact_at(&mut buffer[..n], absolute)?;
            full_hasher.update(&buffer[..n]);
            let chunk_start = absolute;
            let chunk_end = absolute
                .checked_add(n as u64)
                .context("legacy scan position overflow")?;
            while span_index < spans.len() {
                let (payload_start, payload_end, ordinal, _) = spans[span_index];
                if chunk_start >= payload_end {
                    payload_digests[ordinal] = *section_hasher.finalize().as_bytes();
                    section_hasher = blake3::Hasher::new();
                    span_index += 1;
                    continue;
                }
                if chunk_end <= payload_start {
                    break;
                }
                let overlap_start = chunk_start.max(payload_start);
                let overlap_end = chunk_end.min(payload_end);
                if overlap_start < overlap_end {
                    let begin = usize::try_from(overlap_start - chunk_start).unwrap();
                    let end = usize::try_from(overlap_end - chunk_start).unwrap();
                    section_hasher.update(&buffer[begin..end]);
                }
                if chunk_end >= payload_end {
                    payload_digests[ordinal] = *section_hasher.finalize().as_bytes();
                    section_hasher = blake3::Hasher::new();
                    span_index += 1;
                    continue;
                }
                break;
            }
            absolute = chunk_end;
        }
        while span_index < spans.len() {
            let (payload_start, payload_end, ordinal, name) = spans[span_index];
            if payload_end > file_len || payload_start != payload_end {
                bail!("legacy scan did not reach section {name} payload");
            }
            payload_digests[ordinal] = *section_hasher.finalize().as_bytes();
            section_hasher = blake3::Hasher::new();
            span_index += 1;
        }
        self.record_read(file_len);
        Ok(LegacyIdentityScan {
            full_file_blake3: *full_hasher.finalize().as_bytes(),
            encoded_sections_blake3: encoded_content_identity(&self.directory, &payload_digests)?,
            bytes_read: file_len,
        })
    }

    /// Logical file bytes fetched by this reader, including directory/header validation and
    /// compressed section payloads. This is deterministic cache-independent I/O accounting.
    pub fn bytes_read(&self) -> u64 {
        self.bytes_read.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Metadata for the same open file description used by all subsequent section reads.
    pub fn file_metadata(&self) -> Result<std::fs::Metadata> {
        Ok(self.file.metadata()?)
    }

    /// Clone the already-open archive file description without resolving its pathname again.
    /// Callers can use this to copy or hash the exact inode whose directory was validated.
    pub fn try_clone_file(&self) -> Result<std::fs::File> {
        Ok(self.file.try_clone()?)
    }

    fn record_read(&self, bytes: u64) {
        let _ = self.bytes_read.fetch_update(
            std::sync::atomic::Ordering::Relaxed,
            std::sync::atomic::Ordering::Relaxed,
            |current| Some(current.saturating_add(bytes)),
        );
    }

    fn read_exact_at(&self, buffer: &mut [u8], offset: u64) -> Result<()> {
        #[cfg(unix)]
        {
            use std::os::unix::fs::FileExt;
            self.file.read_exact_at(buffer, offset)?;
            Ok(())
        }
        #[cfg(not(unix))]
        {
            use std::io::{Seek, SeekFrom};
            let mut file = self.file.try_clone()?;
            file.seek(SeekFrom::Start(offset))?;
            file.read_exact(buffer)?;
            Ok(())
        }
    }

    fn entry_index(&self, name: &str) -> Result<usize> {
        self.directory
            .iter()
            .position(|(candidate, _, _, _)| candidate == name)
            .with_context(|| format!("archive missing section {name}"))
    }

    fn verify_compressed(&self, index: usize, name: &str, compressed: &[u8]) -> Result<()> {
        if let Some(expected) = self.payload_digests.get(index).copied().flatten() {
            if blake3::hash(compressed).as_bytes() != &expected {
                bail!("archive section {name} compressed-payload digest mismatch");
            }
        }
        Ok(())
    }

    /// Decompress one section. Seeks past the section's inline header to its zstd payload.
    pub fn read(&mut self, name: &str) -> Result<Vec<u8>> {
        use std::io::{Seek, SeekFrom};
        let index = self.entry_index(name)?;
        let (_, offset, raw_len, comp_len) = self.directory[index].clone();
        // Inline layout at `offset`: name_len u8, name, raw u64, comp u64, payload.
        let skip = 1 + name.len() as u64 + 16;
        self.file.seek(SeekFrom::Start(offset + skip))?;
        let comp_len = usize::try_from(comp_len).context("compressed section is too large")?;
        let raw_len_usize = usize::try_from(raw_len).context("raw section is too large")?;
        let mut comp = vec![0u8; comp_len];
        self.file.read_exact(&mut comp)?;
        self.record_read(comp_len as u64);
        self.verify_compressed(index, name, &comp)?;
        let raw = decompress(&comp, raw_len_usize)
            .with_context(|| format!("section {name}"))?;
        Ok(raw)
    }

    pub fn has(&self, name: &str) -> bool {
        self.directory.iter().any(|(n, _, _, _)| n == name)
    }

    /// Like [`Self::read_compressed`], but position-independent (`pread`) and `&self`, so many
    /// sections can be read from worker threads concurrently.
    #[cfg(unix)]
    pub fn read_compressed_at(&self, name: &str) -> Result<(Vec<u8>, usize)> {
        use std::os::unix::fs::FileExt;
        let index = self.entry_index(name)?;
        let (_, offset, raw_len, comp_len) = self.directory[index].clone();
        let skip = 1 + name.len() as u64 + 16;
        let comp_len = usize::try_from(comp_len).context("compressed section is too large")?;
        let raw_len = usize::try_from(raw_len).context("raw section is too large")?;
        let mut comp = vec![0u8; comp_len];
        self.file.read_exact_at(&mut comp, offset + skip)?;
        self.record_read(comp_len as u64);
        self.verify_compressed(index, name, &comp)?;
        Ok((comp, raw_len))
    }

    /// Read one section's zstd payload without decompressing — pair with `decompress` to spread
    /// decode work across threads. Returns (compressed bytes, raw length).
    pub fn read_compressed(&mut self, name: &str) -> Result<(Vec<u8>, usize)> {
        use std::io::{Seek, SeekFrom};
        let index = self.entry_index(name)?;
        let (_, offset, raw_len, comp_len) = self.directory[index].clone();
        let skip = 1 + name.len() as u64 + 16;
        self.file.seek(SeekFrom::Start(offset + skip))?;
        let comp_len = usize::try_from(comp_len).context("compressed section is too large")?;
        let raw_len = usize::try_from(raw_len).context("raw section is too large")?;
        let mut comp = vec![0u8; comp_len];
        self.file.read_exact(&mut comp)?;
        self.record_read(comp_len as u64);
        self.verify_compressed(index, name, &comp)?;
        Ok((comp, raw_len))
    }

    /// Verify every compressed payload against the rooted v2 directory and decode it to validate
    /// its declared raw length. This is the explicit whole-content audit; ordinary reads verify
    /// and decode only selected sections.
    /// Returns the number of compressed payload bytes read.
    pub fn verify_all_payloads(&self) -> Result<u64> {
        if self.root.is_none() {
            bail!("archive v1 has no internal payload commitment");
        }
        let mut total = 0u64;
        for (index, (name, offset, raw_len, comp_len)) in self.directory.iter().enumerate() {
            let skip = 1 + name.len() as u64 + 16;
            let position = offset
                .checked_add(skip)
                .context("verified payload offset overflow")?;
            let compressed_len = usize::try_from(*comp_len)
                .context("verified compressed length exceeds usize")?;
            let mut compressed = vec![0u8; compressed_len];
            self.read_exact_at(&mut compressed, position)?;
            self.record_read(*comp_len);
            total = total
                .checked_add(*comp_len)
                .context("verified payload byte count overflow")?;
            let expected = self.payload_digests[index]
                .context("rooted archive section lacks a payload digest")?;
            if blake3::hash(&compressed).as_bytes() != &expected {
                bail!("archive section {name} compressed-payload digest mismatch");
            }
            let raw_len = usize::try_from(*raw_len).context("verified raw length exceeds usize")?;
            decompress(&compressed, raw_len)
                .with_context(|| format!("verifying archive section {name}"))?;
        }
        Ok(total)
    }
}

/// Varint cursor for decoding streams.
pub struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    pub fn new(buf: &'a [u8]) -> Self {
        Cursor { buf, pos: 0 }
    }
    pub fn varint(&mut self) -> Result<u64> {
        let (mut v, mut shift) = (0u64, 0);
        loop {
            let Some(&b) = self.buf.get(self.pos) else { bail!("varint past end") };
            self.pos += 1;
            if shift == 63 && b & 0x7e != 0 {
                bail!("varint overflows u64");
            }
            v |= ((b & 0x7f) as u64) << shift;
            if b & 0x80 == 0 {
                return Ok(v);
            }
            if shift >= 63 {
                bail!("varint exceeds 10 bytes");
            }
            shift += 7;
        }
    }
    pub fn svarint(&mut self) -> Result<i64> {
        let v = self.varint()?;
        Ok(((v >> 1) as i64) ^ -((v & 1) as i64))
    }
    pub fn byte(&mut self) -> Result<u8> {
        let Some(&b) = self.buf.get(self.pos) else { bail!("byte past end") };
        self.pos += 1;
        Ok(b)
    }
    pub fn is_empty(&self) -> bool {
        self.pos >= self.buf.len()
    }
    pub fn position(&self) -> usize {
        self.pos
    }
    pub fn set_position(&mut self, pos: usize) {
        self.pos = pos;
    }
    pub fn take(&mut self, len: usize) -> Result<&'a [u8]> {
        let end = self.pos.checked_add(len).context("cursor length overflow")?;
        let Some(slice) = self.buf.get(self.pos..end) else { bail!("slice past end") };
        self.pos = end;
        Ok(slice)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archive::{put_svarint, put_varint};

    #[derive(Clone, Debug)]
    struct TestEntryLayout {
        directory_start: usize,
        directory_end: usize,
        name_start: usize,
        offset_start: usize,
        raw_len_start: usize,
        compressed_len_start: usize,
        digest_start: usize,
        inline_start: usize,
        payload_start: usize,
        payload_end: usize,
    }

    fn test_path(tag: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        std::env::temp_dir().join(format!(
            "aie-{tag}-{}-{}.aie",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn write_v2(path: &Path, sections: &[(&str, &[u8])]) {
        let mut writer = SectionWriter::create(path, 1).unwrap();
        for (name, raw) in sections {
            writer.section(name, raw).unwrap();
        }
        writer.finish().unwrap();
    }

    fn test_v2_layout(bytes: &[u8]) -> (usize, usize, Vec<TestEntryLayout>) {
        let footer = bytes.len() - V2_FOOTER_LEN as usize;
        let directory =
            u64::from_le_bytes(bytes[footer..footer + 8].try_into().unwrap()) as usize;
        let count = u32::from_le_bytes(bytes[directory..directory + 4].try_into().unwrap());
        let mut cursor = directory + 4;
        let mut entries = Vec::new();
        for _ in 0..count {
            let entry_start = cursor;
            let name_len = bytes[cursor] as usize;
            cursor += 1;
            let name_start = cursor;
            cursor += name_len;
            let offset_start = cursor;
            let offset =
                u64::from_le_bytes(bytes[cursor..cursor + 8].try_into().unwrap()) as usize;
            cursor += 8;
            let raw_len_start = cursor;
            cursor += 8;
            let compressed_len_start = cursor;
            let compressed_len =
                u64::from_le_bytes(bytes[cursor..cursor + 8].try_into().unwrap()) as usize;
            cursor += 8;
            let digest_start = cursor;
            cursor += 32;
            let payload_start = offset + 1 + name_len + 16;
            entries.push(TestEntryLayout {
                directory_start: entry_start,
                directory_end: cursor,
                name_start,
                offset_start,
                raw_len_start,
                compressed_len_start,
                digest_start,
                inline_start: offset,
                payload_start,
                payload_end: payload_start + compressed_len,
            });
        }
        assert_eq!(cursor, footer);
        (directory, footer, entries)
    }

    fn assert_open_rejected(tag: &str, bytes: &[u8]) {
        let path = test_path(tag);
        std::fs::write(&path, bytes).unwrap();
        assert!(SectionReader::open(&path).is_err(), "{tag} was accepted");
        std::fs::remove_file(path).ok();
    }

    fn test_hex(digest: [u8; 32]) -> String {
        digest.iter().map(|byte| format!("{byte:02x}")).collect()
    }

    fn write_v1(path: &Path, sections: &[(&str, &[u8])]) {
        use std::io::{Seek, Write};
        let mut file = std::fs::File::create(path).unwrap();
        file.write_all(MAGIC).unwrap();
        file.write_all(&SEEKABLE_VERSION.to_le_bytes()).unwrap();
        let mut directory = Vec::new();
        for (name, raw) in sections {
            let compressed = compress(raw, 1).unwrap();
            let offset = file.stream_position().unwrap();
            file.write_all(&[name.len() as u8]).unwrap();
            file.write_all(name.as_bytes()).unwrap();
            file.write_all(&(raw.len() as u64).to_le_bytes()).unwrap();
            file.write_all(&(compressed.len() as u64).to_le_bytes()).unwrap();
            file.write_all(&compressed).unwrap();
            directory.push((name.to_string(), offset, raw.len() as u64, compressed.len() as u64));
        }
        file.write_all(&[0]).unwrap();
        let directory_offset = file.stream_position().unwrap();
        file.write_all(&(directory.len() as u32).to_le_bytes()).unwrap();
        for (name, offset, raw_len, comp_len) in directory {
            file.write_all(&[name.len() as u8]).unwrap();
            file.write_all(name.as_bytes()).unwrap();
            file.write_all(&offset.to_le_bytes()).unwrap();
            file.write_all(&raw_len.to_le_bytes()).unwrap();
            file.write_all(&comp_len.to_le_bytes()).unwrap();
        }
        file.write_all(&directory_offset.to_le_bytes()).unwrap();
        file.write_all(b"AIED").unwrap();
    }

    fn resign_v2(bytes: &mut [u8]) {
        let footer = bytes.len() - V2_FOOTER_LEN as usize;
        let directory_offset =
            u64::from_le_bytes(bytes[footer..footer + 8].try_into().unwrap()) as usize;
        let mut hasher = root_hasher(directory_offset as u64);
        hasher.update(&bytes[directory_offset..footer]);
        bytes[footer + 8..footer + 40].copy_from_slice(hasher.finalize().as_bytes());
    }

    #[test]
    fn sections_roundtrip() {
        let dir = std::env::temp_dir().join(format!("aie-fmt-{}.aie", std::process::id()));
        let mut w = SectionWriter::create(&dir, 3).unwrap();
        w.section("alpha", b"hello world hello world").unwrap();
        w.section("beta", &[0u8; 1000]).unwrap();
        let acct = w.finish().unwrap();
        assert_eq!(acct.len(), 2);
        let back = read_sections(&dir).unwrap();
        assert_eq!(back[0].0, "alpha");
        assert_eq!(back[0].1, b"hello world hello world");
        assert_eq!(back[1].1, vec![0u8; 1000]);
        let mut seekable = SectionReader::open(&dir).unwrap();
        assert_eq!(seekable.archive_version(), VERSION);
        assert!(seekable.content_commitment().is_some());
        assert!(seekable.encoded_content_identity().unwrap().is_some());
        assert_eq!(seekable.read("alpha").unwrap(), b"hello world hello world");
        assert!(seekable.verify_all_payloads().unwrap() > 0);
        std::fs::remove_file(&dir).ok();
    }

    #[test]
    fn finish_with_file_retains_produced_archive_after_path_replacement() {
        let path = test_path("finish-held-path");
        let displaced = test_path("finish-held-displaced");
        let mut writer = SectionWriter::create_new(&path, 1).unwrap();
        writer.section("alpha", b"held archive payload").unwrap();
        let (accounting, file, commitment) = writer.finish_with_file().unwrap();
        assert_eq!(accounting.len(), 1);

        std::fs::rename(&path, &displaced).unwrap();
        std::fs::write(&path, b"concurrent replacement").unwrap();
        let mut reader = SectionReader::from_file(file.try_clone().unwrap()).unwrap();
        assert_eq!(reader.content_commitment(), Some(commitment));
        assert_eq!(reader.read("alpha").unwrap(), b"held archive payload");
        assert_eq!(std::fs::read(&path).unwrap(), b"concurrent replacement");

        std::fs::remove_file(path).ok();
        std::fs::remove_file(displaced).ok();
    }

    #[test]
    fn rooted_empty_many_unknown_and_maximum_name_sections_roundtrip_deterministically() {
        let empty = test_path("empty");
        SectionWriter::create(&empty, 1).unwrap().finish().unwrap();
        let empty_reader = SectionReader::open(&empty).unwrap();
        assert_eq!(empty_reader.section_metadata().len(), 0);
        assert_eq!(empty_reader.verify_all_payloads().unwrap(), 0);

        let first = test_path("surface-a");
        let second = test_path("surface-b");
        let maximum_name = "m".repeat(u8::MAX as usize);
        let mut owned = Vec::new();
        owned.push((maximum_name, b"maximum-name".to_vec()));
        for index in 0..128 {
            owned.push((format!("unknown.optional.{index:03}"), vec![index as u8; index + 1]));
        }
        for path in [&first, &second] {
            let mut writer = SectionWriter::create(path, 1).unwrap();
            for (name, raw) in &owned {
                writer.section(name, raw).unwrap();
            }
            writer.finish().unwrap();
        }
        assert_eq!(std::fs::read(&first).unwrap(), std::fs::read(&second).unwrap());
        let mut reader = SectionReader::open(&first).unwrap();
        assert_eq!(reader.section_metadata().len(), owned.len());
        assert_eq!(reader.read(&owned[0].0).unwrap(), owned[0].1);
        assert_eq!(reader.read("unknown.optional.127").unwrap(), owned[128].1);
        assert!(reader
            .section_metadata()
            .all(|entry| entry.compressed_blake3.is_some()));

        for path in [empty, first, second] {
            std::fs::remove_file(path).ok();
        }
    }

    #[test]
    fn rooted_declared_length_boundaries_are_fixed_and_bounded() {
        let minimum_compressed_for_maximum_raw =
            (MAX_V2_RAW_SECTION_BYTES - V2_COMPRESSION_SLACK)
                .div_ceil(MAX_V2_COMPRESSION_RATIO);
        validate_section_lengths(
            VERSION,
            "maximum",
            MAX_V2_RAW_SECTION_BYTES,
            minimum_compressed_for_maximum_raw,
        )
        .unwrap();
        validate_section_lengths(
            VERSION,
            "maximum-compressed",
            0,
            MAX_V2_COMPRESSED_SECTION_BYTES,
        )
        .unwrap();
        assert!(validate_section_lengths(
            VERSION,
            "raw-over",
            MAX_V2_RAW_SECTION_BYTES + 1,
            1
        )
        .is_err());
        assert!(validate_section_lengths(
            VERSION,
            "compressed-over",
            0,
            MAX_V2_COMPRESSED_SECTION_BYTES + 1
        )
        .is_err());
        assert!(validate_section_lengths(
            VERSION,
            "ratio-over",
            minimum_compressed_for_maximum_raw * MAX_V2_COMPRESSION_RATIO
                + V2_COMPRESSION_SLACK
                + 1,
            minimum_compressed_for_maximum_raw,
        )
        .is_err());
    }

    #[test]
    fn rooted_frozen_commitment_vectors_match_independent_implementation() {
        // These constants were produced by a separately implemented parser. Fixed
        // precompressed bytes keep the vector independent of the linked zstd encoder version.
        let path = test_path("root-vector");
        let mut writer = SectionWriter::create(&path, 1).unwrap();
        writer.section_precompressed("x", 0, &[1, 2]).unwrap();
        writer.finish().unwrap();
        let reader = SectionReader::open(&path).unwrap();
        assert_eq!(reader.file_metadata().unwrap().len(), 135);
        assert_eq!(
            reader.content_commitment().unwrap().to_hex(),
            "feacf981cb0770e7f41616bd23bc72ae45944a74376b7c1579cabf2be8cb88bf"
        );
        assert_eq!(
            test_hex(reader.encoded_content_identity().unwrap().unwrap()),
            "c298fff5a5c2dbe8a18a1098585fc6fd8007257b244dda70d6fcf77e9fc4159f"
        );
        let metadata: Vec<_> = reader.section_metadata().collect();
        assert_eq!(metadata[0].offset, 8);
        assert_eq!(metadata[0].compressed_len, 2);
        assert_eq!(
            test_hex(metadata[0].compressed_blake3.unwrap()),
            "b7d770040f780e9deff6bc038abea66e108b88d098d16d24cd7486eb671060b2"
        );
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn rooted_commitment_covers_header_directory_layout_and_footer() {
        let path = test_path("commitment-coverage");
        write_v2(&path, &[("a", b"alpha"), ("b", b"bravo")]);
        let original = std::fs::read(&path).unwrap();
        let (directory, footer, entries) = test_v2_layout(&original);

        // Each exact root preimage component is fail-closed.  Mutating any directory or root byte
        // without recomputing the root can never reach a section read.
        for index in 0..4 {
            let mut changed = original.clone();
            changed[index] ^= 1;
            assert_open_rejected(&format!("magic-{index}"), &changed);
        }
        for index in 4..8 {
            let mut changed = original.clone();
            changed[index] ^= 1;
            assert_open_rejected(&format!("version-{index}"), &changed);
        }
        for index in directory..footer {
            let mut changed = original.clone();
            changed[index] ^= 1;
            assert_open_rejected(&format!("directory-{index}"), &changed);
        }
        for index in footer..footer + 8 {
            let mut changed = original.clone();
            changed[index] ^= 1;
            assert_open_rejected(&format!("directory-offset-{index}"), &changed);
        }
        for index in footer + 8..footer + 40 {
            let mut changed = original.clone();
            changed[index] ^= 1;
            assert_open_rejected(&format!("root-{index}"), &changed);
        }
        for index in footer + 40..footer + 44 {
            let mut changed = original.clone();
            changed[index] ^= 1;
            assert_open_rejected(&format!("footer-{index}"), &changed);
        }

        // Inline headers are not duplicated in the root preimage, but every byte must agree with
        // the committed directory before open succeeds.
        for entry in &entries {
            for index in entry.inline_start..entry.payload_start {
                let mut changed = original.clone();
                changed[index] ^= 1;
                assert_open_rejected(&format!("inline-{index}"), &changed);
            }
        }
        let mut terminator = original.clone();
        terminator[directory - 1] = 1;
        assert_open_rejected("terminator", &terminator);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn rooted_authenticated_but_noncanonical_directory_edits_are_rejected() {
        let path = test_path("noncanonical-matrix");
        write_v2(&path, &[("a", b"alpha"), ("b", b"bravo")]);
        let original = std::fs::read(&path).unwrap();
        let (_, _, entries) = test_v2_layout(&original);
        assert_eq!(entries.len(), 2);

        let mut swapped = original.clone();
        let left = original[entries[0].directory_start..entries[0].directory_end].to_vec();
        let right = original[entries[1].directory_start..entries[1].directory_end].to_vec();
        assert_eq!(left.len(), right.len());
        swapped[entries[0].directory_start..entries[0].directory_end].copy_from_slice(&right);
        swapped[entries[1].directory_start..entries[1].directory_end].copy_from_slice(&left);
        resign_v2(&mut swapped);
        assert_open_rejected("entry-order", &swapped);

        let mut renamed = original.clone();
        renamed[entries[0].name_start] = b'z';
        resign_v2(&mut renamed);
        assert_open_rejected("entry-name", &renamed);

        for (tag, second_offset) in [
            ("entry-overlap", entries[1].inline_start as u64 - 1),
            ("entry-gap", entries[1].inline_start as u64 + 1),
        ] {
            let mut changed = original.clone();
            changed[entries[1].offset_start..entries[1].offset_start + 8]
                .copy_from_slice(&second_offset.to_le_bytes());
            resign_v2(&mut changed);
            assert_open_rejected(tag, &changed);
        }

        let mut raw_length = original.clone();
        let raw = u64::from_le_bytes(
            raw_length[entries[0].raw_len_start..entries[0].raw_len_start + 8]
                .try_into()
                .unwrap(),
        );
        raw_length[entries[0].raw_len_start..entries[0].raw_len_start + 8]
            .copy_from_slice(&(raw + 1).to_le_bytes());
        resign_v2(&mut raw_length);
        assert_open_rejected("entry-raw-length", &raw_length);

        let mut compressed_length = original.clone();
        let compressed = u64::from_le_bytes(
            compressed_length
                [entries[0].compressed_len_start..entries[0].compressed_len_start + 8]
                .try_into()
                .unwrap(),
        );
        compressed_length
            [entries[0].compressed_len_start..entries[0].compressed_len_start + 8]
            .copy_from_slice(&(compressed + 1).to_le_bytes());
        resign_v2(&mut compressed_length);
        assert_open_rejected("entry-compressed-length", &compressed_length);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn rooted_lazy_and_full_verification_cover_every_payload_byte_and_digest() {
        let path = test_path("payload-coverage");
        write_v2(&path, &[("a", b"aaaaaaaaaaaa"), ("b", b"bbbbbbbbbbbb")]);
        let original = std::fs::read(&path).unwrap();
        let (_, _, entries) = test_v2_layout(&original);

        for (entry_index, entry) in entries.iter().enumerate() {
            for byte_index in entry.payload_start..entry.payload_end {
                let mut changed = original.clone();
                changed[byte_index] ^= 1;
                std::fs::write(&path, &changed).unwrap();
                let mut reader = SectionReader::open(&path).unwrap();
                let selected = if entry_index == 0 { "a" } else { "b" };
                assert!(reader.read(selected).is_err(), "payload byte {byte_index} was accepted");
                assert!(reader.verify_all_payloads().is_err());
            }
        }

        let mut digest = original.clone();
        digest[entries[0].digest_start] ^= 1;
        resign_v2(&mut digest);
        std::fs::write(&path, digest).unwrap();
        let mut reader = SectionReader::open(&path).unwrap();
        assert!(reader.read("a").unwrap_err().to_string().contains("digest mismatch"));

        assert_eq!(entries[0].payload_end - entries[0].payload_start,
                   entries[1].payload_end - entries[1].payload_start);
        let mut swapped_payloads = original.clone();
        let left = original[entries[0].payload_start..entries[0].payload_end].to_vec();
        let right = original[entries[1].payload_start..entries[1].payload_end].to_vec();
        swapped_payloads[entries[0].payload_start..entries[0].payload_end]
            .copy_from_slice(&right);
        swapped_payloads[entries[1].payload_start..entries[1].payload_end]
            .copy_from_slice(&left);
        std::fs::write(&path, swapped_payloads).unwrap();
        let mut reader = SectionReader::open(&path).unwrap();
        assert!(reader.read("a").is_err());
        assert!(reader.read("b").is_err());
        assert!(reader.verify_all_payloads().is_err());

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn rooted_declared_decompressed_length_is_verified_after_authentication() {
        let path = test_path("decompressed-length");
        write_v2(&path, &[("a", b"alpha")]);
        let mut bytes = std::fs::read(&path).unwrap();
        let (_, _, entries) = test_v2_layout(&bytes);
        let entry = &entries[0];
        let declared = u64::from_le_bytes(
            bytes[entry.raw_len_start..entry.raw_len_start + 8].try_into().unwrap(),
        ) + 1;
        bytes[entry.raw_len_start..entry.raw_len_start + 8]
            .copy_from_slice(&declared.to_le_bytes());
        let inline_raw = entry.inline_start + 1 + 1;
        bytes[inline_raw..inline_raw + 8].copy_from_slice(&declared.to_le_bytes());
        resign_v2(&mut bytes);
        std::fs::write(&path, bytes).unwrap();
        let mut reader = SectionReader::open(&path).unwrap();
        let error = reader.read("a").unwrap_err();
        assert!(format!("{error:#}").contains("decompressed"));
        assert!(reader.verify_all_payloads().is_err());
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn rooted_truncation_and_trailing_bytes_are_rejected() {
        let path = test_path("extent");
        write_v2(&path, &[("a", b"alpha"), ("b", b"bravo")]);
        let original = std::fs::read(&path).unwrap();
        for keep in [0, 1, 7, 8, original.len() - 1] {
            assert_open_rejected(&format!("truncate-{keep}"), &original[..keep]);
        }
        let mut trailing = original.clone();
        trailing.extend_from_slice(b"trailing");
        assert_open_rejected("trailing", &trailing);
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn rooted_bytes_read_accounting_distinguishes_metadata_lazy_and_full_reads() {
        let path = test_path("read-accounting");
        write_v2(&path, &[("a", b"alpha"), ("b", b"bravo")]);
        let bytes = std::fs::read(&path).unwrap();
        let (directory, footer, entries) = test_v2_layout(&bytes);
        let expected_open = 8u64
            + V2_FOOTER_LEN
            + (footer - directory) as u64
            + 1
            + entries
                .iter()
                .map(|entry| (entry.payload_start - entry.inline_start) as u64)
                .sum::<u64>();
        let mut reader = SectionReader::open(&path).unwrap();
        assert_eq!(reader.bytes_read(), expected_open);
        let metadata: Vec<_> = reader
            .section_metadata()
            .map(|entry| (entry.compressed_len, entry.compressed_blake3))
            .collect();
        assert_eq!(reader.bytes_read(), expected_open);
        assert_eq!(metadata.len(), 2);
        assert!(metadata.iter().all(|entry| entry.1.is_some()));
        reader.encoded_content_identity().unwrap().unwrap();
        assert_eq!(reader.bytes_read(), expected_open);
        reader.read("a").unwrap();
        assert_eq!(reader.bytes_read(), expected_open + metadata[0].0);
        let verified = reader.verify_all_payloads().unwrap();
        assert_eq!(verified, metadata.iter().map(|entry| entry.0).sum::<u64>());
        assert_eq!(
            reader.bytes_read(),
            expected_open + metadata[0].0 + verified
        );
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn v1_remains_readable_and_sealed_content_identity_is_stable() {
        let v1 = std::env::temp_dir().join(format!("aie-v1-{}.aie", std::process::id()));
        let v2 = std::env::temp_dir().join(format!("aie-v2-{}.aie", std::process::id()));
        let sections: [(&str, &[u8]); 2] = [("alpha", b"one"), ("beta", b"two two")];
        write_v1(&v1, &sections);
        let mut legacy = SectionReader::open(&v1).unwrap();
        assert_eq!(legacy.archive_version(), SEEKABLE_VERSION);
        assert!(legacy.content_commitment().is_none());
        assert_eq!(legacy.read("beta").unwrap(), b"two two");
        let legacy_id = legacy.scan_legacy_identities().unwrap();

        let mut writer = SectionWriter::create(&v2, 1).unwrap();
        for (name, _, raw_len, _) in legacy.entries().to_vec() {
            let (compressed, observed_raw) = legacy.read_compressed(&name).unwrap();
            assert_eq!(observed_raw as u64, raw_len);
            writer.section_precompressed(&name, raw_len, &compressed).unwrap();
        }
        writer.finish().unwrap();
        let rooted = SectionReader::open(&v2).unwrap();
        assert_eq!(
            rooted.encoded_content_identity().unwrap().unwrap(),
            legacy_id.encoded_sections_blake3
        );
        std::fs::remove_file(v1).ok();
        std::fs::remove_file(v2).ok();
    }

    #[test]
    fn v1_seekable_paths_enforce_the_same_count_length_and_ratio_bounds() {
        let source = test_path("v1-bounds-source");
        write_v1(&source, &[("x", b"payload")]);
        let original = std::fs::read(&source).unwrap();
        let footer = original.len() - V1_FOOTER_LEN as usize;
        let directory =
            u64::from_le_bytes(original[footer..footer + 8].try_into().unwrap()) as usize;

        let mut count = original.clone();
        count[directory..directory + 4]
            .copy_from_slice(&((MAX_V2_SECTIONS as u32) + 1).to_le_bytes());
        assert_open_rejected("v1-count-bound", &count);

        let mut raw = original.clone();
        let bad_raw = MAX_V2_RAW_SECTION_BYTES + 1;
        // One-byte name: inline raw starts at 10; directory raw starts after count/name/offset.
        raw[10..18].copy_from_slice(&bad_raw.to_le_bytes());
        raw[directory + 14..directory + 22].copy_from_slice(&bad_raw.to_le_bytes());
        assert_open_rejected("v1-raw-bound", &raw);
        let raw_path = test_path("v1-raw-linear-api");
        std::fs::write(&raw_path, raw).unwrap();
        assert!(read_sections(&raw_path).is_err());
        std::fs::remove_file(raw_path).ok();

        let mut compressed = original.clone();
        let bad_compressed = MAX_V2_COMPRESSED_SECTION_BYTES + 1;
        compressed[18..26].copy_from_slice(&bad_compressed.to_le_bytes());
        compressed[directory + 22..directory + 30]
            .copy_from_slice(&bad_compressed.to_le_bytes());
        assert_open_rejected("v1-compressed-bound", &compressed);

        let mut ratio = original.clone();
        let compressed_len = u64::from_le_bytes(ratio[18..26].try_into().unwrap());
        let bad_ratio = compressed_len * MAX_V2_COMPRESSION_RATIO + V2_COMPRESSION_SLACK + 1;
        ratio[10..18].copy_from_slice(&bad_ratio.to_le_bytes());
        ratio[directory + 14..directory + 22].copy_from_slice(&bad_ratio.to_le_bytes());
        assert_open_rejected("v1-ratio-bound", &ratio);

        std::fs::remove_file(source).ok();
    }

    #[test]
    fn cursor_roundtrips_varints() {
        let mut buf = Vec::new();
        for v in [0u64, 1, 127, 128, 16384, u32::MAX as u64] {
            put_varint(&mut buf, v);
        }
        for v in [0i64, -1, 1, -300, 1 << 40] {
            put_svarint(&mut buf, v);
        }
        let mut c = Cursor::new(&buf);
        for v in [0u64, 1, 127, 128, 16384, u32::MAX as u64] {
            assert_eq!(c.varint().unwrap(), v);
        }
        for v in [0i64, -1, 1, -300, 1 << 40] {
            assert_eq!(c.svarint().unwrap(), v);
        }
        assert!(c.is_empty());
    }

    #[test]
    fn bad_magic_is_rejected() {
        let p = std::env::temp_dir().join(format!("aie-bad-{}.aie", std::process::id()));
        std::fs::write(&p, b"NOPE").unwrap();
        assert!(read_sections(&p).is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn seekable_reader_rejects_future_versions() {
        let p = std::env::temp_dir().join(format!("aie-future-{}.aie", std::process::id()));
        let mut w = SectionWriter::create(&p, 1).unwrap();
        w.section("x", b"x").unwrap();
        w.finish().unwrap();
        let mut bytes = std::fs::read(&p).unwrap();
        bytes[4..8].copy_from_slice(&(VERSION + 1).to_le_bytes());
        std::fs::write(&p, bytes).unwrap();
        assert!(SectionReader::open(&p).err().unwrap().to_string().contains("newer"));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn seekable_reader_rejects_directory_extent_outside_file() {
        let p = std::env::temp_dir().join(format!("aie-extent-{}.aie", std::process::id()));
        let mut w = SectionWriter::create(&p, 1).unwrap();
        w.section("x", b"payload").unwrap();
        w.finish().unwrap();
        let mut bytes = std::fs::read(&p).unwrap();
        let footer = bytes.len() - V2_FOOTER_LEN as usize;
        let dir = u64::from_le_bytes(bytes[footer..footer + 8].try_into().unwrap()) as usize;
        // Directory entry: count (4), name length (1), name (1), offset (8), raw (8), comp (8).
        bytes[dir + 22..dir + 30].copy_from_slice(&u64::MAX.to_le_bytes());
        std::fs::write(&p, bytes).unwrap();
        assert!(SectionReader::open(&p).is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn rooted_reader_rejects_noncanonical_extents_even_with_a_recomputed_root() {
        let p = std::env::temp_dir().join(format!("aie-gap-{}.aie", std::process::id()));
        let mut w = SectionWriter::create(&p, 1).unwrap();
        w.section("x", b"payload").unwrap();
        w.finish().unwrap();
        let mut bytes = std::fs::read(&p).unwrap();
        let footer = bytes.len() - V2_FOOTER_LEN as usize;
        let dir = u64::from_le_bytes(bytes[footer..footer + 8].try_into().unwrap()) as usize;
        bytes[dir + 6..dir + 14].copy_from_slice(&9u64.to_le_bytes());
        resign_v2(&mut bytes);
        std::fs::write(&p, bytes).unwrap();
        assert!(SectionReader::open(&p)
            .err()
            .unwrap()
            .to_string()
            .contains("contiguous"));
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn rooted_selected_and_full_reads_reject_payload_corruption() {
        let p = std::env::temp_dir().join(format!("aie-payload-{}.aie", std::process::id()));
        let mut w = SectionWriter::create(&p, 1).unwrap();
        w.section("x", b"payload").unwrap();
        w.finish().unwrap();
        let reader = SectionReader::open(&p).unwrap();
        let (_, offset, _, comp_len) = reader.entries()[0].clone();
        assert!(comp_len > 0);
        let payload_offset = offset + 1 + 1 + 16;
        let mut bytes = std::fs::read(&p).unwrap();
        bytes[payload_offset as usize] ^= 0x01;
        std::fs::write(&p, bytes).unwrap();
        let mut reader = SectionReader::open(&p).unwrap();
        assert!(reader.read("x").unwrap_err().to_string().contains("digest mismatch"));
        assert!(reader
            .read_compressed("x")
            .unwrap_err()
            .to_string()
            .contains("digest mismatch"));
        #[cfg(unix)]
        assert!(reader
            .read_compressed_at("x")
            .unwrap_err()
            .to_string()
            .contains("digest mismatch"));
        assert!(reader
            .verify_all_payloads()
            .unwrap_err()
            .to_string()
            .contains("digest mismatch"));
        assert!(read_sections(&p).is_err());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn rooted_reader_rejects_declared_allocation_and_compression_bombs() {
        let p = std::env::temp_dir().join(format!("aie-bomb-{}.aie", std::process::id()));
        let mut w = SectionWriter::create(&p, 1).unwrap();
        w.section("x", b"payload").unwrap();
        w.finish().unwrap();
        let mut bytes = std::fs::read(&p).unwrap();
        let footer = bytes.len() - V2_FOOTER_LEN as usize;
        let dir = u64::from_le_bytes(bytes[footer..footer + 8].try_into().unwrap()) as usize;
        let bad_raw = MAX_V2_RAW_SECTION_BYTES + 1;
        // Mutate both the committed directory and inline header, then recompute the root: the
        // parser must reject the bound itself rather than relying on authentication failure.
        bytes[dir + 14..dir + 22].copy_from_slice(&bad_raw.to_le_bytes());
        bytes[10..18].copy_from_slice(&bad_raw.to_le_bytes());
        resign_v2(&mut bytes);
        std::fs::write(&p, bytes).unwrap();
        assert!(SectionReader::open(&p)
            .err()
            .unwrap()
            .to_string()
            .contains("safety limit"));
        std::fs::remove_file(p).ok();

        let p = std::env::temp_dir().join(format!("aie-comp-bomb-{}.aie", std::process::id()));
        let mut w = SectionWriter::create(&p, 1).unwrap();
        w.section("x", b"payload").unwrap();
        w.finish().unwrap();
        let mut bytes = std::fs::read(&p).unwrap();
        let footer = bytes.len() - V2_FOOTER_LEN as usize;
        let dir = u64::from_le_bytes(bytes[footer..footer + 8].try_into().unwrap()) as usize;
        let bad_comp = MAX_V2_COMPRESSED_SECTION_BYTES + 1;
        bytes[dir + 22..dir + 30].copy_from_slice(&bad_comp.to_le_bytes());
        bytes[18..26].copy_from_slice(&bad_comp.to_le_bytes());
        resign_v2(&mut bytes);
        std::fs::write(&p, bytes).unwrap();
        assert!(SectionReader::open(&p)
            .err()
            .unwrap()
            .to_string()
            .contains("safety limit"));
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn writer_rejects_duplicate_and_oversized_names() {
        let p = std::env::temp_dir().join(format!("aie-names-{}.aie", std::process::id()));
        let mut w = SectionWriter::create(&p, 1).unwrap();
        w.section("x", b"one").unwrap();
        assert!(w.section("x", b"two").is_err());
        assert!(w.section(&"z".repeat(256), b"three").is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn cursor_rejects_overlong_and_overflowing_varints() {
        assert!(Cursor::new(&[0x80; 10]).varint().is_err());
        let mut overflow = [0xff; 10];
        overflow[9] = 0x02;
        assert!(Cursor::new(&overflow).varint().is_err());
    }
}
