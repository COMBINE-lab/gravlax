//! The annotation compiler: GTF → transcript/gene models arranged for molecule-compatibility
//! queries.
//!
//! This is the half of replay that changes per annotation; the archive side never changes. The
//! model deliberately contains only what STAR's own assignment consults: each transcript's strand
//! and exon structure, and each gene's identity. It contains no sequence-derived fields.

pub mod assign;
pub mod intent;

use anyhow::{bail, Context, Result};
use hashbrown::HashMap as FxHashMap;
use std::io::{Read, Write};
use std::path::Path;

const COMPILED_MAGIC: [u8; 8] = *b"GRVLXAIC";
const COMPILED_VERSION: u32 = 2;
const COMPILED_VERSION_WITHOUT_IDENTIFIERS: u32 = 1;
const COMPILED_HEADER_LEN: usize = 8 + 4 + 8 + 32;

/// Hash exactly the bytes consumed by an annotation parser. This keeps content binding on the
/// same open file description without adding a separate full-file read before parsing.
struct DigestingReader<R> {
    inner: R,
    hasher: blake3::Hasher,
}

impl<R> DigestingReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            hasher: blake3::Hasher::new(),
        }
    }

    fn digest(self) -> String {
        format!("blake3:{}", self.hasher.finalize().to_hex())
    }
}

impl<R: Read> Read for DigestingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.hasher.update(&buffer[..read]);
        Ok(read)
    }
}

/// Half-open exon `[start, end)`, 0-based, matching `evidence_io::Block` conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Exon {
    pub start: u32,
    pub end: u32,
}

/// One exon record exactly as it appeared in the source annotation. Assignment uses the merged
/// exons in [`Transcript::exons`]; identifier resolution uses these unmerged records so an
/// `exon_id` never acquires the union of a touching or overlapping neighbour's coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceExon {
    pub interval: Exon,
    pub id: Option<String>,
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub gene: u32,
    pub chrom: u32,
    pub strand_rev: bool,
    /// Exons sorted by start; genomic order regardless of transcript orientation.
    pub exons: Vec<Exon>,
}

impl Transcript {
    pub fn span(&self) -> (u32, u32) {
        (
            self.exons.first().map(|e| e.start).unwrap_or(0),
            self.exons.last().map(|e| e.end).unwrap_or(0),
        )
    }
}

pub struct Annotation {
    pub gene_ids: Vec<String>,
    pub gene_names: Vec<String>,
    /// Stable transcript identifiers, parallel to [`Self::transcripts`]. Legacy compiled v1
    /// annotations contain `None` because that format did not retain transcript identity.
    pub transcript_ids: Vec<Option<String>>,
    /// Exact, unmerged source exon records, parallel to [`Self::transcripts`]. Legacy compiled v1
    /// annotations contain empty vectors because that format did not retain exon identity.
    pub source_exons: Vec<Vec<SourceExon>>,
    pub transcripts: Vec<Transcript>,
    /// Chromosome name → id, as encountered in the GTF. Callers must translate BAM reference ids
    /// through names, never assume the orders coincide.
    pub chrom_ids: FxHashMap<String, u32>,
    /// Per chromosome: transcript indices sorted by span start, with a running maximum of span
    /// ends alongside — the classic augmented array for stabbing queries without a tree.
    index: FxHashMap<u32, ChromIndex>,
}

#[derive(Debug, PartialEq, Eq)]
struct ChromIndex {
    starts: Vec<u32>,
    /// Span end per entry, parallel to `starts` — the overlap test reads these flat arrays
    /// instead of chasing into each Transcript (a cache miss per candidate in the hot walk).
    ends: Vec<u32>,
    max_end_prefix: Vec<u32>,
    tx: Vec<u32>,
}

fn build_index(transcripts: &[Transcript]) -> FxHashMap<u32, ChromIndex> {
    let mut index: FxHashMap<u32, ChromIndex> = FxHashMap::new();
    let mut by_chrom: FxHashMap<u32, Vec<u32>> = FxHashMap::new();
    for (i, t) in transcripts.iter().enumerate() {
        by_chrom.entry(t.chrom).or_default().push(i as u32);
    }
    for (chrom, mut txs) in by_chrom {
        txs.sort_by_key(|&i| transcripts[i as usize].span().0);
        let starts: Vec<u32> = txs.iter().map(|&i| transcripts[i as usize].span().0).collect();
        let ends: Vec<u32> = txs.iter().map(|&i| transcripts[i as usize].span().1).collect();
        let mut max_end_prefix = Vec::with_capacity(txs.len());
        let mut m = 0u32;
        for &i in &txs {
            m = m.max(transcripts[i as usize].span().1);
            max_end_prefix.push(m);
        }
        index.insert(chrom, ChromIndex { starts, ends, max_end_prefix, tx: txs });
    }
    index
}

fn push_u32(out: &mut Vec<u8>, value: usize, label: &str) -> Result<()> {
    out.extend_from_slice(
        &u32::try_from(value)
            .with_context(|| format!("{label} exceeds compiled-annotation u32 limit"))?
            .to_le_bytes(),
    );
    Ok(())
}

fn push_string(out: &mut Vec<u8>, value: &str, label: &str) -> Result<()> {
    push_u32(out, value.len(), label)?;
    out.extend_from_slice(value.as_bytes());
    Ok(())
}

struct CompiledCursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> CompiledCursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, at: 0 }
    }

    fn take(&mut self, n: usize, label: &str) -> Result<&'a [u8]> {
        let end = self.at.checked_add(n).context("compiled annotation offset overflow")?;
        if end > self.bytes.len() {
            bail!("truncated compiled annotation while reading {label}");
        }
        let out = &self.bytes[self.at..end];
        self.at = end;
        Ok(out)
    }

    fn u8(&mut self, label: &str) -> Result<u8> {
        Ok(self.take(1, label)?[0])
    }

    fn u32(&mut self, label: &str) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4, label)?.try_into().unwrap()))
    }

    fn count(&mut self, label: &str, minimum_bytes_per_item: usize) -> Result<usize> {
        let count = self.u32(label)? as usize;
        let remaining = self.bytes.len() - self.at;
        if count > remaining / minimum_bytes_per_item {
            bail!(
                "compiled annotation {label} {count} cannot fit in the remaining {remaining} bytes"
            );
        }
        Ok(count)
    }

    fn string(&mut self, label: &str) -> Result<String> {
        let n = self.u32(label)? as usize;
        let bytes = self.take(n, label)?;
        Ok(std::str::from_utf8(bytes)
            .with_context(|| format!("compiled annotation {label} is not UTF-8"))?
            .to_owned())
    }

    fn finish(self) -> Result<()> {
        if self.at != self.bytes.len() {
            bail!("compiled annotation payload has {} trailing bytes", self.bytes.len() - self.at);
        }
        Ok(())
    }
}

#[cfg_attr(not(test), allow(dead_code))]
fn attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    // str shim over the byte scanner, kept for the tests.
    attr_b(attrs.as_bytes(), key.as_bytes()).and_then(|v| std::str::from_utf8(v).ok())
}

/// Decimal u32 from ASCII bytes; `None` on anything else (mirrors `str::parse::<u32>` failures).
fn parse_u32(b: &[u8]) -> Option<u32> {
    if b.is_empty() {
        return None;
    }
    let mut v: u32 = 0;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((c - b'0') as u32)?;
    }
    Some(v)
}

/// GTF attributes: `key "value"; key "value";`. A byte scan beats both a regex and per-line
/// UTF-8 validation here — this parser sees ~1.7M exon lines per annotation (GTF is ASCII).
fn attr_b<'a>(attrs: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    let mut rest = attrs;
    while let Some(p) = rest.windows(key.len()).position(|w| w == key) {
        let after = &rest[p + key.len()..];
        if p == 0 || rest[p - 1] == b';' || rest[p - 1].is_ascii_whitespace() {
            let after = &after[after.iter().take_while(|b| b.is_ascii_whitespace()).count()..];
            if let Some(stripped) = after.strip_prefix(b"\"") {
                if let Some(q) = stripped.iter().position(|&b| b == b'"') {
                    return Some(&stripped[..q]);
                }
            }
        }
        rest = &rest[p + key.len()..];
    }
    None
}

impl Annotation {
    /// Load either an ordinary uncompressed GTF or a guarded Gravlax compiled annotation.
    /// Existing `--gtf` CLI arguments call this method so compiled artifacts are transparent.
    pub fn from_path(path: &Path) -> Result<Annotation> {
        let file =
            std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        Self::from_open_file(file, path)
    }

    /// Parse an annotation from an already-open file. Keeping format detection and parsing on
    /// this descriptor lets callers bind a digest to the exact snapshot that is parsed, even if
    /// the pathname is atomically replaced while loading.
    pub(crate) fn from_open_file(mut file: std::fs::File, path: &Path) -> Result<Annotation> {
        use std::io::{Read, Seek, SeekFrom};
        let mut magic = [0u8; 8];
        let n = file.read(&mut magic)?;
        file.seek(SeekFrom::Start(0))?;
        if n == magic.len() && magic == COMPILED_MAGIC {
            return Self::from_compiled_file(file, path);
        }
        if path.extension().is_some_and(|ext| ext == "aic") {
            bail!("{} does not have Gravlax compiled-annotation magic", path.display());
        }
        Self::from_gtf_file(file)
    }

    /// Parse and content-bind one annotation snapshot in a single pass over the already-open
    /// descriptor. The returned digest uses the public `blake3:<hex>` identity syntax.
    pub(crate) fn from_open_file_with_digest(
        mut file: std::fs::File,
        path: &Path,
    ) -> Result<(Annotation, String)> {
        use std::io::{Seek, SeekFrom};

        let mut magic = [0u8; 8];
        let n = file.read(&mut magic)?;
        file.seek(SeekFrom::Start(0))?;
        let mut reader = DigestingReader::new(file);
        let annotation = if n == magic.len() && magic == COMPILED_MAGIC {
            Self::from_compiled_reader(&mut reader, path)?
        } else {
            if path.extension().is_some_and(|ext| ext == "aic") {
                bail!("{} does not have Gravlax compiled-annotation magic", path.display());
            }
            Self::from_gtf_reader(&mut reader)?
        };
        Ok((annotation, reader.digest()))
    }

    fn compiled_payload(&self, include_identifiers: bool) -> Result<Vec<u8>> {
        let mut payload = Vec::new();
        if self.gene_ids.len() != self.gene_names.len() {
            bail!("annotation has different gene-id and gene-name counts");
        }
        push_u32(&mut payload, self.gene_ids.len(), "gene count")?;
        for id in &self.gene_ids {
            push_string(&mut payload, id, "gene id")?;
        }
        for name in &self.gene_names {
            push_string(&mut payload, name, "gene name")?;
        }

        let mut chrom_names = vec![None; self.chrom_ids.len()];
        for (name, &id) in &self.chrom_ids {
            let slot = chrom_names
                .get_mut(id as usize)
                .with_context(|| format!("chromosome id {id} is outside dictionary"))?;
            if slot.replace(name.as_str()).is_some() {
                bail!("duplicate compiled chromosome id {id}");
            }
        }
        push_u32(&mut payload, chrom_names.len(), "chromosome count")?;
        for (id, name) in chrom_names.into_iter().enumerate() {
            push_string(
                &mut payload,
                name.with_context(|| format!("missing chromosome id {id}"))?,
                "chromosome name",
            )?;
        }

        push_u32(&mut payload, self.transcripts.len(), "transcript count")?;
        for transcript in &self.transcripts {
            payload.extend_from_slice(&transcript.gene.to_le_bytes());
            payload.extend_from_slice(&transcript.chrom.to_le_bytes());
            payload.push(u8::from(transcript.strand_rev));
            push_u32(&mut payload, transcript.exons.len(), "exon count")?;
            for exon in &transcript.exons {
                payload.extend_from_slice(&exon.start.to_le_bytes());
                payload.extend_from_slice(&exon.end.to_le_bytes());
            }
        }

        let mut indexed_chroms: Vec<u32> = self.index.keys().copied().collect();
        indexed_chroms.sort_unstable();
        push_u32(&mut payload, indexed_chroms.len(), "indexed chromosome count")?;
        for chrom in indexed_chroms {
            payload.extend_from_slice(&chrom.to_le_bytes());
            let ix = &self.index[&chrom];
            let n = ix.tx.len();
            if ix.starts.len() != n || ix.ends.len() != n || ix.max_end_prefix.len() != n {
                bail!("annotation chromosome {chrom} has inconsistent overlap-index arrays");
            }
            push_u32(&mut payload, n, "chromosome index length")?;
            for k in 0..n {
                payload.extend_from_slice(&ix.starts[k].to_le_bytes());
                payload.extend_from_slice(&ix.ends[k].to_le_bytes());
                payload.extend_from_slice(&ix.max_end_prefix[k].to_le_bytes());
                payload.extend_from_slice(&ix.tx[k].to_le_bytes());
            }
        }

        if include_identifiers {
            if self.transcript_ids.len() != self.transcripts.len()
                || self.source_exons.len() != self.transcripts.len()
            {
                bail!("annotation identifier metadata is not parallel to its transcripts");
            }
            push_u32(
                &mut payload,
                self.transcripts.len(),
                "transcript metadata count",
            )?;
            for (tx_index, ((transcript, transcript_id), source_exons)) in self
                .transcripts
                .iter()
                .zip(&self.transcript_ids)
                .zip(&self.source_exons)
                .enumerate()
            {
                match transcript_id {
                    Some(id) => {
                        payload.push(1);
                        push_string(&mut payload, id, "transcript id")?;
                    }
                    None => payload.push(0),
                }
                push_u32(&mut payload, source_exons.len(), "source exon count")?;
                for source in source_exons {
                    if source.interval.start >= source.interval.end {
                        bail!("annotation transcript {tx_index} has an invalid source exon");
                    }
                    if !transcript.exons.iter().any(|assignment| {
                        assignment.start <= source.interval.start
                            && source.interval.end <= assignment.end
                    }) {
                        bail!(
                            "annotation transcript {tx_index} has a source exon outside its assignment exons"
                        );
                    }
                    payload.extend_from_slice(&source.interval.start.to_le_bytes());
                    payload.extend_from_slice(&source.interval.end.to_le_bytes());
                    match &source.id {
                        Some(id) => {
                            payload.push(1);
                            push_string(&mut payload, id, "exon id")?;
                        }
                        None => payload.push(0),
                    }
                }
            }
        }
        Ok(payload)
    }

    /// Serialize the fully compiled annotation model and overlap index. The payload checksum and
    /// structural validation make the artifact safe to reuse across replay/query commands.
    /// Version 2 additionally retains transcript and exon identifiers for scientific-intent
    /// resolution. The reader remains compatible with version 1 artifacts.
    pub fn write_compiled(&self, path: &Path) -> Result<()> {
        self.write_compiled_with_identity(path).map(|_| ())
    }

    /// Serialize a compiled annotation and return the payload identity already embedded in its
    /// deterministic header. This avoids rereading a potentially large artifact just to report
    /// the bytes that were written.
    pub fn write_compiled_with_identity(&self, path: &Path) -> Result<String> {
        let mut writer = std::io::BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(path)
                .with_context(|| format!("creating {}", path.display()))?,
        );
        let (identity, _) = self.write_compiled_to(&mut writer)?;
        writer.flush()?;
        Ok(identity)
    }

    /// Write a compiled annotation to an already-owned sink and return its embedded payload
    /// identity and exact byte length. Callers that publish atomically can keep their staging
    /// descriptor open through the no-clobber commit instead of reopening a replaceable pathname.
    pub fn write_compiled_to<W: Write>(&self, writer: &mut W) -> Result<(String, u64)> {
        let payload = self.compiled_payload(true)?;
        let payload_len =
            u64::try_from(payload.len()).context("compiled annotation exceeds u64")?;
        let checksum = blake3::hash(&payload);
        writer.write_all(&COMPILED_MAGIC)?;
        writer.write_all(&COMPILED_VERSION.to_le_bytes())?;
        writer.write_all(&payload_len.to_le_bytes())?;
        writer.write_all(checksum.as_bytes())?;
        writer.write_all(&payload)?;
        let bytes = (COMPILED_HEADER_LEN as u64)
            .checked_add(payload_len)
            .context("compiled annotation byte count overflow")?;
        Ok((
            format!(
                "aic-v{COMPILED_VERSION}-payload-blake3:{}",
                checksum.to_hex()
            ),
            bytes,
        ))
    }

    pub fn from_compiled(path: &Path) -> Result<Annotation> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("opening compiled annotation {}", path.display()))?;
        Self::from_compiled_file(file, path)
    }

    fn from_compiled_file(mut file: std::fs::File, path: &Path) -> Result<Annotation> {
        Self::from_compiled_reader(&mut file, path)
    }

    fn from_compiled_reader<R: Read>(mut reader: R, path: &Path) -> Result<Annotation> {
        let mut bytes = Vec::new();
        reader
            .read_to_end(&mut bytes)
            .with_context(|| format!("reading compiled annotation {}", path.display()))?;
        if bytes.len() < COMPILED_HEADER_LEN {
            bail!("truncated compiled annotation header");
        }
        if bytes[..8] != COMPILED_MAGIC {
            bail!("invalid compiled annotation magic");
        }
        let version = u32::from_le_bytes(bytes[8..12].try_into().unwrap());
        if version != COMPILED_VERSION && version != COMPILED_VERSION_WITHOUT_IDENTIFIERS {
            bail!(
                "unsupported compiled annotation version {version}; expected {} or {COMPILED_VERSION}",
                COMPILED_VERSION_WITHOUT_IDENTIFIERS
            );
        }
        let payload_len = usize::try_from(u64::from_le_bytes(bytes[12..20].try_into().unwrap()))
            .context("compiled annotation payload length exceeds usize")?;
        let expected_len = COMPILED_HEADER_LEN
            .checked_add(payload_len)
            .context("compiled annotation length overflow")?;
        if bytes.len() != expected_len {
            bail!(
                "compiled annotation length mismatch: header declares {expected_len} bytes, file has {}",
                bytes.len()
            );
        }
        let payload = &bytes[COMPILED_HEADER_LEN..];
        if blake3::hash(payload).as_bytes() != &bytes[20..52] {
            bail!("compiled annotation payload checksum mismatch");
        }

        let mut cursor = CompiledCursor::new(payload);
        // Validate all serialized counts against a conservative minimum record
        // size before using them as allocation capacities. The checksum catches
        // accidental damage; these checks also make checksum-valid hostile files
        // fail without attempting attacker-sized allocations.
        let gene_count = cursor.count("gene count", 8)?;
        let mut gene_ids = Vec::with_capacity(gene_count);
        for _ in 0..gene_count {
            gene_ids.push(cursor.string("gene id")?);
        }
        let mut gene_names = Vec::with_capacity(gene_count);
        for _ in 0..gene_count {
            gene_names.push(cursor.string("gene name")?);
        }
        let chrom_count = cursor.count("chromosome count", 4)?;
        let mut chrom_ids = FxHashMap::with_capacity(chrom_count);
        for id in 0..chrom_count {
            let name = cursor.string("chromosome name")?;
            if chrom_ids.insert(name, id as u32).is_some() {
                bail!("compiled annotation has a duplicate chromosome name");
            }
        }
        let transcript_count = cursor.count("transcript count", 13)?;
        let mut transcripts = Vec::with_capacity(transcript_count);
        for _ in 0..transcript_count {
            let gene = cursor.u32("transcript gene")?;
            let chrom = cursor.u32("transcript chromosome")?;
            if gene as usize >= gene_count || chrom as usize >= chrom_count {
                bail!("compiled annotation transcript references an invalid dictionary id");
            }
            let strand_rev = match cursor.u8("transcript strand")? {
                0 => false,
                1 => true,
                value => bail!("compiled annotation has invalid strand byte {value}"),
            };
            let exon_count = cursor.count("exon count", 8)?;
            if exon_count == 0 {
                bail!("compiled annotation transcript has no exons");
            }
            let mut exons = Vec::with_capacity(exon_count);
            for _ in 0..exon_count {
                let exon = Exon {
                    start: cursor.u32("exon start")?,
                    end: cursor.u32("exon end")?,
                };
                if exon.start >= exon.end
                    || exons.last().is_some_and(|previous: &Exon| previous.end >= exon.start)
                {
                    bail!("compiled annotation has invalid or unsorted exons");
                }
                exons.push(exon);
            }
            transcripts.push(Transcript { gene, chrom, strand_rev, exons });
        }

        let indexed_chrom_count = cursor.count("indexed chromosome count", 8)?;
        let mut index = FxHashMap::with_capacity(indexed_chrom_count);
        for _ in 0..indexed_chrom_count {
            let chrom = cursor.u32("indexed chromosome")?;
            if chrom as usize >= chrom_count || index.contains_key(&chrom) {
                bail!("compiled annotation has an invalid or duplicate indexed chromosome");
            }
            let n = cursor.count("chromosome index length", 16)?;
            let mut ix = ChromIndex {
                starts: Vec::with_capacity(n),
                ends: Vec::with_capacity(n),
                max_end_prefix: Vec::with_capacity(n),
                tx: Vec::with_capacity(n),
            };
            for _ in 0..n {
                ix.starts.push(cursor.u32("index start")?);
                ix.ends.push(cursor.u32("index end")?);
                ix.max_end_prefix.push(cursor.u32("index prefix maximum")?);
                ix.tx.push(cursor.u32("index transcript")?);
            }
            index.insert(chrom, ix);
        }
        let (transcript_ids, source_exons) = if version == COMPILED_VERSION {
            let metadata_count = cursor.count("transcript metadata count", 5)?;
            if metadata_count != transcript_count {
                bail!(
                    "compiled annotation has {transcript_count} transcripts but {metadata_count} transcript metadata records"
                );
            }
            let mut transcript_ids = Vec::with_capacity(transcript_count);
            let mut source_exons = Vec::with_capacity(transcript_count);
            for (tx_index, transcript) in transcripts.iter().enumerate() {
                let transcript_id = match cursor.u8("transcript id presence")? {
                    0 => None,
                    1 => Some(cursor.string("transcript id")?),
                    value => {
                        bail!("compiled annotation has invalid transcript-id presence byte {value}")
                    }
                };
                let source_count = cursor.count("source exon count", 9)?;
                let mut transcript_sources = Vec::with_capacity(source_count);
                for _ in 0..source_count {
                    let interval = Exon {
                        start: cursor.u32("source exon start")?,
                        end: cursor.u32("source exon end")?,
                    };
                    if interval.start >= interval.end
                        || !transcript.exons.iter().any(|assignment| {
                            assignment.start <= interval.start && interval.end <= assignment.end
                        })
                    {
                        bail!(
                            "compiled annotation transcript {tx_index} has an invalid source exon"
                        );
                    }
                    let id = match cursor.u8("source exon id presence")? {
                        0 => None,
                        1 => Some(cursor.string("exon id")?),
                        value => {
                            bail!("compiled annotation has invalid exon-id presence byte {value}")
                        }
                    };
                    transcript_sources.push(SourceExon { interval, id });
                }
                transcript_ids.push(transcript_id);
                source_exons.push(transcript_sources);
            }
            (transcript_ids, source_exons)
        } else {
            // AIC v1 deliberately omitted transcript identity. Preserve that absence explicitly;
            // resolvers can still answer gene queries and can explain why transcript/exon lookup
            // is unavailable instead of fabricating identifiers.
            (
                vec![None; transcript_count],
                vec![Vec::new(); transcript_count],
            )
        };
        cursor.finish()?;
        if index != build_index(&transcripts) {
            bail!("compiled annotation overlap index is inconsistent with its transcripts");
        }
        Ok(Annotation {
            gene_ids,
            gene_names,
            transcript_ids,
            source_exons,
            transcripts,
            chrom_ids,
            index,
        })
    }

    /// Parse a (possibly gzip-compressed is NOT handled — pass the uncompressed GTF) annotation.
    /// Only `exon` lines are consulted; transcript and gene records are derived from them, so a
    /// GTF whose feature hierarchy is unusual still compiles as long as its exons are attributed.
    pub fn from_gtf(path: &Path) -> Result<Annotation> {
        let file =
            std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        Self::from_gtf_file(file)
    }

    /// Parse a GTF and bind the model to the exact bytes consumed from the same open descriptor.
    /// The returned identity uses the public `blake3:<hex>` syntax.
    pub fn from_gtf_with_digest(path: &Path) -> Result<(Annotation, String)> {
        let file =
            std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut reader = DigestingReader::new(file);
        let annotation = Self::from_gtf_reader(&mut reader)?;
        Ok((annotation, reader.digest()))
    }

    fn from_gtf_file(file: std::fs::File) -> Result<Annotation> {
        Self::from_gtf_reader(file)
    }

    fn from_gtf_reader<R: Read>(reader: R) -> Result<Annotation> {
        // Two-phase parallel parse. Phase 1 splits the file into line-aligned ranges and scans
        // them in parallel, extracting per-exon field slices (chrom/gene/transcript/attrs +
        // coordinates) — all the byte crunching. Phase 2 walks the records serially IN FILE
        // ORDER, interning ids by first occurrence — identical id assignment to the old
        // line-by-line loop at a fraction of the wall (the serial parse was ~1.5 s of every
        // --gtf command). GTF is ASCII; strings are validated only when first interned.
        struct Rec<'a> {
            chrom: &'a [u8],
            gid: &'a [u8],
            tid: &'a [u8],
            attrs: &'a [u8],
            start: u32,
            end: u32,
            rev: bool,
        }

        let mut gene_ids: Vec<String> = Vec::new();
        let mut gene_names: Vec<String> = Vec::new();
        // Owning only distinct identifiers lets completed GTF batches be dropped; the former
        // borrowed-key maps pinned the whole 3.1 GB GENCODE file in memory.
        let mut gene_lookup: FxHashMap<Vec<u8>, u32> = FxHashMap::new();
        let mut chrom_lookup: FxHashMap<Vec<u8>, u32> = FxHashMap::new();
        let mut tx_lookup: FxHashMap<Vec<u8>, u32> = FxHashMap::new();
        let mut transcript_ids: Vec<Option<String>> = Vec::new();
        let mut source_exons: Vec<Vec<SourceExon>> = Vec::new();
        let mut transcripts: Vec<Transcript> = Vec::new();

        use rayon::prelude::*;
        use std::io::BufRead;
        const INPUT_BATCH: usize = 192 << 20;
        let mut reader = std::io::BufReader::with_capacity(8 << 20, reader);
        loop {
            let mut data = Vec::with_capacity(INPUT_BATCH + (1 << 20));
            let n = reader
                .by_ref()
                .take(INPUT_BATCH as u64)
                .read_to_end(&mut data)?;
            if n == 0 {
                break;
            }
            if data.last() != Some(&b'\n') {
                reader.read_until(b'\n', &mut data)?;
            }

            // Parse line-aligned subranges in parallel, then intern strictly in file order.
            let nchunks = (data.len() / (8 << 20)).clamp(1, 512);
            let mut bounds = vec![0usize];
            for k in 1..nchunks {
                let mut p = k * data.len() / nchunks;
                if p <= *bounds.last().unwrap() {
                    continue;
                }
                while p < data.len() && data[p] != b'\n' {
                    p += 1
                }
                p = (p + 1).min(data.len());
                bounds.push(p);
            }
            bounds.push(data.len());
            let parts: Vec<Vec<Rec>> = bounds
                .par_windows(2)
                .map(|w| -> Result<Vec<Rec>> {
                    let mut out = Vec::new();
                    for line in data[w[0]..w[1]].split(|&b| b == b'\n') {
                        if line.is_empty() || line.first() == Some(&b'#') {
                            continue;
                        }
                        let lossy = || String::from_utf8_lossy(line).into_owned();
                        let mut f = line.split(|&b| b == b'\t');
                        let (Some(chrom), Some(_src), Some(feat)) = (f.next(), f.next(), f.next())
                        else {
                            bail!("malformed GTF line (fewer than 3 fields): {}", lossy())
                        };
                        if feat != b"exon" {
                            continue;
                        }
                        let (Some(start), Some(end), _, Some(strand), _, Some(attrs)) =
                            (f.next(), f.next(), f.next(), f.next(), f.next(), f.next())
                        else {
                            bail!("malformed exon line (expected 9 fields): {}", lossy())
                        };
                        let Some(gid) = attr_b(attrs, b"gene_id") else {
                            bail!("exon without gene_id: {}", lossy())
                        };
                        let Some(tid) = attr_b(attrs, b"transcript_id") else {
                            bail!("exon without transcript_id: {}", lossy())
                        };
                        let start_1: u32 = parse_u32(start).context("invalid exon start")?;
                        let e: u32 = parse_u32(end).context("exon end")?;
                        if start_1 == 0 || e < start_1 {
                            bail!("invalid exon interval {start_1}..{e}: {}", lossy());
                        }
                        let rev = match strand {
                            b"+" => false,
                            b"-" => true,
                            _ => bail!(
                                "invalid exon strand {:?}: {}",
                                String::from_utf8_lossy(strand),
                                lossy()
                            ),
                        };
                        out.push(Rec {
                            chrom,
                            gid,
                            tid,
                            attrs,
                            start: start_1 - 1,
                            end: e,
                            rev,
                        });
                    }
                    Ok(out)
                })
                .collect::<Result<_>>()?;

            for r in parts.iter().flatten() {
                let chrom_id = match chrom_lookup.get(r.chrom) {
                    Some(&id) => id,
                    None => {
                        std::str::from_utf8(r.chrom).context("GTF chromosome name is not UTF-8")?;
                        let id = chrom_lookup.len() as u32;
                        chrom_lookup.insert(r.chrom.to_vec(), id);
                        id
                    }
                };
                let gene_name_bytes = attr_b(r.attrs, b"gene_name").unwrap_or(r.gid);
                let gene_name =
                    std::str::from_utf8(gene_name_bytes).context("GTF gene_name is not UTF-8")?;
                let gene = match gene_lookup.get(r.gid) {
                    Some(&g) => g,
                    None => {
                        let gene_id =
                            std::str::from_utf8(r.gid).context("GTF gene_id is not UTF-8")?;
                        let g = gene_ids.len() as u32;
                        gene_lookup.insert(r.gid.to_vec(), g);
                        gene_ids.push(gene_id.to_owned());
                        gene_names.push(gene_name.to_owned());
                        g
                    }
                };
                let txi = match tx_lookup.get(r.tid) {
                    Some(&t) => {
                        let transcript = &transcripts[t as usize];
                        if (transcript.gene, transcript.chrom, transcript.strand_rev)
                            != (gene, chrom_id, r.rev)
                        {
                            let transcript_id = transcript_ids[t as usize]
                                .as_deref()
                                .unwrap_or("<unavailable>");
                            bail!(
                                "transcript_id {transcript_id:?} is reused with a different gene, chromosome, or strand"
                            );
                        }
                        t
                    }
                    None => {
                        let transcript_id =
                            std::str::from_utf8(r.tid).context("GTF transcript_id is not UTF-8")?;
                        let t = transcripts.len() as u32;
                        tx_lookup.insert(r.tid.to_vec(), t);
                        transcript_ids.push(Some(transcript_id.to_owned()));
                        source_exons.push(Vec::new());
                        transcripts.push(Transcript {
                            gene,
                            chrom: chrom_id,
                            strand_rev: r.rev,
                            exons: Vec::new(),
                        });
                        t
                    }
                };
                let exon_id = attr_b(r.attrs, b"exon_id")
                    .map(|id| {
                        std::str::from_utf8(id)
                            .context("GTF exon_id is not UTF-8")
                            .map(str::to_owned)
                    })
                    .transpose()?;
                transcripts[txi as usize].exons.push(Exon {
                    start: r.start,
                    end: r.end,
                });
                source_exons[txi as usize].push(SourceExon {
                    interval: Exon {
                        start: r.start,
                        end: r.end,
                    },
                    id: exon_id,
                });
            }
        }

        let chrom_ids: FxHashMap<String, u32> = chrom_lookup
            .into_iter()
            .map(|(name, id)| {
                String::from_utf8(name)
                    .context("GTF chromosome name is not UTF-8")
                    .map(|name| (name, id))
            })
            .collect::<Result<_>>()?;

        for (tx_index, t) in transcripts.iter_mut().enumerate() {
            source_exons[tx_index].sort_by(|left, right| {
                (left.interval.start, left.interval.end, left.id.as_deref()).cmp(&(
                    right.interval.start,
                    right.interval.end,
                    right.id.as_deref(),
                ))
            });
            let mut assignment_exons = std::mem::take(&mut t.exons);
            assignment_exons.sort_by_key(|exon| (exon.start, exon.end));
            // Merge touching/overlapping exons: some GENCODE transcripts carry duplicate exon
            // records, and the assignment walk assumes disjoint sorted exons.
            let mut merged: Vec<Exon> = Vec::with_capacity(assignment_exons.len());
            for e in assignment_exons {
                match merged.last_mut() {
                    Some(last) if e.start <= last.end => {
                        last.end = last.end.max(e.end);
                    }
                    _ => merged.push(e),
                }
            }
            t.exons = merged;
        }

        let index = build_index(&transcripts);

        Ok(Annotation {
            gene_ids,
            gene_names,
            transcript_ids,
            source_exons,
            transcripts,
            chrom_ids,
            index,
        })
    }

    /// Transcript indices whose span overlaps `[qs, qe)` on `chrom`. The prefix-max lets the
    /// leftward walk stop at the first position where nothing further left can still reach `qs`.
    pub fn overlapping(&self, chrom: u32, qs: u32, qe: u32) -> Vec<u32> {
        let mut out = Vec::new();
        self.overlapping_into(chrom, qs, qe, &mut out);
        out
    }

    /// Allocation-free variant of [`Self::overlapping`]: clears `out` and appends the hits in the
    /// same order. The per-row classifiers call this millions of times per replay.
    pub fn overlapping_into(&self, chrom: u32, qs: u32, qe: u32, out: &mut Vec<u32>) {
        out.clear();
        let Some(ix) = self.index.get(&chrom) else {
            return;
        };
        let hi = ix.starts.partition_point(|&s| s < qe);
        for k in (0..hi).rev() {
            if ix.max_end_prefix[k] <= qs {
                break;
            }
            // starts[k] < qe holds by construction of hi; only the end needs testing.
            if ix.ends[k] > qs {
                out.push(ix.tx[k]);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_gtf(lines: &[&str]) -> temppath::TempPath {
        let mut f = temppath::NamedTemp::new();
        for l in lines {
            writeln!(f.file, "{l}").unwrap();
        }
        f.into_path()
    }

    // Minimal in-crate temp-file helper so tests need no extra dependency.
    mod temppath {
        pub struct NamedTemp {
            pub file: std::fs::File,
            path: std::path::PathBuf,
        }
        pub struct TempPath(pub std::path::PathBuf);
        impl TempPath {
            pub fn path(&self) -> &std::path::Path {
                &self.0
            }
        }
        impl Drop for TempPath {
            fn drop(&mut self) {
                let _ = std::fs::remove_file(&self.0);
            }
        }
        impl NamedTemp {
            pub fn new() -> NamedTemp {
                let path = std::env::temp_dir().join(format!(
                    "anno-test-{}-{:x}.gtf",
                    std::process::id(),
                    std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()
                ));
                NamedTemp { file: std::fs::File::create(&path).unwrap(), path }
            }
            pub fn into_path(self) -> TempPath {
                TempPath(self.path)
            }
        }
    }

    const A: &str = r#"chr1	X	exon	101	200	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "Alpha";"#;
    const B: &str = r#"chr1	X	exon	301	400	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "Alpha";"#;
    const C: &str = r#"chr1	X	exon	150	260	.	-	.	gene_id "G2"; transcript_id "T2"; gene_name "Beta";"#;

    #[test]
    fn parses_exons_into_transcripts_with_correct_coordinates() {
        let p = write_gtf(&[A, B, C]);
        let a = Annotation::from_gtf(p.path()).unwrap();
        assert_eq!(a.gene_ids, vec!["G1", "G2"]);
        assert_eq!(a.transcript_ids, vec![Some("T1".into()), Some("T2".into())]);
        assert_eq!(a.transcripts.len(), 2);
        let t1 = &a.transcripts[0];
        // 1-based inclusive 101..200 becomes 0-based half-open [100, 200).
        assert_eq!(t1.exons, vec![Exon { start: 100, end: 200 }, Exon { start: 300, end: 400 }]);
        assert!(!t1.strand_rev);
        assert!(a.transcripts[1].strand_rev);
    }

    #[test]
    fn overlap_query_respects_span_not_exons() {
        let p = write_gtf(&[A, B, C]);
        let a = Annotation::from_gtf(p.path()).unwrap();
        let chr1 = a.chrom_ids["chr1"];
        // Query inside T1's intron still returns T1 (span overlap); assignment decides exon/intron.
        let hits = a.overlapping(chr1, 250, 260);
        assert!(hits.iter().any(|&i| a.transcripts[i as usize].gene == 0));
        assert!(hits.iter().any(|&i| a.transcripts[i as usize].gene == 1));
        // Far outside everything.
        assert!(a.overlapping(chr1, 5000, 5100).is_empty());
    }

    #[test]
    fn duplicate_exon_records_are_merged() {
        let p = write_gtf(&[A, A, B]);
        let a = Annotation::from_gtf(p.path()).unwrap();
        assert_eq!(a.transcripts[0].exons.len(), 2);
    }

    #[test]
    fn attr_parser_handles_adjacent_keys() {
        let s = r#"gene_id "G1"; transcript_id "T1"; gene_name "N";"#;
        assert_eq!(attr(s, "gene_id"), Some("G1"));
        assert_eq!(attr(s, "transcript_id"), Some("T1"));
        assert_eq!(attr(s, "gene_name"), Some("N"));
        assert_eq!(attr(s, "absent"), None);
    }

    #[test]
    fn rejects_zero_based_or_reversed_gtf_intervals() {
        let zero = write_gtf(&[
            "chr1\tX\texon\t0\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";",
        ]);
        assert!(Annotation::from_gtf(zero.path()).is_err());

        let reversed = write_gtf(&[
            "chr1\tX\texon\t200\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";",
        ]);
        assert!(Annotation::from_gtf(reversed.path()).is_err());
    }

    #[test]
    fn rejects_truncated_exons_and_invalid_strands() {
        let truncated = write_gtf(&["chr1\tX\texon\t1\t100"]);
        assert!(Annotation::from_gtf(truncated.path()).is_err());

        let unstranded = write_gtf(&[
            "chr1\tX\texon\t1\t100\t.\t.\t.\tgene_id \"G1\"; transcript_id \"T1\";",
        ]);
        assert!(Annotation::from_gtf(unstranded.path()).is_err());
    }

    #[test]
    fn rejects_reused_transcript_ids_with_conflicting_context() {
        let conflicting_gene = write_gtf(&[
            A,
            "chr1\tX\texon\t301\t400\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T1\"; gene_name \"Beta\";",
        ]);
        let error = Annotation::from_gtf(conflicting_gene.path())
            .err()
            .unwrap()
            .to_string();
        assert!(error.contains("transcript_id") && error.contains("reused"));

        let conflicting_strand = write_gtf(&[
            A,
            "chr1\tX\texon\t301\t400\t.\t-\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Alpha\";",
        ]);
        assert!(Annotation::from_gtf(conflicting_strand.path()).is_err());
    }

    #[test]
    fn rejects_non_utf8_annotation_identifiers() {
        for invalid_field in ["chromosome", "gene", "transcript", "exon", "name"] {
            let mut line = b"chr1\tX\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_id \"E1\"; gene_name \"Alpha\";\n".to_vec();
            let needle: &[u8] = match invalid_field {
                "chromosome" => b"chr1",
                "gene" => b"G1",
                "transcript" => b"T1",
                "exon" => b"E1",
                "name" => b"Alpha",
                _ => unreachable!(),
            };
            let offset = line
                .windows(needle.len())
                .position(|value| value == needle)
                .unwrap();
            line[offset] = 0xff;
            let mut file = temppath::NamedTemp::new();
            file.file.write_all(&line).unwrap();
            let path = file.into_path();
            let error = Annotation::from_gtf(path.path()).err().unwrap().to_string();
            assert!(error.contains("UTF-8"), "{invalid_field}: {error}");
        }
    }

    fn compiled_path(label: &str) -> temppath::TempPath {
        let path = std::env::temp_dir().join(format!(
            "anno-test-{}-{:x}-{label}.aic",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        temppath::TempPath(path)
    }

    #[test]
    fn compiled_annotation_is_deterministic_and_semantically_identical() {
        let gtf = write_gtf(&[A, B, C]);
        let annotation = Annotation::from_gtf(gtf.path()).unwrap();
        let first = compiled_path("first");
        let second = compiled_path("second");
        annotation.write_compiled(first.path()).unwrap();
        annotation.write_compiled(second.path()).unwrap();
        assert_eq!(std::fs::read(first.path()).unwrap(), std::fs::read(second.path()).unwrap());

        let restored = Annotation::from_path(first.path()).unwrap();
        assert_eq!(restored.gene_ids, annotation.gene_ids);
        assert_eq!(restored.gene_names, annotation.gene_names);
        assert_eq!(restored.transcript_ids, annotation.transcript_ids);
        assert_eq!(restored.source_exons, annotation.source_exons);
        assert_eq!(restored.chrom_ids, annotation.chrom_ids);
        assert_eq!(restored.transcripts.len(), annotation.transcripts.len());
        for (actual, expected) in restored.transcripts.iter().zip(&annotation.transcripts) {
            assert_eq!((actual.gene, actual.chrom, actual.strand_rev),
                (expected.gene, expected.chrom, expected.strand_rev));
            assert_eq!(actual.exons, expected.exons);
        }
        let chr1 = restored.chrom_ids["chr1"];
        assert_eq!(
            restored.overlapping(chr1, 180, 320),
            annotation.overlapping(chr1, 180, 320)
        );
    }

    #[test]
    fn legacy_v1_compiled_annotation_loads_with_explicitly_unavailable_identifiers() {
        let gtf = write_gtf(&[A, B, C]);
        let annotation = Annotation::from_gtf(gtf.path()).unwrap();
        let legacy = compiled_path("legacy-v1");
        let payload = annotation.compiled_payload(false).unwrap();
        let mut writer = std::io::BufWriter::new(std::fs::File::create(legacy.path()).unwrap());
        writer.write_all(&COMPILED_MAGIC).unwrap();
        writer
            .write_all(&COMPILED_VERSION_WITHOUT_IDENTIFIERS.to_le_bytes())
            .unwrap();
        writer
            .write_all(&(payload.len() as u64).to_le_bytes())
            .unwrap();
        writer.write_all(blake3::hash(&payload).as_bytes()).unwrap();
        writer.write_all(&payload).unwrap();
        writer.flush().unwrap();

        let restored = Annotation::from_compiled(legacy.path()).unwrap();
        assert_eq!(restored.gene_ids, annotation.gene_ids);
        assert_eq!(
            restored.transcript_ids,
            vec![None; restored.transcripts.len()]
        );
        assert!(restored.source_exons.iter().all(Vec::is_empty));

        let resolver = intent::IntentResolver::from_annotation(
            &restored,
            intent::AnnotationIdentity::new("GRCh38", "legacy test").unwrap(),
        )
        .unwrap();
        assert_eq!(resolver.resolve_str("gene:G1").unwrap().stable_id, "G1");
        assert!(matches!(
            resolver.resolve_str("transcript:T1"),
            Err(intent::ResolutionError::IdentifierMetadataUnavailable { .. })
        ));
    }

    #[test]
    fn compiled_annotation_rejects_bad_version_length_checksum_and_structure() {
        let gtf = write_gtf(&[A, B, C]);
        let annotation = Annotation::from_gtf(gtf.path()).unwrap();
        let valid = compiled_path("valid");
        annotation.write_compiled(valid.path()).unwrap();
        let bytes = std::fs::read(valid.path()).unwrap();

        let future = compiled_path("future");
        let mut changed = bytes.clone();
        changed[8..12].copy_from_slice(&(COMPILED_VERSION + 1).to_le_bytes());
        std::fs::write(future.path(), changed).unwrap();
        assert!(Annotation::from_compiled(future.path()).err().unwrap().to_string()
            .contains("unsupported compiled annotation version"));

        let truncated = compiled_path("truncated");
        std::fs::write(truncated.path(), &bytes[..bytes.len() - 1]).unwrap();
        assert!(Annotation::from_compiled(truncated.path()).is_err());

        let trailing = compiled_path("trailing");
        let mut changed = bytes.clone();
        changed.push(0);
        std::fs::write(trailing.path(), changed).unwrap();
        assert!(Annotation::from_compiled(trailing.path()).is_err());

        let corrupt = compiled_path("corrupt");
        let mut changed = bytes.clone();
        *changed.last_mut().unwrap() ^= 1;
        std::fs::write(corrupt.path(), changed).unwrap();
        assert!(Annotation::from_compiled(corrupt.path()).err().unwrap().to_string()
            .contains("checksum mismatch"));

        let hostile_count = compiled_path("hostile-count");
        let mut changed = bytes.clone();
        changed[COMPILED_HEADER_LEN..COMPILED_HEADER_LEN + 4]
            .copy_from_slice(&u32::MAX.to_le_bytes());
        let checksum = blake3::hash(&changed[COMPILED_HEADER_LEN..]);
        changed[20..52].copy_from_slice(checksum.as_bytes());
        std::fs::write(hostile_count.path(), changed).unwrap();
        assert!(Annotation::from_compiled(hostile_count.path())
            .err()
            .unwrap()
            .to_string()
            .contains("cannot fit"));

        let invalid = compiled_path("invalid");
        let mut structurally_bad = Annotation::from_gtf(gtf.path()).unwrap();
        structurally_bad.transcripts[0].gene = u32::MAX;
        structurally_bad.write_compiled(invalid.path()).unwrap();
        assert!(Annotation::from_compiled(invalid.path()).err().unwrap().to_string()
            .contains("invalid dictionary id"));
    }
}
