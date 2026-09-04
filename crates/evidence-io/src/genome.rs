//! Genome identity: a signature of the reference currently bound to an archive.
//!
//! The archive stores coordinates, never sequence; any consumer that *does* consult sequence
//! (internal-priming filters and extension decisions) must use the same bound reference, or its
//! answers are silently wrong. The signature makes that checkable: per-contig BLAKE3 over the
//! uppercased bases (whitespace ignored, so it is invariant to line wrapping and gzip framing),
//! plus one combined digest over the sorted (name, length, hash) triples. It lives in the archive's
//! `meta` section under `"genome_sig"`; readers that predate it ignore the key. A signature proves
//! sequence identity, not that the reference generated the original alignment; logical-v2
//! archives record that relationship separately as an explicit caller declaration.
//!
//! Validation is per-contig on purpose: a windowed query loads one chromosome and can verify just
//! that contig's hash instead of re-hashing 3 GB, and a dev archive restricted to one chromosome
//! validates against a whole-genome signature (archive contigs ⊆ signature contigs).

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub const GENOME_SIG_ALGO: &str = "aie-genome-blake3-v1";

fn finalize_reset(h: &mut blake3::Hasher) -> String {
    let out = h.finalize().to_hex().to_string();
    h.reset();
    out
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct ContigSig {
    pub name: String,
    pub len: u64,
    pub blake3: String,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GenomeSig {
    pub algo: String,
    /// BLAKE3 over "name\tlen\tblake3\n" lines, contigs sorted by name.
    pub digest: String,
    pub contigs: Vec<ContigSig>,
}

impl GenomeSig {
    pub fn contig(&self, name: &str) -> Option<&ContigSig> {
        self.contigs.iter().find(|c| c.name == name)
    }

    pub fn combined_digest(contigs: &[ContigSig]) -> String {
        let mut sorted: Vec<&ContigSig> = contigs.iter().collect();
        sorted.sort_by(|a, b| a.name.cmp(&b.name));
        let mut h = blake3::Hasher::new();
        for c in sorted {
            h.update(format!("{}\t{}\t{}\n", c.name, c.len, c.blake3).as_bytes());
        }
        h.finalize().to_hex().to_string()
    }

    /// Validate a stored normalized signature before it participates in a provenance binding.
    pub fn validate(&self) -> Result<()> {
        if self.algo != GENOME_SIG_ALGO || self.contigs.is_empty() {
            bail!("invalid genome signature algorithm or empty contig set");
        }
        let mut names = std::collections::BTreeSet::new();
        for contig in &self.contigs {
            if contig.name.is_empty()
                || !names.insert(&contig.name)
                || contig.blake3.len() != 64
                || !contig.blake3.bytes().all(|byte| byte.is_ascii_hexdigit())
            {
                bail!("genome signature contains an invalid contig record");
            }
        }
        if self.digest != Self::combined_digest(&self.contigs) {
            bail!("genome signature combined digest is inconsistent with its contigs");
        }
        Ok(())
    }
}

/// A FASTA line source, transparently gunzipped when the file starts with the gzip magic.
fn open_fasta(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    open_fasta_file(file)
}

fn open_fasta_file(f: File) -> Result<Box<dyn BufRead>> {
    let mut head = [0u8; 2];
    use std::io::Seek;
    let mut probe = f.try_clone()?;
    let n = probe.read(&mut head)?;
    probe.seek(std::io::SeekFrom::Start(0))?;
    drop(probe);
    if n == 2 && head == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::with_capacity(1 << 20, flate2::bufread::MultiGzDecoder::new(BufReader::with_capacity(1 << 20, f)))))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, f)))
    }
}

fn contig_name(header: &str) -> String {
    header[1..].split_whitespace().next().unwrap_or("").to_string()
}

/// Hash every contig of a FASTA (plain or gzipped) without keeping sequence in memory.
pub fn sig_from_fasta(path: &Path) -> Result<GenomeSig> {
    sig_from_reader(open_fasta(path)?, &path.display().to_string())
}

/// Hash every contig from an already-open FASTA. Callers that also report byte identity can open
/// and bind all handles before work starts, avoiding a second pathname lookup race.
pub fn sig_from_fasta_file(file: File) -> Result<GenomeSig> {
    sig_from_reader(open_fasta_file(file)?, "open FASTA input")
}

fn sig_from_reader(mut r: Box<dyn BufRead>, source: &str) -> Result<GenomeSig> {
    let mut contigs: Vec<ContigSig> = Vec::new();
    let (mut cur, mut len, mut h) = (None::<String>, 0u64, blake3::Hasher::new());
    let mut line = String::new();
    loop {
        line.clear();
        if r.read_line(&mut line)? == 0 {
            break;
        }
        if line.starts_with('>') {
            if let Some(name) = cur.take() {
                contigs.push(ContigSig { name, len, blake3: finalize_reset(&mut h) });
            }
            cur = Some(contig_name(line.trim_end()));
            len = 0;
        } else if cur.is_some() {
            let bases = line.trim_end().as_bytes();
            let up: Vec<u8> = bases.iter().map(|b| b.to_ascii_uppercase()).collect();
            len += up.len() as u64;
            h.update(&up);
        }
    }
    if let Some(name) = cur.take() {
        contigs.push(ContigSig { name, len, blake3: finalize_reset(&mut h) });
    }
    if contigs.is_empty() {
        bail!("{source}: no FASTA records");
    }
    let digest = GenomeSig::combined_digest(&contigs);
    Ok(GenomeSig { algo: GENOME_SIG_ALGO.into(), digest, contigs })
}

/// Load one contig's uppercased sequence, verifying it against `sig` when one is present.
/// Streams the FASTA up to the contig; for gzipped genomes this costs a scan, so callers doing
/// many contigs should iterate in file order via `for_each_contig`.
pub fn load_contig(path: &Path, name: &str, sig: Option<&GenomeSig>) -> Result<Vec<u8>> {
    let mut found = None;
    for_each_contig(path, |cname, seq| {
        if cname == name {
            found = Some(seq.to_vec());
            false // stop
        } else {
            true
        }
    })?;
    let seq = found.with_context(|| format!("contig {name} not found in {}", path.display()))?;
    if let Some(sig) = sig {
        verify_contig(sig, name, &seq)?;
    }
    Ok(seq)
}

/// Verify a loaded, uppercased contig sequence against the archive's stored signature.
pub fn verify_contig(sig: &GenomeSig, name: &str, seq: &[u8]) -> Result<()> {
    let Some(c) = sig.contig(name) else {
        bail!(
            "genome mismatch: contig {name} is not in the archive's genome signature \
             (signature algo {}, {} contigs)",
            sig.algo,
            sig.contigs.len()
        );
    };
    if c.len != seq.len() as u64 {
        bail!(
            "genome mismatch: contig {name} is {} bp in the supplied FASTA but {} bp in the \
             genome the archive was built from",
            seq.len(),
            c.len
        );
    }
    let got = blake3::hash(seq).to_hex().to_string();
    if got != c.blake3 {
        bail!(
            "genome mismatch: contig {name} sequence differs from the genome the archive was \
             built from (blake3 {got} vs stamped {})",
            c.blake3
        );
    }
    Ok(())
}

/// Stream a FASTA, invoking `f(name, uppercased_sequence)` per contig in file order.
/// `f` returns false to stop early.
pub fn for_each_contig(path: &Path, mut f: impl FnMut(&str, &[u8]) -> bool) -> Result<()> {
    let mut r = open_fasta(path)?;
    let mut cur: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    let mut line = String::new();
    loop {
        line.clear();
        let eof = r.read_line(&mut line)? == 0;
        if eof || line.starts_with('>') {
            if let Some(name) = cur.take() {
                if !f(&name, &seq) {
                    return Ok(());
                }
            }
            if eof {
                return Ok(());
            }
            cur = Some(contig_name(line.trim_end()));
            seq.clear();
        } else if cur.is_some() {
            seq.extend(line.trim_end().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_fasta(tag: &str, content: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("aie-genome-{tag}-{}.fa", std::process::id()));
        std::fs::File::create(&p).unwrap().write_all(content.as_bytes()).unwrap();
        p
    }

    #[test]
    fn sig_is_invariant_to_case_and_wrapping() {
        let a = temp_fasta("wrap", ">chrA desc\nACGT\nACGT\n>chrB\nnnnACG\n");
        let sig1 = sig_from_fasta(&a).unwrap();
        std::fs::write(&a, ">chrA other-desc\nacgtacgt\n>chrB\nNNNacg\n").unwrap();
        let sig2 = sig_from_fasta(&a).unwrap();
        assert_eq!(sig1.digest, sig2.digest);
        assert_eq!(sig1.contigs.len(), 2);
        assert_eq!(sig1.contig("chrA").unwrap().len, 8);
        sig1.validate().unwrap();
        let seq = load_contig(&a, "chrB", Some(&sig1)).unwrap();
        assert_eq!(seq, b"NNNACG");
        std::fs::remove_file(&a).ok();
    }

    #[test]
    fn mismatch_is_detected() {
        let a = temp_fasta("mismatch", ">chrA\nACGTACGT\n");
        let sig = sig_from_fasta(&a).unwrap();
        assert!(verify_contig(&sig, "chrA", b"ACGTACGA").is_err());
        assert!(verify_contig(&sig, "chrA", b"ACGT").is_err());
        assert!(verify_contig(&sig, "chrZ", b"ACGTACGT").is_err());
        assert!(verify_contig(&sig, "chrA", b"ACGTACGT").is_ok());
        std::fs::remove_file(&a).ok();
    }

    #[test]
    fn malformed_stored_signature_is_rejected() {
        let a = temp_fasta("invalid", ">chrA\nACGTACGT\n");
        let mut sig = sig_from_fasta(&a).unwrap();
        let replacement = if sig.digest.starts_with('0') {
            "1"
        } else {
            "0"
        };
        sig.digest.replace_range(..1, replacement);
        assert!(sig.validate().is_err());
        let mut duplicate = sig_from_fasta(&a).unwrap();
        duplicate.contigs.push(duplicate.contigs[0].clone());
        duplicate.digest = GenomeSig::combined_digest(&duplicate.contigs);
        assert!(duplicate.validate().is_err());
        std::fs::remove_file(&a).ok();
    }
}
