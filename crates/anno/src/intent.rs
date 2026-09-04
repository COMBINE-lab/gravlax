//! Resolve biological identifiers to explicit, versioned genomic intent.
//!
//! Coordinate-only archive queries are reproducible but awkward for people. This module is the
//! boundary between names such as `TP53`, `ENSG...`, `ENST...`, or `ENSE...` and those archive
//! coordinates. Resolution always carries caller-supplied assembly and annotation identity;
//! identifiers are never silently interpreted against an unnamed reference.

use crate::Annotation;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnnotationIdentity {
    /// Reference assembly, including patch/haplotype qualification when relevant (for example
    /// `GRCh38.p14`).
    pub assembly: String,
    /// Annotation release or other immutable label (for example `GENCODE 49`).
    pub annotation: String,
    /// Optional expected content digest supplied by the caller or project manifest. A resolver
    /// constructed from a path always fills this with the observed digest of the parsed file.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub digest: Option<String>,
}

impl AnnotationIdentity {
    pub fn new(
        assembly: impl Into<String>,
        annotation: impl Into<String>,
    ) -> Result<Self, ResolverBuildError> {
        let identity = Self {
            assembly: assembly.into(),
            annotation: annotation.into(),
            digest: None,
        };
        identity.validate()?;
        Ok(identity)
    }

    pub fn with_digest(mut self, digest: impl Into<String>) -> Result<Self, ResolverBuildError> {
        self.digest = Some(digest.into());
        self.validate()?;
        Ok(self)
    }

    fn validate(&self) -> Result<(), ResolverBuildError> {
        if self.assembly.trim().is_empty() {
            return Err(ResolverBuildError::InvalidIdentity(
                "assembly must not be empty".into(),
            ));
        }
        if self.annotation.trim().is_empty() {
            return Err(ResolverBuildError::InvalidIdentity(
                "annotation release must not be empty".into(),
            ));
        }
        if self
            .digest
            .as_ref()
            .is_some_and(|value| value.trim().is_empty())
        {
            return Err(ResolverBuildError::InvalidIdentity(
                "annotation digest must not be empty when supplied".into(),
            ));
        }
        if let Some(digest) = &self.digest {
            let Some(encoded) = digest.strip_prefix("blake3:") else {
                return Err(ResolverBuildError::InvalidIdentity(
                    "annotation digest must use the blake3:<64 lowercase hex characters> form"
                        .into(),
                ));
            };
            if encoded.len() != 64
                || !encoded
                    .bytes()
                    .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
            {
                return Err(ResolverBuildError::InvalidIdentity(
                    "annotation digest must use the blake3:<64 lowercase hex characters> form"
                        .into(),
                ));
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum QueryKind {
    Auto,
    Gene,
    Transcript,
    Exon,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct IdentifierQuery {
    pub kind: QueryKind,
    pub value: String,
}

#[expect(
    clippy::result_large_err,
    reason = "resolution errors retain the complete query and annotation identity for callers"
)]
impl IdentifierQuery {
    pub fn new(kind: QueryKind, value: impl Into<String>) -> Result<Self, ResolutionError> {
        let value = value.into();
        let value = value.trim();
        if value.is_empty() {
            return Err(ResolutionError::EmptyIdentifier);
        }
        Ok(Self {
            kind,
            value: value.to_owned(),
        })
    }

    /// Parse an identifier with an optional `gene:`, `transcript:`, or `exon:` prefix. An
    /// unprefixed value uses ambiguity-detecting automatic resolution.
    pub fn parse(value: &str) -> Result<Self, ResolutionError> {
        for (prefix, kind) in [
            ("gene:", QueryKind::Gene),
            ("transcript:", QueryKind::Transcript),
            ("exon:", QueryKind::Exon),
        ] {
            if let Some(identifier) = value.strip_prefix(prefix) {
                return Self::new(kind, identifier);
            }
        }
        Self::new(QueryKind::Auto, value)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FeatureKind {
    Gene,
    Transcript,
    Exon,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MatchBasis {
    StableId,
    StableIdWithoutVersion,
    GeneSymbol,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct GenomicLocus {
    pub contig: String,
    /// Zero-based, inclusive start.
    pub start: u32,
    /// Zero-based, exclusive end.
    pub end: u32,
    pub strand: Strand,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ResolvedFeature {
    pub kind: FeatureKind,
    pub requested: String,
    pub stable_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_name: Option<String>,
    pub matched_by: MatchBasis,
    pub gene_ids: Vec<String>,
    pub transcript_ids: Vec<String>,
    pub loci: Vec<GenomicLocus>,
    pub identity: AnnotationIdentity,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ResolutionCandidate {
    pub kind: FeatureKind,
    pub stable_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub display_name: Option<String>,
    pub loci: Vec<GenomicLocus>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ResolverBuildError {
    InvalidIdentity(String),
    InvalidAnnotation(String),
    Load { path: String, message: String },
}

impl fmt::Display for ResolverBuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidIdentity(message) => write!(f, "invalid annotation identity: {message}"),
            Self::InvalidAnnotation(message) => write!(f, "invalid annotation model: {message}"),
            Self::Load { path, message } => write!(f, "loading annotation {path}: {message}"),
        }
    }
}

impl std::error::Error for ResolverBuildError {}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ResolutionError {
    EmptyIdentifier,
    NotFound {
        query: IdentifierQuery,
        identity: AnnotationIdentity,
    },
    Ambiguous {
        query: IdentifierQuery,
        candidates: Vec<ResolutionCandidate>,
    },
    IdentifierMetadataUnavailable {
        query: IdentifierQuery,
        identity: AnnotationIdentity,
        missing: Vec<FeatureKind>,
    },
}

impl fmt::Display for ResolutionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyIdentifier => write!(f, "identifier must not be empty"),
            Self::NotFound { query, identity } => write!(
                f,
                "{} was not found in {} on {}",
                query.value, identity.annotation, identity.assembly
            ),
            Self::Ambiguous { query, candidates } => write!(
                f,
                "{} is ambiguous in the selected annotation ({} candidates); use a typed prefix or stable ID",
                query.value,
                candidates.len()
            ),
            Self::IdentifierMetadataUnavailable { query, identity, missing } => write!(
                f,
                "cannot resolve {:?} identifier {} in {} on {}: {:?} identifier metadata is unavailable (legacy AIC v1 or source GTF omitted it)",
                query.kind, query.value, identity.annotation, identity.assembly, missing
            ),
        }
    }
}

impl std::error::Error for ResolutionError {}

#[derive(Clone, Debug)]
struct GeneFeature {
    stable_id: String,
    name: String,
    loci: Vec<GenomicLocus>,
}

#[derive(Clone, Debug)]
struct TranscriptFeature {
    stable_id: String,
    gene_id: String,
    gene_name: String,
    locus: GenomicLocus,
}

#[derive(Clone, Debug)]
struct ExonFeature {
    stable_id: String,
    gene_ids: Vec<String>,
    transcript_ids: Vec<String>,
    loci: Vec<GenomicLocus>,
}

#[derive(Clone, Debug, Default)]
struct IdIndex {
    exact: HashMap<String, Vec<usize>>,
    without_version: HashMap<String, Vec<usize>>,
}

impl IdIndex {
    fn insert(&mut self, identifier: &str, index: usize) {
        self.exact
            .entry(identifier.to_owned())
            .or_default()
            .push(index);
        self.without_version
            .entry(without_numeric_version(identifier).to_owned())
            .or_default()
            .push(index);
    }

    fn matches(&self, identifier: &str) -> (Vec<usize>, MatchBasis) {
        if let Some(matches) = self.exact.get(identifier) {
            return (deduplicated(matches), MatchBasis::StableId);
        }
        // A caller who supplied a version asked for that exact biological version. Only an
        // unversioned query may expand to the version carried by the selected annotation.
        if without_numeric_version(identifier) != identifier {
            return (Vec::new(), MatchBasis::StableIdWithoutVersion);
        }
        let base = without_numeric_version(identifier);
        (
            self.without_version
                .get(base)
                .map_or_else(Vec::new, |matches| deduplicated(matches)),
            MatchBasis::StableIdWithoutVersion,
        )
    }
}

fn deduplicated(values: &[usize]) -> Vec<usize> {
    let mut out = values.to_vec();
    out.sort_unstable();
    out.dedup();
    out
}

fn without_numeric_version(identifier: &str) -> &str {
    let Some((base, version)) = identifier.rsplit_once('.') else {
        return identifier;
    };
    if !base.is_empty() && !version.is_empty() && version.bytes().all(|byte| byte.is_ascii_digit())
    {
        base
    } else {
        identifier
    }
}

#[cfg(unix)]
fn same_open_file_snapshot(before: &std::fs::Metadata, after: &std::fs::Metadata) -> bool {
    use std::os::unix::fs::MetadataExt;

    before.dev() == after.dev()
        && before.ino() == after.ino()
        && before.len() == after.len()
        && before.mtime() == after.mtime()
        && before.mtime_nsec() == after.mtime_nsec()
        && before.ctime() == after.ctime()
        && before.ctime_nsec() == after.ctime_nsec()
}

#[cfg(not(unix))]
fn same_open_file_snapshot(before: &std::fs::Metadata, after: &std::fs::Metadata) -> bool {
    before.len() == after.len() && before.modified().ok() == after.modified().ok()
}

/// An annotation model bound to the digest and scientific identity of the exact file snapshot
/// from which it was parsed.
///
/// Commands that both interpret evidence and report annotation provenance should use this loader
/// instead of hashing and reopening a pathname independently.  Digest verification, parsing, and
/// the before/after metadata guard all operate on one open file description.
pub struct BoundAnnotation {
    annotation: Annotation,
    identity: AnnotationIdentity,
}

impl BoundAnnotation {
    pub fn from_path(
        path: &Path,
        identity: AnnotationIdentity,
    ) -> Result<Self, ResolverBuildError> {
        identity.validate()?;
        let file = std::fs::File::open(path).map_err(|error| ResolverBuildError::Load {
            path: path.display().to_string(),
            message: format!("opening annotation: {error}"),
        })?;
        Self::from_open_file(path, file, identity)
    }

    fn from_open_file(
        path: &Path,
        file: std::fs::File,
        mut identity: AnnotationIdentity,
    ) -> Result<Self, ResolverBuildError> {
        identity.validate()?;
        let snapshot = file.try_clone().map_err(|error| ResolverBuildError::Load {
            path: path.display().to_string(),
            message: format!("capturing annotation file descriptor: {error}"),
        })?;
        let before = snapshot
            .metadata()
            .map_err(|error| ResolverBuildError::Load {
                path: path.display().to_string(),
                message: format!("reading annotation metadata: {error}"),
            })?;
        let (annotation, actual) =
            Annotation::from_open_file_with_digest(file, path).map_err(|error| {
                ResolverBuildError::Load {
                    path: path.display().to_string(),
                    message: error.to_string(),
                }
            })?;
        if let Some(expected) = &identity.digest {
            if actual != *expected {
                return Err(ResolverBuildError::InvalidIdentity(format!(
                    "annotation digest does not match {} (expected {expected}, observed {actual})",
                    path.display()
                )));
            }
        }
        identity.digest = Some(actual);
        let after = snapshot.metadata().map_err(|error| ResolverBuildError::Load {
            path: path.display().to_string(),
            message: format!("re-reading annotation metadata: {error}"),
        })?;
        if !same_open_file_snapshot(&before, &after) {
            return Err(ResolverBuildError::Load {
                path: path.display().to_string(),
                message: "annotation changed while its digest and model were being loaded".into(),
            });
        }
        Ok(Self {
            annotation,
            identity,
        })
    }

    pub fn annotation(&self) -> &Annotation {
        &self.annotation
    }

    pub fn identity(&self) -> &AnnotationIdentity {
        &self.identity
    }

    pub fn into_parts(self) -> (Annotation, AnnotationIdentity) {
        (self.annotation, self.identity)
    }
}

/// Immutable lookup index suitable for reuse across a batch or interactive session.
pub struct IntentResolver {
    identity: AnnotationIdentity,
    genes: Vec<GeneFeature>,
    transcripts: Vec<TranscriptFeature>,
    exons: Vec<ExonFeature>,
    gene_ids: IdIndex,
    gene_symbols: HashMap<String, Vec<usize>>,
    transcript_ids: IdIndex,
    exon_ids: IdIndex,
    transcript_metadata_available: bool,
    exon_metadata_available: bool,
}

#[expect(
    clippy::result_large_err,
    reason = "resolution errors retain structured candidates and annotation identity for callers"
)]
impl IntentResolver {
    pub fn from_path(
        path: &Path,
        identity: AnnotationIdentity,
    ) -> Result<Self, ResolverBuildError> {
        let bound = BoundAnnotation::from_path(path, identity)?;
        Self::from_validated_annotation(bound.annotation(), bound.identity().clone())
    }

    /// Build the identifier index over an annotation snapshot that was already loaded and bound
    /// to its observed content digest. This avoids reparsing the annotation when a command needs
    /// both the assignment model and biological-identifier resolution.
    pub fn from_bound_annotation(bound: &BoundAnnotation) -> Result<Self, ResolverBuildError> {
        Self::from_validated_annotation(bound.annotation(), bound.identity().clone())
    }

    #[cfg(test)]
    fn from_open_file(
        path: &Path,
        file: std::fs::File,
        identity: AnnotationIdentity,
    ) -> Result<Self, ResolverBuildError> {
        let bound = BoundAnnotation::from_open_file(path, file, identity)?;
        Self::from_validated_annotation(bound.annotation(), bound.identity().clone())
    }

    pub fn from_annotation(
        annotation: &Annotation,
        identity: AnnotationIdentity,
    ) -> Result<Self, ResolverBuildError> {
        identity.validate()?;
        if identity.digest.is_some() {
            return Err(ResolverBuildError::InvalidIdentity(
                "an annotation digest can only be bound when constructing from_path".into(),
            ));
        }
        Self::from_validated_annotation(annotation, identity)
    }

    fn from_validated_annotation(
        annotation: &Annotation,
        identity: AnnotationIdentity,
    ) -> Result<Self, ResolverBuildError> {
        if annotation.gene_ids.len() != annotation.gene_names.len() {
            return Err(ResolverBuildError::InvalidAnnotation(
                "gene id/name dictionaries have different lengths".into(),
            ));
        }
        if annotation.transcript_ids.len() != annotation.transcripts.len()
            || annotation.source_exons.len() != annotation.transcripts.len()
        {
            return Err(ResolverBuildError::InvalidAnnotation(
                "identifier metadata is not parallel to the transcript table".into(),
            ));
        }

        let mut contigs = vec![None; annotation.chrom_ids.len()];
        for (name, &id) in &annotation.chrom_ids {
            let Some(slot) = contigs.get_mut(id as usize) else {
                return Err(ResolverBuildError::InvalidAnnotation(format!(
                    "contig id {id} is outside its dictionary"
                )));
            };
            if slot.replace(name.clone()).is_some() {
                return Err(ResolverBuildError::InvalidAnnotation(format!(
                    "duplicate contig id {id}"
                )));
            }
        }
        if let Some(index) = contigs.iter().position(Option::is_none) {
            return Err(ResolverBuildError::InvalidAnnotation(format!(
                "contig id {index} is absent from its dictionary"
            )));
        }
        let contigs: Vec<String> = contigs.into_iter().map(Option::unwrap).collect();

        let mut gene_loci: Vec<BTreeMap<(u32, bool), (u32, u32)>> =
            vec![BTreeMap::new(); annotation.gene_ids.len()];
        let mut transcripts = Vec::new();
        let mut exon_aggregates: BTreeMap<String, ExonFeature> = BTreeMap::new();
        let transcript_metadata_available = annotation.transcript_ids.iter().all(Option::is_some);
        let exon_metadata_available = annotation.transcripts.is_empty()
            || annotation.source_exons.iter().all(|source_exons| {
                !source_exons.is_empty()
                    && source_exons
                        .iter()
                        .all(|source_exon| source_exon.id.is_some())
            });

        for (tx_index, transcript) in annotation.transcripts.iter().enumerate() {
            let gene = transcript.gene as usize;
            let chrom = transcript.chrom as usize;
            if gene >= annotation.gene_ids.len() || chrom >= contigs.len() {
                return Err(ResolverBuildError::InvalidAnnotation(format!(
                    "transcript {tx_index} references an invalid gene or contig"
                )));
            }
            let (start, end) = transcript.span();
            let locus = GenomicLocus {
                contig: contigs[chrom].clone(),
                start,
                end,
                strand: if transcript.strand_rev {
                    Strand::Reverse
                } else {
                    Strand::Forward
                },
            };
            gene_loci[gene]
                .entry((transcript.chrom, transcript.strand_rev))
                .and_modify(|span| {
                    span.0 = span.0.min(start);
                    span.1 = span.1.max(end);
                })
                .or_insert((start, end));

            let transcript_id = annotation.transcript_ids[tx_index].as_deref();
            if let Some(transcript_id) = transcript_id {
                transcripts.push(TranscriptFeature {
                    stable_id: transcript_id.to_owned(),
                    gene_id: annotation.gene_ids[gene].clone(),
                    gene_name: annotation.gene_names[gene].clone(),
                    locus: locus.clone(),
                });
            }

            for source_exon in &annotation.source_exons[tx_index] {
                if source_exon.interval.start >= source_exon.interval.end
                    || !transcript.exons.iter().any(|assignment| {
                        assignment.start <= source_exon.interval.start
                            && source_exon.interval.end <= assignment.end
                    })
                {
                    return Err(ResolverBuildError::InvalidAnnotation(format!(
                        "transcript {tx_index} has a source exon outside its assignment exons"
                    )));
                }
                if let Some(exon_id) = &source_exon.id {
                    let exon_locus = GenomicLocus {
                        contig: contigs[chrom].clone(),
                        start: source_exon.interval.start,
                        end: source_exon.interval.end,
                        strand: locus.strand,
                    };
                    let aggregate =
                        exon_aggregates
                            .entry(exon_id.clone())
                            .or_insert_with(|| ExonFeature {
                                stable_id: exon_id.clone(),
                                gene_ids: Vec::new(),
                                transcript_ids: Vec::new(),
                                loci: Vec::new(),
                            });
                    if !aggregate.gene_ids.contains(&annotation.gene_ids[gene]) {
                        aggregate.gene_ids.push(annotation.gene_ids[gene].clone());
                    }
                    if let Some(transcript_id) = transcript_id {
                        if !aggregate
                            .transcript_ids
                            .iter()
                            .any(|id| id == transcript_id)
                        {
                            aggregate.transcript_ids.push(transcript_id.to_owned());
                        }
                    }
                    if !aggregate.loci.contains(&exon_locus) {
                        aggregate.loci.push(exon_locus);
                    }
                }
            }
        }

        let genes: Vec<GeneFeature> = annotation
            .gene_ids
            .iter()
            .zip(&annotation.gene_names)
            .enumerate()
            .map(|(gene, (stable_id, name))| {
                let mut loci: Vec<GenomicLocus> = gene_loci[gene]
                    .iter()
                    .map(|(&(chrom, strand_rev), &(start, end))| GenomicLocus {
                        contig: contigs[chrom as usize].clone(),
                        start,
                        end,
                        strand: if strand_rev {
                            Strand::Reverse
                        } else {
                            Strand::Forward
                        },
                    })
                    .collect();
                loci.sort();
                GeneFeature {
                    stable_id: stable_id.clone(),
                    name: name.clone(),
                    loci,
                }
            })
            .collect();
        let mut exons: Vec<ExonFeature> = exon_aggregates.into_values().collect();
        for exon in &mut exons {
            exon.gene_ids.sort();
            exon.transcript_ids.sort();
            exon.loci.sort();
        }

        let mut gene_ids = IdIndex::default();
        let mut gene_symbols: HashMap<String, Vec<usize>> = HashMap::new();
        for (index, gene) in genes.iter().enumerate() {
            gene_ids.insert(&gene.stable_id, index);
            gene_symbols
                .entry(gene.name.clone())
                .or_default()
                .push(index);
        }
        let mut transcript_ids = IdIndex::default();
        for (index, transcript) in transcripts.iter().enumerate() {
            transcript_ids.insert(&transcript.stable_id, index);
        }
        let mut exon_ids = IdIndex::default();
        for (index, exon) in exons.iter().enumerate() {
            exon_ids.insert(&exon.stable_id, index);
        }

        Ok(Self {
            identity,
            genes,
            transcripts,
            exons,
            gene_ids,
            gene_symbols,
            transcript_ids,
            exon_ids,
            transcript_metadata_available,
            exon_metadata_available,
        })
    }

    pub fn identity(&self) -> &AnnotationIdentity {
        &self.identity
    }

    pub fn resolve_str(&self, query: &str) -> Result<ResolvedFeature, ResolutionError> {
        self.resolve(&IdentifierQuery::parse(query)?)
    }

    pub fn resolve(&self, query: &IdentifierQuery) -> Result<ResolvedFeature, ResolutionError> {
        if query.value.trim().is_empty() {
            return Err(ResolutionError::EmptyIdentifier);
        }
        match query.kind {
            QueryKind::Gene => self.resolve_gene(query),
            QueryKind::Transcript => self.resolve_transcript(query),
            QueryKind::Exon => self.resolve_exon(query),
            QueryKind::Auto => self.resolve_auto(query),
        }
    }

    fn resolve_auto(&self, query: &IdentifierQuery) -> Result<ResolvedFeature, ResolutionError> {
        let mut matches = Vec::new();
        if let Some(value) = self.try_gene(query)? {
            matches.push(value);
        }
        if let Some(value) = self.try_transcript(query)? {
            matches.push(value);
        }
        if let Some(value) = self.try_exon(query)? {
            matches.push(value);
        }
        if matches.len() == 1 {
            return Ok(matches.pop().unwrap());
        }
        if matches.len() > 1 {
            return Err(ResolutionError::Ambiguous {
                query: query.clone(),
                candidates: matches.iter().map(candidate_from_resolved).collect(),
            });
        }
        let mut missing = Vec::new();
        if !self.transcript_metadata_available {
            missing.push(FeatureKind::Transcript);
        }
        if !self.exon_metadata_available {
            missing.push(FeatureKind::Exon);
        }
        if !missing.is_empty() {
            Err(self.metadata_unavailable(query, missing))
        } else {
            Err(self.not_found(query))
        }
    }

    fn resolve_gene(&self, query: &IdentifierQuery) -> Result<ResolvedFeature, ResolutionError> {
        self.try_gene(query)?.ok_or_else(|| self.not_found(query))
    }

    fn try_gene(
        &self,
        query: &IdentifierQuery,
    ) -> Result<Option<ResolvedFeature>, ResolutionError> {
        let (stable_matches, stable_basis) = self.gene_ids.matches(&query.value);
        let symbol_matches = self
            .gene_symbols
            .get(&query.value)
            .map_or_else(Vec::new, |matches| deduplicated(matches));
        let mut combined: BTreeMap<usize, MatchBasis> = BTreeMap::new();
        for index in stable_matches {
            combined.insert(index, stable_basis);
        }
        for index in symbol_matches {
            combined.entry(index).or_insert(MatchBasis::GeneSymbol);
        }
        self.select_gene(query, combined)
    }

    fn select_gene(
        &self,
        query: &IdentifierQuery,
        matches: BTreeMap<usize, MatchBasis>,
    ) -> Result<Option<ResolvedFeature>, ResolutionError> {
        if matches.len() > 1 {
            return Err(ResolutionError::Ambiguous {
                query: query.clone(),
                candidates: matches
                    .keys()
                    .map(|&index| self.gene_candidate(index))
                    .collect(),
            });
        }
        let Some((index, basis)) = matches.into_iter().next() else {
            return Ok(None);
        };
        let gene = &self.genes[index];
        Ok(Some(ResolvedFeature {
            kind: FeatureKind::Gene,
            requested: query.value.clone(),
            stable_id: gene.stable_id.clone(),
            display_name: Some(gene.name.clone()),
            matched_by: basis,
            gene_ids: vec![gene.stable_id.clone()],
            transcript_ids: Vec::new(),
            loci: gene.loci.clone(),
            identity: self.identity.clone(),
        }))
    }

    fn resolve_transcript(
        &self,
        query: &IdentifierQuery,
    ) -> Result<ResolvedFeature, ResolutionError> {
        match self.try_transcript(query)? {
            Some(transcript) => Ok(transcript),
            None if !self.transcript_metadata_available => {
                Err(self.metadata_unavailable(query, vec![FeatureKind::Transcript]))
            }
            None => Err(self.not_found(query)),
        }
    }

    fn try_transcript(
        &self,
        query: &IdentifierQuery,
    ) -> Result<Option<ResolvedFeature>, ResolutionError> {
        let (matches, basis) = self.transcript_ids.matches(&query.value);
        if matches.len() > 1 {
            return Err(ResolutionError::Ambiguous {
                query: query.clone(),
                candidates: matches
                    .into_iter()
                    .map(|index| self.transcript_candidate(index))
                    .collect(),
            });
        }
        let Some(index) = matches.into_iter().next() else {
            return Ok(None);
        };
        let transcript = &self.transcripts[index];
        Ok(Some(ResolvedFeature {
            kind: FeatureKind::Transcript,
            requested: query.value.clone(),
            stable_id: transcript.stable_id.clone(),
            display_name: Some(transcript.gene_name.clone()),
            matched_by: basis,
            gene_ids: vec![transcript.gene_id.clone()],
            transcript_ids: vec![transcript.stable_id.clone()],
            loci: vec![transcript.locus.clone()],
            identity: self.identity.clone(),
        }))
    }

    fn resolve_exon(&self, query: &IdentifierQuery) -> Result<ResolvedFeature, ResolutionError> {
        match self.try_exon(query)? {
            Some(exon) => Ok(exon),
            None if !self.exon_metadata_available => {
                Err(self.metadata_unavailable(query, vec![FeatureKind::Exon]))
            }
            None => Err(self.not_found(query)),
        }
    }

    fn try_exon(
        &self,
        query: &IdentifierQuery,
    ) -> Result<Option<ResolvedFeature>, ResolutionError> {
        let (matches, basis) = self.exon_ids.matches(&query.value);
        if matches.len() > 1 {
            return Err(ResolutionError::Ambiguous {
                query: query.clone(),
                candidates: matches
                    .into_iter()
                    .map(|index| self.exon_candidate(index))
                    .collect(),
            });
        }
        let Some(index) = matches.into_iter().next() else {
            return Ok(None);
        };
        let exon = &self.exons[index];
        Ok(Some(ResolvedFeature {
            kind: FeatureKind::Exon,
            requested: query.value.clone(),
            stable_id: exon.stable_id.clone(),
            display_name: None,
            matched_by: basis,
            gene_ids: exon.gene_ids.clone(),
            transcript_ids: exon.transcript_ids.clone(),
            loci: exon.loci.clone(),
            identity: self.identity.clone(),
        }))
    }

    fn gene_candidate(&self, index: usize) -> ResolutionCandidate {
        let gene = &self.genes[index];
        ResolutionCandidate {
            kind: FeatureKind::Gene,
            stable_id: gene.stable_id.clone(),
            display_name: Some(gene.name.clone()),
            loci: gene.loci.clone(),
        }
    }

    fn transcript_candidate(&self, index: usize) -> ResolutionCandidate {
        let transcript = &self.transcripts[index];
        ResolutionCandidate {
            kind: FeatureKind::Transcript,
            stable_id: transcript.stable_id.clone(),
            display_name: Some(transcript.gene_name.clone()),
            loci: vec![transcript.locus.clone()],
        }
    }

    fn exon_candidate(&self, index: usize) -> ResolutionCandidate {
        let exon = &self.exons[index];
        ResolutionCandidate {
            kind: FeatureKind::Exon,
            stable_id: exon.stable_id.clone(),
            display_name: None,
            loci: exon.loci.clone(),
        }
    }

    fn not_found(&self, query: &IdentifierQuery) -> ResolutionError {
        ResolutionError::NotFound {
            query: query.clone(),
            identity: self.identity.clone(),
        }
    }

    fn metadata_unavailable(
        &self,
        query: &IdentifierQuery,
        missing: Vec<FeatureKind>,
    ) -> ResolutionError {
        ResolutionError::IdentifierMetadataUnavailable {
            query: query.clone(),
            identity: self.identity.clone(),
            missing,
        }
    }
}

fn candidate_from_resolved(feature: &ResolvedFeature) -> ResolutionCandidate {
    ResolutionCandidate {
        kind: feature.kind,
        stable_id: feature.stable_id.clone(),
        display_name: feature.display_name.clone(),
        loci: feature.loci.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;
    use std::io::Write;

    struct TempPath(std::path::PathBuf);

    impl TempPath {
        fn new(extension: &str) -> Self {
            let path = std::env::temp_dir().join(format!(
                "gravlax-intent-{}-{:x}.{extension}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            Self(path)
        }
    }

    impl Drop for TempPath {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    fn annotation() -> (TempPath, Annotation) {
        let path = TempPath::new("gtf");
        let mut file = std::fs::File::create(&path.0).unwrap();
        writeln!(file, "chr1\tX\texon\t101\t200\t.\t+\t.\tgene_id \"ENSG1.4\"; transcript_id \"ENST1.2\"; exon_id \"ENSE1.7\"; gene_name \"ALPHA\";").unwrap();
        writeln!(file, "chr1\tX\texon\t301\t400\t.\t+\t.\tgene_id \"ENSG1.4\"; transcript_id \"ENST1.2\"; exon_id \"ENSE2.1\"; gene_name \"ALPHA\";").unwrap();
        writeln!(file, "chr2\tX\texon\t501\t600\t.\t-\t.\tgene_id \"ENSG2.1\"; transcript_id \"ENST2.8\"; exon_id \"ENSE3.1\"; gene_name \"DUP\";").unwrap();
        writeln!(file, "chr3\tX\texon\t701\t800\t.\t+\t.\tgene_id \"ENSG3.1\"; transcript_id \"ENST3.1\"; exon_id \"ENSE4.1\"; gene_name \"DUP\";").unwrap();
        drop(file);
        let annotation = Annotation::from_gtf(&path.0).unwrap();
        (path, annotation)
    }

    fn identity() -> AnnotationIdentity {
        AnnotationIdentity::new("GRCh38.p14", "GENCODE 49").unwrap()
    }

    #[test]
    fn resolves_symbol_and_versionless_stable_identifiers() {
        let (_path, annotation) = annotation();
        let resolver = IntentResolver::from_annotation(&annotation, identity()).unwrap();
        let gene = resolver.resolve_str("ALPHA").unwrap();
        assert_eq!(gene.kind, FeatureKind::Gene);
        assert_eq!(gene.stable_id, "ENSG1.4");
        assert_eq!(gene.matched_by, MatchBasis::GeneSymbol);
        assert_eq!(gene.identity.assembly, "GRCh38.p14");

        let transcript = resolver.resolve_str("transcript:ENST1").unwrap();
        assert_eq!(transcript.stable_id, "ENST1.2");
        assert_eq!(transcript.matched_by, MatchBasis::StableIdWithoutVersion);
        assert_eq!(transcript.loci[0].start, 100);
        assert_eq!(transcript.loci[0].end, 400);
        assert!(resolver.resolve_str("transcript:ENST1.99").is_err());

        let exon = resolver.resolve_str("exon:ENSE1").unwrap();
        assert_eq!(exon.stable_id, "ENSE1.7");
        assert_eq!(exon.transcript_ids, vec!["ENST1.2"]);
        assert_eq!(exon.loci[0].start, 100);
    }

    #[test]
    fn duplicate_gene_symbols_are_explicitly_ambiguous() {
        let (_path, annotation) = annotation();
        let resolver = IntentResolver::from_annotation(&annotation, identity()).unwrap();
        let error = resolver.resolve_str("gene:DUP").unwrap_err();
        let ResolutionError::Ambiguous { candidates, .. } = error else {
            panic!("expected ambiguity")
        };
        assert_eq!(candidates.len(), 2);
        assert_eq!(
            candidates
                .iter()
                .map(|candidate| candidate.stable_id.as_str())
                .collect::<BTreeSet<_>>(),
            BTreeSet::from(["ENSG2.1", "ENSG3.1"])
        );
    }

    #[test]
    fn compiled_v2_preserves_intent_identifiers() {
        let (_path, annotation) = annotation();
        let compiled = TempPath::new("aic");
        annotation.write_compiled(&compiled.0).unwrap();
        let resolver = IntentResolver::from_path(&compiled.0, identity()).unwrap();
        assert_eq!(
            resolver.resolve_str("ENST2.8").unwrap().kind,
            FeatureKind::Transcript
        );
        assert_eq!(
            resolver.resolve_str("ENSE4.1").unwrap().kind,
            FeatureKind::Exon
        );
    }

    #[test]
    fn overlapping_and_touching_source_exons_keep_exact_identifier_loci() {
        let path = TempPath::new("gtf");
        let mut file = std::fs::File::create(&path.0).unwrap();
        writeln!(file, "chr1\tX\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_id \"E1\";").unwrap();
        writeln!(file, "chr1\tX\texon\t151\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_id \"E2\";").unwrap();
        writeln!(file, "chr1\tX\texon\t251\t300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_id \"E3\";").unwrap();
        drop(file);

        let annotation = Annotation::from_gtf(&path.0).unwrap();
        assert_eq!(
            annotation.transcripts[0].exons,
            vec![crate::Exon {
                start: 100,
                end: 300
            }]
        );
        let expected = [("E1", 100, 200), ("E2", 150, 250), ("E3", 250, 300)];
        let resolver = IntentResolver::from_annotation(&annotation, identity()).unwrap();
        for (id, start, end) in expected {
            let resolved = resolver.resolve_str(&format!("exon:{id}")).unwrap();
            assert_eq!((resolved.loci[0].start, resolved.loci[0].end), (start, end));
        }

        let compiled = TempPath::new("aic");
        annotation.write_compiled(&compiled.0).unwrap();
        let restored = Annotation::from_compiled(&compiled.0).unwrap();
        assert_eq!(
            restored.transcripts[0].exons,
            annotation.transcripts[0].exons
        );
        assert_eq!(restored.source_exons, annotation.source_exons);
        let resolver = IntentResolver::from_annotation(&restored, identity()).unwrap();
        for (id, start, end) in expected {
            let resolved = resolver.resolve_str(&format!("exon:{id}")).unwrap();
            assert_eq!((resolved.loci[0].start, resolved.loci[0].end), (start, end));
        }
    }

    #[test]
    fn identity_is_mandatory_and_prefixes_are_strict() {
        assert!(AnnotationIdentity::new("", "GENCODE 49").is_err());
        assert!(AnnotationIdentity::new("GRCh38", "  ").is_err());
        assert!(identity().with_digest("not-a-content-digest").is_err());
        assert_eq!(
            IdentifierQuery::parse("gene:TP53").unwrap().kind,
            QueryKind::Gene
        );
        assert_eq!(
            IdentifierQuery::parse("TP53").unwrap().kind,
            QueryKind::Auto
        );
        assert!(IdentifierQuery::parse("exon:").is_err());
    }

    #[test]
    fn path_construction_verifies_annotation_content_digest() {
        let (path, annotation) = annotation();
        let digest = format!(
            "blake3:{}",
            blake3::hash(&std::fs::read(&path.0).unwrap()).to_hex()
        );
        assert_eq!(
            IntentResolver::from_path(&path.0, identity())
                .unwrap()
                .identity()
                .digest
                .as_deref(),
            Some(digest.as_str())
        );
        let loaded = BoundAnnotation::from_path(&path.0, identity()).unwrap();
        assert_eq!(loaded.identity().digest.as_deref(), Some(digest.as_str()));
        assert_eq!(loaded.annotation().gene_ids, annotation.gene_ids);
        assert_eq!(
            IntentResolver::from_bound_annotation(&loaded)
                .unwrap()
                .resolve_str("gene:ENSG1.4")
                .unwrap()
                .stable_id,
            "ENSG1.4"
        );
        let bound = identity().with_digest(&digest).unwrap();
        assert_eq!(
            IntentResolver::from_path(&path.0, bound.clone())
                .unwrap()
                .identity()
                .digest
                .as_deref(),
            Some(digest.as_str())
        );

        let wrong = identity()
            .with_digest(format!("blake3:{}", "0".repeat(64)))
            .unwrap();
        assert!(matches!(
            IntentResolver::from_path(&path.0, wrong),
            Err(ResolverBuildError::InvalidIdentity(_))
        ));
        assert!(matches!(
            IntentResolver::from_annotation(&annotation, bound),
            Err(ResolverBuildError::InvalidIdentity(_))
        ));
    }

    #[cfg(unix)]
    #[test]
    fn digest_and_parser_use_the_same_open_annotation_snapshot() {
        let (path, _annotation) = annotation();
        let original = std::fs::read(&path.0).unwrap();
        let open_original = std::fs::File::open(&path.0).unwrap();
        let bound = identity()
            .with_digest(format!("blake3:{}", blake3::hash(&original).to_hex()))
            .unwrap();

        let replacement = TempPath::new("gtf");
        std::fs::write(
            &replacement.0,
            b"chr9\tX\texon\t901\t999\t.\t+\t.\tgene_id \"REPLACEMENT\"; transcript_id \"REPLACEMENT_TX\"; exon_id \"REPLACEMENT_EXON\";\n",
        )
        .unwrap();
        std::fs::rename(&replacement.0, &path.0).unwrap();

        let resolver = IntentResolver::from_open_file(&path.0, open_original, bound).unwrap();
        assert_eq!(
            resolver.resolve_str("gene:ENSG1.4").unwrap().stable_id,
            "ENSG1.4"
        );
        assert!(matches!(
            resolver.resolve_str("gene:REPLACEMENT"),
            Err(ResolutionError::NotFound { .. })
        ));
    }

    #[test]
    fn partial_identifier_metadata_never_turns_unknown_into_not_found() {
        let (_path, mut annotation) = annotation();
        annotation.transcript_ids[1] = None;
        annotation.source_exons[0][1].id = None;
        let resolver = IntentResolver::from_annotation(&annotation, identity()).unwrap();

        assert_eq!(
            resolver
                .resolve_str("transcript:ENST1.2")
                .unwrap()
                .stable_id,
            "ENST1.2"
        );
        assert_eq!(
            resolver.resolve_str("exon:ENSE1.7").unwrap().stable_id,
            "ENSE1.7"
        );
        assert!(matches!(
            resolver.resolve_str("transcript:UNKNOWN"),
            Err(ResolutionError::IdentifierMetadataUnavailable { .. })
        ));
        assert!(matches!(
            resolver.resolve_str("exon:UNKNOWN"),
            Err(ResolutionError::IdentifierMetadataUnavailable { .. })
        ));
    }
}
