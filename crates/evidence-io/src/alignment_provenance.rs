//! Root-bound provenance for the alignments consumed by archive ingest.

use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};

use crate::genome::GenomeSig;
use crate::terminal_tail::TERMINAL_TAIL_RULE;

pub const ALIGNMENT_PROVENANCE_SCHEMA: &str = "gravlax.alignment-provenance.v1";
pub const MOLECULAR_EVIDENCE_SCHEMA: &str = "gravlax.molecular-evidence.v2";
pub const ALIGNMENT_PROVENANCE_SECTION: &str = "alignment.provenance";
pub const JUNCTION_CATALOGUE_SECTION: &str = "alignment.junction-catalogue";
pub const GENOME_REFERENCE_BINDING_SCHEMA: &str = "gravlax.genome-reference-binding.v1";

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum JunctionDiscoveryMode {
    Unspecified,
    OnePass,
    PerLibraryTwoPass,
    FrozenCatalogue,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ProvenanceStatus {
    VerifiedFromConsumedBytes,
    VerifiedBamHeader,
    DeclaredByCaller,
    Unspecified,
}

/// Exact identity of a held input that was read and stability-checked by ingest.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VerifiedFileIdentity {
    pub status: ProvenanceStatus,
    pub scheme: String,
    pub blake3: String,
    pub bytes: u64,
}

impl VerifiedFileIdentity {
    pub fn full_file_blake3(blake3: String, bytes: u64) -> Self {
        Self {
            status: ProvenanceStatus::VerifiedFromConsumedBytes,
            scheme: "full-file-blake3-v1".into(),
            blake3,
            bytes,
        }
    }

    fn validate(&self) -> Result<()> {
        if self.status != ProvenanceStatus::VerifiedFromConsumedBytes
            || self.scheme != "full-file-blake3-v1"
            || self.blake3.len() != 64
            || !self.blake3.bytes().all(|byte| byte.is_ascii_hexdigit())
        {
            bail!("alignment provenance contains an invalid verified file identity");
        }
        Ok(())
    }
}

/// One `@PG` record, in BAM-header order. These values are observations from the exact BAM whose
/// full-file digest is in the same manifest; they are not used to infer the declared junction
/// mode.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BamProgram {
    pub status: ProvenanceStatus,
    pub id: String,
    pub name: Option<String>,
    pub version: Option<String>,
    pub command_line: Option<String>,
    pub previous_program_id: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JunctionCatalogue {
    pub relationship_status: ProvenanceStatus,
    pub role: JunctionCatalogueRole,
    /// Root-bound section containing the exact catalogue bytes.
    pub section: String,
    pub identity: VerifiedFileIdentity,
    /// Nonempty, non-comment tabular records parsed from the exact bytes above.
    pub data_rows: u64,
}

/// Exact bytes Gravlax verified, plus the caller declaration that those bytes participated in
/// alignment. Verification of the file does not independently verify that relationship.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DeclaredAlignmentFile {
    pub relationship_status: ProvenanceStatus,
    pub locator: String,
    pub identity: VerifiedFileIdentity,
}

impl DeclaredAlignmentFile {
    fn validate(&self) -> Result<()> {
        if self.relationship_status != ProvenanceStatus::DeclaredByCaller || self.locator.is_empty()
        {
            bail!("alignment input has an invalid caller declaration");
        }
        self.identity.validate()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum JunctionCatalogueRole {
    PerLibraryPass1,
    FrozenExternal,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlignmentDeclaration {
    pub status: ProvenanceStatus,
    pub junction_discovery: JunctionDiscoveryMode,
    pub programs: Vec<BamProgram>,
    pub junction_catalogue: Option<JunctionCatalogue>,
    /// Exact annotation bytes supplied to the aligner, e.g. STAR `--sjdbGTFfile`, when declared.
    pub alignment_annotation: Option<DeclaredAlignmentFile>,
    /// Ordered source-read or other aligner inputs declared by the caller. Their bytes are
    /// verified by Gravlax; their role in alignment is a caller assertion.
    pub ordered_inputs: Vec<DeclaredAlignmentFile>,
    /// Optional aligner log containing resolved defaults, when retained by the caller.
    pub alignment_log: Option<DeclaredAlignmentFile>,
    pub chemistry: Option<String>,
    pub chemistry_status: ProvenanceStatus,
    /// Content identity or reproducible locator for the aligner index. This is caller-declared;
    /// Gravlax does not hash an opaque index directory implicitly.
    pub index_identity: Option<String>,
    pub index_identity_status: ProvenanceStatus,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlignmentInputs {
    pub bam: VerifiedFileIdentity,
    pub whitelist: VerifiedFileIdentity,
    pub genome_fasta: Option<VerifiedFileIdentity>,
    pub genome_signature: Option<GenomeSig>,
    /// The genome-to-alignment relationship is supplied by the caller; hashing the FASTA proves
    /// its bytes, not that those bytes generated the BAM.
    pub genome_relationship_status: ProvenanceStatus,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum GenomeBindingAction {
    IngestArchive,
    StampGenome,
}

/// Current reference bound to an archive for sequence-consulting queries. This is separate from
/// original alignment provenance because a later stamp cannot prove which reference made a BAM.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenomeReferenceBinding {
    pub schema: String,
    pub bound_by: GenomeBindingAction,
    pub relationship_status: ProvenanceStatus,
    pub identity: VerifiedFileIdentity,
    pub signature: GenomeSig,
    pub compatibility_check: String,
}

impl GenomeReferenceBinding {
    pub fn new(
        bound_by: GenomeBindingAction,
        identity: VerifiedFileIdentity,
        signature: GenomeSig,
    ) -> Self {
        Self {
            schema: GENOME_REFERENCE_BINDING_SCHEMA.into(),
            bound_by,
            relationship_status: ProvenanceStatus::DeclaredByCaller,
            identity,
            signature,
            compatibility_check: "all-archive-contig-names-present-v1".into(),
        }
    }

    pub fn validate(&self) -> Result<()> {
        if self.schema != GENOME_REFERENCE_BINDING_SCHEMA
            || self.relationship_status != ProvenanceStatus::DeclaredByCaller
            || self.compatibility_check != "all-archive-contig-names-present-v1"
        {
            bail!("invalid genome reference binding declaration");
        }
        self.identity.validate()?;
        self.signature.validate()
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IngestProvenance {
    pub program: String,
    pub version: String,
    pub locus_gap: u32,
    pub chunk_bp: u32,
    pub zstd_level: i32,
    pub molecule_chunk_streams: u32,
    pub molecule_codec: String,
    pub barcode_correction: String,
    pub umi_classes: String,
    pub unique_chain_reduction: String,
    pub multimapper_reduction: String,
    pub terminal_tail_rule: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlignmentProvenanceManifest {
    pub schema: String,
    pub molecular_evidence_schema: String,
    pub alignment: AlignmentDeclaration,
    pub inputs: AlignmentInputs,
    pub ingest: IngestProvenance,
}

impl AlignmentProvenanceManifest {
    pub fn validate(&self) -> Result<()> {
        if self.schema != ALIGNMENT_PROVENANCE_SCHEMA
            || self.molecular_evidence_schema != MOLECULAR_EVIDENCE_SCHEMA
        {
            bail!("unsupported alignment-provenance schema");
        }
        self.inputs.bam.validate()?;
        self.inputs.whitelist.validate()?;
        if let Some(identity) = &self.inputs.genome_fasta {
            identity.validate()?;
        }
        if self.inputs.genome_fasta.is_some() != self.inputs.genome_signature.is_some() {
            bail!(
                "alignment provenance genome identity and normalized signature must occur together"
            );
        }
        if let Some(signature) = &self.inputs.genome_signature {
            signature.validate()?;
        }
        match (
            self.inputs.genome_fasta.as_ref(),
            self.inputs.genome_relationship_status,
        ) {
            (Some(_), ProvenanceStatus::DeclaredByCaller) => {}
            (None, ProvenanceStatus::Unspecified) => {}
            _ => bail!("alignment provenance genome relationship has inconsistent status"),
        }
        if self.alignment.status
            != match self.alignment.junction_discovery {
                JunctionDiscoveryMode::Unspecified => ProvenanceStatus::Unspecified,
                _ => ProvenanceStatus::DeclaredByCaller,
            }
        {
            bail!("junction-discovery declaration has inconsistent status");
        }
        for program in &self.alignment.programs {
            if program.status != ProvenanceStatus::VerifiedBamHeader || program.id.is_empty() {
                bail!("alignment provenance contains an invalid BAM program record");
            }
        }
        match (
            self.alignment.junction_discovery,
            self.alignment.junction_catalogue.as_ref(),
        ) {
            (JunctionDiscoveryMode::PerLibraryTwoPass, Some(catalogue))
                if catalogue.role == JunctionCatalogueRole::PerLibraryPass1 => {}
            (JunctionDiscoveryMode::FrozenCatalogue, Some(catalogue))
                if catalogue.role == JunctionCatalogueRole::FrozenExternal => {}
            (JunctionDiscoveryMode::OnePass | JunctionDiscoveryMode::Unspecified, None) => {}
            _ => bail!("junction-discovery mode and catalogue role are inconsistent"),
        }
        if let Some(catalogue) = &self.alignment.junction_catalogue {
            if catalogue.relationship_status != ProvenanceStatus::DeclaredByCaller
                || catalogue.section != JUNCTION_CATALOGUE_SECTION
            {
                bail!("junction catalogue lacks a caller relationship declaration");
            }
            catalogue.identity.validate()?;
        }
        if let Some(annotation) = &self.alignment.alignment_annotation {
            annotation.validate()?;
        }
        for input in &self.alignment.ordered_inputs {
            input.validate()?;
        }
        if let Some(log) = &self.alignment.alignment_log {
            log.validate()?;
        }
        match (&self.alignment.chemistry, self.alignment.chemistry_status) {
            (Some(value), ProvenanceStatus::DeclaredByCaller) if !value.trim().is_empty() => {}
            (None, ProvenanceStatus::Unspecified) => {}
            _ => bail!("alignment chemistry has inconsistent declaration status"),
        }
        match (
            &self.alignment.index_identity,
            self.alignment.index_identity_status,
        ) {
            (Some(value), ProvenanceStatus::DeclaredByCaller) if !value.trim().is_empty() => {}
            (None, ProvenanceStatus::Unspecified) => {}
            _ => bail!("alignment-index identity has inconsistent declaration status"),
        }
        if self.ingest.program != "aie"
            || self.ingest.version.is_empty()
            || self.ingest.chunk_bp == 0
            || self.ingest.molecule_chunk_streams != 10
            || self.ingest.molecule_codec != "rans2"
            || self.ingest.barcode_correction != "unique-hamming1-quality-pseudocount-v1"
            || self.ingest.umi_classes != "global-cell-umi-equivalence-with-1mm-edges-v1"
            || self.ingest.unique_chain_reduction != "junction-chain-span-extremes-v1"
            || self.ingest.multimapper_reduction != "primary-relative-placement-pattern-v1"
        {
            bail!("alignment provenance contains invalid ingest parameters");
        }
        if self
            .ingest
            .terminal_tail_rule
            .as_deref()
            .is_some_and(|rule| rule != TERMINAL_TAIL_RULE)
        {
            bail!("alignment provenance declares an unsupported terminal-tail extraction rule");
        }
        Ok(())
    }

    /// Struct field order plus compact serde encoding defines the canonical section bytes.
    pub fn to_canonical_json(&self) -> Result<Vec<u8>> {
        self.validate()?;
        Ok(serde_json::to_vec(self)?)
    }

    pub fn from_json(bytes: &[u8]) -> Result<Self> {
        let manifest: Self = serde_json::from_slice(bytes)?;
        manifest.validate()?;
        if manifest.to_canonical_json()? != bytes {
            bail!("alignment-provenance JSON is not in canonical compact form");
        }
        Ok(manifest)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn file(byte: u8) -> VerifiedFileIdentity {
        VerifiedFileIdentity::full_file_blake3(format!("{byte:02x}").repeat(32), 7)
    }

    fn manifest() -> AlignmentProvenanceManifest {
        AlignmentProvenanceManifest {
            schema: ALIGNMENT_PROVENANCE_SCHEMA.into(),
            molecular_evidence_schema: MOLECULAR_EVIDENCE_SCHEMA.into(),
            alignment: AlignmentDeclaration {
                status: ProvenanceStatus::DeclaredByCaller,
                junction_discovery: JunctionDiscoveryMode::PerLibraryTwoPass,
                programs: vec![BamProgram {
                    status: ProvenanceStatus::VerifiedBamHeader,
                    id: "STAR".into(),
                    name: Some("STAR".into()),
                    version: Some("2.7".into()),
                    command_line: Some("STAR --twopassMode Basic".into()),
                    previous_program_id: None,
                }],
                junction_catalogue: Some(JunctionCatalogue {
                    relationship_status: ProvenanceStatus::DeclaredByCaller,
                    role: JunctionCatalogueRole::PerLibraryPass1,
                    section: JUNCTION_CATALOGUE_SECTION.into(),
                    identity: file(0xab),
                    data_rows: 12,
                }),
                alignment_annotation: None,
                ordered_inputs: Vec::new(),
                alignment_log: None,
                chemistry: Some("10x-3p-v3".into()),
                chemistry_status: ProvenanceStatus::DeclaredByCaller,
                index_identity: Some("blake3:abc".into()),
                index_identity_status: ProvenanceStatus::DeclaredByCaller,
            },
            inputs: AlignmentInputs {
                bam: file(0xcd),
                whitelist: file(0xef),
                genome_fasta: None,
                genome_signature: None,
                genome_relationship_status: ProvenanceStatus::Unspecified,
            },
            ingest: IngestProvenance {
                program: "aie".into(),
                version: "0.1.4".into(),
                locus_gap: 2_000,
                chunk_bp: 4_000_000,
                zstd_level: 19,
                molecule_chunk_streams: 10,
                molecule_codec: "rans2".into(),
                barcode_correction: "unique-hamming1-quality-pseudocount-v1".into(),
                umi_classes: "global-cell-umi-equivalence-with-1mm-edges-v1".into(),
                unique_chain_reduction: "junction-chain-span-extremes-v1".into(),
                multimapper_reduction: "primary-relative-placement-pattern-v1".into(),
                terminal_tail_rule: None,
            },
        }
    }

    #[test]
    fn canonical_manifest_roundtrips_and_rejects_noncanonical_json() {
        let manifest = manifest();
        let bytes = manifest.to_canonical_json().unwrap();
        assert_eq!(
            AlignmentProvenanceManifest::from_json(&bytes).unwrap(),
            manifest
        );
        let pretty = serde_json::to_vec_pretty(&manifest).unwrap();
        assert!(AlignmentProvenanceManifest::from_json(&pretty).is_err());
    }

    #[test]
    fn declared_mode_requires_matching_verified_catalogue_role() {
        let mut value = manifest();
        value.alignment.junction_discovery = JunctionDiscoveryMode::FrozenCatalogue;
        assert!(value.validate().is_err());
        value.alignment.junction_catalogue.as_mut().unwrap().role =
            JunctionCatalogueRole::FrozenExternal;
        assert!(value.validate().is_ok());
    }

    #[test]
    fn unknown_terminal_tail_rule_is_rejected() {
        let mut value = manifest();
        value.ingest.terminal_tail_rule = Some("unversioned-tail-guess".into());
        assert!(value.validate().is_err());
        value.ingest.terminal_tail_rule = Some(TERMINAL_TAIL_RULE.into());
        assert!(value.validate().is_ok());
    }
}
