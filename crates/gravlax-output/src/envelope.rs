use crate::table::valid_schema_id;
use crate::OutputError;
use serde::{de::Error as _, Deserialize, Deserializer, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const ENVELOPE_SCHEMA: &str = "gravlax.result-envelope.v1";

/// One explicitly named annotation input. Paired/counterfactual analyses use roles such as
/// `before` and `after`; ordinary single-annotation commands may continue to use the legacy
/// singular provenance fields.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnnotationProvenance {
    pub role: String,
    pub assembly: String,
    pub annotation: String,
    pub digest: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Producer {
    pub name: String,
    pub version: String,
}

impl Default for Producer {
    fn default() -> Self {
        Self {
            name: "gravlax".into(),
            version: env!("CARGO_PKG_VERSION").into(),
        }
    }
}

impl Producer {
    fn validate(&self) -> Result<(), OutputError> {
        if self.name.trim().is_empty() || self.version.trim().is_empty() {
            return Err(OutputError::InvalidSchema(
                "result producer name and version must not be empty".into(),
            ));
        }
        Ok(())
    }
}

/// Reproducibility identity shared by every output format. No wall-clock timestamp is inserted,
/// so otherwise identical queries can remain byte-identical.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Provenance {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub archives: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub annotation: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub annotation_digest: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub annotations: Vec<AnnotationProvenance>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub parameters: BTreeMap<String, serde_json::Value>,
}

impl Provenance {
    fn validate(&self) -> Result<(), OutputError> {
        let mut archives = BTreeSet::new();
        for archive in &self.archives {
            if archive.trim().is_empty() {
                return Err(OutputError::InvalidSchema(
                    "provenance archive identities must not be empty".into(),
                ));
            }
            if !archives.insert(archive) {
                return Err(OutputError::InvalidSchema(format!(
                    "duplicate provenance archive identity {archive:?}"
                )));
            }
        }
        for (label, value) in [
            ("assembly", self.assembly.as_deref()),
            ("annotation", self.annotation.as_deref()),
            ("annotation_digest", self.annotation_digest.as_deref()),
        ] {
            if value.is_some_and(|value| value.trim().is_empty()) {
                return Err(OutputError::InvalidSchema(format!(
                    "provenance {label} must not be empty when supplied"
                )));
            }
        }
        let mut roles = BTreeSet::new();
        for annotation in &self.annotations {
            if annotation.role.trim().is_empty()
                || annotation.assembly.trim().is_empty()
                || annotation.annotation.trim().is_empty()
                || annotation.digest.trim().is_empty()
            {
                return Err(OutputError::InvalidSchema(
                    "annotation provenance role, assembly, label, and digest must not be empty"
                        .into(),
                ));
            }
            if !roles.insert(&annotation.role) {
                return Err(OutputError::InvalidSchema(format!(
                    "duplicate annotation provenance role {:?}",
                    annotation.role
                )));
            }
        }
        if self.parameters.keys().any(|name| name.trim().is_empty()) {
            return Err(OutputError::InvalidSchema(
                "provenance parameter names must not be empty".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ResultContext {
    #[serde(default)]
    pub producer: Producer,
    #[serde(default)]
    pub provenance: Provenance,
    #[serde(default)]
    pub warnings: Vec<String>,
}

impl ResultContext {
    pub fn validate(&self) -> Result<(), OutputError> {
        self.producer.validate()?;
        self.provenance.validate()?;
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ResultEnvelope<T> {
    #[serde(rename = "$schema")]
    pub envelope_schema: String,
    pub result_schema: String,
    pub producer: Producer,
    pub provenance: Provenance,
    #[serde(default)]
    pub warnings: Vec<String>,
    pub data: T,
}

impl<T> ResultEnvelope<T> {
    pub fn new(
        result_schema: impl Into<String>,
        context: ResultContext,
        data: T,
    ) -> Result<Self, OutputError> {
        let result_schema = result_schema.into();
        if !valid_schema_id(&result_schema) {
            return Err(OutputError::InvalidSchema(
                "result schema id must be non-empty ASCII letters, digits, '.', '_' or '-'".into(),
            ));
        }
        context.validate()?;
        Ok(Self {
            envelope_schema: ENVELOPE_SCHEMA.into(),
            result_schema,
            producer: context.producer,
            provenance: context.provenance,
            warnings: context.warnings,
            data,
        })
    }
}

impl<'de, T> Deserialize<'de> for ResultEnvelope<T>
where
    T: Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct WireEnvelope<T> {
            #[serde(rename = "$schema")]
            envelope_schema: String,
            result_schema: String,
            producer: Producer,
            provenance: Provenance,
            #[serde(default)]
            warnings: Vec<String>,
            data: T,
        }

        let wire = WireEnvelope::deserialize(deserializer)?;
        if wire.envelope_schema != ENVELOPE_SCHEMA {
            return Err(D::Error::custom(format!(
                "unsupported result envelope schema {:?}; expected {ENVELOPE_SCHEMA}",
                wire.envelope_schema
            )));
        }
        let context = ResultContext {
            producer: wire.producer,
            provenance: wire.provenance,
            warnings: wire.warnings,
        };
        Self::new(wire.result_schema, context, wire.data).map_err(D::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn envelope_round_trip_validates_its_contract() {
        let envelope = ResultEnvelope::new(
            "gravlax.test.envelope.v1",
            ResultContext::default(),
            serde_json::json!({"ok": true}),
        )
        .unwrap();
        let encoded = serde_json::to_vec(&envelope).unwrap();
        let decoded: ResultEnvelope<serde_json::Value> = serde_json::from_slice(&encoded).unwrap();
        assert_eq!(decoded, envelope);

        let mut value = serde_json::to_value(&envelope).unwrap();
        value["$schema"] = serde_json::json!("gravlax.result-envelope.v999");
        assert!(serde_json::from_value::<ResultEnvelope<serde_json::Value>>(value).is_err());
    }

    #[test]
    fn envelope_rejects_invalid_identity_fields() {
        let invalid = serde_json::json!({
            "$schema": ENVELOPE_SCHEMA,
            "result_schema": "not a schema",
            "producer": {"name": "", "version": "1"},
            "provenance": {},
            "warnings": [],
            "data": null
        });
        assert!(serde_json::from_value::<ResultEnvelope<serde_json::Value>>(invalid).is_err());

        let mut context = ResultContext::default();
        context.provenance.archives = vec!["same".into(), "same".into()];
        assert!(ResultEnvelope::new("gravlax.test.v1", context, ()).is_err());

        let paired = AnnotationProvenance {
            role: "before".into(),
            assembly: "GRCh38".into(),
            annotation: "GENCODE 48".into(),
            digest: format!("blake3:{}", "a".repeat(64)),
        };
        let mut context = ResultContext::default();
        context.provenance.annotations = vec![paired.clone(), paired];
        assert!(ResultEnvelope::new("gravlax.test.v1", context, ()).is_err());
    }
}
