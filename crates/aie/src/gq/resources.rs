//! Explicit, content-bound inputs. No network evaluation or implicit assembly aliases.
use super::model::{Truth, Type, Value};
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Member {
    pub sample: String,
    pub path: PathBuf,
    #[serde(default)]
    pub metadata: BTreeMap<String, serde_json::Value>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Federation {
    schema_version: u32,
    assembly: String,
    archives: Vec<Member>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Column {
    pub r#type: Type,
    #[serde(default)]
    pub optional: bool,
}
#[derive(Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Metadata {
    pub columns: BTreeMap<String, Column>,
    /// Namespaced sample -> barcode -> field/value object.
    #[serde(default)]
    pub cells: BTreeMap<String, BTreeMap<String, BTreeMap<String, serde_json::Value>>>,
}
pub struct Resources {
    pub sources: BTreeMap<String, Vec<Member>>,
    pub metadata: Metadata,
    pub digests: BTreeMap<String, String>,
    pub annotations: BTreeMap<String, anno::Annotation>,
}
impl Resources {
    pub fn open(
        bindings: &[String],
        metadata: Option<&Path>,
        project: Option<&Path>,
        assembly: &str,
    ) -> Result<Self> {
        let mut sources = BTreeMap::new();
        let mut digests = BTreeMap::new();
        let mut annotations = BTreeMap::new();
        if let Some(path) = project {
            let context = crate::projectcmd::find_project(Some(path))?;
            digests.insert("project".into(), context.manifest_digest);
            for (alias, r) in context.manifest.resources {
                if r.kind == crate::projectcmd::ResourceKind::Annotation {
                    let identity = r.annotation_identity.context(
                        "annotation binding requires assembly and immutable annotation label",
                    )?;
                    if identity.assembly != assembly || identity.annotation.is_empty() {
                        bail!("annotation {alias} has incompatible assembly or missing label");
                    }
                    let path = if r.external {
                        r.path
                    } else {
                        context.root.join(r.path)
                    };
                    let (annotation, digest) = anno::Annotation::from_open_file_with_digest(
                        std::fs::File::open(&path)?,
                        &path,
                    )?;
                    digests.insert(
                        format!(
                            "annotation:{alias}:{}:{}",
                            identity.assembly, identity.annotation
                        ),
                        digest,
                    );
                    annotations.insert(alias, annotation);
                    continue;
                }
                if r.kind == crate::projectcmd::ResourceKind::Archive {
                    if r.assembly.as_deref() != Some(assembly) {
                        bail!("project archive {alias} must explicitly declare matching assembly {assembly}");
                    }
                    let path = if r.external {
                        r.path
                    } else {
                        context.root.join(r.path)
                    };
                    sources.insert(
                        alias.clone(),
                        vec![Member {
                            sample: alias,
                            path,
                            metadata: BTreeMap::new(),
                        }],
                    );
                }
            }
        }
        for binding in bindings {
            let (alias, path) = binding
                .split_once('=')
                .context("--bind requires alias=archive.aie or alias=federation.json")?;
            if alias.is_empty() || sources.contains_key(alias) {
                bail!("empty or duplicate binding {alias}");
            }
            let path = PathBuf::from(path);
            let members = if path.extension().is_some_and(|e| e == "json") {
                let bytes = std::fs::read(&path)?;
                digests.insert(
                    format!("binding:{alias}"),
                    blake3::hash(&bytes).to_hex().to_string(),
                );
                let f: Federation = serde_json::from_slice(&bytes)?;
                if f.schema_version != 1 || f.assembly != assembly {
                    bail!("federation version or assembly mismatch for {alias}");
                }
                let base = path.parent().unwrap_or(Path::new("."));
                f.archives
                    .into_iter()
                    .map(|mut m| {
                        if m.path.is_relative() {
                            m.path = base.join(m.path);
                        }
                        m
                    })
                    .collect()
            } else {
                vec![Member {
                    sample: alias.into(),
                    path,
                    metadata: BTreeMap::new(),
                }]
            };
            sources.insert(alias.into(), members);
        }
        let metadata = if let Some(path) = metadata {
            let bytes = std::fs::read(path)?;
            digests.insert("metadata".into(), blake3::hash(&bytes).to_hex().to_string());
            serde_json::from_slice(&bytes)?
        } else {
            Metadata::default()
        };
        let out = Self {
            sources,
            metadata,
            digests,
            annotations,
        };
        out.fields()?;
        Ok(out)
    }
    pub fn members(&self, alias: &str) -> Result<&[Member]> {
        let members = self
            .sources
            .get(alias)
            .with_context(|| format!("unbound source @{alias}; use --bind or --project"))?;
        if members.is_empty() {
            bail!("empty federation {alias}");
        }
        let mut names = BTreeSet::new();
        let mut paths = BTreeSet::new();
        for m in members {
            if m.sample.is_empty()
                || !names.insert(&m.sample)
                || !paths.insert(std::fs::canonicalize(&m.path)?)
            {
                bail!("duplicate or empty sample/archive in @{alias}");
            }
        }
        Ok(members)
    }
    pub fn fields(&self) -> Result<BTreeMap<String, Type>> {
        let mut fields = BTreeMap::from([
            ("sample".into(), Type::String),
            ("cell".into(), Type::String),
            ("cell.id".into(), Type::String),
            ("unit.id".into(), Type::String),
        ]);
        for (name, c) in &self.metadata.columns {
            if fields.contains_key(name) || name.starts_with("reads.") {
                bail!("metadata cannot shadow reserved field {name}");
            }
            if !matches!(
                c.r#type,
                Type::String | Type::Number | Type::Truth | Type::List(_)
            ) {
                bail!("invalid metadata type for {name}");
            }
            fields.insert(name.clone(), c.r#type.clone());
        }
        Ok(fields)
    }
    pub fn row(&self, member: &Member, barcode: &str) -> Result<BTreeMap<String, Value>> {
        let cell_id = format!("{}:{barcode}", member.sample);
        let mut fields = BTreeMap::from([
            ("sample".into(), Value::String(member.sample.clone())),
            ("cell".into(), Value::String(cell_id)),
            ("cell.id".into(), Value::String(barcode.into())),
        ]);
        let cell = self
            .metadata
            .cells
            .get(&member.sample)
            .and_then(|rows| rows.get(barcode));
        for (name, column) in &self.metadata.columns {
            let raw = cell
                .and_then(|m| m.get(name))
                .or_else(|| member.metadata.get(name));
            let value=match raw{Some(v)=>typed(v,&column.r#type)?,None if column.optional=>Value::Null,None=>bail!("required metadata column {name} missing for {}:{barcode}; declare optional to permit missing schema entries",member.sample)};
            fields.insert(name.clone(), value);
        }
        Ok(fields)
    }
    pub fn feature(
        &self,
        alias: &str,
        name: &str,
        transcript: bool,
        library: Option<&str>,
    ) -> Result<Value> {
        use super::model::{Pattern, Region};
        let annotation = self.annotations.get(alias).with_context(|| {
            format!("annotation @{alias} must be a pinned project annotation resource")
        })?;
        let opposite=library.context("annotation-derived strand is transcript-frame; bind header library = \"same\" or \"opposite\"")?=="opposite";
        let indices: Vec<_> = if transcript {
            let matches: Vec<_> = annotation
                .transcript_ids
                .iter()
                .enumerate()
                .filter(|(_, id)| id.as_deref() == Some(name))
                .map(|(i, _)| i)
                .collect();
            if matches.len() != 1 {
                bail!(
                    "transcript {name:?} is missing or ambiguous ({} matches)",
                    matches.len()
                );
            }
            matches
        } else {
            let genes: Vec<_> = annotation
                .gene_ids
                .iter()
                .zip(&annotation.gene_names)
                .enumerate()
                .filter(|(_, (id, n))| id.as_str() == name || n.as_str() == name)
                .map(|(i, _)| i as u32)
                .collect();
            if genes.len() != 1 {
                bail!(
                    "gene {name:?} is missing or ambiguous ({} matches); use a stable ID",
                    genes.len()
                );
            }
            annotation
                .transcripts
                .iter()
                .enumerate()
                .filter(|(_, t)| t.gene == genes[0])
                .map(|(i, _)| i)
                .collect()
        };
        let first =
            &annotation.transcripts[*indices.first().context("feature has no transcripts")?];
        let chrom = annotation
            .chrom_ids
            .iter()
            .find(|(_, id)| **id == first.chrom)
            .context("annotation chromosome missing")?
            .0
            .clone();
        let reverse = Some(first.strand_rev ^ opposite);
        let mut intervals = Vec::new();
        let mut junctions = BTreeSet::new();
        for i in indices {
            let t = &annotation.transcripts[i];
            if t.chrom != first.chrom || t.strand_rev != first.strand_rev {
                bail!("feature spans incompatible chromosomes or strands");
            }
            intervals.extend(t.exons.iter().map(|e| (e.start, e.end)));
            junctions.extend(t.exons.windows(2).map(|w| (w[0].end, w[1].start)));
        }
        let mut exons = Region {
            chrom: chrom.clone(),
            intervals,
            reverse,
        };
        exons.union();
        let start = exons.intervals.first().context("feature has no exons")?.0;
        let end = exons.intervals.last().unwrap().1;
        let span = Region {
            chrom: chrom.clone(),
            intervals: vec![(start, end)],
            reverse,
        };
        let make = |js| Pattern {
            chrom: chrom.clone(),
            reverse,
            junctions: js,
            subsequence: false,
            left: 0,
            right: 0,
        };
        let junction_path = if transcript && !junctions.is_empty() {
            Some(make(junctions.iter().copied().collect()))
        } else {
            None
        };
        let junctions = junctions.into_iter().map(|j| make(vec![j])).collect();
        Ok(Value::Feature(Box::new(super::model::Feature {
            transcript,
            span,
            exons,
            junctions,
            junction_path,
        })))
    }
}
fn typed(v: &serde_json::Value, t: &Type) -> Result<Value> {
    if v.is_null() {
        return Ok(Value::Null);
    }
    Ok(match t {
        Type::String => Value::String(v.as_str().context("metadata value must be string")?.into()),
        Type::Number => Value::Number(v.as_f64().context("metadata value must be numeric")?),
        Type::Truth => Value::Truth(Truth::from_bool(
            v.as_bool().context("metadata value must be Boolean")?,
        )),
        Type::List(t) => Value::List(
            v.as_array()
                .context("metadata value must be list")?
                .iter()
                .map(|v| typed(v, t))
                .collect::<Result<_>>()?,
        ),
        _ => bail!("unsupported metadata type"),
    })
}
