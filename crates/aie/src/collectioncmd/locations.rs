//! Locations are untrusted hints; expected identities always come from the collection.
use super::*;
use serde::Deserialize;

const MAX_LOCATION_BYTES: u64 = 4 * 1024 * 1024;
pub(super) const COLLECTION_SCHEME: &str = "aicollection-directory-root-v1";

#[derive(Default)]
pub(super) struct Locations {
    paths: BTreeMap<String, PathBuf>,
    pub provenance: Option<serde_json::Value>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    schema_version: u32,
    locations: Vec<Entry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Entry {
    identity: String,
    path: PathBuf,
}

impl Locations {
    pub fn load(path: Option<&Path>) -> Result<Self> {
        let Some(path) = path else {
            return Ok(Self::default());
        };
        let mut bytes = Vec::new();
        std::fs::File::open(path)
            .with_context(|| format!("opening locations {}", path.display()))?
            .take(MAX_LOCATION_BYTES + 1)
            .read_to_end(&mut bytes)?;
        if bytes.len() as u64 > MAX_LOCATION_BYTES {
            bail!("locations manifest exceeds 4 MiB limit");
        }
        let manifest: Manifest =
            serde_json::from_slice(&bytes).context("parsing locations manifest")?;
        if manifest.schema_version != 1 {
            bail!(
                "unsupported locations schema_version {}; expected 1",
                manifest.schema_version
            );
        }
        if manifest.locations.len() > MAX_ARCHIVES_PER_LAYER {
            bail!("too many location entries");
        }
        let base = path.parent().unwrap_or_else(|| Path::new("."));
        let mut paths = BTreeMap::new();
        for entry in manifest.locations {
            let (scheme, digest) = entry
                .identity
                .split_once(':')
                .context("location identity must be SCHEME:HEX")?;
            if !matches!(
                scheme,
                ROOTED_DIRECTORY_SCHEME | LEGACY_FULL_FILE_SCHEME | COLLECTION_SCHEME
            ) {
                bail!("unsupported location identity scheme {scheme}");
            }
            digest_from_hex(digest, "location identity")?;
            let key = format!("{scheme}:{}", digest.to_ascii_lowercase());
            if entry.path.as_os_str().is_empty() {
                bail!("empty location path for {key}");
            }
            let resolved = if entry.path.is_absolute() {
                entry.path
            } else {
                base.join(entry.path)
            };
            if paths.insert(key.clone(), resolved).is_some() {
                bail!("duplicate location identity {key}");
            }
        }
        Ok(Self {
            paths,
            provenance: Some(
                json!({"path": path, "blake3": blake3::hash(&bytes).to_hex().to_string()}),
            ),
        })
    }

    pub fn resolve(&self, scheme: &str, digest: &str, fallback: &Path, base: &Path) -> PathBuf {
        self.paths
            .get(&format!("{scheme}:{}", digest.to_ascii_lowercase()))
            .cloned()
            .unwrap_or_else(|| {
                if fallback.is_absolute() {
                    fallback.to_path_buf()
                } else {
                    base.join(fallback)
                }
            })
    }
}
