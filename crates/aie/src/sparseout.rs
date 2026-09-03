use anyhow::{bail, Context, Result};
use serde_json::{json, Value};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

type CompressedWriter = zstd::Encoder<'static, BufWriter<File>>;

/// Streaming, zero-reconstructible output for large cohort event queries. Dimensions and
/// catalogue presence are separate from nonzero count facts, so a missing fact is an exact zero
/// rather than an unknown observation.
pub struct SparseCohortWriter {
    out_dir: PathBuf,
    events: CompressedWriter,
    presence: CompressedWriter,
    counts: CompressedWriter,
    event_rows: usize,
    presence_rows: usize,
    count_rows: usize,
}

fn checked_field<'a>(value: &'a str, label: &str) -> Result<&'a str> {
    if value.contains(['\t', '\n', '\r']) {
        bail!("sparse cohort {label} contains a tab or line break");
    }
    Ok(value)
}

fn encoder(path: &Path) -> Result<CompressedWriter> {
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    zstd::Encoder::new(BufWriter::new(file), 6)
        .with_context(|| format!("creating zstd stream {}", path.display()))
}

impl SparseCohortWriter {
    pub fn create(out_dir: &Path) -> Result<Self> {
        if out_dir.exists() {
            bail!(
                "refusing to overwrite existing sparse cohort output {}",
                out_dir.display()
            );
        }
        std::fs::create_dir(out_dir)
            .with_context(|| format!("creating sparse cohort output {}", out_dir.display()))?;
        let mut events = encoder(&out_dir.join("events.tsv.zst"))?;
        let mut presence = encoder(&out_dir.join("presence.tsv.zst"))?;
        let mut counts = encoder(&out_dir.join("counts.tsv.zst"))?;
        writeln!(
            events,
            "event_id\tevent_type\tchrom\tinclusion_junctions\texclusion_junctions\tannotation_genes_json\tstrand\tfully_annotated"
        )?;
        writeln!(presence, "event_id\tsample")?;
        writeln!(
            counts,
            "event_id\tsample\taggregation\tgroup\tinclude_only\texclude_only\tboth\tcells\tselected_cells"
        )?;
        Ok(Self {
            out_dir: out_dir.to_owned(),
            events,
            presence,
            counts,
            event_rows: 0,
            presence_rows: 0,
            count_rows: 0,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn event(
        &mut self,
        event_id: &str,
        event_type: &str,
        chrom: &str,
        inclusion_junctions: &str,
        exclusion_junctions: &str,
        annotation_genes_json: &str,
        strand: &str,
        fully_annotated: &str,
    ) -> Result<()> {
        for (value, label) in [
            (event_id, "event ID"),
            (event_type, "event type"),
            (chrom, "chromosome"),
            (inclusion_junctions, "inclusion junctions"),
            (exclusion_junctions, "exclusion junctions"),
            (annotation_genes_json, "annotation genes"),
            (strand, "strand"),
            (fully_annotated, "annotation flag"),
        ] {
            checked_field(value, label)?;
        }
        writeln!(
            self.events,
            "{event_id}\t{event_type}\t{chrom}\t{inclusion_junctions}\t{exclusion_junctions}\t{annotation_genes_json}\t{strand}\t{fully_annotated}"
        )?;
        self.event_rows += 1;
        Ok(())
    }

    pub fn present(&mut self, event_id: &str, sample: &str) -> Result<()> {
        checked_field(event_id, "event ID")?;
        checked_field(sample, "sample ID")?;
        writeln!(self.presence, "{event_id}\t{sample}")?;
        self.presence_rows += 1;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn count(
        &mut self,
        event_id: &str,
        sample: &str,
        aggregation: &str,
        group: &str,
        include_only: usize,
        exclude_only: usize,
        both: usize,
        cells: usize,
        selected_cells: usize,
    ) -> Result<()> {
        if include_only == 0 && exclude_only == 0 && both == 0 && cells == 0 {
            return Ok(());
        }
        for (value, label) in [
            (event_id, "event ID"),
            (sample, "sample ID"),
            (aggregation, "aggregation"),
            (group, "group"),
        ] {
            checked_field(value, label)?;
        }
        writeln!(
            self.counts,
            "{event_id}\t{sample}\t{aggregation}\t{group}\t{include_only}\t{exclude_only}\t{both}\t{cells}\t{selected_cells}"
        )?;
        self.count_rows += 1;
        Ok(())
    }

    pub fn finish(mut self, mut metadata: Value) -> Result<Value> {
        self.events.flush()?;
        self.presence.flush()?;
        self.counts.flush()?;
        self.events.finish()?.flush()?;
        self.presence.finish()?.flush()?;
        self.counts.finish()?.flush()?;
        let files = ["events.tsv.zst", "presence.tsv.zst", "counts.tsv.zst"];
        let compressed_bytes: u64 = files
            .iter()
            .map(|name| std::fs::metadata(self.out_dir.join(name)).map(|value| value.len()))
            .collect::<std::io::Result<Vec<_>>>()?
            .into_iter()
            .sum();
        let output = json!({
            "files": {
                "events": "events.tsv.zst",
                "presence": "presence.tsv.zst",
                "counts": "counts.tsv.zst",
            },
            "rows": {
                "events": self.event_rows,
                "presence": self.presence_rows,
                "nonzero_counts": self.count_rows,
            },
            "compressed_table_bytes": compressed_bytes,
        });
        metadata
            .as_object_mut()
            .context("sparse cohort metadata must be a JSON object")?
            .insert("output".into(), output);
        let temporary = self.out_dir.join("metadata.json.tmp");
        let final_path = self.out_dir.join("metadata.json");
        std::fs::write(&temporary, serde_json::to_string_pretty(&metadata)? + "\n")?;
        std::fs::rename(&temporary, &final_path)?;
        Ok(metadata)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_fields_reject_ambiguous_tsv_values() {
        assert!(checked_field("safe", "test").is_ok());
        assert!(checked_field("not\tsafe", "test").is_err());
        assert!(checked_field("not\nsafe", "test").is_err());
    }
}
