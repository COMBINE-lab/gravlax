use crate::{ResultContext, ENVELOPE_SCHEMA};
use serde::{de::Error as _, Deserialize, Deserializer, Serialize};
use std::collections::BTreeSet;
use std::fmt;
use std::io::Write;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DataType {
    String,
    Int64,
    #[serde(rename = "uint64")]
    UInt64,
    Float64,
    Boolean,
    /// Arbitrary JSON. Arrow adapters map this to canonical JSON in UTF-8 with the
    /// `gravlax.logical_type=json` field metadata key.
    Json,
}

impl DataType {
    pub fn name(self) -> &'static str {
        match self {
            Self::String => "string",
            Self::Int64 => "int64",
            Self::UInt64 => "uint64",
            Self::Float64 => "float64",
            Self::Boolean => "boolean",
            Self::Json => "json",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Field {
    pub name: String,
    pub data_type: DataType,
    #[serde(default)]
    pub nullable: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
}

impl Field {
    pub fn new(name: impl Into<String>, data_type: DataType) -> Self {
        Self {
            name: name.into(),
            data_type,
            nullable: false,
            description: None,
        }
    }

    pub fn nullable(mut self) -> Self {
        self.nullable = true;
        self
    }
}

/// The logical meaning of a table's rows. This is independent of the order in which a producer
/// happens to emit them.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RowSemantics {
    Set,
    Multiset,
    Sequence,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SortDirection {
    Ascending,
    Descending,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OrderKey {
    pub field: String,
    pub direction: SortDirection,
}

impl OrderKey {
    pub fn ascending(field: impl Into<String>) -> Self {
        Self {
            field: field.into(),
            direction: SortDirection::Ascending,
        }
    }

    pub fn descending(field: impl Into<String>) -> Self {
        Self {
            field: field.into(),
            direction: SortDirection::Descending,
        }
    }
}

/// Logical row metadata. An absent `ordered_by` means that row order is unspecified; writers must
/// not sort merely to make an output deterministic. `key` identifies rows when one exists, while
/// `ordered_by` describes an ordering that the producer already guarantees.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TableSemantics {
    pub row_semantics: RowSemantics,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub key: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ordered_by: Option<Vec<OrderKey>>,
}

impl TableSemantics {
    pub fn new(row_semantics: RowSemantics) -> Self {
        Self {
            row_semantics,
            key: None,
            ordered_by: None,
        }
    }

    pub fn with_key<I, S>(mut self, fields: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.key = Some(fields.into_iter().map(Into::into).collect());
        self
    }

    pub fn ordered_by<I>(mut self, fields: I) -> Self
    where
        I: IntoIterator<Item = OrderKey>,
    {
        self.ordered_by = Some(fields.into_iter().collect());
        self
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct TableSchema {
    /// Versioned, command-specific identifier such as `gravlax.query.region.v1`.
    pub id: String,
    pub fields: Vec<Field>,
    /// Absent only for schemas created before row semantics became part of the uniform contract.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub semantics: Option<TableSemantics>,
}

impl TableSchema {
    pub fn new(id: impl Into<String>, fields: Vec<Field>) -> Result<Self, OutputError> {
        let schema = Self {
            id: id.into(),
            fields,
            semantics: None,
        };
        schema.validate()?;
        Ok(schema)
    }

    pub fn with_semantics(mut self, semantics: TableSemantics) -> Result<Self, OutputError> {
        self.semantics = Some(semantics);
        self.validate()?;
        Ok(self)
    }

    pub fn validate(&self) -> Result<(), OutputError> {
        if !valid_schema_id(&self.id) {
            return Err(OutputError::InvalidSchema(
                "schema id must be non-empty ASCII letters, digits, '.', '_' or '-'".into(),
            ));
        }
        if self.fields.is_empty() {
            return Err(OutputError::InvalidSchema(
                "table must have at least one field".into(),
            ));
        }
        let mut names = BTreeSet::new();
        for field in &self.fields {
            if field.name.is_empty() || field.name.contains(['\t', '\r', '\n']) {
                return Err(OutputError::InvalidSchema(format!(
                    "invalid field name {:?}",
                    field.name
                )));
            }
            if !names.insert(&field.name) {
                return Err(OutputError::InvalidSchema(format!(
                    "duplicate field name {:?}",
                    field.name
                )));
            }
        }
        if let Some(semantics) = &self.semantics {
            validate_semantic_fields("key", semantics.key.as_deref(), &names)?;
            let ordered_by = semantics.ordered_by.as_deref().map(|keys| {
                keys.iter()
                    .map(|key| key.field.as_str())
                    .collect::<Vec<_>>()
            });
            validate_semantic_fields("ordered_by", ordered_by.as_deref(), &names)?;
        }
        Ok(())
    }
}

fn validate_semantic_fields<S>(
    label: &str,
    fields: Option<&[S]>,
    schema_fields: &BTreeSet<&String>,
) -> Result<(), OutputError>
where
    S: AsRef<str>,
{
    let Some(fields) = fields else {
        return Ok(());
    };
    if fields.is_empty() {
        return Err(OutputError::InvalidSchema(format!(
            "{label} must be absent rather than empty"
        )));
    }
    let mut seen = BTreeSet::new();
    for field in fields {
        let field = field.as_ref();
        if !schema_fields
            .iter()
            .any(|candidate| candidate.as_str() == field)
        {
            return Err(OutputError::InvalidSchema(format!(
                "{label} references unknown field {field:?}"
            )));
        }
        if !seen.insert(field) {
            return Err(OutputError::InvalidSchema(format!(
                "{label} repeats field {field:?}"
            )));
        }
    }
    Ok(())
}

pub(crate) fn valid_schema_id(id: &str) -> bool {
    !id.trim().is_empty()
        && id
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
}

#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(untagged)]
pub enum ScalarValue {
    Null,
    String(String),
    Int64(i64),
    UInt64(u64),
    Float64(f64),
    Boolean(bool),
    Json(serde_json::Value),
}

impl ScalarValue {
    fn compatible_with(&self, field: &Field) -> bool {
        match self {
            Self::Null => field.nullable,
            Self::String(_) => field.data_type == DataType::String,
            Self::Int64(_) => field.data_type == DataType::Int64,
            Self::UInt64(_) => field.data_type == DataType::UInt64,
            Self::Float64(value) => field.data_type == DataType::Float64 && value.is_finite(),
            Self::Boolean(_) => field.data_type == DataType::Boolean,
            // A table null has one representation. Treating `Json(null)` as a second spelling
            // would make a nullable JSON cell impossible to round-trip through ordinary JSON.
            Self::Json(value) => field.data_type == DataType::Json && !value.is_null(),
        }
    }
}

impl From<&str> for ScalarValue {
    fn from(value: &str) -> Self {
        Self::String(value.into())
    }
}

impl From<String> for ScalarValue {
    fn from(value: String) -> Self {
        Self::String(value)
    }
}

impl From<i64> for ScalarValue {
    fn from(value: i64) -> Self {
        Self::Int64(value)
    }
}

impl From<u64> for ScalarValue {
    fn from(value: u64) -> Self {
        Self::UInt64(value)
    }
}

impl From<f64> for ScalarValue {
    fn from(value: f64) -> Self {
        Self::Float64(value)
    }
}

impl From<bool> for ScalarValue {
    fn from(value: bool) -> Self {
        Self::Boolean(value)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct TypedTable {
    pub schema: TableSchema,
    pub rows: Vec<Vec<ScalarValue>>,
}

impl<'de> Deserialize<'de> for TypedTable {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct WireTable {
            schema: TableSchema,
            rows: Vec<Vec<serde_json::Value>>,
        }

        let wire = WireTable::deserialize(deserializer)?;
        wire.schema.validate().map_err(D::Error::custom)?;
        let mut rows = Vec::with_capacity(wire.rows.len());
        for (row_index, row) in wire.rows.into_iter().enumerate() {
            if row.len() != wire.schema.fields.len() {
                return Err(D::Error::custom(format!(
                    "invalid output row: row {row_index} has {} values for {} fields",
                    row.len(),
                    wire.schema.fields.len()
                )));
            }
            let mut typed = Vec::with_capacity(row.len());
            for (column, (value, field)) in row.into_iter().zip(&wire.schema.fields).enumerate() {
                let value = scalar_from_json(value, field).map_err(|message| {
                    D::Error::custom(format!(
                        "invalid output row: row {row_index}, column {column} ({}): {message}",
                        field.name
                    ))
                })?;
                typed.push(value);
            }
            rows.push(typed);
        }
        Self::new(wire.schema, rows).map_err(D::Error::custom)
    }
}

fn scalar_from_json(value: serde_json::Value, field: &Field) -> Result<ScalarValue, String> {
    if value.is_null() {
        return if field.nullable {
            Ok(ScalarValue::Null)
        } else {
            Err(format!(
                "null does not match non-nullable {}",
                field.data_type.name()
            ))
        };
    }
    let mismatch = || format!("value does not match {}", field.data_type.name());
    match field.data_type {
        DataType::String => value
            .as_str()
            .map(|value| ScalarValue::String(value.to_owned()))
            .ok_or_else(mismatch),
        DataType::Int64 => value.as_i64().map(ScalarValue::Int64).ok_or_else(mismatch),
        DataType::UInt64 => value.as_u64().map(ScalarValue::UInt64).ok_or_else(mismatch),
        DataType::Float64 => match value {
            serde_json::Value::Number(number) => number
                .as_f64()
                .filter(|value| value.is_finite())
                .map(ScalarValue::Float64)
                .ok_or_else(mismatch),
            _ => Err(mismatch()),
        },
        DataType::Boolean => value
            .as_bool()
            .map(ScalarValue::Boolean)
            .ok_or_else(mismatch),
        DataType::Json => Ok(ScalarValue::Json(value)),
    }
}

impl TypedTable {
    pub fn new(schema: TableSchema, rows: Vec<Vec<ScalarValue>>) -> Result<Self, OutputError> {
        schema.validate()?;
        for (index, row) in rows.iter().enumerate() {
            validate_row(&schema, index, row)?;
        }
        Ok(Self { schema, rows })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    Text,
    Tsv,
    Json,
    Arrow,
    Mex,
}

impl FromStr for OutputFormat {
    type Err = OutputError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "text" => Ok(Self::Text),
            "tsv" => Ok(Self::Tsv),
            "json" => Ok(Self::Json),
            "arrow" | "arrow-ipc" => Ok(Self::Arrow),
            "mex" | "mtx" => Ok(Self::Mex),
            _ => Err(OutputError::UnsupportedFormat(format!(
                "unknown output format {value:?}; expected text, tsv, json, arrow, or mex"
            ))),
        }
    }
}

#[derive(Debug)]
pub enum OutputError {
    InvalidSchema(String),
    InvalidRow(String),
    InvalidMatrix(String),
    UnsupportedFormat(String),
    Sink(String),
    Io(std::io::Error),
    Json(serde_json::Error),
}

impl fmt::Display for OutputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSchema(message) => write!(f, "invalid output schema: {message}"),
            Self::InvalidRow(message) => write!(f, "invalid output row: {message}"),
            Self::InvalidMatrix(message) => write!(f, "invalid sparse matrix: {message}"),
            Self::UnsupportedFormat(message) => write!(f, "unsupported output format: {message}"),
            Self::Sink(message) => write!(f, "output sink failed: {message}"),
            Self::Io(error) => error.fmt(f),
            Self::Json(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for OutputError {}

impl From<std::io::Error> for OutputError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<serde_json::Error> for OutputError {
    fn from(error: serde_json::Error) -> Self {
        Self::Json(error)
    }
}

pub(crate) fn validate_row(
    schema: &TableSchema,
    index: usize,
    row: &[ScalarValue],
) -> Result<(), OutputError> {
    if row.len() != schema.fields.len() {
        return Err(OutputError::InvalidRow(format!(
            "row {index} has {} values for {} fields",
            row.len(),
            schema.fields.len()
        )));
    }
    for (column, (value, field)) in row.iter().zip(&schema.fields).enumerate() {
        if !value.compatible_with(field) {
            return Err(OutputError::InvalidRow(format!(
                "row {index}, column {column} ({}) does not match {}{}",
                field.name,
                field.data_type.name(),
                if field.nullable { " or null" } else { "" }
            )));
        }
    }
    Ok(())
}

/// Write a table from an iterator, validating each row before it is emitted. JSON rows are arrays
/// whose meaning is fixed by `data.schema.fields`; this keeps the stream compact and typed.
pub fn write_table<W, I>(
    mut writer: W,
    schema: &TableSchema,
    rows: I,
    format: OutputFormat,
    context: &ResultContext,
) -> Result<(), OutputError>
where
    W: Write,
    I: IntoIterator<Item = Vec<ScalarValue>>,
{
    schema.validate()?;
    context.validate()?;
    match format {
        OutputFormat::Text => write_delimited(&mut writer, schema, rows, context, true),
        OutputFormat::Tsv => write_delimited(&mut writer, schema, rows, context, false),
        OutputFormat::Json => write_json(&mut writer, schema, rows, context),
        OutputFormat::Arrow => Err(OutputError::UnsupportedFormat(
            "Arrow IPC requires an ArrowBatchSink adapter; use stream_arrow".into(),
        )),
        OutputFormat::Mex => Err(OutputError::UnsupportedFormat(
            "MEX requires sparse matrix dimensions and labels; use write_mex_new".into(),
        )),
    }
}

fn write_delimited<W, I>(
    writer: &mut W,
    schema: &TableSchema,
    rows: I,
    context: &ResultContext,
    human: bool,
) -> Result<(), OutputError>
where
    W: Write,
    I: IntoIterator<Item = Vec<ScalarValue>>,
{
    if human {
        writeln!(writer, "result: {}", schema.id)?;
        writeln!(
            writer,
            "producer: {}@{}",
            context.producer.name, context.producer.version
        )?;
        for archive in &context.provenance.archives {
            writeln!(writer, "archive: {archive}")?;
        }
        if let Some(assembly) = &context.provenance.assembly {
            writeln!(writer, "assembly: {assembly}")?;
        }
        if let Some(annotation) = &context.provenance.annotation {
            writeln!(writer, "annotation: {annotation}")?;
        }
        if let Some(digest) = &context.provenance.annotation_digest {
            writeln!(writer, "annotation digest: {digest}")?;
        }
        for annotation in &context.provenance.annotations {
            writeln!(
                writer,
                "annotation {}: {} on {} ({})",
                annotation.role, annotation.annotation, annotation.assembly, annotation.digest
            )?;
        }
        for (name, value) in &context.provenance.parameters {
            writeln!(
                writer,
                "parameter {name}: {}",
                serde_json::to_string(value)?
            )?;
        }
        for warning in &context.warnings {
            writeln!(writer, "warning: {warning}")?;
        }
        for (column, field) in schema.fields.iter().enumerate() {
            if column > 0 {
                write!(writer, " | ")?;
            }
            write!(writer, "{} ({})", field.name, field.data_type.name())?;
        }
        writeln!(writer)?;
    } else {
        writeln!(writer, "# envelope_schema={ENVELOPE_SCHEMA}")?;
        writeln!(writer, "# result_schema={}", schema.id)?;
        writeln!(
            writer,
            "# producer={}@{}",
            context.producer.name, context.producer.version
        )?;
        for archive in &context.provenance.archives {
            writeln!(writer, "# archive={}", escape_tsv(archive))?;
        }
        if let Some(assembly) = &context.provenance.assembly {
            writeln!(writer, "# assembly={}", escape_tsv(assembly))?;
        }
        if let Some(annotation) = &context.provenance.annotation {
            writeln!(writer, "# annotation={}", escape_tsv(annotation))?;
        }
        if let Some(digest) = &context.provenance.annotation_digest {
            writeln!(writer, "# annotation_digest={}", escape_tsv(digest))?;
        }
        for annotation in &context.provenance.annotations {
            writeln!(
                writer,
                "# annotation_input={}",
                escape_tsv(&serde_json::to_string(annotation)?)
            )?;
        }
        for (name, value) in &context.provenance.parameters {
            writeln!(
                writer,
                "# parameter.{}={}",
                escape_tsv(name),
                escape_tsv(&serde_json::to_string(value)?)
            )?;
        }
        for warning in &context.warnings {
            writeln!(writer, "# warning={}", escape_tsv(warning))?;
        }
        for (index, field) in schema.fields.iter().enumerate() {
            if index > 0 {
                writer.write_all(b"\t")?;
            }
            writer.write_all(field.name.as_bytes())?;
        }
        writeln!(writer)?;
    }

    for (row_index, row) in rows.into_iter().enumerate() {
        validate_row(schema, row_index, &row)?;
        for (column, value) in row.iter().enumerate() {
            if column > 0 {
                writer.write_all(if human { b" | " } else { b"\t" })?;
            }
            let rendered = if human {
                render_text(value)?
            } else {
                render_tsv(value)?
            };
            writer.write_all(rendered.as_bytes())?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

fn write_json<W, I>(
    writer: &mut W,
    schema: &TableSchema,
    rows: I,
    context: &ResultContext,
) -> Result<(), OutputError>
where
    W: Write,
    I: IntoIterator<Item = Vec<ScalarValue>>,
{
    writer.write_all(b"{\"$schema\":")?;
    serde_json::to_writer(&mut *writer, ENVELOPE_SCHEMA)?;
    writer.write_all(b",\"result_schema\":")?;
    serde_json::to_writer(&mut *writer, &schema.id)?;
    writer.write_all(b",\"producer\":")?;
    serde_json::to_writer(&mut *writer, &context.producer)?;
    writer.write_all(b",\"provenance\":")?;
    serde_json::to_writer(&mut *writer, &context.provenance)?;
    writer.write_all(b",\"warnings\":")?;
    serde_json::to_writer(&mut *writer, &context.warnings)?;
    writer.write_all(b",\"data\":{\"schema\":")?;
    serde_json::to_writer(&mut *writer, schema)?;
    writer.write_all(b",\"rows\":[")?;
    let mut first = true;
    for (row_index, row) in rows.into_iter().enumerate() {
        validate_row(schema, row_index, &row)?;
        if !first {
            writer.write_all(b",")?;
        }
        serde_json::to_writer(&mut *writer, &row)?;
        first = false;
    }
    writer.write_all(b"]}}\n")?;
    Ok(())
}

fn render_text(value: &ScalarValue) -> Result<String, OutputError> {
    match value {
        ScalarValue::Null => Ok("—".into()),
        ScalarValue::String(value) => Ok(value.replace('|', "\\|").replace('\n', "\\n")),
        ScalarValue::Int64(value) => Ok(value.to_string()),
        ScalarValue::UInt64(value) => Ok(value.to_string()),
        ScalarValue::Float64(value) => Ok(value.to_string()),
        ScalarValue::Boolean(value) => Ok(value.to_string()),
        ScalarValue::Json(value) => Ok(serde_json::to_string(value)?),
    }
}

fn render_tsv(value: &ScalarValue) -> Result<String, OutputError> {
    match value {
        ScalarValue::Null => Ok("\\N".into()),
        ScalarValue::String(value) => Ok(escape_tsv(value)),
        ScalarValue::Int64(value) => Ok(value.to_string()),
        ScalarValue::UInt64(value) => Ok(value.to_string()),
        ScalarValue::Float64(value) => Ok(value.to_string()),
        ScalarValue::Boolean(value) => Ok(value.to_string()),
        ScalarValue::Json(value) => Ok(escape_tsv(&serde_json::to_string(value)?)),
    }
}

pub(crate) fn escape_tsv(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Producer, Provenance};

    fn schema() -> TableSchema {
        TableSchema::new(
            "gravlax.test.rows.v1",
            vec![
                Field::new("gene", DataType::String),
                Field::new("umis", DataType::UInt64),
                Field::new("score", DataType::Float64).nullable(),
            ],
        )
        .unwrap()
    }

    fn context() -> ResultContext {
        ResultContext {
            producer: Producer {
                name: "aie".into(),
                version: "test".into(),
            },
            provenance: Provenance {
                assembly: Some("GRCh38".into()),
                annotation: Some("GENCODE 49".into()),
                ..Default::default()
            },
            warnings: vec![],
        }
    }

    fn rows() -> Vec<Vec<ScalarValue>> {
        vec![vec!["TP53\talias".into(), 9_u64.into(), ScalarValue::Null]]
    }

    #[test]
    fn streams_valid_json_envelope_with_typed_schema() {
        let mut output = Vec::new();
        write_table(
            &mut output,
            &schema(),
            rows(),
            OutputFormat::Json,
            &context(),
        )
        .unwrap();
        let value: serde_json::Value = serde_json::from_slice(&output).unwrap();
        assert_eq!(value["$schema"], ENVELOPE_SCHEMA);
        assert_eq!(value["result_schema"], "gravlax.test.rows.v1");
        assert_eq!(value["provenance"]["assembly"], "GRCh38");
        assert_eq!(value["data"]["rows"][0][1], 9);
    }

    #[test]
    fn tsv_has_contract_metadata_header_and_reversible_escapes() {
        let mut output = Vec::new();
        write_table(
            &mut output,
            &schema(),
            rows(),
            OutputFormat::Tsv,
            &context(),
        )
        .unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(output.starts_with("# envelope_schema=gravlax.result-envelope.v1\n"));
        assert!(output.contains("gene\tumis\tscore\n"));
        assert!(output.contains("TP53\\talias\t9\t\\N\n"));
    }

    #[test]
    fn rejects_type_mismatch_before_emitting_that_row() {
        let bad = vec![vec![
            "TP53".into(),
            ScalarValue::Int64(9),
            ScalarValue::Null,
        ]];
        let error =
            write_table(Vec::new(), &schema(), bad, OutputFormat::Json, &context()).unwrap_err();
        assert!(error.to_string().contains("column 1 (umis)"));
    }

    #[test]
    fn format_names_are_uniform() {
        assert_eq!(
            "arrow-ipc".parse::<OutputFormat>().unwrap(),
            OutputFormat::Arrow
        );
        assert_eq!("mtx".parse::<OutputFormat>().unwrap(), OutputFormat::Mex);
        assert!("csv".parse::<OutputFormat>().is_err());
    }

    #[test]
    fn typed_table_round_trip_preserves_logical_scalar_types() {
        let schema = TableSchema::new(
            "gravlax.test.roundtrip.v1",
            vec![
                Field::new("small_unsigned", DataType::UInt64),
                Field::new("large_unsigned", DataType::UInt64),
                Field::new("positive_signed", DataType::Int64),
                Field::new("fraction", DataType::Float64),
                Field::new("optional", DataType::String).nullable(),
                Field::new("payload", DataType::Json),
            ],
        )
        .unwrap();
        let table = TypedTable::new(
            schema,
            vec![vec![
                ScalarValue::UInt64(1),
                ScalarValue::UInt64(u64::MAX),
                ScalarValue::Int64(7),
                ScalarValue::Float64(1.0),
                ScalarValue::Null,
                ScalarValue::Json(serde_json::json!({"ok": true})),
            ]],
        )
        .unwrap();
        let encoded = serde_json::to_vec(&table).unwrap();
        let decoded: TypedTable = serde_json::from_slice(&encoded).unwrap();
        assert_eq!(decoded, table);
    }

    #[test]
    fn typed_table_deserialization_rejects_invalid_rows_and_schema() {
        let wrong_type = serde_json::json!({
            "schema": {
                "id": "gravlax.test.decode.v1",
                "fields": [{"name": "count", "data_type": "uint64"}]
            },
            "rows": [[-1]]
        });
        let error = serde_json::from_value::<TypedTable>(wrong_type).unwrap_err();
        assert!(error.to_string().contains("column 0 (count)"), "{error}");

        let wrong_width = serde_json::json!({
            "schema": {
                "id": "gravlax.test.decode.v1",
                "fields": [{"name": "count", "data_type": "uint64"}]
            },
            "rows": [[]]
        });
        assert!(serde_json::from_value::<TypedTable>(wrong_width).is_err());

        let invalid_schema = serde_json::json!({
            "schema": {"id": "", "fields": []},
            "rows": []
        });
        assert!(serde_json::from_value::<TypedTable>(invalid_schema).is_err());

        let noncanonical_type = serde_json::json!({
            "schema": {
                "id": "gravlax.test.decode.v1",
                "fields": [{"name": "count", "data_type": "u_int64"}]
            },
            "rows": [[1]]
        });
        assert!(serde_json::from_value::<TypedTable>(noncanonical_type).is_err());

        let invalid_field = serde_json::json!({
            "schema": {
                "id": "gravlax.test.decode.v1",
                "fields": [{"name": "bad\tname", "data_type": "uint64"}]
            },
            "rows": [[1]]
        });
        assert!(serde_json::from_value::<TypedTable>(invalid_field).is_err());
    }

    #[test]
    fn integer_json_numbers_are_valid_float64_cells() {
        let decoded: TypedTable = serde_json::from_value(serde_json::json!({
            "schema": {
                "id": "gravlax.test.float.v1",
                "fields": [{"name": "score", "data_type": "float64"}]
            },
            "rows": [[1]]
        }))
        .unwrap();
        assert_eq!(decoded.rows, vec![vec![ScalarValue::Float64(1.0)]]);
    }

    #[test]
    fn json_null_has_only_the_table_null_representation() {
        let schema = TableSchema::new(
            "gravlax.test.json-null.v1",
            vec![Field::new("payload", DataType::Json).nullable()],
        )
        .unwrap();
        assert!(TypedTable::new(
            schema.clone(),
            vec![vec![ScalarValue::Json(serde_json::Value::Null)]],
        )
        .is_err());
        let decoded: TypedTable = serde_json::from_value(serde_json::json!({
            "schema": schema,
            "rows": [[null]]
        }))
        .unwrap();
        assert_eq!(decoded.rows, vec![vec![ScalarValue::Null]]);
    }

    #[test]
    fn legacy_schema_wire_shape_remains_unchanged_without_semantics() {
        let schema = TableSchema::new(
            "gravlax.test.legacy-schema.v1",
            vec![Field::new("id", DataType::String)],
        )
        .unwrap();
        assert_eq!(
            serde_json::to_value(&schema).unwrap(),
            serde_json::json!({
                "id": "gravlax.test.legacy-schema.v1",
                "fields": [{"name": "id", "data_type": "string", "nullable": false}]
            })
        );
        let decoded: TableSchema = serde_json::from_value(serde_json::json!({
            "id": "gravlax.test.legacy-schema.v1",
            "fields": [{"name": "id", "data_type": "string"}]
        }))
        .unwrap();
        assert_eq!(decoded.semantics, None);
        decoded.validate().unwrap();
    }

    #[test]
    fn semantic_keys_must_reference_distinct_schema_fields() {
        let base = || {
            TableSchema::new(
                "gravlax.test.semantics.v1",
                vec![
                    Field::new("id", DataType::String),
                    Field::new("count", DataType::UInt64),
                ],
            )
            .unwrap()
        };
        base()
            .with_semantics(
                TableSemantics::new(RowSemantics::Set)
                    .with_key(["id"])
                    .ordered_by([OrderKey::descending("count")]),
            )
            .unwrap();

        let error = base()
            .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["missing"]))
            .unwrap_err();
        assert!(error.to_string().contains("unknown field \"missing\""));

        let error = base()
            .with_semantics(TableSemantics::new(RowSemantics::Set).with_key(["id", "id"]))
            .unwrap_err();
        assert!(error.to_string().contains("repeats field \"id\""));

        let error = base()
            .with_semantics(
                TableSemantics::new(RowSemantics::Set).ordered_by([OrderKey::ascending("missing")]),
            )
            .unwrap_err();
        assert!(error.to_string().contains("unknown field \"missing\""));
    }
}
