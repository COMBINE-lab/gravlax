use crate::table::{escape_tsv, valid_schema_id};
use crate::{
    DataType, Field, OutputError, OutputFormat, ResultContext, RowSemantics, ScalarValue,
    TableSchema, ENVELOPE_SCHEMA,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::io::Write;

/// A borrowed, typed cell used by the allocation-light output path. Numeric and Boolean values are
/// copied; strings and JSON are serialized directly from their source values.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CellValueRef<'a> {
    Null,
    String(&'a str),
    Int64(i64),
    UInt64(u64),
    Float64(f64),
    Boolean(bool),
    Json(&'a serde_json::Value),
}

impl<'a> From<&'a str> for CellValueRef<'a> {
    fn from(value: &'a str) -> Self {
        Self::String(value)
    }
}

impl From<i64> for CellValueRef<'_> {
    fn from(value: i64) -> Self {
        Self::Int64(value)
    }
}

impl From<u64> for CellValueRef<'_> {
    fn from(value: u64) -> Self {
        Self::UInt64(value)
    }
}

impl From<f64> for CellValueRef<'_> {
    fn from(value: f64) -> Self {
        Self::Float64(value)
    }
}

impl From<bool> for CellValueRef<'_> {
    fn from(value: bool) -> Self {
        Self::Boolean(value)
    }
}

impl<'a> From<&'a ScalarValue> for CellValueRef<'a> {
    fn from(value: &'a ScalarValue) -> Self {
        match value {
            ScalarValue::Null => Self::Null,
            ScalarValue::String(value) => Self::String(value),
            ScalarValue::Int64(value) => Self::Int64(*value),
            ScalarValue::UInt64(value) => Self::UInt64(*value),
            ScalarValue::Float64(value) => Self::Float64(*value),
            ScalarValue::Boolean(value) => Self::Boolean(*value),
            ScalarValue::Json(value) => Self::Json(value),
        }
    }
}

impl CellValueRef<'_> {
    fn compatible_with(self, field: &Field) -> bool {
        match self {
            Self::Null => field.nullable,
            Self::String(_) => field.data_type == DataType::String,
            Self::Int64(_) => field.data_type == DataType::Int64,
            Self::UInt64(_) => field.data_type == DataType::UInt64,
            Self::Float64(value) => field.data_type == DataType::Float64 && value.is_finite(),
            Self::Boolean(_) => field.data_type == DataType::Boolean,
            Self::Json(value) => field.data_type == DataType::Json && !value.is_null(),
        }
    }
}

/// Describes an explicitly bounded selection. This is result-instance metadata rather than schema
/// metadata. Unbounded complete tables should omit it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SelectionSummary {
    pub available_rows: u64,
    pub emitted_rows: u64,
    pub truncated: bool,
}

impl SelectionSummary {
    pub fn complete(rows: u64) -> Self {
        Self {
            available_rows: rows,
            emitted_rows: rows,
            truncated: false,
        }
    }

    pub fn selected(available_rows: u64, emitted_rows: u64) -> Result<Self, OutputError> {
        let summary = Self {
            available_rows,
            emitted_rows,
            truncated: emitted_rows < available_rows,
        };
        summary.validate()?;
        Ok(summary)
    }

    pub fn validate(self) -> Result<(), OutputError> {
        if self.emitted_rows > self.available_rows {
            return Err(OutputError::InvalidSchema(format!(
                "selection emits {} rows although only {} are available",
                self.emitted_rows, self.available_rows
            )));
        }
        if self.truncated != (self.emitted_rows < self.available_rows) {
            return Err(OutputError::InvalidSchema(
                "selection truncated must equal emitted_rows < available_rows".into(),
            ));
        }
        Ok(())
    }
}

/// Finish-time knowledge for a capped, genuinely one-pass producer. Unlike [`SelectionSummary`],
/// this does not require the total number of available rows before the table header is written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SelectionCompletion {
    /// The producer counted all available rows. Truncation is derived from the emitted row count.
    Exact { available_rows: u64 },
    /// The producer stopped without counting the full input. `truncated` is `Some(true)` when it
    /// observed at least one additional row, and `None` when it did not peek beyond the cap.
    AvailabilityUnknown { truncated: Option<bool> },
}

impl SelectionCompletion {
    pub fn exact(available_rows: u64) -> Self {
        Self::Exact { available_rows }
    }

    pub fn availability_unknown() -> Self {
        Self::AvailabilityUnknown { truncated: None }
    }

    pub fn availability_unknown_truncated() -> Self {
        Self::AvailabilityUnknown {
            truncated: Some(true),
        }
    }

    fn finalize(self, emitted_rows: u64) -> Result<SelectionOutcome, OutputError> {
        match self {
            Self::Exact { available_rows } => {
                if emitted_rows > available_rows {
                    return Err(OutputError::InvalidRow(format!(
                        "selection emitted {emitted_rows} rows although only {available_rows} are available"
                    )));
                }
                Ok(SelectionOutcome {
                    available_rows: Some(available_rows),
                    emitted_rows,
                    truncated: Some(emitted_rows < available_rows),
                })
            }
            Self::AvailabilityUnknown {
                truncated: Some(false),
            } => Err(OutputError::InvalidSchema(
                "an untruncated completed stream has exact availability; use SelectionCompletion::exact"
                    .into(),
            )),
            Self::AvailabilityUnknown { truncated } => Ok(SelectionOutcome {
                available_rows: None,
                emitted_rows,
                truncated,
            }),
        }
    }
}

/// The serialized selection metadata produced at stream completion. Null availability and
/// truncation mean that the producer deliberately did not scan beyond its output cap.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SelectionOutcome {
    pub available_rows: Option<u64>,
    pub emitted_rows: u64,
    pub truncated: Option<bool>,
}

/// Streams one typed table without constructing `Vec<ScalarValue>` rows. A schema and result
/// context are validated once at construction; each cell is then type-checked as it is encoded.
#[must_use = "a streaming table must be finished to close the result envelope"]
pub struct StreamingTableWriter<'schema, W> {
    writer: W,
    schema: &'schema TableSchema,
    format: OutputFormat,
    selection: Option<SelectionSummary>,
    row_count: u64,
}

impl<'schema, W> StreamingTableWriter<'schema, W>
where
    W: Write,
{
    pub fn new(
        mut writer: W,
        schema: &'schema TableSchema,
        format: OutputFormat,
        context: &ResultContext,
        selection: Option<&SelectionSummary>,
    ) -> Result<Self, OutputError> {
        schema.validate()?;
        context.validate()?;
        let selection = selection.copied();
        if let Some(selection) = selection {
            selection.validate()?;
        }
        reject_non_tabular_format(format)?;
        write_single_table_header(&mut writer, schema, format, context, selection.as_ref())?;
        Ok(Self {
            writer,
            schema,
            format,
            selection,
            row_count: 0,
        })
    }

    /// Write borrowed cells from a slice. No cell strings or row vector are cloned.
    pub fn write_row(&mut self, cells: &[CellValueRef<'_>]) -> Result<(), OutputError> {
        self.write_row_iter(cells.iter().copied())
    }

    /// Write borrowed cells from an iterator, including a stack-allocated array.
    pub fn write_row_iter<'cell, I>(&mut self, cells: I) -> Result<(), OutputError>
    where
        I: IntoIterator<Item = CellValueRef<'cell>>,
    {
        self.write_row_with(|row| {
            for cell in cells {
                row.value(cell)?;
            }
            Ok(())
        })
    }

    /// Encode a row directly from its source object. Validation and serialization are fused, and
    /// the closure may fail without manufacturing a placeholder row.
    pub fn write_row_with<F>(&mut self, encode: F) -> Result<(), OutputError>
    where
        F: FnOnce(&mut RowEncoder<'_, 'schema, W>) -> Result<(), OutputError>,
    {
        begin_row(&mut self.writer, self.format, self.row_count)?;
        let mut row = RowEncoder {
            writer: &mut self.writer,
            schema: self.schema,
            format: self.format,
            row_index: self.row_count,
            column: 0,
        };
        encode(&mut row)?;
        row.finish()?;
        self.row_count = self.row_count.checked_add(1).ok_or_else(|| {
            OutputError::InvalidRow("row count exceeds the uniform interface limit".into())
        })?;
        Ok(())
    }

    pub fn rows_written(&self) -> u64 {
        self.row_count
    }

    pub fn finish(mut self) -> Result<W, OutputError> {
        validate_emitted_rows(self.selection.as_ref(), self.row_count)?;
        finish_single_table(&mut self.writer, self.format)?;
        self.writer.flush()?;
        Ok(self.writer)
    }

    /// Finish with selection metadata that could not be known before the one-pass row producer
    /// ran. This is valid only when construction omitted the exact `SelectionSummary`.
    pub fn finish_with_selection(
        mut self,
        completion: SelectionCompletion,
    ) -> Result<W, OutputError> {
        if self.selection.is_some() {
            return Err(OutputError::InvalidSchema(
                "selection was supplied both at table start and finish".into(),
            ));
        }
        let selection = completion.finalize(self.row_count)?;
        finish_single_table_with_selection(&mut self.writer, self.format, &selection)?;
        self.writer.flush()?;
        Ok(self.writer)
    }
}

/// The per-row encoder exposed by [`StreamingTableWriter::write_row_with`] and bundle tables.
/// Methods return errors at the exact column that violates the schema.
pub struct RowEncoder<'writer, 'schema, W> {
    writer: &'writer mut W,
    schema: &'schema TableSchema,
    format: OutputFormat,
    row_index: u64,
    column: usize,
}

impl<W> RowEncoder<'_, '_, W>
where
    W: Write,
{
    pub fn value(&mut self, value: CellValueRef<'_>) -> Result<(), OutputError> {
        let Some(field) = self.schema.fields.get(self.column) else {
            return Err(OutputError::InvalidRow(format!(
                "row {} has more than {} values",
                self.row_index,
                self.schema.fields.len()
            )));
        };
        if !value.compatible_with(field) {
            return Err(type_mismatch(self.row_index, self.column, field));
        }
        if self.column > 0 {
            self.writer.write_all(match self.format {
                OutputFormat::Text => b" | ",
                OutputFormat::Tsv => b"\t",
                OutputFormat::Json => b",",
                OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
            })?;
        }
        write_cell(self.writer, self.format, value)?;
        self.column += 1;
        Ok(())
    }

    pub fn null(&mut self) -> Result<(), OutputError> {
        self.value(CellValueRef::Null)
    }

    pub fn string(&mut self, value: &str) -> Result<(), OutputError> {
        self.value(CellValueRef::String(value))
    }

    pub fn int64(&mut self, value: i64) -> Result<(), OutputError> {
        self.value(CellValueRef::Int64(value))
    }

    pub fn uint64(&mut self, value: u64) -> Result<(), OutputError> {
        self.value(CellValueRef::UInt64(value))
    }

    pub fn float64(&mut self, value: f64) -> Result<(), OutputError> {
        self.value(CellValueRef::Float64(value))
    }

    pub fn boolean(&mut self, value: bool) -> Result<(), OutputError> {
        self.value(CellValueRef::Boolean(value))
    }

    pub fn json(&mut self, value: &serde_json::Value) -> Result<(), OutputError> {
        self.value(CellValueRef::Json(value))
    }

    fn finish(&mut self) -> Result<(), OutputError> {
        if self.column != self.schema.fields.len() {
            return Err(OutputError::InvalidRow(format!(
                "row {} has {} values for {} fields",
                self.row_index,
                self.column,
                self.schema.fields.len()
            )));
        }
        match self.format {
            OutputFormat::Json => self.writer.write_all(b"]")?,
            OutputFormat::Text | OutputFormat::Tsv => self.writer.write_all(b"\n")?,
            OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
        }
        Ok(())
    }
}

/// Convenience boundary for a producer that wants construction, production, and envelope closure
/// to be one fallible operation.
pub fn write_streaming_table<'schema, W, F>(
    writer: W,
    schema: &'schema TableSchema,
    format: OutputFormat,
    context: &ResultContext,
    selection: Option<&SelectionSummary>,
    produce: F,
) -> Result<W, OutputError>
where
    W: Write,
    F: FnOnce(&mut StreamingTableWriter<'schema, W>) -> Result<(), OutputError>,
{
    let mut stream = StreamingTableWriter::new(writer, schema, format, context, selection)?;
    produce(&mut stream)?;
    stream.finish()
}

/// One-operation counterpart to [`write_streaming_table`] for producers that learn selection
/// completeness only while streaming.
pub fn write_streaming_table_with_deferred_selection<'schema, W, F>(
    writer: W,
    schema: &'schema TableSchema,
    format: OutputFormat,
    context: &ResultContext,
    produce: F,
) -> Result<W, OutputError>
where
    W: Write,
    F: FnOnce(&mut StreamingTableWriter<'schema, W>) -> Result<SelectionCompletion, OutputError>,
{
    let mut stream = StreamingTableWriter::new(writer, schema, format, context, None)?;
    let completion = produce(&mut stream)?;
    stream.finish_with_selection(completion)
}

/// Streams a named, multi-table result. Tables are completed sequentially, so memory use is
/// independent of the total number of rows. TSV and text use explicit table sections; JSON uses a
/// `data.tables` array. Table names must be unique.
#[must_use = "a streaming bundle must be finished to close the result envelope"]
pub struct StreamingBundleWriter<W> {
    writer: W,
    result_schema: String,
    format: OutputFormat,
    first_table: bool,
    table_names: BTreeSet<String>,
}

impl<W> StreamingBundleWriter<W>
where
    W: Write,
{
    pub fn new(
        mut writer: W,
        result_schema: impl Into<String>,
        format: OutputFormat,
        context: &ResultContext,
    ) -> Result<Self, OutputError> {
        let result_schema = result_schema.into();
        if !valid_schema_id(&result_schema) {
            return Err(OutputError::InvalidSchema(
                "bundle result schema id must be non-empty ASCII letters, digits, '.', '_' or '-'"
                    .into(),
            ));
        }
        context.validate()?;
        reject_non_tabular_format(format)?;
        write_bundle_header(&mut writer, &result_schema, format, context)?;
        Ok(Self {
            writer,
            result_schema,
            format,
            first_table: true,
            table_names: BTreeSet::new(),
        })
    }

    /// Begin a bundle whose `data` object contains a typed scientific summary in addition to its
    /// streamed tables. The summary is result data, not provenance, and is serialized only once.
    pub fn new_with_summary<S>(
        mut writer: W,
        result_schema: impl Into<String>,
        format: OutputFormat,
        context: &ResultContext,
        summary: &S,
    ) -> Result<Self, OutputError>
    where
        S: Serialize,
    {
        let result_schema = result_schema.into();
        if !valid_schema_id(&result_schema) {
            return Err(OutputError::InvalidSchema(
                "bundle result schema id must be non-empty ASCII letters, digits, '.', '_' or '-'"
                    .into(),
            ));
        }
        context.validate()?;
        reject_non_tabular_format(format)?;
        match format {
            OutputFormat::Json => {
                write_json_envelope_prefix(&mut writer, &result_schema, context)?;
                writer.write_all(b"{\"summary\":")?;
                serde_json::to_writer(&mut writer, summary)?;
                writer.write_all(b",\"tables\":[")?;
            }
            OutputFormat::Text => {
                write_text_context(&mut writer, &result_schema, context)?;
                writeln!(writer, "summary: {}", serde_json::to_string(summary)?)?;
            }
            OutputFormat::Tsv => {
                write_tsv_context(&mut writer, &result_schema, context)?;
                writeln!(
                    writer,
                    "# summary={}",
                    escape_tsv(&serde_json::to_string(summary)?)
                )?;
            }
            OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
        }
        Ok(Self {
            writer,
            result_schema,
            format,
            first_table: true,
            table_names: BTreeSet::new(),
        })
    }

    pub fn write_table<'schema, F>(
        &mut self,
        name: &str,
        schema: &'schema TableSchema,
        selection: Option<&SelectionSummary>,
        produce: F,
    ) -> Result<(), OutputError>
    where
        F: FnOnce(&mut StreamingBundleTableWriter<'_, 'schema, W>) -> Result<(), OutputError>,
    {
        if !valid_schema_id(name) {
            return Err(OutputError::InvalidSchema(format!(
                "invalid bundle table name {name:?}"
            )));
        }
        if !self.table_names.insert(name.to_owned()) {
            return Err(OutputError::InvalidSchema(format!(
                "duplicate bundle table name {name:?}"
            )));
        }
        schema.validate()?;
        let selection = selection.copied();
        if let Some(selection) = selection {
            selection.validate()?;
        }
        begin_bundle_table(
            &mut self.writer,
            self.format,
            self.first_table,
            name,
            schema,
            selection.as_ref(),
        )?;
        let mut table = StreamingBundleTableWriter {
            writer: &mut self.writer,
            schema,
            format: self.format,
            selection,
            row_count: 0,
        };
        produce(&mut table)?;
        table.finish()?;
        self.first_table = false;
        Ok(())
    }

    /// Stream a table whose availability/truncation metadata is learned only after its one-pass
    /// producer runs. The closure's completion value is written as a table footer.
    pub fn write_table_with_deferred_selection<'schema, F>(
        &mut self,
        name: &str,
        schema: &'schema TableSchema,
        produce: F,
    ) -> Result<(), OutputError>
    where
        F: FnOnce(
            &mut StreamingBundleTableWriter<'_, 'schema, W>,
        ) -> Result<SelectionCompletion, OutputError>,
    {
        if !valid_schema_id(name) {
            return Err(OutputError::InvalidSchema(format!(
                "invalid bundle table name {name:?}"
            )));
        }
        if !self.table_names.insert(name.to_owned()) {
            return Err(OutputError::InvalidSchema(format!(
                "duplicate bundle table name {name:?}"
            )));
        }
        schema.validate()?;
        begin_bundle_table(
            &mut self.writer,
            self.format,
            self.first_table,
            name,
            schema,
            None,
        )?;
        let mut table = StreamingBundleTableWriter {
            writer: &mut self.writer,
            schema,
            format: self.format,
            selection: None,
            row_count: 0,
        };
        let completion = produce(&mut table)?;
        table.finish_with_selection(completion)?;
        self.first_table = false;
        Ok(())
    }

    pub fn finish(mut self) -> Result<W, OutputError> {
        match self.format {
            OutputFormat::Json => self.writer.write_all(b"]}}\n")?,
            OutputFormat::Text | OutputFormat::Tsv => {}
            OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
        }
        self.writer.flush()?;
        Ok(self.writer)
    }

    pub fn result_schema(&self) -> &str {
        &self.result_schema
    }
}

pub struct StreamingBundleTableWriter<'writer, 'schema, W> {
    writer: &'writer mut W,
    schema: &'schema TableSchema,
    format: OutputFormat,
    selection: Option<SelectionSummary>,
    row_count: u64,
}

impl<'schema, W> StreamingBundleTableWriter<'_, 'schema, W>
where
    W: Write,
{
    pub fn write_row(&mut self, cells: &[CellValueRef<'_>]) -> Result<(), OutputError> {
        self.write_row_iter(cells.iter().copied())
    }

    pub fn write_row_iter<'cell, I>(&mut self, cells: I) -> Result<(), OutputError>
    where
        I: IntoIterator<Item = CellValueRef<'cell>>,
    {
        self.write_row_with(|row| {
            for cell in cells {
                row.value(cell)?;
            }
            Ok(())
        })
    }

    pub fn write_row_with<F>(&mut self, encode: F) -> Result<(), OutputError>
    where
        F: FnOnce(&mut RowEncoder<'_, 'schema, W>) -> Result<(), OutputError>,
    {
        begin_row(self.writer, self.format, self.row_count)?;
        let mut row = RowEncoder {
            writer: self.writer,
            schema: self.schema,
            format: self.format,
            row_index: self.row_count,
            column: 0,
        };
        encode(&mut row)?;
        row.finish()?;
        self.row_count = self.row_count.checked_add(1).ok_or_else(|| {
            OutputError::InvalidRow("row count exceeds the uniform interface limit".into())
        })?;
        Ok(())
    }

    pub fn rows_written(&self) -> u64 {
        self.row_count
    }

    fn finish(&mut self) -> Result<(), OutputError> {
        validate_emitted_rows(self.selection.as_ref(), self.row_count)?;
        match self.format {
            OutputFormat::Json => self.writer.write_all(b"]}")?,
            OutputFormat::Text => self.writer.write_all(b"\n")?,
            OutputFormat::Tsv => self.writer.write_all(b"# end_table\n")?,
            OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
        }
        Ok(())
    }

    fn finish_with_selection(
        &mut self,
        completion: SelectionCompletion,
    ) -> Result<(), OutputError> {
        if self.selection.is_some() {
            return Err(OutputError::InvalidSchema(
                "selection was supplied both at table start and finish".into(),
            ));
        }
        let selection = completion.finalize(self.row_count)?;
        match self.format {
            OutputFormat::Json => {
                self.writer.write_all(b"],\"selection\":")?;
                serde_json::to_writer(&mut *self.writer, &selection)?;
                self.writer.write_all(b"}")?;
            }
            OutputFormat::Text => {
                write_text_selection_outcome(self.writer, &selection)?;
                self.writer.write_all(b"\n")?;
            }
            OutputFormat::Tsv => {
                write_tsv_selection_outcome(self.writer, &selection)?;
                self.writer.write_all(b"# end_table\n")?;
            }
            OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
        }
        Ok(())
    }
}

fn reject_non_tabular_format(format: OutputFormat) -> Result<(), OutputError> {
    match format {
        OutputFormat::Text | OutputFormat::Tsv | OutputFormat::Json => Ok(()),
        OutputFormat::Arrow => Err(OutputError::UnsupportedFormat(
            "Arrow IPC requires an ArrowBatchSink adapter; use stream_arrow".into(),
        )),
        OutputFormat::Mex => Err(OutputError::UnsupportedFormat(
            "MEX requires sparse matrix dimensions and labels; use write_mex_new".into(),
        )),
    }
}

fn write_single_table_header<W: Write>(
    writer: &mut W,
    schema: &TableSchema,
    format: OutputFormat,
    context: &ResultContext,
    selection: Option<&SelectionSummary>,
) -> Result<(), OutputError> {
    match format {
        OutputFormat::Json => {
            write_json_envelope_prefix(writer, &schema.id, context)?;
            writer.write_all(b"{\"schema\":")?;
            serde_json::to_writer(&mut *writer, schema)?;
            if let Some(selection) = selection {
                writer.write_all(b",\"selection\":")?;
                serde_json::to_writer(&mut *writer, selection)?;
            }
            writer.write_all(b",\"rows\":[")?;
        }
        OutputFormat::Text => {
            write_text_context(writer, &schema.id, context)?;
            write_text_table_metadata(writer, schema, selection)?;
            write_field_header(writer, schema, true)?;
        }
        OutputFormat::Tsv => {
            write_tsv_context(writer, &schema.id, context)?;
            write_tsv_table_metadata(writer, schema, selection)?;
            write_field_header(writer, schema, false)?;
        }
        OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
    }
    Ok(())
}

fn finish_single_table<W: Write>(writer: &mut W, format: OutputFormat) -> Result<(), OutputError> {
    match format {
        OutputFormat::Json => writer.write_all(b"]}}\n")?,
        OutputFormat::Text | OutputFormat::Tsv => {}
        OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
    }
    Ok(())
}

fn finish_single_table_with_selection<W: Write>(
    writer: &mut W,
    format: OutputFormat,
    selection: &SelectionOutcome,
) -> Result<(), OutputError> {
    match format {
        OutputFormat::Json => {
            writer.write_all(b"],\"selection\":")?;
            serde_json::to_writer(&mut *writer, selection)?;
            writer.write_all(b"}}\n")?;
        }
        OutputFormat::Text => write_text_selection_outcome(writer, selection)?,
        OutputFormat::Tsv => write_tsv_selection_outcome(writer, selection)?,
        OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
    }
    Ok(())
}

fn write_bundle_header<W: Write>(
    writer: &mut W,
    result_schema: &str,
    format: OutputFormat,
    context: &ResultContext,
) -> Result<(), OutputError> {
    match format {
        OutputFormat::Json => {
            write_json_envelope_prefix(writer, result_schema, context)?;
            writer.write_all(b"{\"tables\":[")?;
        }
        OutputFormat::Text => write_text_context(writer, result_schema, context)?,
        OutputFormat::Tsv => write_tsv_context(writer, result_schema, context)?,
        OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
    }
    Ok(())
}

fn begin_bundle_table<W: Write>(
    writer: &mut W,
    format: OutputFormat,
    first: bool,
    name: &str,
    schema: &TableSchema,
    selection: Option<&SelectionSummary>,
) -> Result<(), OutputError> {
    match format {
        OutputFormat::Json => {
            if !first {
                writer.write_all(b",")?;
            }
            writer.write_all(b"{\"name\":")?;
            serde_json::to_writer(&mut *writer, name)?;
            writer.write_all(b",\"schema\":")?;
            serde_json::to_writer(&mut *writer, schema)?;
            if let Some(selection) = selection {
                writer.write_all(b",\"selection\":")?;
                serde_json::to_writer(&mut *writer, selection)?;
            }
            writer.write_all(b",\"rows\":[")?;
        }
        OutputFormat::Text => {
            if !first {
                writer.write_all(b"\n")?;
            }
            writeln!(writer, "table: {name}")?;
            write_text_table_metadata(writer, schema, selection)?;
            write_field_header(writer, schema, true)?;
        }
        OutputFormat::Tsv => {
            if !first {
                writer.write_all(b"\n")?;
            }
            writeln!(writer, "# table={}", escape_tsv(name))?;
            write_tsv_table_metadata(writer, schema, selection)?;
            write_field_header(writer, schema, false)?;
        }
        OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
    }
    Ok(())
}

fn begin_row<W: Write>(
    writer: &mut W,
    format: OutputFormat,
    row_index: u64,
) -> Result<(), OutputError> {
    if format == OutputFormat::Json {
        if row_index > 0 {
            writer.write_all(b",")?;
        }
        writer.write_all(b"[")?;
    }
    Ok(())
}

fn type_mismatch(row: u64, column: usize, field: &Field) -> OutputError {
    OutputError::InvalidRow(format!(
        "row {row}, column {column} ({}) does not match {}{}",
        field.name,
        field.data_type.name(),
        if field.nullable { " or null" } else { "" }
    ))
}

fn validate_emitted_rows(
    selection: Option<&SelectionSummary>,
    actual: u64,
) -> Result<(), OutputError> {
    if let Some(selection) = selection {
        if selection.emitted_rows != actual {
            return Err(OutputError::InvalidRow(format!(
                "selection declares {} emitted rows but producer wrote {actual}",
                selection.emitted_rows
            )));
        }
    }
    Ok(())
}

fn write_cell<W: Write>(
    writer: &mut W,
    format: OutputFormat,
    value: CellValueRef<'_>,
) -> Result<(), OutputError> {
    match format {
        OutputFormat::Json => match value {
            CellValueRef::Null => writer.write_all(b"null")?,
            CellValueRef::String(value) => serde_json::to_writer(&mut *writer, value)?,
            CellValueRef::Int64(value) => serde_json::to_writer(&mut *writer, &value)?,
            CellValueRef::UInt64(value) => serde_json::to_writer(&mut *writer, &value)?,
            CellValueRef::Float64(value) => serde_json::to_writer(&mut *writer, &value)?,
            CellValueRef::Boolean(value) => serde_json::to_writer(&mut *writer, &value)?,
            CellValueRef::Json(value) => serde_json::to_writer(&mut *writer, value)?,
        },
        OutputFormat::Tsv => match value {
            CellValueRef::Null => writer.write_all(b"\\N")?,
            CellValueRef::String(value) => write_escaped_tsv(writer, value)?,
            CellValueRef::Int64(value) => write!(writer, "{value}")?,
            CellValueRef::UInt64(value) => write!(writer, "{value}")?,
            CellValueRef::Float64(value) => write!(writer, "{value}")?,
            CellValueRef::Boolean(value) => write!(writer, "{value}")?,
            CellValueRef::Json(value) => {
                let encoded = serde_json::to_string(value)?;
                write_escaped_tsv(writer, &encoded)?;
            }
        },
        OutputFormat::Text => match value {
            CellValueRef::Null => writer.write_all("—".as_bytes())?,
            CellValueRef::String(value) => write_escaped_text(writer, value)?,
            CellValueRef::Int64(value) => write!(writer, "{value}")?,
            CellValueRef::UInt64(value) => write!(writer, "{value}")?,
            CellValueRef::Float64(value) => write!(writer, "{value}")?,
            CellValueRef::Boolean(value) => write!(writer, "{value}")?,
            CellValueRef::Json(value) => serde_json::to_writer(&mut *writer, value)?,
        },
        OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
    }
    Ok(())
}

fn write_escaped_tsv<W: Write>(writer: &mut W, value: &str) -> std::io::Result<()> {
    write_escaped(writer, value, |byte| match byte {
        b'\\' => Some(b"\\\\"),
        b'\t' => Some(b"\\t"),
        b'\r' => Some(b"\\r"),
        b'\n' => Some(b"\\n"),
        _ => None,
    })
}

fn write_escaped_text<W: Write>(writer: &mut W, value: &str) -> std::io::Result<()> {
    write_escaped(writer, value, |byte| match byte {
        b'|' => Some(b"\\|"),
        b'\n' => Some(b"\\n"),
        _ => None,
    })
}

fn write_escaped<W, F>(writer: &mut W, value: &str, replacement: F) -> std::io::Result<()>
where
    W: Write,
    F: Fn(u8) -> Option<&'static [u8]>,
{
    let bytes = value.as_bytes();
    let mut start = 0;
    for (index, byte) in bytes.iter().copied().enumerate() {
        if let Some(replacement) = replacement(byte) {
            writer.write_all(&bytes[start..index])?;
            writer.write_all(replacement)?;
            start = index + 1;
        }
    }
    writer.write_all(&bytes[start..])
}

fn write_json_envelope_prefix<W: Write>(
    writer: &mut W,
    result_schema: &str,
    context: &ResultContext,
) -> Result<(), OutputError> {
    writer.write_all(b"{\"$schema\":")?;
    serde_json::to_writer(&mut *writer, ENVELOPE_SCHEMA)?;
    writer.write_all(b",\"result_schema\":")?;
    serde_json::to_writer(&mut *writer, result_schema)?;
    writer.write_all(b",\"producer\":")?;
    serde_json::to_writer(&mut *writer, &context.producer)?;
    writer.write_all(b",\"provenance\":")?;
    serde_json::to_writer(&mut *writer, &context.provenance)?;
    writer.write_all(b",\"warnings\":")?;
    serde_json::to_writer(&mut *writer, &context.warnings)?;
    writer.write_all(b",\"data\":")?;
    Ok(())
}

fn write_text_context<W: Write>(
    writer: &mut W,
    result_schema: &str,
    context: &ResultContext,
) -> Result<(), OutputError> {
    writeln!(writer, "result: {result_schema}")?;
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
    Ok(())
}

fn write_tsv_context<W: Write>(
    writer: &mut W,
    result_schema: &str,
    context: &ResultContext,
) -> Result<(), OutputError> {
    writeln!(writer, "# envelope_schema={ENVELOPE_SCHEMA}")?;
    writeln!(writer, "# result_schema={result_schema}")?;
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
    Ok(())
}

fn write_text_table_metadata<W: Write>(
    writer: &mut W,
    schema: &TableSchema,
    selection: Option<&SelectionSummary>,
) -> Result<(), OutputError> {
    if let Some(semantics) = &schema.semantics {
        writeln!(
            writer,
            "row semantics: {}",
            row_semantics_name(semantics.row_semantics)
        )?;
        if let Some(key) = &semantics.key {
            writeln!(writer, "key: {}", key.join(", "))?;
        }
        match &semantics.ordered_by {
            Some(order) => writeln!(writer, "ordered by: {}", format_order(order))?,
            None => writeln!(writer, "ordering: unspecified")?,
        }
    }
    if let Some(selection) = selection {
        writeln!(writer, "available rows: {}", selection.available_rows)?;
        writeln!(writer, "emitted rows: {}", selection.emitted_rows)?;
        writeln!(writer, "truncated: {}", selection.truncated)?;
    }
    Ok(())
}

fn write_tsv_table_metadata<W: Write>(
    writer: &mut W,
    schema: &TableSchema,
    selection: Option<&SelectionSummary>,
) -> Result<(), OutputError> {
    writeln!(writer, "# table_schema={}", schema.id)?;
    if let Some(semantics) = &schema.semantics {
        writeln!(
            writer,
            "# table_semantics={}",
            escape_tsv(&serde_json::to_string(semantics)?)
        )?;
    }
    if let Some(selection) = selection {
        writeln!(writer, "# available_rows={}", selection.available_rows)?;
        writeln!(writer, "# emitted_rows={}", selection.emitted_rows)?;
        writeln!(writer, "# truncated={}", selection.truncated)?;
    }
    Ok(())
}

fn write_text_selection_outcome<W: Write>(
    writer: &mut W,
    selection: &SelectionOutcome,
) -> Result<(), OutputError> {
    writeln!(
        writer,
        "available rows: {}",
        selection
            .available_rows
            .map_or_else(|| "unknown".into(), |rows| rows.to_string())
    )?;
    writeln!(writer, "emitted rows: {}", selection.emitted_rows)?;
    writeln!(
        writer,
        "truncated: {}",
        selection
            .truncated
            .map_or("unknown", |truncated| if truncated {
                "true"
            } else {
                "false"
            })
    )?;
    Ok(())
}

fn write_tsv_selection_outcome<W: Write>(
    writer: &mut W,
    selection: &SelectionOutcome,
) -> Result<(), OutputError> {
    match selection.available_rows {
        Some(rows) => writeln!(writer, "# available_rows={rows}")?,
        None => writeln!(writer, "# available_rows=unknown")?,
    }
    writeln!(writer, "# emitted_rows={}", selection.emitted_rows)?;
    match selection.truncated {
        Some(truncated) => writeln!(writer, "# truncated={truncated}")?,
        None => writeln!(writer, "# truncated=unknown")?,
    }
    Ok(())
}

fn write_field_header<W: Write>(
    writer: &mut W,
    schema: &TableSchema,
    human: bool,
) -> Result<(), OutputError> {
    for (index, field) in schema.fields.iter().enumerate() {
        if index > 0 {
            writer.write_all(if human { b" | " } else { b"\t" })?;
        }
        writer.write_all(field.name.as_bytes())?;
        if human {
            write!(writer, " ({})", field.data_type.name())?;
        }
    }
    writer.write_all(b"\n")?;
    Ok(())
}

fn row_semantics_name(semantics: RowSemantics) -> &'static str {
    match semantics {
        RowSemantics::Set => "set",
        RowSemantics::Multiset => "multiset",
        RowSemantics::Sequence => "sequence",
    }
}

fn format_order(order: &[crate::OrderKey]) -> String {
    order
        .iter()
        .map(|key| {
            format!(
                "{} {}",
                key.field,
                match key.direction {
                    crate::SortDirection::Ascending => "ascending",
                    crate::SortDirection::Descending => "descending",
                }
            )
        })
        .collect::<Vec<_>>()
        .join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Field, OrderKey, Producer, Provenance, SortDirection, TableSemantics};
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;

    fn schema() -> TableSchema {
        TableSchema::new(
            "gravlax.test.streaming.v1",
            vec![
                Field::new("gene", DataType::String),
                Field::new("umis", DataType::UInt64),
                Field::new("score", DataType::Float64).nullable(),
            ],
        )
        .unwrap()
        .with_semantics(
            TableSemantics::new(RowSemantics::Set)
                .with_key(["gene"])
                .ordered_by([OrderKey {
                    field: "umis".into(),
                    direction: SortDirection::Descending,
                }]),
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
                ..Default::default()
            },
            warnings: vec![],
        }
    }

    #[test]
    fn borrowed_stream_is_a_typed_json_envelope() {
        let selection = SelectionSummary::complete(2);
        let output = write_streaming_table(
            Vec::new(),
            &schema(),
            OutputFormat::Json,
            &context(),
            Some(&selection),
            |rows| {
                rows.write_row_iter([
                    CellValueRef::String("TP53"),
                    CellValueRef::UInt64(9),
                    CellValueRef::Null,
                ])?;
                rows.write_row_with(|row| {
                    row.string("BRCA1")?;
                    row.uint64(7)?;
                    row.float64(1.5)
                })
            },
        )
        .unwrap();
        let value: serde_json::Value = serde_json::from_slice(&output).unwrap();
        assert_eq!(value["data"]["rows"][0][0], "TP53");
        assert_eq!(value["data"]["rows"][1][2], 1.5);
        assert_eq!(value["data"]["selection"]["emitted_rows"], 2);
        assert_eq!(value["data"]["schema"]["semantics"]["row_semantics"], "set");
        assert_eq!(
            value["data"]["schema"]["semantics"]["ordered_by"][0]["direction"],
            "descending"
        );
    }

    #[test]
    fn stream_rejects_types_width_nonfinite_and_false_selection_claims() {
        let error = write_streaming_table(
            Vec::new(),
            &schema(),
            OutputFormat::Json,
            &context(),
            None,
            |rows| {
                rows.write_row_iter([
                    CellValueRef::String("TP53"),
                    CellValueRef::Int64(9),
                    CellValueRef::Null,
                ])
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains("column 1 (umis)"));

        let error = write_streaming_table(
            Vec::new(),
            &schema(),
            OutputFormat::Tsv,
            &context(),
            None,
            |rows| rows.write_row_iter([CellValueRef::String("TP53"), CellValueRef::UInt64(9)]),
        )
        .unwrap_err();
        assert!(error.to_string().contains("2 values for 3 fields"));

        let selection = SelectionSummary::complete(2);
        let error = write_streaming_table(
            Vec::new(),
            &schema(),
            OutputFormat::Json,
            &context(),
            Some(&selection),
            |rows| {
                rows.write_row_iter([
                    CellValueRef::String("TP53"),
                    CellValueRef::UInt64(9),
                    CellValueRef::Null,
                ])
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains("declares 2 emitted rows"));
    }

    #[test]
    fn tsv_escapes_borrowed_strings_without_changing_semantics() {
        let output = write_streaming_table(
            Vec::new(),
            &schema(),
            OutputFormat::Tsv,
            &context(),
            None,
            |rows| {
                rows.write_row_iter([
                    CellValueRef::String("TP53\talias"),
                    CellValueRef::UInt64(9),
                    CellValueRef::Null,
                ])
            },
        )
        .unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(output.contains("# table_semantics={"));
        assert!(output.ends_with("TP53\\talias\t9\t\\N\n"));
    }

    #[test]
    fn bundle_streams_named_tables_sequentially() {
        let mut bundle = StreamingBundleWriter::new(
            Vec::new(),
            "gravlax.test.bundle.v1",
            OutputFormat::Json,
            &context(),
        )
        .unwrap();
        bundle
            .write_table("counts", &schema(), None, |rows| {
                rows.write_row_iter([
                    CellValueRef::String("TP53"),
                    CellValueRef::UInt64(9),
                    CellValueRef::Null,
                ])
            })
            .unwrap();
        bundle
            .write_table(
                "empty",
                &schema(),
                Some(&SelectionSummary::complete(0)),
                |_| Ok(()),
            )
            .unwrap();
        let output = bundle.finish().unwrap();
        let value: serde_json::Value = serde_json::from_slice(&output).unwrap();
        assert_eq!(value["result_schema"], "gravlax.test.bundle.v1");
        assert_eq!(value["data"]["tables"][0]["name"], "counts");
        assert_eq!(value["data"]["tables"][1]["name"], "empty");
        assert_eq!(value["data"]["tables"][1]["rows"], serde_json::json!([]));
    }

    #[test]
    fn duplicate_bundle_names_fail_closed() {
        let mut bundle = StreamingBundleWriter::new(
            Vec::new(),
            "gravlax.test.bundle.v1",
            OutputFormat::Json,
            &context(),
        )
        .unwrap();
        bundle
            .write_table("rows", &schema(), None, |_| Ok(()))
            .unwrap();
        assert!(bundle
            .write_table("rows", &schema(), None, |_| Ok(()))
            .is_err());
    }

    struct CountedOnePass {
        next: u64,
        polls: Arc<AtomicUsize>,
    }

    impl Iterator for CountedOnePass {
        type Item = u64;

        fn next(&mut self) -> Option<Self::Item> {
            self.polls.fetch_add(1, Ordering::Relaxed);
            let value = self.next;
            self.next += 1;
            Some(value)
        }
    }

    #[test]
    fn capped_one_pass_bundle_needs_no_prepass_for_selection_metadata() {
        let schema = TableSchema::new(
            "gravlax.test.one-pass.v1",
            vec![Field::new("ordinal", DataType::UInt64)],
        )
        .unwrap()
        .with_semantics(TableSemantics::new(RowSemantics::Sequence))
        .unwrap();
        for format in [OutputFormat::Json, OutputFormat::Tsv, OutputFormat::Text] {
            let polls = Arc::new(AtomicUsize::new(0));
            let source = CountedOnePass {
                next: 0,
                polls: Arc::clone(&polls),
            };
            let mut bundle = StreamingBundleWriter::new_with_summary(
                Vec::new(),
                "gravlax.test.one-pass-result.v1",
                format,
                &context(),
                &serde_json::json!({"cap": 3}),
            )
            .unwrap();
            bundle
                .write_table_with_deferred_selection("rows", &schema, |rows| {
                    for value in source.take(3) {
                        rows.write_row_iter([CellValueRef::UInt64(value)])?;
                    }
                    Ok(SelectionCompletion::availability_unknown())
                })
                .unwrap();
            let output = bundle.finish().unwrap();
            // Exactly the three emitted items were requested. There was no count prepass and no
            // look-ahead solely to populate metadata.
            assert_eq!(polls.load(Ordering::Relaxed), 3);
            match format {
                OutputFormat::Json => {
                    let value: serde_json::Value = serde_json::from_slice(&output).unwrap();
                    let table = &value["data"]["tables"][0];
                    assert_eq!(table["rows"].as_array().unwrap().len(), 3);
                    assert!(table["selection"]["available_rows"].is_null());
                    assert_eq!(table["selection"]["emitted_rows"], 3);
                    assert!(table["selection"]["truncated"].is_null());
                }
                OutputFormat::Tsv => {
                    let output = String::from_utf8(output).unwrap();
                    assert!(output.contains("# available_rows=unknown\n"));
                    assert!(output.contains("# emitted_rows=3\n"));
                    assert!(output.contains("# truncated=unknown\n"));
                }
                OutputFormat::Text => {
                    let output = String::from_utf8(output).unwrap();
                    assert!(output.contains("available rows: unknown\n"));
                    assert!(output.contains("emitted rows: 3\n"));
                    assert!(output.contains("truncated: unknown\n"));
                }
                OutputFormat::Arrow | OutputFormat::Mex => unreachable!(),
            }
        }
    }

    #[test]
    fn finish_time_exact_selection_uses_actual_emitted_count() {
        let schema = TableSchema::new(
            "gravlax.test.finish-selection.v1",
            vec![Field::new("ordinal", DataType::UInt64)],
        )
        .unwrap();
        let mut stream =
            StreamingTableWriter::new(Vec::new(), &schema, OutputFormat::Json, &context(), None)
                .unwrap();
        stream.write_row_iter([CellValueRef::UInt64(0)]).unwrap();
        let output = stream
            .finish_with_selection(SelectionCompletion::exact(2))
            .unwrap();
        let value: serde_json::Value = serde_json::from_slice(&output).unwrap();
        assert_eq!(value["data"]["selection"]["available_rows"], 2);
        assert_eq!(value["data"]["selection"]["emitted_rows"], 1);
        assert_eq!(value["data"]["selection"]["truncated"], true);
    }
}
