use crate::table::validate_row;
use crate::{DataType, OutputError, ResultContext, ScalarValue, TableSchema};
use serde::{Deserialize, Serialize};

/// Dependency-neutral column arrays with a one-to-one Arrow mapping. `Json` is canonical compact
/// JSON carried as UTF-8; adapters must attach `gravlax.logical_type=json` to that Arrow field.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ColumnData {
    String(Vec<Option<String>>),
    Int64(Vec<Option<i64>>),
    UInt64(Vec<Option<u64>>),
    Float64(Vec<Option<f64>>),
    Boolean(Vec<Option<bool>>),
    Json(Vec<Option<String>>),
}

impl ColumnData {
    pub fn len(&self) -> usize {
        match self {
            Self::String(values) | Self::Json(values) => values.len(),
            Self::Int64(values) => values.len(),
            Self::UInt64(values) => values.len(),
            Self::Float64(values) => values.len(),
            Self::Boolean(values) => values.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn data_type(&self) -> DataType {
        match self {
            Self::String(_) => DataType::String,
            Self::Int64(_) => DataType::Int64,
            Self::UInt64(_) => DataType::UInt64,
            Self::Float64(_) => DataType::Float64,
            Self::Boolean(_) => DataType::Boolean,
            Self::Json(_) => DataType::Json,
        }
    }

    fn has_null(&self) -> bool {
        match self {
            Self::String(values) | Self::Json(values) => values.iter().any(Option::is_none),
            Self::Int64(values) => values.iter().any(Option::is_none),
            Self::UInt64(values) => values.iter().any(Option::is_none),
            Self::Float64(values) => values.iter().any(Option::is_none),
            Self::Boolean(values) => values.iter().any(Option::is_none),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ColumnarBatch {
    pub row_count: usize,
    pub columns: Vec<ColumnData>,
}

impl ColumnarBatch {
    pub fn new(
        schema: &TableSchema,
        row_count: usize,
        mut columns: Vec<ColumnData>,
    ) -> Result<Self, OutputError> {
        for (index, column) in columns.iter_mut().enumerate() {
            if let ColumnData::Json(values) = column {
                for value in values.iter_mut().flatten() {
                    *value = deterministic_compact_json(value).map_err(|error| {
                        OutputError::InvalidRow(format!(
                            "column {index} contains invalid JSON: {error}"
                        ))
                    })?;
                }
            }
        }
        let batch = Self { row_count, columns };
        batch.validate(schema)?;
        Ok(batch)
    }

    pub fn validate(&self, schema: &TableSchema) -> Result<(), OutputError> {
        if self.columns.len() != schema.fields.len() {
            return Err(OutputError::InvalidRow(format!(
                "columnar batch has {} columns for {} fields",
                self.columns.len(),
                schema.fields.len()
            )));
        }
        for (index, (column, field)) in self.columns.iter().zip(&schema.fields).enumerate() {
            if column.len() != self.row_count {
                return Err(OutputError::InvalidRow(format!(
                    "column {index} ({}) has {} values for {} rows",
                    field.name,
                    column.len(),
                    self.row_count
                )));
            }
            if column.data_type() != field.data_type {
                return Err(OutputError::InvalidRow(format!(
                    "column {index} ({}) is {} but schema requires {}",
                    field.name,
                    column.data_type().name(),
                    field.data_type.name()
                )));
            }
            if column.has_null() && !field.nullable {
                return Err(OutputError::InvalidRow(format!(
                    "non-nullable column {index} ({}) contains null",
                    field.name
                )));
            }
            if let ColumnData::Float64(values) = column {
                if values.iter().flatten().any(|value| !value.is_finite()) {
                    return Err(OutputError::InvalidRow(format!(
                        "column {index} ({}) contains a non-finite float",
                        field.name
                    )));
                }
            }
            if let ColumnData::Json(values) = column {
                for value in values.iter().flatten() {
                    let normalized = deterministic_compact_json(value).map_err(|error| {
                        OutputError::InvalidRow(format!(
                            "column {index} ({}) contains invalid JSON: {error}",
                            field.name
                        ))
                    })?;
                    if normalized != *value {
                        return Err(OutputError::InvalidRow(format!(
                            "column {index} ({}) contains non-deterministic JSON; object keys must be sorted and insignificant whitespace removed",
                            field.name
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    pub fn from_rows(schema: &TableSchema, rows: &[Vec<ScalarValue>]) -> Result<Self, OutputError> {
        for (index, row) in rows.iter().enumerate() {
            validate_row(schema, index, row)?;
        }
        let mut columns: Vec<ColumnData> = schema
            .fields
            .iter()
            .map(|field| match field.data_type {
                DataType::String => ColumnData::String(Vec::with_capacity(rows.len())),
                DataType::Int64 => ColumnData::Int64(Vec::with_capacity(rows.len())),
                DataType::UInt64 => ColumnData::UInt64(Vec::with_capacity(rows.len())),
                DataType::Float64 => ColumnData::Float64(Vec::with_capacity(rows.len())),
                DataType::Boolean => ColumnData::Boolean(Vec::with_capacity(rows.len())),
                DataType::Json => ColumnData::Json(Vec::with_capacity(rows.len())),
            })
            .collect();
        for row in rows {
            for (value, column) in row.iter().zip(&mut columns) {
                match (value, column) {
                    (ScalarValue::Null, ColumnData::String(values))
                    | (ScalarValue::Null, ColumnData::Json(values)) => values.push(None),
                    (ScalarValue::Null, ColumnData::Int64(values)) => values.push(None),
                    (ScalarValue::Null, ColumnData::UInt64(values)) => values.push(None),
                    (ScalarValue::Null, ColumnData::Float64(values)) => values.push(None),
                    (ScalarValue::Null, ColumnData::Boolean(values)) => values.push(None),
                    (ScalarValue::String(value), ColumnData::String(values)) => {
                        values.push(Some(value.clone()))
                    }
                    (ScalarValue::Int64(value), ColumnData::Int64(values)) => {
                        values.push(Some(*value))
                    }
                    (ScalarValue::UInt64(value), ColumnData::UInt64(values)) => {
                        values.push(Some(*value))
                    }
                    (ScalarValue::Float64(value), ColumnData::Float64(values)) => {
                        values.push(Some(*value))
                    }
                    (ScalarValue::Boolean(value), ColumnData::Boolean(values)) => {
                        values.push(Some(*value))
                    }
                    (ScalarValue::Json(value), ColumnData::Json(values)) => {
                        values.push(Some(deterministic_compact_json_value(value)?))
                    }
                    _ => unreachable!("rows were validated against the same schema"),
                }
            }
        }
        Self::new(schema, rows.len(), columns)
    }
}

fn deterministic_compact_json(input: &str) -> Result<String, serde_json::Error> {
    let mut value: serde_json::Value = serde_json::from_str(input)?;
    sort_json_objects(&mut value);
    serde_json::to_string(&value)
}

fn deterministic_compact_json_value(value: &serde_json::Value) -> Result<String, OutputError> {
    let mut value = value.clone();
    sort_json_objects(&mut value);
    Ok(serde_json::to_string(&value)?)
}

fn sort_json_objects(value: &mut serde_json::Value) {
    match value {
        serde_json::Value::Array(values) => {
            for value in values {
                sort_json_objects(value);
            }
        }
        serde_json::Value::Object(object) => {
            let mut entries: Vec<_> = std::mem::take(object).into_iter().collect();
            for (_, value) in &mut entries {
                sort_json_objects(value);
            }
            entries.sort_by(|left, right| left.0.cmp(&right.0));
            object.extend(entries);
        }
        _ => {}
    }
}

/// Adapter point for Arrow IPC, PyArrow, Polars, or another columnar consumer. Implementations
/// receive exactly one schema and context followed by bounded batches and a finish notification.
pub trait ArrowBatchSink {
    fn begin(&mut self, schema: &TableSchema, context: &ResultContext) -> Result<(), OutputError>;
    fn write_batch(&mut self, batch: &ColumnarBatch) -> Result<(), OutputError>;
    fn finish(&mut self) -> Result<(), OutputError>;
}

pub fn stream_arrow<S, I>(
    sink: &mut S,
    schema: &TableSchema,
    batches: I,
    context: &ResultContext,
) -> Result<(), OutputError>
where
    S: ArrowBatchSink,
    I: IntoIterator<Item = ColumnarBatch>,
{
    schema.validate()?;
    context.validate()?;
    sink.begin(schema, context)?;
    for batch in batches {
        // Re-run structural validation at the trust boundary, including for batches built by
        // deserialization rather than ColumnarBatch::new.
        batch.validate(schema)?;
        sink.write_batch(&batch)?;
    }
    sink.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Field, ResultContext};

    #[derive(Default)]
    struct CountingSink {
        began: bool,
        rows: usize,
        finished: bool,
    }

    impl ArrowBatchSink for CountingSink {
        fn begin(
            &mut self,
            _schema: &TableSchema,
            _context: &ResultContext,
        ) -> Result<(), OutputError> {
            self.began = true;
            Ok(())
        }

        fn write_batch(&mut self, batch: &ColumnarBatch) -> Result<(), OutputError> {
            self.rows += batch.row_count;
            Ok(())
        }

        fn finish(&mut self) -> Result<(), OutputError> {
            self.finished = true;
            Ok(())
        }
    }

    #[test]
    fn streams_bounded_columnar_batches_through_adapter() {
        let schema = TableSchema::new(
            "gravlax.test.arrow.v1",
            vec![Field::new("count", DataType::UInt64)],
        )
        .unwrap();
        let batch = ColumnarBatch::from_rows(
            &schema,
            &[vec![ScalarValue::UInt64(1)], vec![ScalarValue::UInt64(2)]],
        )
        .unwrap();
        let mut sink = CountingSink::default();
        stream_arrow(&mut sink, &schema, [batch], &ResultContext::default()).unwrap();
        assert!(sink.began && sink.finished);
        assert_eq!(sink.rows, 2);
    }

    #[test]
    fn rejects_wrong_column_types_and_nonfinite_values() {
        let schema = TableSchema::new(
            "gravlax.test.arrow.v1",
            vec![Field::new("score", DataType::Float64)],
        )
        .unwrap();
        assert!(ColumnarBatch::new(&schema, 1, vec![ColumnData::UInt64(vec![Some(1)])]).is_err());
        assert!(
            ColumnarBatch::new(&schema, 1, vec![ColumnData::Float64(vec![Some(f64::NAN)])])
                .is_err()
        );
    }

    #[test]
    fn json_columns_are_sorted_and_compacted_on_construction() {
        let schema = TableSchema::new(
            "gravlax.test.arrow-json.v1",
            vec![Field::new("payload", DataType::Json)],
        )
        .unwrap();
        let batch = ColumnarBatch::new(
            &schema,
            1,
            vec![ColumnData::Json(vec![Some(
                r#" { "z": 2, "a": [{"z": 1, "a": 0}] } "#.into(),
            )])],
        )
        .unwrap();
        assert_eq!(
            batch.columns,
            vec![ColumnData::Json(vec![Some(
                r#"{"a":[{"a":0,"z":1}],"z":2}"#.into()
            )])]
        );
    }

    #[test]
    fn validation_rejects_invalid_or_non_compact_json_from_deserialization() {
        let schema = TableSchema::new(
            "gravlax.test.arrow-json.v1",
            vec![Field::new("payload", DataType::Json)],
        )
        .unwrap();
        let invalid = ColumnarBatch {
            row_count: 1,
            columns: vec![ColumnData::Json(vec![Some("{".into())])],
        };
        assert!(invalid.validate(&schema).is_err());
        let non_compact = ColumnarBatch {
            row_count: 1,
            columns: vec![ColumnData::Json(vec![Some(r#"{"z":2,"a":1}"#.into())])],
        };
        assert!(non_compact.validate(&schema).is_err());
    }
}
