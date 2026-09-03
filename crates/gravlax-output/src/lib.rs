//! Stable, typed result envelopes and output boundaries for Gravlax commands.
//!
//! Text, TSV, and JSON are implemented without buffering the full result. Sparse matrices can be
//! written as a guarded MEX directory. Arrow IPC is deliberately kept behind [`ArrowBatchSink`],
//! which lets a CLI or Python adapter provide `arrow-rs` without imposing that dependency on every
//! Gravlax binary.

mod arrow;
mod envelope;
mod mex;
mod publish;
mod stream;
mod table;

pub use arrow::{stream_arrow, ArrowBatchSink, ColumnData, ColumnarBatch};
pub use envelope::{
    AnnotationProvenance, Producer, Provenance, ResultContext, ResultEnvelope, ENVELOPE_SCHEMA,
};
pub use mex::{write_mex_new, MatrixEntry, MatrixValue, MexFeature, MexManifest, SparseMatrix};
pub use publish::{
    canonical_destination_key, install_open_file_no_clobber, publish_file_no_clobber,
    reported_output_path, Durability, PublicationOutcome, LOGICAL_OUTPUT_MAP_ENV,
};
pub use stream::{
    write_streaming_table, write_streaming_table_with_deferred_selection, CellValueRef, RowEncoder,
    SelectionCompletion, SelectionOutcome, SelectionSummary, StreamingBundleTableWriter,
    StreamingBundleWriter, StreamingTableWriter,
};
pub use table::{
    write_table, DataType, Field, OrderKey, OutputError, OutputFormat, RowSemantics, ScalarValue,
    SortDirection, TableSchema, TableSemantics, TypedTable,
};
