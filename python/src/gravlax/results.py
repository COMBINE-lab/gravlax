"""Parser and optional converters for ``gravlax.result-envelope.v1``."""

from __future__ import annotations

import importlib
import json
import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Optional, Union

from ._json import (
    integer,
    mapping,
    nonempty_string,
    required,
    sequence,
    string,
    strings,
    validate_json_value,
)
from ._json import loads_document
from .exceptions import OptionalDependencyError, ProtocolError

ENVELOPE_SCHEMA = "gravlax.result-envelope.v1"
_DATA_TYPES = {"string", "int64", "uint64", "float64", "boolean", "json"}
_MAX_UINT64 = 2**64 - 1
_MAX_UINT32 = 2**32 - 1
_MISSING = object()


def _nonblank_string(value: Any, where: str) -> str:
    result = string(value, where)
    if not result.strip():
        raise ProtocolError(f"{where} must not be empty or whitespace")
    return result


def _schema_id(value: Any, where: str) -> str:
    result = string(value, where)
    valid = bool(result) and all(
        character.isascii()
        and (character.isalnum() or character in {".", "_", "-"})
        for character in result
    )
    if not valid:
        raise ProtocolError(
            f"{where} must contain only non-empty ASCII letters, digits, '.', '_' or '-'"
        )
    return result


def _bounded_integer(value: Any, where: str, *, maximum: int) -> int:
    result = integer(value, where, minimum=0)
    if result > maximum:
        raise ProtocolError(f"{where} must not exceed {maximum}")
    return result


def _uint64(value: Any, where: str) -> int:
    return _bounded_integer(value, where, maximum=_MAX_UINT64)


def _uint32(value: Any, where: str) -> int:
    return _bounded_integer(value, where, maximum=_MAX_UINT32)


def _optional(module: str, extra: str) -> Any:
    try:
        return importlib.import_module(module)
    except ImportError as error:
        raise OptionalDependencyError(
            f"{module} is required for this conversion; "
            f"install 'gravlax-client[{extra}]'"
        ) from error


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    )


@dataclass(frozen=True)
class Producer:
    name: str
    version: str

    @classmethod
    def from_mapping(cls, value: Any, where: str = "producer") -> "Producer":
        document = mapping(value, where)
        return cls(
            name=_nonblank_string(required(document, "name", where), f"{where}.name"),
            version=_nonblank_string(
                required(document, "version", where), f"{where}.version"
            ),
        )


@dataclass(frozen=True)
class AnnotationProvenance:
    role: str
    assembly: str
    annotation: str
    digest: str

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "AnnotationProvenance":
        document = mapping(value, where)
        return cls(
            role=_nonblank_string(required(document, "role", where), f"{where}.role"),
            assembly=_nonblank_string(
                required(document, "assembly", where), f"{where}.assembly"
            ),
            annotation=_nonblank_string(
                required(document, "annotation", where), f"{where}.annotation"
            ),
            digest=_nonblank_string(
                required(document, "digest", where), f"{where}.digest"
            ),
        )

    def as_dict(self) -> dict[str, str]:
        return {
            "role": self.role,
            "assembly": self.assembly,
            "annotation": self.annotation,
            "digest": self.digest,
        }


@dataclass(frozen=True)
class Provenance:
    archives: tuple[str, ...]
    assembly: Optional[str]
    annotation: Optional[str]
    annotation_digest: Optional[str]
    annotations: tuple[AnnotationProvenance, ...]
    parameters: Mapping[str, Any]

    @classmethod
    def from_mapping(cls, value: Any, where: str = "provenance") -> "Provenance":
        document = mapping(value, where)
        assembly = document.get("assembly")
        annotation = document.get("annotation")
        digest = document.get("annotation_digest")
        parameters = mapping(document.get("parameters", {}), f"{where}.parameters")
        archives = strings(document.get("archives", []), f"{where}.archives")
        annotation_values = sequence(
            document.get("annotations", []), f"{where}.annotations"
        )
        annotations = tuple(
            AnnotationProvenance.from_mapping(
                item, f"{where}.annotations[{index}]"
            )
            for index, item in enumerate(annotation_values)
        )
        roles = [item.role for item in annotations]
        if len(roles) != len(set(roles)):
            raise ProtocolError(f"{where}.annotations contains duplicate roles")
        seen_archives: set[str] = set()
        for index, archive in enumerate(archives):
            if not archive.strip():
                raise ProtocolError(f"{where}.archives[{index}] must not be empty")
            if archive in seen_archives:
                raise ProtocolError(
                    f"{where}.archives contains duplicate identity {archive!r}"
                )
            seen_archives.add(archive)
        for name in parameters:
            if not name.strip():
                raise ProtocolError(f"{where}.parameters names must not be empty")
        return cls(
            archives=archives,
            assembly=(
                None
                if assembly is None
                else _nonblank_string(assembly, f"{where}.assembly")
            ),
            annotation=(
                None
                if annotation is None
                else _nonblank_string(annotation, f"{where}.annotation")
            ),
            annotation_digest=(
                None
                if digest is None
                else _nonblank_string(digest, f"{where}.annotation_digest")
            ),
            annotations=annotations,
            parameters=dict(parameters),
        )

    def as_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "archives": list(self.archives),
            "parameters": dict(self.parameters),
        }
        if self.assembly is not None:
            result["assembly"] = self.assembly
        if self.annotation is not None:
            result["annotation"] = self.annotation
        if self.annotation_digest is not None:
            result["annotation_digest"] = self.annotation_digest
        if self.annotations:
            result["annotations"] = [item.as_dict() for item in self.annotations]
        return result


@dataclass(frozen=True)
class Field:
    name: str
    data_type: str
    nullable: bool
    description: Optional[str] = None

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "Field":
        document = mapping(value, where)
        data_type = nonempty_string(
            required(document, "data_type", where), f"{where}.data_type"
        )
        if data_type not in _DATA_TYPES:
            raise ProtocolError(f"{where}.data_type has unknown value {data_type!r}")
        nullable = document.get("nullable", False)
        if not isinstance(nullable, bool):
            raise ProtocolError(f"{where}.nullable must be a boolean")
        description = document.get("description")
        name = nonempty_string(required(document, "name", where), f"{where}.name")
        if any(character in name for character in "\t\r\n"):
            raise ProtocolError(f"{where}.name contains a tab or line break")
        return cls(
            name=name,
            data_type=data_type,
            nullable=nullable,
            description=(
                None
                if description is None
                else string(description, f"{where}.description")
            ),
        )


class RowSemantics(str, Enum):
    """Logical meaning of a table's rows, independent of physical row order."""

    SET = "set"
    MULTISET = "multiset"
    SEQUENCE = "sequence"


class SortDirection(str, Enum):
    ASCENDING = "ascending"
    DESCENDING = "descending"


@dataclass(frozen=True)
class OrderKey:
    field: str
    direction: SortDirection

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "OrderKey":
        document = mapping(value, where)
        field = nonempty_string(required(document, "field", where), f"{where}.field")
        try:
            direction = SortDirection(
                string(required(document, "direction", where), f"{where}.direction")
            )
        except ValueError as error:
            raise ProtocolError(
                f"{where}.direction must be 'ascending' or 'descending'"
            ) from error
        return cls(field=field, direction=direction)


@dataclass(frozen=True)
class TableSemantics:
    row_semantics: RowSemantics
    key: Optional[tuple[str, ...]] = None
    ordered_by: Optional[tuple[OrderKey, ...]] = None

    @classmethod
    def from_mapping(
        cls,
        value: Any,
        fields: tuple[Field, ...],
        where: str,
    ) -> "TableSemantics":
        document = mapping(value, where)
        try:
            row_semantics = RowSemantics(
                string(
                    required(document, "row_semantics", where),
                    f"{where}.row_semantics",
                )
            )
        except ValueError as error:
            raise ProtocolError(
                f"{where}.row_semantics must be 'set', 'multiset', or 'sequence'"
            ) from error

        key_value = document.get("key")
        key = None if key_value is None else strings(key_value, f"{where}.key")
        if key == ():
            raise ProtocolError(f"{where}.key must be absent rather than empty")

        ordered_value = document.get("ordered_by")
        ordered_by = (
            None
            if ordered_value is None
            else tuple(
                OrderKey.from_mapping(item, f"{where}.ordered_by[{index}]")
                for index, item in enumerate(
                    sequence(ordered_value, f"{where}.ordered_by")
                )
            )
        )
        if ordered_by == ():
            raise ProtocolError(
                f"{where}.ordered_by must be absent rather than empty"
            )

        field_names = {field.name for field in fields}
        for label, references in (
            ("key", key),
            (
                "ordered_by",
                None
                if ordered_by is None
                else tuple(order.field for order in ordered_by),
            ),
        ):
            if references is None:
                continue
            if len(references) != len(set(references)):
                raise ProtocolError(f"{where}.{label} repeats a field")
            unknown = [name for name in references if name not in field_names]
            if unknown:
                raise ProtocolError(
                    f"{where}.{label} references unknown field {unknown[0]!r}"
                )
        return cls(
            row_semantics=row_semantics,
            key=key,
            ordered_by=ordered_by,
        )


def _validate_cell(value: Any, field: Field, where: str) -> None:
    if value is None:
        if not field.nullable:
            raise ProtocolError(f"{where} is null but {field.name!r} is not nullable")
        return
    valid = False
    if field.data_type == "string":
        valid = isinstance(value, str)
    elif field.data_type == "int64":
        valid = (
            isinstance(value, int)
            and not isinstance(value, bool)
            and -(2**63) <= value < 2**63
        )
    elif field.data_type == "uint64":
        valid = (
            isinstance(value, int)
            and not isinstance(value, bool)
            and 0 <= value < 2**64
        )
    elif field.data_type == "float64":
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            try:
                valid = math.isfinite(float(value))
            except OverflowError:
                valid = False
    elif field.data_type == "boolean":
        valid = isinstance(value, bool)
    elif field.data_type == "json":
        try:
            _canonical_json(value)
            valid = True
        except (TypeError, ValueError):
            valid = False
    if not valid:
        raise ProtocolError(
            f"{where} does not match declared type {field.data_type!r}"
        )


@dataclass(frozen=True)
class Table:
    schema_id: str
    fields: tuple[Field, ...]
    rows: tuple[tuple[Any, ...], ...]
    semantics: Optional[TableSemantics] = None

    @classmethod
    def from_mapping(
        cls, value: Any, result_schema: str, where: str = "result.data"
    ) -> "Table":
        document = mapping(value, where)
        schema_document = mapping(
            required(document, "schema", where), f"{where}.schema"
        )
        schema_id = _schema_id(
            required(schema_document, "id", f"{where}.schema"),
            f"{where}.schema.id",
        )
        if schema_id != result_schema:
            raise ProtocolError(
                f"table schema {schema_id!r} does not match envelope result_schema "
                f"{result_schema!r}"
            )
        field_values = sequence(
            required(schema_document, "fields", f"{where}.schema"),
            f"{where}.schema.fields",
        )
        fields = tuple(
            Field.from_mapping(item, f"{where}.schema.fields[{index}]")
            for index, item in enumerate(field_values)
        )
        if not fields:
            raise ProtocolError("result table must declare at least one field")
        names = [field.name for field in fields]
        if len(names) != len(set(names)):
            raise ProtocolError("result table contains duplicate field names")
        semantics_value = schema_document.get("semantics")
        semantics = (
            None
            if semantics_value is None
            else TableSemantics.from_mapping(
                semantics_value,
                fields,
                f"{where}.schema.semantics",
            )
        )
        row_values = sequence(required(document, "rows", where), f"{where}.rows")
        rows: list[tuple[Any, ...]] = []
        for row_index, value in enumerate(row_values):
            cells = tuple(sequence(value, f"{where}.rows[{row_index}]"))
            if len(cells) != len(fields):
                raise ProtocolError(
                    f"{where}.rows[{row_index}] has {len(cells)} values for "
                    f"{len(fields)} fields"
                )
            for column, (cell, field) in enumerate(zip(cells, fields)):
                _validate_cell(cell, field, f"{where}.rows[{row_index}][{column}]")
            rows.append(cells)
        result = cls(
            schema_id=schema_id,
            fields=fields,
            rows=tuple(rows),
            semantics=semantics,
        )
        result._validate_keys(where)
        return result

    def _validate_keys(self, where: str) -> None:
        if self.semantics is None:
            return
        key = self.semantics.key
        if key is None and self.semantics.row_semantics is not RowSemantics.SET:
            return
        if key is None:
            columns = tuple(range(len(self.fields)))
        else:
            positions = {field.name: index for index, field in enumerate(self.fields)}
            columns = tuple(positions[name] for name in key)
        seen: set[str] = set()
        for row_index, row in enumerate(self.rows):
            encoded = _canonical_json([row[column] for column in columns])
            if encoded in seen:
                label = "declared key" if key is not None else "set row"
                raise ProtocolError(
                    f"{where}.rows[{row_index}] duplicates a {label}"
                )
            seen.add(encoded)

    @property
    def columns(self) -> tuple[str, ...]:
        return tuple(field.name for field in self.fields)

    def records(self) -> list[dict[str, Any]]:
        """Return rows as ordinary dictionaries, without optional dependencies."""

        return [dict(zip(self.columns, row)) for row in self.rows]

    def to_pandas(self) -> Any:
        """Return a pandas DataFrame with nullable logical dtypes where needed."""

        pandas = _optional("pandas", "pandas")
        frame = pandas.DataFrame.from_records(self.records(), columns=self.columns)
        for field in self.fields:
            if field.data_type == "json":
                continue
            if field.data_type == "string":
                dtype = "string"
            elif field.data_type == "int64":
                dtype = "Int64" if field.nullable else "int64"
            elif field.data_type == "uint64":
                dtype = "UInt64" if field.nullable else "uint64"
            elif field.data_type == "float64":
                dtype = "Float64" if field.nullable else "float64"
            else:
                dtype = "boolean" if field.nullable else "bool"
            frame[field.name] = frame[field.name].astype(dtype)
        return frame

    def to_arrow(self) -> Any:
        """Return a PyArrow Table preserving Gravlax logical types."""

        pyarrow = _optional("pyarrow", "arrow")
        arrow_types = {
            "string": pyarrow.string(),
            "int64": pyarrow.int64(),
            "uint64": pyarrow.uint64(),
            "float64": pyarrow.float64(),
            "boolean": pyarrow.bool_(),
            "json": pyarrow.string(),
        }
        arrays = []
        arrow_fields = []
        for column, field in enumerate(self.fields):
            values = [row[column] for row in self.rows]
            metadata = None
            if field.data_type == "json":
                values = [None if value is None else _canonical_json(value) for value in values]
                metadata = {b"gravlax.logical_type": b"json"}
            arrays.append(pyarrow.array(values, type=arrow_types[field.data_type]))
            arrow_fields.append(
                pyarrow.field(
                    field.name,
                    arrow_types[field.data_type],
                    nullable=field.nullable,
                    metadata=metadata,
                )
            )
        return pyarrow.Table.from_arrays(arrays, schema=pyarrow.schema(arrow_fields))


@dataclass(frozen=True)
class ResultEnvelope:
    envelope_schema: str
    result_schema: str
    producer: Producer
    provenance: Provenance
    warnings: tuple[str, ...]
    data: Any

    @classmethod
    def from_mapping(cls, value: Any) -> "ResultEnvelope":
        where = "result envelope"
        validate_json_value(value, where)
        document = mapping(value, where)
        schema = nonempty_string(
            required(document, "$schema", where), f"{where}.$schema"
        )
        if schema != ENVELOPE_SCHEMA:
            raise ProtocolError(
                f"unsupported result envelope {schema!r}; expected {ENVELOPE_SCHEMA!r}"
            )
        return cls(
            envelope_schema=schema,
            result_schema=_schema_id(
                required(document, "result_schema", where), f"{where}.result_schema"
            ),
            producer=Producer.from_mapping(
                required(document, "producer", where), f"{where}.producer"
            ),
            provenance=Provenance.from_mapping(
                required(document, "provenance", where), f"{where}.provenance"
            ),
            warnings=strings(document.get("warnings", []), f"{where}.warnings"),
            data=required(document, "data", where),
        )

    @classmethod
    def from_json(cls, source: str | bytes) -> "ResultEnvelope":
        return cls.from_mapping(loads_document(source, "result response"))

    @classmethod
    def from_file(cls, path: str | Path) -> "ResultEnvelope":
        try:
            source = Path(path).read_bytes()
        except OSError as error:
            raise ProtocolError(f"cannot read result envelope {path}: {error}") from error
        return cls.from_json(source)

    @property
    def table(self) -> Table:
        """Return the envelope's typed table or raise for a non-table result."""

        document = mapping(self.data, "result.data")
        if "schema" not in document or "rows" not in document:
            raise ProtocolError(
                f"result {self.result_schema!r} is not a typed tabular result"
            )
        return Table.from_mapping(document, self.result_schema)

    def to_pandas(self) -> Any:
        return self.table.to_pandas()

    def to_arrow(self) -> Any:
        table = self.table.to_arrow()
        metadata = {
            b"gravlax.envelope_schema": self.envelope_schema.encode("utf-8"),
            b"gravlax.result_schema": self.result_schema.encode("utf-8"),
            b"gravlax.producer": _canonical_json(
                {"name": self.producer.name, "version": self.producer.version}
            ).encode("utf-8"),
            b"gravlax.provenance": _canonical_json(self.provenance.as_dict()).encode(
                "utf-8"
            ),
            b"gravlax.warnings": _canonical_json(list(self.warnings)).encode("utf-8"),
        }
        return table.replace_schema_metadata(metadata)

    def to_anndata(self, *, obs_names: str | None = None) -> Any:
        """Represent a typed row result as AnnData observations with no inferred matrix."""

        anndata = _optional("anndata", "anndata")
        frame = self.to_pandas()
        for field in self.table.fields:
            if field.data_type == "json":
                frame[field.name] = frame[field.name].map(
                    lambda value: None if value is None else _canonical_json(value)
                )
        if obs_names is not None:
            if obs_names not in frame.columns:
                raise KeyError(f"observation-name column {obs_names!r} is not present")
            frame.index = frame[obs_names].astype(str)
            frame.index.name = obs_names
        result = anndata.AnnData(obs=frame)
        result.uns["gravlax"] = self.metadata()
        return result

    def metadata(self) -> dict[str, Any]:
        return {
            "envelope_schema": self.envelope_schema,
            "result_schema": self.result_schema,
            "producer": {
                "name": self.producer.name,
                "version": self.producer.version,
            },
            "provenance": self.provenance.as_dict(),
            "warnings": list(self.warnings),
        }

    @property
    def bundle(self) -> "UniformResultBundle":
        """Return named streamed tables from a uniform multi-table result."""

        return UniformResultBundle.from_envelope(self)


@dataclass(frozen=True)
class ExactSelection:
    """Selection metadata whose complete available-row count is known."""

    available_rows: int
    emitted_rows: int
    truncated: bool


@dataclass(frozen=True)
class DeferredSelection:
    """One-pass selection for which total availability was deliberately not scanned."""

    available_rows: None
    emitted_rows: int
    truncated: Optional[bool]


TableSelection = Union[ExactSelection, DeferredSelection]


def _parse_selection(value: Any, where: str) -> TableSelection:
    document = mapping(value, where)
    available_value = required(document, "available_rows", where)
    emitted_rows = _uint64(
        required(document, "emitted_rows", where), f"{where}.emitted_rows"
    )
    truncated_value = required(document, "truncated", where)

    if available_value is None:
        if truncated_value is not None and not isinstance(truncated_value, bool):
            raise ProtocolError(f"{where}.truncated must be a boolean or null")
        if truncated_value is False:
            raise ProtocolError(
                f"{where} cannot claim complete selection while available_rows is null"
            )
        return DeferredSelection(
            available_rows=None,
            emitted_rows=emitted_rows,
            truncated=truncated_value,
        )

    available_rows = _uint64(available_value, f"{where}.available_rows")
    if not isinstance(truncated_value, bool):
        raise ProtocolError(
            f"{where}.truncated must be a boolean when available_rows is known"
        )
    if emitted_rows > available_rows:
        raise ProtocolError(
            f"{where}.emitted_rows exceeds {where}.available_rows"
        )
    expected_truncation = emitted_rows < available_rows
    if truncated_value != expected_truncation:
        raise ProtocolError(
            f"{where}.truncated must equal emitted_rows < available_rows"
        )
    return ExactSelection(
        available_rows=available_rows,
        emitted_rows=emitted_rows,
        truncated=truncated_value,
    )


@dataclass(frozen=True)
class BundleSummary:
    """Validated command-specific JSON summary for a future or unknown bundle schema."""

    values: Mapping[str, Any]

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "BundleSummary":
        document = mapping(value, where)
        validate_json_value(document, where)
        return cls(values=dict(document))

    def as_dict(self) -> dict[str, Any]:
        return dict(self.values)

    def __getitem__(self, name: str) -> Any:
        return self.values[name]

    def get(self, name: str, default: Any = None) -> Any:
        return self.values.get(name, default)


@dataclass(frozen=True)
class RegionQuerySummary:
    coordinates: str
    anchor_semantics: bool
    chrom: str
    start: int
    end: int
    molecules: int
    umis: int
    cells: int
    chunks_decoded: int

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "RegionQuerySummary":
        document = mapping(value, where)
        coordinates = string(
            required(document, "coordinates", where), f"{where}.coordinates"
        )
        if coordinates != "0-based half-open":
            raise ProtocolError(
                f"{where}.coordinates must be '0-based half-open'"
            )
        anchor_semantics = required(document, "anchor_semantics", where)
        if anchor_semantics is not True:
            raise ProtocolError(f"{where}.anchor_semantics must be true")
        start = _uint32(required(document, "start", where), f"{where}.start")
        end = _uint32(required(document, "end", where), f"{where}.end")
        if start >= end:
            raise ProtocolError(f"{where} must have start < end")
        return cls(
            coordinates=coordinates,
            anchor_semantics=True,
            chrom=_nonblank_string(
                required(document, "chrom", where), f"{where}.chrom"
            ),
            start=start,
            end=end,
            molecules=_uint64(
                required(document, "molecules", where), f"{where}.molecules"
            ),
            umis=_uint64(required(document, "umis", where), f"{where}.umis"),
            cells=_uint64(required(document, "cells", where), f"{where}.cells"),
            chunks_decoded=_uint64(
                required(document, "chunks_decoded", where),
                f"{where}.chunks_decoded",
            ),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "coordinates": self.coordinates,
            "anchor_semantics": self.anchor_semantics,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "molecules": self.molecules,
            "umis": self.umis,
            "cells": self.cells,
            "chunks_decoded": self.chunks_decoded,
        }


@dataclass(frozen=True)
class JunctionQuerySummary:
    coordinates: str
    chrom: str
    donor: int
    acceptor: int
    archive_supporting_children: int
    archive_posting_chunks: int
    umis: int
    cells: int

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "JunctionQuerySummary":
        document = mapping(value, where)
        coordinates = string(
            required(document, "coordinates", where), f"{where}.coordinates"
        )
        if coordinates != "0-based half-open junction boundaries":
            raise ProtocolError(
                f"{where}.coordinates must be '0-based half-open junction boundaries'"
            )
        donor = _uint32(required(document, "donor", where), f"{where}.donor")
        acceptor = _uint32(
            required(document, "acceptor", where), f"{where}.acceptor"
        )
        if donor >= acceptor:
            raise ProtocolError(f"{where} must have donor < acceptor")
        return cls(
            coordinates=coordinates,
            chrom=_nonblank_string(
                required(document, "chrom", where), f"{where}.chrom"
            ),
            donor=donor,
            acceptor=acceptor,
            archive_supporting_children=_uint64(
                required(document, "archive_supporting_children", where),
                f"{where}.archive_supporting_children",
            ),
            archive_posting_chunks=_uint64(
                required(document, "archive_posting_chunks", where),
                f"{where}.archive_posting_chunks",
            ),
            umis=_uint64(required(document, "umis", where), f"{where}.umis"),
            cells=_uint64(required(document, "cells", where), f"{where}.cells"),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "coordinates": self.coordinates,
            "chrom": self.chrom,
            "donor": self.donor,
            "acceptor": self.acceptor,
            "archive_supporting_children": self.archive_supporting_children,
            "archive_posting_chunks": self.archive_posting_chunks,
            "umis": self.umis,
            "cells": self.cells,
        }


UniformBundleSummary = Union[
    BundleSummary,
    RegionQuerySummary,
    JunctionQuerySummary,
]


@dataclass(frozen=True)
class NamedTable:
    name: str
    table: Table
    selection: Optional[TableSelection]

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "NamedTable":
        document = mapping(value, where)
        name = _schema_id(required(document, "name", where), f"{where}.name")
        schema_document = mapping(
            required(document, "schema", where), f"{where}.schema"
        )
        schema_id = _schema_id(
            required(schema_document, "id", f"{where}.schema"),
            f"{where}.schema.id",
        )
        table = Table.from_mapping(document, schema_id, where)
        if table.semantics is None:
            raise ProtocolError(
                f"{where}.schema.semantics is required for a uniform bundle table"
            )
        selection_value = document.get("selection", _MISSING)
        selection = (
            None
            if selection_value is _MISSING
            else _parse_selection(selection_value, f"{where}.selection")
        )
        if selection is not None and selection.emitted_rows != len(table.rows):
            raise ProtocolError(
                f"{where}.selection.emitted_rows disagrees with the number of rows"
            )
        return cls(name=name, table=table, selection=selection)

    @property
    def schema_id(self) -> str:
        return self.table.schema_id

    @property
    def fields(self) -> tuple[Field, ...]:
        return self.table.fields

    @property
    def rows(self) -> tuple[tuple[Any, ...], ...]:
        return self.table.rows

    @property
    def semantics(self) -> TableSemantics:
        assert self.table.semantics is not None
        return self.table.semantics

    def records(self) -> list[dict[str, Any]]:
        return self.table.records()

    def to_pandas(self) -> Any:
        return self.table.to_pandas()

    def to_arrow(self) -> Any:
        return self.table.to_arrow()


REGION_UNIFORM_RESULT_SCHEMA = "gravlax.query.region.result.v1"
REGION_UNIFORM_COUNTS_SCHEMA = "gravlax.query.region.counts.v1"
JUNCTION_UNIFORM_RESULT_SCHEMA = "gravlax.query.junction.result.v1"
JUNCTION_UNIFORM_COUNTS_SCHEMA = "gravlax.query.junction.counts.v1"


@dataclass(frozen=True)
class UniformResultBundle:
    """Strict named-table view of a streamed uniform JSON result bundle."""

    envelope: ResultEnvelope
    summary: Optional[UniformBundleSummary]
    tables: tuple[NamedTable, ...]

    @classmethod
    def from_envelope(cls, envelope: ResultEnvelope) -> "UniformResultBundle":
        data = mapping(envelope.data, "uniform result data")
        table_values = sequence(
            required(data, "tables", "uniform result data"),
            "uniform result data.tables",
        )
        tables = tuple(
            NamedTable.from_mapping(item, f"uniform result data.tables[{index}]")
            for index, item in enumerate(table_values)
        )
        names = [table.name for table in tables]
        if len(names) != len(set(names)):
            raise ProtocolError("uniform result data.tables contains duplicate names")

        summary_value = data.get("summary", _MISSING)
        if summary_value is _MISSING:
            summary: Optional[UniformBundleSummary] = None
        elif envelope.result_schema == REGION_UNIFORM_RESULT_SCHEMA:
            summary = RegionQuerySummary.from_mapping(
                summary_value, "uniform result data.summary"
            )
        elif envelope.result_schema == JUNCTION_UNIFORM_RESULT_SCHEMA:
            summary = JunctionQuerySummary.from_mapping(
                summary_value, "uniform result data.summary"
            )
        else:
            summary = BundleSummary.from_mapping(
                summary_value, "uniform result data.summary"
            )

        result = cls(envelope=envelope, summary=summary, tables=tables)
        result._validate_known_query_bundle()
        return result

    @classmethod
    def from_mapping(cls, value: Any) -> "UniformResultBundle":
        return cls.from_envelope(ResultEnvelope.from_mapping(value))

    @classmethod
    def from_json(cls, source: str | bytes) -> "UniformResultBundle":
        return cls.from_envelope(ResultEnvelope.from_json(source))

    @classmethod
    def from_file(cls, path: str | Path) -> "UniformResultBundle":
        try:
            source = Path(path).read_bytes()
        except OSError as error:
            raise ProtocolError(f"cannot read uniform result bundle {path}: {error}") from error
        return cls.from_json(source)

    @property
    def result_schema(self) -> str:
        return self.envelope.result_schema

    @property
    def provenance(self) -> Provenance:
        return self.envelope.provenance

    @property
    def warnings(self) -> tuple[str, ...]:
        return self.envelope.warnings

    @property
    def table_names(self) -> tuple[str, ...]:
        return tuple(table.name for table in self.tables)

    def table(self, name: str) -> NamedTable:
        for table in self.tables:
            if table.name == name:
                return table
        raise KeyError(f"uniform result has no table {name!r}")

    def _validate_known_query_bundle(self) -> None:
        expected_table_schema: Optional[str]
        expected_summary: type[Any]
        if self.result_schema == REGION_UNIFORM_RESULT_SCHEMA:
            expected_table_schema = REGION_UNIFORM_COUNTS_SCHEMA
            expected_summary = RegionQuerySummary
        elif self.result_schema == JUNCTION_UNIFORM_RESULT_SCHEMA:
            expected_table_schema = JUNCTION_UNIFORM_COUNTS_SCHEMA
            expected_summary = JunctionQuerySummary
        else:
            return

        if not isinstance(self.summary, expected_summary):
            raise ProtocolError(
                f"uniform result {self.result_schema!r} requires a typed summary"
            )
        if len(self.tables) != 1 or self.tables[0].name != "counts":
            raise ProtocolError(
                f"uniform result {self.result_schema!r} requires exactly the 'counts' table"
            )
        named = self.tables[0]
        if named.schema_id != expected_table_schema:
            raise ProtocolError(
                f"uniform result {self.result_schema!r} requires table schema "
                f"{expected_table_schema!r}"
            )
        _require_columns(
            named.table,
            ("aggregation", "entity", "umis", "cells", "selected_cells"),
            "uniform query counts table",
        )
        _require_types(
            named.table,
            (
                ("string", False),
                ("string", False),
                ("uint64", False),
                ("uint64", True),
                ("uint64", True),
            ),
            "uniform query counts table",
        )
        if (
            named.semantics.row_semantics is not RowSemantics.SET
            or named.semantics.key != ("aggregation", "entity")
            or named.semantics.ordered_by is not None
        ):
            raise ProtocolError(
                "uniform query counts table must be an unordered set keyed by "
                "aggregation and entity"
            )
        if not isinstance(named.selection, ExactSelection):
            raise ProtocolError(
                "uniform query counts table requires exact selection metadata"
            )

        aggregation = self.provenance.parameters.get("aggregation")
        if aggregation not in {"cell", "group", "bulk"}:
            raise ProtocolError(
                "uniform query provenance.parameters.aggregation must be cell, group, or bulk"
            )
        summary = self.summary
        assert isinstance(summary, (RegionQuerySummary, JunctionQuerySummary))
        for row_index, row in enumerate(named.rows):
            if row[0] != aggregation:
                raise ProtocolError(
                    f"uniform query counts row {row_index} disagrees with provenance aggregation"
                )
            if not row[1].strip():
                raise ProtocolError(
                    f"uniform query counts row {row_index} has an empty entity"
                )
            if aggregation == "cell" and (row[3] is not None or row[4] is not None):
                raise ProtocolError(
                    "cell count rows must leave cells and selected_cells null"
                )
            if aggregation in {"group", "bulk"}:
                if row[3] is None or row[4] is None or row[3] > row[4]:
                    raise ProtocolError(
                        f"uniform query counts row {row_index} has invalid cell totals"
                    )

        selection = named.selection
        if aggregation == "cell" and selection.available_rows != summary.cells:
            raise ProtocolError(
                "uniform query cell selection availability disagrees with summary.cells"
            )
        if aggregation != "cell" and selection.truncated:
            raise ProtocolError("group and bulk query count tables must not be truncated")
        if aggregation == "bulk":
            if len(named.rows) != 1 or named.rows[0][1] != "bulk":
                raise ProtocolError("bulk query counts must contain exactly the bulk row")
            if named.rows[0][3] != summary.cells:
                raise ProtocolError("bulk query row cells disagrees with summary.cells")
        if aggregation == "group" and sum(row[3] for row in named.rows) != summary.cells:
            raise ProtocolError("group query row cells disagree with summary.cells")
        if not selection.truncated and sum(row[2] for row in named.rows) != summary.umis:
            raise ProtocolError(
                "complete uniform query count rows disagree with summary.umis"
            )


ANNOTATION_COMPARISON_SCHEMA = "gravlax.annotation.compare.v1"
ANNOTATION_COUNT_DELTAS_SCHEMA = "gravlax.annotation.compare.count-deltas.v1"
ANNOTATION_CLASS_TRANSITIONS_SCHEMA = (
    "gravlax.annotation.compare.class-transitions.v1"
)
ANNOTATION_CONTRIBUTING_CAUSES_SCHEMA = (
    "gravlax.annotation.compare.contributing-causes.v1"
)
ANNOTATION_WITNESSES_SCHEMA = "gravlax.annotation.compare.witnesses.v1"

TRANSCRIPT_EQUIVALENCE_SCHEMA = "gravlax.query.transcript-ecs.v1"
TRANSCRIPT_EQUIVALENCE_CATALOG_SCHEMA = (
    "gravlax.query.transcript-ecs.catalog.v1"
)
TRANSCRIPT_EQUIVALENCE_COUNTS_SCHEMA = "gravlax.query.transcript-ecs.counts.v1"
TRANSCRIPT_EQUIVALENCE_MEMBERSHIP_SCHEMA = (
    "gravlax.query.transcript-ecs.membership.v1"
)


def _require_columns(table: Table, expected: tuple[str, ...], where: str) -> None:
    if table.columns != expected:
        raise ProtocolError(
            f"{where} columns {table.columns!r} do not match schema {expected!r}"
        )


def _require_types(
    table: Table, expected: tuple[tuple[str, bool], ...], where: str
) -> None:
    actual = tuple((field.data_type, field.nullable) for field in table.fields)
    if actual != expected:
        raise ProtocolError(
            f"{where} field types {actual!r} do not match schema {expected!r}"
        )


@dataclass(frozen=True)
class AnnotationComparisonResult:
    """Typed multi-table annotation comparison; transition rows are not delta attribution."""

    envelope: ResultEnvelope
    summary: Mapping[str, Any]
    semantics: Mapping[str, Any]
    count_deltas: Table
    class_transitions: Table
    contributing_causes: Table
    witnesses: Table

    @classmethod
    def from_envelope(cls, envelope: ResultEnvelope) -> "AnnotationComparisonResult":
        if envelope.result_schema != ANNOTATION_COMPARISON_SCHEMA:
            raise ProtocolError(
                f"expected {ANNOTATION_COMPARISON_SCHEMA!r}, got {envelope.result_schema!r}"
            )
        roles = {item.role for item in envelope.provenance.annotations}
        if roles != {"before", "after"}:
            raise ProtocolError(
                "annotation comparison provenance requires exactly 'before' and 'after' annotation roles"
            )
        if len(envelope.provenance.archives) != 1:
            raise ProtocolError(
                "annotation comparison provenance requires exactly one archive identity"
            )
        provenance_assembly = envelope.provenance.assembly
        annotation_assemblies = {
            item.assembly for item in envelope.provenance.annotations
        }
        if len(annotation_assemblies) != 1 or (
            provenance_assembly is not None
            and annotation_assemblies != {provenance_assembly}
        ):
            raise ProtocolError(
                "annotation comparison provenance assemblies must agree"
            )
        data = mapping(envelope.data, "annotation comparison data")
        summary = mapping(
            required(data, "summary", "annotation comparison data"),
            "annotation comparison data.summary",
        )
        semantics = mapping(
            required(data, "semantics", "annotation comparison data"),
            "annotation comparison data.semantics",
        )
        semantic_requirements = {
            "final_count_deltas_are_exact": True,
            "class_transition_ledger_is_complete": True,
            "contributing_causes_are_nonexclusive": True,
            "contributing_causes_are_additive_attributions": False,
            "annotation_order_tie_break_is_biological_change": False,
            "molecule_witnesses_are_bounded": True,
        }
        for name, expected in semantic_requirements.items():
            if semantics.get(name) is not expected:
                raise ProtocolError(
                    f"annotation comparison data.semantics.{name} must be {expected!r}"
                )
        count_deltas = Table.from_mapping(
            required(data, "count_deltas", "annotation comparison data"),
            ANNOTATION_COUNT_DELTAS_SCHEMA,
            "annotation comparison data.count_deltas",
        )
        class_transitions = Table.from_mapping(
            required(data, "class_transitions", "annotation comparison data"),
            ANNOTATION_CLASS_TRANSITIONS_SCHEMA,
            "annotation comparison data.class_transitions",
        )
        contributing_causes = Table.from_mapping(
            required(data, "contributing_causes", "annotation comparison data"),
            ANNOTATION_CONTRIBUTING_CAUSES_SCHEMA,
            "annotation comparison data.contributing_causes",
        )
        witnesses = Table.from_mapping(
            required(data, "witnesses", "annotation comparison data"),
            ANNOTATION_WITNESSES_SCHEMA,
            "annotation comparison data.witnesses",
        )
        _require_columns(
            count_deltas,
            (
                "cell",
                "cell_barcode",
                "comparison_gene_id",
                "annotation_a_gene_id",
                "annotation_b_gene_id",
                "annotation_a_count",
                "annotation_b_count",
                "signed_delta_b_minus_a",
            ),
            "annotation comparison count-deltas table",
        )
        _require_types(
            count_deltas,
            (
                ("uint64", False),
                ("string", False),
                ("string", False),
                ("string", True),
                ("string", True),
                ("uint64", False),
                ("uint64", False),
                ("int64", False),
            ),
            "annotation comparison count-deltas table",
        )
        _require_columns(
            contributing_causes,
            (
                "cell",
                "cell_barcode",
                "umi_class",
                "transition_kind",
                "contributing_cause",
                "nonexclusive",
                "additive_count_attribution",
            ),
            "annotation comparison contributing-causes table",
        )
        _require_types(
            contributing_causes,
            (
                ("uint64", False),
                ("string", False),
                ("uint64", False),
                ("string", False),
                ("string", False),
                ("boolean", False),
                ("boolean", False),
            ),
            "annotation comparison contributing-causes table",
        )
        _require_columns(
            class_transitions,
            (
                "cell",
                "cell_barcode",
                "umi_class",
                "transition_kind",
                "molecule_records",
                "evidence_rows",
                "changed_evidence_rows",
                "annotation_a_selected_comparison_gene_id",
                "annotation_a_selected_gene_id",
                "annotation_a_selected_weight",
                "annotation_a_counted",
                "annotation_a_canonical_class",
                "annotation_a_gene_support",
                "annotation_a_same_gene_neighbors",
                "annotation_b_selected_comparison_gene_id",
                "annotation_b_selected_gene_id",
                "annotation_b_selected_weight",
                "annotation_b_counted",
                "annotation_b_canonical_class",
                "annotation_b_gene_support",
                "annotation_b_same_gene_neighbors",
                "contributing_cause_count",
                "molecule_witnesses",
                "omitted_molecule_witnesses",
                "changed_row_witnesses",
                "omitted_changed_row_witnesses",
            ),
            "annotation comparison class-transitions table",
        )
        _require_types(
            class_transitions,
            (
                ("uint64", False),
                ("string", False),
                ("uint64", False),
                ("string", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("string", True),
                ("string", True),
                ("uint64", False),
                ("boolean", False),
                ("uint64", True),
                ("json", False),
                ("json", False),
                ("string", True),
                ("string", True),
                ("uint64", False),
                ("boolean", False),
                ("uint64", True),
                ("json", False),
                ("json", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
            ),
            "annotation comparison class-transitions table",
        )
        _require_columns(
            witnesses,
            (
                "archive_ordinal",
                "cell",
                "cell_barcode",
                "umi_class",
                "chrom",
                "anchor",
                "evidence_rows",
                "changed_rows_total",
                "changed_rows_omitted",
                "annotation_a_selected_comparison_gene_id",
                "annotation_a_selected_gene_id",
                "annotation_a_counted",
                "annotation_a_canonical_class",
                "annotation_b_selected_comparison_gene_id",
                "annotation_b_selected_gene_id",
                "annotation_b_counted",
                "annotation_b_canonical_class",
                "contributing_causes",
                "changed_row_witnesses",
            ),
            "annotation comparison witnesses table",
        )
        _require_types(
            witnesses,
            (
                ("uint64", False),
                ("uint64", False),
                ("string", False),
                ("uint64", False),
                ("string", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("string", True),
                ("string", True),
                ("boolean", False),
                ("uint64", True),
                ("string", True),
                ("string", True),
                ("boolean", False),
                ("uint64", True),
                ("json", False),
                ("json", False),
            ),
            "annotation comparison witnesses table",
        )
        expected_rows = {
            "count_delta_rows": len(count_deltas.rows),
            "class_transition_rows": len(class_transitions.rows),
            "molecule_witness_rows": len(witnesses.rows),
        }
        for name, expected in expected_rows.items():
            actual = integer(
                required(summary, name, "annotation comparison data.summary"),
                f"annotation comparison data.summary.{name}",
                minimum=0,
            )
            if actual != expected:
                raise ProtocolError(
                    f"annotation comparison summary {name} disagrees with its nested table"
                )
        for row_index, row in enumerate(count_deltas.rows):
            if row[6] - row[5] != row[7]:
                raise ProtocolError(
                    "annotation comparison count-deltas row "
                    f"{row_index} violates signed B-minus-A arithmetic"
                )
        for row_index, row in enumerate(contributing_causes.rows):
            if row[5] is not True or row[6] is not False:
                raise ProtocolError(
                    "annotation comparison contributing-causes row "
                    f"{row_index} must remain nonexclusive and non-additive"
                )
        return cls(
            envelope=envelope,
            summary=dict(summary),
            semantics=dict(semantics),
            count_deltas=count_deltas,
            class_transitions=class_transitions,
            contributing_causes=contributing_causes,
            witnesses=witnesses,
        )

    @classmethod
    def from_json(cls, source: str | bytes) -> "AnnotationComparisonResult":
        return cls.from_envelope(ResultEnvelope.from_json(source))


def parse_annotation_comparison(
    source: str | bytes | Mapping[str, Any],
) -> AnnotationComparisonResult:
    envelope = (
        ResultEnvelope.from_mapping(source)
        if isinstance(source, Mapping)
        else ResultEnvelope.from_json(source)
    )
    return AnnotationComparisonResult.from_envelope(envelope)


@dataclass(frozen=True)
class TranscriptEquivalenceResult:
    """Annotation-conditional compatibility classes, never abundance or phasing calls."""

    envelope: ResultEnvelope
    scope: Mapping[str, Any]
    semantics: Mapping[str, Any]
    summary: Mapping[str, Any]
    catalog: Table
    counts: Table
    membership: Optional[Table]

    @classmethod
    def from_envelope(cls, envelope: ResultEnvelope) -> "TranscriptEquivalenceResult":
        if envelope.result_schema != TRANSCRIPT_EQUIVALENCE_SCHEMA:
            raise ProtocolError(
                f"expected {TRANSCRIPT_EQUIVALENCE_SCHEMA!r}, got {envelope.result_schema!r}"
            )
        provenance = envelope.provenance
        if len(provenance.archives) > 1:
            raise ProtocolError(
                "transcript equivalence provenance permits at most one archive identity"
            )
        if (
            provenance.assembly is None
            or provenance.annotation is None
            or provenance.annotation_digest is None
        ):
            raise ProtocolError(
                "transcript equivalence provenance requires assembly, annotation, and annotation_digest"
            )
        data = mapping(envelope.data, "transcript equivalence data")
        scope = mapping(
            required(data, "scope", "transcript equivalence data"),
            "transcript equivalence data.scope",
        )
        semantics = mapping(
            required(data, "semantics", "transcript equivalence data"),
            "transcript equivalence data.semantics",
        )
        compatibility = mapping(
            required(
                semantics,
                "compatibility",
                "transcript equivalence data.semantics",
            ),
            "transcript equivalence data.semantics.compatibility",
        )
        for name in ("abundance_inferred", "full_isoform_phasing_claimed"):
            if compatibility.get(name) is not False:
                raise ProtocolError(
                    "transcript equivalence data.semantics.compatibility."
                    f"{name} must be False"
                )
        summary = mapping(
            required(data, "summary", "transcript equivalence data"),
            "transcript equivalence data.summary",
        )
        catalog = Table.from_mapping(
            required(data, "catalog", "transcript equivalence data"),
            TRANSCRIPT_EQUIVALENCE_CATALOG_SCHEMA,
            "transcript equivalence data.catalog",
        )
        counts = Table.from_mapping(
            required(data, "counts", "transcript equivalence data"),
            TRANSCRIPT_EQUIVALENCE_COUNTS_SCHEMA,
            "transcript equivalence data.counts",
        )
        membership_value = data.get("membership")
        membership = (
            None
            if membership_value is None
            else Table.from_mapping(
                membership_value,
                TRANSCRIPT_EQUIVALENCE_MEMBERSHIP_SCHEMA,
                "transcript equivalence data.membership",
            )
        )
        _require_columns(
            catalog,
            (
                "ec_id",
                "transcript_ids",
                "gene_ids",
                "ambiguous",
                "archived_umi_class_count",
                "cell_count",
                "complete_umi_class_count",
            ),
            "transcript equivalence catalog table",
        )
        _require_types(
            catalog,
            (
                ("string", False),
                ("json", False),
                ("json", False),
                ("boolean", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
            ),
            "transcript equivalence catalog table",
        )
        _require_columns(
            counts,
            (
                "aggregation",
                "key",
                "cell_id",
                "ec_id",
                "archived_umi_class_count",
                "ambiguous_umi_class_count",
                "no_compatible_transcript_umi_class_count",
                "conflicting_umi_class_count",
                "complete_umi_class_count",
                "incomplete_umi_class_count",
                "retained_record_count",
                "represented_alignment_count",
            ),
            "transcript equivalence counts table",
        )
        _require_types(
            counts,
            (
                ("string", False),
                ("string", False),
                ("uint64", True),
                ("string", True),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
                ("uint64", False),
            ),
            "transcript equivalence counts table",
        )
        if membership is not None:
            _require_columns(
                membership,
                (
                    "umi_class",
                    "cell_id",
                    "barcode",
                    "aggregation",
                    "key",
                    "ec_id",
                    "retained_record_count",
                    "represented_alignment_count",
                    "compatible_record_count",
                    "unmatched_record_count",
                    "ambiguous",
                    "no_compatible_transcript",
                    "conflict",
                    "complete_within_archive_quotient",
                    "retained_representatives_complete",
                ),
                "transcript equivalence membership table",
            )
            _require_types(
                membership,
                (
                    ("uint64", False),
                    ("uint64", False),
                    ("string", False),
                    ("string", False),
                    ("string", False),
                    ("string", True),
                    ("uint64", False),
                    ("uint64", False),
                    ("uint64", False),
                    ("uint64", False),
                    ("boolean", False),
                    ("boolean", False),
                    ("boolean", False),
                    ("boolean", False),
                    ("boolean", False),
                ),
                "transcript equivalence membership table",
            )
        summary_counts = {
            name: integer(
                required(summary, name, "transcript equivalence data.summary"),
                f"transcript equivalence data.summary.{name}",
                minimum=0,
            )
            for name in (
                "scoped_umi_classes",
                "transcript_ecs",
                "count_rows",
                "membership_rows",
                "assigned_umi_classes",
                "unassigned_umi_classes",
                "complete_umi_classes",
                "incomplete_umi_classes",
            )
        }
        if (
            summary_counts["assigned_umi_classes"]
            + summary_counts["unassigned_umi_classes"]
            != summary_counts["scoped_umi_classes"]
            or summary_counts["complete_umi_classes"]
            + summary_counts["incomplete_umi_classes"]
            != summary_counts["scoped_umi_classes"]
        ):
            raise ProtocolError(
                "transcript equivalence summary violates UMI-class conservation"
            )
        if (
            len(catalog.rows) != summary_counts["transcript_ecs"]
            or len(counts.rows) != summary_counts["count_rows"]
            or (0 if membership is None else len(membership.rows))
            != summary_counts["membership_rows"]
        ):
            raise ProtocolError(
                "transcript equivalence summary row counts disagree with nested tables"
            )
        emit_membership = provenance.parameters.get("emit_membership")
        if not isinstance(emit_membership, bool):
            raise ProtocolError(
                "transcript equivalence provenance.parameters.emit_membership must be a boolean"
            )
        if emit_membership != (membership is not None):
            raise ProtocolError(
                "transcript equivalence membership table disagrees with emit_membership provenance"
            )
        caps = {
            name: integer(
                required(
                    provenance.parameters,
                    name,
                    "transcript equivalence provenance.parameters",
                ),
                f"transcript equivalence provenance.parameters.{name}",
                minimum=1,
            )
            for name in ("max_ecs", "max_memberships", "max_count_rows")
        }
        if (
            summary_counts["transcript_ecs"] > caps["max_ecs"]
            or summary_counts["count_rows"] > caps["max_count_rows"]
            or summary_counts["membership_rows"] > caps["max_memberships"]
        ):
            raise ProtocolError(
                "transcript equivalence result exceeds a disclosed fail-closed row cap"
            )
        for row_index, row in enumerate(catalog.rows):
            if row[6] > row[4]:
                raise ProtocolError(
                    "transcript equivalence catalog row "
                    f"{row_index} has more complete than archived UMI classes"
                )
        for row_index, row in enumerate(counts.rows):
            if row[8] + row[9] != row[4]:
                raise ProtocolError(
                    "transcript equivalence counts row "
                    f"{row_index} violates complete/incomplete conservation"
                )
        return cls(
            envelope=envelope,
            scope=dict(scope),
            semantics=dict(semantics),
            summary=dict(summary),
            catalog=catalog,
            counts=counts,
            membership=membership,
        )

    @classmethod
    def from_json(cls, source: str | bytes) -> "TranscriptEquivalenceResult":
        return cls.from_envelope(ResultEnvelope.from_json(source))


def parse_transcript_equivalence(
    source: str | bytes | Mapping[str, Any],
) -> TranscriptEquivalenceResult:
    envelope = (
        ResultEnvelope.from_mapping(source)
        if isinstance(source, Mapping)
        else ResultEnvelope.from_json(source)
    )
    return TranscriptEquivalenceResult.from_envelope(envelope)


def parse_result(source: str | bytes | Mapping[str, Any]) -> ResultEnvelope:
    """Parse a result envelope from JSON text, bytes, or an already-decoded mapping."""

    if isinstance(source, Mapping):
        return ResultEnvelope.from_mapping(source)
    return ResultEnvelope.from_json(source)


def parse_uniform_bundle(
    source: str | bytes | Mapping[str, Any],
) -> UniformResultBundle:
    """Parse and validate a named-table uniform JSON result bundle."""

    if isinstance(source, Mapping):
        return UniformResultBundle.from_mapping(source)
    return UniformResultBundle.from_json(source)
