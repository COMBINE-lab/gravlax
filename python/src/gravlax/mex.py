"""Guarded reader and AnnData adapter for Gravlax MEX result directories."""

from __future__ import annotations

import math
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, BinaryIO

from ._json import integer, mapping, nonempty_string, required
from .exceptions import ProtocolError
from .results import ResultEnvelope, _optional

FileIdentity = tuple[int, str]


def _open_file_identity(path: Path) -> FileIdentity:
    try:
        with path.open("rb") as handle:
            return _handle_identity(handle)
    except OSError as error:
        raise ProtocolError(f"cannot read MEX component {path}: {error}") from error


def _handle_identity(handle: BinaryIO) -> FileIdentity:
    hasher = hashlib.sha256()
    total = 0
    handle.seek(0)
    while True:
        chunk = handle.read(1024 * 1024)
        if not chunk:
            break
        total += len(chunk)
        hasher.update(chunk)
    handle.seek(0)
    return total, hasher.hexdigest()


def _safe_result_file(root: Path, value: Any, where: str) -> Path:
    name = nonempty_string(value, where)
    relative = Path(name)
    if relative.is_absolute() or len(relative.parts) != 1 or relative.name in {".", ".."}:
        raise ProtocolError(f"{where} must name one file directly inside the MEX directory")
    candidate = root / relative
    try:
        resolved = candidate.resolve(strict=True)
    except OSError as error:
        raise ProtocolError(f"cannot resolve {where} file {candidate}: {error}") from error
    if resolved.parent != root:
        raise ProtocolError(f"{where} escapes the MEX directory through a symlink")
    if not resolved.is_file():
        raise ProtocolError(f"{where} is not a file: {resolved}")
    return resolved


def _unescape(value: str, where: str) -> str:
    result: list[str] = []
    index = 0
    escapes = {"\\": "\\", "t": "\t", "r": "\r", "n": "\n"}
    while index < len(value):
        char = value[index]
        if char != "\\":
            result.append(char)
            index += 1
            continue
        if index + 1 >= len(value) or value[index + 1] not in escapes:
            raise ProtocolError(f"{where} contains an invalid TSV escape")
        result.append(escapes[value[index + 1]])
        index += 2
    return "".join(result)


@dataclass(frozen=True)
class MexFeature:
    id: str
    name: str
    feature_type: str


@dataclass(frozen=True)
class MexDataset:
    directory: Path
    envelope: ResultEnvelope
    matrix_path: Path
    features_path: Path
    barcodes_path: Path
    feature_count: int
    barcode_count: int
    nonzero_count: int
    value_type: str
    _matrix_identity: FileIdentity = field(repr=False)
    _features_identity: FileIdentity = field(repr=False)
    _barcodes_identity: FileIdentity = field(repr=False)

    @classmethod
    def open(cls, directory: str | Path) -> "MexDataset":
        try:
            root = Path(directory).resolve(strict=True)
        except OSError as error:
            raise ProtocolError(f"cannot open MEX directory {directory}: {error}") from error
        if not root.is_dir():
            raise ProtocolError(f"MEX path is not a directory: {root}")
        metadata_path = root / "metadata.json"
        if not metadata_path.is_file():
            raise ProtocolError(
                f"MEX directory is incomplete: completion marker {metadata_path} is missing"
            )
        metadata_path = _safe_result_file(root, "metadata.json", "MEX completion marker")
        envelope = ResultEnvelope.from_file(metadata_path)
        data = mapping(envelope.data, "MEX metadata.data")
        format_name = nonempty_string(
            required(data, "format", "MEX metadata.data"), "MEX metadata.data.format"
        )
        if format_name != "matrix_market_coordinate":
            raise ProtocolError(f"unsupported MEX format {format_name!r}")
        index_base = integer(
            required(data, "index_base", "MEX metadata.data"),
            "MEX metadata.data.index_base",
        )
        if index_base != 1:
            raise ProtocolError(f"unsupported MEX index base {index_base}; expected 1")
        value_type = nonempty_string(
            required(data, "value_type", "MEX metadata.data"),
            "MEX metadata.data.value_type",
        )
        if value_type not in {"integer", "real"}:
            raise ProtocolError(f"unsupported MEX value type {value_type!r}")
        matrix_path = _safe_result_file(
            root,
            required(data, "matrix", "MEX metadata.data"),
            "MEX metadata.data.matrix",
        )
        features_path = _safe_result_file(
            root,
            required(data, "features", "MEX metadata.data"),
            "MEX metadata.data.features",
        )
        barcodes_path = _safe_result_file(
            root,
            required(data, "barcodes", "MEX metadata.data"),
            "MEX metadata.data.barcodes",
        )
        dataset = cls(
            directory=root,
            envelope=envelope,
            matrix_path=matrix_path,
            features_path=features_path,
            barcodes_path=barcodes_path,
            feature_count=integer(
                required(data, "feature_count", "MEX metadata.data"),
                "MEX metadata.data.feature_count",
                minimum=0,
            ),
            barcode_count=integer(
                required(data, "barcode_count", "MEX metadata.data"),
                "MEX metadata.data.barcode_count",
                minimum=0,
            ),
            nonzero_count=integer(
                required(data, "nonzero_count", "MEX metadata.data"),
                "MEX metadata.data.nonzero_count",
                minimum=0,
            ),
            value_type=value_type,
            _matrix_identity=_open_file_identity(matrix_path),
            _features_identity=_open_file_identity(features_path),
            _barcodes_identity=_open_file_identity(barcodes_path),
        )
        dataset._validate_matrix()
        if len(dataset.features()) != dataset.feature_count:
            raise ProtocolError("features.tsv row count does not match MEX metadata")
        if len(dataset.barcodes()) != dataset.barcode_count:
            raise ProtocolError("barcodes.tsv row count does not match MEX metadata")
        return dataset

    def _validate_matrix(self) -> None:
        """Validate dimensions and every coordinate without SciPy.

        Gravlax emits coordinates in strict feature/barcode order. Enforcing
        that canonical order detects duplicates with constant memory while
        also rejecting a reordered file that could not have come from the
        completed writer represented by ``metadata.json``.
        """

        try:
            with self.matrix_path.open("rb") as handle:
                def line_text(raw: bytes) -> str:
                    try:
                        return raw.decode("utf-8")
                    except UnicodeDecodeError as error:
                        raise ProtocolError("matrix.mtx is not UTF-8 text") from error

                header = line_text(handle.readline()).strip().split()
                expected = [
                    "%%MatrixMarket",
                    "matrix",
                    "coordinate",
                    self.value_type,
                    "general",
                ]
                if header != expected:
                    raise ProtocolError(
                        "matrix.mtx header does not match the guarded MEX metadata"
                    )
                for raw_line in handle:
                    line = line_text(raw_line)
                    if line.startswith("%") or not line.strip():
                        continue
                    dimensions = line.split()
                    break
                else:
                    raise ProtocolError("matrix.mtx has no dimensions line")
                if len(dimensions) != 3:
                    raise ProtocolError(
                        "matrix.mtx dimensions line must have three integers"
                    )
                try:
                    rows, columns, declared_entries = (
                        int(value, 10) for value in dimensions
                    )
                except ValueError as error:
                    raise ProtocolError("matrix.mtx dimensions are not integers") from error
                if (rows, columns, declared_entries) != (
                    self.feature_count,
                    self.barcode_count,
                    self.nonzero_count,
                ):
                    raise ProtocolError("matrix.mtx dimensions do not match MEX metadata")

                previous: tuple[int, int] | None = None
                observed_entries = 0
                for line_number, raw_line in enumerate(handle, start=3):
                    line = line_text(raw_line)
                    if line.startswith("%") or not line.strip():
                        continue
                    values = line.split()
                    if len(values) != 3:
                        raise ProtocolError(
                            f"matrix.mtx line {line_number} must have row, column, and value"
                        )
                    try:
                        feature = int(values[0], 10)
                        barcode = int(values[1], 10)
                    except ValueError as error:
                        raise ProtocolError(
                            f"matrix.mtx line {line_number} has a non-integer coordinate"
                        ) from error
                    if not 1 <= feature <= self.feature_count:
                        raise ProtocolError(
                            f"matrix.mtx line {line_number} feature index is out of bounds"
                        )
                    if not 1 <= barcode <= self.barcode_count:
                        raise ProtocolError(
                            f"matrix.mtx line {line_number} barcode index is out of bounds"
                        )
                    coordinate = (feature, barcode)
                    if previous is not None and coordinate <= previous:
                        if coordinate == previous:
                            raise ProtocolError(
                                f"matrix.mtx line {line_number} duplicates a coordinate"
                            )
                        raise ProtocolError(
                            "matrix.mtx coordinates are not in canonical feature/barcode order"
                        )
                    previous = coordinate

                    if self.value_type == "integer":
                        try:
                            value = int(values[2], 10)
                        except ValueError as error:
                            raise ProtocolError(
                                f"matrix.mtx line {line_number} has a non-integer value"
                            ) from error
                        if value == 0:
                            raise ProtocolError(
                                f"matrix.mtx line {line_number} stores an explicit zero"
                            )
                    else:
                        try:
                            value = float(values[2])
                        except (ValueError, OverflowError) as error:
                            raise ProtocolError(
                                f"matrix.mtx line {line_number} has an invalid real value"
                            ) from error
                        if value == 0.0 or not math.isfinite(value):
                            raise ProtocolError(
                                f"matrix.mtx line {line_number} stores a zero or non-finite value"
                            )
                    observed_entries += 1
                    if observed_entries > self.nonzero_count:
                        raise ProtocolError(
                            "matrix.mtx contains more coordinates than MEX metadata declares"
                        )
                if _handle_identity(handle) != self._matrix_identity:
                    raise ProtocolError(
                        "matrix.mtx changed after the MEX dataset was opened"
                    )
        except OSError as error:
            raise ProtocolError(f"cannot read matrix.mtx: {error}") from error
        if observed_entries != self.nonzero_count:
            raise ProtocolError(
                "matrix.mtx is truncated: coordinate count does not match MEX metadata"
            )

    def features(self) -> tuple[MexFeature, ...]:
        values: list[MexFeature] = []
        seen_ids: set[str] = set()
        try:
            with self.features_path.open("rb") as handle:
                for index, raw_line in enumerate(handle):
                    try:
                        line = raw_line.decode("utf-8")
                    except UnicodeDecodeError as error:
                        raise ProtocolError("features.tsv is not UTF-8 text") from error
                    columns = line.rstrip("\r\n").split("\t")
                    if len(columns) != 3:
                        raise ProtocolError(
                            f"features.tsv line {index + 1} does not have three columns"
                        )
                    feature = MexFeature(
                        id=_unescape(columns[0], f"features.tsv line {index + 1}"),
                        name=_unescape(columns[1], f"features.tsv line {index + 1}"),
                        feature_type=_unescape(
                            columns[2], f"features.tsv line {index + 1}"
                        ),
                    )
                    if not feature.id:
                        raise ProtocolError(
                            f"features.tsv line {index + 1} has an empty feature id"
                        )
                    if feature.id in seen_ids:
                        raise ProtocolError(
                            f"features.tsv line {index + 1} has duplicate feature id {feature.id!r}"
                        )
                    seen_ids.add(feature.id)
                    values.append(feature)
                if _handle_identity(handle) != self._features_identity:
                    raise ProtocolError(
                        "features.tsv changed after the MEX dataset was opened"
                    )
        except OSError as error:
            raise ProtocolError(f"cannot read features.tsv: {error}") from error
        return tuple(values)

    def barcodes(self) -> tuple[str, ...]:
        values: list[str] = []
        seen_barcodes: set[str] = set()
        try:
            with self.barcodes_path.open("rb") as handle:
                for index, raw_line in enumerate(handle):
                    try:
                        line = raw_line.decode("utf-8")
                    except UnicodeDecodeError as error:
                        raise ProtocolError("barcodes.tsv is not UTF-8 text") from error
                    encoded = line.rstrip("\r\n")
                    if "\t" in encoded:
                        raise ProtocolError(
                            f"barcodes.tsv line {index + 1} has more than one column"
                        )
                    barcode = _unescape(encoded, f"barcodes.tsv line {index + 1}")
                    if not barcode:
                        raise ProtocolError(
                            f"barcodes.tsv line {index + 1} has an empty barcode"
                        )
                    if barcode in seen_barcodes:
                        raise ProtocolError(
                            f"barcodes.tsv line {index + 1} has duplicate barcode {barcode!r}"
                        )
                    seen_barcodes.add(barcode)
                    values.append(barcode)
                if _handle_identity(handle) != self._barcodes_identity:
                    raise ProtocolError(
                        "barcodes.tsv changed after the MEX dataset was opened"
                    )
        except OSError as error:
            raise ProtocolError(f"cannot read barcodes.tsv: {error}") from error
        return tuple(values)

    def to_scipy(self) -> Any:
        """Return the exact feature-by-barcode matrix as SciPy CSR."""

        scipy_io = _optional("scipy.io", "anndata")
        try:
            with self.matrix_path.open("rb") as handle:
                if _handle_identity(handle) != self._matrix_identity:
                    raise ProtocolError(
                        "matrix.mtx changed after the MEX dataset was opened"
                    )
                decoded = scipy_io.mmread(handle)
                if _handle_identity(handle) != self._matrix_identity:
                    raise ProtocolError(
                        "matrix.mtx changed while it was being decoded"
                    )
        except OSError as error:
            raise ProtocolError(f"cannot read matrix.mtx: {error}") from error
        if decoded.shape != (self.feature_count, self.barcode_count):
            raise ProtocolError("decoded matrix shape does not match MEX metadata")
        if decoded.nnz != self.nonzero_count:
            raise ProtocolError("decoded matrix nonzero count does not match MEX metadata")
        return decoded.tocsr()

    def to_anndata(self) -> Any:
        """Return cells by features, with labels and provenance preserved."""

        anndata = _optional("anndata", "anndata")
        pandas = _optional("pandas", "anndata")
        features = self.features()
        barcodes = self.barcodes()
        matrix = self.to_scipy().transpose().tocsr()
        obs = pandas.DataFrame(index=pandas.Index(barcodes, name="barcode"))
        var = pandas.DataFrame(
            {
                "feature_name": [feature.name for feature in features],
                "feature_type": [feature.feature_type for feature in features],
            },
            index=pandas.Index([feature.id for feature in features], name="feature_id"),
        )
        result = anndata.AnnData(X=matrix, obs=obs, var=var)
        result.uns["gravlax"] = self.envelope.metadata()
        return result


def read_mex(directory: str | Path) -> MexDataset:
    """Open and validate a completed Gravlax MEX result directory."""

    return MexDataset.open(directory)
