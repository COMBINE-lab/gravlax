"""Safe subprocess client for Gravlax's versioned JSON protocols."""

from __future__ import annotations

import ctypes
import errno
import os
import secrets
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Union

from ._json import loads_document
from .exceptions import (
    CommandError,
    CommandTimeoutError,
    ExecutableNotFoundError,
    ProtocolError,
)
from .models import DoctorReport, Project, ResolvedPlan, StepCompletion
from .results import (
    AnnotationComparisonResult,
    ResultEnvelope,
    TranscriptEquivalenceResult,
    UniformResultBundle,
)

PathToken = Union[str, os.PathLike[str]]
BinaryCommand = Union[PathToken, Sequence[PathToken]]
_STALE_FILE_ERROR = getattr(errno, "ESTALE", errno.EIO)


def _token(value: PathToken, where: str) -> str:
    try:
        result = os.fspath(value)
    except TypeError as error:
        raise TypeError(f"{where} must be a string or path-like value") from error
    if isinstance(result, bytes):
        result = os.fsdecode(result)
    if not result:
        raise ValueError(f"{where} must not be empty")
    if "\0" in result:
        raise ValueError(f"{where} must not contain a NUL byte")
    return result


def _binary_prefix(value: BinaryCommand) -> tuple[str, ...]:
    if isinstance(value, (str, os.PathLike)):
        return (_token(value, "binary"),)
    if not isinstance(value, Sequence) or not value:
        raise ValueError("binary command must contain at least one argument")
    return tuple(_token(item, f"binary[{index}]") for index, item in enumerate(value))


def _args(values: Sequence[PathToken]) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, os.PathLike)):
        raise TypeError("command arguments must be a sequence of tokens, not a command string")
    return tuple(_token(value, f"args[{index}]") for index, value in enumerate(values))


def _option(name: str, value: PathToken) -> str:
    """Use ``--name=value`` so a leading hyphen in a path stays data."""

    return f"--{name}={_token(value, name)}"


def _choice(value: str, allowed: set[str], where: str) -> str:
    if value not in allowed:
        choices = ", ".join(sorted(allowed))
        raise ValueError(f"{where} must be one of {choices}")
    return value


def _positive(value: int, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{where} must be a positive integer")
    return value


def _nonnegative(value: int, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{where} must be a nonnegative integer")
    return value


def _uniform_count_query_args(
    archive: PathToken,
    command: str,
    locus: str,
    *,
    top: int,
    cells: Optional[PathToken],
    groups: Optional[PathToken],
    aggregation: str,
) -> list[PathToken]:
    if cells is not None and groups is not None:
        raise ValueError("cells and groups are mutually exclusive")
    aggregation = _choice(
        aggregation, {"auto", "cell", "group", "bulk"}, "aggregation"
    )
    if aggregation == "group" and groups is None:
        raise ValueError("group aggregation requires groups")
    top = _nonnegative(top, "top")
    args: list[PathToken] = [
        "query",
        _token(archive, "archive"),
        command,
        _token(locus, "locus"),
        _option("top", str(top)),
        "--format=json",
    ]
    if cells is not None:
        args.append(_option("cells", cells))
    if groups is not None:
        args.append(_option("groups", groups))
    if aggregation != "auto":
        args.append(_option("agg", aggregation))
    return args


def _cooccurrence_query_args(
    archive: PathToken,
    predicates: Mapping[str, str],
    expression: str,
    universe: str,
    *,
    unit: str,
    region_match: str,
    placements: str,
    allow_full_scan: bool,
    cells: Optional[PathToken],
    groups: Optional[PathToken],
    aggregation: str,
    emit_membership: bool,
    max_memberships: int,
    max_pattern_rows: int,
    max_chunks: int,
    max_evidence_records: int,
    max_terminal_events: int,
) -> list[PathToken]:
    if not predicates:
        raise ValueError("predicates must contain at least one named predicate")
    if not isinstance(expression, str) or not expression.strip():
        raise ValueError("expression must not be empty")
    if not isinstance(universe, str) or not universe.strip():
        raise ValueError("universe must not be empty")
    if universe not in predicates:
        raise ValueError("universe must name one of predicates")
    if cells is not None and groups is not None:
        raise ValueError("cells and groups are mutually exclusive")
    aggregation = _choice(
        aggregation, {"auto", "cell", "group", "bulk"}, "aggregation"
    )
    if aggregation == "group" and groups is None:
        raise ValueError("group aggregation requires groups")
    unit = _choice(unit, {"molecule-record", "umi-class"}, "unit")
    region_match = _choice(
        region_match, {"anchor", "aligned-block"}, "region_match"
    )
    placements = _choice(placements, {"unique", "direct", "all"}, "placements")
    if unit == "umi-class" and not allow_full_scan:
        raise ValueError("unit='umi-class' requires allow_full_scan=True")
    max_memberships = _positive(max_memberships, "max_memberships")
    max_pattern_rows = _positive(max_pattern_rows, "max_pattern_rows")
    max_chunks = _positive(max_chunks, "max_chunks")
    max_evidence_records = _positive(
        max_evidence_records, "max_evidence_records"
    )
    max_terminal_events = _nonnegative(
        max_terminal_events, "max_terminal_events"
    )
    args: list[PathToken] = [
        "query",
        _token(archive, "archive"),
        "cooccur",
    ]
    for name, descriptor in predicates.items():
        if not isinstance(name, str) or not name:
            raise ValueError("predicate names must be non-empty strings")
        if not isinstance(descriptor, str) or not descriptor:
            raise ValueError(f"predicate {name!r} must have a non-empty descriptor")
        args.append(_option("predicate", f"{name}={descriptor}"))
    args.extend(
        [
            _option("where", expression),
            _option("universe", universe),
            _option("unit", unit),
            _option("region-match", region_match),
            _option("placements", placements),
            _option("max-pattern-rows", str(max_pattern_rows)),
            _option("max-chunks", str(max_chunks)),
            _option("max-evidence-records", str(max_evidence_records)),
            _option("max-terminal-events", str(max_terminal_events)),
            "--format=json",
        ]
    )
    if allow_full_scan:
        args.append("--allow-full-scan")
    if emit_membership:
        args.append("--emit-membership")
        args.append(_option("max-memberships", str(max_memberships)))
    if cells is not None:
        args.append(_option("cells", cells))
    if groups is not None:
        args.append(_option("groups", groups))
    if aggregation != "auto":
        args.append(_option("agg", aggregation))
    return args


def _collection_find_events_args(
    collection: PathToken,
    *,
    kinds: Sequence[str],
    design: Optional[PathToken],
    groups: Optional[PathToken],
    require_groups: Sequence[str],
    min_group_umi_classes: int,
    min_donors: int,
    min_samples: int,
    min_umi_classes: int,
    min_side_umi_classes: int,
    min_support: int,
    terminal_cluster_bp: int,
    max_terminal_events: int,
    annotation: Optional[PathToken],
    assembly: Optional[str],
    annotation_label: Optional[str],
    annotation_digest: Optional[str],
    novel_only: bool,
    solo_strand: str,
    max_candidates: int,
    max_candidates_considered: int,
    max_routed_entries: int,
    max_exact_match_attempts: int,
    max_annotation_comparisons: int,
    verify_content: bool,
) -> list[PathToken]:
    selected_kinds = tuple(
        _choice(
            kind,
            {
                "junction",
                "alt-acceptor",
                "alt-donor",
                "cassette",
                "terminal-tail",
            },
            "kind",
        )
        for kind in kinds
    )
    if len(selected_kinds) != len(set(selected_kinds)):
        raise ValueError("kinds must not contain duplicates")
    selected_groups = _args(require_groups)
    if selected_groups and groups is None:
        raise ValueError("require_groups requires groups")
    if len(selected_groups) != len(set(selected_groups)):
        raise ValueError("require_groups must not contain duplicates")
    if novel_only and annotation is None:
        raise ValueError("novel_only requires annotation")
    annotation_identity = (assembly, annotation_label, annotation_digest)
    if annotation is None and any(value is not None for value in annotation_identity):
        raise ValueError(
            "assembly, annotation_label, and annotation_digest require annotation"
        )
    if annotation is not None:
        if not isinstance(assembly, str) or not assembly.strip():
            raise ValueError("annotation requires a non-empty assembly")
        if not isinstance(annotation_label, str) or not annotation_label.strip():
            raise ValueError("annotation requires a non-empty annotation_label")
    solo_strand = _choice(
        solo_strand,
        {"forward", "reverse", "unstranded"},
        "solo_strand",
    )
    terminal_requested = not selected_kinds or "terminal-tail" in selected_kinds
    if terminal_requested and solo_strand != "forward":
        raise ValueError("terminal-tail discovery requires solo_strand='forward'")
    min_group_umi_classes = _positive(
        min_group_umi_classes, "min_group_umi_classes"
    )
    min_donors = _positive(min_donors, "min_donors")
    min_samples = _positive(min_samples, "min_samples")
    min_umi_classes = _positive(min_umi_classes, "min_umi_classes")
    min_side_umi_classes = _positive(
        min_side_umi_classes, "min_side_umi_classes"
    )
    min_support = _nonnegative(min_support, "min_support")
    terminal_cluster_bp = _nonnegative(
        terminal_cluster_bp, "terminal_cluster_bp"
    )
    max_terminal_events = _nonnegative(
        max_terminal_events, "max_terminal_events"
    )
    max_candidates = _positive(max_candidates, "max_candidates")
    max_candidates_considered = _positive(
        max_candidates_considered, "max_candidates_considered"
    )
    max_routed_entries = _positive(max_routed_entries, "max_routed_entries")
    max_exact_match_attempts = _positive(
        max_exact_match_attempts, "max_exact_match_attempts"
    )
    max_annotation_comparisons = _positive(
        max_annotation_comparisons, "max_annotation_comparisons"
    )
    args: list[PathToken] = [
        "collection",
        "find-events",
        _token(collection, "collection"),
    ]
    for kind in selected_kinds:
        args.append(_option("kind", kind))
    if design is not None:
        args.append(_option("design", design))
    if groups is not None:
        args.append(_option("groups", groups))
    for group in selected_groups:
        args.append(_option("require-group", group))
    if selected_groups:
        args.append(
            _option("min-group-umi-classes", str(min_group_umi_classes))
        )
    args.extend(
        [
            _option("min-donors", str(min_donors)),
            _option("min-samples", str(min_samples)),
            _option("min-umi-classes", str(min_umi_classes)),
            _option("min-side-umi-classes", str(min_side_umi_classes)),
            _option("min-support", str(min_support)),
            _option("terminal-cluster-bp", str(terminal_cluster_bp)),
            _option("max-terminal-events", str(max_terminal_events)),
            _option("solo-strand", solo_strand),
            _option("max-candidates", str(max_candidates)),
            _option("max-candidates-considered", str(max_candidates_considered)),
            _option("max-routed-entries", str(max_routed_entries)),
            _option("max-exact-match-attempts", str(max_exact_match_attempts)),
            _option("max-annotation-comparisons", str(max_annotation_comparisons)),
        ]
    )
    if annotation is not None:
        args.append(_option("annotation", annotation))
        args.append(_option("assembly", assembly))
        args.append(_option("annotation-label", annotation_label))
        if annotation_digest is not None:
            args.append(_option("annotation-digest", annotation_digest))
    if novel_only:
        args.append("--novel-only")
    if verify_content:
        args.append("--verify-content")
    args.append("--format=json")
    return args


def _same_file_identity(left: os.stat_result, right: os.stat_result) -> bool:
    return (left.st_dev, left.st_ino) == (right.st_dev, right.st_ino)


def _path_matches_file(
    path: Path, expected: os.stat_result
) -> bool:
    try:
        observed = os.stat(path, follow_symlinks=False)
    except FileNotFoundError:
        return False
    return _same_file_identity(observed, expected)


def _verify_path_matches_descriptor(
    path: Path, descriptor: int, where: str
) -> os.stat_result:
    expected = os.fstat(descriptor)
    try:
        observed = os.stat(path, follow_symlinks=False)
    except FileNotFoundError as error:
        raise OSError(
            _STALE_FILE_ERROR,
            f"{where} disappeared before its identity could be verified",
            path,
        ) from error
    if not _same_file_identity(observed, expected):
        raise OSError(
            _STALE_FILE_ERROR,
            f"{where} does not refer to the completed output file",
            path,
        )
    return observed


def _link_descriptor_exact(descriptor: int, destination: Path) -> None:
    """Link the inode held by ``descriptor`` without trusting its staging name."""

    if not sys.platform.startswith("linux") or not Path("/proc/self/fd").is_dir():
        raise OSError(errno.ENOSYS, "descriptor linking is unavailable")

    # os.link(..., follow_symlinks=True) may still call link(2), which treats the
    # procfs magic link as a different filesystem. linkat(2) with
    # AT_SYMLINK_FOLLOW reliably dereferences /proc/self/fd/N to the held inode.
    at_fdcwd = -100
    at_symlink_follow = 0x400
    try:
        linkat = ctypes.CDLL(None, use_errno=True).linkat
    except (AttributeError, OSError) as error:
        raise OSError(errno.ENOSYS, "linkat is unavailable") from error
    linkat.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
    ]
    linkat.restype = ctypes.c_int
    source = os.fsencode(f"/proc/self/fd/{descriptor}")
    target = os.fsencode(destination)
    if linkat(
        at_fdcwd,
        source,
        at_fdcwd,
        target,
        at_symlink_follow,
    ) != 0:
        error_number = ctypes.get_errno()
        raise OSError(error_number, os.strerror(error_number), destination)


def _link_held_descriptor(
    descriptor: int, staging: Path, destination: Path
) -> None:
    """Create a hard link to an open output, falling back only with detection."""

    unsupported = {
        errno.EACCES,
        errno.EINVAL,
        errno.ENOENT,
        errno.ENOSYS,
        errno.EPERM,
        errno.EXDEV,
    }
    if hasattr(errno, "EOPNOTSUPP"):
        unsupported.add(errno.EOPNOTSUPP)
    try:
        _link_descriptor_exact(descriptor, destination)
    except FileExistsError:
        raise
    except OSError as error:
        if error.errno not in unsupported:
            raise
        expected = os.fstat(descriptor)
        if not _path_matches_file(staging, expected):
            raise OSError(
                _STALE_FILE_ERROR,
                "staging path changed and exact descriptor linking is unavailable",
                staging,
            ) from error
        try:
            os.link(staging, destination, follow_symlinks=False)
        except (TypeError, NotImplementedError):
            # Older platforms may not expose follow_symlinks. The post-link
            # inode check below detects a source-name substitution.
            os.link(staging, destination)
    _verify_path_matches_descriptor(
        destination, descriptor, "newly linked output"
    )


def _link_descriptor_to_unique_path(
    descriptor: int, staging: Path, destination: Path
) -> Path:
    for _ in range(16):
        candidate = destination.parent / (
            f".gravlax-{secrets.token_hex(16)}.publish"
        )
        try:
            _link_held_descriptor(descriptor, staging, candidate)
        except FileExistsError:
            continue
        return candidate
    raise FileExistsError(
        f"cannot reserve a private publication name beside {destination}"
    )


def _unlink_if_same_file(path: Path, expected: os.stat_result) -> None:
    """Remove only a name still referring directly to the held staging inode."""

    if not _path_matches_file(path, expected):
        return
    try:
        path.unlink()
    except FileNotFoundError:
        pass


@dataclass(frozen=True)
class CommandResult:
    argv: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


@dataclass(frozen=True)
class FileCommandResult:
    argv: tuple[str, ...]
    returncode: int
    output_path: Path
    bytes: int
    stderr: str


class Client:
    """Invoke one configured ``aie`` binary without a shell.

    ``binary`` may be a path or an argument prefix such as
    ``(sys.executable, "/path/to/fake-aie.py")``. Every user value remains one
    subprocess argument; command strings are intentionally not accepted.
    """

    def __init__(
        self,
        binary: BinaryCommand = "aie",
        *,
        cwd: Optional[PathToken] = None,
        env: Optional[Mapping[str, str]] = None,
        timeout: Optional[float] = None,
    ) -> None:
        self._binary = _binary_prefix(binary)
        self._cwd = None if cwd is None else Path(_token(cwd, "cwd"))
        if timeout is not None and timeout <= 0:
            raise ValueError("timeout must be positive")
        self._timeout = timeout
        if env is None:
            self._env = None
        else:
            merged = dict(os.environ)
            for key, value in env.items():
                if not isinstance(key, str) or not key or "=" in key or "\0" in key:
                    raise ValueError(f"invalid environment variable name {key!r}")
                if not isinstance(value, str) or "\0" in value:
                    raise ValueError(f"invalid value for environment variable {key!r}")
                merged[key] = value
            self._env = merged

    @property
    def binary(self) -> tuple[str, ...]:
        return self._binary

    def run(
        self,
        args: Sequence[PathToken],
        *,
        check: bool = True,
        input_text: Optional[str] = None,
        timeout: Optional[float] = None,
    ) -> CommandResult:
        """Run an argv sequence with captured UTF-8 output and no shell."""

        argv = self._binary + _args(args)
        selected_timeout = self._timeout if timeout is None else timeout
        if selected_timeout is not None and selected_timeout <= 0:
            raise ValueError("timeout must be positive")
        try:
            completed = subprocess.run(
                argv,
                input=input_text,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self._cwd,
                env=self._env,
                timeout=selected_timeout,
                check=False,
                shell=False,
                text=True,
                encoding="utf-8",
            )
        except FileNotFoundError as error:
            raise ExecutableNotFoundError(
                f"cannot find Gravlax executable {self._binary[0]!r}; "
                "install aie or pass Client(binary=...)"
            ) from error
        except subprocess.TimeoutExpired as error:
            raise CommandTimeoutError(argv, selected_timeout) from error
        except OSError as error:
            raise ExecutableNotFoundError(
                f"cannot start Gravlax executable {self._binary[0]!r}: {error}"
            ) from error
        result = CommandResult(
            argv=argv,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
        )
        if check and result.returncode != 0:
            raise CommandError(result)
        return result

    def run_to_file(
        self,
        args: Sequence[PathToken],
        output: PathToken,
        *,
        replace: bool = False,
        input_text: Optional[str] = None,
        timeout: Optional[float] = None,
    ) -> FileCommandResult:
        """Stream stdout to a completed file without retaining it in memory.

        Output first goes to a temporary file in the destination directory. On
        success it is installed atomically; an existing destination is refused
        unless ``replace=True`` was explicit. Failed commands leave no partial
        destination.
        """

        argv = self._binary + _args(args)
        destination = Path(_token(output, "output"))
        if not destination.is_absolute():
            base = self._cwd if self._cwd is not None else Path.cwd()
            destination = base / destination
        destination = Path(os.path.abspath(destination))
        if not destination.name:
            raise ValueError("output must name a file")
        if not destination.parent.is_dir():
            raise FileNotFoundError(
                f"output directory does not exist: {destination.parent}"
            )
        if destination.exists() and not replace:
            raise FileExistsError(f"refusing to overwrite output {destination}")

        selected_timeout = self._timeout if timeout is None else timeout
        if selected_timeout is not None and selected_timeout <= 0:
            raise ValueError("timeout must be positive")
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{destination.name}.",
            suffix=".gravlax-tmp",
            dir=destination.parent,
        )
        temporary = Path(temporary_name)
        staging_identity = os.fstat(descriptor)
        publication_link: Optional[Path] = None
        try:
            # Keep the original descriptor alive through publication. A
            # directory writer can replace `temporary`, but cannot change
            # which inode this descriptor and the child process refer to.
            with os.fdopen(descriptor, "wb", closefd=False) as output_handle:
                try:
                    completed = subprocess.run(
                        argv,
                        input=input_text,
                        stdout=output_handle,
                        stderr=subprocess.PIPE,
                        cwd=self._cwd,
                        env=self._env,
                        timeout=selected_timeout,
                        check=False,
                        shell=False,
                        text=True,
                        encoding="utf-8",
                    )
                except FileNotFoundError as error:
                    raise ExecutableNotFoundError(
                        f"cannot find Gravlax executable {self._binary[0]!r}; "
                        "install aie or pass Client(binary=...)"
                    ) from error
                except subprocess.TimeoutExpired as error:
                    raise CommandTimeoutError(argv, selected_timeout) from error
                except OSError as error:
                    raise ExecutableNotFoundError(
                        f"cannot start Gravlax executable {self._binary[0]!r}: {error}"
                    ) from error
                output_handle.flush()
            if completed.returncode != 0:
                raise CommandError(
                    CommandResult(
                        argv=argv,
                        returncode=completed.returncode,
                        stdout="",
                        stderr=completed.stderr,
                    )
                )
            if replace:
                # rename/replace accepts only a name, so first give the held
                # inode a private hard-link name. A staging-name substitution
                # cannot affect it; a publication-name substitution is caught
                # by the destination identity check below.
                publication_link = _link_descriptor_to_unique_path(
                    descriptor, temporary, destination
                )
                os.replace(publication_link, destination)
                publication_link = None
            else:
                try:
                    # Linking the open descriptor provides an atomic
                    # create-new operation without trusting the staging name.
                    _link_held_descriptor(descriptor, temporary, destination)
                except FileExistsError as error:
                    raise FileExistsError(
                        f"refusing to overwrite output {destination}"
                    ) from error
            published = _verify_path_matches_descriptor(
                destination, descriptor, "published output"
            )
            return FileCommandResult(
                argv=argv,
                returncode=completed.returncode,
                output_path=destination,
                bytes=published.st_size,
                stderr=completed.stderr,
            )
        finally:
            # A hostile or concurrent writer may have swapped either private
            # name. Never remove a replacement merely because it reused ours.
            try:
                if publication_link is not None:
                    _unlink_if_same_file(publication_link, staging_identity)
                _unlink_if_same_file(temporary, staging_identity)
            finally:
                os.close(descriptor)

    def run_json(self, args: Sequence[PathToken]) -> Any:
        """Run a successful command and decode its strict JSON response."""

        result = self.run(args)
        return loads_document(result.stdout, "aie JSON response")

    def result_raw(self, args: Sequence[PathToken]) -> Any:
        """Decode a legacy command-specific JSON result without reinterpreting it."""

        return self.run_json(args)

    def result_envelope(self, args: Sequence[PathToken]) -> ResultEnvelope:
        """Parse JSON only when the command documents the shared result envelope."""

        result = self.run(args)
        return ResultEnvelope.from_json(result.stdout)

    def result_bundle(self, args: Sequence[PathToken]) -> UniformResultBundle:
        """Run an explicit JSON result API and validate its named uniform tables.

        This method intentionally captures the response because it returns materialized
        Python rows. Use :meth:`result_bundle_to_file` for output whose size is not known
        to be modest.
        """

        result = self.run(args)
        return UniformResultBundle.from_json(result.stdout)

    def result_bundle_from_file(self, path: PathToken) -> UniformResultBundle:
        """Parse a completed uniform JSON bundle from a local file."""

        selected = Path(_token(path, "result bundle path"))
        if not selected.is_absolute() and self._cwd is not None:
            selected = self._cwd / selected
        return UniformResultBundle.from_file(selected)

    def result_bundle_to_file(
        self,
        args: Sequence[PathToken],
        output: PathToken,
        *,
        replace: bool = False,
        timeout: Optional[float] = None,
    ) -> FileCommandResult:
        """Stream a command's uniform JSON stdout to an atomic completed file.

        The response is not loaded into Python memory. ``args`` must explicitly request
        the command's JSON uniform format; the resulting file can later be validated with
        :meth:`result_bundle_from_file` when materialization is desired.
        """

        return self.run_to_file(args, output, replace=replace, timeout=timeout)

    def query_region(
        self,
        archive: PathToken,
        locus: str,
        *,
        top: int = 20,
        cells: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        aggregation: str = "auto",
    ) -> UniformResultBundle:
        """Return typed uniform region counts; ``top=0`` requests every row."""

        return self.result_bundle(
            _uniform_count_query_args(
                archive,
                "region",
                locus,
                top=top,
                cells=cells,
                groups=groups,
                aggregation=aggregation,
            )
        )

    def query_region_to_file(
        self,
        archive: PathToken,
        locus: str,
        output: PathToken,
        *,
        top: int = 20,
        cells: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        aggregation: str = "auto",
        replace: bool = False,
        timeout: Optional[float] = None,
    ) -> FileCommandResult:
        """Stream typed uniform region counts to a completed JSON file."""

        args = _uniform_count_query_args(
            archive,
            "region",
            locus,
            top=top,
            cells=cells,
            groups=groups,
            aggregation=aggregation,
        )
        return self.result_bundle_to_file(
            args, output, replace=replace, timeout=timeout
        )

    def query_junction(
        self,
        archive: PathToken,
        locus: str,
        *,
        top: int = 20,
        cells: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        aggregation: str = "auto",
    ) -> UniformResultBundle:
        """Return typed uniform exact-junction counts; ``top=0`` requests every row."""

        return self.result_bundle(
            _uniform_count_query_args(
                archive,
                "junction",
                locus,
                top=top,
                cells=cells,
                groups=groups,
                aggregation=aggregation,
            )
        )

    def query_junction_to_file(
        self,
        archive: PathToken,
        locus: str,
        output: PathToken,
        *,
        top: int = 20,
        cells: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        aggregation: str = "auto",
        replace: bool = False,
        timeout: Optional[float] = None,
    ) -> FileCommandResult:
        """Stream typed uniform exact-junction counts to a completed JSON file."""

        args = _uniform_count_query_args(
            archive,
            "junction",
            locus,
            top=top,
            cells=cells,
            groups=groups,
            aggregation=aggregation,
        )
        return self.result_bundle_to_file(
            args, output, replace=replace, timeout=timeout
        )

    def query_cooccurrence(
        self,
        archive: PathToken,
        predicates: Mapping[str, str],
        expression: str,
        *,
        universe: str,
        unit: str = "molecule-record",
        region_match: str = "anchor",
        placements: str = "unique",
        allow_full_scan: bool = False,
        cells: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        aggregation: str = "auto",
        emit_membership: bool = False,
        max_memberships: int = 1_000_000,
        max_pattern_rows: int = 1_000_000,
        max_chunks: int = 100_000,
        max_evidence_records: int = 10_000_000,
        max_terminal_events: int = 10_000_000,
    ) -> UniformResultBundle:
        """Evaluate Boolean predicates on the same retained archive evidence unit.

        ``universe`` names the positive predicate that bounds candidate records. A
        negative expression term means only "not observed in retained evidence";
        selection is unknown when representative reduction prevents that absence
        from being established.
        ``unit='umi-class'`` explicitly merges records with the same
        barcode-corrected cell and exact raw UMI value across the archive. It
        does not collapse one-mismatch UMI edges and may combine distinct
        physical molecules after a cell-local UMI collision.
        ``placements='unique'`` excludes multimappers. The explicit ``direct``
        and ``all`` modes are diagnostics and do not prove physical-molecule
        co-occurrence.
        """

        return self.result_bundle(
            _cooccurrence_query_args(
                archive,
                predicates,
                expression,
                universe,
                unit=unit,
                region_match=region_match,
                placements=placements,
                allow_full_scan=allow_full_scan,
                cells=cells,
                groups=groups,
                aggregation=aggregation,
                emit_membership=emit_membership,
                max_memberships=max_memberships,
                max_pattern_rows=max_pattern_rows,
                max_chunks=max_chunks,
                max_evidence_records=max_evidence_records,
                max_terminal_events=max_terminal_events,
            )
        )

    def query_cooccurrence_to_file(
        self,
        archive: PathToken,
        predicates: Mapping[str, str],
        expression: str,
        output: PathToken,
        *,
        universe: str,
        unit: str = "molecule-record",
        region_match: str = "anchor",
        placements: str = "unique",
        allow_full_scan: bool = False,
        cells: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        aggregation: str = "auto",
        emit_membership: bool = False,
        max_memberships: int = 1_000_000,
        max_pattern_rows: int = 1_000_000,
        max_chunks: int = 100_000,
        max_evidence_records: int = 10_000_000,
        max_terminal_events: int = 10_000_000,
        replace: bool = False,
        timeout: Optional[float] = None,
    ) -> FileCommandResult:
        """Stream a Boolean co-occurrence result to a completed JSON file."""

        args = _cooccurrence_query_args(
            archive,
            predicates,
            expression,
            universe,
            unit=unit,
            region_match=region_match,
            placements=placements,
            allow_full_scan=allow_full_scan,
            cells=cells,
            groups=groups,
            aggregation=aggregation,
            emit_membership=emit_membership,
            max_memberships=max_memberships,
            max_pattern_rows=max_pattern_rows,
            max_chunks=max_chunks,
            max_evidence_records=max_evidence_records,
            max_terminal_events=max_terminal_events,
        )
        return self.result_bundle_to_file(
            args, output, replace=replace, timeout=timeout
        )

    def collection_find_events(
        self,
        collection: PathToken,
        *,
        kinds: Sequence[str] = (),
        design: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        require_groups: Sequence[str] = (),
        min_group_umi_classes: int = 1,
        min_donors: int = 1,
        min_samples: int = 1,
        min_umi_classes: int = 1,
        min_side_umi_classes: int = 1,
        min_support: int = 2,
        terminal_cluster_bp: int = 25,
        max_terminal_events: int = 10_000_000,
        annotation: Optional[PathToken] = None,
        assembly: Optional[str] = None,
        annotation_label: Optional[str] = None,
        annotation_digest: Optional[str] = None,
        novel_only: bool = False,
        solo_strand: str = "forward",
        max_candidates: int = 100_000,
        max_candidates_considered: int = 1_000_000,
        max_routed_entries: int = 10_000_000,
        max_exact_match_attempts: int = 25_000_000,
        max_annotation_comparisons: int = 10_000_000,
        verify_content: bool = False,
    ) -> UniformResultBundle:
        """Reverse-search unique-chain events across a collection.

        Candidate coordinates come from collection metadata; final sample, donor,
        and group thresholds are applied to exact raw-UMI-value classes in the
        rooted source archives. One-mismatch UMI edges and multimappers are not
        included. Annotation assembly and label values are caller declarations.
        """

        return self.result_bundle(
            _collection_find_events_args(
                collection,
                kinds=kinds,
                design=design,
                groups=groups,
                require_groups=require_groups,
                min_group_umi_classes=min_group_umi_classes,
                min_donors=min_donors,
                min_samples=min_samples,
                min_umi_classes=min_umi_classes,
                min_side_umi_classes=min_side_umi_classes,
                min_support=min_support,
                terminal_cluster_bp=terminal_cluster_bp,
                max_terminal_events=max_terminal_events,
                annotation=annotation,
                assembly=assembly,
                annotation_label=annotation_label,
                annotation_digest=annotation_digest,
                novel_only=novel_only,
                solo_strand=solo_strand,
                max_candidates=max_candidates,
                max_candidates_considered=max_candidates_considered,
                max_routed_entries=max_routed_entries,
                max_exact_match_attempts=max_exact_match_attempts,
                max_annotation_comparisons=max_annotation_comparisons,
                verify_content=verify_content,
            )
        )

    def collection_find_events_to_file(
        self,
        collection: PathToken,
        output: PathToken,
        *,
        kinds: Sequence[str] = (),
        design: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        require_groups: Sequence[str] = (),
        min_group_umi_classes: int = 1,
        min_donors: int = 1,
        min_samples: int = 1,
        min_umi_classes: int = 1,
        min_side_umi_classes: int = 1,
        min_support: int = 2,
        terminal_cluster_bp: int = 25,
        max_terminal_events: int = 10_000_000,
        annotation: Optional[PathToken] = None,
        assembly: Optional[str] = None,
        annotation_label: Optional[str] = None,
        annotation_digest: Optional[str] = None,
        novel_only: bool = False,
        solo_strand: str = "forward",
        max_candidates: int = 100_000,
        max_candidates_considered: int = 1_000_000,
        max_routed_entries: int = 10_000_000,
        max_exact_match_attempts: int = 25_000_000,
        max_annotation_comparisons: int = 10_000_000,
        verify_content: bool = False,
        replace: bool = False,
        timeout: Optional[float] = None,
    ) -> FileCommandResult:
        """Stream collection event discovery to a completed uniform JSON file."""

        args = _collection_find_events_args(
            collection,
            kinds=kinds,
            design=design,
            groups=groups,
            require_groups=require_groups,
            min_group_umi_classes=min_group_umi_classes,
            min_donors=min_donors,
            min_samples=min_samples,
            min_umi_classes=min_umi_classes,
            min_side_umi_classes=min_side_umi_classes,
            min_support=min_support,
            terminal_cluster_bp=terminal_cluster_bp,
            max_terminal_events=max_terminal_events,
            annotation=annotation,
            assembly=assembly,
            annotation_label=annotation_label,
            annotation_digest=annotation_digest,
            novel_only=novel_only,
            solo_strand=solo_strand,
            max_candidates=max_candidates,
            max_candidates_considered=max_candidates_considered,
            max_routed_entries=max_routed_entries,
            max_exact_match_attempts=max_exact_match_attempts,
            max_annotation_comparisons=max_annotation_comparisons,
            verify_content=verify_content,
        )
        return self.result_bundle_to_file(
            args, output, replace=replace, timeout=timeout
        )

    def resolve(
        self,
        annotation_file: PathToken,
        identifiers: Sequence[PathToken],
        *,
        assembly: str,
        annotation: str,
        annotation_digest: Optional[str] = None,
    ) -> ResultEnvelope:
        """Resolve biological identifiers with explicit reference identity.

        The command uses Gravlax's shared typed result envelope. An ambiguous or
        missing identifier fails the entire call instead of returning a partial
        batch.
        """

        selected_identifiers = _args(identifiers)
        if not selected_identifiers:
            raise ValueError("identifiers must contain at least one value")
        args: list[PathToken] = [
            "resolve",
            _token(annotation_file, "annotation_file"),
            _option("assembly", assembly),
            _option("annotation", annotation),
            "--format=json",
        ]
        if annotation_digest is not None:
            args.append(_option("annotation-digest", annotation_digest))
        args.append("--")
        args.extend(selected_identifiers)
        return self.result_envelope(args)

    def compare_annotations(
        self,
        archive: PathToken,
        annotation_a: PathToken,
        annotation_b: PathToken,
        *,
        assembly: str,
        annotation_a_label: str,
        annotation_b_label: str,
        annotation_a_digest: Optional[str] = None,
        annotation_b_digest: Optional[str] = None,
        gene_key: str = "unversioned",
        solo_strand: str = "forward",
        max_molecule_witnesses: int = 10_000,
        max_row_transitions_per_molecule: int = 32,
        allow_identical: bool = False,
    ) -> AnnotationComparisonResult:
        """Run two independent annotation reductions and parse all explanation tables.

        Final count deltas are exact B-minus-A results. Class transitions and
        contributing causes are deliberately not exposed as additive delta attribution.
        """

        gene_key = _choice(gene_key, {"unversioned", "exact"}, "gene_key")
        solo_strand = _choice(
            solo_strand,
            {"forward", "reverse", "unstranded"},
            "solo_strand",
        )
        max_molecule_witnesses = _nonnegative(
            max_molecule_witnesses, "max_molecule_witnesses"
        )
        max_row_transitions_per_molecule = _nonnegative(
            max_row_transitions_per_molecule,
            "max_row_transitions_per_molecule",
        )
        args: list[PathToken] = [
            "compare-annotations",
            _token(archive, "archive"),
            _option("annotation-a", annotation_a),
            _option("annotation-b", annotation_b),
            _option("assembly", assembly),
            _option("annotation-a-label", annotation_a_label),
            _option("annotation-b-label", annotation_b_label),
            _option("gene-key", gene_key),
            _option("solo-strand", solo_strand),
            _option("max-molecule-witnesses", str(max_molecule_witnesses)),
            _option(
                "max-row-transitions-per-molecule",
                str(max_row_transitions_per_molecule),
            ),
            "--format=json",
        ]
        if annotation_a_digest is not None:
            args.append(_option("annotation-a-digest", annotation_a_digest))
        if annotation_b_digest is not None:
            args.append(_option("annotation-b-digest", annotation_b_digest))
        if allow_identical:
            args.append("--allow-identical")
        result = self.run(args)
        return AnnotationComparisonResult.from_json(result.stdout)

    def transcript_ecs(
        self,
        archive: PathToken,
        annotation_file: PathToken,
        *,
        assembly: str,
        annotation_label: str,
        feature: Optional[str] = None,
        locus: Optional[str] = None,
        annotation_digest: Optional[str] = None,
        solo_strand: str = "forward",
        cells: Optional[PathToken] = None,
        groups: Optional[PathToken] = None,
        aggregation: str = "auto",
        emit_membership: bool = False,
        max_ecs: int = 100_000,
        max_memberships: int = 1_000_000,
    ) -> TranscriptEquivalenceResult:
        """Return annotation-conditional transcript compatibility sets for UMI classes.

        This call performs no abundance estimation, isoform selection, or phasing.
        """

        if (feature is None) == (locus is None):
            raise ValueError("exactly one of feature or locus must be supplied")
        if cells is not None and groups is not None:
            raise ValueError("cells and groups are mutually exclusive")
        aggregation = _choice(
            aggregation, {"auto", "cell", "group", "bulk"}, "aggregation"
        )
        if aggregation == "group" and groups is None:
            raise ValueError("group aggregation requires groups")
        solo_strand = _choice(
            solo_strand,
            {"forward", "reverse", "unstranded"},
            "solo_strand",
        )
        max_ecs = _positive(max_ecs, "max_ecs")
        max_memberships = _positive(max_memberships, "max_memberships")
        args: list[PathToken] = [
            "query",
            _token(archive, "archive"),
            "transcript-ecs",
            _option("annotation-file", annotation_file),
            _option("assembly", assembly),
            _option("annotation-label", annotation_label),
            _option("solo-strand", solo_strand),
            _option("max-ecs", str(max_ecs)),
            _option("max-memberships", str(max_memberships)),
            "--format=json",
        ]
        if annotation_digest is not None:
            args.append(_option("annotation-digest", annotation_digest))
        if feature is not None:
            args.append(_option("feature", feature))
        else:
            assert locus is not None
            args.append(_option("locus", locus))
        if cells is not None:
            args.append(_option("cells", cells))
        if groups is not None:
            args.append(_option("groups", groups))
        if aggregation != "auto":
            args.append(_option("agg", aggregation))
        if emit_membership:
            args.append("--emit-membership")
        result = self.run(args)
        return TranscriptEquivalenceResult.from_json(result.stdout)

    def project_init(
        self, directory: PathToken = ".", *, name: Optional[str] = None
    ) -> CommandResult:
        args: list[PathToken] = ["project", "init"]
        if name is not None:
            args.append(_option("name", name))
        args.extend(["--", directory])
        return self.run(args)

    def project_add(
        self,
        name: str,
        path: PathToken,
        *,
        kind: Optional[str] = None,
        project: Optional[PathToken] = None,
        replace: bool = False,
        external: bool = False,
        assembly: Optional[str] = None,
        annotation_label: Optional[str] = None,
    ) -> CommandResult:
        args: list[PathToken] = ["project", "add"]
        if kind is not None:
            args.append(_option("kind", kind))
        if project is not None:
            args.append(_option("project", project))
        if replace:
            args.append("--replace")
        if external:
            args.append("--external")
        if assembly is not None:
            args.append(_option("assembly", assembly))
        if annotation_label is not None:
            args.append(_option("annotation-label", annotation_label))
        args.extend(["--", name, path])
        return self.run(args)

    def project_show(self, *, project: Optional[PathToken] = None) -> Project:
        args: list[PathToken] = ["project", "show", "--json"]
        if project is not None:
            args.append(_option("project", project))
        return Project.from_mapping(self.run_json(args))

    def plan_check(
        self,
        plan: PathToken,
        *,
        project: Optional[PathToken] = None,
        explain: bool = False,
    ) -> ResolvedPlan:
        args: list[PathToken] = ["plan", "check", "--json"]
        if project is not None:
            args.append(_option("project", project))
        if explain:
            args.append("--explain")
        args.extend(["--", plan])
        result = self.run(args)
        document = loads_document(result.stdout, "resolved-plan response")
        return ResolvedPlan.from_mapping(document, explanation_text=result.stderr)

    def plan_run(
        self,
        plan: PathToken,
        *,
        project: Optional[PathToken] = None,
        explain: bool = False,
        dry_run: bool = False,
        resume: bool = False,
    ) -> CommandResult:
        args: list[PathToken] = ["plan", "run"]
        if project is not None:
            args.append(_option("project", project))
        if explain:
            args.append("--explain")
        if dry_run:
            args.append("--dry-run")
        if resume:
            args.append("--resume")
        args.extend(["--", plan])
        return self.run(args)

    def step_completion(self, path: PathToken) -> StepCompletion:
        """Parse a plan-run completion record from the local project."""

        selected = Path(_token(path, "completion path"))
        if not selected.is_absolute() and self._cwd is not None:
            selected = self._cwd / selected
        try:
            source = selected.read_bytes()
        except OSError as error:
            raise ProtocolError(f"cannot read step completion {selected}: {error}") from error
        return StepCompletion.from_mapping(
            loads_document(source, f"step completion {selected}")
        )

    def doctor(
        self,
        paths: Sequence[PathToken] = (),
        *,
        project: Optional[PathToken] = None,
        verify_content: bool = False,
        strict: bool = False,
    ) -> DoctorReport:
        args: list[PathToken] = ["doctor", "--json"]
        if project is not None:
            args.append(_option("project", project))
        if verify_content:
            args.append("--verify-content")
        if strict:
            args.append("--strict")
        path_args = _args(paths)
        if path_args:
            args.append("--")
            args.extend(path_args)
        result = self.run(args, check=False)
        try:
            document = loads_document(result.stdout, "doctor response")
            return DoctorReport.from_mapping(
                document, exit_code=result.returncode, stderr=result.stderr
            )
        except ProtocolError:
            if result.returncode != 0:
                raise CommandError(result)
            raise
