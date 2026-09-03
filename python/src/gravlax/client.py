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
