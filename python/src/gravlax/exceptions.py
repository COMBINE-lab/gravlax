"""Exceptions raised by the Gravlax Python client."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .client import CommandResult


class GravlaxError(Exception):
    """Base class for client and protocol errors."""


class ExecutableNotFoundError(GravlaxError):
    """The configured ``aie`` executable could not be started."""


class CommandError(GravlaxError):
    """An ``aie`` invocation returned a non-zero exit status."""

    def __init__(self, result: "CommandResult") -> None:
        self.result = result
        detail = result.stderr.strip() or result.stdout.strip() or "no diagnostic output"
        super().__init__(
            f"aie command failed with exit status {result.returncode}: {detail}"
        )


class CommandTimeoutError(GravlaxError):
    """An ``aie`` invocation exceeded the configured timeout."""

    def __init__(self, argv: tuple[str, ...], timeout: float | None) -> None:
        self.argv = argv
        self.timeout = timeout
        super().__init__(f"aie command timed out after {timeout} seconds")


class ProtocolError(GravlaxError, ValueError):
    """A JSON response does not match the advertised Gravlax protocol."""


class OptionalDependencyError(ImportError):
    """A requested converter's optional Python dependency is unavailable."""
