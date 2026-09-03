"""Strict JSON and small protocol-validation helpers."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from typing import Any

from .exceptions import ProtocolError


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number {value}")


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def loads_document(source: str | bytes, where: str = "response") -> Any:
    """Parse one strict JSON document, rejecting duplicate keys and NaN values."""

    try:
        return json.loads(
            source,
            parse_constant=_reject_constant,
            object_pairs_hook=_unique_object,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ProtocolError(f"{where} is not valid strict JSON: {error}") from error


def mapping(value: Any, where: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ProtocolError(f"{where} must be a JSON object")
    if not all(isinstance(key, str) for key in value):
        raise ProtocolError(f"{where} contains a non-string key")
    return value


def sequence(value: Any, where: str) -> Sequence[Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise ProtocolError(f"{where} must be a JSON array")
    return value


def string(value: Any, where: str) -> str:
    if not isinstance(value, str):
        raise ProtocolError(f"{where} must be a string")
    return value


def nonempty_string(value: Any, where: str) -> str:
    result = string(value, where)
    if not result:
        raise ProtocolError(f"{where} must not be empty")
    return result


def integer(value: Any, where: str, *, minimum: int | None = None) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ProtocolError(f"{where} must be an integer")
    if minimum is not None and value < minimum:
        raise ProtocolError(f"{where} must be at least {minimum}")
    return value


def boolean(value: Any, where: str) -> bool:
    if not isinstance(value, bool):
        raise ProtocolError(f"{where} must be a boolean")
    return value


def strings(value: Any, where: str) -> tuple[str, ...]:
    return tuple(string(item, f"{where}[{index}]") for index, item in enumerate(sequence(value, where)))


def required(document: Mapping[str, Any], key: str, where: str) -> Any:
    if key not in document:
        raise ProtocolError(f"{where} is missing required field {key!r}")
    return document[key]


def validate_json_value(value: Any, where: str = "value") -> None:
    """Validate an already-decoded value as finite, ordinary JSON data."""

    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ProtocolError(f"{where} contains a non-finite number")
        return
    if isinstance(value, Mapping):
        for key, child in value.items():
            if not isinstance(key, str):
                raise ProtocolError(f"{where} contains a non-string object key")
            validate_json_value(child, f"{where}.{key}")
        return
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for index, child in enumerate(value):
            validate_json_value(child, f"{where}[{index}]")
        return
    raise ProtocolError(f"{where} contains a non-JSON value of type {type(value).__name__}")
