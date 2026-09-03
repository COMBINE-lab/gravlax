"""Typed views of the project, plan, and doctor JSON protocols."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional

from ._json import (
    boolean,
    integer,
    mapping,
    nonempty_string,
    required,
    sequence,
    string,
    strings,
    validate_json_value,
)
from .exceptions import ProtocolError


@dataclass(frozen=True)
class AnnotationIdentityMetadata:
    assembly: str
    annotation: str

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "AnnotationIdentityMetadata":
        document = mapping(value, where)
        return cls(
            assembly=nonempty_string(
                required(document, "assembly", where), f"{where}.assembly"
            ),
            annotation=nonempty_string(
                required(document, "annotation", where), f"{where}.annotation"
            ),
        )


@dataclass(frozen=True)
class ProjectResource:
    name: str
    kind: str
    path: Path
    external: bool
    resolved_path: Optional[Path]
    status: str
    error: Optional[str] = None
    assembly: Optional[str] = None
    annotation_identity: Optional[AnnotationIdentityMetadata] = None

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ProjectResource":
        document = mapping(value, where)
        resolved = document.get("resolved_path")
        error = document.get("error")
        assembly = document.get("assembly")
        annotation_identity = document.get("annotation_identity")
        return cls(
            name=nonempty_string(required(document, "name", where), f"{where}.name"),
            kind=nonempty_string(required(document, "kind", where), f"{where}.kind"),
            path=Path(nonempty_string(required(document, "path", where), f"{where}.path")),
            external=boolean(
                required(document, "external", where), f"{where}.external"
            ),
            resolved_path=(
                None
                if resolved is None
                else Path(nonempty_string(resolved, f"{where}.resolved_path"))
            ),
            status=nonempty_string(
                required(document, "status", where), f"{where}.status"
            ),
            error=None if error is None else string(error, f"{where}.error"),
            assembly=(
                None
                if assembly is None
                else nonempty_string(assembly, f"{where}.assembly")
            ),
            annotation_identity=(
                None
                if annotation_identity is None
                else AnnotationIdentityMetadata.from_mapping(
                    annotation_identity, f"{where}.annotation_identity"
                )
            ),
        )


@dataclass(frozen=True)
class Project:
    schema_version: int
    name: str
    root: Path
    manifest: Path
    resources: tuple[ProjectResource, ...]

    @classmethod
    def from_mapping(cls, value: Any) -> "Project":
        where = "project response"
        validate_json_value(value, where)
        document = mapping(value, where)
        version = integer(
            required(document, "schema_version", where),
            f"{where}.schema_version",
            minimum=1,
        )
        if version != 1:
            raise ProtocolError(f"unsupported project schema version {version}; expected 1")
        items = sequence(required(document, "resources", where), f"{where}.resources")
        return cls(
            schema_version=version,
            name=nonempty_string(required(document, "name", where), f"{where}.name"),
            root=Path(nonempty_string(required(document, "root", where), f"{where}.root")),
            manifest=Path(
                nonempty_string(
                    required(document, "manifest", where), f"{where}.manifest"
                )
            ),
            resources=tuple(
                ProjectResource.from_mapping(item, f"{where}.resources[{index}]")
                for index, item in enumerate(items)
            ),
        )


@dataclass(frozen=True)
class ResolvedResource:
    kind: str
    path: Path
    external: bool
    bytes: int
    identity: "ContentIdentity"
    assembly: Optional[str] = None
    annotation_identity: Optional[AnnotationIdentityMetadata] = None

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedResource":
        document = mapping(value, where)
        assembly = document.get("assembly")
        annotation_identity = document.get("annotation_identity")
        return cls(
            kind=nonempty_string(required(document, "kind", where), f"{where}.kind"),
            path=Path(nonempty_string(required(document, "path", where), f"{where}.path")),
            external=boolean(
                required(document, "external", where), f"{where}.external"
            ),
            bytes=integer(
                required(document, "bytes", where), f"{where}.bytes", minimum=0
            ),
            identity=ContentIdentity.from_mapping(
                required(document, "identity", where), f"{where}.identity"
            ),
            assembly=(
                None
                if assembly is None
                else nonempty_string(assembly, f"{where}.assembly")
            ),
            annotation_identity=(
                None
                if annotation_identity is None
                else AnnotationIdentityMetadata.from_mapping(
                    annotation_identity, f"{where}.annotation_identity"
                )
            ),
        )


@dataclass(frozen=True)
class ContentIdentity:
    scheme: str
    digest: str

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ContentIdentity":
        document = mapping(value, where)
        return cls(
            scheme=nonempty_string(
                required(document, "scheme", where), f"{where}.scheme"
            ),
            digest=nonempty_string(
                required(document, "digest", where), f"{where}.digest"
            ),
        )


@dataclass(frozen=True)
class ResolvedProducer:
    name: str
    version: str
    plan_engine: str
    executable_identity: ContentIdentity

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedProducer":
        document = mapping(value, where)
        return cls(
            name=nonempty_string(required(document, "name", where), f"{where}.name"),
            version=nonempty_string(
                required(document, "version", where), f"{where}.version"
            ),
            plan_engine=nonempty_string(
                required(document, "plan_engine", where), f"{where}.plan_engine"
            ),
            executable_identity=ContentIdentity.from_mapping(
                required(document, "executable_identity", where),
                f"{where}.executable_identity",
            ),
        )


@dataclass(frozen=True)
class ResolvedEmbeddedResource:
    owner_resource: str
    sample: str
    role: str
    declared_path: str
    kind: str
    path: Path
    external: bool
    bytes: int
    identity: ContentIdentity

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedEmbeddedResource":
        document = mapping(value, where)
        return cls(
            owner_resource=nonempty_string(
                required(document, "owner_resource", where),
                f"{where}.owner_resource",
            ),
            sample=nonempty_string(
                required(document, "sample", where), f"{where}.sample"
            ),
            role=nonempty_string(required(document, "role", where), f"{where}.role"),
            declared_path=nonempty_string(
                required(document, "declared_path", where), f"{where}.declared_path"
            ),
            kind=nonempty_string(required(document, "kind", where), f"{where}.kind"),
            path=Path(nonempty_string(required(document, "path", where), f"{where}.path")),
            external=boolean(
                required(document, "external", where), f"{where}.external"
            ),
            bytes=integer(
                required(document, "bytes", where), f"{where}.bytes", minimum=0
            ),
            identity=ContentIdentity.from_mapping(
                required(document, "identity", where), f"{where}.identity"
            ),
        )


@dataclass(frozen=True)
class ResolvedAnnotationSemantics:
    assembly: str
    annotation: str
    source_path: Path
    source_identity: ContentIdentity

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedAnnotationSemantics":
        document = mapping(value, where)
        return cls(
            assembly=nonempty_string(
                required(document, "assembly", where), f"{where}.assembly"
            ),
            annotation=nonempty_string(
                required(document, "annotation", where), f"{where}.annotation"
            ),
            source_path=Path(
                nonempty_string(
                    required(document, "source_path", where), f"{where}.source_path"
                )
            ),
            source_identity=ContentIdentity.from_mapping(
                required(document, "source_identity", where),
                f"{where}.source_identity",
            ),
        )


@dataclass(frozen=True)
class ResolvedOutput:
    name: str
    path: Path
    kind: str
    resource_kind: str
    staging_path: Path
    annotation_semantics: Optional[ResolvedAnnotationSemantics] = None
    assembly: Optional[str] = None

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedOutput":
        document = mapping(value, where)
        kind = nonempty_string(required(document, "kind", where), f"{where}.kind")
        if kind not in {"file", "directory"}:
            raise ProtocolError(f"{where}.kind must be 'file' or 'directory'")
        annotation_semantics = document.get("annotation_semantics")
        assembly = document.get("assembly")
        return cls(
            name=nonempty_string(required(document, "name", where), f"{where}.name"),
            path=Path(nonempty_string(required(document, "path", where), f"{where}.path")),
            kind=kind,
            resource_kind=nonempty_string(
                required(document, "resource_kind", where), f"{where}.resource_kind"
            ),
            staging_path=Path(
                nonempty_string(
                    required(document, "staging_path", where), f"{where}.staging_path"
                )
            ),
            annotation_semantics=(
                None
                if annotation_semantics is None
                else ResolvedAnnotationSemantics.from_mapping(
                    annotation_semantics, f"{where}.annotation_semantics"
                )
            ),
            assembly=(
                None
                if assembly is None
                else nonempty_string(assembly, f"{where}.assembly")
            ),
        )


@dataclass(frozen=True)
class ResolvedStepInput:
    reference: str
    producer_step: str
    output_name: str
    path: Path
    kind: str
    resource_kind: str

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedStepInput":
        document = mapping(value, where)
        kind = nonempty_string(required(document, "kind", where), f"{where}.kind")
        if kind not in {"file", "directory"}:
            raise ProtocolError(f"{where}.kind must be 'file' or 'directory'")
        return cls(
            reference=nonempty_string(
                required(document, "reference", where), f"{where}.reference"
            ),
            producer_step=nonempty_string(
                required(document, "producer_step", where), f"{where}.producer_step"
            ),
            output_name=nonempty_string(
                required(document, "output_name", where), f"{where}.output_name"
            ),
            path=Path(nonempty_string(required(document, "path", where), f"{where}.path")),
            kind=kind,
            resource_kind=nonempty_string(
                required(document, "resource_kind", where), f"{where}.resource_kind"
            ),
        )


@dataclass(frozen=True)
class ResolvedPreparedInput:
    path: Path
    bytes: int
    identity: ContentIdentity
    content: str

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedPreparedInput":
        document = mapping(value, where)
        return cls(
            path=Path(nonempty_string(required(document, "path", where), f"{where}.path")),
            bytes=integer(
                required(document, "bytes", where), f"{where}.bytes", minimum=0
            ),
            identity=ContentIdentity.from_mapping(
                required(document, "identity", where), f"{where}.identity"
            ),
            content=string(required(document, "content", where), f"{where}.content"),
        )


@dataclass(frozen=True)
class ResolvedAssemblyCompatibility:
    resource: str
    kind: str
    status: str
    declared_assembly: Optional[str]
    chromosome_digest: Optional[str]
    genome_digest: Optional[str]
    note: str

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedAssemblyCompatibility":
        document = mapping(value, where)
        status = nonempty_string(
            required(document, "status", where), f"{where}.status"
        )
        if status not in {"verified", "unverified"}:
            raise ProtocolError(f"{where}.status has unknown value {status!r}")

        def optional_text(name: str) -> Optional[str]:
            value = document.get(name)
            return None if value is None else nonempty_string(value, f"{where}.{name}")

        return cls(
            resource=nonempty_string(
                required(document, "resource", where), f"{where}.resource"
            ),
            kind=nonempty_string(required(document, "kind", where), f"{where}.kind"),
            status=status,
            declared_assembly=optional_text("declared_assembly"),
            chromosome_digest=optional_text("chromosome_digest"),
            genome_digest=optional_text("genome_digest"),
            note=nonempty_string(required(document, "note", where), f"{where}.note"),
        )


@dataclass(frozen=True)
class ResolvedAnnotationInput:
    role: str
    resource: str
    annotation_path: Path
    source_path: Path
    assembly: str
    annotation: str
    source_identity: ContentIdentity
    expected_command_identity: Optional[ContentIdentity]
    compatibility: tuple[ResolvedAssemblyCompatibility, ...]

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedAnnotationInput":
        document = mapping(value, where)
        role = nonempty_string(required(document, "role", where), f"{where}.role")
        if role not in {"a", "b", "query"}:
            raise ProtocolError(f"{where}.role has unknown value {role!r}")
        expected = document.get("expected_command_identity")
        compatibility_values = sequence(
            required(document, "compatibility", where), f"{where}.compatibility"
        )
        return cls(
            role=role,
            resource=nonempty_string(
                required(document, "resource", where), f"{where}.resource"
            ),
            annotation_path=Path(
                nonempty_string(
                    required(document, "annotation_path", where),
                    f"{where}.annotation_path",
                )
            ),
            source_path=Path(
                nonempty_string(
                    required(document, "source_path", where), f"{where}.source_path"
                )
            ),
            assembly=nonempty_string(
                required(document, "assembly", where), f"{where}.assembly"
            ),
            annotation=nonempty_string(
                required(document, "annotation", where), f"{where}.annotation"
            ),
            source_identity=ContentIdentity.from_mapping(
                required(document, "source_identity", where),
                f"{where}.source_identity",
            ),
            expected_command_identity=(
                None
                if expected is None
                else ContentIdentity.from_mapping(
                    expected, f"{where}.expected_command_identity"
                )
            ),
            compatibility=tuple(
                ResolvedAssemblyCompatibility.from_mapping(
                    item, f"{where}.compatibility[{index}]"
                )
                for index, item in enumerate(compatibility_values)
            ),
        )


@dataclass(frozen=True)
class ResolvedAnnotationComparisonIntent:
    annotation_a_resource: str
    annotation_b_resource: str
    assembly: str
    gene_key: str
    solo_strand: str
    final_count_delta_semantics: str
    transition_evidence_semantics: str

    @classmethod
    def from_mapping(
        cls, value: Any, where: str
    ) -> "ResolvedAnnotationComparisonIntent":
        document = mapping(value, where)
        gene_key = nonempty_string(
            required(document, "gene_key", where), f"{where}.gene_key"
        )
        if gene_key not in {"unversioned", "exact"}:
            raise ProtocolError(f"{where}.gene_key has unknown value {gene_key!r}")
        solo_strand = nonempty_string(
            required(document, "solo_strand", where), f"{where}.solo_strand"
        )
        if solo_strand not in {"forward", "reverse", "unstranded"}:
            raise ProtocolError(
                f"{where}.solo_strand has unknown value {solo_strand!r}"
            )
        return cls(
            annotation_a_resource=nonempty_string(
                required(document, "annotation_a_resource", where),
                f"{where}.annotation_a_resource",
            ),
            annotation_b_resource=nonempty_string(
                required(document, "annotation_b_resource", where),
                f"{where}.annotation_b_resource",
            ),
            assembly=nonempty_string(
                required(document, "assembly", where), f"{where}.assembly"
            ),
            gene_key=gene_key,
            solo_strand=solo_strand,
            final_count_delta_semantics=nonempty_string(
                required(document, "final_count_delta_semantics", where),
                f"{where}.final_count_delta_semantics",
            ),
            transition_evidence_semantics=nonempty_string(
                required(document, "transition_evidence_semantics", where),
                f"{where}.transition_evidence_semantics",
            ),
        )


@dataclass(frozen=True)
class ResolvedBiologicalIntent:
    requested: str
    resolved_kind: str
    stable_id: str
    display_name: Optional[str]
    matched_by: str
    gene_ids: tuple[str, ...]
    transcript_ids: tuple[str, ...]
    annotation_resource: str
    annotation_path: Path
    assembly: str
    annotation: str
    annotation_digest: str
    contig: str
    start: int
    end: int
    strand: str
    locus: str
    compatibility: tuple[ResolvedAssemblyCompatibility, ...]

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedBiologicalIntent":
        document = mapping(value, where)
        resolved_kind = nonempty_string(
            required(document, "resolved_kind", where), f"{where}.resolved_kind"
        )
        if resolved_kind not in {"gene", "transcript", "exon"}:
            raise ProtocolError(
                f"{where}.resolved_kind has unknown value {resolved_kind!r}"
            )
        matched_by = nonempty_string(
            required(document, "matched_by", where), f"{where}.matched_by"
        )
        if matched_by not in {"stable_id", "stable_id_without_version", "gene_symbol"}:
            raise ProtocolError(f"{where}.matched_by has unknown value {matched_by!r}")
        strand = nonempty_string(
            required(document, "strand", where), f"{where}.strand"
        )
        if strand not in {"forward", "reverse"}:
            raise ProtocolError(f"{where}.strand has unknown value {strand!r}")
        start = integer(required(document, "start", where), f"{where}.start", minimum=0)
        end = integer(required(document, "end", where), f"{where}.end", minimum=0)
        if start >= end:
            raise ProtocolError(f"{where} must have start < end")
        display_name = document.get("display_name")
        compatibility = sequence(
            required(document, "compatibility", where), f"{where}.compatibility"
        )
        return cls(
            requested=nonempty_string(
                required(document, "requested", where), f"{where}.requested"
            ),
            resolved_kind=resolved_kind,
            stable_id=nonempty_string(
                required(document, "stable_id", where), f"{where}.stable_id"
            ),
            display_name=(
                None
                if display_name is None
                else nonempty_string(display_name, f"{where}.display_name")
            ),
            matched_by=matched_by,
            gene_ids=strings(required(document, "gene_ids", where), f"{where}.gene_ids"),
            transcript_ids=strings(
                required(document, "transcript_ids", where),
                f"{where}.transcript_ids",
            ),
            annotation_resource=nonempty_string(
                required(document, "annotation_resource", where),
                f"{where}.annotation_resource",
            ),
            annotation_path=Path(
                nonempty_string(
                    required(document, "annotation_path", where),
                    f"{where}.annotation_path",
                )
            ),
            assembly=nonempty_string(
                required(document, "assembly", where), f"{where}.assembly"
            ),
            annotation=nonempty_string(
                required(document, "annotation", where), f"{where}.annotation"
            ),
            annotation_digest=nonempty_string(
                required(document, "annotation_digest", where),
                f"{where}.annotation_digest",
            ),
            contig=nonempty_string(
                required(document, "contig", where), f"{where}.contig"
            ),
            start=start,
            end=end,
            strand=strand,
            locus=nonempty_string(required(document, "locus", where), f"{where}.locus"),
            compatibility=tuple(
                ResolvedAssemblyCompatibility.from_mapping(
                    item, f"{where}.compatibility[{index}]"
                )
                for index, item in enumerate(compatibility)
            ),
        )


@dataclass(frozen=True)
class ResolvedIoEstimate:
    known_selected_input_bytes: int
    known_selected_input_files: int
    unknown_prior_step_outputs: int
    read_bytes_lower_bound: int
    read_bytes_upper_bound: Optional[int]
    bound: str
    note: str

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedIoEstimate":
        document = mapping(value, where)
        bound = nonempty_string(required(document, "bound", where), f"{where}.bound")
        if bound not in {"whole-selected-files", "known-inputs-only"}:
            raise ProtocolError(f"{where}.bound has unknown value {bound!r}")
        upper = document.get("read_bytes_upper_bound")
        return cls(
            known_selected_input_bytes=integer(
                required(document, "known_selected_input_bytes", where),
                f"{where}.known_selected_input_bytes",
                minimum=0,
            ),
            known_selected_input_files=integer(
                required(document, "known_selected_input_files", where),
                f"{where}.known_selected_input_files",
                minimum=0,
            ),
            unknown_prior_step_outputs=integer(
                required(document, "unknown_prior_step_outputs", where),
                f"{where}.unknown_prior_step_outputs",
                minimum=0,
            ),
            read_bytes_lower_bound=integer(
                required(document, "read_bytes_lower_bound", where),
                f"{where}.read_bytes_lower_bound",
                minimum=0,
            ),
            read_bytes_upper_bound=(
                None
                if upper is None
                else integer(upper, f"{where}.read_bytes_upper_bound", minimum=0)
            ),
            bound=bound,
            note=nonempty_string(required(document, "note", where), f"{where}.note"),
        )


@dataclass(frozen=True)
class ResolvedUniformIo:
    kind: str
    format: str
    publication: str
    output: Optional[Path]

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "ResolvedUniformIo":
        document = mapping(value, where)
        kind = nonempty_string(required(document, "kind", where), f"{where}.kind")
        if kind not in {"result", "report"}:
            raise ProtocolError(f"{where}.kind has unknown value {kind!r}")
        format_name = nonempty_string(
            required(document, "format", where), f"{where}.format"
        )
        if format_name not in {"text", "tsv", "json"}:
            raise ProtocolError(f"{where}.format has unknown value {format_name!r}")
        publication = nonempty_string(
            required(document, "publication", where), f"{where}.publication"
        )
        if publication not in {"stdout", "atomic-no-clobber-file"}:
            raise ProtocolError(
                f"{where}.publication has unknown value {publication!r}"
            )
        raw_output = document.get("output")
        output = (
            None
            if raw_output is None
            else Path(nonempty_string(raw_output, f"{where}.output"))
        )
        if (publication == "stdout") != (output is None):
            raise ProtocolError(
                f"{where}.output must be absent exactly when publication is stdout"
            )
        return cls(
            kind=kind,
            format=format_name,
            publication=publication,
            output=output,
        )


@dataclass(frozen=True)
class ResolvedStep:
    id: str
    kind: str
    args: tuple[str, ...]
    stdout: Optional[Path]
    uniform_io: Optional[ResolvedUniformIo]
    outputs: tuple[ResolvedOutput, ...]
    input_resources: tuple[str, ...]
    step_inputs: tuple[ResolvedStepInput, ...]
    embedded_resources: tuple[str, ...]
    prepared_inputs: tuple[ResolvedPreparedInput, ...]
    biological_intent: Optional[ResolvedBiologicalIntent]
    annotation_inputs: tuple[ResolvedAnnotationInput, ...]
    annotation_comparison: Optional[ResolvedAnnotationComparisonIntent]
    output_schema_ids: tuple[str, ...]
    io_estimate: Optional[ResolvedIoEstimate]
    explanation: tuple[str, ...]

    @classmethod
    def from_mapping(
        cls, value: Any, where: str, *, resolved_plan_version: int = 3
    ) -> "ResolvedStep":
        document = mapping(value, where)
        stdout = document.get("stdout")
        output_values = sequence(
            required(document, "outputs", where), f"{where}.outputs"
        )
        prepared_values = sequence(
            document.get("prepared_inputs", []), f"{where}.prepared_inputs"
        )
        step_input_values = sequence(
            document.get("step_inputs", []), f"{where}.step_inputs"
        )
        biological_intent = document.get("biological_intent")
        annotation_input_values = sequence(
            document.get("annotation_inputs", []), f"{where}.annotation_inputs"
        )
        annotation_inputs = tuple(
            ResolvedAnnotationInput.from_mapping(
                item, f"{where}.annotation_inputs[{index}]"
            )
            for index, item in enumerate(annotation_input_values)
        )
        roles = [item.role for item in annotation_inputs]
        if len(roles) != len(set(roles)):
            raise ProtocolError(f"{where}.annotation_inputs contains duplicate roles")
        annotation_comparison_value = document.get("annotation_comparison")
        annotation_comparison = (
            None
            if annotation_comparison_value is None
            else ResolvedAnnotationComparisonIntent.from_mapping(
                annotation_comparison_value, f"{where}.annotation_comparison"
            )
        )
        io_estimate = document.get("io_estimate")
        if resolved_plan_version >= 4:
            output_schema_value = required(document, "output_schema_ids", where)
            io_estimate = required(document, "io_estimate", where)
        else:
            output_schema_value = document.get("output_schema_ids", [])
        uniform_io_value = document.get("uniform_io")
        step = cls(
            id=nonempty_string(required(document, "id", where), f"{where}.id"),
            kind=nonempty_string(required(document, "kind", where), f"{where}.kind"),
            args=strings(required(document, "args", where), f"{where}.args"),
            stdout=(
                None
                if stdout is None
                else Path(nonempty_string(stdout, f"{where}.stdout"))
            ),
            uniform_io=(
                None
                if uniform_io_value is None
                else ResolvedUniformIo.from_mapping(
                    uniform_io_value, f"{where}.uniform_io"
                )
            ),
            outputs=tuple(
                ResolvedOutput.from_mapping(item, f"{where}.outputs[{index}]")
                for index, item in enumerate(output_values)
            ),
            input_resources=strings(
                document.get("input_resources", []), f"{where}.input_resources"
            ),
            step_inputs=tuple(
                ResolvedStepInput.from_mapping(item, f"{where}.step_inputs[{index}]")
                for index, item in enumerate(step_input_values)
            ),
            embedded_resources=strings(
                document.get("embedded_resources", []),
                f"{where}.embedded_resources",
            ),
            prepared_inputs=tuple(
                ResolvedPreparedInput.from_mapping(
                    item, f"{where}.prepared_inputs[{index}]"
                )
                for index, item in enumerate(prepared_values)
            ),
            biological_intent=(
                None
                if biological_intent is None
                else ResolvedBiologicalIntent.from_mapping(
                    biological_intent, f"{where}.biological_intent"
                )
            ),
            annotation_inputs=annotation_inputs,
            annotation_comparison=annotation_comparison,
            output_schema_ids=strings(
                output_schema_value, f"{where}.output_schema_ids"
            ),
            io_estimate=(
                None
                if io_estimate is None
                else ResolvedIoEstimate.from_mapping(
                    io_estimate, f"{where}.io_estimate"
                )
            ),
            explanation=strings(
                required(document, "explanation", where), f"{where}.explanation"
            ),
        )
        if len(step.output_schema_ids) != len(set(step.output_schema_ids)):
            raise ProtocolError(f"{where}.output_schema_ids contains duplicates")
        if resolved_plan_version >= 6 and step.uniform_io is not None:
            if step.stdout is not None:
                raise ProtocolError(
                    f"{where}.stdout must be absent when uniform_io is present"
                )
            if step.uniform_io.output is not None and not any(
                output.path == step.uniform_io.output for output in step.outputs
            ):
                raise ProtocolError(
                    f"{where}.uniform_io output is missing from the declared outputs"
                )
        if resolved_plan_version >= 5 and step.kind == "compare-annotations":
            if step.annotation_comparison is None or set(roles) != {"a", "b"}:
                raise ProtocolError(
                    f"{where} compare-annotations requires annotation roles 'a' and 'b' and annotation_comparison"
                )
            by_role = {item.role: item for item in step.annotation_inputs}
            comparison = step.annotation_comparison
            if (
                by_role["a"].resource != comparison.annotation_a_resource
                or by_role["b"].resource != comparison.annotation_b_resource
                or by_role["a"].assembly != comparison.assembly
                or by_role["b"].assembly != comparison.assembly
            ):
                raise ProtocolError(
                    f"{where}.annotation_comparison disagrees with its role-named annotation inputs"
                )
            if step.biological_intent is not None:
                raise ProtocolError(
                    f"{where} compare-annotations must not declare a locus biological intent"
                )
            compare_schemas = {
                "gravlax.annotation.compare.v1",
                "gravlax.annotation.compare.count-deltas.v1",
                "gravlax.annotation.compare.class-transitions.v1",
                "gravlax.annotation.compare.contributing-causes.v1",
                "gravlax.annotation.compare.witnesses.v1",
            }
            if not step.output_schema_ids or not set(step.output_schema_ids).issubset(
                compare_schemas
            ):
                raise ProtocolError(
                    f"{where}.output_schema_ids are not an annotation-comparison contract"
                )
        if resolved_plan_version >= 5 and step.kind == "query-transcript-ecs":
            if roles != ["query"]:
                raise ProtocolError(
                    f"{where} query-transcript-ecs requires exactly one 'query' annotation input"
                )
            if step.annotation_comparison is not None:
                raise ProtocolError(
                    f"{where} query-transcript-ecs must not declare annotation_comparison"
                )
            annotation_input = step.annotation_inputs[0]
            if step.biological_intent is not None:
                intent = step.biological_intent
                if intent.resolved_kind != "gene":
                    raise ProtocolError(
                        f"{where} query-transcript-ecs biological intent must resolve to a gene"
                    )
                if (
                    intent.annotation_resource != annotation_input.resource
                    or intent.annotation_path != annotation_input.annotation_path
                    or intent.assembly != annotation_input.assembly
                    or intent.annotation != annotation_input.annotation
                    or intent.annotation_digest
                    != f"blake3:{annotation_input.source_identity.digest}"
                ):
                    raise ProtocolError(
                        f"{where}.biological_intent disagrees with its query annotation input"
                    )
            transcript_schemas = {
                "gravlax.query.transcript-ecs.v1",
                "gravlax.query.transcript-ecs.catalog.v1",
                "gravlax.query.transcript-ecs.counts.v1",
                "gravlax.query.transcript-ecs.membership.v1",
            }
            if not step.output_schema_ids or not set(step.output_schema_ids).issubset(
                transcript_schemas
            ):
                raise ProtocolError(
                    f"{where}.output_schema_ids are not a transcript-equivalence contract"
                )
        if resolved_plan_version >= 5 and step.kind in {
            "compare-annotations",
            "query-transcript-ecs",
        }:
            for index, annotation_input in enumerate(step.annotation_inputs):
                if annotation_input.source_identity.scheme != "full-file-blake3-v1":
                    raise ProtocolError(
                        f"{where}.annotation_inputs[{index}].source_identity must use full-file-blake3-v1"
                    )
                expected = annotation_input.expected_command_identity
                if expected is not None and expected != annotation_input.source_identity:
                    raise ProtocolError(
                        f"{where}.annotation_inputs[{index}].expected_command_identity disagrees with its source identity"
                    )
                if not annotation_input.compatibility:
                    raise ProtocolError(
                        f"{where}.annotation_inputs[{index}].compatibility must disclose coordinate compatibility"
                    )
        return step


@dataclass(frozen=True)
class ResolvedPlan:
    schema_version: int
    plan_schema_version: int
    producer: ResolvedProducer
    name: str
    source_path: Path
    source_digest: str
    project_name: str
    project_root: Path
    project_manifest: Path
    project_manifest_digest: str
    resources: Mapping[str, ResolvedResource]
    embedded_resources: Mapping[str, ResolvedEmbeddedResource]
    steps: tuple[ResolvedStep, ...]
    explanation_text: str = field(default="", compare=False)

    @classmethod
    def from_mapping(cls, value: Any, *, explanation_text: str = "") -> "ResolvedPlan":
        where = "resolved plan"
        validate_json_value(value, where)
        document = mapping(value, where)
        version = integer(
            required(document, "schema_version", where),
            f"{where}.schema_version",
            minimum=1,
        )
        plan_version = integer(
            required(document, "plan_schema_version", where),
            f"{where}.plan_schema_version",
            minimum=1,
        )
        if version not in {3, 4, 5, 6} or plan_version != 1:
            raise ProtocolError(
                "unsupported resolved-plan protocol: "
                f"schema_version={version}, plan_schema_version={plan_version}; "
                "supported resolved-plan versions are 3, 4, 5, and 6 with plan schema 1"
            )
        resources_document = mapping(
            required(document, "resources", where), f"{where}.resources"
        )
        embedded_document = mapping(
            document.get("embedded_resources", {}), f"{where}.embedded_resources"
        )
        step_values = sequence(required(document, "steps", where), f"{where}.steps")
        return cls(
            schema_version=version,
            plan_schema_version=plan_version,
            producer=ResolvedProducer.from_mapping(
                required(document, "producer", where), f"{where}.producer"
            ),
            name=nonempty_string(required(document, "name", where), f"{where}.name"),
            source_path=Path(
                nonempty_string(
                    required(document, "source_path", where), f"{where}.source_path"
                )
            ),
            source_digest=nonempty_string(
                required(document, "source_digest", where), f"{where}.source_digest"
            ),
            project_name=nonempty_string(
                required(document, "project_name", where), f"{where}.project_name"
            ),
            project_root=Path(
                nonempty_string(
                    required(document, "project_root", where), f"{where}.project_root"
                )
            ),
            project_manifest=Path(
                nonempty_string(
                    required(document, "project_manifest", where),
                    f"{where}.project_manifest",
                )
            ),
            project_manifest_digest=nonempty_string(
                required(document, "project_manifest_digest", where),
                f"{where}.project_manifest_digest",
            ),
            resources={
                nonempty_string(name, f"{where}.resources key"): ResolvedResource.from_mapping(
                    resource, f"{where}.resources[{name!r}]"
                )
                for name, resource in resources_document.items()
            },
            embedded_resources={
                nonempty_string(
                    name, f"{where}.embedded_resources key"
                ): ResolvedEmbeddedResource.from_mapping(
                    resource, f"{where}.embedded_resources[{name!r}]"
                )
                for name, resource in embedded_document.items()
            },
            steps=tuple(
                ResolvedStep.from_mapping(
                    step,
                    f"{where}.steps[{index}]",
                    resolved_plan_version=version,
                )
                for index, step in enumerate(step_values)
            ),
            explanation_text=explanation_text,
        )


@dataclass(frozen=True)
class DoctorSummary:
    passed: int
    warnings: int
    failures: int

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "DoctorSummary":
        document = mapping(value, where)
        return cls(
            passed=integer(
                required(document, "passed", where), f"{where}.passed", minimum=0
            ),
            warnings=integer(
                required(document, "warnings", where), f"{where}.warnings", minimum=0
            ),
            failures=integer(
                required(document, "failures", where), f"{where}.failures", minimum=0
            ),
        )


@dataclass(frozen=True)
class DoctorCheck:
    id: str
    status: str
    summary: str
    detail: Optional[str]
    remedy: Optional[str]
    data: Any

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "DoctorCheck":
        document = mapping(value, where)
        status = nonempty_string(
            required(document, "status", where), f"{where}.status"
        )
        if status not in {"pass", "warn", "fail"}:
            raise ProtocolError(f"{where}.status has unknown value {status!r}")
        detail = document.get("detail")
        remedy = document.get("remedy")
        return cls(
            id=nonempty_string(required(document, "id", where), f"{where}.id"),
            status=status,
            summary=nonempty_string(
                required(document, "summary", where), f"{where}.summary"
            ),
            detail=None if detail is None else string(detail, f"{where}.detail"),
            remedy=None if remedy is None else string(remedy, f"{where}.remedy"),
            data=document.get("data"),
        )


@dataclass(frozen=True)
class DoctorReport:
    schema: str
    version: str
    ok: bool
    strict: bool
    summary: DoctorSummary
    checks: tuple[DoctorCheck, ...]
    exit_code: int = field(default=0, compare=False)
    stderr: str = field(default="", compare=False)

    @classmethod
    def from_mapping(
        cls, value: Any, *, exit_code: int = 0, stderr: str = ""
    ) -> "DoctorReport":
        where = "doctor response"
        if isinstance(exit_code, bool) or not isinstance(exit_code, int):
            raise ProtocolError("doctor process exit status must be an integer")
        validate_json_value(value, where)
        document = mapping(value, where)
        schema = nonempty_string(
            required(document, "schema", where), f"{where}.schema"
        )
        if schema != "gravlax.doctor.v1":
            raise ProtocolError(
                f"unsupported doctor response schema {schema!r}; expected 'gravlax.doctor.v1'"
            )
        check_values = sequence(required(document, "checks", where), f"{where}.checks")
        report = cls(
            schema=schema,
            version=nonempty_string(
                required(document, "version", where), f"{where}.version"
            ),
            ok=boolean(required(document, "ok", where), f"{where}.ok"),
            strict=boolean(required(document, "strict", where), f"{where}.strict"),
            summary=DoctorSummary.from_mapping(
                required(document, "summary", where), f"{where}.summary"
            ),
            checks=tuple(
                DoctorCheck.from_mapping(check, f"{where}.checks[{index}]")
                for index, check in enumerate(check_values)
            ),
            exit_code=exit_code,
            stderr=stderr,
        )
        actual = {
            "pass": sum(check.status == "pass" for check in report.checks),
            "warn": sum(check.status == "warn" for check in report.checks),
            "fail": sum(check.status == "fail" for check in report.checks),
        }
        if (
            report.summary.passed != actual["pass"]
            or report.summary.warnings != actual["warn"]
            or report.summary.failures != actual["fail"]
        ):
            raise ProtocolError("doctor summary does not match its check statuses")
        expected_ok = report.summary.failures == 0 and (
            not report.strict or report.summary.warnings == 0
        )
        if report.ok != expected_ok:
            raise ProtocolError(
                "doctor ok flag is inconsistent with strict mode and its summary"
            )
        if (report.exit_code == 0) != report.ok:
            raise ProtocolError(
                "doctor process exit status is inconsistent with its ok flag"
            )
        return report


@dataclass(frozen=True)
class CompletedOutput:
    path: Path
    kind: str
    bytes: int
    identity: ContentIdentity

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "CompletedOutput":
        document = mapping(value, where)
        kind = nonempty_string(required(document, "kind", where), f"{where}.kind")
        if kind not in {"file", "directory"}:
            raise ProtocolError(f"{where}.kind must be 'file' or 'directory'")
        return cls(
            path=Path(nonempty_string(required(document, "path", where), f"{where}.path")),
            kind=kind,
            bytes=integer(
                required(document, "bytes", where), f"{where}.bytes", minimum=0
            ),
            identity=ContentIdentity.from_mapping(
                required(document, "identity", where), f"{where}.identity"
            ),
        )


@dataclass(frozen=True)
class CompletedStepInput:
    producer_step: str
    output_name: str
    path: Path
    kind: str
    bytes: int
    identity: ContentIdentity

    @classmethod
    def from_mapping(cls, value: Any, where: str) -> "CompletedStepInput":
        document = mapping(value, where)
        completed = CompletedOutput.from_mapping(value, where)
        return cls(
            producer_step=nonempty_string(
                required(document, "producer_step", where), f"{where}.producer_step"
            ),
            output_name=nonempty_string(
                required(document, "output_name", where), f"{where}.output_name"
            ),
            path=completed.path,
            kind=completed.kind,
            bytes=completed.bytes,
            identity=completed.identity,
        )


@dataclass(frozen=True)
class StepCompletion:
    schema_version: int
    resolved_plan_digest: str
    step_id: str
    step_digest: str
    inputs: tuple[CompletedStepInput, ...]
    outputs: tuple[CompletedOutput, ...]

    @classmethod
    def from_mapping(cls, value: Any) -> "StepCompletion":
        where = "step completion"
        validate_json_value(value, where)
        document = mapping(value, where)
        version = integer(
            required(document, "schema_version", where),
            f"{where}.schema_version",
            minimum=1,
        )
        if version != 2:
            raise ProtocolError(
                f"unsupported step-completion schema version {version}; expected 2"
            )
        input_values = sequence(required(document, "inputs", where), f"{where}.inputs")
        output_values = sequence(
            required(document, "outputs", where), f"{where}.outputs"
        )
        return cls(
            schema_version=version,
            resolved_plan_digest=nonempty_string(
                required(document, "resolved_plan_digest", where),
                f"{where}.resolved_plan_digest",
            ),
            step_id=nonempty_string(
                required(document, "step_id", where), f"{where}.step_id"
            ),
            step_digest=nonempty_string(
                required(document, "step_digest", where), f"{where}.step_digest"
            ),
            inputs=tuple(
                CompletedStepInput.from_mapping(item, f"{where}.inputs[{index}]")
                for index, item in enumerate(input_values)
            ),
            outputs=tuple(
                CompletedOutput.from_mapping(output, f"{where}.outputs[{index}]")
                for index, output in enumerate(output_values)
            ),
        )
