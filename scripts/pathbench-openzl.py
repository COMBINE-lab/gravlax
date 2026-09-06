#!/usr/bin/env python3
"""Ablate typed OpenZL graphs against byte zstd on path-bench exported streams.

Writes a NEW output directory. This is a backend-only benchmark, not an archive
codec. Values are pre-entropy u64s (signed residuals already zigzag transformed).
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess


def varints(raw):
    out = bytearray()
    for (value,) in struct.iter_unpack("<Q", raw):
        while value >= 128:
            out.append((value & 127) | 128)
            value >>= 7
        out.append(value)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--helper", type=Path, required=True)
    a = p.parse_args()
    inputs = sorted(a.source.glob("typed-*.u64le"))
    if not inputs:
        p.error("no typed streams; run aie dev path-bench first")
    a.out.mkdir(exist_ok=False)
    rows = []
    for source in inputs:
        raw = source.read_bytes()
        vpath = a.out / (source.stem + ".varints")
        vpath.write_bytes(varints(raw))
        candidates = []
        for name, mode, level, path in [
            ("varint-zstd19", 3, 19, vpath),
            ("u64-zstd19", 3, 19, source),
            ("numeric-generic9", 0, 9, source),
            ("numeric-generic19", 0, 19, source),
            ("numeric-delta9", 1, 9, source),
            ("numeric-fieldlz9", 2, 9, source),
            ("numeric-range-generic9", 4, 9, source),
            ("numeric-range-fse9", 5, 9, source),
            ("numeric-range-zstd19", 6, 19, source),
            ("numeric-flatpack9", 7, 9, source),
            ("numeric-token-fse9", 8, 9, source),
            ("numeric-bitpack9", 9, 9, source),
        ]:
            output = a.out / (source.stem + "." + name + ".frame")
            run = subprocess.run([str(a.helper), str(mode), str(level), str(path), str(output)], capture_output=True, text=True)
            if run.returncode:
                raise RuntimeError(f"{source.name}/{name}: {run.stderr}")
            candidates.append(dict(name=name, **json.loads(run.stdout)))
        rows.append(dict(stream=source.name, input_sha256=hashlib.sha256(raw).hexdigest(), candidates=candidates))
    names = [c["name"] for c in rows[0]["candidates"]]
    totals = []
    for name in names + ["openzl-best-per-stream", "hybrid-best-per-stream"]:
        selected = []
        for row in rows:
            if name.endswith("best-per-stream"):
                eligible = [c for c in row["candidates"] if c["name"].startswith("numeric-") or (name.startswith("hybrid") and c["name"] == "varint-zstd19")]
                selected.append(min(eligible, key=lambda c: c["bytes"]))
            else:
                selected.append(next(c for c in row["candidates"] if c["name"] == name))
        encode = sum(c["encode_seconds"] for c in selected)
        if name.endswith("best-per-stream"):
            encode = sum(c["encode_seconds"] for r in rows for c in r["candidates"] if c["name"].startswith("numeric-") or (name.startswith("hybrid") and c["name"] == "varint-zstd19"))
        totals.append(dict(name=name, frame_bytes=sum(c["bytes"] for c in selected), encode_seconds=encode, decode_seconds=sum(c["decode_seconds"] for c in selected)))
    report = dict(schema="gravlax.path-openzl-backend.v1", openzl_release="v0.2.0", openzl_commit="3dceb64867840201fb8f57a29d179995f700c9b8", streams=len(rows), totals=totals, rows=rows, limitations="Actual persisted backend frames, every stream reopened and exactly verified. No outer archive/directory/index integration; frame totals are NOT whole archive sizes. Timings exclude conversion to/from u64 and geometry reconstruction. Simple standard graphs, no trained clustering or exhaustive graph search; best-per-stream pays all candidate encode costs. Empty streams included.")
    report["source_root"] = json.loads((a.source / "report.json").read_text())["source_root"]
    report["source_benchmark"] = str(a.source)
    (a.out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=2))


if __name__ == "__main__":
    main()
