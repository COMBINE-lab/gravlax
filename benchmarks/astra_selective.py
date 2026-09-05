#!/usr/bin/env python3
"""Compare actual archive bytes and exact cooccur results using two pinned binaries.

Artifacts include raw timings/RSS, identities and result JSON. This locus fixture
measures selective-access mechanics; it is not an atlas throughput benchmark.
"""
import argparse
import hashlib
import json
import pathlib
import statistics
import subprocess
import time


def identity(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser()
    for name in ("baseline", "candidate", "bam", "whitelist", "gtf", "barcodes", "out"):
        p.add_argument("--" + name, required=True, type=pathlib.Path)
    args = p.parse_args()
    args.out.mkdir(exist_ok=False)
    results = {"schema": "gravlax.astra-selective-benchmark.v1", "inputs": {
        k: {"path": str(getattr(args, k).resolve()), "sha256": identity(getattr(args, k))}
        for k in ("baseline", "candidate", "bam", "whitelist", "gtf", "barcodes")}, "configurations": {}}

    def run(command, label):
        rss = args.out / (label + ".rss")
        start = time.perf_counter()
        proc = subprocess.run(["/usr/bin/time", "-f", "%M", "-o", str(rss), *map(str, command)], capture_output=True)
        elapsed = time.perf_counter() - start
        (args.out / (label + ".stderr")).write_bytes(proc.stderr)
        if proc.returncode:
            raise RuntimeError(f"{command}: {proc.stderr.decode()}")
        return proc.stdout, {"seconds": elapsed, "peak_rss_kib": int(rss.read_text().strip())}

    expected = {}
    expected_mex = {}
    for name, binary, flags in [("baseline", args.baseline, []), ("reader", args.candidate, []),
                                ("indexed", args.candidate, ["--access-index", "--chunk-records", "4096"]),
                                ("compressed", args.candidate, ["--compression-tuning"]),
                                ("fidelity", args.candidate, ["--geometry-fidelity", "--access-index", "--chunk-records", "4096"]),
                                ("all", args.candidate, ["--geometry-fidelity", "--access-index", "--chunk-records", "4096", "--compression-tuning"]),
                                ("compressed-l3", args.candidate, ["--compression-tuning"]),
                                ("compressed-l9", args.candidate, ["--compression-tuning"])]:
        archive = args.out / (name + ".aie")
        report, ingest = run([binary, "ingest-archive", args.bam, "--whitelist", args.whitelist,
                             "--out", archive, "--zstd-level", name.split("-l")[-1] if "-l" in name else "19", "--report-format", "json", *flags], name + "-ingest")
        (args.out / (name + "-ingest.json")).write_bytes(report)
        run([binary, "doctor", archive, "--verify-content", "--json"], name + "-doctor")
        times = []
        for rep in range(7):
            raw, measurement = run([binary, "query", archive, "cooccur", "--predicate", "a=region:1:198696711-198696909:-",
                "--predicate", "j=junction:1:198692373-198696711:-", "--where", "a & j", "--universe", "a",
                "--unit", "umi-class", "--region-match", "aligned-block", "--placements", "unique",
                "--emit-membership", "--format", "json", *([] if "--access-index" in flags else ["--allow-full-scan"])], f"{name}-query-{rep}")
            result = json.loads(raw)
            tables = result["data"]["tables"]
            group = "fidelity" if "--geometry-fidelity" in flags else "compact"
            if group not in expected:
                expected[group] = tables
            assert tables == expected[group], "Scientific result tables changed within geometry mode"
            times.append(measurement)
            if rep == 0:
                (args.out / (name + "-query.json")).write_bytes(raw)
        replay_times = []
        for rep in range(3):
            mex = args.out / f"{name}-mex-{rep}"
            _, measurement = run([binary, "replay-rows", archive, "--gtf", args.gtf,
                "--barcodes", args.barcodes, "--out-dir", mex, "--solo-strand", "unstranded"], f"{name}-replay-{rep}")
            hashes = {f: identity(mex / f) for f in ("matrix.mtx", "features.tsv", "barcodes.tsv")}
            if group not in expected_mex:
                expected_mex[group] = hashes
            assert hashes == expected_mex[group], "Replay matrix changed within geometry mode"
            replay_times.append(measurement)
        results["configurations"][name] = {"archive_bytes": archive.stat().st_size, "archive_sha256": identity(archive),
            "ingest": ingest, "query": times, "query_median_seconds": statistics.median(t["seconds"] for t in times),
            "query_median_rss_kib": statistics.median(t["peak_rss_kib"] for t in times),
            "replay": replay_times, "replay_median_seconds": statistics.median(t["seconds"] for t in replay_times),
            "replay_median_rss_kib": statistics.median(t["peak_rss_kib"] for t in replay_times), "mex_sha256":hashes}
    configs = results["configurations"]
    assert configs["baseline"]["archive_sha256"] == configs["reader"]["archive_sha256"]
    assert configs["compressed"]["archive_bytes"] <= configs["reader"]["archive_bytes"]
    assert configs["all"]["archive_bytes"] <= configs["fidelity"]["archive_bytes"]
    (args.out / "benchmark.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
