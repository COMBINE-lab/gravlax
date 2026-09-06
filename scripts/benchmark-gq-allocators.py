#!/usr/bin/env python3
"""Interleaved allocator comparison: post-parse GQ execution and child peak RSS."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import tempfile


def semantics(document):
    data = copy.deepcopy(document['data'])
    data['summary'].pop('timings_ms', None)
    return data


def digest_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', action='append', required=True, help='allocator-label=/absolute/binary')
    parser.add_argument('--archive', type=Path, action='append', required=True)
    parser.add_argument('--repeats', type=int, default=9)
    parser.add_argument('--cpu', type=int, help='Pin harness and child processes to one allowed CPU')
    parser.add_argument('--engine', choices=['auto', 'reference'], default='auto')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1 or args.output.exists():
        parser.error('use positive repetitions and a new output filename')
    binaries = {}
    for binding in args.binary:
        name, separator, path = binding.partition('=')
        if not separator or not name or name in binaries:
            parser.error('--binary requires distinct labels and paths')
        binaries[name] = str(Path(path).resolve(strict=True))
    if len(binaries) < 2:
        parser.error('compare at least two binaries')
    if args.cpu is not None:
        if args.cpu not in os.sched_getaffinity(0):
            parser.error('--cpu is outside the allowed affinity mask')
        os.sched_setaffinity(0, {args.cpu})
    load_start = os.getloadavg()
    env = dict(os.environ, RAYON_NUM_THREADS='4')
    rng = random.Random(20260905)
    cases = {
        'geometry_tally': 'from @x.records |> within U |> derive {a=any unique {overlaps(R,min:12bp)},b=any unique {J},joint=any unique {overlaps(R,min:12bp) & J}} |> tally {a,b,joint} by {sample}',
        'truth_read_summary': 'from @x.records |> within U |> derive {q=any unique {J}} |> summarize {n=count(),t=count_true(q),reads=sum(reads.total)} by {sample}',
        'class_summary': 'from @x.classes |> within U |> derive {q=any record {any unique {J}}} |> summarize {n=count(),t=count_true(q),reads=sum(reads.total)} by {sample}',
        'distinct_counts': 'from @x.records |> within U |> summarize {cells=count_distinct(cell.id)} by {sample}',
        'cell_grouped_summary': 'from @x.records |> within U |> summarize {n=count(),reads=sum(reads.total)} by {cell.id}',
        'read_support_top_k': 'from @x.records |> within U |> derive {s=support_reads(overlaps(R,min:12bp))} |> select {unit.id,s} |> sort {unit.id} |> take 20',
        'junction_enumeration': 'from junctions(@x,within:U) |> support(unit:class,by:{sample})',
    }
    experiments = []
    with tempfile.TemporaryDirectory(prefix='gravlax-gq-allocator-') as scratch:
        for archive in args.archive:
            archive = str(archive.resolve(strict=True))
            for name, body in cases.items():
                source = ('header {gq=1,assembly="GRCh38"}\n'
                          'let U=g[1:198600000..198800000]\n'
                          'let R=g[1:198696711..198696909:-]\n'
                          'let J=j[1:198692373..198703297:-]\n' + body)
                path = Path(scratch) / (name + '.gq')
                path.write_text(source)
                samples = {label: [] for label in binaries}
                expected = None
                for repetition in range(args.repeats + 1):
                    order = list(binaries)
                    rng.shuffle(order)
                    for label in order:
                        command = ['/usr/bin/time', '-f', 'gq_rss_kib=%M', binaries[label],
                                   'gq', 'run', str(path), '--bind', f'x={archive}', '--profile',
                                   '--engine', args.engine, '--allow-full-scan']
                        completed = subprocess.run(command, env=env, capture_output=True, check=True)
                        data = semantics(json.loads(completed.stdout))
                        if expected is None:
                            expected = data
                        if data != expected:
                            raise RuntimeError((archive, name, label, 'semantic or work-counter mismatch'))
                        lines = completed.stderr.decode().splitlines()
                        profile = next(json.loads(line.split('=', 1)[1]) for line in lines if line.startswith('gq_profile='))
                        rss = next(int(line.split('=', 1)[1]) for line in lines if line.startswith('gq_rss_kib='))
                        if repetition:
                            samples[label].append(dict(profile, rss_kib=rss))
                medians = {label: {key: statistics.median(sample[key] for sample in values)
                                   for key in ['execution_ms', 'post_parse_compute_ms', 'serialization_ms', 'rss_kib']}
                           for label, values in samples.items()}
                experiments.append({'case': name, 'archive': archive, 'source': source,
                                    'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
                                    'median': medians, 'samples': samples,
                                    'verified_result_sha256': hashlib.sha256(json.dumps(expected, sort_keys=True).encode()).hexdigest(),
                                    'verified_summary': expected['summary']})
                print(name, Path(archive).name, json.dumps(medians), flush=True)
    report = {'method': 'One warmup, randomized allocator interleaving, warm cache, internal execution clock excluding parse/startup/binding/serialization; whole child peak RSS from GNU time; exact tables and semantic/work summaries compared every run. No allocator environment tuning or native malloc interposition.',
              'host': platform.node(), 'platform': platform.platform(), 'engine': args.engine,
              'affinity': sorted(os.sched_getaffinity(0)), 'rayon_threads': 4,
              'load_start': load_start, 'load_end': os.getloadavg(), 'repeats': args.repeats,
              'binaries': {label: {'path': path, 'sha256': digest_file(path), 'bytes': Path(path).stat().st_size} for label, path in binaries.items()},
              'experiments': experiments}
    with args.output.open('x') as output:
        json.dump(report, output, indent=2)


if __name__ == '__main__':
    main()
