#!/usr/bin/env python3
"""Interleaved pre/post GQ benchmark, with matched built-in commands where available."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import runpy
import statistics
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent
semantic_result = runpy.run_path(str(HERE / 'benchmark-gq-execution.py'))['semantic_result']
states = runpy.run_path(str(HERE / 'benchmark-gq.py'))['states']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--after', type=Path, required=True)
    parser.add_argument('--archive', type=Path, action='append', required=True)
    parser.add_argument('--repeats', type=int, default=11)
    parser.add_argument('--cpus', help='Comma-separated allowed CPU IDs')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1 or args.output.exists():
        parser.error('use positive repetitions and a new output path')
    if args.cpus:
        cpus = {int(c) for c in args.cpus.split(',')}
        if not cpus <= os.sched_getaffinity(0):
            parser.error('requested CPUs outside allowed affinity')
        os.sched_setaffinity(0, cpus)
    binaries = {k: str(p.resolve(strict=True)) for k, p in [('before', args.before), ('after', args.after)]}
    env = dict(os.environ, RAYON_NUM_THREADS='4')
    rng = random.Random(6026)
    prefix = ('header {gq=1,assembly="GRCh38"}\n'
              'let U=g[1:198600000..198800000]\n'
              'let R=g[1:198696711..198696909:-]\n'
              'let A=j[1:198692373..198696711:-]\n'
              'let B=j[1:198692373..198703297:-]\n')
    cases = {
        'record_marginals': ('records', 'where any stored {overlaps(U)} |> derive {q=any unique {A} & any unique {B}} |> tally {q}', ('record', 'A & B')),
        'same_observation': ('records', 'where any stored {overlaps(U)} |> derive {q=any unique {A | B}} |> tally {q}', ('any-placement', 'A | B')),
        'universal': ('records', 'where any stored {overlaps(U)} |> derive {q=nonempty(unique) & all unique {A | B}} |> tally {q}', ('all-placements', 'A | B')),
        'geometry_tally': ('records', 'derive {a=any unique {overlaps(R,min:12bp)},b=any unique {B},joint=any unique {overlaps(R,min:12bp) & B}} |> tally {a,b,joint} by {sample}', None),
        'truth_read_summary': ('records', 'derive {q=any unique {B}} |> summarize {n=count(),t=count_true(q),reads=sum(reads.total)} by {sample}', None),
        'class_summary': ('classes', 'derive {q=any record {any unique {B}}} |> summarize {n=count(),t=count_true(q),reads=sum(reads.total)} by {sample}', None),
        'distinct_counts': ('records', 'summarize {cells=count_distinct(cell.id)} by {sample}', None),
        'cell_grouped_summary': ('records', 'summarize {n=count(),reads=sum(reads.total)} by {cell.id}', None),
        'support_fallback': ('records', 'derive {s=support_reads(overlaps(R,min:12bp))} |> select {unit.id,s} |> sort {unit.id} |> take 20', None),
        'junction_enumeration': (None, 'from junctions(@x,within:U) |> support(unit:class,by:{sample})', None),
    }
    report = {'method': 'Warm cache, one warmup, randomized interleaving. GQ internal post-parse clock excludes source reading, parsing, process startup and publication. Built-in complete process wall time includes startup, argument parsing and output: this comparison favors GQ. Pre/post GQ results and work counters checked exactly, built-in state totals checked where matched. RSS is whole-process peak from GNU time.',
              'load_start': os.getloadavg(), 'affinity': sorted(os.sched_getaffinity(0)),
              'rayon_threads': 4, 'repeats': args.repeats,
              'binaries': {k: {'path': p, 'sha256': hashlib.sha256(Path(p).read_bytes()).hexdigest()} for k, p in binaries.items()},
              'experiments': []}
    with tempfile.TemporaryDirectory(prefix='gravlax-native-bench-') as tmp:
        for archive in args.archive:
            archive = str(archive.resolve(strict=True))
            for name, (unit, body, builtin) in cases.items():
                source = prefix + (f'from @x.{unit} |> within U |> {body}' if unit else body)
                path = Path(tmp) / 'query.gq'
                path.write_text(source)
                commands = {k: [p, 'gq', 'run', str(path), '--bind', f'x={archive}', '--allow-full-scan', '--profile'] for k, p in binaries.items()}
                if builtin:
                    commands['builtin'] = [binaries['after'], 'query', archive, 'cooccur',
                        '--predicate', 'u=region:1:198600000-198800000',
                        '--predicate', 'A=junction:1:198692373-198696711:-',
                        '--predicate', 'B=junction:1:198692373-198703297:-',
                        '--universe', 'u', '--where', builtin[1], '--match-within', builtin[0],
                        '--placements', 'unique', '--region-match', 'aligned-block',
                        '--allow-full-scan', '--agg', 'bulk', '--format', 'json']
                expected = expected_states = None
                samples = {k: [] for k in commands}
                for repeat in range(args.repeats + 1):
                    order = list(commands)
                    rng.shuffle(order)
                    for label in order:
                        start = time.perf_counter_ns()
                        completed = subprocess.run(['/usr/bin/time', '-f', 'peak_rss=%M'] + commands[label], env=env, capture_output=True, check=True)
                        sample = {'process_ms': (time.perf_counter_ns() - start) / 1e6}
                        document = json.loads(completed.stdout)
                        lines = completed.stderr.decode().splitlines()
                        sample['rss_kib'] = next(int(line.split('=', 1)[1]) for line in lines if line.startswith('peak_rss='))
                        if label != 'builtin':
                            data = semantic_result(document)
                            if expected is None:
                                expected = data
                            assert expected == data, (archive, name, label, 'semantic/work mismatch')
                            sample.update(next(json.loads(line.split('=', 1)[1]) for line in lines if line.startswith('gq_profile=')))
                        if builtin:
                            observed = states(document, label != 'builtin')
                            if expected_states is None:
                                expected_states = observed
                            assert observed == expected_states, (archive, name, label, observed, expected_states)
                        if repeat:
                            samples[label].append(sample)
                medians = {k: {field: statistics.median(s[field] for s in runs) for field in runs[0] if isinstance(runs[0][field], (float, int))} for k, runs in samples.items()}
                phases = {k: {p: statistics.median(s['phases_ms'][p] for s in runs) for p in runs[0]['phases_ms']} for k, runs in samples.items() if k != 'builtin'}
                report['experiments'].append({'archive': archive, 'case': name, 'source': source, 'median': medians, 'phases_ms': phases, 'samples': samples,
                    'states': expected_states, 'sources': expected['summary']['explain']['sources'],
                    'verified_result_sha256': hashlib.sha256(json.dumps(expected, sort_keys=True).encode()).hexdigest()})
                print(Path(archive).name, name, json.dumps(medians), flush=True)
    report['load_end'] = os.getloadavg()
    with args.output.open('x') as output:
        json.dump(report, output, indent=2)


if __name__ == '__main__':
    main()
