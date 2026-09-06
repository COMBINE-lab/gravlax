#!/usr/bin/env python3
"""End-to-end, warm-cache, interleaved query comparison; fail on row disagreement."""
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import random
import statistics
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--archive', action='append', required=True)
    parser.add_argument('--junctions', required=True, help='CD45 junction catalogue JSON')
    parser.add_argument('--output', required=True)
    parser.add_argument('--repeats', type=int, default=11)
    args = parser.parse_args()
    junctions = sorted(json.loads(Path(args.junctions).read_text())['junctions'], key=lambda x: -x['umis'])
    assert len({(j['donor'], j['acceptor']) for j in junctions}) >= 63, 'need at least 63 distinct junctions; do not inflate panels with repeated atoms'
    env = dict(os.environ, RAYON_NUM_THREADS='4')
    rng = random.Random(271828)
    results = []
    identities = []
    def binary_digest(path):
        digest = hashlib.sha256()
        with open(path, 'rb') as source:
            for chunk in iter(lambda: source.read(1024 * 1024), b''):
                digest.update(chunk)
        return digest.hexdigest()
    provenance = {'host': platform.node(), 'platform': platform.platform(),
                  'before_binary_sha256': binary_digest(args.before),
                  'after_binary_sha256': binary_digest(args.after), 'archives': identities}
    with tempfile.TemporaryDirectory(prefix='gravlax-query-bench-') as scratch:
        for archive in args.archive:
            base = ['query', archive, 'cooccur', '--predicate', 'u=region:1:198600000-198800000',
                    '--where', 'u', '--universe', 'u', '--format', 'json', '--emit-membership']
            population = json.loads(subprocess.check_output([args.before, *base], env=env))
            identities.append({'path': archive, 'content_roots': population['provenance']['archives']})
            memberships = next(t['rows'] for t in population['data']['tables'] if t['name'] == 'memberships')
            barcodes = sorted({row[2] for row in memberships})
            scope = Path(scratch) / 'cells.txt'
            scope.write_text('\n'.join(barcodes[::10]) + '\n')
            for scoped in (False, True):
                for count in (2, 8, 32, 64):
                    predicates = ['u=region:1:198600000-198800000']
                    for i in range(count - 1):
                        j = junctions[i]
                        predicates.append(f'j{i}=junction:1:{j["donor"]}-{j["acceptor"]}')
                    command = ['query', archive, 'cooccur', '--where', ' | '.join(f'j{i}' for i in range(count - 1)),
                               '--universe', 'u', '--format', 'json', '--agg', 'bulk']
                    for predicate in predicates:
                        command += ['--predicate', predicate]
                    if scoped:
                        command += ['--cells', str(scope)]
                    variants = {'before': [args.before], 'scalar': [args.after],
                                'compiled': [args.after], 'auto': [args.after]}
                    runs = {name: [] for name in variants}
                    expected = None
                    for repetition in range(args.repeats + 1):
                        order = list(variants)
                        rng.shuffle(order)
                        for name in order:
                            invocation = [*variants[name], *command]
                            if name != 'before':
                                invocation += ['--engine', name]
                            start = time.perf_counter_ns()
                            completed = subprocess.run(['/usr/bin/time', '-f', 'GRAVLAX_BENCH_RSS=%M', *invocation],
                                env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                            elapsed = (time.perf_counter_ns() - start) / 1e6
                            data = json.loads(completed.stdout)['data']
                            logical = {'tables': {table['name']: table['rows'] for table in data['tables']},
                                       'counts': {k: data['summary'][k] for k in ('candidate_units', 'selected_units', 'indeterminate_units', 'candidate_cells', 'selected_cells')}}
                            digest = hashlib.sha256(json.dumps(logical, sort_keys=True).encode()).hexdigest()
                            if expected is None:
                                expected = digest
                            assert digest == expected, (archive, scoped, count, name, 'logical output mismatch')
                            rss = int(completed.stderr.decode().rsplit('GRAVLAX_BENCH_RSS=', 1)[1].strip())
                            if repetition:
                                runs[name].append({'ms': elapsed, 'rss_kib': rss})
                    row = {'archive': archive, 'scoped': scoped, 'scope_cells': len(barcodes[::10]) if scoped else len(barcodes),
                           'predicates': count, 'logical_sha256': expected, 'variants': {name: {
                               'median_ms': statistics.median(r['ms'] for r in samples),
                               'median_rss_kib': statistics.median(r['rss_kib'] for r in samples), 'samples': samples}
                               for name, samples in runs.items()}}
                    results.append(row)
                    print(json.dumps({k: v for k, v in row.items() if k != 'variants'}),
                          {name: round(data['median_ms'], 2) for name, data in row['variants'].items()}, flush=True)
    Path(args.output).write_text(json.dumps({'method': 'warm filesystem cache; one warmup; interleaved randomized order; 4 Rayon threads; complete process including JSON output; all logical rows/counts checked',
                                          'provenance': provenance, 'repeats': args.repeats, 'results': results}, indent=2) + '\n')


if __name__ == '__main__':
    main()
