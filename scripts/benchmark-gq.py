#!/usr/bin/env python3
"""Differential, interleaved warm-cache GQ/cooccur benchmark (complete processes)."""
import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import tempfile
import time


def table(document, name):
    selected = next(t for t in document['data']['tables'] if t['name'] == name)
    columns = [f['name'] for f in selected['schema']['fields']]
    return [dict(zip(columns, row)) for row in selected['rows']]


def states(document, gq):
    counts = collections.Counter()
    for row in table(document, 'results' if gq else 'patterns'):
        counts[row['q_state' if gq else 'selection_state']] += row['count' if gq else 'evidence_units']
    return dict(sorted(counts.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--archive', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    env = dict(os.environ, RAYON_NUM_THREADS='4')
    rng = random.Random(19)
    experiments = []
    with tempfile.TemporaryDirectory(prefix='gravlax-gq-bench-') as tmp:
        for archive in args.archive:
            archive = str(archive.resolve())
            for name, cli_mode, expression, gq in [
                ('record_marginals', 'record', 'A & B', 'any unique {A} & any unique {B}'),
                ('same_observation', 'any-placement', 'A | B', 'any unique {A | B}'),
                ('nonvacuous_universal', 'all-placements', 'A | B', 'nonempty(unique) & all unique {A | B}'),
            ]:
                source = ('header {gq=1,assembly="GRCh38"}\n'
                          'let U=g[1:198600000..198800000]\n'
                          'let A=j[1:198692373..198696711:-]\n'
                          'let B=j[1:198692373..198703297:-]\n'
                          'from @x.records |> within U\n'
                          '|> where any stored {overlaps(U)}\n'
                          f'|> derive {{q={gq}}} |> tally {{q}}\n')
                path = Path(tmp) / f'{name}.gq'
                path.write_text(source)
                commands = {
                    'gq': [binary, 'gq', 'run', str(path), '--bind', f'x={archive}', '--allow-full-scan'],
                    'cooccur': [binary, 'query', archive, 'cooccur',
                                '--predicate', 'u=region:1:198600000-198800000',
                                '--predicate', 'A=junction:1:198692373-198696711:-',
                                '--predicate', 'B=junction:1:198692373-198703297:-',
                                '--universe', 'u', '--where', expression,
                                '--match-within', cli_mode, '--placements', 'unique',
                                '--region-match', 'aligned-block', '--allow-full-scan',
                                '--agg', 'bulk', '--format', 'json'],
                }
                durations = collections.defaultdict(list)
                expected = None
                work = {}
                for repeat in range(args.repeats + 1):
                    order = list(commands)
                    rng.shuffle(order)
                    for implementation in order:
                        start = time.perf_counter_ns()
                        result = subprocess.run(commands[implementation], env=env, capture_output=True, check=True)
                        elapsed = (time.perf_counter_ns() - start) / 1e6
                        document = json.loads(result.stdout)
                        observed = states(document, implementation == 'gq')
                        if expected is None:
                            expected = observed
                        assert observed == expected, (name, implementation, observed, expected)
                        work[implementation] = document['data']['summary']
                        if repeat:
                            durations[implementation].append(elapsed)
                experiments.append({'archive': archive, 'case': name, 'states': expected,
                                    'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
                                    'sources': work['gq']['explain']['sources'],
                                    'median_ms': {k: statistics.median(v) for k, v in durations.items()},
                                    'runs_ms': dict(durations),
                                    'gq_work': {k: work['gq'][k] for k in ('decoded_chunks', 'decoded_records', 'expression_steps')}})
    report = {'method': 'one warmup; randomized interleaving; warm filesystem cache; four Rayon threads; complete processes and JSON output; exact Truth count equivalence required on every run',
              'binary_sha256': hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
              'repeats': args.repeats, 'experiments': experiments}
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
    print(json.dumps([{k: e[k] for k in ('archive', 'case', 'states', 'median_ms')} for e in experiments], indent=2))


if __name__ == '__main__':
    main()
