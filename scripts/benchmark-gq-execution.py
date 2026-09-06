#!/usr/bin/env python3
"""Interleaved GQ engine comparison using internal post-parse phase timings."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import tempfile


def semantic_result(document):
    data = copy.deepcopy(document['data'])
    data['summary'].pop('timings_ms', None)
    data['summary']['explain'].pop('physical_strategy', None)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--archive', type=Path, action='append', required=True)
    parser.add_argument('--repeats', type=int, default=7)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    binary = str(args.binary.resolve())
    env = dict(os.environ, RAYON_NUM_THREADS='4')
    rng = random.Random(37)
    cases = {
        'record_marginals': ('records', 'derive {q=any unique {A} & any unique {B}} |> tally {q}'),
        'same_observation': ('records', 'derive {q=any unique {A | B}} |> tally {q}'),
        'universal': ('records', 'derive {q=nonempty(unique) & all unique {A | B}} |> tally {q}'),
        'geometry_joint': ('records', 'derive {ra=any unique {overlaps(R,min:12bp)},ro=any unique {B},joint=any unique {overlaps(R,min:12bp) & B}} |> tally {ra,ro,joint} by {sample}'),
        'truth_and_read_sums': ('records', 'derive {q=any unique {A | B}} |> summarize {n=count(),t=count_true(q),f=count_false(q),u=count_unknown(q),reads=sum(reads.total)} by {sample}'),
        'class_closure': ('classes', 'derive {q=any record {any unique {A | B}}} |> summarize {n=count(),t=count_true(q),reads=sum(reads.total)} by {sample}'),
        'distinct_counts': ('records', 'summarize {cells=count_distinct(cell.id)} by {sample}'),
        'support_fallback': ('records', 'derive {s=support_reads(overlaps(R,min:12bp))} |> select {unit.id,s} |> sort {unit.id} |> take 20'),
        'junction_enumeration': (None, 'from junctions(@x,within:U) |> support(unit:class,by:{sample})'),
    }
    experiments = []
    with tempfile.TemporaryDirectory(prefix='gravlax-gq-execution-') as tmp:
        for archive in args.archive:
            archive = str(archive.resolve())
            for name, (unit, pipeline) in cases.items():
                source = ('header {gq=1,assembly="GRCh38"}\n'
                          'let U=g[1:198600000..198800000]\n'
                          'let R=g[1:198696711..198696909:-]\n'
                          'let A=j[1:198692373..198696711:-]\n'
                          'let B=j[1:198692373..198703297:-]\n'
                          + (f'from @x.{unit} |> within U |> {pipeline}\n' if unit else pipeline + '\n'))
                path = Path(tmp) / (name + '.gq')
                path.write_text(source)
                runs = {'auto': [], 'reference': []}
                expected = None
                strategies = {}
                for repeat in range(args.repeats + 1):
                    order = list(runs)
                    rng.shuffle(order)
                    for engine in order:
                        command = [binary, 'gq', 'run', str(path), '--bind', f'x={archive}',
                                   '--engine', engine, '--profile', '--allow-full-scan']
                        result = subprocess.run(command, env=env, capture_output=True)
                        if result.returncode:
                            raise RuntimeError((name, engine, result.stderr.decode()))
                        document = json.loads(result.stdout)
                        observed = semantic_result(document)
                        if expected is None:
                            expected = observed
                        assert observed == expected, (archive, name, engine, 'semantic mismatch')
                        profile = next(json.loads(line.split('=', 1)[1]) for line in result.stderr.decode().splitlines() if line.startswith('gq_profile='))
                        strategies[engine] = document['data']['summary']['explain']['physical_strategy']
                        if repeat:
                            runs[engine].append(profile)
                medians = {engine: {metric: statistics.median(run[metric] for run in values)
                                    for metric in ('preparation_ms', 'execution_ms', 'serialization_ms', 'post_parse_compute_ms', 'publication_ms')}
                           for engine, values in runs.items()}
                phases = {engine: {phase: statistics.median(run['phases_ms'][phase] for run in values)
                                   for phase in values[0]['phases_ms']} for engine, values in runs.items()}
                experiments.append({'archive': archive, 'case': name, 'source': source,
                                    'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
                                    'strategies': strategies, 'median_ms': medians, 'phase_medians_ms': phases,
                                    'speedup': medians['reference']['execution_ms'] / medians['auto']['execution_ms'],
                                    'runs': runs, 'verified_result': expected})
                print(name, Path(archive).name, json.dumps(medians), flush=True)
    report = {'method': 'One warmup, randomized engine interleaving, warm filesystem cache, four Rayon threads. Internal execution timing excludes startup, source reading, parsing, binding/type checking, serialization and publication; preparation and output costs reported separately. Complete tables and semantic/work summaries must match on every run. Fused evaluation/aggregation measured together, not attributed artificially to separate phases.',
              'binary_sha256': hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
              'repeats': args.repeats, 'experiments': experiments}
    with args.output.open('x') as output:
        json.dump(report, output, indent=2)


if __name__ == '__main__':
    main()
