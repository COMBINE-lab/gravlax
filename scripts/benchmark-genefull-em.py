#!/usr/bin/env python3
"""Interleaved complete-command GeneFull benchmarks with result equality gates.

--binaries contains aie-baseline, aie-features, aie-memory and aie-index.
These isolate mode selection, actual-support spilling, and the block/strand index.
Input project files are read only; --out must be new. Python standard library only.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def scientific(data):
    return {m['name']: m for m in data['modes']}, data['evaluation_counts']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--binaries', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--repetitions', type=int, default=3)
    args = parser.parse_args()
    root, binaries, out = args.project.resolve(), args.binaries.resolve(), args.out.resolve()
    if args.repetitions < 1:
        parser.error('repetitions must be positive')
    out.mkdir(parents=True, exist_ok=False)
    (out / 'tmp').mkdir()
    archive = root / 'runs/archive/d2p.aie'
    annotation = root / 'runs/genefull-20260910/gencode-v49.aic'
    called = root / 'runs/oracle/d2p/v49/Solo.out/Gene/filtered/barcodes.tsv'
    barcodes = root / 'runs/oracle/d2p/v49/Solo.out/Gene/raw/barcodes.tsv'
    groups = root / 'runs/genefull-paper-20260910/fixed-nuclei-groups.tsv'
    env = dict(os.environ, RAYON_NUM_THREADS='16', TMPDIR=str(out / 'tmp'))
    cases = [('nine-baseline', 'baseline', 'nine'), ('nine-index', 'index', 'nine'),
             ('four-features', 'features', 'four'), ('four-memory', 'memory', 'four'),
             ('four-index', 'index', 'four'), ('pooled-index', 'index', 'pooled'),
             ('replay-baseline', 'baseline', 'replay'), ('replay-index', 'index', 'replay')]
    manifest = {'threads': 16, 'hostname': os.uname().nodename, 'archive': str(archive),
                'annotation': str(annotation), 'called_nuclei': len(called.read_text().splitlines()),
                'binary_sha256': {v: digest(binaries / f'aie-{v}') for v in ['baseline', 'features', 'memory', 'index']},
                'scope': 'EM scores fixed 6460 nuclei; priors use all archive barcodes. Replay emits all raw barcodes.',
                'case_order': cases, 'repetitions': args.repetitions, 'runs': []}
    # Read-through warmup; no root cache eviction and no changes to archive data.
    manifest['archive_sha256'] = digest(archive)
    reference = None
    replay_reference = None
    for rep in range(args.repetitions):
        order = cases[rep:] + cases[:rep]
        for name, version, mode in order:
            label = f'{rep+1}-{name}'
            metrics = out / f'{label}.json'
            binary = binaries / f'aie-{version}'
            if mode == 'replay':
                argv = [binary, 'replay-rows', archive, '--gtf', annotation, '--gene-full',
                        '--barcodes', barcodes, '--out-dir', out / label]
            else:
                argv = [binary, 'dev', 'em', archive, '--gtf', annotation, '--gene-full',
                        '--mask', '0.2', '--seed', '7', '--alpha', '20', '--metrics-json', metrics]
                argv += ['--groups', groups] if mode == 'nine' else ['--eval-barcodes', called]
                if mode == 'pooled':
                    argv += ['--modes', 'pooled']
            argv = list(map(str, argv))
            (out / f'{label}-command.json').write_text(json.dumps(argv, indent=2) + '\n')
            timing = out / f'{label}.time'
            started = time.time()
            with (out / f'{label}.log').open('w') as stream:
                subprocess.run(['/usr/bin/time', '-f', '%e\t%M\t%U\t%S', '-o', str(timing), *argv],
                               stdout=stream, stderr=subprocess.STDOUT, env=env, check=True)
            elapsed, rss, user, system = map(float, timing.read_text().split())
            result = {'case': name, 'rep': rep + 1, 'started_unix': started,
                      'elapsed_seconds': elapsed, 'peak_rss_kib': int(rss), 'user_seconds': user,
                      'system_seconds': system, 'command_file': f'{label}-command.json'}
            if mode == 'replay':
                observed = {f: digest(out / label / f) for f in ['matrix.mtx', 'features.tsv', 'barcodes.tsv']}
                if replay_reference is None:
                    replay_reference = observed
                assert observed == replay_reference, (label, observed, replay_reference)
                result['matrix_sha256'] = observed
            else:
                data = json.loads(metrics.read_text())
                observed, coverage = scientific(data)
                if reference is None:
                    assert mode == 'nine'
                    reference = (observed, coverage)
                assert coverage == reference[1], (label, coverage, reference[1])
                assert observed == {name: reference[0][name] for name in observed}, label
                expected = 9 if mode == 'nine' else 4 if mode == 'four' else 1
                assert len(observed) == expected, label
                result['modes_verified'] = list(observed)
                result['support_storage'] = data.get('support_storage')
            manifest['runs'].append(result)
            (out / 'results.json').write_text(json.dumps(manifest, indent=2) + '\n')
            print(label, f'{elapsed:.2f}s', f'{rss / 1024**2:.3f} GiB', 'equality passed', flush=True)
    manifest['summary'] = {name: {
        'median_seconds': statistics.median(r['elapsed_seconds'] for r in manifest['runs'] if r['case'] == name),
        'median_peak_rss_kib': statistics.median(r['peak_rss_kib'] for r in manifest['runs'] if r['case'] == name),
    } for name, _, _ in cases}
    manifest['passed'] = True
    (out / 'results.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
