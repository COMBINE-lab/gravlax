#!/usr/bin/env python3
"""Repeat release benchmark GQ sources through all execution paths (not timed)."""
import argparse
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import tempfile
from gq_report_summary import compact_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--benchmark', type=Path, required=True)
    parser.add_argument('--demo-dir', type=Path, required=True)
    parser.add_argument('--archive-dir', type=Path, required=True)
    parser.add_argument('--stress-archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must not exist')
    semantics = runpy.run_path(str(Path(__file__).with_name('benchmark-gq-execution.py')))['semantic_result']
    manifest = json.loads((args.demo_dir / 'demo-manifest.json').read_text())
    binary = str(args.binary.resolve(strict=True))
    checked = []
    for key, entry in manifest['resources'].items():
        if not key.startswith('archive_'):
            continue
        path = args.demo_dir / entry['filename']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256']
        result = subprocess.run([binary, 'inspect-archive', str(path), '--verify-content', '--json'], capture_output=True, check=True)
        checked.append(json.loads(result.stdout))
    results = []
    with tempfile.TemporaryDirectory(prefix='gravlax-gq-acceptance-') as directory:
        tmp = Path(directory)
        federation = tmp / 'cohort.json'
        federation.write_text(json.dumps({'schema_version': 1, 'assembly': 'GRCh38', 'archives': [
            {'sample': key, 'path': str((args.demo_dir / entry['filename']).resolve())}
            for key, entry in manifest['resources'].items() if key.startswith('archive_')]}))
        bindings = {'cd45_compact': args.archive_dir / 'indexed.aie', 'cd45_fidelity': args.archive_dir / 'fidelity.aie',
                    'sez_g': args.demo_dir / 'sez-donor-g.demo.aie', 'sez_eight': federation, 'synthetic_scale': args.stress_archive}
        report = json.loads(args.benchmark.read_text())
        for case in report['experiments']:
            if not case.get('query'):
                continue
            name = case['case'] if 'case' in case else case['name']
            source = tmp / 'query.gq'
            source.write_text(case['query'])
            base = [binary, 'gq', 'run', str(source), '--bind', f"x={bindings[name.split('/')[0]].resolve()}",
                    '--allow-full-scan', '--max-records', '5000000', '--max-rows', '1000000', '--max-steps', '500000000']
            expected = None
            for mode, flags in [('auto', []), ('parallel', ['--parallel-decode']), ('reference', ['--engine', 'reference'])]:
                proc = subprocess.run(base + flags, capture_output=True, check=True)
                data = semantics(json.loads(proc.stdout))
                if expected is None:
                    expected = data
                assert data == expected, (name, mode, 'result/work mismatch')
            results.append({'case': name, 'modes': ['auto', 'parallel', 'reference'],
                            'result_sha256': hashlib.sha256(json.dumps(expected, sort_keys=True).encode()).hexdigest(),
                            'summary': compact_summary(expected['summary'])})
            print(name, 'all execution paths agree', flush=True)
    with args.output.open('x') as stream:
        json.dump({'binary_sha256': hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
                   'verified_archives': checked, 'queries': results}, stream, indent=2)


if __name__ == '__main__':
    main()
