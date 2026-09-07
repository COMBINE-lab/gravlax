#!/usr/bin/env python3
"""Offline executable acceptance: relocate an unchanged layered collection bundle."""
import argparse
import base64
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--relocation-dir', type=Path, help='Optional second filesystem for the copied bundle')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    if args.report and args.report.exists():
        parser.error('report must not exist')
    binary = str(args.binary.resolve(strict=True))
    env = dict(os.environ, RAYON_NUM_THREADS='4')
    def run(command, succeeds=True):
        p = subprocess.run([binary, *map(str, command)], env=env, capture_output=True, text=True)
        if (p.returncode == 0) != succeeds:
            raise RuntimeError((command, p.returncode, p.stdout, p.stderr))
        return json.loads(p.stdout) if succeeds else p.stderr
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory(prefix='gravlax-relocation-original-') as first, tempfile.TemporaryDirectory(prefix='gravlax-relocation-copy-', dir=args.relocation_dir) as second:
        root, moved = Path(first) / 'bundle', Path(second) / 'bundle'
        root.mkdir()
        fixtures = Path(__file__).with_name('fixtures')
        for sample, filename in [('a', 'gq-smoke.aie.b64'), ('b', 'gq-smoke-coarse.aie.b64'), ('c', 'gq-smoke-middle.aie.b64')]:
            (root / f'{sample}.aie').write_bytes(base64.b64decode(b''.join((fixtures / filename).read_bytes().split()), validate=True))
        run(['collection', 'build', '--sample', f'a={root / "a.aie"}', '--sample', f'b={root / "b.aie"}', '--allow-unstamped', '--shape-routes', '--out', root / 'base.aicollection', '--json'])
        run(['collection', 'build', '--base', root / 'base.aicollection', '--sample', f'c={root / "c.aie"}', '--allow-unstamped', '--shape-routes', '--out', root / 'atlas.aicollection', '--json'])
        inspection = run(['collection', 'inspect', root / 'atlas.aicollection'])
        entries = [{'identity': f'{a["native_identity"]["scheme"]}:{a["native_identity"]["blake3"]}', 'path': f'{a["id"]}.aie'} for a in inspection['archives']]
        entries += [{'identity': f'aicollection-directory-root-v1:{layer["root_digest"]}', 'path': Path(layer['path']).name} for layer in inspection['layers']]
        (root / 'locations.json').write_text(json.dumps({'schema_version': 1, 'locations': entries}))
        commands = {
            'inspect': ['inspect', '{collection}'],
            'verify_routes': ['inspect', '{collection}', '--verify-routes'],
            'verify_content': ['inspect', '{collection}', '--verify-content'],
            'region': ['region', '{collection}', '1:0-1000', '--top', '0', '--format', 'json'],
            'junction': ['junction', '{collection}', '1:160-180', '--top', '0', '--format', 'json'],
            'jset': ['jset', '{collection}', '--include', '1:160-180', '--exclude', '1:999-1000', '--top', '0', '--format', 'json'],
            'find_events': ['find-events', '{collection}', '--kind', 'junction', '--min-support', '0', '--format', 'json'],
        }
        def execute(where, with_locations):
            result, timings = {}, {}
            for name, cmd in commands.items():
                command = ['collection', *[str(where / 'atlas.aicollection') if x == '{collection}' else x for x in cmd]]
                if with_locations:
                    command += ['--locations', str(where / 'locations.json')]
                start = time.perf_counter()
                result[name] = run(command)
                timings[name] = (time.perf_counter() - start) * 1000
            return result, timings
        before, before_ms = execute(root, False)
        hashes = {p.name: sha(p) for p in root.iterdir() if p.suffix != '.json'}
        shutil.copytree(root, moved)
        different_filesystems = root.stat().st_dev != moved.stat().st_dev
        # Original paths become unavailable. Keep the original bytes for the hash comparison.
        root.rename(Path(first) / 'original-unavailable')
        run(['collection', 'inspect', moved / 'atlas.aicollection'], succeeds=False)
        after, after_ms = execute(moved, True)
        def normalized(value):
            if isinstance(value, dict):
                return {k: normalized(v) for k, v in value.items() if not k.endswith('_seconds') and k not in ('timings_ms', 'locations_manifest')}
            if isinstance(value, list):
                return [normalized(v) for v in value]
            if isinstance(value, str):
                return value.replace(str(root), '<bundle>').replace(str(moved), '<bundle>')
            return value
        for name in commands:
            if normalized(before[name]) != normalized(after[name]):
                raise RuntimeError((name, 'relocation changed results or I/O', normalized(before[name]), normalized(after[name])))
        if hashes != {name: sha(moved / name) for name in hashes}:
            raise RuntimeError('relocation changed archive/collection bytes')
        # An incremental build can reuse a relocated base without rescanning its indexes.
        extended = moved / 'extended.aicollection'
        run(['collection', 'build', '--base', moved / 'base.aicollection', '--locations', moved / 'locations.json',
             '--sample', f'c={moved / "c.aie"}', '--allow-unstamped', '--shape-routes', '--out', extended, '--json'])
        extension = run(['collection', 'inspect', extended, '--locations', moved / 'locations.json'])
        if len(extension['archives']) != 3:
            raise RuntimeError('relocated incremental build lost sources')
        run(['collection', 'build', '--sample', f'a={moved / "a.aie"}', '--sample', f'b={moved / "b.aie"}',
             '--locations', moved / 'locations.json', '--allow-unstamped', '--out', moved / 'unused-map.aicollection', '--json'], succeeds=False)
        if (moved / 'unused-map.aicollection').exists():
            raise RuntimeError('invalid location option published a collection')
        # Existing path plus an identical replacement inode must work without a mapping too.
        standalone = moved / 'standalone.aicollection'
        run(['collection', 'build', '--sample', f'a={moved / "a.aie"}', '--sample', f'b={moved / "b.aie"}', '--allow-unstamped', '--out', standalone, '--json'])
        replacement = moved / 'replacement.aie'
        shutil.copyfile(moved / 'a.aie', replacement)
        os.replace(replacement, moved / 'a.aie')
        run(['collection', 'inspect', standalone])
        # Wrong hints are authoritative failures, never silently ignored in favor of old paths.
        wrong = {'schema_version': 1, 'locations': [{**entries[0], 'path': 'b.aie'}]}
        (moved / 'wrong.json').write_text(json.dumps(wrong))
        run(['collection', 'inspect', standalone, '--locations', moved / 'wrong.json'], succeeds=False)
        report = {'cases': list(commands), 'different_filesystems': different_filesystems,
                  'archive_and_collection_sha256': hashes, 'before_ms': before_ms, 'after_ms': after_ms,
                  'note': 'Single-pass acceptance timings, not a statistical speed benchmark. Complete results and I/O counters agree after removing timings and locator provenance.',
                  'verified_results': normalized(after)}
        if args.report:
            with args.report.open('x') as out:
                json.dump(report, out, indent=2)
    print('Collection relocation acceptance passed: unchanged bytes, layered roots, query/route parity, identical archive I/O, wrong-root rejection')


if __name__ == '__main__':
    main()
