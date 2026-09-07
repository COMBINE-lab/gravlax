#!/usr/bin/env python3
"""Paired original/relocated collection benchmark on the real eight-donor capsule."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--demo-dir', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=9)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 3 or args.output.exists():
        parser.error('require at least three repeats and a new report path')
    binary = str(args.binary.resolve(strict=True))
    env = dict(os.environ, RAYON_NUM_THREADS='4')
    rng = random.Random(907)
    def run(command):
        result = subprocess.run([binary, *map(str, command)], capture_output=True, env=env, check=True)
        return json.loads(result.stdout)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = json.loads((args.demo_dir / 'demo-manifest.json').read_text())
    archives = sorted((name, entry) for name, entry in manifest['resources'].items() if name.startswith('archive_'))
    if len(archives) != 8:
        raise RuntimeError('expected the eight-donor demo capsule')
    report = {'method': 'One warmup and nine (or --repeats) randomized paired repetitions; warm cache; four Rayon threads; same filesystem. Full process wall time includes startup, argument parsing and output. Complete science/work results agree after excluding timings and locator provenance. No archive payload verification or decoding is added by relocation.',
              'repeats': args.repeats, 'binary_sha256': sha(Path(binary)), 'archive_sha256': {}, 'cases': []}
    with tempfile.TemporaryDirectory(prefix='gravlax-locations-benchmark-') as temporary:
        tmp = Path(temporary); original = tmp / 'original'; original.mkdir()
        for name, entry in archives:
            source = args.demo_dir / entry['filename']
            if sha(source) != entry['sha256']:
                raise RuntimeError(('capsule checksum mismatch', source))
            shutil.copyfile(source, original / f'{name}.aie')
            report['archive_sha256'][name] = entry['sha256']
        for filename, members, base in [('base.aicollection', archives[:4], []), ('atlas.aicollection', archives[4:], ['--base', original / 'base.aicollection'])]:
            command = ['collection', 'build', *base, '--shape-routes', '--out', original / filename, '--json']
            for name, _ in members:
                command += ['--sample', f'{name}={original / (name + ".aie")}']
            run(command)
        inspected = run(['collection', 'inspect', original / 'atlas.aicollection'])
        discovered = run(['collection', 'find-events', original / 'atlas.aicollection', '--kind', 'junction', '--min-support', '0', '--format', 'json'])
        components = next(t for t in discovered['data']['tables'] if t['name'] == 'components')
        fields = [f['name'] for f in components['schema']['fields']]
        observed = [dict(zip(fields, row)) for row in components['rows']]
        observed.sort(key=lambda row: row['exact_umi_classes'], reverse=True)
        junctions = list(dict.fromkeys(f'{r["chrom"]}:{r["donor"]}-{r["acceptor"]}' for r in observed if r['donor'] is not None))
        if len(junctions) < 2:
            raise RuntimeError('benchmark needs two observed junctions')
        entries = [{'identity': f'{a["native_identity"]["scheme"]}:{a["native_identity"]["blake3"]}', 'path': f'{a["id"]}.aie'} for a in inspected['archives']]
        entries += [{'identity': f'aicollection-directory-root-v1:{l["root_digest"]}', 'path': Path(l['path']).name} for l in inspected['layers']]
        (original / 'locations.json').write_text(json.dumps({'schema_version': 1, 'locations': entries}))
        relocated = tmp / 'relocated'; shutil.copytree(original, relocated)
        report['collection_sha256'] = {p.name: sha(p) for p in original.glob('*.aicollection')}
        for name, digest in report['collection_sha256'].items():
            if sha(relocated / name) != digest:
                raise RuntimeError('collection copy differs')
        def clean(value):
            if isinstance(value, dict):
                return {k: clean(v) for k, v in value.items() if not k.endswith('_seconds') and k not in ('locations_manifest', 'timings_ms')}
            if isinstance(value, list): return [clean(v) for v in value]
            if isinstance(value, str): return value.replace(str(original), '<bundle>').replace(str(relocated), '<bundle>')
            return value
        for case, tail in [
            ('region_chr1', ['region', 'chr1:70122500-70142500']),
            ('region_chr12', ['region', 'chr12:18714000-18854000']),
            ('junction_observed', ['junction', junctions[0]]),
            ('jset_observed', ['jset', '--include', junctions[0], '--exclude', junctions[1]]),
        ]:
            values = {'original': [], 'relocated': []}; expected = None
            for repeat in range(args.repeats + 1):
                order = list(values); rng.shuffle(order)
                for mode in order:
                    root = original if mode == 'original' else relocated
                    command = ['collection', tail[0], root / 'atlas.aicollection', *tail[1:], '--top', '0', '--format', 'json']
                    if mode == 'relocated': command += ['--locations', relocated / 'locations.json']
                    start = time.perf_counter_ns(); result = run(command); elapsed = (time.perf_counter_ns() - start) / 1e6
                    observed = clean(result)
                    if expected is None: expected = observed
                    if expected != observed: raise RuntimeError((case, mode, 'result or I/O mismatch'))
                    if repeat: values[mode].append(elapsed)
            ratios = [b / a for a, b in zip(values['original'], values['relocated'])]
            boot = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(4000))
            report['cases'].append({'case': case, 'query': tail, 'process_ms': values,
                'median_ms': {k: statistics.median(v) for k, v in values.items()},
                'paired_relocated_over_original': statistics.median(ratios), 'paired_ratio_ci95': [boot[100], boot[3899]],
                'verified_result_sha256': hashlib.sha256(json.dumps(expected, sort_keys=True).encode()).hexdigest(),
                'summary': expected['data']['summary']})
            print(case, report['cases'][-1]['median_ms'], flush=True)
    with args.output.open('x') as output:
        json.dump(report, output, indent=2)


if __name__ == '__main__':
    main()
