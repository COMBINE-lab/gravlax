#!/usr/bin/env python3
"""Offline acceptance of an actual installed GQ executable, including a federation."""
import argparse
import base64
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    args = parser.parse_args()
    binary = str(args.binary.resolve(strict=True))
    encoded = Path(__file__).with_name('fixtures') / 'gq-smoke.aie.b64'
    with tempfile.TemporaryDirectory(prefix='gravlax-gq-smoke-') as tmp:
        tmp = Path(tmp)
        archive = tmp / 'fixture.aie'
        archive.write_bytes(base64.b64decode(b''.join(encoded.read_bytes().split()), validate=True))
        second = tmp / 'second.aie'
        second.write_bytes(base64.b64decode(b''.join(encoded.with_name('gq-smoke-coarse.aie.b64').read_bytes().split()), validate=True))
        federation = tmp / 'cohort.json'
        federation.write_text(json.dumps({'schema_version': 1, 'assembly': 'test',
            'archives': [{'sample': 'A', 'path': 'fixture.aie'}, {'sample': 'B', 'path': 'second.aie'}]}))
        query = tmp / 'count.gq'
        query.write_text('header {gq=1,assembly="test"}\nfrom @x.records |> within all '
                         '|> derive {j=any unique{j[1:160..180]}} '
                         '|> summarize {records=count(),positive=count_true(j),reads=sum(reads.total)}')
        env = dict(os.environ, RAYON_NUM_THREADS='4')
        subprocess.run([binary, 'gq', 'validate', str(query)], env=env, check=True, capture_output=True)
        for source, multiple in [(archive, 1), (federation, 2)]:
            for flags in [[], ['--parallel-decode'], ['--engine', 'reference']]:
                output = subprocess.run([binary, 'gq', 'run', str(query), '--bind', f'x={source}',
                    '--allow-full-scan', *flags], env=env, check=True, capture_output=True)
                document = json.loads(output.stdout)
                table = next(t for t in document['data']['tables'] if t['name'] == 'results')
                columns = [f['name'] for f in table['schema']['fields']]
                rows = [dict(zip(columns, row)) for row in table['rows']]
                expected = [{'records': 4 * multiple, 'positive': multiple, 'reads': 17 * multiple}]
                if rows != expected:
                    raise RuntimeError((source, flags, rows, expected))
        output = tmp / 'result.json'
        command = [binary, 'gq', 'run', str(query), '--bind', f'x={archive}', '--allow-full-scan', '--output', str(output)]
        subprocess.run(command, env=env, check=True, capture_output=True)
        original = output.read_bytes()
        again = subprocess.run(command, env=env, capture_output=True)
        if again.returncode == 0 or output.read_bytes() != original:
            raise RuntimeError('GQ publication replaced an existing result')
    print('GQ executable acceptance passed: archive/federation, serial/parallel/reference, no-clobber publication')


if __name__ == '__main__':
    main()
