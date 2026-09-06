#!/usr/bin/env python3
"""Exercise new geometry semantics on the existing CD45 archives and save exact commands."""
import argparse
import json
import os
from pathlib import Path
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True)
    p.add_argument('--archive', action='append', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    atoms = ['u=region:1:198600000-198800000',
             'a=junction:1:198692373-198696711:-',
             'b=junction:1:198696909-198699563:-',
             'p=path:1:198692373-198696711,198696909-198699563:-',
             'n=junction-near:1:198692373-198696711:-/1',
             'o=overlap:1:198696711-198696909:-/12',
             's=start:1:198696711-198696712:-',
             'e=end:1:198696909-198696910:-']
    cases = [('junction_pair_record', 'a & b', 'record'),
             ('junction_pair_one_placement', 'a & b', 'any-placement'),
             ('junction_pair_all_placements', 'a | b', 'all-placements'),
             ('consecutive_path', 'p', 'record'),
             ('exact_junction', 'a', 'record'),
             ('near_junction_1bp', 'n', 'record'),
             ('overlap_12bp', 'o', 'record'),
             ('left_start_boundary', 's', 'record'),
             ('right_end_boundary', 'e', 'record')]
    results = []
    for archive in args.archive:
        for name, expression, within in cases:
            command = [args.binary, 'query', archive, 'cooccur', '--where', expression,
                       '--universe', 'u', '--format', 'json', '--agg', 'bulk', '--match-within', within]
            for atom in atoms:
                command += ['--predicate', atom]
            expected = None
            for engine in ['scalar', 'compiled', 'auto']:
                run = subprocess.run([*command, '--engine', engine], env=dict(os.environ, RAYON_NUM_THREADS='4'),
                                     check=True, capture_output=True)
                data = json.loads(run.stdout)['data']
                if expected is None:
                    expected = data
                assert data == expected, (archive, name, engine)
            s = expected['summary']
            row = {'archive': archive, 'case': name, 'command': command,
                   **{key: s[key] for key in ['candidate_units', 'selected_units', 'indeterminate_units', 'chunks_read']}}
            results.append(row)
            print(json.dumps({k: v for k, v in row.items() if k != 'command'}), flush=True)
    Path(args.output).write_text(json.dumps({'results': results, 'all_engines_identical': True}, indent=2) + '\n')


if __name__ == '__main__':
    main()
