#!/usr/bin/env python3
"""Release allocator assessment: matched GQ, native query, replay and ingest."""
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
from gq_report_summary import compact_summary

HERE = Path(__file__).resolve().parent
semantics = runpy.run_path(str(HERE / 'benchmark-gq-execution.py'))['semantic_result']


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--system', type=Path, required=True)
    p.add_argument('--mimalloc', type=Path, required=True)
    p.add_argument('--demo-dir', type=Path, required=True)
    p.add_argument('--cd45-dir', type=Path, required=True)
    p.add_argument('--archive-dir', type=Path, required=True)
    p.add_argument('--stress-archive', type=Path)
    p.add_argument('--case-filter', action='append', default=[])
    p.add_argument('--mimalloc-env', action='append', default=[], help='Explicit MIMALLOC_NAME=value tuning, recorded in the report')
    p.add_argument('--repeats', type=int, default=9)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--cpus', default='1,2,7,9')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.repeats < 3 or args.output.exists():
        p.error('use at least three repetitions and a new output file')
    os.sched_setaffinity(0, {int(n) for n in args.cpus.split(',')})
    binaries = {k: str(v.resolve(strict=True)) for k, v in [('system', args.system), ('mimalloc', args.mimalloc)]}
    demo = args.demo_dir.resolve(strict=True)
    manifest = json.loads((demo / 'demo-manifest.json').read_text())
    roots = {}
    for key, entry in manifest['resources'].items():
        if key.startswith('archive_'):
            archive = demo / entry['filename']
            if digest(archive) != entry['sha256']:
                raise RuntimeError(f'capsule checksum mismatch: {archive}')
            roots[str(archive)] = entry
    env = dict(os.environ, RAYON_NUM_THREADS=str(args.threads))
    tuning = dict(value.split('=', 1) for value in args.mimalloc_env)
    if any(not key.startswith('MIMALLOC_') for key in tuning):
        p.error('only MIMALLOC_ environment options are accepted')
    rng = random.Random(7069)
    report = {'method': 'One warmup, randomized paired allocator order, warm filesystem cache. GQ post_parse_compute_ms excludes parsing/startup/publication. Native commands use complete process wall time. Peak RSS is whole child process. Exact GQ semantic/work results, native data tables, replay files and ingest bytes are checked. Each pair supplies one speed ratio; percentile bootstrap CI uses resampling of paired median ratios, not independent samples. Serial GQ decoding is the release default; allocator tuning is explicitly recorded in mimalloc_tuning; no C malloc interposition.',
        'threads': args.threads, 'mimalloc_tuning': tuning, 'affinity': sorted(os.sched_getaffinity(0)), 'load_start': os.getloadavg(), 'repeats': args.repeats,
        'binaries': {k: {'path': v, 'sha256': digest(v)} for k, v in binaries.items()}, 'capsule_archives': roots, 'experiments': []}
    with tempfile.TemporaryDirectory(prefix='gravlax-release-allocator-') as tmp:
        tmp = Path(tmp)
        federation = tmp / 'cohort.json'
        federation.write_text(json.dumps({'schema_version': 1, 'assembly': 'GRCh38', 'archives': [
            {'sample': key, 'path': str(demo / entry['filename'])} for key, entry in manifest['resources'].items() if key.startswith('archive_')]}))
        cases = []
        for label, archive, region in [
            ('cd45_compact', args.archive_dir / 'indexed.aie', '1:198696711..198696909:-'),
            ('cd45_fidelity', args.archive_dir / 'fidelity.aie', '1:198696711..198696909:-'),
            ('sez_g', demo / 'sez-donor-g.demo.aie', 'chr1:70122500..70142500'),
            ('sez_eight', federation, 'chr1:70122500..70142500'),
            *([('synthetic_scale', args.stress_archive, '1:0..2000000')] if args.stress_archive else [])]:
            for name, body in {
                'geometry_summary': 'from @x.records |> within all |> derive {q=any unique{overlaps(R)},m=any multimap{any alternative{overlaps(R)}}} |> summarize {n=count(),q=count_true(q),m=count_true(m),reads=sum(reads.total)} by {sample}',
                'class_summary': 'from @x.classes |> within all |> derive {q=any record{any unique{overlaps(R)}}} |> summarize {n=count(),q=count_true(q),reads=sum(reads.total)} by {sample}',
                'cell_summary': 'from @x.records |> within all |> summarize {n=count(),reads=sum(reads.total)} by {cell.id,sample}',
                'support_top_k': 'from @x.records |> within all |> derive {s=support_reads(overlaps(R))} |> select {unit.id,s} |> sort {unit.id} |> take 20',
            }.items():
                source = f'header {{gq=1,assembly="GRCh38"}} let R=g[{region}] ' + body
                query = tmp / f'{label}-{name}.gq'; query.write_text(source)
                cases.append((f'{label}/{name}', 'gq', ['gq','run',str(query),'--bind',f'x={archive}','--allow-full-scan','--max-records','5000000','--max-rows','1000000','--max-steps','500000000','--profile'], source))
        for label, archive, locus in [
            ('cd45', args.archive_dir / 'indexed.aie', '1:198696711-198696909:-'),
            ('sez_g', demo / 'sez-donor-g.demo.aie', 'chr1:70122500-70142500')]:
            cases.append((label + '/cooccur', 'native', ['query',str(archive),'cooccur','--predicate',f'a=region:{locus}','--where','a','--universe','a','--placements','unique','--region-match','aligned-block','--allow-full-scan','--agg','bulk','--format','json'], None))
        cases.append(('cd45/replay','replay',['replay-rows',str(args.archive_dir / 'indexed.aie'),'--gtf',str(args.cd45_dir / 'CD45_exons_nochr.hg38.txt'),'--barcodes',str(args.cd45_dir / '737K-august-2016.txt'),'--solo-strand','unstranded','--out-dir','{out}'],None))
        cases.append(('cd45/ingest','ingest',['ingest-archive',str(args.cd45_dir / 'ptprc-grch38-complete-multiplicity.bam'),'--whitelist',str(args.cd45_dir / '737K-august-2016.txt'),'--out','{out}'],None))
        for name, kind, command, source in cases:
            if args.case_filter and not any(part in name for part in args.case_filter):
                continue
            samples = {k: [] for k in binaries}; expected = None; summary = None
            for repeat in range(args.repeats + 1):
                order = list(binaries); rng.shuffle(order)
                for label in order:
                    out = tmp / f'out-{name.replace("/","-")}-{repeat}-{label}'
                    cmd = [str(out) if a == '{out}' else a for a in command]
                    start = time.perf_counter_ns()
                    proc = subprocess.run(['/usr/bin/time','-f','allocator_rss=%M',binaries[label],*cmd],env=dict(env, **(tuning if label == 'mimalloc' else {})),capture_output=True)
                    if proc.returncode:
                        raise RuntimeError((name,label,proc.stderr.decode()))
                    wall = (time.perf_counter_ns() - start) / 1e6
                    lines = proc.stderr.decode().splitlines()
                    sample = {'process_ms': wall, 'rss_kib': next(int(s.split('=',1)[1]) for s in lines if s.startswith('allocator_rss='))}
                    if kind == 'gq':
                        document = json.loads(proc.stdout); data = semantics(document); summary = data['summary']
                        sample.update(next(json.loads(s.split('=',1)[1]) for s in lines if s.startswith('gq_profile=')))
                    elif kind == 'native':
                        document = json.loads(proc.stdout); data = document['data']; summary = data['summary']
                    elif kind == 'ingest':
                        data = digest(out)
                    else:
                        data = {str(f.relative_to(out)): digest(f) for f in sorted(out.rglob('*')) if f.is_file()}
                    if expected is None: expected = data
                    if data != expected: raise RuntimeError((name,label,'semantic/artifact mismatch'))
                    if repeat: samples[label].append(sample)
            clock = 'post_parse_compute_ms' if kind == 'gq' else 'process_ms'
            median = {k: {f: statistics.median(s[f] for s in values) for f in [clock,'rss_kib']} for k,values in samples.items()}
            ratios = [a[clock]/b[clock] for a,b in zip(samples['system'],samples['mimalloc'])]
            boot = sorted(statistics.median(rng.choices(ratios,k=len(ratios))) for _ in range(4000))
            experiment = {'case':name,'kind':kind,'query':source,'median':median,'paired_speedup':statistics.median(ratios),'paired_speedup_ci95':[boot[100],boot[3899]],'samples':samples,'verified_result_sha256':hashlib.sha256(json.dumps(expected,sort_keys=True).encode()).hexdigest(),'summary':summary}
            experiment['summary'] = compact_summary(summary)
            report['experiments'].append(experiment)
            print(name,json.dumps({k:experiment[k] for k in ['median','paired_speedup','paired_speedup_ci95']}),flush=True)
    report['load_end'] = os.getloadavg()
    with args.output.open('x') as output: json.dump(report,output,indent=2)


if __name__ == '__main__': main()
