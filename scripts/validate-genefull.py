#!/usr/bin/env python3
"""Synthetic STARsolo and representation checks; writes only a new output directory.

Requires aie and STAR 2.7.11b. Uses deterministic synthetic sequence,
not project data. Compare identical annotation-aware alignments first; the
biological comparison separately measures annotation-free alignment replay.
"""
import argparse
import json
from pathlib import Path
import random
import subprocess


def run(argv, log):
    with log.open('w') as f:
        subprocess.run(list(map(str, argv)), stdout=f, stderr=subprocess.STDOUT, check=True)


def matrix(path):
    genes=[line.split('\t')[0] for line in (path/'features.tsv').read_text().splitlines()]
    barcodes=(path/'barcodes.tsv').read_text().splitlines()
    rows=[line for line in (path/'matrix.mtx').read_text().splitlines() if not line.startswith('%')]
    shape=list(map(int,rows[0].split()))
    assert shape[:2]==[len(genes),len(barcodes)]
    result={}
    for row in rows[1:]:
        g,c,n=map(int,row.split())
        if n: result[(genes[g-1],barcodes[c-1])]=n
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--aie',type=Path,required=True)
    ap.add_argument('--star',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    rng=random.Random(991)
    genome=''.join(rng.choice('ACGT') for _ in range(6000))
    # Canonical splice boundaries for A and C; D has disjoint isoforms.
    for start,end in [(200,400),(700,900)]:
        genome=genome[:start]+'GT'+genome[start+2:end-2]+'AG'+genome[end:]
    fasta=out/'genome.fa';fasta.write_text('>chr1\n'+genome+'\n')
    annotation=[('A','A1','+',[(100,200),(400,500)]),('B','B1','+',[(250,300)]),
                ('C','C1','-',[(600,700),(900,1000)]),
                ('D','D1','+',[(1200,1250)]),('D','D2','+',[(1500,1550)]),
                ('E','E1','+',[(2000,2100)])]
    gtf=out/'genes.gtf'
    gtf.write_text(''.join(f'chr1\tfixture\texon\t{s+1}\t{e}\t.\t{strand}\t.\tgene_id "{g}"; transcript_id "{t}"; gene_name "{g}";\n'
                           for g,t,strand,exons in annotation for s,e in exons))
    cb='ACGTACGTACGTACGT';wl=out/'barcodes.tsv';wl.write_text(cb+'\n')
    # Separate molecular tags (minimum Hamming distance >=3), plus one intended 1MM edge.
    umis=[]
    while len(umis)<30:
        u=''.join(rng.choice('ACGT') for _ in range(12))
        if all(sum(a!=b for a,b in zip(u,v))>=3 for v in umis): umis.append(u)
    specs=[(110,50,False),(210,30,False),(260,30,False),(480,50,False),
           (1000,50,False),(1300,50,False),(710,50,True),(620,50,True),
           (2050,50,False),(350,30,False)]
    reads=[(genome[s:s+n],umis[i],rev) for i,(s,n,rev) in enumerate(specs)]
    reads.extend([(genome[150:200]+genome[400:450],umis[10],False),
                  (genome[140:190]+genome[410:460],umis[11],False)])
    # Repeated geometries on one molecule and a less abundant adjacent UMI.
    reads.extend([(genome[110:160],umis[0],False)]*3)
    last='A' if umis[0][-1]!='A' else 'C'
    reads.append((genome[120:170],umis[0][:-1]+last,False))
    r1=out/'r1.fq';r2=out/'r2.fq'
    complement=str.maketrans('ACGT','TGCA')
    r1.write_text(''.join(f'@r{i}\n{cb}{u}\n+\n'+('I'*28)+'\n' for i,(_,u,_) in enumerate(reads)))
    r2.write_text(''.join(f'@r{i}\n{seq.translate(complement)[::-1] if rev else seq}\n+\n'+('I'*len(seq))+'\n' for i,(seq,_,rev) in enumerate(reads)))
    index=out/'index';index.mkdir()
    run([args.star,'--runMode','genomeGenerate','--runThreadN','2','--genomeDir',index,
         '--genomeFastaFiles',fasta,'--sjdbGTFfile',gtf,'--sjdbOverhang','49',
         '--genomeSAindexNbases','3','--genomeChrBinNbits','12','--outFileNamePrefix',str(out/'index-')],out/'index.log')
    checks=[]
    for strand in ['Forward','Reverse','Unstranded']:
        star=out/strand;star.mkdir()
        run([args.star,'--genomeDir',index,'--runThreadN','2','--readFilesIn',r2,r1,
             '--soloType','CB_UMI_Simple','--soloCBwhitelist',wl,'--soloCBlen','16','--soloUMIstart','17','--soloUMIlen','12',
             '--soloBarcodeReadLength','0','--soloFeatures','Gene','GeneFull','--soloStrand',strand,
             '--soloUMIdedup','1MM_CR','--soloUMIfiltering','MultiGeneUMI_CR','--soloCellFilter','None',
             '--outSAMtype','BAM','SortedByCoordinate','--outSAMattributes','NH','HI','AS','nM','CR','CY','UR','UY',
             '--outFileNamePrefix',str(star)+'/', '--outFilterMatchNmin','20',
             '--outFilterScoreMinOverLread','0','--outFilterMatchNminOverLread','0'],star/'command.log')
        bam=star/'Aligned.sortedByCoord.out.bam';archive=star/'fixture.aie'
        run([args.aie,'ingest-archive',bam,'--whitelist',wl,'--out',archive,'--zstd-level','1'],star/'ingest.log')
        compiled=star/'genes.aic'
        run([args.aie,'compile-annotation',gtf,'--out',compiled],star/'compile.log')
        for model in ['Gene','GeneFull']:
            expected=matrix(star/'Solo.out'/model/'raw')
            variants=[('stream',archive,gtf,[]),('eager',archive,gtf,['--eager']),
                      ('compiled',archive,compiled,[]),('bam',bam,gtf,['--from-bam','--whitelist',wl])]
            outputs=[]
            for name,source,anno,extra in variants:
                dest=star/f'{model}-{name}';report=star/f'{model}-{name}.json'
                run([args.aie,'replay-rows',source,'--gtf',anno,'--barcodes',wl,'--out-dir',dest,
                     '--solo-strand',strand.lower(),'--report-format','json','--report-output',report,
                     *(['--gene-full'] if model=='GeneFull' else []),*extra],star/f'{model}-{name}.log')
                observed=matrix(dest)
                assert observed==expected,(strand,model,name,observed,expected)
                outputs.append(dest)
                payload=json.loads(report.read_text());summary=payload['data']['summary']
                stats=summary['assignment_statistics']
                assert stats['assigned_molecule_records']<=stats['molecule_records']
                assert stats['assigned_representative_rows']<=stats['representative_rows']
                assert stats['assigned_umi_classes']<=stats['umi_classes']
                assert summary['assigned_molecules']==stats['assigned_molecule_records']
                metadata=json.loads((dest/'metadata.json').read_text())
                assert metadata['provenance']['parameters']['counting_model']==model
                checks.append({'strand':strand,'model':model,'input':name,'umis':sum(observed.values()),'statistics':stats})
            for dest in outputs[1:]:
                for file in ['matrix.mtx','features.tsv','barcodes.tsv']:
                    assert (dest/file).read_bytes()==(outputs[0]/file).read_bytes()
        for incompatible in ['--velocity','--audit-multigene']:
            result=subprocess.run([str(args.aie),'replay-rows',str(archive),'--gtf',str(gtf),'--barcodes',str(wl),
                                   '--out-dir',str(star/'must-not-exist'),'--gene-full',incompatible],capture_output=True)
            assert result.returncode!=0 and b'cannot be used with' in result.stderr
            assert not (star/'must-not-exist').exists()
    (out/'summary.json').write_text(json.dumps({'checks':checks,'passed':True},indent=2)+'\n')
    print(f'{len(checks)} STARsolo/representation comparisons passed')

if __name__=='__main__': main()
