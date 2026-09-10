#!/usr/bin/env python3
"""Matched Gene/GeneFull comparisons with explicit units and nucleus-set contrasts.

Inputs are preserved. Outputs include integer mass differences (no assumption
that half the L1 distance counts physically relocated molecules), matched
per-nucleus summaries, and a deterministic shared PCA/k-means sensitivity analysis.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import io, sparse
from scipy.sparse.linalg import svds
from scipy.cluster.vq import kmeans2


def load_selected(path, selected):
    genes=[line.split('\t')[0] for line in (path/'features.tsv').read_text().splitlines()]
    names=[line.split('\t')[1] for line in (path/'features.tsv').read_text().splitlines()]
    assert len(set(genes))==len(genes)
    ids={b:i for i,b in enumerate(selected)}
    columns=np.fromiter((ids.get(line.strip(),-1) for line in (path/'barcodes.tsv').open()),dtype=np.int32)
    assert sum(columns>=0)==len(selected),'missing or repeated selected barcodes'
    m=io.mmread(path/'matrix.mtx').tocoo()
    assert m.shape==(len(genes),len(columns))
    assert np.all(m.data>=0) and np.all(m.data==np.floor(m.data))
    raw_total=int(m.data.sum()); raw_nnz=int(np.count_nonzero(m.data))
    cols=columns[m.col];keep=cols>=0
    result=sparse.coo_matrix((m.data[keep].astype(np.int64),(m.row[keep],cols[keep])),shape=(len(genes),len(selected))).tocsr()
    result.sum_duplicates();result.eliminate_zeros()
    return result,genes,names,{'raw_collapsed_umis':raw_total,'raw_matrix_nonzero_entries':raw_nnz,'raw_barcode_columns':len(columns)}


def distribution(values):
    return {'min':float(np.min(values)),'median':float(np.median(values)),
            'p10':float(np.quantile(values,.1)),'p90':float(np.quantile(values,.9)),
            'max':float(np.max(values)),'mean':float(np.mean(values))}


def compare(a,b):
    d=(b-a).tocoo();gain=int(d.data[d.data>0].sum());loss=int(-d.data[d.data<0].sum())
    total_a=int(a.sum());total_b=int(b.sum())
    ac=np.asarray(a.sum(axis=0)).ravel();bc=np.asarray(b.sum(axis=0)).ravel()
    l1=np.asarray(abs(b-a).sum(axis=0)).ravel()
    return {'nuclei':a.shape[1],'gene_rows':a.shape[0],
            'a_collapsed_umis':total_a,'b_collapsed_umis':total_b,
            'positive_count_difference_umis':gain,'negative_count_difference_umis':loss,
            'net_count_difference_umis':total_b-total_a,'l1_count_difference_umis':gain+loss,
            'l1_over_a_umis':(gain+loss)/max(1,total_a),
            'half_l1_over_a_umis':(gain+loss)/(2*max(1,total_a)),
            'changed_gene_by_nucleus_entries':int(np.count_nonzero(d.data)),
            'possible_gene_by_nucleus_entries':int(a.shape[0]*a.shape[1]),
            'changed_nuclei':int(np.count_nonzero(l1)),
            'a_umis_per_nucleus':distribution(ac),'b_umis_per_nucleus':distribution(bc),
            'a_genes_per_nucleus':distribution(np.asarray((a>0).sum(axis=0)).ravel()),
            'b_genes_per_nucleus':distribution(np.asarray((b>0).sum(axis=0)).ravel()),
            'paired_umi_ratio_b_over_a':distribution(bc/np.maximum(ac,1))}


def ari(a,b):
    _,aa=np.unique(a,return_inverse=True);_,bb=np.unique(b,return_inverse=True)
    c=np.zeros((aa.max()+1,bb.max()+1),dtype=np.int64);np.add.at(c,(aa,bb),1)
    choose=lambda x:np.sum(x*(x-1)/2)
    n=choose(c);sa=choose(c.sum(axis=1));sb=choose(c.sum(axis=0));pairs=len(a)*(len(a)-1)/2
    expected=sa*sb/pairs;denom=(sa+sb)/2-expected
    return float((n-expected)/denom) if denom else 1.


def downstream(matrices,genes,names,barcodes,out):
    labels=list(matrices)
    lengths=[m.shape[1] for m in matrices.values()]
    normalized=[]
    for m in matrices.values():
        x=m.T.astype(np.float64).tocsr();tot=np.asarray(x.sum(axis=1)).ravel()
        x=sparse.diags(1e4/np.maximum(tot,1))@x;x.data=np.log1p(x.data)
        normalized.append(x)
    pooled=sparse.vstack(normalized).tocsr()
    mean=np.asarray(pooled.mean(axis=0)).ravel()
    var=np.asarray(pooled.power(2).mean(axis=0)).ravel()-mean**2
    selected=np.argsort(var,kind='stable')[-2000:]
    x=pooled[:,selected].toarray();x-=x.mean(axis=0);x/=np.maximum(x.std(axis=0),1e-8)
    u,s,_=svds(x,k=30,random_state=314159)
    scores=u[:,::-1]*s[::-1]
    _,pooled_cluster=kmeans2(scores,15,iter=100,minit='++',seed=2718)
    assignments={};start=0
    for label,n in zip(labels,lengths):
        assignments[label]=pooled_cluster[start:start+n]
        start+=n
    result={'method':'joint log1p(10000 library-normalized counts), pooled 2000 highest-variance genes, gene z-scores, 30 joint PCs, shared pooled k-means k=15, seed 2718',
            'scope':'fixed Gene-called nuclei; exploratory sensitivity, not cell-type accuracy',
            'pairwise_adjusted_rand_index':{a+' vs '+b:ari(assignments[a],assignments[b]) for i,a in enumerate(labels) for b in labels[i+1:]}}
    (out/'downstream.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'cluster-assignments.tsv').open('w') as f:
        f.write('barcode\t'+'\t'.join(labels)+'\n')
        for i,bc in enumerate(barcodes):f.write(bc+'\t'+'\t'.join(str(assignments[l][i]) for l in labels)+'\n')
    with (out/'cluster-markers.tsv').open('w') as f:
        f.write('matrix\tcluster\tnuclei\tgene_id\tgene_name\tmean_log1p_in\tmean_log1p_out\tdifference\n')
        for label,norm in zip(labels,normalized):
            cluster=assignments[label]
            for k in range(15):
                mask=cluster==k
                if not mask.any() or mask.all():continue
                inside=np.asarray(norm[mask].mean(axis=0)).ravel();outside=np.asarray(norm[~mask].mean(axis=0)).ravel()
                difference=inside-outside
                for g in np.argsort(difference)[-10:][::-1]:
                    f.write(f'{label}\t{k}\t{mask.sum()}\t{genes[g]}\t{names[g]}\t{inside[g]}\t{outside[g]}\t{difference[g]}\n')
    panel=['SLC17A7','GAD1','GAD2','AQP4','GFAP','MBP','PLP1','PDGFRA','CSF1R','PTPRC','CLDN5','FLT1']
    with (out/'marker-detection.tsv').open('w') as f:
        f.write('matrix\tgene_id\tgene_name\tdetected_nuclei\tdenominator_nuclei\ttotal_umis\n')
        for label,m in matrices.items():
            for name in panel:
                for g in [i for i,n in enumerate(names) if n==name]:
                    f.write(f'{label}\t{genes[g]}\t{name}\t{m[g].count_nonzero()}\t{m.shape[1]}\t{int(m[g].sum())}\n')
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reference',type=Path,required=True,help='STAR Solo.out directory')
    ap.add_argument('--gene',type=Path,required=True)
    ap.add_argument('--gene-full',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--skip-downstream',action='store_true')
    args=ap.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=False)
    called={m:set((args.reference/m/'filtered/barcodes.tsv').read_text().splitlines()) for m in ['Gene','GeneFull']}
    union=sorted(called['Gene']|called['GeneFull']);idx={b:i for i,b in enumerate(union)}
    sets={'fixed_Gene_called':sorted(called['Gene']),'intersection_called':sorted(called['Gene']&called['GeneFull']),
          'GeneFull_called':sorted(called['GeneFull']),'Gene_only_called':sorted(called['Gene']-called['GeneFull']),
          'GeneFull_only_called':sorted(called['GeneFull']-called['Gene'])}
    for name,bc in sets.items():(out/(name+'.barcodes.tsv')).write_text(''.join(b+'\n' for b in bc))
    paths={'STAR_Gene':args.reference/'Gene/raw','STAR_GeneFull':args.reference/'GeneFull/raw',
           'replay_Gene':args.gene,'replay_GeneFull':args.gene_full}
    matrices={};raw={};gene_order=None;gene_names=None
    for label,path in paths.items():
        m,genes,names,stats=load_selected(path,union)
        if gene_order is None:gene_order=genes;gene_names=names
        elif genes!=gene_order:
            assert set(genes)==set(gene_order)
            positions={g:i for i,g in enumerate(genes)};m=m[[positions[g] for g in gene_order]]
        matrices[label]=m;raw[label]=stats;sparse.save_npz(out/(label+'.npz'),m)
    (out/'genes.tsv').write_text(''.join(g+'\t'+n+'\n' for g,n in zip(gene_order,gene_names)))
    (out/'union.barcodes.tsv').write_text(''.join(b+'\n' for b in union))
    results={'counting_units':{'matrix':'collapsed gene-specific UMI counts','nucleus':'distinct called barcode','mass_difference':'differences between integer matrix entries; not tracked molecular identities'},
             'calling':{'Gene_called':len(called['Gene']),'GeneFull_called':len(called['GeneFull']),
                        'intersection':len(sets['intersection_called']),'union':len(union),
                        'Gene_only':len(sets['Gene_only_called']),'GeneFull_only':len(sets['GeneFull_only_called']),
                        'jaccard_intersection_over_union':len(sets['intersection_called'])/len(union)},'raw':raw,'comparisons':{}}
    for scope,bc in sets.items():
        if not bc:continue
        cols=[idx[b] for b in bc]
        selected={label:m[:,cols] for label,m in matrices.items()}
        for a,b in [('STAR_Gene','replay_Gene'),('STAR_GeneFull','replay_GeneFull'),('replay_Gene','replay_GeneFull'),('STAR_Gene','STAR_GeneFull')]:
            results['comparisons'][scope+'/'+a+'_to_'+b]=compare(selected[a],selected[b])
        if scope=='fixed_Gene_called':
            if not args.skip_downstream:results['downstream']=downstream(selected,gene_order,gene_names,bc,out)
            with (out/'per-nucleus.tsv').open('w') as f:
                f.write('barcode\t'+'\t'.join(label+'_umis\t'+label+'_genes' for label in selected)+'\n')
                arrays=[(np.asarray(m.sum(axis=0)).ravel(),np.asarray((m>0).sum(axis=0)).ravel()) for m in selected.values()]
                for i,barcode in enumerate(bc):f.write(barcode+'\t'+'\t'.join(str(int(x[i])) for pair in arrays for x in pair)+'\n')
    (out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps({'calling':results['calling'],'fixed_nucleus_comparisons':{k:v for k,v in results['comparisons'].items() if k.startswith('fixed')}},indent=2))

if __name__=='__main__':main()
