#!/usr/bin/env python3
"""Summarize saved nucleus calling and annotation perturbations without overwrites."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import numpy as np
from scipy import io,sparse

spec=importlib.util.spec_from_file_location('comparison',Path(__file__).with_name('compare-genefull-brain.py'))
comparison=importlib.util.module_from_spec(spec);spec.loader.exec_module(comparison)


def stable_id(gene):
    return re.sub(r'^(ENSG\d+)\.\d+(_PAR_Y)?$',lambda m:m[1]+(m[2] or ''),gene)


def align(matrix,genes,union):
    assert len(set(genes))==len(genes),'duplicate normalized gene identifiers'
    mapping={g:i for i,g in enumerate(union)};coo=matrix.tocoo()
    row=np.array([mapping[g] for g in genes],dtype=np.int32)[coo.row]
    return sparse.coo_matrix((coo.data,(row,coo.col)),shape=(len(union),matrix.shape[1])).tocsr()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--reference',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=False);run=args.run
    fixed=(args.reference/'Gene/filtered/barcodes.tsv').read_text().splitlines()
    called={};sets={}
    for label,path in [(f'STAR_{m}',args.reference/m/'filtered') for m in ['Gene','GeneFull']]+[(f'replay_{m}',run/f'calling-replay-{m}') for m in ['Gene','GeneFull']]:
        bc=path.joinpath('barcodes.tsv').read_text().splitlines();sets[label]=set(bc)
        m=io.mmread(path/'matrix.mtx').tocsr()
        called[label]={'called_nuclei':len(bc),'collapsed_umis_in_called_nuclei':int(m.sum()),
                       'umis_per_called_nucleus':comparison.distribution(np.asarray(m.sum(axis=0)).ravel()),
                       'genes_per_called_nucleus':comparison.distribution(np.asarray((m>0).sum(axis=0)).ravel())}
    comparisons={}
    for a,b in [('replay_Gene','replay_GeneFull'),('STAR_Gene','STAR_GeneFull'),('STAR_Gene','replay_Gene'),('STAR_GeneFull','replay_GeneFull')]:
        x,y=sets[a],sets[b];comparisons[a+' to '+b]={'a_nuclei':len(x),'b_nuclei':len(y),'intersection':len(x&y),'union':len(x|y),
            'a_only':len(x-y),'b_only':len(y-x),'jaccard':len(x&y)/len(x|y)}
    (out/'nucleus-calling.json').write_text(json.dumps({'method':'STAR 2.7.11b EmptyDrops_CR default parameters; input each raw matrix separately','called':called,'comparisons':comparisons},indent=2)+'\n')
    effects={}
    for model in ['Gene','GeneFull']:
        base,basegenes,_,_=comparison.load_selected(run/f'replay-v49-{model}',fixed)
        for condition,otherpath in [('v32_to_v49',run/f'replay-v32-{model}-raw'),('v49_to_extended',run/f'replay-extended-{model}-raw')]:
            other,othergenes,_,_=comparison.load_selected(otherpath,fixed)
            ga=list(map(stable_id,basegenes));gb=list(map(stable_id,othergenes));union=sorted(set(ga)|set(gb))
            ba=align(base,ga,union);ot=align(other,gb,union)
            a,b=(ot,ba) if condition=='v32_to_v49' else (ba,ot)
            effects[model+'/'+condition]=comparison.compare(a,b)
            effects[model+'/'+condition]['gene_matching']='Ensembl accession without version, preserving _PAR_Y suffix'
            effects[model+'/'+condition]['scope']='same 6460 Gene-called nuclei; replay of fixed annotation-free alignments'
    (out/'annotation-effects.json').write_text(json.dumps(effects,indent=2)+'\n')
    print(json.dumps({'calling':comparisons,'annotation_effects':{k:{name:v[name] for name in ['a_collapsed_umis','b_collapsed_umis','positive_count_difference_umis','negative_count_difference_umis','half_l1_over_a_umis']} for k,v in effects.items()}},indent=2))

if __name__=='__main__':main()
