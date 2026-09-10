#!/usr/bin/env python3
"""Render a standalone scientific summary from completed comparison artifacts."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser();ap.add_argument('--comparison',type=Path,required=True);ap.add_argument('--calling',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
s=json.loads((a.comparison/'summary.json').read_text());c=json.loads(a.calling.read_text())
p=np.genfromtxt(a.comparison/'per-nucleus.tsv',names=True,delimiter='\t',dtype=None,encoding='utf-8')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(2,2,figsize=(10,8),layout='constrained')
for ax,unit,title in [(axs[0,0],'umis','A  Collapsed UMIs per fixed nucleus'),(axs[0,1],'genes','B  Detected genes per fixed nucleus')]:
 x=p['replay_Gene_'+unit];y=p['replay_GeneFull_'+unit]
 ax.hexbin(x,y,gridsize=45,mincnt=1,xscale='log',yscale='log',cmap='Blues')
 lo=min(x.min(),y.min());hi=max(x.max(),y.max());ax.plot([lo,hi],[lo,hi],color='0.4',ls='--',lw=1)
 ax.set(xlabel='Gene',ylabel='GeneFull',title=title)
 ax.text(.04,.95,'Same 6,460 nuclei',transform=ax.transAxes,va='top')
ax=axs[1,0]
values=[100*s['comparisons']['fixed_Gene_called/STAR_'+m+'_to_replay_'+m]['half_l1_over_a_umis'] for m in ['Gene','GeneFull']]
ax.bar(['Gene','GeneFull'],values,color=['#b7683a','#326b92'],width=.55)
for i,v in enumerate(values):ax.text(i,v+.015,f'{v:.3f}%',ha='center')
ax.set(ylabel='100 × half-L1 / matched STARsolo UMI total',ylim=(0,1),title='C  End-to-end matrix difference')
ax=axs[1,1];v=c['comparisons']['replay_Gene to replay_GeneFull']
ax.barh(['Gene only','Shared','GeneFull only'],[v['a_only'],v['intersection'],v['b_only']],color=['#b7683a','#71816d','#326b92'])
ax.set(xscale='log',xlabel='Number of called barcodes (log scale)',title='D  Nucleus calling on each replay matrix')
for i,n in enumerate([v['a_only'],v['intersection'],v['b_only']]):ax.text(n*1.08,i,f'{n:,}',va='center')
ax.set_xlim(1,15000)
fig.suptitle('GeneFull replay in brain nuclei',fontsize=16)
for ext in ['png','pdf']:fig.savefig(a.out/f'genefull-summary.{ext}',dpi=180)
