"""Every number and figure derives from saved experiment records."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from biomotiflab.report import logo_figure

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'results'; F=R/'figures';F.mkdir(exist_ok=True)
df=pd.read_csv(R/'raw.csv'); main=df[df.group=='main']
summary=main.groupby(['heterogeneity','mutation','method'])[['nucleotide_f1','exact_start','fit_seconds']].agg(['mean','std'])
summary.to_csv(R/'summary.csv')
fig,axes=plt.subplots(1,2,figsize=(10,3.7),sharey=True)
for ax,(h,g) in zip(axes,main.groupby('heterogeneity')):
    for method,mg in g.groupby('method'):
        v=mg.groupby('mutation').nucleotide_f1.agg(['mean','std'])
        ax.errorbar(v.index,v['mean'],yerr=v['std'],marker='o',capsize=3,label=method)
    ax.set(title=f'GC heterogeneity = {h:g}',xlabel='Per-base mutation probability',ylim=(-.05,1.05)); ax.grid(alpha=.2)
axes[0].set_ylabel('Held-out nucleotide F1');axes[1].legend();fig.tight_layout();fig.savefig(F/'noise.png',dpi=190);plt.close(fig)
base=main[(main.heterogeneity==1)&(main.mutation==.15)&(main.method.isin(['global','adaptive']))].copy()
base['variant']=base.method
ab=df[df.group=='ablation'].copy();ab['variant']=np.where(ab.seed_strategy=='random','adaptive/random seeds','adaptive/no shrinkage')
ab=pd.concat([base,ab]); ab.groupby('variant')[['nucleotide_f1','fit_seconds']].agg(['mean','std']).to_csv(R/'ablation_summary.csv')
fig,ax=plt.subplots(figsize=(8,3.5));v=ab.groupby('variant').nucleotide_f1.agg(['mean','std']);ax.barh(v.index,v['mean'],xerr=v['std'],capsize=3);ax.set(xlabel='Held-out nucleotide F1',xlim=(0,1.1));fig.tight_layout();fig.savefig(F/'ablation.png',dpi=190);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(11,3.5))
for ax,(group,x) in zip(axes,[('tau','shrinkage'),('width','width'),('scale_length','length')]):
    subset=pd.concat([df[df.group==group],base[base.method=='adaptive']]);v=subset.groupby(x).nucleotide_f1.agg(['mean','std']);ax.errorbar(v.index,v['mean'],yerr=v['std'],marker='o',capsize=3);ax.set(xlabel=x,ylabel='Held-out nucleotide F1',ylim=(0,1.1));ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(F/'sensitivity.png',dpi=190);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(8,3.5))
for ax,(g,x) in zip(axes,[('scale_n','n_train'),('scale_length','length')]):
    sub=pd.concat([df[df.group==g],base[base.method=='adaptive']]);v=sub.groupby(x).fit_seconds.agg(['mean','std']);ax.errorbar(v.index,v['mean'],yerr=v['std'],marker='o',capsize=3);ax.set(xlabel=x,ylabel='Fit time (seconds)');ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(F/'runtime.png',dpi=190);plt.close(fig)
sig=pd.read_csv(R/'significance.csv');sr=[]
for (seed,method),g in sig.groupby(['seed','method']):
    for threshold in ['p_value','q_value']:
        pred=g[threshold]<=.05; truth=g.positive;tp=int((pred&truth).sum());fp=int((pred&~truth).sum());fn=int((~pred&truth).sum());tn=int((~pred&~truth).sum())
        sr.append(dict(seed=int(seed),method=method,threshold=threshold,precision=tp/max(1,tp+fp),recall=tp/max(1,tp+fn),f1=2*tp/max(1,2*tp+fp+fn),false_positive_rate=fp/max(1,fp+tn),tp=tp,fp=fp,fn=fn,tn=tn))
pd.DataFrame(sr).to_csv(R/'classification.csv',index=False)
pd.DataFrame(sr).groupby(['method','threshold'])[['precision','recall','f1','false_positive_rate']].agg(['mean','std']).to_csv(R/'classification_summary.csv')
real=json.loads((R/'real.json').read_text()); rr=[]
for item in real:
    held=item['heldout']; rr.append(dict(seed=item['seed'],method=item['method'],consensus=item['model']['consensus'],fit_seconds=item['model']['runtime_seconds'],p_le_05=sum(x['p_value']<=.05 for x in held),q_le_05=sum(x['q_value']<=.05 for x in held),n_test=len(held)))
pd.DataFrame(rr).to_csv(R/'real_summary.csv',index=False)
for item in real:
    if item['seed']==11 and item['method']=='adaptive':logo_figure(item['model']['pwm']).savefig(F/'crp_logo.png',dpi=190)
details=json.loads((R/'details.json').read_text()); candidates=[x for x in details if x['config']['group']=='main' and x['config']['method']=='adaptive']; cases={'success':max(candidates,key=lambda x:x['config']['nucleotide_f1']),'failure':min(candidates,key=lambda x:x['config']['nucleotide_f1'])}
(R/'cases.json').write_text(json.dumps(cases,indent=2))
for label,item in cases.items():logo_figure(item['model']['pwm']).savefig(F/f'{label}_logo.png',dpi=190)
print(summary.to_string());print('\nAblation\n',ab.groupby('variant').nucleotide_f1.agg(['mean','std']).to_string());print('\nReal\n',pd.DataFrame(rr).groupby('method')[['p_le_05','q_le_05']].agg(['mean','std']).to_string())
