"""Run real, reproducible experiments. No cached or invented measurements."""
import argparse, json, platform, time, resource, sys
from pathlib import Path
import numpy as np
import pandas as pd
from biomotiflab.core import discover, score_sequences
from biomotiflab.stats import calibrate
from biomotiflab.synthetic import generate_dataset
from biomotiflab.io import read_fasta, write_fasta

ROOT = Path(__file__).resolve().parents[1]


def localization(predictions,truth,width):
    overlaps, exact = [], []
    for p,t in zip(predictions,truth):
        if t is None: continue
        inter = max(0,min(p['start']+width,t['end'])-max(p['start'],t['start']))
        overlaps.append(2*inter/(width+t['end']-t['start']))
        exact.append(p['start']==t['start'])
    return dict(nucleotide_f1=float(np.mean(overlaps)),exact_start=float(np.mean(exact)))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default=str(ROOT/'results'))
    parser.add_argument('--seeds',default='11,22,33,44,55');a=parser.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    seeds=[int(x) for x in a.seeds.split(',')]
    rows=[]; details=[]
    config=dict(seeds=seeds,n_train=30,n_test=30,length=80,width=10,starts=6,iterations=50,
        mutation=[0.,.15,.3],heterogeneity=[0.,1.],methods=['kmer','uniform','global','adaptive'],
        significance_shuffles=199,split='train seed s; test seed s+10000; independent generation',
        resource_metric='process high-water RSS, macOS bytes / Linux KiB; not incremental method allocation')
    (out/'config.json').write_text(json.dumps(config,indent=2))
    def run(group,seed,noise=.15,het=1.,method='adaptive',tau=40.,strategy='enriched',n=30,length=80,fit_width=10):
        train=generate_dataset(n=n,length=length,width=10,mutation=noise,heterogeneity=het,seed=seed)
        test=generate_dataset(n=30,length=length,width=10,mutation=noise,heterogeneity=het,seed=seed+10000)
        t=time.perf_counter()
        r=discover(train['sequences'],fit_width,method=method,seed=seed,n_starts=6,max_iter=50,shrinkage=tau,seed_strategy=strategy)
        fit_time=time.perf_counter()-t
        sites=score_sequences(test['sequences'],r['pwm'],method='global' if method=='kmer' else method,shrinkage=tau,global_background=r['global_background'])
        row=dict(group=group,seed=seed,mutation=noise,heterogeneity=het,method=method,shrinkage=tau,seed_strategy=strategy,n_train=n,length=length,width=fit_width,
                 fit_seconds=fit_time,**localization(sites,test['sites'],fit_width))
        rows.append(row); details.append(dict(config=row,model=r,test_predictions=sites,test_truth=test['sites']))
        pd.DataFrame(rows).to_csv(out/'raw.csv',index=False)
        print(len(rows),group,seed,method,round(row['nucleotide_f1'],3),round(fit_time,3),flush=True)
    for het in [0.,1.]:
        for noise in [0.,.15,.3]:
            for seed in seeds:
                for method in config['methods']: run('main',seed,noise,het,method)
    for seed in seeds:
        run('ablation',seed,tau=0.)
        run('ablation',seed,strategy='random')
        for tau in [5.,200.,1000.]: run('tau',seed,tau=tau)
        for width in [8,12]: run('width',seed,fit_width=width)
        for n in [15,60]: run('scale_n',seed,n=n)
        for length in [160,320]: run('scale_length',seed,length=length)
    # Independent train/test fixed-model classification, with true null negatives.
    sig=[]
    for seed in seeds:
        train=generate_dataset(n=30,length=80,mutation=.15,heterogeneity=1.,seed=seed)
        test=generate_dataset(n=40,length=80,mutation=.15,heterogeneity=1.,seed=seed+20000,positive_fraction=.5)
        for method in ['global','adaptive']:
            r=discover(train['sequences'],10,method=method,seed=seed,n_starts=6,max_iter=50)
            records=calibrate(test['sequences'],r['pwm'],method=method,n_shuffles=199,seed=seed+30000,global_background=r['global_background'])
            for rec,truth in zip(records,test['sites']):sig.append(dict(seed=seed,method=method,positive=truth is not None,**rec))
    pd.DataFrame(sig).to_csv(out/'significance.csv',index=False)
    # CRP: no validated location labels; six held-out sequences per split.
    real=read_fasta(ROOT/'data/real/crp.fasta'); real_results=[]
    for seed in seeds:
        order=np.random.default_rng(seed).permutation(len(real)); train=[real[i] for i in order[:-6]]; test=[real[i] for i in order[-6:]]
        for method in ['kmer','uniform','global','adaptive']:
            r=discover(train,20,method=method,seed=seed,n_starts=8,max_iter=60)
            held=calibrate(test,r['pwm'],method='global' if method=='kmer' else method,n_shuffles=199,seed=seed,global_background=r['global_background'])
            real_results.append(dict(seed=seed,method=method,train_indices=order[:-6].tolist(),test_indices=order[-6:].tolist(),model=r,heldout=held))
    (out/'real.json').write_text(json.dumps(real_results,indent=2))
    (out/'details.json').write_text(json.dumps(details,indent=2))
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (out/'environment.json').write_text(json.dumps(dict(python=sys.version,platform=platform.platform(),processor=platform.processor(),numpy=np.__version__,peak_rss_mb=rss/(1024**2 if sys.platform=='darwin' else 1024),measure='cumulative whole-experiment process peak RSS; not per-algorithm memory'),indent=2))
    example=generate_dataset(n=24,length=80,mutation=.05,seed=101)
    write_fasta(example['sequences'],ROOT/'data/example/train.fasta')
    (ROOT/'data/example/truth.json').write_text(json.dumps(example,indent=2))
    write_fasta(generate_dataset(n=12,length=80,mutation=.05,seed=102)['sequences'],ROOT/'data/example/test.fasta')

if __name__=='__main__': main()
