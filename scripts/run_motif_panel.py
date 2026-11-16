"""Post-hoc exploratory motif panel; all prelisted panel conditions are retained.

This extension was designed after the fixed-motif main study was inspected. It
is explicitly exploratory and is not a preregistered confirmation. The panel
includes mixed, GC-rich, AT-rich and approximately palindromic motifs; no motif
is removed on account of performance. Same protocol as main held-out study.
"""
import argparse
import json
from pathlib import Path
import time
import pandas as pd
from biomotiflab.core import discover, score_sequences
from biomotiflab.synthetic import generate_dataset
from run_experiments import localization

ROOT = Path(__file__).resolve().parents[1]
MOTIFS = ['ACGTCAGTCA', 'GCGCGGCCGC', 'ATATTAATAT', 'TGACGTCAAG']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default=str(ROOT / 'results/motif_panel.csv'))
    args = parser.parse_args()
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows, details = [], []
    for motif in MOTIFS:
        for seed in [101, 202, 303]:
            for het in [0., 1.]:
                train = generate_dataset(n=30, length=80, width=10, mutation=.3, heterogeneity=het, seed=seed, motif=motif)
                test = generate_dataset(n=30, length=80, width=10, mutation=.3, heterogeneity=het, seed=seed+10000, motif=motif)
                for method in ['global', 'adaptive']:
                    t = time.perf_counter()
                    model = discover(train['sequences'],10,method=method,seed=seed,n_starts=6,max_iter=50)
                    seconds = time.perf_counter()-t
                    predictions = score_sequences(test['sequences'],model['pwm'],method=method,global_background=model['global_background'])
                    row = dict(stage='posthoc_exploratory', motif=motif, seed=seed, heterogeneity=het,
                               mutation=.3, n_train=30, n_test=30, length=80, width=10,
                               method=method, shrinkage=40., starts=6, iterations=50,
                               fit_seconds=seconds, **localization(predictions,test['sites'],10))
                    rows.append(row)
                    details.append(dict(config=row,model=model,test_predictions=predictions,test_truth=test['sites']))
                    pd.DataFrame(rows).to_csv(path,index=False)
                    print(len(rows), motif, seed, het, method, round(row['nucleotide_f1'],3), flush=True)
    pd.DataFrame(rows).groupby(['motif', 'heterogeneity', 'method']).nucleotide_f1.agg(['mean', 'std']).to_csv(path.with_name(path.stem + '_summary.csv'))
    path.with_suffix('.json').write_text(json.dumps(dict(
        protocol='Post-hoc exploratory panel, designed after inspecting main single-motif results; retain all failures.',
        motifs=MOTIFS,seeds=[101,202,303],records=details),indent=2))


if __name__ == '__main__':
    main()
