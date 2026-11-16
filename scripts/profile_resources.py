"""Measure each algorithm in a fresh interpreter; includes library baseline RSS."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--worker', choices=['kmer', 'uniform', 'global', 'adaptive'])
    parser.add_argument('--out', default=str(ROOT / 'results/resources.csv'))
    parser.add_argument('--n', type=int, default=30)
    parser.add_argument('--length', type=int, default=80)
    parser.add_argument('--seed', type=int, default=11)
    args = parser.parse_args()
    if args.worker:
        import resource
        import time
        from biomotiflab.core import discover
        from biomotiflab.synthetic import generate_dataset
        data = generate_dataset(n=args.n, length=args.length, mutation=.15, heterogeneity=1., seed=args.seed)
        divisor = 1024 ** 2 if sys.platform == 'darwin' else 1024
        baseline = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / divisor
        t = time.perf_counter()
        model = discover(data['sequences'], 10, method=args.worker, seed=args.seed, n_starts=6, max_iter=50)
        seconds = time.perf_counter() - t
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / divisor
        print(json.dumps(dict(method=args.worker, n=args.n, length=args.length, seed=args.seed,
                              width=10, mutation=.15, heterogeneity=1., starts=6, max_iterations=50,
                              fit_seconds=seconds, baseline_rss_mb=baseline, peak_rss_mb=peak,
                              peak_above_baseline_mb=max(0., peak-baseline),
                              metric='fresh subprocess ru_maxrss; includes interpreter and libraries; high-water difference is not total allocation',
                              python=sys.version.split()[0], platform=sys.platform)))
        return
    rows = []
    for method in ['kmer', 'uniform', 'global', 'adaptive']:
        run = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker', method,
                              '--n', str(args.n), '--length', str(args.length), '--seed', str(args.seed)],
                             cwd=ROOT, text=True, capture_output=True, check=True)
        rows.append(json.loads(run.stdout))
        print(method, rows[-1]['peak_rss_mb'], 'MiB', flush=True)
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
