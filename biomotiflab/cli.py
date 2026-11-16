import argparse
from pathlib import Path
from .io import read_fasta
from .core import discover, reverse_complement
from .stats import calibrate
from .report import json_text, markdown_report, logo_figure


def main():
    p = argparse.ArgumentParser(description='Discover an OOPS DNA motif; optionally test independent sequences.')
    p.add_argument('fasta')
    p.add_argument('--test-fasta')
    p.add_argument('--width', type=int, default=10)
    p.add_argument('--method', choices=['uniform','global','adaptive','kmer'], default='adaptive')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--starts', type=int, default=8)
    p.add_argument('--iterations', type=int, default=60)
    p.add_argument('--shrinkage', type=float, default=40.)
    p.add_argument('--shuffles', type=int, default=199)
    p.add_argument('--out', default='analysis')
    a = p.parse_args()
    try:
        seqs = read_fasta(a.fasta)
        result = discover(seqs,a.width,method=a.method,seed=a.seed,n_starts=a.starts,max_iter=a.iterations,shrinkage=a.shrinkage)
        result['n_train'] = len(seqs)
        result['configuration'] = vars(a)
        if a.test_fasta:
            test = read_fasta(a.test_fasta)
            if {min(s, reverse_complement(s)) for s in seqs} & {min(s, reverse_complement(s)) for s in test}:
                raise ValueError('Training and test sequences overlap (including reverse complements); independent test sequences required.')
            result['heldout'] = calibrate(test,result['pwm'],method=a.method if a.method != 'kmer' else 'global',n_shuffles=a.shuffles,seed=a.seed,shrinkage=a.shrinkage,global_background=result['global_background'])
        out = Path(a.out); out.mkdir(parents=True,exist_ok=True)
        (out/'results.json').write_text(json_text(result),encoding='utf-8')
        (out/'report.md').write_text(markdown_report(result),encoding='utf-8')
        logo_figure(result['pwm']).savefig(out/'logo.png',dpi=180)
        print(f"{result['consensus']} → {out.resolve()}")
    except (ValueError, OSError) as exc:
        p.error(str(exc))

if __name__ == '__main__':
    main()
