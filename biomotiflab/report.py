"""Serializable reports and honest sequence logos."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logomaker


def logo_figure(pwm):
    p = np.asarray(pwm, dtype=float)
    information = np.maximum(0, 2 + (p * np.log2(np.maximum(p, 1e-12))).sum(axis=1))
    fig, ax = plt.subplots(figsize=(max(6, len(p)*.4), 2.8))
    logomaker.Logo(pd.DataFrame(p * information[:, None], columns=list('ACGT')), ax=ax)
    ax.set(xlabel='Motif position (0-based)', ylabel='Bits', ylim=(0, 2))
    fig.tight_layout()
    return fig


def json_text(result):
    return json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)


def markdown_report(result):
    lines = ['# BioMotifLab analysis', '', f"Method: {result['method']}",
             f"Consensus: {result['consensus']}",
             f"Training sequences: {result.get('n_train', 'see configuration')}",
             '', 'Coordinates: 0-based starts on original forward input; strand refers to motif orientation.',
             'OOPS assumes exactly one motif per sequence. Scores are log likelihood ratios, not binding probabilities.',
             'Permutation p-values, when supplied, concern a fixed motif on held-out sequences under an IID composition-preserving shuffle null.',
             'They are not experimental proof of TF binding and do not model dinucleotide dependence.', '',
             '## Configuration and results', '```json', json_text(result), '```']
    return '\n'.join(lines)
