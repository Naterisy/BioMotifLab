"""Held-out, fixed-motif conditional randomization tests and BH adjustment.

The null randomly permutes ACGT bases within each sequence (N locations fixed),
conditioning on mononucleotide counts. It does NOT preserve dinucleotides. Use
new held-out sequences: testing the motif-training set would be selective and
anti-conservative. A sequence-level test uses its maximum over both strands and
all valid windows, automatically including that within-sequence search.
"""
import numpy as np
from .core import score_sequences, estimate_background
from .io import normalize_sequence


def benjamini_hochberg(p_values):
    """Monotone BH adjusted p-values; FDR control requires suitable dependence."""
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError('p_values must be a finite one-dimensional vector in [0, 1]')
    if not len(p):
        return []
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    result = np.empty_like(q)
    result[order] = np.minimum(q, 1)
    return result.tolist()


def calibrate(sequences, pwm, method='adaptive', n_shuffles=199, seed=0,
              shrinkage=40.0, global_background=None):
    """Fixed-motif empirical p=(1+null>=observed)/(B+1); BH across sequences.

    ``global_background`` should normally be the discovery training background.
    If omitted it is estimated from all supplied sequences and frozen; because
    each permutation preserves counts this remains a conditional randomization
    test. q-values do not correct across repeated choices of PWM/width/method.
    """
    if not isinstance(n_shuffles, (int, np.integer)) or n_shuffles < 1:
        raise ValueError('n_shuffles must be a positive integer')
    sequences = [normalize_sequence(s) for s in sequences]
    if global_background is None:
        _, global_background = estimate_background(sequences, method, shrinkage)
    observed = score_sequences(sequences, pwm, method, shrinkage, global_background)
    rng = np.random.default_rng(seed)
    p_values = []
    for sequence, site in zip(sequences, observed):
        chars = np.array(list(sequence))
        mask = chars != 'N'
        letters = chars[mask].copy()
        shuffled = []
        for _ in range(n_shuffles):
            chars[mask] = rng.permutation(letters)
            shuffled.append(''.join(chars))
        null = score_sequences(shuffled, pwm, method, shrinkage, global_background)
        exceedances = sum(s['score'] >= site['score'] - 1e-12 for s in null)
        p_values.append((exceedances + 1) / (n_shuffles + 1))
        site.update({'p_value': p_values[-1], 'null_exceedances': exceedances,
                     'n_shuffles': int(n_shuffles), 'null_model': 'within-sequence mononucleotide permutation',
                     'inference_scope': 'fixed PWM on held-out sequences; sequence maximum over both strands'})
    for site, q in zip(observed, benjamini_hochberg(p_values)):
        site['q_value'] = q
    return observed


bh_adjust = benjamini_hochberg
