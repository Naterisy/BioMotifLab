"""Reproducible DNA generation with exact embedded-site ground truth."""
import numpy as np
from .core import BASES, reverse_complement
from .io import normalize_sequence


def generate_dataset(n=40, length=100, width=10, mutation=0.1,
                     heterogeneity=0.0, seed=0, motif=None, positive_fraction=1.0):
    """Generate IID backgrounds, optionally bimodal GC, and zero/one site each.

    ``heterogeneity`` in [0,1] selects per-sequence GC 0.5 +/- 0.4*h. ``mutation``
    is the probability a planted consensus base changes to one of the other 3.
    The default motif is fixed across dataset seeds (seed 1729), and is never
    derived from the sampled backgrounds. Truth has None for negative sequences.
    """
    if not isinstance(n, (int, np.integer)) or n < 1 or not isinstance(length, (int, np.integer)) or length < 2:
        raise ValueError('n and length must be positive integers (length >= 2)')
    if not isinstance(width, (int, np.integer)) or not 2 <= width <= min(length, 50):
        raise ValueError('width must be between 2 and min(length,50)')
    if any(not np.isfinite(x) or not 0 <= x <= 1 for x in (mutation, heterogeneity, positive_fraction)):
        raise ValueError('mutation, heterogeneity and positive_fraction must lie in [0,1]')
    rng = np.random.default_rng(seed)
    if motif is None:
        motif = ''.join(np.random.default_rng(1729).choice(list(BASES), width))
    else:
        motif = normalize_sequence(motif)
        if len(motif) != width or 'N' in motif:
            raise ValueError('motif must be unambiguous and match width')
    pwm = np.full((width, 4), mutation / 3)
    for j, base in enumerate(motif):
        pwm[j, BASES.index(base)] = 1 - mutation
    sequences, sites, gc_values = [], [], []
    n_positive = round(n * positive_fraction)
    positive_indices = set(rng.permutation(n)[:n_positive].tolist())
    for i in range(n):
        gc = 0.5 + rng.choice([-1, 1]) * 0.4 * heterogeneity
        background = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        seq = rng.choice(list(BASES), size=length, p=background)
        truth = None
        if i in positive_indices:
            planted = ''.join(rng.choice(list(BASES), p=pwm[j]) for j in range(width))
            start = int(rng.integers(length - width + 1))
            strand = str(rng.choice(['+', '-']))
            oriented = planted if strand == '+' else reverse_complement(planted)
            seq[start:start + width] = list(oriented)
            truth = {'sequence_index': i, 'start': start, 'end': start + width,
                     'strand': strand, 'sequence': planted}
        sequences.append(''.join(seq))
        sites.append(truth)
        gc_values.append(float(gc))
    return {'sequences': sequences, 'sites': sites, 'motif': motif, 'pwm': pwm.tolist(),
            'background_gc': gc_values,
            'config': {'n': n, 'length': length, 'width': width, 'mutation': mutation,
                       'heterogeneity': heterogeneity, 'seed': seed, 'motif': motif,
                       'positive_fraction': positive_fraction}}
