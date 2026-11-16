"""Transparent two-strand OOPS EM motif discovery and composition-aware scoring.

Backgrounds are strand-symmetrized zero-order models. ``adaptive`` shrinks each
sequence's counts toward a pooled background; it is a project hypothesis, not a
claim to invent EM, adaptive backgrounds, or motif enrichment. PWM shape: W x 4
in A,C,G,T order. Coordinates: zero-based start on the original forward sequence.
"""
from collections import Counter
from time import perf_counter
import numpy as np
from scipy.special import logsumexp
from .io import normalize_sequence

BASES = 'ACGT'
COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')
METHODS = ('uniform', 'global', 'adaptive', 'kmer')


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def _counts(sequence):
    c = np.array([sequence.count(base) for base in BASES], dtype=float)
    return (c + c[::-1]) / 2


def estimate_background(sequences, method='adaptive', shrinkage=40.0, global_background=None):
    """Return (per-sequence 4-vectors, pooled 4-vector), with Jeffreys smoothing."""
    if method not in METHODS:
        raise ValueError(f'Unknown method: {method}')
    if not np.isfinite(shrinkage) or shrinkage < 0:
        raise ValueError('shrinkage must be finite and nonnegative')
    counts = np.array([_counts(s) for s in sequences])
    pooled = counts.sum(axis=0) + 0.5
    pooled /= pooled.sum()
    if global_background is not None:
        pooled = np.array(global_background, dtype=float)
        if pooled.shape != (4,) or not np.isfinite(pooled).all() or np.any(pooled <= 0):
            raise ValueError('global_background must have four positive finite entries')
        pooled = (pooled + pooled[::-1]) / 2
        pooled /= pooled.sum()
    if method == 'uniform':
        backgrounds = np.full((len(sequences), 4), 0.25)
        pooled = np.full(4, 0.25)
    elif method in ('global', 'kmer'):
        backgrounds = np.tile(pooled, (len(sequences), 1))
    else:
        backgrounds = counts + shrinkage * pooled + 0.5
        backgrounds /= backgrounds.sum(axis=1, keepdims=True)
    return backgrounds, pooled


def _prepare(sequences, width):
    if isinstance(sequences, str):
        raise ValueError('sequences must be a collection of DNA strings')
    sequences = [normalize_sequence(s) for s in sequences]
    if not sequences:
        raise ValueError('At least one sequence is required')
    if not isinstance(width, (int, np.integer)) or width < 2 or width > 50:
        raise ValueError('width must be an integer between 2 and 50')
    windows, starts = [], []
    for i, sequence in enumerate(sequences):
        encoded = np.array([BASES.find(c) for c in sequence], dtype=np.int8)
        if len(encoded) < width:
            raise ValueError(f'Sequence {i} is shorter than motif width {width}')
        all_windows = np.lib.stride_tricks.sliding_window_view(encoded, width)
        good = np.all(all_windows >= 0, axis=1)
        forward = all_windows[good]
        if not len(forward):
            raise ValueError(f'Sequence {i} has no unambiguous window of width {width}')
        windows.append(np.concatenate([forward, 3 - forward[:, ::-1]], axis=0))
        starts.append(np.flatnonzero(good))
    return sequences, windows, starts


def _pwm(pwm):
    p = np.array(pwm, dtype=float)
    if p.ndim != 2 or p.shape[1] != 4 or not np.isfinite(p).all() or np.any(p < 0) or np.any(p.sum(axis=1) <= 0):
        raise ValueError('PWM must be W x 4, finite, nonnegative, with positive row sums')
    p = np.maximum(p, 1e-12)
    return p / p.sum(axis=1, keepdims=True)


def _scores(windows, pwm, backgrounds):
    columns = np.arange(len(pwm))[None, :]
    lp = np.log(pwm)
    return [(lp[columns, w] - np.log(b)[w]).sum(axis=1) for w, b in zip(windows, backgrounds)]


def _sites(sequences, windows, starts, scores):
    answer = []
    for i, (w, pos, score) in enumerate(zip(windows, starts, scores)):
        j = int(np.argmax(score))
        start = int(pos[j % len(pos)])
        answer.append({'sequence_index': i, 'start': start, 'end': start + w.shape[1],
                       'strand': '+' if j < len(pos) else '-', 'score': float(score[j]),
                       'sequence': ''.join(BASES[x] for x in w[j])})
    return answer


def score_sequences(sequences, pwm, method='adaptive', shrinkage=40.0, global_background=None):
    """Return each sequence's maximum natural-log likelihood-ratio site."""
    pwm = _pwm(pwm)
    sequences, windows, starts = _prepare(sequences, len(pwm))
    backgrounds, _ = estimate_background(sequences, method, shrinkage, global_background)
    return _sites(sequences, windows, starts, _scores(windows, pwm, backgrounds))


def _seed_words(windows, backgrounds, n_starts, strategy, rng):
    if strategy == 'random':
        chosen = []
        for _ in range(n_starts):
            w = windows[int(rng.integers(len(windows)))]
            chosen.append(tuple(w[int(rng.integers(len(w)))]))
        return chosen
    counts = Counter()
    for w in windows:
        # A sequence contributes once, so repeats cannot dominate seed ranking.
        words = {min(tuple(row), tuple(3 - row[::-1])) for row in w[:len(w)//2]}
        counts.update(words)
    # Batch background expectations to avoid a Python loop over words x sequences.
    ranked = []
    words = list(counts)
    lengths = np.array([len(w) / 2 for w in windows])
    log_backgrounds = np.log(backgrounds).T
    for offset in range(0, len(words), 1024):
        chunk = words[offset:offset + 1024]
        encoded = np.array(chunk, dtype=np.int8)
        base_counts = np.array([(encoded == base).sum(axis=1) for base in range(4)]).T
        multiplicity = np.where(np.all(encoded == 3 - encoded[:, ::-1], axis=1), 1, 2)
        rate = np.exp(base_counts @ log_backgrounds) * multiplicity[:, None] * lengths[None, :]
        expected = (-np.expm1(-rate)).sum(axis=1)
        for word, e in zip(chunk, expected):
            observed = counts[word]
            ranked.append(((observed - e) / np.sqrt(e + 1), observed, word))
    ranked.sort(reverse=True)
    return [item[2] for item in ranked[:n_starts]]


def discover(sequences, width, method='adaptive', seed=0, n_starts=8,
             max_iter=60, shrinkage=40.0, seed_strategy='enriched', tol=1e-6):
    """Discover one OOPS motif. ``kmer`` reports the best seed without EM.

    All EM methods use identical seed strategy, pseudocount and stopping rule.
    ``seed_strategy='random'`` is the candidate-screening ablation. Objective is
    sum_i log(mean_j exp(LLR_ij)); the background-only likelihood is omitted.
    """
    start_time = perf_counter()
    if method not in METHODS or seed_strategy not in ('enriched', 'random'):
        raise ValueError('Unknown method or seed_strategy')
    if not isinstance(n_starts, (int, np.integer)) or n_starts < 1 or not isinstance(max_iter, (int, np.integer)) or max_iter < 1:
        raise ValueError('n_starts and max_iter must be positive integers')
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError('tol must be positive and finite')
    sequences, windows, starts = _prepare(sequences, width)
    backgrounds, pooled = estimate_background(sequences, method, shrinkage)
    seeds = _seed_words(windows, backgrounds, n_starts, seed_strategy, np.random.default_rng(seed))
    best = None
    columns = np.arange(width)[None, :]
    for word in seeds[:1] if method == 'kmer' else seeds:
        pwm = np.full((width, 4), 0.05)
        pwm[np.arange(width), word] = 0.85
        iterations = 0
        if method != 'kmer':
            for iteration in range(max_iter):
                scores = _scores(windows, pwm, backgrounds)
                counts = np.full((width, 4), 0.25)
                for w, s in zip(windows, scores):
                    responsibilities = np.exp(s - logsumexp(s))
                    for base in range(4):
                        counts[:, base] += (responsibilities[:, None] * (w == base)).sum(axis=0)
                updated = counts / counts.sum(axis=1, keepdims=True)
                iterations = iteration + 1
                change = np.max(np.abs(updated - pwm))
                pwm = updated
                if change < tol:
                    break
        scores = _scores(windows, pwm, backgrounds)
        objective = float(sum(logsumexp(s) - np.log(len(s)) for s in scores))
        if best is None or objective > best['log_likelihood']:
            best = {'pwm': pwm.tolist(), 'sites': _sites(sequences, windows, starts, scores),
                    'log_likelihood': objective, 'iterations': iterations,
                    'consensus': ''.join(BASES[x] for x in pwm.argmax(axis=1))}
    best.update({'method': method, 'width': width, 'seed': int(seed), 'n_starts': len(seeds),
                 'seed_strategy': seed_strategy, 'shrinkage': float(shrinkage),
                 'background': pooled.tolist(), 'global_background': pooled.tolist(), 'sequence_backgrounds': backgrounds.tolist(),
                 'runtime_seconds': perf_counter() - start_time,
                 'model': 'two-strand OOPS; zero-order strand-symmetric background',
                 'coordinate_system': 'zero-based, half-open, original forward sequence'})
    return best
