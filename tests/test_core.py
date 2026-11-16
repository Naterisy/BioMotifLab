import numpy as np
import pytest
from biomotiflab.io import read_fasta, write_fasta, normalize_sequence
from biomotiflab.core import discover, score_sequences, estimate_background, reverse_complement
from biomotiflab.stats import benjamini_hochberg, calibrate
from biomotiflab.synthetic import generate_dataset


def test_fasta_round_trip_and_ambiguity(tmp_path):
    seqs = read_fasta('>a\nacgtRy\n>b\nTTGG\n', is_text=True)
    assert seqs == ['ACGTNN', 'TTGG']
    path = tmp_path / 'x.fa'
    write_fasta(seqs, path)
    assert read_fasta(path) == seqs
    with pytest.raises(ValueError):
        read_fasta('ACGT', is_text=True)
    with pytest.raises(ValueError):
        read_fasta('>a\n>b\nACGT', is_text=True)
    with pytest.raises(ValueError):
        normalize_sequence('ACG!')


def test_reverse_strand_coordinates_independent_score():
    motif = 'ACGACCTG'
    pwm = np.full((len(motif), 4), .01)
    for j, b in enumerate(motif):
        pwm[j, 'ACGT'.index(b)] = .97
    seq = 'NNNN' + reverse_complement(motif) + 'NNN'
    site = score_sequences([seq], pwm, method='uniform')[0]
    assert (site['start'], site['end'], site['strand'], site['sequence']) == (4, 12, '-', motif)
    assert site['score'] == pytest.approx(8 * np.log(.97 / .25))


def test_background_shrinkage_and_strand_symmetry():
    seqs = ['A' * 80 + 'G' * 20, 'C' * 80 + 'T' * 20]
    b, global_b = estimate_background(seqs, 'adaptive', 40)
    expected_count = np.array([40, 10, 10, 40], float)
    expected = (expected_count + 40 * global_b + .5) / 142
    np.testing.assert_allclose(b[0], expected)
    np.testing.assert_allclose(b.sum(axis=1), 1)
    np.testing.assert_allclose(b[:, ::-1], b)
    near_global, _ = estimate_background(seqs, 'adaptive', 1e10)
    np.testing.assert_allclose(near_global[0], global_b, atol=1e-7)


def test_bh_known_example():
    assert benjamini_hochberg([.01, .04, .03, .002]) == pytest.approx([.02, .04, .04, .008])
    assert benjamini_hochberg([]) == []
    with pytest.raises(ValueError):
        benjamini_hochberg([1.1])


def test_strong_motif_recovery_and_normalization():
    d = generate_dataset(n=35, length=70, width=10, mutation=0, seed=4, motif='ACGACCTGTA')
    result = discover(d['sequences'], 10, method='global', n_starts=4)
    np.testing.assert_allclose(np.sum(result['pwm'], axis=1), 1)
    accuracy = np.mean([a['start'] == b['start'] for a, b in zip(result['sites'], d['sites'])])
    assert accuracy >= .9
    assert result['consensus'] in [d['motif'], reverse_complement(d['motif'])]


def test_random_seed_reproducibility():
    d = generate_dataset(n=8, length=30, width=6, seed=1)
    a = discover(d['sequences'], 6, seed_strategy='random', seed=23, n_starts=2, max_iter=8)
    b = discover(d['sequences'], 6, seed_strategy='random', seed=23, n_starts=2, max_iter=8)
    assert a['pwm'] == b['pwm']
    assert a['sites'] == b['sites']


def test_permutation_inference_bounds_and_reproducibility():
    d = generate_dataset(n=4, length=30, width=6, mutation=.05, seed=12)
    a = calibrate(d['sequences'], d['pwm'], n_shuffles=19, seed=20)
    b = calibrate(d['sequences'], d['pwm'], n_shuffles=19, seed=20)
    assert a == b
    assert all(.05 <= site['p_value'] <= site['q_value'] <= 1 for site in a)
    # An uninformative PWM produces identical scores under a uniform background.
    flat = calibrate(['ACGTACGTACGT'], np.ones((4, 4)) / 4, method='uniform', n_shuffles=9)
    assert flat[0]['p_value'] == flat[0]['q_value'] == 1


def test_validation():
    for sequences, width in [([], 6), (['AAAA'], 6), (['NNNNNN'], 4), (['ACGT'], 0)]:
        with pytest.raises(ValueError):
            discover(sequences, width)
    with pytest.raises(ValueError):
        discover(['ACGTACGT'], 4, method='made_up')
    with pytest.raises(ValueError):
        score_sequences(['ACGTACGT'], [[1, -1, 1, 1]] * 4)
    with pytest.raises(ValueError):
        discover(['ACGTACGT'], 4, shrinkage=-1)


def test_synthetic_truth_matches_sequence_both_strands():
    d = generate_dataset(n=30, length=50, width=8, seed=9, mutation=.2, positive_fraction=.5)
    assert sum(site is not None for site in d['sites']) == 15
    for seq, site in zip(d['sequences'], d['sites']):
        if site:
            actual = seq[site['start']:site['end']]
            assert actual == (site['sequence'] if site['strand'] == '+' else reverse_complement(site['sequence']))
