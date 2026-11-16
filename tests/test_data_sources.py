"""Offline data-integrity checks; the optional download test needs explicit opt-in."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'real'
spec = importlib.util.spec_from_file_location('biomotiflab_download_data', ROOT / 'scripts' / 'download_data.py')
downloader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(downloader)


def test_bundled_jaspar_matrix_identity_and_integrity():
    """Matrix identity and count content, not mutable API metadata, are the pin."""
    manifest = json.loads((DATA / 'manifest.json').read_text())['datasets'][1]
    matrix = json.loads((DATA / manifest['file']).read_text())
    assert matrix['matrix_id'] == 'MA0139.1'
    assert matrix['name'] == 'CTCF'
    assert set(matrix['pfm']) == set('ACGT')
    assert {len(row) for row in matrix['pfm'].values()} == {19}
    assert all(value >= 0 for row in matrix['pfm'].values() for value in row)
    canonical = json.dumps(matrix['pfm'], sort_keys=True, separators=(',', ':')).encode()
    assert hashlib.sha256(canonical).hexdigest() == manifest['canonical_pfm_sha256']
    assert manifest['license'] == 'CC BY 4.0'


def test_local_crp_file_matches_official_pin():
    """Restricted source FASTA is optional in an offline source checkout."""
    path = DATA / 'crp.fasta'
    if not path.exists():
        pytest.skip('Run scripts/download_data.py to obtain the official CRP data')
    from biomotiflab.io import read_fasta
    assert hashlib.sha256(path.read_bytes()).hexdigest() == downloader.CRP_SHA256
    sequences = read_fasta(path)
    assert len(sequences) == 18
    assert {len(s) for s in sequences} == {105}


def test_changed_source_is_rejected_before_data_write(tmp_path, monkeypatch):
    """A changed upstream file must not silently replace the benchmark."""
    monkeypatch.setattr(downloader, 'download', lambda *_: b'>altered\nACGT\n')
    monkeypatch.setattr(sys, 'argv', ['download_data.py', '--output', str(tmp_path)])
    with pytest.raises(ValueError, match='CRP source changed'):
        downloader.main()
    assert not list(tmp_path.iterdir())


@pytest.mark.skipif(os.environ.get('BIOMOTIFLAB_ONLINE_TEST') != '1', reason='Opt in to network validation with BIOMOTIFLAB_ONLINE_TEST=1')
def test_online_official_sources(tmp_path, monkeypatch):
    args = ['download_data.py', '--output', str(tmp_path)]
    if os.environ.get('SSL_CERT_FILE'):
        args += ['--cafile', os.environ['SSL_CERT_FILE']]
    monkeypatch.setattr(sys, 'argv', args)
    downloader.main()
    assert hashlib.sha256((tmp_path / 'crp.fasta').read_bytes()).hexdigest() == downloader.CRP_SHA256
    assert (tmp_path / 'MEME-COPYRIGHT.html').stat().st_size > 100
    assert json.loads((tmp_path / 'MA0139.1.json').read_text())['matrix_id'] == 'MA0139.1'
