"""Exercise the CLI and Streamlit's actual analysis callbacks.

Streamlit AppTest has no native file-uploader setter. Only that widget's return
is supplied as BytesIO; FASTA decoding, validation, computation, plotting,
comparison and export creation execute the production application code.
"""
from io import BytesIO
import json
from pathlib import Path
import subprocess
import sys
import pytest
from biomotiflab.io import write_fasta
from biomotiflab.synthetic import generate_dataset

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def sample_files(tmp_path):
    train = generate_dataset(n=12, length=40, width=10, mutation=.02, seed=63)
    test = generate_dataset(n=4, length=40, width=10, mutation=.02, seed=65)
    a, b = tmp_path / 'train.fa', tmp_path / 'test.fa'
    write_fasta(train['sequences'], a)
    write_fasta(test['sequences'], b)
    return a, b


def test_cli_outputs_and_rejects_overlapping_test(sample_files, tmp_path):
    train, test = sample_files
    out = tmp_path / 'result'
    args = [sys.executable, '-m', 'biomotiflab.cli', str(train), '--test-fasta', str(test),
            '--width', '10', '--starts', '2', '--iterations', '15', '--shuffles', '9', '--out', str(out)]
    run = subprocess.run(args, cwd=ROOT, text=True, capture_output=True, check=False)
    assert run.returncode == 0, run.stderr
    data = json.loads((out / 'results.json').read_text())
    assert len(data['sites']) == 12
    assert len(data['heldout']) == 4
    assert data['configuration']['shuffles'] == 9
    assert (out / 'logo.png').stat().st_size > 1000
    assert 'fixed motif' in (out / 'report.md').read_text()
    args[args.index('--test-fasta') + 1] = str(train)
    rejected = subprocess.run(args, cwd=ROOT, text=True, capture_output=True, check=False)
    assert rejected.returncode == 2
    assert 'overlap' in rejected.stderr


def test_streamlit_real_compute_and_exports(monkeypatch, sample_files):
    import streamlit as st
    from streamlit.testing.v1 import AppTest
    train, test = sample_files
    def upload(label, **kwargs):
        return BytesIO((train if label.startswith('上传发现集') else test).read_bytes())
    monkeypatch.setattr(st, 'file_uploader', upload)
    app = AppTest.from_file(str(ROOT / 'app.py'), default_timeout=30).run()
    assert not app.exception
    app.sidebar.multiselect[0].set_value(['uniform', 'adaptive'])
    app.sidebar.slider[1].set_value(2)
    app.sidebar.slider[2].set_value(10)
    app.sidebar.select_slider[0].set_value(49)
    app.button[0].click().run()
    assert not app.exception
    assert not app.error
    results = app.session_state['results']
    assert [r['method'] for r in results] == ['uniform', 'adaptive']
    assert all(len(r['sites']) == 12 and len(r['heldout']) == 4 for r in results)
    assert all(r['configuration']['n_starts'] == 2 for r in results)
    assert all(r['iterations'] <= 10 for r in results)
    assert len(app.get('download_button')) == 4
    assert len(app.dataframe) == 7  # comparison plus PWM/sites/heldout for two methods


def test_streamlit_invalid_upload_is_visible_error(monkeypatch):
    import streamlit as st
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(st, 'file_uploader', lambda label, **kw: BytesIO(b'>bad\nACGT!\n') if label.startswith('上传发现集') else None)
    app = AppTest.from_file(str(ROOT / 'app.py'), default_timeout=30).run()
    app.button[0].click().run()
    assert not app.exception
    assert len(app.error) == 1
    assert 'Invalid DNA' in app.error[0].value
    assert 'results' not in app.session_state


def test_cli_rejects_reverse_complement_overlap(sample_files, tmp_path):
    from biomotiflab.core import reverse_complement
    from biomotiflab.io import read_fasta
    train, _ = sample_files
    test = tmp_path / 'rc.fa'
    write_fasta([reverse_complement(read_fasta(train)[0])], test)
    run = subprocess.run([sys.executable, '-m', 'biomotiflab.cli', str(train), '--test-fasta', str(test),
                          '--starts', '1', '--iterations', '2', '--out', str(tmp_path / 'rc-result')],
                         cwd=ROOT, text=True, capture_output=True, check=False)
    assert run.returncode == 2
    assert 'reverse complements' in run.stderr


def test_streamlit_rejects_reverse_complement_overlap(monkeypatch, sample_files):
    import streamlit as st
    from streamlit.testing.v1 import AppTest
    from biomotiflab.core import reverse_complement
    from biomotiflab.io import read_fasta
    train, _ = sample_files
    rc = ('>rc\n' + reverse_complement(read_fasta(train)[0]) + '\n').encode()
    monkeypatch.setattr(st, 'file_uploader', lambda label, **kw: BytesIO(train.read_bytes() if label.startswith('上传发现集') else rc))
    app = AppTest.from_file(str(ROOT / 'app.py'), default_timeout=30).run()
    app.button[0].click().run()
    assert not app.exception
    assert len(app.error) == 1
    assert '反向互补' in app.error[0].value
