#!/usr/bin/env python3
"""Fetch pinned public biological examples; never silently accept changed bytes.

Only standard-library Python is needed. The CRP FASTA and MEME license are kept
local (see data/real/.gitignore); JASPAR's reference matrix is CC BY 4.0.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import ssl
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
CRP_URL = 'https://meme-suite.org/meme-software/meme_example_inputs/crp0.s'
CRP_SHA256 = 'dc81a5200ba8bbcb4f12d6ae1dc073f9f3ebeba26e43833bd43478c151972098'
JASPAR_URL = 'https://jaspar.elixir.no/api/v1/matrix/MA0139.1/'
LICENSE_URL = 'https://meme-suite.org/meme/doc/copyright.html'

def download(url: str, cafile: str | None) -> bytes:
    context = ssl.create_default_context(cafile=cafile)
    request = urllib.request.Request(url, headers={'User-Agent': 'BioMotifLab/0.1 (academic reproducibility)'})
    with urllib.request.urlopen(request, timeout=60, context=context) as response:
        data = response.read(2_000_001)
    if len(data) > 2_000_000:
        raise ValueError(f'Refusing unexpectedly large download: {url}')
    return data

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=ROOT / 'data' / 'real')
    p.add_argument('--cafile', help='Optional trusted CA bundle; TLS verification is never disabled.')
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    crp = download(CRP_URL, args.cafile)
    digest = hashlib.sha256(crp).hexdigest()
    if digest != CRP_SHA256:
        raise ValueError(f'CRP source changed: expected {CRP_SHA256}, received {digest}; inspect before updating the pin')
    # Fetch the governing notice before making the restricted example available.
    notice = download(LICENSE_URL, args.cafile)
    if b'educational' not in notice or b'Copyright' not in notice:
        raise ValueError('Unexpected MEME copyright response')
    pfm_bytes = download(JASPAR_URL, args.cafile)
    pfm = json.loads(pfm_bytes)
    if pfm.get('matrix_id') != 'MA0139.1' or set(pfm.get('pfm', {})) != set('ACGT'):
        raise ValueError('Unexpected JASPAR matrix response')
    # Metadata URLs may change; pin biological count matrix separately.
    canonical = json.dumps(pfm['pfm'], sort_keys=True, separators=(',', ':')).encode()
    expected = '0083d788dc856815962516190f51fe525d53c5789fd31c3c8cd7d862d783df96'
    if hashlib.sha256(canonical).hexdigest() != expected:
        raise ValueError('JASPAR MA0139.1 counts changed; inspect before updating the pin')
    for name, data in [('crp.fasta', crp), ('MEME-COPYRIGHT.html', notice), ('MA0139.1.json', pfm_bytes)]:
        dest = args.output / name
        tmp = dest.with_suffix(dest.suffix + '.tmp')
        tmp.write_bytes(data)
        tmp.replace(dest)
        print(f'{name}: {len(data)} bytes; SHA256={hashlib.sha256(data).hexdigest()}')

if __name__ == '__main__':
    main()
