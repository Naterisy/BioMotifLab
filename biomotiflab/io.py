"""Strict FASTA parsing. Ambiguous DNA bases become N; N-overlapping windows are skipped."""
from pathlib import Path

DNA = frozenset('ACGTRYSWKMBDHVN')


def normalize_sequence(sequence):
    sequence = ''.join(str(sequence).split()).upper()
    if not sequence:
        raise ValueError('DNA sequence is empty')
    invalid = set(sequence) - DNA
    if invalid:
        raise ValueError('Invalid DNA character(s): ' + ', '.join(sorted(invalid)))
    return ''.join(c if c in 'ACGT' else 'N' for c in sequence)


def read_fasta(path_or_text, is_text=False, return_names=False):
    """Return sequence strings (or names, sequences); reject malformed/empty records."""
    text = str(path_or_text) if is_text else Path(path_or_text).read_text(encoding='utf-8-sig')
    names, sequences, chunks = [], [], []
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(';'):
            continue
        if line.startswith('>'):
            if names:
                sequences.append(normalize_sequence(''.join(chunks)))
            name = line[1:].strip()
            if not name:
                raise ValueError(f'Empty FASTA header at line {number}')
            names.append(name)
            chunks = []
        else:
            if not names:
                raise ValueError(f'Expected FASTA header before line {number}')
            chunks.append(line)
    if not names:
        raise ValueError('No FASTA records found')
    sequences.append(normalize_sequence(''.join(chunks)))
    return (names, sequences) if return_names else sequences


def write_fasta(sequences, path, names=None):
    seqs = [normalize_sequence(s) for s in sequences]
    names = list(names) if names is not None else [f'seq_{i + 1}' for i in range(len(seqs))]
    if len(names) != len(seqs):
        raise ValueError('Number of names and sequences must match')
    if any('\n' in n or '\r' in n or not n for n in names):
        raise ValueError('FASTA names must be nonempty single lines')
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(''.join(f'>{name}\n' + '\n'.join(s[i:i + 80] for i in range(0, len(s), 80)) + '\n' for name, s in zip(names, seqs)), encoding='utf-8')
