#!/usr/bin/env python3
"""Validate FASTA, strip terminal stops, mask internal stops as X."""
from pathlib import Path
import sys
seen = set()
def records(path):
    name, sequence = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(sequence)
            name, sequence = line[1:].split()[0], []
        elif line.strip():
            if name is None:
                raise ValueError('Sequence before first FASTA header')
            sequence.append(line.strip())
    if name is not None:
        yield name, ''.join(sequence)
with open(sys.argv[2], 'w') as output:
    for name, sequence in records(sys.argv[1]):
        sequence = sequence.upper().rstrip('*').replace('*', 'X')
        if name in seen or not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
            raise ValueError(f'Duplicate ID, empty sequence or invalid protein characters: {name}')
        seen.add(name)
        output.write(f'>{name}\n{sequence}\n')
if not seen:
    raise ValueError('No proteins found')
