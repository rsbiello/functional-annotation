#!/usr/bin/env python3
"""Conservative cleanup of merged names/products; never change feature structure."""
import argparse
from pathlib import Path
from urllib.parse import unquote


def clean_attributes(text):
    lines = []
    in_fasta = False
    for line in text.splitlines(keepends=True):
        if line.startswith('##FASTA'):
            in_fasta = True
        if in_fasta or line.startswith('#') or not line.strip():
            lines.append(line)
            continue
        ending = '\r\n' if line.endswith('\r\n') else '\n' if line.endswith('\n') else ''
        fields = line.rstrip('\r\n').split('\t')
        if len(fields) != 9:
            raise ValueError('Expected nine tab-separated GFF columns')
        parts = fields[8].split(';')
        attrs = dict(p.split('=', 1) for p in parts if '=' in p)
        replacements = {}
        if fields[2] == 'gene' and attrs.get('Name'):
            names = list(dict.fromkeys(attrs['Name'].split(',')))
            locus = unquote(attrs.get('locus_tag', ''))
            # Remove a locus placeholder only when exactly one distinct symbol remains.
            symbols = [n for n in names if unquote(n) != locus]
            if locus and len(names) > 1 and len(symbols) == 1:
                replacements['Name'] = symbols[0]
        if 'product' in attrs:
            seen, products = set(), []
            for value in attrs['product'].split(','):
                key = ' '.join(unquote(value).split()).casefold()
                if key not in seen:
                    products.append(value)
                    seen.add(key)
            if len(products) < len(attrs['product'].split(',')):
                replacements['product'] = ','.join(products)
        if replacements:
            parts = [p.split('=', 1)[0] + '=' + replacements[p.split('=', 1)[0]]
                     if '=' in p and p.split('=', 1)[0] in replacements else p for p in parts]
            fields[8] = ';'.join(parts)
            line = '\t'.join(fields) + ending
        lines.append(line)
    return ''.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input')
    parser.add_argument('output', help='A new output file; the input is never overwritten')
    args = parser.parse_args()
    if Path(args.input).resolve() == Path(args.output).resolve():
        parser.error('Choose a different output path')
    text = clean_attributes(Path(args.input).read_text())
    with open(args.output, 'x') as handle:
        handle.write(text)


if __name__ == '__main__':
    main()
