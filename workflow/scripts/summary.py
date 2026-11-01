"""Report protein-level annotation coverage, including valid zero-hit outputs."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path


def rows(path):
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith('#'):
                yield line.rstrip('\n').split('\t')


def main(sm):
    with open(sm.input.proteins) as handle:
        proteins = {line[1:].split()[0] for line in handle if line.startswith('>')}
    if not proteins:
        raise ValueError('Cannot summarize an empty protein set')
    categories = {name: set() for name in (
        'swissprot_hit', 'named_final', 'interproscan_any', 'interpro_entry', 'pfam',
        'go_interproscan', 'eggnog_hit', 'go_eggnog', 'kegg_ko')}
    categories['swissprot_hit'] = {r[0] for r in rows(sm.input.blast)}
    for row in rows(sm.input.gff):
        if len(row) == 9 and row[2] == 'mRNA':
            attr = dict(x.split('=', 1) for x in row[8].split(';') if '=' in x)
            if attr.get('Name') and attr.get('ID'):
                categories['named_final'].add(attr['ID'])
    for row in rows(sm.input.ipr):
        categories['interproscan_any'].add(row[0])
        if len(row) > 11 and row[11] not in ('', '-'):
            categories['interpro_entry'].add(row[0])
        if len(row) > 3 and row[3] == 'Pfam':
            categories['pfam'].add(row[0])
        if len(row) > 13 and row[13] not in ('', '-'):
            categories['go_interproscan'].add(row[0])
    for row in rows(sm.input.eggnog):
        categories['eggnog_hit'].add(row[0])
        if len(row) > 9 and row[9] not in ('', '-'):
            categories['go_eggnog'].add(row[0])
        if len(row) > 11 and row[11] not in ('', '-'):
            categories['kegg_ko'].add(row[0])
    categories['go_any'] = categories['go_interproscan'] | categories['go_eggnog']
    categories['any_annotation'] = categories['swissprot_hit'] | categories['interpro_entry'] | categories['eggnog_hit']
    if sm.input.rbh:
        with open(sm.input.rbh[0]) as handle:
            categories['ref_ortholog'] = {r['query_id'] for r in csv.DictReader(handle, delimiter='\t')}
        with open(sm.input.comparison) as handle:
            comparisons = list(csv.DictReader(handle, delimiter='\t'))
        for status in ('same', 'different'):
            categories['ref_name_' + status] = {r['mRNA_id'] for r in comparisons if r['comparison'] == status}
    with open(sm.output.summary, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['category', 'n_genes', 'percent'])
        for name, identifiers in categories.items():
            count = len(identifiers & proteins)
            writer.writerow([name, count, f'{100 * count / len(proteins):.1f}'])
        writer.writerow(['total_genes', len(proteins), '100.0'])
    metadata = []
    for filename in sm.input:
        path = Path(filename).resolve()
        stat = path.stat()
        metadata.append({'path': str(path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns})
    # Snakemake stores execution/code/Conda metadata in the working directory.
    report = {'completed_utc': datetime.now(timezone.utc).isoformat(),
              'config': dict(sm.params.settings), 'input_metadata': metadata,
              'note': 'Execution metadata and deployed environments: .snakemake in the working directory. Database versions must be recorded by the operator.'}
    Path(sm.output.provenance).write_text(json.dumps(report, indent=2) + '\n')
    Path(sm.log[0]).write_text(f'Summarized {len(proteins)} proteins.\n')


main(snakemake)
