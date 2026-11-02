import importlib.util
from pathlib import Path
import unittest
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('cleaner', root/'scripts/clean_gff_attributes.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
clean = module.clean_attributes


def row(kind, attrs):
    return f'chr1\tGnomon\t{kind}\t1\t402\t.\t+\t0\t{attrs}\n'


class CleanupTests(unittest.TestCase):
    def test_placeholder(self):
        old = row('gene', 'ID=gene-SagSer_016351;Name=SagSer_016351,FAU;locus_tag=SagSer_016351')
        self.assertEqual(clean(old), old.replace('Name=SagSer_016351,FAU', 'Name=FAU'))

    def test_products_and_encoded_comma(self):
        old = row('mRNA', 'ID=t1;Parent=g1;product=Fusion%2C protein,fusion%2C protein,other protein')
        self.assertEqual(clean(old), old.replace(',fusion%2C protein', ''))

    def test_conflicts_and_locus_only_unchanged(self):
        for names in ('g1,FAU,OTHER', 'g1', 'FAU,OTHER'):
            old = row('gene', f'ID=g1;Name={names};locus_tag=g1')
            self.assertEqual(clean(old), old)

    def test_cds_and_fasta_unchanged(self):
        old = row('CDS', 'ID=c1;Name=g1,FAU;locus_tag=g1;Parent=t1') + '##FASTA\n>chr1\nACGT\n'
        self.assertEqual(clean(old), old)

    def test_ontology_and_reference_quote(self):
        old = row('mRNA', 'ID=t1;Ontology_term=GO:0003735,GO:0003735,GO:0005762;ref_product="mitochondrial ribosomal protein L49')
        new = clean(old)
        self.assertEqual(new, old.replace('GO:0003735,GO:0003735', 'GO:0003735').replace('ref_product="', 'ref_product='))
        self.assertEqual(clean(new), new)

    def test_quotes_and_escapes(self):
        for value in ('protein "domain"', 'protein%3Bname%2Cextra', 'protein', '"domain" protein'):
            old = row('mRNA', 'ID=t1;ref_product=' + value)
            self.assertEqual(clean(old), old)
        old = row('mRNA', 'ID=t1;ref_product="protein"')
        self.assertEqual(clean(old), old.replace('"protein"', 'protein'))

    def test_reference_parser(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / 'genes.tsv'
            fasta = ('>p1 gene:g1 gene_symbol:MRPL49 description:"mitochondrial ribosomal protein L49 [Source:NCBI;Acc:1]"\nMAA\n'
                     '>p2 gene:g2 description:unquoted protein [Source:NCBI;Acc:2]\nMAA\n'
                     '>p3 gene:g3 description:"protein without source"\nMAA\n')
            subprocess.run(['awk', '-v', f'tsv={table}', '-f', str(root/'scripts/prepare_reference.awk')],
                           input=fasta, text=True, capture_output=True, check=True)
            descriptions = {r.split('\t')[0]: r.split('\t')[2] for r in table.read_text().splitlines()}
            self.assertEqual(descriptions, {'g1': 'mitochondrial ribosomal protein L49',
                                           'g2': 'unquoted protein', 'g3': 'protein without source'})

    def test_idempotent(self):
        old = row('gene', 'ID=g1;Name=g1,FAU;locus_tag=g1;product=Protein,protein')
        self.assertEqual(clean(clean(old)), clean(old))


if __name__ == '__main__':
    unittest.main()
