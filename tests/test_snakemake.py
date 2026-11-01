from pathlib import Path
import tempfile, subprocess, os, json
root=Path(__file__).resolve().parents[1]
import sys
SNK=os.environ.get('SNAKEMAKE', 'snakemake')
t=Path(tempfile.mkdtemp(prefix='annotation-test-')); b=t/'bin'; b.mkdir()
gff='##gff-version 3\nchr1\ttest\tgene\t1\t9\t.\t+\t.\tID=g1\nchr1\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n'
mock='''#!/usr/bin/env python3
import sys,os
from pathlib import Path
name=Path(sys.argv[0]).name; a=sys.argv[1:]
def arg(k): return a[a.index(k)+1]
with open(os.environ['CALLS'],'a') as f:f.write(name+' '+str(a)+'\\n')
if '--version' in a or '--help' in a or '-version' in a or a==['version']: print('test version'); sys.exit()
if name=='agat_sp_keep_longest_isoform.pl': Path(arg('--output')).write_text(Path(arg('--gff')).read_text())
elif name=='agat_sp_extract_sequences.pl': Path(arg('--output')).write_text('>t1 description\\nMAA*\\n')
elif name=='diamond':
 if a[0]=='blastp': Path(arg('-o')).write_text('')
 elif a[0]=='makedb': Path(arg('-d')+'.dmnd').write_text('index')
elif name=='interproscan.sh':
 Path(arg('-o')).write_text('t1\\tmd5\\t3\\tPfam\\tPF0001\\tdesc\\t1\\t3\\t1e-9\\tT\\tdate\\tIPR0001\\tdesc\\tGO:0003677\\n' if os.environ.get('MATCHES') else '')
 if os.environ.get('FAIL_IPR'): sys.exit(7)
elif name=='emapper.py': Path(arg('--output_dir'),arg('-o')+'.emapper.annotations').write_text('t1\\tseed\\t1e-9\\t100\\tog\\tlevel\\tS\\tdesc\\tabc\\tGO:0003677\\t-\\tko:K00001\\n' if os.environ.get('MATCHES') else '# empty hits\\n')
elif name=='agat_sp_manage_functional_annotation.pl':
 p=Path(arg('--output'));p.mkdir();(p/'result.gff3').write_text(Path(arg('--gff')).read_text())
'''
for name in ['agat_sp_keep_longest_isoform.pl','agat_sp_extract_sequences.pl','diamond','interproscan.sh','emapper.py','agat_sp_manage_functional_annotation.pl']:
 p=b/name;p.write_text(mock);p.chmod(0o755)
for name in ['swiss','egg']: (t/name).mkdir()
for name in ['swiss/uniprot_sprot.fasta','swiss/uniprot_sprot.dmnd','egg/eggnog.db','egg/eggnog.taxa.db','egg/eggnog_proteins.dmnd','genome.fa']:(t/name).write_text('fixture')
(t/'genes.gff3').write_text(gff)

c=t/'config.yaml'
settings={'project': 'project_a', 'species': 'test_species', 'genome': str(t/'genome.fa'),
          'gff': str(t/'genes.gff3'), 'swissprot_fasta': str(t/'swiss/uniprot_sprot.fasta'),
          'swissprot_db': str(t/'swiss/uniprot_sprot.dmnd'), 'eggnog_db': str(t/'egg'),
          'outdir': str(t/'results with spaces'), 'threads': 2}
c.write_text(json.dumps(settings))
env=dict(os.environ, PATH=str(b)+':'+str(Path(sys.executable).parent)+':'+os.environ['PATH'],
         CALLS=str(t/'calls'), XDG_CACHE_HOME=str(t/'cache'))
def run(ok=True, flags=(), **extra):
 r=subprocess.run([SNK, '--snakefile', str(root/'Snakefile'), '--directory', str(t),
                   '--configfile', str(c), '--cores', '2', '--rerun-incomplete', *flags],
                  env=dict(env, **extra), capture_output=True, text=True)
 assert (r.returncode==0)==ok, r.stdout+r.stderr
 return r
run(flags=('--dry-run',))
run(False, FAIL_IPR='yes')
run()
out=Path(settings['outdir'])
assert (out/'06_final/test_species.functional.gff3').read_text()==gff
assert 'total_genes\t1' in (out/'06_final/test_species.functional_annotation_summary.tsv').read_text()
before=(t/'calls').read_text(); run(); assert (t/'calls').read_text()==before
# A postprocessing setting change must not rerun expensive search rules.
settings['postprocess']={'name_case': 'upper'};c.write_text(json.dumps(settings));run()
assert (t/'calls').read_text()==before
# Changing a search threshold must rerun DIAMOND and its dependent merge.
settings['search']={'evalue': 1e-8};c.write_text(json.dumps(settings));run()
assert 'blastp' in (t/'calls').read_text()[len(before):]
# An optional generic reference with zero reciprocal hits remains valid.
ref=t/'reference.fa';ref.write_text('>r1\nMAAA\n')
settings['reference']={'species':'reference', 'proteome':str(ref), 'name_mode':'fill'}
c.write_text(json.dumps(settings));run()
assert (out/'05_reference/rbh.tsv').read_text().startswith('query_id\t')
assert 'ref_ortholog\t0\t0.0' in (out/'06_final/test_species.functional_annotation_summary.tsv').read_text()
# A second assembly config uses an isolated output directory.
settings['project']='project_b';settings['outdir']=str(t/'project_b');c.write_text(json.dumps(settings));run()
assert (t/'project_b/06_final/test_species.functional.gff3').exists()
print('PASS: Snakemake dry-run, failure recovery, zero hits, no-op resume, selective config reruns, optional reference, project isolation, paths with spaces')
print(t)

settings['project']='with_hits';settings['outdir']=str(t/'with_hits');settings['reference']={}
settings['postprocess']={'name_case':'keep'};c.write_text(json.dumps(settings));run(MATCHES='yes')
final=(t/'with_hits/06_final/test_species.functional.gff3').read_text()
assert final.count('Name=abc')==2,final
summary=(t/'with_hits/06_final/test_species.functional_annotation_summary.tsv').read_text()
assert 'pfam\t1\t100.0' in summary and 'kegg_ko\t1\t100.0' in summary,summary
print('PASS: positive domain/GO/KEGG hits and eggNOG names propagated to transcript and gene')
