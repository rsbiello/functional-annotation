from pathlib import Path
import tempfile, subprocess, os, json
root=Path(__file__).resolve().parents[1]
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
 Path(arg('-b')+'.tsv').write_text('')
 if os.environ.get('FAIL_IPR'): sys.exit(7)
elif name=='emapper.py': Path(arg('--output_dir'),arg('-o')+'.emapper.annotations').write_text('# empty hits\\n')
elif name=='agat_sp_manage_functional_annotation.pl':
 p=Path(arg('--output'));p.mkdir();(p/'result.gff3').write_text(Path(arg('--gff')).read_text())
'''
for name in ['agat_sp_keep_longest_isoform.pl','agat_sp_extract_sequences.pl','diamond','interproscan.sh','emapper.py','agat_sp_manage_functional_annotation.pl']:
 p=b/name;p.write_text(mock);p.chmod(0o755)
for name in ['swiss','egg']: (t/name).mkdir()
for name in ['swiss/uniprot_sprot.fasta','swiss/uniprot_sprot.dmnd','egg/eggnog.db','egg/eggnog.taxa.db','egg/eggnog_proteins.dmnd','genome.fa']:(t/name).write_text('fixture')
(t/'genes.gff3').write_text(gff)
base='SPECIES="test_species"\nPROJECT="project_a"\nGENOME="genome.fa"\nGFF="genes.gff3"\nDB_DIR="swiss"\nEGGNOG_DB="egg"\nOUTPUT_ROOT="results"\nCONDA_ENV=""\n'
c=t/'config.sh';c.write_text(base)
env=dict(os.environ,PATH=str(b)+':'+os.environ['PATH'],CALLS=str(t/'calls'))
def run(ok=True,**extra):
 r=subprocess.run(['bash',str(root/'functional_annotation.sh'),str(c)],env=dict(env,**extra),capture_output=True,text=True,cwd='/tmp')
 assert (r.returncode==0)==ok, r.stdout+r.stderr
 return r
run(False,FAIL_IPR='yes');run()
out=t/'results/project_a/test_species'
assert (out/'06_final/test_species.functional.gff3').read_text()==gff
before=(t/'calls').read_text();run();after=(t/'calls').read_text()[len(before):]
assert 'blastp' not in after and '--output' not in after,after
assert 'total_genes\t1' in (out/'06_final/test_species.functional_annotation_summary.tsv').read_text()
c.write_text(base+'NAME_CASE="upper"\n');assert 'changed' in run(False).stderr
c.write_text(base.replace('project_a','project_b'));run();assert (t/'results/project_b/test_species/06_final/test_species.functional.gff3').exists()
c.write_text(base);(out/'.running').mkdir();assert 'locked' in run(False).stderr;(out/'.running').rmdir()
# Empty eggNOG input and multiple passes preserve all GFF records.
empty=t/'empty';empty.touch()
r=subprocess.run(['awk','-F\t','-v','OFS=\t','-f',str(root/'scripts/add_eggnog_names.awk'),str(empty),str(t/'genes.gff3'),str(t/'genes.gff3')],capture_output=True,text=True,check=True)
assert r.stdout==gff,r.stdout
print('PASS: failure recovery, zero-hit summary, resume, changed-config rejection, project isolation, locking, empty eggNOG input')
print(t)
# Generic references retain IDs, Ensembl isoforms collapse to longest gene record.
ref=t/'reference.fa';ref.write_text('>generic.1\nMAA\n>p1 gene:g2.1 gene_symbol:abc\nMA\n>p2 gene:g2.1 gene_symbol:abc\nMAAA\n')
r=subprocess.run(['awk','-v',f'tsv={t}/ref.tsv','-f',str(root/'scripts/prepare_reference.awk'),str(ref)],capture_output=True,text=True,check=True)
assert '>generic.1\nMAA' in r.stdout and '>g2\nMAAA' in r.stdout
rbh=t/'rbh.tsv';rbh.write_text('query_id\tref_gene_id\tpident\tevalue\tref_symbol\tref_description\nt1\tr1\t99\t1e-9\tabc\tprotein;description\n')
r=subprocess.run(['awk','-F\t','-v','OFS=\t','-v','mode=fill','-v','upper=keep','-v','ref=plant','-f',str(root/'scripts/apply_reference_names.awk'),str(rbh),str(t/'genes.gff3'),str(t/'genes.gff3')],capture_output=True,text=True,check=True)
assert r.stdout.count('Name=abc')==2 and 'protein%3Bdescription' in r.stdout
raw=t/'raw.fa';raw.write_text('>t1 details\nMA*AA*\n')
subprocess.run(['python3',str(root/'scripts/clean_proteins.py'),str(raw),str(t/'clean.fa')],check=True)
assert (t/'clean.fa').read_text()=='>t1\nMAXAA\n'
print('PASS: generic and Ensembl references, reference name propagation/escaping, internal and terminal stops')
