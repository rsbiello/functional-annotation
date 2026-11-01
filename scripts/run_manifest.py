#!/usr/bin/env python3
"""Reject stale results when configuration, code or input metadata changes."""
import hashlib
import json
import os
from pathlib import Path
import sys

out, config, code, *inputs = sys.argv[1:]
root = Path(code)
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def metadata(path):
    p = Path(path).resolve()
    st = p.stat()
    return [str(p), st.st_size, st.st_mtime_ns]
files = [root / 'functional_annotation.sh', *sorted((root / 'scripts').glob('*'))]
manifest = {'config': digest(Path(config)), 'code': {f.name: digest(f) for f in files if f.is_file()},
            'inputs': [], 'environment': {k: os.environ.get(k) for k in ('CONDA_PREFIX',)}}
for value in inputs:
    if not value:
        continue
    p = Path(value)
    if p.is_dir():
        manifest['inputs'].extend(metadata(f) for f in sorted(p.iterdir()) if f.is_file())
    else:
        manifest['inputs'].append(metadata(p))
path = Path(out) / 'run_manifest.json'
if path.exists():
    if json.loads(path.read_text()) != manifest:
        sys.exit('Inputs, databases, configuration or pipeline changed. Use a new OUTDIR (or a new PROJECT); existing results were not modified.')
else:
    if any((Path(out) / x).exists() for x in ('01_proteins', '02_swissprot', '06_final')):
        sys.exit('Untracked previous results exist. Use a new OUTDIR; legacy files cannot be resumed safely.')
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(manifest, indent=2) + '\n')
    temporary.replace(path)
