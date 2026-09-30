#!/usr/bin/env python3
"""Validate the distributed files, input checksums, panel map and local documentation links."""
from pathlib import Path
import csv
import hashlib
import gzip
import zipfile
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
errors = []
def require(condition, message):
    if not condition:
        errors.append(message)

def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest() if hasattr(hashlib, 'file_digest') else hashlib.sha256(stream.read()).hexdigest()

def runtime_path(path):
    parts = path.relative_to(ROOT).parts
    return any(x in parts for x in ('.git', 'output', '__pycache__', 'renv')) or (path.suffix == '.rds' and 'raw' in parts) or path.name in ('.DS_Store', 'Rplots.pdf')

for module in range(1, 6):
    base = ROOT / 'figures' / f'figure{module}'
    for name in ('.here', '.Rprofile', 'renv.lock', 'scripts/run.R', 'scripts/00_check_input_data.R'):
        require((base / name).is_file(), f'Missing module file: {base.relative_to(ROOT)}/{name}')
    require(not (base / '.git').exists(), f'Nested Git history: figure{module}')

# Prefer the exact tracked distribution when Git metadata is available.
if (ROOT / '.git').exists():
    names = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
    files = [ROOT / name for name in names if name]
    if not files:
        files = [p for p in ROOT.rglob('*') if p.is_file() and not runtime_path(p)]
else:
    files = [p for p in ROOT.rglob('*') if p.is_file() and not runtime_path(p)]

sensitive = re.compile(r'(?:/(?:Users|home|Volumes)/[^\s"<>]+|gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|AKIA[0-9A-Z]{16}|-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----)')
for path in files:
    rel = path.relative_to(ROOT)
    require(path.is_file(), f'Missing distributed file: {rel}')
    if not path.is_file():
        continue
    require(not path.is_symlink(), f'Symlink in distribution: {rel}')
    require(path.stat().st_size < 50_000_000, f'Large file in Git distribution: {rel}')
    require(path.suffix.lower() not in ('.docx', '.pptx', '.pdf', '.png', '.jpg', '.rds', '.nd2', '.tif', '.tiff', '.fcs'), f'Unexpected artifact in source distribution: {rel}')
    require(path.name not in ('AGENTS.md', 'CLAUDE.md', '.env'), f'Unexpected local configuration: {rel}')
    if path.suffix.lower() == '.xlsx':
        with zipfile.ZipFile(path) as archive:
            texts = [(name, archive.read(name).decode(errors='replace')) for name in archive.namelist() if name.endswith(('.xml', '.rels'))]
    elif path.suffix.lower() == '.gz':
        with gzip.open(path, 'rt', errors='replace') as stream:
            texts = [('compressed contents', stream.read())]
    else:
        texts = [('', path.read_text(errors='replace'))]
    for part, text in texts:
        require(sensitive.search(text) is None, f'Check machine path or credential-like text: {rel} {part}')

with (ROOT / 'inputs.tsv').open() as stream:
    inputs = list(csv.DictReader(stream, delimiter='\t'))
require(bool(inputs), 'Input manifest is empty')
require(len({r['path'] for r in inputs}) == len(inputs), 'Duplicate input path')
for row in inputs:
    path = ROOT / row['path']
    require(path.is_file(), f'Missing input: {row["path"]}')
    if path.is_file():
        require(path.stat().st_size == int(row['bytes']), f'Input size mismatch: {row["path"]}')
        require(sha256(path) == row['sha256'], f'Input checksum mismatch: {row["path"]}')

with (ROOT / 'panels.tsv').open() as stream:
    panels = list(csv.DictReader(stream, delimiter='\t'))
expected = set()
for figure, main, supp in [(1,10,5),(2,9,17),(3,11,10),(4,8,8),(5,10,11)]:
    expected.update(f'{figure}{chr(65+i)}' for i in range(main))
    expected.update(f'S{figure}{chr(65+i)}' for i in range(supp))
require(len({r['panel'] for r in panels}) == len(panels), 'Duplicate panel label')
require({r['panel'] for r in panels} == expected, 'Panel coverage does not match the manuscript inventory')
for row in panels:
    require(row['kind'] in ('code', 'mixed', 'external'), f'Unknown panel type: {row["panel"]}')
    scripts = [p for p in row['scripts'].split(';') if p]
    require(bool(scripts) == (row['kind'] != 'external'), f'Inconsistent panel mapping: {row["panel"]}')
    for script in scripts:
        require((ROOT / script).is_file(), f'Missing panel script: {script}')

for path in files:
    if path.suffix != '.md' or not path.exists():
        continue
    for link in re.findall(r'\]\(([^)]+)\)', path.read_text()):
        if '://' in link or link.startswith(('#', 'mailto:')):
            continue
        target = link.split('#', 1)[0]
        require((path.parent / target).exists(), f'Broken local link in {path.relative_to(ROOT)}: {target}')

if errors:
    print('\n'.join('ERROR: ' + message for message in errors), file=sys.stderr)
    sys.exit(1)
print(f'PASS: {len(files)} source files, {len(inputs)} checksummed inputs, {len(panels)} manuscript panels; local links valid.')
