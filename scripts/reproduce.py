#!/usr/bin/env python3
"""Run one or more manuscript figure workflows in isolated R processes."""
from pathlib import Path
import argparse, subprocess, sys
ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--figures', type=int, nargs='+', choices=range(1,6), default=list(range(1,6)))
p.add_argument('--check-inputs', action='store_true')
p.add_argument('--light', action='store_true', help='Figure 3 only: screen time courses and editing efficiency; does not build transcriptomic panels')
a = p.parse_args()
if a.light and a.figures != [3]: p.error('--light requires --figures 3')
for n in a.figures:
    module = ROOT / 'figures' / f'figure{n}'
    script = 'scripts/00_check_input_data.R' if a.check_inputs else 'scripts/run.R'
    cmd = ['Rscript', script] + (['--light'] if a.light else [])
    print(f'Figure {n}: {" ".join(cmd)}', flush=True)
    subprocess.run(cmd, cwd=module, check=True)
