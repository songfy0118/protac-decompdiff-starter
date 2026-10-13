"""Run the conditional DeepTernary smoke experiment with explicit local paths."""
import json
import os
import argparse
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser()
parser.add_argument('--seeds', type=int, default=1)
parser.add_argument('--inputs',type=Path)
parser.add_argument('--out',type=Path)
parser.add_argument('--name',default='linker5_conditional')
parser.add_argument('--upstream',type=Path)
args = parser.parse_args()
upstream = (args.upstream or root/'work/deepternary-upstream').resolve()
inputs = (args.inputs or root/'outputs/deepternary-linker5-inputs-2026-10-12').resolve()
output = (args.out or root/f'outputs/deepternary-linker5-{args.seeds}seeds-2026-10-12').resolve()
output.mkdir(exist_ok=True)
command = [sys.executable, 'predict.py', '--task', 'PROTAC', '--name', args.name,
           '--lig', str(inputs/'ligand.sdf'), '--p1', str(inputs/'protein1.pdb'),
           '--p2', str(inputs/'protein2.pdb'), '--unbound-lig1', str(inputs/'lig1.pdb'),
           '--unbound-lig2', str(inputs/'lig2.pdb'), '--lig1-mask', str(inputs/'lig1.pdb'),
           '--lig2-mask', str(inputs/'lig2.pdb'), '--outdir', str(output), '--seeds', str(args.seeds),
           '--device', 'cpu', '--config', 'deepternary/configs/protac.py',
           '--checkpoint', 'output/checkpoints/PROTAC/checkpoint.pth']
(output/'command.json').write_text(json.dumps(command, indent=2))
with (output/'run.log').open('w') as log:
    temp = root/'work/deepternary-tmp'
    temp.mkdir(exist_ok=True)
    env = dict(os.environ, DGLBACKEND='pytorch', TEMP=str(temp), TMP=str(temp))
    result = subprocess.run(command, cwd=upstream, env=env, stdout=log, stderr=subprocess.STDOUT)
print('Exit code:', result.returncode)
print((output/'run.log').read_text(errors='replace')[-6000:])
sys.exit(result.returncode)
