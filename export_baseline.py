"""Export locally produced baseline results; do not use with untrusted .pt files."""
import json
import argparse
from pathlib import Path
import torch
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdDepictor
from run_checks import export_ready_mol

parser=argparse.ArgumentParser()
parser.add_argument('--run-dir',type=Path,required=True)
args=parser.parse_args()
root=args.run_dir.resolve()
manifest=json.loads((root/'run_manifest.json').read_text(encoding='utf8'))
if manifest['status']!='completed': raise ValueError('Only completed runs can be exported')
files=list(root.glob('*/result.pt'))
if len(files)!=1: raise RuntimeError(f'Expected one baseline result, found {len(files)}')
# This file was produced locally by run_baseline, not supplied externally.
results=torch.load(files[0],map_location='cpu',weights_only=False)
if len(results)!=manifest['requested']: raise ValueError('Sample count does not match run manifest')
writer=Chem.SDWriter(str(root/'generated.sdf')); rows=[]; mols=[]
for i,r in enumerate(results):
    mol=r['mol']
    if mol is None:
        rows.append({'sample':i,'valid':False}); continue
    Chem.SanitizeMol(mol); mol=export_ready_mol(mol); writer.write(mol); mols.append(mol)
    rows.append({'sample':i,'valid':True,'connected':len(Chem.GetMolFrags(mol))==1,
                 'smiles':Chem.MolToSmiles(mol),'heavy_atoms':mol.GetNumHeavyAtoms(),
                 'molecular_weight':Descriptors.MolWt(mol)})
writer.close()
if mols:
    drawings=[Chem.Mol(m) for m in mols]
    for mol in drawings: rdDepictor.Compute2DCoords(mol)
    Draw.MolsToGridImage(drawings,molsPerRow=1,subImgSize=(700,420)).save(str(root/'generated.png'))
report={'method':'Official DecompDiff reference prior, pretrained weights',
        'data_id':manifest['data_id'],'steps':manifest['steps'],
        'checkpoint_sha256':manifest['checkpoint_sha256'],'samples':len(results),'results':rows,
        'stereochemistry':'Inferred by RDKit from generated 3D coordinates; unvalidated.',
        'limits':'A single runtime test; no affinity, degradation, novelty benchmark or PROTAC fine-tuning claim.'}
(root/'report.json').write_text(json.dumps(report,indent=2),encoding='utf8')
print(json.dumps(report,indent=2))
