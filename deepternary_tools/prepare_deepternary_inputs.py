"""Prepare a conditional pose experiment, retaining the transferred-pose caveat."""
import hashlib
import argparse
import json
from pathlib import Path
import shutil
from rdkit import Chem

parser = argparse.ArgumentParser()
parser.add_argument('--case',type=Path,default=Path('outputs/friend-linker-cases-2026-10-12/5_'))
parser.add_argument('--ligand',type=Path)
parser.add_argument('--protein1',type=Path)
parser.add_argument('--out',type=Path,default=Path('outputs/deepternary-linker5-inputs-2026-10-12'))
args = parser.parse_args()
case,out = args.case,args.out
out.mkdir(exist_ok=True)
ligand = args.ligand or case/'ak1690_observed.sdf'
mol = next(x for x in Chem.SDMolSupplier(str(ligand)) if x)
tasks = json.loads((case/'generation_tasks.json').read_text())
groups = {t['role']: t['generate_indices'] for t in tasks}
if set(groups)!= {'warhead','e3_ligand'}:
    raise ValueError('Expected exactly one task for each endpoint')
if any(not ids or len(ids)!=len(set(ids)) or any(i<0 or i>=mol.GetNumAtoms() for i in ids) for ids in groups.values()):
    raise ValueError('Endpoint indices are empty, duplicated, or outside the ligand atom range')
if set(groups['warhead']) & set(groups['e3_ligand']):
    raise ValueError('Endpoint index sets overlap')
for role, filename in [('warhead', 'lig2.pdb'), ('e3_ligand', 'lig1.pdb')]:
    if role not in groups:
        raise ValueError(f'Missing role {role}; found {list(groups)}')
    keep = set(groups[role])
    fragment = Chem.RWMol(mol)
    for i in reversed(range(mol.GetNumAtoms())):
        if i not in keep:
            fragment.RemoveAtom(i)
    fragment = fragment.GetMol()
    Chem.SanitizeMol(fragment)
    Chem.MolToPDBFile(fragment, str(out/filename))
    reread = Chem.MolFromPDBFile(str(out/filename))
    assert reread is not None and reread.GetNumAtoms() == len(keep)
    print(role, len(keep))
for source, target in [('ak1690_observed.sdf', 'ligand.sdf'),
                       ('crbn_transferred_hypothesis.pdb', 'protein1.pdb'),
                       ('stat6_observed.pdb', 'protein2.pdb')]:
    origin = ligand if target=='ligand.sdf' else args.protein1 if target=='protein1.pdb' and args.protein1 else case/source
    shutil.copyfile(origin, out/target)
manifest = {
    'status': 'conditional_pose_experiment_inputs',
    'source_case': str(case), 'source_ligand':str(ligand), 'source_protein1':str(args.protein1 or case/'crbn_transferred_hypothesis.pdb'),
    'limitations': ['Ligand bond orders are reconstructed hypotheses.',
                    'Protein1 is a transferred structure; see source audit for species and alignment.',
                    'Endpoint coordinates come from the input ligand, not new experimentally observed bound ligands.',
                    'Predicted RMSD is a model ranking score, not measured structural accuracy.'],
    'files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},
}
(out/'manifest.json').write_text(json.dumps(manifest, indent=2))
print(out)
