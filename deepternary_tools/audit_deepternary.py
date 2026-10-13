import csv
import argparse
import json
from pathlib import Path
import numpy as np
from rdkit import Chem

parser = argparse.ArgumentParser()
parser.add_argument('--out',type=Path,default=Path('outputs/deepternary-linker5-40seeds-2026-10-12'))
parser.add_argument('--reference',type=Path,default=Path('outputs/deepternary-linker5-inputs-2026-10-12/ligand.sdf'))
args = parser.parse_args()
out = args.out
reference = next(iter(Chem.SDMolSupplier(str(args.reference))))
reference_smiles = Chem.MolToSmiles(reference, isomericSmiles=True)
summaries = list(out.glob('summary_*.csv'))
assert len(summaries)==1
rows = list(csv.DictReader(summaries[0].open()))
records = []
ptable = Chem.GetPeriodicTable()
for row in rows:
    mol = next(iter(Chem.SDMolSupplier(str(out/f"ligand_pred_{row['seed']}.sdf"))))
    assert mol is not None and mol.GetNumAtoms() == reference.GetNumAtoms()
    xyz = mol.GetConformer().GetPositions()
    assert np.isfinite(xyz).all()
    proteins = np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in Path(row['complex_pred']).read_text().splitlines() if l.startswith('ATOM  ') and l[76:78].strip() != 'H'])
    external = np.linalg.norm(xyz[:,None]-proteins[None,:],axis=2)
    distance = np.linalg.norm(xyz[:,None]-xyz[None,:],axis=2)
    vdw = np.array([ptable.GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    overlaps = np.triu((Chem.GetDistanceMatrix(mol)>2)&(distance<.75*(vdw[:,None]+vdw[None,:])),1)
    extreme = []
    for bond in mol.GetBonds():
        a,b = bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        ratio = distance[a,b]/sum(ptable.GetRcovalent(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in [a,b])
        if ratio<.65 or ratio>1.4:
            extreme.append([a,b,float(distance[a,b])])
    stored = Chem.MolToSmiles(mol,isomericSmiles=True)==reference_smiles
    Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
    stereo = Chem.MolToSmiles(mol,isomericSmiles=True)==reference_smiles
    records.append(dict(seed=int(row['seed']), surrogate_score=float(row['pred_p2_rmsd']),
        complex_pred=row['complex_pred'],ca_clash_ratio=float(row['clash_ratio']), protein_pairs_below_1_5_A=int((external<1.5).sum()),
        ligand_nonbonded_overlaps=int(overlaps.sum()), extreme_bonds=extreme,
        stored_chemistry_matches=stored, geometry_stereochemistry_matches=stereo))
clean = [r for r in records if not r['protein_pairs_below_1_5_A'] and not r['ligand_nonbonded_overlaps'] and not r['extreme_bonds'] and r['geometry_stereochemistry_matches']]
report = dict(predictions=len(records), geometry_screen_pass=len(clean), best_screened=clean[0] if clean else None,
              records=records, limits='Conditional transferred-protein experiment; scores are surrogate predictions, not measured RMSD or affinity.')
(out/'geometry_audit.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
