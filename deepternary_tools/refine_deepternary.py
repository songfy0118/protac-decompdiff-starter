"""Local geometry repair only; neither docking validation nor binding-energy estimation."""
import json
import argparse
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.optimize import minimize

parser = argparse.ArgumentParser()
parser.add_argument('--out',type=Path,default=Path('outputs/deepternary-linker5-40seeds-2026-10-12'))
parser.add_argument('--reference',type=Path,default=Path('outputs/deepternary-linker5-inputs-2026-10-12/ligand.sdf'))
parser.add_argument('--maxiter',type=int,default=1200)
args = parser.parse_args()
out = args.out
audit = json.loads((out/'geometry_audit.json').read_text())
reference = next(iter(Chem.SDMolSupplier(str(args.reference))))
target_smiles = Chem.MolToSmiles(reference,isomericSmiles=True)
results = []
# Only refine poses whose observed 3D stereochemistry agrees with the input.
eligible = [r for r in audit['records'] if r['geometry_stereochemistry_matches']][:5]
for record in eligible:
    seed = record['seed']
    mol = next(iter(Chem.SDMolSupplier(str(out/f'ligand_pred_{seed}.sdf'))))
    Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
    mol = Chem.AddHs(mol,addCoords=True)
    nheavy = reference.GetNumAtoms()
    start = mol.GetConformer().GetPositions().copy()
    complex_path = Path(record.get('complex_pred',out/f'complex_pred_linker5_conditional_{seed}.pdb'))
    protein = np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in complex_path.read_text().splitlines() if l.startswith('ATOM  ') and l[76:78].strip()!='H'])
    # Limit to atoms near the starting ligand, while retaining a 6-A margin.
    protein = protein[np.linalg.norm(protein[:,None]-start[None,:nheavy],axis=2).min(1)<8]
    props = AllChem.MMFFGetMoleculeProperties(mol,mmffVariant='MMFF94s')
    force = AllChem.MMFFGetMoleculeForceField(mol,props)
    def objective(flat):
        coords = flat.reshape(-1,3)
        energy = force.CalcEnergy(tuple(flat))
        grad = np.array(force.CalcGrad(tuple(flat))).reshape(-1,3)
        delta = coords[:nheavy,None,:]-protein[None,:,:]
        distance = np.linalg.norm(delta,axis=2)
        deficit = np.maximum(2.2-distance,0)
        energy += 50*float((deficit**2).sum())
        grad[:nheavy] += (-100*deficit[:,:,None]*delta/np.maximum(distance[:,:,None],1e-8)).sum(1)
        displacement = coords[:nheavy]-start[:nheavy]
        energy += 2*float((displacement**2).sum())
        grad[:nheavy] += 4*displacement
        return energy,grad.ravel()
    optimized = minimize(objective,start.ravel(),jac=True,method='L-BFGS-B',options={'maxiter':args.maxiter,'ftol':1e-8})
    xyz = optimized.x.reshape(-1,3)
    for i,p in enumerate(xyz):
        mol.GetConformer().SetAtomPosition(i,tuple(float(x) for x in p))
    mol = Chem.RemoveHs(mol)
    Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
    stereo = Chem.MolToSmiles(mol,isomericSmiles=True)==target_smiles
    with Chem.SDWriter(str(out/f'ligand_refined_{seed}.sdf')) as writer:
        writer.write(mol)
    heavy = xyz[:nheavy]
    protein_pairs = int((np.linalg.norm(heavy[:,None]-protein[None,:],axis=2)<1.5).sum())
    ptable = Chem.GetPeriodicTable()
    vdw = np.array([ptable.GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    dist = np.linalg.norm(heavy[:,None]-heavy[None,:],axis=2)
    internal = int(np.triu((Chem.GetDistanceMatrix(mol)>2)&(dist<.75*(vdw[:,None]+vdw[None,:])),1).sum())
    results.append(dict(seed=seed,optimizer_converged=bool(optimized.success),iterations=int(optimized.nit),
        heavy_atom_displacement_rmsd_A=float(np.sqrt(((heavy-start[:nheavy])**2).sum(1).mean())),
        protein_pairs_below_1_5_A=protein_pairs,ligand_nonbonded_overlaps=internal,
        geometry_stereochemistry_matches=stereo))
report = dict(results=results,limits='Local MMFF94s geometry repair with protein repulsion and positional restraints. Protein coordinates fixed. No binding-energy or activity claims.')
(out/'refinement_audit.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
