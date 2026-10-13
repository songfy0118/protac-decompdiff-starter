"""Geometry and bond-inference robustness checks, not affinity scoring."""
import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors
from audit_linker_xyz import read_xyz


def graph_from_xyz(elements, positions):
    from openbabel import openbabel as ob
    conversion = ob.OBConversion(); conversion.SetInFormat('xyz')
    text = str(len(elements))+'\n\n'+'\n'.join(
        f'{e} {p[0]:.9f} {p[1]:.9f} {p[2]:.9f}' for e,p in zip(elements, positions))
    molecule = ob.OBMol()
    if not conversion.ReadString(molecule, text):
        raise ValueError('XYZ inference failed')
    return {(min(b.GetBeginAtomIdx()-1,b.GetEndAtomIdx()-1),
             max(b.GetBeginAtomIdx()-1,b.GetEndAtomIdx()-1)):b.GetBondOrder()
            for b in ob.OBMolBondIter(molecule)}


def inference_stability(path, trials=30, sigma=.03, seed=2026):
    elements, positions = read_xyz(path.read_text())
    original = {p:t for p,t in graph_from_xyz(elements,positions).items() if max(p)>=63}
    random = np.random.default_rng(seed)
    same_edges = same_orders = 0; variants = {}
    for _ in range(trials):
        noisy = positions.copy(); noisy[63:] += random.normal(0,sigma,noisy[63:].shape)
        inferred = {p:t for p,t in graph_from_xyz(elements,noisy).items() if max(p)>=63}
        same_edges += inferred.keys()==original.keys()
        same_orders += inferred==original
        key = str(sorted(inferred.items())); variants[key] = variants.get(key,0)+1
    return {'trials':trials,'coordinate_jitter_sigma_A':sigma,'seed':seed,
            'linker_edge_agreement':same_edges/trials,'linker_bond_order_agreement':same_orders/trials,
            'distinct_linker_graphs':len(variants),'graph_hypotheses_counts':variants,
            'limits':'Numerical robustness of one inference method, not proof of unique/correct chemistry.'}


def assess(molecule, protein_path=None):
    molecule = Chem.RemoveHs(molecule)
    if not molecule.GetNumConformers():
        raise ValueError('Geometry assessment requires coordinates')
    positions = molecule.GetConformer().GetPositions()
    if not np.isfinite(positions).all():
        raise ValueError('Nonfinite molecular coordinates')
    ptable = Chem.GetPeriodicTable()
    distances = np.linalg.norm(positions[:,None,:]-positions[None,:,:],axis=2)
    topological = Chem.GetDistanceMatrix(molecule)
    vdw = np.array([ptable.GetRvdw(a.GetAtomicNum()) for a in molecule.GetAtoms()])
    nonlocal_pairs = np.triu((topological>2)&(distances<.75*(vdw[:,None]+vdw[None,:])),1)
    stretched = []
    for bond in molecule.GetBonds():
        a,b = bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        expected = sum(ptable.GetRcovalent(molecule.GetAtomWithIdx(i).GetAtomicNum()) for i in (a,b))
        ratio = float(distances[a,b]/expected)
        if ratio<.65 or ratio>1.4:
            stretched.append({'atoms':[a,b],'length_A':float(distances[a,b]),'covalent_radius_ratio':ratio})
    result = {'heavy_atoms':molecule.GetNumAtoms(),'connected':len(Chem.GetMolFrags(molecule))==1,
        'molecular_weight':Descriptors.MolWt(molecule),'cLogP':Crippen.MolLogP(molecule),
        'TPSA':rdMolDescriptors.CalcTPSA(molecule),'rotatable_bonds':rdMolDescriptors.CalcNumRotatableBonds(molecule),
        'nonbonded_overlap_pairs':int(nonlocal_pairs.sum()),'extreme_bond_lengths':stretched,
        'lipinski_filter_applied':False}
    if protein_path:
        protein_positions = []
        for line in protein_path.read_text().splitlines():
            if line.startswith('ATOM  ') and line[76:78].strip()!='H':
                protein_positions.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        external = np.linalg.norm(positions[:,None,:]-np.array(protein_positions)[None,:,:],axis=2)
        result['protein_pairs_below_1_5_A'] = int((external<1.5).sum())
        result['minimum_protein_distance_A'] = float(external.min())
    # Optimize a copy only. Energy changes are comparable within a molecule,
    # not across chemically different candidates. Preserve stereo hypothesis.
    relaxed = Chem.AddHs(Chem.Mol(molecule),addCoords=True)
    properties = AllChem.MMFFGetMoleculeProperties(relaxed,mmffVariant='MMFF94s')
    if properties is None:
        result['relaxation_status'] = 'MMFF parameters unavailable'
    else:
        forcefield = AllChem.MMFFGetMoleculeForceField(relaxed,properties)
        before = forcefield.CalcEnergy(); status = 1; iteration_budget = 0
        while status and iteration_budget<3000:
            status = forcefield.Minimize(maxIts=500); iteration_budget += 500
        after = forcefield.CalcEnergy()
        result['relaxation'] = {'converged':status==0,'iteration_budget':iteration_budget,'initial_energy_kcal_mol':before,
            'relaxed_energy_kcal_mol':after,'energy_drop_kcal_mol':before-after,
            'limits':'Unrestrained isolated-molecule relaxation; not binding energy, and coordinates not used to replace the supplied pose.'}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',type=Path,required=True)
    parser.add_argument('--xyz',type=Path,required=True)
    parser.add_argument('--warhead-run',type=Path,required=True)
    parser.add_argument('--e3-run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args(); results = []
    for case in sorted(args.cases.iterdir()):
        if not case.is_dir(): continue
        reference = case/'ak1690_observed.sdf'
        if not reference.exists(): continue
        mol = next(iter(Chem.SDMolSupplier(str(reference))))
        row = {'case':case.name,'kind':'inferred_linker_reference',**assess(mol,case/'stat6_observed.pdb')}
        row['inference_stability'] = inference_stability(args.xyz/(case.name+'.xyz'))
        results.append(row)
    for name,run in [('new_warhead',args.warhead_run),('sequential_e3',args.e3_run)]:
        molecules = [m for m in Chem.SDMolSupplier(str(run/'accepted.sdf')) if m is not None]
        for index,mol in enumerate(molecules):
            results.append({'case':name,'sample':index,'kind':'generated_candidate',**assess(mol,run/'pocket.pdb')})
    report = {'results':results,'limits':'Reconstructed chemistry and transferred protein poses remain hypothetical; no affinity, degradation, or synthesis validation.'}
    with args.output.open('x',encoding='utf-8') as stream: json.dump(report,stream,indent=2)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
