"""Construct and audit a transferred CRBN pose; never an experimental ternary pose."""
import argparse
import json
from pathlib import Path
import zipfile

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--source-pdb', type=Path, help='Optional 5YJ1 public R-thalidomide mouse CRBN domain')
    args = parser.parse_args()
    case = args.case
    if not args.source_pdb:
        with zipfile.ZipFile(args.archive) as archive:
            ligand_block = archive.read('data/BTK_CRBN/Y70.pdb').decode()
        protein_block = (case/'unbound_CRBN_4CI3_B.pdb').read_text()
    source_description = 'Archive 4CI3/Y70 reference'
    if args.source_pdb:
        lines = args.source_pdb.read_text().splitlines()
        if not lines or not lines[0].startswith('HEADER') or '5YJ1' not in lines[0]:
            raise ValueError('The optional source route expects the verified public 5YJ1 structure')
        ligand_block = '\n'.join(l for l in lines if l.startswith('HETATM') and l[17:20]=='6EL' and l[21]=='A')+'\nEND\n'
        protein_block = '\n'.join(l for l in lines if l.startswith('ATOM  ') and l[21]=='A')
        source_description = '5YJ1 chain A; R-thalidomide; mouse CRBN binding domain, not human full-length CRBN'
    reference = next(iter(Chem.SDMolSupplier(str(case/'ak1690_observed.sdf'))))
    source = Chem.MolFromPDBBlock(ligand_block, sanitize=False, removeHs=True)
    # Y70 is pomalidomide; bond orders restored from this explicit template.
    template = Chem.MolFromSmiles('O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1' if args.source_pdb
                                else 'O=C1CCC(N2C(=O)c3cccc(N)c3C2=O)C(=O)N1')
    source = AllChem.AssignBondOrdersFromTemplate(template, source)
    Chem.AssignStereochemistryFrom3D(source)
    pairs = []
    for stereo in ('@', '@@'):
        query = Chem.MolFromSmarts(f'O=C1CC[C{stereo}H](N)C(=O)N1')
        pairs.extend((a, b) for a in source.GetSubstructMatches(query, useChirality=True)
                     for b in reference.GetSubstructMatches(query, useChirality=True))
    if not pairs:
        # Absolute stereochemistry must be compatible; do not silently reflect.
        raise ValueError('No stereochemically compatible glutarimide match')
    best = None
    for source_match, target_match in pairs:
        x = source.GetConformer().GetPositions()[list(source_match)]
        y = reference.GetConformer().GetPositions()[list(target_match)]
        u, _, vt = np.linalg.svd((x-x.mean(0)).T@(y-y.mean(0)))
        correction = np.eye(3); correction[-1, -1] = np.linalg.det(u@vt)
        rotation = u@correction@vt
        translation = y.mean(0)-x.mean(0)@rotation
        rmsd = float(np.sqrt(np.mean(np.sum((x@rotation+translation-y)**2, axis=1))))
        if best is None or rmsd < best[0]:
            best = (rmsd, rotation, translation, source_match, target_match)
    rmsd, rotation, translation, a, b = best
    lines = []; protein_coords = []
    for line in protein_block.splitlines():
        if line.startswith('ATOM  '):
            point = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            point = point@rotation+translation
            protein_coords.append(point)
            line = line[:30]+''.join(f'{v:8.3f}' for v in point)+line[54:]
            lines.append(line)
    coords = reference.GetConformer().GetPositions()
    distances = np.linalg.norm(np.array(protein_coords)[:, None, :]-coords[None, :, :], axis=2)
    stat6 = []
    for line in (case/'stat6_observed.pdb').read_text().splitlines():
        if line.startswith('ATOM  '):
            stat6.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    interprotein = np.linalg.norm(np.array(protein_coords)[:, None, :]-np.array(stat6)[None, :, :], axis=2)
    report = {'status': 'unvalidated_transferred_pose', 'source': source_description,
              'matched_atoms': len(a), 'core_rmsd_A': rmsd,
              'source_indices': a, 'reference_indices': b,
              'crbn_ligand_pairs_below_1_5_A': int((distances<1.5).sum()),
              'crbn_stat6_pairs_below_1_5_A': int((interprotein<1.5).sum()),
              'limits': 'Ligand-core alignment from another CRBN complex; not experimentally validated, not DeepTernary prediction.'}
    with (case/'crbn_transferred_hypothesis.pdb').open('x') as handle:
        handle.write('\n'.join(lines)+'\nEND\n')
    with (case/'crbn_pose_audit.json').open('x') as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report))


if __name__ == '__main__':
    main()
