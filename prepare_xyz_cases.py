"""Prepare explicitly hypothetical cases from template-restored XYZ candidates."""
import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem


def fit(x, y):
    u, _, vt = np.linalg.svd((x-x.mean(0)).T@(y-y.mean(0)))
    d = np.eye(3); d[-1, -1] = np.linalg.det(u@vt)
    rotation = u@d@vt; translation = y.mean(0)-x.mean(0)@rotation
    rmsd = float(np.sqrt(np.mean(np.sum((x@rotation+translation-y)**2, axis=1))))
    return rotation, translation, rmsd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--parent-case', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Use a new output directory')
    audit = json.loads((args.audit/'xyz_audit.json').read_text())
    reference = next(iter(Chem.SDMolSupplier(str(args.parent_case/'ak1690_observed.sdf'))))
    mapping = audit['fixed_atom_reference_indices']
    warhead = [i for i in range(46) if i != 2]
    e3 = [2]+list(range(46, 63))
    args.output.mkdir(parents=True)
    summary = []
    for row in audit['results']:
        if not row.get('usable_as_exploratory_input'):
            continue
        stem = Path(row['file']).stem
        source = args.audit/(stem+'_template_restored_unverified.sdf')
        mol = next(iter(Chem.SDMolSupplier(str(source))))
        output = args.output/stem; output.mkdir()
        (output/'ak1690_observed.sdf').write_bytes(source.read_bytes())
        tasks = []
        for role, generated in [('warhead', warhead), ('e3_ligand', e3)]:
            fixed = sorted(set(range(mol.GetNumAtoms()))-set(generated))
            boundaries = [{'indices':[b.GetBeginAtomIdx(), b.GetEndAtomIdx()]}
                for b in mol.GetBonds() if (b.GetBeginAtomIdx() in generated)!=(b.GetEndAtomIdx() in generated)]
            tasks.append({'role':role, 'generate_indices':generated, 'fixed_indices':fixed,
                          'attachment_bonds':boundaries, 'reference_chemistry_inferred':True})
        (output/'generation_tasks.json').write_text(json.dumps(tasks, indent=2))
        r, t, rmsd = fit(reference.GetConformer().GetPositions()[np.array(mapping)[warhead]],
                         mol.GetConformer().GetPositions()[warhead])
        protein_positions = []; lines = []
        for line in (args.parent_case/'stat6_observed.pdb').read_text().splitlines():
            if line.startswith('ATOM  '):
                pos = np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])@r+t
                protein_positions.append(pos)
                line = line[:30]+''.join(f'{v:8.3f}' for v in pos)+line[54:]
                lines.append(line)
        (output/'stat6_observed.pdb').write_text('\n'.join(lines)+'\nEND\n')
        distances = np.linalg.norm(mol.GetConformer().GetPositions()[:,None,:]-np.array(protein_positions)[None,:,:],axis=2)
        provenance = {'reference_is_generated':True, 'source_xyz':row['file'],
            'archive_sha256':audit['archive_sha256'], 'bond_source':'Known fixed bonds from 9BIG; linker bonds inferred from XYZ',
            'pose_source':'Rigid transfer of STAT6 using the warhead atoms; unvalidated',
            'warhead_alignment_RMSD_A':rmsd, 'pairs_below_1_5_A':int((distances<1.5).sum()),
            'minimum_protein_ligand_distance_A':float(distances.min()),
            'limits':'Exploratory reconstructed chemistry and transferred binding pose; not a verified experimental complex.'}
        (output/'provenance.json').write_text(json.dumps(provenance,indent=2))
        summary.append(provenance)
    (args.output/'case_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
