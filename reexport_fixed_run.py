"""Re-export saved diffusion arrays after stereo handling fixes; preserve originals."""
import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).parent/'DecompDiff-main'))
from utils import reconstruct as recon
from utils.transforms import get_atomic_number_from_index
from run_checks import run_record, sha256, export_ready_mol
from validate_fixed_candidate import validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--atom-mode', choices=['basic','add_aromatic','full'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = json.loads((args.run/'report.json').read_text())
    reference_path = args.case/'ak1690_observed.sdf'
    if report['reference_sha256'] != sha256(reference_path):
        raise ValueError('Reference differs from original sampling input')
    reference = next(iter(Chem.SDMolSupplier(str(reference_path))))
    task = next(t for t in json.loads((args.case/'generation_tasks.json').read_text()) if t['role']==report['role'])
    metadata = {'method':'Stereo-corrected export of existing sampling arrays; no new diffusion run',
                'source_report_sha256':sha256(args.run/'report.json'), 'atom_encoding_mode':args.atom_mode,
                'reference_sha256':sha256(reference_path)}
    with run_record(args.output, metadata) as output:
        rows = []
        with Chem.SDWriter(str(output/'candidates.sdf')) as all_writer, Chem.SDWriter(str(output/'accepted.sdf')) as accepted_writer:
            for original in report['results']:
                i = original['sample']; source = args.run/f'sample_{i}_raw.npz'
                with np.load(source, allow_pickle=False) as raw:
                    mol = recon.reconstruct_from_generated_with_bond(raw['positions'].astype(np.float64),
                        get_atomic_number_from_index(torch.tensor(raw['atom_types']), mode=args.atom_mode),
                        raw['bond_index'].tolist(), raw['bond_types'])
                mol.GetConformer().Set3D(True)
                mol = export_ready_mol(mol)
                result = dict(sample=i, **validate(reference,mol,task))
                result['smiles'] = Chem.MolToSmiles(mol)
                result['trajectory_constraints_preserved'] = original['trajectory_constraints_preserved']
                mol.SetProp('sample_id',str(i)); mol.SetProp('role',report['role'])
                mol.SetProp('status','stereo_corrected_export_activity_unvalidated')
                all_writer.write(mol)
                if result['accepted']:
                    accepted_writer.write(mol)
                shutil.copyfile(source,output/source.name)
                rows.append(result)
        corrected = dict(report, results=rows, accepted=sum(r['accepted'] for r in rows),
            accepted_changed=sum(r['accepted'] and r.get('changed',False) for r in rows),
            export_correction=metadata)
        (output/'report.json').write_text(json.dumps(corrected,indent=2))
        print(json.dumps({'accepted':corrected['accepted'],'accepted_changed':corrected['accepted_changed']}))


if __name__=='__main__':
    main()
