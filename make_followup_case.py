"""Use an accepted generated endpoint as the fixed endpoint of the next stage."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from validate_fixed_candidate import validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-case', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--sample', type=int, default=0)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Use a new followup case directory')
    reference_path = args.parent_case/'ak1690_observed.sdf'
    reference = next(iter(Chem.SDMolSupplier(str(reference_path))))
    report = json.loads((args.run/'report.json').read_text())
    if report['reference_sha256'] != hashlib.sha256(reference_path.read_bytes()).hexdigest():
        raise ValueError('Parent reference differs from the reference used during sampling')
    tasks = json.loads((args.parent_case/'generation_tasks.json').read_text())
    source_task = next(t for t in tasks if t['role'] == report['role'])
    candidates = [m for m in Chem.SDMolSupplier(str(args.run/'accepted.sdf'))
                  if m is not None and m.GetProp('sample_id') == str(args.sample)]
    if len(candidates) != 1:
        raise ValueError('Select exactly one accepted sample')
    candidate = candidates[0]
    result = validate(reference, candidate, source_task)
    if not result['accepted'] or not result['changed']:
        raise ValueError('Followup stage requires a changed candidate passing the source constraints')
    # Index reuse is allowed only after checking the exported SDF preserves raw
    # sampled node coordinates/order. Do not infer identity from SMILES ordering.
    with np.load(args.run/f'sample_{args.sample}_raw.npz', allow_pickle=False) as raw:
        positions = raw['positions']
    if candidate.GetNumAtoms() != reference.GetNumAtoms() or positions.shape != (reference.GetNumAtoms(), 3):
        raise ValueError('Node-count changes require a new atom mapping')
    if not np.allclose(candidate.GetConformer().GetPositions(), positions, atol=0.001, rtol=0):
        raise ValueError('Exported atom order differs from raw sampling; remap before continuing')
    opposite = 'e3_ligand' if report['role'] == 'warhead' else 'warhead'
    next_tasks = [t for t in tasks if t['role'] == opposite]
    if len(next_tasks) != 1:
        raise ValueError('Parent case must contain exactly one opposite endpoint task')
    validate(candidate, candidate, next_tasks[0])
    args.output.mkdir(parents=True)
    with Chem.SDWriter(str(args.output/'ak1690_observed.sdf')) as writer:
        writer.write(candidate)
    for name in ['stat6_observed.pdb', 'crbn_transferred_hypothesis.pdb', 'crbn_pose_audit.json']:
        source = args.parent_case/name
        if source.exists():
            (args.output/name).write_bytes(source.read_bytes())
    (args.output/'generation_tasks.json').write_text(json.dumps(tasks, indent=2))
    provenance = {'reference_is_generated': True, 'next_role': opposite,
                  'parent_reference_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(),
                  'source_report_sha256': hashlib.sha256((args.run/'report.json').read_bytes()).hexdigest(),
                  'source_role': report['role'], 'source_sample': args.sample,
                  'limits': 'Generated reference, not an experimental AK-1690 structure; protein pose limitations carry over.'}
    if (args.parent_case/'provenance.json').exists():
        provenance['parent_reference_provenance'] = json.loads((args.parent_case/'provenance.json').read_text())
    (args.output/'provenance.json').write_text(json.dumps(provenance, indent=2))
    print(json.dumps(provenance))


if __name__ == '__main__':
    main()
