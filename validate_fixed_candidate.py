"""Check chemical preservation independently of atom order and coordinates.

This is an acceptance gate, not a constrained generator or activity predictor.
"""
import argparse
import json
from pathlib import Path

from rdkit import Chem


def signature(atom):
    return (atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetIsotope(),
            atom.GetIsAromatic(), atom.GetTotalNumHs())


def validate(reference, candidate, task):
    if reference is None or any(a.GetAtomicNum() == 1 for a in reference.GetAtoms()):
        raise ValueError('Reference must be a readable hydrogen-removed molecule')
    fixed_input = task.get('fixed_indices', [])
    generated = task.get('generate_indices', [])
    indices = fixed_input + generated
    if (not fixed_input or not generated or
            any(type(i) is not int for i in indices) or
            len(indices) != len(set(indices)) or
            set(indices) != set(range(reference.GetNumAtoms()))):
        raise ValueError('Fixed/generated indices must form a nonempty disjoint complete partition')
    if candidate is None:
        return {'accepted': False, 'reason': 'unreadable molecule'}
    try:
        candidate = Chem.RemoveHs(candidate)
        Chem.SanitizeMol(candidate)
    except (ValueError, RuntimeError):
        return {'accepted': False, 'reason': 'invalid chemistry'}
    if len(Chem.GetMolFrags(candidate)) != 1:
        return {'accepted': False, 'reason': 'disconnected molecule'}
    fixed = sorted(task['fixed_indices'])
    fixed_set = set(fixed)
    editable = Chem.RWMol(reference)
    for index in reversed(range(reference.GetNumAtoms())):
        if index not in fixed_set:
            editable.RemoveAtom(index)
    query = editable.GetMol()
    # Removing the generated side changes valence; retain the original H counts
    # for explicit comparison against the final, attached candidate below.
    Chem.SanitizeMol(query)
    boundaries = {}
    for bond in reference.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (a in fixed_set) != (b in fixed_set):
            inside = a if a in fixed_set else b
            boundaries.setdefault(inside, []).append(bond.GetBondType())
    for match in candidate.GetSubstructMatches(query, useChirality=True, maxMatches=10000):
        mapped = set(match)
        valid = True
        for source, target in zip(fixed, match):
            original = reference.GetAtomWithIdx(source)
            atom = candidate.GetAtomWithIdx(target)
            if signature(original) != signature(atom):
                valid = False
                break
            external = [bond.GetBondType() for bond in atom.GetBonds()
                        if bond.GetOtherAtomIdx(target) not in mapped]
            if sorted(map(str, external)) != sorted(map(str, boundaries.get(source, []))):
                valid = False
                break
        # An induced subgraph must not contain extra bonds absent in the query.
        expected = query.GetNumBonds()
        actual = sum(b.GetBeginAtomIdx() in mapped and b.GetEndAtomIdx() in mapped
                     for b in candidate.GetBonds())
        if valid and actual == expected:
            return {'accepted': True, 'fixed_atom_mapping': dict(zip(fixed, match)),
                    'fixed_atoms': len(fixed), 'candidate_atoms': candidate.GetNumAtoms(),
                    'changed': Chem.MolToSmiles(reference) != Chem.MolToSmiles(candidate)}
    return {'accepted': False, 'reason': 'fixed chemistry, stereochemistry or attachment changed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--candidates', type=Path, required=True)
    parser.add_argument('--role', choices=['warhead', 'e3_ligand'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    references = list(Chem.SDMolSupplier(str(args.case/'ak1690_observed.sdf')))
    if len(references) != 1 or references[0] is None:
        raise ValueError('Case must contain exactly one valid reference molecule')
    reference = references[0]
    tasks = json.loads((args.case/'generation_tasks.json').read_text())
    task = next(t for t in tasks if t['role'] == args.role)
    results = [dict(index=i, **validate(reference, molecule, task))
               for i, molecule in enumerate(Chem.SDMolSupplier(str(args.candidates)))]
    if not results:
        raise ValueError('Candidate SDF contains no records')
    report = {'role': args.role, 'results': results,
              'accepted': sum(r['accepted'] for r in results),
              'accepted_changed': sum(r['accepted'] and r.get('changed', False) for r in results),
              'total': len(results), 'limits': 'Chemical preservation only; no binding or activity validation.'}
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps({'accepted': report['accepted'], 'total': len(results)}))


if __name__ == '__main__':
    main()
