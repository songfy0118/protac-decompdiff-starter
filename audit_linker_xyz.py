"""Audit generated XYZs; inferred bonds are hypotheses, never source bond data."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path, PurePosixPath

import numpy as np
from rdkit import Chem


def read_xyz(text):
    lines = text.splitlines()
    count = int(lines[0])
    rows = [line.split() for line in lines[2:] if line.strip()]
    if len(rows) != count or any(len(row) != 4 for row in rows):
        raise ValueError('XYZ count or coordinate columns invalid')
    elements = [row[0] for row in rows]
    positions = np.array([[float(v) for v in row[1:]] for row in rows])
    if not np.isfinite(positions).all():
        raise ValueError('Nonfinite coordinates')
    return elements, positions


def main():
    from openbabel import openbabel as ob
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Use a fresh output directory')
    with zipfile.ZipFile(args.archive) as archive:
        names = [i.filename for i in archive.infolist() if not i.is_dir()]
        if any('..' in PurePosixPath(n).parts or PurePosixPath(n).is_absolute() for n in names):
            raise ValueError('Unsafe archive path')
        texts = {Path(n).name: archive.read(n).decode('utf-8-sig') for n in names if n.endswith('.xyz')}
    fixed_elements, fixed_positions = read_xyz(texts['frag_.xyz'])
    reference = next(iter(Chem.SDMolSupplier(str(args.reference))))
    # Explicit case-specific proposal, checked against every element below.
    original_indices = list(range(46))+list(range(54, 71))
    if fixed_elements != [reference.GetAtomWithIdx(i).GetSymbol() for i in original_indices]:
        raise ValueError('Fixed element sequence does not match the proposed 9BIG mapping')
    original_to_xyz = dict(zip(original_indices, range(63)))
    expected = {tuple(sorted((original_to_xyz[b.GetBeginAtomIdx()], original_to_xyz[b.GetEndAtomIdx()])))
                for b in reference.GetBonds()
                if b.GetBeginAtomIdx() in original_to_xyz and b.GetEndAtomIdx() in original_to_xyz}
    args.output.mkdir(parents=True)
    converter = ob.OBConversion(); converter.SetInAndOutFormats('xyz', 'sdf')
    reports = []
    for name, text in texts.items():
        (args.output/name).write_text(text)
        elements, positions = read_xyz(text)
        distances = np.linalg.norm(positions[:, None, :]-positions[None, :, :], axis=2)
        np.fill_diagonal(distances, np.inf)
        row = {'file': name, 'atoms': len(elements), 'minimum_atom_distance_A': float(distances.min()),
               'overlapping_pairs_below_0_3_A': int(np.triu(distances<0.3, 1).sum())}
        if name not in ('true_.xyz', 'frag_.xyz'):
            row['fixed_prefix_preserved'] = (elements[:63] == fixed_elements and
                                            np.allclose(positions[:63], fixed_positions, atol=1e-7, rtol=0))
            row['generated_linker_atoms'] = len(elements)-63
            molecule = ob.OBMol()
            if not converter.ReadString(molecule, text):
                raise ValueError('OpenBabel XYZ parsing failed')
            block = converter.WriteString(molecule)
            sdf_path = args.output/(Path(name).stem+'_inferred.sdf')
            sdf_path.write_text(block)
            mol = Chem.MolFromMolBlock(block, sanitize=True, removeHs=True)
            inferred = {tuple(sorted((b.GetBeginAtomIdx()-1, b.GetEndAtomIdx()-1))) for b in ob.OBMolBondIter(molecule)}
            fixed_inferred = {pair for pair in inferred if max(pair)<63}
            row['missing_fixed_edges_under_inference'] = sorted(expected-fixed_inferred)
            row['extra_fixed_edges_under_inference'] = sorted(fixed_inferred-expected)
            row['inferred_boundary_edges'] = sorted(pair for pair in inferred if min(pair)<63<=max(pair))
            row['rdkit_sanitized'] = mol is not None
            graph = {i: set() for i in range(len(elements))}
            for a, b in inferred:
                graph[a].add(b); graph[b].add(a)
            seen = {0}; queue = [0]
            while queue:
                current = queue.pop()
                for neighbor in graph[current]-seen:
                    seen.add(neighbor); queue.append(neighbor)
            row['connected_under_geometric_inference'] = len(seen) == len(elements)
            row['expected_attachment_sites_preserved_under_inference'] = sorted(
                min(pair) for pair in inferred if min(pair)<63<=max(pair)) == [44, 46]
            hybrid = Chem.RWMol(Chem.MolFromMolBlock(block, sanitize=False, removeHs=False))
            # XYZ is heavy-atom-only; OpenBabel's inferred explicit valence flags
            # must not suppress the H counts recalculated after template repair.
            for atom in hybrid.GetAtoms():
                atom.SetNoImplicit(False)
                atom.SetNumExplicitHs(0)
            for source, target in original_to_xyz.items():
                atom = reference.GetAtomWithIdx(source)
                hybrid.GetAtomWithIdx(target).SetFormalCharge(atom.GetFormalCharge())
            for bond in reference.GetBonds():
                a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if a in original_to_xyz and b in original_to_xyz:
                    observed = hybrid.GetBondBetweenAtoms(original_to_xyz[a], original_to_xyz[b])
                    if observed is not None:
                        observed.SetBondType(bond.GetBondType())
                        observed.SetIsAromatic(bond.GetIsAromatic())
            try:
                hybrid = hybrid.GetMol(); hybrid.GetConformer().Set3D(True); Chem.SanitizeMol(hybrid)
                Chem.AssignStereochemistryFrom3D(hybrid)
                row['template_restored_sanitized'] = True
                with Chem.SDWriter(str(args.output/(Path(name).stem+'_template_restored_unverified.sdf'))) as writer:
                    writer.write(hybrid)
                row['template_restored_smiles'] = Chem.MolToSmiles(hybrid)
            except (ValueError, RuntimeError):
                row['template_restored_sanitized'] = False
            if mol is not None:
                row['inferred_smiles'] = Chem.MolToSmiles(mol)
            row['usable_as_verified_model_input'] = False
            row['usable_as_exploratory_input'] = bool(row['template_restored_sanitized'] and
                row['connected_under_geometric_inference'] and row['fixed_prefix_preserved'] and
                row['expected_attachment_sites_preserved_under_inference'] and
                row['overlapping_pairs_below_0_3_A'] == 0 and not row['missing_fixed_edges_under_inference'] and
                not row['extra_fixed_edges_under_inference'])
        reports.append(row)
    report = {'archive_sha256': hashlib.sha256(args.archive.read_bytes()).hexdigest(),
              'fixed_atom_reference_indices': original_indices, 'fixed_atoms': 63,
              'mapping_status': 'Element-sequence proposal, geometry and stereochemistry require confirmation',
              'results': reports,
              'limits': 'XYZ contains no bonds, bond orders, formal charges or explicit stereochemistry; inferred SDF is unverified.'}
    (args.output/'xyz_audit.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
