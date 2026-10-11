"""Read-only source inspection and atom-name-based PROTAC case preparation.

Does not sample molecules or infer missing experimental coordinates.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
import zipfile

import gemmi
import numpy as np
from rdkit import Chem


def records(block, category):
    table = block.find_mmcif_category(category)
    keys = [str(t).split('.', 1)[1] for t in table.tags]
    return [dict(zip(keys, [gemmi.cif.as_string(str(v)) for v in row])) for row in table]


def prepare(cif, archive, output):
    if output.exists():
        raise ValueError('Use a new output directory; existing results are preserved')
    block = gemmi.cif.read_file(str(cif)).sole_block()
    atoms = [r for r in records(block, '_chem_comp_atom.') if r['comp_id'] == 'A1AQQ']
    bonds = [r for r in records(block, '_chem_comp_bond.') if r['comp_id'] == 'A1AQQ']
    sites = records(block, '_atom_site.')
    entities = records(block, '_entity.')
    rw = Chem.RWMol()
    indices = {}
    for row in atoms:
        atom = Chem.Atom(row['type_symbol'])
        atom.SetProp('source_atom_name', row['atom_id'])
        indices[row['atom_id']] = rw.AddAtom(atom)
    bond_types = {'SING': Chem.BondType.SINGLE, 'DOUB': Chem.BondType.DOUBLE,
                  'TRIP': Chem.BondType.TRIPLE, 'AROM': Chem.BondType.AROMATIC}
    for row in bonds:
        rw.AddBond(indices[row['atom_id_1']], indices[row['atom_id_2']], bond_types[row['value_order'].upper()])
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    # CCD stereochemical labels are declarations, not inferred from missing coordinates.
    desired = {r['atom_id']: r['pdbx_stereo_config'] for r in atoms if r['pdbx_stereo_config'] in ('R', 'S')}
    for atom in mol.GetAtoms():
        if atom.GetProp('source_atom_name') in desired:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for atom in mol.GetAtoms():
        target = desired.get(atom.GetProp('source_atom_name'))
        if target and atom.HasProp('_CIPCode') and atom.GetProp('_CIPCode') != target:
            atom.InvertChirality()
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for atom in mol.GetAtoms():
        target = desired.get(atom.GetProp('source_atom_name'))
        if target and (not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode') != target):
            raise ValueError('Cannot reproduce CCD stereochemistry')
    names = [a.GetProp('source_atom_name') for a in mol.GetAtoms()]
    observed = defaultdict(dict)
    for row in sites:
        if row['label_comp_id'] != 'A1AQQ' or row['type_symbol'] == 'H':
            continue
        group = (row['pdbx_PDB_model_num'], row['label_asym_id'])
        name = row['label_atom_id']
        occupancy = float(row['occupancy'])
        if name not in observed[group] or occupancy > observed[group][name]['occupancy']:
            observed[group][name] = {'xyz': [float(row[k]) for k in ('Cartn_x','Cartn_y','Cartn_z')], 'occupancy': occupancy}
    if not observed:
        raise ValueError('No observed AK-1690 atoms')
    selected = max(observed, key=lambda key: len(observed[key]))
    coords = observed[selected]
    missing = [n for n in names if n not in coords]
    # For AK-1690 only: explicit seven-carbon heptynyl bridge between O and aryl C.
    query = Chem.MolFromSmarts('[O]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[C]#[C]-[c]')
    matches = mol.GetSubstructMatches(query)
    if len(matches) != 1:
        raise ValueError(f'Expected one AK-1690 linker match, found {len(matches)}')
    match = matches[0]
    linker = set(match[1:-1])
    cuts = [(match[0], match[1]), (match[-2], match[-1])]
    adjacency = {a.GetIdx(): set() for a in mol.GetAtoms()}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (a,b) not in cuts and (b,a) not in cuts:
            adjacency[a].add(b); adjacency[b].add(a)
    components = []
    unseen = set(adjacency)
    while unseen:
        stack = [min(unseen)]; group = set()
        while stack:
            atom = stack.pop()
            if atom in group:
                continue
            group.add(atom); stack.extend(adjacency[atom] - group)
        unseen -= group; components.append(group)
    warhead = next(c for c in components if any(mol.GetAtomWithIdx(i).GetSymbol() == 'P' for i in c))
    e3 = next(c for c in components if c != linker and c != warhead)
    if len(components) != 3 or linker not in components or len(warhead|linker|e3) != len(names):
        raise ValueError('Invalid or incomplete decomposition')
    part = {i: role for role, group in [('warhead',warhead),('linker',linker),('e3_ligand',e3)] for i in group}
    tasks = []
    for role, generated in [('warhead',warhead),('e3_ligand',e3)]:
        tasks.append({'role': role, 'generate_indices': sorted(generated),
                      'fixed_indices': sorted(set(range(len(names))) - generated),
                      'attachment_bonds': [{'indices': [a,b], 'atom_names':[names[a],names[b]]} for a,b in cuts if (a in generated) != (b in generated)],
                      'generation_supported_by_current_runner': False,
                      'fixed_policy': 'preserve atom identity and chemical bonds; allow conformation changes'})
    protein_positions = np.array([[float(r[k]) for k in ('Cartn_x','Cartn_y','Cartn_z')] for r in sites
                                 if r['group_PDB']=='ATOM' and r['type_symbol']!='H' and r['pdbx_PDB_model_num']==selected[0]])
    distances = []
    for name in names:
        if name in coords:
            distances.append(float(np.linalg.norm(protein_positions-np.array(coords[name]['xyz']), axis=1).min()))
    with zipfile.ZipFile(archive) as z:
        files = [i for i in z.infolist() if not i.is_dir() and i.file_size]
        listed = set(i.filename for i in files)
        unsafe = [i.filename for i in z.infolist() if '..' in PurePosixPath(i.filename).parts or PurePosixPath(i.filename).is_absolute() or '\\' in i.filename]
        if unsafe:
            raise ValueError('Unsafe archive paths')
        missing_refs = []
        with z.open('data/Affinity/test_list.csv') as f:
            rows = list(csv.DictReader(f.read().decode('utf-8-sig').splitlines()))
        for row in rows:
            for field,value in row.items():
                if field != 'name' and value and 'data/Affinity/'+value not in listed:
                    missing_refs.append({'sample':row['name'],'field':field,'path':value})
        crbn_files = [n for n in listed if n.startswith('data/BTK_CRBN/') and 'unbound_CRBN' in n]
        audit = {'nonempty_files':len(files), 'extensions':dict(Counter(Path(i.filename).suffix for i in files)),
                 'affinity_samples':len(rows), 'missing_affinity_references':missing_refs,
                 'crbn_reference_files': sorted(crbn_files)}
        output.mkdir(parents=True)
        for name in crbn_files:
            (output/Path(name).name).write_bytes(z.read(name))
    (output/'9BIG.cif').write_bytes(cif.read_bytes())
    protein_structure = gemmi.read_structure(str(cif))
    protein_structure.remove_ligands_and_waters()
    protein_structure.write_pdb(str(output/'stat6_observed.pdb'))
    (output/'ak1690_reference.smiles').write_text(Chem.MolToSmiles(mol)+'\n',encoding='utf-8')
    # Write 3D SDF only when every atom has an experimental coordinate.
    if not missing:
        conf = Chem.Conformer(len(names)); conf.Set3D(True)
        for i,name in enumerate(names): conf.SetAtomPosition(i,coords[name]['xyz'])
        mol.AddConformer(conf)
        with Chem.SDWriter(str(output/'ak1690_observed.sdf')) as writer: writer.write(mol)
    else:
        mol.RemoveAllConformers()
        with Chem.SDWriter(str(output/'ak1690_connectivity.sdf')) as writer: writer.write(mol)
    with (output/'atom_mapping.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f);writer.writerow(['index','atom_name','element','part','observed','x','y','z'])
        for i,name in enumerate(names):
            writer.writerow([i,name,mol.GetAtomWithIdx(i).GetSymbol(),part[i],name in coords,*(coords[name]['xyz'] if name in coords else ['','',''])])
    report = {'source_sha256':hashlib.sha256(cif.read_bytes()).hexdigest(),
              'protein_entities':[r['pdbx_description'] for r in entities if r['type']=='polymer'],
              'crystal_resolution':block.find_value('_refine.ls_d_res_high'),
              'ligand_heavy_atoms':len(names),'observed_heavy_atoms':len(set(names)&set(coords)),
              'missing_atom_names':missing,'part_sizes':dict(Counter(part.values())),
              'minimum_observed_ligand_protein_distance':min(distances),
              'archive':audit,'tasks':tasks,
              'limits':['9BIG is a binary STAT6-ligand complex, not an experimental ternary complex.',
                        'CRBN from BTK example is a starting reference; compatibility and pose alignment unverified.',
                        'Seven-carbon linker boundary is an explicit case-specific assumption.',
                        'No diffusion generation, affinity or degradation validation performed.']}
    (output/'inspection.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    (output/'generation_tasks.json').write_text(json.dumps(tasks,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('ligand_heavy_atoms','observed_heavy_atoms','missing_atom_names','part_sizes')},ensure_ascii=False))
    return report


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--cif',type=Path,required=True)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    prepare(args.cif,args.archive,args.output)
