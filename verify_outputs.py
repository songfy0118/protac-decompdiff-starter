"""Verify exported records and molecules without executing serialized models."""
import csv
import json
from collections import defaultdict, Counter
from pathlib import Path
from rdkit import Chem
import argparse
import hashlib
import math
from prepare_fragments import normalize, split_for

def verify_fragment_tables(rows, provenance, summary):
    identifiers=[r['record_id'] for r in rows]
    origin_ids=[r['record_id'] for r in provenance]
    if len(set(identifiers))!=len(identifiers) or len(set(origin_ids))!=len(origin_ids):
        raise ValueError('Duplicate fragment or provenance ID')
    if set(identifiers)!=set(origin_ids): raise ValueError('Missing provenance link')
    origins={r['record_id']:r['origins'] for r in provenance}
    groups=defaultdict(set)
    for row in rows:
        if row['role'] not in ('warhead','e3_ligand'): raise ValueError('Unknown fragment role')
        canonical,parent,markers=normalize(row['smiles'])
        if canonical!=row['smiles'] or parent!=row['parent_smiles'] or markers!=int(row['attachment_count']):
            raise ValueError('Canonical structure or attachment metadata mismatch')
        identity=hashlib.sha256((row['role']+'|'+canonical).encode()).hexdigest()[:20]
        if identity!=row['record_id']: raise ValueError('Incorrect fragment record ID')
        if row['split']!=split_for(parent): raise ValueError('Incorrect deterministic split')
        linked=origins[identity]
        if not linked or len(linked)!=int(row['source_count']): raise ValueError('Source count mismatch')
        if any(p['parent']!=parent or p['attachment_count']!=markers for p in linked):
            raise ValueError('Provenance chemistry metadata mismatch')
        groups[parent].add(row['split'])
    if any(len(s)!=1 for s in groups.values()): raise ValueError('Parent structure leaks across splits')
    if dict(Counter(r['role'] for r in rows))!=summary['fragments_by_role']:
        raise ValueError('Summary role counts mismatch')
    if dict(Counter(r['role']+'/'+r['split'] for r in rows))!=summary['splits']:
        raise ValueError('Summary split counts mismatch')
    return len(groups)

def verify_molecules(directory, sdf_name):
    report=json.loads((directory/'report.json').read_text(encoding='utf8'))
    rows=report['results']
    if len(rows)!=report.get('requested',report.get('samples')):
        raise ValueError('Report sample count mismatch')
    if [r.get('sample') for r in rows]!=list(range(len(rows))):
        raise ValueError('Sample IDs must be unique and ordered')
    manifest_path=directory/'run_manifest.json'
    if manifest_path.exists():
        manifest=json.loads(manifest_path.read_text(encoding='utf8'))
        if manifest['status']!='completed': raise ValueError('Run has not completed')
        if manifest['requested']!=len(rows): raise ValueError('Manifest sample count mismatch')
        for key in ('steps','checkpoint_sha256'):
            if manifest[key]!=report.get(key): raise ValueError(f'Manifest {key} mismatch')
        for key in ('seed','data_id','role','protein_sha256','reference_sha256'):
            if key in manifest and manifest[key]!=report.get(key):
                raise ValueError(f'Manifest {key} mismatch')
    for row in rows:
        if row.get('smiles') is None and row.get('connected',False):
            raise ValueError('Failed reconstruction cannot be connected')
    valid=[r for r in rows if r.get('smiles') is not None]
    sdf=directory/sdf_name
    mols=list(Chem.SDMolSupplier(str(sdf))) if sdf.read_text().strip() else []
    if len(mols)!=len(valid) or any(m is None for m in mols):
        raise ValueError('SDF count/readability differs from report')
    if 'reconstructed' in report and report['reconstructed']!=len(valid):
        raise ValueError('Reconstructed count differs from report')
    for mol,row in zip(mols,valid):
        if mol.HasProp('sample_id') and mol.GetProp('sample_id')!=str(row['sample']):
            raise ValueError('SDF sample ID differs from report')
        if not all(math.isfinite(float(x)) for xyz in mol.GetConformer().GetPositions() for x in xyz):
            raise ValueError('SDF has nonfinite coordinates')
        if Chem.MolToSmiles(mol)!=row['smiles']:
            raise ValueError('SDF structure differs from report SMILES')
        if (len(Chem.GetMolFrags(mol))==1)!=row['connected']:
            raise ValueError('SDF connectivity differs from report')
    return {'readable_molecules':len(mols),'connected':sum(len(Chem.GetMolFrags(m))==1 for m in mols)}

def main():
    root=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser()
    parser.add_argument('--baseline',type=Path,default=root/'sample_results/baseline')
    parser.add_argument('--endpoint',type=Path,default=root/'sample_results/endpoint')
    parser.add_argument('--data-dir',type=Path,help='Optional processed tables; omit to check shipped molecule results only')
    parser.add_argument('--report',type=Path,default=root/'runs/verification.json')
    args=parser.parse_args()
    rows=list(csv.DictReader((args.data_dir/'fragments.csv').open(encoding='utf-8-sig'))) if args.data_dir else []
    parent_count=None; pairs=[]
    if args.data_dir:
        provenance=[json.loads(line) for line in (args.data_dir/'provenance.jsonl').open(encoding='utf8')]
        summary=json.loads((args.data_dir/'summary.json').read_text(encoding='utf8'))
        parent_count=verify_fragment_tables(rows,provenance,summary)
        pairs=list(csv.DictReader((args.data_dir/'paired_protacs.csv').open(encoding='utf-8-sig')))
        if len(pairs)!=summary['paired']['rows']: raise ValueError('Pair count differs from summary')
        errors=list(csv.DictReader((args.data_dir/'validation_errors.csv').open(encoding='utf-8-sig')))
        if len(errors)!=summary['invalid_fields']: raise ValueError('Error count differs from summary')
    quarantined=0
    for row in pairs:
        expected='valid' if all(row[k] for k in ('warhead','e3_ligand','protac')) else 'quarantined'
        if row['status']!=expected: raise ValueError('Incorrect pair quarantine status')
        quarantined+=expected=='quarantined'
    mol_checks={}
    for name,directory,file in [('baseline',args.baseline,'generated.sdf'),('endpoint_demo',args.endpoint,'candidates.sdf')]:
        mol_checks[name]=verify_molecules(directory,file)
    report={'tables_checked':args.data_dir is not None,'fragment_rows':len(rows) if args.data_dir else None,'unique_parent_groups':parent_count,
            'cross_split_parent_overlap':0 if args.data_dir else None,'pair_rows':len(pairs) if args.data_dir else None,'quarantined_pairs':quarantined if args.data_dir else None,
            'provenance_links_valid':True if args.data_dir else None,'molecule_exports':mol_checks}
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report,indent=2))

if __name__ == '__main__': main()
