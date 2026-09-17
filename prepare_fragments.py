"""Validate supplied PROTAC tables and preserve attachment-marked fragments."""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import re
from rdkit import Chem, RDLogger
from openpyxl import load_workbook
from run_checks import fresh_output

SOURCES = {
    'e3_ligand': 'e3_ligand.xlsx',
    'warhead': 'warhead.xlsx',
    'protac': 'protac.xlsx',
    'paired': 'protacDB_smiles.csv',
}

def normalize(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError('SMILES must be nonempty text')
    text = re.sub(r'\[R(\d+)\]', lambda m: '[*:'+m[1]+']', smiles)
    mol = Chem.MolFromSmiles(text)
    if mol is None or mol.GetNumAtoms()==0: raise ValueError('SMILES parse/sanitize failed')
    attachments = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if any(a.GetDegree() != 1 for a in attachments):
        raise ValueError('Attachment marker must have exactly one neighbor')
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    # Cap dummy atoms with H to group identical parent structures across exit vectors.
    capped = Chem.Mol(mol)
    for atom in capped.GetAtoms():
        atom.SetAtomMapNum(0)
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1); atom.SetIsotope(0)
    capped = Chem.RemoveHs(capped)
    Chem.SanitizeMol(capped)
    parent = Chem.MolToSmiles(capped, isomericSmiles=True)
    return canonical, parent, len(attachments)

def split_for(parent):
    value = int(hashlib.sha256(parent.encode()).hexdigest()[:8],16) % 100
    return 'train' if value < 80 else 'validation' if value < 90 else 'test'

def write_csv(path, rows, fields):
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=fields); writer.writeheader(); writer.writerows(rows)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--source-dir',type=Path,required=True)
    p.add_argument('--source-map',type=Path,help='JSON mapping e3_ligand/warhead/protac/paired to input filenames')
    p.add_argument('--output-dir',type=Path,default=Path(__file__).parent/'data/processed')
    args=p.parse_args()
    sources=json.loads(args.source_map.read_text(encoding='utf-8-sig')) if args.source_map else SOURCES
    if not isinstance(sources,dict) or set(sources)!=set(SOURCES):
        p.error('source map must define exactly e3_ligand, warhead, protac, paired')
    for filename in sources.values():
        if not isinstance(filename,str) or not (args.source_dir/filename).is_file():
            p.error(f'Input file not found: {filename}')
    out=fresh_output(args.output_dir)
    RDLogger.DisableLog('rdApp.error')
    errors=[]; summary={}; fragments=defaultdict(list); pairs=[]
    for role,filename in sources.items():
        source=args.source_dir/filename
        if source.suffix=='.csv':
            with source.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
        else:
            wb=load_workbook(source,read_only=True,data_only=True)
            values=iter(wb.active.values); headers=next(values)
            rows=[dict(zip(headers,r)) for r in values if any(x is not None for x in r)]; wb.close()
        required={'id_protac','linker_canonical','warhead_canonical','e3_ligand_canonical','smiles_canonical'} if role=='paired' else {'Compound ID','Uniprot','Target','Smiles'}
        if not rows or required-set(rows[0]):
            raise ValueError(f'{filename}: empty table or missing columns {sorted(required-set(rows[0] if rows else {}))}')
        summary[role]={'rows':len(rows),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
        for line,row in enumerate(rows,2):
            fields = [('warhead','warhead_canonical'),('e3_ligand','e3_ligand_canonical'),('protac','smiles_canonical')] if role=='paired' else [(role,'Smiles')]
            pair={'source_row':line,'protac_id':row.get('id_protac')}
            for molecule_role,column in fields:
                try:
                    canonical,parent,markers=normalize(row[column])
                    if role=='paired' and molecule_role != 'protac' and markers != 1:
                        raise ValueError(f'Expected one attachment marker, found {markers}')
                except (ValueError,TypeError,RuntimeError) as exc:
                    errors.append({'source':filename,'row':line,'field':column,'error':str(exc),'smiles':row[column]}); continue
                if role=='paired': pair[molecule_role]=canonical
                if molecule_role in ('warhead','e3_ligand'):
                    fragments[(molecule_role,canonical)].append({
                        'source':filename,'row':line,'compound_id':row.get('Compound ID',row.get('id_protac')),
                        'target':row.get('Target'),'uniprot':row.get('Uniprot'),
                        'parent':parent,'attachment_count':markers})
            if role=='paired':
                pair['linker_raw']=row['linker_canonical']
                pair['status']='valid' if all(k in pair for k in ('warhead','e3_ligand','protac')) else 'quarantined'
                pairs.append(pair)
    records=[]; provenance=[]
    for (role,canonical),origins in sorted(fragments.items()):
        parent=origins[0]['parent']
        record_id=hashlib.sha256((role+'|'+canonical).encode()).hexdigest()[:20]
        provenance.append({'record_id':record_id,'origins':origins})
        records.append({'role':role,'smiles':canonical,'parent_smiles':parent,
                        'attachment_count':origins[0]['attachment_count'],'split':split_for(parent),
                        'source_count':len(origins),'record_id':record_id})
    write_csv(out/'fragments.csv',records,['role','smiles','parent_smiles','attachment_count','split','source_count','record_id'])
    with (out/'provenance.jsonl').open('w',encoding='utf8') as f:
        for item in provenance: f.write(json.dumps(item,ensure_ascii=False)+'\n')
    write_csv(out/'paired_protacs.csv',pairs,['source_row','protac_id','warhead','e3_ligand','protac','linker_raw','status'])
    write_csv(out/'validation_errors.csv',errors,['source','row','field','error','smiles'])
    summary['fragments_by_role']=dict(Counter(r['role'] for r in records))
    summary['splits']=dict(Counter(r['role']+'/'+r['split'] for r in records))
    summary['invalid_fields']=len(errors)
    summary['notes']=['Splits are per parent molecule for separate endpoint models; paired_protacs.csv has no pair-level split.',
                      'SMILES validation does not establish affinity, activity, or 3D pocket alignment.',
                      'Attachment-free table records and attachment-marked CSV records remain distinct.',
                      'No cross-file join by Compound ID is assumed.']
    (out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf8')
    print(json.dumps(summary,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
