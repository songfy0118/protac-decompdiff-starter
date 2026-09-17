import unittest
from prepare_fragments import normalize, split_for, read_table
from input_validation import validate_groups, validate_coordinates
from run_checks import positive_int, fresh_output, export_ready_mol, run_record, load_endpoint_spec
from rdkit import Chem
from rdkit.Chem import AllChem
from tempfile import TemporaryDirectory
from pathlib import Path
import numpy as np
import json
import subprocess
import sys
try:
    import torch
except ImportError:
    torch=None
from verify_outputs import verify_molecules, verify_fragment_tables
from package_project import build_archive, verify_archive
import copy
import hashlib
import zipfile

class FragmentTests(unittest.TestCase):
    def test_excel_source_rows_survive_blank_rows(self):
        from openpyxl import Workbook
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'input.xlsx'
            workbook=Workbook(); sheet=workbook.active
            sheet.append(['Smiles']); sheet.append(['CC']); sheet.append([None]); sheet.append(['CO'])
            workbook.save(path); workbook.close()
            self.assertEqual(read_table(path,{'Smiles'}),[(2,{'Smiles':'CC'}),(4,{'Smiles':'CO'})])

    def test_csv_physical_rows_include_multiline_and_blank_lines(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'input.CSV'
            path.write_text('Smiles,Note\nCC,"first\nsecond"\n\nCO,last\n',encoding='utf8')
            rows=read_table(path,{'Smiles'})
            self.assertEqual([line for line,row in rows],[2,5])
            self.assertEqual(rows[1][1]['Smiles'],'CO')

    def test_ambiguous_and_empty_tables_rejected(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'input.csv'
            for text in ('Smiles,Smiles\nCC,CO\n','Smiles\nCC,CO\n','Smiles\n','Other\nCC\n'):
                path.write_text(text)
                with self.assertRaises(ValueError): read_table(path,{'Smiles'})

    def test_attachment_is_preserved(self):
        canonical,parent,count=normalize('[R1]CCO')
        self.assertIn('[*:1]',canonical)
        self.assertEqual(parent,'CCO')
        self.assertEqual(count,1)

    def test_parent_split_matches_without_attachment(self):
        _,parent,_=normalize('[R2]CCO')
        _,other,_=normalize('OCC')
        self.assertEqual(split_for(parent),split_for(other))

    def test_invalid_or_unconnected_marker_rejected(self):
        for text in ('not_smiles','[R1].CC','C[R1]C'):
            with self.assertRaises(ValueError): normalize(text)

    def test_stereochemistry_preserved(self):
        a=normalize('N[C@@H](C)C(=O)O')[1]
        b=normalize('N[C@H](C)C(=O)O')[1]
        self.assertNotEqual(a,b)

    def test_atom_groups_reject_bad_length_and_degenerate_prior(self):
        for groups,count in [([0,0],3),([0,1,1],3),([1,1],2),([-1,-1],2)]:
            with self.assertRaises(ValueError): validate_groups(groups,count)
        self.assertEqual(validate_groups([-1,-1,0,0],4),1)

    def test_empty_smiles_rejected(self):
        for value in ('', '  ', None):
            with self.assertRaises(ValueError): normalize(value)

    def test_bad_group_types_rejected(self):
        for groups in (None, '00', [True,True], [0.0,0.0]):
            with self.assertRaises(ValueError): validate_groups(groups,2)

    def test_invalid_geometry_rejected(self):
        valid=np.zeros((2,3))
        for coords in (np.empty((0,3)), [[np.nan,0,0]], [[np.inf,0,0]], [0,0,0]):
            with self.assertRaises(ValueError): validate_coordinates(coords,valid)
            with self.assertRaises(ValueError): validate_coordinates(valid,coords)

    def test_nonpositive_counts_rejected(self):
        for n in (0,-1):
            with self.assertRaises(ValueError): positive_int(n)

    def test_existing_output_preserved(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'result'; fresh_output(path)
            (path/'keep.txt').write_text('original')
            with self.assertRaises(FileExistsError): fresh_output(path)
            self.assertEqual((path/'keep.txt').read_text(),'original')

    def test_sdf_report_stereochemistry_matches(self):
        mol=Chem.AddHs(Chem.MolFromSmiles('CC(O)F'))
        self.assertEqual(AllChem.EmbedMolecule(mol,randomSeed=2021),0)
        mol=Chem.RemoveHs(mol)
        exported=export_ready_mol(mol)
        reread=Chem.MolFromMolBlock(Chem.MolToMolBlock(exported))
        self.assertEqual(Chem.MolToSmiles(exported),Chem.MolToSmiles(reread))
        self.assertIn('@',Chem.MolToSmiles(exported))

    def test_report_structure_and_count_mismatches_rejected(self):
        with TemporaryDirectory() as tmp:
            directory=Path(tmp); writer=Chem.SDWriter(str(directory/'m.sdf'))
            writer.write(Chem.MolFromSmiles('CCO')); writer.close()
            for report in [
                {'samples':2,'results':[{'sample':0,'smiles':'CCO','connected':True}]},
                {'samples':1,'results':[{'sample':0,'smiles':'CCN','connected':True}]},
                {'samples':1,'results':[{'sample':0,'smiles':'CCO','connected':False}]}]:
                (directory/'report.json').write_text(json.dumps(report))
                with self.assertRaises(ValueError): verify_molecules(directory,'m.sdf')

    @unittest.skipIf(torch is None, 'PyTorch optional: install inference dependencies to test .pt export')
    def test_export_uses_run_metadata(self):
        with TemporaryDirectory() as tmp:
            directory=Path(tmp); (directory/'case').mkdir()
            torch.save([{'mol':Chem.MolFromSmiles('CCO')}],directory/'case/result.pt')
            manifest={'status':'completed','requested':1,'data_id':17,'steps':987,'checkpoint_sha256':'test_fixture'}
            (directory/'run_manifest.json').write_text(json.dumps(manifest))
            subprocess.run([sys.executable,str(Path(__file__).parent/'export_baseline.py'),
                            '--run-dir',str(directory)],check=True,capture_output=True,text=True)
            report=json.loads((directory/'report.json').read_text())
            self.assertEqual((report['data_id'],report['steps']),(17,987))

    def test_run_record_success(self):
        with TemporaryDirectory() as tmp:
            output=Path(tmp)/'run'
            with run_record(output,{'requested':1}): pass
            report=json.loads((output/'run_manifest.json').read_text())
            self.assertEqual(report['status'],'completed')
            self.assertIn('finished_at',report)

    def test_run_record_failure_and_interrupt_propagate(self):
        with TemporaryDirectory() as tmp:
            for error,status in [(ValueError('fixture failure'),'failed'),(KeyboardInterrupt(),'interrupted')]:
                output=Path(tmp)/status
                with self.assertRaises(type(error)):
                    with run_record(output,{}): raise error
                report=json.loads((output/'run_manifest.json').read_text())
                self.assertEqual(report['status'],status)
                self.assertEqual(report['error_type'],type(error).__name__)

    def test_spec_schema_and_file_validation(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'spec.json'
            for spec in ([],{}, {'role':'bad','protein_pdb':'p','reference_sdf':'s','atom_groups':[]}):
                path.write_text(json.dumps(spec))
                with self.assertRaises(ValueError): load_endpoint_spec(path)
            spec={'role':'test_demo','protein_pdb':'p','reference_sdf':'s','atom_groups':[0,0]}
            path.write_text(json.dumps(spec))
            with self.assertRaises(FileNotFoundError): load_endpoint_spec(path)
            (Path(tmp)/'p').touch(); (Path(tmp)/'s').touch()
            self.assertEqual(load_endpoint_spec(path),spec)

    def test_output_file_is_rejected(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'file'; path.write_text('keep')
            with self.assertRaises(ValueError): fresh_output(path)

    def test_fragment_metadata_corruption_rejected(self):
        canonical,parent,count=normalize('[R1]CCO')
        identity=hashlib.sha256(('warhead|'+canonical).encode()).hexdigest()[:20]
        row={'role':'warhead','smiles':canonical,'parent_smiles':parent,'attachment_count':str(count),
             'record_id':identity,'source_count':'1','split':split_for(parent)}
        provenance=[{'record_id':identity,'origins':[{'parent':parent,'attachment_count':count}]}]
        summary={'fragments_by_role':{'warhead':1},'splits':{'warhead/'+row['split']:1}}
        self.assertEqual(verify_fragment_tables([row],provenance,summary),1)
        for key,value in [('source_count','2'),('split','wrong'),('record_id','wrong'),('parent_smiles','C')]:
            altered=dict(row,**{key:value})
            with self.assertRaises(ValueError): verify_fragment_tables([altered],provenance,summary)
        with self.assertRaises(ValueError): verify_fragment_tables([row],provenance*2,summary)
        altered=copy.deepcopy(summary); altered['fragments_by_role']['warhead']=2
        with self.assertRaises(ValueError): verify_fragment_tables([row],provenance,altered)

    def test_archive_excludes_local_assets_and_secrets(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp)/'project'; root.mkdir()
            (root/'README.md').write_text('public')
            for name in ('data/private.csv','assets/model.pt','runs/result.json','.venv/pyvenv.cfg','.env','.env.local','cache.lmdb-lock'):
                path=root/name; path.parent.mkdir(parents=True,exist_ok=True); path.write_text('private')
            output=Path(tmp)/'bundle.zip'
            self.assertEqual(build_archive(root,output),1)
            with zipfile.ZipFile(output) as archive:
                self.assertEqual(set(archive.namelist()),{'project/README.md','PACKAGE_MANIFEST.json'})

    def test_archive_checksum_rejects_tampering(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp)/'project'; root.mkdir(); (root/'README.md').write_text('original')
            output=Path(tmp)/'bundle.zip'
            self.assertEqual(build_archive(root,output),1)
            self.assertEqual(verify_archive(output),1)
            with zipfile.ZipFile(output) as z: contents={n:z.read(n) for n in z.namelist()}
            contents['project/README.md']=b'changed'
            with zipfile.ZipFile(output,'w') as z:
                for name,data in contents.items(): z.writestr(name,data)
            with self.assertRaises(ValueError): verify_archive(output)

    def test_run_manifest_mismatch_and_unfinished_run_rejected(self):
        with TemporaryDirectory() as tmp:
            directory=Path(tmp); writer=Chem.SDWriter(str(directory/'m.sdf'))
            writer.write(Chem.MolFromSmiles('CCO')); writer.close()
            manifest={'status':'completed','requested':1,'steps':1000,'seed':1,'checkpoint_sha256':'fixture'}
            report={'samples':1,'steps':1000,'seed':2,'checkpoint_sha256':'fixture',
                    'results':[{'sample':0,'smiles':'CCO','connected':True}]}
            (directory/'run_manifest.json').write_text(json.dumps(manifest))
            (directory/'report.json').write_text(json.dumps(report))
            with self.assertRaises(ValueError): verify_molecules(directory,'m.sdf')
            report['seed']=1; manifest['status']='running'
            (directory/'run_manifest.json').write_text(json.dumps(manifest))
            (directory/'report.json').write_text(json.dumps(report))
            with self.assertRaises(ValueError): verify_molecules(directory,'m.sdf')

if __name__=='__main__': unittest.main()
