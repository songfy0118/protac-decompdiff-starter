"""Behavioral regression for stereo-preserving upstream conformer preparation."""
import ast
import os
from pathlib import Path
import unittest
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit.DistanceGeometry import DoTriangleSmoothing


class StereoEmbeddingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source=Path(os.environ.get('DEEPTERNARY_SOURCE','work/deepternary-upstream/predict.py'))
        if not source.exists():
            raise unittest.SkipTest('Optional upstream DeepTernary checkout unavailable')
        tree=ast.parse(source.read_text())
        selected=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in {'adjust_bounds','get_rdkit_coords_with_fixed_ends'}]
        namespace=dict(Chem=Chem,rdDistGeom=rdDistGeom,DoTriangleSmoothing=DoTriangleSmoothing)
        exec(compile(ast.Module(body=selected,type_ignores=[]),str(source),'exec'),namespace)
        cls.embed=staticmethod(namespace['get_rdkit_coords_with_fixed_ends'])

    def test_multicenter_stereo_survives_multiple_seeds(self):
        molecule=Chem.AddHs(Chem.MolFromSmiles('C[C@H](O)CC[C@@H](N)C(=O)O'))
        self.assertEqual(AllChem.EmbedMolecule(molecule,randomSeed=19),0)
        molecule=Chem.RemoveHs(molecule)
        expected=Chem.MolToSmiles(molecule,isomericSmiles=True)
        for seed in range(6):
            result,coords=self.embed(Chem.Mol(molecule),Chem.MolFromSmiles('CC(O)'),Chem.MolFromSmiles('NC(C)C(=O)O'),seed)
            Chem.AssignStereochemistryFrom3D(result,replaceExistingTags=True)
            self.assertEqual(Chem.MolToSmiles(result,isomericSmiles=True),expected,f'seed {seed}')
            self.assertEqual(coords.shape,(molecule.GetNumAtoms(),3))
            self.assertEqual([a.GetAtomicNum() for a in result.GetAtoms()],[a.GetAtomicNum() for a in molecule.GetAtoms()])

    def test_inverted_embedding_is_rejected(self):
        molecule=Chem.AddHs(Chem.MolFromSmiles('C[C@H](O)CC'))
        self.assertEqual(AllChem.EmbedMolecule(molecule,randomSeed=23),0)
        molecule=Chem.RemoveHs(molecule)
        def inverted_embed(mol,params):
            conf=mol.GetConformer()
            for i,xyz in enumerate(conf.GetPositions()):
                conf.SetAtomPosition(i,(float(xyz[0]),float(xyz[1]),-float(xyz[2])))
            return 0
        class ForcedInversion:
            ETKDGv3=staticmethod(rdDistGeom.ETKDGv3)
            EmbedMolecule=staticmethod(inverted_embed)
        namespace=self.embed.__globals__
        original=namespace['rdDistGeom']
        try:
            namespace['rdDistGeom']=ForcedInversion
            with self.assertRaisesRegex(RuntimeError,'reference stereochemistry'):
                self.embed(molecule,Chem.MolFromSmiles('[Xe]'),Chem.MolFromSmiles('[Xe]'),1)
        finally:
            namespace['rdDistGeom']=original


if __name__=='__main__':
    unittest.main()
