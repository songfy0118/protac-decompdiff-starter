import unittest
import numpy as np
from rdkit import Chem
from run_checks import export_ready_mol
from audit_linker_xyz import read_xyz
from prepare_xyz_cases import fit


class XYZTests(unittest.TestCase):
    def test_valid_xyz(self):
        elements, positions = read_xyz('2\ncomment\nC 0 1 2\nO 3 4 5\n')
        self.assertEqual(elements, ['C', 'O'])
        np.testing.assert_equal(positions, [[0, 1, 2], [3, 4, 5]])

    def test_count_mismatch(self):
        with self.assertRaises(ValueError):
            read_xyz('2\n\nC 0 0 0\n')

    def test_bad_columns(self):
        with self.assertRaises(ValueError):
            read_xyz('1\n\nC 0 0 0 extra\n')

    def test_nonfinite(self):
        with self.assertRaises(ValueError):
            read_xyz('1\n\nC nan 0 0\n')

    def test_fit_preserves_proper_rotation(self):
        x = np.array([[0.,0,0],[1,0,0],[0,2,0],[0,0,3]])
        rotation = np.array([[0.,1,0],[-1,0,0],[0,0,1]])
        y = x@rotation+np.array([4.,5,6])
        r, t, rmsd = fit(x, y)
        self.assertLess(rmsd, 1e-10)
        np.testing.assert_allclose(x@r+t, y, atol=1e-10)
        self.assertAlmostEqual(np.linalg.det(r), 1.)

    def test_fit_does_not_reflect_stereochemistry(self):
        x = np.array([[0.,0,0],[1,0,0],[0,2,0],[0,0,3]])
        y = x*np.array([-1.,1,1])
        r, t, rmsd = fit(x, y)
        self.assertAlmostEqual(np.linalg.det(r), 1.)
        self.assertGreater(rmsd, .1)

    def test_export_assigns_alkene_stereo_from_3d(self):
        molecule = Chem.MolFromSmiles('CC=CC')
        conf = Chem.Conformer(4); conf.Set3D(True)
        for i, point in enumerate([(-1.,1.,0.),(0.,0.,0.),(1.34,0.,0.),(2.34,-1.,0.)]):
            conf.SetAtomPosition(i, point)
        molecule.AddConformer(conf)
        result = export_ready_mol(molecule)
        self.assertEqual(result.GetBondBetweenAtoms(1,2).GetStereo(), Chem.BondStereo.STEREOE)


if __name__ == '__main__':
    unittest.main()
