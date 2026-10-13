import unittest
from rdkit import Chem
from assess_exploratory_candidates import assess


def positioned(smiles, coordinates):
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(len(coordinates)); conf.Set3D(True)
    for i, point in enumerate(coordinates):
        conf.SetAtomPosition(i,point)
    mol.AddConformer(conf)
    return mol


class GeometryTests(unittest.TestCase):
    def test_missing_coordinates_rejected(self):
        with self.assertRaises(ValueError):
            assess(Chem.MolFromSmiles('CC'))

    def test_extreme_bond_flagged(self):
        result = assess(positioned('CC',[(0.,0.,0.),(5.,0.,0.)]))
        self.assertEqual(len(result['extreme_bond_lengths']),1)

    def test_nonlocal_overlap_flagged(self):
        result = assess(positioned('CCCC',[(0.,0.,0.),(1.54,0.,0.),(1.54,1.54,0.),(.1,0.,0.)]))
        self.assertGreater(result['nonbonded_overlap_pairs'],0)

    def test_bonded_neighbors_not_counted_as_steric_overlap(self):
        result = assess(positioned('CC',[(0.,0.,0.),(1.54,0.,0.)]))
        self.assertEqual(result['nonbonded_overlap_pairs'],0)
        self.assertFalse(result['lipinski_filter_applied'])


if __name__=='__main__':
    unittest.main()
