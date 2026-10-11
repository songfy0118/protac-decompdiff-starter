import json
import os
import unittest
from pathlib import Path

from rdkit import Chem
from validate_fixed_candidate import validate


@unittest.skipUnless(os.environ.get('PROTAC_CASE'), 'Optional local 9BIG case; no private data required')
class FixedCandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        case = Path(os.environ['PROTAC_CASE'])
        cls.reference = next(iter(Chem.SDMolSupplier(str(case/'ak1690_observed.sdf'))))
        cls.tasks = json.loads((case/'generation_tasks.json').read_text())

    def test_reference_and_reordering(self):
        reordered = Chem.RenumberAtoms(self.reference, list(reversed(range(71))))
        for task in self.tasks:
            self.assertTrue(validate(self.reference, reordered, task)['accepted'])
            self.assertFalse(validate(self.reference, reordered, task)['changed'])

    def test_conformation_can_change(self):
        candidate = Chem.Mol(self.reference)
        candidate.RemoveAllConformers()
        for task in self.tasks:
            self.assertTrue(validate(self.reference, candidate, task)['accepted'])

    def test_fixed_atom_change_rejected(self):
        for task in self.tasks:
            candidate = Chem.RWMol(self.reference)
            index = next(i for i in task['fixed_indices']
                         if candidate.GetAtomWithIdx(i).GetSymbol() == 'F') if task['role']=='e3_ligand' else 70
            candidate.GetAtomWithIdx(index).SetIsotope(13)
            self.assertFalse(validate(self.reference, candidate.GetMol(), task)['accepted'])

    def test_attachment_removed_rejected(self):
        for task in self.tasks:
            candidate = Chem.RWMol(self.reference)
            a, b = task['attachment_bonds'][0]['indices']
            candidate.RemoveBond(a, b)
            Chem.SanitizeMol(candidate)
            self.assertFalse(validate(self.reference, candidate.GetMol(), task)['accepted'])

    def test_editable_isotope_change_accepted(self):
        for task in self.tasks:
            candidate = Chem.Mol(self.reference)
            candidate.GetAtomWithIdx(task['generate_indices'][0]).SetIsotope(13)
            result = validate(self.reference, candidate, task)
            self.assertTrue(result['accepted'])
            self.assertTrue(result['changed'])

    def test_fixed_stereochemistry_change_rejected(self):
        for task in self.tasks:
            candidate = Chem.Mol(self.reference)
            center = next(i for i in task['fixed_indices']
                          if candidate.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)
            candidate.GetAtomWithIdx(center).InvertChirality()
            self.assertFalse(validate(self.reference, candidate, task)['accepted'])

    def test_warhead_substitution_allowed_only_for_warhead_task(self):
        candidate = Chem.Mol(self.reference)
        index = next(a.GetIdx() for a in candidate.GetAtoms() if a.GetSymbol() == 'F')
        candidate.GetAtomWithIdx(index).SetAtomicNum(17)
        Chem.SanitizeMol(candidate)
        for task in self.tasks:
            self.assertEqual(validate(self.reference, candidate, task)['accepted'], task['role']=='warhead')


class PortableTests(unittest.TestCase):
    def setUp(self):
        # Synthetic chemistry, independently runnable without the research archive.
        self.reference = Chem.MolFromSmiles('F[C@H](C)COCCN')
        self.task = {'fixed_indices': list(range(5)), 'generate_indices': list(range(5, 8))}

    def test_reference_and_reordering(self):
        candidate = Chem.RenumberAtoms(self.reference, list(reversed(range(8))))
        self.assertTrue(validate(self.reference, candidate, self.task)['accepted'])

    def test_changed_generated_end(self):
        candidate = Chem.MolFromSmiles('F[C@H](C)COCCO')
        self.assertTrue(validate(self.reference, candidate, self.task)['accepted'])

    def test_changed_fixed_end(self):
        candidate = Chem.MolFromSmiles('Cl[C@H](C)COCCN')
        self.assertFalse(validate(self.reference, candidate, self.task)['accepted'])

    def test_stereochemistry(self):
        candidate = Chem.MolFromSmiles('F[C@@H](C)COCCN')
        self.assertFalse(validate(self.reference, candidate, self.task)['accepted'])

    def test_disconnected(self):
        candidate = Chem.MolFromSmiles('F[C@H](C)COCCN.C')
        self.assertFalse(validate(self.reference, candidate, self.task)['accepted'])

    def test_invalid_task(self):
        for fixed, generated in [([], list(range(8))), ([0, 0], list(range(1, 8))),
                                 ([0], [2]), ([False], list(range(1, 8)))]:
            with self.assertRaises(ValueError):
                validate(self.reference, self.reference,
                         {'fixed_indices': fixed, 'generate_indices': generated})

    def test_unreadable_candidate(self):
        self.assertFalse(validate(self.reference, None, self.task)['accepted'])


if __name__ == '__main__':
    unittest.main()
