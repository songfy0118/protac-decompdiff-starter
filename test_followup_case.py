import hashlib
import json
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from make_followup_case import main


class FollowupTests(unittest.TestCase):
    def prepare(self, root, changed=True):
        parent = root/'parent'; run = root/'run'
        parent.mkdir(); run.mkdir()
        reference = Chem.MolFromSmiles('F[C@H](C)COCCN')
        candidate = Chem.MolFromSmiles('F[C@H](C)COCCO' if changed else 'F[C@H](C)COCCN')
        for molecule in (reference, candidate):
            AllChem.Compute2DCoords(molecule)
        with Chem.SDWriter(str(parent/'ak1690_observed.sdf')) as writer:
            writer.write(reference)
        tasks = [{'role': 'warhead', 'fixed_indices': list(range(5)), 'generate_indices': list(range(5, 8))},
                 {'role': 'e3_ligand', 'fixed_indices': list(range(5, 8)), 'generate_indices': list(range(5))}]
        (parent/'generation_tasks.json').write_text(json.dumps(tasks))
        reference_hash = hashlib.sha256((parent/'ak1690_observed.sdf').read_bytes()).hexdigest()
        (run/'report.json').write_text(json.dumps({'role': 'warhead', 'reference_sha256': reference_hash}))
        candidate.SetProp('sample_id', '0')
        with Chem.SDWriter(str(run/'accepted.sdf')) as writer:
            writer.write(candidate)
        np.savez(run/'sample_0_raw.npz', positions=candidate.GetConformer().GetPositions())
        return parent, run, root/'child'

    def invoke(self, parent, run, child):
        with patch.object(sys, 'argv', ['make_followup_case', '--parent-case', str(parent),
                                       '--run', str(run), '--output', str(child)]):
            main()

    def test_valid_generated_reference_has_provenance(self):
        with TemporaryDirectory() as tmp:
            parent, run, child = self.prepare(Path(tmp))
            self.invoke(parent, run, child)
            provenance = json.loads((child/'provenance.json').read_text())
            self.assertTrue(provenance['reference_is_generated'])
            self.assertEqual(provenance['next_role'], 'e3_ligand')

    def test_wrong_parent_hash_rejected(self):
        with TemporaryDirectory() as tmp:
            parent, run, child = self.prepare(Path(tmp))
            (run/'report.json').write_text(json.dumps({'role':'warhead', 'reference_sha256':'wrong'}))
            with self.assertRaises(ValueError):
                self.invoke(parent, run, child)
            self.assertFalse(child.exists())

    def test_reordered_raw_coordinates_rejected(self):
        with TemporaryDirectory() as tmp:
            parent, run, child = self.prepare(Path(tmp))
            with np.load(run/'sample_0_raw.npz') as raw:
                positions = raw['positions'][::-1]
            np.savez(run/'sample_0_raw.npz', positions=positions)
            with self.assertRaises(ValueError):
                self.invoke(parent, run, child)
            self.assertFalse(child.exists())

    def test_unchanged_reference_rejected(self):
        with TemporaryDirectory() as tmp:
            parent, run, child = self.prepare(Path(tmp), changed=False)
            with self.assertRaises(ValueError):
                self.invoke(parent, run, child)


if __name__ == '__main__':
    unittest.main()
