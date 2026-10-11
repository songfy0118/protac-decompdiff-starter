"""Sampler invariants without a checkpoint or research input."""
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

try:
    import torch
    sys.path.insert(0, str(Path(__file__).parent/'DecompDiff-main'))
    from models.decompdiff import DecompScorePosNet3D
    from run_fixed_endpoint import constraints
except ImportError:
    torch = None


@unittest.skipIf(torch is None, 'Requires the isolated model environment')
class ConstraintTests(unittest.TestCase):
    def test_exit_bond_and_absent_cross_bonds(self):
        data = SimpleNamespace(ligand_pos=torch.zeros(3, 3),
            ligand_atom_feature_full=torch.tensor([1, 2, 3]),
            ligand_bond_index=torch.tensor([[0, 1, 1, 2], [1, 0, 2, 1]]),
            ligand_bond_type=torch.tensor([1, 1, 2, 2]),
            ligand_fc_bond_index=torch.tensor([[0, 0, 1, 1, 2, 2], [1, 2, 0, 2, 0, 1]]))
        types, bonds = constraints(data, [0])
        self.assertEqual(types.tolist(), [1, -1, -1])
        self.assertEqual(bonds.tolist(), [1, 0, 1, -1, 0, -1])

    def test_missing_exit_bond_rejected(self):
        data = SimpleNamespace(ligand_pos=torch.zeros(2, 3),
            ligand_atom_feature_full=torch.tensor([1, 2]),
            ligand_bond_index=torch.tensor([[0, 1], [1, 0]]),
            ligand_bond_type=torch.tensor([1, 1]),
            ligand_fc_bond_index=torch.empty(2, 0, dtype=torch.long))
        with self.assertRaises(ValueError):
            constraints(data, [0])

    def run_toy(self, anchor=False, start_timestep=None):
        class Transition:
            def q_v_posterior(self, prediction, *args):
                return prediction

        class Toy(torch.nn.Module):
            num_timesteps = 3
            num_classes = 2
            num_bond_classes = 5
            model_mean_type = 'C0'
            bond_diffusion = True
            atom_type_trans = Transition()
            bond_type_trans = Transition()
            posterior_logvar = torch.full((3,), -20.)

            def forward(self, **inputs):
                return {'pred_ligand_pos': inputs['init_ligand_pos']+1,
                        'pred_ligand_v': torch.zeros(2, 2),
                        'pred_bond': torch.zeros(2, 5)}

            def q_pos_posterior(self, x0, **kwargs):
                return x0

        zeros = torch.zeros(2, dtype=torch.long)
        result = DecompScorePosNet3D.sample_diffusion(Toy(),
            protein_pos=torch.zeros(2, 3), protein_v=torch.zeros(2, 1),
            batch_protein=zeros, protein_group_idx=zeros,
            init_ligand_pos=torch.zeros(2, 3), init_ligand_v=zeros,
            ligand_v_aux=torch.zeros(2, 1), batch_ligand=zeros, ligand_group_idx=zeros,
            prior_centers=torch.zeros(1, 3), prior_stds=torch.ones(1, 3),
            prior_num_atoms=torch.tensor([2]), batch_prior=torch.tensor([0]), prior_group_idx=torch.tensor([0]),
            ligand_decomp_batch=zeros, ligand_decomp_index=zeros,
            ligand_fc_bond_index=torch.tensor([[0, 1], [1, 0]]),
            init_ligand_fc_bond_type=zeros, batch_ligand_bond=zeros,
            center_pos_mode='protein', fixed_atom_types=torch.tensor([1, -1]),
            fixed_bond_types=torch.tensor([3, 3]),
            fixed_position_mask=torch.tensor([True, False]) if anchor else None,
            start_timestep=start_timestep)
        return result

    def test_every_step_clamps_chemistry_without_freezing_coordinates(self):
        result = self.run_toy()
        for types, bonds in zip(result['v_traj'], result['bond_traj']):
            self.assertEqual(types[0].item(), 1)
            self.assertEqual(bonds.tolist(), [3, 3])
        self.assertTrue((result['pos'][0] > 2.9).all())

    def test_optional_anchor_freezes_only_fixed_coordinates(self):
        result = self.run_toy(anchor=True)
        for positions in result['pos_traj']:
            self.assertTrue(torch.equal(positions[0], torch.zeros(3)))
        self.assertTrue((result['pos'][1] > 2.9).all())

    def test_reference_local_schedule_reaches_zero(self):
        result = self.run_toy(start_timestep=1)
        self.assertEqual(len(result['pos_traj']), 2)
        self.assertTrue((result['pos'][0] > 1.9).all())

    def test_out_of_range_timestep_rejected(self):
        with self.assertRaises(ValueError):
            self.run_toy(start_timestep=3)


if __name__ == '__main__':
    unittest.main()
