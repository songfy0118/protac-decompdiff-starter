"""Experimental chemistry-constrained sampling; not a PROTAC-trained model."""
import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import torch
from easydict import EasyDict
from rdkit import Chem

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'DecompDiff-main'))
from models.decompdiff import DecompScorePosNet3D, log_sample_categorical
from models.transitions import index_to_log_onehot
from datasets.pl_data import torchify_dict
from utils.data import PDBProtein, ProteinLigandData, parse_sdf_file
from utils.prior import compute_golden_prior_from_data
from utils import transforms as trans
from utils import reconstruct as recon
from run_checks import positive_int, run_record, sha256, export_ready_mol
from validate_fixed_candidate import validate


def constraints(data, fixed_indices):
    """Keep all bonds/nonbonds touching fixed atoms, including the exit bond."""
    n = data.ligand_pos.size(0)
    fixed = torch.zeros(n, dtype=torch.bool, device=data.ligand_pos.device)
    fixed[fixed_indices] = True
    types = torch.full_like(data.ligand_atom_feature_full, -1)
    types[fixed] = data.ligand_atom_feature_full[fixed]
    matrix = torch.zeros((n, n), dtype=torch.long, device=fixed.device)
    matrix[data.ligand_bond_index[0], data.ligand_bond_index[1]] = data.ligand_bond_type
    a, b = data.ligand_fc_bond_index
    bonds = torch.full((len(a),), -1, dtype=torch.long, device=fixed.device)
    protected = fixed[a] | fixed[b]
    bonds[protected] = matrix[a[protected], b[protected]]
    # Sparse bond modes must contain every protected reference bond.
    represented = set(zip(a.tolist(), b.tolist()))
    for x, y in zip(*data.ligand_bond_index.tolist()):
        if (fixed[x] or fixed[y]) and (x, y) not in represented:
            raise ValueError('Model bond graph omits a required fixed/exit bond')
    return types, bonds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    parser.add_argument('--role', choices=['warhead', 'e3_ligand'], required=True)
    parser.add_argument('--protein', type=Path, help='Same-frame binding protein; defaults to STAT6')
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--samples', type=positive_int, default=1)
    parser.add_argument('--steps', type=positive_int, help='Short runs test plumbing only')
    parser.add_argument('--refine-timestep', type=int,
                        help='Reference-local editing: forward-noise to this timestep, then denoise to zero')
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--anchor-fixed-coordinates', action='store_true',
                        help='Stricter experimental mode: freeze the fixed fragment coordinates too')
    args = parser.parse_args()
    if args.role == 'e3_ligand' and args.protein is None:
        raise ValueError('E3 generation needs an aligned E3 binding protein; STAT6 alone is not an E3 pocket')
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    torch.set_num_threads(4)
    reference_path = args.case/'ak1690_observed.sdf'
    reference = next(iter(Chem.SDMolSupplier(str(reference_path))))
    tasks = json.loads((args.case/'generation_tasks.json').read_text())
    task = next(t for t in tasks if t['role'] == args.role)
    validate(reference, reference, task)  # Validate the task partition before loading a model.
    protein_path = args.protein or args.case/'stat6_observed.pdb'
    metadata = {'role': args.role, 'requested': args.samples, 'seed': args.seed,
                'checkpoint_sha256': sha256(args.checkpoint),
                'reference_sha256': sha256(reference_path), 'protein_sha256': sha256(protein_path),
                'anchor_fixed_coordinates': args.anchor_fixed_coordinates,
                'method': 'Experimental per-step atom/bond clamping; final chemistry must pass acceptance gate'}
    if (args.case/'provenance.json').exists():
        metadata['reference_provenance'] = json.loads((args.case/'provenance.json').read_text())
    if args.role == 'e3_ligand' and (args.case/'crbn_pose_audit.json').exists():
        metadata['e3_pose_audit'] = json.loads((args.case/'crbn_pose_audit.json').read_text())
    with run_record(args.output, metadata) as output:
        with torch.serialization.safe_globals([EasyDict]):
            checkpoint = torch.load(args.checkpoint, map_location='cpu', weights_only=True)
        cfg = checkpoint['config']; transform = cfg.data.transform
        if args.refine_timestep is not None and args.steps is not None:
            raise ValueError('Use either reference-local refinement or a step-limited smoke test')
        if args.refine_timestep is not None and not 0 <= args.refine_timestep < cfg.model.num_diffusion_timesteps:
            raise ValueError('Refinement timestep outside the trained schedule')
        steps = args.refine_timestep+1 if args.refine_timestep is not None else (args.steps or cfg.model.num_diffusion_timesteps)
        if steps > cfg.model.num_diffusion_timesteps:
            raise ValueError('Steps exceed the trained diffusion schedule')
        protein_parser = PDBProtein(str(protein_path))
        ligand = parse_sdf_file(str(reference_path))
        # Crop by whole residues, retaining the supplied coordinate frame.
        residues = protein_parser.query_residues_ligand(ligand, radius=10)
        pocket = protein_parser.residues_to_pdb_block(residues)
        (output/'pocket.pdb').write_text(pocket)
        protein = PDBProtein(pocket, mode='block').to_dict_atom()
        data = ProteinLigandData.from_protein_ligand_dicts(torchify_dict(protein), torchify_dict(ligand))
        groups = torch.full((reference.GetNumAtoms(),), -1, dtype=torch.long)
        groups[task['generate_indices']] = 0
        data.ligand_atom_mask = groups
        data.num_arms = 1; data.num_scaffold = 1
        compute_golden_prior_from_data(data)
        data.pocket_atom_masks = data.pocket_prior_masks[:1]
        if not data.pocket_atom_masks.any():
            raise ValueError('Generated endpoint has no nearby protein atoms; check binding pose')
        protein_features = trans.FeaturizeProteinAtom()
        ligand_features = trans.FeaturizeLigandAtom(transform.ligand_atom_mode,
                                                   prior_types=cfg.model.get('prior_types', False))
        indicator = trans.AddDecompIndicator(max_num_arms=transform.max_num_arms,
            global_prior_index=ligand_features.ligand_feature_dim, add_ord_feat=transform.get('add_ord_feat', True))
        data = protein_features(data); data = ligand_features(data)
        data = trans.ComputeLigandAtomNoiseDist(version='ref_prior')(data)
        data = indicator(data)
        data = trans.FeaturizeLigandBond(mode=transform.ligand_bond_mode, set_bond_type=True)(data)
        data = data.to(args.device)
        model = DecompScorePosNet3D(cfg.model,
            protein_atom_feature_dim=protein_features.protein_feature_dim+indicator.protein_feature_dim,
            ligand_atom_feature_dim=ligand_features.ligand_feature_dim+indicator.ligand_feature_dim,
            num_classes=ligand_features.ligand_feature_dim,
            prior_atom_types=ligand_features.atom_types_prob,
            prior_bond_types=ligand_features.bond_types_prob).to(args.device)
        model.load_state_dict(checkpoint['model'], strict=True); model.eval()
        fixed_types, fixed_bonds = constraints(data, task['fixed_indices'])
        results = []
        with Chem.SDWriter(str(output/'candidates.sdf')) as writer, \
                Chem.SDWriter(str(output/'accepted.sdf')) as accepted_writer:
            for sample in range(args.samples):
                n = reference.GetNumAtoms()
                centers = data.ligand_decomp_centers[data.ligand_decomp_mask]
                stds = data.ligand_decomp_stds[data.ligand_decomp_mask]
                pos = centers + torch.randn_like(data.ligand_pos)*stds
                # Fixed chemistry starts in its experimentally observed conformation.
                pos[task['fixed_indices']] = data.ligand_pos[task['fixed_indices']]
                types = log_sample_categorical(torch.zeros(n, model.num_classes, device=args.device))
                bonds = log_sample_categorical(torch.zeros(len(fixed_bonds), model.num_bond_classes, device=args.device))
                if args.refine_timestep is not None:
                    timestep = torch.tensor([args.refine_timestep], device=args.device)
                    alpha = model.alphas_cumprod[args.refine_timestep]
                    pos = alpha.sqrt()*(data.ligand_pos-centers)+(1-alpha).sqrt()*torch.randn_like(pos)*stds+centers
                    pos[task['fixed_indices']] = data.ligand_pos[task['fixed_indices']]
                    types, _ = model.atom_type_trans.q_v_sample(
                        index_to_log_onehot(data.ligand_atom_feature_full, model.num_classes), timestep,
                        torch.zeros(n, dtype=torch.long, device=args.device))
                    bonds, _ = model.bond_type_trans.q_v_sample(
                        index_to_log_onehot(data.ligand_fc_bond_type, model.num_bond_classes), timestep,
                        torch.zeros(len(bonds), dtype=torch.long, device=args.device))
                r = model.sample_diffusion(
                    protein_pos=data.protein_pos, protein_v=data.protein_atom_feature.float(),
                    batch_protein=torch.zeros(len(data.protein_pos), dtype=torch.long, device=args.device),
                    protein_group_idx=data.protein_decomp_group_idx,
                    init_ligand_pos=pos, init_ligand_v=types,
                    ligand_v_aux=data.ligand_atom_aux_feature.float(),
                    batch_ligand=torch.zeros(n, dtype=torch.long, device=args.device),
                    ligand_group_idx=data.ligand_decomp_group_idx,
                    fixed_position_mask=(fixed_types >= 0) if args.anchor_fixed_coordinates else None,
                    prior_centers=data.ligand_decomp_centers, prior_stds=data.ligand_decomp_stds,
                    prior_num_atoms=data.ligand_decomp_num_atoms,
                    batch_prior=torch.zeros(len(data.ligand_decomp_centers), dtype=torch.long, device=args.device),
                    prior_group_idx=data.prior_group_idx,
                    ligand_decomp_batch=data.ligand_decomp_mask, ligand_decomp_index=data.ligand_atom_mask,
                    ligand_fc_bond_index=data.ligand_fc_bond_index, init_ligand_fc_bond_type=bonds,
                    batch_ligand_bond=torch.zeros(len(bonds), dtype=torch.long, device=args.device),
                    num_steps=steps, center_pos_mode='protein', fixed_atom_types=fixed_types,
                    fixed_bond_types=fixed_bonds, start_timestep=args.refine_timestep)
                protected_atoms = fixed_types >= 0; protected_bonds = fixed_bonds >= 0
                invariant = all(torch.equal(v.to(args.device)[protected_atoms], fixed_types[protected_atoms])
                                for v in r['v_traj']) and all(
                    torch.equal(b.to(args.device)[protected_bonds], fixed_bonds[protected_bonds])
                    for b in r['bond_traj'])
                if not invariant:
                    raise RuntimeError('Chemical constraint violated during diffusion')
                np.savez_compressed(output/f'sample_{sample}_raw.npz',
                    positions=r['pos'].cpu().numpy(), atom_types=r['v'].cpu().numpy(),
                    bond_types=r['bond'].cpu().numpy(), bond_index=data.ligand_fc_bond_index.cpu().numpy())
                row = {'sample': sample, 'trajectory_constraints_preserved': invariant}
                try:
                    mol = recon.reconstruct_from_generated_with_bond(r['pos'].cpu().numpy().astype(np.float64),
                        trans.get_atomic_number_from_index(r['v'].cpu(), mode=transform.ligand_atom_mode),
                        data.ligand_fc_bond_index.cpu().tolist(), r['bond'].cpu().numpy())
                    mol = export_ready_mol(mol)
                    row.update(validate(reference, mol, task)); row['smiles'] = Chem.MolToSmiles(mol)
                    mol.SetProp('sample_id', str(sample)); mol.SetProp('fixed_check_passed', str(row['accepted']))
                    mol.SetProp('role', args.role)
                    mol.SetProp('status', 'chemistry_accepted_activity_unvalidated' if row['accepted'] else 'constraint_rejected')
                    writer.write(mol)
                    if row['accepted']:
                        accepted_writer.write(mol)
                except (recon.MolReconsError, ValueError, RuntimeError) as exc:
                    row.update(accepted=False, reason='reconstruction_failed', error=str(exc))
                results.append(row)
                print(json.dumps(row), flush=True)
        report = dict(metadata, steps=steps, full_schedule=steps==cfg.model.num_diffusion_timesteps,
            refine_timestep=args.refine_timestep,
            results=results, accepted=sum(r['accepted'] for r in results),
            accepted_changed=sum(r['accepted'] and r.get('changed', False) for r in results),
            limits='Untrained constraint adaptation; fixed reference node count; H predictions may change exported heavy-atom count; no affinity, degradation or synthesis validation.')
        (output/'report.json').write_text(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
