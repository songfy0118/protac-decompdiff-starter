"""Pocket-conditioned endpoint sampling with a reference decomposition.

The role is output metadata, not a newly trained conditioning token. Reference
coordinates must already be aligned to the supplied protein. This does not
predict PROTAC activity or select an attachment site.
"""
import argparse
import json
import logging
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
import torch
from easydict import EasyDict
from rdkit import Chem
from torch_geometric.transforms import Compose

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'DecompDiff-main'))
from utils.data import PDBProtein, ProteinLigandData, parse_sdf_file
from datasets.pl_data import torchify_dict
from utils.prior import compute_golden_prior_from_data
import utils.transforms as trans
from models.decompdiff import DecompScorePosNet3D
import scripts.sample_diffusion_decomp as sampling
from run_checks import positive_int, sha256, export_ready_mol, run_record, load_endpoint_spec

from input_validation import validate_groups, validate_coordinates

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--spec',type=Path,required=True)
    p.add_argument('--checkpoint',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--samples',type=positive_int,default=3)
    p.add_argument('--device',default='cuda:0')
    p.add_argument('--seed',type=int,default=2021)
    args=p.parse_args(); spec=load_endpoint_spec(args.spec)
    def path(value): return str((args.spec.resolve().parent/value).resolve())
    protein=PDBProtein(path(spec['protein_pdb'])).to_dict_atom()
    ligand=parse_sdf_file(path(spec['reference_sdf']))
    validate_coordinates(protein['pos'],ligand['pos'])
    groups=spec['atom_groups']; arms=validate_groups(groups,len(ligand['element']))
    data=ProteinLigandData.from_protein_ligand_dicts(torchify_dict(protein),torchify_dict(ligand))
    data.ligand_atom_mask=torch.tensor(groups,dtype=torch.long)
    data.num_arms=arms; data.num_scaffold=int(-1 in groups)
    compute_golden_prior_from_data(data)
    data.pocket_atom_masks=data.pocket_prior_masks[:arms]
    if not data.pocket_atom_masks.any(dim=1).all():
        raise ValueError('At least one arm has no nearby pocket atoms; check coordinate alignment')
    for prior in data.arms_prior+data.scaffold_prior:
        if not np.isfinite(prior[2]).all() or np.trace(prior[2])<=0:
            raise ValueError('Reference group has invalid or zero coordinate variance')
    with torch.serialization.safe_globals([EasyDict]):
        ckpt=torch.load(args.checkpoint,map_location='cpu',weights_only=True)
    cfg=ckpt['config']; transform=cfg.data.transform
    protein_features=trans.FeaturizeProteinAtom()
    ligand_features=trans.FeaturizeLigandAtom(transform.ligand_atom_mode,prior_types=cfg.model.get('prior_types',False))
    indicator=trans.AddDecompIndicator(max_num_arms=transform.max_num_arms,
        global_prior_index=ligand_features.ligand_feature_dim,add_ord_feat=transform.get('add_ord_feat',True))
    initial=Compose([trans.ComputeLigandAtomNoiseDist(version='ref_prior'),indicator,
                     trans.FeaturizeLigandBond(mode=transform.ligand_bond_mode,set_bond_type=False)])
    metadata={'role':spec['role'],'requested':args.samples,'seed':args.seed,
              'steps':cfg.model.num_diffusion_timesteps,
              'checkpoint_sha256':sha256(args.checkpoint),
              'protein_sha256':sha256(path(spec['protein_pdb'])),
              'reference_sha256':sha256(path(spec['reference_sdf']))}
    with run_record(args.output,metadata) as output:
        args.output=output
        model=DecompScorePosNet3D(cfg.model,
            protein_atom_feature_dim=protein_features.protein_feature_dim+indicator.protein_feature_dim,
            ligand_atom_feature_dim=ligand_features.ligand_feature_dim+indicator.ligand_feature_dim,
            num_classes=ligand_features.ligand_feature_dim,
            prior_atom_types=ligand_features.atom_types_prob,prior_bond_types=ligand_features.bond_types_prob).to(args.device)
        model.load_state_dict(ckpt['model'],strict=True); model.eval()
        logging.basicConfig(level=logging.INFO)
        sampling.logger=logging.getLogger('endpoint')
        sampling.full_protein_pos=data.protein_pos.clone()
        sampling.args=SimpleNamespace(recon_with_bond=True)
        sampling.misc.seed_all(args.seed)
        results=sampling.sample_diffusion_ligand_decomp(model,protein_features(data),initial,
            num_samples=args.samples,batch_size=1,device=args.device,prior_mode='ref_prior',
            num_steps=cfg.model.num_diffusion_timesteps,center_pos_mode='protein',num_atoms_mode='ref',
            atom_enc_mode=transform.ligand_atom_mode,bond_fc_mode=transform.ligand_bond_mode,
            atom_prior_probs=ligand_features.atom_types_prob,bond_prior_probs=ligand_features.bond_types_prob,
            energy_drift_opt=[EasyDict(type='armsca_prox',min_d=1.2,max_d=1.9),EasyDict(type='clash',sigma=2,gamma=4)])
        writer=Chem.SDWriter(str(args.output/'candidates.sdf'))
        rows=[]
        for i,result in enumerate(results):
            mol=result['mol']
            if mol is None:
                rows.append({'sample':i,'smiles':None,'connected':False,'status':'reconstruction_failed'})
                continue
            mol=export_ready_mol(mol)
            smiles=Chem.MolToSmiles(mol)
            mol.SetProp('sample_id',str(i))
            mol.SetProp('role',spec['role']); mol.SetProp('status','unvalidated_candidate')
            writer.write(mol)
            rows.append({'sample':i,'smiles':smiles,'connected':len(Chem.GetMolFrags(mol))==1})
        writer.close()
        report={'role':spec['role'],'requested':args.samples,
                'steps':cfg.model.num_diffusion_timesteps,'seed':args.seed,'checkpoint_sha256':sha256(args.checkpoint),
                'spec':spec,'protein_sha256':sha256(path(spec['protein_pdb'])),
                'reference_sha256':sha256(path(spec['reference_sdf'])),
                'reconstructed':sum(r['smiles'] is not None for r in rows),'results':rows,
                'method':'Unmodified pretrained DecompDiff with reference-derived priors',
                'stereochemistry':'Inferred by RDKit from generated 3D coordinates; unvalidated.',
                'limits':'No PROTAC fine-tuning, attachment-site guarantee, affinity or degradation validation.'}
        (args.output/'report.json').write_text(json.dumps(report,indent=2),encoding='utf8')
        print(json.dumps(report,indent=2))

if __name__=='__main__': main()
