"""Run the supplied checkpoint against the official reference-prior test data."""
import argparse
import os
from pathlib import Path
import runpy
import sys
import json
import torch
import yaml
from easydict import EasyDict
from run_checks import positive_int, sha256, run_record

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--assets', type=Path, required=True)
    parser.add_argument('--samples', type=positive_int, default=1)
    parser.add_argument('--data-id', type=int, default=0)
    parser.add_argument('--steps', type=positive_int, default=1000)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--device', default='cuda:0')
    args = parser.parse_args()
    assets = args.assets.resolve()
    repo = ROOT / 'DecompDiff-main'
    if args.data_id < 0: parser.error('--data-id must be nonnegative')
    with torch.serialization.safe_globals([EasyDict]):
        checkpoint = torch.load(assets/'decompdiff.pt', map_location='cpu', weights_only=True)
    if args.steps != checkpoint['config'].model.num_diffusion_timesteps:
        parser.error('--steps must match the trained diffusion horizon; truncated reverse steps are not accelerated full sampling')
    manifest = {'data_id':args.data_id,'steps':args.steps,'requested':args.samples,
                'checkpoint_sha256':sha256(assets/'decompdiff.pt'),'status':'running'}
    with run_record(args.output,manifest) as out:
        def plain(x):
            if isinstance(x, dict): return {k:plain(v) for k,v in x.items()}
            if isinstance(x, (tuple,list)): return [plain(v) for v in x]
            return x
        train = plain(checkpoint['config'])
        train['data']['path'] = str(assets/'data/crossdocked_v1.1_rmsd1.0_processed')
        train['data']['split'] = str(assets/'data/split_by_name.pt')
        train_path = out/'training_local.yml'
        train_path.write_text(yaml.safe_dump(train), encoding='utf8')
        config = yaml.safe_load((repo/'configs/sampling_drift.yml').read_text())
        config['model']['checkpoint'] = str(assets/'decompdiff.pt')
        config['model']['train_config'] = str(train_path)
        config['sample']['num_samples'] = args.samples
        config['sample']['num_steps'] = args.steps
        config_path = out/'sampling_local.yml'
        config_path.write_text(yaml.safe_dump(config), encoding='utf8')
        del checkpoint
        old_cwd,old_argv,old_path=os.getcwd(),sys.argv[:],sys.path[:]
        sys.path.insert(0,str(repo))
        os.chdir(repo)
        sys.argv = ['sample_diffusion_decomp.py',str(config_path),'--outdir',str(out),
                    '--ori_data_path',str(assets/'data/test_set'),
                    '--index_path',str(assets/'data/test_index.pkl'),
                    '--batch_size','1','--data_id',str(args.data_id),'--device',args.device]
        try:
            runpy.run_path(str(repo/'scripts/sample_diffusion_decomp.py'),run_name='__main__')
        finally:
            os.chdir(old_cwd)
            sys.argv=old_argv
            sys.path[:]=old_path

if __name__ == '__main__': main()
