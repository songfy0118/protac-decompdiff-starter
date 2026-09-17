"""Report installed capabilities without exposing local paths or credentials."""
import argparse
import importlib
import json
import platform

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inference',action='store_true')
    parser.add_argument('--device',choices=['cpu','cuda'],default='cuda')
    args=parser.parse_args()
    packages=['numpy','rdkit','openpyxl']
    if args.inference:
        packages+=['torch','torch_geometric','torch_scatter','torch_sparse','torch_cluster','openbabel','lmdb','easydict','scipy','sklearn','yaml','tqdm']
    report={'python':platform.python_version(),'mode':'inference' if args.inference else 'data','checks':{}}
    for package in packages:
        try:
            module=importlib.import_module(package)
            report['checks'][package]={'ok':True,'version':getattr(module,'__version__',None)}
        except Exception as exc:
            report['checks'][package]={'ok':False,'error':f'{type(exc).__name__}: {exc}'}
    if args.inference and report['checks']['torch']['ok']:
        import torch
        try:
            if args.device=='cuda' and not torch.cuda.is_available(): raise RuntimeError('CUDA is unavailable')
            value=torch.ones(2,device=args.device).sum().item()
            report['checks']['device']={'ok':value==2,'device':args.device,'cuda_build':torch.version.cuda}
        except Exception as exc:
            report['checks']['device']={'ok':False,'error':str(exc)}
    report['ok']=all(check['ok'] for check in report['checks'].values())
    print(json.dumps(report,indent=2))
    return 0 if report['ok'] else 1

if __name__=='__main__': raise SystemExit(main())
