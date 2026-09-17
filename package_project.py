"""Create and verify a portable source/data bundle using only the standard library."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile

def verify_archive(path):
    with zipfile.ZipFile(path) as archive:
        names=archive.namelist()
        if len(names)!=len(set(names)): raise ValueError('Duplicate archive entry')
        for name in names:
            if PurePosixPath(name).is_absolute() or '..' in PurePosixPath(name).parts or '\\' in name or ':' in name:
                raise ValueError('Unsafe archive path')
        manifest=json.loads(archive.read('PACKAGE_MANIFEST.json'))
        if set(names)-{'PACKAGE_MANIFEST.json'}!=set(manifest['files']):
            raise ValueError('Archive contents differ from manifest')
        for name,digest in manifest['files'].items():
            if hashlib.sha256(archive.read(name)).hexdigest()!=digest:
                raise ValueError(f'Checksum mismatch: {name}')
        return len(manifest['files'])

def build_archive(root, output):
    root=Path(root).resolve(); output=Path(output).resolve()
    if output.is_relative_to(root): raise ValueError('Put the archive outside the project folder')
    if output.exists(): raise FileExistsError(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    manifest={'format_version':1,'files':{},'excluded':'Local data, environments, weights, runs, caches, secrets and optional SA score asset'}
    temporary=output.with_suffix('.partial')
    with zipfile.ZipFile(temporary,'x',zipfile.ZIP_DEFLATED) as archive:
        for file in sorted(root.rglob('*')):
            if file.is_symlink(): raise ValueError(f'Symlinks are not packaged: {file}')
            parts=file.relative_to(root).parts
            if not file.is_file() or any(part in ('__pycache__','.git','assets','data','runs','results','work') or part.startswith('.venv') for part in parts): continue
            if file.suffix in ('.pt','.ckpt','.pyc','.zip','.partial') or '.lmdb' in file.name: continue
            if file.name in ('assets_local.json','fpscores.pkl.gz') or file.name=='.env' or file.name.startswith('.env.'): continue
            name=(PurePosixPath(root.name)/file.relative_to(root).as_posix()).as_posix()
            data=file.read_bytes(); manifest['files'][name]=hashlib.sha256(data).hexdigest()
            archive.writestr(name,data)
        archive.writestr('PACKAGE_MANIFEST.json',json.dumps(manifest,indent=2))
    count=verify_archive(temporary)
    temporary.replace(output)
    return count

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--output',type=Path)
    action.add_argument('--verify',type=Path)
    args=parser.parse_args()
    count=verify_archive(args.verify) if args.verify else build_archive(Path(__file__).parent,args.output)
    print(f'Verified {count} files')
