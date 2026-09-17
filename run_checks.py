"""Shared checks for reproducible local sampling runs."""
import hashlib
import json
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from rdkit import Chem

def export_ready_mol(mol):
    # Derive report SMILES from the same serialization a downstream SDF reader sees.
    exported = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), sanitize=True, removeHs=True)
    if exported is None: raise ValueError('Generated molecule failed SDF round trip')
    return exported

def positive_int(value):
    value = int(value)
    if value <= 0:
        raise ValueError('sample count and step count must be positive')
    return value

def fresh_output(path):
    path = Path(path).resolve()
    if path.exists() and not path.is_dir():
        raise ValueError(f'Output path must be a directory: {path}')
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f'Output directory is not empty; choose a new directory: {path}')
    path.mkdir(parents=True, exist_ok=True)
    return path

@contextmanager
def run_record(output, metadata):
    """Record failures and interruptions as well as successful completion."""
    output=fresh_output(output)
    record=dict(metadata, status='running', started_at=datetime.now(timezone.utc).isoformat())
    destination=output/'run_manifest.json'
    # Exclusive creation also prevents two processes claiming the same empty folder.
    with destination.open('x',encoding='utf8') as stream:
        json.dump(record,stream,indent=2)
    try:
        yield output
    except BaseException as exc:
        record.update(status='interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed',
                      error_type=type(exc).__name__,error=str(exc))
        raise
    else:
        record['status']='completed'
    finally:
        record['finished_at']=datetime.now(timezone.utc).isoformat()
        temporary=destination.with_suffix('.tmp')
        temporary.write_text(json.dumps(record,indent=2),encoding='utf8')
        temporary.replace(destination)

def load_endpoint_spec(path):
    path=Path(path).resolve()
    spec=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(spec,dict): raise ValueError('Input specification must be a JSON object')
    required={'role','protein_pdb','reference_sdf','atom_groups'}
    if required-set(spec): raise ValueError(f'Missing specification fields: {sorted(required-set(spec))}')
    if spec['role'] not in ('warhead','e3_ligand','test_demo'):
        raise ValueError('role must be warhead, e3_ligand, or test_demo')
    for key in ('protein_pdb','reference_sdf'):
        if not isinstance(spec[key],str) or not spec[key].strip():
            raise ValueError(f'{key} must be a nonempty file path')
        if not (path.parent/spec[key]).is_file(): raise FileNotFoundError(f'{key}: {spec[key]}')
    return spec

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()
