"""Transfer human CRBN using homologous backbone atoms, never mirror ligand stereo."""
import json
import argparse
from pathlib import Path
import numpy as np
from Bio import Align
from Bio.SeqUtils import seq1

def residues(path,chain):
    data=[]
    for line in path.read_text().splitlines():
        if line.startswith('ATOM  ') and line[21]==chain and line[12:16].strip()=='CA' and line[16] in (' ','A'):
            data.append((seq1(line[17:20]),np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]),line[22:27].strip()))
    return data

parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,default=Path('work/4TZ4.pdb'))
parser.add_argument('--target',type=Path,default=Path('outputs/friend-linker-cases-2026-10-12/5_/crbn_transferred_hypothesis.pdb'))
parser.add_argument('--out',type=Path,default=Path('outputs/human-crbn-transfer-2026-10-12'))
args=parser.parse_args()
human,mouse,out=args.source,args.target,args.out
header=human.read_text()[:12000]
assert '4TZ4' in header and 'HOMO SAPIENS' in header and 'PROTEIN CEREBLON' in header
out.mkdir(exist_ok=True)
source=residues(human,'C'); target=residues(mouse,'A')
aligner=Align.PairwiseAligner(mode='local',match_score=2,mismatch_score=-1,open_gap_score=-5,extend_gap_score=-.5)
alignment=aligner.align(''.join(x[0] for x in source),''.join(x[0] for x in target))[0]
pairs=[]
for (sa,sb),(ta,tb) in zip(*alignment.aligned):
    pairs.extend((i,j) for i,j in zip(range(sa,sb),range(ta,tb)))
x=np.array([source[i][1] for i,j in pairs]); y=np.array([target[j][1] for i,j in pairs])
u,_,vt=np.linalg.svd((x-x.mean(0)).T@(y-y.mean(0)))
correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
rotation=u@correction@vt; translation=y.mean(0)-x.mean(0)@rotation
rmsd=float(np.sqrt(((x@rotation+translation-y)**2).sum(1).mean()))
print('Matched CA / initial RMSD:',len(pairs),rmsd,flush=True)
initial_rmsd=rmsd
keep=np.ones(len(pairs),dtype=bool)
for iteration in range(5):
    residual=np.linalg.norm(x@rotation+translation-y,axis=1)
    updated=residual<3
    if np.array_equal(updated,keep):
        break
    keep=updated
    assert keep.sum()>=60
    xx,yy=x[keep],y[keep]
    u,_,vt=np.linalg.svd((xx-xx.mean(0)).T@(yy-yy.mean(0)))
    correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
    rotation=u@correction@vt;translation=yy.mean(0)-xx.mean(0)@rotation
rmsd=float(np.sqrt(((x[keep]@rotation+translation-y[keep])**2).sum(1).mean()))
assert len(pairs)>=80 and keep.sum()>=60 and rmsd<2 and np.linalg.det(rotation)>.999
lines=[]
for line in human.read_text().splitlines():
    if line.startswith('ATOM  ') and line[21]=='C' and line[16] in (' ','A'):
        point=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])@rotation+translation
        lines.append(line[:30]+''.join(f'{v:8.3f}' for v in point)+line[54:])
(out/'human_crbn_transferred.pdb').write_text('\n'.join(lines)+'\nEND\n')
report=dict(source='4TZ4 human CRBN chain C; experimentally lenalidomide-bound',
    target='Prior 5YJ1 mouse-domain transfer into linker5 input coordinate frame',
    matched_CA=len(pairs),sequence_identity=sum(source[i][0]==target[j][0] for i,j in pairs)/len(pairs),
    initial_all_matched_rmsd_A=initial_rmsd,fit_CA=int(keep.sum()),excluded_CA=int((~keep).sum()),
    outlier_threshold_A=3,backbone_rmsd_A=rmsd,rotation_determinant=float(np.linalg.det(rotation)),
    status='human_protein_homology_transferred_pose',
    limits='Human CRBN protein coordinates are experimental; candidate binding pose is a homology-transfer hypothesis. No experimental R-candidate binding pose is claimed. DDB1 excluded.')
(out/'audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
