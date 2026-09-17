"""Input checks that do not require the GPU/model environment."""
import numpy as np

def validate_groups(groups, count):
    if not isinstance(groups,list) or len(groups)!=count or any(type(x) is not int for x in groups):
        raise ValueError('groups must contain one integer per hydrogen-removed SDF atom')
    arms=sorted(set(groups)-{-1})
    if not arms or arms!=list(range(len(arms))) or len(arms)>10:
        raise ValueError('arm groups must be contiguous 0..N-1, with 1 <= N <= 10')
    if any(groups.count(g)<2 for g in set(groups)):
        raise ValueError('Every group must have at least 2 atoms for a nondegenerate prior')
    return len(arms)

def validate_coordinates(protein, ligand):
    for name,coords in [('protein',protein),('ligand',ligand)]:
        coords=np.asarray(coords)
        if coords.ndim!=2 or coords.shape[1]!=3 or len(coords)==0 or not np.isfinite(coords).all():
            raise ValueError(f'{name} must contain finite nonempty N x 3 coordinates')
