# DeepTernary Windows CPU compatibility

Validated against upstream commit `827821dccca31a5918bd0355e2d6bf70c072b6dd` on Python 3.12 with a separate environment. Do not replace the existing DecompDiff environment.

Install CPU torch 2.2.1 and torchvision 0.17.1 from the official CPU wheel index, then install `requirements-deepternary-windows.txt`. Explicitly pin torch during dependency resolution to prevent unintended upgrades. These versions differ from the upstream Linux CUDA environment; a successful inference verifies execution, not equivalence of every numerical result.

Apply `patches/deepternary-sdf-export.patch` to that upstream checkout. The documented SDF input is otherwise passed to a PDB-only reader at export. The patch exports the already-loaded RDKit molecule, retains atom order, and additionally saves predicted ligand SDF files. It does not change model inference.

Set process environment `DGLBACKEND=pytorch`, and set `TEMP` and `TMP` to a writable local project directory in sandboxed environments. Use explicit `--config` and `--checkpoint` paths rather than relying on the checkpoint archive's original machine paths.

Official resources:
- https://github.com/youqingxiaozhua/DeepTernary
- https://github.com/youqingxiaozhua/DeepTernary/releases/tag/v1.0.0
- PROTAC checkpoint SHA256: `4ac2b4c67d05705e5c1907f878b6865b8b788cbd9111578965e01fd6c1b4e683`

Validation: checkpoint load and a 40-seed CPU prediction completed, exporting 40 complex PDBs and 40 ligand SDFs. Private input structures, predictions, and case-level metrics are not included in this repository.

Interpretation: the predicted P2 RMSD is a surrogate ranking score, not RMSD measured against an experimental ternary structure. The upstream C-alpha clash ratio does not replace all-heavy-atom ligand/protein clash checks. Independently audit bond lengths, ligand self-overlaps, stereochemistry derived from coordinates, and ligand/protein contacts. Positional-restraint MMFF relaxation can repair local geometry but does not establish binding affinity, degradation activity, or biological relevance of a transferred protein template.
