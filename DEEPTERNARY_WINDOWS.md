# DeepTernary Windows CPU compatibility

Validated against upstream commit `827821dccca31a5918bd0355e2d6bf70c072b6dd` on Python 3.12 with a separate environment. Do not replace the existing DecompDiff environment.

Install CPU torch 2.2.1 and torchvision 0.17.1 from the official CPU wheel index, then install `requirements-deepternary-windows.txt`. Explicitly pin torch during dependency resolution to prevent unintended upgrades. These versions differ from the upstream Linux CUDA environment; a successful inference verifies execution, not equivalence of every numerical result.

Apply `patches/deepternary-inference.patch` to a clean checkout of that upstream commit. This comprehensive patch supersedes `deepternary-sdf-export.patch`; do not apply both. It includes SDF export and stereo-preserving conformer preparation. Model weights are unchanged, but starting conformers and corrected outputs change.

The documented SDF input is otherwise passed to a PDB-only reader at export. The patch exports the already-loaded RDKit molecule, retains atom order, and additionally saves predicted ligand SDF files.

Conformers are embedded with explicit hydrogens and enforced chirality, then converted back to the original heavy-atom graph. Coordinate-derived stereochemistry must match the reference before inference. A failed conformer raises an error instead of automatically downloading a different molecule. Correction coordinates are cloned before the model mutates graph coordinates in place. Local geometry repair uses a 1200-iteration budget and still records convergence separately from geometry checks.

Set process environment `DGLBACKEND=pytorch`, and set `TEMP` and `TMP` to a writable local project directory in sandboxed environments. Use explicit `--config` and `--checkpoint` paths rather than relying on the checkpoint archive's original machine paths.

Official resources:
- https://github.com/youqingxiaozhua/DeepTernary
- https://github.com/youqingxiaozhua/DeepTernary/releases/tag/v1.0.0
- PROTAC checkpoint SHA256: `4ac2b4c67d05705e5c1907f878b6865b8b788cbd9111578965e01fd6c1b4e683`

Validation: checkpoint load and a 40-seed CPU prediction completed, exporting 40 complex PDBs and 40 ligand SDFs. Private input structures, predictions, and case-level metrics are not included in this repository.

Interpretation: the predicted P2 RMSD is a surrogate ranking score, not RMSD measured against an experimental ternary structure. The upstream C-alpha clash ratio does not replace all-heavy-atom ligand/protein clash checks. Independently audit bond lengths, ligand self-overlaps, stereochemistry derived from coordinates, and ligand/protein contacts. Positional-restraint MMFF relaxation can repair local geometry but does not establish binding affinity, degradation activity, or biological relevance of a transferred protein template.

## Local evaluation helpers

`deepternary_tools/` provides scripts for private local cases. Run from a directory containing your case data; pass explicit paths rather than relying on the development-case defaults. None of the research inputs are bundled.

1. `transfer_human_crbn.py --source 4TZ4.pdb --target transferred_mouse_crbn.pdb --out human_transfer` transfers the human chain-C backbone using homologous CA atoms. It records sequence identity, excluded outliers, RMSD, and proper-rotation determinant. This uses an experimental human protein reference, but the candidate pose remains a homology-transfer hypothesis. The original reference is lenalidomide-bound, not an experimental structure of the candidate. DDB1 is excluded.
2. `prepare_deepternary_inputs.py --case CASE --ligand GENERATED.sdf --protein1 human_transfer/human_crbn_transferred.pdb --out INPUTS` uses `generation_tasks.json` endpoint indices and `stat6_observed.pdb` from CASE. Endpoints are extracted from the selected molecule in its atom order, not silently substituted with the original molecule. Atom counts and index partitions must be compatible with that case.
3. Using the DeepTernary environment, run `run_deepternary_linker5.py --upstream UPSTREAM --inputs INPUTS --out PREDICTIONS --name SAMPLE --seeds 40`.
4. `audit_deepternary.py --out PREDICTIONS --reference INPUTS/ligand.sdf` checks each exported ligand against its reference, including stereochemistry reassigned from coordinates. Unchanged stored stereo labels alone are insufficient.
5. `refine_deepternary.py --out PREDICTIONS --reference INPUTS/ligand.sdf` locally relaxes at most five highest-ranked stereochemically compatible predictions. It records convergence, displacement, clashes, and stereo. Nonconverged outputs remain exploratory. Recheck full-protein contacts and bond lengths before selecting any result.

Evaluation validation: original and newly generated ligands each completed 40-seed CPU runs against the transferred human reference. Geometry audits and local refinement completed; the existing four geometry-quality regression tests passed. Successful execution does not imply every output is valid or that the generated ligand outperforms its parent.

Stereo regression: `test_deepternary_stereo.py` checks multiple seeds of a molecule with multiple stereocenters, unchanged heavy-atom order, and rejection of forced inverted coordinates. It reads only the two pure conformer functions from the optional checkout, avoiding model imports. Set `DEEPTERNARY_SOURCE` to the patched `predict.py` if running outside the development workspace. These two tests and the four geometry-quality tests passed. Both private 40-seed prediction sets were rerun after the correction; structures and case metrics remain local.
