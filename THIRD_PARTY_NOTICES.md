# Third-party notices

`DecompDiff-main/` is adapted from ByteDance's official DecompDiff implementation:
https://github.com/bytedance/DecompDiff

Paper: Jiaqi Guan, Xiangxin Zhou, Yuwei Yang, Yu Bao, Jian Peng, Jianzhu Ma, Qiang Liu, Liang Wang, Quanquan Gu. *DecompDiff: Diffusion Models with Decomposed Priors for Structure-Based Drug Design*. ICML 2023, PMLR 202.

The upstream LICENSE is retained (CC BY-NC 4.0). Compatibility edits include explicit Python package initialization, safe PyTorch checkpoint loading, and NumPy 2 integer types/input validation. The model architecture and training objective are unchanged.

Preserve file-specific notices as well: upstream `utils/reconstruct.py` identifies code originating from liGAN fitting code and carries a GPL v2 notice. The root license does not erase that notice. Dependencies retain their respective licenses. Model weights and full datasets are not redistributed here.

The new project wrappers and documentation are supplied under CC BY-NC 4.0; third-party files retain their original notices and applicable terms. No endorsement by the upstream authors is implied.

`examples/test_pocket.pdb` and `examples/test_reference.sdf` derive from the supplied official CrossDocked test case 0 (BSD_ASPTE; structure identifier 2z3h); atom groups come from official processed record 99165. They demonstrate the pipeline, not PROTAC biological activity.

`examples/table_inputs/` contains synthetic table fixtures only. Private source tables and derived full fragment libraries are not included.
