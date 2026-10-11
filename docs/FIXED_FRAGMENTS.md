# 固定端点：当前完成程度

完成 9BIG 结构准备与候选化学约束检查，并新增实验性的固定端点采样入口。它是对原预训练模型的约束适配，不是已训练或验证过的 PROTAC 专用模型。

## 本地结构准备

```bash
python -m pip install -r requirements-structure.txt
python inspect_9big.py --cif /local/9BIG.cif --archive /local/data.zip --output runs/9big-case
```

原始压缩包和输出不提交。程序导出同坐标系的 STAT6 PDB、AK-1690 SDF、原子映射、两个方向的固定/生成任务和数据审计。
仅适用于本例 A1AQQ 化学组分，七碳 linker 的边界是明确假设：氧端留在 warhead、芳香端留在 E3。实际生成 linker 尚需另行替换。
9BIG 是二元结构；压缩包里的 CRBN 是另外体系的参考，尚未对齐，不能当作实验三元结构。

## 候选检查

```bash
python validate_fixed_candidate.py --case runs/9big-case --candidates /local/candidates.sdf --role warhead --output runs/check-warhead.json
```

另一方向使用 `--role e3_ligand`。检查固定原子的元素、电荷、同位素、氢数、芳香性、立体化学、固定内部键以及跨边界出口键。允许构象和原子顺序变化，拒绝额外不连通片段。任务分区必须完整、无重复且两侧非空。
`accepted` 表示满足固定条件；`accepted_changed` 额外排除与参考相同的分子。通过检查不能证明结合、活性、合成可行性。子结构匹配最多枚举 10,000 个匹配，高对称性候选可能被保守拒绝。

## 验证

```bash
python -m unittest -v test_fixed_candidate
```

7 项合成分子测试无需私有资料。设置 `PROTAC_CASE` 为本地结构准备目录后，额外运行 7 项 AK-1690 回归测试。模型推理环境与结构检查环境分开。

## 实验性采样

```bash
python run_fixed_endpoint.py --case runs/9big-case --role warhead --checkpoint assets/decompdiff.pt --output runs/fixed-warhead --samples 1
```

默认运行完整扩散日程。在每一步固定原子类型、所有涉及固定原子的键及非键，防止出现新的连接位点。坐标默认可以变化；`--anchor-fixed-coordinates` 可额外冻结固定片段的参考坐标。采样节点数暂与参考相同，预测为氢的节点可能影响导出的重原子数，不能据此主张支持任意端点大小或一般化生成。

`--refine-timestep 99` 使用参考结构加噪后，从第 99 步完整去噪到 0，适合参考附近的局部编辑。它与从纯噪声运行完整日程不同，必须分别报告。`--steps 5` 则仅运行日程开头的少数步骤，不能用于判断生成质量。

`report.json` 记录重建和最终化学检查；`sample_*_raw.npz` 保存原始采样结果。`candidates.sdf` 包含所有可重建分子，`accepted.sdf` 仅保存通过固定条件的分子。没有通过者时该文件为空。原模型可能重建失败或被筛除，采样完成不能代替候选通过。

顺序任务先生成 warhead，再把通过检查的新候选作为下一阶段的固定参考：

```bash
python make_followup_case.py --parent-case runs/9big-case --run runs/fixed-warhead --sample 0 --output runs/followup-case
python run_fixed_endpoint.py --case runs/followup-case --role e3_ligand --protein /local/aligned-crbn.pdb --checkpoint assets/decompdiff.pt --output runs/fixed-e3 --samples 1
```

跟进案例会核对参考哈希、样本通过情况及原始节点与 SDF 坐标顺序。节点数或顺序改变时拒绝复用原分区。`provenance.json` 明确标记这是生成的参考，不是实验 AK-1690。

E3 方向必须额外提供 `--protein`，其坐标需要与完整参考配体一致；程序不会把 STAT6 当作 E3 结合口袋。`prepare_crbn_hypothesis.py` 可转移配体核心匹配的 CRBN 姿态并记录 RMSD/碰撞。默认归档 S 构型配体与本例 R 构型不匹配，会拒绝；可用 `--source-pdb` 指定公开 5YJ1 R 构型参考，但它是小鼠结合域，这种转移仅供技术试跑，不是人源体系或可靠三元结构。

模型环境验证路线：Python 3.12、PyTorch 2.8.0+cu128、仓库 requirements.txt、与 torch-2.8.0+cu128 对应的 PyG 扩展 wheel；与纯结构环境分开。公开作者权重 SHA-256 仍为 docs/INFERENCE.md 所列值。

```bash
python -m unittest -v test_pipeline test_fixed_candidate test_diffusion_constraints
```

采样约束测试无需 checkpoint，包括出口键保持、禁止额外跨边界键，以及每一步保持化学信息但允许坐标变化。DeepTernary 仍属于后续三元构象预测/评估路线。
