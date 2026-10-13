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

### 仅有 DiffPROTACs XYZ 时

官方 `Fenglei104/DiffPROTACs` 的采样路径预测元素与坐标，使用 `visualizer.save_xyz_file` 导出 XYZ。`dataset.parse_xyz` 在输入固定片段后追加重叠的碳占位点；`trainer.pred` 将这一输入也保存为 `true_.xyz`。自定义采样中的该文件不是实验参考分子，因此不能用于真实结构恢复率计算。

本项目可对本例的63个固定原子，恢复9BIG已知片段的内部键，并用OpenBabel推断linker键。新增脚本不会把推断SDF标记为验证过的化学输入：

```bash
python audit_linker_xyz.py --archive /local/linker-xyz.zip --reference runs/9big-case/ak1690_observed.sdf --output runs/xyz-audit
python prepare_xyz_cases.py --audit runs/xyz-audit --parent-case runs/9big-case --output runs/xyz-cases
```

这是本例的原子顺序/片段分界假设，不能用于任意XYZ。明确区分：推断连通性、价态检查、固定片段匹配、退出位点保持、立体化学与坐标。氢数在固定键修复后重新计算，否则无氢XYZ的OpenBabel价态标记会丢失手性信息。

`prepare_xyz_cases.py` 导出可探索的参考、两个方向的任务分区及临时STAT6姿态。对齐来自刚性warhead模板迁移，报告RMSD与碰撞，不能替代对接或已验证的结合姿态。所有假设保存到 `provenance.json` 并由采样报告继承。

导出生成分子前，需要从3D坐标明确分配双键E/Z与原子手性，否则MolBlock可能把已有的双键几何写成 `STEREOANY`。本项目的共享导出函数已补上这一处理，并有烯烃回归测试。旧采样可用 `reexport_fixed_run.py --case ... --run ... --atom-mode basic --output ...` 从保存的NPZ重新导出；原运行不覆盖，新报告记录修正来源，不把重新导出当成新采样。

### 候选质量与重建稳健性

`assess_exploratory_candidates.py` 检查连通性、异常键长、非键原子重叠及蛋白碰撞，并报告常规分子描述符，不对PROTAC套用硬性Lipinski过滤。
对linker坐标做小幅扰动，比较OpenBabel推断的边和键级，保存不同图假设及频数。这是单一推断方法的数值稳健性，不是正确键级概率或结合成功率。
单分子MMFF94s放松在拷贝上进行，最多3000步，记录收敛与同一分子的能量下降；不覆盖原坐标，也不能把不同分子的绝对能量当作亲和力排序。

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
