# 固定端点：当前完成程度

现阶段完成 **9BIG 结构准备与候选化学约束检查**，尚未完成约束扩散采样，也没有新候选的活性或降解验证。

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

## 下一阶段

采样入口仍不能保持指定端点的化学结构。上游 atom mask 同时冻结类型和坐标，却未约束固定化学键；必须分别设计化学固定、出口键保持和构象采样，再在可复现的权重/运行环境中验证。DeepTernary 用于完整候选的后续三元构象评估。
