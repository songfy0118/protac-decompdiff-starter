# PROTAC DecompDiff Starter

基于 [DecompDiff](https://github.com/bytedance/DecompDiff) 的可复现项目起点：整理 PROTAC 两端片段数据、运行蛋白口袋条件生成、导出并核验分子结果。

**这是工程复现与数据准备工具，不是已经训练好的 PROTAC 专用模型。** `warhead` / `e3_ligand` 是任务标签；当前网络仍是作者的预训练 DecompDiff。尚未实现 linker 连接位点保持、PROTAC 微调或活性验证。

## 从哪里开始

| 你想做什么 | 要求 | 入口 |
|---|---|---|
| 先确认项目能跑 | Python 3.10，无需 GPU | 下方五分钟演示 |
| 整理自己的表格 | Excel/CSV | `prepare_fragments.py` |
| 生成一个口袋中的候选配体 | 模型依赖、权重、对齐的 PDB/SDF | `run_endpoint.py` |
| 复现官方参考先验流程 | 完整模型环境和官方分解数据 | [基线说明](docs/INFERENCE.md) |

## 五分钟数据演示（无需 GPU）

从仓库根目录执行。Python 3.10 已测试。

```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# Linux/macOS 使用: source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-data.txt
python doctor.py
python -m unittest -v test_pipeline
python prepare_fragments.py --source-dir examples/table_inputs --source-map examples/table_inputs/sources.json --output-dir runs/demo
python verify_outputs.py --data-dir runs/demo
```

不便激活环境时，Windows 可将 `python` 换成 `.venv\Scripts\python.exe`。
无 PyTorch 的环境会明确跳过一个 `.pt` 导出测试，其余检查照常执行。

演示输入是人为构造的玩具分子，**不是已知活性 PROTAC**。预期输出为 4 条片段记录、2 个母体分组、1 条配对记录、0 条隔离记录；完整研究数据没有随仓库上传。
重复运行请换一个新的 `--output-dir`，已有非空目录不会被覆盖。

## 数据与结果

`prepare_fragments.py` 接受标准列名的 Excel/CSV，将 `[R1]` / `[R2]` 转成保留映射的 `[*:1]` / `[*:2]`，按封端母体进行确定性划分，保留来源和异常记录。[输入格式](docs/DATA.md)

`sample_results/` 保存两次实际 GPU 运行的输出及摘要。每次使用一个样本、完整 1,000 步，均得到一个可读取的连通分子。这仅验证运行链路，不是质量或成功率基准。SDF 的手性由生成的三维坐标推断，未做实验验证。

![实际生成示例](sample_results/baseline/generated.png)

```bash
python verify_outputs.py
```

该命令核验示例 SDF、SMILES、数量、连通性及运行清单；加入 `--data-dir` 才会检查处理后的表格，未检查项不会标为通过。

## 完整模型运行

见 [安装与推理指南](docs/INFERENCE.md)。权重和大型 LMDB 不在 Git 中，下载位置、期望文件名和目录结构均列在指南里。

```bash
python run_endpoint.py --spec examples/test_spec.json --checkpoint assets/decompdiff.pt --output runs/endpoint-001 --samples 1
```

自定义输入需要蛋白 PDB、处于同一坐标系的参考配体 SDF，以及明确的原子分组。不能把 SMILES 的任意三维构象直接放进蛋白口袋。`arms/scaffold` 不等同于 `warhead/E3`。

## 项目结构

- `DecompDiff-main/`：保留上游实现及少量兼容性修改。
- `prepare_fragments.py`：输入校验、标准化、去重、划分和来源记录。
- `run_endpoint.py` / `run_baseline.py`：自定义口袋与官方基线入口。
- `run_checks.py`：运行记录、输入检查和一致的分子导出。
- `verify_outputs.py`：数据与结果的交叉核验。
- `doctor.py`：环境检查；`test_pipeline.py`：回归测试。
- `package_project.py`：可选打包及逐文件哈希校验。

## 研究边界与后续工作

下一步需明确一个靶蛋白和一种 E3、参考结合结构及 linker 连接位点，再构建对应的三维训练/评估数据、实现位点保持和微调。现有代码没有声称这些步骤已经完成。
数据分割目前按母体分子分组，不是按骨架或靶点留出，不能由此主张跨骨架或跨靶点泛化。

## 来源与许可

原论文：[DecompDiff: Diffusion Models with Decomposed Priors for Structure-Based Drug Design](https://proceedings.mlr.press/v202/guan23a.html)，Guan 等，ICML 2023。
上游以 CC BY-NC 4.0 发布；本仓库保留原许可和文件级声明，见 [LICENSE](LICENSE) 与 [第三方说明](THIRD_PARTY_NOTICES.md)。本项目是改编和复现工具，不代表原作者背书。
