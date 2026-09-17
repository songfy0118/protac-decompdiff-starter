# 完整模型环境与运行

本仓库的支持范围是数据准备和两个推理入口。上游的完整评估脚本需要额外的 docking 工具和评分资源；可选 SA 评分文件 `DecompDiff-main/utils/evaluation/fpscores.pkl.gz` 未随此仓库发布，需要评估时从原作者仓库获取。这里的分子导出和一致性核验不依赖该文件。

## 安装

Python 3.10。以下命令在项目根目录运行；与轻量数据环境分开安装也可以。

```bash
python -m venv .venv-model
# Windows: .venv-model\Scripts\Activate.ps1
# Linux: source .venv-model/bin/activate
python -m pip install --upgrade pip
python -m pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu118
python -m pip install -r requirements.txt
python -m pip install torch-scatter torch-sparse torch-cluster --only-binary=:all: -f https://data.pyg.org/whl/torch-2.7.0+cu118.html
python doctor.py --inference
python -m unittest -v test_pipeline
```

先升级 pip，旧版本可能误判官方 wheel 的包名元数据。扩展 wheel 必须匹配 PyTorch/CUDA/Python；不要混用不同构建。CPU 推理可传 `--device cpu`，但未做性能验证；macOS 模型环境未验证。

## 权重

从[作者提供的权重目录](https://drive.google.com/drive/folders/1JAB5pp25rEM5Wt-i373_rrAyTsLvAACZ)获取 checkpoint，保存为 `assets/decompdiff.pt`。
本项目测试过的权重 SHA-256 为 `9d3daf75ff37fbc314f2dd43f899c7c104a8f68992cad5c0e8446574fac39c17`。其他权重需重新验证兼容性。

```bash
python run_endpoint.py --spec examples/test_spec.json --checkpoint assets/decompdiff.pt --output runs/endpoint-001 --samples 1
python verify_outputs.py --endpoint runs/endpoint-001
```

输出包括 `candidates.sdf`、`report.json` 和 `run_manifest.json`。失败和中断保留记录；非空输出目录不会被覆盖。
示例结构来自官方测试编号 0，而非特定 PROTAC 项目靶点。

## 自定义结构

复制 `examples/test_spec.json`：`protein_pdb` 和 `reference_sdf` 相对于 JSON 所在目录；`atom_groups` 按 SDF 去氢后的原子顺序填写。`-1` 是 scaffold，`0..N-1` 是 arms，每组至少两个原子。保留参考坐标与蛋白的对齐关系。
角色可为 `warhead` / `e3_ligand`，但这个字段仅标记结果，未改变预训练模型的条件输入。生成并不自动保持 linker 出口位点。

## 官方基线（可选，大数据）

从[作者的预处理数据目录](https://drive.google.com/drive/folders/1z74dKcDKQbwpo8Uf8EJpGi12T4GCD8_Z)获取匹配的分解 LMDB、name2id 和 split。测试结构/索引按上游 README 获取。

```text
assets/
  decompdiff.pt
  data/
    crossdocked_v1.1_rmsd1.0_processed_full_ref_prior_aromatic.lmdb
    crossdocked_v1.1_rmsd1.0_processed_full_ref_prior_aromatic_name2id.pt
    split_by_name.pt
    test_index.pkl
    test_set/...
```

```bash
python run_baseline.py --assets assets --output runs/baseline-001 --samples 1
python export_baseline.py --run-dir runs/baseline-001
python verify_outputs.py --baseline runs/baseline-001
```

`crossdocked_pocket10_mol.lmdb` 是另一种预处理格式，不能直接改名替代本项目所需 LMDB。基线依赖可信的作者序列化数据；`export_baseline.py` 只用于自己刚生成的 `.pt`，不要加载来历不明的 pickle 文件。
