# 输入表格与输出

默认输入文件名：`e3_ligand.xlsx`、`warhead.xlsx`、`protac.xlsx`、`protacDB_smiles.csv`。
通过 `--source-map sources.json` 可以映射任意文件名，不依赖聊天软件生成的哈希名称。四种角色必须齐全，CSV 和 XLSX 均可读取。

```json
{"e3_ligand":"e3.xlsx","warhead":"warhead.xlsx","protac":"protac.xlsx","paired":"pairs.csv"}
```

前三份表包含 `Compound ID`, `Uniprot`, `Target`, `Smiles`；完整 PROTAC 表可额外提供 `E3 ligase`。
配对表包含 `id_protac`, `linker_canonical`, `warhead_canonical`, `e3_ligand_canonical`, `smiles_canonical`。
配对片段目前按一个连接标记处理；多标记片段会被隔离，但不据此判断化学结构无效。linker 原字符串保留，并未验证整个分子可由三个片段正确拼接。

输出：

| 文件 | 内容 |
|---|---|
| fragments.csv | 标准化片段、连接标记、母体、分组、记录 ID |
| provenance.jsonl | 通过记录 ID 对应到源文件、行号和靶点 |
| paired_protacs.csv | 配对记录与 valid/quarantined 状态 |
| validation_errors.csv | 具体失败字段、原值和原因 |
| summary.json | 数量、源文件摘要和方法说明 |

同一母体的不同连接位置落在同一个训练/验证/测试集合；比例是确定性哈希的约 80/10/10，不强制小样本满足比例。配对表没有用于完整 PROTAC 模型的独立划分。记录 ID 不能代替活性或三维结构标识。

个人数据放在 Git 忽略的 `data/`，处理输出建议放 `runs/`。将同学的数据通过适当渠道单独共享，不需要修改源码中的文件名。

## 用自己的文件运行

文件名可以不改。在项目根目录新建 `sources.json`，填入上面的四项映射，再运行：

```bash
python prepare_fragments.py --source-dir data/raw --source-map sources.json --output-dir runs/my-data-01
python verify_outputs.py --data-dir runs/my-data-01
```

映射中的文件路径相对于 `--source-dir`；映射文件本身的位置不影响它。XLSX 读取当前活动工作表，第一行必须是列名。CSV 使用 UTF-8（支持 BOM）。列名区分大小写，例如 `Smiles` 和 `smiles` 不相同。

空行会跳过，但来源行号仍对应原文件：Excel 是工作表行号，CSV 是该记录开始的物理行号。重复列名、缺少必需列、空表，以及 CSV 一行的数据多于表头列数，都会停止处理并指出文件。

## 怎样看错误和配对状态

`validation_errors.csv` 记录的是字段错误，所以同一行可能出现多次。它包含文件、行号、列名、原始字符串和原因。先修正原表，再用新的输出目录重跑；不要只改生成的 CSV，否则来源和摘要会失去对应关系。

`paired_protacs.csv` 中的 `valid` 只表示两端片段和完整分子的 SMILES 通过当前检查，且两端各有一个连接标记。linker 只保留原字符串，尚未检查三部分能否拼回完整分子；因此 `valid` 不等于可训练、可合成或有活性。`quarantined` 表示至少一个必需结构未通过检查，原记录仍保留。

部分字段出错时，其他通过检查的片段仍可能进入 `fragments.csv`，即使它所在的配对记录被隔离。这适合分别整理两端片段；如果要训练成对模型，必须另行筛选配对，并为配对数据设计划分。
