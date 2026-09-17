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
