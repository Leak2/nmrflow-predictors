# 评估套件（tools/）

这个目录随便携包发布。目的只有一个：让 README §5 的留出法数字从"作者自报"
变成**你能自己跑出来**的东西。

涉及的数字：¹H 0.39 / ¹³C 3.41（Tier A+B）、HOSE 自身 0.25 / 3.16、
¹⁹F 8.76（有机 C–F 7.32）/ ³¹P 19.71。

---

## 1. 先自检，再跑评估

```powershell
.\runtime\python\python.exe app\tools\check_eval_kit.py
```

9 项检查：脚本齐备、索引齐备、逐字节哈希、四核齐全、训练/杂核索引的核种隔离、
¹³C 训练/全量行数比（应约 0.82）、以及**训练索引的每个 (核, HOSE 码) 键都能在
全量索引里找到**（说明留出集确实没漏进训练索引）。

全部 PASS 再往下走；任何一项 FAIL，后面跑出来的数字都不算数。

```
  [PASS] 评估脚本齐备            13 个
  [PASS] 索引齐备              full / train / hetero
  [PASS] 哈希一致              3 条
  [PASS] 全量索引四核齐全          {'13C': 238295, '19F': 1727, '1H': 77188, '31P': 11781}
  [PASS] 训练索引只含 ¹H/¹³C     {'13C': 195316, '1H': 64223}
  [PASS] 杂核索引只含 ¹⁹F/³¹P    {'19F': 1397, '31P': 9559}
  [PASS] ¹³C 训练/全量行数比      0.820（期望 0.72~0.9，留出 20%）
  [PASS] 训练索引的键全部来自全量索引  259539 键，多出 0
  [PASS] 杂核索引的键全部来自全量索引  10956 键，多出 0
```

---

## 2. 需要你自己准备的一样东西：原始谱库

评估的输入是 NMRShiftDB2 的 `nmrshiftdb2withsignals.sd`（151 MB，约 5 万条实验谱）。

**它不随包**，原因有两条：

1. 体积是便携包的 1/4，而绝大多数用户只想要"输入 SMILES 出谱图"；
2. 它是**整个** NMRShiftDB2 数据库。再分发整库意味着这个包同时也是数据库分发者，
   要连带一整套 ODbL 要求的声明与追踪（许可允许，但我们只分发按需求解的聚合索引，
   见 `app/nmrflow/data/HOSE_INDEX_NOTICE.txt`）。

获取方式：SourceForge 的 `nmrshiftdb2/data/` 目录，无需认证，也不需要注册。
0.2.8 核对时用的是这一版，建议先比对：

```text
bytes    158,569,118
records  58,187
sha256   0E86688360E23C88CCF0EB82A1251315FA57EC6F5B376DC8886F3311793A0AFE
```

（该记录来自 2026-10-30 的一次核对；仓库没有保留原始文件，所以本机无法复核。）

需要的 Python 依赖只有 `rdkit` 与 `numpy`，都已在 `runtime\python` 里：

```
python 3.13.14 / rdkit 2026.09.1 / numpy 2.5.4
```

---

## 3. 随包的两份训练索引（这是能复算的关键）

留出法的做法是：**80% 的化合物建索引，20% 留出来评估**。如果用随包的全量索引去评估，
就是拿答案考自己——数字会漂亮得离谱。所以另带两份训练索引：

| 文件 | 字节 | 内容 | 行数 |
|---|---|---|---|
| `app/build/hose_index_train.sqlite` | 25,202,688 | ¹H + ¹³C | 259,539 |
| `app/build/hose_index_hetero.sqlite` | 1,523,712 | ¹⁹F + ³¹P | 10,956 |

sha256（与 `app/build/EVAL_KIT_SHA256.txt` 同一份记录）：

```
29738486456DC55D078B0BA0737C21E1ADD122F1EA0C6D90330354F5475F2D3A  hose_index_train.sqlite
263EB4C6519C6E833B85EA51F32884A801D34B0E976C6B89AE953E4FCBEAE02A  hose_index_hetero.sqlite
9581BD5BD9BFD9D5F2337CEC40AF2908E9E60EDBD753BB9D581CF255402ACAA6  data/hose_index.sqlite（全量，运行时用）
```

这两份就是生成 README 数字时用的那两份，逐字节相同。所以你**不需要**重建索引
（从 151 MB 的 SD 建一次约 6 分钟），直接拿它们评估即可。

### 留出法怎么划的

不是随机抽行，也不是按原子抽，而是按**化合物**：

```python
key = InChIKey(分子)                                  # 同一化合物去重
h   = blake2b(f"nmrflow-hose-v1:{key}", digest_size=8).digest()
留出 = int.from_bytes(h, "big") % 1_000_000 < 200_000   # 20%
```

`salt` 写成常数而不是随机种子，是为了让"同一份数据、同一份代码"永远划出同一批分子。
换 `--salt` 会换一批留出分子，MAE 会动零点几——这也是为什么这个字符串要写在文档里。

---

## 4. 命令 → 数字

在包根目录（`nmrflow-0.2.7-portable\`）执行，`<SD>` 换成你的 SD 路径。

### ① ¹H 0.39 / ¹³C 3.41（合并后，默认配置）

```powershell
.\runtime\python\python.exe app\tools\eval_tier_ab.py --sd <SD> --nuclei 1H,13C `
    --train-index app\build\hose_index_train.sqlite
```

看输出最后的「合并策略」表：`仅 Tier A` 一行是 0.59 / 6.94，
`≥1球` 一行是 **0.39 / 3.41**，紧跟的 P90 是 0.85 / 8.60。

### ② HOSE 自身 0.25 / 3.16，以及 1~4 层球的 MAE 表

```powershell
.\runtime\python\python.exe app\tools\eval_hose_index.py --sd <SD> --nuclei 1H,13C
```

看 `MAE by sphere` 表：1 层 0.53/0.73/6.33/8.48，逐层降到 4 层 0.10/0.12/0.92/1.27。

### ③ ¹⁹F 8.76 / ³¹P 19.71，以及兜底率与物理合理性过滤

```powershell
.\runtime\python\python.exe app\tools\eval_tier_ab.py --sd <SD> --nuclei 19F,31P `
    --train-index app\build\hose_index_hetero.sqlite
```

三张表各有用途，别只看第一张：

- 「合并策略」→ 19.71 / 8.76（
  同一张表里 `仅 Tier A` 的 40.40 / 23.62 是兜底常数造成的，不是化学预测）
- 「Tier A 走到兜底常数时…」→ ³¹P 兜底率 56%、¹⁹F 24%，以及逐球层能不能替掉兜底
  （¹⁹F 1 层球 96.06 比兜底还差，所以 `min_level=2` 的依据就在这一行）
- 「实验值的物理合理性过滤」→ `有机C–F` 一行给出 7.32（与 FluoBase 对照用的口径）

### ④ 其它诊断（都要 SD，都是"当时怎么发现问题的"）

| 脚本 | 用途 |
|---|---|
| `diag_sp2.py` | sp² 碳逐原子诊断（`--set holdout` 切留出集） |
| `diag_sp3.py` | sp³ 碳：环烷 / 烯丙位 / 苄位 / 羰基 α 位 |
| `diag_bicyclic.py` | 稠环桥环按「环体系签名」统计 Tier A vs HOSE（0.2.6 撤特例的依据） |
| `diag_hetero.py` | ¹⁹F/³¹P 离群排查：库内分布 + 最差样本的 (A, B, 实测) 三联 |
| `mine_hetero.py` | 按邻居签名挖未覆盖环境，输出「中位数 vs 兜底」的收益 |
| `ml_baseline_19f.py` | 岭回归 vs HOSE 查表的同信息量对照（只依赖 numpy） |
| `bench_spin.py` | 分块对角化的耗时/内存基准 |
| `check_reference_records.py` | 核对 SD 里十氢萘/降冰片烷的记录与文献值（顺/反异构、C7 归属），避免用错记录改模型 |
| `build_hose_index.py` | SD → 索引（`--nuclei` 选核，`--spheres 4`） |

---

## 5. 数字会飘的地方

| 来源 | 影响 |
|---|---|
| RDKit 版本（本项目 2026.09.1） | 芳香化/去重/InChIKey 生成若有变，留出批次会变，MAE 动 0.1~0.5 |
| SD 文件版本 | 库在持续增补，换一版数据 MAE 必然变 |
| `--salt` | 换划分批次，MAE 动零点几 |

同数据同版本下，两位小数应当稳定复现；**不要**指望跨 RDKit 版本逐位相同。

---

## 6. 仍然没做到的

- **SD 只有外来的校验值**：0.2.8 记下了 2026-10-30 那版的 size/records/sha256，
  但仓库没有保留原文件，无法在本机复核。
- **参考记录的判据是文献值，不是"标准答案"**：`check_reference_records.py` 用的是
  Lippmaa 1970（降冰片烷）与 trans/cis 桥头 44.0/36.0（十氢萘）。文献本身有 ±0.5
  的分布，所以容差给到 2 ppm。
- 留出法只跑了一次划分（不是 k 折），所以 3.41 本身也带划分噪声。
- `tests/` 里的 230 条 pytest 用例没有全部随包：便携 runtime 不带 pytest，
  包内只跑 `tests\run_tests.py`（stdlib unittest，20 条）。
  完整测试仍在源码仓库。
