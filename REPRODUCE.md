# 可复现性说明（0.3.3）

这个便携包现在能复现**两部分**：回归行为，和 README §5 的留出法误差。
还剩下一样东西需要你自己下载（原始谱库），理由见 §2。

---

## 1. 包内直接可复现的

### 便携包回归测试（36 项，stdlib unittest，无需 pytest）

```powershell
cd nmrflow-0.3.3-portable
.\runtime\python\python.exe tests\run_tests.py
```

覆盖：

- HRMS 只有 Calculated；没有 Found；CLI 不再有 `--found` / `--seed`；界面不得出现该开关
- 乙醇 CH₂ 四重峰（¹H 线宽修复）
- 胆固醇 ¹H 无负峰、峰区间上限、积分 46.0、>10 自旋 warning
- ¹⁹F 负位移保留
- 带电物种 [M]⁺/[M]⁻、盐的对离子过滤
- 十氢萘 27.1 / 34.7 / 44.2、降冰片烷 29.6 / 36.3 / 38.3（0.2.6 起走 HOSE，
  不再有"小双环特例"）
- 磷酸三甲酯 ³¹P 保留 curated 0.0（不被浅层 HOSE 覆盖）
- 咖啡因 ¹³C 不再被 Tier A 分歧降成 low
- 水不输出空 ¹³C 块、未知溶剂返回 1
- **bat 字节约束**：全 CRLF；含中文的必须是 GBK 且不带 `chcp 65001`
  （0.2.8 用的是 LF + UTF-8，双击卸载向导会把行拆成碎片执行）
- **自检能在 GBK 控制台跑完**：`¹`(U+00B9) 在 GBK 里没有，0.2.8 上直接崩
- **评估套件本身**：`app/tools/*.py` 能 import（`--help` 跑通）、
  `check_eval_kit --quick` 九项全 PASS、以及把训练索引指向全量索引时自检必须报错

GUI 的批量模式 / 结构预览 / ZIP 下载在开发仓库的 `tests/test_app.py` 里用 Streamlit AppTest 验证；
仓库全量测试是 `pytest tests/` 的 278 条，便携包只带上面 36 条交付级回归。

### 索引完整性

```powershell
.\runtime\python\python.exe app\tools\check_eval_kit.py
```

```
  [PASS] 评估脚本齐备            12 个
  [PASS] 索引齐备              full / train / hetero
  [PASS] 哈希一致              3 条
  [PASS] 全量索引四核齐全          {'13C': 238295, '19F': 1727, '1H': 77188, '31P': 11781}
  [PASS] 训练索引只含 ¹H/¹³C     {'13C': 195316, '1H': 64223}
  [PASS] 杂核索引只含 ¹⁹F/³¹P    {'19F': 1397, '31P': 9559}
  [PASS] ¹³C 训练/全量行数比      0.820（期望 0.72~0.9，留出 20%）
  [PASS] 训练索引的键全部来自全量索引  259539 键，多出 0
  [PASS] 杂核索引的键全部来自全量索引  10956 键，多出 0
```

三份索引的摘要（同一份记录在 `app/build/EVAL_KIT_SHA256.txt`）：

```text
9581BD5BD9BFD9D5F2337CEC40AF2908E9E60EDBD753BB9D581CF255402ACAA6  nmrflow/data/hose_index.sqlite   （全量，32,538,624 B / 328,991 行）
29738486456DC55D078B0BA0737C21E1ADD122F1EA0C6D90330354F5475F2D3A  build/hose_index_train.sqlite   （¹H+¹³C，25,202,688 B / 259,539 行）
263EB4C6519C6E833B85EA51F32884A801D34B0E976C6B89AE953E4FCBEAE02A  build/hose_index_hetero.sqlite  （¹⁹F+³¹P，1,523,712 B / 10,956 行）
```

索引只存 `(nucleus, code, n, Σδ, Σδ²)`，不存原始结构/谱图；许可与来源声明见
`nmrflow/data/nmrshiftdb2_database_license.txt` 与 `HOSE_INDEX_NOTICE.txt`。

### 基本自检

```powershell
.\runtime\python\python.exe -c "import nmrflow; print(nmrflow.__version__)"
.\runtime\python\python.exe -m nmrflow.cli "CCO" -o out
# 期望：HRMS 只有 Calculated；没有 Found
.\runtime\python\python.exe -m nmrflow.cli "CCO" -o out
# 期望：report.txt 末尾有 Predicted by nmrflow ... - not experimental data.
```

---

## 2. 留出法误差：包内可复算，但要你补一个 SD 文件

```powershell
# ① 自检（不用 SD）
.\runtime\python\python.exe app\tools\check_eval_kit.py

# ② ¹H 0.39 / ¹³C 3.41
.\runtime\python\python.exe app\tools\eval_tier_ab.py --sd <SD> --nuclei 1H,13C `
    --train-index app\build\hose_index_train.sqlite

# ③ ¹⁹F 8.76 / ³¹P 19.71
.\runtime\python\python.exe app\tools\eval_tier_ab.py --sd <SD> --nuclei 19F,31P `
    --train-index app\build\hose_index_hetero.sqlite

# ④ HOSE 自身 0.25 / 3.16
.\runtime\python\python.exe app\tools\eval_hose_index.py --sd <SD> --nuclei 1H,13C
```

`<SD>` = NMRShiftDB2 的 `nmrshiftdb2withsignals.sd`，SourceForge 免费下载，不需要注册。
0.2.8 核对时用的是这一版，**建议先比对校验值**（这也是 0.2.7 遗留的那个缺口）：

```text
bytes    158,569,118
records  58,187 raw / 56,786 with signals
sha256   0E86688360E23C88CCF0EB82A1251315FA57EC6F5B376DC8886F3311793A0AFE
```

（该记录来自 2026-10-30 的一次完整核对；本仓库没有保留原始文件，所以无法在本机复核。
但两份**训练索引**仍然逐字节随包，`check_eval_kit.py` 会核对它们，不受 SD 版本影响。）

它不随包：体积与便携包同量级，而它是**整个**谱库，再分发整库意味着这个包同时也是
数据库分发者（许可允许，但要连带全套声明）。我们只分发按需求解的聚合索引。

命令与数字的逐条对应、留出法的确切定义、以及三个会让数字飘的变量
（RDKit 版本 / SD 版本 / `--salt`），都在 **`app\tools\README.md`**。

### 剩下的缺口（如实写）

- **SD 的校验值有了，但仍无法在本机复核**：0.2.8 记录下了 2026-10-30 那版的
  size / records / sha256（见上），数据本身在持续增补，换一版 MAE 必然变。
- **十氢萘/降冰片烷的参考记录**：SD 里同一骨架有多条记录（cis/trans 共用
  `InChIKey` 前 14 位）。`app\tools\check_reference_records.py` 会把它们与文献值
  逐条比对，跑一次就知道手上的 SD 是否可用。
- 只跑了一次 20% 划分，不是 k 折；3.41 这个数本身带划分噪声。
- 仓库全量测试是 `pytest tests/` 的 278 条；便携 runtime 不带 pytest，包内只跑
  36 条 stdlib unittest，两者分工明确。

---

## 3. 0.3.3 之后仍然已知的薄弱区

- **HOSE 在纯烃稠环上是 ±2 ppm 的量级，不是 0.5**。十氢萘现在 42.5 / 32.6 / 25.9
  对实测 44.2 / 34.7 / 27.1；这属于"没有专门的补丁"，不是"算错了"。
  想要纯规则路径用 `--no-tier-b`。
- **大于 10 个自旋的强耦合块**仍退回一级近似，只在 `report.warnings` 里提示。
- **³¹P** 的 TMP 负收益已用 curated 保护修掉，但整体 MAE 19.71 / P90 51
  没有变；它仍是四核里最弱的一个（兜底率 56%）。
- **¹⁹F/³¹P 尾部**（P90 21 / 27 ppm，最大数百 ppm）没有变。
- **稠合杂芳环/大环**仍走粗缩放；HOSE 的环闭合信息缺失是结构性问题，
  真正的解法是重建带环闭合扩展的索引，不是再加常数。
- **置信度**是"球层 + n=1 折扣 + σ 降级"，比 0.2.2 合理，但仍不是逐分子的误差棒。
- **exe 仍然没有代码签名**（没有证书），SmartScreen 照样会拦一次；但 0.2.7 起
  有版本资源了：右键 → 属性 → 详细信息能看到 `0.3.3.0`，`NMR-Predictor.exe --version`
  也会弹窗报版本（窗体程序没有 stdout，所以是弹窗而不是打印）。

---

*0.3.3｜2026-11-13*
