# nmrflow CHANGELOG

## 0.3.0 (2026-11-03) — 防误用收口 + 发布元数据

这是 0.2.9 之后的**行为版本**：应用、exe、包名统一为 0.3.0；核心预测逻辑不变。

### 防误用

- **删除模拟 Found**：CLI `--found` / `--no-found` / `--seed`、`RunOptions(ms_found_mode=...)`、
  `AdductResult.found` / `found_is_simulated`、`SIMULATED_FOUND_NOTE` 全部移除；HRMS 只给 Calculated。
- **默认预测页脚**：所有 `report.txt` 末尾带
  `Predicted by nmrflow 0.3.0 - not experimental data.`
- **PNG 水印/metadata**：右下角 `PREDICTED - not experimental`；PNG 写
  `Software=nmrflow 0.3.0`、`Description=Predicted NMR spectrum - not experimental data`。
- **JSON provenance**：顶层 `schema_version` / `data_type` / `generator` / `generated_at` /
  `input` / `options`。
- **manifest.json**：CLI / GUI / batch 每次运行都写，含输出 SHA256。
- **ETHICS.md**：预测≠实验、禁止用途、合法用途、输出标记、代码/数据许可区别。

### 发布元数据

- **`CITATION.cff`** + **`.zenodo.json`**：Zenodo 归档后把 DOI 填进 CITATION。
- **README 许可 FAQ**：代码 BSD-3-Clause；HOSE 索引是 NMRShiftDB2 衍生数据，按
  nmrshiftdb2 Database License 分发。
- **README 示例图**：`assets/example-predicted-13c.png`，右下角带 `PREDICTED` 水印。
- **版本**：应用 `0.3.0`、exe `0.3.0.0`、包名 `nmrflow-0.3.0-portable`。
- **测试**：仓库 258 passed；包内 **36 项**。

## 0.2.9 (2026-11-03) — 修 0.2.8 的两个交付缺陷（bat 换行 / 自检编码）

0.2.8 的**功能**没问题，但**交付物**有两个只有用户双击才会发现的坑。
2026-11-02 的外部交付（`update_26.11.03`）把它们指了出来，本版修掉，并补上能拦住它们的检查。

### 1. 所有 .bat 是 LF-only（严重：卸载向导不可用）

`make_release` 把 bat 当普通文本写（`text.encode()`），于是换行全是 LF。
cmd.exe 读 LF-only 的批处理时，一遇到多行 `if (...)`、标签、`%~1` 子串就把行拼错：
实测 `卸载.bat` 被拆成 `'cho'`、`'HOICE"=="3" goto cancel'` 这样的碎片**去执行**，
窗口刷一堆 `不是内部或外部命令`，选 2 完全没反应。

修法不是改某一个文件，而是把生成器改对：新增 `bat_bytes()`，写出前统一转 CRLF
（`\r\n`），任何 bat 都走这条路径。

### 2. 含中文的 bat 还必须是 GBK，且不能带 `chcp 65001`

外部版本只改了换行，我这边仍能复现两处错位（`'会' is not recognized`）。
原因是 cmd 按**当前代码页**读批处理，而 `chcp 65001` 是第 2 行才生效：
之前读进缓冲区的 UTF-8 汉字按 CP936 解码，三字节字符被当成 1.5 个字符，行从汉字中间断。

三组实测（同一份 `卸载.bat 3`，输出重定向到文件后按 GBK 解码）：

| 编码 / 换行 | 解析错误 |
|---|---|
| UTF-8 + LF | 7 处以上碎片 |
| UTF-8 + CRLF + `chcp 65001` | 仍有 2 处 |
| **GBK + CRLF，无 chcp** | **0 处** |

所以：含非 ASCII 的 bat 存 GBK 并去掉 `chcp`；纯 ASCII 的（启动/停止/自检/安装）
保留 `chcp 65001`——那是给下面 Python 的输出用的，没有解码风险。
`自检.bat` 另加 `set PYTHONUTF8=1`。

### 3. `check_eval_kit.py` 在 GBK 控制台直接崩

`¹H / ¹³C` 里的 `¹`(U+00B9) 在 GBK 里没有，`print` 抛
`UnicodeEncodeError: 'gbk' codec can't encode character '\xb9'`——
用户双击 `自检.bat` 看到的第一件事就是崩溃。

修法：核名一律用 ASCII（`1H / 13C / 19F / 31P`），并给 stdout/stderr 加
`reconfigure(errors="backslashreplace")` 兜底。诊断工具里可读性优先于排版。

### 4. 打包环节加字节级检查（防第三次）

`make_release.check_bat_files()`：每个 bat 必须全 CRLF；非 ASCII 的必须是 GBK
（能按 GBK 解、**不能**按 UTF-8 解）、且不含 `chcp 65001`。不合格直接中止打包。
`tests/test_regressions.py` 增加同名的两条断言，包内跑测试就能发现。

### 5. 顺带

- **zip 压缩级别 1 → 6**：条目数不变（12,935），包从 246 MB 降到约 220 MB。
- **版本号 0.2.8 → 0.2.9 并重新冻结 exe**（`0.2.9.0`）。
  0.2.8 已经发过一版有问题的包，再叫 0.2.8 会让"手上这份到底是哪版"无法判断。
- UI 右上角版本角标随之显示 `v0.2.9`（沿用外部交付的实现）。
- 测试：仓库 **258 passed**，包内 **33 项**。
- **删除模拟 Found**：CLI `--found` / `--no-found` / `--seed`、`RunOptions(ms_found_mode=...)`、`AdductResult.found` / `found_is_simulated` 全部移除；HRMS 只给 Calculated。
- **防误用 P0**：`report.txt` 默认带 `Predicted by nmrflow ... - not experimental data.` 页脚；PNG 加 `PREDICTED - not experimental` 水印与 `Software`/`Description` 元数据；JSON 顶层加 `schema_version` / `data_type` / `generator` / `generated_at` / `input` / `options`；新增 3 条回归测试。
- **P1 manifest / ETHICS**：每次运行写 `manifest.json`（输入、选项、输出 SHA256）；新增 `ETHICS.md` 和 README 使用边界；新增 manifest/`--found` 拒绝/ETHICS 回归。

## 0.2.8 (2026-11-02) — 吸收外部交付的界面/脚本改进，拒绝「恢复小双环守卫」

2026-10-30 收到一份外部交付（0.2.8 / 0.2.8-r2 / 0.2.9 三个包 + 补丁）。逐项复核后
**拆成两半处理**：界面与脚本的改进照收，位移层的改动退回。

### 采纳

- **关闭所有页面后自动停止**（`nmrflow/autostop.py`）：每个会话每 30 s 心跳，
  守护线程在最后心跳超过 **180 s** 后 `os._exit(0)`。侧栏可勾选，**默认关闭**；
  不用本地 HTTP 端点，因此没有"别的网页能关你的服务"的攻击面。
  代价写在文档里：刷新/切标签/休眠不会误停，但关页面 3 分钟后才退出，
  正在跑批量任务时可能被中断。
- **`停止 NMR.bat` 重写**：可带端口参数、`/quiet`（供卸载向导调用）、
  `taskkill /T` 带走子进程、结束后复查端口是否真的释放、失败返回 1。
- **卸载向导**（`卸载.bat`）：三选一。选 2 时先把 `results` 移到桌面、再停服务、
  最后删除目录；**服务停不掉就中止删除**，不留删一半的目录。
- **`自检.bat`**：双击跑 `check_eval_kit.py`，用户自己能确认包没坏。
- **SD 的校验值**：`158,569,118 B / 58,187 条 / sha256 0E86…0AFE`。
  这正好补上 0.2.7 REPRODUCE 里承认的缺口。

### 拒绝：恢复 (5,5)/(6,6) 窄守卫

外部版本认为「0.2.6 撤特例的依据是 SD 错记录」，于是把特例加回来并收窄到 (5,5)/(6,6)。
核对文献后**没有采纳**，理由是那份"文献值"本身错了：

| 分子 | 文献值 | 外部版本用的"文献值" |
|---|---|---|
| 降冰片烷 | C7 = **39.2**（全分子最去屏蔽） | C7 = 24.5（弄反了） |
| 十氢萘 trans | 桥头 **44.0**（cis 才是 36.0，差 8 ppm） | 「44.2 不合理」 |

两处根因：① 把 `InChIKey` 前 14 位（骨架哈希）当成完整键，cis/trans 十氢萘本就共用它，
于是"同一分子两条记录冲突"被当成脏数据；② 拿 24.5 当 C7 的文献值，而 C7 恰恰是
降冰片烷里最去屏蔽的碳。实测对比（同 SMILES，默认路径）：

| 分子 | 文献 | 本版（HOSE） | 外部 0.2.9（窄守卫） |
|---|---|---|---|
| 十氢萘 桥头 | 44.0 | 42.5（1.5） | 34.7（**9.3**） |
| 十氢萘 CH₂ | 34.7 | 32.6（2.1） | 26.4（**8.3**） |
| 降冰片烷 C7 | 39.2 | 38.2（1.0） | 24.9（**14.3**） |

外部版本同时把 ¹³C 合并 MAE 从 3.41 推到 3.46 并解释为"正确参考的代价"——方向读反了。
本版保持 3.41 / P90 8.60。

### 新增与修正

- **`tools/check_reference_records.py`**（自写，未采用外部的版本）：把 SD 里十氢萘/
  降冰片烷的每条记录聚成 `(位移, 碳数)`，与文献值逐项比对，判定 **trans / cis / 无法解释**。
  外部版本把 `34.5/26.9/26.3` 和 `36.8/28.5/24.5` 写成期望值，会对正确的 SD **永远报 FAIL**——
  一个反向的守卫，比没有守卫更糟。本版把文献出处写进 docstring，容差 2 ppm。
- 卸载向导增加**安全闸**：必须先看到 `NMR-Predictor.exe` 与 `app\nmrflow\app.py`
  才允许删除；拒绝在磁盘根目录执行。外部版本直接 `rmdir /s /q "%~dp0"`，
  单独把这个 bat 拷到别的文件夹双击就会删掉无关目录。
- `停止 / 卸载 / 自检` 三个 bat 改为 `packaging/bat/` 下的文件模板（`%` 与引号太多，
  内联成 Python 字符串无法维护），打包时只做 `{{PORT}}` 替换。
- 测试回退 `tempfile.mkdtemp()`（外部版本改成固定 `tests/.work`，包内原地跑会
  PermissionError，还会留下目录）。
- **exe 重新冻结**：版本资源 `0.2.8.0`，与应用版本、包名、CHANGELOG 对齐
  （外部三个包的 exe 都还停在 0.2.7.0）。
- 测试：仓库 **254 passed**，包内 **27 项**（新增 4 条 autostop、3 条参考记录）。

## 0.2.7 (2026-10-28) — 评估套件随包（可复核）+ exe 版本资源

评审原文（0.2.1 快照）：「0.39 / 3.41 / 7.32 等留出法数字……包里没有 `tools/eval_*.py`，
也没有原始 SD 文件或索引构建参数……**无法从交付物复核**。」这一版把这条补上。

- **`tools/` 随包**，放在 `app/tools/`（12 个文件）。放这里而不是包根目录，是因为
  脚本开头就是 `sys.path.insert(0, Path(__file__).parent.parent)`——只有让 `tools/`
  与 `nmrflow/` 同级，`import nmrflow` 才成立。好处是**源码仓库与便携包共用同一份
  文件，不给随包版本打补丁**。
- **两份训练索引随包**：`app/build/hose_index_train.sqlite`（¹H+¹³C，25.2 MB）与
  `hose_index_hetero.sqlite`（¹⁹F+³¹P，1.5 MB）。运行时那份是全量索引，拿它评估
  等于拿答案考自己；这两份是留出 20% 之后的 80%，与出数字时逐字节相同，
  也省掉从 151 MB 的 SD 重建的 6 分钟。
- **新增 `tools/check_eval_kit.py`**（stdlib，不依赖 rdkit）。九项检查，把"文件在"
  升级成"文件是那一份"：逐字节哈希、四核齐全、训练索引的核种隔离、¹³C 训练/全量
  行数比 0.82、以及**训练索引的每个 (核, HOSE 码) 键都能在全量索引里找到**
  （证明留出集确实没漏进训练索引）。任何一项不过，后面跑出来的数字都不算数。
- **新增 `tools/README.md`**：命令 ↔ README 数字的对应表（含 `--nuclei`、
  `--train-index` 该指哪个文件）、留出法的确切定义（`blake2b("nmrflow-hose-v1:"+InChIKey)`）、
  以及会说谎的三个地方（RDKit 版本、SD 版本、`--salt`）。
- `SHA256SUMS.txt` 增列两份训练索引；`make_release.py` 新增 `copy_eval_kit()`，
  随包清单写死在 `EVAL_TOOLS` / `EVAL_INDICES` 里，`check_layout()` 逐项核对。
- 回归从 17 条扩到 **20 条**：`--help` 能跑通（验证 `app/` 的路径假设）、
  `check_eval_kit --quick` 全 PASS、以及把训练索引指向全量时**自检必须报错**
  （守住"守卫本身不是走过场"）。
- 仍然没做到：**SD 文件没有校验值**。当初那份 151 MB 文件没留 sha256，
  无法确认复算者下到的是同一版——这一条写在 `tools/README.md` §6 与 REPRODUCE.md，
  没有藏起来。

**exe 版本资源（评审 P2 遗留）**

- `NMR-Predictor.exe` 以前没有 VS_VERSIONINFO：右键属性是空白，任务管理器只显示文件名，
  用户报问题时"你用的哪一版"只能靠猜。
- 新增 `tools/version_info.py`：从 `nmrflow/__init__.py` 读 `__version__`（单一口径），
  生成 PyInstaller 的 `--version-file`。带 zh-CN(0x0804) 与 en-US(0x0409) 两个语言块，
  属性页显示 `FileVersion 0.2.7.0` / `ProductName nmrflow` /
  `FileDescription NMR 谱图预测 (SMILES → NMR + HRMS)` / `LegalCopyright BSD 3-Clause License`。
- 两个容易踩的坑写在那个文件的 docstring 里：① PyInstaller 读该文件的方式变过，
  非 ASCII 直接写有被按 cp1252/cp936 解码的风险，所以文件**保持纯 ASCII**、
  中文一律写成 `\uXXXX`；② `--version-file` / `--icon` 的路径是**相对 `--specpath`**
  解析的，传相对路径会拼成 `build/pyinstaller/build/pyinstaller/...`。
- 打包时新增两道核对：`make_release.VERSION` 必须等于 `nmrflow.__version__`；
  exe 产出后用 `pefile` 把版本资源读回来比对（`--version-file` 拼错时 PyInstaller
  可能只警告不报错，产出照样能跑的"无版本" exe）。
- `--skip-exe` 现在会**校验被回收的 exe 的版本号**：沿用上一版 exe 会让属性页写着旧版本，
  而包名/CHANGELOG/SHA256SUMS 都是新的——比没有版本号更糟。对不上就要求重新冻结。
- 启动器新增 `--version`（弹窗，因为窗体程序没有 stdout）与 `--selftest` 里的
  `version` / `built` 两行；支持排错时先问版本。

## 0.2.6 (2026-10-28) — 界面修复 + 评审遗留项清理

- **修「点一次下载就跳回首页、第二张图下不了」**。`st.download_button` 被点击同样会
  触发一次脚本重跑，而 `st.button` 在重跑时返回 `False`；此前结果只活在
  「点了生成」那一支里，于是任何一次重跑（下载、打开目录、改任意控件）都会把
  整屏结果丢掉，页面上只剩「输入 SMILES 后点生成」的提示。
  现在生成结果与已算好的下载内容一起存进 `st.session_state`，重跑只做渲染。
- 顺带的好处：ZIP 与 report 的字节只在生成时算一次，重跑不再重复打包；
  下载到的内容也不再随重跑变化。
- 结果区顶部新增「当前显示：<SMILES> · 溶剂 · 频率（生成于 hh:mm:ss）」标注。
  持久化之后「输入框改了、下面还是旧图」是正常状态，不标注就成了误导。
- 结果区顶部新增「清空结果」按钮。放在主区而非侧栏：侧栏先于生成渲染，
  放侧栏会等到下一次交互才出现。
- 批量模式同样常驻；「超过 30 个只取前 30」的提示随结果一起保留，不再重跑后消失。
- **「示例分子」下拉默认改为「（自定义）」**。原先默认选中「乙酸乙酯」，而输入框里
  显示的正是它的 SMILES——下拉与输入框看起来对不上，用户容易以为选错了。
  改为默认「（自定义）」+ 输入框留空 + 占位提示 `例如 CCO、c1ccccc1`。
  同时修掉这条路径上的一个坑：旧写法 `st.text_input(value=示例的 SMILES)`，
  选「（自定义）」时 `value` 变成 `"CCO"`，Streamlit 会重建控件、把用户手打的 SMILES
  清掉。改成显式 `key` + 只在「换了另一个示例」时写入。
- 输入框为空时点生成，给的是「请先输入 SMILES」的提示，而不是一句 SMILES 解析失败。
- 测试：`tests/test_app.py` 从 2 条扩到 11 条，覆盖连续重跑、换分子替换、清空、
  非法 SMILES 不残留、批量常驻、下拉默认值、示例填充与自定义保留手输内容、
  以及"界面不得再出现模拟 Found 开关"。
- **「模拟 Found」从界面撤下（功能保留）**。它是内测期的排版演示开关
  （Found = Calculated + 随机抖动）。正式版不再在侧栏出现，界面固定走
  `MS_FOUND_MODE = "none"`，只给 Calculated。功能本身保留在 CLI `--found`
  与 `RunOptions(ms_found_mode="simulate")`，导出物照旧带 `[SIMULATED]` 标注，
  所以即便绕道开启也不会被误当成实测值。加了一条护栏测试钉住"界面不能再出现该开关"。
  （0.2.9 已彻底删除 `--found`、`ms_found_mode` 与所有 Found 字段。）
- README 里三处写着 `Found` 的样例输出已改回默认形态（只有 Calculated），
  并在格式说明、精度表、便携包《安装使用说明》里同步注明"默认不给 Found"。

**稠环/桥环：撤掉「小双环特例」（评审 P1-1 的收尾）**

- 0.2.3 加过一条特例：判定「恰好 2 个环 + 环原子 ≤10 + 有桥头」就放弃 HOSE、
  改用经验规则 + 三个硬编码常数（`(6,6)桥头 −7.0`、`(5,5)桥头 −2.2`、`(5,5)桥碳 −9.0`）。
  本轮整段删除（`shift_hose` 的 ¹³C 拒绝、`shift_13c` 的桥头补丁、`pipeline` 的提示、
  `prepare.is_small_aliphatic_bicycle`）。
- **原因是标定用的参考值错了**：旧常数照「34.5 / 26.9 / 26.3」标定，而 NMRShiftDB2 里
  十氢萘的两条独立谱都是 **27.1 / 34.7 / 44.2**。照错值改出来的 34.7/26.4/25.2
  反而比 HOSE 原始预测 42.5/32.6/25.9 更远离实测（9.5 vs 2.1）。
  降冰片烷同理：旧补丁 36.8/27.4/24.9 对实测 38.3/36.3/29.6 差 8.9，
  HOSE 给 38.2/36.3/29.6 几乎完全对上。
- 量化（留出集，262 个脂肪族多环分子 / 2452 个 ¹³C 信号）：
  旧特例 MAE 7.19 / P90 16.00，关掉补丁 7.15 / 15.47，**放开 HOSE 4.07 / 9.40**。
  逐分子：227 变好、24 变差（2~5 ppm）、11 持平，平均 +3.16 ppm。
- 全库不变：仅 Tier A 6.94 / 合并 3.41（P90 8.60），与上一版一致——特例只影响那 262 个分子。
- 新增 `tools/diag_bicyclic.py`：按「环体系签名」统计多环体系的 A/B 误差，
  以后要再动这类规则，先用它出数。
- 回归测试改为断言**库内实测值**（十氢萘 27.1/34.7/44.2、降冰片烷 29.6/36.3/38.3），
  并检查这些分子现在走 Tier B。原来的测试断言的正是那组错的参考值。

**清理评审遗留项（第二批）**

- **删除 `cn-pdf` 排版档**。它是为一次性对照交付做的（复刻中文 PDF 复制稿的全角 ＝、
  整数与小数点间插空格、U+2011 连字符），实际使用频率极低，却会让输出出现
  `3 .64`、`J＝` 这类东西，容易被误当成产品缺陷。整档移除：`RunOptions.style`、
  CLI `--style`、界面「报道值排版」单选、`formatter` 里的 `pdf` 分支与对应测试。
- **`M+.` → `M⁺·`**（EI 的自由基阳离子）。原先输出 `... 46.0413 M+.`，既不是
  `M⁺·` 也不是 `[M]+`。`M⁺·` 自带结尾符号，所以行尾不再叠句点。
  （`HRMS (EI)` 这个名字保留：高分辨 EI 在文献里就是这么写的。）
- **CLI 部分失败返回非零**。原先 `return 0 if ok else 1`，批量里 3 个成功 1 个失败
  仍返回 0，脚本会拿一份不完整的结果继续跑。现在有任何失败就返回 1，
  并在 stderr 打一行 `[汇总] 成功 n，失败 m`。
- **依赖分家**：`requirements.txt` 去掉 `pytest`，新增 `requirements-dev.txt`
  （`-r requirements.txt` + `pytest`）。运行时不再需要测试框架。
- **删死代码**：`Signal.second_order`（从未赋值，`Peak.apparent` 因此恒为 False）、
  `Peak.apparent`、`magnetic_equivalence_ok`（定义后从未调用）。
  `hose_codes_with_hydrogens` **保留**——它被 `tests/test_hose.py` 用来钉住
  "氢作为节点参与展开"以及记录"非对映异位氢同码"这条已知局限，不是死代码。
- README 流程图里「Tier B 指纹检索（可选，未启用）」已过期，改为「默认开启」；
  删掉 `--style cn-pdf` 示例。
- 测试：**229 passed**（删掉 1 条 cn-pdf 用例，其余不变）。

## 0.2.5 (2026-10-28) — 回灌 0.2.2–0.2.4，并修 AX/AB 丢线

- 把 0.2.2（5 个 P0）、0.2.3（P1）、0.2.4（P2 可做部分）的改动合并回**源码仓库**。
  此前这些改动只存在于便携包与补丁里，从仓库打包会退回 0.2.1。
- **修 `exact_transitions` 的耦合预算**：0.2.4 用全局 `Σ|J|/2` 判断“谱线离所属自旋
  太远就丢”，而 AX / AB 体系的四条线正好落在 ±J/2 附近，二阶修正后略超即被丢弃——
  AX、AB 都只剩 2 条线，强度又被重新归一化，总积分仍然正确，所以不报警。
  改为**每个自旋自己的** `Σ|J|`。
- **强度归一化移到所有过滤之后**：此前被过滤掉的线会带走强度（胆固醇 ¹H 46.0 → 45.99）。
- 仓库测试扩到 **222 条**，其中包内那 17 条（`tests/test_regressions.py`）现在也由 pytest 执行，
  避免再出现“包里测不出、仓库测得到”的分离；新增 AB 体系 4 条线的回归。

## 0.2.4 (2026-10-27) — P2 交付/交互

- GUI 单分子模式增加 **2D 结构预览**（RDKit），便于发现 SMILES 输错。
- GUI 增加 **批量模式**：每行一个 SMILES，一次最多 30 个；结果写到
  `results/batch_时间戳/`，提供汇总表与整包 ZIP。
- GUI 增加 **下载全部结果（ZIP）**（单分子也支持）。
- 新增 `停止 NMR.bat`：结束监听 8901 的 nmrflow 服务，替代“去任务管理器找进程”。
- CLI 新增 `--version`。
- 新增 `CHANGELOG.md`、`SHA256SUMS.txt`。
- `NMR-Predictor.exe` 仍未签名（无代码签名证书）；README/说明里明确提示核对 SHA256
  或改用 `启动 NMR.bat`。exe 的版本资源仍未补（需要重新冻结 exe）。

## 0.2.3 (2026-10-27) — P1 科学/置信度/可复现边界

- 置信度改为“球层 + n=1 折扣”；与 Tier A 分歧只提示、不降 HOSE 置信度。
- 十氢萘/降冰片烷等小双环脂肪体系不再采信 HOSE，改走经验规则 + 桥头/桥碳修正；
  双环[2.2.2]辛烷等仍未标定。
- ³¹P 默认 `protect_curated`：磷酸三甲酯回到 0.00（不再被 2 层球 HOSE 带到 14.35）；
  ³¹P 整体 MAE 19.71 / P90 51 未变。
- 强耦合伪线增加逐跃迁归属距离过滤；合并峰增加跨度上限；胆固醇 ¹H 最大峰跨度 0.14 ppm。
- 包内新增 17 项回归（stdlib unittest）、`REPRODUCE.md`、索引 SHA256。
- 完整 202 条测试 / `tools/` 评估脚本 / 原始 SD 仍不随包。

## 0.2.2 (2026-10-26) — P0 修复

- 默认不生成 Found；`--found` 显式开启，导出物带 `[SIMULATED]`（0.2.9 已删除该参数）。
- ¹H 低置信度线宽按核限幅（2 Hz），乙醇 CH₂ 恢复四重峰。
- 强耦合不再输出负 ppm 镜像线；>10 自旋 warning 接回报告。
- 带电物种按 [M]⁺/[M]⁻ 报；Na⁺/Cl⁻ 不再出 HRMS。
- GUI 把 PNG/TXT 落到 `results/` 并提供下载。

## 0.2.1 (2026-10-26)

- HOSE 索引、四核预测、FluoBase 对照、便携包 0.2.1。
