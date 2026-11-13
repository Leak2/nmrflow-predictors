# nmrflow 使用边界（ETHICS）

**一句话**：nmrflow 输出的是**预测数据**，不是实验数据。  
可以拿它做实验设计、谱图指认、教学和方法比较；**不要**把它包装成实测谱或论文里的 “Found”。

---

## 1. 预测不是实验

- `Calculated` 是模型算出来的；没有 `Found`，也没有模拟 `Found`。
- 预测误差在 README 的精度表里写得很清楚：常见分子 ¹H ~0.2–0.4 ppm、¹³C ~1–4 ppm；复杂/新骨架会更大。
- 预测谱不能替代 FID、原始谱图、元素分析或 HRMS 实测。

## 2. 明确禁止

- 把预测谱、预测 HRMS、预测位移写成 “实测/Found/experimental”。
- 去掉 `PREDICTED` 水印、`Description` metadata、`report.txt` 页脚、`manifest.json` 后当实验数据使用。
- 用本工具批量生成“论文补充材料”里的假谱图、假 Found、假纯度/收率证据。
- 在司法、临床、法规或安全相关场景中把预测当实测结论。

## 3. 鼓励的用途

- 实验前预判：该不该做、大概长什么样、哪里可能重叠。
- 谱图指认：候选结构对比、峰归属检查。
- 教学：讲化学位移、耦合、对称性。
- 方法研究：Tier A / HOSE / 自旋哈密顿量的误差分析。
- 流水线测试：在没有标准品时验证解析、绘图、导出格式。

## 4. 输出里已经带的标记

- `report.txt` 末尾：`Predicted by nmrflow 0.3.2 - not experimental data.`
- PNG 右下角：`PREDICTED - not experimental`
- PNG metadata：`Software=nmrflow 0.3.2`、`Description=Predicted NMR spectrum - not experimental data`
- JSON 顶层：`data_type: "predicted"`、`generator.version`、`generated_at`
- `manifest.json`：输入、选项、输出 SHA256、`data_type`

这些标记是给读者和审稿人看的，不是装饰。**不要移除它们。**

## 5. 许可证与可执行性

- 代码是 **BSD-3-Clause**：法律上允许使用、修改、再分发，包括商业用途。
- 本文件是**社区规范**，不是 BSD 附加条款，也不试图限制合法用途。
- 如果你需要“法律上可执行的禁止伪造”条款，那就不能同时声称是 OSI 开源；应改用非 OSI 的伦理/使用限制许可证，并明确说明。
- HOSE 索引来自 NMRShiftDB2，数据部分按 nmrshiftdb2 Database License 分发；代码与数据许可分开，详见 README §9。

## 6. 发现误用

如果你看到有人把 nmrflow 输出当实测数据，欢迎在 issue 里指出（附截图/链接）。  
项目不提供“帮你把预测做得更像实验”的支持；这类请求可以被拒绝。

---

*最后更新：2026-11-03*
