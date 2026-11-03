# TRACE 修订稿 Overleaf 项目

更新日期：2026-09-29。

## 最新正文修订：第 3.3 节的期望近似与交替更新

- 从最优源权重及其对应的最优变分因子出发，先展示精确的 log-odds 方程，再介绍相关期望与四个完整残差近似。
- 删除两个额外的更新函数记号，直接写出最优概率的近似关系和实际迭代公式。最优关系使用近似号，实际更新使用等号，左侧均保留 logit。
- 本节统一使用 expectations 描述近似，并明确先更新 inclusion、再更新 borrowing、最后更新源权重的顺序。
- 初始化依次使用 target-only inclusion、初始系数估计、Huber 源权重和第一轮 borrowing。
- 将系数估计与入选概率加权预测分开，算法框中的输出记号同步更新。
- 理论和实验的前向说明具体指向固定维度一致性与经验表现。

本次更新主文 TeX 和预览 PDF。已完整编译并检查第 3.3 节公式及算法框，无未解析引用或超宽公式提示。已有 Overleaf 项目只需替换 `draft_revised_crossrefs.tex`，图片无需更新。

## 此前附录修订：K_0 的引入位置

附录 B2 仅保留噪声和最小信号强度条件，K_0 改在 B4 首次引入。大模型引理中 K_0 的来源引用同步改为 B4。本次仅修改附录，替换 `supplement_information_revised_crossrefs.tex` 即可，图片无需更新。

## 理论记号介绍的文字修订

理论部分开头按“真实模型与索引方式、源信息与重要集合、算法与概率记号、范数与渐近记号”重新组织，并进一步压缩重复解释和过渡句。英文说明由 594 词减至 386 词，约缩短 35%（不计公式及标题）。全部记号定义、适用范围和五组独立公式均保留，数学含义及后续假设、定理和证明保持不变。该次记号修订仅更新主文文字。

## 本次 Simulation 修订

- 正文保留三个 setting，引言按“比较方法、实验设计、补充实验”组织。
- 各 setting 先解释源的构造逻辑、源与目标的关系及对应模型结构，再介绍公式和参数。
- 结果以文字分析为主，依次讨论 TRACE 的相对优势、迁移信息增益、污染变化、坐标借用筛选与源权重自适应更新，并保留图中显示的困难情况。
- 补充解释 Settings 1、3 在大支持集、强相关设计下的 recall 下降，以及 Setting 2 I 中额外误选造成的 precision 和 F1 损失；附录补充初始化和借用先验的诊断依据，未将可能机制写成已验证的单一原因。
- 正文共三幅组合图、八个 F1 面板：Setting 1 的 s=10/15/30；Setting 2 的 TRADER I/III；Setting 3 的 s=10/15/30。Settings 1、3 均使用 rho=0.6、SNR=3。
- 组合图统一颜色、线型、纵轴与共享图例；Setting 2 图标题不含 SNR。图中 264 个均值和 Monte Carlo 标准误均与原始结果逐项核对，没有重新运行或改变模拟数据。
- 原 Setting 4 移至附录 D.4，标题为 “Sensitivity to source-weight estimation”。该实验固定真实 alpha，通过 TRACE-EM 与 Oracle TRACE 的比较考察权重估计对模型选择的影响，并补充源权重误差和 F1 差异表。
- 附录保留 Settings 1、3 的完整 F1 网格；Setting 2 的两张 F1 图集中放入正文，附录增加 precision、recall 和 prediction MSE 的汇总表。
- 已修正 Setting 2 关于 CONCERT 是否参与比较的旧说明。理论和真实数据分析内容保持原样。

## 上传与编译

在 Overleaf 中选择 “New Project / Upload Project”，上传 `TRACE_revised_Overleaf.zip`。
编译器选择 **pdfLaTeX**。正文主文件为 `draft_revised_crossrefs.tex`；补充材料主文件为 `supplement_information_revised_crossrefs.tex`。两份文件独立编译，不需要额外的 `.bib` 文件。

也可以直接上传两份 `.tex` 和完整 `figures/` 文件夹。请保留相对目录和文件名。

本次正文已重新编译，共 26 页。附录沿用现有版本。正文交叉引用完整，第 3.3 节公式及算法框已检查排版。

### 已有 Overleaf 项目的增量上传

若已同步此前修订，本次只需替换 `draft_revised_crossrefs.tex`。若尚未同步此前的附录和 Simulation 修订，请同时上传当前附录及以下三幅组合图：

- `setting01_main_f1.pdf`：Setting 1 的 s=10、15、30。
- `setting02_main_f1.pdf`：Setting 2 的 I、III。
- `setting03_main_f1.pdf`：Setting 3 的 s=10、15、30。

当前两份文稿实际需要以下八个图文件；原有文件若已上传，无需重复上传：

| 用途 | `figures/` 下的文件 |
| --- | --- |
| 正文 Setting 1 | `setting01_main_f1.pdf` |
| 正文 Setting 2 | `setting02_main_f1.pdf` |
| 正文 Setting 3 | `setting03_main_f1.pdf` |
| 正文真实数据 | `real_data_amygdala_mse.pdf` |
| 正文真实数据 | `real_data_hippocampus_mse.pdf` |
| 附录 Setting 1 完整网格 | `setting01_all_f1.pdf` |
| 附录 Setting 3 完整网格 | `setting03_all_f1.pdf` |
| 附录 alpha 敏感性实验 | `setting04_all_f1.pdf` |

旧的 `setting02_I_f1.pdf`、`setting02_III_f1.pdf` 和各 `panel...` 文件目前未被引用，不是编译依赖。

## 图表与复现

正文新组合图为：

- `figures/setting01_main_f1.pdf`
- `figures/setting02_main_f1.pdf`
- `figures/setting03_main_f1.pdf`

`figures/simulation_f1_data.csv` 保存三幅组合图的 264 个 F1 均值、标准误、样本数及来源。
`figures/reproduce_simulation_figures.R` 可仅用 R 基础绘图重建这三幅图：

```sh
Rscript figures/reproduce_simulation_figures.R
```

`figure_sources.json` 列出原图来源和新组合图的重绘来源。原始单面板图和完整网格图继续保留，因此目录和上传包包含 15 个图 PDF；当前两份文稿实际引用其中 8 个文件。

simulation 来源根目录：

`/Users/sivanfia/Desktop/Research_process/CVB/draft/simulation_and_real_data_application/TRACE_four_settings_with_CONCERT_server_v1/results_50rep_with_CONCERT`

real_data 来源根目录：

`/Users/sivanfia/Desktop/Research_process/CVB/draft/simulation_and_real_data_application/TRACE_real_data_with_CONCERT_server/results_50rep`

项目编译与新组合图的复现仅依赖项目内文件，不依赖这些来源绝对路径。

## 文件同步

本文件夹中的两份 `.tex`、两份预览 PDF、`overleaf_project/` 和 `TRACE_revised_Overleaf.zip` 均与本次修订同步。上传压缩包包括两份 TeX、图文件、绘图数据与复现文件、本说明及来源清单，不包含编译中间文件。
