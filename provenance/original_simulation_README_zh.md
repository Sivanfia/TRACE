# 四组 simulation：修正 T0，保留对比方法结果

本包可独立上传服务器运行。默认只重跑 **TRACE-EM、Oracle TRACE**，然后自动合并原有 ACVB、Trans-Lasso、TRADER、CONCERT 结果，重新汇总并画图。原来的 TRACE 结果不参与合并。

## 开始运行

解压后进入本文件夹。继续使用生成旧结果时的服务器 R 环境，建议不要升级已经安装的 R 包。旧结果的 R 版本为 4.4.1。

若依赖尚未安装，先运行一次：

```bash
Rscript install_dependencies.R
```

默认每个参数组合运行 50 次重复，使用 48 个并行进程，每个进程的 BLAS/OpenMP 线程数为 1。四个 Setting 依次运行，不会同时启动四组 48 进程。

```bash
nohup bash run_48cores.sh > run.log 2>&1 &
```

查看进度：

```bash
tail -f run.log
```

仅核对配置而不开始拟合：

```bash
SIM_CONFIG_ONLY=TRUE bash run_all_50rep.sh
```

也可以单独运行一组，例如：

```bash
bash run_setting.sh 03
```

用较少进程：

```bash
SIM_NCORES=24 bash run_all_50rep.sh
```

## 实验设置

| Setting | 参数组合数 | 保留的设置 |
|---|---:|---|
| 1 | 48 | s=10,15,30；rho=0.3,0.6；SNR=3,4；污染比例=0,0.2,0.4,0.6；每个源500个样本 |
| 2 | 20 | TRADER I、III，各自 Ka=1,...,10；I 的源样本量为200，III 为500 |
| 3 | 48 | 与原 Setting 3 相同的网格；每个源500个样本 |
| 4 | 21 | s=10,15,30；rho=0.3；sigma=3；污染比例=0,0.1,...,0.6；每个源500个样本 |

合计 137 个参数组合，50 次重复对应 6,850 个数据重复、13,700 次 TRACE 拟合。先验、初始化及随机种子按原结果的运行清单保留。Setting 2 文件中后来改变过的 `z_prior` 和 `gamma_prior` 已恢复为运行清单实际使用的 `dimension_mild`、`constant_half`，避免混入先验变化。

本包已保留 26,900 条对比方法记录，所有已有方法/参数组合均有50次重复。Setting 2 的 CONCERT 只保留原有 Design I 的500条记录，Design III 不添加未运行过的 CONCERT 结果。

## 结果与图的位置

默认输出到本文件夹内的 `results_T0_corrected_trace/`。每组包含：

- `trace_refit/`：本次重跑的 TRACE 原始结果、迭代诊断、断点记录。
- `comparison/`：新 TRACE 与保留基线合并后的结果、均值、标准误及全部比较图。
- `comparison/result_provenance.csv`：新算与保留结果的记录数。

四组主文图还会集中复制到：

```text
results_T0_corrected_trace/paper_figures/
    setting01_main_f1.pdf / .png
    setting02_main_f1.pdf / .png
    setting03_main_f1.pdf / .png
    setting04_main_f1.pdf / .png
```

这些 PDF 的文件名与当前论文一致。Setting 1 使用 SNR=4 的两行三列图，Setting 3 使用 SNR=3 的两行三列图，Setting 2 为 I/III 两个面板，Setting 4 为一行三个面板。其余完整网格和指标也保留在每组 `comparison/figures/`。

需要重新画图而不重新拟合时运行：

```bash
Rscript plot_all.R
```

如使用自定义输出位置：

```bash
SIM_OUTPUT_BASE=/你的路径/TRACE_T0_new_results bash run_all_50rep.sh
Rscript plot_all.R /你的路径/TRACE_T0_new_results
```

## 保留结果与输入一致性

- `retained_results/` 保存原对比方法的逐重复结果，自动绘图不需要再提供旧文件。
- Setting 1、3、4 直接复用打包的 5,850 份源摘要，避免重复估计带来变化。Setting 2 原实验没有保存源摘要，仍按原代码和原种子估计。
- 本包保存了目标初始化和源权重初始化的数值对照。新旧输入或初始化不能对齐时，该重复会报错，不会静默与旧基线混合。请沿用原服务器的 R 和包环境，并保留错误日志用于排查。
- 只合并已经完成的新 TRACE 重复及其同一参数组合、重复编号、种子的基线记录，不会用旧 TRACE 补缺。
- 不覆盖任何旧实验目录。代码或参数改变后，需要新的输出目录，不能混合断点结果。
- 默认会继续本包同版本已经完成的重复。中断后再次运行同一命令即可续跑。保留基线仅有50次重复，因此本包的合并重跑支持 `SIM_NREP` 为1至50。

## 可选的全方法模式

默认不需要重算对比方法。如确实要运行各组原配置中的全部方法：

```bash
SIM_MODE=all bash run_all_50rep.sh
```

全方法模式另存到 `results_T0_corrected_all/`。Setting 1、3、4 会运行六种方法，Setting 2 按原 I/III 包运行五种方法，并继续保留已有的 CONCERT Design I 记录。相同方法已有的新结果会替换旧行，不会重复计数。

## 验证

```bash
bash validate_all.sh
```

验证使用独立的小规模临时数据，不改动正式实验参数或保留结果。其中 Setting 2 的输入路由检查使用 oracle 源摘要以缩短测试时间，正式实验仍使用原 spike-and-slab 摘要。

本包保留了原有第三方方法及其来源说明。校验文件为 `SHA256SUMS`。
