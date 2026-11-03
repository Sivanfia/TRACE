# GTEx：修正 T0 后仅重跑 TRACE-EM

这个包保留 ACVB、TRADER、Trans-Lasso 和 CONCERT 的原结果，仅重新拟合 TRACE-EM。没有 Oracle TRACE：真实数据中没有已知的真实源权重。

## 默认运行

解压并进入这个文件夹，使用原服务器的 R 环境：

```bash
nohup bash run_48cores.sh > run.log 2>&1 &
```

默认 48 个并行任务、六种目标组织各 50 次划分，共 300 次新 TRACE-EM 拟合。每个任务的 BLAS/OpenMP 只用一个线程。

查看日志：

```bash
tail -f run.log
```

Ctrl+C 仅停止查看日志，后台计算继续。中断后运行相同启动命令会跳过已完成的新 TRACE 拟合；失败任务会保留记录并重新尝试。

如果缺少依赖：

```bash
Rscript scripts/install_dependencies.R
```

只需要 glmnet、quadprog，不需要重新编译或运行其它比较方法。安装脚本只安装缺少的包，不主动升级已有包。

## 保留与重算的范围

- 原结果：`TRACE_real_data_with_CONCERT_server/results_50rep`，每个比较方法 300 条，共 1200 条逐次指标；同时保留预测、基因选择及 CONCERT 的迁移概率。
- 原文表格对应这次运行的配置为 `rho_z=0.15`、`rho_c=0.3`。包内采用该运行清单，未沿用后来编辑的 `rho_c=0.1`。`config/default.R` 从 `retained_results/original_manifest.rds` 直接读取原配置，避免数值转成文本时的舍入差异；`config/original_values.txt` 提供可读副本。
- 旧 TRACE 的指标、预测、最终权重和最终选择不会参与合并。只保留与 T0 修正无关的目标初始化作为输入一致性核对。
- 六个脑组织、固定五个来源、原 80/20 供体划分、筛选后 250 个基因、训练标准化、源系数摘要和 bootstrap 精度校准全部复用。包内含所需的六个 GTEx 模块数据文件，无需联网下载数据。
- TRACE 的变量选择方式和原预测评估不变。选择模型后，用目标训练数据的供体级五折交叉验证选择 ridge 惩罚；同一划分的交叉验证种子保持为原来的 `seed+1`，不因现在只运行一种方法而改变。
- T0 按删除当前坐标后的模型计算，T1 的对应补偿项同步修改。源权重更新、先验、停止规则和其它方法的结果保持原分析方式。

## 输出

正式结果写到新目录 `results_T0_corrected/`，不会覆盖包内旧结果。

- `trace_refit/jobs/`：新 TRACE 的逐次拟合、权重和初始化核对。
- `comparison/jobs/`：每次新 TRACE 与同一次划分的四种原比较方法合并。
- `comparison/tables/gtex_mse_table.tex`、`.csv`：六组织 MSE 表，均值和标准差，数值最低者在 TeX 表中加粗。
- `comparison/tables/model_size_table.tex`、`.csv`：选择变量数量表。
- `comparison/tables/metrics_per_rep.csv`、`metrics_summary.csv`：逐次指标及汇总，含原方法预测诊断。
- `comparison/tables/paired_mse_differences.csv`：同划分比较。
- `comparison/tables/trace_alpha_per_rep.csv`、`gene_results_per_rep.csv`、`gene_stability.csv`：新权重及选择结果。
- `comparison/tables/result_provenance.csv`：每条指标是重算还是保留的说明。
- `comparison/paper_figures/real_data_amygdala_mse.pdf`、`.png`，`real_data_hippocampus_mse.pdf`、`.png`：原稿使用的两张图。
- `comparison/paper_figures/real_data_all_tissues_mse.pdf`、`.png`：六组织合图，另有每个组织的单图。
- `comparison/figures/`：所有方法的指标图、TRACE 权重和基因稳定性、CONCERT 原迁移概率等附加图。
- `job_status.csv`：任务失败原因。未成功的新 TRACE 不会与旧比较结果混合成一个有效划分。

只重画已有结果，不重新拟合：

```bash
Rscript plot_all.R
```

修改线程数或新输出位置：

```bash
WORKERS=24 REAL_OUTPUT_DIR=results_T0_corrected_24cores bash run_48cores.sh
```

更改配置、代码、数据或 R/依赖版本后，请使用新的结果目录。程序会拒绝混用不同配置的检查点。

## 本地小规模检查

```bash
bash validate_all.sh
```

运行矩阵公式检查、保存划分的重建核对、合并保护测试，以及六种组织各一次的真实数据拟合和重新画表画图。检查输出放在临时目录，正式 300 次结果未预先放入包中。

原计算环境是 R 4.4.1 / Linux。优先沿用原服务器环境。重算前会核对原始数据、供体、基因顺序、标准化、源摘要、精度和目标初始化，一致性不满足时记录错误，不合并对应划分。

`SHA256SUMS` 用于验证上传完整性：

```bash
sha256sum -c SHA256SUMS
```
