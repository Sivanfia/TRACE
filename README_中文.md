# TRACE 源码版使用说明

本目录适合放在 GitHub 的 `main` 分支。算法、模拟、真实数据分析及绘图源码均保留，完整数据、逐次结果和缓存放在 Release 的 `simulation_and_real_data_analysis.zip` 中。

## 恢复数据并重画图表

1. 下载 **Releases 中手动上传的完整代码包**，不是 GitHub 自动生成的 “Source code” 压缩包。
2. 将完整包解压到另一个位置。
3. 在本源码目录中运行：

```bash
Rscript scripts/restore_release_data.R /完整包解压后/simulation_and_real_data_analysis
bash run_analysis.sh figures
```

路径中有空格时，请加英文双引号。恢复脚本检查源码版本，只复制所需的外部输入和结果，不覆盖源码。如果发现内容不同的本地数据，会停止复制。

新图表位于 `outputs/reproduced/`。这一步不重新拟合模型，只需支持 PDF、PNG 的 R。恢复的数据与新输出已由 `.gitignore` 排除，不会自动进入下一次 Git 提交。

## 重新拟合

安装缺少的 R 包后，在服务器上运行：

```bash
Rscript scripts/install_dependencies.R
nohup bash run_48cores.sh all > run_all.log 2>&1 &
```

这会依次重跑四组模拟及真实数据的全部方法。只重跑 TRACE 并保留其他方法时：

```bash
CORES=48 bash run_analysis.sh trace
```

模拟中重跑 TRACE-EM 和 Oracle TRACE，真实数据中只重跑 TRACE-EM。后者没有已知的真实源权重，因此不设置 Oracle。保留其他方法的运行方式会检查输入及初始化是否与原实验一致，应使用原实验的 R 和包版本。

默认使用 50 次重复。正式运行结果写入 `runs/`，原始结果保留。完整参数、环境记录和检验命令见英文 `README.md`。

原创代码暂未指定开放许可，第三方代码和数据保留原有条款。
