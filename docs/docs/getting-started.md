# 上手指南

这份文档写给刚接触这个仓库、计算机基础还不太熟的新人。目标不是让人一天内完全看懂所有代码，而是先知道：这个项目在做什么、主线文件有哪些、哪些目录可以先不管，以及看到一个输出文件时该怎么判断它来自哪一步。

如果你是在 GitHub 上看的，很多数据和外部软件目录不会出现，这是正常的。`data/`、大型数据库、paired MSA、中间日志、`DPAM/` 和 `RoseTTAFold2-PPI/` 这类内容通常只在服务器或作者本地环境中存在。GitHub 上主要放代码、配置和文档。

## 先理解项目在做什么

这个仓库的目标是筛选 piRNA 相关蛋白之间可能存在的蛋白-蛋白相互作用（PPI）。简单说，它先把较长的蛋白按结构域切成更适合建模的小片段，再把片段两两组合，给每个组合生成 paired MSA，然后交给 RoseTTAFold2-PPI 预测互作概率，最后把高概率结果筛出来。

可以把主流程记成这一条线：

```text
蛋白序列和结构结果
  -> DPAM/InterProScan 结构域注释
  -> 蛋白分段
  -> 生成蛋白/片段组合
  -> 生成 paired MSA
  -> 修正 A3M
  -> RoseTTAFold2-PPI 预测
  -> 解析和筛选概率
```

## 第一次打开仓库先看什么

建议按这个顺序看，不要一上来从所有 Python 文件开始读。

1. 看顶层 `README.zh-CN.md`，先了解项目目标、环境、数据目录和主流程。
2. 看 `notebooks/piRNA_PPI.ipynb` 的 Markdown 标题，了解实际分析大概按哪些步骤做过。
3. 看 `pirna_ppi/config.py`，弄清楚 `data/raw`、`data/interim`、`data/processed` 这些路径是什么意思。
4. 看 `pirna_ppi/protein_segmenter/README_zh.md`，理解为什么要切分大蛋白。
5. 再看主线脚本：`protein_pairs.py`、`generate_paired_msa.py`、`fix_msa.py`、`parse_rf2ppi_results.py`、`filter_msa_by_probability.py`。

前几步只需要 GitHub 上的代码和文档就能读懂。真正运行流程时，才需要服务器上的 `data/`、外部数据库和工具容器。

## 环境和常用命令

这个项目用 Pixi 管环境。第一次使用时先安装依赖：

```bash
pixi install
```

进入环境：

```bash
pixi shell
```

如果只是想检查代码格式：

```bash
make lint
```

如果想看某个命令行脚本有哪些参数，可以加 `--help`。例如：

```bash
pixi run python -m pirna_ppi.protein_segmenter.cli --help
pixi run python -m pirna_ppi.protein_pairs --help
pixi run python -m pirna_ppi.generate_paired_msa --help
```

## 主线脚本怎么读

`pirna_ppi/protein_segmenter/`

这是蛋白切分模块。它读入蛋白序列、DPAM 结构域注释、InterProScan 注释和结构文件，然后输出四类 FASTA：`segmented.fasta`、`unsegmented.fasta`、`oversized_segments.fasta`、`no_domain.fasta`。

`pirna_ppi/protein_pairs.py`

这一步把上面得到的片段或蛋白记录两两组合。它会跳过来自同一个原始蛋白的片段组合，并按长度阈值筛掉太长的记录。

`pirna_ppi/generate_paired_msa.py`

这一步最重，依赖 MMseqs2 和 ColabFold 数据库。它根据 pair FASTA 生成 RoseTTAFold2-PPI 需要的 paired MSA 文件和 index。

`pirna_ppi/fix_msa.py`

RoseTTAFold2-PPI 对 A3M 格式比较敏感。这一步会去掉 A3M 里的小写插入状态，并可同步生成新的 index 文件。

`pirna_ppi/parse_rf2ppi_results.py`

RoseTTAFold2-PPI 跑完会产生日志。这一步把日志整理成片段级结果和蛋白对级结果。

`pirna_ppi/filter_msa_by_probability.py`

这一步根据互作概率阈值复制高概率结果对应的 MSA 文件，方便后续用 ColabFold 或其他方法继续建模。

## 本地结果怎么看

如果你拿到了服务器或作者本地的中间结果，它们通常在 `data/interim/piRNA_PPI_20260210/`。GitHub clone 里默认没有这些文件。几个最有用的 summary 文件是：

- `proteins_segment_output/summary.tsv`：蛋白切分结果统计。作者某次本地运行中有 687 条片段/蛋白记录，来自 276 个蛋白。
- `proteins_pair_output/summary.tsv`：蛋白或片段组合统计。作者某次本地运行中长度筛选后得到 220673 个有效组合。
- `pairs_msa_generation_output/summary.tsv`：paired MSA 生成统计。作者某次本地运行中生成了 220673 个 MSA 文件。
- `pairs_msa_generation_output/index_fixed.tsv`：RoseTTAFold2-PPI 使用的输入 index，每行通常是 `MSA路径<TAB>第一条蛋白长度`。

如果你不知道一个文件来自哪一步，先看它所在目录名。比如 `proteins_segment_output` 多半来自蛋白切分，`proteins_pair_output` 多半来自组合生成，`pairs_msa_generation_output` 多半来自 MSA 生成和 RF2-PPI 前后的处理。

## 哪些可以先不看

这些内容一开始可以跳过：

- `dataset.py`、`code_scaffold.py`、`plots.py`、`modeling/`：这些是模板或还没有接入当前主流程的文件。
- `run_dca.py`：这是实验性旁路，不是当前最主要的 RF2-PPI 流程。
- `DPAM/` 和 `RoseTTAFold2-PPI/`：这是外部工具代码。GitHub clone 里可能没有；即使服务器上有，也不需要一开始就读它们的源码。
- `data/interim/DPAM_input/step*` 和大量 `.log`：这些是运行产物，GitHub clone 里默认没有，主要用于服务器环境中的追踪和排错。
- `references/zhang_sm.txt`：这是论文补充材料文本，适合在理解方法背景时查阅，不是运行代码的入口。

## 建议的学习节奏

第一天先把路径和文件类型搞清楚：FASTA 是序列，A3M 是多序列比对，TSV 是表格，PDB 是结构，JSON 里常放 PAE 矩阵。

第二步读蛋白切分模块，因为它决定后面到底在预测“整条蛋白对”还是“片段对”。

第三步读 `protein_pairs.py`。如果你拿到了本地 summary，再结合 summary 确认为什么会得到这么多组合；如果只是在 GitHub 上看，先理解代码逻辑即可。

最后再看 paired MSA 和 RoseTTAFold2-PPI。它们依赖数据库、GPU、容器和大文件，先理解输入输出比直接重跑更重要。

## 小心事项

在服务器或作者本地环境中，不要随便删除 `data/interim` 下的结果；很多步骤计算量大，重新跑可能很久。GitHub 上默认没有这个目录。

不要直接修改 `DPAM/` 和 `RoseTTAFold2-PPI/` 的代码，除非你明确知道是在调外部工具。

如果某个命令看起来会跑很久，先用 `--help` 看参数，再检查输入文件是否存在。大型步骤最好在 `tmux` 或服务器任务系统里运行。
