# piRNA PPI

[English](README.md) | 中文

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

本仓库包含一套 piRNA 相关蛋白-蛋白相互作用（PPI）筛选流程。该流程结合蛋白结构域切分、paired MSA 构建、RoseTTAFold2-PPI 打分以及后处理工具，用于排序候选互作蛋白对。

项目整体结构参考 Cookiecutter Data Science，并使用 Pixi 管理 Python 和命令行环境。

## 当前状态

当前可用的工作流程围绕以下阶段组织：

1. 准备蛋白序列和结构预测衍生输入。
2. 解析 DPAM 和 InterProScan 结构域注释。
3. 使用 `pirna_ppi.protein_segmenter` 切分大蛋白。
4. 使用 `pirna_ppi.protein_pairs` 生成片段或蛋白组合。
5. 使用 `pirna_ppi.generate_paired_msa` 生成 paired MSA。
6. 使用 `pirna_ppi.fix_msa` 去除 A3M 中的小写插入状态。
7. 在外部运行 RoseTTAFold2-PPI，对 paired MSA index 进行预测。
8. 使用 `pirna_ppi.parse_rf2ppi_results` 和 `pirna_ppi.filter_msa_by_probability` 汇总并筛选预测结果。

仓库中也保留了探索性代码或脚手架代码。`dataset.py`、`code_scaffold.py`、`plots.py` 和 `modeling/` 不属于当前主流程。`run_dca.py` 是实验性旁路。`DPAM/` 和 `RoseTTAFold2-PPI/` 目录如果在本地出现，属于外部工具副本，不应视为本项目的一手源码。

## 快速开始

创建项目环境：

```bash
make requirements
# 或
pixi install
```

激活环境：

```bash
pixi shell
```

常用代码质量检查：

```bash
make lint
make format
```

构建文档：

```bash
pixi run mkdocs build -f docs/mkdocs.yml
```

## 环境与依赖

环境定义在 `pixi.toml` 中，并由 `pixi.lock` 锁定。

- Python：`~=3.10`
- 开发工具：`ruff`、`tqdm`、`typer`、`loguru`
- 生物信息和命令行工具：`mmseqs2`、`hmmer`、`aria2`
- Python 包：`python-dotenv`、`mkdocs`、`biopython (<1.86)`、`polars`、`numpy`、`pyreadr`、`joblib`
- 外部 Python 依赖：来自 GitHub 的 `colabfold`

部分 GPU 相关依赖，尤其是 ColabFold 工作流中可能用到的 JAX/JAXLIB，可能需要根据本机 CUDA/cuDNN 环境手动安装。建议在 Pixi 环境中使用 `python -m pip ...`，并避免混用用户级包目录（`~/.local`）中的包。

示例：

```bash
pixi shell
python -m pip install \
  jaxlib==0.3.25+cuda11.cudnn82 \
  -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

## 数据目录

默认路径定义在 `pirna_ppi/config.py` 中。`data/` 目录不会提交到 GitHub，因为原始输入、paired MSA、中间日志和序列数据库通常体积很大，并且依赖具体机器环境。克隆仓库后，需要另外准备或获取这些文件，才能运行完整流程。

- `data/raw`：原始输入，例如蛋白 FASTA 文件。
- `data/interim`：中间产物，包括蛋白切分输出、paired MSA、index、日志和外部工具工作文件。
- `data/processed`：下游汇总表和最终分析结果。
- `data/external`：第三方数据库，例如 ColabFold/MMseqs2 数据库。
- `models`：后续若产生项目级模型文件，可放在此处。
- `reports/figures`：分析报告中使用的生成图件。

大型数据和外部软件目录已被 Git 忽略。服务器或工作站上的本地 checkout 可能保留这些文件，便于复现和检查；但 GitHub clone 预期只包含代码、配置和文档。

## 主要模块

- `pirna_ppi/config.py`：集中定义路径。
- `pirna_ppi/utils.py`：分析 notebook 中使用的小型辅助函数。
- `pirna_ppi/protein_segmenter/`：基于结构域信息的蛋白切分包。
- `pirna_ppi/protein_pairs.py`：构建有效的片段或蛋白组合，并输出 pair FASTA 和汇总表。
- `pirna_ppi/generate_paired_msa.py`：以 MMseqs2/ColabFold 风格生成 RoseTTAFold2-PPI 所需的 paired MSA。
- `pirna_ppi/fix_msa.py`：移除 A3M 小写插入状态，并重写 MSA index 路径。
- `pirna_ppi/extract_known_ppi_index.py`：从 MSA index 中筛选已知 PPI 对应条目。
- `pirna_ppi/parse_rf2ppi_results.py`：将 RF2-PPI 日志转为片段级和蛋白对级汇总表。
- `pirna_ppi/filter_msa_by_probability.py`：将高概率互作对应的 MSA 文件复制到筛选后的输出目录。
- `pirna_ppi/repeat_rf2ppi.py`：重复运行 RF2-PPI，用于预测波动性分析。

## 外部工具

当前流程依赖若干外部系统：

- ColabFold 或 AlphaFold 输出提供结构预测文件和 PAE JSON，用于 DPAM 和蛋白切分前处理。
- DPAM 提供来自预测结构的结构域注释。
- InterProScan 提供备用结构域注释。
- ColabFold/MMseqs2 数据库用于 paired MSA 生成。
- RoseTTAFold2-PPI 通常通过 Singularity 在外部运行，输入为本仓库生成的 paired MSA index。

`DPAM/` 和 `RoseTTAFold2-PPI/` 目录如果存在，是上游工具的本地副本。它们可能不会出现在 GitHub clone 中，应按上游项目说明单独安装或 checkout。其上游 README 不作为本项目主文档。

## 预期数据与本地输出

本仓库预期本地存在一个 `data/` 目录，用于放置原始输入、中间文件和生成结果。该目录不会上传到 GitHub。以作者的一次本地运行为例，`data/interim/piRNA_PPI_20260210/` 中曾包含：

- 蛋白切分输出：来自 276 个蛋白的 687 条片段/蛋白记录。
- 蛋白组合输出：长度筛选后得到 220673 个有效组合。
- paired MSA 生成汇总：生成 220673 个 MSA 文件。

这些数字只是一次本地分析状态的示例；新读者克隆仓库后不应预期马上看到这些文件。如果输入 FASTA、切分阈值、数据库或筛选阈值发生变化，结果也会随之改变。

## 文档

项目文档存放在 `docs/docs/` 中，并使用 MkDocs 构建：

```bash
pixi run mkdocs build -f docs/mkdocs.yml
pixi run mkdocs serve -f docs/mkdocs.yml
```

构建文档不需要 `data/` 或外部工具目录。

`docs/docs/getting-started.md` 是写给不太熟悉编程和命令行的新成员的中文学习指南。

## 许可协议

MIT License，详见 `LICENSE`。
