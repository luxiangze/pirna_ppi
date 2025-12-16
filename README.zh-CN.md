# piRNA PPI

[English](README.md) | 中文

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

本仓库包含一个 **piRNA-蛋白相互作用 / 蛋白-蛋白相互作用（PPI）** 项目的代码。
项目整体结构参考 Cookiecutter Data Science，并使用 **Pixi** 管理开发环境。

## 当前状态

`protein_pairs.py` 已实现，可用于将蛋白 FASTA 生成两两组合的蛋白对 FASTA。

## 快速开始

### 1) 创建环境

本项目使用 Pixi（`pixi.toml`、`pixi.lock`）。

```bash
make requirements
# 或
pixi install
```

激活环境：

```bash
pixi shell
```

## 环境与依赖

环境定义在 `pixi.toml` 中。

- Python：`~=3.10`
- 工具：`ruff`、`tqdm`、`typer`、`loguru`
- 生信/命令行依赖：`mmseqs2`、`hmmer`、`aria2`
- PyPI 依赖：`python-dotenv`、`mkdocs`、`biopython (<1.86)`、`colabfold`（Git 依赖）

## 手动依赖配置（可选）

有些依赖（尤其是 **GPU 版本的 JAX/JAXLIB**）通常需要你根据机器上的 CUDA/cuDNN 情况手动用 `pip` 安装。

建议：

- 进入 pixi 环境后，用 `python -m pip ...` 替代直接执行 `pip`。
- 尽量避免把 `~/.local` 的用户级包混进环境里。

示例（CUDA 轮子，版本请按你的机器调整）：

```bash
pixi shell
python -m pip install \
  jaxlib==0.3.25+cuda11.cudnn82 \
  -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

## 代码质量

```bash
make lint
make format
```

## 数据目录

默认数据目录在 `pirna_ppi/config.py` 中定义：

- `data/raw`：原始输入（尽量保持不可变）
- `data/interim`：中间产物
- `data/processed`：模型可直接使用的最终数据
- `data/external`：第三方数据

## 项目结构

```text
├── LICENSE
├── Makefile
├── README.md
├── README.zh-CN.md
├── pixi.toml
├── pixi.lock
├── pyproject.toml
├── data/
├── docs/
├── models/
├── notebooks/
├── references/
├── reports/
└── pirna_ppi/
    ├── __init__.py
    ├── config.py
    ├── dataset.py
    ├── features.py
    ├── protein_pairs.py
    ├── utils.py
    ├── plots.py
    └── modeling/
        ├── train.py
        └── predict.py
```

## 许可协议

MIT License，详见 `LICENSE`。
