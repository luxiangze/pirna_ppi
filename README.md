 # piRNA PPI

 English | [中文](README.zh-CN.md)

 <a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
     <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
 </a>

 This repository contains a **piRNA–protein interaction / protein–protein interaction (PPI)** project.
 The codebase follows a Cookiecutter Data Science-like structure, and uses **Pixi** to manage the development
 environment.

 ## Status

 **placeholders** and should be replaced with project-specific logic.

 The `protein_pairs.py` utility is implemented and can generate all-vs-all protein pair FASTA records.

 ## Quickstart

 ### 1) Create the environment

 This project uses Pixi (`pixi.toml`, `pixi.lock`).

 ```bash
 make requirements
 # or
 pixi install
 ```

 Activate the environment:

 ```bash
 pixi shell
 ```

## Environment & Dependencies

The environment is declared in `pixi.toml`.

- Python: `~=3.10`
- Tools: `ruff`, `tqdm`, `typer`, `loguru`
- Bio/CLI deps: `mmseqs2`, `hmmer`, `aria2`
- PyPI deps: `python-dotenv`, `mkdocs`, `biopython (<1.86)`, `colabfold` (Git dependency)

## Manual dependency setup (optional)

Some dependencies (notably **GPU-enabled JAX/JAXLIB**) are commonly installed manually with `pip` to match your
CUDA/cuDNN stack.

Recommended practices:

 - Use `python -m pip ...` (inside the pixi environment) instead of calling `pip` directly.
 - Avoid mixing user-site packages (`~/.local`) into the environment.

 Example (CUDA wheel, adjust versions to your machine):

 ```bash
 pixi shell
 python -m pip install \
   jaxlib==0.3.25+cuda11.cudnn82 \
   -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
 ```

 ## Code Quality

 ```bash
 make lint
 make format
 ```

 ## Data Layout

 Default data directories are defined in `pirna_ppi/config.py`:

 - `data/raw`: immutable raw inputs
 - `data/interim`: intermediate artifacts
 - `data/processed`: canonical model-ready datasets
 - `data/external`: third-party data

 ## Project Organization

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

 ## License

 MIT License. See `LICENSE`.

