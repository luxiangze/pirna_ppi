# piRNA PPI

English | [中文](README.zh-CN.md)

<a target="_blank" href="https://cookiecutter-data-science.drivendata.org/">
    <img src="https://img.shields.io/badge/CCDS-Project%20template-328F97?logo=cookiecutter" />
</a>

This repository contains a piRNA-related protein-protein interaction (PPI)
screening workflow. It combines protein domain segmentation, paired MSA
construction, RoseTTAFold2-PPI scoring, and post-processing utilities for
ranking candidate interacting protein pairs.

The project follows a Cookiecutter Data Science-style layout and uses Pixi to
manage the Python and command-line environment.

## Status

The current working pipeline is organized around the following stages:

1. Prepare protein sequences and structure-derived inputs.
2. Parse DPAM and InterProScan domain annotations.
3. Split large proteins with `pirna_ppi.protein_segmenter`.
4. Generate segment or protein pairs with `pirna_ppi.protein_pairs`.
5. Generate paired MSAs with `pirna_ppi.generate_paired_msa`.
6. Remove A3M lowercase insertion states with `pirna_ppi.fix_msa`.
7. Run RoseTTAFold2-PPI externally on the paired MSA index.
8. Summarize and filter predictions with `pirna_ppi.parse_rf2ppi_results` and
   `pirna_ppi.filter_msa_by_probability`.

The repository also contains exploratory or scaffold code. `dataset.py`,
`code_scaffold.py`, `plots.py`, and `modeling/` are not part of the active
pipeline. `run_dca.py` is an experimental side path. The `DPAM/` and
`RoseTTAFold2-PPI/` directories, when present locally, are external tool
checkouts and should not be treated as first-party project code.

## Quickstart

Create the project environment:

```bash
make requirements
# or
pixi install
```

Activate the environment:

```bash
pixi shell
```

Common quality checks:

```bash
make lint
make format
```

Build the documentation:

```bash
pixi run mkdocs build -f docs/mkdocs.yml
```

## Environment and Dependencies

The environment is declared in `pixi.toml` and locked in `pixi.lock`.

- Python: `~=3.10`
- Development tools: `ruff`, `tqdm`, `typer`, `loguru`
- Bioinformatics and CLI tools: `mmseqs2`, `hmmer`, `aria2`
- Python packages: `python-dotenv`, `mkdocs`, `biopython (<1.86)`,
  `polars`, `numpy`, `pyreadr`, `joblib`
- External Python dependency: `colabfold` from GitHub

Some GPU-specific packages, especially JAX/JAXLIB used by ColabFold workflows,
may need to be installed manually to match the local CUDA/cuDNN stack. When
doing so, run `python -m pip ...` inside the Pixi environment and avoid mixing
packages from the user site (`~/.local`).

Example:

```bash
pixi shell
python -m pip install \
  jaxlib==0.3.25+cuda11.cudnn82 \
  -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

## Data Layout

Default paths are defined in `pirna_ppi/config.py`. The `data/` directory is
intentionally not committed to GitHub because raw inputs, paired MSAs,
intermediate logs, and sequence databases can be large and machine-specific.
After cloning the repository, prepare or obtain these files separately before
running the full workflow.

- `data/raw`: original inputs, such as protein FASTA files.
- `data/interim`: intermediate artifacts, including segmentation outputs,
  paired MSAs, indices, logs, and external-tool working files.
- `data/processed`: downstream summaries and finalized analysis results.
- `data/external`: third-party databases, such as ColabFold/MMseqs2 databases.
- `models`: project-level model artifacts, if any are produced later.
- `reports/figures`: generated figures for analysis reports.

Large data and external software directories are ignored by Git. A local server
or workstation checkout may contain these files for reproducibility and
inspection, but a GitHub clone is expected to contain only code, configuration,
and documentation.

## Main Modules

- `pirna_ppi/config.py`: central path definitions.
- `pirna_ppi/utils.py`: small helpers used by the analysis notebook.
- `pirna_ppi/protein_segmenter/`: domain-aware protein segmentation package.
- `pirna_ppi/protein_pairs.py`: builds all valid segment or protein pairs and
  writes pair FASTA plus summary tables.
- `pirna_ppi/generate_paired_msa.py`: MMseqs2/ColabFold-style paired MSA
  generation for RoseTTAFold2-PPI.
- `pirna_ppi/fix_msa.py`: removes lowercase A3M insertion states and rewrites
  MSA index paths.
- `pirna_ppi/extract_known_ppi_index.py`: filters MSA index entries for known
  PPI pairs.
- `pirna_ppi/parse_rf2ppi_results.py`: converts RF2-PPI log files into
  segment-level and protein-pair-level summary tables.
- `pirna_ppi/filter_msa_by_probability.py`: copies high-probability MSA files
  to a filtered output directory.
- `pirna_ppi/repeat_rf2ppi.py`: repeats RF2-PPI runs for variability analysis.

## External Tools

The active workflow depends on several external systems:

- ColabFold or AlphaFold outputs provide predicted structures and PAE JSON
  files used before DPAM and segmentation.
- DPAM provides domain annotations from predicted structures.
- InterProScan provides fallback domain annotations.
- ColabFold/MMseqs2 databases are required for paired MSA generation.
- RoseTTAFold2-PPI is run externally, usually through Singularity, using the
  paired MSA index generated by this repository.

The directories `DPAM/` and `RoseTTAFold2-PPI/` are local copies of upstream
tools when they are present. They may be absent from a GitHub clone and should
be installed or checked out separately according to the upstream projects. Their
upstream README files are not the main project documentation.

## Expected Data and Local Outputs

The repository expects a local `data/` tree with raw inputs, intermediate files,
and generated outputs. This tree is not uploaded to GitHub. In the author's
local run, `data/interim/piRNA_PPI_20260210/` included:

- Protein segmentation output: 687 segment/protein entries from 276 proteins.
- Pair generation output: 220673 valid pairs after the length filter.
- Paired MSA generation summary: 220673 generated MSA files.

These numbers are examples from one local analysis state; they are not files a
new reader should expect to see immediately after cloning the repository. They
may change if the input FASTA, segmentation threshold, databases, or filtering
thresholds change.

## Documentation

Project documentation is stored in `docs/docs/` and built with MkDocs:

```bash
pixi run mkdocs build -f docs/mkdocs.yml
pixi run mkdocs serve -f docs/mkdocs.yml
```

The documentation build does not require `data/` or external tool directories.

`docs/docs/getting-started.md` contains a Chinese learning guide for new lab
members who are less familiar with programming and command-line workflows.

## License

MIT License. See `LICENSE`.
