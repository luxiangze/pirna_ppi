# Documentation

This directory contains the MkDocs documentation for the piRNA PPI project.
It is intended to build from a normal GitHub clone, even when large data files
and external tool checkouts are not present.

## Source Files

- `docs/mkdocs.yml`: MkDocs configuration.
- `docs/docs/index.md`: documentation landing page.
- `docs/docs/getting-started.md`: Chinese onboarding guide for new lab members.

Top-level project summaries live outside this directory:

- `README.md`: English peer-facing project overview.
- `README.zh-CN.md`: Chinese version with matching structure and content.

## Build Locally

From the repository root:

```bash
pixi run mkdocs build -f docs/mkdocs.yml
```

The generated site is written to `docs/site/`, which is ignored by Git.
The build does not require `data/`, sequence databases, DPAM outputs, or
RoseTTAFold2-PPI outputs.

## Preview Locally

From the repository root:

```bash
pixi run mkdocs serve -f docs/mkdocs.yml
```

MkDocs will print a local URL for previewing the documentation in a browser.
