# piRNA PPI Documentation

This documentation describes the project-specific workflow for piRNA-related
protein-protein interaction screening. The top-level `README.md` and
`README.zh-CN.md` are the peer-facing project summaries; this MkDocs site keeps
the operational notes and onboarding material closer to the code.

The docs are written for both GitHub readers and local/server users. A GitHub
clone is expected to contain code, configuration, and documentation, while large
data directories, external databases, and external tool checkouts are prepared
separately in the local analysis environment.

## Pages

- [Getting started](getting-started.md): Chinese onboarding guide for new lab
  members, including what to read first, which scripts form the active
  pipeline, and which files can be skipped at the beginning.

## Workflow Summary

The active workflow is:

1. Prepare protein FASTA and structure-derived inputs.
2. Use DPAM and InterProScan annotations for domain-aware protein segmentation.
3. Generate valid segment or protein pairs.
4. Build paired MSAs for RoseTTAFold2-PPI.
5. Run RoseTTAFold2-PPI externally.
6. Parse, summarize, and filter interaction probabilities.

When the external tool directories `DPAM/` and `RoseTTAFold2-PPI/` are present
locally, they should be treated as upstream tool checkouts. The
project-specific logic is mainly under `pirna_ppi/`.

The MkDocs site can be built without `data/` or external tool directories.
