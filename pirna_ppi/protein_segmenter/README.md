# Protein Segmenter

`pirna_ppi.protein_segmenter` splits large proteins into domain-aware segments
for downstream PPI screening. It is designed for cases where full-length
proteins are too long for efficient paired-MSA generation or RF2-PPI
prediction.

## How It Works

1. Load protein sequences from a FASTA file.
2. Load domain annotations from DPAM and InterProScan.
3. Assign each protein a domain source in this order: DPAM, InterProScan, none.
4. If a long protein has only one assigned domain but InterProScan provides
   multiple domains, fall back to InterProScan to enable segmentation.
5. Decide whether and how to split each protein:
   - `length <= threshold`: keep one trimmed segment and write it to
     `unsegmented.fasta`.
   - no domain annotations: keep the full sequence and write it to
     `no_domain.fasta`.
   - `length > threshold` with structure: use C-alpha contacts and the PAE
     matrix to decide whether adjacent domains should be merged.
   - `length > threshold` without structure: treat domains as independent
     segment candidates.
6. Write segment FASTA files plus summary and detail tables.

The structure-aware merge criterion is:

```text
N_contact / L_min + N_contact / 100 > Avg_PAE / 8
```

## Module Structure

```text
protein_segmenter/
├── cli.py                     # CLI entry point and run_segmentation()
├── config.py                  # SegmenterConfig dataclass
├── models.py                  # Domain, Protein, Segment, ElementPairMetrics
├── test_segmenter.py          # Lightweight tests
├── core/
│   ├── segmenter.py           # Core segmentation algorithm
│   └── contact_calculator.py  # Inter-domain contact calculation
└── parsers/
    ├── domain_parser.py       # DPAM domain parser
    ├── fasta_parser.py        # FASTA parser
    ├── iprscan_parser.py      # InterProScan TSV parser
    ├── pae_parser.py          # AlphaFold/ColabFold PAE JSON parser
    └── pdb_parser.py          # PDB C-alpha coordinate parser
```

## Inputs

The paths below are expected local analysis paths. They are not guaranteed to
exist in a fresh GitHub clone because `data/` is intentionally kept outside
version control. Prepare these files before running the segmenter.

| Input | Default path | Description |
| --- | --- | --- |
| Protein FASTA | `data/raw/piRNA_PPI_protein.fasta` | Protein sequences to segment |
| DPAM domains | `data/interim/DPAM_input/piRNA_domains` | Structure-based domain annotations |
| InterProScan TSV | `data/interim/iprscan5-output.tsv` | Fallback sequence/domain annotations |
| Structure directory | `data/interim/DPAM_input/piRNA` | PDB files and PAE JSON files named by protein ID |

For a protein with ID `PROTEIN_ID`, the structure-aware path expects
`PROTEIN_ID.pdb` and `PROTEIN_ID.json` under the structure directory.

## Usage

### Command Line

Prepare the input FASTA, domain annotations, and optional structure files first,
then run:

```bash
pixi run python -m pirna_ppi.protein_segmenter.cli \
  --fasta-path data/raw/piRNA_PPI_protein.fasta \
  --dpam-domains-path data/interim/DPAM_input/piRNA_domains \
  --iprscan-path data/interim/iprscan5-output.tsv \
  --structure-dir data/interim/DPAM_input/piRNA \
  --output-dir data/interim/piRNA_PPI_20260210/proteins_segment_output \
  --max-segment-length 600
```

Use `--help` to inspect the current CLI options:

```bash
pixi run python -m pirna_ppi.protein_segmenter.cli --help
```

### Python API

The Python API uses the same expected local files:

```python
from pathlib import Path

from pirna_ppi.protein_segmenter.config import SegmenterConfig
from pirna_ppi.protein_segmenter.cli import run_segmentation

config = SegmenterConfig(
    fasta_path=Path("data/raw/piRNA_PPI_protein.fasta"),
    dpam_domains_path=Path("data/interim/DPAM_input/piRNA_domains"),
    iprscan_path=Path("data/interim/iprscan5-output.tsv"),
    structure_dir=Path("data/interim/DPAM_input/piRNA"),
    output_dir=Path("data/interim/piRNA_PPI_20260210/proteins_segment_output"),
    max_segment_length=600,
)
run_segmentation(config)
```

## Outputs

The output directory contains:

| File | Description |
| --- | --- |
| `segmented.fasta` | Valid segments from proteins split into multiple parts |
| `unsegmented.fasta` | Proteins or domain-trimmed records kept as one segment |
| `oversized_segments.fasta` | Segments still longer than the threshold |
| `no_domain.fasta` | Proteins without usable domain annotations |
| `segment_detail.tsv` | Per-protein segmentation decisions and segment ranges |
| `summary.tsv` | Counts and length statistics for each FASTA output |
| `segmenter.log` | Detailed processing log |

FASTA headers use this format:

```text
>PROTEIN_ID-seg1 range=1-450 domains=nD1,nD2
MSILKAYR...
```

## Configuration

| Parameter | Default | Description |
| --- | --- | --- |
| `max_segment_length` | `750` | Length threshold used to decide whether segmentation is needed |
| `min_sequence_separation` | `10` | Minimum sequence separation for contact counting |
| `max_contact_distance` | `8.0` | Maximum C-alpha distance for residue contact counting |
| `iprscan_priority` | `["Pfam", "Gene3D"]` | Preferred InterProScan annotation sources |

The current piRNA run used a 600 aa threshold for pair generation consistency.

## Testing

```bash
pixi run python -m pirna_ppi.protein_segmenter.test_segmenter
```
