# Protein Segmenter

`pirna_ppi.protein_segmenter` 用于将较长蛋白按结构域信息切分为更适合下游 PPI 筛选的小片段。它主要用于全长蛋白过长、不利于 paired MSA 生成或 RF2-PPI 预测的情况。

## 工作原理

1. 从 FASTA 文件读取蛋白序列。
2. 读取 DPAM 和 InterProScan 结构域注释。
3. 按 DPAM、InterProScan、无注释的顺序为每个蛋白分配结构域来源。
4. 如果某个长蛋白只分到一个结构域，但 InterProScan 提供了多个结构域，则回退到 InterProScan，以支持切分。
5. 判断每个蛋白是否需要切分以及如何切分：
   - `length <= threshold`：保留一个按结构域范围裁剪的片段，写入 `unsegmented.fasta`。
   - 无结构域注释：保留完整序列，写入 `no_domain.fasta`。
   - `length > threshold` 且有结构：使用 C-alpha 接触和 PAE 矩阵判断相邻结构域是否应合并。
   - `length > threshold` 且无结构：将结构域作为独立片段候选。
6. 输出片段 FASTA 文件，以及 summary 和 detail 表格。

结构信息参与时，结构域合并判定公式为：

```text
N_contact / L_min + N_contact / 100 > Avg_PAE / 8
```

## 模块结构

```text
protein_segmenter/
├── cli.py                     # CLI 入口和 run_segmentation()
├── config.py                  # SegmenterConfig 配置 dataclass
├── models.py                  # Domain, Protein, Segment, ElementPairMetrics
├── test_segmenter.py          # 轻量测试
├── core/
│   ├── segmenter.py           # 核心切分算法
│   └── contact_calculator.py  # 结构域间接触计算
└── parsers/
    ├── domain_parser.py       # DPAM 结构域解析
    ├── fasta_parser.py        # FASTA 解析
    ├── iprscan_parser.py      # InterProScan TSV 解析
    ├── pae_parser.py          # AlphaFold/ColabFold PAE JSON 解析
    └── pdb_parser.py          # PDB C-alpha 坐标解析
```

## 输入

下面的路径是本地分析环境中的预期路径。由于 `data/` 不进入版本控制，刚从 GitHub clone 下来的仓库里不一定存在这些文件。运行切分程序前，需要先准备好这些输入。

| 输入 | 默认路径 | 说明 |
| --- | --- | --- |
| 蛋白 FASTA | `data/raw/piRNA_PPI_protein.fasta` | 待切分的蛋白序列 |
| DPAM domains | `data/interim/DPAM_input/piRNA_domains` | 基于结构的结构域注释 |
| InterProScan TSV | `data/interim/iprscan5-output.tsv` | 备用的序列/结构域注释 |
| 结构目录 | `data/interim/DPAM_input/piRNA` | 按蛋白 ID 命名的 PDB 文件和 PAE JSON 文件 |

对于 ID 为 `PROTEIN_ID` 的蛋白，结构辅助切分会在结构目录下查找 `PROTEIN_ID.pdb` 和 `PROTEIN_ID.json`。

## 使用方式

### 命令行

先准备输入 FASTA、结构域注释和可选结构文件，然后运行：

```bash
pixi run python -m pirna_ppi.protein_segmenter.cli \
  --fasta-path data/raw/piRNA_PPI_protein.fasta \
  --dpam-domains-path data/interim/DPAM_input/piRNA_domains \
  --iprscan-path data/interim/iprscan5-output.tsv \
  --structure-dir data/interim/DPAM_input/piRNA \
  --output-dir data/interim/piRNA_PPI_20260210/proteins_segment_output \
  --max-segment-length 600
```

可用 `--help` 查看当前 CLI 参数：

```bash
pixi run python -m pirna_ppi.protein_segmenter.cli --help
```

### Python API

Python API 使用同样的本地输入文件：

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

## 输出

输出目录包含：

| 文件 | 说明 |
| --- | --- |
| `segmented.fasta` | 多段切分后符合长度要求的片段 |
| `unsegmented.fasta` | 未切分、或只按结构域范围裁剪后保留为单段的记录 |
| `oversized_segments.fasta` | 切分后仍超过长度阈值的片段 |
| `no_domain.fasta` | 没有可用结构域注释的蛋白 |
| `segment_detail.tsv` | 每个蛋白的切分决策和片段范围 |
| `summary.tsv` | 每类 FASTA 输出的数量和长度统计 |
| `segmenter.log` | 详细运行日志 |

FASTA header 格式如下：

```text
>PROTEIN_ID-seg1 range=1-450 domains=nD1,nD2
MSILKAYR...
```

## 配置

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `max_segment_length` | `750` | 判断是否需要切分的长度阈值 |
| `min_sequence_separation` | `10` | 接触计数所需的最小序列间隔 |
| `max_contact_distance` | `8.0` | 计算残基接触时允许的最大 C-alpha 距离 |
| `iprscan_priority` | `["Pfam", "Gene3D"]` | 优先使用的 InterProScan 注释来源 |

当前 piRNA 运行中，为了与后续组合生成保持一致，使用过 600 aa 阈值。

## 测试

```bash
pixi run python -m pirna_ppi.protein_segmenter.test_segmenter
```
