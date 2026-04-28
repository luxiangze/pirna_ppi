import os
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from loguru import logger
from tqdm import tqdm
import typer

from pirna_ppi.config import INTERIM_DATA_DIR

app = typer.Typer()

DEFAULT_INPUT_DIR = INTERIM_DATA_DIR / "paired_msa"
DEFAULT_OUTPUT_DIR = INTERIM_DATA_DIR / "paired_msa_fixed"

# Pattern to match lowercase amino acid insertion states in A3M format
LOWERCASE_PATTERN = re.compile(r"[acdefghiklmnpqrstvwy]+")


def _fix_single_file(args: tuple[Path, Path]) -> None:
    """Process a single A3M file."""
    input_file, output_file = args
    content = input_file.read_text()
    fixed_content = LOWERCASE_PATTERN.sub("", content)
    output_file.write_text(fixed_content)


@app.command()
def main(
    input_dir: Path = typer.Option(
        DEFAULT_INPUT_DIR,
        help="Input directory containing A3M files",
    ),
    output_dir: Path = typer.Option(
        DEFAULT_OUTPUT_DIR,
        help="Output directory for fixed A3M files",
    ),
    input_index: Path = typer.Option(
        None,
        help="Input index TSV file (path<TAB>length). If provided, generates a new index with updated paths.",
    ),
    output_index: Path = typer.Option(
        None,
        help="Output index TSV file. Defaults to <output_dir>_index.tsv if input_index is given.",
    ),
    workers: int = typer.Option(
        0,
        help="Number of worker processes (0 = auto, use all CPUs minus 2)",
    ),
    resume: bool = typer.Option(
        True,
        "--resume/--no-resume",
        help="Skip files that already exist in output directory",
    ),
):
    """Remove lowercase insertion states from A3M files for RoseTTAFold2-PPI compatibility."""
    output_dir.mkdir(parents=True, exist_ok=True)

    input_files = list(input_dir.glob("*.a3m"))
    logger.info("Found {} A3M files in {}", len(input_files), input_dir)

    # Filter out already processed files if resume is enabled
    if resume:
        existing = {f.name for f in output_dir.glob("*.a3m")}
        tasks = [(f, output_dir / f.name) for f in input_files if f.name not in existing]
        skipped = len(input_files) - len(tasks)
        if skipped > 0:
            logger.info("Resuming: skipping {} already processed files", skipped)
    else:
        tasks = [(f, output_dir / f.name) for f in input_files]

    if not tasks:
        logger.success("All files already processed!")
        return

    # Determine worker count
    if workers <= 0:
        workers = max(1, os.cpu_count() - 2)
    logger.info("Using {} worker processes for {} files", workers, len(tasks))

    with ProcessPoolExecutor(max_workers=workers) as executor:
        list(tqdm(executor.map(_fix_single_file, tasks), total=len(tasks), desc="Fixing A3M files"))

    logger.success("Fixed {} files, saved to {}", len(tasks), output_dir)

    # Generate new index file if input_index is provided
    if input_index is not None:
        if output_index is None:
            output_index = output_dir.parent / f"{output_dir.name}_index.tsv"
        _generate_index(input_index, input_dir, output_dir, output_index)


def _generate_index(input_index: Path, input_dir: Path, output_dir: Path, output_index: Path) -> None:
    """Generate a new index file with paths pointing to the output directory."""
    count = 0
    with input_index.open() as fin, output_index.open("w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            old_path = Path(parts[0])
            new_path = output_dir / old_path.name
            parts[0] = str(new_path)
            fout.write("\t".join(parts) + "\n")
            count += 1
    logger.success("Wrote index file: {} ({} entries)", output_index, count)


if __name__ == "__main__":
    app()
