"""
Generate Paired MSA for RoseTTAFold2-PPI.

Logic adapted from: ColabFold (sokrypton/ColabFold)
Purpose: Concatenate sequences from the same UniRef cluster for complex prediction.
Compatible with ColabFold MSA server output format.
"""

import csv
import math
import os
import shutil
from pathlib import Path
import subprocess

from loguru import logger
from tqdm import tqdm
import typer

from pirna_ppi.config import EXTERNAL_DATA_DIR, INTERIM_DATA_DIR

# MMseqs2 module output positions for skip detection
MODULE_OUTPUT_POS = {
    "align": 4,
    "convertalis": 4,
    "expandaln": 5,
    "filterresult": 4,
    "lndb": 2,
    "mergedbs": 2,
    "mvdb": 2,
    "pairaln": 4,
    "result2msa": 4,
    "search": 3,
}

app = typer.Typer()


def run_mmseqs(mmseqs: Path, params: list[str | Path]) -> None:
    """Run MMseqs2 command with skip detection for existing outputs."""
    module = params[0]
    if module in MODULE_OUTPUT_POS:
        output_pos = MODULE_OUTPUT_POS[module]
        output_path = Path(params[output_pos]).with_suffix(".dbtype")
        if output_path.exists():
            logger.info(f"Skipping {module} because {output_path} already exists")
            return

    params_log = " ".join(str(i) for i in params)
    logger.info(f"Running {mmseqs} {params_log}")
    os.environ["MMSEQS_CALL_DEPTH"] = "1"
    subprocess.check_call([str(mmseqs)] + [str(p) for p in params])


def mmseqs_search_monomer(
    dbbase: Path,
    base: Path,
    mmseqs: Path = Path("mmseqs"),
    uniref_db: str = "uniref30_2302_db",
    metagenomic_db: str = "colabfold_envdb_202108_db",
    use_env: bool = True,
    filter_msa: bool = True,
    expand_eval: float = math.inf,
    align_eval: int = 10,
    diff: int = 3000,
    qsc: float = -20.0,
    max_accept: int = 1000000,
    prefilter_mode: int = 0,
    s: float = 8,
    db_load_mode: int = 2,
    threads: int = 32,
    gpu: int = 0,
) -> None:
    """Run MMseqs2 monomer search (uniref + optional envdb)."""
    if filter_msa:
        align_eval = 10
        qsc = 0.8
        max_accept = 100000

    # Check database index
    if (
        not dbbase.joinpath(f"{uniref_db}.idx").is_file()
        and not dbbase.joinpath(f"{uniref_db}.idx.index").is_file()
    ) or os.environ.get("MMSEQS_IGNORE_INDEX", False):
        logger.info("Search does not use index")
        db_load_mode = 0
        dbSuffix1 = "_seq"
        dbSuffix2 = "_aln"
    else:
        dbSuffix1 = ".idx"
        dbSuffix2 = ".idx"

    search_param = [
        "--num-iterations", "3",
        "--db-load-mode", str(db_load_mode),
        "-a", "-e", "0.1",
        "--max-seqs", "10000",
    ]
    if gpu:
        search_param += ["--gpu", str(gpu), "--prefilter-mode", "1"]
    else:
        search_param += ["--prefilter-mode", str(prefilter_mode)]
        if s is not None:
            search_param += ["-s", f"{s:.1f}"]

    filter_param = [
        "--filter-msa", str(int(filter_msa)),
        "--filter-min-enable", "1000",
        "--diff", str(diff),
        "--qid", "0.0,0.2,0.4,0.6,0.8,1.0",
        "--qsc", "0",
        "--max-seq-id", "0.95",
    ]
    expand_param = [
        "--expansion-mode", "0",
        "-e", str(expand_eval),
        "--expand-filter-clusters", str(int(filter_msa)),
        "--max-seq-id", "0.95",
    ]

    # UniRef search
    if not base.joinpath("uniref.a3m.dbtype").exists():
        # Clean up incomplete tmp directory if exists
        if (base / "tmp").exists():
            shutil.rmtree(base / "tmp")
        run_mmseqs(mmseqs, ["search", base / "qdb", dbbase / uniref_db, base / "res", base / "tmp", "--threads", str(threads)] + search_param)
        run_mmseqs(mmseqs, ["mvdb", base / "tmp/latest/profile_1", base / "prof_res"])
        run_mmseqs(mmseqs, ["lndb", base / "qdb_h", base / "prof_res_h"])
        run_mmseqs(mmseqs, ["expandaln", base / "qdb", dbbase / f"{uniref_db}{dbSuffix1}", base / "res", dbbase / f"{uniref_db}{dbSuffix2}", base / "res_exp", "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + expand_param)
        run_mmseqs(mmseqs, ["align", base / "prof_res", dbbase / f"{uniref_db}{dbSuffix1}", base / "res_exp", base / "res_exp_realign", "--db-load-mode", str(db_load_mode), "-e", str(align_eval), "--max-accept", str(max_accept), "--threads", str(threads), "--alt-ali", "10", "-a"])
        run_mmseqs(mmseqs, ["filterresult", base / "qdb", dbbase / f"{uniref_db}{dbSuffix1}", base / "res_exp_realign", base / "res_exp_realign_filter", "--db-load-mode", str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0", "--threads", str(threads), "--max-seq-id", "1.0", "--filter-min-enable", "100"])
        run_mmseqs(mmseqs, ["result2msa", base / "qdb", dbbase / f"{uniref_db}{dbSuffix1}", base / "res_exp_realign_filter", base / "uniref.a3m", "--msa-format-mode", "6", "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + filter_param)
        # Cleanup
        for db in ["res_exp_realign_filter", "res_exp_realign", "res_exp", "res"]:
            run_mmseqs(mmseqs, ["rmdb", base / db])
    else:
        logger.info("Skipping uniref search because uniref.a3m already exists")

    # Environmental DB search
    if use_env and not base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m.dbtype").exists():
        run_mmseqs(mmseqs, ["search", base / "prof_res", dbbase / metagenomic_db, base / "res_env", base / "tmp3", "--threads", str(threads)] + search_param)
        run_mmseqs(mmseqs, ["expandaln", base / "prof_res", dbbase / f"{metagenomic_db}{dbSuffix1}", base / "res_env", dbbase / f"{metagenomic_db}{dbSuffix2}", base / "res_env_exp", "-e", str(expand_eval), "--expansion-mode", "0", "--db-load-mode", str(db_load_mode), "--threads", str(threads)])
        run_mmseqs(mmseqs, ["align", base / "tmp3/latest/profile_1", dbbase / f"{metagenomic_db}{dbSuffix1}", base / "res_env_exp", base / "res_env_exp_realign", "--db-load-mode", str(db_load_mode), "-e", str(align_eval), "--max-accept", str(max_accept), "--threads", str(threads), "--alt-ali", "10", "-a"])
        run_mmseqs(mmseqs, ["filterresult", base / "qdb", dbbase / f"{metagenomic_db}{dbSuffix1}", base / "res_env_exp_realign", base / "res_env_exp_realign_filter", "--db-load-mode", str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0", "--max-seq-id", "1.0", "--threads", str(threads), "--filter-min-enable", "100"])
        run_mmseqs(mmseqs, ["result2msa", base / "qdb", dbbase / f"{metagenomic_db}{dbSuffix1}", base / "res_env_exp_realign_filter", base / "bfd.mgnify30.metaeuk30.smag30.a3m", "--msa-format-mode", "6", "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + filter_param)
        for db in ["res_env_exp_realign_filter", "res_env_exp_realign", "res_env_exp", "res_env"]:
            run_mmseqs(mmseqs, ["rmdb", base / db])
    elif use_env:
        logger.info("Skipping envdb search because bfd.mgnify30.metaeuk30.smag30.a3m already exists")

    # Merge results
    if use_env:
        run_mmseqs(mmseqs, ["mergedbs", base / "qdb", base / "final.a3m", base / "uniref.a3m", base / "bfd.mgnify30.metaeuk30.smag30.a3m"])
    else:
        run_mmseqs(mmseqs, ["mvdb", base / "uniref.a3m", base / "final.a3m"])


def mmseqs_search_pair(
    dbbase: Path,
    base: Path,
    mmseqs: Path = Path("mmseqs"),
    uniref_db: str = "uniref30_2302_db",
    prefilter_mode: int = 0,
    s: float = 8,
    threads: int = 64,
    gpu: int = 0,
    db_load_mode: int = 2,
    pairing_strategy: int = 0,
) -> None:
    """Run MMseqs2 pairaln for complex pairing."""
    if (
        not dbbase.joinpath(f"{uniref_db}.idx").is_file()
        and not dbbase.joinpath(f"{uniref_db}.idx.index").is_file()
    ) or os.environ.get("MMSEQS_IGNORE_INDEX", False):
        db_load_mode = 0
        dbSuffix1 = "_seq"
        dbSuffix2 = "_aln"
    else:
        dbSuffix1 = ".idx"
        dbSuffix2 = ".idx"

    search_param = [
        "--num-iterations", "3",
        "--db-load-mode", str(db_load_mode),
        "-a", "-e", "0.1",
        "--max-seqs", "10000",
    ]
    if gpu:
        search_param += ["--gpu", str(gpu), "--prefilter-mode", "1"]
    else:
        search_param += ["--prefilter-mode", str(prefilter_mode)]
        if s is not None:
            search_param += ["-s", f"{s:.1f}"]

    expand_param = ["--expansion-mode", "0", "-e", "inf", "--expand-filter-clusters", "0", "--max-seq-id", "0.95"]

    # Pairing search pipeline
    run_mmseqs(mmseqs, ["search", base / "qdb", dbbase / uniref_db, base / "res_pair", base / "tmp_pair", "--threads", str(threads)] + search_param)
    run_mmseqs(mmseqs, ["mvdb", base / "tmp_pair/latest/profile_1", base / "prof_res_pair"])
    run_mmseqs(mmseqs, ["lndb", base / "qdb_h", base / "prof_res_pair_h"])
    run_mmseqs(mmseqs, ["expandaln", base / "qdb", dbbase / f"{uniref_db}{dbSuffix1}", base / "res_pair", dbbase / f"{uniref_db}{dbSuffix2}", base / "res_pair_exp", "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + expand_param)
    run_mmseqs(mmseqs, ["align", base / "prof_res_pair", dbbase / f"{uniref_db}{dbSuffix1}", base / "res_pair_exp", base / "res_pair_exp_realign", "--db-load-mode", str(db_load_mode), "-e", "0.001", "--max-accept", "1000000", "--threads", str(threads)])
    run_mmseqs(mmseqs, ["pairaln", base / "qdb", dbbase / uniref_db, base / "res_pair_exp_realign", base / "res_pair_aln", "--db-load-mode", str(db_load_mode), "--pairing-mode", str(pairing_strategy), "--pairing-dummy-mode", "0", "--threads", str(threads)])
    run_mmseqs(mmseqs, ["align", base / "prof_res_pair", dbbase / f"{uniref_db}{dbSuffix1}", base / "res_pair_aln", base / "res_pair_aln_bt", "--db-load-mode", str(db_load_mode), "-e", "inf", "-a", "--threads", str(threads)])
    run_mmseqs(mmseqs, ["pairaln", base / "qdb", dbbase / uniref_db, base / "res_pair_aln_bt", base / "res_pair_final", "--db-load-mode", str(db_load_mode), "--pairing-mode", str(pairing_strategy), "--pairing-dummy-mode", "1", "--threads", str(threads)])
    run_mmseqs(mmseqs, ["result2msa", base / "qdb", dbbase / f"{uniref_db}{dbSuffix1}", base / "res_pair_final", base / "pair.a3m", "--db-load-mode", str(db_load_mode), "--msa-format-mode", "5", "--threads", str(threads)])

    # Cleanup
    for db in ["res_pair", "res_pair_exp", "res_pair_exp_realign", "res_pair_aln", "res_pair_aln_bt", "res_pair_final", "prof_res_pair", "prof_res_pair_h"]:
        if (base / db).with_suffix(".dbtype").exists():
            run_mmseqs(mmseqs, ["rmdb", base / db])
    if (base / "tmp_pair").exists():
        shutil.rmtree(base / "tmp_pair")


def parse_a3m(a3m_string: str) -> tuple[list[str], list[str]]:
    """Parse a3m format string and return sequences with headers.

    Compatible with ColabFold MSA server output.
    Preserves lowercase letters (insertion states in a3m format).

    Returns:
        A tuple of (sequences, headers).
    """
    sequences: list[str] = []
    headers: list[str] = []
    current_seq: list[str] = []

    for line in a3m_string.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            if current_seq:
                sequences.append("".join(current_seq))
            headers.append(line[1:])
            current_seq = []
        else:
            current_seq.append(line)

    if current_seq:
        sequences.append("".join(current_seq))

    return sequences, headers


def get_query_len(seq: str) -> int:
    """Get query length excluding lowercase insertion states."""
    return sum(1 for c in seq if c.isupper() or c == "-")


def load_pairs_and_sequences_from_fasta(
    pairs_file: Path,
) -> tuple[list[tuple[str, str]], dict[str, str]]:
    pairs: list[tuple[str, str]] = []
    protein_to_seq: dict[str, str] = {}

    current_pair: tuple[str, str] | None = None
    current_seq_lines: list[str] = []

    def _flush_current() -> None:
        nonlocal current_pair, current_seq_lines
        if current_pair is None:
            return
        seq_line = "".join(current_seq_lines).strip()
        if not seq_line:
            raise typer.BadParameter(f"Missing sequence for pair: {current_pair[0]}_{current_pair[1]}")
        if ":" not in seq_line:
            raise typer.BadParameter(
                f"Invalid sequence line for pair: {current_pair[0]}_{current_pair[1]}. Expected 'SEQ1:SEQ2'"
            )
        seq1, seq2 = seq_line.split(":", 1)
        seq1 = seq1.strip().upper()
        seq2 = seq2.strip().upper()
        if not seq1 or not seq2:
            raise typer.BadParameter(f"Empty sequence in pair: {current_pair[0]}_{current_pair[1]}")

        p1, p2 = current_pair
        pairs.append((p1, p2))

        if p1 in protein_to_seq and protein_to_seq[p1] != seq1:
            raise typer.BadParameter(f"Inconsistent sequence for protein {p1}")
        if p2 in protein_to_seq and protein_to_seq[p2] != seq2:
            raise typer.BadParameter(f"Inconsistent sequence for protein {p2}")
        protein_to_seq[p1] = seq1
        protein_to_seq[p2] = seq2

        current_pair = None
        current_seq_lines = []

    with pairs_file.open() as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _flush_current()
                header = line[1:]
                if "_" not in header:
                    raise typer.BadParameter(f"Invalid pair header: {header}")
                p1, p2 = header.split("_", 1)
                current_pair = (p1, p2)
            else:
                current_seq_lines.append(line)

    _flush_current()

    return pairs, protein_to_seq


# Output file names
LOG_FILE = "generate_paired_msa.log"
INDEX_TSV = "index.tsv"
SUMMARY_TSV = "summary.tsv"
MSA_SUBDIR = "paired_msa"


def _write_summary(
    output_dir: Path,
    n_pairs: int,
    n_proteins: int,
    n_unique_seqs: int,
    success_count: int,
    low_depth_count: int,
) -> None:
    """Write a TSV summary table."""
    path = output_dir / SUMMARY_TSV
    fieldnames = ["category", "count"]
    rows = [
        {"category": "total_pairs", "count": n_pairs},
        {"category": "unique_proteins", "count": n_proteins},
        {"category": "unique_sequences", "count": n_unique_seqs},
        {"category": "success_msa", "count": success_count},
        {"category": "low_depth_msa (<10 paired)", "count": low_depth_count},
        {"category": "failed_msa", "count": n_pairs - success_count},
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@app.command()
def local_pairaln(
    pairs_file: Path = typer.Option(
        INTERIM_DATA_DIR / "piRNA_PPI_20260210" / "protein_pairs_output" / "protein_pairs.fasta",
        help="FASTA file with pairs (>ID1_ID2) and sequences (SEQ1:SEQ2)",
    ),
    dbbase: Path = typer.Option(
        EXTERNAL_DATA_DIR / "ColabFold_db",
        help="Directory containing ColabFold MMseqs2 databases",
    ),
    work_dir: Path = typer.Option(
        INTERIM_DATA_DIR / "piRNA_PPI_20260210" / "paired_msa_work",
        help="Working directory for MMseqs2 intermediate files",
    ),
    output_dir: Path = typer.Option(
        INTERIM_DATA_DIR / "piRNA_PPI_20260210" / "paired_msa_output",
        help="Output directory (contains log, index, summary, paired_msa/)",
    ),
    mmseqs: Path = typer.Option(Path("mmseqs"), help="Path to mmseqs binary"),
    threads: int = typer.Option(64, help="Threads for MMseqs2"),
    gpu: int = typer.Option(0, help="Use GPU for MMseqs2 (0/1)"),
    use_env: bool = typer.Option(True, "--use-env/--no-use-env", help="Use environmental DB"),
    db_load_mode: int = typer.Option(2, help="MMseqs2 db-load-mode"),
    pairing_strategy: int = typer.Option(0, help="pairaln pairing strategy (0: greedy, 1: all)"),
    skip_existing: bool = typer.Option(True, "--skip-existing/--no-skip-existing", help="Skip if outputs exist"),
) -> None:
    """Generate paired MSAs using ColabFold's pairaln strategy.

    This command:
    1. Extracts unique proteins from pairs FASTA
    2. Creates MMseqs2 query database with proper lookup for pairing
    3. Runs monomer search (uniref + optional envdb)
    4. Runs pairaln to generate paired MSAs
    5. Unpacks and organizes output files
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    msa_dir = output_dir / MSA_SUBDIR
    msa_dir.mkdir(parents=True, exist_ok=True)
    index_path = output_dir / INDEX_TSV

    # Setup file logger
    log_id = logger.add(
        output_dir / LOG_FILE,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        level="DEBUG",
    )

    try:
        logger.info(f"Pairs file  : {pairs_file}")
        logger.info(f"Database dir: {dbbase}")
        logger.info(f"Work dir    : {work_dir}")
        logger.info(f"Output dir  : {output_dir}")

        # Load pairs and sequences
        pairs, protein_to_seq = load_pairs_and_sequences_from_fasta(pairs_file)
        logger.info("Loaded {} pairs and {} unique proteins", len(pairs), len(protein_to_seq))

        # Create query FASTA with unique sequences (deduplicated)
        unique_seqs: dict[str, list[str]] = {}  # seq -> list of protein IDs with that seq
        for pid, seq in protein_to_seq.items():
            if seq not in unique_seqs:
                unique_seqs[seq] = []
            unique_seqs[seq].append(pid)

        # Map: unique_seq_idx -> (seq, protein_ids)
        seq_idx_map: list[tuple[str, list[str]]] = list(unique_seqs.items())
        # Map: protein_id -> unique_seq_idx
        pid_to_seq_idx = {pid: idx for idx, (_, pids) in enumerate(seq_idx_map) for pid in pids}

        logger.info("Found {} unique sequences from {} proteins", len(seq_idx_map), len(protein_to_seq))

        # Create query.fas with unique sequences
        query_file = work_dir / "query.fas"
        if not query_file.exists() or not skip_existing:
            with query_file.open("w") as f:
                for idx, (seq, _) in enumerate(seq_idx_map):
                    f.write(f">{101 + idx}\n{seq}\n")
            logger.info("Wrote {} unique sequences to {}", len(seq_idx_map), query_file)

        # Create qdb
        qdb = work_dir / "qdb"
        if not qdb.with_suffix(".dbtype").exists():
            run_mmseqs(mmseqs, ["createdb", query_file, qdb, "--shuffle", "0", "--dbtype", "1"])

        # Create qdb.lookup for pairing
        # Format: id<TAB>pair_name<TAB>file_number
        # Each pair needs two entries (one for each protein)
        lookup_file = work_dir / "qdb.lookup"
        if not lookup_file.exists() or not skip_existing:
            with lookup_file.open("w") as f:
                for file_number, (p1, p2) in enumerate(pairs):
                    idx1 = pid_to_seq_idx[p1]
                    idx2 = pid_to_seq_idx[p2]
                    pair_name = f"{p1}_{p2}"
                    f.write(f"{idx1}\t{pair_name}\t{file_number}\n")
                    f.write(f"{idx2}\t{pair_name}\t{file_number}\n")
            logger.info("Wrote lookup file with {} pairs", len(pairs))

        # Run monomer search
        logger.info("Running monomer search...")
        mmseqs_search_monomer(
            dbbase=dbbase,
            base=work_dir,
            mmseqs=mmseqs,
            use_env=use_env,
            db_load_mode=db_load_mode,
            threads=threads,
            gpu=gpu,
        )

        # Unpack monomer MSAs
        if not (work_dir / "0.a3m").exists():
            run_mmseqs(mmseqs, ["unpackdb", work_dir / "final.a3m", work_dir, "--unpack-name-mode", "0", "--unpack-suffix", ".a3m"])

        # Run pairaln
        logger.info("Running pairaln for complex pairing...")
        mmseqs_search_pair(
            dbbase=dbbase,
            base=work_dir,
            mmseqs=mmseqs,
            prefilter_mode=0,
            s=8,
            threads=threads,
            gpu=gpu,
            db_load_mode=db_load_mode,
            pairing_strategy=pairing_strategy,
        )

        # Load pair.a3m index: each pair produces 2 entries (one per protein),
        # stored sequentially as row 2*file_number and 2*file_number+1.
        # Bypass mmseqs unpackdb which crashes (SIGSEGV) when multiple index
        # entries share the same seq_id and are written concurrently.
        pair_index: list[tuple[int, int]] = []  # (offset, length) per index row
        pair_index_file = work_dir / "pair.a3m.index"
        with pair_index_file.open() as _f:
            for _line in _f:
                parts = _line.split()
                pair_index.append((int(parts[1]), int(parts[2])))
        logger.info("Loaded {} entries from pair.a3m index", len(pair_index))

        # Combine unpaired + paired MSAs and write to output directory
        logger.info("Combining unpaired and paired MSAs...")
        success_count = 0
        low_depth_count = 0

        with index_path.open("w") as idx_f, open(work_dir / "pair.a3m", "rb") as pair_fh:
            for file_number, (p1, p2) in enumerate(tqdm(pairs, desc="Generating final paired MSAs")):
                idx1 = pid_to_seq_idx[p1]
                idx2 = pid_to_seq_idx[p2]

                # Read monomer MSAs
                mono1_file = work_dir / f"{idx1}.a3m"
                mono2_file = work_dir / f"{idx2}.a3m"

                out_file = msa_dir / f"{p1}_vs_{p2}.a3m"
                len1 = get_query_len(protein_to_seq[p1])
                idx_f.write(f"{out_file.resolve()}\t{len1}\n")

                if skip_existing and out_file.exists():
                    success_count += 1
                    continue

                # Combine MSAs using msa_to_str logic
                try:
                    unpaired_msa = []
                    if mono1_file.exists():
                        unpaired_msa.append(mono1_file.read_text())
                    if mono2_file.exists():
                        unpaired_msa.append(mono2_file.read_text())

                    # Read paired MSA directly from pair.a3m via index
                    # Row 2*file_number -> protein1, row 2*file_number+1 -> protein2
                    paired_msa = []
                    row1 = file_number * 2
                    row2 = file_number * 2 + 1
                    for row in (row1, row2):
                        if row < len(pair_index):
                            offset, length = pair_index[row]
                            pair_fh.seek(offset)
                            paired_msa.append(pair_fh.read(length).decode("utf-8", errors="replace"))

                    # Write combined MSA
                    query_seqs = [protein_to_seq[p1], protein_to_seq[p2]]
                    msa_content = _msa_to_str(unpaired_msa, paired_msa, query_seqs)
                    out_file.write_text(msa_content)

                    # Count paired sequences
                    num_paired = sum(1 for line in msa_content.splitlines() if line.startswith(">") and "\t" in line) - 1
                    success_count += 1
                    if num_paired < 10:
                        low_depth_count += 1
                except Exception as e:
                    logger.warning(f"Failed to process {p1}_{p2}: {e}")

        # Write summary
        _write_summary(
            output_dir=output_dir,
            n_pairs=len(pairs),
            n_proteins=len(protein_to_seq),
            n_unique_seqs=len(seq_idx_map),
            success_count=success_count,
            low_depth_count=low_depth_count,
        )
        logger.info(f"Written: {SUMMARY_TSV}")

        if low_depth_count > 0:
            logger.warning(f"{low_depth_count} pairs have <10 paired sequences (low MSA depth)")
        logger.success(f"Generated {success_count} paired MSA files in {msa_dir}")
        logger.success(f"Wrote index file: {index_path}")
    finally:
        logger.remove(log_id)


def _msa_to_str(
    unpaired_msa: list[str],
    paired_msa: list[str],
    query_seqs: list[str],
) -> str:
    """Combine unpaired and paired MSAs into ColabFold format.

    Adapted from ColabFold's msa_to_str function.
    """
    # Parse all MSAs
    unpaired_seqs = []
    for msa_str in unpaired_msa:
        seqs, _ = parse_a3m(msa_str)
        unpaired_seqs.append(seqs)

    paired_seqs = []
    for msa_str in paired_msa:
        seqs, headers = parse_a3m(msa_str)
        paired_seqs.append((seqs, headers))

    # Calculate lengths
    Ls = [get_query_len(seq) for seq in query_seqs]

    lines = []
    # NOTE: Removed ColabFold-style header line (e.g., "#L1,L2\t1,1")
    # because RoseTTAFold2-PPI does not skip # lines and treats them as sequences

    # Query header and sequence
    query_header = "\t".join([str(101 + i) for i in range(len(query_seqs))])
    lines.append(f">{query_header}")
    lines.append("".join(query_seqs))

    # Paired sequences (with tab-separated headers)
    for seqs, headers in paired_seqs:
        for seq, header in zip(seqs[1:], headers[1:]):  # Skip query
            if "\t" in header:
                lines.append(f">{header}")
                lines.append(seq)

    # Unpaired sequences (gap-padded)
    for chain_idx, seqs in enumerate(unpaired_seqs):
        gap_before = sum(Ls[:chain_idx]) * "-"
        gap_after = sum(Ls[chain_idx + 1:]) * "-"
        for seq in seqs[1:]:  # Skip query
            lines.append(f">{101 + chain_idx}")
            lines.append(f"{gap_before}{seq}{gap_after}")

    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    app()