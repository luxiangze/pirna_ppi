from pathlib import Path
import subprocess
import os
from Bio import SeqIO

def run_cmd(cmd: str, outputs: list[Path], force: bool = False):
    for o in outputs:
        if o.exists() and not force:
            print(f"⏭️ {o} exists, skip")
            return
        else:
            o.parent.mkdir(parents=True, exist_ok=True)

    subprocess.run(cmd, shell=True, check=True)

def get_len1(first_gene_name:str, protein_seqs: Path) -> int:
    protein_seqs_len = protein_seqs.with_name(f"{protein_seqs.stem}_length.tsv")
    if not protein_seqs_len.exists():
        # Generate the length file if it doesn't exist
        with open(protein_seqs, 'r') as handle, open(protein_seqs_len, 'w') as out_handle:
            for record in SeqIO.parse(handle, "fasta"):
                out_handle.write(f"{record.id}\t{len(record.seq)}\n")

    found_any = False
    with open(protein_seqs_len, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            found_any = True
            gene, length_str = line.split("\t", 1)
            if gene == first_gene_name:
                return int(length_str)

    if not found_any:
        raise ValueError(f"Length file {protein_seqs_len} is empty")
    raise ValueError(f"Gene {first_gene_name} not found in length file {protein_seqs_len}")
