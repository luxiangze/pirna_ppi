from pathlib import Path
import subprocess

def run_cmd(cmd: str, outputs: list[Path]):
    outputs = [Path(o) for o in outputs]

    if all(o.exists() for o in outputs):
        print("⏭️ output exists, skip")
        return

    subprocess.run(cmd, shell=True, check=True)
