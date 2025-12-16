from pathlib import Path
import subprocess

def run_cmd(cmd: str, outputs: list[Path]):
    for o in outputs:
        if o.exists():
            print(f"⏭️ {o} exists, skip")
            return
        else:
            o.parent.mkdir(parents=True, exist_ok=True)

    subprocess.run(cmd, shell=True, check=True)
