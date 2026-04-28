#!/usr/bin/env bash
# Run 4 ColabFold jobs (2 input variants x 2 networks) to verify
# whether the `#L1,L2\t1,1` header drives ColabFold to add a 200 aa
# residue-index gap between the two chains.
#
# Params mirror notebook cell 22.
set -euo pipefail

SIF=/home/gyk/sifs/colabfold_1.5.5-cuda12.2.2.sif
CACHE=/home/gyk/alphafold/cache
WORK=/home/gyk/projects/pirna_ppi/test/colabfold_chain_split

run_one() {
  local variant="$1"   # baseline | fixed
  local model_type="$2"  # alphafold2_ptm | alphafold2_multimer_v3
  local tag="$3"       # AF2 | AFmm
  local out="$WORK/outputs/${variant}_${tag}"
  mkdir -p "$out"
  echo "[run_test] $(date '+%H:%M:%S') START variant=$variant tag=$tag"
  singularity run --nv \
    -B "$CACHE":/cache \
    -B "$WORK":/work \
    "$SIF" \
    colabfold_batch \
      --model-order 3,5,1,2,4 \
      --num-models 5 \
      --num-recycle 3 \
      --save-all \
      --model-type "$model_type" \
      "/work/inputs/${variant}/" "/work/outputs/${variant}_${tag}/"
  echo "[run_test] $(date '+%H:%M:%S') DONE  variant=$variant tag=$tag"
}

run_one baseline alphafold2_ptm         AF2
run_one baseline alphafold2_multimer_v3 AFmm
run_one fixed    alphafold2_ptm         AF2
run_one fixed    alphafold2_multimer_v3 AFmm

echo "[run_test] All 4 runs complete."
