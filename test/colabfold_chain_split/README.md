# ColabFold Chain-Split Test (200 aa gap)

Verifies whether prepending `#L1,L2\t1,1` to the input a3m makes ColabFold treat
the paired MSA as a true two-chain complex (and so insert a 200 aa residue-index
gap), instead of a single concatenated monomer.

## What it tests

A 2x2 controlled experiment on a single small pair
(`KWMTBOMO00189-seg1_vs_KWMTBOMO15589-seg2`, L1=155, L2=25):

| Variant      | Network              | Expectation                                                         |
|--------------|----------------------|---------------------------------------------------------------------|
| `baseline`   | `alphafold2_ptm`     | Header `#180\t1`; pickle len 189 (9-G linker); border contact ~ 1   |
| `baseline`   | `alphafold2_multimer_v3` | Header `#180\t1`; chain split absent; border contact still high |
| `fixed`      | `alphafold2_ptm`     | Header `#155,25\t1,1`; pickle len 189; border contact lower         |
| `fixed`      | `alphafold2_multimer_v3` | Header `#155,25\t1,1`; PDB has chains A,B; chain B residue index gap = 200 |

## Run

```bash
# Run all 4 ColabFold jobs sequentially (~20 min on a 4090):
bash run_test.sh

# Compare results:
pixi run python analyze.py
```

`analyze.py` prints a Polars table covering the columns above and writes
`summary.json`.

## Pass/fail criteria

The fix is considered working if all four hold:

1. `fixed_*` runs have output a3m header `#155,25\t1,1`.
2. `fixed_AFmm` produces a PDB with two chains (`A`, `B`).
3. Chain B residue index gap (`B[0] - A[-1]`) is 200 in `fixed_AFmm`.
4. `border_(L1-1, B0)` contact prob in `fixed_AFmm` is materially lower than in
   `baseline_AFmm` (e.g. drops below ~ 0.5).

## File layout

```text
inputs/{baseline,fixed}/<pair>.a3m   <- 2 input variants
run_test.sh                          <- 4 ColabFold invocations (singularity)
analyze.py                           <- comparison script
outputs/{baseline,fixed}_{AF2,AFmm}/ <- ColabFold outputs (created at runtime)
summary.json                         <- analyze.py result
```
