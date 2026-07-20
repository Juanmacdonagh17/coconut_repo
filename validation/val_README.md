# Validation

Numerical checks that Coconut reproduces the codon-usage metrics of established,
independent tools. All checks run on one example human transcript
(HBB, `../examples/example.fasta`) using the codon tables shipped in `../examples/`.

> These validate the **correctness of the descriptive metrics** against reference
> implementations, on an example gene. They are a per-metric correctness floor —
> not a performance benchmark.

## Prerequisites

Compile Coconut in the repository root first (the notebooks call `../coconut`):

```bash
gcc -o coconut coconut.c -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lcurl -lm
```

- `benchmark.Rmd` — R with **seqinr** (`install.packages("seqinr")`).
- `benchmark_mm.ipynb` — Python with **CodonTransformer** (`pip install CodonTransformer`),
  plus `numpy` and `matplotlib`.

## Notebooks

### `benchmark.Rmd` — CU, RSCU, CAI

Computes Codon Usage, RSCU and CAI with **seqinr** (`uco`, `cai`) and compares
against Coconut (`-all`, `-cai`) on the same sequence and the same codon table.

Result on HBB: CU and RSCU match to rounding (max|Δ| ≈ 5e-4); CAI matches
(seqinr 0.8698 vs Coconut 0.8700).

### `benchmark_mm.ipynb` — %MinMax and RRT

- **%MinMax** vs. CodonTransformer's `get_min_max_percentage` (same codon table,
  window 18): the positional profiles match to Coconut's output precision
  (Pearson r = 1.0000). It follows the %MinMax framework (Rodriguez et al. 2018).

## Notes

- The `*_coconut*.csv` and `codon_metrics*.csv` files are produced by the notebooks
  at run time and do not need to be committed.
- Everything is validated on a single gene; scaling the checks over many CDS is a
  natural next step.

## References

- seqinr — Charif & Lobry (2007).
- CodonTransformer — https://github.com/adibvafa/CodonTransformer
- %MinMax / null model — Rodriguez, Wright, Emrich & Clark (2018), *Protein Science*.
  https://doi.org/10.1002/pro.3336