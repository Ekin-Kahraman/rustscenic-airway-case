# Contributing

Contributions are welcome when they improve reproducibility, validation, figure quality, or scientific interpretation.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/validate_outputs.py
python -m py_compile scripts/*.py
```

To rerun the full case study, provide the Ziegler h5ad:

```bash
export ZIEGLER_H5AD=/path/to/ziegler2021_nasopharyngeal.h5ad
python scripts/01_grn_aucell_celltype_regulons.py 2>&1 | tee results/run.log
python scripts/02_make_figures.py
python scripts/03_headtohead_pyscenic_aucell.py
python scripts/04_comparison_figures.py
python scripts/05_covid_differential_regulons.py
python scripts/validate_outputs.py
```

The pySCENIC head-to-head script needs a separate environment with `pyscenic` and `ctxcore` installed.

## Pull request checklist

- `python scripts/validate_outputs.py` passes.
- `python -m py_compile scripts/*.py` passes.
- Changed figures are regenerated intentionally.
- Any changed headline metric is reflected in `README.md` and `CASE_STUDY.md`.
- Validation claims include command, version, dataset, hardware where relevant, and caveats.

## Scientific correctness

Do not change headline values silently. If a change affects canonical TF hits, Pearson agreement, runtime, or COVID differential regulons, include the before and after values in the pull request.

## Licence

By contributing, you agree that your contributions are licensed under the MIT licence.
