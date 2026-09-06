# Gene regulation in COVID-19 airway cells

[![CI](https://github.com/Ekin-Kahraman/rustscenic-airway-case/actions/workflows/ci.yml/badge.svg)](https://github.com/Ekin-Kahraman/rustscenic-airway-case/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20230540.svg)](https://doi.org/10.5281/zenodo.20230540)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This [RustScenic](https://github.com/Ekin-Kahraman/rustscenic) case study combines
software validation with an exploratory analysis of gene regulation in a
published airway atlas (Ziegler et al. 2021, 58 donors).

1. **Do the activity scores agree with pySCENIC?** Compare both tools using the same expression data and candidate gene sets.
2. **Which gene programmes differ with COVID-19 status?** Compare activity within airway cell types, extending the [cell-type proportion analysis](https://github.com/Ekin-Kahraman/covid-airway-deconvolution).

The source h5ad is not committed. To rerun the analysis, place the Ziegler file locally and set `ZIEGLER_H5AD=/path/to/ziegler2021_nasopharyngeal.h5ad`, or keep the sibling-repo default at `../covid-airway-deconvolution/data/ziegler2021_nasopharyngeal.h5ad`.

## Tool validation

Against pyscenic on the same 31,602 cells × 59 regulons on the same machine:

A regulon is a candidate set of genes associated with a transcription factor.
The two pySCENIC settings use either equal gene weights or network-derived weights.

| | rustscenic | pyscenic-unit | pyscenic-weighted |
|---|---:|---:|---:|
| Per-cell Pearson vs rustscenic | 1.000 | **0.984** | 0.949 |
| % cells > 0.95 | - | 91.7 % | 71.6 % |
| Canonical TF hits (of 14) | 8 | 8 | 9 |
| AUCell wall-time | **0.25 s** | 6.81 s | 5.29 s |

The methods share five missed expected transcription factors. This comparison
checks activity scoring on fixed gene sets, not agreement of complete network
inference or proof that every biological prediction is correct. See
[`CASE_STUDY.md`](CASE_STUDY.md) for the full comparison.

## Biological analysis

The analysis compares COVID-positive and negative cells within 11 cell types:

- Higher interferon-associated activity (IRF7 gene set) in 7 of 11 cell types.
- Lower AP-1 and stress-associated activity in squamous cells.
- Higher activity of several gene sets associated with regeneration in secretory cells.

These are exploratory associations. Tests compare individual cells and adjust
for multiple comparisons within each cell type; they do **not** model donor
dependence. Donor-level confirmation is needed before interpreting the reported
significance as evidence of disease effects or regulatory mechanisms.

See [`figures/fig6_covid_differential.png`](figures/fig6_covid_differential.png).

## Structure

- `scripts/01_grn_aucell_celltype_regulons.py` - full rustscenic pipeline on Ziegler
- `scripts/02_make_figures.py` - single-tool figures (heatmap + canonical TF benchmark)
- `scripts/03_headtohead_pyscenic_aucell.py` - head-to-head vs pyscenic.aucell (uses scenic-env with pyscenic installed)
- `scripts/04_comparison_figures.py` - 3-way comparison + runtime + Pearson distribution
- `scripts/05_covid_differential_regulons.py` - per-celltype COVID+/− Wilcoxon + BH-FDR
- `CASE_STUDY.md` - long-form writeup linking validation + biology
- `friction_log.md` - 9 first-user observations fed back to rustscenic v0.2
- `figures/` - 6 PNGs (committed; small)
- `results/` - intermediate parquets/CSVs (gitignored; reproducible from scripts)

## Reproduction

Requires:

1. The Ziegler h5ad from the [covid-airway-deconvolution](https://github.com/Ekin-Kahraman/covid-airway-deconvolution) repo
2. rustscenic ≥ 0.4.4 (`pip install -r requirements.txt`)
3. For head-to-head: the pinned pySCENIC reference env in [`reference/`](reference/)

Run in order:

```bash
export ZIEGLER_H5AD=/path/to/ziegler2021_nasopharyngeal.h5ad
python scripts/01_grn_aucell_celltype_regulons.py 2>&1 | tee results/run.log
python scripts/02_make_figures.py
python scripts/03_headtohead_pyscenic_aucell.py
python scripts/04_comparison_figures.py
python scripts/05_covid_differential_regulons.py
python scripts/validate_outputs.py
```

To verify only the cross-tool AUCell environment:

```bash
python -m pip install -r reference/requirements-pyscenic-aucell.txt
python reference/smoke_pyscenic_reference.py
```

## Cross-references

- Tool: [`Ekin-Kahraman/rustscenic`](https://github.com/Ekin-Kahraman/rustscenic) - where the Ziegler numbers + figures also live, as `validation/ziegler_headtohead_2026-04-19.md`
- Source dataset: [`Ekin-Kahraman/covid-airway-deconvolution`](https://github.com/Ekin-Kahraman/covid-airway-deconvolution) - cell-type deconvolution of GSE152075 using Ziegler reference
- Upstream biology: Ziegler et al. 2021 *Cell* - the scRNA-seq atlas this uses

## Licence and citation

MIT. Citation metadata is in [`CITATION.cff`](CITATION.cff).

Zenodo archive:

- Version DOI for `v0.1.1`: [`10.5281/zenodo.20230541`](https://doi.org/10.5281/zenodo.20230541)
- Concept DOI for all versions: [`10.5281/zenodo.20230540`](https://doi.org/10.5281/zenodo.20230540)
