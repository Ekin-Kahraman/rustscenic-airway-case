# rustscenic on Ziegler 2021 nasopharyngeal atlas - validation + biology

**Companion repo to** [`Ekin-Kahraman/rustscenic`](https://github.com/Ekin-Kahraman/rustscenic) - not a demo. Two deliverables:

1. **Tool validation** - head-to-head of [rustscenic](https://github.com/Ekin-Kahraman/rustscenic) against pyscenic on a published 58-donor, 32,588-cell airway atlas (Ziegler et al. 2021 *Cell*). Same input, same regulons, same env, isolated AUCell kernel. Answers: "does rustscenic produce pyscenic's numbers on real atlas-scale data?"
2. **Biology** - COVID+ vs COVID− differential regulon analysis, extending the [covid-airway-deconvolution](https://github.com/Ekin-Kahraman/covid-airway-deconvolution) project from "which cells are perturbed" to "which regulatory programmes rewire during SARS-CoV-2 infection". Candidate for a standalone paper.

The source h5ad is not committed. To rerun the analysis, place the Ziegler file locally and set `ZIEGLER_H5AD=/path/to/ziegler2021_nasopharyngeal.h5ad`, or keep the sibling-repo default at `../covid-airway-deconvolution/data/ziegler2021_nasopharyngeal.h5ad`.

## Headline (tool validation)

Against pyscenic on the same 31,602 cells × 59 regulons on the same machine:

| | rustscenic | pyscenic-unit | pyscenic-weighted |
|---|---:|---:|---:|
| Per-cell Pearson vs rustscenic | 1.000 | **0.984** | 0.949 |
| % cells > 0.95 | - | 91.7 % | 71.6 % |
| Canonical TF hits (of 14) | 8 | 8 | 9 |
| AUCell wall-time | **0.25 s** | 6.81 s | 5.29 s |

**All three tools miss the same 5 TFs** (STAT1, MYB, IRF7, SOX2, PAX5) - the tool-to-tool variation is smaller than the dataset-inherent noise. See [`CASE_STUDY.md`](CASE_STUDY.md) for the full interpretation.

## Headline (biology)

Wilcoxon rank-sum (BH-FDR per cell type) on 11 cell types with ≥100 cells per COVID arm:

- **Type I IFN programme (IRF7 regulon) ↑ in COVID+** across 7/11 cell types - strongest in Ionocytes (+1.34 log₂ FC, q = 6e-17). Canonical antiviral response.
- **AP-1 / stress-response programme ↓ in squamous cells** (JUN, JUNB, NR4A1, XBP1, all log₂ FC < −0.5, q < 1e-90). Consistent with the SARS-CoV-2-induced squamous metaplasia Ziegler reported - the canonical squamous stress programme is suppressed.
- **WNT / regenerative programme ↑ in secretory cells** (TCF7 +1.10, LEF1 +1.12, EOMES +1.34). Consistent with regenerative response to epithelial damage.

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
