# rustscenic validation on Ziegler 2021 nasopharyngeal scRNA-seq

Eat-your-own-dogfood validation of rustscenic on a real atlas-scale dataset
(32,588 cells × 32,871 genes, 18 airway cell types, 58 donors, balanced
COVID+/- cohort from Ziegler et al. 2021 *Cell*).

This repo is **private, research-only**. Its purpose is (a) to flush out v0.1
user-experience issues before rustscenic ships publicly, (b) to produce a
biological result we can cite alongside the rustscenic preprint, and (c) to
extend the covid-airway-deconvolution project from 'which cells are present
in COVID samples' to 'which TFs regulate those cells'.

Dataset source: Ziegler et al. 2021, GEO-equivalent BCH/Broad/UMMC release.
Local copy at ../covid-airway-deconvolution/data/ziegler2021_nasopharyngeal.h5ad.

## Structure
- `scripts/` — the actual rustscenic pipeline
- `results/` — outputs (gitignored; reproducible from scripts)
- `friction_log.md` — user-experience issues logged as we hit them. Each
  entry links to a candidate rustscenic issue.
- `figures/` — plots for the case study

