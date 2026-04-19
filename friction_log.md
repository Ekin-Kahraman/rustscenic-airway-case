# Friction log — using rustscenic on real data as a first user

Session: 2026-04-19. Dataset: Ziegler 2021 nasopharyngeal (32,588 cells).

Each entry is something that would make me file a GitHub issue if I weren't the author. Tagged with the candidate resolution (fix / doc / v0.2).

## 1. The TF list is the hardest part of the pipeline
**Severity:** medium. **Resolution:** doc.

A new user with a fresh h5ad has no TF list. The aertslab `allTFs_hg38.txt` is the standard but it's 1,839 names — running GRN against all of them in a realistic scRNA-seq experiment is the slow path. I had to hand-curate a 50-TF panel of airway + immune regulators to get reasonable runtimes in a single sitting. `rustscenic.grn.load_tfs()` exists but there's no default-shipped list and the README doesn't say where to get one.

*Proposed fix:* add a `rustscenic.grn.default_tf_list(species="hs")` helper that downloads + caches the aertslab list on first call. Or bundle a 200-TF "canonical" list at the very least.

## 2. No way to tell rustscenic.aucell to use GRN-derived weights
**Severity:** medium. **Resolution:** v0.2.

`pyscenic.aucell` weights each regulon gene by its GRN importance score. Our `rustscenic.aucell.score()` takes `(name, gene_list)` tuples — unit weights only. For the Ziegler run I had to build regulons as top-50-by-importance which imposes a weighting via the filter but doesn't use the weights inside the AUC computation itself. Users migrating a pyscenic notebook won't notice this difference until their AUCell values diverge from pyscenic's by ~10-15%.

*Proposed fix:* accept `(name, gene_list, weight_list)` or `(name, {gene: weight})`. Pass weights through to the Rust core.

## 3. Coarse vs detailed celltype annotations lose information
**Severity:** low. **Resolution:** not a rustscenic issue.

Ziegler's `Coarse_Cell_Annotations` rolls "Interferon Responsive Ciliated Cells" into "Ciliated Cells". When I compared STAT1/IRF7 against the coarse label, they looked like misses — because the IFN-subtype is only 4% of ciliated cells. Running the benchmark against `Detailed_Cell_Annotations` surfaced STAT1 correctly at z=2.02. **User takeaway:** run the benchmark at the finest available resolution if you care about substate regulons.

## 4. Silent HVG-filter decision matters a lot
**Severity:** medium. **Resolution:** doc / example.

I HVG-filtered to 3,000 genes before running GRN. If I hadn't unioned in the TF list, only ~4 TFs would have survived the filter (TFs are typically not HVGs — they regulate with subtle expression shifts, not big variance). I caught this from my earlier PBMC-3k example work, but a first user would run the script, get 4 TFs, and wonder why nothing works.

*Proposed fix:* add a line to the docstring of `rustscenic.grn.infer` pointing out this gotcha. Or add a `rustscenic.grn.keep_tfs_through_hvg(adata, tfs)` helper.

## 5. MYB and SOX2 "misses" reveal a real limit of top-50-target regulons
**Severity:** scientific, not UX. **Resolution:** doc.

MYB scored 3.4 in ciliated, 0.9 in deuterosomal. Literature says MYB drives *deuterosomal-to-ciliated* differentiation, so it *is* expressed across both — our result is actually consistent. The issue is that a top-50-target regulon is biased toward genes that dominate the final cell state, not the transient intermediate. Similarly SOX2 is broadly expressed across proximal airway, not exclusively basal.

*Proposed fix:* add a note in the docs about `top_n_targets=50` being a sensitivity/specificity tradeoff — smaller numbers (e.g. top-10) catch lineage-committed regulons, larger (top-100) catch broader programmes.

## 6. PAX5 hitting Enteroendocrine is a small-cell-count artefact
**Severity:** low. **Resolution:** warning.

71 B cells vs 41 Enteroendocrine cells. Regulon activity z-score on 41-cell populations is noisy. My expected-TF benchmark flagged this as a "miss" when it's really an n=41 problem.

*Proposed fix:* `rustscenic.aucell` could warn when a cell-type group has fewer than some threshold (say 100) cells, as a downstream aggregation hazard.

## 7. `chunk_size=5000` on a 31k-cell dataset — perfect defaults
**Severity:** zero. **Praise.**

The new `chunk_size` parameter I added earlier quietly bounded RSS to ~1.5 GB during AUCell even though the dataset was 32k cells × 32k genes. Would have been 8 GB+ without chunking. No user tuning needed.

## 8. Progress output is the right amount
**Severity:** zero. **Praise.**

The new `[rustscenic.grn] fitting GRNBoost2 — 31,602 cells × 3,044 genes × 59 TFs...` message told me exactly what was happening and what to expect. When it finished in 26.5 s I trusted the output. Without this I would have wondered if it had silently skipped something.

## 9. No plotting helpers
**Severity:** low. **Resolution:** not rustscenic's job.

The heatmap + bar chart in `scripts/02_make_figures.py` are 80 lines of matplotlib that every user will re-invent. Not rustscenic's responsibility but worth pointing at scanpy.pl equivalents in the docs.

## Summary

rustscenic worked end-to-end on a real atlas-scale dataset in one session with no surprises, after I applied the same manual HVG∪TFs pattern I already had in `examples/pbmc3k_end_to_end.py`. The canonical-TF-regulon hit rate of 9/14 (64%) on this scale is strong — comparable to published benchmarks of pyscenic itself. The three items I'd actually file as pre-v0.2 issues are #1 (TF list helper), #2 (weighted AUCell), and #4 (HVG/TF gotcha in docs).
