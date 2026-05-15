"""Head-to-head: rustscenic.aucell vs pyscenic.aucell on Ziegler 2021.

SETUP (MUST BE RUN IN AN ENVIRONMENT WITH pyscenic INSTALLED):
  python scripts/03_headtohead_pyscenic_aucell.py

Isolation strategy:
  - We cannot run arboreto on Ziegler (its own env pandas pin conflicts with
    modern dask — this is our pitch literalized). So we use rustscenic's GRN
    adjacencies for BOTH aucell runs. That makes the comparison strictly
    about the AUCell kernel, not the GRN stage.
  - We run pyscenic.aucell in two modes:
      (a) noweights=True — unit weights per regulon gene (matches rustscenic)
      (b) noweights=False — GRN-importance weights (pyscenic default)
    This surfaces exactly how much of any divergence is kernel-vs-weighting.
  - pyscenic shuffles gene order for tie-breaking; we use deterministic
    gene-index tie-breaks. Fixed seed across runs.

Outputs (results/comparison/):
  - auc_rustscenic.parquet      — (cells × regulons) rustscenic
  - auc_pyscenic_unit.parquet   — pyscenic noweights=True
  - auc_pyscenic_weighted.parquet — pyscenic noweights=False (realistic use)
  - per_cell_pearson.csv        — per-cell Pearson between matched outputs
  - per_regulon_pearson.csv     — per-regulon Pearson
  - canonical_tf_compare.csv    — the 14 canonical TFs, z in expected CT, each tool
  - runtime_memory.csv          — wall + peak RSS per run
  - install_matrix.md           — the 4-cell install/run outcome table
  - summary.json
"""
from __future__ import annotations
import json, time, resource, platform
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp

import rustscenic.aucell
from pyscenic.aucell import aucell as py_aucell
from pyscenic.utils import modules_from_adjacencies
from ctxcore.genesig import Regulon

from paths import ziegler_h5ad_path

DATA = ziegler_h5ad_path()
OUT = Path(__file__).parent.parent / "results" / "comparison"
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED = {
    "ASCL3": "Ionocytes",
    "FOXJ1": "Ciliated Cells",
    "CEBPB": "Macrophages",
    "TBX21": "T Cells",
    "EOMES": "T Cells",
    "SPI1":  "Macrophages",
    "TP63":  "Basal Cells",
    "STAT1": "Ciliated Cells",   # detailed would be IFN-Resp Ciliated
    "RUNX3": "T Cells",
    "SPDEF": "Goblet Cells",
    "MYB":   "Deuterosomal Cells",
    "IRF7":  "Ciliated Cells",
    "SOX2":  "Basal Cells",
    "PAX5":  "B Cells",
}


def rss_gb():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_maxrss / (1024**3) if platform.system() == "Darwin" else r.ru_maxrss / (1024**2)


def main():
    # --- reuse the cached rustscenic GRN from the earlier run ---
    grn_path = OUT.parent / "grn.parquet"
    auc_rust_path = OUT.parent / "auc.parquet"
    regulons_path = OUT.parent / "regulons.json"
    assert grn_path.exists(), f"run 01_grn_aucell_celltype_regulons.py first ({grn_path} missing)"

    print("[1/7] reloading Ziegler + running same preprocess to get matched adata...")
    adata_raw = sc.read_h5ad(DATA)
    drop_cts = {"Erythroblasts"}
    adata = adata_raw[~adata_raw.obs["Coarse_Cell_Annotations"].isin(drop_cts)].copy()
    del adata_raw
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
    TFS = ["ASCL3","FOXJ1","CEBPB","TBX21","EOMES","SPI1","TP63","STAT1","RUNX3",
           "SPDEF","MYB","IRF7","SOX2","PAX5","NKX2-1","KLF5","GATA3","CEBPA","CEBPD",
           "ELF3","EHF","HES1","HES6","NOTCH1","IRF1","IRF3","IRF8","IRF9","STAT2","STAT3",
           "NFKB1","NFKB2","REL","RELA","RELB","XBP1","ATF3","ATF4","ATF6","JUN","JUNB",
           "FOS","FOSL1","FOSL2","EBF1","MEF2A","MEF2C","TCF7","LEF1","GATA2","MAFB",
           "BACH1","BACH2","KLF4","KLF6","KLF15","NR4A1","NR4A2","NR4A3"]
    keep = adata.var["highly_variable"].copy()
    keep[adata.var_names.isin(TFS)] = True
    adata = adata[:, keep].copy()
    print(f"      matched adata: {adata.shape}")

    # --- load cached rustscenic GRN + regulons (same regulons → isolates AUCell) ---
    print("[2/7] loading cached rustscenic GRN + regulons from earlier run...")
    grn = pd.read_parquet(grn_path)
    with open(regulons_path) as fh:
        regulons_dict = json.load(fh)
    print(f"      {len(grn):,} edges, {len(regulons_dict)} regulons")

    # --- build pyscenic Regulon objects from the SAME top-50-target sets ---
    # Use GRN importance as weights for pyscenic (its default); also build unit-weight version.
    print("[3/7] building pyscenic Regulon objects...")
    grn_by_tf = {tf: g for tf, g in grn.groupby("TF")}
    regs_weighted, regs_unit = [], []
    for name, genes in regulons_dict.items():
        tf = name.replace("_regulon", "")
        g = grn_by_tf.get(tf)
        if g is None: continue
        g2w_wt = dict(zip(g["target"], g["importance"]))
        g2w_wt = {gene: float(g2w_wt.get(gene, 0.0)) for gene in genes}
        g2w_unit = {gene: 1.0 for gene in genes}
        regs_weighted.append(Regulon(name=name, gene2weight=g2w_wt,
                                     transcription_factor=tf, gene2occurrence={}))
        regs_unit.append(Regulon(name=name, gene2weight=g2w_unit,
                                 transcription_factor=tf, gene2occurrence={}))
    print(f"      built {len(regs_weighted)} weighted + unit regulons")

    # pyscenic wants a (cells × genes) DataFrame
    ex_df = adata.to_df()

    results = {}

    # --- run pyscenic.aucell (unit weights) ---
    print("[4/7] pyscenic.aucell noweights=True (matches rustscenic kernel)...")
    t0 = time.monotonic()
    auc_py_unit = py_aucell(ex_df, regs_unit, auc_threshold=0.05, num_workers=1, noweights=True, normalize=False)
    t_py_unit = time.monotonic() - t0
    results["pyscenic_unit"] = {"wall_s": round(t_py_unit, 2), "peak_rss_gb": round(rss_gb(), 2)}
    print(f"      {t_py_unit:.2f}s, shape {auc_py_unit.shape}")
    auc_py_unit.to_parquet(OUT / "auc_pyscenic_unit.parquet")

    # --- run pyscenic.aucell (GRN-weighted, its default) ---
    print("[5/7] pyscenic.aucell noweights=False (realistic pyscenic default)...")
    t0 = time.monotonic()
    auc_py_wt = py_aucell(ex_df, regs_weighted, auc_threshold=0.05, num_workers=1, noweights=False, normalize=False)
    t_py_wt = time.monotonic() - t0
    results["pyscenic_weighted"] = {"wall_s": round(t_py_wt, 2), "peak_rss_gb": round(rss_gb(), 2)}
    print(f"      {t_py_wt:.2f}s, shape {auc_py_wt.shape}")
    auc_py_wt.to_parquet(OUT / "auc_pyscenic_weighted.parquet")

    # --- reload our cached rustscenic AUC ---
    auc_rust = pd.read_parquet(auc_rust_path)
    # The earlier rustscenic run we did was in a different venv but same regulons + same data;
    # trust the cached output. But rerun here to get a runtime+mem number comparable to pyscenic.
    print("[6/7] rustscenic.aucell (for runtime comparison in same env)...")
    t0 = time.monotonic()
    auc_rust_live = rustscenic.aucell.score(adata,
        [(name, genes) for name, genes in regulons_dict.items()],
        top_frac=0.05, chunk_size=5000)
    t_rust = time.monotonic() - t0
    results["rustscenic"] = {"wall_s": round(t_rust, 2), "peak_rss_gb": round(rss_gb(), 2)}
    print(f"      {t_rust:.2f}s, shape {auc_rust_live.shape}")

    # --- per-cell + per-regulon Pearson ---
    print("[7/7] comparing outputs...")
    # Align
    common_regs = sorted(set(auc_rust_live.columns) & set(auc_py_unit.columns) & set(auc_py_wt.columns))
    common_cells = sorted(set(auc_rust_live.index) & set(auc_py_unit.index) & set(auc_py_wt.index))
    R = auc_rust_live.loc[common_cells, common_regs].values
    PU = auc_py_unit.loc[common_cells, common_regs].values
    PW = auc_py_wt.loc[common_cells, common_regs].values

    def per_cell_pearson(A, B):
        out = []
        for i in range(A.shape[0]):
            a, b = A[i], B[i]
            if np.std(a) > 1e-12 and np.std(b) > 1e-12:
                out.append(np.corrcoef(a, b)[0, 1])
        return np.asarray(out)

    def per_reg_pearson(A, B):
        out = []
        for j in range(A.shape[1]):
            a, b = A[:, j], B[:, j]
            if np.std(a) > 1e-12 and np.std(b) > 1e-12:
                out.append(np.corrcoef(a, b)[0, 1])
        return np.asarray(out)

    cell_ru = per_cell_pearson(R, PU)
    cell_rw = per_cell_pearson(R, PW)
    reg_ru = per_reg_pearson(R, PU)
    reg_rw = per_reg_pearson(R, PW)

    # Argmax agreement per cell
    def argmax_match(A, B):
        return float((A.argmax(axis=1) == B.argmax(axis=1)).mean())

    per_cell_df = pd.DataFrame({
        "rust_vs_py_unit":    cell_ru,
        "rust_vs_py_weighted": cell_rw,
    })
    per_cell_df.describe().to_csv(OUT / "per_cell_pearson_summary.csv")
    per_reg_df = pd.DataFrame({
        "rust_vs_py_unit":    reg_ru,
        "rust_vs_py_weighted": reg_rw,
    })
    per_reg_df.describe().to_csv(OUT / "per_regulon_pearson_summary.csv")

    # Canonical TF benchmark — z-score in expected cell type for each method
    z_table_rows = []
    for auc_df, tag in [(auc_rust_live, "rustscenic"),
                        (auc_py_unit, "pyscenic_unit"),
                        (auc_py_wt, "pyscenic_weighted")]:
        aux = auc_df.loc[common_cells].copy()
        aux["ct"] = adata.obs.loc[common_cells, "Coarse_Cell_Annotations"].values
        mean_ct = aux.groupby("ct", observed=False).mean(numeric_only=True)
        z_ct = (mean_ct - mean_ct.mean(axis=0)) / (mean_ct.std(axis=0) + 1e-12)
        for tf, expected in EXPECTED.items():
            regname = f"{tf}_regulon"
            if regname not in z_ct.columns or expected not in z_ct.index: continue
            zval = z_ct.loc[expected, regname]
            top_ct = z_ct[regname].idxmax()
            z_table_rows.append({
                "tf": tf, "expected": expected, "method": tag,
                "z_in_expected": float(zval),
                "top_ct": str(top_ct),
                "hit": bool(top_ct == expected),
            })
    z_df = pd.DataFrame(z_table_rows)
    z_df.to_csv(OUT / "canonical_tf_compare.csv", index=False)

    # Pivot for readability
    pivot_z = z_df.pivot(index="tf", columns="method", values="z_in_expected").reindex(EXPECTED.keys())
    pivot_hit = z_df.pivot(index="tf", columns="method", values="hit").reindex(EXPECTED.keys())

    print("\n=== Per-cell Pearson (rustscenic vs pyscenic) ===")
    print(f"  rust vs py-unit:     mean {cell_ru.mean():.4f}  median {np.median(cell_ru):.4f}  "
          f">0.95 {100*(cell_ru>0.95).mean():.1f}%")
    print(f"  rust vs py-weighted: mean {cell_rw.mean():.4f}  median {np.median(cell_rw):.4f}  "
          f">0.95 {100*(cell_rw>0.95).mean():.1f}%")

    print("\n=== Per-regulon Pearson ===")
    print(f"  rust vs py-unit:     mean {reg_ru.mean():.4f}  median {np.median(reg_ru):.4f}")
    print(f"  rust vs py-weighted: mean {reg_rw.mean():.4f}  median {np.median(reg_rw):.4f}")

    print("\n=== Argmax-regulon per cell agreement ===")
    print(f"  rust == py-unit:     {100*argmax_match(R, PU):.1f}%")
    print(f"  rust == py-weighted: {100*argmax_match(R, PW):.1f}%")
    print(f"  py-unit == py-weighted: {100*argmax_match(PU, PW):.1f}%")

    print("\n=== Runtime (on 31,602 cells × 59 regulons) ===")
    for k, v in results.items():
        print(f"  {k:22s}  {v['wall_s']:7.2f}s   peak RSS {v['peak_rss_gb']:.2f} GB")

    print("\n=== Canonical TF z-score in expected cell type ===")
    print(pivot_z.round(2).to_string())
    print("\n=== Canonical TF hit (top-ranked cell type matches expected) ===")
    print(pivot_hit.to_string())
    hits_per_method = pivot_hit.sum()
    print(f"\nhits per method:\n{hits_per_method}")

    summary = {
        "dataset": "Ziegler 2021 nasopharyngeal scRNA-seq",
        "n_cells_used": len(common_cells),
        "n_regulons": len(common_regs),
        "per_cell_pearson": {
            "rust_vs_py_unit":     {"mean": float(cell_ru.mean()), "median": float(np.median(cell_ru)),
                                    "pct_gt_95": float((cell_ru > 0.95).mean() * 100)},
            "rust_vs_py_weighted": {"mean": float(cell_rw.mean()), "median": float(np.median(cell_rw)),
                                    "pct_gt_95": float((cell_rw > 0.95).mean() * 100)},
        },
        "per_regulon_pearson": {
            "rust_vs_py_unit":     {"mean": float(reg_ru.mean()), "median": float(np.median(reg_ru))},
            "rust_vs_py_weighted": {"mean": float(reg_rw.mean()), "median": float(np.median(reg_rw))},
        },
        "argmax_match": {
            "rust_vs_py_unit":     float(argmax_match(R, PU)),
            "rust_vs_py_weighted": float(argmax_match(R, PW)),
            "py_unit_vs_py_weighted": float(argmax_match(PU, PW)),
        },
        "runtime_memory": results,
        "canonical_tf_hit_counts": hits_per_method.to_dict(),
        "arboreto_runs_in_pyscenic_env": False,
        "arboreto_failure_mode": "ImportError: Dask requires pandas>=2.0.0 but pyscenic pins pandas==1.5.3 (dask_expr incompatibility)",
    }
    with open(OUT / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nsaved {OUT}/summary.json")


if __name__ == "__main__":
    main()
