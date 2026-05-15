"""Stage 1: GRN + AUCell on Ziegler 2021 airway scRNA-seq (32,588 cells).

Scientific question: which transcription factors drive each airway cell type?
Known benchmarks from the airway biology literature we'll check against:
  - FOXJ1 — master regulator of multiciliated cells
  - SPDEF — master regulator of goblet cells / mucus secretion
  - TP63 — basal cell identity
  - SCGB1A1 (CC10) target network — secretory club cells (SCGB1A1 is a
    marker, but CEBPA/CEBPB drive it)
  - SOX2 — general airway epithelial progenitor
  - ASCL3 — ionocyte specification
  - MYB — deuterosomal / ciliogenesis intermediate
  - STAT1/IRF1/IRF7/ISGs — interferon-responsive ciliated cells
  - RUNX3/TBX21/EOMES — T cell cytotoxic programme
  - SPI1/CEBPB — myeloid (macrophages, DCs)
  - PAX5/EBF1 — B cells

Pipeline:
  1. Load raw Ziegler h5ad (~32.5k cells, 18 coarse cell types)
  2. Subset to epithelial + immune compartments (drop erythroblasts for clarity)
  3. HVG-filter to ~3000 genes UNION with the TF candidate list
  4. rustscenic.grn.infer on HVG × TF subset
  5. Build top-50-target regulons per TF
  6. rustscenic.aucell on the full-cell matrix
  7. Mean regulon activity per celltype
  8. Check: does each expected-TF regulon top the expected cell type?
  9. Write per-celltype regulon rankings + the friction log

Run:
  python scripts/01_grn_aucell_celltype_regulons.py
Outputs in results/
"""
from __future__ import annotations
import json, time
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import rustscenic.grn
import rustscenic.aucell

from paths import ziegler_h5ad_path

DATA = ziegler_h5ad_path()
OUT = Path(__file__).parent.parent / "results"
OUT.mkdir(exist_ok=True)

# --- Expected-TF truth set (from airway literature) ----------------------
EXPECTED = {
    # TF: (expected top-ranked cell type, lit reference)
    "FOXJ1":  ("Ciliated Cells", "You & Brody 2007 Physiol Rev — FOXJ1 master regulator of motile cilia"),
    "SPDEF":  ("Goblet Cells",   "Park et al 2007 JCI — SPDEF drives goblet metaplasia"),
    "TP63":   ("Basal Cells",    "Daniely et al 2004 AJP — p63 airway basal cell marker"),
    "SOX2":   ("Basal Cells",    "Que et al 2009 Development — SOX2 maintains proximal airway identity"),
    "MYB":    ("Deuterosomal Cells", "Pan et al 2014 Nature — MYB drives multiciliogenesis"),
    "ASCL3":  ("Ionocytes",      "Plasschaert et al 2018 Nature — ASCL3 ionocyte lineage"),
    "STAT1":  ("Ciliated Cells",  "Canonical IFN-I response"),
    "IRF7":   ("Ciliated Cells",  "Canonical IFN-I response"),
    "SPI1":   ("Macrophages",    "Scott et al 1994 Science — PU.1/SPI1 myeloid master"),
    "CEBPB":  ("Macrophages",    "Akira et al — myeloid/inflammatory"),
    "PAX5":   ("B Cells",        "Urbánek et al 1994 — PAX5 B-lymphoid"),
    "TBX21":  ("T Cells",        "T-bet Th1/NK cytotoxic programme"),  # T Cells collapsed
    "EOMES":  ("T Cells",        "Eomes CD8/NK cytotoxic"),
    "RUNX3":  ("T Cells",        "Egawa et al 2007 — CD8 T cell specification"),
}

# TF candidate list — union of expected TFs + broader airway + immune panel.
TFS_CANDIDATE = sorted(set(list(EXPECTED.keys()) + [
    "NKX2-1", "KLF5", "GATA3", "CEBPA", "CEBPD", "ELF3", "EHF",
    "HES1", "HES6", "NOTCH1", "SCGB3A2",
    "IRF1", "IRF3", "IRF8", "IRF9", "STAT2", "STAT3",
    "NFKB1", "NFKB2", "REL", "RELA", "RELB",
    "XBP1", "ATF3", "ATF4", "ATF6", "JUN", "JUNB", "FOS", "FOSL1", "FOSL2",
    "EBF1", "MEF2A", "MEF2C",
    "TCF7", "LEF1", "GATA2", "MAFB", "BACH1", "BACH2",
    "KLF4", "KLF6", "KLF15",
    "NR4A1", "NR4A2", "NR4A3",
]))


def main():
    print(f"[1/8] loading {DATA.name}...")
    t0 = time.monotonic()
    adata = sc.read_h5ad(DATA)
    print(f"     shape: {adata.shape}  loaded in {time.monotonic()-t0:.1f}s")

    print("[2/8] keeping epithelial + T + Macrophage + B + ionocyte (drop erythroid for clarity)...")
    drop_cts = {"Erythroblasts"}  # nucleated but terminal; noise for TF discovery
    mask = ~adata.obs["Coarse_Cell_Annotations"].isin(drop_cts)
    adata = adata[mask].copy()
    print(f"     after filter: {adata.shape}")

    print("[3/8] scanpy preprocess + HVG∪TFs...")
    # X is raw counts (sparse). Normalise + log.
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
    keep = adata.var["highly_variable"].copy()
    keep[adata.var_names.isin(TFS_CANDIDATE)] = True
    adata = adata[:, keep].copy()
    print(f"     after HVG∪TFs: {adata.shape}")

    tfs_in = [t for t in TFS_CANDIDATE if t in adata.var_names]
    print(f"     TFs in HVG-filtered matrix: {len(tfs_in)}/{len(TFS_CANDIDATE)}")

    print("[4/8] rustscenic.grn.infer...")
    t0 = time.monotonic()
    grn = rustscenic.grn.infer(adata, tf_names=tfs_in, n_estimators=500, seed=777)
    t_grn = time.monotonic() - t0
    print(f"     {len(grn):,} edges in {t_grn:.1f}s")
    grn.to_parquet(OUT / "grn.parquet")

    print("[5/8] building top-50-target regulons per TF...")
    regulons = {}
    for tf in grn["TF"].unique():
        top = grn[grn["TF"] == tf].nlargest(50, "importance")["target"].tolist()
        if len(top) >= 10:
            regulons[f"{tf}_regulon"] = top
    print(f"     regulons: {len(regulons)}")
    with open(OUT / "regulons.json", "w") as fh:
        json.dump(regulons, fh, indent=2)

    print("[6/8] rustscenic.aucell per cell...")
    reg_list = list(regulons.items())
    t0 = time.monotonic()
    auc = rustscenic.aucell.score(adata, reg_list, top_frac=0.05, chunk_size=5000)
    t_auc = time.monotonic() - t0
    print(f"     activity matrix: {auc.shape}  in {t_auc:.1f}s")
    auc.to_parquet(OUT / "auc.parquet")

    print("[7/8] mean regulon activity per cell type...")
    auc_df = auc.copy()
    auc_df["celltype"] = adata.obs["Coarse_Cell_Annotations"].values
    auc_df["detailed"] = adata.obs["Detailed_Cell_Annotations"].values
    mean_per_ct = auc_df.groupby("celltype", observed=False).mean(numeric_only=True)
    mean_per_detailed = auc_df.groupby("detailed", observed=False).mean(numeric_only=True)
    mean_per_ct.to_csv(OUT / "mean_activity_per_coarse_celltype.csv")
    mean_per_detailed.to_csv(OUT / "mean_activity_per_detailed_celltype.csv")

    # z-score across celltypes so relative-specificity surfaces
    z_coarse = (mean_per_ct - mean_per_ct.mean(axis=0)) / (mean_per_ct.std(axis=0) + 1e-12)
    z_coarse.to_csv(OUT / "z_activity_per_coarse_celltype.csv")
    z_detailed = (mean_per_detailed - mean_per_detailed.mean(axis=0)) / (mean_per_detailed.std(axis=0) + 1e-12)

    print("\n[8/8] checking expected-TF benchmarks against lit:")
    print(f"     (z-score ≥ 1 means the TF's regulon is 1+ sigma above the mean celltype activity)\n")
    rows = []
    for tf, (expected_ct, ref) in EXPECTED.items():
        reg_name = f"{tf}_regulon"
        if reg_name not in z_coarse.columns:
            rows.append((tf, expected_ct, "(no regulon — TF not in GRN output)", None, None, ref))
            continue
        if expected_ct not in z_coarse.index:
            if expected_ct in z_detailed.index:
                zval = z_detailed.loc[expected_ct, reg_name]
                top_ct = z_detailed[reg_name].idxmax()
                top_z = z_detailed[reg_name].max()
            else:
                rows.append((tf, expected_ct, "(expected ct not in data)", None, None, ref))
                continue
            rows.append((tf, expected_ct, top_ct, zval, top_z, ref))
            continue
        zval = z_coarse.loc[expected_ct, reg_name]
        top_ct = z_coarse[reg_name].idxmax()
        rows.append((tf, expected_ct, top_ct, zval, z_coarse[reg_name].max(), ref))

    bench = pd.DataFrame(rows, columns=["TF", "expected_celltype", "observed_top_celltype",
                                         "z_in_expected", "z_of_top", "lit_reference"])
    bench.to_csv(OUT / "expected_tf_benchmark.csv", index=False)
    print(bench[["TF", "expected_celltype", "observed_top_celltype",
                 "z_in_expected", "z_of_top"]].to_string(index=False))

    # Summary
    hits = (bench["expected_celltype"] == bench["observed_top_celltype"]).sum()
    total = bench.dropna(subset=["z_in_expected"]).shape[0]
    summary = {
        "dataset": "Ziegler 2021 nasopharyngeal scRNA-seq (BCH/Broad/UMMC)",
        "n_cells_used": int(adata.n_obs),
        "n_genes_used": int(adata.n_vars),
        "n_donors": int(adata.obs["donor_id"].nunique()),
        "n_tfs_in_grn": len(tfs_in),
        "n_regulons": len(regulons),
        "grn_wall_s": round(t_grn, 1),
        "aucell_wall_s": round(t_auc, 1),
        "expected_tf_hit_rate": f"{hits}/{total}",
        "expected_tf_hit_pct": round(100 * hits / max(total, 1), 1),
    }
    with open(OUT / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    print(f"\n=== SUMMARY ===")
    print(json.dumps(summary, indent=2))
    try:
        out_display = OUT.relative_to(Path.cwd())
    except ValueError:
        out_display = OUT
    print(f"\nartefacts in {out_display}")


if __name__ == "__main__":
    main()
