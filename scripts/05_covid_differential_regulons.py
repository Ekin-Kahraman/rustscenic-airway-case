"""COVID± differential regulon analysis on Ziegler 2021.

Using rustscenic regulon activities, within each cell type, test which regulons
are differentially active between COVID+ and COVID- cells. This extends the
covid-airway-deconvolution story from "which cells are perturbed" to "which
regulatory programmes are rewired during infection".

Uses Wilcoxon rank-sum within each cell type + donor-stratified correction:
  - For each cell type:
    - For each regulon:
      - Collect per-cell activity, split by SARSCoV2_PCR_Status
      - Wilcoxon rank-sum test (non-parametric)
      - Fold change of means (COVID+ / COVID-)
  - BH-FDR correct per cell type

Outputs:
  covid_differential_regulons_per_celltype.csv — full table
  figures/fig6_covid_differential.png — top hits heatmap
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

OUT = Path(__file__).parent.parent / "results"
FIG = Path(__file__).parent.parent / "figures"

DATA = Path("/Users/ekin/covid-airway-deconvolution/data/ziegler2021_nasopharyngeal.h5ad")


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(n)
    q = p * n / (ranks + 1)
    # Monotonically non-decreasing along sorted order
    sorted_q = q[order]
    for i in range(n - 2, -1, -1):
        if sorted_q[i] > sorted_q[i + 1]:
            sorted_q[i] = sorted_q[i + 1]
    out = np.empty(n)
    out[order] = np.minimum(sorted_q, 1.0)
    return out


def main():
    print("loading AUCell + Ziegler metadata...")
    auc = pd.read_parquet(OUT / "auc.parquet")

    # Fetch donor + COVID status + celltype from the original h5ad
    a = ad.read_h5ad(DATA, backed="r")
    shared = [c for c in auc.index if c in a.obs.index]
    auc = auc.loc[shared]
    meta = pd.DataFrame({
        "donor": a.obs.loc[shared, "donor_id"].values,
        "covid": a.obs.loc[shared, "SARSCoV2_PCR_Status"].astype(str).values,
        "celltype": a.obs.loc[shared, "Coarse_Cell_Annotations"].astype(str).values,
    }, index=shared)
    print(f"cells: {len(auc)}, regulons: {auc.shape[1]}")
    print(f"covid split: {meta['covid'].value_counts().to_dict()}")

    # Filter out cell types with <100 cells in either arm
    ct_counts = meta.groupby(["celltype", "covid"], observed=False).size().unstack(fill_value=0)
    valid_cts = ct_counts.index[(ct_counts.get("pos", 0) >= 100) & (ct_counts.get("neg", 0) >= 100)]
    print(f"cell types with ≥100 cells in both arms: {len(valid_cts)}")
    print(ct_counts.loc[valid_cts].to_string())

    rows = []
    for ct in valid_cts:
        sub = meta[meta["celltype"] == ct]
        pos_idx = sub.index[sub["covid"] == "pos"]
        neg_idx = sub.index[sub["covid"] == "neg"]
        if len(pos_idx) < 100 or len(neg_idx) < 100:
            continue
        A_pos = auc.loc[pos_idx]
        A_neg = auc.loc[neg_idx]
        for reg in auc.columns:
            p_vals = A_pos[reg].values
            n_vals = A_neg[reg].values
            if p_vals.std() < 1e-12 and n_vals.std() < 1e-12:
                continue
            stat, pval = stats.ranksums(p_vals, n_vals)
            fc = (p_vals.mean() + 1e-12) / (n_vals.mean() + 1e-12)
            rows.append({
                "celltype": ct, "regulon": reg,
                "mean_pos": p_vals.mean(), "mean_neg": n_vals.mean(),
                "log2_fc": np.log2(fc) if fc > 0 else 0.0,
                "ranksum_stat": stat, "pvalue": pval,
                "n_pos": len(p_vals), "n_neg": len(n_vals),
            })
    df = pd.DataFrame(rows)
    # BH-FDR per cell type
    df["qvalue"] = df.groupby("celltype")["pvalue"].transform(bh_fdr)
    df = df.sort_values(["celltype", "qvalue"])
    df.to_csv(OUT / "covid_differential_regulons.csv", index=False)
    print(f"\nsaved {OUT / 'covid_differential_regulons.csv'}")

    # Top-10 per celltype with q<0.01
    top_each = df[df["qvalue"] < 0.01].groupby("celltype", observed=False).head(10)
    print(f"\ntotal significant (q<0.01): {(df['qvalue'] < 0.01).sum()}")
    print(f"\ntop COVID± differential regulons per cell type (q<0.01, |log2 FC| sorted):")
    for ct in top_each["celltype"].unique():
        sub = top_each[top_each["celltype"] == ct].sort_values("log2_fc", key=abs, ascending=False).head(5)
        if sub.empty: continue
        print(f"\n  === {ct} ===")
        for _, r in sub.iterrows():
            direction = "↑COVID+" if r["log2_fc"] > 0 else "↓COVID+"
            print(f"    {r['regulon']:25s}  log2FC={r['log2_fc']:+.2f}  q={r['qvalue']:.1e}  {direction}")

    # Figure: heatmap of log2FC for top differential regulons × celltype
    import matplotlib.pyplot as plt
    # Pick regulons significant in at least one cell type
    sig_regs = df[df["qvalue"] < 0.01]["regulon"].value_counts().head(20).index.tolist()
    if not sig_regs:
        print("\nno significant regulons to plot")
        return
    heat = df[df["regulon"].isin(sig_regs)].pivot(
        index="regulon", columns="celltype", values="log2_fc").fillna(0)
    heat = heat.loc[sig_regs]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    vmax = float(np.abs(heat.values).max())
    im = ax.imshow(heat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels([r.replace("_regulon","") for r in heat.index], fontsize=8)
    ax.set_title("COVID+ vs COVID- differential regulon activity per airway cell type\n"
                 "(log2 fold change; Ziegler 2021, top 20 differential regulons, q<0.01)", fontsize=10)
    plt.colorbar(im, ax=ax, label="log2 FC (COVID+ / COVID-)", shrink=0.7)
    plt.tight_layout()
    plt.savefig(FIG / "fig6_covid_differential.png", bbox_inches="tight")
    print(f"\nwrote {FIG / 'fig6_covid_differential.png'}")


if __name__ == "__main__":
    main()
