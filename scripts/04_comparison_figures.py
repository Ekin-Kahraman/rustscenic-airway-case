"""Figures for the head-to-head comparison.

Generates:
  fig3_canonical_tf_3way.png — 14 TFs × 3 methods (rustscenic / pyscenic-unit /
    pyscenic-weighted) side-by-side bar chart of z in expected cell type
  fig4_per_cell_pearson_hist.png — distribution of per-cell Pearson (rust vs py)
  fig5_runtime_comparison.png — wall-time bar chart
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({"figure.dpi": 120, "savefig.dpi": 200, "font.size": 9})

OUT = Path(__file__).parent.parent / "results" / "comparison"
FIG = Path(__file__).parent.parent / "figures"


def fig3_canonical_3way():
    z = pd.read_csv(OUT / "canonical_tf_compare.csv")
    # Keep the ordering from EXPECTED
    order = ["ASCL3","FOXJ1","CEBPB","TBX21","EOMES","SPI1","TP63","STAT1","RUNX3",
             "SPDEF","MYB","IRF7","SOX2","PAX5"]
    pivot_z = z.pivot(index="tf", columns="method", values="z_in_expected").loc[order]
    pivot_hit = z.pivot(index="tf", columns="method", values="hit").loc[order]

    fig, ax = plt.subplots(figsize=(10, 6))
    methods = ["rustscenic", "pyscenic_unit", "pyscenic_weighted"]
    colors = {"rustscenic": "#2a9d8f", "pyscenic_unit": "#e9c46a", "pyscenic_weighted": "#e76f51"}
    n = len(order)
    width = 0.27
    xs = np.arange(n)
    for i, m in enumerate(methods):
        vals = pivot_z[m].values
        hits = pivot_hit[m].values
        # Edge highlight for hits
        bars = ax.bar(xs + (i - 1) * width, vals, width=width * 0.9,
                      color=colors[m], alpha=0.85, label=m.replace("_", " "))
        for bar, hit in zip(bars, hits):
            if hit:
                bar.set_edgecolor("black")
                bar.set_linewidth(1.3)
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axhline(1.0, color="grey", linestyle=":", linewidth=0.8)
    ax.set_xticks(xs); ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("z-score of regulon activity in expected cell type")
    ax.set_title("Canonical airway TFs: head-to-head on Ziegler 2021 (n=31,602 cells)\n"
                 "black outline = cell type with highest regulon activity matches lit-expected", fontsize=10)
    ax.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG / "fig3_canonical_tf_3way.png", bbox_inches="tight")
    plt.close()
    print(f"  wrote {FIG / 'fig3_canonical_tf_3way.png'}")


def fig4_per_cell_hist():
    s = json.load(open(OUT / "summary.json"))
    pct = s["per_cell_pearson"]

    # We need the actual distributions — re-derive from the parquet outputs
    auc_rust = pd.read_parquet(OUT.parent / "auc.parquet")
    auc_py_u = pd.read_parquet(OUT / "auc_pyscenic_unit.parquet")
    auc_py_w = pd.read_parquet(OUT / "auc_pyscenic_weighted.parquet")
    common_regs = sorted(set(auc_rust.columns) & set(auc_py_u.columns) & set(auc_py_w.columns))
    common_cells = sorted(set(auc_rust.index) & set(auc_py_u.index) & set(auc_py_w.index))
    R = auc_rust.loc[common_cells, common_regs].values
    PU = auc_py_u.loc[common_cells, common_regs].values
    PW = auc_py_w.loc[common_cells, common_regs].values

    def per_cell(A, B):
        out = []
        for i in range(A.shape[0]):
            a, b = A[i], B[i]
            if np.std(a) > 1e-12 and np.std(b) > 1e-12:
                out.append(np.corrcoef(a, b)[0, 1])
        return np.asarray(out)
    rho_u = per_cell(R, PU)
    rho_w = per_cell(R, PW)

    fig, ax = plt.subplots(figsize=(8, 4.2))
    bins = np.linspace(0.4, 1.0, 50)
    ax.hist(rho_u, bins=bins, alpha=0.7, color="#2a9d8f",
            label=f"rust vs pyscenic-unit (mean {rho_u.mean():.3f}, >0.95 {100*(rho_u>0.95).mean():.1f}%)",
            edgecolor="black", linewidth=0.3)
    ax.hist(rho_w, bins=bins, alpha=0.6, color="#e76f51",
            label=f"rust vs pyscenic-weighted (mean {rho_w.mean():.3f}, >0.95 {100*(rho_w>0.95).mean():.1f}%)",
            edgecolor="black", linewidth=0.3)
    ax.axvline(0.95, linestyle=":", color="grey")
    ax.set_xlabel("Per-cell Pearson of regulon-activity profile (rustscenic vs pyscenic)")
    ax.set_ylabel(f"number of cells (n = {len(rho_u):,})")
    ax.set_title("Agreement between rustscenic and pyscenic on Ziegler 2021 nasopharyngeal\n"
                 "(31,602 cells, 59 regulons, identical GRN adjacencies on both sides)", fontsize=10)
    ax.legend(loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG / "fig4_per_cell_pearson_hist.png", bbox_inches="tight")
    plt.close()
    print(f"  wrote {FIG / 'fig4_per_cell_pearson_hist.png'}")


def fig5_runtime():
    s = json.load(open(OUT / "summary.json"))
    rt = s["runtime_memory"]
    methods = ["rustscenic", "pyscenic_unit", "pyscenic_weighted"]
    walls = [rt[m]["wall_s"] for m in methods]
    speedups = [rt["pyscenic_unit"]["wall_s"] / rt["rustscenic"]["wall_s"],
                rt["pyscenic_weighted"]["wall_s"] / rt["rustscenic"]["wall_s"]]

    fig, ax = plt.subplots(figsize=(6.5, 4))
    colors = ["#2a9d8f", "#e9c46a", "#e76f51"]
    bars = ax.bar([m.replace("_", " ") for m in methods], walls, color=colors, alpha=0.9,
                  edgecolor="black", linewidth=0.5)
    for bar, wall in zip(bars, walls):
        ax.text(bar.get_x() + bar.get_width()/2, wall + max(walls)*0.02,
                f"{wall:.2f}s", ha="center", fontsize=10, fontweight="bold")
    ax.set_ylabel("wall-clock (s)")
    ax.set_title(f"AUCell runtime on Ziegler 2021 — 31,602 cells × 59 regulons\n"
                 f"rustscenic is {speedups[0]:.0f}× faster than pyscenic-unit, "
                 f"{speedups[1]:.0f}× faster than pyscenic-weighted", fontsize=10)
    plt.tight_layout()
    plt.savefig(FIG / "fig5_runtime_comparison.png", bbox_inches="tight")
    plt.close()
    print(f"  wrote {FIG / 'fig5_runtime_comparison.png'}")


def main():
    fig3_canonical_3way()
    fig4_per_cell_hist()
    fig5_runtime()


if __name__ == "__main__":
    main()
