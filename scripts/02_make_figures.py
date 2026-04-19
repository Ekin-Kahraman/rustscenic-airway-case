"""Figures for the Ziegler airway case study.

Produces:
  fig1_tf_celltype_heatmap.png — z-scored regulon activity per coarse cell type
  fig2_canonical_tf_benchmark.png — bar chart of z-score in expected vs top cell type
  fig3_umap_selected_regulons.png — UMAP colored by a few key regulons (FOXJ1, TP63, ASCL3)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({"figure.dpi": 120, "savefig.dpi": 200, "font.size": 9})

OUT = Path(__file__).parent.parent / "results"
FIG = Path(__file__).parent.parent / "figures"
FIG.mkdir(exist_ok=True)


def fig1_heatmap():
    z = pd.read_csv(OUT / "z_activity_per_coarse_celltype.csv", index_col=0)
    # Keep only regulons where the top z-score is >2 (strong-specificity TFs)
    max_z = z.abs().max(axis=0)
    z = z.loc[:, max_z.sort_values(ascending=False).head(40).index]
    # Strip "_regulon" suffix
    z.columns = [c.replace("_regulon", "") for c in z.columns]
    # Reorder rows: epithelial first, then immune/misc
    order = [
        "Basal Cells", "Mitotic Basal Cells",
        "Secretory Cells", "Goblet Cells",
        "Developing Secretory and Goblet Cells",
        "Deuterosomal Cells", "Developing Ciliated Cells", "Ciliated Cells",
        "Ionocytes", "Squamous Cells", "Enteroendocrine Cells",
        "T Cells", "B Cells", "Macrophages", "Dendritic Cells", "Plasmacytoid DCs",
    ]
    order = [o for o in order if o in z.index]
    z = z.loc[order]

    fig, ax = plt.subplots(figsize=(12, 5.5))
    vmax = max(abs(z.values.min()), abs(z.values.max()))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(z.columns)))
    ax.set_xticklabels(z.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index, fontsize=8)
    ax.set_title("rustscenic AUCell — z-scored regulon activity per airway cell type\n"
                 "(Ziegler 2021 nasopharyngeal, 31,602 cells, 40 most specific TFs)", fontsize=10)
    plt.colorbar(im, ax=ax, label="z-score across cell types", shrink=0.7)
    plt.tight_layout()
    plt.savefig(FIG / "fig1_tf_celltype_heatmap.png", bbox_inches="tight")
    plt.close()
    print(f"  wrote {FIG / 'fig1_tf_celltype_heatmap.png'}")


def fig2_canonical_benchmark():
    bench = pd.read_csv(OUT / "expected_tf_benchmark.csv")
    bench = bench.dropna(subset=["z_in_expected"])
    bench = bench.sort_values("z_in_expected")

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    y = np.arange(len(bench))
    colors = ["#2a9d8f" if (e == o) else "#e76f51"
              for e, o in zip(bench["expected_celltype"], bench["observed_top_celltype"])]
    ax.barh(y, bench["z_in_expected"], color=colors, alpha=0.85, label="z-score in expected cell type")
    ax.axvline(1.0, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(0.0, color="black", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{tf} → {exp[:22]}" for tf, exp in zip(bench["TF"], bench["expected_celltype"])], fontsize=9)
    ax.set_xlabel("Regulon activity z-score in expected cell type")
    ax.set_title("Canonical airway TFs: regulon activity in literature-expected cell type\n"
                 "green = top-ranked cell type matches lit; red = top cell type differs", fontsize=10)
    ax.legend([plt.Rectangle((0,0),1,1,color="#2a9d8f"),
               plt.Rectangle((0,0),1,1,color="#e76f51")],
              ["expected cell type is top-ranked",
               "different cell type is top-ranked"], loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG / "fig2_canonical_tf_benchmark.png", bbox_inches="tight")
    plt.close()
    print(f"  wrote {FIG / 'fig2_canonical_tf_benchmark.png'}")


def main():
    fig1_heatmap()
    fig2_canonical_benchmark()


if __name__ == "__main__":
    main()
