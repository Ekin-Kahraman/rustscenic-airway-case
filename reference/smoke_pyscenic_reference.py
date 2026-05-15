"""Smoke-test the pinned RustScenic/pySCENIC AUCell reference environment."""

from __future__ import annotations

import warnings

warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")

import numpy as np
import pandas as pd
from ctxcore.genesig import Regulon
from pyscenic.aucell import aucell as pyscenic_aucell

import rustscenic.aucell


EXPR = pd.DataFrame(
    [
        [10, 5, 0, 0, 1],
        [8, 4, 0, 1, 0],
        [0, 1, 9, 7, 3],
        [1, 0, 8, 6, 2],
    ],
    index=["cell_a", "cell_b", "cell_c", "cell_d"],
    columns=["GATA3", "KRT5", "FOXJ1", "TP63", "ACTB"],
)

REGULONS = {
    "basal": ["GATA3", "KRT5"],
    "ciliated": ["FOXJ1", "TP63"],
}


def _pyscenic_regulons() -> list[Regulon]:
    return [
        Regulon(
            name="basal",
            gene2weight={"GATA3": 1.0, "KRT5": 1.0},
            transcription_factor="GATA3",
            gene2occurrence={},
        ),
        Regulon(
            name="ciliated",
            gene2weight={"FOXJ1": 1.0, "TP63": 1.0},
            transcription_factor="FOXJ1",
            gene2occurrence={},
        ),
    ]


def main() -> None:
    py_auc = pyscenic_aucell(
        EXPR,
        _pyscenic_regulons(),
        auc_threshold=0.4,
        num_workers=1,
        noweights=True,
        normalize=False,
    )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="top_frac=0.4 is unusually high")
        rust_auc = rustscenic.aucell.score(EXPR, REGULONS, top_frac=0.4)

    py_auc.index.name = None
    py_auc.columns.name = None
    pd.testing.assert_index_equal(rust_auc.index, py_auc.index)
    pd.testing.assert_index_equal(rust_auc.columns, py_auc.columns)
    np.testing.assert_allclose(rust_auc.to_numpy(), py_auc.to_numpy(), rtol=1e-12, atol=1e-12)

    assert rust_auc.loc["cell_a", "basal"] > rust_auc.loc["cell_a", "ciliated"]
    assert rust_auc.loc["cell_c", "ciliated"] > rust_auc.loc["cell_c", "basal"]
    print("pySCENIC/RustScenic AUCell reference smoke passed")


if __name__ == "__main__":
    main()
