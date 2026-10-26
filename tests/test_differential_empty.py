"""Exercise empty differential results using synthetic cells only."""

from contextlib import redirect_stdout
import importlib.util
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
RESULT_COLUMNS = [
    "celltype", "regulon", "mean_pos", "mean_neg", "log2_fc",
    "ranksum_stat", "pvalue", "n_pos", "n_neg", "qvalue",
]


class EmptyDifferentialTests(unittest.TestCase):
    def check_empty_result(
        self, n_per_arm: int, constant: bool, statuses: tuple[str, str] = ("pos", "neg")
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cells = [f"cell_{i}" for i in range(n_per_arm * 2)]
            obs = pd.DataFrame(
                {
                    "donor_id": ["synthetic_donor"] * len(cells),
                    "SARSCoV2_PCR_Status": [statuses[0]] * n_per_arm + [statuses[1]] * n_per_arm,
                    "Coarse_Cell_Annotations": ["synthetic_type"] * len(cells),
                },
                index=cells,
            )
            dataset = root / "synthetic.h5ad"
            ad.AnnData(np.zeros((len(cells), 1)), obs=obs).write_h5ad(dataset)
            results = root / "results"
            results.mkdir()
            activity = np.zeros(len(cells)) if constant else np.linspace(0, 1, len(cells))
            pd.DataFrame({"TF_regulon": activity}, index=cells).to_parquet(results / "auc.parquet")

            spec = importlib.util.spec_from_file_location(
                "differential_analysis", SCRIPTS / "05_covid_differential_regulons.py"
            )
            module = importlib.util.module_from_spec(spec)
            with patch.dict(os.environ, {"ZIEGLER_H5AD": str(dataset)}), patch.object(
                sys, "path", [str(SCRIPTS), *sys.path]
            ):
                spec.loader.exec_module(module)
            module.OUT = results
            module.FIG = root / "figures"
            output = io.StringIO()
            with redirect_stdout(output):
                module.main()

            table = pd.read_csv(results / "covid_differential_regulons.csv")
            self.assertTrue(table.empty)
            self.assertEqual(table.columns.tolist(), RESULT_COLUMNS)
            self.assertIn("no significant regulons to plot", output.getvalue())
            self.assertFalse(module.FIG.exists())

    def test_no_cell_type_meets_minimum_group_size(self) -> None:
        self.check_empty_result(n_per_arm=60, constant=False)

    def test_all_regulons_have_constant_activity(self) -> None:
        self.check_empty_result(n_per_arm=100, constant=True)

    def test_neither_expected_covid_status_is_present(self) -> None:
        self.check_empty_result(n_per_arm=100, constant=False, statuses=("unknown", "unknown"))

    def test_one_expected_covid_status_is_missing(self) -> None:
        for status in ("pos", "neg"):
            with self.subTest(status=status):
                self.check_empty_result(n_per_arm=100, constant=False, statuses=(status, "unknown"))


if __name__ == "__main__":
    unittest.main()
