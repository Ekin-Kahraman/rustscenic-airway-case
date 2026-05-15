"""Validate committed case-study outputs without requiring the raw atlas."""

from __future__ import annotations

import csv
import json
import struct
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
COMPARISON = RESULTS / "comparison"
FIGURES = ROOT / "figures"
FIGURE_NAMES = [
    "tf_celltype_heatmap",
    "canonical_tf_benchmark",
    "canonical_tf_3way",
    "per_cell_pearson_hist",
    "runtime_comparison",
    "covid_differential",
]


def require(path: Path, min_bytes: int = 1) -> None:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.stat().st_size < min_bytes:
        raise AssertionError(f"{path} is unexpectedly small")


def load_json(path: Path) -> dict:
    require(path)
    return json.loads(path.read_text(encoding="utf-8"))


def count_csv_rows(path: Path) -> tuple[int, list[str]]:
    require(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = sum(1 for _ in reader)
        return rows, list(reader.fieldnames or [])


def validate_png(path: Path) -> None:
    require(path, min_bytes=1_000)
    with path.open("rb") as handle:
        signature = handle.read(8)
    if signature != b"\x89PNG\r\n\x1a\n":
        raise AssertionError(f"{path} is not a PNG")
    with path.open("rb") as handle:
        handle.read(16)
        width, height = struct.unpack(">II", handle.read(8))
    if width < 800 or height < 500:
        raise AssertionError(f"{path} is too small for publication review: {width}x{height}")


def validate_run_log(summary: dict) -> None:
    run_log = RESULTS / "run.log"
    require(run_log, min_bytes=100)
    text = run_log.read_text(encoding="utf-8")
    required_snippets = [
        f'"n_cells_used": {summary["n_cells_used"]}',
        f'"n_regulons": {summary["n_regulons"]}',
        f'"grn_wall_s": {summary["grn_wall_s"]}',
        f'"expected_tf_hit_rate": "{summary["expected_tf_hit_rate"]}"',
    ]
    for snippet in required_snippets:
        if snippet not in text:
            raise AssertionError(f"{run_log} does not match summary.json: missing {snippet}")


def main() -> None:
    summary = load_json(RESULTS / "summary.json")
    assert summary["n_cells_used"] == 31602
    assert summary["n_donors"] == 58
    assert summary["n_regulons"] >= 50
    assert summary["expected_tf_hit_rate"].startswith("8/")
    validate_run_log(summary)

    comparison = load_json(COMPARISON / "summary.json")
    assert comparison["n_cells_used"] == 31602
    assert comparison["n_regulons"] >= 50
    assert comparison["per_cell_pearson"]["rust_vs_py_unit"]["mean"] > 0.95

    required_results = [
        RESULTS / "regulons.json",
        RESULTS / "z_activity_per_coarse_celltype.csv",
        RESULTS / "mean_activity_per_coarse_celltype.csv",
        RESULTS / "mean_activity_per_detailed_celltype.csv",
        RESULTS / "comparison_run.log",
        RESULTS / "run.log",
        COMPARISON / "per_cell_pearson_summary.csv",
        COMPARISON / "per_regulon_pearson_summary.csv",
    ]
    for path in required_results:
        require(path, min_bytes=100)

    covid_rows, covid_fields = count_csv_rows(RESULTS / "covid_differential_regulons.csv")
    assert covid_rows > 100
    assert {"celltype", "regulon", "log2_fc", "qvalue"}.issubset(covid_fields)

    bench_rows, bench_fields = count_csv_rows(RESULTS / "expected_tf_benchmark.csv")
    assert bench_rows >= 14
    assert {"TF", "expected_celltype", "observed_top_celltype"}.issubset(bench_fields)

    canonical_rows, canonical_fields = count_csv_rows(COMPARISON / "canonical_tf_compare.csv")
    assert canonical_rows >= 40
    assert {"tf", "expected", "method", "z_in_expected", "hit"}.issubset(canonical_fields)

    for i, name in enumerate(FIGURE_NAMES, start=1):
        validate_png(FIGURES / f"fig{i}_{name}.png")


if __name__ == "__main__":
    main()
