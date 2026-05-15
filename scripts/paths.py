"""Shared local paths for the case-study scripts."""

from __future__ import annotations

import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def ziegler_h5ad_path() -> Path:
    """Return the Ziegler h5ad path from env or sibling repo default."""
    env_path = os.environ.get("ZIEGLER_H5AD")
    path = (
        Path(env_path).expanduser()
        if env_path
        else ROOT.parent / "covid-airway-deconvolution" / "data" / "ziegler2021_nasopharyngeal.h5ad"
    )
    if not path.exists():
        raise FileNotFoundError(
            f"Ziegler h5ad not found at {path}. Set ZIEGLER_H5AD=/path/to/ziegler2021_nasopharyngeal.h5ad"
        )
    return path
