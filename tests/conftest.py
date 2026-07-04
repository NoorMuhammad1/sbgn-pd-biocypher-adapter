"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


@pytest.fixture()
def upper_glycolysis_path() -> Path:
    """Path to the upper-glycolysis sample SBGN-ML file."""
    return DATA_DIR / "glycolysis_upper.sbgn"


@pytest.fixture()
def lower_glycolysis_path() -> Path:
    """Path to the lower-glycolysis sample SBGN-ML file."""
    return DATA_DIR / "glycolysis_lower.sbgn"


@pytest.fixture()
def both_glycolysis_paths(upper_glycolysis_path: Path, lower_glycolysis_path: Path) -> list[Path]:
    """Both sample files together (upper first)."""
    return [upper_glycolysis_path, lower_glycolysis_path]
