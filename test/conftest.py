"""Shared fixtures for the desidescsn test suite.

The suite exercises the pure / offline functions of the efficiency, surveys and
footprint modules against synthetic dataframes — no external catalogs, survey
files or network access are needed. ``catalogs`` (which needs ``GCRCatalogs``)
is only compile-checked / import-skipped.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

TEST_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = TEST_DIR.parent
for _p in (str(PACKAGE_ROOT), str(TEST_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)


@pytest.fixture
def host_catalog():
    """A small synthetic host-galaxy catalog with mass / SFR / redshift."""
    return pd.DataFrame(
        {
            "stellar_mass": [1e10, 2e10, 3e10, 4e10],
            "sfr": [1.0, 2.0, 0.0, 5.0],
            "redshift_true": [0.1, 0.3, 0.6, 0.9],
            "ra": [10.0, 10.1, 10.2, 10.3],
            "dec": [-1.0, -0.9, -0.8, -0.7],
        }
    )
