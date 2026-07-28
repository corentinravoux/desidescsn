"""Smoke tests: every shipped module byte-compiles and imports (deps allowing).

``catalogs`` imports ``GCRCatalogs`` (a heavy DESC dependency) at module load,
so its import is skipped when that package is unavailable while its syntax is
still checked.
"""

import importlib
import py_compile

import pytest

from conftest import PACKAGE_ROOT

CORE_MODULES = ["surveys", "efficiency", "footprint"]

ALL_PY = sorted((PACKAGE_ROOT / "desidescsn").glob("*.py"))


@pytest.mark.parametrize("path", ALL_PY, ids=[p.name for p in ALL_PY])
def test_module_compiles(path):
    """Every shipped .py file byte-compiles (no syntax errors)."""
    py_compile.compile(str(path), doraise=True)


@pytest.mark.parametrize("module", CORE_MODULES)
def test_core_module_imports(module):
    """Core modules import (they need only numpy/scipy/astropy/healpy/desimodel)."""
    importlib.import_module(f"desidescsn.{module}")


def test_catalogs_imports_or_skips():
    """``catalogs`` imports when GCRCatalogs is present, else is skipped."""
    try:
        importlib.import_module("desidescsn.catalogs")
    except Exception as exc:  # GCRCatalogs missing
        pytest.skip(f"GCRCatalogs unavailable: {exc!r}")
