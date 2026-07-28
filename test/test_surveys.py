"""Unit tests for :mod:`desidescsn.surveys` (selection cuts + efficiencies)."""

import numpy as np
import pandas as pd
import pytest

from desidescsn import surveys


# --- constant host efficiencies --------------------------------------------


def test_simple_host_efficiencies_are_products():
    assert surveys.host_efficiency_DESIBGS_simple() == pytest.approx(
        surveys.fiber_assignement_efficiency_DESIBGS
        * surveys.redshift_success_rate_DESIBGS
    )
    assert surveys.host_efficiency_DESILRG_simple() == pytest.approx(0.989)
    for eff in (
        surveys.host_efficiency_DESI2_simple(),
        surveys.host_efficiency_CRSBG_simple(),
        surveys.host_efficiency_CRSLRG_simple(),
        surveys.host_efficiency_4HS_simple(),
    ):
        assert 0.0 < eff <= 1.0


# --- magnitude selections --------------------------------------------------


def test_mask_magnitude_desibgs():
    galaxies = pd.DataFrame({"mag_r": [19.0, 20.5, 21.0, 20.175]})
    mask = surveys.mask_magnitude_DESIBGS(galaxies, {"r": "mag_r"})
    # strictly below 20.175 passes
    np.testing.assert_array_equal(mask.values, [True, False, False, False])


def test_mask_magnitude_4hs_color_cut():
    # J < 18 and (J - K) < 0.45
    galaxies = pd.DataFrame({"J": [17.0, 17.0, 19.0], "K": [17.0, 16.4, 15.0]})
    hashing = {"J": "J", "K": "K"}
    # row0: J<18 & J-K=0<0.45 pass; row1: J-K=0.6 fails colour; row2: J>18 fails mag
    mask = surveys.mask_magnitude_4HS(galaxies, hashing, cut_color=True)
    np.testing.assert_array_equal(mask.values, [True, False, False])
    # without the colour cut, only the magnitude limit applies
    mask_nocolor = surveys.mask_magnitude_4HS(galaxies, hashing, cut_color=False)
    np.testing.assert_array_equal(mask_nocolor.values, [True, True, False])


def test_mask_magnitude_crslrg():
    galaxies = pd.DataFrame(
        {"J": [18.5, 17.0, 20.0], "r": [19.0, 19.0, 19.0], "W1": [17.0, 17.0, 17.0]}
    )
    hashing = {"J": "J", "r": "r", "W1": "W1"}
    mask = surveys.mask_magnitude_CRSLRG(galaxies, hashing, cut_color=False)
    # 18 < J < 19.5
    np.testing.assert_array_equal(mask.values, [True, False, False])
