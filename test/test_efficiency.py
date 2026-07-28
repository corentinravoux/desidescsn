"""Unit tests for :mod:`desidescsn.efficiency`."""

import numpy as np
import pandas as pd
import pytest

from desidescsn import efficiency


# --- parametric n(z) -------------------------------------------------------


def test_n_z_function_positive_and_shape():
    z = np.linspace(0.05, 1.5, 50)
    n_z = efficiency.n_z_function(z, A=1.0, z0=0.3, beta=2.0, d=2.0)
    assert n_z.shape == z.shape
    assert np.all(n_z >= 0)
    # a scalar input gives a scalar output
    assert np.isscalar(efficiency.n_z_function(0.3, 1.0, 0.3, 2.0, 2.0))


def test_fit_n_z_recovers_parameters():
    z = np.linspace(0.05, 1.5, 80)
    true = [2.0, 0.35, 2.0, 2.5]
    y = efficiency.n_z_function(z, *true)
    popt = efficiency.fit_n_z(z, y)
    np.testing.assert_allclose(popt, true, rtol=1e-2)


# --- SN weight models ------------------------------------------------------


def test_return_weight_model_noweights(host_catalog):
    w = efficiency.return_weight_model(host_catalog, "noweigths", {})
    np.testing.assert_array_equal(w, np.ones(len(host_catalog)))


def test_return_weight_model_snia(host_catalog):
    w = efficiency.return_weight_model(host_catalog, "snia", {"A": 1.0, "B": 2.0})
    expected = 1.0 * host_catalog["stellar_mass"] + 2.0 * host_catalog["sfr"]
    np.testing.assert_allclose(w, expected)


def test_return_weight_model_sncc_sSFR_cut(host_catalog):
    # only galaxies above the sSFR cut get a non-zero weight
    w = efficiency.return_weight_model(
        host_catalog, "sncc", {"C": 0.36, "sSFR_cut": -11.5}
    )
    ssfr = np.log10(host_catalog["sfr"] / host_catalog["stellar_mass"])
    assert np.all(w[ssfr <= -11.5] == 0)
    assert np.all(w[ssfr > -11.5] > 0)


# --- SN efficiency ---------------------------------------------------------


def test_compute_sn_efficiency_range():
    sn = pd.DataFrame({"redshift_true": np.linspace(0.1, 0.9, 60)})
    mask = (sn["redshift_true"] < 0.5).values
    z, eff = efficiency.compute_sn_efficiency(sn, mask, N_z=7)
    assert z.shape == (6,)
    # efficiency is a fraction, dropping to ~0 at high redshift where the mask is False
    assert np.nanmin(eff) == pytest.approx(0.0, abs=1e-9)
    assert np.nanmax(eff) <= 1.0 + 1e-9


# --- host efficiency ratio -------------------------------------------------


def test_compute_host_efficiency_is_callable():
    z = np.linspace(0.05, 1.5, 60)
    n_survey = efficiency.n_z_function(z, 1.0, 0.3, 2.0, 2.0)
    n_simu = efficiency.n_z_function(z, 1.0, 0.4, 2.0, 2.0)
    eff_fn = efficiency.compute_host_efficiency(z, n_survey, z, n_simu)
    assert callable(eff_fn)
    val = eff_fn(0.3)
    assert np.isfinite(val)


def test_compute_host_efficiency_capped():
    z = np.linspace(0.05, 1.5, 60)
    n_survey = efficiency.n_z_function(z, 2.0, 0.3, 2.0, 2.0)
    n_simu = efficiency.n_z_function(z, 1.0, 0.3, 2.0, 2.0)
    eff_fn = efficiency.compute_host_efficiency(
        z, n_survey, z, n_simu, maximum_redshift_efficiency=1.0
    )
    # the cap bounds the ratio from above
    assert eff_fn(0.3) <= 1.0
