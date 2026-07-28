"""Per-survey target selection, host efficiency and ``n(z)`` for SN forecasts.

Encodes the approximate spectroscopic target-selection cuts, fibre-assignment x
redshift-success efficiencies and observed redshift distributions of the surveys
considered for DESI/DESC SN host follow-up: DESI (BGS, LRG), DESI2, 4MOST CRS
(BG, LRG) and 4HS. Module-level constants hold each survey's magnitude / colour
thresholds; the ``mask_magnitude_*`` functions apply them to a galaxy sample,
``host_efficiency_*_simple`` give the constant host efficiency, and the
``n_z_*`` functions load the survey ``n(z)``. Band names are mapped to catalog
columns through a ``hashing_table`` dict.
"""

import fitsio
import numpy as np
import pandas as pd
import yaml

from desidescsn import efficiency

### DESI

# BGS Hypothesis:
# - No faint rfiber cut to remove faint galaxies with low fiber fluxes (++)
# - No differenciation between BGS Bright and Faint (-)
# - No very bright limit cut (---)
# - No extreme colot cut (---)
# - No bright rfiber cut to remove ‘shredded’ galaxies (---)
# - No star-galaxy separation (---)


magnitude_cut_DESIBGS = 20.175
band_DESIBGS = "r"
fiber_assignement_efficiency_DESIBGS = 0.85
redshift_success_rate_DESIBGS = 0.95

# LRG Hypothesis:
# - zfiber = z (+++)
# - W1 (3400 nm) = ROMAN_obs_K213 (2125 nm) (+)
# - No star-galaxy separation (---)


magnitude_cut_DESILRG = 21.60
color_cut_DESILRG_1 = -0.6
color_cut_DESILRG_2 = 2.9
color_cut_DESILRG_3 = 1.8
color_cut_DESILRG_4 = -30.852
color_cut_DESILRG_5 = -16.33
color_cut_DESILRG_6 = 3.3
band_DESILRG_1 = "zfiber"
band_DESILRG_2 = "z"
band_DESILRG_3 = "W1"
band_DESILRG_4 = "r"
band_DESILRG_5 = "g"

fiber_assignement_efficiency_DESILRG = 1.0
redshift_success_rate_DESILRG = 0.989


### DESI2

# Hypothesis:
# - zfiber = z (+++)
# - No stellar rejection (+)
# - No lower rfibertot and zfibertot cuts (--)
# - No star-galaxy separation (---)


color_cut_DESI2_1 = -0.6
color_cut_DESI2_2 = 2.9
color_cut_DESI2_3 = 1.8
band_DESI2_1 = "zfiber"
band_DESI2_2 = "z"
band_DESI2_3 = "W1"
band_DESI2_4 = "r"
band_DESI2_5 = "g"
fiber_assignement_efficiency_DESI2 = 0.8
redshift_success_rate_DESI2 = 0.95


### 4MOST CRS

# BG Hypothesis:
# - W1 (3400 nm) = ROMAN_obs_K213 (2125 nm) (+)
# - VHS = ROMAN (--)

magnitude_cut_CRSBG_1 = 16
magnitude_cut_CRSBG_2 = 18
color_cut_CRSBG_1 = -1.6
color_cut_CRSBG_2 = -0.5
color_cut_CRSBG_3 = 0.1
color_cut_CRSBG_4 = 0.1

band_CRSBG_1 = "J"
band_CRSBG_2 = "K"
band_CRSBG_3 = "W1"

fiber_assignement_efficiency_CRSBG = 1.0
redshift_success_rate_CRSBG = 0.95

# LRG Hypothesis:
# - Only one color cut (+++)
# - W1 (3400 nm) = ROMAN_obs_K213 (2125 nm) (+)
# - VHS = ROMAN (--)

magnitude_cut_CRSLRG_1 = 18
magnitude_cut_CRSLRG_2 = 19.5
color_cut_CRSLRG_1 = -7.7
band_CRSLRG_1 = "J"
band_CRSLRG_2 = "r"
band_CRSLRG_3 = "W1"

fiber_assignement_efficiency_CRSLRG = 1.0
redshift_success_rate_CRSLRG = 0.75


### 4MOST 4HS

# Hypothesis:
# - No PV sub-sample included (+)
# - No 2MASS extended source catalog (+)
# - VHS = ROMAN (--)

magnitude_cut_4HS = 18
color_cut_4HS = 0.45
band_4HS_1 = "J"
band_4HS_2 = "K"

fiber_assignement_efficiency_4HS = 1.0
redshift_success_rate_4HS = 0.98


# DESI functions


def mask_magnitude_DESIBGS(
    file_sn,
    hashing_table,
    cut_color=True,
):
    """Apply the DESI BGS magnitude selection to a galaxy sample.

    Args:
        file_sn (pandas.DataFrame): Galaxy sample.
        hashing_table (dict): Band-name -> column-name mapping.
        cut_color (bool, optional): Unused for BGS (magnitude cut only).

    Returns:
        numpy.ndarray: Boolean mask of galaxies passing the r-band cut.
    """
    key_r = hashing_table[band_DESIBGS]
    mask_magnitude = file_sn[key_r] < magnitude_cut_DESIBGS

    return mask_magnitude


def host_efficiency_DESIBGS_simple():
    """DESI BGS constant host efficiency (fibre assignment x redshift success).

    Returns:
        float: The efficiency product.
    """
    return fiber_assignement_efficiency_DESIBGS * redshift_success_rate_DESIBGS


def mask_magnitude_DESILRG(
    file_sn,
    hashing_table,
    cut_color=True,
):
    """Apply the DESI LRG magnitude and colour selection to a galaxy sample.

    Args:
        file_sn (pandas.DataFrame): Galaxy sample.
        hashing_table (dict): Band-name -> column-name mapping (uses
            ``zfiber``/``z``/``W1``/``r``/``g``).
        cut_color (bool, optional): Also apply the LRG colour cuts.

    Returns:
        numpy.ndarray: Boolean mask of galaxies passing the selection.
    """
    key_zfiber = hashing_table[band_DESILRG_1]
    key_z = hashing_table[band_DESILRG_2]
    key_W1 = hashing_table[band_DESILRG_3]
    key_r = hashing_table[band_DESILRG_4]
    key_g = hashing_table[band_DESILRG_5]

    color_to_cut_1 = (
        file_sn[key_z] - file_sn[key_W1] - 0.8 * (file_sn[key_r] - file_sn[key_z])
    )
    color_to_cut_2 = file_sn[key_g] - file_sn[key_W1]
    color_to_cut_3 = file_sn[key_r] - file_sn[key_W1]
    color_to_cut_4 = file_sn[key_r] - file_sn[key_W1] - 1.8 * file_sn[key_W1]
    color_to_cut_5 = file_sn[key_r] - 2 * file_sn[key_W1]
    color_to_cut_6 = file_sn[key_r] - file_sn[key_W1]

    mask_magnitude = file_sn[key_zfiber] < magnitude_cut_DESILRG  # 1a
    if cut_color:
        mask_magnitude &= color_to_cut_1 > color_cut_DESILRG_1  # 1b
        mask_magnitude &= (color_to_cut_2 > color_cut_DESILRG_2) | (
            color_to_cut_3 > color_cut_DESILRG_3
        )  # 1c
        mask_magnitude &= (
            (color_to_cut_4 > color_cut_DESILRG_4)
            & (color_to_cut_5 > color_cut_DESILRG_5)
        ) | (
            color_to_cut_6 > color_cut_DESILRG_6
        )  # 1d

    return mask_magnitude


def host_efficiency_DESILRG_simple():
    """DESI LRG constant host efficiency (fibre assignment x redshift success).

    Returns:
        float: The efficiency product.
    """
    return fiber_assignement_efficiency_DESILRG * redshift_success_rate_DESILRG


def n_z_DESI_from_target_selection(
    targets,
    n_z_file,
    redshift_efficiency_file,
):
    """Build the DESI ``n(z)`` from a target-selection density file.

    Weights the normalised redshift histogram by the good-redshift number
    density (observed targets x per-tracer success rate from the yaml file).

    Args:
        targets (str): ``"BGS"`` or ``"LRG"``.
        n_z_file (str): Text file with redshift bin edges and counts.
        redshift_efficiency_file (str): Yaml with per-tracer nobs/success rates.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: ``(z_centers, n_z)`` (per deg^2 per
        unit redshift).
    """
    with open(redshift_efficiency_file, "r") as file:
        redshift_efficiency = yaml.safe_load(file)
    if targets == "BGS":
        n_good_z = (
            redshift_efficiency["nobs_bgs_faint"]
            * redshift_efficiency["success_bgs_faint"]
            + redshift_efficiency["nobs_bgs_bright"]
            * redshift_efficiency["success_bgs_bright"]
        )
    elif targets == "LRG":
        n_good_z = redshift_efficiency["nobs_lrg"] * redshift_efficiency["success_lrg"]

    n_redshift = np.loadtxt(n_z_file, skiprows=3)
    n_z_normalized = n_redshift[:, 2] / n_redshift.sum()
    z_centers = (n_redshift[:, 0] + n_redshift[:, 1]) / 2
    delta_z = np.mean(z_centers[1:] - z_centers[:-1])
    n_z = n_z_normalized * n_good_z / delta_z
    return z_centers, n_z


def n_z_DESI_from_y1(
    targets,
    n_z_file,
    redshift_efficiency_file,
    zmin=0.0,
    zmax=4.0,
    bins=200,
):
    """Build the DESI ``n(z)`` from a DESI Y1 redshift catalog (private data).

    Like :func:`n_z_DESI_from_target_selection` but histograms the observed
    ``Z_not4clus`` redshifts from a Y1 FITS clustering catalog.

    Args:
        targets (str): ``"BGS"``, ``"LRG"``, ``"ELG"`` or ``"QSO"``.
        n_z_file (str): Y1 clustering FITS catalog.
        redshift_efficiency_file (str): Yaml with per-tracer nobs/success rates.
        zmin (float, optional): Histogram lower edge.
        zmax (float, optional): Histogram upper edge.
        bins (int, optional): Number of redshift bins.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: ``(z_centers, n_z)``.
    """
    with open(redshift_efficiency_file, "r") as file:
        redshift_efficiency = yaml.safe_load(file)
    if targets == "BGS":
        n_good_z = (
            redshift_efficiency["nobs_bgs_faint"]
            * redshift_efficiency["success_bgs_faint"]
            + redshift_efficiency["nobs_bgs_bright"]
            * redshift_efficiency["success_bgs_bright"]
        )
    elif targets == "LRG":
        n_good_z = redshift_efficiency["nobs_lrg"] * redshift_efficiency["success_lrg"]
    elif targets == "ELG":
        n_good_z = redshift_efficiency["nobs_elg"] * redshift_efficiency["success_elg"]
    elif targets == "QSO":
        n_good_z = redshift_efficiency["nobs_qso"] * redshift_efficiency["success_qso"]

    file = fitsio.FITS(n_z_file)[1]
    file = file[["Z_not4clus"]]
    data = pd.DataFrame(file.read().byteswap().newbyteorder())
    hist, edges = np.histogram(data["Z_not4clus"], bins=bins, range=(zmin, zmax))
    n_z_normalized = hist / np.sum(hist)
    z_centers = (edges[:-1] + edges[1:]) / 2
    delta_z = np.mean(z_centers[1:] - z_centers[:-1])
    n_z = n_z_normalized * n_good_z / delta_z
    return z_centers, n_z


# DESI2 functions


def mask_magnitude_DESI2(
    file_sn,
    hashing_table,
    strategy_file=None,
    strategy_index=None,
    cut_color=True,
):
    """Apply the DESI2 magnitude and colour selection to a galaxy sample.

    The z-fibre magnitude threshold is read from a survey-strategy FITS file.

    Args:
        file_sn (pandas.DataFrame): Galaxy sample.
        hashing_table (dict): Band-name -> column-name mapping.
        strategy_file (str): Strategy FITS file with ``zfibermax_cut``.
        strategy_index (int): Strategy row index.
        cut_color (bool, optional): Also apply the colour cuts.

    Returns:
        numpy.ndarray: Boolean mask of galaxies passing the selection.

    Raises:
        ValueError: If ``strategy_file`` is None.
    """
    if strategy_file is None:
        raise ValueError("No strategy file for DESI2 magnitude masking")

    key_zfiber = hashing_table[band_DESI2_1]
    key_z = hashing_table[band_DESI2_2]
    key_W1 = hashing_table[band_DESI2_3]
    key_r = hashing_table[band_DESI2_4]
    key_g = hashing_table[band_DESI2_5]

    color_to_cut_1 = (
        file_sn[key_z] - file_sn[key_W1] - 0.8 * (file_sn[key_r] - file_sn[key_z])
    )
    color_to_cut_2 = file_sn[key_g] - file_sn[key_W1]
    color_to_cut_3 = file_sn[key_r] - file_sn[key_W1]

    mag_cut = fitsio.FITS(strategy_file)[1]["zfibermax_cut"][strategy_index]
    mask_magnitude = file_sn[key_zfiber] < mag_cut

    if cut_color:
        mask_magnitude &= color_to_cut_1 > color_cut_DESI2_1  # LRG non-stellar cut
        mask_magnitude &= (color_to_cut_2 > color_cut_DESI2_2) | (
            color_to_cut_3 > color_cut_DESI2_3  # low-z cuts
        )

    return mask_magnitude


def host_efficiency_DESI2_simple():
    """DESI2 constant host efficiency (fibre assignment x redshift success).

    Returns:
        float: The efficiency product.
    """
    return fiber_assignement_efficiency_DESI2 * redshift_success_rate_DESI2


def n_z_DESI2_from_target_selection(
    n_z_file,
    strategy_index,
):
    """Build the DESI2 ``n(z)`` from a survey-strategy FITS file.

    Evaluates the fitted :func:`desidescsn.efficiency.n_z_function` with the
    strategy's parameters, scaled by its target density.

    Args:
        n_z_file (str): Strategy FITS file (``Nspec/sqdeg``, ``nzfit_params_ccut``).
        strategy_index (int): Strategy row index.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: ``(z_centers, n_z)``.
    """
    desi2_strategy_file = fitsio.FITS(n_z_file)[1]
    density = desi2_strategy_file["Nspec/sqdeg"][strategy_index]
    n_z_params = desi2_strategy_file["nzfit_params_ccut"][strategy_index]

    z_centers = np.linspace(0, 1.2, 100)
    n_z = efficiency.n_z_function(z_centers, 1.0, *n_z_params) * density
    return z_centers, n_z


# CRS functions


def mask_magnitude_CRSBG(
    file_sn,
    hashing_table,
    cut_color=True,
):
    """Apply the 4MOST CRS bright-galaxy magnitude and colour selection.

    Args:
        file_sn (pandas.DataFrame): Galaxy sample.
        hashing_table (dict): Band-name -> column-name mapping (``J``/``K``/``W1``).
        cut_color (bool, optional): Also apply the colour cuts.

    Returns:
        numpy.ndarray: Boolean mask of galaxies passing the selection.
    """
    key_J = hashing_table[band_CRSBG_1]
    key_K = hashing_table[band_CRSBG_2]
    key_W1 = hashing_table[band_CRSBG_3]

    color_to_cut_1 = file_sn[key_J] - file_sn[key_K]
    color_to_cut_2 = file_sn[key_J] - file_sn[key_W1]

    mask_magnitude = file_sn[key_J] > magnitude_cut_CRSBG_1
    mask_magnitude &= file_sn[key_K] < magnitude_cut_CRSBG_2
    if cut_color:
        mask_magnitude &= color_to_cut_1 - 1.6 * color_to_cut_2 > color_cut_CRSBG_1
        mask_magnitude &= color_to_cut_1 - 1.6 * color_to_cut_2 < color_cut_CRSBG_2
        mask_magnitude &= color_to_cut_1 + 2.5 * color_to_cut_2 > color_cut_CRSBG_3
        mask_magnitude &= color_to_cut_1 + 0.5 * color_to_cut_2 > color_cut_CRSBG_4

    return mask_magnitude


def host_efficiency_CRSBG_simple():
    """4MOST CRS BG constant host efficiency (fibre assignment x z success).

    Returns:
        float: The efficiency product.
    """
    return fiber_assignement_efficiency_CRSBG * redshift_success_rate_CRSBG


def mask_magnitude_CRSLRG(
    file_sn,
    hashing_table,
    cut_color=True,
):
    """Apply the 4MOST CRS LRG magnitude and colour selection.

    Args:
        file_sn (pandas.DataFrame): Galaxy sample.
        hashing_table (dict): Band-name -> column-name mapping (``J``/``r``/``W1``).
        cut_color (bool, optional): Also apply the colour cut.

    Returns:
        numpy.ndarray: Boolean mask of galaxies passing the selection.
    """
    key_J = hashing_table[band_CRSLRG_1]
    key_r = hashing_table[band_CRSLRG_2]
    key_W1 = hashing_table[band_CRSLRG_3]
    color_to_cut_1 = file_sn[key_r] - file_sn[key_W1]

    mask_magnitude = file_sn[key_J] > magnitude_cut_CRSLRG_1
    mask_magnitude &= file_sn[key_J] < magnitude_cut_CRSLRG_2
    if cut_color:
        mask_magnitude &= color_to_cut_1 - 0.5 * file_sn[key_W1] > color_cut_CRSLRG_1

    return mask_magnitude


def host_efficiency_CRSLRG_simple():
    """4MOST CRS LRG constant host efficiency (fibre assignment x z success).

    Returns:
        float: The efficiency product.
    """
    return fiber_assignement_efficiency_CRSLRG * redshift_success_rate_CRSLRG


def n_z_CRS_from_target_selection(
    n_z_file,
    zmax=1.2,
):
    """Load the 4MOST CRS ``n(z)`` from a two-column text file.

    Args:
        n_z_file (str): Text file with redshift and density columns.
        zmax (float, optional): Upper redshift cut.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: ``(z_centers, n_z)`` below ``zmax``.
    """
    n_redshift = np.loadtxt(n_z_file)
    mask = n_redshift[:, 0] < zmax
    z_centers = n_redshift[:, 0][mask]
    n_z = n_redshift[:, 1][mask]
    return z_centers, n_z


# 4HS functions


def mask_magnitude_4HS(
    file_sn,
    hashing_table,
    cut_color=True,
):
    """Apply the 4HS magnitude and colour selection to a galaxy sample.

    Args:
        file_sn (pandas.DataFrame): Galaxy sample.
        hashing_table (dict): Band-name -> column-name mapping (``J``/``K``).
        cut_color (bool, optional): Also apply the J-K colour cut.

    Returns:
        numpy.ndarray: Boolean mask of galaxies passing the selection.
    """
    key_J = hashing_table[band_4HS_1]
    key_K = hashing_table[band_4HS_2]

    color_to_cut = file_sn[key_J] - file_sn[key_K]
    mask_magnitude = file_sn[key_J] < magnitude_cut_4HS
    if cut_color:
        mask_magnitude &= color_to_cut < color_cut_4HS

    return mask_magnitude


def host_efficiency_4HS_simple():
    """4HS constant host efficiency (fibre assignment x redshift success).

    Returns:
        float: The efficiency product.
    """
    return fiber_assignement_efficiency_4HS * redshift_success_rate_4HS


def n_z_4HS_from_target_selection(
    n_z_file,
    redshift_name="Z_HELIO",
    zmax=1.2,
    bins=50,
    area_gama=180,
):
    """Build the 4HS ``n(z)`` from a GAMA redshift catalog.

    Args:
        n_z_file (str): GAMA FITS catalog.
        redshift_name (str, optional): Redshift column name.
        zmax (float, optional): Histogram upper edge.
        bins (int, optional): Number of redshift bins.
        area_gama (float, optional): GAMA area (deg^2) for normalisation.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: ``(z_centers, n_z)`` (per deg^2 per
        unit redshift).
    """
    z_gama = fitsio.FITS(n_z_file)[1][redshift_name][:]
    hist, edges = np.histogram(z_gama, bins=bins, range=(0.0, zmax))
    z_centers = (edges[:-1] + edges[1:]) / 2
    delta_z = np.mean(z_centers[1:] - z_centers[:-1])
    n_z = hist / (area_gama * delta_z)

    return z_centers, n_z


### General functions:


def get_n_z_survey(
    method,
    survey,
    n_z_file=None,
    redshift_efficiency_file=None,
    strategy_index=0,
):
    """Dispatch to the right ``n(z)`` loader for a survey and method.

    Args:
        method (str): ``"target_selection"`` (public predictions) or ``"y1"``
            (private DESI Y1 data).
        survey (str): Survey identifier (e.g. ``"DESIBGS"``, ``"CRSLRG"``).
        n_z_file (str, optional): Survey ``n(z)`` file.
        redshift_efficiency_file (str, optional): Redshift-efficiency yaml file.
        strategy_index (int, optional): Strategy row index (DESI2).

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: ``(redshift_survey, n_z_survey)``.

    Raises:
        ValueError: If the method is not available for the survey.
    """
    if method == "target_selection":
        if survey == "DESIBGS":
            redshift_survey, n_z_survey = n_z_DESI_from_target_selection(
                "BGS",
                n_z_file,
                redshift_efficiency_file,
            )
        elif survey == "DESILRG":
            redshift_survey, n_z_survey = n_z_DESI_from_target_selection(
                "LRG",
                n_z_file,
                redshift_efficiency_file,
            )
        elif survey == "DESI2":
            redshift_survey, n_z_survey = n_z_DESI2_from_target_selection(
                n_z_file,
                strategy_index,
            )
        elif (survey == "CRSBG") | (survey == "CRSLRG"):
            redshift_survey, n_z_survey = n_z_CRS_from_target_selection(
                n_z_file,
            )
        elif survey == "4HS":
            redshift_survey, n_z_survey = n_z_4HS_from_target_selection(
                n_z_file,
            )

    elif method == "y1":
        print("WARNING: This is private data, be sure to be part of DESI")
        if survey == "DESIBGS":
            redshift_survey, n_z_survey = n_z_DESI_from_y1(
                "BGS",
                n_z_file,
                redshift_efficiency_file,
            )
        elif survey == "DESILRG":
            redshift_survey, n_z_survey = n_z_DESI_from_y1(
                "LRG",
                n_z_file,
                redshift_efficiency_file,
            )
        else:
            raise ValueError(
                f"The method {method} is not available for the survey {survey}"
            )

    return redshift_survey, n_z_survey
