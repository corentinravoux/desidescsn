"""Host- and SN-redshift efficiency modelling for DESI/DESC SN forecasts.

Two efficiency ingredients are combined to forecast how many type-Ia supernova
hosts a spectroscopic survey recovers a redshift for:

* the **host redshift efficiency** — how the target-selection cuts of a survey
  reweight the simulated host-galaxy redshift distribution towards the survey's
  observed ``n(z)`` (via a fitted parametric ``n(z)`` model);
* the **SN redshift efficiency** — on a simulated SN-host sample, the fraction
  passing a survey's magnitude/colour cuts as a function of redshift.

Survey-specific cuts and ``n(z)`` come from :mod:`desidescsn.surveys`.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma

from desidescsn import surveys

### Host redshift efficiency: Normalize a sample of galaxy obtain from cuts to an input n(z)


def n_z_function(z, A, z0, beta, d):
    """Parametric redshift distribution ``n(z)`` (Smail-type).

    Args:
        z (float | array-like): Redshift.
        A (float): Overall amplitude.
        z0 (float): Characteristic redshift.
        beta (float): High-redshift fall-off exponent.
        d (float): Low-redshift rise exponent.

    Returns:
        float | array-like: The number density at ``z``.
    """
    return (
        A
        * (beta / gamma(d / beta))
        * (z ** (d - 1) / z0**d)
        * np.exp(-((z / z0) ** beta))
    )


def fit_n_z(xdata, ydata):
    """Fit the :func:`n_z_function` parameters to a measured ``n(z)``.

    Args:
        xdata (numpy.ndarray): Redshift bin centres.
        ydata (numpy.ndarray): Measured number density per bin.

    Returns:
        numpy.ndarray: Best-fit ``[A, z0, beta, d]`` parameters.
    """
    p0 = [np.max(ydata), 0.3, 2.0, 2.0]
    popt, _ = curve_fit(
        n_z_function,
        xdata=xdata,
        ydata=ydata,
        p0=p0,
        bounds=([0.0, 0.0, 0.0, 0.0], [np.inf, np.inf, np.inf, np.inf]),
        maxfev=10000,
    )
    return popt


def compute_host_efficiency(
    redshift_survey,
    n_z_survey,
    redshift_simu,
    n_z_simu,
    maximum_redshift_efficiency=None,
):
    """Build the host redshift-efficiency function from two ``n(z)`` samples.

    Fits :func:`n_z_function` to the survey and simulation redshift
    distributions and returns their ratio ``n_survey(z) / n_simu(z)`` as a
    callable, optionally capped at a maximum value.

    Args:
        redshift_survey (numpy.ndarray): Survey redshift bin centres.
        n_z_survey (numpy.ndarray): Survey number density.
        redshift_simu (numpy.ndarray): Simulation redshift bin centres.
        n_z_simu (numpy.ndarray): Simulation number density.
        maximum_redshift_efficiency (float, optional): Upper cap on the ratio.

    Returns:
        callable: ``efficiency(z)`` giving the host redshift efficiency.
    """
    n_z_parameter_survey = fit_n_z(redshift_survey, n_z_survey)
    n_z_parameter_simu = fit_n_z(redshift_simu, n_z_simu)

    redshift_efficiency = lambda z: n_z_function(
        z, *n_z_parameter_survey
    ) / n_z_function(z, *n_z_parameter_simu)
    if maximum_redshift_efficiency is not None:
        redshift_efficiency_maximized = np.vectorize(
            lambda z: np.min([maximum_redshift_efficiency, redshift_efficiency(z)])
        )
        return redshift_efficiency_maximized
    else:
        return redshift_efficiency


def compute_host_efficiency_survey(
    method,
    survey,
    n_z_file=None,
    redshift_efficiency_file=None,
    redshift_array_simu=None,
    area_simu=None,
    maximum_redshift_efficiency=None,
    strategy_index=0,
):
    """Return the host redshift-efficiency function for a survey.

    For the ``"simple"`` method this is the constant fibre-assignment x
    redshift-success product of the survey; otherwise the simulated host ``n(z)``
    is histogrammed and combined with the survey ``n(z)`` (via
    :func:`compute_host_efficiency`).

    Args:
        method (str): ``"simple"`` or a data-driven method
            (``"target_selection"`` / ``"y1"``, see
            :func:`desidescsn.surveys.get_n_z_survey`).
        survey (str): Survey identifier (e.g. ``"DESIBGS"``).
        n_z_file (str, optional): Survey ``n(z)`` file.
        redshift_efficiency_file (str, optional): Redshift-efficiency yaml file.
        redshift_array_simu (numpy.ndarray, optional): Simulated host redshifts.
        area_simu (float, optional): Simulated area (deg^2) for normalisation.
        maximum_redshift_efficiency (float, optional): Upper cap on the ratio.
        strategy_index (int, optional): Survey-strategy row index.

    Returns:
        callable: ``efficiency(z)`` giving the host redshift efficiency.
    """
    if method == "simple":
        return lambda z: eval(f"surveys.host_efficiency_{survey}_simple()")

    else:
        n_z_simu, redshift_simu = np.histogram(redshift_array_simu, bins=30)
        redshift_simu_centers = (redshift_simu[:-1] + redshift_simu[1:]) / 2
        delta_z = np.mean(redshift_simu_centers[1:] - redshift_simu_centers[:-1])
        n_z_simu = n_z_simu / (area_simu * delta_z)

        (
            redshift_survey_centers,
            n_z_survey,
        ) = surveys.get_n_z_survey(
            method,
            survey,
            n_z_file=n_z_file,
            redshift_efficiency_file=redshift_efficiency_file,
            strategy_index=strategy_index,
        )
        redshift_efficiency = compute_host_efficiency(
            redshift_survey_centers,
            n_z_survey,
            redshift_simu_centers,
            n_z_simu,
            maximum_redshift_efficiency=maximum_redshift_efficiency,
        )
        return redshift_efficiency


### SN redshift efficiency: On a given sample, efficiency of getting the SNIa host redshift.


def get_host_eff_one_band_cut(
    file,
    seed,
    model_weights,
    parameters_weigths,
    magnitude_cut,
    band,
    number_years,
    rate,
    cosmo,
    N_z=30,
):
    """Compute the SN redshift efficiency after a single-band magnitude cut.

    Draws SN explosions on a host sample, applies a magnitude cut in one band
    and returns the recovered fraction vs redshift.

    Args:
        file (pandas.DataFrame): Host-galaxy sample.
        seed (int): RNG seed for the SN draw.
        model_weights (str): Weight model (see :func:`return_weight_model`).
        parameters_weigths (dict): Weight-model parameters.
        magnitude_cut (float): Magnitude threshold.
        band (str): Column name of the magnitude band to cut on.
        number_years (float): Survey duration (years).
        rate (float): SN rate (per volume per year).
        cosmo: Astropy cosmology (for comoving volume).
        N_z (int, optional): Number of redshift bins.

    Returns:
        tuple: ``(redshift, efficiency, file_sn)`` — bin centres, recovered
        fraction and the exploded-SN sub-sample.
    """
    mask_sn_explosion = sn_explosion(
        file,
        seed,
        model_weights,
        parameters_weigths,
        number_years,
        rate,
        cosmo,
    )
    file_sn = file[mask_sn_explosion]

    mask_magnitude = file_sn[band] < magnitude_cut

    redshift, efficiency = compute_sn_efficiency(
        file_sn,
        mask_magnitude,
        N_z=N_z,
    )

    return redshift, efficiency, file_sn


def sn_explosion(
    file,
    seed,
    model_weights,
    parameters_weigths,
    number_years,
    rate,
    cosmo,
):
    """Draw which host galaxies host an observable SN Ia over a survey window.

    Estimates the expected number of SNe from the comoving volume, rate and
    duration, then randomly selects host galaxies weighted by the chosen SN
    weight model.

    Note:
        Only works for a rectangular RA/Dec area.

    Args:
        file (pandas.DataFrame): Host-galaxy sample (needs ``ra``/``dec``/
            ``redshift_true`` and the columns used by the weight model).
        seed (int): RNG seed.
        model_weights (str): Weight model (see :func:`return_weight_model`).
        parameters_weigths (dict): Weight-model parameters.
        number_years (float): Survey duration (years).
        rate (float): SN rate (per volume per year).
        cosmo: Astropy cosmology (for comoving distance).

    Returns:
        numpy.ndarray: Boolean mask of host galaxies with a drawn SN.
    """
    # CR - only works for rectangular area.
    area = (file["ra"].max() - file["ra"].min()) * (
        file["dec"].max() - file["dec"].min()
    )
    vol = (
        area
        / (360.0 * 360.0 / np.pi)
        * 4.0
        * np.pi
        / 3.0
        * (
            cosmo.comoving_distance(file["redshift_true"].max()).value ** 3
            - cosmo.comoving_distance(file["redshift_true"].min()).value ** 3
        )
    )
    N_sn = int(rate * number_years * vol)

    # Define weights for A + B model, remove passive galaxies,
    weights = return_weight_model(file, model_weights, parameters_weigths)
    weights = np.array(weights).astype("float64")
    weights = weights / np.sum(weights)

    np.random.seed(seed)
    indexes = np.random.choice(file.index, p=weights, size=N_sn)
    mask_sn_explosion = np.in1d(file.index, indexes)

    return mask_sn_explosion


def return_weight_model(
    file,
    model,
    parameters,
):
    """Return per-galaxy SN-host weights for a chosen SN model.

    Args:
        file (pandas.DataFrame): Host-galaxy sample (needs ``stellar_mass`` and
            ``sfr`` for the non-trivial models).
        model (str): One of ``"noweigths"`` (uniform), ``"snia"``
            (A*M + B*SFR), ``"snia_pec"`` (A*M + B*SFR above an sSFR cut) or
            ``"sncc"`` (M^C above an sSFR cut).
        parameters (dict): Model parameters (``A``/``B``/``C``/``sSFR_cut``).

    Returns:
        numpy.ndarray: Per-galaxy weights (unnormalised).
    """
    if model == "noweigths":
        weights = np.ones(file.shape[0])
        
    elif model == "snia":
        A = parameters["A"]
        B = parameters["B"]
        weights = A * file["stellar_mass"] + B * file["sfr"]
        
    elif model == "snia_pec":
        A = parameters["A"]
        B = parameters["B"]
        sSFR_cut = parameters["sSFR_cut"] # -11.5 from Vincenxi et al. 2020 https://academic.oup.com/mnras/article/505/2/2819/6284776
        weights = np.zeros(len(file["stellar_mass"]))
        mask = np.log10(file['sfr']/file['stellar_mass']) >  sSFR_cut 
        weights[mask] =  A * file["stellar_mass"][mask] + B * file["sfr"][mask]
        
    elif model == "sncc":
        C = parameters['C'] # 0.16 for SNeII , 0.36 for SNIb/c from Vincenxi et al. 2020 https://academic.oup.com/mnras/article/505/2/2819/6284776
        sSFR_cut = parameters["sSFR_cut"] # -11.5 from Vincenxi et al. 2020 https://academic.oup.com/mnras/article/505/2/2819/6284776
        weights = np.zeros(len(file["stellar_mass"]))
        mask = np.log10(file['sfr']/file['stellar_mass']) > sSFR_cut 
        weights[mask] = file["stellar_mass"][mask] **C
    
    return weights


def compute_sn_efficiency(
    file_sn,
    mask_magnitude,
    N_z=30,
    redshift_range=None,
):
    """Fraction of SN hosts passing a magnitude mask, as a function of redshift.

    Args:
        file_sn (pandas.DataFrame): SN-host sample (needs ``redshift_true``).
        mask_magnitude (numpy.ndarray): Boolean mask of hosts passing the cut.
        N_z (int, optional): Number of redshift bins.
        redshift_range (tuple[float, float], optional): Restrict the masked
            hosts to this redshift window.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: ``(z_centers, efficiency)`` — the
        recovered fraction per redshift bin.
    """
    if type(mask_magnitude) == list:
        print("Several masks not implemented for now")

    else:
        bins_z = np.linspace(
            file_sn["redshift_true"].min(), file_sn["redshift_true"].max(), N_z
        )
        if redshift_range is not None:
            mask_redshift = (
                file_sn["redshift_true"][mask_magnitude] > redshift_range[0]
            ) & (file_sn["redshift_true"][mask_magnitude] < redshift_range[1])
        else:
            mask_redshift = np.full(
                file_sn["redshift_true"][mask_magnitude].shape, True
            )
        bins_z_centers = (bins_z[1:] + bins_z[:-1]) / 2
        count_z, _ = np.histogram(file_sn["redshift_true"], bins=bins_z)
        count_z_masked, _ = np.histogram(
            file_sn["redshift_true"][mask_magnitude][mask_redshift], bins=bins_z
        )
    return bins_z_centers, count_z_masked / count_z


def compute_full_efficiency(
    survey,
    hashing_table,
    file_simu_sn,
    file_simu,
    host_efficiency_estimator,
    n_z_file=None,
    redshift_efficiency_file=None,
    area_simu=None,
    maximum_redshift_efficiency=None,
    N_z=30,
    cut_color=True,
    strategy_file=None,
    strategy_index=0,
):
    """Combine SN and host redshift efficiencies into the full survey efficiency.

    Applies the survey magnitude/colour mask to both the SN-host and full host
    samples, computes the SN efficiency vs redshift and multiplies it by the
    host redshift efficiency.

    Args:
        survey (str): Survey identifier (selects ``surveys.mask_magnitude_*``).
        hashing_table (dict): Band-name -> dataframe-column mapping.
        file_simu_sn (pandas.DataFrame): Simulated SN-host sample.
        file_simu (pandas.DataFrame): Simulated full host sample.
        host_efficiency_estimator (str): Host-efficiency method (``"simple"`` ...).
        n_z_file (str, optional): Survey ``n(z)`` file.
        redshift_efficiency_file (str, optional): Redshift-efficiency yaml file.
        area_simu (float, optional): Simulated area (deg^2).
        maximum_redshift_efficiency (float, optional): Cap on the host efficiency.
        N_z (int, optional): Number of redshift bins.
        cut_color (bool, optional): Apply the colour cuts as well as magnitude.
        strategy_file (str, optional): DESI2-style strategy file.
        strategy_index (int, optional): Strategy row index.

    Returns:
        tuple: ``(redshift, normalized_efficiency, sn_efficiency,
        host_efficiency)``.
    """
    if strategy_file is not None:
        mask_magnitude_sn = eval(f"surveys.mask_magnitude_{survey}")(
            file_simu_sn,
            hashing_table,
            cut_color=cut_color,
            strategy_file=strategy_file,
            strategy_index=strategy_index,
        )

        mask_magnitude_full = eval(f"surveys.mask_magnitude_{survey}")(
            file_simu,
            hashing_table,
            cut_color=cut_color,
            strategy_file=strategy_file,
            strategy_index=strategy_index,
        )
    else:
        mask_magnitude_sn = eval(f"surveys.mask_magnitude_{survey}")(
            file_simu_sn,
            hashing_table,
            cut_color=cut_color,
        )

        mask_magnitude_full = eval(f"surveys.mask_magnitude_{survey}")(
            file_simu,
            hashing_table,
            cut_color=cut_color,
        )

    redshift, sn_efficiency = compute_sn_efficiency(
        file_simu_sn,
        mask_magnitude_sn,
        N_z=N_z,
    )

    host_efficiency = compute_host_efficiency_survey(
        host_efficiency_estimator,
        survey,
        n_z_file=n_z_file,
        redshift_efficiency_file=redshift_efficiency_file,
        redshift_array_simu=file_simu["redshift_true"][mask_magnitude_full],
        area_simu=area_simu,
        maximum_redshift_efficiency=maximum_redshift_efficiency,
        strategy_index=strategy_index,
    )

    normalized_efficiency = sn_efficiency * host_efficiency(redshift)

    return redshift, normalized_efficiency, sn_efficiency, host_efficiency(redshift)
