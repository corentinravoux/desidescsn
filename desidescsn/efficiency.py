import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma

from desidescsn import surveys

### Host redshift efficiency: Normalize a sample of galaxy obtain from cuts to an input n(z)


def n_z_function(z, A, z0, beta, d):
    return (
        A
        * (beta / gamma(d / beta))
        * (z ** (d - 1) / z0**d)
        * np.exp(-((z / z0) ** beta))
    )


def fit_n_z(xdata, ydata):
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
    if model == "noweigths":
        weights = np.ones(file.shape[0])
    elif model == "ab":
        A = parameters["A"]
        B = parameters["B"]
        weights = A * file["stellar_mass"] + B * file["sfr"]

    return weights


def compute_sn_efficiency(
    file_sn,
    mask_magnitude,
    N_z=30,
    redshift_range=None,
):
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
