import numpy as np


def get_redshift_eff_one_band_cut(
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

    redshift, efficiency = compute_efficiency(
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


def return_weight_model(file, model, parameters):
    if model == "noweigths":
        weights = np.ones(file.shape[0])
    elif model == "ab":
        A = parameters["A"]
        B = parameters["B"]
        weights = A * file["stellar_mass"] + B * file["sfr"]

    return weights


def compute_efficiency(file_sn, mask_magnitude, N_z=30, redshift_range=None):
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
            mask_redshift = np.full(file_sn["redshift_true"].shape, True)
        bins_z_centers = (bins_z[1:] + bins_z[:-1]) / 2
        count_z, _ = np.histogram(file_sn["redshift_true"], bins=bins_z)
        count_z_masked, _ = np.histogram(
            file_sn["redshift_true"][mask_magnitude][mask_redshift], bins=bins_z
        )
    return bins_z_centers, count_z_masked / count_z
