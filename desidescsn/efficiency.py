import numpy as np


def get_redshift_eff_one_band_cut(
    file,
    seed,
    A,
    B,
    magnitude_cut,
    band,
    number_years,
    rate,
    cosmo,
    N_z=30,
):

    zmin, zmax = file["redshift_true"].min(), file["redshift_true"].max()
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
            cosmo.comoving_distance(zmax).value ** 3
            - cosmo.comoving_distance(zmin).value ** 3
        )
    )
    N_sn = int(rate * number_years * vol)

    # Define weights for A + B model, remove passive galaxies,
    weights = A * file["stellar_mass"] + B * file["sfr"]
    weights = np.array(weights).astype("float64")
    weights = weights / np.sum(weights)

    np.random.seed(seed)
    indexes = np.random.choice(file.index, p=weights, size=N_sn)
    mask_indexes = np.in1d(file.index, indexes)
    file_sn = file[mask_indexes]

    mask_magnitude = file_sn[band] < magnitude_cut

    bins_z = np.linspace(zmin, zmax, N_z)
    bins_z_centers = (bins_z[1:] + bins_z[:-1]) / 2
    count_z, _ = np.histogram(file_sn["redshift_true"], bins=bins_z)
    count_z_masked, _ = np.histogram(
        file_sn["redshift_true"][mask_magnitude], bins=bins_z
    )

    return bins_z_centers, count_z_masked / count_z, file_sn
