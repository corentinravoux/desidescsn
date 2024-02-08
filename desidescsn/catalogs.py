import GCRCatalogs
import numpy as np
import pandas as pd

feature_list_roman_rubin = np.array(
    [
        "galaxy_id",
        "ra",
        "dec",
        "redshift",
        "redshift_true",
        "mag_true_u_lsst",
        "mag_true_g_lsst",
        "mag_true_r_lsst",
        "mag_true_i_lsst",
        "mag_true_z_lsst",
        "mag_true_Y_lsst",
        "mag_true_u_sdss",
        "mag_true_g_sdss",
        "mag_true_r_sdss",
        "mag_true_i_sdss",
        "mag_true_z_sdss",
        "ROMAN_obs_F184",
        "ROMAN_obs_H158",
        "ROMAN_obs_J129",
        "ROMAN_obs_K213",
        "ROMAN_obs_R062",
        "ROMAN_obs_W146",
        "ROMAN_obs_Y106",
        "ROMAN_obs_Z087",
        "sfr",
        "stellar_mass",
    ]
)


def extract_dc2_properties(name_catalog, feature_list, path_out):
    gc = GCRCatalogs.load_catalog(name_catalog)
    dict_quantities = gc.get_quantities(feature_list)
    df = pd.DataFrame(dict_quantities)
    df.to_parquet(path_out)
