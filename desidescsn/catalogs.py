"""Extract host-galaxy properties from DESC DC2 / Roman-Rubin mock catalogs.

Thin wrapper over ``GCRCatalogs`` that pulls a fixed feature list (photometry in
LSST/SDSS/Roman bands, redshift, star-formation rate and stellar mass) from a
DESC mock catalog into a pandas dataframe and writes it to parquet, ready for
the efficiency forecasts in :mod:`desidescsn.efficiency`.
"""

import GCRCatalogs
import numpy as np
import pandas as pd

name_catalog_roman_rubin = "roman_rubin_2023_v1.1.1_elais"

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

feature_list_roman_rubin_light = np.array(
    [
        "ra",
        "dec",
        "redshift_true",
        "mag_true_g_sdss",
        "mag_true_r_sdss",
        "mag_true_i_sdss",
        "mag_true_z_sdss",
        "ROMAN_obs_J129",
        "ROMAN_obs_K213",
        "sfr",
        "stellar_mass",
    ]
)


def extract_dc2_properties(name_catalog, feature_list, path_out, redshift_cut=None):
    """Extract a feature list from a DESC mock catalog and write it to parquet.

    Args:
        name_catalog (str): GCR catalog name (e.g. ``name_catalog_roman_rubin``).
        feature_list (Sequence[str]): Quantities to extract (e.g.
            ``feature_list_roman_rubin``).
        path_out (str): Output parquet file path.
        redshift_cut (float, optional): Keep only ``redshift_true`` below this.
    """
    gc = GCRCatalogs.load_catalog(name_catalog)
    dict_quantities = gc.get_quantities(feature_list)
    df = pd.DataFrame(dict_quantities)
    if redshift_cut is not None:
        df = df[df["redshift_true"] < redshift_cut]
    df.to_parquet(path_out)
