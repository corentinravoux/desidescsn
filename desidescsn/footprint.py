import healpy as hp
import numpy as np
from desimodel.footprint import is_point_in_desi


def IndexToDeclRa(index, nside):
    theta, phi = hp.pixelfunc.pix2ang(nside, index)
    return np.degrees(np.pi * 2.0 - phi), -np.degrees(theta - np.pi / 2.0)


def DeclRaToIndex(RA, DEC, nside):
    return hp.pixelfunc.ang2pix(nside, np.radians(-DEC + 90.0), np.radians(360.0 - RA))


def plot_property(property_name, tiles_selected, npix, nside, radius):
    property_map = np.full(npix, np.nan)
    for i in range(len(tiles_selected)):
        ind = DeclRaToIndex(tiles_selected["RA"][i], tiles_selected["DEC"][i], nside)
        vec = hp.pix2vec(nside, ind)
        ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(radius))
        property_map[ipix_disc] = tiles_selected[property_name][i]

    hp.mollview(
        (property_map), hold=True, nest=False, notext=True, fig=1, title=property_name
    )
    hp.graticule()


def get_desi_footprint(nside, desi_footprint_file):
    npix = hp.nside2npix(nside)
    index = np.arange(npix)
    ra, dec = IndexToDeclRa(index, nside)
    mask = is_point_in_desi(tiles_selected, ra, dec)
