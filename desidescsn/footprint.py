import fitsio
import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from desimodel.footprint import is_point_in_desi


def get_desi_mask(
    desi_footprint_file,
    ra,
    dec,
    program="BRIGHT",
):
    tiles = Table.read(desi_footprint_file)
    selection = (tiles["IN_DESI"] == 1) & (tiles["PROGRAM"] == program)
    tiles_selected = tiles[selection]
    return is_point_in_desi(tiles_selected, ra, dec)


def get_desiext_mask(
    desiext_footprint_file,
    nside,
    field="DESI1B_BRIGHT",
):
    desiext_footprint = fitsio.FITS(desiext_footprint_file)[1][field][:]

    desiext_footprint_resolved = hp.pixelfunc.ud_grade(
        desiext_footprint,
        nside,
        pess=False,
        order_in="NEST",
        order_out=None,
    )
    mask = desiext_footprint_resolved != 0.0
    return mask


def get_desi2_mask(
    desi2_footprint_file,
    nside,
    field="IBIS_WIDE",
):
    desi2_footprint = fitsio.FITS(desi2_footprint_file)[1][field][:]

    desi2_footprint_resolved = hp.pixelfunc.ud_grade(
        desi2_footprint,
        nside,
        pess=False,
        order_in="NEST",
        order_out=None,
    )
    mask = desi2_footprint_resolved != 0.0
    return mask


def get_crs_mask(crs_footprint_file, nside):
    lsm_map = fitsio.FITS(crs_footprint_file)[1]["LSM"][:]
    lsm_map_resolved = hp.pixelfunc.ud_grade(
        lsm_map,
        nside,
        pess=False,
        order_in="NEST",
        order_out=None,
    )
    mask = lsm_map_resolved != 0.0
    return mask


def get_4hs_mask(
    ra,
    dec,
):
    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame="icrs").galactic
    b = np.array([c[i].b.value for i in range(len(c))])

    mask = dec < 0
    mask &= np.abs(b) > 20.0

    return mask


def create_healpix_map(nside):
    npix = hp.nside2npix(nside)
    ra, dec = hp.pix2ang(
        nside,
        range(npix),
        nest=True,
        lonlat=True,
    )
    hpx_map = np.zeros(npix)
    return ra, dec, hpx_map
