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


def get_desi2_mask(
    ra,
    dec,
    variant=0,
):
    mask = (ra < 50) & (-15 < dec) & (dec < 3)
    mask |= (ra > -50 + 360) & (-15 < dec) & (dec < 3)
    mask |= (ra > 130) & (ra < 220) & (-15 < dec) & (dec <= 0)
    mask |= (ra > 120) & (ra < 250) & (0 <= dec) & (dec < 15)
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
