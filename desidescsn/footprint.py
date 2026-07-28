"""Survey sky-footprint masks for the DESI/DESC SN forecasts.

Build boolean sky masks for the survey footprints considered (DESI, DESI
extended, DESI2, 4MOST CRS, 4HS), either from tile files via ``desimodel`` or by
degrading HEALPix footprint maps, plus a helper to lay out an empty HEALPix map.
Coordinates are RA/Dec in degrees; HEALPix maps use the NESTED scheme.
"""

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
    """Whether points fall inside the DESI tiles of a given program.

    Args:
        desi_footprint_file (str): DESI tiles file.
        ra (array-like): Right ascension (deg).
        dec (array-like): Declination (deg).
        program (str, optional): DESI program (e.g. ``"BRIGHT"``).

    Returns:
        numpy.ndarray: Boolean in-footprint mask.
    """
    tiles = Table.read(desi_footprint_file)
    selection = (tiles["IN_DESI"] == 1) & (tiles["PROGRAM"] == program)
    tiles_selected = tiles[selection]
    return is_point_in_desi(tiles_selected, ra, dec)


def get_desiext_mask(
    desiext_footprint_file,
    nside,
    field="DESI1B_BRIGHT",
):
    """Build a mask from a DESI-extended HEALPix footprint map.

    Args:
        desiext_footprint_file (str): FITS footprint file.
        nside (int): Target HEALPix NSIDE (the map is degraded to it).
        field (str, optional): Footprint column name.

    Returns:
        numpy.ndarray: Boolean mask (True where the footprint is non-zero).
    """
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
    """Build a mask from a DESI2 HEALPix footprint map.

    Args:
        desi2_footprint_file (str): FITS footprint file.
        nside (int): Target HEALPix NSIDE (the map is degraded to it).
        field (str, optional): Footprint column name.

    Returns:
        numpy.ndarray: Boolean mask (True where the footprint is non-zero).
    """
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
    """Build a mask from the 4MOST CRS low-surface-material HEALPix map.

    Args:
        crs_footprint_file (str): FITS footprint file (``LSM`` column).
        nside (int): Target HEALPix NSIDE (the map is degraded to it).

    Returns:
        numpy.ndarray: Boolean mask (True where the map is non-zero).
    """
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
    """Build the 4HS footprint mask (southern sky, away from the galactic plane).

    Selects points with ``dec < 0`` and galactic latitude ``|b| > 20 deg``.

    Args:
        ra (array-like): Right ascension (deg).
        dec (array-like): Declination (deg).

    Returns:
        numpy.ndarray: Boolean in-footprint mask.
    """
    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame="icrs").galactic
    b = np.array([c[i].b.value for i in range(len(c))])

    mask = dec < 0
    mask &= np.abs(b) > 20.0

    return mask


def create_healpix_map(nside):
    """Create an empty HEALPix map and its per-pixel RA/Dec (NESTED scheme).

    Args:
        nside (int): HEALPix NSIDE.

    Returns:
        tuple: ``(ra, dec, hpx_map)`` — per-pixel RA/Dec (deg) and a zero map of
        length ``npix``.
    """
    npix = hp.nside2npix(nside)
    ra, dec = hp.pix2ang(
        nside,
        range(npix),
        nest=True,
        lonlat=True,
    )
    hpx_map = np.zeros(npix)
    return ra, dec, hpx_map
