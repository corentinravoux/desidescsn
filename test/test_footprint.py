"""Unit tests for the offline parts of :mod:`desidescsn.footprint`."""

import numpy as np

from desidescsn import footprint


def test_create_healpix_map():
    nside = 16
    ra, dec, hpx = footprint.create_healpix_map(nside)
    npix = 12 * nside**2
    assert len(ra) == len(dec) == len(hpx) == npix
    assert np.all(hpx == 0)
    # RA in [0, 360), Dec in [-90, 90]
    assert ra.min() >= 0.0 and ra.max() < 360.0
    assert dec.min() >= -90.0 and dec.max() <= 90.0


def test_get_4hs_mask():
    ra = np.array([10.0, 10.0, 200.0, 200.0])
    dec = np.array([-30.0, 30.0, -40.0, -5.0])
    mask = footprint.get_4hs_mask(ra, dec)
    # only southern points far from the galactic plane (|b| > 20) are kept
    assert mask[1] == False  # northern
    assert mask.dtype == bool
    assert mask.shape == ra.shape


def test_get_4hs_mask_rejects_northern():
    ra = np.array([45.0, 45.0])
    dec = np.array([-50.0, 10.0])
    mask = footprint.get_4hs_mask(ra, dec)
    assert mask[1] == False  # dec > 0 is always rejected
