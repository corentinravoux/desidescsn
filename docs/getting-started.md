# Getting started

## Installation

```bash
git clone https://github.com/corentinravoux/desidescsn.git
cd desidescsn
pip install .            # or: pip install -e ".[dev]" for tests
```

Core dependencies (`pandas`, `numpy`, `scipy`, `matplotlib`, `astropy`) install
automatically. Some modules need extra packages that are **not** required to use
the rest of the library:

- [`footprint`](api/footprint.md) uses `healpy`, `fitsio` and `desimodel`;
- [`surveys`](api/surveys.md) uses `fitsio` and `pyyaml`;
- [`catalogs`](api/catalogs.md) uses `GCRCatalogs` (DESC) to read mock catalogs.

## The redshift distribution model

The parametric `n(z)` used throughout lives in [`efficiency`](api/efficiency.md):

```python
import numpy as np
from desidescsn import efficiency

z = np.linspace(0.01, 1.5, 100)
n_z = efficiency.n_z_function(z, A=1.0, z0=0.3, beta=2.0, d=2.0)

# recover the parameters from a measured n(z)
popt = efficiency.fit_n_z(z, n_z)
```

## Applying a survey's target selection

Each survey exposes a `mask_magnitude_*` selection and a
`host_efficiency_*_simple` constant. Band names are mapped to catalog columns
with a `hashing_table`:

```python
import pandas as pd
from desidescsn import surveys

galaxies = pd.DataFrame({"mag_r": [19.0, 20.5, 21.0]})
hashing_table = {"r": "mag_r"}

mask = surveys.mask_magnitude_DESIBGS(galaxies, hashing_table)
eff = surveys.host_efficiency_DESIBGS_simple()   # fibre-assign x z-success
```

## SN population weights and efficiency

```python
from desidescsn import efficiency

weights = efficiency.return_weight_model(
    host_catalog, model="snia", parameters={"A": 1.0, "B": 1.0},
)

redshift, sn_eff = efficiency.compute_sn_efficiency(sn_host_sample, magnitude_mask)
```

## Sky footprints

```python
from desidescsn import footprint

ra, dec, hpx = footprint.create_healpix_map(nside=64)
mask_4hs = footprint.get_4hs_mask(ra, dec)   # southern, |b| > 20 deg
```
