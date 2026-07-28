# desidescsn

Forecast how efficiently spectroscopic surveys recover the **host-galaxy
redshifts of type-Ia supernovae** — the key ingredient for combining **DESI**
(and DESI2 / 4MOST) spectroscopy with **DESC / LSST** supernova samples.

📖 **Full documentation (API reference + user guide):** **<https://corentinravoux.github.io/desidescsn/>**

---

## What it does

```
DESC mock catalog ─▶ catalogs    host photometry / z / mass / SFR
                  ─▶ surveys     target-selection cuts + n(z)   ─┐
                  ─▶ efficiency  SN population + magnitude cuts  ─┤
                                                                 ├─▶ full efficiency(z)
                                 host efficiency × SN efficiency ─┘
footprint ─▶ sky masks (DESI / DESI2 / 4MOST CRS / 4HS)
```

## Install

```bash
pip install .              # core deps: pandas, numpy, scipy, matplotlib, astropy
pip install -e ".[dev]"    # + pytest, to run the test suite
```

Some modules need extra packages: `healpy` + `fitsio` + `desimodel`
([`footprint`](desidescsn/footprint.py)), `pyyaml` + `fitsio`
([`surveys`](desidescsn/surveys.py)), and `GCRCatalogs`
([`catalogs`](desidescsn/catalogs.py), DESC mock catalogs).

## Quick start

```python
import numpy as np, pandas as pd
from desidescsn import efficiency, surveys

# parametric redshift distribution
z = np.linspace(0.01, 1.5, 100)
n_z = efficiency.n_z_function(z, A=1.0, z0=0.3, beta=2.0, d=2.0)

# a survey's target selection + constant host efficiency
galaxies = pd.DataFrame({"mag_r": [19.0, 20.5, 21.0]})
mask = surveys.mask_magnitude_DESIBGS(galaxies, {"r": "mag_r"})
eff  = surveys.host_efficiency_DESIBGS_simple()
```

## Package map

| Module | Role |
|---|---|
| `surveys` | Per-survey target-selection cuts, constant host efficiencies and `n(z)` loaders — DESI BGS/LRG, DESI2, 4MOST CRS BG/LRG, 4HS |
| `efficiency` | Parametric `n(z)` model + fit; host and SN redshift efficiencies and their combination |
| `footprint` | Boolean sky-footprint masks (DESI via `desimodel`, others via HEALPix maps / cuts) |
| `catalogs` | Extract host-galaxy properties from a DESC DC2 / Roman-Rubin mock catalog |

## Conventions

- Redshift distributions `n(z)` are per deg² per unit redshift; efficiencies are
  dimensionless fractions.
- Sky coordinates are RA/Dec in **degrees**; HEALPix maps use the **NESTED**
  scheme.
- Band names are mapped to catalog columns through a `hashing_table` dict.

## Documentation

📖 **Online: <https://corentinravoux.github.io/desidescsn/>**

Built from the in-code docstrings with
[MkDocs Material](https://squidfunk.github.io/mkdocs-material/) +
[mkdocstrings](https://mkdocstrings.github.io/) and deployed to GitHub Pages on
every push to `main`. Build locally with:

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

## Tests

```bash
pip install -e ".[dev]"
pytest -q
```

See [`test/README.md`](test/README.md) for the layout of the suite.
