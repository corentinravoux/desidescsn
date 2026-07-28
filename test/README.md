# desidescsn test suite

Run from the package root:

```bash
pip install -e ".[dev]"   # installs pytest alongside the runtime deps
pytest -q                 # or: pytest test/
```

The suite exercises the **pure / offline** functions against synthetic
dataframes — no external catalogs, survey files or network access are required.

| File | Scope | What it covers |
|---|---|---|
| `test_efficiency.py` | unit | the parametric `n_z_function` and its `fit_n_z` recovery; the `return_weight_model` SN weight models (`noweigths`/`snia`/`sncc`, incl. the sSFR cut); `compute_sn_efficiency`; and the `compute_host_efficiency` ratio callable (with cap) |
| `test_surveys.py` | unit | the constant `host_efficiency_*_simple` products; and the `mask_magnitude_*` magnitude/colour selections for DESIBGS, 4HS and CRSLRG |
| `test_footprint.py` | unit | `create_healpix_map` (NESTED, correct npix and RA/Dec ranges) and the `get_4hs_mask` southern + galactic-latitude cut |
| `test_imports_smoke.py` | smoke | every shipped module byte-compiles; `surveys`/`efficiency`/`footprint` import; `catalogs` imports when `GCRCatalogs` is present, else skips |

## Note on optional dependencies

- **`catalogs`** imports `GCRCatalogs` (a heavy DESC dependency) at module load;
  its import test `skip`s gracefully when that package is absent, while its
  syntax is still checked.
- The file-reading functions (`n_z_*` loaders, footprint mask builders that read
  FITS/HEALPix files) need survey data products and are documented but not
  exercised here.
