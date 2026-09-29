# Finite-sample bias in ordinal complexity–entropy planes

Code repository for the manuscript:

**A. M. S. Borin Jr., D. Sierra-Porta, R. H. Rosa, R. R. Ferreira**
*Finite-sample bias in ordinal complexity–entropy planes: exact identities under overlapping patterns*
(submitted to Physica Scripta)

## Structure

```
finite_sample_ordinal_complexity_reproduction.ipynb   main reproducibility notebook
figures/                                              all paper figures (PDF + PNG)
real_data_applications/                               general-delay / general-dimension extensions
                                                       and the two real-data applications
```

## Main notebook

`finite_sample_ordinal_complexity_reproduction.ipynb` reproduces the exact
finite-sample identities, the Monte Carlo bias/RMSE study, the bootstrap
coverage study, and the calibrated rank test. Running it end to end (`Run
All`) regenerates:

| Output | Paper reference |
|---|---|
| Table 1 | Section 2.3, exact overlap terms $g(h)$ |
| `figure_1.pdf` (Fig. 1) | Section 3.1, sign reversal of the distinct-pairs bias |
| `figure_2.pdf` (Fig. 2) | Section 3.2, bias–variance trade-offs |
| `figure_3.pdf` (Fig. 3) | Section 3.3–3.4, gap estimator residual and ties |
| `figure_4.pdf` (Fig. 4) | Section 3.5, bootstrap coverage |
| `figure_5.pdf` (Fig. 5) | Section 3.6, calibrated rank test |

## Real-data applications and extensions (`real_data_applications/`)

Contributed by David Sierra-Porta: the general-delay overlap identities, the
general-dimension AR(1) reference, and the two real-data applications
(heliospheric neutron-monitor counts during two Forbush decreases, and
cardiac RR intervals). This module is independent of the main notebook code
(it only reuses the same estimator definitions) and produces:

| Output | Paper reference |
|---|---|
| Table 2 (`validate_extensions.py`, stdout) | Section 2.7–2.8, validation of the two extensions |
| `figure_forbush.pdf` (Fig. 6), Table 3, `forbush_summary.csv` | Section 3.7, heliospheric application |
| `figure_rr.pdf` (Fig. 7), `rr_summary.csv` | Section 3.8, physiological application |

See `real_data_applications/README.md` for how to run each script and for
the data provenance (NMDB, PhysioNet/MIT-BIH).

## Requirements

`numpy`, `scipy`, `pandas`, `matplotlib`; `real_data_applications/rr_extract.py`
additionally needs `wfdb`. No other dependencies.

## Data availability

No proprietary data. The heliospheric records are public exports from the
NMDB database (Rome and Oulu neutron-monitor stations); the physiological
record is public record 14134 of the MIT-BIH Long-Term ECG Database on
PhysioNet. Both are redistributed here under `real_data_applications/RealData/`
together with the scripts that reproduce them from source.
