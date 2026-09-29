# Real-data applications and general-delay / general-dimension extensions

Code for the parts of the paper contributed by David Sierra-Porta: the exact
overlap identities for a general embedding delay τ, the AR(1) reference for a
general embedding dimension m, and the two real-data applications (a
heliospheric record and a physiological record). Everything here builds on
the estimators defined in the main notebook (`../finite_sample_ordinal_complexity_reproduction.ipynb`) —
same permutation encoding, same H/D/distinct-pairs/Jensen–Shannon definitions.

## Requirements

`numpy`, `scipy`, `pandas`, `matplotlib`. `rr_extract.py` additionally needs
`wfdb` (only for that one script; not required to reproduce the figures from
the already-extracted `RealData/RR_14134.csv`).

## Files and what they reproduce

| File | Produces | Paper reference |
|---|---|---|
| `fsb_extensions.py` | Core math: `g_tau(h)` overlap identity, AR(1) orthant-probability reference for general m | Section 2.7–2.8, Eqs. (12)–(16) |
| `apply_diagnostics.py` | Shared estimator/diagnostic functions (symbol encoding, H/D/U/gap, rank test) used by the three scripts below | — (library only, not run directly) |
| `validate_extensions.py` | Numerical validation printed to stdout: τ>1 white-noise bias (exact vs Monte Carlo) and AR(1) orthant vs direct simulation, for m=3,4,5 | Table 2 |
| `forbush_analysis.py` | `figure_forbush.(png\|pdf)` (copied into `../figures/`), `forbush_summary.csv` | Figure 6, Table 3, Section 3.7 |
| `rr_extract.py` | `RealData/RR_14134.csv` from the WFDB source (`RealData/14134.{dat,hea,atr}`) | provenance for Figure 7 / Section 3.8 |
| `rr_analysis.py` | `figure_rr.(png\|pdf)` (copied into `../figures/`), `rr_summary.csv` | Figure 7, Section 3.8 |

## Data

`RealData/` holds the inputs:
- `FD240510_*.txt`, `FD240324_*.txt`: NMDB neutron-monitor exports (Rome and Oulu stations) for the two Forbush decreases analysed in the paper.
- `14134.dat`, `14134.hea`, `14134.atr`: WFDB source files for MIT-BIH Long-Term ECG Database record 14134 (PhysioNet), the physical provenance of the RR series.
- `RR_14134.csv`: the extracted RR-interval series (one column, integer samples), regenerated from the WFDB source by `rr_extract.py`.

## Running

```bash
pip install numpy scipy pandas matplotlib wfdb

python3 validate_extensions.py      # Table 2 numbers, printed to stdout
python3 forbush_analysis.py         # figure_forbush.{png,pdf}, forbush_summary.csv
python3 rr_extract.py               # regenerates RealData/RR_14134.csv (optional; already provided)
python3 rr_analysis.py              # figure_rr.{png,pdf}, rr_summary.csv
```

All four scripts are self-contained (no absolute paths; each resolves its own
directory at runtime) and run in a few seconds to a couple of minutes each.
`apply_diagnostics.py` is a library imported by the others and is not run on
its own.
