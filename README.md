# Finite-sample bias in ordinal complexity–entropy planes

Reproducibility repository for the manuscript:

**Finite-sample bias in ordinal complexity–entropy planes: exact identities under overlapping patterns**

Airton M. S. Borin Jr., Robson H. Rosa, and Rafael R. Ferreira  
Federal Institute of Triângulo Mineiro (IFTM), Uberaba, MG, Brazil

## Overview

This repository contains the computational material required to reproduce the analytical benchmarks, Monte Carlo experiments, and figures reported in the manuscript.

The main result concerns finite-sample bias in ordinal complexity–entropy analysis when ordinal patterns are extracted from overlapping windows. The notebook implements the exact overlap identities, reference-process calculations, simulation experiments, bootstrap coverage study, tie/quantization experiment, and Monte Carlo rank test used in the paper.

No empirical data are used.

## Repository structure

```text
finite-sample-ordinal-complexity/
├── README.md
├── requirements.txt
├── finite_sample_ordinal_complexity_reproduction.ipynb
└── figures/
    ├── figure_1.png
    ├── figure_2.png
    ├── figure_3.png
    ├── figure_4.png
    └── figure_5.png
```

The figures directory is created automatically by the notebook if it does not already exist.

## Reproducibility

The notebook is self-contained: no precomputed simulation tables are required. All numerical results and figures are regenerated directly from the implemented analytical calculations and simulations using fixed random seeds.

Two execution modes are available:

```python
FULL_REPRODUCTION = False
```

runs a reduced version for quick inspection, whereas

```python
FULL_REPRODUCTION = True
```

uses the manuscript-scale simulation settings.

For the full reproduction, the principal settings are:

- 5000 independent replications for bias and RMSE calculations;
- 500 outer replications for the bootstrap coverage experiment;
- 499 bootstrap resamples per outer replication;
- 9999 null simulations for calibration of the Monte Carlo rank test;
- raw lengths \(N = 50, 100, 200, 500, 1000, 2000\).

## Running the notebook

Create a Python environment, install the dependencies,

```bash
pip install -r requirements.txt
```

start Jupyter,

```bash
jupyter notebook
```

and open:

```text
finite_sample_ordinal_complexity_reproduction.ipynb
```

Set `FULL_REPRODUCTION = True` for the manuscript-scale calculation and use **Run All**.

The full calculation may take several minutes depending on the machine.

## Numerical audit

The final section of the notebook performs automatic checks against numerical values reported in the manuscript. Analytical targets and Monte Carlo targets are treated separately because independent Monte Carlo realizations are not expected to reproduce identical final digits.

The audit covers, among other quantities:

- exact IID overlap collision terms \(g(h)\) for \(m=3,4,5\);
- analytical logistic-map entropies and allowed-pattern counts;
- analytical AR(1) disequilibria;
- exact finite-sample plug-in and distinct-pairs biases;
- entropy bias and the process-specific second-order correction;
- quantization/tie effects;
- bootstrap coverage;
- size and power of the Monte Carlo rank test;
- the 89 structural zeros for the logistic map at \(m=5\).

## Reference processes

The computational experiments include:

- IID Gaussian white noise;
- stationary Gaussian AR(1) processes;
- non-overlapping white-noise ordinal patterns where required;
- the logistic map at parameter 4.

The AR(1) simulations are initialized from their stationary Gaussian distribution.

For the logistic map, analytical reference probabilities are used. In finite-precision floating-point simulation, exceptionally rare trajectories can develop exact repeated values. The notebook detects such trajectories and regenerates them so that numerical floating-point ties do not alter the continuous-map ordinal experiment. A dedicated audit checks the generated trajectories for this condition.

## Figures

Running the notebook regenerates Figures 1–5 in the local `figures/` directory. The figures supplied in this repository correspond to the manuscript-scale calculations.

The plotting code uses English labels and publication-oriented formatting. Curves are distinguished by line style and marker as well as color.

## Software

The reproduction was developed and audited with Python 3.11 and the packages listed in `requirements.txt`.

The notebook uses only standard numerical/scientific Python dependencies; no external datasets or proprietary software are required.

## Data availability

No empirical data were used in this study. The analytical calculations, simulation code, fixed random seeds, numerical checks, and figure-generation workflow are contained in the reproducibility notebook.

A permanent archival DOI can be added here after the repository is released and archived.

## Citation

If you use this code, please cite the associated manuscript. The final journal citation and archival DOI will be added after publication/repository archiving.

## License

Unless a separate license file is supplied with the repository, the manuscript and source code remain subject to their respective authors' rights. A repository license should be selected before public release if reuse permissions are to be granted explicitly.
