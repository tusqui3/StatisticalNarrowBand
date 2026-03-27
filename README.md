# Narrow Band Radiation Model

A python implementation to compute narrow band parameters for thermal radiation in gas mixtures using the Malkmus model for gaseous species in high temperature applications. Based on RADIS and the following reference,

 *Updated band model parameters for H2O, CO2, CH4 and CO radiation at high temperature*, Philippe Rivière, Anouar Soufiani, 2012

The code returns the statistical narrow band parameter table for a single species as a function of wavenumber band and temperature.

**Outputs:**

| Parameter | Units | Description |
|-----------|-------|-------------|
| `k` | `cm⁻¹ bar⁻¹` | Band-mean absorption coefficient per unit partial pressure |
| `delta` | `cm⁻¹` | Optimal mean line spacing (Malkmus parameter) |

Results are written to a plain-text file with the following structure as in the paper
- First block (one row per wavenumber band): `wavenumber  k(T1)  k(T2)  ...  k(Tn)`
- Second block (one row per wavenumber band): `wavenumber  delta(T1)  delta(T2)  ...  delta(Tn)`

## Overview

The Statistical narrow band model is a simplified approach for computing spectral transmission coefficients through gas mixtures on a band-by-band basis. It is considered one of the most accurate engineering radiation models, though it is difficult to incorporate scattering within the narrow band framework and is expensive in comparison with other gray models. thee method works by least-squares fitting of the Malkmus equation to transmittances computed from line-by-line (LBL) spectra — generated here via RADIS — to find the parameters that minimize the error with respect to the LBL reference.

## Inputs

The following parameters must be specified inside the `main()` function:

| Parameter | Description | Example |
|-----------|-------------|---------|
| `molecule` | Target species in capital letters | `'CO'`, `'CO2'`, `'H2O'`, `'CH4'` |
| `isotope` | Isotopologue selection | `'all'` |
| `databankPath` | Spectroscopic database to use | `'hitran'`, `'hitemp'` |
| `pressure` | Total pressure (bar) | `1.0` |
| `moleFraction` | estimation of approximately how much of this specie will be present in the mixture | `0.01` |
| `numberOfFittingPoints` | Number of transmittance points used in the Malkmus fit (higher = more accurate) | `20` |
| `delta_wn` | Narrow band width (cm⁻¹) | `25` |
| `wn_range` | Array of band-center wavenumbers (cm⁻¹) | `np.arange(2000, 2100, 25)` |
| `Tgas_range` | Array of temperatures K | `np.arange(200, 700, 50)` |
| `wavenumberStep` | Spectral resolution for LBL calculation (cm⁻¹) | `0.01` |
| `cutoff` | Line strength cutoff threshold (cm⁻¹) | `1e-27` |
| `broadeningMethod` | Line profile: `'voigt'`, `'lorentzian'`, or `'gaussian'` | `'voigt'` |
| `truncation` | Line wing truncation distance (cm⁻¹) | `50` |

> **Note:**  averaged Lorentz half-widths (`γ`) are implemented for `H2O`, `CO2`, `CH4`, and `CO`. The code automatically switch between these, To use any other species, a custom `gamma_*` function must be added to `utilities.py`.

## Basic Usage
```python
from main import main

# Configure parameters inside main() then run:
main()
```
