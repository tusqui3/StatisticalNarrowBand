# Narrow Band Radiation Model

A python implementation to compute narrow band parameters for thermal radiation in gas mixtures using the Malkmus model for gaseous species in high temperature applications such as combustion chambers. Based on RADIS and the following reference,

 *Updated band model parameters for H2O, CO2, CH4 and CO radiation at high temperature*, Philippe Rivière, Anouar Soufiani, 2012

The code returns the statistical narrow band parameter table for a single species as a function of wavenumber band and temperature.

Outputs:

| Parameter | Units | Description |
|-----------|-------|-------------|
| `k` | `cm⁻¹ bar⁻¹` | Band-mean absorption coefficient per unit partial pressure |
| `delta` | `cm⁻¹` | Optimal mean line spacing (Malkmus parameter) |

Results are written to a plain-text file with the following structure as in the paper
- First block (one row per wavenumber band): `wavenumber  k(T1)  k(T2)  ...  k(Tn)`
- Second block (one row per wavenumber band): `wavenumber  delta(T1)  delta(T2)  ...  delta(Tn)`

## summary

The RTE is monochromatic, so requires integrating over all frequencies, however each local κ_ν in a real gas is a the contribution of millions broadened lines, making full line-by-line integration expensive. Gas mixtures add further complexity because species lines overlap non-additively. The narrow band approach models this by dividing the spectrum into windows where line statistics are constant and fitting the Malkmus equation to LBL reference data, giving a analytic transmittance per band.

the method works by least-squares fitting of the Malkmus equation to transmittances computed from line-by-line (LBL) spectra — generated here via RADIS — to find the parameters that minimize the error with respect to the LBL reference. The narrow-band model is built around computing transmittances along a path, which implicitly assumes a purely absorbing/emitting media, that does not allow scattering.

## input parameters

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

## use
```python
from main import main

# Configure parameters inside main() then run:
main()
```
