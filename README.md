# Narrow Band Radiation Model

A Python implementation for computing narrow band radiation model parameters using the Malkmus model for gaseous species.

## Overview

This code computes absorption coefficients and optimal delta values for narrow band radiation modeling. It uses the RADIS library for spectral calculations and fits the Malkmus model to obtain temperature and pressure dependent parameters.

## Features

- Supports multiple molecular species (CO, CO2, H2O)
- Temperature and pressure dependent gamma functions
- Malkmus model parameter fitting
- Efficient batch processing across wavenumber and temperature ranges
- Results saved in NumPy format for further analysis

## Requirements

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```python
from build_narrow_band import main

# Run with default parameters
main()
```

### Custom Configuration

Edit the configuration parameters in `build_narrow_band.py`:

```python
molecule = 'CO'  # or 'CO2', 'H2O'
pressure = 1.0   # atm
moleFraction = 0.1
wn_range = np.arange(2000, 2100, 25)  # cm-1
Tgas_range = np.arange(200, 1000 + 1, 50)  # K
```

### Using Different Species

For different molecular species, update the gamma function in the main loop:

```python
# For CO2
from utilities import gamma_co2
gamma_val = gamma_co2(T, pressure, moleFraction)

# For H2O
from utilities import gamma_h2o
gamma_val = gamma_h2o(T, pressure, moleFraction)
```

## Output

The script generates the following files:

- `k_results.npy`: Absorption coefficients [wavenumber × temperature]
- `delta_opt_results.npy`: Optimal delta values [wavenumber × temperature]
- `wn_range.npy`: Wavenumber range used
- `Tgas_range.npy`: Temperature range used

## Functions

### Core Functions

- `init_spectrum()`: Initialize RADIS SpectrumFactory
- `transmittance()`: Calculate transmittance for given path length
- `abs_coeff()`: Calculate absorption coefficient
- `tau_malkmus()`: Malkmus model transmittance calculation
- `find_pathlength_for_transmittance()`: Find path length for target transmittance

### Gamma Functions

Temperature and pressure dependent gamma functions for different species:

- `gamma_co()`: Carbon monoxide
- `gamma_co2()`: Carbon dioxide  
- `gamma_h2o()`: Water vapor

## References

Gamma function formulations based on:
*International Journal of Heat and Mass Transfer*, 2012, DOI: 10.1016/j.ijheatmasstransfer.2012.03.019

## License

MIT License