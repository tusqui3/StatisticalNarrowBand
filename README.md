# Narrow Band Radiation Model

A Python implementation for computing narrow band radiation model in gas mixtures using the Malkmus model for gaseous species in high temperature applications. Based on Radis and *International Journal of Heat and Mass Transfer*, 2012, DOI: 10.1016/j.ijheatmasstransfer.2012.03.019. The code returns the statistical narrow band parameters as a function of the wavenumber band and the temperature, for each specie. 

## Overview

Statistical Narrow band is a simplified model for modelling the transmission coefficients through gas mixtures in a non-gray basis, band-by-band. It is considered one of the most precise models, although is very difficult to introduce scattering when using the Narrow Band approach. The method is just based on a least squares fitting of the malkmus equation to find the parameters thaat minimize the error compared with line by line calculations. 

## Inputs

Given a certain specie, the following parameters have to be specified:

- Supports virtually any species, althoy

## Basic Usage

```python
from build_narrow_band import main

# Run with default parameters
main()
```

## Output

The script generates the following files:

- `k_results.npy`: Absorption coefficients [wavenumber × temperature]
- `delta_opt_results.npy`: Optimal delta values [wavenumber × temperature]
- `wn_range.npy`: Wavenumber range used
- `Tgas_range.npy`: Temperature range used
