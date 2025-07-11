"""
Build narrow band radiation model parameters for gaseous species.

This script computes absorption coefficients and optimal delta values
for narrow band radiation modeling using the Malkmus model.
"""

import numpy as np
from scipy.optimize import curve_fit
from utilities import (
    init_spectrum, transmittance, abs_coeff, 
    find_pathlength_for_transmittance, tau_malkmus, gamma_co
)

def main():
    # Configuration parameters
    molecule = 'CO'
    isotope = 'all'
    wavenumberStep = 0.01  # cm-1
    cutoff = 25  # cm-1
    broadeningMethod = 'voigt'
    diluent = 'air'
    truncation = 25  # cm-1
    medium = 'air'
    databankPath = 'hitran'
    databankType = 'hitran'
    
    pressure = 1.0
    moleFraction = 0.1
    numberOfFittingPoints = 20
    
    delta_wn = 25  # cm-1
    wn_range = np.arange(2000, 2100, delta_wn)
    Tgas_range = np.arange(200, 1000 + 1, 50)  # K
    
    # Initialize result arrays
    k_results = np.zeros((len(wn_range), len(Tgas_range)))
    delta_opt_results = np.zeros((len(wn_range), len(Tgas_range)))
    tau_values = np.zeros(numberOfFittingPoints)
    
    for i, wn in enumerate(wn_range):
        sf = init_spectrum(
            molecule=molecule,
            wavenumberCenter=wn,
            wavenumberRange=delta_wn * 2,
            isotope=isotope,
            pressure=pressure,
            moleFraction=moleFraction,
            wavenumberStep=wavenumberStep,
            cutoff=cutoff,
            broadeningMethod=broadeningMethod,
            diluent=diluent,
            truncation=truncation,
            medium=medium,
            databankPath=databankPath,
            databankType=databankType
        )
        
        print(f"Processing wavenumber band {i+1}/{len(wn_range)}: {wn} cm-1")
        
        for j, T in enumerate(Tgas_range):
            try:
                s = sf.eq_spectrum(Tgas=T)
                
                l1 = find_pathlength_for_transmittance(s, 0.02, delta_wn)
                l2 = find_pathlength_for_transmittance(s, 0.95, delta_wn)
                l_values = np.linspace(l2, l1, numberOfFittingPoints)
                
                for idx, l in enumerate(l_values):
                    tau_values[idx] = transmittance(s, l, delta_wn)
                
                k = abs_coeff(s, delta_wn, moleFraction, pressure)
                gamma_val = gamma_co(T, pressure, moleFraction)
                
                popt, _ = curve_fit(
                    lambda l, delta: tau_malkmus(l, k, pressure, delta, moleFraction, gamma_val),
                    l_values,
                    tau_values,
                    p0=[0.01],
                    bounds=(0, np.inf)
                )
                
                delta_opt_results[i, j] = popt[0]
                k_results[i, j] = k
                
            except Exception as e:
                print(f"Error at wn={wn}, T={T}: {e}")
                delta_opt_results[i, j] = np.nan
                k_results[i, j] = np.nan
    
    print("Calculation completed!")
    print(f"k_results shape: {k_results.shape}")
    print(f"delta_opt_results shape: {delta_opt_results.shape}")
    
    # Save results
    np.save('k_results.npy', k_results)
    np.save('delta_opt_results.npy', delta_opt_results)
    np.save('wn_range.npy', wn_range)
    np.save('Tgas_range.npy', Tgas_range)
    
    print("Results saved to .npy files")


if __name__ == "__main__":
    main()