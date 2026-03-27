"""
Input parameters can be configured in main() function.

ONLY WORKS WITH H2O, CO2, CH4 and CO.
IF YOU WANT TO USE ANOTHER SPECIE, YOU HAVE TO SPECIFY A customized avg Lorentz line-widths. in utilities.py

OUTPUTS:
k cm⁻¹ bar⁻¹ Band-mean absorption coefficient per unit partial pressure
delta cm⁻¹Optimal mean line spacing (Malkmus parameter)
"""

import numpy as np
from scipy.optimize import curve_fit
from utilities import (
    init_spectrum, transmittance, abs_coeff, 
    find_pathlength_for_transmittance, tau_malkmus, gamma_co, gamma_ch4, gamma_co2, gamma_h2o
)


def build_narrow_band():
    ####### BEGIN INPUT PARAMETERS
    output_file = "results.txt"
    molecule = 'CO' # must use capital letters!
    isotope = 'all'
    databankPath = 'hitran' # hitemp is a must a high temperature applications.
    databankType = 'hitran'
    pressure = 1.0 # bar
    moleFraction = 0.01 #use a guess of approx how much of the current specie will be in your desired mixture approx! For combustion, for examle: 0.1 for h2o and co2 and 0.01 for co.

    numberOfFittingPoints = 20 # higher better but more cost
    delta_wn = 25  # cm-1 bw of the output model

    # DESIRED WN RANGE AND T DISCRETIZATION
    wn_range = np.arange(2000, 2100, delta_wn)
    Tgas_range = np.arange(200.0, 700.0 + 1, 50.0)  # K

    #optimization = 'none' This is hardcoded later in the code because RADIS/Numba missmatch bug does not allow for optimizations at the moment 03/2026
    wavenumberStep = 0.01  # cm-1. using 'auto' sometimes does not work properly. The smaller the better.
    # CUTOFF: discard linestrengths that are lower that this, to reduce calculation times. 1e-27 is what is generally used to generate databases such as CDSD. If 0, no cutoff. Default 1e-27.
    cutoff =  1e-27  # cm-1.
    # BRAOADENING PROFILE. gaussian is better at low pressure(high altitute) lorentzian at atm pressure is good. Voigt is better since is the convolution of both. Its most accurate.
    broadeningMethod = 'voigt'
    diluent = 'air'
    # The higher the better.
    truncation = 50  # cm-1
    medium = 'air' # or 'vacuum' 4 space applications?
    ####### END INPUT PARAMETERS
    
    k_results = np.zeros((len(wn_range), len(Tgas_range)))
    delta_opt_results = np.zeros((len(wn_range), len(Tgas_range)))
    tau_values = np.zeros(numberOfFittingPoints)
    
    for i, wn in enumerate(wn_range):
        sf = init_spectrum(
            molecule=molecule,
            wavenumberCenter=wn,
            wavenumberRange=delta_wn,
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
        pressureAtm = pressure*1.01325
        
        print(f"Processing wavenumber band {i+1}/{len(wn_range)}: {wn} cm-1")
        
        for j, T in enumerate(Tgas_range):
            try:
                s = sf.eq_spectrum(Tgas=T)

                l1 = find_pathlength_for_transmittance(s, 0.02, delta_wn)
                l2 = find_pathlength_for_transmittance(s, 0.95, delta_wn)
                l_values = np.linspace(l2, l1, numberOfFittingPoints)
                
                

                for idx, l in enumerate(l_values):
                    tau_values[idx] = transmittance(s, l, delta_wn)
                
                k = abs_coeff(s, delta_wn, moleFraction, pressureAtm)
                if molecule == 'CO':
                    gamma_val = gamma_co(T, pressureAtm, moleFraction)
                if molecule == 'CO2':
                    gamma_val = gamma_co2(T, pressureAtm, moleFraction)
                if molecule == 'H2O':
                    gamma_val = gamma_h2o(T, pressureAtm, moleFraction)
                if molecule == 'CH4':
                    gamma_val = gamma_ch4(T, pressureAtm, moleFraction)
                
                

                popt, _ = curve_fit(
                    lambda l, delta: tau_malkmus(l, k, pressureAtm, delta, moleFraction, gamma_val),
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
    
    # Save results
    #np.save('k_results.npy', k_results)
    #np.save('delta_opt_results.npy', delta_opt_results)
    #np.save('wn_range.npy', wn_range)
    #np.save('Tgas_range.npy', Tgas_range)

    # Save to custom TXT file
    
    nWn = k_results.shape[0]   # number of wavenumbers
    nT  = k_results.shape[1]   # number of temperatures

    with open(output_file, "w") as f:
        for i in range(nWn):
            row = "  "+f"{float(wn_range[i]):15.10E} " + "  ".join(f"{k_results[i, j]:15.10E}" for j in range(nT))
            f.write(row + "\n")

        for i in range(nWn):
            row = "  "+ f"{float(wn_range[i]):15.10E} " + "  ".join(f"{delta_opt_results[i, j]:15.10E}" for j in range(nT))
            f.write(row + "\n")

    
    print("Results saved to .txt files")


if __name__ == "__main__":
    build_narrow_band()