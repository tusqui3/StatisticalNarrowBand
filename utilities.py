import numpy as np
from scipy.optimize import brentq
from radis.lbl.factory import SpectrumFactory


def init_spectrum(
    molecule,
    wavenumberCenter,
    wavenumberRange,
    isotope,
    pressure,
    moleFraction,
    wavenumberStep,
    cutoff,
    broadeningMethod,
    diluent,
    truncation,
    medium,
    databankPath,
    databankType,
    verbose=0
):
    """Initialize a SpectrumFactory object for given parameters."""
    sf = SpectrumFactory(
        wavenumberCenter - (wavenumberRange / 2),
        wavenumberCenter + (wavenumberRange / 2),
        molecule=molecule,
        isotope=isotope,
        pressure=pressure,
        mole_fraction=moleFraction,
        wstep=wavenumberStep,
        cutoff=cutoff,
        broadening_method=broadeningMethod,
        diluent=diluent,
        truncation=truncation,
        medium=medium,
        verbose=verbose,
        optimization=None
    )
    try:
        sf.load_databank(
        path=databankPath,
        format=databankType,
        db_use_cached=True
    )
    except Exception:
        # First run — download and cache it
        sf.fetch_databank(databankType)
    return sf


def transmittance(spectrum, x, delta_wn):
    """Calculate transmittance for given path length x (cm)."""
    spectrum.rescale_path_length(x)
    transmittance_noslit = spectrum.get('transmittance_noslit')
    
    if len(transmittance_noslit[0]) > 0:
        integral_transmittance = np.trapezoid(transmittance_noslit[1], transmittance_noslit[0])
        t = integral_transmittance / delta_wn
        return max(0.0, min(1.0, t))
    else:
        print("No lines found, returning transmittance=1")
        return 1.0


def abs_coeff(spectrum, delta_wn, molar_fraction, pressure):
    """Calculate absorption coefficient. Pressure in atm."""
    absorbance = spectrum.get("abscoeff")
    integral_k = np.trapezoid(absorbance[1], absorbance[0])
    return (integral_k / delta_wn) * (1 / (molar_fraction * pressure))


def find_pathlength_for_transmittance(spectrum, targetTransmittance, delta_wn):
    """Find path length that produces target transmittance."""
    def func(x):
        return transmittance(spectrum, x, delta_wn) - targetTransmittance
    
    x_min, x_max = 1e-4, 10.0
    while func(x_max) > 0:
        x_max *= 2
    while func(x_min) < 0:
        x_min /= 2
    return brentq(func, x_min, x_max, xtol=0.01)


def tau_malkmus(l, k, pressure, delta, moleFraction, gamma):
    """
    Calculate transmittance using Malkmus model.
    
    Parameters:
    l: path length (cm)
    pressure: pressure (atm)
    delta: line spacing (cm-1)
    """
    tau = np.exp(
        -(2 * gamma / delta) * (np.sqrt(1 + (moleFraction * pressure * k * l * delta / gamma)) - 1)
    )
    return tau

def gamma_co(temperature, pressure, molar_fraction):
    """Temperature and pressure dependent gamma for CO."""
    Tref = 296
    Pref = 1.0
    return ((Tref / temperature)**0.7) * (pressure / Pref) * 0.06

def gamma_co2(temperature, pressure, molar_fraction):
    """Temperature and pressure dependent gamma for CO2."""
    Tref = 296
    Pref = 1.0
    return ((Tref / temperature)**0.7) * (pressure / Pref) * (0.07*molar_fraction + (1-molar_fraction)*0.058)

def gamma_h2o(temperature, pressure, molar_fraction):
    """Temperature and pressure dependent gamma for H2O."""
    Tref = 296
    Pref = 1.0
    return (pressure / Pref) * (0.462*molar_fraction*((Tref / temperature)) + ((Tref/ temperature)**0.5) * 0.0792)

def gamma_ch4(temperature, pressure, molar_fraction):
    """Temperature and pressure dependent gamma for CH4."""
    Tref = 296
    Pref = 1.0
    return (pressure / Pref) * 0.051* (Tref/ temperature)**0.75