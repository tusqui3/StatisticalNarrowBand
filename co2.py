import numpy as np
import matplotlib.pyplot as plt
from radis.lbl.factory import SpectrumFactory
from radis.spectrum.rescale import rescale_path_length
from scipy.optimize import brentq
from scipy.optimize import curve_fit
import pandas as pd

Tgas_range = np.arange(300,3000+1,100)
delta_wn = 5.0
wn_range = np.arange(250,8300+1, delta_wn)

def init_spectrum(wmed):
  sf = SpectrumFactory(wmed-(delta_wn/2),wmed+(delta_wn/2),
  molecule ='CO2',
  isotope='1,2',
  pressure=1.01325,
  mole_fraction=0.1,
  wstep='auto',
  cutoff=1e-166,
  broadening_method='voigt',
  diluent='air',
  truncation=50,
  medium='air',
  verbose=0
  #//warnings={'HighTemperatureWarning': False}
  )

  sf.load_databank(
  path='/local/discoF/u5012198/hyperion/createNarrowBand/hitran/co2/*',
  format='hitemp',
  db_use_cached=True
  )
  return sf

def transmittance(s,x):
  s.rescale_path_length(x) # beware! x is in cm
  transmittance_noslit = s.get('transmittance_noslit')
  if (len(transmittance_noslit[0])>0):
    integral_transmittance = np.trapz(transmittance_noslit[1], transmittance_noslit[0])
    t = integral_transmittance / (delta_wn)
    return max(0.0, min(1.0, t))
  else:
    print("No lines found, reuturnining transmittance=1")
    return 1.0

def abs_coeff(s):
  absorbance = s.get("abscoeff")
  integral_k = np.trapz(absorbance[1], absorbance[0])
  return (integral_k / (delta_wn)) * (1 / 0.1)

def find_pathlength_for_transmittance(s,target_y):
	def func(x): return transmittance(s,x) - target_y
	x_min, x_max = 1e-4, 10.0
	while func(x_max) > 0: x_max *= 2
	while func(x_min) < 0: x_min /= 2
	return brentq(func, x_min, x_max, xtol=0.01)

def tau_malkmus_co2(l, k, T, deltaCO2):  # p en atm
    gammaCO2 =  (296 / T)**0.7 * (0.07*0.1 + 0.058*(1-0.1))
    tau = np.exp(  -(2 * gammaCO2 / deltaCO2) * (np.sqrt(1 + (0.1 * 1 * k * l * deltaCO2 / gammaCO2)) - 1)  )
    return tau


k_results = np.zeros((len(wn_range), len(Tgas_range)))
delta_opt_results = np.zeros((len(wn_range), len(Tgas_range)))
tau_values = np.zeros(30)
k = np.zeros(30)
for i, wn in enumerate(wn_range):
  sf = init_spectrum(wn)
  print(str(i)+"/"+str(len(wn_range)))
  for j, T in enumerate(Tgas_range):
    s = sf.eq_spectrum(Tgas=T) # CHECK UNITS!
    l1 = find_pathlength_for_transmittance(s,0.02)
    l2 = find_pathlength_for_transmittance(s,0.95)
    l_values = np.linspace(l2, l1, 30)
    for idx, l in enumerate(l_values): tau_values[idx] = transmittance(s,l)
    k = abs_coeff(s)
    popt, _ = curve_fit(
      lambda l, deltaCO2: tau_malkmus_co2(l, k, T, deltaCO2),
      l_values,
      tau_values,
      p0=[0.001],
      bounds=(0, np.inf)
    )
    delta_opt_results[i, j] = 1 / popt[0]
    k_results[i, j] = k



k_df = pd.DataFrame(k_results, index=wn_range, columns=Tgas_range)
delta_opt_df = pd.DataFrame(delta_opt_results, index=wn_range, columns=Tgas_range)

k_df.index.name = 'Wavenumber (cm$^{-1}$)'
k_df.columns.name = 'Gas Temperature (K)'

delta_opt_df.index.name = 'Wavenumber (cm$^{-1}$)'
delta_opt_df.columns.name = 'Gas Temperature (K)'

k_df.to_csv('co2k_results.csv')
delta_opt_df.to_csv('co2delta_opt_results.csv')
