import numpy as np
import matplotlib.pyplot as plt
from radis.lbl.factory import SpectrumFactory
from radis.spectrum.rescale import rescale_path_length
from scipy.optimize import brentq
from scipy.optimize import curve_fit
import pandas as pd

Tgas_range = np.arange(300,3000+1,100)
delta_wn = 5
wn_range = np.arange(1600,6425+1, delta_wn)

print(len(Tgas_range))
print(len(wn_range))

def init_spectrum(wmed):
  sf = SpectrumFactory(wmed-(delta_wn/2),wmed+(delta_wn/2),
  molecule ='CO',
  isotope='all',
  pressure=1.01325,
  mole_fraction=0.01,
  wstep='auto',
  cutoff=2.294e-166,
  broadening_method='voigt',
  diluent='air',
  truncation=500,
  medium='air',
  verbose=0
  #//warnings={'HighTemperatureWarning': False}
  )

  sf.load_databank(
  path='/local/discoF/u5012198/radisinfo/05_HITEMP2019_all.par',
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
  return (integral_k / (delta_wn)) * (1 / 0.01)

def find_pathlength_for_transmittance(s,target_y):
	def func(x): return transmittance(s,x) - target_y
	x_min, x_max = 1e-4, 10.0
	while func(x_max) > 0: x_max *= 2
	while func(x_min) < 0: x_min /= 2
	return brentq(func, x_min, x_max, xtol=0.01)

def tau_malkmus_co(l, k, T, deltaCO):  # p en atm
    gammaCO =  0.06 * (296 / T) ** 0.7
    tau = np.exp(  -(2 * gammaCO / deltaCO) * (np.sqrt(1 + (0.01 * 1 * k * l * deltaCO / gammaCO)) - 1)  )
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
      lambda l, deltaCO: tau_malkmus_co(l, k, T, deltaCO),
      l_values,
      tau_values,
      p0=[0.001],
      bounds=(0, 100000)
    )
    delta_opt_results[i, j] = 1 / popt[0]
    k_results[i, j] = k



k_df = pd.DataFrame(k_results, index=wn_range, columns=Tgas_range)
delta_opt_df = pd.DataFrame(delta_opt_results, index=wn_range, columns=Tgas_range)

k_df.index.name = 'Wavenumber (cm$^{-1}$)'
k_df.columns.name = 'Gas Temperature (K)'

delta_opt_df.index.name = 'Wavenumber (cm$^{-1}$)'
delta_opt_df.columns.name = 'Gas Temperature (K)'

# Save the DataFrames to CSV files
k_df.to_csv('cok_results.csv')
delta_opt_df.to_csv('codelta_opt_results.csv')
