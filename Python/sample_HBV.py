# -*- coding: utf-8 -*-
"""

Sample run for the HBV96


In this folder there are two implementations of the HBV model.

The HBV96.py is the pure python implementation.
The HBV96x.pyx is its equivalent the Cython implementation

If you are planning to run several 

Created on Tue Oct 10 16:36:27 2017

@author: chaco3
"""

import numpy as np
import HBV96
import matplotlib.pyplot as plt

# Load the sample data (Prec, evap, T, Q, LAT)
data = np.loadtxt('SampleData.csv', skiprows=1, delimiter=',')

prec = data[:,0]  #Precipitation
et = data[:,1]  # Evapotranspiration
t = data[:,2]  # Temperature
q = data[:,3]  # Discharge
lat = data[:,4]  # Long term average temperature

# get random parameter set
pars = HBV96.get_random_pars()
p2 = [1.0, 142.0]  # TFAC and Area [Km2]

reload(HBV96)
# Run the simulation with the random paramete set
q_sim, st_sim = HBV96.simulate(prec, t, et, pars, p2, ll_temp=lat)

# Calibrate the model using simplified algorithm
pars, perf = HBV96.calibrate(q, prec, t, et, p2, wu=50)
q_sim, st_sim = HBV96.simulate(prec, t, et, pars, p2, ll_temp=lat, q_0=2.0)

plt.plot(q_sim, label='Sim')
plt.plot(q, label='Rec')
plt.legend()
