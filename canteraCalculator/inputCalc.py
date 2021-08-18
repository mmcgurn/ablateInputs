from math import sqrt

import pandas as pd
import numpy as np

import time

# from https://cantera.org/examples/jupyter/reactors/batch_reactor_ignition_delay_NTC.ipynb.html

import cantera as ct
print('Runnning Cantera version: ' + ct.__version__)

# matplotlib notebook
import matplotlib.pyplot as plt

plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['figure.autolayout'] = True

plt.style.use('ggplot')
plt.style.use('seaborn-pastel')

# Load in the cti mech file
gas = ct.Solution('../ignitionDelayGriMech/grimech30.cti')

# Define inlet density
inletTemperature = 300  # Kelvin
inletPressure = 101325  # Pascals
inletYi = {'N2': 0.769306931, 'O2': 0.230693069}
gas.TPY = inletTemperature, inletPressure, inletYi

inletDensity = gas.density

# add in the blowing conditions
blowTemperature = 1500  # Kelvin
blowPressure = 101325  # Pascals
blowYi = {'CH4': 1.0}
gas.TPY = blowTemperature, blowPressure, blowYi

blowDensity = gas.density

# compute the outlet density
outletTemperature = 300  # Kelvin
outletPressure = 101325  # Pascals
outletYi = {'N2': 0.769306931, 'O2': 0.230693069}
gas.TPY = outletTemperature, outletPressure, outletYi

outletDensity = gas.density

# compute the mass coming in
mass = 0.0;
inletArea = 0.0254**2 # 3D
# inletArea = 0.0276 #2D
mass += inletDensity * 5.0 * inletArea
slabArea = 0.00762*sqrt((.036501-0.025)**2 + (.01146- 0.0)**2) #3D
# slabArea = sqrt((.036501-0.025)**2 + (.01146- 0.0)**2) #2D
mass += blowDensity * sqrt(5**2 + 5**2) * slabArea

outletVelocity = mass /(inletArea*outletDensity)
print("OutletVelocity", outletVelocity)