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

# Define the reactor temperature and pressure
density = 1000
pressure = 101325
yi = {"H2": 0.0, "H2O": 0.0, "N2": .0, "CO": 0.0, "CH4": 0.0, "NO": 0, "CH2": 0, "O2": 1.0}

# compute the reference enthalpy
gas.TPY = 298.15, 101325.0, yi
enthalpyOfFormation = gas.h

# # Set the real t and p
# gas.DPY = density, pressure, yi

# print the Cp
print("rho: ", gas.density)
print("Cp: ", gas.cp_mass)
print("Cv: ", gas.cv_mass)
print("u: ", gas.u - enthalpyOfFormation)
print("h: ", gas.h - enthalpyOfFormation)

