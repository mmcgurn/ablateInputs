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
referenceTemperature = 298.15
temperature = 350.0  # Kelvin
pressure = 101325  # Pascals

# print the enthalpy

# Output information
print("Species: ", gas.species_names)
print("NumberOfSpecies: ", len(gas.species_names))

# march over each species
for sp in gas.species_names:
    # set the gas to a single species
    gas.Y = {sp: 1.0}
    gas.TP = temperature, pressure
    hi = gas.h

    # Get the ref state
    gas.TP = referenceTemperature, pressure
    hf = gas.h
    print(sp, ": ", (hi-hf))

