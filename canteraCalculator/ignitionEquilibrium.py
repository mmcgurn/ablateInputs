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
reactor_temperature = 300  # Kelvin
reactor_pressure = 101325  # Pascals

gas.TP = reactor_temperature, reactor_pressure

# Output information
print("Species: ", gas.species_names)
print("NumberOfSpecies: ", len(gas.species_names))
# Define the fuel, oxidizer and set the stoichiometry
gas.set_equivalence_ratio(phi=1.0, fuel="CH4", oxidizer={"o2": 1.0, "n2": 3.76})

# Print the mass fraction of each species
print("Init State:")
print("rho: ", gas.density)
print("T: ", gas.T)
print("p: ", gas.P)
print("u + hf: ", gas.u)
print("rho*yi: ", np.array2string(gas.Y * gas.density, separator=',', max_line_width=100000))
print("yi: ", np.array2string(gas.Y, separator=',', max_line_width=100000))
print("basis: ", gas.basis)

# equilibrate
gas.equilibrate('HP')

# Print the mass fraction of each species
print("equilibrate State:")
print("rho: ", gas.density)
print("T: ", gas.T)
print("p: ", gas.P)
print("u + hf: ", gas.u)
print("rho*yi: ", np.array2string(gas.Y * gas.density, separator=',', max_line_width=100000))
print("yi: ", np.array2string(gas.Y, separator=',', max_line_width=100000))
print("basis: ", gas.basis)
