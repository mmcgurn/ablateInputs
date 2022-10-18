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
gas = ct.Solution('../mechanisms/grimech30.cti')

# Calculate simple fuel decomposition (paraffin wax), and oxidizer stoichiometric mass coefficient
# Paraffin wax is CnH_(2n+2), taking n = 31...
n = 31
x = n
y = (2*n+2)
nCH = x
nH2 = (y-x)/2.
ntotGas = nCH + nH2

# Calculate the equivalent mole Fractions of CH and H2 in the wax
XCH = nCH/ntotGas
XH2 = nH2/ntotGas

# Set the surface params
Tboil = 643.15

gas.TPX = Tboil, ct.one_atm, {"CH": XCH, "H2": XH2}
gas.equilibrate('TP')

print("equilibrate State:")
print("rho: ", gas.density)
print("T: ", gas.T)
print("p: ", gas.P)
print("u + hf: ", gas.u)
print("rho*yi: ", np.array2string(gas.Y * gas.density, separator=',', max_line_width=100000))
print("yi: ", np.array2string(gas.Y, separator=',', max_line_width=100000))
print("xi: ", np.array2string(gas.X, separator=',', max_line_width=100000))
print("basis: ", gas.basis)
# Print mass fractions if they are meaningful
for (species, yi) in zip(gas.species_names, gas.Y):
    # only print the sigificant values
    if yi > 1E-8:
        print(species , ":", yi)
