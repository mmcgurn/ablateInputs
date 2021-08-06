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
density = 999.9
pressure = 281125963.5
yi = {"H2": .1, "H2O": .2, "N2": .3, "CO": .4}
gas.DPY = density, pressure, yi


# print the Cp
print("Cp: ", gas.cp_mass)



