
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
gas = ct.Solution('../mechanisms/gri30.yaml')

# Define inlet density
inletTemperature = 300  # Kelvin
inletPressure = 101325  # Pascals
inletYi = {'O2': 1.0}
gas.TPY = inletTemperature, inletPressure, inletYi
inletDensity = gas.density

# compute the mass coming in
gAtFuel = 22.19 # kg/m2/s

# set areas
inletArea = 0.027686*(0.0127*2)
fuelArea = 0.0114*(0.00381*2)

# fuel normal area
fuelFlowArea = inletArea - fuelArea;

# compute G at inlet
gAtInlet = gAtFuel*fuelFlowArea/inletArea
inletVel = gAtInlet/inletDensity

print("O2Density", inletDensity)
print("velocityAtFuel", gAtFuel/inletDensity)
print("inletVel", inletVel)
