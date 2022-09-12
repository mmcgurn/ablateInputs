
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

# add in the blowing conditions
blowTemperature = 1500  # Kelvin
blowPressure = 101325  # Pascals
blowYi = {'CH4': 1.0}
gas.TPY = blowTemperature, blowPressure, blowYi
blowDensity = gas.density

# compute the mass coming in
inletG = 22.19 # kg/m2/s
inletArea = 0.027686*(0.0127*2)

inletVel = inletG/inletDensity
inletMassRate = inletG*inletArea # kg/(m2 s) * m2 = kg/s

# compute stoic ratio
o2FuelRatio = 3.989

# compute the fuelInlet
fuelRate = inletMassRate/o2FuelRatio # kg/s
fuelArea = 0.06858*(0.00381*2) #m2
fuelVelocity = fuelRate/(blowDensity*fuelArea)

print("inletVelocity", inletVel)
print("fuelVelocity", fuelVelocity)