from sympy import *
import cantera as ct
print('Runnning Cantera version: ' + ct.__version__)

# matplotlib notebook
import matplotlib.pyplot as plt

# determine the density of gas at the inlet
# # Load in the cti mech file
gas = ct.Solution('../mechanisms/gri30.yaml')

# Define inlet density
inletTemperature = 300  # Kelvin
inletPressure = 101325  # Pascals
inletYi = {'O2': 1.0}
gas.TPY = inletTemperature, inletPressure, inletYi
inletDensity = gas.density

# define the geometry
yo = 0.0;
yh = 0.0254
dia = (yh - yo)
yc = (yo + yh) * 0.5
depth = 0.0254
inletMassFlowKgMin = 0.192 #kg/min
inletMassFlow = inletMassFlowKgMin/60.0

# setup symbols
y = Symbol('y')
f = Symbol('f')
m = Symbol('m')

# compute the mass flux function
# for laminar (https://www.kau.edu.sa/Files/0057863/Subjects/Chapter%208.pdf 8-17)
# velocityProfile = f*(1-((y-yc)**2)/((dia**2)/4))
# for turbulent
velocityProfile = f*(1-sqrt((y-yc)**2)/(dia/2))**(1/7)

# for turblent flow
massFluxProfile = velocityProfile * inletDensity

# Perform a simple 1D integration (assume that z dimension is one)
massFlow = integrate(massFluxProfile, (y, yo, yh))*depth

# Solver for f
conserved = Eq(massFlow, m)
fSolution = solve(conserved, f)[0]
print("VelocityFactor: ", fSolution)
print("VelocityFactor for ", inletMassFlow, " kg/s is ",  fSolution.subs(m, inletMassFlow))

#
# Plot the velocity
fSolution = fSolution.subs(m,inletMassFlow)
velocityProfile = velocityProfile.subs(f,fSolution)
plot(velocityProfile, (y, yo, yh))
