import sympy
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
inletMassFlowKgMin = 0.155 * inletDensity  # kg/min
inletMassFlow = inletMassFlowKgMin / 60.0

# setup symbols
y = Symbol('y')
r = Symbol('r', positive=True)
f = Symbol('f')
m = Symbol('m')
frac = sympy.Rational

# compute the mass flux function
# for laminar (https://www.kau.edu.sa/Files/0057863/Subjects/Chapter%208.pdf 8-17)
# velocityProfile = f*(1-(r**2)/((dia**2)/4))
# for turbulent
velocityProfile = f * (1 - (r / (dia / 2))) ** frac('1/7')

# for turblent flow
# int_0..2Pi int_0_dia/2 density*vel*r dr dtheta (note the extra r for integration)
massFluxProfile = velocityProfile * inletDensity * r * 2 * pi

# this equation does nto change in theta, so scale by 2*pi
massFlow = integrate(massFluxProfile, (r, 0.0, (dia / 2)-1E-10))

# Solver for f
conserved = Eq(massFlow, m)
fSolution = solve(conserved, f)[0]
print("VelocityFactor: ", fSolution)
print("VelocityFactor for ", inletMassFlow, " kg/s is ", fSolution.subs(m, inletMassFlow))

#
# Plot the velocity
fSolution = fSolution.subs(m, inletMassFlow)
velocityProfile = velocityProfile.subs(f, fSolution)
velocityProfile = velocityProfile.subs(r, (y - yc))

plot(velocityProfile, (y, yo, yh))
