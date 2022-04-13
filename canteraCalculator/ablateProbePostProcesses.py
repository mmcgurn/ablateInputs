import h5py
import numpy
import numpy as np
import os, sys
import cantera as ct
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

print('Runnning Cantera version: ' + ct.__version__)

# Load in the cti mech file
gas = ct.Solution('../mechanisms/grimech30.cti')


def computeTemperature(row):
    density = row['euler_rho']
    totalEnergy = row['euler_rhoE'] / density

    kineticEnergy = 0.0
    for c in range(3):
        kineticEnergy += row['euler_rhoVel'+str(c)]/density
    internalEnergy = totalEnergy - 0.5 * kineticEnergy
    specificVolume = 1.0 / density

    yiArray = []
    for species in gas.species_names:
        yi = row['densityYi_'+ species]/density
        yiArray.append(yi)

    # compute the reference enthalpy
    gas.TPY = 298.15, 101325.0, yiArray
    enthalpyOfFormation = gas.h

    gas.UVY = (internalEnergy + enthalpyOfFormation), specificVolume, yiArray

    return gas.T



def main(probeFileName):
    probeData = pd.read_csv(probeFileName)
    probeData['temperature'] = probeData.apply(computeTemperature, axis=1)

    # updated probe name
    probeFilePath = Path(probeFileName)
    probeFileOutputDir =  os.path.join(str(probeFilePath.parent) ,"mod")
    os.makedirs(probeFileOutputDir, exist_ok=True)
    probeFileOutput = os.path.join(probeFileOutputDir, probeFilePath.stem + ".csv")
    probeData.to_csv(probeFileOutput)
    return probeData

if __name__ == '__main__':
    # dir = '/Users/mcgurn/scratch/results/slabBurnerChem.3D.V18-pgs_560x80x80n384/med/'
    # temperatureData = pd.DataFrame()
    # for file in os.listdir(dir):
    #     if file.endswith(".csv"):
    #         print('starting ...', file)
    #         subFile = main(os.path.join(dir, file))
    #
    #         temperatureData[file] = subFile['temperature']
    #         temperatureData['time'] = subFile['time']
    #
    #         print('done with ...', file)
    #
    # temperatureData.to_csv(os.path.join(dir, 'temperature.csv'))

    # plot
    probeData = pd.read_csv('/Users/mcgurn/scratch/results/slabBurnerChem.3D.V18-pgs_560x80x80n384/back/temperature.csv')
    probeData.set_index('time', drop=True)
    print(probeData.head())
    probeData.plot(x='time', y=['back.1','back.2','back.3','back.4','back.5','back.6','back.7','back.8'])

    plt.xlabel('$time(s)$',fontsize=10)
    plt.ylabel('$T$ (K)',fontsize=10)
    plt.grid(True)
    plt.legend(loc=0)

    plt.show()


