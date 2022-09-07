import h5py
import numpy
import numpy as np
import os, sys
import cantera as ct

print('Runnning Cantera version: ' + ct.__version__)

def toStringAtt(string):
    ascii_type = h5py.string_dtype('ascii', len(string))
    return np.array(string.encode("latin-1"), dtype=ascii_type)


def computeYi(fields):
    eulerField = fields['solution_euler']
    densityYiField = fields['solution_densityYi']

    # create the field
    if 'yi' in fields:
        del fields['yi']
    yiField = fields.create_dataset("yi", densityYiField.shape, densityYiField.dtype)

    # copy over the attributes
    for name, value in densityYiField.attrs.items():
      yiField.attrs.create(name, value)

    components = yiField.shape[-1]

    for t in range(len(eulerField)):
        density = eulerField[t, :, 0]
        for c in range(components):
            yiField[t, :, c] = densityYiField[t, :, c]/density

def computeVel(fields):
    eulerField = fields['solution_euler']

    dim = eulerField.shape[-1]-2
    velShape = list(eulerField.shape)
    velShape[-1] = dim

    # create the field
    if 'vel' in fields:
        del fields['vel']
    velField = fields.create_dataset("vel", tuple(velShape), eulerField.dtype)

    # add the required fields
    velField.attrs.create('timestepping', eulerField.attrs['timestepping'])
    velField.attrs.create('componentName0', toStringAtt('vel0'))
    velField.attrs.create('componentName1', toStringAtt('vel1'))
    velField.attrs.create('componentName2', toStringAtt('vel2'))
    velField.attrs.create('vector_field_type', toStringAtt('vector'))

    for t in range(len(eulerField)):
        density = eulerField[t, :, 0]
        for d in range(dim):
            velField[t, :, d] = eulerField[t, :, d+2]/density


def computeTemperature(fields):
    eulerField = fields['solution_euler']
    yiField = fields['yi']
    velField = fields['vel']

    tShape = list(eulerField.shape)
    tShape.pop()
    dim = velField.shape[-1]

    # create the field
    if 'T' in fields:
        del fields['T']
    tField = fields.create_dataset("T", tuple(tShape), eulerField.dtype)

    # add the required fields
    tField.attrs.create('timestepping', eulerField.attrs['timestepping'])
    tField.attrs.create('componentName0', toStringAtt('0'))
    tField.attrs.create('vector_field_type', toStringAtt('scalar'))

    # Load in the cti mech file
    gas = ct.Solution('../mechanisms/grimech30.cti')

    speciesList = []
    for s in range(yiField.shape[-1]):
        attName = 'componentName' + str(s)
        speciesList.append(yiField.attrs[attName])

    for t in range(len(eulerField)):
        density = eulerField[t, :, 0]
        totalEnergy = eulerField[t, :, 1]/density
        for p in range(len(totalEnergy)):
            # compute the kinetic totalEnergy
            kineticEnergy = 0.0
            for c in range(dim):
                kineticEnergy += velField[t, p, c]**2
            internalEnergy = totalEnergy[p] - 0.5*kineticEnergy
            specificVolume = 1.0/density[p]

            # build the yi
            yi = dict(zip(speciesList, yiField[t, p, :]))

            # compute the reference enthalpy
            gas.TPY = 298.15, 101325.0, yi
            enthalpyOfFormation = gas.h

            # set the internal energy
            gas.UVY = (internalEnergy + enthalpyOfFormation), specificVolume, yi
            tField[t, p] = gas.T

def main(hdfFileName):
    h5 = h5py.File(hdfFileName, 'r+')

    # Check each of the cell fields
    fields = h5['cell_fields']

    import time
    start = time.time()
    computeYi(fields)
    computeVel(fields)
    computeTemperature(fields)
    end = time.time()
    print(end - start)
    for field in fields.items():
        print(field)

if __name__ == '__main__':
    hdfFileName = "/Users/mcgurn/scratch/ablateInputs/slabBurnerChem.2D.V17/_2dSlabBurner/domain.hdf5"
    main(hdfFileName)
