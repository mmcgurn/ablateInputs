import h5py
import numpy as np
import os, sys

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


def main(hdfFileName):
    h5 = h5py.File(hdfFileName, 'r+')

    # Check each of the cell fields
    fields = h5['cell_fields']

    computeYi(fields)
    computeVel(fields)

    for field in fields.items():
        print(field)

if __name__ == '__main__':
    hdfFileName = "/Users/mcgurn/scratch/ablateInputs/slabBurnerChem.2D.V17/_2dSlabBurner/domain.hdf5"
    main(hdfFileName)
