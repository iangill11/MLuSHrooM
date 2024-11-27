import numpy as np
import deepxde as dde
import h5py

dde.config.set_default_backend("pytorch")

hdf5_filepath = '/Users/iangill/Downloads/test_MuSHrooM_processed-60024954.h5'

with h5py.File(hdf5_filepath, 'r') as f:
    # Read the datasets
    #tempIon_data = f['tempIon_data'][:]
    #denIon_data = f['denIon_data'][:]
    #tempElc_data = f['tempElc_data'][:]
    #denElc_data = f['denElc_data'][:]
    phi_data = f['phi_data'][:]
    times = f['times'][:]
    grid = f['X'][:]
    # Read the attributes (dimensions of data)
    Nx = f.attrs['Nx']
    Ny = f.attrs['Ny']
    Nt = f.attrs['Nt']

# exchange time dimension with batch dimension
phi_data = np.transpose(phi_data, (2, 1, 0))



