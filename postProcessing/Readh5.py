import h5py
import numpy as np


hdf5_filepath = 'C:/Users/lukas/OneDrive - Yale University/Personal/Studies/Semester 3/S&DS 689/Project/test_datah5.h5'

# Open the HDF5 file in read mode
with h5py.File(hdf5_filepath, 'r') as f:
    # Read the datasets
    tempIon_data = f['tempIon_data'][:]
    denIon_data = f['denIon_data'][:]
    tempElc_data = f['tempElc_data'][:]
    denElc_data = f['denElc_data'][:]
    phi_data = f['phi_data'][:]
    grid = f['X'][:]
    times = f['times'][:]
    

    Nx = f.attrs['Nx']
    Ny = f.attrs['Ny']
    Nt = f.attrs['Nt']
