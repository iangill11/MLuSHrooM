import numpy as np
import torch
import h5py
from torchvision import transforms
import matplotlib.pyplot as plt

# Load the training data
hdf5_filepath = '/Users/rubencarpenter/Desktop/Yale_Submissions/Year 2/S&DS 689/Project/test_MuSHrooM_processed-60024954.h5'

with h5py.File(hdf5_filepath, 'r') as f:
    phi_data = f['phi_data'][:]
    times = f['times'][:]
    grid = f['X'][:]
    Nx = f.attrs['Nx']
    Ny = f.attrs['Ny']
    Nt = f.attrs['Nt']

# Define transform for resizing
input_dim = 64  # Target dimension for resizing
transform = transforms.Resize((input_dim, input_dim))  # Resizing to (64, 64)

# Function to perform downsampling
def downsample_data(data, transform):
    downsampled_data = []
    for t in range(data.shape[2]):
        slice_2d = torch.tensor(data[:, :, t]).unsqueeze(0).unsqueeze(0)
        resized_slice = transform(slice_2d).squeeze().numpy()
        downsampled_data.append(resized_slice)
    return np.stack(downsampled_data, axis=-1)

phi_data_downsampled = downsample_data(phi_data, transform)

# Rescale data to range [0, 1]
def rescale(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

phi_data_downsampled = rescale(phi_data_downsampled)

def perform_pod(data, num_time_slices):
    Nx, Ny, T = data.shape

    # Ensure num_time_slices does not exceed the actual number of time slices available
    num_time_slices = min(num_time_slices, T)

    reshaped_data = data[:, :, :num_time_slices].reshape(Nx * Ny, num_time_slices).T

    mean_data = np.mean(reshaped_data, axis=0)
    centered_data = reshaped_data - mean_data

    covariance_matrix = np.dot(centered_data.T, centered_data) / (num_time_slices - 1)
    U, S, Vt = np.linalg.svd(covariance_matrix)

    return U, S, mean_data

num_time_slices = 544  # adjust to have code run faster
pod_modes, singular_values, mean_data = perform_pod(phi_data_downsampled, num_time_slices)

# Number of modes to plot
num_modes_to_plot = 10

# Plot singular values
plt.figure(figsize=(10, 5))
plt.plot(singular_values[:num_modes_to_plot], marker='o')
plt.title('Singular Values')
plt.xlabel('Mode Number')
plt.ylabel('Singular Value')
plt.grid()
plt.show()

# Plot cumulative energy
cumulative_energy = np.cumsum(singular_values) / np.sum(singular_values)
plt.figure(figsize=(10, 5))
plt.plot(cumulative_energy[:num_modes_to_plot], marker='o')
plt.title('Cumulative Energy')
plt.xlabel('Mode Number')
plt.ylabel('Cumulative Energy')
plt.grid()
plt.show()

# Plot first few POD modes
fig, axes = plt.subplots(2, 5, figsize=(15, 6))
for i in range(num_modes_to_plot):
    mode = pod_modes[:, i].reshape(input_dim, input_dim)
    ax = axes.flat[i]
    im = ax.imshow(mode, cmap='jet')
    ax.set_title(f'POD Mode {i+1}')
    plt.colorbar(im, ax=ax, orientation='vertical')
plt.tight_layout()
plt.show()
