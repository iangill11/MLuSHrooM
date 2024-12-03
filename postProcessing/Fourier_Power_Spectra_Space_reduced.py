import h5py
import numpy as np
import matplotlib.pyplot as plt
from torchvision import transforms
import torch
hdf5_filepath = '/Users/iangill/Downloads/test_MuSHrooM_processed-60024954.h5'


with h5py.File(hdf5_filepath, 'r') as f:
    phi_data = f['phi_data'][:]
    grid = f['X'][:]
    times = f['times'][:]

    Nx = f.attrs['Nx']
    Ny = f.attrs['Ny']
    Nt = f.attrs['Nt']
phi_data = np.transpose(phi_data, (2, 1, 0))

def compute_fourier_power_spectrum(data):
    Nt, Nx, Ny = data.shape  # Unpacking dimensions based on your data shape
    Ex = np.zeros_like(data)
    Ey = np.zeros_like(data)

    for t in range(Nt):
        # Compute gradients in the spatial directions
        Ex[t, :, :], Ey[t, :, :] = np.gradient(data[t, :, :], dx, dy)

    # Compute the magnitude of the electric field
    E_magnitude = np.sqrt(Ex**2 + Ey**2)
    
    # Compute the average magnitude of the electric field across spatial dimensions
    E_avg = np.mean(E_magnitude, axis=(1, 2))  # Average over x and y
    delta_E = E_magnitude - E_avg[:, np.newaxis, np.newaxis]  # Subtract the average

    # Perform 2D Fourier Transform (on x and y axes) for each time slice
    delta_E_k = np.fft.fft2(delta_E, axes=(1, 2))
    
    # Compute the power spectrum
    power_spectrum = np.abs(delta_E_k)**2
    
    # Average the power spectrum over time to get a single value for each spatial frequency
    power_spectrum_avg = np.mean(power_spectrum, axis=0)

    return power_spectrum_avg

# Compute dx and dy
input_dim = 64
transform = transforms.Resize((input_dim, input_dim))
x_coords = torch.tensor(grid[0, :, :]).unsqueeze(0)
y_coords = torch.tensor(grid[1, :, :]).unsqueeze(0)
resized_x_coords = transform(x_coords).squeeze().numpy()
resized_y_coords = transform(y_coords).squeeze().numpy()
dx = resized_x_coords[ 1, 0] - resized_x_coords[ 0, 0]
dy = resized_y_coords[ 0, 1] - resized_y_coords[ 0, 0]
phi_data = torch.tensor(phi_data.copy())
phi_data = transform(phi_data).numpy().astype(np.float32)
# Compute Fourier Power Spectrum for the training data
input_power_spectrum_avg = compute_fourier_power_spectrum(phi_data)

kx = np.fft.fftfreq(64, d=dx) * 2 * np.pi  # Frequencies in x-direction
ky = np.fft.fftfreq(64, d=dy) * 2 * np.pi  # Frequencies in y-direction

kx_grid, ky_grid = np.meshgrid(kx, ky)  # Create a mesh grid for kx, ky
k = np.sqrt(kx_grid**2 + ky_grid**2)  # Magnitude of k

# Flatten the power spectrum and k arrays
k_flat = k.flatten()
input_power_spectrum_avg_flat = input_power_spectrum_avg.flatten()

# Sort the k values and corresponding power spectrum values
sorted_indices = np.argsort(k_flat)
k_sorted = k_flat[sorted_indices]
input_power_spectrum_avg_sorted = input_power_spectrum_avg_flat[sorted_indices]

# Plot the Fourier Power Spectrum for the training data
plt.figure(figsize=(8, 6))
plt.plot(k_sorted, input_power_spectrum_avg_sorted, label='Training Data Power Spectrum', linestyle='-', color='blue')
plt.xlim(0, 6)  # Adjust the x-axis as needed
plt.xlabel(r'$k$')
plt.ylabel(r'$|\delta \hat{E}(k)|$')
plt.legend()
plt.title('Fourier Power Spectrum of Training Data')
plt.show()
