import h5py
import numpy as np
import matplotlib.pyplot as plt

hdf5_filepath = '/Users/rubencarpenter/Desktop/Yale Submissions/Year 2/S&DS 689/Project/test_MuSHrooM_processed-60024954.h5'

with h5py.File(hdf5_filepath, 'r') as f:
    phi_data = f['phi_data'][:]
    grid = f['X'][:]
    times = f['times'][:]

    Nx = f.attrs['Nx']
    Ny = f.attrs['Ny']
    Nt = f.attrs['Nt']

Ex = np.zeros_like(phi_data)
Ey = np.zeros_like(phi_data)

dx = grid[0, 1, 0] - grid[0, 0, 0]
dy = grid[1, 0, 1] - grid[1, 0, 0]

for t in range(0, Nt):
    Ex[:, :, t], Ey[:, :, t] = np.gradient(phi_data[:, :, t], dx, dy)

E_magnitude = np.sqrt(Ex**2 + Ey**2)
E_avg = np.mean(E_magnitude, axis=(0, 1))
delta_E = E_magnitude - E_avg[np.newaxis, np.newaxis, :]

delta_E_k = np.fft.fft2(delta_E, axes=(0, 1))

power_spectrum = np.abs(delta_E_k)**2
power_spectrum_avg = np.mean(power_spectrum, axis=2)

kx = np.fft.fftfreq(Nx, d=dx) * 2 * np.pi
ky = np.fft.fftfreq(Ny, d=dy) * 2 * np.pi

kx_grid, ky_grid = np.meshgrid(kx, ky)
k = np.sqrt(kx_grid**2 + ky_grid**2)

k_flat = k.flatten()
power_spectrum_avg_flat = power_spectrum_avg.flatten()

sorted_indices = np.argsort(k_flat)
k_sorted = k_flat[sorted_indices]
power_spectrum_avg_sorted = power_spectrum_avg_flat[sorted_indices]

plt.figure(figsize=(8, 6))
plt.plot(k_sorted, power_spectrum_avg_sorted, label='Can I get this to work?', linestyle='-', color='blue')
plt.xlim(0, 6)
plt.xlabel(r'$k$')
plt.ylabel(r'$|\delta \hat{E}(k)|$')
plt.legend()
plt.title('Fourier Power Spectrum')
plt.show()
