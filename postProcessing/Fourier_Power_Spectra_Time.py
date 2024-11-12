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

for t in range(0, 300):
    print(t)
    Ex[:, :, t], Ey[:, :, t] = np.gradient(phi_data[:, :, t], dx, dy)

E_magnitude = np.sqrt(Ex**2 + Ey**2)

E_avg = np.mean(E_magnitude, axis=(0, 1))

delta_E = E_magnitude - E_avg[np.newaxis, np.newaxis, :]

delta_E_k = np.fft.fft(delta_E, axis=2)

power_spectrum_t = np.abs(delta_E_k)

power_spectrum_avg_t = np.mean(power_spectrum_t, axis=(0, 1))

frequencies = np.fft.fftfreq(Nt, d=times[1] - times[0])

plt.figure(figsize=(8, 6))
plt.plot(frequencies[:Nt // 2], power_spectrum_avg_t[:Nt // 2], label='Training Data', color='blue')
plt.xlabel('Frequency $\omega$ (Hz)')
plt.ylabel(r'$|\delta \hat{E}(\omega)|$')
plt.legend()
plt.title('Fourier Power Spectrum of Electric Field Turbulence Averaged Over Space')
plt.show()
