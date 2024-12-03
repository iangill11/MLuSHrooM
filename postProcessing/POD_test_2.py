# Perform POD analysis on the data
def perform_pod(data, num_time_slices):
    T, Nx, Ny = data.shape

    num_time_slices = min(num_time_slices, T)

    # Reshaping the data to match the dimensions for POD (flatten the spatial dimensions)
    reshaped_data = data[:num_time_slices, :, :].reshape(num_time_slices, Nx * Ny)

    # Center the data by subtracting the mean
    mean_data = np.mean(reshaped_data, axis=0)
    centered_data = reshaped_data - mean_data

    # Calculate the covariance matrix
    covariance_matrix = np.dot(centered_data.T, centered_data) / (num_time_slices - 1)

    # Perform Singular Value Decomposition
    U, S, Vt = np.linalg.svd(covariance_matrix)

    return U, S, mean_data

# Function to plot POD results
def plot_pod_analysis(singular_values, pod_modes, Nx, Ny, title_prefix):
    num_modes_to_plot = 8

    # Plot singular values
    plt.figure(figsize=(10, 5))
    plt.plot(singular_values[:num_modes_to_plot], marker='o')
    plt.title(f'{title_prefix} Singular Values')
    plt.xlabel('Mode Number')
    plt.ylabel('Singular Value')
    plt.grid()
    plt.show()

    # Plot cumulative energy
    cumulative_energy = np.cumsum(singular_values) / np.sum(singular_values)
    plt.figure(figsize=(10, 5))
    plt.plot(cumulative_energy[:num_modes_to_plot], marker='o')
    plt.title(f'{title_prefix} Cumulative Energy')
    plt.xlabel('Mode Number')
    plt.ylabel('Cumulative Energy')
    plt.grid()
    plt.show()

    # Plot first few POD modes
    fig, axes = plt.subplots(2, 4, figsize=(15, 6))
    for i in range(num_modes_to_plot):
        mode = pod_modes[:, i].reshape(Nx, Ny)  # Reshape to the original spatial dimensions
        ax = axes.flat[i]
        im = ax.imshow(mode, cmap='jet')
        ax.set_title(f'{title_prefix} POD Mode {i+1}')
        plt.colorbar(im, ax=ax, orientation='vertical')
    plt.tight_layout()
    plt.show()

# Function to plot comparison of original and reconstructed POD modes
def plot_pod_modes_comparison(original_U, reconstructed_U, Nx, Ny, title_prefix):
    num_modes_to_plot = 8
    fig, axes = plt.subplots(num_modes_to_plot, 2, figsize=(15, 24))
    for i in range(num_modes_to_plot):
        # Reshape the modes to the original spatial dimensions
        original_mode = original_U[:, i].reshape(Nx, Ny)
        reconstructed_mode = reconstructed_U[:, i].reshape(Nx, Ny)

        # Plot original mode
        ax = axes[i, 0]
        im = ax.imshow(original_mode, cmap='jet')
        ax.set_title(f'Original Mode {i+1}')
        plt.colorbar(im, ax=ax, orientation='vertical')

        # Plot reconstructed mode
        ax = axes[i, 1]
        im = ax.imshow(reconstructed_mode, cmap='jet')
        ax.set_title(f'Reconstructed Mode {i+1}')
        plt.colorbar(im, ax=ax, orientation='vertical')

    plt.tight_layout()
    plt.show()

# Initialize parameters
batch_size_cvae = 32
num_time_slices = 544  # Adjust to the desired number of time slices

# Split data into train and test
N_train = batch_size_cvae * 16

model.eval()
original_data_all = []
reconstructed_data_all = []

with torch.no_grad():
    for batch_idx, data_batch in enumerate(phi_train_loader):
        batch_data = data_batch  # (batch_size, 1, 64, 64)
        
        # Perform reconstruction
        recon_batch, _, _, _ = model(batch_data)

        # Collect the original and reconstructed batches
        original_data_all.append(batch_data.squeeze(1).numpy())  # Remove channel dimension: (batch_size, 64, 64)
        reconstructed_data_all.append(recon_batch.squeeze(1).numpy())  # Remove channel dimension

    # Convert all batches to numpy arrays
    original_data_all = np.concatenate(original_data_all, axis=0)
    reconstructed_data_all = np.concatenate(reconstructed_data_all, axis=0)


    # Plot the first time slice of original data and reconstructed data
    time_slice_idx = 400  # Select the desired time slice

    # Get the corresponding slices from original and reconstructed data
    original_slice = original_data_all[time_slice_idx]
    reconstructed_slice = reconstructed_data_all[time_slice_idx]

    # Plotting original vs reconstructed slice for the selected time slice
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot original data
    axes[0].imshow(original_slice, cmap='jet')
    axes[0].set_title(f'Original Time Slice {time_slice_idx}')
    axes[0].axis('off')

    # Plot reconstructed data
    axes[1].imshow(reconstructed_slice, cmap='jet')
    axes[1].set_title(f'Reconstructed Time Slice {time_slice_idx}')
    axes[1].axis('off')

    plt.tight_layout()
    plt.show()


    # Perform POD analysis on the entire dataset
    original_U, original_S, original_mean = perform_pod(original_data_all, num_time_slices)
    reconstructed_U, reconstructed_S, reconstructed_mean = perform_pod(reconstructed_data_all, num_time_slices)

    # Plot the comparison of the first 8 modes
    plot_pod_modes_comparison(original_U, reconstructed_U, 64, 64, 'POD Modes Comparison')

    # Optionally, plot singular values and cumulative energy for both original and reconstructed data
    plt.figure(figsize=(10, 5))
    plt.plot(original_S[:10], marker='o', label='Original Singular Values')
    plt.plot(reconstructed_S[:10], marker='o', label='Reconstructed Singular Values')
    plt.title('Comparison of Singular Values (First 10)')
    plt.xlabel('Mode Number')
    plt.ylabel('Singular Value')
    plt.legend()
    plt.grid()
    plt.show()

    # Plot cumulative energy
    cumulative_energy_original = np.cumsum(original_S) / np.sum(original_S)
    cumulative_energy_reconstructed = np.cumsum(reconstructed_S) / np.sum(reconstructed_S)
    
    plt.figure(figsize=(10, 5))
    plt.plot(cumulative_energy_original[:10], marker='o', label='Original Cumulative Energy')
    plt.plot(cumulative_energy_reconstructed[:10], marker='o', label='Reconstructed Cumulative Energy')
    plt.title('Comparison of Cumulative Energy (First 10 Modes)')
    plt.xlabel('Mode Number')
    plt.ylabel('Cumulative Energy')
    plt.legend()
    plt.grid()
    plt.show()
  
# Function to calculate relative L2 difference between original and reconstructed modes
def relative_l2_difference(original_modes, reconstructed_modes, num_modes_to_plot=20):
    l2_differences = []
    
    for i in range(num_modes_to_plot):
        original_mode = original_modes[:, i].reshape(64, 64)
        reconstructed_mode = reconstructed_modes[:, i].reshape(64, 64)

        l2_diff = np.linalg.norm(original_mode - reconstructed_mode) / np.linalg.norm(original_mode)
        l2_differences.append(l2_diff)
    
    return l2_differences

# Add the calculation and plot of the relative L2 difference
l2_diff = relative_l2_difference(original_U, reconstructed_U, num_modes_to_plot=20)

# Plot the relative L2 difference for the first 20 modes
plt.figure(figsize=(10, 5))
plt.plot(l2_diff, marker='o', color='b')
plt.title('Relative L2 Difference Between Original and Reconstructed POD Modes (First 20 Modes)')
plt.xlabel('Mode Number')
plt.ylabel('Relative L2 Difference')
plt.grid()
plt.show()
