# general imports
import matplotlib.pyplot as plt
import numpy as np
import os

# torch imports
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

# data loading and processing
import h5py
from torchvision import transforms
from torch.utils.data import Dataset
from torch.utils.data import DataLoader


# Set Script ID
ID = 1

# make directory for saving files
dir_name = 'CVAE_{}'.format(ID)

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
training_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CVAE Training")
model_dir = os.path.join(training_dir, dir_name)

os.chdir(training_dir + "") # go to training directory
os.makedirs(dir_name, exist_ok=True) # make directory for model files
os.chdir(script_dir) # go back to script directory

""" Set Hyperparameters """
# CVAE
# model training
lr_cvae = 0.001
num_epochs_cvae = 1000


# model data
batch_size_cvae = 32
input_dim = 64 # for square images
latent_dim = 128

w_phi = 1.0
w_grad = 10.0
w_kl = 0.01

# DNN
# model training
lr_dnn = 0.001
num_epochs_dnn = 10000

# model data
batch_size_dnn = 32



""" Load Test and Train Data """
# load data
phi_data = torch.load('downscaled_data{}old.pt'.format(input_dim))
N_train_batches = 38


# define class for dataset
class PhiDataset(Dataset):
    def __init__(self, phi_data, transform=None):
        self.data = phi_data
        self.transform = transform

    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, idx):
        image = torch.tensor(self.data[idx,:,:], dtype=torch.float32)
        # Add a channel dimension to make it (1, 256, 256)
        image = image.unsqueeze(0)  # Shape becomes (1, 256, 256)

        if self.transform:
            image = self.transform(image)

        return image
    
N_train = batch_size_cvae*N_train_batches

# train data loader
phi_train = PhiDataset(phi_data[0:N_train, :, :])
phi_train_loader = DataLoader(phi_train, batch_size=batch_size_cvae, shuffle=True)

# test data loader
phi_test = PhiDataset(phi_data[N_train:, :, :])
phi_test_loader = DataLoader(phi_test, batch_size=batch_size_cvae, shuffle=True)



""" Define the CVAE model """
# Encoder Network
class Encoder(nn.Module):
    def __init__(self, latent_dim):
        super(Encoder, self).__init__()
        
        # Convolutional layers
        self.conv1 = nn.Conv2d(1, 32, kernel_size=4, stride=2, padding=1)  # 64x64 -> 32x32
        self.bn1 = nn.BatchNorm2d(32)
        
        self.conv2 = nn.Conv2d(32, 64, kernel_size=4, stride=2, padding=1)  # 32x32 -> 16x16
        self.bn2 = nn.BatchNorm2d(64)
        
        self.conv3 = nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1)  # 16x16 -> 16x16
        self.bn3 = nn.BatchNorm2d(128)
        
        self.conv4 = nn.Conv2d(128, 256, kernel_size=4, stride=2, padding=1)  # 16x16 -> 8x8
        self.bn4 = nn.BatchNorm2d(256)
        
        self.conv5 = nn.Conv2d(256, 256, kernel_size=3, stride=1, padding=1)  # 8x8 -> 8x8
        self.bn5 = nn.BatchNorm2d(256)
        
        # Fully connected layers
        self.fc_mu = nn.Linear(256 * int((input_dim / 8))**2, latent_dim)
        self.fc_logvar = nn.Linear(256 * int((input_dim / 8))**2, latent_dim)

    def forward(self, x):
        # Apply convolution -> batch normalization -> activation
        x = F.relu(self.bn1(self.conv1(x)))
        x = F.relu(self.bn2(self.conv2(x)))
        x = F.relu(self.bn3(self.conv3(x)))
        x = F.relu(self.bn4(self.conv4(x)))
        x = F.relu(self.bn5(self.conv5(x)))
        
        # Flatten before fully connected layers
        x = x.view(x.size(0), -1)
        
        # Compute mu and logvar for the latent space
        mu = self.fc_mu(x)
        logvar = self.fc_logvar(x)
        return mu, logvar
    
    
# Decoder Network
class Decoder(nn.Module):
    def __init__(self, latent_dim):
        super(Decoder, self).__init__()
        
        # Fully connected layer to project latent space to feature map
        self.fc = nn.Linear(latent_dim, 256 * int((input_dim / 8))**2)
        
        # Deconvolutional layers
        self.deconv1 = nn.ConvTranspose2d(256, 256, kernel_size=3, stride=1, padding=1)  # 8x8 -> 8x8
        self.deconv2 = nn.ConvTranspose2d(256, 128, kernel_size=4, stride=2, padding=1)  # 8x8 -> 16x16
        self.deconv3 = nn.ConvTranspose2d(128, 64, kernel_size=3, stride=1, padding=1)   # 16x16 -> 16x16
        self.deconv4 = nn.ConvTranspose2d(64, 32, kernel_size=4, stride=2, padding=1)    # 16x16 -> 32x32
        self.deconv5 = nn.ConvTranspose2d(32, 1, kernel_size=4, stride=2, padding=1)     # 32x32 -> 64x64

    def forward(self, z):
        # Project latent space to feature map
        z = F.relu(self.fc(z))
        z = z.view(z.size(0), 256, int((input_dim / 8)), int((input_dim / 8)))  # Reshape to (batch_size, 256, 8, 8)
        
        # Apply deconvolution -> activation
        z = F.relu(self.deconv1(z))
        z = F.relu(self.deconv2(z))
        z = F.relu(self.deconv3(z))
        z = F.relu(self.deconv4(z))
        return torch.sigmoid(self.deconv5(z))  # Output scaled between 0 and 1


class CVAE(nn.Module):
    def __init__(self, latent_dim=latent_dim):
        super(CVAE, self).__init__()
        self.encoder = Encoder(latent_dim)
        self.decoder = Decoder(latent_dim)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        epsilon = torch.randn_like(std) # std just gives dimension of tensor to give back
        return mu + epsilon * std

    def forward(self, x):
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decoder(z)
        return recon_x, mu, logvar, z # return reconstructed data, mu, logvar, and the latent state z
    
    
# define function to get gradient loss
def gradient_loss(recon_x, x):
    # Ensure gradients can be computed
    x.requires_grad_(True)
    recon_x.requires_grad_(True)

    # Compute gradient of input with respect to spatial dimensions
    grad_input = torch.autograd.grad(outputs=x.sum(), inputs=x, create_graph=True)[0]

    # Compute gradient of output (reconstructed) with respect to spatial dimensions
    grad_output = torch.autograd.grad(outputs=recon_x.sum(), inputs=recon_x, create_graph=True)[0]

    # Compute L2 norm of the difference between the gradients
    grad_loss = F.mse_loss(grad_output, grad_input, reduction='sum')
    
    return grad_loss


def loss_function(recon_x, x, mu, logvar):
    # Reconstruction loss (binary cross-entropy)
    recon_loss = F.mse_loss(recon_x, x, reduction='sum') # might need to change to L2 norm
    
    recon_grad_loss = gradient_loss(recon_x, x)
    
    # KL divergence - regularizes the distribution of the latent space to be close to a standard normal distribution
    kl_divergence = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    
    return w_phi*recon_loss + w_grad*recon_grad_loss + w_kl*kl_divergence


""" Train the CVAE model """	
os.chdir(model_dir) # go to model directory to save files

model = CVAE(latent_dim=latent_dim)
optimizer = optim.Adam(model.parameters(), lr=lr_cvae)

# Training loop
losses = []

for epoch in range(num_epochs_cvae):
    model.train()
    train_loss = 0
    for i, batch in enumerate(phi_train_loader):
        optimizer.zero_grad()
        recon_batch, mu, logvar, z = model(batch)
        loss = loss_function(recon_batch, batch, mu, logvar)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()
    
    avg_loss = train_loss / len(phi_train_loader.dataset)
    
    losses.append(avg_loss)
    # open file to write loss
    with open('CVAE_loss_{}.txt'.format(ID), 'a') as f:
        f.write('{},{}\n'.format(epoch, avg_loss))
        
# plot loss and save image
plt.figure(figsize=(10, 8))
plt.plot(losses, color='blue')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('CVAE Loss')
plt.grid()
plt.savefig('CVAE_loss{}.png'.format(ID))

    
    
# Set the model to evaluation mode
model.eval()
with torch.no_grad():
    # Get a batch of test data
    test_batch = next(iter(phi_test_loader))
    recon_batch, mu, logvar, z = model(test_batch)
    
    # Select the 27th image from the batch for display
    original_image = test_batch[12].squeeze().numpy()
    reconstructed_image = recon_batch[12].squeeze().numpy()

    # Create a single figure with two vertically-stacked subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))

    # Display the original image in the first subplot
    cax1 = ax1.imshow(original_image, cmap='jet', interpolation='nearest')
    cbar1 = fig.colorbar(cax1, ax=ax1, format='%.1e')
    cbar1.set_label(r'[e$L_n$ / ($\rho_s T_e$)] $\phi(x, y)$', rotation=270, labelpad=20)
    ax1.set_title('Original Image')
    ax1.set_xlabel(r'$x / \rho_s$')
    ax1.set_ylabel(r'$y / \rho_s$')

    # Display the reconstructed image in the second subplot
    cax2 = ax2.imshow(reconstructed_image, cmap='jet', interpolation='nearest')
    cbar2 = fig.colorbar(cax2, ax=ax2, format='%.1e')
    cbar2.set_label(r'[e$L_n$ / ($\rho_s T_e$)] $\phi(x, y)$', rotation=270, labelpad=20)
    ax2.set_title('Reconstructed Image')
    ax2.set_xlabel(r'$x / \rho_s$')
    ax2.set_ylabel(r'$y / \rho_s$')

    # Adjust layout to prevent overlapping
    plt.tight_layout()
    plt.savefig('CVAE_reconstruction{}.png'.format(ID))

# save model
torch.save(model.state_dict(), 'CVAE_model_{}.pth'.format(ID))



""" DNN Training """
# make directory for saving files
dir_name = 'CVAE_{}'.format(ID)

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
training_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CVAE Training")
model_dir = os.path.join(training_dir, dir_name)

# change directory to script directory
os.chdir(script_dir)


# define class for dataset
class PhiDataset(Dataset):
    def __init__(self, phi_data, transform=None):
        self.data = phi_data
        self.transform = transform

    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, idx):
        image = torch.tensor(self.data[idx,:,:], dtype=torch.float32)
        # Add a channel dimension to make it (1, 256, 256)
        image = image.unsqueeze(0)  # Shape becomes (1, 256, 256)

        if self.transform:
            image = self.transform(image)

        return image
    
N_train = batch_size_cvae*N_train_batches

# train data loader
phi_train = PhiDataset(phi_data[0:N_train, :, :])
phi_train_loader = DataLoader(phi_train, batch_size=batch_size_cvae, shuffle=True)

# test data loader
phi_test = PhiDataset(phi_data[N_train:, :, :])
phi_test_loader = DataLoader(phi_test, batch_size=batch_size_cvae, shuffle=True)



""" Define and Load the CVAE model """
# Encoder Network
class Encoder(nn.Module):
    def __init__(self, latent_dim):
        super(Encoder, self).__init__()
        
        # Convolutional layers
        self.conv1 = nn.Conv2d(1, 32, kernel_size=4, stride=2, padding=1)  # 64x64 -> 32x32
        self.bn1 = nn.BatchNorm2d(32)
        
        self.conv2 = nn.Conv2d(32, 64, kernel_size=4, stride=2, padding=1)  # 32x32 -> 16x16
        self.bn2 = nn.BatchNorm2d(64)
        
        self.conv3 = nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1)  # 16x16 -> 16x16
        self.bn3 = nn.BatchNorm2d(128)
        
        self.conv4 = nn.Conv2d(128, 256, kernel_size=4, stride=2, padding=1)  # 16x16 -> 8x8
        self.bn4 = nn.BatchNorm2d(256)
        
        self.conv5 = nn.Conv2d(256, 256, kernel_size=3, stride=1, padding=1)  # 8x8 -> 8x8
        self.bn5 = nn.BatchNorm2d(256)
        
        # Fully connected layers
        self.fc_mu = nn.Linear(256 * int((input_dim / 8))**2, latent_dim)
        self.fc_logvar = nn.Linear(256 * int((input_dim / 8))**2, latent_dim)

    def forward(self, x):
        # Apply convolution -> batch normalization -> activation
        x = F.relu(self.bn1(self.conv1(x)))
        x = F.relu(self.bn2(self.conv2(x)))
        x = F.relu(self.bn3(self.conv3(x)))
        x = F.relu(self.bn4(self.conv4(x)))
        x = F.relu(self.bn5(self.conv5(x)))
        
        # Flatten before fully connected layers
        x = x.view(x.size(0), -1)
        
        # Compute mu and logvar for the latent space
        mu = self.fc_mu(x)
        logvar = self.fc_logvar(x)
        return mu, logvar
    
    
# Decoder Network
class Decoder(nn.Module):
    def __init__(self, latent_dim):
        super(Decoder, self).__init__()
        
        # Fully connected layer to project latent space to feature map
        self.fc = nn.Linear(latent_dim, 256 * int((input_dim / 8))**2)
        
        # Deconvolutional layers
        self.deconv1 = nn.ConvTranspose2d(256, 256, kernel_size=3, stride=1, padding=1)  # 8x8 -> 8x8
        self.deconv2 = nn.ConvTranspose2d(256, 128, kernel_size=4, stride=2, padding=1)  # 8x8 -> 16x16
        self.deconv3 = nn.ConvTranspose2d(128, 64, kernel_size=3, stride=1, padding=1)   # 16x16 -> 16x16
        self.deconv4 = nn.ConvTranspose2d(64, 32, kernel_size=4, stride=2, padding=1)    # 16x16 -> 32x32
        self.deconv5 = nn.ConvTranspose2d(32, 1, kernel_size=4, stride=2, padding=1)     # 32x32 -> 64x64

    def forward(self, z):
        # Project latent space to feature map
        z = F.relu(self.fc(z))
        z = z.view(z.size(0), 256, int((input_dim / 8)), int((input_dim / 8)))  # Reshape to (batch_size, 256, 8, 8)
        
        # Apply deconvolution -> activation
        z = F.relu(self.deconv1(z))
        z = F.relu(self.deconv2(z))
        z = F.relu(self.deconv3(z))
        z = F.relu(self.deconv4(z))
        return torch.sigmoid(self.deconv5(z))  # Output scaled between 0 and 1


class CVAE(nn.Module):
    def __init__(self, latent_dim=latent_dim):
        super(CVAE, self).__init__()
        self.encoder = Encoder(latent_dim)
        self.decoder = Decoder(latent_dim)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        epsilon = torch.randn_like(std) # std just gives dimension of tensor to give back
        return mu + epsilon * std

    def forward(self, x):
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decoder(z)
        return recon_x, mu, logvar, z # return reconstructed data, mu, logvar, and the latent state z
    
# Load the model
model = CVAE(latent_dim=latent_dim)
model.load_state_dict(torch.load('CVAE Training/CVAE_3/CVAE_model_3.pth'))

""" Prepare Test and Train data for the DNN model """
# get length of training and test set
train_len = len(phi_train_loader.dataset)
test_len = len(phi_test_loader.dataset)

# initialize empty tensors to store z-vectors
z_train_data = torch.empty(train_len, latent_dim)
z_test_data = torch.empty(test_len, latent_dim)

# get z-vectors for training and test set
model.eval()
with torch.no_grad():
    for i, batch in enumerate(phi_train_loader):
        _, _, _, z = model(batch)
        z_train_data[i*batch_size_cvae:(i+1)*batch_size_cvae, :] = z
        
    for i, batch in enumerate(phi_test_loader):
        _, _, _, z = model(batch)
        z_test_data[i*batch_size_cvae:(i+1)*batch_size_cvae, :] = z
        
# print shapes of data
print("Training data shape: ", z_train_data.shape)
print("Test data shape: ", z_test_data.shape)


# define class for dataset
class ZDataset(Dataset):
    def __init__(self, z_data):
        self.data = z_data

    def __len__(self):
        return (self.data.shape[0] - 1) # the sataset consists of samples and targets, but there are only n-1 targets

    def __getitem__(self, idx):
        sample = self.data[idx]
        target = self.data[idx + 1]
        return sample, target
    
# train data loader
z_train = (z_train_data)
z_train_loader = DataLoader(z_train, batch_size=batch_size_dnn, shuffle=True)

# test data loader
z_test = ZDataset(z_test_data)
z_test_loader = DataLoader(z_test, batch_size=batch_size_dnn, shuffle=True)


""" Define the DNN model """
# DNN Model
class DNN(nn.Module):
    def __init__(self, latent_dim):
        super(DNN, self).__init__()
        self.fc1 = nn.Linear(latent_dim, 128)
        self.fc2 = nn.Linear(128, 256)
        self.fc3 = nn.Linear(256, 256)
        self.fc4 = nn.Linear(256, 128)
        self.fc5 = nn.Linear(128, latent_dim)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        x = F.relu(self.fc4(x))
        return self.fc5(x)

# define loss function for DNN
def loss_function_dnn(output, z):
    return F.mse_loss(output, z, reduction='sum')


""" Train the DNN model """
os.chdir(model_dir) # go to model directory to save files

dnnmodel = DNN(latent_dim)
optimizer = optim.Adam(dnnmodel.parameters(), lr=lr_dnn)

# create loss list to visualize loss
losses = []

# Training loop
for epoch in range(num_epochs_dnn):
    dnnmodel.train()
    train_loss = 0
    for i, batch in enumerate(z_train_loader): # get sample and target batches
        sample, target = batch
        optimizer.zero_grad()
        output = dnnmodel(sample)
        loss = loss_function_dnn(output, target)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()
    
    avg_loss = train_loss / len(z_train_loader.dataset)
    # open file to write loss
    with open('DNN_loss_{}.txt'.format(ID), 'a') as f:
        f.write('{},{}\n'.format(epoch, avg_loss))
    losses.append(avg_loss)
    
# plot loss and save image
plt.figure(figsize=(10, 8))
plt.plot(losses, color='blue')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('DNN Loss')
plt.grid()
plt.savefig('DNN_loss{}.png'.format(ID))

# save the model
torch.save(dnnmodel.state_dict(), 'DNN_model_{}.pth'.format(ID))


""" Define the GAIT model"""
# combine the CVAE and DNN models
class GAIT(nn.Module):
    def __init__(self, cvae, dnn):
        super(GAIT, self).__init__()
        self.cvae = cvae
        self.dnn = dnn

    def forward(self, x):
        mu, logvar = self.cvae.encoder(x)
        z = self.cvae.reparameterize(mu, logvar)
        z_pred = self.dnn(z)
        return self.cvae.decoder(z_pred), mu, logvar, z
    
# initialize the GAIT model
gait = GAIT(model, dnnmodel)

# make prediction
first_batch = next(iter(phi_test_loader))

it = 20
gait.eval()
with torch.no_grad():
    pred, mu, logvar, z = gait(first_batch)

    # compare original and predicted on same subplot
    fig, ax = plt.subplots(1, 3)
    # set size to 10x5
    fig.set_size_inches(10, 5)
    ax[0].imshow(first_batch[it].squeeze(), cmap='jet')
    ax[0].set_title('Original (i)')
    ax[1].imshow(pred[it].squeeze().detach().numpy(), cmap='jet')
    ax[1].set_title('Predicted (i+1)')
    ax[2].imshow(first_batch[it+1].squeeze(), cmap='jet')
    ax[2].set_title('Actual (i+1)')	
    # set tight layout
    plt.tight_layout()
    plt.savefig('GAIT_prediction{}.png'.format(ID))
    
# save the model
torch.save(gait.state_dict(), 'GAIT_model_{}.pth'.format(ID))