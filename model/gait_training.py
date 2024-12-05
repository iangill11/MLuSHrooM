# general imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
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
ID = 2

# make directory for saving files
dir_name = 'CVAE_{}'.format(ID)

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
training_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CVAE Training")
model_dir = os.path.join(training_dir, dir_name)

# change directory to script directory
os.chdir(script_dir)


""" Set Hyperparameters """
# model training
lr_dnn = 0.001
num_epochs_dnn = 1

# model data
batch_size_cvae = 32
batch_size_dnn = 32
input_dim = 64 # for square images
latent_dim = 64

""" Manual Parameters (Careful with these)"""
N_train_batches = 38


""" Load Data """
# load data
phi_data = torch.load('downscaled_data{}old.pt'.format(input_dim))

# get training and test data lengths
total_len = phi_data.shape[0]
train_len = N_train_batches * batch_size_cvae
test_len = total_len - train_len
# print lengths
print("Total data length: ", total_len)
print("Training data length: ", train_len)
print("Test data length: ", test_len)

# print data shape
print("Data shape: ", phi_data.shape)


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
#model.load_state_dict(torch.load('CVAE Training/CVAE_{}/CVAE_model_{}.pth'.format(ID, ID)))
model.load_state_dict(torch.load('CVAE Training/CVAE_3/CVAE_model_3.pth'))

""" Prepare Test and Train data for the DNN model """
# initialize empty tensors to store z-vectors
z_train_data = torch.empty(train_len, latent_dim)
z_test_data = torch.empty(test_len, latent_dim)

# split phi_train and phi_test data
phi_train = torch.tensor(phi_data[:train_len])
phi_test = torch.tensor(phi_data[train_len:])

print("phi_train shape: ", phi_train.shape)
print("phi_test shape: ", phi_test.shape)

model.eval()
with torch.no_grad():
    for i in range(train_len):
        input = phi_train[i].unsqueeze(0).unsqueeze(0)
        z_train_data[i] = model.encoder(input)[0]
        
    for i in range(test_len):
        input = phi_test[i].unsqueeze(0).unsqueeze(0)
        z_test_data[i] = model.encoder(input)[0]

# print shapes of data
print("Training data shape: ", z_train_data.shape)
print("Test data shape: ", z_test_data.shape)


# print first 10 z-vectors
print("First 10 z-vectors: ", z_train_data[:10])
 
z_train_numpy = z_train_data.numpy()
z_test_numpy = z_test_data.numpy()


# define class for dataset
class ZDataset(Dataset):
    def __init__(self, z_data):
        self.data = z_data

    def __len__(self):
        return (self.data.shape[0] - 1) # the dataset consists of samples and targets, but there are only n-1 targets

    def __getitem__(self, idx):
        sample = self.data[idx, :]
        target = self.data[idx + 1, :]
        return sample, target
    
# train data loader
z_train = ZDataset(z_train_data)
z_train_loader = DataLoader(z_train, batch_size=batch_size_dnn, shuffle=True)

# test data loader
z_test = ZDataset(z_test_data)
z_test_loader = DataLoader(z_test, batch_size=batch_size_dnn,  shuffle=True)


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
        
        # # get one sample vector and target vector
        # if i == 0:
        #     sample_arr = sample[:3].detach().numpy()
        #     target_arr = target[:3].detach().numpy()
        #     output_arr = output[:3].detach().numpy()
            
        #     print("Sample vector: ", sample_arr[0, 1:10])
        #     print("Target vector: ", target_arr[0, 1:10])
        #     print("Output vector: ", output_arr[0, 1:10])

    
    avg_loss = train_loss / len(z_train_loader.dataset)
    # open file to write loss
    with open('DNN_loss_{}.txt'.format(ID), 'a') as f:
        f.write('{},{}\n'.format(epoch, avg_loss))
    losses.append(avg_loss)
    
# plot loss and save image
plt.figure(figsize=(10, 8))
plt.plot(losses[10:], color='blue')
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
    
# save the model
torch.save(gait.state_dict(), 'GAIT_model_{}.pth'.format(ID))