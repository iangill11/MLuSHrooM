import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import time
from timeit import default_timer
from Utilities import *
import time

# ====================================
#  Calculate transient response
# ====================================  
class Transient(nn.Module):
    def __init__(self, in_channels, out_channels, modes11, modes12):
        super(Transient, self).__init__()

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes11 = modes11
        self.modes12 = modes12
        self.scale = (1 / (in_channels*out_channels))
        self.weights_pole1 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes11,  dtype=torch.cfloat))
        self.weights_pole2 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes12, dtype=torch.cfloat))
        self.weights_residue = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes11,  self.modes12, dtype=torch.cfloat))

    def output_PR(self, lambda1, lambda2, alpha, weights_pole1, weights_pole2, weights_residue):  
        term1=torch.div(1,torch.einsum("pbix,qbik->pqbixk",torch.sub(weights_pole1,lambda1),torch.sub(weights_pole2,lambda2)))
        Pk=torch.einsum("bixk,pqbixk->pqbixk",weights_residue,term1) 
        output_residue1=torch.einsum("biox,oxikpq->bkpq", alpha, Pk)
        return output_residue1    

    def forward(self, x):
        ## 1 Computing Transient part
        #Compute input poles and resudes by FFT
        x = x[:,:,::2,::2]
        tx=T[:,::2].cuda()
        ty=X[:,::2].cuda()
        dty=(ty[0,1]-ty[0,0]).item()  # location interval
        dtx=(tx[0,1]-tx[0,0]).item()  # time interval
        alpha = torch.fft.fft2(x, dim=[-2,-1])
        omega1=torch.fft.fftfreq(ty.shape[1], dty)*2*np.pi*1j   # location frequency
        omega2=torch.fft.fftfreq(tx.shape[1], dtx)*2*np.pi*1j   # time frequency
        omega1=omega1.unsqueeze(-1).unsqueeze(-1).unsqueeze(-1)
        omega2=omega2.unsqueeze(-1).unsqueeze(-1).unsqueeze(-1)
        lambda1=omega1.cuda()
        lambda2=omega2.cuda()    
        
        # Obtain output poles and residues for transient part and steady-state part
        output_residue1 = self.output_PR(lambda1, lambda2, alpha, self.weights_pole1, self.weights_pole2, self.weights_residue)
        
        # Obtain time histories of transient response and steady-state response
        tx_out=T.cuda()
        ty_out=X.cuda()
        term1=torch.einsum("bip,kz->bipz", self.weights_pole1, ty_out.type(torch.complex64).reshape(1,-1))
        term2=torch.einsum("biq,kx->biqx", self.weights_pole2, tx_out.type(torch.complex64).reshape(1,-1))
        term3=torch.einsum("bipz,biqx->bipqzx", torch.exp(term1),torch.exp(term2))
        x1=torch.einsum("kbpq,bipqzx->kizx", output_residue1,term3)
        x1=torch.real(x1)
        x1=x1/x.size(-1)/x.size(-2)
        return x1

 # ====================================
#  Calculate steady-state response, which is same to the Fourier layer
# ====================================     
class Steady(nn.Module):
    def __init__(self, in_channels, out_channels, modes21, modes22):
        super(Steady, self).__init__()

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes21 = modes21
        self.modes22 = modes22
        self.scale = (1 / (in_channels*out_channels))
        self.weights1 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes21, self.modes22, dtype=torch.cfloat))
        self.weights2 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes21, self.modes22, dtype=torch.cfloat))

    def forward(self, x):
        ## 2 Computing Steady-state part
        x = x[:,:,::2,::2]
        tx=T[:,::2].cuda()
        ty=X[:,::2].cuda()
        dty=(ty[0,1]-ty[0,0]).item()  # location interval
        dtx=(tx[0,1]-tx[0,0]).item()  # time interval
        FW1 = torch.fft.rfft2(x)
        omega1=torch.fft.fftfreq(ty.shape[1], dty)*2*np.pi*1j   # location frequency
        omega2=torch.fft.rfftfreq(tx.shape[1], dtx)*2*np.pi*1j   # time frequency
        lambda1=omega1.cuda()
        lambda2=omega2.cuda() 

        out_2 = torch.zeros(x.shape[0], self.out_channels,  x.size(-2), x.size(-1)//2 + 1, dtype=torch.cfloat, device=x.device)
        out_2[:, :, :self.modes21, :self.modes22] = \
            torch.einsum("bixy,ioxy->boxy", FW1[:, :, :self.modes21, :self.modes22], self.weights1)
        out_2[:, :, -self.modes21:, :self.modes22] = \
            torch.einsum("bixy,ioxy->boxy", FW1[:, :, -self.modes21:, :self.modes22], self.weights2)
        
        #Return to physical space
        tx_out=T.cuda()
        ty_out=X.cuda()
        term4=torch.einsum("xk,kz->xz", lambda1.reshape(-1,1), ty_out.type(torch.complex64).reshape(1,-1))
        term5=torch.einsum("bk,ki->bi", lambda2.reshape(-1,1), tx_out.type(torch.complex64).reshape(1,-1))
        term6=torch.einsum("xz,bi->xbzi", torch.exp(term4),torch.exp(term5))
        x2=torch.einsum("pqxb,xbzi->pqzi", out_2,term6)
        x2=torch.real(x2)
        x2=x2/x.size(-1)/x.size(-2)
        return x2

class LNO2d(nn.Module):
    def __init__(self, width1,width2,width3,modes11,modes12,modes21,modes22):
        super(LNO2d, self).__init__()

        self.width1 = width1
        self.width2 = width2
        self.width3 = width3
        self.modes11 = modes11
        self.modes12 = modes12
        self.modes21 = modes21
        self.modes22 = modes22

        self.fc2 = nn.Linear(3, self.width2) 
        self.conv_s0 = Steady(self.width2, self.width2, modes21,modes22)
        self.conv_s1 = Steady(self.width2, self.width2, modes21,modes22)
        self.conv_s2 = Steady(self.width2, self.width2, modes21,modes22)
        self.conv_s3 = Steady(self.width2, self.width2, modes21,modes22)
        self.w0 = nn.Conv2d(self.width3, self.width3, 1)
        self.w1 = nn.Conv2d(self.width3, self.width3, 1)
        self.w2 = nn.Conv2d(self.width3, self.width3, 1)
        self.w3 = nn.Conv2d(self.width3, self.width3, 1)
        self.norm2 = nn.InstanceNorm2d(self.width2)
        self.fc5 = nn.Linear(self.width2, 128)
        self.fc8 = nn.Linear(128, 1)
        self.fc1 = nn.Linear(3, self.width1) 
        self.conv_t0 = Transient(self.width1, self.width1, modes11,modes12)
        self.conv_t1 = Transient(self.width1, self.width1, modes11,modes12)
        self.conv_t2 = Transient(self.width1, self.width1, modes11,modes12)
        self.conv_t3 = Transient(self.width1, self.width1, modes11,modes12)
        self.norm1 = nn.InstanceNorm2d(self.width1)
        self.fc4 = nn.Linear(self.width1, 128)
        self.fc7 = nn.Linear(128, 1)


# ====================================
#  The total output includes three parts: transient part, steady-state part and W part. 
#  One can choose how to combine them in this part according to the problem
# ====================================  
    def forward(self,f):
        grid = self.get_grid(f.shape, f.device)
        f = torch.cat((f, grid), dim=-1)

        # # # # # Transient part
        f1 = self.fc1(f)
        f1 = f1.permute(0, 3, 1, 2)
        x1 = self.norm1(self.conv_t0(self.norm1(f1)))
        x1 =  torch.sin(x1)
        x1 = self.norm1(self.conv_t1(self.norm1(x1)))
        x1 =  torch.sin(x1)
        x1 = self.norm1(self.conv_t2(self.norm1(x1)))
        x1 =  torch.sin(x1)
        x1 = self.norm1(self.conv_t3(self.norm1(x1)))
        x1 = x1.permute(0, 2, 3, 1)
        x1 = self.fc4(x1)
        x1 =  torch.sin(x1)
        x1 = self.fc7(x1)

        # Steady-state part + W part
        f2 = self.fc2(f)
        f2 = f2.permute(0, 3, 1, 2)
        x2 = self.norm2(self.conv_s0(self.norm2(f2)))
        x23=x2
        x23 =  torch.sin(x23)

        x2 = self.norm2(self.conv_s1(self.norm2(x23)))
        x3 = self.w1(x23)
        x23=x2+x3
        x23 =  torch.sin(x23)

        x2 = self.norm2(self.conv_s2(self.norm2(x23)))
        x3 = self.w2(x23)
        x23=x2+x3
        x23 =  torch.sin(x23)

        x2 = self.norm2(self.conv_s3(self.norm2(x23)))
        x3 = self.w3(x23)
        x23=x2+x3

        x23 = x23.permute(0, 2, 3, 1)
        x23 = self.fc5(x23)
        x23 =  torch.sin(x23)
        x23 = self.fc8(x23)

        # # # Transient part+W part
        # f1 = self.fc1(f)
        # f1 = f1.permute(0, 2, 1)
        # x1 = self.conv_t0(f1)
        # x3 = self.w0(f1)
        # x13=x1+x3
        
        # x13 = x13.permute(0, 2, 1)
        # x13 = self.fc4(x13)
        # x13 =  torch.sin(x13)
        # x13 = self.fc7(x13)

        # # Steady-state part
        # f2 = self.fc2(f)
        # f2 = f2.permute(0, 2, 1)
        # x2 = self.conv_s0(f2)
        # x2 =  torch.sin(x2)

        # x2 = self.conv_s1(x2)

        # x2 = x2.permute(0, 2, 1)
        # x2 = self.fc5(x2)
        # x2 =  torch.sin(x2)
        # x2 = self.fc7(x2)

        # # W part
        # f3 = self.fc3(f)
        # f3 = f3.permute(0, 2, 1)
        # x3 = self.w0(f3)
        # x3 =  torch.sin(x3)

        # x3 = self.w1(x3)

        # x3 = x3.permute(0, 2, 1)
        # x3 = self.fc6(x3)
        # x3 =  torch.sin(x3)
        # x3 = self.fc7(x3)

        #return x23
        return x1+x23

    def get_grid(self, shape, device):
        batchsize, size_x, size_y = shape[0], shape[1], shape[2]
        gridx = torch.tensor(np.linspace(0, 1, size_x), dtype=torch.float)
        gridx = gridx.reshape(1, size_x, 1, 1).repeat([batchsize, 1, size_y, 1])
        gridy = torch.tensor(np.linspace(0, 1, size_y), dtype=torch.float)
        gridy = gridy.reshape(1, 1, size_y, 1).repeat([batchsize, size_x, 1, 1])
        return torch.cat((gridx, gridy), dim=-1).to(device)
    
 
# ====================================
#  Define parameters and Load data
# ====================================
     
ntrain = 800
nvali = 100
ntest = 100

batch_size_train = 10
batch_size_vali = 10

learning_rate = 0.01
epochs = 1000
step_size = 100
gamma = 0.5

modes11 =4
modes12 = 4
modes21 = 4
modes22 = 4
width1 = 16
width2 = 16
width3 = 16

reader = MatReader('Data/data.mat')
T = reader.read_field('t')
X = reader.read_field('x')

x_train = reader.read_field('f_train')
y_train = reader.read_field('u_train')
for idx1 in range(x_train.shape[1]):
    for idx2 in range(x_train.shape[2]):
        if idx1 % 2 != 0 and idx2 % 2 != 0:
            x_train[:, idx1, idx2] = 0

x_vali = reader.read_field('f_vali')
y_vali = reader.read_field('u_vali')
for idx1 in range(x_vali.shape[1]):
    for idx2 in range(x_vali.shape[2]):
        if idx1 % 2 != 0 and idx2 % 2 != 0:
            x_vali[:, idx1, idx2] = 0

x_test = reader.read_field('f_test')
y_test = reader.read_field('u_test')
for idx1 in range(x_test.shape[1]):
    for idx2 in range(x_test.shape[2]):
        if idx1 % 2 != 0 and idx2 % 2 != 0:
            x_test[:, idx1, idx2] = 0

x_train = x_train.reshape(ntrain,x_train.shape[1],x_train.shape[2],1)
x_vali = x_vali.reshape(nvali,x_vali.shape[1],x_vali.shape[2],1)
x_test = x_test.reshape(ntest,x_test.shape[1],x_test.shape[2],1)

train_loader = torch.utils.data.DataLoader(torch.utils.data.TensorDataset(x_train, y_train), batch_size=batch_size_train, shuffle=True)
vali_loader = torch.utils.data.DataLoader(torch.utils.data.TensorDataset(x_vali, y_vali), batch_size=batch_size_vali, shuffle=True)
# model
model = LNO2d(width1,width2,width3,modes11,modes12,modes21,modes22).cuda()
