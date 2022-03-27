import os
import torch
import torch.nn as nn
import torchvision 
from torch.nn import Linear
from torch import optim
import pandas as pd
import numpy as np
from torch.utils.data.dataloader import DataLoader
from torch.utils.data import random_split
import matplotlib.pyplot as plt
from torchvision.transforms import ToTensor
from torch.utils.data.dataloader import DataLoader
from torch.utils.data import random_split
from mpl_toolkits import mplot3d
torch.set_printoptions(threshold=50,edgeitems=300,linewidth=100)
import math
import statistics 
import optuna
import plotly
torch.set_default_dtype(torch.float64)
import torchvision.transforms
torch.cuda.empty_cache()
import jovian

start = torch.cuda.Event(enable_timing=True)
end = torch.cuda.Event(enable_timing=True)
start.record() 



def device_common():
    if torch.cuda.is_available():
        device = torch.device("cuda:0")
        print('device is assigned as CUDA')
        torch.cuda.empty_cache()

    else:
        device = torch.device("cpu")   
        print('device is assigned as CPU')
    
    torch.cuda.device_count()
    return device

random.seed(manualSeed)

class My_Model(nn.Module):
    def __init__(self, in_size, h,s,outsize, p=0):
        
        super(My_Model,self). __init__()
        self.in_size = in_size
        self.outsize = outsize
        self.h = h
        self.s = s
        self.hidden=nn.ModuleList()
        self.drop=nn.Dropout(p=p)
        self.Layers = [in_size, *[self.s for i in range(h)], outsize]
        print(self.Layers)
        self.norm = nn.BatchNorm1d(s)
        for input_size,output_size in zip(self.Layers, self.Layers[1:]):
            self.hidden.append(nn.Linear(input_size, output_size, bias=False))
            
    def forward(self, activation):
        L=len(self.hidden)
        for(l,linear_transform) in zip(range(L),self.hidden):
            if l<L-1:
            
                activation = linear_transform(activation)
                activation=torch.relu(self.norm(activation))
                activation=self.drop(activation)
            else:
                activation= linear_transform(activation)
        return activation



def train(model,train_loader, val_loader, optimizer, epochs):
    model.to(device)
    total_loss = {'training_loss':[],'validation_loss':[]}
    for epoch in range(epochs):
        batch_loss=[]
        for i,(x,y) in enumerate(train_loader):
            x = x.to(device)
            y = y.to(device)
            yhat = model(x)
            loss = MSE(y,yhat)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            batch_loss.append(torch.sqrt(loss).item())
        total_loss['training_loss'].append(torch.Tensor(batch_loss).mean().item())
#         print(total_loss['training_loss'])
  
        with torch.no_grad():
            model.eval()
            batch_loss=[]        
            for x,y in val_loader:
                x = x.to(device)
                y = y.to(device)
                yhat = model(x)
                loss = MSE(y,yhat)
                batch_loss.append(torch.sqrt(loss).item())
            total_loss['validation_loss'].append(torch.Tensor(batch_loss).mean().item()) 
#             print(total_loss['validation_loss'])

    return total_loss


def loss_top():
    loss = torch.sum(torch.div(Jelem,nn_rho**penal))/obj0; # compliance
    volConstraint =((torch.mean(nn_rho)/desiredVolumeFraction) - 1.0);
    return loss, volConstraint




#initialisation of parameters
seed_list = np.array([9,89,232,5688,87843,216848,9867985,65468935,135878968,9999999999])
epochs = 30
in_size = 3 
outsize = 14943
s = 520

if torch.cuda.is_available():
    torch.cuda.empty_cache()
    
for j, h in enumerate(range(1,8)):
    for k,seed in enumerate(seed_list):
        torch.manual_seed(seed)
        print("\n hidden layer = {} seed = {}" .format(h,seed))
        MSE = nn.MSELoss()
        model = My_Model(in_size, h,s,outsize,  p=0.30).to(device)
        lrfder = 0.00609 #lrfinder(model)
        optimizer = optim.RMSprop(model.parameters(), lr=lrfder, weight_decay= 1e-3)
        total_loss= train(model, train_loader, val_loader, optimizer,epochs)
        for training_loss,validation_loss in zip(total_loss['training_loss'],total_loss['validation_loss']):
            print('training_loss:{:.2f},     validation_loss:{:.2f}'.format(training_loss,validation_loss))
        r2_score[j,k] = r2score(model)
        torch.save(model.state_dict(), '/media/sneha/E0DCBE20DCBDF140/Users/sneha/Downloads/files/Hemispherical/Model properties/model_{}_{}.ckpt'.format(h,seed))
        del model
        torch.cuda.empty_cache()

torch.cuda.synchronize()
print("time taken to run this program is {:0.4f} seconds" .format(start.elapsed_time(end)))
torch.cuda.empty_cache()