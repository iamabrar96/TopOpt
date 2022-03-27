import numpy as np
import torch
import torch.nn as nn
from common import GMSH_helper, Topology_viz, device_common

from sklearn.preprocessing import MinMaxScaler
from parameters import Parameters

from topopt import FE_solver

class NN_Optimizer:
    def __init__(self, params) -> None:
        self.params = params
        self.solver= FE_solver(params)
        self.helper= GMSH_helper(params)
        self.topology_visualizer= Topology_viz(params)

        self.device= device_common()
        self.model= TopOpt_NN(params).to(self.device)
        self.inputs= torch.from_numpy(self.helper.get_element_centers()).to(self.device)

        self.criterion= loss_fn
    def train(self):
        train_loss = []
        valid_loss = []

        for epoch in range(epochs):
            model.train()
            batch_loss = []

            optimizer.zero_grad()

            density_hat = model(self.inputs)

            loss = self.criterion(density,density_hat)

            loss.backward()
            optimizer.step()

            batch_loss.append(torch.log10(loss).item())

            if (i+1) % 1 == 0:
                print (f'Epoch [{epoch+1}/{epochs}], Step [{i+1}/{n_total_steps}], Train-Step-Loss: {torch.log10(loss.item()):.4f}')

            train_loss.append(torch.Tensor(batch_loss).mean().item())

            with torch.no_grad():
                model.eval()
                batch_loss = []
                
                for self.inputs,density in val_loader:
                    
                    self.inputs = self.inputs.to(device)
                    density = density.to(device)
                    density_hat = model(self.inputs.float())
                    
                    loss = self.criterion(density,density_hat)
                    
                    batch_loss.append(torch.log10(loss).item())
                    
                valid_loss.append(torch.Tensor(batch_loss).mean().item()) 
            print(f'Epoch [{epoch+1}/{epochs}], Train-Loss: {train_loss[epoch]:.4f}, Valid-Loss: {valid_loss[epoch]:.4f}')
        return train_loss, valid_loss





class TopOpt_NN(nn.Module):
    def __init__(self, params:Parameters):
        
        super(My_Model,self). __init__()
        self.in_size = in_size
        self.outsize = outsize
        self.n_hidden = n_hidden
        self.size_hidden = size_hidden
        self.hidden=nn.ModuleList()
        self.drop=nn.Dropout(p=p)
        self.Layers = [in_size, *[self.size_hidden for i in range(n_hidden)], outsize]
        print(self.Layers)
        self.norm = nn.BatchNorm1d(s)
        for input_size,output_size in zip(self.Layers, self.Layers[1:]):
            self.hidden.append(nn.Linear(input_size, output_size, bias=False))
            
    def forward(self, x):
        L=len(self.hidden)
        for(l,linear_transform) in zip(range(L),self.hidden):
            if l<L-1:
            
                x = linear_transform(x)
                x=torch.relu(self.norm(x))
                x=self.drop(x)
            else:
                x= linear_transform(x)
        return x


def loss_fn():
    loss = torch.sum(torch.div(Jelem,nn_rho**penal))/obj0; # compliance
    volConstraint =((torch.mean(nn_rho)/desiredVolumeFraction) - 1.0);
    return loss, volConstraint


#initialisation of parameters
seed_list = np.array([9,89,232,5688,87843,216848,9867985,65468935,135878968,9999999999])
epochs = 30
in_size = parameters.n_dim
outsize = parameters.out_dim
size_hidden = 520

if torch.cuda.is_available():
    torch.cuda.empty_cache()
    
for j, h in enumerate(range(1,8)):
    for k,seed in enumerate(seed_list):
        torch.manual_seed(seed)
        print("\n hidden layer = {} seed = {}" .format(h,seed))
        MSE = nn.MSELoss()
        model = My_Model(in_size, n_hidden, size_hidden, outsize,  p=0.30).to(device)
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