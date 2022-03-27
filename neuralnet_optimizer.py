import numpy as np
import torch
import torch.nn as nn
from common import GMSH_helper, Topology_viz, device_common, loss_fn
from parameters import Parameters

from topopt import FE_solver

class NN_Optimizer:
    def __init__(self, params: Parameters) -> None:
        self.params = params
        self.helper= GMSH_helper(params)
        self.topology_visualizer= Topology_viz(params)
        self.solver= FE_solver(params)

        self.device= device_common()
        self.model= TopOpt_NN(params).to(self.device)
        self.inputs= torch.from_numpy(self.helper.get_element_centers()).to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.params.lrfder)
        self.convergenceHistory = []
        
    def train(self):

        train_loss = []
        alphaMax = 100*self.params.volfrac
        alpha = params.alphaIncrement
        for epoch in range(params.epochs):
            self.optimizer.zero_grad()
            density_predicted = self.model(self.inputs).to(self.device)
            self.model.train()
            density_predicted_FE = density_predicted.cpu().detach().numpy()
            U, Jelem, _ = self.solver.solve(density_predicted)
            if(epoch == 0):
                self.inital_obj = ( self.params.E0*(density_predicted_FE**self.params.p)*Jelem).sum()
            # For sensitivity analysis, exponentiate by 2p here and divide by p in the loss func hence getting -ve sign
            Jelem = torch.tensor(self.params.E0*(density_predicted_FE**(2*self.params.p))*Jelem).view(-1).float().to(self.device)
            loss, volConstraint = self.loss_fn()
            currentVolumeFraction = np.average(density_predicted_FE)
            loss = self.loss+ alpha*pow(volConstraint,2)
            alpha = min(alphaMax, alpha + self.params.alphaIncrement)
            loss.backward()
            self.optimizer.step()
            self.convergenceHistory.append([loss.item(), currentVolumeFraction])
            self.params.p = min(4.0,self.params.p + 0.01)
            if(epoch % 20 == 0):
                self.topology_visualizer.add_view(density_predicted)
                print("{:3d} J: {:.3F}; loss: {:.3F}; ".format(epoch, currentVolumeFraction,loss.item()))

            self.topology_visualizer.add_view(density_predicted)
            print("{:3d} J: {:.3F}; loss: {:.3F}; ".format(epoch, currentVolumeFraction,loss.item()))

            train_loss.append(torch.Tensor(loss).mean().item())

        return self.convergenceHistory


def loss_fn(self):
    loss = torch.sum(torch.div(self.Jelem,self.density_predicted**self.params.p))/self.inital_obj
    volConstraint =((torch.mean(self.density_predicted)/self.params.volfrac) - 1.0);
    return loss, volConstraint


class TopOpt_NN(nn.Module):
    def __init__(self, params:Parameters):
        
        super(TopOpt_NN,self). __init__()
        self.params = params
        self.hidden=nn.ModuleList()
        self.drop=nn.Dropout(p=params.dropout)
        self.Layers = [params.n_dim, *[self.params.size_hidden for i in range(params.num_hidden_layer)], params.out_dim]
        print(self.Layers)
        self.norm = nn.BatchNorm1d(self.params.size_hidden)
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



if torch.cuda.is_available():
    torch.cuda.empty_cache()
params = Parameters()
torch.manual_seed(params.seed_item)

train_obj = NN_Optimizer(params)
train_obj.train()

del train_obj
torch.cuda.empty_cache()

torch.cuda.synchronize()
print("time taken to run this program is {:0.4f} seconds" .format(start.elapsed_time(end)))
torch.cuda.empty_cache()