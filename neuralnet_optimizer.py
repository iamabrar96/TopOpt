import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import parameters
import common

from sklearn.preprocessing import MinMaxScaler


class Dataset():
    def __init__(self):
        super(Dataset,self).__init__()
        # read and define inputs
        self.data_name=pd.read_csv('Experiments_1_without_Error_Term.csv', sep=';', decimal=',') # small dataset
        self.x = self.data_name[['F_press','fs','h_BT']].values
        # normilize inputs between 0 and 1
        scaler_x = MinMaxScaler(feature_range=(0,1))
        self.x = scaler_x.fit_transform(self.x)
        self.x = torch.tensor(self.x, dtype=torch.float64)
        print(f'size input: {self.x.size()}, min input: {self.x.min()}, max input: {self.x.max()}')
        # read and define outputs
        disp_all = np.load("disp_all_81000001.npy") # small dataset
        self.y = np.reshape(disp_all, (disp_all.shape[0],-1))
        # normalize displacement results between -1 and 1
        scaler_y = MinMaxScaler(feature_range=(0,1))
        self.y = scaler_y.fit_transform(self.y)
        self.y = torch.tensor(self.y, dtype=torch.float64)        
        print(f'size output: {self.y.size()}, min output: {self.y.min()}, max output: {self.y.max()}')   
        self.n_samples = self.x.shape[0]
        print('size dataset:', self.n_samples)
        
    def __getitem__(self, index):
        return self.x[index],self.y[index]  
    
    def __len__(self):
        return self.n_samples

# class Model(nn.Module):
#     def __init__(self, in_size, s, outsize):
#         super(Model,self). __init__()

#         self.l1 = nn.Linear(in_size, s, bias=True)
#         self.relu = nn.ReLU()
#         self.l2 = nn.Linear(s, outsize, bias=True)

#     def forward(self, x):
#         out = self.l1(x)
#         out = self.relu(out)
#         out = self.l2(out)
#         return out

def train(model, train_loader, val_loader, optimizer, epochs, device, criterion):
    model.to(device)
    train_loss = []
    valid_loss = []
    n_total_steps = len(train_loader)
    
    for epoch in range(epochs):
        model.train()
        batch_loss = []
        for i,(x,y) in enumerate(train_loader):
            optimizer.zero_grad()
            x = x.to(device)
            y = y.to(device)
            yhat = model(x.float())
            loss = criterion(y,yhat)
            loss.backward()
            optimizer.step()
            batch_loss.append(torch.log10(loss).item())

            if (i+1) % 1 == 0:
                print (f'Epoch [{epoch+1}/{epochs}], Step [{i+1}/{n_total_steps}], Train-Step-Loss: {torch.log10(loss.item()):.4f}')

        train_loss.append(torch.Tensor(batch_loss).mean().item())

        with torch.no_grad():
            model.eval()
            batch_loss = []
            for x,y in val_loader:
                x = x.to(device)
                y = y.to(device)
                yhat = model(x.float())
                loss = criterion(y,yhat)
                batch_loss.append(torch.log10(loss).item())
            valid_loss.append(torch.Tensor(batch_loss).mean().item()) 
        print(f'Epoch [{epoch+1}/{epochs}], Train-Loss: {train_loss[epoch]:.4f}, Valid-Loss: {valid_loss[epoch]:.4f}')
    return train_loss, valid_loss


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

class Model(nn.Module):
    def __init__(self, in_size, n_hidden,size_hidden,outsize, p=0):
        
        super(My_Model,self). __init__()
        self.in_size = in_size
        self.outsize = outsize
        self.n_hidden = n_hidden
        self.size_hidden = size_hideen
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

def loss_top(y, yhat):
    loss = torch.sum(torch.div(Jelem,nn_rho**penal))/obj0; # compliance
    volConstraint =((torch.mean(nn_rho)/desiredVolumeFraction) - 1.0);
    return loss, volConstraint


if torch.cuda.is_available():
    torch.cuda.empty_cache()
    
torch.manual_seed(params.seed_item)
MSE = loss_top(y, yhat)
model = My_Model(parameters : params).to(device)
lrfder = 0.00609 #lrfinder(model)
optimizer = getattr(torch.optim, optimizer_name)(model.parameters(), lr=lrfder, weight_decay= 1e-3)

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