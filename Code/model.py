# -*- coding: utf-8 -*-
# ************************************************************ #
#   FileName      : model.py
#   Author        : Mengchen
#   Email         : mengchenpu@gmail.com
#   Create on     : 06-21-2022
#   Last modified : 06-21-2022 17:14
#   Version       : V1.0
#   Description   : 
# ************************************************************ #
import torch
import torchvision
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
from torch.utils.data import Dataset

class ANNmodel(nn.Module):
    def __init__(self, input_dim=200, hidden_dim=128, output_dim=2):
        super(ANNmodel, self).__init__()
 
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.relu1 = nn.ReLU()
 
        self.fc2 = nn.Linear(hidden_dim, hidden_dim)
        self.relu2 = nn.ReLU()
 
        self.fc3 = nn.Linear(hidden_dim, output_dim)
 
    def forward(self, x):
        out = self.fc1(x)
        out = self.relu1(out)
        out = self.fc2(out)
        out = self.relu2(out)
        out = self.fc3(out)
        return out

class ANNmodel2(nn.Module):
    def __init__(self, input_dim=200, hidden_dim=128, output_dim=2):
        super(ANNmodel, self).__init__()

        self.fc1 = nn.Linear(input_dim, 784)
        self.drop1 = nn.Dropout(0.6)
        self.relu1 = nn.ReLU()

        self.fc2 = nn.Linear(784, 512)
        self.drop2 = nn.Dropout(0.6)
        self.relu2 = nn.ReLU()

        self.fc3 = nn.Linear(512, hidden_dim)
        self.drop3 = nn.Dropout(0.6)
        self.relu3 = nn.ReLU()

        self.fc4 = nn.Linear(hidden_dim, hidden_dim)
        self.drop4 = nn.Dropout(0.6)
        self.relu4 = nn.ReLU()

        self.fc5 = nn.Linear(hidden_dim, hidden_dim)
        self.drop5 = nn.Dropout(0.4)
        self.relu5 = nn.ReLU()

        self.fc6 = nn.Linear(hidden_dim, hidden_dim)
        self.drop6 = nn.Dropout(0.2)
        self.relu6 = nn.ReLU()
       
        self.fc7 = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu1(self.drop1(out))
        out = self.fc2(out)
        out = self.relu2(self.drop2(out))
        out = self.fc3(out)
        out = self.relu3(self.drop3(out))
        out = self.fc4(out)
        out = self.relu4(self.drop4(out))
        out = self.fc5(out)
        out = self.relu5(self.drop5(out))
        out = self.fc6(out)
        out = self.relu6(self.drop6(out))
        out = self.fc7(out)
        return out
