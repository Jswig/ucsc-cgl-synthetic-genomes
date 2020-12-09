import torch
import torch.nn as nn
import torch.nn.functional as F

class ResidualBlock(nn.Module):
    def __init__(self, d_in, d_out, relu): #pick nn.ReLU for generator, nn.LeakyRelu for discriminator
        super().__init__()
        self.relu_1 = relu()
        self.conv_1 = nn.Conv1d(in_channels=d_in, out_channels=d_out, kernel_size=5, padding = 2)
        self.relu_2 = relu()
        self.conv_2 = nn.Conv1d(in_channels=d_in, out_channels=d_out, kernel_size=5, padding = 2)
        # bias not needed, already present by default in Conv1d()
        
    def forward(self, x):
        y = self.relu_1(x)
        y = self.conv_1(y)
        y = self.relu_2(y)
        y = self.conv_2(y)
        return x + (0.3 * y)

class Generator(nn.Module):
    def __init__(self, seq_len, batch_size):
        super(Generator, self).__init__()
        self.linear_in = nn.Linear(in_features=100, out_features=100*seq_len)
        self.resblocks = []
        for _ in range(5):
            self.resblocks.append(ResidualBlock(100, 100, relu=nn.ReLU))
        self.conv_out = nn.Conv1d(in_channels=100, out_channels=4, kernel_size=1)
        
        self.batch_size = batch_size
        self.seq_len = seq_len
    
    def forward(self,z):
        y = self.linear_in(z)
        y = torch.reshape(y, (self.batch_size, 100, self.seq_len))
        for i in range(5):
            y = self.resblocks[i](y)
        y = self.conv_out(y)
        y = F.softmax(y, dim=1)
        return y

        
class Discriminator(nn.Module):
    def __init__(self, seq_len, batch_size):
        super().__init__()
        self.conv_in = nn.Conv1d(in_channels=4, out_channels=100, kernel_size=1)
        self.resblocks = []
        for _ in range(5):
            self.resblocks.append(ResidualBlock(100, 100, relu=nn.LeakyReLU))
        self.linear_out = nn.Linear(in_features=100*seq_len, out_features=1)
        
        self.seq_len = seq_len
        self.batch_size = batch_size
        
    def forward(self, x):
        y = self.conv_in(x)
        for i in range(5):
            y = self.resblocks[i](y)
        y = torch.reshape(y, (self.batch_size, 100*self.seq_len))
        y = self.linear_out(y)
        return y