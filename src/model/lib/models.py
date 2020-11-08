import torch
import torch.nn as nn

class ResidualBlock(nn.module):
    def __init__(self, d_in, d_out):
        super().__init__()
        self.relu_1 = nn.ReLU()
        self.conv_1 = nn.Conv1d(in_channels=d_in, out_channels=d_out, kernel_size=5, padding=4)
        self.relu_2 = nn.ReLU()
        self.conv_2 = nn.Conv1d(in_channels=d_in, out_channels=d_out, kernel_size=5, padding=4)

    def forward(self, x):
        y = self.relu_1(x)
        y = self.conv_1(y)
        y = self.relu_2(y)
        y = self.conv_2(y)
        return x + (0.3 * y)