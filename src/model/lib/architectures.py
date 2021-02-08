class ResidualBlock(nn.Module):
    def __init__(self, d_in, d_out, relu): #pick nn.ReLU for generator, nn.LeakyRelu for discriminator
        super().__init__()
        self.bnorm_1 = nn.BatchNorm1d(d_in)
        self.relu_1 = relu()
        self.conv_1 = nn.Conv1d(in_channels=d_in, out_channels=d_out, kernel_size=5, padding = 2)
        self.relu_2 = relu()
        self.bnorm_2 = nn.BatchNorm1d(d_in)
        self.conv_2 = nn.Conv1d(in_channels=d_in, out_channels=d_out, kernel_size=5, padding = 2)
        # bias not needed, already present by default in Conv1d()
        
    def forward(self, x):
        y = self.bnorm_1(x)
        y = self.relu_1(y)
        y = self.conv_1(y)
        y = self.bnorm_1(y)
        y = self.relu_2(y)
        y = self.conv_2(y)
        return x + (0.3 * y)


class Generator(nn.Module):
    def __init__(self, seq_len, batch_size, latent_dim=100):
        super(Generator, self).__init__()
        # layers
        self.bnorm_in = nn.BatchNorm1d(latent_dim)
        self.linear_in = nn.Linear(in_features=latent_dim, out_features=latent_dim*seq_len)
        self.resblocks = nn.Sequential(
            *[ResidualBlock(latent_dim, latent_dim, relu=nn.LeakyReLU) for _ in range(5)]
        )
        self.conv_out = nn.Conv1d(in_channels=latent_dim, out_channels=4, kernel_size=1)
        # hyper-parameters
        self.seq_len = seq_len
        self.latent_dim = latent_dim
        self.batch_size = batch_size

    def forward(self,z):
        y = self.bnorm_in(z)
        y = self.linear_in(y)    
        y = torch.reshape(y, (self.batch_size, self.latent_dim, self.seq_len))
        y = self.resblocks(y)
        y = self.conv_out(y)
        y = F.softmax(y, dim=1)
        return y

        
class Discriminator(nn.Module):
    def __init__(self, seq_len, batch_size, latent_dim=100):
        super().__init__()
        # layers
        self.conv_in = nn.Conv1d(in_channels=4, out_channels=latent_dim, kernel_size=1)
        self.resblocks = nn.Sequential(
            *[ResidualBlock(latent_dim, latent_dim, relu=nn.LeakyReLU) for _ in range(5)]
        )
        self.bnorm_out = nn.BatchNorm1d(latent_dim*seq_len)
        self.linear_out = nn.Linear(in_features=latent_dim*seq_len, out_features=1)
        # hyper-parameters
        self.seq_len = seq_len
        self.latent_dim = latent_dim
        self.batch_size = batch_size
        
    def forward(self, x):
        y = self.conv_in(x)
        y = self.resblocks(y)
        y = torch.reshape(y, (self.batch_size, self.latent_dim*self.seq_len))
        y = self.bnorm_out(y)
        y = self.linear_out(y)
        return y