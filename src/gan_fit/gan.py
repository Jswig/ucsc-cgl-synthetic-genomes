import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
class ResidualBlock(nn.Module):
	def __init__(self, d_in: int, d_out: int, seq_len: int, relu: nn.Module): #pick nn.ReLU for generator, nn.LeakyRelu for discriminator
		super().__init__()
		self.lnorm_1 = nn.LayerNorm([d_in, seq_len])
		self.relu_1 = relu()
		self.conv_1 = nn.Conv1d(in_channels=d_in, out_channels=d_out, kernel_size=5, padding = 2)
		self.relu_2 = relu()
		self.lnorm_2 = nn.LayerNorm([d_in, seq_len])
		self.conv_2 = nn.Conv1d(in_channels=d_in, out_channels=d_out, kernel_size=5, padding = 2)
		# bias not needed, already present by default in Conv1d()
		
	def forward(self, x: torch.Tensor):
		x_0 = x
		y = self.lnorm_1(x_0)
		x = self.relu_1(x)
		x = self.conv_1(x)
		y = self.lnorm_1(y)
		x = self.relu_2(x)
		x = self.conv_2(x)
		return x_0 + (0.3 * x)


class Generator(nn.Module):
	def __init__(self, seq_len: int, batch_size: int, latent_dim: int=100):
		super(Generator, self).__init__()
		# layers
		self.lnorm_in = nn.LayerNorm(latent_dim)
		self.linear_in = nn.Linear(in_features=latent_dim, out_features=latent_dim*seq_len)
		self.resblocks = nn.Sequential(
			*[ResidualBlock(
				latent_dim, latent_dim, seq_len=seq_len, relu=nn.LeakyReLU
			) for _ in range(5)]
		)
		self.conv_out = nn.Conv1d(in_channels=latent_dim, out_channels=4, kernel_size=1)
		# hyper-parameters
		self.seq_len = seq_len
		self.latent_dim = latent_dim
		self.batch_size = batch_size

	def forward(self, x: torch.Tensor):
		x = self.lnorm_in(x)
		x = self.linear_in(x)    
		x = torch.reshape(x, (self.batch_size, self.latent_dim, self.seq_len))
		x = self.resblocks(x)
		x = self.conv_out(x)
		y = F.softmax(x, dim=1)
		return y

		
class Discriminator(nn.Module):
	def __init__(self, seq_len: int, batch_size: int, latent_dim: int=100):
		super().__init__()
		# layers
		self.conv_in = nn.Conv1d(in_channels=4, out_channels=latent_dim, kernel_size=1)
		self.resblocks = nn.Sequential(
			*[ResidualBlock(
				latent_dim, latent_dim, seq_len=seq_len, relu=nn.LeakyReLU
			) for _ in range(5)]
		)
		self.lnorm_out = nn.LayerNorm(latent_dim*seq_len)
		self.linear_out = nn.Linear(in_features=latent_dim*seq_len, out_features=1)
		# hyper-parameter
		self.seq_len = seq_len
		self.latent_dim = latent_dim
		self.batch_size = batch_size
		
	def forward(self, x: torch.Tensor):
		x = self.conv_in(x)
		x = self.resblocks(x)
		x = torch.reshape(x, (self.batch_size, self.latent_dim*self.seq_len))
		x = self.lnorm_out(x)
		y = self.linear_out(x)
		return y

class BRCADataset(Dataset):
	"""Dataset for one of the BRCA genes"""

	def __init__(self, feather_file: str):
		""" 
			feather_file (str): Path to the feather file with the dataset
		"""
		self.sequences = pd.read_feather(feather_file)
		self.seq_len = self.sequences.shape[1]

	def __len__(self):
		return(len(self.sequences))

	def __getitem__(self, idx: int):
		item = (pd
			.get_dummies(self.sequences.iloc[idx,:])
			.values 
			.astype(np.float32)
			.reshape((4,self.seq_len))
		)
		return torch.from_numpy(item)


def compute_gradient_penalty(discriminator: nn.Module, real_samples: torch.Tensor, fake_samples: torch.Tensor):
	"""
	Calculates the gradient penalty loss for WGAN GP
	"""

	batch_size = real_samples.shape[0]
	# Random weight term for interpolation between real and fake samples
	alpha = torch.randn((batch_size, 1, 1), device=device)
	# Get random interpolation between real and fake samples
	interpolates = (alpha * real_samples + ((1 - alpha) * fake_samples)).requires_grad_(True)
	d_interpolates = discriminator(interpolates)
	fake = torch.autograd.Variable(torch.ones((batch_size, 1), device=device), requires_grad=False)
	# Get gradient w.r.t. interpolates
	gradients = torch.autograd.grad(
		outputs=d_interpolates,
		inputs=interpolates,
		grad_outputs=fake,
		create_graph=True,
		retain_graph=True,
		only_inputs=True,
	)[0]
	gradients = gradients.view(gradients.size(0), -1)
	gradient_penalty = ((gradients.norm(2, dim=1) - 1) ** 2).mean()
	return gradient_penalty