import argparse
import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from gan import Generator, Discriminator, BRCADataset, compute_gradient_penalty


device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
argparse = argparse.ArgumentParser(
	description='Train GAN with given parameters'
)
argparse.add_argument('input', type=str, help='Input DNA sequences in .feather format')
argparse.add_argument('checkpoint', type=str, help='Checkpoint path')
argparse.add_argument(
	'--epochs', default=100, type=int, help='Number of training epochs' 
)
argparse.add_argument(
	'--batch_size', default=16, type=int, help='Training batch size'
)
argparse.add_argument(
	'--latent_dim', default=25, type=int, help='Dimension of latent space'
)
argparse.add_argument('--lr', default=0.001, type=float, help='Learning rate')
argparse.add_argument('--gp', default=10, type=floar, help='Gradient penalty')


def train(
	dataloader: DataLoader,
	generator: nn.Module,
	discriminator: nn.Module,
	checkpt_path: str, 
	epochs: int, 
	batch_size: int, 
	lr: float, 
	lambda_gp: float,
):
	dataset = BRCADataset(input_sequences)
	dataloader = DataLoader(
		dataset, 
		batch_size=batch_size, 
		shuffle=True,
		drop_last=True,
	)
	optimizer_G = torch.optim.AdamW(
		generator.parameters(), 
		lr=lr, 
		weight_decay=0.5,
	)
	optimizer_D = torch.optim.AdamW(
		discriminator.parameters(), 
		lr=lr,
		weight_decay=0.5,
	)   
	for k in range(epochs):
		for i, seq in enumerate(dataloader):
			seq = seq.to(device)

			optimizer_G.zero_grad()            
			
			z = torch.autograd.variable(torch.randn((batch_size, latent_dim)).to(device))
			fake_seq = generator(z)

			real_score = discriminator(seq)
			fake_score = discriminator(fake_seq)
			gp = compute_gradient_penalty(discriminator, seq, fake_seq)

			D_loss = -torch.mean(real_score) + torch.mean(fake_score) + lambda_gp * gp
			D_loss.backward()
			optimizer_D.step()

			if i % 5 == 0:
				# generator update
				fake_seq = generator(z)
				fake_score = discriminator(fake_seq)
				G_loss = -torch.mean(fake_score)
				G_loss.backward()
				optimizer_G.step()
				print(
					f"[epoch : {k}] [sample : {i} ] [ generator loss : {G_loss}] [ discriminator loss: {D_loss}]"
				)

		if i % 10 == 0:
			# checkpoint model
			torch.save({
				'epoch': k,
				'discriminator_state': discriminator.state_dict(),
				'generator_state': generator.state_dict(),  
				'optimizer_D_state': optimizer_D.state_dict(),
				'optimizer_G_state': optimizer_G.state_dict(),
				'D_loss': D_loss,
				'G_loss': G_loss,
			}, os.path.join(checkpt_path, 'checkpt_{}.pt'.format(k)))
			
			
	torch.save(generator.state_dict(), os.path.join(checkpt_path, 'generator_final.pt'))
	torch.save(discriminator.state_dict(), os.path.join(checkpt_path, 'discriminator_final.pt'))

if __name__ == "__main__":
	args = parser.parse_args()

	dataset = BRCADataset(args.input)
	dataloader = DataLoader(
		dataset, 
		batch_size=args.batch_size, 
		shuffle=True,
		drop_last=True,
	)
	generator = Generator(
		seq_len=dataset.seq_len, 
		batch_size=args.batch_size,
		latent_dim=args.latent_dim,
	).to(device)
	discriminator = Discriminator(
		seq_len=dataset.seq_len, 
		batch_size=args.batch_size,
		latent_dim=args.latent_dim,
	).to(device)

	train(
		dataloader=dataloader,
		generator=generator,
		discriminator=discriminator,
		checkpt_path=args.checkpoint,
		epochs=args.epochs,
		batch_size=args.batch_size,
		lr=args.lr,
		lambda_gp=args.gp,
	)