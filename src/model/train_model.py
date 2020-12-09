import argparse
import torch
from torch.utils.data import DataLoader

from lib.architectures import Generator, Discriminator 
from lib.dataset import BRCADataset 
from lib.optim import compute_gradient_penalty

cuda = True if torch.cuda.is_available() else False
Tensor = torch.cuda.FloatTensor if cuda else torch.FloatTensor

def train_model(input_sequences: str, epochs: int, batch_size: int, lr: float, lambda_gp: float):

    dataset = BRCADataset(input_sequences)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    generator = Generator(seq_len=dataset.seq_len, batch_size=batch_size)
    discriminator = Discriminator(seq_len=dataset.seq_len, batch_size=batch_size)
    optimizer_G = torch.optim.AdamW(
        generator.parameters(), 
        lr=lr, 
        lr_decay=0.5, 
        weight_decay=0.9,
    )
    optimizer_D = torch.optim.AdamW(
        discriminator.parameters(), 
        lr=lr,
        lr_decay=0.5,
        weight_decay=0.9,
    )   

    if cuda:
        generator.cuda()
        discriminator.cuda()

    for k in range(epochs):
        for i, seq in enumerate(dataloader):
            optimizer_G.zero_grad()            
            
            z = torch.autograd.variable(Tensor(torch.randn(100)))
            fake_seq = generator(z)

            real_score = discriminator(seq)
            fake_score = discriminator(fake_seq)
            gp = compute_gradient_penalty(discriminator, seq, fake_seq)

            D_loss = -torch.mean(real_score) + torch.mean(fake_score) + lambda_gp * gp
            D_loss.backward()
            optimizer_D.step()

            if i % 5 == 0:
                # discriminator update
                fake_seq = generator(z)
                fake_score = discriminator(fake_seq)
                G_loss = -torch.mean(fake_score)
                G_loss.backward()
                optimizer_G.step()

                print(
                    "[epoch : {}] [sample : {} ] [ generator loss : {}] [ discriminator loss: {}]".format(
                        k, i, G_loss, D_loss
                ))

if __name__ == "__main__":
    pass