import argparse
import os
import torch
from torch.utils.data import DataLoader

from lib.architectures import Generator, Discriminator 
from lib.dataset import BRCADataset 
from lib.optim import compute_gradient_penalty

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')

def train(
    input_sequences: str,
    checkpt_path: str, 
    epochs: int, 
    batch_size: int, 
    lr: float, 
    lambda_gp: float,
):

    dataset = BRCADataset(input_sequences)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    generator = Generator(
        seq_len=dataset.seq_len, 
        batch_size=batch_size
    ).to(device)
    discriminator = Discriminator(
        seq_len=dataset.seq_len, 
        batch_size=batch_size
    ).to(device)
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
            
            z = torch.autograd.variable(torch.randn(100).to(device))
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
                    "[epoch : {}] [sample : {} ] [ generator loss : {}] [ discriminator loss: {}]".format(
                        k, i, G_loss, D_loss
                ))

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
    train(
        epochs=1,
        input_sequences='data/processed/brca2_seqs.feather',
        checkpt_path='output',
        batch_size=1,
        lr=0.0001,
        lambda_gp=10,
    )