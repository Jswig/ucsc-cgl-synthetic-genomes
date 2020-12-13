import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset

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

    def __getitem__(self, idx):
        item = (pd
            .get_dummies(self.sequences.iloc[idx,:])
            .values 
            .astype(np.float32)
            .reshape((4,self.seq_len))
        )
        return torch.from_numpy(item)

