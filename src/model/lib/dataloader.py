import pandas as pd
import numpy as np
from torch.utils.data import Dataset, DataLoader

class BRCADataset(Dataset):
    """Dataset for one of the BRCA genes"""

    def __init__(self, feather_file: str):
        """ 
            feather_file (str): Path to the feather file with the dataset
        """
        self.sequences = pd.read_feather(feather_file)
        self.seq_length = self.sequences.shape[1]

    def __len__(self):
        return(len(self.sequences))

    def __getitem__(self, idx):
        item = pd.get_dummies(self.sequences.iloc[idx,:])
        return np.reshape(item.values, newshape=(1,4,len(item)))

