import pandas as pd
import numpy as np
from torch.utils.data import DataLoader

class BRCADataset(Dataset):
    """Dataset for one of the BRCA genes"""

    def __init__(self, feather_file: str):
        """ 
            feather_file (str): Path to the feather file with the dataset
        """
        self.sequences = pd.read_feather(feather_file)    

    def __len__(self):
        pass 

    def __getitem__(self, idx):
        pass

