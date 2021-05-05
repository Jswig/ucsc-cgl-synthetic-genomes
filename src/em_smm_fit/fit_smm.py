import allel
import argparse
import numpy as np
import tqdm
from numba import jit 

parser = argparse.ArgumentParser(
    description='fit symmetric mixture model'
)
parser.add_argument('input', help='Input VCF file')

if __name__ == '__main__':
    args = parser.parse_args()