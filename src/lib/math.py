import numba
import numpy as np
from numba import jit 

@numba.jit
def logsum_approx(summands: np.array):
    n = len(summands)
    tot = 0 
    while(n > 0):
        tot += np.log(
            1 + np.exp(summands[n] - np.log(np.sum(summands[:n])))
        )
        n -= 1
    tot += np.log(summands[n])

    return tot 
