import numpy as np

def isNaN(num):
    return num != num
    # raise ValueError('is nan')


def isinf(num):
    return num == np.inf
    # raise ValueError('is nan')
