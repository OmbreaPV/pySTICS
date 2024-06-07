import numpy as np


def plant_height(hautmax, hautbase, lai, laisen):
    '''
    This module computes canopy height.
    '''

    hauteur = (
        hautmax * (1 - np.exp(-0.7 * (lai + laisen)))
        + hautbase
    ) 

    return hauteur