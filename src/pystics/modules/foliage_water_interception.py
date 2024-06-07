import numpy as np

def foliage_water(codlocirrig, stemflowmax, kstemflow, airg, trr, mouillabil, lai_prev):
    '''
    This module water intercpted by foliage.
    '''

    # Water runoff along the stem

    if codlocirrig == 1:
        stemflow = stemflowmax * (trr + airg) * (1 - np.exp(-kstemflow * lai_prev))
    else:
        stemflow = stemflowmax * trr * (1 - np.exp(-kstemflow * lai_prev))

    # Water intercepted by foliage
    if codlocirrig== 1:
        mouill = (trr + airg - stemflow) * mouillabil
    else:
        mouill = (trr - stemflow) * mouillabil

    return stemflow, mouill