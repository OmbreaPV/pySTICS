import numpy as np


def soil_temperature(temp_min, tcultmax, tcult, depth, diftherm, tsol_i_prev):
    '''
    This module computes soil temperature for every soil layers (cm by cm).
    See Section 10.2 of STICS book.
    '''
    
    # Daily thermal amplitude
    amplsurf = tcultmax - temp_min

    # Daily thermal amplitude for each soil layer
    amplz_i = np.empty(depth)
    amplz_i[:] = np.nan
    for z in range(depth):
        amplz_i[z] = amplsurf * np.exp(-z * (7.272*0.00001 / (2*diftherm))**(1/2))

    # Soil temperature for each soil layer
    tsol_i = np.empty(depth)
    tsol_i[:] = np.nan
    for z in range(depth):
        tsol_i[z] = tsol_i_prev[z] - amplz_i[z] * (tcult - temp_min) / amplsurf + 0.1 * (tcult - tsol_i_prev[z]) + amplz_i[z] / 2 

    return amplsurf, amplz_i, tsol_i