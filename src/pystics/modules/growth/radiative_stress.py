import numpy as np


def radiative_stress(df, coefb):
    '''
    This module computes the radiative factor affecting biomass production based on photoresponse curve defined by ebmax and coefb.
    It is defined as the ratio of the carbon assimilation from raint(i) on the carbon assimilation in optimal conditions (optimal radiation).
    '''

    def carbon_assimilation(r, ebmax):
        return (-coefb) * r**2 + ebmax * r
    
    # Optimal radiation
    df['optimal_radiation'] = - df['ebmax'] / (2 * (-coefb))
    
    # Radiative stress
    radiative_stress_list = np.zeros(len(df['raint']))
    radiative_stress_list[df.raint < 1e-3] = 0
    radiative_stress_list[df.raint >= 1e-3] = df.loc[df.raint >= 1e-3,:].apply(lambda row:carbon_assimilation(row['raint'], row['ebmax']), axis=1) / df.loc[df.raint >= 1e-3].apply(lambda row:carbon_assimilation(row['optimal_radiation'], row['ebmax']), axis=1)

    return radiative_stress_list
