import numpy as np

def tvar(x):
    """
    This function computes the saturation vapor pressure from air temperature.
    """
    return 6.1070 * (1 + 2 ** (1 / 2) * np.sin(0.017453293 * x / 3)) ** (8.827)

def potential_etp(trg, temp, wind, tpm, gamma, codeetp, alphapt, fracinsol):
    '''
    This modules commputes potential evapotranspiration (etp) according to Penman equation (codeetp = 2) or Priestley & Taylor (codeetp = 4).
    '''

    deltat = tvar(temp + 0.5) - tvar(temp - 0.5)

    if (codeetp == 2):
        # Penman calculation
        dsat = tvar(temp) - tpm
        L = (2500840 - 2358.6 * temp) * 1e-6
        rglo = 4.9e-9 * ((temp + 273.16) ** 4)* (0.1 + 0.9 * fracinsol) * (0.56 - 0.08 * (tpm ** (1/2)))
        rnet_PE = (1 - 0.2)  * trg - rglo
        etp =  rnet_PE / L * deltat / (deltat + gamma) + (gamma / (deltat + gamma))  * 0.26 * (1 + 0.54 * wind) * dsat
        etp = max(etp,0)

    elif (codeetp == 4):
        # Priestley & Taylor calculation
        rnetpt = 0.8 * 0.72 * trg - 0.9504
        etp = alphapt * deltat * rnetpt / (2.5* (deltat + gamma))

    return etp