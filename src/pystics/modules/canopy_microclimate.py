import numpy as np
from pystics.exceptions import pysticsException

def radiation_interception(parsurrg, trg, extin, lai):
    '''
    This module computes intercepted radiation, based on Beer's law.
    See section 9.2.1.1 of STICS book.
    '''
    
    # Intercepted radiation
    raint = (
        0.95
        * parsurrg
        * trg
        * (1 - np.exp(-extin * lai))
    )

    return raint

def net_radiation(albedo, hur_i0, hminf_1, hccf_1, albveg, lai, trg, tcult, temp, tpm, fracinsol, codernet):
    '''
    This modules computes soil and vegetation albedos, long wave radiation and net radiation.
    See Section 9.3.1 of STICS book.

    Simplifications compared to STICS : 
        - No modelling of mulch.

    '''

    # Soil albedo
    albsol = albedo * (
        1
        - 0.517
        * (hur_i0 - hminf_1)
        / (hccf_1 - hminf_1)
    )
    albsolhum = 0.483 * albedo
    albsol = max(min(albsol, albedo), albsolhum)

    # Surface albedo = soil + crop
    albedolai = albveg - (albveg - albsol) * np.exp(-0.75 * lai)


    if codernet == 1: # Brunt method

        # Long wave radiations with STICS method (FA0 method is implemented but not used, see below)
        rglo = 4.903e-9 * ((tcult + 273) ** 4)* (0.1 + 0.9 * fracinsol) * (0.56 - 0.08 * (tpm ** 1 / 2))
        # rglo = -(
        #     0.000000004903
        #     * ((tcult + 273) ** 4)
        #     * (0.1 + 0.9 * fracinsol)
        #     * (0.34 - 0.14 * (tpm ** 1 / 2)))
    
    elif codernet == 2: # Brutsaert method

        eabrut = 1.24 * (tpm / (temp + 273.15))**(1.0 / 7.0)
        emissa = eabrut + (1.0 - fracinsol) * (1.0 - eabrut) * (1.0 - 4.0 * 11.0 / (273.15 + temp))
        ratm = 5.67e-8 * emissa * (temp + 273.15)**4.0
        ratm = ratm * 3600.0 *24.0 * 1e-6

        # Long wave radiations
        rsolglo = 5.67e-8 * (tcult + 273.15)**4
        rglo = -(ratm - (rsolglo * 3600.0 *24.0 *1e-6)) # j'ai rajouté un - pr pouvoir faire rnet en dehors du if


    # Net radiation
    rnet = (
        1 - albedolai
    ) * trg - rglo

    rnet = max(rnet, 0.01)

    
    return rnet, rglo, albedolai, albsol






def crop_temperature(lev_i_prev, temp_max, rnet, et, temp_min, temp, z0):
    ''''
    This module computes crop surface temperature with the empirical approach.
    See Section 9.3.2.1 of STICS book.
    '''

    if lev_i_prev > 0:

        # Max crop surface temperature
        tcultmax = temp_max + (
            rnet / 2.46 - et - 1.27
        ) / (1.68 / np.log(1 / z0))

        tcultmax = max(tcultmax, temp_max)

        # Crop temperature (tcultmin = tmin)
        tcult = (
            tcultmax + temp_min
        ) / 2

    else: # tcult before emergence does not make sense, but is necessary for soil temperature calculation
        tcult = temp
        tcultmax = temp_max

    return tcult, tcultmax


def wind_profile(zosolnu, hauteur):
    '''
    This module computes the characteristics of wind profile.
    See Section 9.3.2.1 of STICS book.
    '''
    
    # Crop roughness
    z0 = max(zosolnu, 0.13 * hauteur)
    
    # Displacement height
    dh = 6.6 * z0

    return dh, z0


def iterative_calculation(temp, lev_i_prev, temp_max, temp_min, et, z0, albedo, hur_i0, hminf_1, hccf_1, albveg, lai, trg, tpm, fracinsol, codernet, ):
    '''
    This module performs the iterative calculation of crop temperature and net radiation.
    See section 9.3.2.3 of STICS book.
    '''
    tcult_tmp = temp

    j = 0
    while j < 5:

        # Net radiation
        rnet, rglo, albedolai, albsol = net_radiation(albedo, hur_i0, hminf_1, hccf_1, albveg, lai, trg, tcult_tmp, temp, tpm, fracinsol, codernet)

        # Crop temperature
        tcult, tcultmax = crop_temperature(lev_i_prev, temp_max, rnet, et, temp_min, temp, z0)

        if np.isnan(tcult):
            raise pysticsException(
            detail="tcult is nan",
        )

        if tcult - tcult_tmp < 0.5:
            converge = True
            break
        tcult_tmp = tcult
        j += 1
    
    return rnet, rglo, albedolai, albsol, tcult, tcultmax, converge