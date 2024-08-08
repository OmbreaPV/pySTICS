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

def net_radiation(albedo, hur_i0, hminf_1, hccf_1, albveg, lai, trg, tcult, temp, tpm, fracinsol, codernet,
                  raint, parsurrg):
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

    ### Available energy and soil/plant distribution ###
    fapar = raint / (
    parsurrg * trg
    )

    # Available energy for the plant
    rnetp1 = (
        0.83 * fapar * rnet
    )  # eq 9.26

    # Available energy for the soil
    rnets = rnet - rnetp1

    return rnet, rglo, albedolai, albsol, ratm, rnets


def crop_temperature(lev_i_prev, temp_max, rnet, et, temp_min, temp, z0, codecaltemp,
    ratm, tcultmin, wind, lai, trg, tcultmax, albedolai, daylen, rnets,
    zr, lai_prev, z0solnu, hauteur):
    ''''
    This module computes crop surface temperature with the empirical approach or energy balance approach.
    See Section 9.3.2 of STICS book.
    '''
    if lev_i_prev > 0:
    
        if codecaltemp == 1: # Empirical method

            # Max crop surface temperature
            tcultmax = temp_max + (
                rnet / 2.46 - et - 1.27
            ) / (1.68 / np.log(1 / z0))

            tcultmax = max(tcultmax, temp_max)
            tcultmin = temp_min

        elif codecaltemp == 2: # Energy balance method (P. Cellier)
            # hyp : atm radiation is constant over day
    
            # min & max aerodynamic resistances
            wind_min = max(0.02,wind)*0.5
            raamin, _ = calraero(zr, wind_min, lai_prev, z0solnu, hauteur)
            wind_max = max(0.02,wind)*1.5
            raamax, _ = calraero(zr, wind_max, lai_prev, z0solnu, hauteur)
            
            ratmh = ratm / 24 * 1e6 / 3600

            sigma = 5.67e-8
            rsolglo = sigma * (tcultmin + 273.15)**4
            rglo = ratmh - rsolglo

            # wind reduction
            v = max(wind * np.exp(-0.96 * lai), 0.2)
            rnetmin = rglo
            gmin = gsol(6.0,v*0.5,rnetmin,0.0)
            tcultmin = (rnetmin - gmin) * raamin / 1200 + temp_min

            etmax = et * 3.14 / 2 / daylen * 2.46 * 1e6 / 3600
            rgmax = trg * 3.14 / 2 / daylen * 1e6 / 3600

            rsolglo = sigma * (tcultmax + 273.15)**4
            rglo = ratmh - rsolglo
            rnetmax = (1 - albedolai) * rgmax + rglo

            rnetmax = max(rnetmax, rnetmin)

            rnetsmax = rnets / rnet * rnetmax
            gmax = 0.25 * rnetsmax
            tcultmax = (rnetmax - gmax - etmax) * raamax / 1200 + temp_max

        # Crop temperature 
        tcult = (
            tcultmin + tcultmax
        ) / 2

    else: # tcult before emergence does not make sense, but is necessary for soil temperature calculation
        tcult = temp

    return tcult, tcultmax, tcultmin

def gsol(heure,v,rnet,hsens):
    '''
    This module computes the minimum daily soil heat flux.
    '''
    # soil radiation
    omega = 2*3.14/24.0
    if (heure <= 7.0 or heure >= 18):
        Pgh = 1
    else:
        Pgh = np.cos(omega*(heure-12.0+1.5)) / np.cos(omega*(heure-12.0-1.0))

    if (rnet > 0.0):
        gsol = 1.36/np.sqrt(v)*hsens*Pgh
    elif (v < 1.0):
        gsol = 0.9*rnet
    else:
        gsol = max(0.3,(0.9+0.1*(v-1)))*rnet
        
    return gsol

def calraero(zr, wind, lai_prev, z0solnu, hauteur):
    '''
    This module computes the resistances to diffusion.
    '''

    karm = 0.41
    nconv = 2.5
    hauteur = max(hauteur, z0solnu / 0.10)
    d  = 0.66 * hauteur
    z0 = 0.10 * hauteur

    raalim = np.log((zr - d) / z0) / (karm**2 * wind) \
                * (np.log((zr - d) / (hauteur - d)) + hauteur / (nconv *(hauteur - d)) * 
                (np.exp(nconv * (1.0 - (d + z0) / hauteur)) - 1.0))
     
    raslim = np.log((zr - d) / z0) / (karm**2 * wind) * hauteur \
                / (nconv * (hauteur - d)) * (np.exp(nconv) - np.exp(nconv * (1.0 - (d + z0) / hauteur)))
    
    raszero = np.log(zr / z0solnu) * np.log((d + z0) / z0solnu) / (karm**2 * wind)
    raazero = np.log(zr / z0solnu)**2 / (karm**2 * wind) - raszero

    if (lai_prev < 4.0):
        raa = (0.25 * lai_prev * raalim) + (0.25 * (4 - lai_prev) * raazero)
        ras = (0.25 * lai_prev * raslim) + (0.25 * (4 - lai_prev) * raszero)
    else:
        raa = raalim
        ras = raslim

    return raa, ras


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


def iterative_calculation(temp, lev_i_prev, temp_max, temp_min, et, z0, albedo, hur_i0, hminf_1, hccf_1, albveg, lai, trg, tpm, fracinsol, codernet, codecaltemp,
                          raint, parsurrg, ratm, tcultmin, wind, tcultmax, daylen, zr, lai_prev, z0solnu, hauteur):
    '''
    This module performs the iterative calculation of crop temperature and net radiation.
    See section 9.3.2.3 of STICS book.
    '''

    tcult_tmp = temp
    tcultmin = temp_min
    tcultmax = temp_max

    j = 0
    while j < 5:

        # Net radiation
        rnet, rglo, albedolai, albsol, ratm, rnets = net_radiation(albedo, hur_i0, hminf_1, hccf_1, albveg, lai, trg, tcult_tmp, temp, tpm, fracinsol, codernet,
                                                            raint, parsurrg)

        # Crop temperature
        tcult, tcultmax, tcultmin = crop_temperature(lev_i_prev, temp_max, rnet, et, temp_min, temp, z0, codecaltemp,
                                ratm, tcultmin, wind, lai, trg, tcultmax, albedolai, daylen, rnets, zr, lai_prev, z0solnu, hauteur)

        if np.isnan(tcult):
            raise pysticsException(
            detail="tcult is nan",
        )

        if tcult - tcult_tmp < 0.5:
            converge = True
            break
        tcult_tmp = tcult
        j += 1
    
    return rnet, rglo, albedolai, albsol, tcult, tcultmax, converge, tcultmin