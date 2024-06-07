import numpy as np

def biomass_production(masec_prev, raint, lev, amf_prev, drp_prev, efcroijuv, efcroirepro, efcroiveg,
                             coefb, swfac_prev, fco2, ftemp, deltai, slamax):
    '''
    This module computes the daily aerial biomass production according to phenological stage, radiation, thermal and water stress, and CO2 effect.
    '''

    # Maximal radiation use efficiency (RUE)
    ebmax = 0
    if drp_prev == 1:
        ebmax = round(efcroirepro / 100, 5)
    elif amf_prev == 1:
        ebmax = round(efcroiveg / 100, 5)
    elif lev == 1:
        ebmax = round(efcroijuv / 100, 5)

    # Daily biomass production
    dltams = shoot_biomass_production(
        raint,
        ebmax,
        coefb,
        ftemp,
        swfac_prev,
        fco2,
    )

    # Green leaf biomass production
    dltafv = min(deltai / slamax * 100, dltams)

    # Cumulated aboveground biomass
    masec = masec_prev + dltams

    return dltams, masec, ebmax, dltafv

def shoot_biomass_production(radiation, ebmax, coefb, ftemp, swfac_prev, fco2):
    """
    This function computes the aerial biomass production.
    """
    rue = (ebmax - coefb * radiation)
    dltams = radiation * rue * ftemp * swfac_prev * fco2
    return dltams


def fallen_biomass(abscission, dltamsen, mafeuiltombe_prev, mafeuiljaune_prev):
    '''
    This module computes the biomass of fallen leaves and removes it from the biomass of yellow leaves.
    '''
    dltamstombe = abscission * dltamsen
    mafeuiltombe = mafeuiltombe_prev + dltamstombe
    mafeuiljaune = mafeuiljaune_prev + dltamsen - dltamstombe
    
    return dltamstombe, mafeuiltombe, mafeuiljaune
