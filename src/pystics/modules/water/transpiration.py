import numpy as np
from datetime import datetime

def potential_transpiration_coef(lai, etp, kmax, codeintercept, eos, esol, extin, mouill):
    """
    This module computes potential (eo) and maximal (eop) plant transpirations with crop coefficient approach (codebeso = 1).
    """

    # Potential plant transpiration
    eo = etp * (
        1
        + (kmax - 1)
        / (1 + np.exp(-1.5 * (lai - 3)))
    ) 

    if codeintercept == 1:
        emd = foliage_water_evap_lai(codeintercept, mouill, eo, etp, lai, extin)
    else:
        emd = 0

    # Maximal plant transpiration
    edirectm = min(eos + emd, eo)
    edirect = esol + emd
    
    if edirectm > 0:
        eop = (eo - edirectm) * (1.39999998 + (1-1.39999998) * (edirect / edirectm))
    else:
        eop = eo
    eop = max(eop, 0)
        
    return eo, eop, emd, edirect, edirectm

def foliage_water_evap_lai(codeintercept, mouill, eo, etp, lai_prev, extin):
    """
    This module computes evaporation (emd) from water intercepted by leaves.
    """

    if codeintercept == 2:
        emd = 0
    elif codeintercept == 1:
        emd = min(
            mouill,
            eo
            - etp
            * np.exp(-(extin - 0.2) * lai_prev),
        )

    return emd

def actual_transpiration(swfac_prev, eop):
    '''
    This module computes actual plant transpiration (ep) from potential transpiration and water stress index.
    '''

    ep = swfac_prev * eop

    return ep

def soil_contribution_to_transpiration(ep, epz_i, depth, hur_i, hmin, lracz_i, lev_i_prev, profsem, zrac):
    '''
    This module computes soil layers contributions (epz) to actual transpiration.
    See Section 11.5.3 of STICS book.
    '''

    if (ep != 0) & (lev_i_prev == 1):
        
        # Transpiration distribution to soil profile
        water_to_distribute = ep.copy()
        h = np.maximum(0, hur_i - hmin)

        while water_to_distribute > 1e-15:
            water_distributed = 0
            
            # Total soil contribution 
            total_soil_contribution = (lracz_i * (h - epz_i)).sum()

            if np.abs(total_soil_contribution) < 1.e-8:
                break
            
            for z_index in range(max(0,int(profsem)-1), min(round(zrac)+1, depth)): # same range as root density calculation

                # Soil layer (1 cm) contribution / total soil contribution
                soil_layer_contribution = water_to_distribute * lracz_i[z_index] * (h[z_index] - epz_i[z_index]) / total_soil_contribution

                # Check if soil layer contribution does not exceed available water
                if (soil_layer_contribution + epz_i[z_index]) > h[z_index]:
                    soil_layer_contribution = h[z_index] - epz_i[z_index]
                    epz_i[z_index] = h[z_index].copy()
                else:
                    epz_i[z_index] = epz_i[z_index] + soil_layer_contribution
                
                water_distributed = water_distributed + soil_layer_contribution
            
            water_to_distribute = water_to_distribute - water_distributed
            
    return epz_i