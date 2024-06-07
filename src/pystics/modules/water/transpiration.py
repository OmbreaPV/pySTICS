import numpy as np


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


def potential_transpiration_energy_budget(lai_prev, rnetp1, emd, deltat, gamma, rc, rac, dos, L):
    '''
    This module computes maximal transpiration (eop) with resistance approach (codebeso = 2).
    '''
    if lai_prev == 0:
        eop = 0  # C'est pas qu'il n'y a pas d'énergie dispo pr transpi, mais que rac n'est pas calculabele
    else:
        # Calcul de l'énergie dispo pour la transpi de la plante (= énergie dispo pour la plante - ce qui est utilisé pr évap eau sur plante)
        rnetp2 = (rnetp1 - emd)

        # Calcul de la transpi plante maximale (potentiel ce flux --> réel va dépendre de eau dispo dans sol)
        # on l'appelle eop alors qu'on parle encore du potentiel. Dans Benj / miro j'avais mis eo = potentielle et eop = réel. Stics : eop = transpi max
        eop = (
            deltat * rnetp2
            + 105.03 * dos / rac
        ) / (
            L
            * (
                deltat
                + gamma * (1 + rc / rac)
            )
        )
        eop = max(eop, 0) # Max evaporation can only be positive

    return eop, rnetp2

def actual_transpiration(swfac_prev, eop):
    '''
    This module computes actual plant transpiration (ep) from potential transpiration and water stress index.
    '''

    ep = swfac_prev * eop

    return ep

def soil_contribution_to_transpiration(ep, epz_i, depth, hur_i, hmin, lracz_i, lev_i_prev):
    '''
    This module computes soil layers contributions (epz) to actual transpiration.
    See Section 11.5.3 of STICS book.
    '''


    if (ep == 0) | (lev_i_prev == 0):
        epz_i[:] = [0 for z in range(depth)]
    
    elif (lev_i_prev == 1):
        
        # Transpiration distribution to soil profile
        water_to_distribute = ep.copy()
        epz_i[0:depth] = 0
        h = np.maximum(0, hur_i - hmin)

        while water_to_distribute > 1e-15:
            water_distributed = 0
            total_soil_contribution = 0
            
            # Total soil contribution 
            for z_index in range(depth):
                total_soil_contribution = total_soil_contribution + np.where(np.isnan(lracz_i[z_index]),0,lracz_i[z_index]) * (h[z_index]- np.where(np.isnan(epz_i[z_index]),0,epz_i[z_index]))

            if np.abs(total_soil_contribution) < 1.e-8:
                break

            for z_index in range(depth):

                # Soil layer (1 cm) contribution / total soil contribution
                soil_layer_contribution = water_to_distribute * np.where(np.isnan(lracz_i[z_index]),0,lracz_i[z_index]) * (h[z_index] - np.where(np.isnan(epz_i[z_index]),0,epz_i[z_index])) / total_soil_contribution # le dernier terme sert à gérer boucle while = là on veut considérer que la partie de cette itération

                # Check if soil layer contribution does not exceed available water
                if (soil_layer_contribution + epz_i[z_index]) > h[z_index]:
                    soil_layer_contribution = h[z_index] - epz_i[z_index]
                    epz_i[z_index] = h[z_index].copy()
                else:
                    epz_i[z_index] = np.where(np.isnan(epz_i[z_index]), 0, epz_i[z_index]) + soil_layer_contribution
                
                water_distributed = water_distributed + soil_layer_contribution

            
            water_to_distribute = water_to_distribute - water_distributed

    return epz_i