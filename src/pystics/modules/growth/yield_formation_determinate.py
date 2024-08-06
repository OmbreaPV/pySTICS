import numpy as np

def harvested_organs_number(dltams_list, ind_drp, nbjgrain, cgrainv0, cgrain, nbgrmax, nbgrmin):
    '''
    This module computes the number of harvested organs.
    See section 8.1.1 of STICS book.
    '''
    
    if ind_drp == 0:
        vitmoy = 0
        nbgrains = 0
    else:
        # Au stade 69 ('fin floraison, début formation fruit visible') --> Calcul du nb de grains.
        # 71 = 10% des fruits ont leur taille max ou les fruits ont 10% de leur taille max

        # Average biomass production on the grains/fruits formation period (before filling start)
        vitmoy = (dltams_list[ind_drp - nbjgrain + 1 : ind_drp+1].sum() / nbjgrain)*100

        # Grains/fruits number
        nbgrains = (cgrainv0 + cgrain * vitmoy) * nbgrmax
        nbgrains = min(nbgrains, nbgrmax) 
        nbgrains = max(nbgrains, nbgrmin)

    return vitmoy, nbgrains

def frost_reduction_harvested_organs(i, ind_drp, nbgrains_prev, fgelflo, pgrain_list, nbgraingel_list):
    '''
    This module computes the mass of harvestable organs destroyed by frost (pgraingel) during grains/fruits filling.
    See section 8.1.1 of STICS book.
    '''

    nbgrains = nbgrains_prev.copy()
    pgraingel = 0
    if (i > ind_drp) & (ind_drp != 0):
        nbgraingel_list[i]  = nbgrains_prev * (1 - fgelflo)
        nbgrains =  (nbgrains_prev - nbgraingel_list[i])

        pgraingel = (pgrain_list.loc[ind_drp:i-1] * nbgraingel_list[ind_drp+1:i+1]).sum()

    return nbgrains, nbgraingel_list, pgraingel


def carbon_harvest_index(ircarb_list, ind_mat, ind_drp, codeir, vitircarb, vitircarbt, irmax, sum_upvt_post_lev_list):
    """
    This module computes the harvest index (ircarb) during fruit filling.
    See section 8.1.1 of STICS book.
    """
    if ind_drp != 0:
        index = np.array([i for i in range(len(ircarb_list))])
        if codeir == 1:
            ircarb_list[(index > ind_drp)] = np.minimum(vitircarb * (index[(index > ind_drp)] - ind_drp +1), irmax)
        else:
            ircarb_list[(index > ind_drp)] = vitircarbt * (sum_upvt_post_lev_list[(index > ind_drp)] - sum_upvt_post_lev_list[ind_drp])
        if ind_mat != 0:
            ircarb_list[index > ind_mat] = ircarb_list[ind_mat]

    return ircarb_list



def harvested_organs_mass(i, ircarb, masec, ircarb_prev, masec_prev, ftempremp, pgrainmaxi, nbgrains, ind_mat, pgraingel, deltags_list, ind_drp, mafruit_prev, pgrain_prev, nbgraingel):
    '''
    This module computes the mass of harvested organs (deltags) between fruit filling start (drp) and physiological maturity (mat), reduced by the mass of frozen grains.
    See section 8.1.1 of STICS book.
    '''

    # Daily grain/fruit mass increase
    deltags_list[i] = (ircarb * masec - ircarb_prev * masec_prev) * ftempremp
    if (ind_mat != 0) & (i >= ind_mat):
        deltags_list[i] = 0

    # Cumulated grains/fruits mass, reduced with frozen grains
    mafruit = mafruit_prev + deltags_list[i] - pgrain_prev * nbgraingel

    if mafruit > (pgrainmaxi * nbgrains / 100):
        mafruit = pgrainmaxi * nbgrains / 100
        deltags_list[i] = 0

    # Average mass per grain/fruit
    if nbgrains > 0:
        pgrain = np.minimum(mafruit / nbgrains * 100, pgrainmaxi)
    else:
        pgrain = 0

    return mafruit, deltags_list, pgrain


def grains_water_content(tcult_list, temp_list, ind_debdes, h2ofrvert, deshydbase, tempdeshyd):
    '''
    This module computes grains/fruits water content (teaugrain).
    See section 8.2.1 of STICS book.
    '''

    teaugrain_list = np.array([h2ofrvert for i in range(len(tcult_list))])
    index = np.array([i for i in range(len(tcult_list))])
    teaugrain_list[index >= ind_debdes] = h2ofrvert - deshydbase * (index[index >= ind_debdes] - ind_debdes + 1) - tempdeshyd * (tcult_list[index >= ind_debdes] - temp_list[index >= ind_debdes]).cumsum()
    teaugrain_list = np.maximum(teaugrain_list, 0)

    return teaugrain_list