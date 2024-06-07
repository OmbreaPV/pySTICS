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
        vitmoy = (dltams_list[ind_drp - nbjgrain + 1 : ind_drp].sum() / nbjgrain)

        # Grains/fruits number
        nbgrains = (cgrainv0 + cgrain * vitmoy) * nbgrmax
        nbgrains = min(nbgrains, nbgrmax) 
        nbgrains = max(nbgrains, nbgrmin)

    return vitmoy, nbgrains


def frost_reduction_harvested_organs(i, ind_drp, nbgrains_prev, fgelflo):
    '''
    This module computes the number of harvestable organs destroyed by frost (nbgraingel) during grains/fruits filling.
    See section 8.1.1 of STICS book.
    '''

    if (i > ind_drp) & (ind_drp != 0): 
        nbgraingel = nbgrains_prev * (1 - fgelflo)
        nbgrains = (nbgrains_prev - nbgraingel)

    return nbgrains, nbgraingel


def carbon_harvest_index(ircarb_list, somcourdrp_list, nbgrains_list, ind_mat, ind_drp, codeir, vitircarb, vitircarbt, irmax, sum_upvt_post_lev_list):
    """
    This module computes the harvest index (ircarb) during fruit filling.
    See section 8.1.1 of STICS book.
    """
    if ind_drp != 0:
        index = np.array([i for i in range(len(nbgrains_list))])
        if codeir == 1:
            ircarb_list[(index > ind_drp)] = np.minimum(vitircarb * (index[(index > ind_drp)] - ind_drp +1), irmax)
        else:
            ircarb_list[(index > ind_drp)] = vitircarbt * (sum_upvt_post_lev_list[(index > ind_drp)] - sum_upvt_post_lev_list[ind_drp])
        if ind_mat != 0:
            ircarb_list[index > ind_mat] = ircarb_list[ind_mat]

    return ircarb_list, somcourdrp_list


def harvested_organs_mass(ircarb_list, masec_list, ircarb_prev_list, masec_prev_list, ftempremp_list, pgraingel_list, pgrainmaxi, nbgrains_list, ind_mat):
    '''
    This module computes the mass of harvested organs (deltags) between fruit filling start (drp) and physiological maturity (mat).
    See section 8.1.1 of STICS book.
    '''

    # Daily grain/fruit mass increase
    deltags_list = (ircarb_list * masec_list - ircarb_prev_list * masec_prev_list) * ftempremp_list
    index = np.array([i for i in range(len(ircarb_list))])
    if ind_mat != 0:
        deltags_list[index >= ind_mat] = 0

    # Cumulated grains/fruits mass
    mafruit_list = np.minimum((deltags_list - pgraingel_list / 100).cumsum(), pgrainmaxi * nbgrains_list)

    # Average mass per grain/fruit
    pgrain_list = np.minimum(mafruit_list / nbgrains_list * 100, pgrainmaxi)

    return mafruit_list, deltags_list, pgrain_list


def grains_water_content(tcult_list, temp_list, ind_debdes, h2ofrvert, deshydbase, tempdeshyd):
    '''
    This module computes grains/fruits water content (teaugrain).
    See section 8.2.1 of STICS book.
    '''

    teaugrain_list = np.array([h2ofrvert for i in range(len(tcult_list))])
    index = np.array([i for i in range(len(tcult_list))])
    teaugrain_list[index >= ind_debdes] = h2ofrvert - deshydbase * (index[index >= ind_debdes] - ind_debdes + 1) - tempdeshyd * (tcult_list[index >= ind_debdes] - temp_list[index >= ind_debdes]).cumsum()


    return teaugrain_list