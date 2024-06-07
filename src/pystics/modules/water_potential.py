import numpy as np

def plant_soil_water_potential(lev_i_prev, zrac0, depth, psihucc, psihumin, daf, hmin, hcc, potgermi, profsem, zrac, lracz_i, hur_i, racinepsi_i):
    '''
    This module computes the soil potential (psisol) and the predawn leaf water potential (psibase).
    '''

    # Soil potential
    wsat = np.array([1 for i in range(depth)]) - daf / 2.66
    bpsisol = np.log(psihucc / psihumin) / np.log(hmin / hcc)
    psisols = psihumin * np.power(hmin / (wsat*10 ), bpsisol)
    psisol_i = psisols * np.power(hur_i / wsat, -bpsisol)
    humpotsol_i = wsat * 10 * np.power(np.array([potgermi for i in range(depth)]) / psisols, -1/bpsisol) # dans Divers_develop.f90 mais pas dans bible.

    # Predawn leaf water potential
    if ((lev_i_prev == 1) | (zrac0 != 0)) & (zrac != 0):
        # print('zrac',zrac)
        # print('psisol_i > -1.5',psisol_i > -1.5)
        # print('lracz_i[psisol_i > -1.5]',lracz_i[psisol_i > -1.5])
        
        racinepsi_i[psisol_i > -1.5] = lracz_i[psisol_i > -1.5]
        # print(' racinepsi_i', racinepsi_i)
        # print('racinepsi_i[max(0,int(profsem)-1):round(zrac)+1].sum()',racinepsi_i[max(0,int(profsem)-1):round(zrac)+1].sum())
        psibase = (psisol_i[max(0,int(profsem)-1):int(zrac)+1] * psisol_i[max(0,int(profsem)-1):int(zrac)+1]).sum() / racinepsi_i[max(0,int(profsem)-1):round(zrac)+1].sum()
    else:
        psibase = 0


    return humpotsol_i, psisol_i, racinepsi_i, psibase