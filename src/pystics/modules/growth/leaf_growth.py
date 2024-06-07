import numpy as np


def leaf_growth(lev_i, lax_i, sum_upvt, stlevamf, vlaimax, stamflax, udlaimax, dlaimax, pentlaimax,
                    tcult_prev, tcxstop, tcmax, tcmin, adens, bdens, densite, turfac_prev, phoi, phoi_prev, tigefeuil,
                    phobase, rfpi, dlaimin, lai_prev, laicomp, stopfeuille, vmax_prev,
                    codlainet, dltaisenat_prev, fstressgel_prev, laisen_prev, lax, slamin, slamax, dltamsen_prev, dltaisen, sen, lan):
    '''
    This module computes the leaf growth (deltai) and its different components (ulai, deltai_dev, deltai_t, deltai_stress).
    See section 4.1 of STICS book.
    '''

    deltai, deltai_dev, deltaidens, deltai_t, ulai, deltai_stress, efdensite = 0, 0, 0, 0, 0, 0, 1
    vmax = vmax_prev.copy()
    dltaisenat = dltaisenat_prev.copy()

    # Leaf growth stop criteria
    compute_leaf_growth = True
    if ((stopfeuille == 'LAX') & (lax_i == 1)): # determinate growth crops
        compute_leaf_growth = False

    if compute_leaf_growth & (lev_i == 1):
        
        # Development effect on leaf growth
        ulai = np.where(
            sum_upvt
            < stlevamf,
            1 + (vlaimax - 1) * (sum_upvt / stlevamf),
            np.where(
                sum_upvt
                < stlevamf
                + stamflax,
                vlaimax
                + (3 - vlaimax)
                * (sum_upvt - stlevamf)
                / stamflax,
                0,
            ),
        )
        
        # Interplant competition effect on leaf growth
        efdensite = 1
        if ulai > 1:
            if lai_prev < laicomp:
                efdensite = 1
            else:
                if densite == 0:
                    efdensite = 0
                else:
                    efdensite = min(1, np.exp(np.log(densite / bdens) * (adens)))
        else:
            if densite == 0:
                efdensite = 0
            else:
                efdensite = min(1, np.exp(np.log(densite / bdens) * (adens)))
        
        if ulai == 1:
            efdensite = 1
        
        # Development and interplant competition components of leaf growth
        if (udlaimax == 3) | (ulai <= udlaimax):
            deltai_dev = dlaimax * (dlaimin + (1 - dlaimin)  / (1 + np.exp(pentlaimax * (vlaimax - ulai))))
            deltaidens = densite * efdensite
            vmax = deltai_dev * deltaidens

        else:    
            deltai_dev = vmax * (1 - (ulai - udlaimax) / (3 - udlaimax))**2 
            deltaidens = 1 
        
        # Thermal component of leaf growth
        deltai_t = max(0, tcult_prev - tcmin)
        if tcxstop >= 100:
            if tcult_prev > tcmax:
                deltai_t = tcmax - tcmin
        else:   
            if tcult_prev > tcmax:
                deltai_t = max(0,(tcmax - tcmin) * (tcult_prev - tcxstop) / (tcmax - tcxstop))

        # Water stress affecting leaf growth
        deltai_stress = turfac_prev

        # Leaf growth
        deltai = (
            deltai_dev
            * deltaidens
            * deltai_t
            * deltai_stress
        )

        # Negative effect of photoperiod during days shortening period
        if phoi < phoi_prev:
            deltai = (deltai * rfpi)
        if phoi < phobase:
            rfpi = 0
            deltai = (deltai * rfpi)

    elif lev_i == 1:
        ulai = 3
    
    # Cumulated outputs
    if codlainet == 1:
        lai = lai_prev + deltai - dltaisenat

        if (sen == 1) & (lan == 0):
            dltaisenat = lai_prev - lai
            if lai <= 0:
                lan = 0
                lai = 0
    else:
        lai = lai_prev + deltai - dltaisen

    mafeuilverte = lai / slamax * 100
    
    
    if (codlainet == 1) & (lai > 0):
        if (tcult_prev < tcmin) | (fstressgel_prev < 1):
            dltaisen = dltamsen_prev /100 * slamin # tustress not computed
            lai = lai - dltaisen
            dltaisen = dltaisen + dltaisenat
            lai = max(lai, 0)

        else:
            dltaisen = dltaisenat
    
    # Cumulated outputs
    laisen = laisen_prev + dltaisen # dltaisen i-1 if codlainet = 2, i if codlainet = 1

    if (lai <= 0) & (lan == 0) & (lax == 1):
        sen = 1
        lan = 1


    return deltai, deltai_dev, deltaidens, deltai_t, ulai, deltai_stress, efdensite, vmax, lai, mafeuilverte, dltaisen, dltaisenat, laisen, lan, sen