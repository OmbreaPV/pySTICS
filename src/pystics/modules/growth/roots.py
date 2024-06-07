import numpy as np


def root_growth(ger_i, lax_i, codeperenne, findorm_i, profsem, zrac, codetemprac,
                tcult_prev, tcmax, tcmin, tgmin, croirac, hur_i, tsol_i, hmin, sensrsec, depth, daseuilbas, daseuilhaut, contrdamax, daf, lev_i, hcc, herbaceous, znonli):
    """
    This module computes daily root growth (deltaz), maximum root depth (zrac) and maximum root depth without physical or phenological stop (znonli).

    Simplifications compared to STICS :
        - Growth stop criteria not implemented : when perennial plant and STOPRAC == 'sen' (stoprac = LAX in pySTICS) or when OBSTARAC reached.
        - Water excess stress index affecting root growth (izrac) not implemented.
    """

    deltaz, deltaz_t, deltaz_stress, humirac_mean = 0, 0 , 0, 1

    # Bulk density effect on root growth
    efda = compute_efda(daf[np.clip(int(zrac), 0, depth-1)], daseuilbas, daseuilhaut, contrdamax)

    # Phenological root growth stop
    compute_root_growth = True
    if lax_i == 1: # # TODO : add criteria on STOPRAC ('LAX' or 'SEN'). pySTICS : we consider that STOPRAC = LAX for all species
        compute_root_growth = False


    if ((
        ger_i == 1
        & ((codeperenne == 1) | herbaceous)
    ) or (
        findorm_i == 1
        and codeperenne == 2
    )):
        
        # Root apex depths
        ap = range(max(0,min(int(zrac)-1,depth-1)),min(int(zrac)+2,depth))
        
        # Thermal component of root growth
        if codetemprac == 1:
            deltaz_t = max(0,(min(tcult_prev,tcmax) - tcmin))
        elif codetemprac == 2:
            deltaz_t = max(0,(min(tsol_i[ap].mean(),tcmax) - tgmin))
        

        # Water stress affecting root growth
        hur_mean = hur_i[ap].mean()
        humin_mean = hmin[ap].mean()
        if lev_i == 0:
            hcc_mean = hcc[ap].mean()
            if hur_mean > humin_mean:
                x = (hur_mean - humin_mean) / (hcc_mean - humin_mean)
                humirac_ap = np.clip(sensrsec + (1 - sensrsec) * x, 0, 1)
            elif hur_mean <= humin_mean:
                humirac_ap = np.clip(hur_mean * sensrsec / humin_mean, 0, 1)
        else:
            if hur_mean >= humin_mean:
                humirac_ap = 1
            elif hur_mean < humin_mean:
                humirac_ap = np.clip(hur_mean * sensrsec / humin_mean, 0, 1)

        # Water stress on whole depth
        humirac_list = []
        if lev_i == 0:
            for z_index in range(max(0,int(profsem)-1), max(int(profsem)+1,int(min(round(zrac), depth)))):
                if hur_i[z_index] > hmin[z_index]:
                    x = (hur_i[z_index] - hmin[z_index]) / (hcc[z_index] - hmin[z_index])
                    humirac_list.append(np.clip(sensrsec + (1 - sensrsec) * x, 0, 1))
                elif hur_i[z_index] <= hmin[z_index]:
                    humirac_list.append(np.clip(hur_i[z_index] * sensrsec / hmin[z_index], 0, 1))
        else:
            for z_index in range(max(0,int(profsem)-1), int(min(round(zrac), depth))):
                if hur_i[z_index] >= hmin[z_index]:
                    humirac_list.append(1)
                elif hur_i[z_index] < hmin[z_index]:
                    humirac_list.append(np.clip(hur_i[z_index] * sensrsec / hmin[z_index], 0, 1))
        humirac_mean = np.mean(humirac_list)

        # Stress component of root growth
        deltaz_stress = humirac_ap * efda

        # Root growth
        deltaz = croirac * deltaz_t * deltaz_stress

        # Maximum root depths
        zrac = min(deltaz * compute_root_growth + zrac, depth)  
        znonli = znonli + deltaz
        deltaz = deltaz * compute_root_growth

    return zrac, deltaz, deltaz_t, deltaz_stress, efda, humirac_mean, znonli

def root_density(lracz_i, zrac, znonli, depth, zprlim, zpente, ger_i, codeperenne, lvopt, s, profsem,
                 hur_i, hmin, humirac_i):
    '''
    This module computes root density with standard profile method (coderacine = 2).
    
    Simplifications compared to STICS :
        - True density approach (coderacine=1) is not implemented.
    '''

    humirac_mean = 1

    # Necessary root depth to absorb 20% of water
    zdemi = max(znonli - zprlim + zpente, (np.log(4) / s))

    ZRac_Ceil = min(round(zrac)+1, depth)

    if ((ger_i == 1) or (codeperenne == 2)) & (zrac != 0):
        humirac_mean = 0
        for z_index in range(max(0,int(profsem)-1), round(ZRac_Ceil)): 
            
            # Water stress index affecting root density
            if hur_i[z_index] >= hmin[z_index]:  
                humirac_i[z_index] = 1
            else:
                humirac_i[z_index] = min(1, max(0, 0 * hur_i[z_index] / hmin[z_index])) 
            
            # Root density
            lracz_i[z_index] = lvopt * humirac_i[z_index] / (1 + np.exp(-s * ((z_index + z_index+1)/2 - zdemi)))
            humirac_mean = humirac_mean + humirac_i[z_index]
        
        humirac_mean = humirac_mean / (round(ZRac_Ceil) - max(0,int(profsem)-1))
                                       
    # Cumulated root length density
    cumlracz = np.nansum(lracz_i)

    return lracz_i, cumlracz, zdemi, humirac_i


def compute_efda(daf, daseuilbas, daseuilhaut, contrdamax):
    '''
    This function computes bulk density effect on root growth.
    '''
         
    if (daf <= daseuilbas):
        efda = 1
        return efda
    
    if (daf >= daseuilhaut):
        efda = contrdamax
        return efda
   
    dx = daseuilhaut - daseuilbas
    if abs(dx) < 1e-8:
        dx = 1.e-10
    pente = (contrdamax - 1) / dx

    efda = 1 + (pente * (daf - daseuilbas))

    return efda
      