



def lethal_frost(lai, mafeuilverte, laisen_prev, dltaisen, stopfeuille_stage, restemp, codeperenne, amf, lax, sen, lan, drp, mat, sum_upvt_post_lev, stsenlan, stlaxsen, stlevamf, stamflax, stlevdrp, stdrpmat):
    '''
    This modules computes the new pools when a lethal frost has occured.
    '''

    # Leaves and biomass pools
    dltaisen = lai
    dltamsen = mafeuilverte
    laisen = laisen_prev + dltaisen

    # Change phenological stages reached
    if (stopfeuille_stage == 1) & (lan == 0):
        lan = 1
        stsenlan = sum_upvt_post_lev
        if sen == 0:
            sen = 1
            stlaxsen = sum_upvt_post_lev
            stsenlan = 0
    elif (restemp == 0) & (codeperenne == 1):
        if amf == 0:
            amf = 1
            amf = 1
            lax = 1
            sen = 1
            lan = 1
            stlevamf = sum_upvt_post_lev
            stamflax = 0
            stlaxsen = 0
            stsenlan = 0
        if (lax == 0) & (amf == 1):
            lax = 1
            sen = 1
            lan = 1
            stamflax = sum_upvt_post_lev
            stlaxsen = 0
            stsenlan = 0
        if (sen == 0) & (lax == 1):
            sen = 1
            lan = 1
            stlaxsen = sum_upvt_post_lev
            stsenlan = 0    
    
    stlevdrp = sum_upvt_post_lev
    stdrpmat = 0
    drp = 1
    mat = 1
    
    return dltaisen, dltamsen, laisen, amf, lax, sen, lan, drp, mat, stsenlan, stlaxsen, stlevamf, stamflax, stlevdrp, stdrpmat