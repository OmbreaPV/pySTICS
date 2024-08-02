



def cumultative_partitioning(lev, mafeuil, matigestruc, mafruit_prev, maenfruit_prev, codeperenne, masec, restemp0, resplmax, densite, remobres, deltai, slamin, ratiotf, codeindetermin, cumdltaremobil_prev, dltams):
    '''
    This module computes the reserves remobilisation with the cumulative partitioning approach (see chapter 7 of STICS book with code_acti_reserves = 2).
    '''

    #########################
    ### Vegetative organs ###
    #########################

    # Reserve = total biomass - leaves / stem / harvested organs = biomass that can be accumulated or mobilized
    if (lev == 0) and (codeperenne == 2): 
        restemp = restemp0
    elif lev == 1:
        restemp = max(0, masec - mafeuil - matigestruc) #- mafruit_prev - maenfruit_prev
    else:
        restemp = 0
        
    # Reserve limit --> sink/source applied if this limit is reached
    restempmax = 10 * resplmax * densite
    if restemp > restempmax:
        dltams = 0
    
    dltarestemp = max(0, restemp-restempmax)

    ##############################
    ### Reserves remoblisation ###
    ##############################

    # 1. Source / sink ratio

    # Sink strength of vegegative organs
    fpv = (
        deltai
        * 10000
        / (slamin / (1 + ratiotf))
    )

    # Sink strength of reproductive organs (only for indeterminate growth plants)
    if codeindetermin == 1:
        fpft = 0
    else:
        fpft = 0 # TODO

    # Source / sink ratio (first computation)
    sourcepuits1 = min(1,
        dltams * 100
        / (fpv + fpft)
        if fpv + fpft != 0
        else 1
    )

    # 2. Reserves mobilisation if source/sink ratio < 1

    # Remobilised reserved
    if sourcepuits1 < 1:
        remob = (fpv + fpft) / 100 - dltams
    else:
        remob = 0

    dltaremobil = min(remob, remobres * restemp) # remobilisation limited by remobres parameter
    # cumdltaremobil = cumdltaremobil_prev + dltaremobil

    # Sink / source ratio after remobilisation
    sourcepuits = (dltams + dltaremobil) / (fpv + fpft) if fpv + fpft != 0 else 1 # same as sourcepuits1 here ?

    return dltaremobil, restemp, dltams, fpv, sourcepuits, dltarestemp