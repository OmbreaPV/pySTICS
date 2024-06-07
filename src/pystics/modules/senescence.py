import numpy as np
from pystics.modules.growth.thermal_stress_indices import frost_stress



def senescence(i, lev, ulai, somtemp, vlaimax, durviei, durvief, senfac_prev, udevcult, fstressgel,
               dayLAIcreation_list, senstress_list, tdevelop_list, durvie_list, durage_list, deltai_list,
                somsenreste_prev, lai0, ndebsen, dltafv_list, ratiosen, forage, msres_prev, msresiduel, msresjaune_prev, codlainet):
    '''
    This module computes leaf and biomass senescence.
    See section 4.1.2 of STICS book.

    Simplifications compared to STICS :
        - Nitrogen effect on senescence not implemented
    '''

    # Initialize values
    dltaisen, dltamsen, somsenreste, deltamsresen, msresjaune = 0, 0, 0, 0, 0
    msres = msres_prev

    if lev[i] > 0:        
        ind_EMERG = np.where(lev > 0)[0][0]

        # Initialize oldest non-senescent leaf area
        if i == ind_EMERG:
            dayLAIcreation_list[i] = i

        # Stress index affecting senescence
        senstress_list[i] = min(senfac_prev, fstressgel)
        if codlainet == 1:
            senstress_list[i] = 1

        # Update lifespan of all non-senescent leaf areas
        for j in range(int(dayLAIcreation_list[i]),i):
            durvie_list[j] = min(durvie_list[j], durage_list[j] * senstress_list[i])
        
        # Lifespan of leaf area produced on day i
        if (ulai <= vlaimax):
            durage_list[i] = (durviei)
        else: 
            durage_list[i] = durviei + (durvief - durviei) * (ulai - vlaimax) / (3 - vlaimax)

        durvie_list[i] = durage_list[i] * senstress_list[i]

        # Development temperature (to compare to lifespan)
        tdevelop_list[i] = 2**(udevcult/10)

        # Senescence of residual biomass for forage crops
        if forage:
            if msres_prev > 0:
                deltamsresen = msresiduel * ratiosen * 2**(udevcult / 10) / durviei
                msres_tmp = max(0, msres_prev - deltamsresen)
                deltamsresen = msres_prev - msres_tmp
                msres = msres_prev - deltamsresen
            
            # Cumulated senescent residual biomass
            msresjaune = msresjaune_prev + deltamsresen
            
        
        deltai_disappears = True
        nb_deltai_senescent = 0
        

        if somtemp > durvie_list[i]:

            if ndebsen == 0: # first day of senescence
                ndebsen = i
                if codlainet == 2:
                    dltaisen = dltaisen + lai0
                dltamsen =  dltafv_list[i] * ratiosen

            somsen = somsenreste_prev + tdevelop_list[int(dayLAIcreation_list[i]):i+1].sum()
            while deltai_disappears: # as long as deltai become senescent, we check senescence of following day (several deltai can become senescent on day i)
                
                # Oldest non-senescent deltai
                j = int(dayLAIcreation_list[i]) + nb_deltai_senescent

                # Compare development temperature to leaf area lifespan, and trigger senescence
                if somsen >= durvie_list[j]:
                    if codlainet == 2:
                        dltaisen = dltaisen + deltai_list[j]
                    somsen = somsen - durvie_list[j]
                    dltamsen = dltamsen + ratiosen * dltafv_list[j]
                    nb_deltai_senescent +=1
                else:
                    deltai_disappears = False
            
            somsenreste = somsen # fortran je comprends pas le sens de ce reste
        
        # Update oldest non-senescent deltai
        dayLAIcreation_list[i+1] = dayLAIcreation_list[i] + nb_deltai_senescent

    return dayLAIcreation_list, durage_list, senstress_list, tdevelop_list, durvie_list, dltaisen, somsenreste, ndebsen, somtemp, dltamsen, deltamsresen, msres, msresjaune, durvie_list[i]

def senescence_stress(lev_i, ulai, vlaimax, temp_min_prev, tgeljuv10, tgeljuv90, tgelveg10, tgelveg90, tletale, tdebgel, codgeljuv, codgelveg):
    '''
    This module computes frost stress affecting senescence.
    '''

    if lev_i > 0:
        
        # Frost stress
        if (
            ulai < vlaimax
        ):  # si on est avant stade iamf --> fgeljuv
            if codgeljuv == 2:
                fstressgel = frost_stress(temp_min_prev, tgeljuv90, tgeljuv10, tletale, tdebgel)

        else:  # après stade iamf -> fgelveg
            if codgelveg == 2:
                fstressgel = frost_stress(temp_min_prev, tgelveg90, tgelveg10, tletale, tdebgel)
    else:
        fstressgel = 1

    return fstressgel

