import numpy as np
from pystics.modules.growth.thermal_stress_indices import frost_stress
from pystics.modules.water.water_stress import water_stress_on_root_growth

def emergence_macro(i, outputs, crop, soil, manage, tsol, humirac, hur, humpotsol):
    outputs['moist'], outputs.loc[i,'nbjgrauto'], outputs.loc[i,'nbjhumec'], humirac, outputs.loc[i,'somger'], outputs['ger'], outputs.loc[i,'zrac'], outputs.loc[i,'elong'], outputs['lev'], outputs.loc[i,'coeflev'], outputs['densite'], outputs['let'], outputs.loc[i,'udev'], outputs.loc[i,'somfeuille'], outputs.loc[i,'nbfeuille'], outputs.loc[i,'fgellev'], outputs.loc[i,'somelong'] = emergence(i, outputs['densite'].array, outputs['lev'].array, outputs.loc[i-1,'lev'], outputs['ger'].array, outputs.loc[i-1,'ger'], outputs['moist'].array, outputs.loc[i-1,'moist'], outputs['let'].array, manage.PROFSEM, hur[i], humpotsol[i], crop.PROPJGERMIN, crop.NBJGERLIM, crop.TDMAX, crop.TDMIN, crop.TGMIN,
                                         tsol, outputs.loc[i-1,'nbjhumec'], soil.HMIN, crop.SENSRSEC, soil.HCC, outputs.loc[i-1,'somger'], crop.STPLTGER, tsol[i], manage.DENSITESEM, crop.CODEHYPO,
                                         outputs.loc[i,'zrac'], crop.BELONG, crop.CELONG, crop.ELMAX, outputs.loc[i-1,'tcult'], crop.NLEVLIM1, crop.NLEVLIM2, crop.TCXSTOP, outputs.loc[i-1,'somfeuille'], outputs.loc[i,'nbfeuille'], outputs.loc[i-1,'temp_min'],
                                         crop.TGELLEV90, crop.TGELLEV10, crop.TLETALE, crop.TDEBGEL, crop.PHYLLOTHERME, crop.NBFGELLEV, humirac, soil.DEPTH, crop.CODGELLEV, outputs.loc[i-1,'let'], outputs.loc[i-1,'somelong'], crop.CODETEMP, outputs.loc[i,'temp'])

    return outputs, humirac


def emergence(i, densite_list, lev, lev_i_prev, ger, ger_i_prev, moist, moist_i_prev, let, profsem, hur_i, humpotsol_i, propjgermin, nbjgerlim, tdmax, tdmin,tgmin,
              tsol, nbjhumec_prev, hmin, sensrsec, hcc, somger_prev, stpltger, tsol_i, densitesem, codehypo,
              zrac, belong, celong, elmax, tcult_prev, nlevlim1, nlevlim2, tcxstop, somfeuille_prev, nbfeuille, temp_min_prev,
            tgellev90,tgellev10, tletale, tdebgel, phyllotherme, nbfgellev, humirac, depth, codgellev, let_i_prev, somelong_prev, codetemp, temp):
    '''
    This module computes emergence of herbaceous plants : moistening, germination and plantlet growth.
    See Section 3.4.1.3 of STICS book.

    Simplifications compared to STICS :
        - Soil crusting not implemented
    '''

    nbjgrauto, nbjhumec, somger, lev_i, elong, coeflev, udev, somfeuille, fgellev, somelong = 0,0,0,0,0,0,0,0,1, 0

    if ger_i_prev == 0: 
        
        ##################
        ### MOISTENING ###
        ##################

        # Seed bed depth (profsem +/-1)
        sb = range(max(0,int(profsem) - 1),int(profsem) + 2)

        # Moistening
        if moist_i_prev != 1:
            if hur_i[sb].mean() > humpotsol_i[sb].mean():
                moist[i:len(moist)] = 1

        # Number of autotrophy days after moistening
        if moist[i] == 1:
            ind_moist = np.where(moist > 0)[0][0]
            nbjgrauto = max(propjgermin * nbjgerlim, min(nbjgerlim,(1 - propjgermin) / (tdmax - tdmin) * (tsol[ind_moist:i,sb].sum(axis=0).mean() / (i - ind_moist + 1) - tgmin) +1))     
            nbjhumec = nbjhumec_prev + 1

        ###################
        ### GERMINATION ###
        ###################

        # Water stress index affecting germination
        len_sb = len([i for i in sb])
        hur_sb = sum([hur_i[z_index] for z_index in sb]) / len_sb
        hmin_sb = sum([hmin[z_index] for z_index in sb]) / len_sb
        hcc_sb = sum([hcc[z_index] for z_index in sb]) / len_sb

        humirac[i,sb] = water_stress_on_root_growth(hur_sb, hmin_sb, hcc_sb, sensrsec, 2)

        # Growing degree days to reach germination
        somger = somger_prev + max(0, tsol_i[int(profsem)] - tgmin) * humirac[i,sb].mean() # does not depend on sowing day because this module is called when plt=1
        if somger >= stpltger:
            ger[i:len(ger)] = 1
            zrac = profsem
            somelong = somger - stpltger

        # Plant density reduction because of late germination after humectation
        if nbjhumec >  nbjgrauto:
            densite_list[i] = max(densitesem, densitesem * somger / stpltger)
        else:
            densite_list[i] = densitesem
        
    elif lev_i_prev == 0:

        ind_ger = np.where(ger > 0)[0][0]

        #######################
        ### PLANTLET GROWTH ###
        #######################
        if codehypo == 2:
            lev[i] = 1
        
        else:

            sb = range(max(0,int(profsem) - 1),int(profsem) + 2)
            hb = range(sb[0], min(max(sb[-1],int(zrac))+1, depth))

            len_hb = len([i for i in hb])
            hur_hb = sum([hur_i[z_index] for z_index in hb]) / len_hb
            hmin_hb = sum([hmin[z_index] for z_index in hb]) / len_hb
            hcc_hb = sum([hcc[z_index] for z_index in hb]) / len_hb

            somelong = somelong_prev + max(0, tsol_i[int(profsem)] - tgmin) * water_stress_on_root_growth(hur_hb, hmin_hb, hcc_hb, sensrsec, 2)

            # Plantlet elongation
            elong = elmax * (1 - np.exp(-(belong* somelong)**celong))

            # Emergence
            if elong > profsem:
                lev[i:len(lev)] = 1

    if (lev_i_prev == 0) and (lev[i] == 1): # emergence day

        # Germination and emergence dates
        ind_ger = np.where(ger > 0)[0][0]
        ind_lev = np.where(lev > 0)[0][0]

        # Density reduction component
        if (ind_lev - ind_ger < nlevlim1) or (codehypo == 2):
            coeflev = 1
        elif (ind_lev - ind_ger >= nlevlim1) and (ind_lev - ind_ger <= nlevlim2):
            coeflev = (nlevlim2 - (ind_lev - ind_ger)) / (nlevlim2 - nlevlim1) # the later is the emergence, the smaller is coeflev
        else:
            coeflev = 0

        # Density reduction because of late emergence after germination
        densite_list[i] = densite_list[ind_ger] * coeflev

    
    if lev[i] == 1:
        # Thermal time of day i
        if codetemp == 1:
            udev = effective_temperature(temp, tdmax, tdmin, tcxstop)
        elif codetemp == 2:
            udev = effective_temperature(tcult_prev, tdmax, tdmin, tcxstop)

        # Cumulated thermal time
        somfeuille = somfeuille_prev + udev


    if (lev_i_prev == 1) and (let_i_prev == 0): # frost sensitivty period

        # Leaves number
        if somfeuille > phyllotherme:
            nbfeuille = nbfeuille + 1
            somfeuille  = somfeuille - phyllotherme

        # Frost sensitivity period end
        if nbfeuille > nbfgellev:
            let[i:len(let)] = 1
        elif (codgellev == 2):
            # Emergence date
            ind_lev = np.where(lev > 0)[0][0]

            # Frost stress affecting density
            fgellev = frost_stress(
            temp_min_prev,
            tgellev90,
            tgellev10,
            tletale,
            tdebgel,
            )

            # Density reduction
            densite_list[i] = min(densite_list[i], densite_list[ind_lev] * fgellev) # fgellev not applied to densite(i-1) but to densite(lev)

    return moist, nbjgrauto, nbjhumec, humirac, somger, ger, zrac, elong, lev, coeflev, densite_list, let, udev, somfeuille, nbfeuille, fgellev, somelong


def budding(i, codedormance, q10, temp_max_list, temp_min_list, jvc, lev_i_prev, tdmindeb, tdmaxdeb, hourly_temp, gdh_prev,stdordebour):
    '''
    This modules computes budding for ligneous plants : dormancy break and post-dormancy period.
    See section 3.3.4.2 and 3.4.3 of STICS book.
    '''
    

    # Dormancy break
    if codedormance == 1:
        findorm_i = 1 

    else:
        # Degree days based on Q10 and air temperature
        cu = sum(
            [
                q10 ** (-temp_max_list[j] / 10)
                + q10 ** (-temp_min_list[j] / 10)
                for j in range(0, i + 1)
            ]
        )

        # Dormancy break
        findorm_i = np.where(
            cu > jvc, 1, 0
        )

    # Budding
    if (findorm_i == 1) and (
        lev_i_prev == 0
    ): 

        # Hourly temperatures
        thn = [
            0
            if t < tdmindeb
            else (
                tdmaxdeb - tdmindeb
                if t > tdmaxdeb
                else t - tdmindeb
            )
            for t in hourly_temp
        ]

        # Growing degree hours
        gdh = gdh_prev + sum(thn)
        
        # Budding
        lev_i = int(
            gdh > stdordebour
        )

    elif (
        lev_i_prev == 1
    ):
        lev_i = 1
    
    return findorm_i, cu, lev_i_prev, thn, gdh, lev_i


def development_temperature(tcult_prev, temp, tdmax, tdmin, tcxstop, coderetflo, stressdev, turfac_prev, codetemp, somtemp_prev, drp_i_prev):
    '''
    This module computes development temperature for phenology.
    See Section 3.3.2 of STICS book.
    '''
    
    # Air or crop temperature
    if codetemp == 1:
        udevcult = effective_temperature(temp, tdmax, tdmin, tcxstop)
    elif codetemp == 2:
        udevcult = effective_temperature(tcult_prev, tdmax, tdmin, tcxstop)
    
    # Crop sensitive to water stress during vegetative phase
    if (coderetflo == 1) & (drp_i_prev == 0):
        udevcult = udevcult * (stressdev * turfac_prev + 1 - stressdev)

    tdevelop = 2 ** (udevcult / 10)
    somtemp = somtemp_prev + tdevelop

    return udevcult, somtemp


def effective_temperature(temp, tdmax, tdmin, tcxstop):

    # Development temperature
    udevcult = max(0, temp - tdmin)
    if temp > tdmax:
        if tcxstop >= 100:
            udevcult = tdmax - tdmin
        else:
            udevcult = max(0,(tdmax - tdmin) * (temp - tcxstop) / (tdmax - tcxstop))
    
    return udevcult


def photoperiod_effect(herbaceous, lev_i, findorm_i, drp_i_prev, sensiphot, phoi, phosat, phobase):
    '''
    This module computes the photoperiod effect on development.
    It has an effect on a specific period : between emergence (lev) and fruit/grain filling start (drp) for herbaceous, between dormancy break and fruit filling start (drp) for ligneous.
    See section 3.3.3 of STICS book.
    '''

    if (
        herbaceous
        and (lev_i == 1)
        and (drp_i_prev == 0)
    ):
        rfpi = (1 - sensiphot) * (
            phoi - phosat
        ) / (phosat - phobase) + 1

    elif (
        (not herbaceous)
        and (findorm_i == 1)
        and (drp_i_prev == 0)
    ):
        rfpi = (
            phoi - phosat
        ) / (phosat - phobase) + 1
    else:
        rfpi = 1
    
    return rfpi

def vernalisation_effect(herbaceous, codebfroid, ger_i, tfroid, tcult_prev, ampfroid, jvi_list, jvcmini, jvc, findorm_i, doy, julvernal, codeperenne, rfvi, vernalisation_ongoing):
    '''
    This module computes the vernalisation effect on development.
    For herbaceous plants, vernalisation effect is computed from germination, it varies between 0 and 1 until chilling requirements are met.
    For ligneous plants, vernalisation effect is equal to 0 before dormancy break and 1 after. 
    See section 3.3.4 of STICS book.
    '''

    jvi = 0

    # Activate vernalisation on germination for annual, and on julvernal day for perennial crops
    condition_annual = (codeperenne == 1) & (codebfroid == 2) & (ger_i == 1)
    condition_perennial = (codeperenne == 2) & (codebfroid != 1)
    
    if doy == julvernal:
        vernalisation_ongoing = True

    if vernalisation_ongoing & (condition_annual | condition_perennial):
        if herbaceous:
            # Number of vernalising days
            jvi = max(
                (
                    1
                    - (
                        (tfroid- tcult_prev)
                        / ampfroid
                    )
                    ** 2
                ),
                0,
            ) 
            
            # Vernalisation effect = number of vernalising days / number of days needed
            rfvi = (
                max(
                    0,
                    (jvi_list.sum() + jvi - jvcmini)
                    / (jvc - jvcmini),
                )
                if (jvi_list.sum() + jvi - jvcmini)
                / (jvc - jvcmini)
                < 1
                else 1
            ) 
            
        elif findorm_i == 1:
            rfvi = 1
        
        if rfvi == 1:
            vernalisation_ongoing = False

    return rfvi, jvi, vernalisation_ongoing


def phenological_stage(lev_i, udevcult, rfpi, rfvi, sum_upvt_post_lev_prev, stlevamf, stamflax, stlevdrp, stflodrp,
               stdrpdes, codeindetermin, stdrpmat, stdrpnou, codlainet, stlaxsen, stsenlan, lan_i_prev, somcour_prev):
    '''
    This module computes the phenological stage.
    Temperature acts on development from germination for herbaceous plants, and from dormancy break for ligneous plants.
    A stage is reached when development temperature, reduced by vernalisation and photoperiod effect, exceeds the degree days needed.
    '''

    # Development temperature reduced by vernalisation and photoperiod effects
    upvt_post_lev = (
        lev_i
        * udevcult
        * rfpi
        * rfvi
    )

    # Cumulated degree days from emergence (lev)
    sum_upvt_post_lev = sum_upvt_post_lev_prev + upvt_post_lev

    # amf stage : 1 = amf stage reached, 0 else
    amf_i = 1 if sum_upvt_post_lev >= stlevamf else 0
    
    # lax stage : 1 = lax stage reached, 0 else
    lax_i = 1 if sum_upvt_post_lev >= stlevamf + stamflax else 0
    
    # flo stage : 1 = flo stage reached, 0 else
    flo_i = 1 if sum_upvt_post_lev >= stlevdrp - stflodrp else 0

    # drp stage : 1 = drp stage reached, 0 else
    drp_i = 1 if sum_upvt_post_lev >= stlevdrp else 0

    # debdes stage : 1 = debdes stage reached, 0 else
    debdes_i = 1 if sum_upvt_post_lev >= stlevdrp + stdrpdes else 0
    
    if codlainet == 1:
        sen_i = 1 if sum_upvt_post_lev >= stlevamf + stamflax + stlaxsen else 0
    else:
        sen_i = 0
    
    if sen_i == 1:
        somcour = somcour_prev + upvt_post_lev
    else:
        somcour = 0

    if lan_i_prev == 1:
        lan_i = 1
    if codlainet == 1:
        lan_i = 1 if sum_upvt_post_lev >= stlevamf + stamflax + stlaxsen + stsenlan else 0
    else:
        lan_i = 0

    if codeindetermin == 1:
        # mat stage : 1 = mat stage reached, 0 else
        mat_i = 1 if sum_upvt_post_lev >= stlevdrp + stdrpmat else 0
    else:
        mat_i = lax_i

    if codeindetermin == 2:
        # nou stage : 1 = nou stage reached, 0 else
        nou_i = 1 if sum_upvt_post_lev >= stlevdrp + stdrpnou else 0

    if codeindetermin == 1:
        return upvt_post_lev, sum_upvt_post_lev, amf_i, lax_i, flo_i, drp_i, debdes_i, mat_i, sen_i, lan_i, somcour
    elif codeindetermin == 2:
        return upvt_post_lev, sum_upvt_post_lev, amf_i, lax_i, flo_i, drp_i, debdes_i, mat_i, sen_i, lan_i, nou_i, somcour



def phenological_stage_dates(lev, amf, debdes, drp, nou, flo, findorm, mat, lax, codeindetermin, codeperenne):
    '''
    This module retrieves dates (julian days) of phenological stages, and associated BBCH codes.
    '''

    # Emergence day
    ind_lev = np.where(lev > 0)[0][0]

    # AMF day
    ind_amf = np.where(amf > 0)[0][0]

    # End of setting day
    if codeindetermin == 2:
        ind_nou = np.where(nou > 0)[0][0]
    
    # DRP day
    if drp.max() == 1:
        ind_drp = np.where(drp > 0)[0][0]
    else:
        ind_drp = 0

    # DEBDES day
    if debdes.max() == 1:
        ind_debdes = np.where(debdes > 0)[0][0]
    else:
        ind_debdes = 0

    # Maturation day
    if mat.max() == 1:
        ind_mat = np.where(mat > 0)[0][0]
    else:
        ind_mat = 0

    # BBCH code associated to each phenological stage
    bbch_list = np.zeros(len(lev)) -1
    if codeperenne == 2:
        bbch_list[findorm > 0] = -0.5
        bbch_list[lev > 0] = 0
        bbch_list[amf > 0] = 3
        bbch_list[flo > 0] = 6
        if drp.max() == 1:
            bbch_list[drp > 0] = 7
        if codeindetermin == 2:
            bbch_list[nou > 0] = 7.5
        if debdes.max() == 1:
            bbch_list[debdes > 0] = 7.5
        if lax.max() == 1:
            bbch_list[lax > 0] = 7.75
        if mat.max() == 1:
            bbch_list[mat > 0] = 8.9
    else:
        bbch_list[lev > 0] = 0
        bbch_list[amf > 0] = 3
        bbch_list[lax > 0] = 4
        bbch_list[flo > 0] = 6
        if drp.max() == 1:
            bbch_list[drp > 0] = 7
        if codeindetermin == 2:
            bbch_list[nou > 0] = 7.5
        if mat.max() == 1:
            bbch_list[mat > 0] = 8.9
        if debdes.max() == 1:
            bbch_list[debdes > 0] = 9.2

    if codeindetermin == 2:
        return bbch_list, ind_drp, ind_lev, ind_amf, ind_debdes, ind_mat, ind_nou
    elif codeindetermin == 1:
        return bbch_list, ind_drp, ind_lev, ind_amf, ind_mat, ind_debdes