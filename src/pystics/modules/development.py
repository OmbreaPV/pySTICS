import numpy as np
from pystics.modules.growth.thermal_stress_indices import frost_stress

def emergence_macro(i, outputs, crop, soil, manage, tsol, humirac, hur, humpotsol):
    outputs['moist'], outputs.loc[i,'nbjgrauto'], outputs.loc[i,'nbjhumec'], humirac, outputs.loc[i,'somger'], outputs['ger'], outputs.loc[i,'zrac'], outputs.loc[i,'elong'], outputs['lev'], outputs.loc[i,'coeflev'], outputs['densite'], outputs['let'], outputs.loc[i,'udevcult_lev'], outputs.loc[i,'somfeuille'], outputs.loc[i,'nbfeuille'], outputs.loc[i,'fgellev'] = emergence(i, outputs['densite'].array, outputs['lev'].array, outputs.loc[i-1,'lev'], outputs['ger'].array, outputs.loc[i-1,'ger'], outputs['moist'].array, outputs.loc[i-1,'moist'], outputs['let'].array, manage.PROFSEM, hur[i], humpotsol[i], crop.PROPJGERMIN, crop.NBJGERLIM, crop.TDMAX, crop.TDMIN, crop.TGMIN,
                                         tsol, outputs.loc[i-1,'nbjhumec'], soil.HMIN, crop.SENSRSEC, soil.HCC, outputs.loc[i-1,'somger'], crop.STPLTGER, tsol[i], manage.DENSITESEM, crop.CODEHYPO,
                                         outputs.loc[i,'zrac'], crop.BELONG, crop.CELONG, crop.ELMAX, outputs.loc[i-1,'tcult'], crop.NLEVLIM1, crop.NLEVLIM2, crop.TCXSTOP, outputs.loc[i-1,'somfeuille'], outputs.loc[i,'nbfeuille'], outputs.loc[i,'temp_min'],
                                         crop.TGELLEV90, crop.TGELLEV10, crop.TLETALE, crop.TDEBGEL, crop.PHYLLOTHERME, crop.NBFGELLEV, humirac, soil.DEPTH, crop.CODGELLEV, outputs.loc[i-1,'let'])

    return outputs, humirac


def emergence(i, densite_list, lev, lev_i_prev, ger, ger_i_prev, moist, moist_i_prev, let, profsem, hur_i, humpotsol_i, propjgermin, nbjgerlim, tdmax, tdmin,tgmin,
              tsol, nbjhumec_prev, hmin, sensrsec, hcc, somger_prev, stpltger, tsol_i, densitesem, codehypo,
              zrac, belong, celong, elmax, tcult_prev, nlevlim1, nlevlim2, tcxstop, somfeuille_prev, nbfeuille, temp_min,
            tgellev90,tgellev10, tletale, tdebgel, phyllotherme, nbfgellev, humirac, depth, codgellev, let_i_prev):
    '''
    This module computes emergence of herbaceous plants : moistening, germination and plantlet growth.
    See Section 3.4.1.3 of STICS book.

    Simplifications compared to STICS :
        - Soil crusting not implemented
    '''

    nbjgrauto, nbjhumec, somger, lev_i, elong, coeflev, udevcult_lev, somfeuille, fgellev = 0,0,0,0,0,0,0,0,1

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
        for z_index in sb:
            if hur_i[z_index] > hmin[z_index]:  
                humirac[i,z_index] = 1
            else:
                humirac[i,z_index] = min(1, max(0, sensrsec * hur_i[z_index] / hmin[z_index]))
            
        # Growing degree days to reach germination
        somger = somger_prev + max(0, tsol_i[sb].mean() - tgmin) * humirac[i,sb].mean() # does not depend on sowing day because this module is called when plt=1
        if somger > stpltger:
            ger[i:len(ger)] = 1
            zrac = profsem

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
            lev_i = 1
        
        else:

            sb = range(max(0,int(profsem) - 1),int(profsem) + 2)
            hb = range(sb[0], min(max(sb[-1],int(zrac))+1, depth))

            # Water stress index affecting elongation
            for z_index in hb: 
                if hur_i[z_index] > hmin[z_index]:  
                    humirac[i,z_index] = 1
                else:
                    humirac[i,z_index] = min(1, max(0, sensrsec * hur_i[z_index] / hmin[z_index]))

            # Plantlet elongation
            elong = elmax * (1 - np.exp(-(belong* np.nanmean(humirac[ind_ger:i,hb] * np.maximum(0,tsol[ind_ger:i,hb] - tgmin), axis=1).sum())**celong))

            # Emergence
            if elong > profsem:
                lev[i:len(lev)] = 1

    if (lev_i_prev == 0) and (lev_i == 1): # emergence day

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


    if (lev_i_prev == 1) and (let_i_prev == 0): # frost sensitivty period

        # Thermal time of day i
        udevcult_lev = max(0, tcult_prev - tdmin)
        if tcxstop >= 100:
            if tcult_prev > tdmax:
                udevcult_lev = tdmax - tdmin
        else:
            if tcult_prev > tdmax:
                udevcult_lev = max(0,(tdmax - tdmin) * (tcult_prev - tcxstop) / (tdmax - tcxstop))

        # Cumulated thermal time
        somfeuille = somfeuille_prev + udevcult_lev

        # Leaves number
        if somfeuille > phyllotherme:
            nbfeuille = nbfeuille + 1
            somfeuille  = somfeuille - phyllotherme

        # Frost sensitivity period end
        if nbfeuille >= nbfgellev:
            let[i:len(let)] = 1
        elif (codgellev == 2):
            # Emergence date
            ind_lev = np.where(lev > 0)[0][0]

            # Frost stress affecting density
            fgellev = frost_stress(
            temp_min,
            tgellev90,
            tgellev10,
            tletale,
            tdebgel,
            )

            # Density reduction
            densite_list[i] = min(densite_list[i], densite_list[ind_lev] * fgellev) # fgellev not applied to densite(i-1) but to densite(lev)

    return moist, nbjgrauto, nbjhumec, humirac, somger, ger, zrac, elong, lev, coeflev, densite_list, let, udevcult_lev, somfeuille, nbfeuille, fgellev


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


def development_temperature(tcult_prev, temp, tdmax, tdmin, tcxstop, coderetflo, stressdev, turfac_prev, codetemp, somtemp_prev):
    '''
    This module computes development temperature for phenology.
    See Section 3.3.2 of STICS book.
    '''
    
    # Air or crop temperature
    if codetemp == 1:
        temp_consid = temp 
    elif codetemp == 2:
        temp_consid = tcult_prev
    
    # Development temperature
    udevcult = max(0, temp_consid - tdmin)
    if temp_consid > tdmax:
        if tcxstop >= 100:
            udevcult = tdmax - tdmin
        else:
            udevcult = max(0,(tdmax - tdmin) * (temp_consid - tcxstop) / (tdmax - tcxstop))
    
    # Crop sensitive to water stress
    if coderetflo == 1:
        udevcult = udevcult * (stressdev * turfac_prev + 1 - stressdev)

    tdevelop = 2 ** (udevcult / 10)
    somtemp = somtemp_prev + tdevelop

    return udevcult, somtemp


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

def vernalisation_effect(herbaceous, codebfroid, ger_i, tfroid, tcult_prev, ampfroid, jvi_list, jvcmini, jvc, findorm_i):
    '''
    This module computes the vernalisation effect on development.
    For herbaceous plants, vernalisation effect is computed from germination, it varies between 0 and 1 until chilling requirements are met.
    For ligneous plants, vernalisation effect is equal to 0 before dormancy break and 1 after. 
    See section 3.3.4 of STICS book.
    '''

    rfvi, jvi = 1, 0
    if herbaceous & (codebfroid == 2):
        if ger_i == 1:

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
        else:
            rfvi = 0 
    elif (
        not herbaceous
    ):
        if codebfroid != 1:
            if findorm_i == 0:
                rfvi = 0

    return rfvi, jvi


def phenological_stage(lev_i, udevcult, rfpi, rfvi, sum_upvt_post_lev_prev, stlevamf, stamflax, stlevdrp, stflodrp,
               stdrpdes, codeindetermin, stdrpmat, stdrpnou, codlainet, stlaxsen, stsenlan):
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
    amf_i = np.where(
        sum_upvt_post_lev >= stlevamf, 1, 0
    )
    
    # lax stage : 1 = lax stage reached, 0 else
    lax_i = np.where(
        sum_upvt_post_lev
        >= stlevamf + stamflax,
        1,
        0,
    )
    
    # flo stage : 1 = flo stage reached, 0 else
    flo_i = np.where(
        sum_upvt_post_lev
        >= stlevdrp - stflodrp,
        1,
        0,
    )

    # drp stage : 1 = drp stage reached, 0 else
    drp_i = np.where(
        sum_upvt_post_lev >= stlevdrp, 1, 0
    )

    # debdes stage : 1 = debdes stage reached, 0 else
    debdes_i = np.where(sum_upvt_post_lev >= stlevdrp + stdrpdes, 1, 0)
    
    if codlainet == 1:
        sen_i =  np.where(
            sum_upvt_post_lev
            >= stlevamf + stamflax + stlaxsen,
            1,
            0,
        )
    else:
        sen_i = 0

    if codlainet == 1:
        lan_i =  np.where(
            sum_upvt_post_lev
            >= stlevamf + stamflax + stlaxsen + stsenlan,
            1,
            0,
        )
    else:
        lan_i = 0

    if codeindetermin == 1:

        # mat stage : 1 = mat stage reached, 0 else
        mat_i = np.where(
            sum_upvt_post_lev
            >= stlevdrp + stdrpmat,
            1,
            0,
        )
    else:
        mat_i = lax_i

    if codeindetermin == 2:

        # nou stage : 1 = nou stage reached, 0 else
        nou_i = np.where(
            sum_upvt_post_lev >= stlevdrp + stdrpnou,
            1,
            0,
        ) 

    if codeindetermin == 1:
        return upvt_post_lev, sum_upvt_post_lev, amf_i, lax_i, flo_i, drp_i, debdes_i, mat_i, sen_i, lan_i, sum_upvt_post_lev
    elif codeindetermin == 2:
        return upvt_post_lev, sum_upvt_post_lev, amf_i, lax_i, flo_i, drp_i, debdes_i, mat_i, sen_i, lan_i,  sum_upvt_post_lev, nou_i



def phenological_stage_dates(lev, amf, debdes, drp, nou, flo, findorm, mat, lax, sum_upvt_list,
                                stlevdrp, codeindetermin, codeperenne):
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