import numpy as np
from pystics.modules.canopy_microclimate import wind_profile
from pystics.modules.irrigation import automated_irrig
from pystics.modules.foliage_water_interception import foliage_water
from pystics.modules.canopy_microclimate import radiation_interception, iterative_calculation, wind_profile
from pystics.modules.growth.plant_height import plant_height
from pystics.modules.water_potential import plant_soil_water_potential
from pystics.modules.development import development_temperature, emergence, budding, photoperiod_effect, vernalisation_effect, phenological_stage
from pystics.modules.growth.roots import root_density, root_growth
from pystics.modules.growth.leaf_growth import leaf_growth
from pystics.modules.senescence import senescence, senescence_stress
from pystics.modules.growth.thermal_stress_indices import thermal_stress_on_grain_filling, thermal_stress_on_biomass_broduction, frost_stress_on_fruit_number
from pystics.modules.growth.biomass import biomass_production, fallen_biomass
from pystics.modules.growth.yield_formation_forage import yield_formation_forage
from pystics.modules.water.soil_evaporation import potential_soil_evap_beer_law, actual_soil_evaporation, soil_layers_contribution_to_evaporation
from pystics.modules.water.transpiration import potential_transpiration_coef, actual_transpiration, soil_contribution_to_transpiration
from pystics.modules.water.water_stress import water_stress
from pystics.modules.water.water_balance import water_balance

def compute_outputs_day_zero(outputs, crop, soil, manage, constants, initial, station, lracz, hur, wi, esz, epz, tsol, humirac, humpotsol, psisol, racinepsi, amplz, ndebsen, hurlim, aevap, s, fco2s, durviei, gamma, msresiduel, fco2):
    ''''
    This module computes all variables for day zero.
    It is very similar to run_pystics_simulation function which calls the same modules in a daily loop.
    '''

    # Plant density
    outputs.loc[0,'densite'] = manage.DENSITESEM

    ###################
    ### Development
    ###################

    ## 0. Sowing / planting
    if outputs.loc[0,'plt'] == 0:
        if (outputs.loc[0,'doy'] == manage.IPLT0):
            outputs['plt'] = 1

    ## 1. Emergence / budding
    if (crop.CODEPERENNE == 1) | (crop.HERBACEOUS):
        outputs['moist'], outputs.loc[0,'nbjgrauto'], outputs.loc[0,'nbjhumec'], humirac, outputs.loc[0,'somger'], outputs['ger'], outputs.loc[0,'zrac'], outputs.loc[0,'elong'], outputs['lev'], outputs.loc[0,'coeflev'], outputs['densite'], outputs['let'], outputs.loc[0,'udevcult_lev'], outputs.loc[0,'somfeuille'], outputs.loc[0,'nbfeuille'], outputs.loc[0,'fgellev'] = emergence(0, outputs['densite'].array, outputs['lev'].array, outputs.loc[0,'lev'], outputs['ger'].array, outputs.loc[0,'ger'], outputs['moist'].array, outputs.loc[0,'moist'], outputs['let'].array, manage.PROFSEM, hur[0], humpotsol[0], crop.PROPJGERMIN, crop.NBJGERLIM, crop.TDMAX, crop.TDMIN, crop.TGMIN,
                                         tsol, outputs.loc[0,'nbjhumec'], soil.HMIN, crop.SENSRSEC, soil.HCC, outputs.loc[0,'somger'], crop.STPLTGER, tsol[0], manage.DENSITESEM, crop.CODEHYPO,
                                         outputs.loc[0,'zrac'], crop.BELONG, crop.CELONG, crop.ELMAX, outputs.loc[0,'temp'], crop.NLEVLIM1, crop.NLEVLIM2, crop.TCXSTOP, outputs.loc[0,'somfeuille'], outputs.loc[0,'nbfeuille'], outputs.loc[0,'temp_min'],
                                         crop.TGELLEV90, crop.TGELLEV10, crop.TLETALE, crop.TDEBGEL, crop.PHYLLOTHERME, crop.NBFGELLEV, humirac, soil.DEPTH, crop.CODGELLEV, outputs.loc[0,'let'])
    elif crop.CODEPERENNE == 2:
        outputs.loc[0,'findorm'], outputs.loc[0,'cu'], outputs.loc[0,'lev'], outputs.loc[0,'thn'], outputs.loc[0,'gdh'], outputs.loc[0,'lev'] = budding(0, crop.CODEDORMANCE, crop.Q10, outputs['Temp_max'], outputs['Temp_min'], crop.JVC, outputs.loc[0,'lev'], crop.TDMINDEB, crop.TDMAXDEB, outputs.loc[0,'hourly_temp'], outputs.loc[0,'gdh'], crop.STDORDEBOUR)

    ## 2. Development temperature
    if outputs.loc[0,'plt'] == 1:
        outputs.loc[0,'udevcult'], outputs.loc[0,'somtemp'] = development_temperature(outputs.loc[0,'temp'], outputs.loc[0,'temp'], crop.TDMAX, crop.TDMIN, crop.TCXSTOP, crop.CODERETFLO, crop.STRESSDEV, 1, crop.CODETEMP, outputs.loc[0,'somtemp'])

    ## 3. Effect of photoperiod
    if crop.CODEPHOT == 1:
        outputs.loc[0,'rfpi'] = photoperiod_effect(crop.HERBACEOUS, outputs.loc[0,'lev'], outputs.loc[0,'findorm'], outputs.loc[0,'drp'], crop.SENSIPHOT, outputs.loc[0,'phoi'], crop.PHOSAT, crop.PHOBASE)
    else:
        outputs.loc[0,'rfpi'] = 1

    ## 4. Effect of vernalisation
    outputs.loc[0,'rfvi'], outputs.loc[0,'jvi'] = vernalisation_effect(crop.HERBACEOUS, crop.CODEBFROID, outputs.loc[0,'ger'], crop.TFROID, outputs.loc[0,'temp'], crop.AMPFROID, outputs.loc[0,'jvi'], crop.JVCMINI, crop.JVC, outputs.loc[0,'findorm'])

    ## 5. Compute of phenological stage this day
    if crop.CODEINDETERMIN == 1:
        outputs.loc[0,'upvt_post_lev'], outputs.loc[0,'sum_upvt_post_lev'], outputs.loc[0,'amf'], outputs.loc[0,'lax'], outputs.loc[0,'flo'], outputs.loc[0,'drp'], outputs.loc[0,'debdes'], outputs.loc[0,'mat'], outputs.loc[0,'sen'], outputs.loc[0,'lan'], outputs.loc[0,'sum_upvt'] = phenological_stage(outputs.loc[0,'lev'], outputs.loc[0,'udevcult'], outputs.loc[0,'rfpi'], outputs.loc[0,'rfvi'], outputs.loc[0,'sum_upvt_post_lev'], crop.STLEVAMF,
                                                crop.STAMFLAX, crop.STLEVDRP, crop.STFLODRP, crop.STDRPDES, crop.CODEINDETERMIN, crop.STDRPMAT, crop.STDRPNOU, crop.CODLAINET, crop.STLAXSEN, crop.STSENLAN)
    elif crop.CODEINDETERMIN == 2:  
        outputs.loc[0,'upvt_post_lev'], outputs.loc[0,'sum_upvt_post_lev'], outputs.loc[0,'amf'], outputs.loc[0,'lax'], outputs.loc[0,'flo'], outputs.loc[0,'drp'], outputs.loc[0,'debdes'], outputs.loc[0,'mat'], outputs.loc[0,'sen'], outputs.loc[0,'lan'], outputs.loc[0,'sum_upvt'], outputs.loc[0,'nou'] = phenological_stage(outputs.loc[0,'lev'], outputs.loc[0,'udevcult'], outputs.loc[0,'rfpi'], outputs.loc[0,'rfvi'], outputs.loc[0,'sum_upvt_post_lev'], crop.STLEVAMF,
                                                crop.STAMFLAX, crop.STLEVDRP, crop.STFLODRP, crop.STDRPDES, crop.CODEINDETERMIN, crop.STDRPMAT, crop.STDRPNOU, crop.CODLAINET, crop.STLAXSEN, crop.STSENLAN)
    
    #############
    ### Leaf growth
    #############
    outputs.loc[0,'deltai'], outputs.loc[0,'deltai_dev'], outputs.loc[0,'deltaidens'], outputs.loc[0,'deltai_t'], outputs.loc[0,'ulai'], outputs.loc[0,'deltai_stress'], outputs.loc[0,'efdensite'], outputs.loc[0,'vmax'], outputs.loc[0,'lai'], outputs.loc[0,'mafeuilverte'], outputs.loc[0,'dltaisen'], outputs.loc[0,'dltaisenat'], outputs.loc[0,'laisen'], outputs.loc[0,'lan'], outputs.loc[0,'sen'] = leaf_growth(outputs.loc[0,'lev'], outputs.loc[0,'lax'], outputs.loc[0,'sum_upvt'], crop.STLEVAMF, crop.VLAIMAX, crop.STAMFLAX, crop.UDLAIMAX, crop.DLAIMAX, crop.PENTLAIMAX,
                outputs.loc[0,'tcult'], crop.TCXSTOP, crop.TCMAX, crop.TCMIN, crop.ADENS, crop.BDENS, outputs.loc[0,'densite'], outputs.loc[0,'turfac'], outputs.loc[0,'phoi'], outputs.loc[0,'phoi'], crop.TIGEFEUIL,
                crop.PHOBASE, outputs.loc[0,'rfpi'], crop.DLAIMIN, outputs.loc[0,'lai'], crop.LAICOMP, crop.STOPFEUILLE, outputs.loc[0,'vmax'], crop.CODLAINET, outputs.loc[0,'dltaisenat'], outputs.loc[0,'fstressgel'], outputs.loc[0,'laisen'], outputs.loc[0,'lax'], crop.SLAMIN, crop.SLAMAX, outputs.loc[0,'dltamsen'], outputs.loc[0,'dltaisen'], outputs.loc[0,'sen'], outputs.loc[0,'lan'])
    
    if crop.CODLAINET == 1:
        outputs.loc[0,'lai'] = outputs.loc[0,'lai'] + outputs.loc[0,'deltai'] - outputs.loc[0,'dltaisenat']
    else:
        outputs.loc[0,'lai'] = outputs.loc[0,'lai'] + outputs.loc[0,'deltai'] - outputs.loc[0,'dltaisen']
    outputs.loc[0,'mafeuilverte'] = outputs.loc[0, 'lai'] / crop.SLAMAX * 100

   
    #############################
    ### Intercepted radiation ###
    #############################
    outputs.loc[0,'raint'] = radiation_interception(constants.PARSURRG, outputs.loc[0,'trg'], crop.EXTIN, outputs.loc[0,'lai'])

    #############################
    ### Thermal stress on RUE ###
    #############################
    outputs.loc[0,'ftemp'] = thermal_stress_on_biomass_broduction(outputs.loc[0,'temp'], crop.TEMIN, crop.TEMAX, crop.TEOPT, crop.TEOPTBIS)

    ##########################
    ### Biomass production ###
    ##########################
    outputs.loc[0,'dltams'], outputs.loc[0,'masec'], outputs.loc[0,'ebmax'], outputs.loc[0,'dltafv'] = biomass_production(outputs.loc[0,'masec'], outputs.loc[0,'raint'], outputs.loc[0,'lev'], outputs.loc[0,'amf'], outputs.loc[0,'drp'], crop.EFCROIJUV, crop.EFCROIREPRO, crop.EFCROIVEG,
                            constants.COEFB, outputs.loc[0,'swfac'], fco2, outputs.loc[0,'ftemp'], outputs.loc[0, "deltai"], crop.SLAMAX)

    #############
    ### Senescence
    #############
    outputs.loc[0,'fstressgel'] = senescence_stress(outputs.loc[0,'lev'], outputs.loc[0,'ulai'], crop.VLAIMAX, outputs.loc[0,'temp_min'], crop.TGELJUV10, crop.TGELJUV90, crop.TGELVEG10, crop.TGELVEG90, crop.TLETALE, crop.TDEBGEL, crop.CODGELJUV, crop.CODGELVEG)

    outputs['dayLAIcreation'], outputs['durage'], outputs['senstress'], outputs['tdevelop'], outputs['durvie'], outputs.loc[0,'dltaisen'], outputs.loc[0,'somsenreste'], ndebsen, outputs.loc[0,'somtemp'], outputs.loc[0,'dltamsen'], outputs.loc[0,'deltamsresen'], outputs.loc[0,'msres'], outputs.loc[0,'msresjaune'], outputs.loc[0,'durvie_n']   = senescence(0, outputs['lev'], outputs.loc[0,'ulai'], outputs.loc[0,'somtemp'], crop.VLAIMAX, durviei, crop.DURVIEF, 1, outputs.loc[0,'udevcult'], outputs.loc[0,'fstressgel'],
            outputs['dayLAIcreation'].array, outputs['senstress'].array, outputs['tdevelop'].array, outputs['durvie'].array, outputs['durage'].array, outputs['deltai'].array, outputs.loc[0,'somsenreste'], initial.LAI0, ndebsen, outputs['dltafv'], crop.RATIOSEN, crop.forage, outputs.loc[0,'msres'], msresiduel, outputs.loc[0,'msresjaune'], crop.CODLAINET)
    
    outputs.loc[0,'laisen'] = outputs.loc[0,'dltaisen']

    ######################
    ### Fallen biomass ###
    ######################
    outputs.loc[0,'dltamstombe'], outputs.loc[0,'mafeuiltombe'], outputs.loc[0,'mafeuiljaune'] = fallen_biomass(crop.ABSCISSION, outputs.loc[0,'dltamsen'], outputs.loc[0,'mafeuiltombe'], outputs.loc[0,'mafeuiljaune'])

    ###################
    ### Root growth ###
    ###################
    outputs.loc[0,'zrac'], outputs.loc[0,'deltaz'], outputs.loc[0,'deltaz_t'], outputs.loc[0,'deltaz_stress'], outputs.loc[0,'efda'], outputs.loc[0,'humirac_mean'], outputs.loc[0,'znonli'] = root_growth(outputs.loc[0,'ger'], outputs.loc[0,'lax'], crop.CODEPERENNE, outputs.loc[0,'findorm'], manage.PROFSEM, outputs.loc[0,'zrac'], crop.CODETEMPRAC,
            outputs.loc[0,'tcult'], crop.TCMAX, crop.TCMIN, crop.TGMIN, crop.CROIRAC, hur[0], tsol[0], soil.HMIN, crop.SENSRSEC, soil.DEPTH, constants.DASEUILBAS, constants.DASEUILHAUT, crop.CONTRDAMAX, soil.DAF, outputs.loc[0,'lev'], soil.HCC, crop.HERBACEOUS, outputs.loc[0, 'znonli'])

    ####################
    ### Root density ###
    ####################
    lracz[0], outputs.loc[0,'cumlracz'], outputs.loc[0,'zdemi'], humirac[0] = root_density(lracz[0], outputs.loc[0,'zrac'], outputs.loc[0,'znonli'], soil.DEPTH, crop.ZPRLIM, crop.ZPENTE, outputs.loc[0,'ger'], crop.CODEPERENNE, constants.LVOPT, s, manage.PROFSEM, 
                                                                                            hur[0], soil.HMIN, humirac[0])

    #######################
    ### Water potential ###
    #######################
    humpotsol[0], psisol[0], racinepsi[0], outputs.loc[0,'psibase'] = plant_soil_water_potential(outputs.loc[0,'lev'], initial.ZRAC0, soil.DEPTH, constants.PSIHUCC, constants.PSIHUMIN, soil.DAF, soil.HMIN, soil.HCC, crop.POTGERMI, manage.PROFSEM, outputs.loc[0,'zrac'], lracz[0], hur[0], racinepsi[0])

    ############################################################
    ### Compute irrigation when it is automatically computed ###
    ############################################################
    if manage.CODECALIRRIG == 1:
        outputs.loc[0,'airg'] = automated_irrig(manage.CODEDATE_IRRIGAUTO, outputs.loc[0,'compute_airg'], manage.RATIOL, 1, manage.EFFIRR, manage.DOSEIRRIGMIN, manage.DOSIMX, soil.HCC, hur[0], outputs.loc[0,'zrac'], outputs.loc[0,'doy'], manage.IPLT0, constants.IRRLEV) # swfac set to 1 on day 0

    ####################################
    ### Water intercepted by foliage ###
    ####################################
    if crop.CODEINTERCEPT == 1:
        outputs.loc[0,'stemflow'], outputs.loc[0,'mouill'] = foliage_water(manage.CODLOCIRRIG, crop.STEMFLOWMAX, crop.KSTEMFLOW, outputs.loc[0,'airg'], outputs.loc[0,'trr'], crop.MOUILLABIL, outputs.loc[0,'lai'])
    elif crop.CODEINTERCEPT == 2:
        outputs.loc[0,'mouill'] = 0
    
    ########################
    ### Soil evaporation ###
    ########################

    # Potential
    outputs.loc[0,'eos'] = potential_soil_evap_beer_law(outputs.loc[0,'etp'], crop.EXTIN, outputs.loc[0,'lai'])

    # Actual
    outputs.loc[0,'esol'], outputs.loc[0,'sumes0'], outputs.loc[0,'sumes1'], outputs.loc[0,'sumes2'], outputs.loc[0,'sesj0'], outputs.loc[0,'ses2j0'], outputs.loc[0,'smes02'], outputs.loc[0,'nstoc'], outputs.loc[0,'stoc'] = actual_soil_evaporation(outputs.loc[0,'trr'], outputs.loc[0,'airg'], outputs.loc[0,'eos'], aevap, soil.Q0, outputs.loc[0,'sumes0'], outputs.loc[0,'sumes1'], outputs.loc[0,'sumes2'], outputs.loc[0,'sesj0'], outputs.loc[0,'ses2j0'], outputs.loc[0,'smes02'], outputs.loc[0,'nstoc'], outputs.loc[0,'stoc']) 

    # Soil layers contribution
    outputs.loc[0,'esol'], esz[0], wi[0] = soil_layers_contribution_to_evaporation(outputs.loc[0,'esol'], soil.ZESX, hurlim, soil.CFES, soil.HCC, wi[0], hur[0], esz[0])

    ###########################
    ### Plant transpiration ###
    ###########################

    # Potential
    outputs.loc[0,'eo'], outputs.loc[0,'eop'], outputs.loc[0,'emd'], outputs.loc[0,'edirect'], outputs.loc[0,'directm'] = potential_transpiration_coef(outputs.loc[0,'lai'], outputs.loc[0,'etp'], crop.KMAX, crop.CODEINTERCEPT, outputs.loc[0,'eos'], outputs.loc[0,'esol'], crop.EXTIN, outputs.loc[0,'mouill'])

    #######################
    ### Water balance #####
    #######################
    hur[0], outputs.loc[0,'drained_water'] = water_balance(outputs.loc[0,'trr'], outputs.loc[0,'airg'], outputs.loc[0,'mouill'], soil.DEPTH, hur[0], np.array([0 for i in range(len(outputs))]), esz[0], soil.HCC, hurlim)

    ###########################
    ### Plant transpiration ###
    ###########################

    # Actual
    outputs.loc[0,'ep'] = actual_transpiration(1, outputs.loc[0,'eop'])
    epz[0] = soil_contribution_to_transpiration(outputs.loc[0,'ep'], epz[0], soil.DEPTH, hur[0], soil.HMIN, lracz[0], outputs.loc[0,'lev'])


    outputs.loc[0,'et'] = outputs.loc[0,'esol'] + outputs.loc[0,'ep']

    ######################
    ### Water stresses ###
    ######################
    outputs.loc[0,'teta'], outputs.loc[0,'swfac'], outputs.loc[0,'tetstomate'], outputs.loc[0,'turfac'], outputs.loc[0,'teturg'], outputs.loc[0,'senfac'], outputs.loc[0,'tetsen'], outputs.loc[0,'resrac'] = water_stress(outputs.loc[0,'eop'], crop.SWFACMIN, crop.RAPSENTURG, outputs.loc[0,'cumlracz'], crop.PSITURG, crop.PSISTO, crop.RAYON, outputs.loc[0,'zrac'], outputs.loc[0,'lev'], soil.HMIN, hur[0], manage.PROFSEM, soil.DEPTH)

    ####################
    ### Plant height ###
    ####################
    outputs.loc[0,'hauteur'] = plant_height(crop.HAUTMAX, crop.HAUTBASE, outputs.loc[0,'lai'], outputs.loc[0,'laisen'])

    ####################
    ### Wind profile ###
    ####################
    outputs.loc[0,'dh'], outputs.loc[0,'z0'] = wind_profile(soil.Z0SOLNU, outputs.loc[0, "hauteur"])

    ###################################################################
    ### Iterative calculation of crop temperature and net radiation ###
    ###################################################################
    outputs.loc[0,'rnet'], outputs.loc[0,'rglo'], outputs.loc[0,'albedolai'], outputs.loc[0,'albsol'], outputs.loc[0,'tcult'], outputs.loc[0,'tcultmax'], outputs.loc[0,'converge'] = iterative_calculation(outputs.loc[0,'temp'], outputs.loc[0,'lev'], outputs.loc[0,'temp_max'], outputs.loc[0,'temp_min'], outputs.loc[0,'et'], outputs.loc[0,'z0'], soil.ALBEDO, hur[0,0], soil.HMINF_1, soil.HCCF_1, station.ALBVEG, outputs.loc[0,'lai'], outputs.loc[0,'trg'], outputs.loc[0,'tpm'], outputs.loc[0,'fracinsol'], station.CODERNET)

    ########################
    ### Soil temperature ###
    ########################
    outputs.loc[0,'amplsurf'] = outputs.loc[0,'temp_max'] - outputs.loc[0,'temp_min']
    for z in range(soil.DEPTH):
        amplz[0,z] = outputs.loc[0,'amplsurf'] * np.exp(-z * (7.272*0.00001 / (2*constants.DIFTHERM))**(1/2))
        tsol[0,z] = outputs.loc[0, "tcult"]


    ##############################
    ### Thermal stress indices ###
    ##############################
    outputs.loc[0,'ftempremp'] = thermal_stress_on_grain_filling(crop.CODETREMP, outputs.loc[0,'temp_min'], outputs.loc[0,'tcultmax'], crop.TMINREMP, crop.TMAXREMP)
    outputs.loc[0,'fgelflo'] = frost_stress_on_fruit_number(outputs.loc[0,'temp_min'], crop.TGELFLO90, crop.TGELFLO10, crop.TLETALE, crop.TDEBGEL, crop.CODGELFLO)
        

    # For forage crops only
    if crop.forage:

        # Specific outputs
        outputs.loc[0,'msneojaune'] = outputs.loc[0,'mafeuiljaune'].copy()
        outputs.loc[0,'masec'] = outputs.loc[0,'masec'] - outputs.loc[0,'dltamstombe'] # should this also be done for annual plants ?
        outputs.loc[0,'masectot'] = initial.MASEC0
        outputs.loc[0,'masecneo'] = outputs.loc[0,'dltams']
        
        #######################
        ### Yield formation ###
        #######################
        outputs.loc[0,'mafruit'], outputs.loc[0,'msrec_fou'] = yield_formation_forage(outputs.loc[0,'masec'], outputs.loc[0,'msresjaune'], outputs.loc[0,'msneojaune'], msresiduel)

    return outputs, humirac, hur, wi, esz, tsol, amplz, ndebsen

    







