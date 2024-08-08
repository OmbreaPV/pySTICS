import pandas as pd
import numpy as np

from pystics.modules.weather_variables import compute_weather_variables, check_weather_variables
from pystics.params import StationParams, CropParams, ManageParams, SoilParams, Constants, InitialParams
from pystics.modules.initialization import initialize_outputs_df, initialize_soil_matrix, compute_constants
from pystics.modules.day_zero import compute_outputs_day_zero
from pystics.modules.irrigation import irrigation_dates_and_amounts, dates_irrig, automated_irrig
from pystics.modules.foliage_water_interception import foliage_water
from pystics.modules.canopy_microclimate import radiation_interception, iterative_calculation, wind_profile
from pystics.modules.growth.plant_height import plant_height
from pystics.modules.soil_temperature import soil_temperature
from pystics.modules.water_potential import plant_soil_water_potential
from pystics.modules.development import development_temperature, emergence_macro, budding, photoperiod_effect, vernalisation_effect, phenological_stage, phenological_stage_dates
from pystics.modules.growth.roots import root_density, root_growth
from pystics.modules.growth.leaf_growth import leaf_growth
from pystics.modules.senescence import senescence, senescence_stress
from pystics.modules.growth.thermal_stress_indices import thermal_stress_on_grain_filling, thermal_stress_on_biomass_broduction, frost_stress_on_fruit_number
from pystics.modules.growth.biomass import biomass_production, fallen_biomass
from pystics.modules.growth.harvest_criteria import harvest_criteria
from pystics.modules.growth.yield_formation_determinate import harvested_organs_number, frost_reduction_harvested_organs, carbon_harvest_index, harvested_organs_mass, grains_water_content
from pystics.modules.growth.radiative_stress import radiative_stress
from pystics.modules.forage_cutting import cutting 
from pystics.modules.growth.yield_formation_forage import yield_formation_forage
from pystics.modules.water.soil_evaporation import potential_soil_evap_beer_law, actual_soil_evaporation, soil_layers_contribution_to_evaporation
from pystics.modules.water.transpiration import potential_transpiration_coef, actual_transpiration, soil_contribution_to_transpiration
from pystics.modules.water.water_stress import water_stress
from pystics.modules.water.water_balance import water_balance
from pystics.exceptions import pysticsException
from pystics.modules.growth.biomass_partitioning import cumultative_partitioning
from pystics.modules.lethal_frost import lethal_frost


def run_pystics_simulation(weather : pd.DataFrame, crop : CropParams, soil : SoilParams, constants : Constants, manage : ManageParams, station : StationParams, initial : InitialParams) -> pd.DataFrame :
    '''
    This function runs a pySTICS simulation from weather and parameters (crop, soil, constants, practices, station and intial values).
    '''

    if manage.PROFSEM > soil.DEPTH:
        raise pysticsException(detail="Sowing depth is greater than soil depth.")
    
    # In case of annual plants, start of vernalisation can not be before sowing
    if (crop.CODEPERENNE == 1) & (crop.JULVERNAL < manage.IPLT0):
        crop.JULVERNAL = manage.IPLT0

    ################################
    ### Constants for simulation ###
    ################################
    
    hurlim, aevap, s, fco2, fco2s, durviei, gamma, msresiduel, lairesiduel = compute_constants(crop, soil, station, weather.co2[0], manage)

    ###################################
    ### Initialization of variables ###
    ###################################

    # Check weather inputs and consistency with options
    weather = check_weather_variables(weather, station.CODEETP)

    # Irrigation when given as inputs (weather dates must be consistent with manage.IRRIGATION_INTERVENTIONS)
    if (manage.CODECALIRRIG == 2) & (manage.CODEDATEAPPH2O == 2):
        weather = irrigation_dates_and_amounts(weather,manage)

    # Outputs dataframe (1 line = 1 day, 1 col = 1 output) based on weather dataframe
    outputs = initialize_outputs_df(weather, crop, manage, initial)

    # Soil matrixes corresponding to outputs with a depth dimension (1 line = 1 day, 1 col = 1 soil cm)
    outputs, lracz, tsol, amplz, psisol, racinepsi, epz, hur, wi, esz, humirac, humpotsol = initialize_soil_matrix(outputs.shape[0], soil, initial, outputs, manage)

    # Initialize values
    ndebsen = 0
    cut_number = 0

    #########################
    ### Weather variables ###
    #########################
    outputs = compute_weather_variables(outputs, station, gamma)

    ########################################
    ### Irrigation automatic calculation ###
    ########################################
    # When irrigation is automated calculated, we compute here the days where irrigation is 'active' (compute_airg = True) and should be computed.
    if manage.CODECALIRRIG == 1 and manage.CODEDATE_IRRIGAUTO == 1:
        outputs['compute_airg'] = dates_irrig(outputs, manage.DATEDEB_IRRIGATUO, manage.DATEFIN_IRRIGATUO)
    
    #############
    ### Day 0 ###
    #############
    outputs, humirac, hur, wi, esz, tsol, amplz, ndebsen = compute_outputs_day_zero(outputs, crop, soil, manage, constants, initial, station, lracz, hur, wi, esz, epz, tsol, humirac, humpotsol, psisol, racinepsi, amplz, ndebsen, hurlim, aevap, s, fco2s, durviei, gamma, msresiduel, fco2)

    ##################
    ### Daily loop ###
    ##################
    for i in range(1, len(outputs)-1):
        
        #########################################################
        ### Update variables with value from the previous day ###
        #########################################################
        outputs.loc[i,'densite'] = outputs.loc[i-1,'densite'].copy()
        outputs.loc[i,'zrac'] = outputs.loc[i-1,'zrac'].copy()
        outputs.loc[i,'znonli'] = outputs.loc[i-1,'znonli'].copy()
        outputs.loc[i,'nbfeuille'] = outputs.loc[i-1,'nbfeuille'].copy()
        outputs.loc[i,'dayLAIcreation'] = outputs.loc[i-1,'dayLAIcreation'].copy()
        hur[i] = hur[i-1].copy()
        humpotsol[i] = humpotsol[i-1].copy()
        tsol[i] = tsol[i-1].copy()
        outputs.loc[i,'rfvi'] = outputs.loc[i-1,'rfvi'].copy()
        outputs.loc[i,'vernalisation_ongoing'] = outputs.loc[i-1,'vernalisation_ongoing'].copy()
        

        ###################
        ### Development ###
        ###################
    
        ## 0. Sowing / planting
        if outputs.loc[i,'plt'] == 0:
            if (outputs.loc[i,'doy'] == manage.IPLT0):
                outputs.loc[i:,'plt'] = 1

        ## 1. Emergence / budding
        if outputs.loc[i,'plt'] == 1:
            if (crop.CODEPERENNE == 1) | (crop.HERBACEOUS):
                outputs, humirac = emergence_macro(i, outputs, crop, soil, manage, tsol, humirac, hur, humpotsol)
            elif crop.CODEPERENNE == 2:
                outputs.loc[i,'findorm'], outputs.loc[i,'cu'], outputs.loc[i-1,'lev'], outputs.loc[i,'thn'], outputs.loc[i,'gdh'], outputs.loc[i,'lev'] = budding(i, crop.CODEDORMANCE, crop.Q10, outputs['Temp_max'], outputs['Temp_min'], crop.JVC, outputs.loc[i-1,'lev'], crop.TDMINDEB, crop.TDMAXDEB, outputs.loc[i,'hourly_temp'], outputs.loc[i-1,'gdh'], crop.STDORDEBOUR)
        
        ## 2. Development temperature
        if outputs.loc[i,'plt'] == 1:
            outputs.loc[i,'udevcult'], outputs.loc[i,'somtemp'] = development_temperature(outputs.loc[i-1,'tcult'], outputs.loc[i,'temp'], crop.TDMAX, crop.TDMIN, crop.TCXSTOP, crop.CODERETFLO, crop.STRESSDEV, outputs.loc[i-1,'turfac'], crop.CODETEMP, outputs.loc[i-1,'somtemp'], outputs.loc[i-1,'drp'])

        ## 3. Effect of photoperiod
        if crop.CODEPHOT == 1:
            outputs.loc[i,'rfpi'] = photoperiod_effect(crop.HERBACEOUS, outputs.loc[i,'lev'], outputs.loc[i,'findorm'], outputs.loc[i-1,'drp'], crop.SENSIPHOT, outputs.loc[i,'phoi'], crop.PHOSAT, crop.PHOBASE)
        else:
            outputs.loc[i,'rfpi'] = 1

        ## 4. Effect of vernalisation
        if crop.CODEBFROID != 1:
            outputs.loc[i,'rfvi'], outputs.loc[i,'jvi'], outputs.loc[i,'vernalisation_ongoing'] = vernalisation_effect(crop.HERBACEOUS, crop.CODEBFROID, outputs.loc[i,'ger'], crop.TFROID, outputs.loc[i-1,'tcult'], crop.AMPFROID, outputs.loc[0:i,'jvi'], crop.JVCMINI, crop.JVC, outputs.loc[i,'findorm'], outputs.loc[i,'doy'], crop.JULVERNAL, crop.CODEPERENNE, outputs.loc[i,'rfvi'], outputs.loc[i,'vernalisation_ongoing'])
        else:
            outputs.loc[i,'rfvi'] = 1

        ## 5. Compute of phenological stage this day
        if crop.CODEINDETERMIN == 1:
            outputs.loc[i,'upvt_post_lev'], outputs.loc[i,'sum_upvt_post_lev'], outputs.loc[i,'amf'], outputs.loc[i,'lax'], outputs.loc[i,'flo'], outputs.loc[i,'drp'], outputs.loc[i,'debdes'], outputs.loc[i,'mat'], outputs.loc[i,'sen'], outputs.loc[i,'lan'], outputs.loc[i,'somcour'] = phenological_stage(outputs.loc[i,'lev'], outputs.loc[i,'udevcult'], outputs.loc[i,'rfpi'], outputs.loc[i,'rfvi'], outputs.loc[i-1,'sum_upvt_post_lev'], crop.STLEVAMF,
                                                    crop.STAMFLAX, crop.STLEVDRP, crop.STFLODRP, crop.STDRPDES, crop.CODEINDETERMIN, crop.STDRPMAT, crop.STDRPNOU, crop.CODLAINET, crop.STLAXSEN, crop.STSENLAN, outputs.loc[i-1,'lan'], outputs.loc[i-1,'somcour'])
        elif crop.CODEINDETERMIN == 2:  
            outputs.loc[i,'upvt_post_lev'], outputs.loc[i,'sum_upvt_post_lev'], outputs.loc[i,'amf'], outputs.loc[i,'lax'], outputs.loc[i,'flo'], outputs.loc[i,'drp'], outputs.loc[i,'debdes'], outputs.loc[i,'mat'], outputs.loc[i,'sen'], outputs.loc[i,'lan'], outputs.loc[i,'nou'], outputs.loc[i,'somcour'] = phenological_stage(outputs.loc[i,'lev'], outputs.loc[i,'udevcult'], outputs.loc[i,'rfpi'], outputs.loc[i,'rfvi'], outputs.loc[i-1,'sum_upvt_post_lev'], crop.STLEVAMF,
                                                    crop.STAMFLAX, crop.STLEVDRP, crop.STFLODRP, crop.STDRPDES, crop.CODEINDETERMIN, crop.STDRPMAT, crop.STDRPNOU, crop.CODLAINET, crop.STLAXSEN, crop.STSENLAN, outputs.loc[i-1,'lan'], outputs.loc[i-1,'somcour'])

        ###################
        ### Leaf growth ###
        ###################
        outputs.loc[i,'deltai'], outputs.loc[i,'deltai_dev'], outputs.loc[i,'deltai_dens'], outputs.loc[i,'deltai_t'], outputs.loc[i,'ulai'], outputs.loc[i,'deltai_stress'], outputs.loc[i,'efdensite'], outputs.loc[i,'vmax'], outputs.loc[i,'lai'], outputs.loc[i,'mafeuilverte'], outputs.loc[i,'dltaisen'], outputs.loc[i,'dltaisenat'], outputs.loc[i,'laisen'], outputs.loc[i,'lan'], outputs.loc[i,'sen'], outputs.loc[i,'ratiotf'], outputs.loc[i,'stopfeuille_stage'] = leaf_growth(i, outputs.loc[i,'lev'], outputs.loc[i,'lax'], outputs.loc[i,'sum_upvt_post_lev'], crop.STLEVAMF, crop.VLAIMAX, crop.STAMFLAX, crop.UDLAIMAX, crop.DLAIMAX, crop.PENTLAIMAX,
                    outputs.loc[i-1,'tcult'], crop.TCXSTOP, crop.TCMAX, crop.TCMIN, crop.ADENS, crop.BDENS, outputs.loc[i,'densite'], outputs.loc[i-1,'turfac'], outputs.loc[i,'phoi'], outputs.loc[i-1,'phoi'], outputs.loc[i,'ratiotf'],
                    crop.PHOBASE, outputs.loc[i,'rfpi'], crop.DLAIMIN, outputs.loc[i-1,'lai'], crop.LAICOMP, crop.STOPFEUILLE, outputs.loc[i-1,'vmax'], crop.CODLAINET, outputs.loc[i-1,'dltaisenat'], outputs.loc[i-1,'fstressgel'], outputs.loc[i-1,'laisen'], outputs.loc[i,'lax'], crop.SLAMIN, crop.SLAMAX, outputs.loc[i-1,'dltamsen'], outputs.loc[i-1,'dltaisen'], outputs.loc[i,'lan'], crop.CODEPHOT_PART, outputs.loc[i,'amf'], crop.TIGEFEUIL, outputs.loc[i-1,'dltams'], crop.CODEINDETERMIN, outputs.loc[i-1,'sla'], outputs['sen'].array, outputs.loc[i-1,'sen'], outputs.loc[i,'somcour'], crop.STSENLAN, outputs['lai'].array)

        #############################
        ### Intercepted radiation ###
        #############################
        outputs.loc[i,'raint'] = radiation_interception(constants.PARSURRG, outputs.loc[i,'trg'], crop.EXTIN, outputs.loc[i,'lai'])


        #############################
        ### Thermal stress on RUE ###
        #############################
        outputs.loc[i,'ftemp'] = thermal_stress_on_biomass_broduction(outputs.loc[i-1,'tcult'], crop.TEMIN, crop.TEMAX, crop.TEOPT, crop.TEOPTBIS)


        ##########################
        ### Biomass production ###
        ##########################
        outputs.loc[i,'dltams'], outputs.loc[i,'ebmax'], outputs.loc[i,'dltafv'], outputs.loc[i,'pfeuilverte'] = biomass_production(outputs.loc[i,'raint'], outputs.loc[i,'lev'], outputs.loc[i-1,'amf'], outputs.loc[i-1,'drp'], crop.EFCROIJUV, crop.EFCROIREPRO, crop.EFCROIVEG,
                                constants.COEFB, outputs.loc[i-1,'swfac'], fco2, outputs.loc[i,'ftemp'], outputs.loc[i, "deltai"], crop.SLAMAX)

        # Biomass remobilised on previous day
        # outputs.loc[i,'dltams'] = outputs.loc[i,'dltams'] + outputs.loc[i-1,'dltaremobil']

        outputs.loc[i,'mafeuil'] = outputs.loc[i,'mafeuilverte'] + outputs.loc[i-1,'mafeuiljaune']

        # Stem structural part biomass
        # = leaves biomass proportion
        outputs.loc[i,'matigestruc'] = crop.TIGEFEUIL * outputs.loc[i,'mafeuil'] # eq 7.7

        # Biomass remobilisation
        outputs.loc[i,'dltaremobil'], outputs.loc[i,'restemp'], outputs.loc[i,'dltams'], outputs.loc[i,'fpv'], outputs.loc[i,'sourcepuits'], outputs.loc[i,'dltarestemp'] = cumultative_partitioning(outputs.loc[i,'lev'], outputs.loc[i,'mafeuil'], outputs.loc[i,'matigestruc'], outputs.loc[i-1,'mafruit'], outputs.loc[i-1,'maenfruit'], crop.CODEPERENNE, outputs.loc[i,'masec'], initial.RESTEMP0, crop.RESPLMAX, outputs.loc[i,'densite'], crop.REMOBRES, outputs.loc[i,'deltai'], crop.SLAMIN, outputs.loc[i,'ratiotf'], crop.CODEINDETERMIN, outputs.loc[i-1,'cumdltaremobil'], outputs.loc[i,'dltams'])

        outputs.loc[i,'masec'] = outputs.loc[i-1,'masec'] + outputs.loc[i,'dltams']
        outputs.loc[i,'masecnp'] = outputs.loc[i,'matigestruc'] + outputs.loc[i-1,'masecnp'] + outputs.loc[i,'dltams'] - outputs.loc[i,'dltarestemp']

        ##################
        ### Senescence ###
        ##################
        outputs.loc[i,'fstressgel'] = senescence_stress(outputs.loc[i,'lev'], outputs.loc[i,'ulai'], crop.VLAIMAX, outputs.loc[i-1,'temp_min'], crop.TGELJUV10, crop.TGELJUV90, crop.TGELVEG10, crop.TGELVEG90, crop.TLETALE, crop.TDEBGEL, crop.CODGELJUV, crop.CODGELVEG)


        # Plant death because of frost
        if outputs.loc[i,'fstressgel'] == 0:
            outputs.loc[i,'dltaisen'], outputs.loc[i,'dltamsen'], outputs.loc[i,'laisen'], outputs.loc[i,'amf'], outputs.loc[i,'lax'], outputs.loc[i,'sen'], outputs.loc[i,'lan'], outputs.loc[i,'drp'], outputs.loc[i,'mat'], crop.STSENLAN, crop.STLAXSEN, crop.STLEVAMF, crop.STAMFLAX, crop.STLEVDRP, crop.STDRPMAT = lethal_frost(outputs.loc[i,'lai'], outputs.loc[i,'mafeuilverte'], outputs.loc[i-1,'laisen'], outputs.loc[i,'dltaisen'], outputs.loc[i,'stopfeuille_stage'], outputs.loc[i,'restemp'], crop.CODEPERENNE, outputs.loc[i,'amf'], outputs.loc[i,'lax'], outputs.loc[i,'sen'], outputs.loc[i,'lan'], outputs.loc[i,'drp'], outputs.loc[i,'mat'], outputs.loc[i,'sum_upvt_post_lev'], crop.STSENLAN, crop.STLAXSEN, crop.STLEVAMF, crop.STAMFLAX, crop.STLEVDRP, crop.STDRPMAT)
        else:
            outputs['dayLAIcreation'], outputs['durage'], outputs['senstress'], outputs['tdevelop'], outputs['durvie'], outputs.loc[i,'dltaisen'], outputs.loc[i,'somsenreste'], ndebsen, outputs.loc[i,'somtemp'], outputs.loc[i,'dltamsen'], outputs.loc[i,'deltamsresen'], outputs.loc[i,'msres'], outputs.loc[i,'msresjaune'], outputs.loc[i,'durvie_n'] = senescence(i, outputs['lev'], outputs.loc[i,'ulai'], outputs.loc[i,'somtemp'], crop.VLAIMAX, durviei, crop.DURVIEF, outputs.loc[i-1,'senfac'], outputs.loc[i,'udevcult'], outputs.loc[i,'fstressgel'],
                outputs['dayLAIcreation'].array, outputs['senstress'].array, outputs['tdevelop'].array, outputs['durvie'].array, outputs['durage'].array, outputs['deltai'].array, outputs.loc[i-1,'somsenreste'], initial.LAI0, ndebsen, outputs['dltafv'], crop.RATIOSEN, crop.CODEPLANTE, outputs.loc[i-1,'msres'], msresiduel, outputs.loc[i-1,'msresjaune'], crop.CODLAINET, outputs['pfeuilverte'], outputs['dltams'])
        

        ######################
        ### Fallen biomass ###
        ######################
        outputs.loc[i,'dltamstombe'], outputs.loc[i,'mafeuiltombe'], outputs.loc[i,'mafeuiljaune'] = fallen_biomass(crop.ABSCISSION, outputs.loc[i,'dltamsen'], outputs.loc[i-1,'mafeuiltombe'], outputs.loc[i-1,'mafeuiljaune'])

        ###################
        ### Root growth ###
        ###################
        outputs.loc[i,'zrac'], outputs.loc[i,'deltaz'], outputs.loc[i,'deltaz_t'], outputs.loc[i,'deltaz_stress'], outputs.loc[i,'efda'], outputs.loc[i,'znonli'] = root_growth(outputs.loc[i,'ger'], outputs.loc[i,'lax'], outputs.loc[i,'sen'], outputs.loc[i,'rec'], crop.STOPRAC, crop.CODEPERENNE, outputs.loc[i,'findorm'], manage.PROFSEM, outputs.loc[i,'zrac'], crop.CODETEMPRAC,
                outputs.loc[i-1,'tcult'], crop.TCMAX, crop.TCMIN, crop.TGMIN, crop.CROIRAC, hur[i], tsol[i], soil.HMIN, crop.SENSRSEC, soil.DEPTH, constants.DASEUILBAS, constants.DASEUILHAUT, crop.CONTRDAMAX, soil.DAF, outputs.loc[i,'lev'], soil.HCC, crop.HERBACEOUS, outputs.loc[i-1, 'znonli'])

        ####################
        ### Root density ###
        ####################
        lracz[i], outputs.loc[i,'cumlracz'], outputs.loc[i,'zdemi'], humirac[i], outputs.loc[i,'humirac_mean'] = root_density(lracz[i], outputs.loc[i,'zrac'], outputs.loc[i,'znonli'], soil.DEPTH, crop.ZPRLIM, crop.ZPENTE, outputs.loc[i,'ger'], crop.CODEPERENNE, constants.LVOPT, s, manage.PROFSEM, 
                                                                                               hur[i], soil.HMIN, humirac[i])
        
        ######################
        ## Water potential ###
        ######################
        humpotsol[i], psisol[i], racinepsi[i], outputs.loc[i,'psibase'] = plant_soil_water_potential(outputs.loc[i-1,'lev'], initial.ZRAC0, soil.DEPTH, constants.PSIHUCC, constants.PSIHUMIN, soil.DAF, soil.HMIN, soil.HCC, crop.POTGERMI, manage.PROFSEM, outputs.loc[i,'zrac'], lracz[i], hur[i], racinepsi[i])


        ############################################################
        ### Compute irrigation when it is automatically computed ###
        ############################################################
        if manage.CODECALIRRIG == 1:
            outputs.loc[i,'airg'] = automated_irrig(manage.CODEDATE_IRRIGAUTO, outputs.loc[i,'compute_airg'], manage.RATIOL, outputs.loc[i-1,'swfac'], manage.EFFIRR, manage.DOSEIRRIGMIN, manage.DOSIMX, soil.HCC, hur[i-1], outputs.loc[i-1,'zrac'], outputs.loc[i,'doy'], manage.IPLT0, constants.IRRLEV)

        ####################################
        ### Water intercepted by foliage ###
        ####################################
        if crop.CODEINTERCEPT == 1:
            outputs.loc[i,'stemflow'], outputs.loc[i,'mouill'] = foliage_water(manage.CODLOCIRRIG, crop.STEMFLOWMAX, crop.KSTEMFLOW, outputs.loc[i,'airg'], outputs.loc[i,'trr'], crop.MOUILLABIL, outputs.loc[i-1,'lai'])
        elif crop.CODEINTERCEPT == 2:
            outputs.loc[i,'mouill'] = 0

        #######################
        ## Soil evaporation ###
        #######################

        # Potential
        outputs.loc[i,'eos'] = potential_soil_evap_beer_law(outputs.loc[i,'etp'], crop.EXTIN, outputs.loc[i,'lai'])

        # Actual
        outputs.loc[i,'sumes0'] = outputs.loc[i-1,'sumes0'].copy()
        outputs.loc[i,'sumes1'] = outputs.loc[i-1,'sumes1'].copy()
        outputs.loc[i,'sumes2'] = outputs.loc[i-1,'sumes2'].copy()
        outputs.loc[i,'sesj0'] = outputs.loc[i-1,'sesj0'].copy()
        outputs.loc[i,'ses2j0'] = outputs.loc[i-1,'ses2j0'].copy()
        outputs.loc[i,'smes02'] = outputs.loc[i-1,'smes02'].copy()
        outputs.loc[i,'nstoc'] = outputs.loc[i-1,'nstoc'].copy()
        outputs.loc[i,'stoc'] = outputs.loc[i-1,'stoc'].copy()

        outputs.loc[i,'esol'], outputs.loc[i,'sumes0'], outputs.loc[i,'sumes1'], outputs.loc[i,'sumes2'], outputs.loc[i,'sesj0'], outputs.loc[i,'ses2j0'], outputs.loc[i,'smes02'], outputs.loc[i,'nstoc'], outputs.loc[i,'stoc'] = actual_soil_evaporation(outputs.loc[i,'trr'], outputs.loc[i,'airg'], outputs.loc[i,'eos'], aevap, soil.Q0, outputs.loc[i,'sumes0'], outputs.loc[i,'sumes1'], outputs.loc[i,'sumes2'], outputs.loc[i,'sesj0'], outputs.loc[i,'ses2j0'], outputs.loc[i,'smes02'], outputs.loc[i,'nstoc'], outputs.loc[i,'stoc']) 

        # Soil layers contribution
        outputs.loc[i,'esol'], esz[i], wi[i] = soil_layers_contribution_to_evaporation(outputs.loc[i,'esol'], soil.ZESX, hurlim, soil.CFES, soil.HCC, wi[i], hur[i], esz[i])

        ###########################
        ### Plant transpiration ###
        ###########################

        # Potential
        outputs.loc[i,'eo'], outputs.loc[i,'eop'], outputs.loc[i,'emd'], outputs.loc[i,'edirect'], outputs.loc[i,'directm'] = potential_transpiration_coef(outputs.loc[i,'lai'], outputs.loc[i,'etp'], crop.KMAX, crop.CODEINTERCEPT, outputs.loc[i,'eos'], outputs.loc[i,'esol'], crop.EXTIN, outputs.loc[i,'mouill'])
        
        #####################
        ### Water balance ###
        #####################
        hur[i], outputs.loc[i,'drain'] = water_balance(outputs.loc[i,'trr'], outputs.loc[i,'airg'], outputs.loc[i,'mouill'], soil.DEPTH, hur[i].copy(), epz[i-1], esz[i], soil.HCC, hurlim)

        ###########################
        ### Plant transpiration ###
        ###########################

        # Actual
        outputs.loc[i,'ep']  = actual_transpiration(outputs.loc[i-1,'swfac'], outputs.loc[i,'eop'])

        epz[i] =  soil_contribution_to_transpiration(outputs.loc[i,'ep'], epz[i], soil.DEPTH, hur[i], soil.HMIN, lracz[i], outputs.loc[i-1,'lev'], manage.PROFSEM, outputs.loc[i,'zrac'])

        outputs.loc[i,'et'] = outputs.loc[i,'esol'] + outputs.loc[i,'ep']


        ######################
        ### Water stresses ###
        ######################
        outputs.loc[i,'teta'], outputs.loc[i,'swfac'], outputs.loc[i,'tetstomate'], outputs.loc[i,'turfac'], outputs.loc[i,'teturg'], outputs.loc[i,'senfac'], outputs.loc[i,'tetsen'], outputs.loc[i,'resrac'] = water_stress(outputs.loc[i,'eop'], crop.SWFACMIN, crop.RAPSENTURG, outputs.loc[i,'cumlracz'], crop.PSITURG, crop.PSISTO, crop.RAYON, outputs.loc[i,'zrac'], outputs.loc[i-1,'lev'], soil.HMIN, hur[i], manage.PROFSEM, soil.DEPTH)

        ####################
        ### Plant height ###
        ####################
        outputs.loc[i,'hauteur'] = plant_height(crop.HAUTMAX, crop.HAUTBASE, outputs.loc[i,'lai'], outputs.loc[i,'laisen'])

        ####################
        ### Wind profile ###
        ####################
        outputs.loc[i,'dh'], outputs.loc[i,'z0'] = wind_profile(soil.Z0SOLNU, outputs.loc[i,'hauteur'])
        
        ##########################
        ### Specific leaf area ###
        ##########################
        outputs.loc[i,'tursla'] = (outputs.loc[i-1,'tursla'] + outputs.loc[i,'turfac']) / 2
        outputs.loc[i,'sla'] = max(outputs.loc[i,'tursla'] * crop.SLAMAX, crop.SLAMIN)

        ###################################################################
        ### Iterative calculation of crop temperature and net radiation ###
        ###################################################################
        outputs.loc[i,'rnet'], outputs.loc[i,'rglo'], outputs.loc[i,'albedolai'], outputs.loc[i,'albsol'], outputs.loc[i,'tcult'], outputs.loc[i,'tcultmax'], outputs.loc[i,'converge'] = iterative_calculation(outputs.loc[i,'temp'], outputs.loc[i-1,'lev'], outputs.loc[i,'temp_max'], outputs.loc[i,'temp_min'], outputs.loc[i,'et'], outputs.loc[i,'z0'], soil.ALBEDO, hur[i,0], soil.HMINF_1, soil.HCCF_1, station.ALBVEG, outputs.loc[i,'lai'], outputs.loc[i,'trg'], outputs.loc[i,'tpm'], outputs.loc[i,'fracinsol'], station.CODERNET,
            station.CODECALTEMP, outputs.loc[i,'raint'], constants.PARSURRG, outputs.loc[i,'ratm'], outputs.loc[0,'tcultmin'], outputs.loc[i,'wind'], outputs.loc[0,'tcultmax'], outputs.loc[i,'daylen'], station.ZR, outputs.loc[i-1,'lai'], soil.Z0SOLNU, outputs.loc[i,'hauteur'])

        ########################
        ### Soil temperature ###
        #######"################
        outputs.loc[i,'amplsurf'], amplz[i], tsol[i] = soil_temperature(outputs.loc[i,'temp_min'], outputs.loc[i,'tcultmax'], outputs.loc[i,'tcult'], soil.DEPTH, constants.DIFTHERM, tsol[i-1])


        ##############################
        ### Thermal stress indices ###
        ##############################
        outputs.loc[i,'ftempremp'] = thermal_stress_on_grain_filling(crop.CODETREMP, outputs.loc[i-1,'temp_min'], outputs.loc[i-1,'tcultmax'], crop.TMINREMP, crop.TMAXREMP)
        if outputs.loc[i,'flo'] == 1:
            outputs.loc[i,'fgelflo'] = frost_stress_on_fruit_number(outputs.loc[i-1,'temp_min'], crop.TGELFLO90, crop.TGELFLO10, crop.TLETALE, crop.TDEBGEL, crop.CODGELFLO)
        else:
            outputs.loc[i,'fgelflo'] = 1

        # For forage crops only
        if crop.CODEPLANTE == 'FOU':

            # Specific outputs
            outputs.loc[i,'msneojaune'] = outputs.loc[i,'mafeuiljaune'].copy()
            outputs.loc[i,'masec'] = outputs.loc[i,'masec'] - outputs.loc[i,'dltamstombe'] # should this also be done for annual plants ?
            outputs.loc[i,'masectot'] = outputs.loc[i-1,'masectot'].copy()
            outputs.loc[i,'masecneo'] = outputs.loc[i-1,'masecneo'] + outputs.loc[i,'dltams']

            #######################
            ### Yield formation ###
            #######################
            outputs.loc[i,'mafruit'], outputs.loc[i,'msrec_fou'] = yield_formation_forage(outputs.loc[i,'masec'], outputs.loc[i,'msresjaune'], outputs.loc[i,'msneojaune'], msresiduel)

            ###############
            ### Cutting ###
            ###############
            outputs.loc[i,'masectot'], outputs.loc[i,'msrec_fou'], outputs.loc[i,'mafruit'], outputs.loc[i,'masecneo'], outputs.loc[i,'msresjaune'], outputs.loc[i,'msneojaune'], outputs.loc[i,'mafeuiljaune'], outputs.loc[i,'masec'], outputs.loc[i,'msres'], outputs.loc[i,'lai'], cut_number, outputs.loc[i,'lev'], outputs.loc[i,'amf'], outputs.loc[i,'lax'], outputs.loc[i,'flo'], outputs.loc[i,'drp'], outputs.loc[i,'debdes'], outputs.loc[i,'dayLAIcreation'], outputs.loc[i,'somsenreste'], outputs.loc[i,'dltaisen'], outputs.loc[i,'laisen'], outputs.loc[i,'sum_upvt_post_lev'] = cutting(manage.CODEFAUCHE, cut_number, manage.JULFAUCHE, outputs.loc[i,'doy'], manage.CODEMODFAUCHE, msresiduel, outputs.loc[i,'masec'], outputs.loc[i,'masecneo'], lairesiduel, outputs.loc[i,'lai'], outputs.loc[i,'msrec_fou'], outputs.loc[i,'mafruit'], manage.MSCOUPEMINI, outputs.loc[i,'masectot'], outputs.loc[i,'msresjaune'], outputs.loc[i,'msneojaune'], outputs.loc[i,'mafeuiljaune'], outputs.loc[i,'msres'], initial.STADE0, outputs.loc[i,'lev'], outputs.loc[i,'amf'], outputs.loc[i,'lax'], outputs.loc[i,'flo'], outputs.loc[i,'drp'], outputs.loc[i,'debdes'], i, outputs.loc[i,'dayLAIcreation'], outputs.loc[i,'deltai'], outputs.loc[i,'somsenreste'], outputs.loc[i,'dltaisen'], outputs.loc[i,'laisen'], outputs.loc[i,'sum_upvt_post_lev'])

        #########################
        ### End of daily loop ###
        #########################
    
    ###################################################
    ### Dates of phenological stages and BCCH codes ###
    ###################################################
    if crop.CODEINDETERMIN == 2:
        outputs['bbch'], ind_drp, ind_lev, ind_amf, ind_debdes, ind_mat, ind_nou = phenological_stage_dates(outputs['lev'], outputs['amf'], outputs['debdes'], outputs['drp'], outputs['nou'], outputs['flo'], outputs['findorm'], outputs['mat'], outputs['lax'],
                                crop.CODEINDETERMIN, crop.CODEPERENNE)
    elif crop.CODEINDETERMIN == 1:
        outputs['bbch'], ind_drp, ind_lev, ind_amf, ind_mat, ind_debdes = phenological_stage_dates(outputs['lev'], outputs['amf'], outputs['debdes'], outputs['drp'], outputs['nou'], outputs['flo'], outputs['findorm'], outputs['mat'], outputs['lax'],
                                crop.CODEINDETERMIN, crop.CODEPERENNE)
    

    ##################################
    ### Radiative stress indicator ###
    ##################################
    outputs['radiative_stress'] = radiative_stress(outputs[['raint','ebmax']].copy(), constants.COEFB)


    #######################
    ### Yield formation ###
    #######################
    if (crop.CODEINDETERMIN == 1) & (crop.CODEPLANTE != 'FOU'):

        #####################
        ### Grains number ###
        #####################
        outputs['vitmoy'], outputs.loc[ind_drp,'nbgrains'] = harvested_organs_number(outputs['dltams'].array, ind_drp, crop.NBJGRAIN, crop.CGRAINV0, crop.CGRAIN, crop.NBGRMAX, crop.NBGRMIN)

        ############################
        ### Carbon Harvest Index ###
        ############################
        outputs['ircarb'] = carbon_harvest_index(outputs['ircarb'].array, ind_mat, ind_drp, crop.CODEIR, crop.VITIRCARB, crop.VITIRCARBT, crop.IRMAX, outputs['sum_upvt_post_lev'])
    
        for i in range(ind_drp+1, len(outputs)-1):
            ######################################
            ### Grains number reduced by frost ###
            ######################################
            outputs.loc[i,'nbgrains'], outputs['nbgraingel'], outputs.loc[i,'pgraingel'] = frost_reduction_harvested_organs(i, ind_drp, outputs.loc[i-1,'nbgrains'], outputs.loc[i,'fgelflo'], outputs['pgrain'], outputs['nbgraingel'].array)

            ###################
            ### Grains mass ###
            ###################
            outputs.loc[i,'mafruit'], outputs['deltags'], outputs.loc[i,'pgrain'] = harvested_organs_mass(i, outputs.loc[i,'ircarb'], outputs.loc[i,'masec'], outputs.loc[i-1,'ircarb'], outputs.loc[i-1,'masec'], outputs.loc[i,'ftempremp'], crop.PGRAINMAXI, outputs.loc[i,'nbgrains'], ind_mat, outputs.loc[i,'pgraingel'], outputs['deltags'].array, ind_drp, outputs.loc[i-1,'mafruit'], outputs.loc[i-1,'pgrain'], outputs.loc[i,'nbgraingel'])

            # Fruits/grains envelope biomass
            outputs.loc[i,'maenfruit'] = outputs.loc[i,'nbgrains'] * crop.ENVFRUIT * crop.PGRAINMAXI * 0.01

        ############################
        ### Grains water content ###
        ############################
        if ind_drp != 0:
            outputs['teaugrain'] = grains_water_content(outputs['tcult'], outputs['temp'], ind_debdes, crop.H2OFRVERT, crop.DESHYDBASE, crop.TEMPDESHYD)


    ########################
    ### Harvest criteria ###
    ########################
    outputs['maturity_criteria'], outputs['water_content_min_criteria'], outputs['water_content_max_criteria'] = harvest_criteria(outputs['teaugrain'], outputs['mat'], ind_debdes, manage.CODRECOLTE, manage.CODEAUMIN, manage.H2OGRAINMIN, manage.H2OGRAINMAX)
    
    if crop.CODEPLANTE != 'FOU':
        outputs['rec'] = (outputs.water_content_max_criteria | outputs.maturity_criteria | outputs.water_content_min_criteria).astype(int)
        if outputs.rec.max() == 1: # TODO : why harvest not reached ?
            outputs.loc[outputs['rec'] == 1,'mafruit'] = 0
            outputs['mafruit_rec'] = outputs['mafruit'][outputs.rec == 0].values[-1]
        else:
            outputs['mafruit_rec'] = outputs['mafruit'].values[ind_mat]

    # Extra indicators calculation
    outputs['water_stress_day'] = (outputs.swfac < 1).astype(int)
    outputs['water_stress_day_value'] = (outputs.swfac < 1) * outputs.swfac
    outputs.loc[outputs.swfac == 1,'water_stress_day_value'] = np.nan

    outputs['thermal_stress_day'] = (outputs.ftemp < 1).astype(int)
    outputs['thermal_stress_day_value'] = (outputs.ftemp < 1) * outputs.ftemp
    outputs.loc[outputs.ftemp == 1,'thermal_stress_day_value'] = np.nan

    outputs['hur_0_10_cm'] = pd.DataFrame(hur[:,0:10]).mean(axis=1)
    if hur.shape[1] >= 20:
        outputs['hur_10_20_cm'] = pd.DataFrame(hur[:,10:20]).mean(axis=1)
    if hur.shape[1] >= 30:
        outputs['hur_20_30_cm'] = pd.DataFrame(hur[:,20:30]).mean(axis=1)
    if hur.shape[1] >= 40:
        outputs['hur_30_40_cm'] = pd.DataFrame(hur[:,30:40]).mean(axis=1)
    if hur.shape[1] >= 50:
        outputs['hur_40_50_cm'] = pd.DataFrame(hur[:,40:50]).mean(axis=1)
    if hur.shape[1] >= 60:
        outputs['hur_50_60_cm'] = pd.DataFrame(hur[:,50:60]).mean(axis=1)

    outputs['water_use_efficiency'] = [0 if et==0 else dltams*1000/et for et,dltams in zip(outputs.et.shift(1).fillna(0), outputs.dltams)]

    #########################
    ### END OF SIMULATION ###
    #########################

    return outputs, [lracz, tsol, amplz, psisol, racinepsi, epz, hur, wi, esz, humirac, humpotsol]





