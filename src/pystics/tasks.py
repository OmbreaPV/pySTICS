from pystics.stics_daily import *
from pystics.params import StationParams, CropParams, ManageParams, SoilParams, Constants, UserOptions

def run_stics_simulation(crop : CropParams, soil :SoilParams, constants : Constants, meteo_raw : pd.DataFrame, user : UserOptions, manage : ManageParams, station : StationParams) -> pd.DataFrame :
    """
    Run a simulation of crop growth and yield using the STICS (Simulateur mulTIdisciplinaire pour les Cultures Standard) model.

    Parameters
    ----------
    crop : CropParams
        Object containing the parameters for the crop being simulated.
    soil : SoilParams
        Object containing the parameters for the soil being used in the simulation.
    meteo_raw : pd.DataFrame
        A Pandas DataFrame containing the daily raw meteorological data for the simulation
    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame containing the simulated crop growth and yield data
    """


    # Construction du fichier météo
    if crop.CODEPERENNE == 1: # plante annuelle = le 1er jour simulation est la date de semis.
        ministics, meteo_gamma = meteo_df(meteo_raw, station, starting_doy= manage.IPLT0)
    elif crop.CODEPERENNE == 2: # plante pérenne = début autonomne si dat dorm breka calculée, sinon à date dorm break
        if user.CODEDORMANCE == 1: # date dorm break fixée = date de début de simulation
            ministics, meteo_gamma = meteo_df(meteo_raw, station, starting_doy= user.IFINDORM)
        elif user.CODEDORMANCE == 3: # date dorm break calculée voir 3.3.4.2
            ministics, meteo_gamma = meteo_df(meteo_raw, station, starting_doy=250) # dans 3.3.4.2 = dit de prendre automne ou été, peu d'effet sur résultat final.
    
    # Initialisation des df contenant les variables calculées à chaque pas de temps
    ministics = initialise_vars(ministics, soil, crop, station)
    LRAC_MATRIX = initialize_lrac_matrix(ministics, soil, crop)
    if crop.CODEINDETERMIN == 2: # plantes indéterminées
        alpha, beta, dfr = param_croiss_indeter(crop)

# Problème : il faut LAI /phéno/ LRAC matrix pour calculer tcult/transpiration !! et tcult pr calculer ces 3
# Solution choisie ATTENTION est ce que c'est grave ce décalage :
#   - pour calculer la transpi on prend phase phéno j-1, et LRAC_matrix j-1
#   - pour calculer tcult on prend LAI cumulé jusqu'au j-1.
#   - pour calculer rnet, je prends le LAI cumulé jusqu'au j-1
#   - pour calculer EOS (ds evap sol), je prends le lai cumulé jusqu'au j-1
#   - pr le calcul de cumulracZ/TETA1/TETA2/TETA_moy/ep1th/ep2th, je prends LRACMATRIX de j-1
#   - pour calculer TETSTOMATE, je prends ZRAC de j-1
#   - dans crop_temp et net_rad, condiiton sur DOY_EMERG > 1 se fait à j-1
#   - dans rnet, pr calculer albsol, je mets HUR1 j-1
#   - LAI i-1 dans radiation_intercep
# ces approx me semblent pas horribles A VERIFIER PAC
# je préfère calculer tcult d'abord car il est hyper important pr la phéno et l'augmentation du LAI
# j'ai ajouté des if i ==0 pour gérer les premiers jours où le 'j-1' ne sera pas accessible.

# Attention je vais jusqu'à len(mini)-1 pr être sur que ça pose pas de prob avec les j+1. Voir quand tout marchera bien comment on peut calculer sur le dernier jour.

    # Boucle 1 : calcul tcult / phénologie / croiss racinaire / densité racinaire / LAI
    for i in range(len(ministics)-1):

        ministics = radiation_intercep(ministics, crop, constants, i)

        ministics = hauteur(ministics, i, crop, soil)

        ministics = iterative_tcult_rnet(ministics, i, crop, soil, constants, station, LRAC_MATRIX, meteo_gamma, user) # calcul itératif de radiation nette - etp - evap sol - transpi plante - temp culture

        # que vaut tcult quand y'a pas de pousse qui reçoit rayonnement ? vérifier qu'on applique pas la même formule (ou ça dép lai).
        ministics = phenology(ministics, i, crop, user)

        ministics = root_growth(ministics, crop, soil, i)
        LRAC_MATRIX = root_density(ministics, crop, soil, constants, LRAC_MATRIX, i)

        ministics = leaf_area_index(ministics, crop, i)


   # Calcul de dates phénologiques utiles
    if crop.CODEINDETERMIN == 2:
        ministics, ind_Z69, ind_EMERG, ind_AMF, ind_NOU = phenology_dates(ministics, crop)
    elif crop.CODEINDETERMIN == 1:
        ministics, ind_Z69, ind_EMERG, ind_AMF = phenology_dates(ministics, crop)
    
    # Boucle 2 : calcul des stress impactant rendement et RUE / radiation intercept / prod biomass / partitionem biomass / formation rendement
    for i in range(len(ministics)-1):
        ministics = temp_stress(ministics, crop, i, ind_Z69)
        ministics = total_biomass_growth(ministics, crop, constants, i, ind_Z69)

        if crop.CODEINDETERMIN == 1: # plantes déterminées
            # ministics = biomass_partitioning(ministics, crop, constants, i, ind_EMERG)
            ministics = harvested_organs_number_determinate(ministics, i, crop, ind_Z69)
            ministics = yield_formation_determinate(ministics, crop, i, ind_Z69)
        elif crop.CODEINDETERMIN == 2: # plantes indéterminées
            ministics = force_puit_fruits(ministics, i, ind_Z69, crop, alpha, beta, dfr)
            ministics = biomass_partitioning(ministics, crop, constants, i, ind_EMERG)
            ministics = yield_formation_indeterminate(ministics, crop, i, ind_AMF, ind_Z69, ind_NOU, alpha, beta,dfr)
        ministics = water_drainage(ministics, soil, i)


    return ministics, LRAC_MATRIX

