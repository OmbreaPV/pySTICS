import pandas as pd
import numpy as np
from datetime import datetime



def initialize_outputs_df(weather, crop, manage, initial):
    '''
    This function initializes the outputs dataframe (all computed variables each day of simulation).
    It creates a pandas dataframe based on the weather file and initializes zeros columns for each computed variable.
    '''

    # All computed variables initialized to zero
    colonnes_nulles = [
        'deltaz','zrac','ulai','deltai_dev','deltai_stress','deltai',
        'lai','laisen','dayLAIcreation','eos','flagrain','seos','flagphase','ses','es_a','esol_sum','zdemi',
        'cumlracz','ep_sum','dltams','masec',
        'eo','eop','teta','ep','swfac',
        "maxlai_sofar",
        "nbgrains",
        "nbgraingel",
        "pgrain",
        "pgraingel",
        "ircarb",
        "deltags",
        "mafruit",
        "raint",
        "mafeuilverte",
        "mafeuiljaune",
        "mafeuil",
        "matigestruc",
        "masecnp",
        "maenfruit",
        "msrac",
        "restemp",
        "dltamsen",
        "durage",
        "durvie",'durvie_n',
        "somsen",
        "somtemp",
        "dltaisen",
        'dltaisenat',
        "resperenne",
        "resperennemax",
        "resperennestruc",
        "maperenne",
        "tetsen",
        "senfac",
        "turfac",
        "fstressgel",
        'day_lai_creation',
        "dltafv",
        "dltat",
        "dltares",
        "dltaremobsen",
        "masecveg",
        "restempmax",
        "dltarestemp",
        "senstress",
        'sen',
        'tdevelop',
        "fpv",
        "sourcepuits1",
        "sourcepuit2",
        "sourcepuits",
        "remob",
        "dltaremobilbrut",
        "dltaremobil",
        "dltaCO2resperenne",
        "remobilj",
        "fpft",
        "devjour",
        "spfruit",
        "nbinflo_recall",
        "nfruitnou",
        "allocfruit",
        "albsol",
        "albedolai",
        "rglo",
        "rnet",
        "et",
        "hauteur",
        "z0",
        "tcultmax",
        "tcult",
        "lev",
        "flo",
        "nou",
        "drp",
        "udevcult",
        "udevcult_lev",
        "rfpi",
        "jvi",
        "sum_dcdays",
        "plt",
        "ger",
        "jvi_cum",
        "rfvi",
        "upvt_ger_flo",
        "sum_upvt_ger_flo",
        "lax",
        "amf",
        "mat",
        "upvt_post_lev",
        "sum_upvt_post_lev",
        "sum_upvt",
        "ebmax",
        "tetstomate",
        "teturg",
        "gdh",
        "findorm",
        "cu",
        "deltaidens",
        "deltai_t",
        "esol",
        "stemflow",
        "mouill",
        "fapar",
        "rnetp1",
        "rnets",
        "dh",
        "ras0",
        "raa0",
        "rasinf",
        "raainf",
        "ras",
        "raa",
        "deltat",
        "L",
        "dsat",
        "rac",
        "rc",
        "ept",
        "dos",
        "emd",
        "rnetp2",
        "edirectm",
        "edirect",
        'dltamstombe','mafeuiltombe', 'nbjgrauto','amplsurf',
        'terme1','terme2','terme3',
        'psibase', 'moist', 'deltaz_t',
        'teaugrain', 'debdes',
        'drained_water',
        'deltaz_stress',
        'somger', 'densite','elong',
        'nbfeuille', 'somfeuille','nbjhumec','let','fgellev',
        'fgelflo','ftemp','ftempremp',
        'water_content_max_criteria', 'water_content_min_criteria', 'maturity_criteria',
        'coeflev',
        'efdensite', 'fgelveg','fgeljuv',
        'compute_airg','somcourdrp',
        'somsenreste',
        'efda',
        'humirac_mean',
        'znonli',
        'vmax',
        'msrec_fou', 'masectot', 'masecneo', 'msresjaune', 'msneojaune','deltamsresen', 'msres',
        'hur_0_10_cm','hur_10_20_cm','hur_20_30_cm','hur_30_40_cm','hur_40_50_cm','hur_50_60_cm',
        'water_stress_day','water_stress_day_value','thermal_stress_day','thermal_stress_day_value',
        'sumes0','sumes1','sumes2','ses2j0','sesj0','smes02','nstoc','stoc',
        'resrac'
        
    ]

    # Concatenate weather dataframe and zeros columns for all computed variables
    df2 = pd.DataFrame(0., index=weather.index, columns=colonnes_nulles)
    df = pd.concat([weather, df2], axis=1)

    # When irrigation is automatically computed, irrigation column is initialized to zero.
    if manage.CODECALIRRIG == 1:
        df['airg'] = 0

    # Initialize bbch variable at nan because 0 is emergence
    df["bbch"] = np.array([np.nan for i in range(df.shape[0])])
    df['bbch'] = df['bbch'].astype('float')

    # Initialize stem.leaf ratio to initial value
    df["ratiotf"] = crop.TIGEFEUIL

    # For indeterminate growth variables, initialize variables depending on plant parameters
    for K in range(crop.NBOITE):
        df[f"fpft{K}"] = 0
        df[f"potcroifruit{K}"] = 0
        df[f"nfruit{K}"] = 0
        df[f"croifruit{K}"] = 0

    # Initialize iterative calculation (rnet, tcult) convergence
    df["converge"] = False

    # Initialize variables with input initial values. 
    df.loc[0,'lai'] = initial.LAI0
    df.loc[0,'zrac'] = initial.ZRAC0
    df.loc[0,'znonli'] = initial.ZRAC0
    df.loc[0,'masec'] = initial.MASEC0
    df.loc[0,'msres'] = initial.MASEC0
    
    if initial.STADE0 == 'lev':
        if crop.HERBACEOUS:
            df['moist'] = 1
            df['ger'] = 1
        else:
            df['findorm'] = 1
        df['plt'] = 1
        df['lev'] = 1
    elif initial.STADE0 == 'amf':
        if crop.HERBACEOUS:
            df['plt'] = 1
            df['moist'] = 1
            df['ger'] = 1
        else:
            df['findorm'] = 1
        df['plt'] = 1
        df['lev'] = 1
        df['amf'] = 1
        df.loc[0,'sum_upvt_post_lev'] = crop.STLEVAMF

    return df


def initialize_soil_matrix(nb_day, soil, initial, outputs, manage):
    '''
    This function initializes soil matrixes (lines = days, columns = soil depths) with nan values or input initial values.
    '''

    # Root length density profile
    LRACZ_MATRIX = np.empty((nb_day, soil.DEPTH))
    LRACZ_MATRIX[:] = np.nan

    LRACZ_MATRIX[0,0:soil.EPC_1] = initial.DENSINITIAL_1
    if soil.EPC_2 != 0:
        LRACZ_MATRIX[0,soil.EPC_1:soil.EPC_1+soil.EPC_2] = initial.DENSINITIAL_2
    if soil.EPC_3 != 0:
        LRACZ_MATRIX[0,soil.EPC_1+soil.EPC_2:soil.EPC_1+soil.EPC_2+soil.EPC_3] = initial.DENSINITIAL_3
    if soil.EPC_4 != 0:
        LRACZ_MATRIX[0,soil.EPC_1+soil.EPC_2+soil.EPC_3:soil.EPC_1+soil.EPC_2+soil.EPC_3+soil.EPC_4] = initial.DENSINITIAL_4
    if soil.EPC_5 != 0:
        LRACZ_MATRIX[0,soil.EPC_1+soil.EPC_2+soil.EPC_3+soil.EPC_4:soil.EPC_1+soil.EPC_2+soil.EPC_3+soil.EPC_4+soil.EPC_5] = initial.DENSINITIAL_5

    LRACZ_MATRIX[0,int(initial.ZRAC0):] = 0
    LRACZ_MATRIX[0,:int(manage.PROFSEM)] = 0

    # Thermal amplitude profile
    AMPLZ_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    AMPLZ_MATRIX[:] = np.nan

    # Soil temperature profile
    TSOL_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    TSOL_MATRIX[:] = np.nan

    # Soil potential profile
    PSISOL_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    PSISOL_MATRIX[:] = np.nan

    #  Efficient root length densityparticipating in predawn potential profile, located in the moist layers
    RACINEPSI_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    RACINEPSI_MATRIX[:] = np.nan

    # Root water extraction from transpiration profile
    EPZ_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    EPZ_MATRIX[:] = np.nan

    # Microporosity elementary layer soil water content profile
    HUR_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    HUR_MATRIX[:] = np.nan
    HUR_MATRIX[0,0:soil.EPC_1] = initial.HINITF_1 * soil.DAF_1
    if soil.EPC_2 != 0:
        HUR_MATRIX[0,soil.EPC_1:soil.EPC_1+soil.EPC_2] = initial.HINITF_2 * soil.DAF_2
    if soil.EPC_3 != 0:
        HUR_MATRIX[0,soil.EPC_1+soil.EPC_2:soil.EPC_1+soil.EPC_2+soil.EPC_3] = initial.HINITF_3 * soil.DAF_3 
    if soil.EPC_4 != 0:
        HUR_MATRIX[0,soil.EPC_1+soil.EPC_2+soil.EPC_3:soil.EPC_1+soil.EPC_2+soil.EPC_3+soil.EPC_4] = initial.HINITF_4 * soil.DAF_4
    if soil.EPC_5 != 0:
        HUR_MATRIX[0,soil.EPC_1+soil.EPC_2+soil.EPC_3+soil.EPC_4:soil.EPC_1+soil.EPC_2+soil.EPC_3+soil.EPC_4+soil.EPC_5] = initial.HINITF_5 * soil.DAF_5

    # Soil water availability index profile
    WI_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    WI_MATRIX[:] = np.nan

    # Root water extraction from evaporation profile
    ESZ_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    ESZ_MATRIX[:] = np.nan

    # Water stress on root growth and density profile
    HUMIRAC_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    HUMIRAC_MATRIX[:] = np.nan

    # Potential to water content conversion profile
    HUMPOTSOL_MATRIX = np.empty(
        (nb_day, soil.DEPTH)
    )
    HUMPOTSOL_MATRIX[:] = np.nan

    return outputs, LRACZ_MATRIX, AMPLZ_MATRIX, TSOL_MATRIX, PSISOL_MATRIX, RACINEPSI_MATRIX, EPZ_MATRIX, HUR_MATRIX, WI_MATRIX, ESZ_MATRIX, HUMIRAC_MATRIX, HUMPOTSOL_MATRIX


def compute_constants(crop, soil, station, co2, manage):
    '''
    This module computes constants from crop, soil, station and technical parameters.
    '''

    # Soil related constants
    ha = soil.ARGI / 100 / 15 * soil.DAF_1
    hurlim = 10 * ha
    aevap = 0.5 * station.ACLIM * ((0.63 - ha) ** (5 / 3)) * (soil.HCCF_1/10 - ha) 

    # Crop related constants
    s = -np.log(1e-2) / (crop.ZLABOUR - crop.ZPENTE)
    fco2 = 2 - np.exp(np.log(2 - crop.ALPHACO2) * (co2 - 350) / (600 - 350))
    fco2s = 1 / (1 + 0.77 * (1 - fco2 / 2.5) * (1 - co2 / 330))
    durviei = (crop.RATIODURVIEI * crop.DURVIEF)

    # Station related constants
    patm = 1000 # atmospheric pressure
    gamma = 0.65 * patm / 1000

    # Forage crops constants
    msresiduel = crop.COEFMSHAUT * (manage.HAUTCOUPE - crop.HAUTBASE)
    lairesiduel = -1 / crop.KHAUT * np.log(1 - (manage.HAUTCOUPE - crop.HAUTBASE) / crop.HAUTMAX)

    return hurlim, aevap, s, fco2, fco2s, durviei, gamma, msresiduel, lairesiduel