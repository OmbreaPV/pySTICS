import pandas as pd
from datetime import datetime
import numpy as np
from scipy.optimize import fsolve
from math import log
import ast
from pyeto import fao
import math

from pystics.utils import isNaN, isinf

def frost_stress(t, tgel90, tgel10, tletale, tdebgel):
    '''
    Calcul de l'indice de stress de gel
    Voir Fig 4.16, sorte de sigmoide qui dépend de la Tmin du jour
    Idée : pour reconstruire la courbe Figure 4.16, on calcule la pente des droites constituant la courbe
    '''
    if (t >= tgel10) and (tgel10 < tdebgel):
        a = (0.9 - 1.0) / (tgel10 - tdebgel)
        b = 1.0 - (a * tdebgel)
        stress = (a * t) + b
    elif (t < tgel10) and (t >= tgel90):
        a = (0.9 - 0.1) / (tgel10 - tgel90)
        b = 0.9 - (a * tgel10)
        stress = (a * t) + b
    elif (t < tgel90):
        a = (0. - 0.1) / (tletale - tgel90)
        b = 0. - (a * tletale)
        stress = (a * t) + b
    stress = np.clip(stress, 0, 1)
    return stress


def thermal_stress_rue(t, temin, teopt, teoptbis):
    ''' 
    combined responses of photosynthesis and respiration to temperature
    Calcul de l'indice de stress thermique impactant la RUE (ftemp). Eq 4.25
    Si entre les 2 opt, stress vaut 1, sinon il diminue jusqu'à 0 au fur et à mesure que la Temp s'éloigne de cet intervalle
    '''
    if t<=teopt:
        stress = 1 - ((t - teopt) / (temin - teopt))**2
    elif t>=teoptbis:
        stress = 1 - ((t - teoptbis) / (temin - teoptbis))**2
    else:
        stress = 1
    stress = np.clip(stress, 0, 1)
    return stress

def thermal_stress_fruitfilling(t1, t2, tmin, tmax):
    '''
    Calcul de l'indice de stress thermique impactant le remplissage des grains (ftempremp). Fig 4.17
    Attention il faut utiliser Temp_min et Temp_max et pas Temp_moy
    Si entre min et max = stress vaut 1, sinon 0
    '''
    if (t1>tmin) and (t2<tmax):
        stress = 1
    else:
        stress = 0
    return stress


def shoot_biomass_production(radiation, ebmax, COEFB, ftemp, SWFAC, FCO2):
    """
    Calcul de la biomasse aérienne produite 
    Ministics : avait juste multiplié raint par RUE fixée, ici j'ai repris Eq 4.18. sans stress azoté (inns) et sans facteur en cas d'engorgement eau (waterlogging exobiom)
    Prochaine étape : dltaremobil
    """
    RUE = (ebmax - COEFB * radiation)
    DELTAMS = radiation * RUE * ftemp * SWFAC * FCO2
    return DELTAMS

def etp_stics(ministics, i, meteo_gamma, station):
    '''
    Calcul de l'évapotranspiration potentielle (ETP) avec l'équation de Penman (l'une des 3 options dans STICS : CODEETP)
    
    Remarques :
    J'ai remplacé trg par rnet
    Quand codeetp = 1 --> ETP fourni par météo (utilisé pr exemples de STICS)
    '''
    if station.CODEETP != 1:
        ministics.loc[i, "ETP"] = (
            0.408 * ministics.loc[i, "Delta"] * ministics.loc[i, "rnet"]
            + meteo_gamma
            * 900
            / (ministics.loc[i, "Temp"] + 273)
            * ministics.loc[i, "Wind"]
            * (ministics.loc[i, "esat"] - ministics.loc[i, "TPM"])
        ) / (
            ministics.loc[i, "Delta"] + meteo_gamma * (1 + 0.34 * ministics.loc[i, "Wind"])
        )

    return ministics


def iterative_tcult_rnet(ministics, i, CropParams, SoilParams, Constants, station, LRAC_MATRIX, meteo_gamma, user):  
    '''
    Réalistion du processus itératif (partie 9.2.2.3) pour calculer l'évapotranspiration, la température de culture et la radiation nette
    Explication : ces 3 calculs sont imbriqués et nécessitent de partir d'une approximation --> tcult = temp moy.
    '''
    # Processus itératif de calcul de rnet et tcult
    ministics.loc[i, "tcult"] = ministics.loc[i, "Temp"]
    j = 0
    while j < 5:
        ministics = net_radiation(ministics, i, station, SoilParams)

        ministics = evapotranspiration(
            ministics,
            i,
            CropParams,
            SoilParams,
            Constants,
            station,
            meteo_gamma,
            user,
            LRAC_MATRIX,
        )

        ministics = crop_temperature(
            ministics, i, CropParams
        )  # attention ça doit calculer tcult pour i+1 dans cette fonction

        if isNaN(ministics.loc[i + 1, "tcult"]):
            raise ValueError("tcult nan")

        if ministics.loc[i + 1, "tcult"] - ministics.loc[i, "tcult"] < 0.5:
            ministics.loc[i, "tcult"] = ministics.loc[
                i + 1, "tcult"
            ]  # à vérifier compréhension stics, mais c'est en i+1 qu'on calcule
            ministics.loc[
                i, "converge"
            ] = True  # impro juste pr les tests pr voir si ça converge ou pas chaque jour
            break
        ministics.loc[i, "tcult"] = ministics.loc[i + 1, "tcult"]
        j += 1

    # Après avoir obtenu la tcult finale, on recalcule les net radiation / etp / evap / transpi
    # ministics = net_radiation(ministics, i, Constants, SoilParams)
    # ministics = etp_stics(ministics, i, meteo_gamma)
    # ministics = soil_evaporation(ministics, CropParams, SoilParams, i, user)
    # ministics = plant_transpiration(ministics, CropParams, SoilParams, LRAC_MATRIX, i)

    return ministics


def evapotranspiration(
    ministics, i, CropParams, SoilParams, Constants, station, meteo_gamma, user, LRAC_MATRIX
):
    '''
    Calcul de l'évapotranspiration potentielle puis réelle.
    
    Explications :
    L'ETP est calculé à partir de la radiation nette.
    Deux options sont ensuite possibles pour calculer l'évapotranspiration réelle :
        - Beer's Law pour évap sol + Coeff culture pour transpiration
        - Bilan d'énergie pour répartir l'énergie nette entre les différents flux
    '''
    ministics = etp_stics(ministics, i, meteo_gamma, station)

    if CropParams.CODEINTERCEPT == 1:
        ministics = foliage_water(ministics, CropParams, i)
    elif CropParams.CODEINTERCEPT == 2:
        ministics.loc[i, "mouill"] = 0

    # Option Beer law + Coeff culture
    if CropParams.CODEBESO == 1:
        
        # D'abord soil et ça sert dans calcul EOP (potential_transpi_crop_coef eq 11.12)
        ministics = potential_soil_evap_beer_law(ministics, CropParams, i)
        ministics = actual_soil_evap(ministics, i, SoilParams, user, station)

        # Calcul transpi avec coeff culture
        ministics = potential_transpi_crop_coef(ministics, i, CropParams)  # faut avoir EOS avant
        ministics = actual_transpi(ministics, i, CropParams, SoilParams, LRAC_MATRIX)

    # Option Bilan d'énergie
    elif CropParams.CODEBESO == 2:
        # Energy budget pr répartir énergie entre plante et sol : approche résistive
        ministics = energy_budget(ministics, CropParams, i, Constants, user, SoilParams)
        ministics = intermed_sat_deficit(ministics, i, meteo_gamma)

        # Eau sur plante (que si LAI > 0)
        ministics = foliage_water_evap_energy_budget(ministics, i, meteo_gamma, CropParams)

        # Transpi
        ministics = potential_transpi_energy_budget(ministics, i, meteo_gamma)
        ministics = actual_transpi(ministics, i, CropParams, SoilParams, LRAC_MATRIX)

        # Evap eau sol
        ministics = potential_soil_evap_energy_budget(
            ministics, i, meteo_gamma
        )  # besoin dos et rnetS (calculé ds energy budget)

        ministics = actual_soil_evap(ministics, i, SoilParams, user, station)

    # Calcul de l'évapotranspiration journalière (evap sol + transpi plante)
    ministics.loc[i, "et"] = ministics.loc[i, "esol_a"] + ministics.loc[i, "EP"]
    # print('et',ministics.loc[i, "et"])

    return ministics


def net_radiation(ministics, i, station, soil):
    '''
    Calcul de la radiation nette.
    Explication : calcul à partir de l'albédo du sol + végétation et des long wave radiations
    Simplifications : 
        - Pas pris en compte le couvermulch, à voir avec pac si c'est important
        - Pas utilisé le reste équation 9.15 pr modifier HUR1, faut il le faire ?
    '''
    
    if i == 0:
        ministics.loc[i, "albsol"] = soil.ALBEDO * (
            1 - 0.517 * (soil.HUR1_Init - soil.HMINF_1) / (soil.HCCF_1 - soil.HMINF_1)
        )
    else:
        ministics.loc[i, "albsol"] = soil.ALBEDO * (
            1
            - 0.517
            * (ministics.loc[i - 1, "HUR1"] - soil.HMINF_1)
            / (soil.HCCF_1 - soil.HMINF_1)
        )

        if ministics.loc[i - 1, "DOY_2EMERG"] > 0:  # pas de rad nette avant émergence

            # Calcul de l'albédo de surface = albedo sol + végétation
            ministics.loc[i, "albedolai"] = station.ALBVEG - (
                station.ALBVEG - ministics.loc[i, "albsol"]
            ) * np.exp(-0.75 * ministics.loc[i - 1, "LAI"])

            # Calcul de la fraction d'ensoleillement (entre 0 et 1)
            # rgex = radiaiton extraterrestre. Pas de formule dans STICS, renvoie à papier
            # ministics.loc[i,'fracinsol'] = (((ministics.loc[i,'trg'] / ministics.loc[i,'rgex']) - station.AANGST) / station.BANGST)
            ministics.loc[i, "fracinsol"] = (
                1.35 * (ministics.loc[i, "trg"] / ministics.loc[i, "rgex"]) - 0.35
            )  # fracinsol as in pyet0 fao

            # Calcul des long wave radiation
            # tpm = pression de vapeur d'eau. Rel_Hum = pression vapeur eau / pression vapeur saturante
            # ministics.loc[i,'Rglo'] = -0.000000004903 * ((ministics.loc[i,'tcult'] + 273.16)**4) * (0.1 + 0.9 * ministics.loc[i,'fracinsol']) * (0.56 - 0.08 * (ministics.loc[i,'ea'])**(1/2))
            ministics.loc[i, "Rglo"] = -(
                0.000000004903
                * ((ministics.loc[i, "tcult"] + 273) ** 4)
                * ministics.loc[i, "fracinsol"]
                * (0.34 - 0.14 * (ministics.loc[i, "TPM"] ** 1 / 2))
            )  # coeffs de fracinsol et ea as in fao pyet0, car on obtient des bons ordres de grandeurs des longWaveRadiation

            # Calcul de la radiation nette
            ministics.loc[i, "rnet"] = (
                1 - ministics.loc[i, "albedolai"]
            ) * ministics.loc[i, "trg"] + ministics.loc[i, "Rglo"]
            # question : qd on parle de Radiaiton du coup on ne prend pas en cpte Rglo donc là faut l'ajouter ?
        else:
            ministics.loc[i, "rnet"] = 0
        # print('albedolai',ministics.loc[i, "albedolai"])
        # print('trg',ministics.loc[i, "trg"])
        # print('rgex',ministics.loc[i, "rgex"])
        # print('TPM',ministics.loc[i, "TPM"])
        # print('Rglo',ministics.loc[i, "Rglo"])
        # print('rnet',ministics.loc[i, "rnet"])

    return ministics


def hauteur(ministics, i, CropParams, soil):
    '''
    
    Calcul de la hauteur de la canopée (partie 9.2.1.1)
    Explications : on prend LAI i-1 parce qu'on calcule tcult avant le LAI
    Simplifications :
        - laisen pas calculé encore, voir avec module sénescence.
    '''
   
    if (
        i == 0
    ):  # le 1er jour y'a pas de LAI - à corriger pr plantes pérennes qui ont feuilles au début simulation
        ministics.loc[i, "hauteur"] = CropParams.HAUTBASE
    else:
        ministics.loc[i, "hauteur"] = (
            CropParams.HAUTMAX * (1 - np.exp(-0.7 * (ministics.loc[i - 1, "LAI"])))
            + CropParams.HAUTBASE
        )  # à add : + ministics.loc[i,'laisen']
        # 0.7 = khaut = indép plante = coeff extinction reliant LAI et hauteur

    # Calcul de z0 (rugosité culture)
    ministics.loc[i, "z0"] = max(
        soil.Z0SOLNU, 0.10 * ministics.loc[i, "hauteur"]
    )  # 0.1 ou 0.13 selon passage dans STICS

    # Calcul du 'displacement height' en mètres
    ministics.loc[i, "dh"] = (
        6.6 * ministics.loc[i, "z0"]
    )  # Displacement Height = the height at which the logarithm of the wind profile projects to be zero for purposes of computing the surface layer turbulent fluxes.

    return ministics



def crop_temperature(ministics, i, CropParams):
    ''''
    Calcul de la température moyenne journalière de surface de la plante (partie 9.3.2)
    Explications : 2 approches --> empirique et bilan d'énergie. 
    Simplifications :
        - N'est implémentée que l'approche empirique.
    '''

    if i > 0:
        if ministics.loc[i - 1, "DOY_2EMERG"] > 0:  # après émergence
            # Calcul de tcultmax (qui doit être > Tmax)
            ministics.loc[i, "tcultmax"] = ministics.loc[i, "Temp_max"] + (
                ministics.loc[i, "rnet"] / 2.46 - ministics.loc[i, "et"] - 1.27
            ) / (1.68 / log(1 / ministics.loc[i, "z0"]))
            # print('tcultmax',ministics.loc[i, "tcultmax"])
            # et = evapotranspiration (esol + ep) journalière en mm

            # Calcul de tcult
            # On considère tcultmin = tmin
            ministics.loc[i + 1, "tcult"] = (
                ministics.loc[i, "tcultmax"] + ministics.loc[i, "Temp_min"]
            ) / 2

        else:  # Avant émergence : tcult = Temp air. A l'avenir mettre plutôt Tsol, c'est une meilleure approximation.
            ministics.loc[i + 1, "tcult"] = ministics.loc[i, "Temp"]
    else:  # Avant émergence : tcult = Temp air. A l'avenir mettre plutôt Tsol, c'est une meilleure approximation.
        ministics.loc[i + 1, "tcult"] = ministics.loc[i, "Temp"]
    # print('tcult',ministics.loc[i + 1, "tcult"])

    return ministics


def param_croiss_indeter(CropParams):
    """
    Calcul de paramètres utilisées pour la formation du rendement des plantes indéterminées
    Explication : paramètres calculés 1 seule fois par simulation, utilisés dans le module de croissance des plantes indéterminées.
    """

    def func(x):
        return [
            CropParams.DFPF * (1 - np.exp(-CropParams.CFPF * 0))
            + x[0] / (1 + np.exp(-CropParams.BFPF * (0 - CropParams.AFPF)))
            - x[1],
            CropParams.DFPF * (1 - np.exp(-CropParams.CFPF * 1))
            + x[0] / (1 + np.exp(-CropParams.BFPF * (1 - CropParams.AFPF)))
            - x[1]
            - CropParams.PGRAINMAXI,
        ]

    alpha, beta = fsolve(
        func, [1, 1]
    )  # A voir si on peut changer l'ini pour accélérer l'optimisation ??

    # Calcul du stade de dév du fruit = c'est juste num boite / nb boites
    dfr = [K / CropParams.NBOITE for K in range(CropParams.NBOITE)]  # eq 8.13

    return alpha, beta, dfr


def force_puit_fruits(ministics, i, ind_Z69, CropParams, alpha, beta, dfr):
    '''
    Calcul de la force du puit de carbone des fruits de chaque boite - en g/m2/j
    '''

    if i >= ind_Z69:  # calculé à partir première apparition des fruits
        for K in range(CropParams.NBOITE):
            # Calcul du potentiel de croiss des fruits (potcroifruit) ds comparti K = somme de 2 phases de croiss des fruits :
            #   - 1 = division cellul --> fonction expo avec 2 params. Augmente et tend vers DFPF quand les fruits vont vers le dernier comparti
            #   - 2 = expansion des cell formées --> fonction sigmoi avec 2 params
            # Eq 8.11 et 8.12. Alpha et beta déterminés par optimisation pour potcroifruit connu pour DFR = 0 et 1
            ministics.loc[i, f"potcroifruit{K}"] = (
                CropParams.PGRAINMAXI
                * CropParams.DFPF
                * (1 - np.exp(-CropParams.CFPF * dfr[K]))
                + alpha / (1 + np.exp(-CropParams.BFPF * (dfr[K] - CropParams.AFPF)))
                - beta
            )

            # Si on garde que sigmoïde, les param changent bcp forme courbe : sigmoide/liné
            # Rq :  potcroi on s'en sert pas mais c'est le alpha qui ns intéresse j'ia l'impression

            # Calcul du taux de dév journalier
            # ça donne qqch entre 0 et 1 car normalisé par durée max
            ministics.loc[i, "devjour"] = (
                ministics.loc[i, "tcult"] - CropParams.TDMIN
            ) / CropParams.DUREEFRUIT

            # Calcul de la force du puit des fruits du compartiment K - Eq 8.14
            Y = np.exp(-CropParams.BFPF * (dfr[K] - CropParams.AFPF))
            ministics.loc[i, f"fpft{K}"] = (
                CropParams.PGRAINMAXI
                * ministics.loc[i, "devjour"]
                * (
                    CropParams.DFPF
                    * CropParams.CFPF
                    * np.exp(-CropParams.CFPF * dfr[K])
                    + CropParams.BFPF * alpha * Y / (1 + Y) ** 2
                )
            )

        # Comment calculer force du puit global --> faudrait pondérer par nb de fruits par boite mais on ne l'a pas encore calculé le jour j donc je prends j-1.
        # IMPRO
        if (
            i == ind_Z69
        ):  # le 1er jour d'apparition de fruits, c'est le fpft de la boite 0
            ministics.loc[i, "fpft"] = ministics.loc[i, "fpft0"]
        else:  # on fait moy des fpft de chaque boite pondérés par le nb de fruit ds chaque boite
            ministics.loc[i, "fpft"] = sum(
                [
                    ministics.loc[i, f"fpft{K}"] * ministics.loc[i - 1, f"nfruit{K}"]
                    for K in range(CropParams.NBOITE)
                ]
            ) / sum(
                [ministics.loc[i - 1, f"nfruit{K}"] for K in range(CropParams.NBOITE)]
            )

    return ministics


def yield_formation_indeterminate(
    ministics, CropParams, i, ind_AMF, ind_Z69, ind_NOU, alpha, beta, dfr
):
    '''
    Module de formation du rendement pour les plantes indéterminées

    Explications:
    Feuilles des plantes indéter continuent de croitre pdt que organes récoltés sont produits et croissent.
    --> Modélisé par interaction trophique entre groupes d'organes et organes récoltés, approche sources-puits en définissant notion stress trophique
    Méthode 'boxcartrain' :
       - les fruits passent par des boites associées à des âges phénologiques.
       - Le temps passé par boite dépend de la température. Définition d'une durée thermique à passer dans chaque boite (la même pr ttes boites)
       - Dans chaque boite, la croiss du fruit = force du puit * ratio source-puit. Force du puit = dérivée fonction logistique prenant cpte potentiel genetic de croiss
    Attention : si période de calcul du nb de grain est longue, alors il faut bcp de compartiments (fruits à des âges différents) = cohérence entre ces deux param
    Notion d'inflorescence utile quand régulation technique/trophique se fait au niveau de l'inflorescence. Ex : vigne (grappes)

    # Choix du nb de boites : dépend de durée de fructification, de la dynamique de croiss des fruits, et date début dynamique eau dans organes prélevés

    '''

    ### 1. Nb de fruits apparus + répartition dans boites ###
    # Les fruits apparaissent sur la période AMF - NOU
    if (i >= ind_Z69) & (
        i < ind_NOU
    ):  # <= ou < --> à 1 jour près ça peut ajouter des fruits attention et ça joue énormément sur rendement.
        # Calcul du stress trophique impactant le nb de fruit
        # Pour avoir sourcepuits, faut avoir calculé fpft avant mais on le calcule après A VOIR
        # Idée : ratio source-puit * pente courbe avec indice min = 0 et indice max = 1
        a = (1 - 0) / (CropParams.SPFRMAX - CropParams.SPFRMIN)
        b = -a * CropParams.SPFRMIN
        ministics.loc[i, "spfruit"] = min(
            1, max(0, ministics.loc[i, "sourcepuits"] * a + b)
        )  # pq c'est fpft et pas sourcepuits ici ??

        # A faire : ajouter fpft dans calcul sourcepuits pr plantes indéterminées

        # Calcul du nb d'inflorescences par plant (pour vigne : inflo = grappe)
        # Idée : le nb d'inflo est soit fixé (CODCALINFLO=1), soit (CODCALINFLO=2) calculé comme une fonction du 'statut trophique' à un stade phéno précoce
        # Intuition : biomasse aérienne au mom accélération max * réserve pérenne de la saison précéd * facteur diminuant avec densité semis et augmentant avec param réglable
        if CropParams.CODCALINFLO == 1:
            ministics.loc[i, "nbinflo_recall"] = CropParams.NBINFLO
        elif CropParams.CODCALINFLO == 2:
            ministics.loc[i, "nbinflo_recall"] = min(
                CropParams.INFLOMAX,
                CropParams.PENTINFLORES
                / CropParams.DENSITESEM
                * (ministics.loc[ind_AMF, "MASEC"] + ministics.loc[0, "resperenne"]),
            )

        # Calcul du nb de fruits apparus le jour i
        # Intuition : nb max de fruits par inflo * temp dev * nb inflo par plant * densité * stress trophique * stress gel
        ministics.loc[i, "nfruitnou"] = (
            CropParams.AFRUITPOT
            * (ministics.loc[i, "Sum_UPVT"] - ministics.loc[i - 1, "Sum_UPVT"])
            * ministics.loc[i, "nbinflo_recall"]
            * CropParams.DENSITESEM
            * ministics.loc[i, "spfruit"]
            * ministics.loc[i, "fgelflo"]
        )

        # Diff entre upvt et udevcult ? La diff de somme ça marche de sûr ? à reprendre le calcul et comparer à STICS

    if i >= ind_Z69:
        # Répartion des fruits dans les compartiments - IMPROVISE y'a pas d'équation dans STICS
        duree_par_boite = (
            CropParams.DUREEFRUIT / CropParams.NBOITE
        )  # quand somme upvt depuis DRP dépasse cette durée, ça change de boite

        for K in range(
            CropParams.NBOITE
        ):  # pr chaque boite, on va calculer le nb de fruits présents le jour i
            nbfruits = 0
            for j in range(
                ind_Z69, i
            ):  # pr chaque jour j avant le jour i depuis DRP, si la somme thermique depuis ce jour j correspond à la durée de vie associée à la boite K, alors on ajoute les fruits produits le jour j à la boite K du jour i
                if (
                    ministics.loc[i, "Sum_UPVT_DRP"]
                    - ministics.loc[j - 1, "Sum_UPVT_DRP"]
                    >= K * duree_par_boite
                ) and (
                    (
                        ministics.loc[i, "Sum_UPVT_DRP"]
                        - ministics.loc[j - 1, "Sum_UPVT_DRP"]
                        < (K + 1) * duree_par_boite
                    )
                    or K == CropParams.NBOITE - 1
                ):
                    nbfruits = nbfruits + ministics.loc[j, "nfruitnou"]

            ministics.loc[
                i, f"nfruit{K}"
            ] = nbfruits  # on attribue à la boite K du jour i, les fruits produits les jours < i
        ministics.loc[i, "nfruit0"] = (
            ministics.loc[i, "nfruit0"] + ministics.loc[i, "nfruitnou"]
        )  # on ajoute les fruits produits le jour i à la boite 0, car la boucle d'avant ne répartit que tous les fruits produits les jours < i

    ### 2. Remplissage des fruits ###
    # Idée = 1er comparti = boite d'apparition des fruits, et dernier comparti =  pas de croiss de fruit = fruit a atteint maturité

    if i >= ind_Z69:  # A partir du 1er jour d'apparition de fruits
        for K in range(CropParams.NBOITE):
            # Calcul de la croiss du fruit ds le compartiment K - en g/m2
            # Intuition : nb de fruits comparti K * force puit (= croiss par fruit) * stress trophique impactant dev fruit = ratio source-puits * stress thermique
            #  ils assimilent le stress qui est spfruit au ratio source-puit je capte pas ??
            ministics.loc[i, f"croifruit{K}"] = (
                ministics.loc[i, f"nfruit{K}"]
                * ministics.loc[i, "fpft"]
                * ministics.loc[i, "sourcepuits"]
                * ministics.loc[i, "ftempremp"]
            )  # eq 8.10

        # Calcul du ratio d'allocation des assimilats aux fruits et correction des valeurs de croiss si dépassement du seuil
        # Intuition : ça contrôle 2 choses :
        #   - y'a un max d'allocation (en %) de la biom aux fruits à ne pas dépasser, sinon on réduit la croiss des fruits
        #   - mm si allocfrmax = 1 (tte la biomasse prod va dans fruits), on ne veut pas que la masse accumulée par les fruits au jour i soit supérieure à biomasse produite jour i
        # --> pr réduire la croiss des fruits si 1 de ces 2 conditions n'est pas vérifiée, on diminue le facteur sourcepuits

        ministics.loc[i, "allocfruit"] = (
            sum([ministics.loc[i, f"croifruit{K}"] for K in range(CropParams.NBOITE)])
            / ministics.loc[i, "DELTAMS"]
        )

        if ministics.loc[i, "allocfruit"] > CropParams.ALLOCFRMAX:
            ministics.loc[i, "sourcepuits"] = (
                ministics.loc[i, "sourcepuits"]
                * ministics.loc[i, "allocfruit"]
                / CropParams.ALLOCFRMAX
            )

            for K in range(CropParams.NBOITE):
                ministics.loc[i, f"croifruit{K}"] = ministics.loc[
                    i, f"croifruit{K}"
                ] / (
                    ministics.loc[i, "allocfruit"] / CropParams.ALLOCFRMAX
                )  # en gros on diminue chaque croiss par le ratio entre allocfruit et allocfrMAX pour que allocfruit = allocfrMAX
        # vérifier que ça marche bien

        # Calcul improvisé de masse de fruit totale - je prends les 2 mm noms de variable que pr croiss déterminée
        ministics.loc[i, "deltags"] = sum(
            [ministics.loc[i, f"croifruit{K}"] for K in range(CropParams.NBOITE)]
        )
        ministics.loc[i, "mafruit"] = ministics.loc[ind_Z69:i, "deltags"].sum()

    return ministics


# ajout *rfpi dans LAI si photoP diminue d'un jour à l'autre
# Module à lancer après LAI et après prod biomasse



def cumultative_partitioning(ministics, CropParams, Constants, i, ind_EMERG):
    '''
    MODULE EN DEVELOPPEMENT

    Partionnement de la biomasse avec option de cumul des réserves
    Explications :
    utilisé quand code_acti_reserve=2
    Marche pr plantes pér/croiss indéter et pr plantes croiss déter
    '''

    #######################################
    ### Biomasse des organes végétatifs ###
    #######################################

    # 1. Feuilles

    # Calcul de la biomasse de feuille verte
    # Intuition : (m2 de feuille / m2 de sol) / (m2 de feuille / g de feuille) --> donne g de feuille / m2 de sol
    ministics.loc[i, "mafeuilverte"] = (
        ministics.loc[i, "LAI"] / CropParams.SLAMAX * 100
    )  # eq 7.7

    # Calcul de la biomasse de feuilles jaunes depuis émergence, avant = 0
    # FAUT AVOIR calculé sénescence avant
    ministics.loc[i, "mafeuiljaune"] = ministics.loc[ind_EMERG:i, "dltamsen"].sum()

    # Calcul feuilles vertes + jaunes
    ministics.loc[i, "mafeuil"] = ministics.loc[i, "mafeuilverte"] + ministics.loc[i, "mafeuiljaune"]

    # Calcul de la qtté journalière de feuilles tombées et le cumulé
    ministics.loc[i, "dltamstombe"] = CropParams.ABSCISSION * ministics.loc[i,'dltamsen']
    ministics.loc[i, "mafeuiltombe"] = ministics.loc[i, "dltamstombe"]

    # 2. Tiges

    # Calcul de la biomasse de la 'partie structurale' de la tige
    # Intuition = juste une proportion de la masse de la feuille
    ministics.loc[i, "matigestruc"] = (
        CropParams.TIGEFEUIL * ministics.loc[i, "mafeuilverte"]
    )  # eq 7.7

    # 3. Réserves temporaires

    # Calcul des réserves temporaires
    if (i < ind_EMERG) and (CropParams.CODEPERENNE == 2): # Pr les plantes pérennes, on initialise à RESTEMP0
         # condition faite à chaque itération = pas optimisé
        ministics.loc[i, "restemp"] = CropParams.RESTEMP0
    elif i >= ind_EMERG:
    # Faut avoir calculé MASEC, MZFRUIT, MAENFRUIT avant
        ministics.loc[i, "restemp"] = ( ministics.loc[i, "MASEC"] - ministics.loc[i, "mafeuil"] - ministics.loc[i, "matigestruc"])

    # 4. Atteinte de la taille max du compartiment de réserve
    # cf 7.3.1.1.2 : pq on met dltams à 0 ? si toutes réserves temporaires sont pleines, on considère que y'a plus de production de matière ?
    # Attention densité : dépend pas de t
    if ministics.loc[i, "restemp"] > 10 * CropParams.RESPLMAX * CropParams.DENSITESEM:
        ministics.loc[i, "DELTAMS"] = 0 # pas trop compris, c encore utilisé ici ?


    ###################################
    ### Remobilisation des réserves ###
    ###################################
    # Prio réserves pérennes ou temporaires ?
    # Ici on parle que des rés pérennes ?
    # 2 types de réserves : aériennes (restemp) via remobilj et organes de stockage des plantes pérennes (resperenne) via dltaremobil
    #   --> aussi notion de ce qui a été accumulé durant ce cycle = remobilj, et ce qui a été accumulé à un autre cycle = dltaremobil

    # 1. Ratio source/puit = production biomasse / demande par organe végé et reprod

    # Calcul de la demande des organes végé = fpv
    # Intuition : demande biomasse dépend du LAI produit et est d'autant plus petit qu'il y a une forte surface par gramme de biomasse.
    # Question : pq y'a ratioTF ici sachant qu'on parle de 'végétatif' -> on veut séparer tige/feuille ?
    ministics.loc[i, "fpv"] = (
        ministics.loc[i, "DELTAI"]
        * 10000
        / (CropParams.SLAMIN / (1 + ministics.loc[i, "ratioTF"]))
    )

    # Calcul du ratio source/puit : 1er calcul car il faut ça pr la suite, mais on le recalcule ensuite. Sorte de boucle imbriquée faut choisir un 1er calcul à faire
    # Pr les plantes annuelles, fpft = 0
    # Dans le code : y'a un min(1,calcul)
    ministics.loc[i, "sourcepuits1"] = min(1,
        ministics.loc[i, "DELTAMS"]
        / (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"])
        if ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"] != 0
        else 1
    )

    # 2. Mobilisation des réserves
    # on parle que des réserves pérennes dans le 7.4.2.1 : pourtant pr moi dans le formalisme 'cumulative' y'a pas de réserves pérenne, enfin pas séparées des rés temp
    # Diff entre remobilj et dltaremobil ? 1ère du cycle actuel, et 2ème de l'autre cycle ?

    # Calcul de la biomasse remobilisée depuis réserves pérennes (??), avec un max possible
    if ministics.loc[i, "sourcepuits1"] < 1:
        ministics.loc[i, "remob"] = (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"]) / 100 - ministics.loc[i, "DELTAMS"]
    ministics.loc[i, "dltaremobil"] = min(ministics.loc[i, "remob"], CropParams.REMOBRES * ministics.loc[i, "restemp"])
    ministics.loc[i, "culdltaremobil"] = ministics.loc[0:i, "dltaremobil"].sum()



    # Ratio source-puits pr prendre en compte mobilisation des réserves
    # Dans DETLAMS on n'a pas ajouté encore dltaremobil
    # Pas calculé dltaremobil, et remobilj encore
    ministics.loc[i, "sourcepuits"] = (
        (
            ministics.loc[i, "DELTAMS"]
            + ministics.loc[i, "dltaremobil"]
            + ministics.loc[i, "remobilj"]
        )
        / (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"])
        if (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"]) != 0
        else 1
    )





def daily_partitioning(ministics, CropParams, Constants, i, ind_EMERG):
    '''
    MODULE EN DEVELOPPEMENT

    utilisé quand code_acti_reserve=1
    Ce formalisme n'est développé dans STICS que pr les plantes pérennes pr le moment.
    '''
    pass

def biomass_partitioning(ministics, CropParams, Constants, i, ind_EMERG):
    '''
    MODULE EN DEVELOPPEMENT
    '''
    # Distinction des organes pérennes et non pérennes. Aucune allocation aux organes pérennes sauf en cas de remobilisation : réserves temp (non pér) --flux--> réserves (pér)

    if i == 1:
        # Calcul de la biomasse de feuille verte
        # Intuition : (m2 de feuille / m2 de sol) / (m2 de feuille / g de feuille) --> donne g de feuille / m2 de sol
        ministics.loc[i, "mafeuilverte"] = (
            ministics.loc[i, "LAI"] / CropParams.SLAMAX * 100
        )  # eq 7.7

        # Calcul de la biomasse de la 'partie structurale' de la tige
        # Intuition = juste une proportion de la masse de la feuille
        ministics.loc[i, "matigestruc"] = (
            CropParams.TIGEFEUIL * ministics.loc[i, "mafeuilverte"]
        )  # eq 7.7

        # Calcul de la biommasse non pérenne
        # Hypothèse : pas de fruit, ni réserves, ni feuille jaune A CHANGER PR PLANTE PÉRENNE SIMULÉE SANS GRAINE
        ministics.loc[i, "masecnp"] = (
            ministics.loc[i, "mafeuilverte"] + ministics.loc[i, "matigestruc"]
        )  # eq. 7.7

    elif i >= ind_EMERG:
        # -----------------------
        # A. Organes non pérennes
        # -----------------------

        # 1. Feuilles
        # Calcul de la biomasse journalière de la partie structurale des feuilles vertes
        # hyp = l'allocation a comme priorité le lieu où est produit la biomasse = les feuilles
        # = basé sur l'augmentation de LAI journalier
        # Intuition : (m2 feuille / m2 sol) / (m2 feuille / g de feuille)
        ministics.loc[i, "dltafv"] = (
            ministics.loc[i, "DELTAI"] / CropParams.SLAMAX * 100
        )

        # Calcul de la biomasse journalière REELLE de la partie structurale des feuilles verte
        # Prend en cpte le fait que ça ne peut pas être supérieur à la biomasse totale produite, donc on prend le min
        ministics.loc[i, "dltafv"] = min(
            ministics.loc[i, "dltafv"], ministics.loc[i, "DELTAMS"]
        )

        # Eq 7.9 : ça se mort la queue, il faut maj deltai si le 'deltams est insuffisant pour garantir le potentiel de grandiss du deltai' alors que deltams est justement calculé à partir du deltai
        # donc on calcule deltai, puis à partir de ça deltams puis on regarde si ce deltams est suffisant pour deltai c pas possible ??

        # Calcul de la masse de feuille verte totale = j-1 + produciton jour i
        ministics.loc[i, "mafeuilverte"] = (
            ministics.loc[i - 1, "mafeuilverte"] + ministics.loc[i, "dltafv"]
        )

        # Appel du module sénescence avant calcul mafeuiljaune A VOIR
        ministics = senescence_stress(ministics, CropParams, i, ind_EMERG)
        ministics = senescence(ministics, CropParams, i, ind_EMERG)

        # Calcul de la biomasse de feuilles jaunes depuis émergence, avant = 0
        ministics.loc[i, "mafeuiljaune"] = ministics.loc[ind_EMERG:i, "dltamsen"].sum()

        # Calcul de la biomasse de feuille (verte + jaune). Mais on n'a pas calculé biomasse jaune encore ??
        ministics.loc[i, "mafeuil"] = (
            ministics.loc[i, "mafeuilverte"] + ministics.loc[i, "mafeuiljaune"]
        )

        # 2. Tiges
        # Calcul de l'allocation à la tige s'il reste de la biomasse à allouer
        # Intuition : s'il reste de la biomasss après allocation à feuille, la biomasse allouée à la tige est une proportion de la biomasse tige+racine
        # Questions : 1. on prend ratioTF * dltafv mais on enlève pas (1-ratioTF)*dltafv à la biomasse de feuille.
        # 2. pq prendre mafeuil j < mafeuil j-1 -> c'est pas égal à dltafv ? Lien avec feuilles jaunes ?
        # Problème : ratioTF n'est jamais initialisé ni calculé si ? Sauf période de racourciss jours ou quand photoP < photoP base
        if ministics.loc[i, "dltafv"] < ministics.loc[i, "DELTAMS"]:
            ministics.loc[i, "dltat"] = max(
                ministics.loc[i, "ratioTF"]
                * (ministics.loc[i, "mafeuil"] - ministics.loc[i - 1, "mafeuil"]),
                0,
            )
        # On prend le min entre la biomasse calculée et  biomasse totale - biomasse feuille (= biomasse allouée à tige forcément inf à biomasse tot - biomass feuille)
        ministics.loc[i, "dltat"] = min(
            ministics.loc[i, "dltat"],
            ministics.loc[i, "DELTAMS"] - ministics.loc[i, "dltafv"],
        )

        # Calcul de la biomasse totale de tige au jour i = proportion cste * biomasse de feuilles. Pareil que d'ajouter à matigestruc j-1 le dltat
        ministics.loc[i, "matigestruc"] = (
            CropParams.TIGEFEUIL * ministics.loc[i, "mafeuil"]
        )
        # Pas compris paragraphe après Eq 7.4

        # 3. Réserves temporaires
        # Calcul de l'allocation de biomasse restante aux réserves temporaires (contenues dans feuille + tige)
        ministics.loc[i, "dltares"] = (
            ministics.loc[i, "DELTAMS"]
            - ministics.loc[i, "dltafv"]
            - ministics.loc[i, "dltat"]
        )
        ministics.loc[i, "restemp"] = (
            ministics.loc[i - 1, "restemp"] + ministics.loc[i, "dltares"]
        )

        # Calcul de la biomasse remobilisée pdt la sénescence = ce qui a été produit * % de biomasse sénescente * % de feuilles vertes dans biomasse tot non sénesc
        # On calcule cette biomasse remob à partir de dltfv (la biom de fv produite jour i) mais à quel moment faut-il ajouter ça à restemp ? c'est pas au jour i que c'est remobilisé
        # pfeuilverte = % de feuilles vertes parmi la biomasse non sénescente. Où est-il calculé ? J'ai choisi de le calculer comme (dltafv / DELTAMS) donc les DELTAMS s'annulent
        ministics.loc[i, "dltaremobsen"] = (1 - CropParams.RATIOSEN) * ministics.loc[
            ind_EMERG:i, "dltafv"
        ].sum()  # eq 4.12
        # Puis ce qui est remob va dans restemp. Dans stics Eq 7.12 le calcul ça écraserait ce qui a déjà été calculé dans restemp juste avant, donc moi je somme les 2
        ministics.loc[i, "restemp"] = (
            ministics.loc[i, "restemp"] + ministics.loc[i, "dltaremobsen"]
        )

        # 4.Masse sèche végé

        # Calcul de la biomasse végétative
        # Si on soustrait magrain à masecnp, il faut ajouter magrain(i) à masecnp(i) àa la fin de Yield formation du jour i. STICS n'est pas clair.
        # Mon choix : maseveg = feuilles + tige + res temp
        # ministics.loc[i,'masecveg'] = ministics.loc[i,'masecnp'] - ministics.loc[i-1,'magrain']
        ministics.loc[i, "masecveg"] = (
            ministics.loc[i, "mafeuil"]
            + ministics.loc[i, "restemp"]
            + ministics.loc[0:i, "dltat"].sum()
        )

        # max réserves mobilisables depuis organes stockage pérennes * proportion de feuilles vertes parmi fv + fj * masse organes végétatifs. MAIS feuille c'est pas pérenne
        # Calcul je capte pas trop
        ministics.loc[i, "restempmax"] = (
            CropParams.PROPRES
            * ministics.loc[i, "mafeuilverte"]
            / ministics.loc[i, "mafeuil"]
            * ministics.loc[i, "masecveg"]
        )

        # Calcul du surplus qui ne peut pas être stocké dans les réserves temporaires (sera stocké dans organes pérennes) et maj des réserves temporaires en conséquence.
        # Eq 7.15
        if (
            ministics.loc[i, "restemp"] > ministics.loc[i, "restempmax"]
        ):  # si les réserves temp dépassent le seuil max de stockage -> on stocke le surplus dans les organes pérennes
            ministics.loc[i, "dltarestemp"] = (
                ministics.loc[i, "restemp"] - ministics.loc[i, "restempmax"]
            )  # = surplus qui sera stocké dans les organes pérennes
            ministics.loc[i, "restemp"] = ministics.loc[
                i, "restempmax"
            ]  # les réserves = réserves max

        # -----------------------
        # B. Organes pérennes (= réserves ou parties structurelles des organes pérennes)
        # -----------------------
        # Hyp = réserves stockées en priorité dans organes pérennes existant, et quand limite max de stockage atteinte --> stockage dans partie structurelles
        # Autre hyp = les organes pérennes ne sont pas des puits avec une demande en biomasse, ils stockent juste le surplus non alloué aux organes végé aériens et reproducteurs

        # Calcul du seuil de réserve max dans les organes pérennes --> tout ce qui est en surplus est ensuite stocké dans Parties structurelles
        # % max de réserve mobilisable * biomasse pérenne (réserves + parties structurelles organes pérennes). Pq on utilise maperenne et pas resperenne car y'a que ça de mobilisable.
        # maperenne jamais calculé et pq c'est maperenne et dlarestemp comme les équaitons suivantes ?
        ministics.loc[i, "resperennemax"] = (
            CropParams.PROPRESP * ministics.loc[i, "maperenne"]
        )

        # Calcul de l'allocation de biomasse aux réserves Eq 7.15
        # Idée = on ajoute aux réserves pérennes le surplus de réserves temp, dans la mesure du possible de max de réserves pérennes
        # Si le seuil max de réserves pérennes n'est pas dépassé en ajoutant aux réserves pérennes le dltarestemp, on ne fait rien d'autre
        # resperenne jamais calculé ?
        if (
            ministics.loc[i, "resperenne"] + ministics.loc[i, "dltarestemp"]
            < ministics.loc[i, "resperennemax"]
        ):
            ministics.loc[i, "resperenne"] = (
                ministics.loc[i, "resperenne"] + ministics.loc[i, "dltarestemp"]
            )  #

        else:  # si le seuil est dépassé avec dltarestemp, on applique un coeff pour que tout le surplus de rés temp ne soit pas ajouté aux rés pérennes
            ministics.loc[i, "resperenne"] = (
                ministics.loc[i, "resperenne"]
                + CropParams.PROPRESP * ministics.loc[i, "dltarestemp"]
            )
            # Rien sur resperennestruc ?? on met pas ce qui reste si dépasse max ? il faudrait faire resperennestruc = resperennemax - resperenne dans le else nan ?
            # j'improvise
            ministics.loc[i, "resperennestruc"] = (
                ministics.loc[i, "resperennestruc"]
                + (1 - CropParams.PROPRESP) * ministics.loc[i, "dltarestemp"]
            )

        # Maj de la masse pérenne. Comment l'initialiser en début de simulation ?? Si c'est 0, le resperennemax va rester à 0 tt le temps
        ministics.loc[i, "maperenne"] = (
            ministics.loc[i, "resperenne"] + ministics.loc[i, "resperennestruc"]
        )

        # ------------------------------
        # C. Remobilisation des réserves
        # ------------------------------
        # Biomasse remob vient en prio des réserves pér puis réserves temp d'après code, mais l'inverse écrit dans STICS.
        # La remobilisation arrive si le ratio source-puits est inf à 1

        # 1. Ratio source/puit = production biomasse / demande par organe végé et reprod

        # Calcul de la demande des organes végé = fpv
        # Intuition : demande biomasse dépend du LAI produit et est d'autant plus petit qu'il y a une forte surface par gramme de biomasse.
        # Question : pq y'a ratioTF ici sachant qu'on parle de 'végétatif' -> on veut séparer tige/feuille ?
        ministics.loc[i, "fpv"] = (
            ministics.loc[i, "DELTAI"]
            * 10000
            / (CropParams.SLAMIN / (1 + ministics.loc[i, "ratioTF"]))
        )

        # Calcul du ratio source/puit (diff entre 1 et rien ??)
        # On considère que des plantes déterminées = pas de fpft
        ministics.loc[i, "sourcepuits1"] = (
            ministics.loc[i, "DELTAMS"]
            / (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"])
            if ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"] != 0
            else 1
        )

        # 2. Remobilisation des réserves pérennes
        # Principe = calcul des carbohydrates remobilisés chaque jour --> réserves temporaires en prio, puis réserves pérennes ??? pas l'inverse

        # Si demande > source, on calcule l'écart
        if (ministics.loc[i, "sourcepuits1"] < 1) and (
            ministics.loc[i, "resperenne"] > 0
        ):  # si demande > source et que réserves pérennes non nulles
            # Eq 7.22 j'ai enlevé le /100 sinon ça devient négatif la remob mais j'ai peut être pas compris
            ministics.loc[i, "remob"] = (
                ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"]
            ) - ministics.loc[
                i, "DELTAMS"
            ]  # biomasse remob = demande - source
        else:  # si y'a pas de demande ou plus de réserve pérenne = rien n'est remobilisé. Impro
            ministics.loc[i, "remob"] = 0

        # On garde le min entre écart demande-source, le max de ce qui peut être remobilisé chaque jour et les réserves pérennes (à quoi sert dernier terme)
        ministics.loc[i, "dltaremobilbrut"] = min(
            ministics.loc[i, "remob"],
            CropParams.REMOBRES * ministics.loc[i, "resperenne"],
            ministics.loc[i, "resperenne"],
        )

        # On maj les réserves pérennes en y retirant ce qui est mobilisé
        ministics.loc[i, "resperenne"] = (
            ministics.loc[i, "resperenne"] - ministics.loc[i, "dltaremobilbrut"]
        )

        # Calcul du flux net de remobilisation des carbohydrates
        # Idée = une partie de ce qui est remobilisé est perdu par le biais de la respiration demander à PAC précisions.
        ministics.loc[i, "dltaremobil"] = (
            CropParams.EFREMOBIL * ministics.loc[i, "dltaremobilbrut"]
        )
        ministics.loc[i, "dltaCO2resperenne"] = (
            0.4 * (1 - CropParams.EFREMOBIL) * ministics.loc[i, "dltaremobilbrut"]
        )

        # 3. Remobilisation des réserves temporaires (pas prio du coup?)
        # Idée = si biomasse produite + celle remob depuis rés pérennes ne suffit pas à satisfaire demande, les rés temp sont remobilisées

        ministics.loc[i, "sourcepuit2"] = (
            (ministics.loc[i, "DELTAMS"] + ministics.loc[i, "dltaremobil"])
            * 100
            / (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"])
            if (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"]) != 0
            else 1
        )

        if (ministics.loc[i, "sourcepuit2"] < 1) & (
            ministics.loc[i, "restemp"] > 0
        ):  # si demande > source et réserv temp non nulles
            # Eq 7.24 j'ai enlevé le /100 sinon ça devient négatif la remob mais j'ai peut être pas compris
            ministics.loc[i, "remob"] = (
                (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"])
                - ministics.loc[i, "DELTAMS"]
                + ministics.loc[i, "dltaremobil"]
            )  # pq + à la fin

        # On garde le min entre écart demande-source, max mobilisable, et les réserves temp (dernier terme inutile si REMOBRES < 1)
        ministics.loc[i, "remobilj"] = min(
            ministics.loc[i, "remob"],
            CropParams.REMOBRES * ministics.loc[i, "restemp"],
            ministics.loc[i, "restemp"],
        )

        # Maj des réserves temp en retirant ce qui a été mobilisé
        ministics.loc[i, "restemp"] = (
            ministics.loc[i, "restemp"] - ministics.loc[i, "remobilj"]
        )

        # Ratio source-puits pr prendre en compte mobilisation des réserves
        # Dans DETLAMS on n'a pas ajouté encore dltaremobil
        ministics.loc[i, "sourcepuits"] = (
            (
                ministics.loc[i, "DELTAMS"]
                + ministics.loc[i, "dltaremobil"]
                + ministics.loc[i, "remobilj"]
            )
            / (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"])
            if (ministics.loc[i, "fpv"] + ministics.loc[i, "fpft"]) != 0
            else 1
        )

        # Ajout des biomasses remob à la biomasse ? je suis pas sûr
        ministics.loc[i, "DELTAMS"] = (
            ministics.loc[i, "DELTAMS"]
            + ministics.loc[i, "dltaremobil"]
            + ministics.loc[i, "remobilj"]
        )
        ministics.loc[i, "MASEC"] = ministics.DELTAMS[0:i].sum() / 100

        # -------------------
        # D. Organes récoltés
        # -------------------

        # Calcul du nb et masse des grains (chap 8)

        # Calcul masse enveloppe = nb fruit *
        # ministics.loc[i,'maenfruit'] = ministics.loc[i,'nbgrains'] * CropParams.ENVFRUIT * CropParams.PGRAINMAXI * 0.01

        # -------------------
        # E. Racines
        # -------------------
        # A FAIRE : comment identifier le stade i tel que c'est la récolte. C'est maturation
        # Avec option de profil racinaire, calcul de la biomasse racinaire à la récolte.
        # Intuition : la masse racinaire à la récolte est une proportion de la masse aérienne
        # Question : pas du tout basé sur profondeur et densité des racines.
        ministics.loc[i, "msrac"] = (
            ministics.loc[i, "masecnp"] * Constants.PROPRAC + Constants.Y0MSRAC
        )

    return ministics


def temp_stress(ministics, CropParams, i, ind_Z69):
    """
    Calcul des indices de stress thermique agissant sur la RUE et la formation du rendement.
    Améliorations : 
        - Stress basés sur Temp air mais à terme il faudrait replacer par température de la plante (tcult)
    """

    ### STRESS AGISSANT SUR LA FORMATION DU RENDEMENT ###

    # Calcul de l'indice de stress de gel impactant le nb de fruits (fgelflo) pdt la période de remplissage des grains
    if i >= ind_Z69:
        ministics.loc[i, "fgelflo"] = frost_stress(
            ministics.loc[i, "Temp_min"],
            CropParams.TGELFLO90,
            CropParams.TGELFLO10,
            CropParams.TLETALE,
            CropParams.TDEBGEL,
        )

    # Calcul de l'indice de stress thermique impactant rempliss des grains
    ministics.loc[i, "ftempremp"] = thermal_stress_fruitfilling(
        ministics.loc[i, "Temp_min"],
        ministics.loc[i, "tcultmax"],
        CropParams.TMINREMP,
        CropParams.TMAXREMP,
    )

    ### STRESS AGISSANT SUR LA RUE ###
    # Calcul de l'indice de stress thermique impactant la RUE
    ministics.loc[i, "ftemp"] = thermal_stress_rue(
        ministics.loc[i, "tcult"],
        CropParams.TEMIN,
        CropParams.TEOPT,
        CropParams.TEOPTBIS,
    )

    ### STRESS AGISSANT SUR LA SENESCENCE ###

    return ministics


def harvested_organs_number_determinate(ministics, i, CropParams, ind_Z69):
    """'
    Calcul du nombre de fruits pour les plantes à croissance déterminée.
    Explications : 
    Spécificité de STICS : calcul du nb d'organes récoltés AVANT l'allocation de biomasse à ces organes, et joue le rôle de limite génétique du rendement.
    Calcul du nb de grains avant le stade de remplissage des grains, et sa réduction par le gel pdt le remplissage
    """

    # Au stade 69 ('fin floraison, début formation fruit visible') --> Calcul du nb de grains.
    # 71 = 10% des fruits ont leur taille max ou les fruits ont 10% de leur taille max
    if i == ind_Z69:
        # Calcul de la vitesse moy de croiss de la canopée sur la période de calcul des grains jusqu'au début de remplissage des grains
        # Idée : on fixe une durée en jours (nbjgrain) sur laquelle va être calculée le nb de grains, et on calcule la vitesse moyenne de croiss de la canopée sur cette période.
        # Intuition : c'est la moyenne d'augmentation de biomasse sur une période d'intérêt qui correspond à la période où le nb de grains va être déterminé.
        # Question : comment est trouvé NBJGRAIN, explication biologique ??
        vitmoy = (
            ministics.loc[ind_Z69 - CropParams.NBJGRAIN + 1 : ind_Z69, "DELTAMS"].sum()
            / CropParams.NBJGRAIN
        )  # Eq 8.1

        # Calcul du nb de grains (calculé à date ind_Z69) = nb de grains final par plante
        # Intuition : y'a nb de grains de base, qui augmente linéairement avec vit moy de croiss de canopée, et y'a un max atteint = sorte de sigmoïde avec pente droite
        # stress sur période intégrés via vitmoy, et les autres param sont espèce-dpdt = traduisent la sensibilité de l'espèce aux conditions de croiss
        ministics.loc[i, "nbgrains"] = (
            CropParams.CGRAINV0 + CropParams.CGRAIN * vitmoy
        ) * CropParams.NBGRMAX  # Eq 8.2
        ministics.loc[i, "nbgrains"] = min(
            ministics.loc[i, "nbgrains"], CropParams.NBGRMAX
        )  # on garde min entre calculé et max possible
        ministics.loc[i, "nbgrains"] = max(
            ministics.loc[i, "nbgrains"], CropParams.NBGRMIN
        )  # on garde max entre calculé et min possible

    # Remarque : très important de bien choisir la période d'intérêt, parce que c'est la croiss de biomasse moy sur cette période qui va faire qu'on aura bcp de grains ou pas.

    # Après stade 69 --> Calcul du nb grains détruits par le gel et maj du nb de grains, après la date de début de remplissage
    elif i > ind_Z69:  # si on est après le début de rempliss
        ministics.loc[i, "nbgraingel"] = ministics.loc[i - 1, "nbgrains"] * (
            1 - ministics.loc[i, "fgelflo"]
        )
        ministics.loc[i, "nbgrains"] = (
            ministics.loc[i - 1, "nbgrains"] - ministics.loc[i, "nbgraingel"]
        )

        # Calcul masse grains détruits par gel
        ministics.loc[i, "pgraingel"] = (
            ministics.loc[ind_Z69 : i - 1, "pgrain"]
            * ministics.loc[ind_Z69 + 1 : i, "nbgraingel"]
        ).sum()

    return ministics


# Description de la formation du rendement
# yield = masse et qualité des organes récoltés.
# Organes récoltés = organes reprod (grain = déhysdraté, fruit = hydraté), organes végé de stockage comme tige (sucre canne) ou racines (tubercules), ou tout (fourrage)

# 2 Méthodes possibles :

#   - Approche source/puits : --> PLANTES INDETERMINES
# Approche source / puit : chaque organe joue soit rôle de source de C, soit de puit de C, rôle qui peut changer au cours du cycle. Ex: feuille puit puis source
# Force source = nouveaux assimilats + ressource remob depuis org végé
# Force puit = taux de croiss potentiel --> fonction de l'âge de l'organe.
# Limites = force source dure à quantifier, nécesite choisir prio entre organes
# 'Solution' = imposer distribution cste des assimilat par phase phéno
# Indeterminate = espèce avec compétition trophique significative entre organes végétatifs et organes récoltés.

#   - Extension du Harvest Index = ratio grain biomasse / biomasse tot) à l'accumulation de biomasse dans grains --> PLANTES DETERMINEES
#       Idée : on multiplie la biomass totale produite par le harvest index qui varie linéairement avec le temps.
#       Avantage = mutualise 2 sources d'assimilats, moins de param à connaitre
#       Limite = il faut imposer un seuil max pr pas avoir des rembo irréalistes ou dépasser le rempliss max permis par nb de grains et masse max d'1 org.
#       Marche pr céréales, et plein d'autres : pois, vigne.


# Spécificité de STICS : calcul du nb d', et, etes récoltés AVANT l'allocation de biomasse à ces organes. Intérêt : joue le rôle de limite génétique du rendement --> masse fruit limitée à nb grain * poids max d'un grain
# Méthode de calcul du nb grain / fruit varie selon si plante déter / indéter


def yield_formation_determinate(ministics, CropParams, i, ind_Z69):
    """
    Calcul de la biomasse accumulée dans les grains par méthode carbon Harvest index
    """

    # Entre début rempliss grain et maturité physio = on calcule le HI puis la masse de fruit
    if (i >= ind_Z69) & (ministics.loc[i, "DOY_2Z89"] == 0):
        # Calcul du carbone HI entre 0 et 1
        # HI augmente linéairem avec temps entre début remplissage et maturation physio, pente et max fixés par biblio
        if CropParams.CODEIR == 1:
            ministics.loc[i, "ircarb"] = min(
                CropParams.VITIRCARB * (i - ind_Z69), CropParams.IRMAX
            )
        elif CropParams.CODEIR == 2:
            ministics.loc[i, "ircarb"] = min(
                CropParams.VITIRCARBT * (i - ind_Z69), CropParams.IRMAX
            )

        # Calcul de la biomasse accumulée journalière dans les grains avec prise en compte du stress thermique qui vient limiter le remplissage en carbone des grains
        # Intuition : on mutilplie le HI par la MS aérienne totale, et on soustrait ce même calcul au j-1 et j pour obtenir  dans les grains le jour j. HI doit être appliqué à MASEC et pas DELTAMS par définition.
        # Et on prend en compte le stress thermique impactant le remplissage.
        if (
            i == ind_Z69
        ):  # Le 1er jour de remplissage il n'y a pas de soustraction car c'est le premier jour où on ajoute de la masse aux fruits
            ministics.loc[i, "deltags"] = (
                ministics.loc[i, "ircarb"] * ministics.loc[i, "MASEC"]
            ) * ministics.loc[i, "ftempremp"]
        else:
            ministics.loc[i, "deltags"] = (
                ministics.loc[i, "ircarb"] * ministics.loc[i, "MASEC"]
                - ministics.loc[i - 1, "ircarb"] * ministics.loc[i - 1, "MASEC"]
            ) * ministics.loc[i, "ftempremp"]

        # Calcul de la biomasse de fruit accumulée depuis début remplissage des grains
        ministics.loc[i, "mafruit"] = min(
            (
                ministics.loc[ind_Z69:i, "deltags"]
                - ministics.loc[ind_Z69:i, "pgraingel"] / 100
            ).sum(),
            CropParams.PGRAINMAXI * ministics.loc[i, "nbgrains"],
        )  # garde min entre calculé et masse max théorique

        # Calcul de la biomasse moyenne par fruit
        ministics.loc[i, "pgrain"] = min(
            ministics.loc[i, "mafruit"] / ministics.loc[i, "nbgrains"] * 100,
            CropParams.PGRAINMAXI,
        )  # garde min entre calculé et masse max théorique par grain

    return ministics


def photoperiod(zlat,jday):
    ''' 
    Calcul de la photoP
    Exlication : nécessaire pr les exemples de STICS, mais quand on prend météo sur ERA5 on l'a déjà.
    '''

    maxjd = 365
    alat = zlat / 57.296
    zday = float(jday)
    days = [zday - 1.0, zday, zday + 1.0]
    pi = 3.14159

    if days[2] > float(maxjd):
        days[2] = 1.0
    if days[0] < 1.0:
        days[0] = float(maxjd)

    dl = [0.0, 0.0, 0.0]
    photp = [0.0, 0.0, 0.0]

    for i in range(3):
        theta1 = 2.0 * pi * (days[i] - 80.0) / 365.0
        theta2 = 0.034 * (math.sin(2.0 * pi * days[i] / 365.0) - math.sin(2.0 * pi * 80.0 / 365.0))
        theta = theta1 + theta2
        dec = math.asin(0.3978 * math.sin(theta))
        x1 = math.cos(alat) * math.cos(dec)
        b = -0.01454 / x1
        x2 = math.tan(alat) * math.tan(dec)
        h = math.acos(b - x2)
        daylen = 24.0 * h / pi
        dl[i] = daylen
        d = -1.0 * 0.10453 / (math.cos(alat) * math.cos(dec))
        p = d - (math.tan(alat) * math.tan(dec))
        p = math.acos(p)
        photp[i] = 24.0 * p / pi
        if photp[i] > 24.0:
            photp[i] = 24.0

    daylen = dl[1] # utilisé pr quoi dans STICS ?
    photp = photp[1]
    return photp

def meteo_df(meteo, station, starting_doy=295, user_doy_final=600):
    """
    Calcul des paramètres météo nécessaires pour STICS à partir des données météo issues d'ERA5.
    """
    year_2import = datetime.strptime("01/01/" + str(meteo.ANNEE.min()), "%d/%m/%Y")

    meteo_P = 101.3 * ((293 - 0.0065 * station.ZR) / 293) ** (5.26)
    meteo_gamma = 0.665 * 10 ** (-3) * meteo_P

    if 'latRad' not in meteo.columns: # cas où on utilise exemple STICS (qd ça vient d'ERA5 on l'a déjà)
        meteo['latRad'] = station.LATITUDE * np.pi / 180
    
    if 'Temp_moy' not in meteo.columns: # cas où on utilise exemple STICS (qd ça vient d'ERA5 on l'a déjà)
        meteo['Temp_moy'] = (meteo.Temp_min + meteo.Temp_max) / 2
    
    if 'Daily_Temp' not in meteo.columns: # cas où on utilise exemple STICS = on l'a pas mais ça pose prob après de pas avoir la col
        meteo['Daily_Temp'] = 0

    meteo["sol_dec"] = [fao.sol_dec(meteo.DOY[i]) for i in meteo.index]
    meteo["sunset_hour_angle"] = [
        fao.sunset_hour_angle(meteo.latRad[i], meteo.sol_dec[i]) for i in meteo.index
    ]
    meteo["inv_dist_earth_sun"] = [
        fao.inv_rel_dist_earth_sun(meteo.DOY[i]) for i in meteo.index
    ]

    meteo["rgex"] = [
        fao.et_rad(
            meteo.latRad[i],
            meteo.sol_dec[i],
            meteo.sunset_hour_angle[i],
            meteo.inv_dist_earth_sun[i],
        )
        * 0.75
        for i in meteo.index
    ]  # on multiplie par 0.75 pour passer des extra radiations au clear sky radiation - fao cs_rad
    meteo["Date"] = pd.to_datetime(meteo["Date"], format="%d/%m/%Y")
    meteo["esat"] = 0.6108 * np.exp(17.27 * meteo["Temp_moy"] / (meteo["Temp_moy"] + 237.3))
    if 'TPM' in meteo.columns: # quand météo vient d'un exemple STICS
        if meteo.TPM[0] == -999.9:
            meteo['TPM'] = (meteo['Temp_min'] - station.CORECTROSEE).apply(tvar) # voir paragraphe entre éq 9.17 et 9.18
    else: # météo vient de ERA5 donc il y a Rel_Hum --> d'où sort cette formule ? faire comme juste au-dessus sinon ?
        meteo["TPM"] = meteo["Rel_Hum"] / 100 * meteo["esat"]  # pression vapeur eau
    meteo["Delta"] = 4098 * meteo["esat"] / ((meteo["Temp_moy"] - 237.3) ** 2)

    if 'Photoperiod' not in meteo.columns: # cas où on utilise exemple STICS = faut calculer la photopériode
        meteo['Photoperiod'] = meteo.DOY.apply(lambda x:photoperiod(station.LATITUDE, x))


    if type(meteo.Daily_Temp[0]) is str:
        meteo["Daily_Temp"] = [ast.literal_eval(i) for i in meteo.Daily_Temp]
    
    if 'ETP' in meteo.columns: # exemple STICS
        meteo = meteo[
            [
                "Date",
                "Rain",
                "Radiation",
                "Temp_moy",
                "Photoperiod",
                "Temp_min",
                "Temp_max",
                "TPM",
                "rgex",
                "Wind",
                "esat",
                "Delta",
                "Daily_Temp", 'ETP'
            ]
        ]
    else: # meteo ERA5
        meteo = meteo[
            [
                "Date",
                "Rain",
                "Radiation",
                "Temp_moy",
                "Photoperiod",
                "Temp_min",
                "Temp_max",
                "TPM",
                "rgex",
                "Wind",
                "esat",
                "Delta",
                "Daily_Temp"
            ]
        ]


    # Rename Photoperiod-->PHOI
    meteo = meteo.rename(columns={'Photoperiod':'PHOI', 'Temp_moy':'Temp'})

    # Create a range of DOY from User.doy_sowing to User.doy_final
    doy_range = np.arange(starting_doy, user_doy_final + 1)

    # Create a DataFrame with DOY values
    ministics = pd.DataFrame({"doy": doy_range})

    # Add 'year_2import' to 'doy' to calculate 'D_Cal'
    ministics["D_Cal"] = year_2import + pd.to_timedelta(ministics["doy"] - 1, unit="d")

    # Merge 'ministics' with 'meteo' on 'D_Cal' and 'Date' columns, respectively
    ministics = pd.merge(ministics, meteo, left_on="D_Cal", right_on="Date")
    ministics = ministics.rename(columns={"Radiation": "trg"})

    return ministics, meteo_gamma


def phenology_dates(ministics, CropParams):
    '''
    Calcul des dates de changement de phase phénologique
    Explication : cette fonction associe à chaque date son stade phénologique (calculés dans la fonction phenology) et permet de générer des indices facilement utilisables dans d'autres fonctions.
    '''
    ind_EMERG = np.where(ministics["DOY_2EMERG"] > 0)[0][0]  # 1er jour d'émergence
    ind_AMF = np.where(ministics["DOY_2Z30"] > 0)[0][0]  # 1er jour AMF
    if CropParams.CODEINDETERMIN == 2:
        ind_NOU = np.where(ministics["DOY_NOU"] > 0)[0][0]  # jour de 'end of setting'
    ind_Z69 = np.where(ministics["DOY_2Z69"] > 0)[0][
        0
    ]  # 1er jour de remplissage des grains

    # Calcul de la somme UPVT DEPUIS DRP - formation rendement plantes indéterminées
    dico = {
        "Sum_UPVT_DRP": ministics["Sum_UPVT"] - CropParams.STLEVDRP
    }  # gdd cumulé entre DRP et maturation (PAS depuis début)
    ministics2 = pd.concat(dico.values(), axis=1, ignore_index=True)
    ministics2.columns = dico.keys()
    ministics = pd.concat([ministics, ministics2], axis=1)
    # ministics['Sum_UPVT_DRP'] = ministics['Sum_UPVT'] - CropParams.STLEVDRP # gdd cumulé entre DRP et maturation (PAS depuis début)
    ministics.loc[ministics.index < ind_Z69, "Sum_UPVT_DRP"] = 0

    # DELETE UNNECESSARY COLUMNS
    ministics = ministics.drop(["UPVT_GER_FLO", "Sum_UPVT_GER_FLO"], axis=1)
    # Colonne BBCH avec les num de phases phénologiques
    ministics["BBCH"] = -1.
    if CropParams.CODEPERENNE == 2:  # plantes pérennes
        ministics.loc[ministics["DOY_2DORMBREAK"] > 0, "BBCH"] = -0.5  # sortie dormance
        ministics.loc[ministics["DOY_2EMERG"] > 0, "BBCH"] = 0  # émergence
        ministics.loc[
            ministics["DOY_2Z30"] > 0, "BBCH"
        ] = 3  # élongation = après 1er noeud = feuilles se forment JE CONNAIS PAS BBCH
        ministics.loc[ministics["DOY_2Z65"] > 0, "BBCH"] = 6  # floraison
        ministics.loc[ministics["DOY_2Z69"] > 0, "BBCH"] = 7  # début remplissage grains
        if CropParams.CODEINDETERMIN == 2:
            ministics.loc[ministics["DOY_NOU"] > 0, "BBCH"] = 7.5  # end of setting
        ministics.loc[
            ministics["DOY_ILAX"] > 0, "BBCH"
        ] = 7.75  #  LAI max, jsp quel stade BBCH ça correspond
        ministics.loc[
            ministics["DOY_2Z89"] > 0, "BBCH"
        ] = 8.9  # maturation physiologique

    else:  # plantes annuelles
        ministics.loc[ministics["DOY_2EMERG"] > 0, "BBCH"] = 0  # émergence
        ministics.loc[
            ministics["DOY_2Z30"] > 0, "BBCH"
        ] = 3  # élongation = après 1er noeud = feuilles se forment = montaison
        ministics.loc[
            ministics["DOY_ILAX"] > 0, "BBCH"
        ] = 4  # montaison = après LAI max
        ministics.loc[ministics["DOY_2Z65"] > 0, "BBCH"] = 6  # floraison
        ministics.loc[
            ministics["DOY_2Z69"] > 0, "BBCH"
        ] = 7  # remplissage grains = développement des graines, stade aqueux puis laiteux
        if CropParams.CODEINDETERMIN == 2:
            ministics.loc[
                ministics["DOY_NOU"] > 0, "BBCH"
            ] = 7.5  # end of settingQ10JVC
        ministics.loc[
            ministics["DOY_2Z89"] > 0, "BBCH"
        ] = 8.9  # maturation physiologique

    if CropParams.CODEINDETERMIN == 2:
        return ministics, ind_Z69, ind_EMERG, ind_AMF, ind_NOU
    elif CropParams.CODEINDETERMIN == 1:
        return ministics, ind_Z69, ind_EMERG, ind_AMF


def phenology(ministics, i, CropParams, user):
    '''
    Calcul des stades phénologiques.

    Explications :
    Plantes herbacées / annuelles : unité de dev commencent à s'accumuler depuis émerg (ilev ds 3.3.2), en étant réduits par RFVI (entre 0 et 1) qui se rapproche de 1 avec accumulation de froid
    MAIS dans ministics on calcule un udevcult bis UDEVCULT_emerg depuis semis parce qu'on fait calcule date émergence avec ST2EMERGE, et pas avec la méthode stics
    Plantes ligneuses / pérennes : unité de dev ne commencent à s'accumuler que au moment dormancy break (DOY_2DORMBREAK = 1), donc avant on a RFVI=0
    --> traduit par rfvi qui vaut entre 0 et 1 pdt vernalisation des plantes annuelles, mais vaut 0 (donc annule udevcult) pour plantes ligneuses jusqu'à dorm break
    
    Simplifications : émergence EN DEVELOPPEMENT
    '''

    # Calcul de l'émergence / débourrement selon (différents formalismes si plante pérenne ou annuelle)

    if (
        CropParams.CODEPERENNE == 1
    ):  # si la plante est annuelle --> émergence quand seuil > ST2EMERGE (stics plus compliqué)
        # Calcul de UDEVCULT_emerg : parce que dans STICS udevcult doit être calculé à partir émerg, mais Benj s'en sert pour l'émergence
        ministics.loc[i, "UDEVCULT_emerg"] = np.where(
            ministics.loc[i, "tcult"] < CropParams.TDMIN,
            0,
            np.where(
                ministics.loc[i, "tcult"] < CropParams.TDMAX,
                ministics.loc[i, "tcult"] - CropParams.TDMIN,
                np.where(
                    ministics.loc[i, "tcult"] < CropParams.TCXSTOP,
                    (CropParams.TDMAX - CropParams.TDMIN)
                    * (ministics.loc[i, "tcult"] - CropParams.TCXSTOP)
                    / (CropParams.TDMAX - CropParams.TCXSTOP),
                    0,
                ),
            ),
        )

        # Calcul des gdd (gdd - somme degrés nécessaires jusqu'à émergence)
        # V2RIFIER QUE C ZRAC0 parce que dans stics ils parlent de profsemis
        ministics.loc[i, "Sum_DcDays"] = (
            -(CropParams.STPLTGER + CropParams.ST2EMERGE * CropParams.ZRAC0)
            + ministics.loc[0:i, "UDEVCULT_emerg"].cumsum()[i]
        )  # somme degrés jour - sommes théoriques pour atteindre émergence

        ministics.loc[i, "DOY_2GER"] = np.where(
            ministics.loc[i, "Sum_DcDays"]
            >= -(CropParams.ST2EMERGE * CropParams.ZRAC0),
            1,
            0,
        )  # jours après germination = 1, reste = 0
        ministics.loc[i, "DOY_2EMERG"] = np.where(
            ministics.loc[i, "Sum_DcDays"] >= 0, 1, 0
        )  # jours après émergence = 1, reste = 0

    elif (
        CropParams.CODEPERENNE == 2
    ):  # si la plante est pérenne --> dorm break calculé ou fixe puis débouremm calculé avec degrés heure. (choix de l'option codebfroid = 3 tt le temps)
        # Détermination de la date de sortie de dormance
        if CropParams.CODEDORMANCE == 1:  # date dorm break fixée = date de simulaton
            ministics.loc[i, "DOY_2DORMBREAK"] = 1  # 1 si après dorm break, 0 sinon

        else:  # date dorm break calculée (3.3.4.2), début dormance = 1er jour simulation.
            # Somme degrés jours avec Q10 et temp atmosphère (pas tcult)
            ministics.loc[i, "cu"] = sum(
                [
                    CropParams.Q10 ** (-ministics.loc[j, "Temp_max"] / 10)
                    + CropParams.Q10 ** (-ministics.loc[j, "Temp_min"] / 10)
                    for j in range(0, i + 1)
                ]
            )

            # Calcul de si dorm break dépasse ou non
            ministics.loc[i, "DOY_2DORMBREAK"] = np.where(
                ministics.loc[i, "cu"] > CropParams.JVC, 1, 0
            )

        # Calcul de la date de débourrement = degrés heure (on le fait après dorm break, et on ne le fait plus après émergence)
        if i != 0:
            if (ministics.loc[i, "DOY_2DORMBREAK"] == 1) and (
                ministics.loc[i - 1, "DOY_2EMERG"] == 0
            ):  # si on est après le dorm break et avant l'émergence
                # Calcul des T(h,n)

                thn = [
                    0
                    if t < CropParams.TDMINDEB
                    else (
                        CropParams.TDMAXDEB - CropParams.TDMINDEB
                        if t > CropParams.TDMAXDEB
                        else t - CropParams.TDMINDEB
                    )
                    for t in ministics.loc[i, "Daily_Temp"]
                ]
                # Calcul des GDH
                ministics.loc[i, "GDH"] = ministics.loc[i - 1, "GDH"] + sum(thn)

                ministics.loc[i, "DOY_2EMERG"] = int(
                    ministics.loc[i, "GDH"] > CropParams.STDORDEBOUR
                )  # vérif que ça marche

            elif (
                ministics.loc[i - 1, "DOY_2EMERG"] == 1
            ):  # après émergence, condition nécessaire pr pas que ça fasse ça avant dorm break
                ministics.loc[i, "DOY_2EMERG"] = 1

    # Calcul de la température effective de développement de la plante pour les autres stades phéno
    # 3.3.2 : unité de dev s'accumulent depuis émergence (annuelles) ou dorm break (pérennes)
    if (CropParams.CODEPERENNE == 1 and ministics.loc[i, "DOY_2EMERG"] == 1) or (
        CropParams.CODEPERENNE == 2 and ministics.loc[i, "DOY_2DORMBREAK"] == 1
    ):
        ministics.loc[i, "UDEVCULT"] = np.where(
            ministics.loc[i, "tcult"] < CropParams.TDMIN,
            0,
            np.where(
                ministics.loc[i, "tcult"] < CropParams.TDMAX,
                ministics.loc[i, "tcult"] - CropParams.TDMIN,
                np.where(
                    ministics.loc[i, "tcult"] < CropParams.TCXSTOP,
                    (CropParams.TDMAX - CropParams.TDMIN)
                    * (ministics.loc[i, "tcult"] - CropParams.TCXSTOP)
                    / (CropParams.TDMAX - CropParams.TCXSTOP),
                    0,
                ),
            ),
        )

    # Calculs de RFPI = réduction de dvlp liée à la photoP
    # S'applique sur période précise : Herbacée : entre émerg et DRP, pérenne : entre dorm break et DRP
    if CropParams.CODEPHOT:
        if i != 0:  # 1er jour pr le mom c'est RFPI = 1
            # entre émerg et DRP (annuelle), OU  entre dorm break et DRP (pérennes)
            # je mets sensiphot pour annuelles, parce que pérennes il est à -999 tt le temps MAIS PAS EXPLIQUE DANS STICS
            if (
                (CropParams.CODEPERENNE == 1)
                and (ministics.loc[i, "DOY_2EMERG"] == 1)
                and (ministics.loc[i - 1, "DOY_2Z69"] == 0)
            ):
                ministics.loc[i, "RFPI"] = (1 - CropParams.SENSIPHOT) * (
                    ministics.loc[i, "PHOI"] - CropParams.PHOSAT
                ) / (CropParams.PHOSAT - CropParams.PHOBASE) + 1
                # selon plante,
            elif (
                (CropParams.CODEPERENNE == 2)
                and (ministics.loc[i, "DOY_2DORMBREAK"] == 1)
                and (ministics.loc[i - 1, "DOY_2Z69"] == 0)
            ):
                ministics.loc[i, "RFPI"] = (
                    ministics.loc[i, "PHOI"] - CropParams.PHOSAT
                ) / (CropParams.PHOSAT - CropParams.PHOBASE) + 1
            else:
                ministics.loc[i, "RFPI"] = 1
        else:
            ministics.loc[i, "RFPI"] = 1
    else:
        ministics.loc[i, "RFPI"] = 1

    # Calcul de RFVI =  réduction de dvlp liée à la vernalisation.
    if (
        CropParams.CODEPERENNE == 1
    ):  # si la plante est annuelle = RFVI entre 0 et 1 pdt vernalisation
        if CropParams.CODEBFROID:  # si la plante a besoin d'une vernalisation
            if ministics.loc[i, "DOY_2GER"] == 1:
                # on somme jvi depuis jour 0 mais on les calcule que à partir germ donc en fait ça revient à sommer depuis Germ
                ministics.loc[i, "JVI"] = max(
                    (
                        1
                        - (
                            (CropParams.TFROID - ministics.loc[i, "tcult"])
                            / CropParams.AMPFROID
                        )
                        ** 2
                    ),
                    0,
                )  # = entre 0 et 1 = contribution vernalisation de ce jour = proche 1 si temp proche temp optimale de vernalisation
                ministics.loc[i, "RFVI"] = (
                    max(
                        0,
                        (ministics.loc[0:i, "JVI"].cumsum()[i] - CropParams.JVCMINI)
                        / (CropParams.JVC - CropParams.JVCMINI),
                    )
                    if (ministics.loc[0:i, "JVI"].cumsum()[i] - CropParams.JVCMINI)
                    / (CropParams.JVC - CropParams.JVCMINI)
                    < 1
                    else 1
                )  # avancée de la vernalisation = nb jour vernalisant subis / nb de jours nécessaire en théorie. 1 si on a déjà atteint le nb de jours nécessaires. Max(0,value) pas dans doc mais dans code source javastics.

        else:  # si la plante n'a pas besoin de vernalisation
            ministics.loc[i, "RFVI"] = 1

    elif (
        CropParams.CODEPERENNE == 2
    ):  # si la plante est pérenne = RFVI vaut 0 jusqu'à dorm break puis 1 après
        if (
            CropParams.CODEBFROID
        ):  # si la plante a besoin d'une vernalisation. RAJOUTE APRES A VERIF
            if ministics.loc[i, "DOY_2DORMBREAK"] == 0:  # si on est avant le dorm break
                ministics.loc[i, "RFVI"] = 0
            else:  # si on est après le dorm break
                ministics.loc[i, "RFVI"] = 1
        else:
            ministics.loc[i, "RFVI"] = 1

    # on peut conditionner tout ce qui est en-dessous par DOY_2EMERG == 1 non ?

    # Calcul des dates de chaque stade phénologique avec degré jour contraint par RFPI et RFVI : 1 col = 1 stade avec 1 si la plante est à ce stade et 0 sinon.
    # DOY_2EMERG prend en cpte si plante pérenne ou annuelle = 2 calculs diffs jusqu'à émergence, et mtn on peut faire la mm chose pr pér ou annuelle.
    # RFVI vaut 1 pr plantes pér après émerg, mais annuelle c'est peut être possible que mm après émerg il soit tjrs <1 d'où le fait que benj l'a mis là ? Et pour FLO MAT, pq que

    # ancien à garder au cas où
    # ministics.loc[i,'UPVT_GER_FLO'] = ministics.loc[i,'DOY_2EMERG'] * ministics.loc[i,'UDEVCULT'] * ministics.loc[i,'RFPI'] * ministics.loc[i,'RFVI'] # si après émerg, T°plante * rfpi * rfvi, sinon 0.
    # ministics.loc[i,'Sum_UPVT_GER_FLO'] = ministics.loc[0:i,'UPVT_GER_FLO'].cumsum()[i] # gdd après germ, contraint par RFPI et RFVI

    ministics.loc[i, "UPVT_post_EMERG"] = (
        ministics.loc[i, "DOY_2EMERG"]
        * ministics.loc[i, "UDEVCULT"]
        * ministics.loc[i, "RFPI"]
        * ministics.loc[i, "RFVI"]
    )  # RFPI prend déjà en compte le stade où on est
    ministics.loc[i, "Sum_UPVT_post_EMERG"] = ministics.loc[
        0:i, "UPVT_post_EMERG"
    ].cumsum()[
        i
    ]  # gdd après émerg, contraint par RFPI et RFVI

    ministics.loc[i, "DOY_2Z30"] = np.where(
        ministics.loc[i, "Sum_UPVT_post_EMERG"] >= CropParams.STLEVAMF, 1, 0
    )  # 1 si gdd cumulé >= stade 1st node, sinon 0.
    ministics.loc[i, "DOY_ILAX"] = np.where(
        ministics.loc[i, "Sum_UPVT_post_EMERG"]
        >= CropParams.STLEVAMF + CropParams.STAMFLAX,
        1,
        0,
    )  # 1 si gdd cumulé >= stade flag_leaf, sinon 0.
    ministics.loc[i, "DOY_2Z65"] = np.where(
        ministics.loc[i, "Sum_UPVT_post_EMERG"]
        >= CropParams.STLEVDRP - CropParams.STFLODRP,
        1,
        0,
    )  # 1 si gdd cumulé >= stade floraison, sinon 0.

    # ancien à garder au cas où mais --> pour moi, dans flo_mat faut faire* 'DOY_2EMERG' parce que tous les seuils partent de là. Benj fait * 'DOY_2Z65' et du coup son calcul de 'DOY_2Z89' me semble pas logique
    # ministics.loc[i,'UPVT_FLO_MAT'] = ministics.loc[i,'DOY_2Z65'] * ministics.loc[i,'UDEVCULT'] * ministics.loc[i,'RFVI'] # si après floraison/émerg, T° plante * RFVI. Pq pas RFPI ?
    # ministics.loc[i,'Sum_UPVT_FLO_MAT'] = ministics.loc[0:i,'UPVT_FLO_MAT'].cumsum()[i] # gdd entre flo et mat, contraint par RFVI

    ministics.loc[i, "DOY_2Z69"] = np.where(
        ministics.loc[i, "Sum_UPVT_post_EMERG"] >= CropParams.STLEVDRP, 1, 0
    )  # 1 si gdd (cumulé et contraint) >= stade début rempliss grains, sinon 0. (pr module formaiton rendement)

    # J'ai l'impression que STDRPMAT est tjrs à -999 pr plantes indéterminées, donc on le calcule que pour déter.
    if CropParams.CODEINDETERMIN == 1:
        ministics.loc[i, "DOY_2Z89"] = np.where(
            ministics.loc[i, "Sum_UPVT_post_EMERG"]
            >= CropParams.STLEVDRP + CropParams.STDRPMAT,
            1,
            0,
        )  # 1 si gdd (cumulé et contraint) >= stade maturation, sinon 0.
    else:  # D'après 3.2.2, ILAX = début maturation pour la vigne (et attention ILAX après DRP pour vigne).
        # PAC a confirmé c'est ok ILMAT = ILAX.
        ministics.loc[i, "DOY_2Z89"] = ministics.loc[i, "DOY_ILAX"]

    # Pour Sum_UPVT on regarde d'abord quel jour on est avant de le calculer.
    if ministics.loc[i, "DOY_2EMERG"] == 0:
        ministics.loc[i, "Sum_UPVT"] = ministics.loc[i, "Sum_DcDays"]
    else:
        ministics.loc[i, "Sum_UPVT"] = ministics.loc[i, "Sum_UPVT_post_EMERG"]
    # ancien à garder au cas où
    # elif ministics.loc[i,'DOY_2Z65'] == 0:
    #     ministics.loc[i,'Sum_UPVT'] = ministics.loc[i,'Sum_UPVT_GER_FLO']
    # elif (ministics.loc[i,'DOY_2Z65'] > 0) and ((ministics.loc[i-1,'DOY_2Z65'] == 0)): # le jour de floraison --> on change ce qu'on avait mis la veille. Ca va pas changer ce qu'on a fait la veille...
    #     ministics.loc[i-1,'Sum_UPVT'] = ministics.loc[i-1,'Sum_UPVT_GER_FLO'] # on l'a bien changé, mais ça ne change pas les calculs réalisés à i-1...
    #     ministics.loc[i,'Sum_UPVT'] = ministics.loc[i,'Sum_UPVT_FLO_MAT']
    # elif ministics.loc[i,'DOY_2Z65'] > 0:
    #     ministics.loc[i,'Sum_UPVT'] = ministics.loc[i,'Sum_UPVT_FLO_MAT']
    # col Sum_UPVT :
    #   - pour les jours avant émergence : gdd cumulé depuis semis - somme théorique pour atteindre émergence
    #   - pour les jours entre émerg et floraison (et jour juste avant floraison): gdd cumulé contraint par RFPI et RFVI entre émerg et floraison (PAS DEPUIS LE DEBUT)
    #   - pour les jours après floraison (et veille du 1er jour floraison) : gdd cumulé contraint par RFVI entre émerg et maturation (PAS DEPUIS LE DEBUT)

    # Calcul stade 'end of setting' INOU - plantes indéter
    if CropParams.CODEINDETERMIN == 2:
        ministics.loc[i, "DOY_NOU"] = np.where(
            ministics.loc[i, "Sum_UPVT"] >= CropParams.STLEVDRP + CropParams.STDRPNOU,
            1,
            0,
        )  # 1 si gdd (cumulé et contraint) >= stade end of setting, sinon 0.

    return ministics


def initialise_vars(ministics, soil, CropParams, station):
    '''
    Initialisation des param qui vont être calculés chaque jour
    '''
    # Pq initialiser en nan --> on ne sait pas si > NAN donne vrai ou faux c'est trompeur.
    ministics[
        "PFZ1"
    ] = np.nan  # PFZ1 = 1 si point de flétrissement horiz 1 non atteint, 0 sinon.
    ministics[
        "PFZ2"
    ] = np.nan  # PFZ2 = 1 si point de flétrissement horiz 2 non atteint, 0 sinon.
    ministics[
        "DeltaZEff"
    ] = np.nan  # valeur d'approfondeissement du front racinaire ce jour_ci
    ministics[
        "DeltaZ"
    ] = (
        np.nan
    )  # min entre valeur d'approfondissement du front racinaire DeltaZEff et ce qu'il restait entre deltaZ j-1 et prof max du sol (ProfSoil)
    ministics["ZRac"] = np.nan  # prof max atteinte par système racinaire (cm)

    ministics["ULAI"] = np.nan  # = unité de dvp relatif du LAI (entre 0 et 3)
    ministics[
        "DELTAI_dev"
    ] = np.nan  # composante de dvp phasique du DELTAI (m2/plant/gdd)
    ministics[
        "DELTAI"
    ] = np.nan  # augmentation journélière d'index de feuilles par surf de sol (m2/m2/j)
    ministics["DELTAI_stress"] = np.nan  # composante stress du DELTAI
    ministics["DELTAI_sen"] = np.nan
    ministics["DELTAIunstr_sen"] = np.nan
    ministics["LAI"] = np.nan
    ministics["LAI_unstressed"] = np.nan

    ministics["EOS"] = np.nan  # flux évap max journalier (mm/j)
    ministics["FlagRain"] = np.nan
    ministics["SEOS"] = np.nan
    ministics["FlagPhase"] = np.nan
    ministics["SES"] = np.nan
    ministics["ES"] = np.nan
    ministics["ES_a"] = np.nan
    ministics["ESol_sum"] = np.nan
    ministics["HUR1"] = np.nan

    ministics[
        "ZDemi"
    ] = (
        np.nan
    )  # prof racinaire assurant extraction proche surface de 20% de l'eau disponible
    ministics["CUMLRACZ"] = np.nan
    ministics["TETA1"] = np.nan
    ministics["TETA2"] = np.nan
    ministics["EP_sum"] = np.nan
    ministics[
        "HUR1"
    ] = (
        np.nan
    )  # contenu en eau de l'horiz 1. C'est par cm de sol j'ai l'impression parce qu'on multiplie par prof sol à chaque fois, et on divise par ça pour maj HUR1
    ministics[
        "HUR2"
    ] = np.nan  # contenu en eau de l'horiz 2. C'est par cm de sol j'ai l'impression

    ministics["DELTAMS"] = np.nan
    ministics["DELTAMS_unstressed"] = np.nan
    ministics["MASEC"] = np.nan
    ministics["MASEC_unstressed"] = np.nan

    ministics["DRAIN1"] = np.nan
    ministics["HUR1"] = soil.HUR1_Init
    ministics["DRAIN2"] = np.nan
    ministics["HUR2"] = soil.HUR2_Init

    colonnes_nulles = [
        'EO','EOP','TETA','TETA_moy','EP','EP_1_a','EP_2_a','SWFAC','EP_1_th','EP_2_th',
        "MAXLAI_sofar",
        "MAXLAIunstr_sofar",
        "nbgrains",
        "nbgraingel",
        "pgrain",
        "pgraingel",
        "ircarb",
        "deltags",
        "mafruit",
        "raint",
        "raint_unstressed",
        "ebmax",
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
        "durvie",
        "somsen",
        "dltaisen",
        "resperenne",
        "resperennemax",
        "resperennestruc",
        "maperenne",
        "tetsen",
        "senfac",
        "turfac",
        "fstressgel",
        "dltafv",
        "dltat",
        "dltares",
        "dltaremobsen",
        "masecveg",
        "restempmax",
        "dltarestemp",
        "senstress",
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
        "fracinsol",
        "Rglo",
        "rnet",
        "et",
        "hauteur",
        "z0",
        "tcultmax",
        "tcult",
        "DOY_2EMERG",
        "DOY_2Z30",
        "DOY_NOU",
        "DOY_2Z69",
        "UDEVCULT",
        "UDEVCULT_emerg",
        "RFPI",
        "JVI",
        "Sum_DcDays",
        "DOY_2GER",
        "JVI_cum",
        "RFVI",
        "UPVT_GER_FLO",
        "Sum_UPVT_GER_FLO",
        "DOY_ILAX",
        "DOY_2Z65",
        "DOY_2Z89",
        "UPVT_post_EMERG",
        "Sum_UPVT_post_EMERG",
        "Sum_UPVT",
        "ebmax",
        "TETSTOMATE",
        "TETURG",
        "GDH",
        "DOY_2DORMBREAK",
        "cu",
        "DELTAIDENS",
        "DELTAIT",
        "esol_a",
        "esol",
        "stemflow",
        "mouill",
        "fapar",
        "rnetP1",
        "rnetS",
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
        "ePT",
        "dos",
        "Emd",
        "rnetP2",
        "EdirectM",
        "Edirect",
        'dltamstombe','mafeuiltombe'
    ]

    if station.CODEETP != 1: # ETP doit être calculé
        colonnes_nulles.append('ETP')

    dico = {}
    for col in colonnes_nulles:
        dico[col] = pd.Series([0. for i in range(len(ministics))])

    ministics2 = pd.concat(dico.values(), axis=1, ignore_index=True)
    ministics2.columns = dico.keys()

    ministics = pd.concat([ministics, ministics2], axis=1)

    ministics[
        "ratioTF"
    ] = (
        CropParams.TIGEFEUIL
    )  # sinon ça reste à nan quand on n'est pas en période de racourcissement des jours

    ministics["fgelflo"] = np.nan
    ministics["ftemp"] = np.nan
    ministics["ftempremp"] = np.nan
    ministics["BBCH"] = np.nan  # ne pas mettre 0 car 0 = émergence
    ministics['BBCH'] = ministics['BBCH'].astype('float')

    for K in range(CropParams.NBOITE):
        ministics[f"fpft{K}"] = 0
        ministics[f"potcroifruit{K}"] = 0
        ministics[f"nfruit{K}"] = 0
        ministics[f"croifruit{K}"] = 0

    ministics["converge"] = False

    # ind_Z65 = np.where(ministics['DOY_2Z65'] > 0)[0][0] # 1er jour du stade floraison
    # ind_Z89 = np.where(ministics['DOY_2Z89'] > 0)[0][0] # 1er jour du stade maturation
    # Slope_LAI_Sen = ministics.loc[ind_Z89, 'Sum_UPVT'] - ministics.loc[ind_Z65, 'Sum_UPVT'] # Différence entre gdd cumulé au 1er jour maturation (=entre flo et matura) et 1er jour floraison (=gdd ind_Z65-1 + gdd ind_Z65)

    return ministics


def initialize_lrac_matrix(ministics, soil, CropParams):
    '''
    Initialisation de la matrice des densités racinaires par profondeur
    '''
    LRAC_MATRIX = np.empty(
        (ministics.shape[0], soil.ProfSoil)
    )  # profil de densité racinaire pr chaque prof sol : taille =  nb jours simulation * prof sol
    LRAC_MATRIX[:] = np.nan

    if (
        CropParams.CODEPERENNE == 2
    ):  # si plante pérenne : je prends h1 et h2 de vigne_ini (au lieu de 5 faudra peut être adapter)
        # Question en suspens : sert à rien d'initialiser car LRAC du jour j+1 ne dépend pas de LRAC j-1 mais que de ZRAC. Ca va être écrasé tout de suite par le calcul de LRAC au jour 0.
        LRAC_MATRIX[0, 0 : soil.EPC_1] = 0.0044
        LRAC_MATRIX[0, soil.EPC_1 : soil.EPC_1 + soil.EPC_2] = 0.0400

    return LRAC_MATRIX


def root_growth(ministics, CropParams, SoilParams, i):
    """
    Calcul de la croissance racinaire journalière et de la nouvelle profondeur max (ZRac)
    Explications : calcule si on est au pt de flétrissement, dans ce cas = pas de croissance racinaire, sinon croissance journalière = cste * (temp plante - Tmin), avec max de croissance qui ne doit pas dépasser la prof du sol.
    
    Simplification :
        - pas d'effet du stress hydrique sur croissance racinaire, mais la teneur en eau du sol peut limiter la croiss racinaire. ????
        - critères d'arrêt de croiss non implémentés ici : si plante pérenne ET réserves simulées ET 'date = (date of LAI = 0)' (pas compris) OU STOPRAC atteint (crop specific) ou OBSTARAC atteint( soil specific)
    """

    if (
        i == 0
    ):  # 1er jour = on intialise des valeurs, pq ne pas les intialiser dans la cellule précèd ?
        ministics.loc[i, "PFZ1"] = 1
        ministics.loc[i, "PFZ2"] = 1
        ministics.loc[i, "DeltaZEff"] = 0
        ministics.loc[i, "DeltaZ"] = 0
        ministics.loc[i, "ZRac"] = CropParams.ZRAC0
    else:
        # Calcul du dépassement ou non du pt de flétrissement
        ministics.loc[i, "PFZ1"] = (
            1 if ministics.loc[i - 1, "HUR1"] >= SoilParams.HMINF_1 else 0
        )  # PFZ1 = 1 si point de flétrissement horiz 1 non atteint, 0 sinon.
        ministics.loc[i, "PFZ2"] = (
            1 if ministics.loc[i - 1, "HUR2"] >= SoilParams.HMINF_2 else 0
        )  # PFZ2 = 1 si point de flétrissement horiz 2 non atteint, 0 sinon.

        # PAs de croiss racinaire avant germination pr plantes annuelles, et avant dorm break pr plantes pérennes
        if (
            ministics.loc[i, "DOY_2GER"] == 0
            and ministics.loc[i, "DOY_2Z89"] == 0
            and CropParams.CODEPERENNE == 1
        ) or (
            ministics.loc[i, "DOY_2DORMBREAK"] == 0
            and ministics.loc[i, "DOY_2Z89"] == 0
            and CropParams.CODEPERENNE == 2
        ):  # jours avant germ et avant maturation (????)
            ministics.loc[i, "DeltaZEff"] = 0
            ministics.loc[i, "DeltaZ"] = 0  # pas de croiss racinaire
            ministics.loc[i, "ZRac"] = CropParams.ZRAC0  # ZRAC reste ini

        # jours entre germ et floraison pr annuelles (Benj), et entre dorm break et floraison pr pérennes (impro, STICS pas précis)
        elif (
            ministics.loc[i, "DOY_2GER"] == 1
            and CropParams.CODEPERENNE == 1
            and ministics.loc[i, "DOY_2Z65"] == 0
        ) or (
            ministics.loc[i, "DOY_2DORMBREAK"] == 1
            and CropParams.CODEPERENNE == 2
            and ministics.loc[i, "DOY_2Z65"] == 0
        ):
            # Selon si prof max racinaire j-1 dépasse horiz 1 ou pas, on utilise PFZ1 ou PFZ2
            # Retourne 0 si point de flétrissement atteint = les racines ne peuvent pas grandir sans eau, sinon : taux de croiss racinaire * [min(temp plante j-1, temp optimale de croiss) - temp min croiss]
            if ministics.loc[i, "tcult"] <= CropParams.TCMIN:
                ministics.loc[i, "DeltaZEff"] = 0
            else:  # Eq 5.1 et 5.2 (ici on fait direct deltaZ_T * deltaz_stress)
                ministics.loc[i, "DeltaZEff"] = (
                    CropParams.CROIRAC
                    * (
                        max(
                            0,
                            min(ministics.loc[i, "tcult"], CropParams.TCMAX)
                            - CropParams.TCMIN,
                        )
                    )
                    * ministics.loc[i, "PFZ1"]
                    if ministics.loc[i - 1, "ZRac"] < SoilParams.EPC_1
                    else CropParams.CROIRAC
                    * (
                        min(ministics.loc[i, "tcult"], CropParams.TCMAX)
                        - CropParams.TCMIN
                    )
                    * ministics.loc[i, "PFZ2"]
                )

            ministics.loc[i, "DeltaZ"] = min(
                ministics.loc[i, "DeltaZEff"],
                SoilParams.ProfSoil - ministics.loc[i - 1, "ZRac"],
            )  # min entre DeltaZEff calculé et (prof max du sol - prof max racinaire) = pour ne pas dépasser la prof max du sol
            ministics.loc[i, "ZRac"] = (
                ministics.loc[i, "DeltaZ"] + ministics.loc[i - 1, "ZRac"]
            )  # ajout de la croiss à la prof j-1

        # jours après floraison = pas de croissance racinaire --> STICS dit qu'il y a un stade à partir duquel root growth stops, mais il en propose plsrs dont floraison, mat, ilax; Benj a choisi ça.
        elif ministics.loc[i, "DOY_2Z65"] == 1:
            ministics.loc[i, "DeltaZEff"] = 0
            ministics.loc[i, "DeltaZ"] = 0
            ministics.loc[i, "ZRac"] = ministics.loc[i - 1, "ZRac"].copy()

    return ministics


def root_density(ministics, CropParams, SoilParams, Constants, LRAC_MATRIX, i):
    """
    Calcul de la densité racinaire : profil de densité racinaire + biomasse racinaire morte au mom récolte

    Explications :
    Choix : 2 méthodes de calcul dans STICS, miniSTICS = méthode du 'profil standard'
    Principe : détermine profil racinaire efficace en termes d'absorption --> densité racinaire en fonction profondeur = suit sigmoïde, 1 courbe sigmoïde pour chaque prof max ZRAC. Voir Figure 5.4
    Hyp miniSTICS sur ZRAC0 : il y a des racines au-delà de ZRAC0 (mm au mom émergence).
    """
    # Calcul de la prof racinaire nécessaire pour absorber 20% de l'eau (proche surface)
    ministics.loc[i, "ZDemi"] = max(
        ministics.loc[i, "ZRac"] - CropParams.ZPRLIM + CropParams.ZPENTE,
        (np.log(4) / CropParams.S),
    )  # 1.4 / S = seuil pr assurer minimum extraciton de 20% de l'eau proche surface

    ZRac_Ceil = np.ceil(ministics.loc[i, "ZRac"])

    if (
        ministics.loc[i, "DOY_2GER"] == 1 or CropParams.CODEPERENNE == 2
    ):  # si après germ pr plantes annuelles, ou tout le temps pour plantes pérennes (ça overwrite l'ini du coup nan pr pérennes?)
        for c in range(
            int(ZRac_Ceil)
        ):  # pour chaque prof jusqu'à prof max --> on remplit la matrice de densité racinaire pour chaque profondeur
            LRAC_MATRIX[i, c] = Constants.LVOPT / (
                1 + np.exp(-CropParams.S * (c - ministics.loc[i, "ZDemi"]))
            )  # LVOPT = palier de la sigmoïde = max de densité

    # Calcul de somme des densités racinaires de surface jusqu'à fin horizon 2 le jour i
    # en cm/cm3
    ministics.loc[i, "CUMLRACZ"] = np.nansum(
        LRAC_MATRIX[i, 0 : SoilParams.EPC_1 + SoilParams.EPC_2]
    )
    # ministics.loc[i+1,'CUMLRACZ'] = ministics.loc[i,'CUMLRACZ'].copy()

    return LRAC_MATRIX


def leaf_area_index(ministics, CropParams, i):
    '''
    Calcul du LAI

    Explications :
    Principe : LAI augmente de émergence à LAI max, puis plateau jusqu'à floraison, puis décroissance linéaire des gdd jusqu'à maturation
    Phase d'agumentaiton jusqu'au plateau : on calcule l'augmentation journalière de LAI = DELTAI qui suit une sigmoïde avec en abscisse l'avancement relatif entre émergence et stade flag leaf
    Phase de diminution entre floraison et mauration : on calcule le LAI MAX rencontré avant et on soustrait uen proportion de ce max corrélée à l'avancement relatif entre floraison et maturation
    Prend en compte la compétition interplantes par la densité du semis et un facteur de compétition fixé en input
    Prend en compte le stress hydrique stomatique (calcul plus tard) entr e0 et 1 qui vient diminuer le LAI.
    Amélioration :
        - ne calculer le MAX-LAI que si on est après floraison
    '''

    # si stade après emerg et avant flag leaf : ULAI = 1.x ou 2.x = quantifie (entre 0 et 3) l'avancement dans le stade entre émergence et flag leaf.
    # En gros : à émerg -> ULAI = 1, à stade LAI max -> ULAI = 3, reste des stades -> ULAI=0

    if (
        ministics.loc[i, "DOY_2EMERG"] == 1 and ministics.loc[i, "DOY_ILAX"] == 0
    ):  # après ULAI = 0 donc sert à rien de calculer deltai ?
        ministics.loc[i, "ULAI"] = np.where(
            ministics.loc[i, "Sum_UPVT"]
            < CropParams.STLEVAMF,  # si entre émerg et 1st node --> returne 1.x avec x = à quel % de la transition émerg-1st node le jour i est
            1 + (2 - 1) * (ministics.loc[i, "Sum_UPVT"] / CropParams.STLEVAMF),
            np.where(
                ministics.loc[i, "Sum_UPVT"]
                < CropParams.STLEVAMF
                + CropParams.STAMFLAX,  # si entre 1st node et flag leaf --> retourne 2.x avec x = à quel % de la transition 1st node-fla gleaf le jour i est
                2
                + (3 - 2)
                * (ministics.loc[i, "Sum_UPVT"] - CropParams.STLEVAMF)
                / CropParams.STAMFLAX,
                0,
            ),
        )

        # Calcul de l'augmentation de LAI du jour i en m2/pl/gdd
        # DELTAI_dev = fonction sigmoïde de augmentation LAI (DELTAI en m2/pl/gdd) en fonction ULAI. Plateau = DLAIMAXBRUT
        # PENTLAIMAx = pente de la sigmoïde. VLAIMAX = valeur de ULAI au point d'inflexion qui fait chuter le plateau
        # Intuition : entre émerg et flag leaf, l'augmentation de LAI en m2/pl/gdd augmente selon l'avancée vers le stade flag leaf en suivant une sigmoïde. Ici, on calcule où en est le jour i sur cette courbe.
        # Pas de udlaimax car il vaut 3 dans tous les xml.
        ministics.loc[i, "DELTAI_dev"] = np.where(
            ministics.loc[i, "ULAI"] > 0,
            CropParams.DLAIMAXBRUT
            / (
                1
                + np.exp(
                    CropParams.PENTLAIMAX
                    * (CropParams.VLAIMAX - ministics.loc[i, "ULAI"])
                )
            ),
            0,
        )

        # Calcul de l'effet de la composante gdd du DELTAI.
        # On utilise pas directement UDEVCULT car ce n'est pas exactement pareil : ici on utilise tcmin/tcmax, au lieu de tdmin/tdmax
        ministics.loc[i, "DELTAIT"] = np.where(
            ministics.loc[i, "tcult"] < CropParams.TCMIN,
            0,
            np.where(
                ministics.loc[i, "tcult"] < CropParams.TCMAX,
                ministics.loc[i, "tcult"] - CropParams.TCMIN,
                np.where(
                    ministics.loc[i, "tcult"] < CropParams.TCXSTOP,
                    (CropParams.TCMAX - CropParams.TCMIN)
                    * (ministics.loc[i, "tcult"] - CropParams.TCXSTOP)
                    / (CropParams.TCMAX - CropParams.TCXSTOP),
                    0,
                ),
            ),
        )

        # Calcul de l'effet de la compétition interplantes sur le deltai
        ministics.loc[i, "DELTAIDENS"] = CropParams.DENSITESEM * (
            CropParams.DENSITESEM / CropParams.BDENS
        ) ** (
            CropParams.ADENS
        )  # à mettre en ini car constant au cours tps

        # Calcul de l'augmentation de LAI du jour i en m2/pl en prenant en compte COMPETITION
        # DELTAI en m2 feuille / m2 sol / j = composante phasique * gdd (donne qqch en m2 pl) * densité * 1/(rapport densité / densité min pour compét)^fac de compét -
        # Intuition : on multiplie le DELTAI_dev par les gdd et le nb de plante par m2 de sol = vla en m2 pl /m2 de sol, et on réduit l'augmentation obtenue selon la  compétition en déoulant et on soustrait la Tbase
        # Eq 4.1 (deltai_dev * deltai_dens juste)
        ministics.loc[i, "DELTAI"] = (
            ministics.loc[i, "DELTAI_dev"]
            * ministics.loc[i, "DELTAIDENS"]
            * ministics.loc[i, "DELTAIT"]
        )

        # DELTAI_stress = indice de stress hydrique stomatique (entre 0 et 1) du jour précédent, sauf 1er jour = 1
        ministics.loc[i, "DELTAI_stress"] = (
            1 if i == 0 else ministics.loc[i - 1, "SWFAC"]
        )

        # Prise en compte de l'impact négatif de la photoPériode durant la période de racourcissement des jours.
        # Idée = RFPI vient diminuer l'augmentation du LAI
        if i > 0:
            if ministics.loc[i, "PHOI"] < ministics.loc[i - 1, "PHOI"]:
                ministics.loc[i, "DELTAI"] = (
                    ministics.loc[i, "DELTAI"] * ministics.loc[i, "RFPI"]
                )

                ministics.loc[i, "ratioTF"] = CropParams.TIGEFEUIL

        # Pas de croissance si la photoP est inf à la photoB base + on fait diminuer la part de tige dans la masse tige+feuille avec rfpi pq ??
        # Si photoP < PHOBASE -> rfpi = 0 et croiss LAI = 0
        if ministics.loc[i, "PHOI"] < CropParams.PHOBASE:
            ministics.loc[i, "RFPI"] = 0
            ministics.loc[i, "DELTAI"] = (
                ministics.loc[i, "DELTAI"] * ministics.loc[i, "RFPI"]
            )
            ministics.loc[i, "ratioTF"] = (
                CropParams.TIGEFEUIL * ministics.loc[i, "RFPI"]
            )  # je comprends pas l'équation 7.8

        # Calcul du LAI (stressé) max rencontré entre jour 1 et jour i. Sachant que LAI jour i = somme des DELTAI depuis jour 1 pondérés par stress hydrique stomatique
        # Remarque : le LAI ne peut faire que augmenter si on ajoute qqch chaque jour (mm si tout petit), sauf si DELTAI peut être négatif
        if i > 0:
            ministics.loc[i, "MAXLAI_sofar"] = (
                sum(
                    ministics.loc[0:i, "DELTAI"] * ministics.loc[0:i, "DELTAI_stress"]
                )  # somme DELTAI entre jour 1 et i pondérés par leur stress hydrique stomatique journaliers si cette sommme > MAXLAI_sofar, sinon = MAXLAI_sofar
                if sum(
                    ministics.loc[0:i, "DELTAI"] * ministics.loc[0:i, "DELTAI_stress"]
                )
                > ministics.loc[i - 1, "MAXLAI_sofar"]
                else ministics.loc[i - 1, "MAXLAI_sofar"]
            )

            # Calcul du LAI max rencontré sans prise en compte du stress hydrique stomatique.
            ministics.loc[i, "MAXLAIunstr_sofar"] = (
                sum(ministics.loc[0:i, "DELTAI"])
                if sum(ministics.loc[0:i, "DELTAI"])
                > ministics.loc[i - 1, "MAXLAIunstr_sofar"]
                else ministics.loc[i - 1, "MAXLAIunstr_sofar"]
            )

        # Si entre floraison et maturité : (différence entre ggd jour i et gdd jour floraison) / gdd nécessaire entre floraison et maturité * max LAI stressé rencontré / , sinon 0
        # DELTAI_sen = durant transition floraison-maturité, ça va de 0 à MAXLAI_sofar en avançant linéairement selon l'avancement relatif (en %) dans cette transiton selon les gdd requis
        # Pourquoi on avance de cette façon entre floraison et maturité ?

        if ministics.loc[i, "DOY_2Z65"] == 1 and ministics.loc[i, "DOY_2Z89"] == 0:
            ind_Z65 = np.where(ministics["DOY_2Z65"] > 0)[0][
                0
            ]  # 1er jour du stade floraison

            ministics.loc[i, "DELTAI_sen"] = ministics.loc[i, "MAXLAI_sofar"] * min(
                1,
                (ministics.loc[i, "Sum_UPVT"] - ministics.loc[ind_Z65, "Sum_UPVT"])
                / CropParams.STDRPMAT,
            )  # LAI stressé.
            # ministics.loc[i,'DELTAI_sen'] = ministics.loc[i,'MAXLAI_sofar'] * (ministics.loc[i,'Sum_UPVT']-ministics.loc[ind_Z65,'Sum_UPVT']) / CropParams.STDRPMAT # LAI stressé.
            #   Je triche avec le max car ça dépasse un tout petit 1 et donc ça donne un deltai négatif.
            ministics.loc[i, "DELTAIunstr_sen"] = (
                ministics.loc[i, "MAXLAIunstr_sofar"]
                * (ministics.loc[i, "Sum_UPVT"] - ministics.loc[ind_Z65, "Sum_UPVT"])
                / CropParams.STDRPMAT
            )  # Pareil avec LAI max non stressé

        else:
            ministics.loc[i, "DELTAI_sen"] = 0
            ministics.loc[i, "DELTAIunstr_sen"] = 0

    # Si avant émergence ou après stade LAI max (=flag leaf)
    else:
        ministics.loc[i, "ULAI"] = 0
        ministics.loc[i, "DELTAI_dev"] = 0
        ministics.loc[i, "DELTAI"] = 0
        ministics.loc[i, "DELTAI_stress"] = 0
        ministics.loc[i, "DELTAI_sen"] = 0
        ministics.loc[i, "DELTAIunstr_sen"] = 0

    # Calcul du LAI : somme des DELTAI entre jour 1 et i pondérés par stress hydrique stomatique - baisse LAi liée à la sénescence des feuilles durant la transition floraison-maturité, à la maturité on retire le max de LAI
    ministics.loc[i, "LAI"] = (
        sum(ministics.loc[0:i, "DELTAI"] * ministics.loc[0:i, "DELTAI_stress"])
        - ministics.loc[i, "DELTAI_sen"]
        if ministics.loc[i, "DOY_2Z89"] == 0
        else 0
    )

    # Calcul du LAI en omettant le stress hydrique stomatique
    ministics.loc[i, "LAI_unstressed"] = (
        sum(ministics.loc[0:i, "DELTAI"]) - ministics.loc[i, "DELTAIunstr_sen"]
        if ministics.loc[i, "DOY_2Z89"] == 0
        else 0
    )

    return ministics


def tvar(x):
    """
   Calcul de la pression vap saturante - température (équation 9.19 de la Doc)
    """
    return 6.1070 * (1 + 2 ** (1 / 2) * np.sin(0.017453293 * x / 3)) ** 8.827


def foliage_water(ministics, CropParams, i):
    """
    Calul de l'eau retenue par feuilles après pluie ou irrigation de la canopée
    Remarque : tous les XML = pas cette modélisation = pas d'eau sur feuilles. Important en zones tropicales plutôt ptêtre
    """
    if i == 0:  # pas de feuille que ce soit pér ou annuelle
        ministics.loc[
            i, "stemflow"
        ] = 0  # valeur qd LAI nul mais je comprends pas  la logique
    else:
        # Calcul de l'eau s'écoulant le long des tigesfoliage_water_energy_budget
        ministics.loc[i, "stemflow"] = (
            CropParams.STEMFLOWMAX
            * ministics.loc[i, "Rain"]
            * (1 - np.exp(-CropParams.KSTEMFLOW * ministics.loc[i - 1, "LAI"]))
        )

    # Calcul de l'eau retenue par les feuilles : eau non écoulée * facteur de mouillabilité
    # Calcul pas écrit dans STICS, y'a écrit dépend LAI ? ici stemflow dépend lai
    ministics.loc[i, "mouill"] = (
        ministics.loc[i, "Rain"] - ministics.loc[i, "stemflow"]
    ) * CropParams.MOUILLABIL

    return ministics


def energy_budget(ministics, CropParams, i, constants, user, soil):
    """
    Calcul du bilan d'énergie avec la resistive approach.

    Explications : 
    Flux d'énergie liées à la radiation incidente
    Approche résistance prenant en compte tous les flux d'énergie :
        1. Soil evap (EOS)
        2. Plant transpi (EOP)
        3. Evap directe de l'eau interceptée par feuilles (Emd)
        4. 2 types de résistance :
            - résist à la diffusion --> entre canopée/sol (ras) et entre 'cover' (=haut canopée) / niveau ref météo (raa)
            - résist de surface de la canopée (rc) et 'couche limite de canopée' (rac)
    Tous les flux = flux réels avec même formule sauf plant transpi = le flux maximal potentiel
    Diffusion = fait que vapeur d'eau se déplace dans l'air humide pour homogénéiser la compisition du milieu.
    Résistance à diffusion = fait qu'un matériau laisse passer la vapeur pour se diffusion. Résist faible = vapeur peut fuire le matériau.
    """

    ### Energie disponible et distribution plante/sol ###
    # pas calculée encore : Emd

    # Calcul la fraction de PAR interceptée par plante
    ministics.loc[i, "fapar"] = ministics.loc[i, "raint"] / (
        constants.PARSURRG * ministics.loc[i, "trg"]
    )

    # Calcul de l'énergie dispo pour la plante = rnet * fapar * 0.83 = équivalent du coeff d'extinction pr la net radiation
    # = énergie avant évap eau sur plante et avant transpi
    ministics.loc[i, "rnetP1"] = (
        0.83 * ministics.loc[i, "fapar"] * ministics.loc[i, "rnet"]
    )  # eq 9.26

    # Calcul de l'énergie dispo pour le sol = rnet - énergie interceptée/dispo pour la plante
    ministics.loc[i, "rnetS"] = ministics.loc[i, "rnet"] - ministics.loc[i, "rnetP1"]

    ### Résistances à la diffusion ###

    # Pour lai(t) = 0
    # Vérif que Wind en m.s-1
    ministics.loc[i, "ras0"] = (
        np.log(user.ZR / soil.Z0SOLNU)
        * np.log((ministics.loc[i, "z0"] + ministics.loc[i, "dh"]) / soil.Z0SOLNU)
        / (0.41**2 * ministics.loc[i, "Wind"])
    )
    ministics.loc[i, "raa0"] = (
        np.log(user.ZR / soil.Z0SOLNU) ** 2
        / (0.41**2 * ministics.loc[i, "Wind"])
        - ministics.loc[i, "ras0"]
    )

    # Pour lai(t) >= 4
    ministics.loc[i, "rasinf"] = np.log(
        (user.ZR - ministics.loc[i, "dh"]) / ministics.loc[i, "z0"]
    ) / (0.41**2 * ministics.loc[i, "Wind"]) - ministics.loc[i, "hauteur"] / (
        2.5 * (ministics.loc[i, "hauteur"] - ministics.loc[i, "dh"])
    ) * (
        12.18
        - np.exp(
            2.5
            * (
                1
                - (ministics.loc[i, "dh"] + ministics.loc[i, "z0"])
                / ministics.loc[i, "hauteur"]
            )
        )
    )
    ministics.loc[i, "raainf"] = (
        np.log((user.ZR - ministics.loc[i, "dh"]) / ministics.loc[i, "z0"])
        / (0.41**2 * ministics.loc[i, "Wind"])
        * (
            np.log((user.ZR - ministics.loc[i, "dh"]) / ministics.loc[i, "hauteur"])
            + ministics.loc[i, "hauteur"]
            / (2.5 * (ministics.loc[i, "hauteur"] - ministics.loc[i, "dh"]))
            * np.exp(
                2.5
                * (
                    1
                    - (ministics.loc[i, "dh"] + ministics.loc[i, "z0"])
                    / ministics.loc[i, "hauteur"]
                )
                - 1
            )
        )
    )

    # Pour lai entre 0 et 4raa
    if i == 0:
        ministics.loc[i, "ras"] = (4 - 0) * ministics.loc[i, "ras0"]
        ministics.loc[i, "raa"] = (4 - 0) * ministics.loc[i, "raa0"]
    else:
        ministics.loc[i, "ras"] = 0.25 * (
            ministics.loc[i - 1, "LAI"] * ministics.loc[i, "rasinf"]
            + (4 - ministics.loc[i - 1, "LAI"]) * ministics.loc[i, "ras0"]
        )
        ministics.loc[i, "raa"] = 0.25 * (
            ministics.loc[i - 1, "LAI"] * ministics.loc[i, "raainf"]
            + (4 - ministics.loc[i - 1, "LAI"]) * ministics.loc[i, "raa0"]
        )

    ### Autres param nécessaires ###

    # Calcul de la pente de la relation pression vapeur sat - température ###
    # tair en °C    # print('i',i)
    # print('raainf',ministics.loc[i,'raainf'])
    # if isinf(ministics.loc[i,'raainf']):
    #         print('fraction 1',np.log((user.ZR - ministics.loc[i,'dh'])/ministics.loc[i,'z0']) / (0.41**2 * ministics.loc[i,'Wind']))
    #         print('ln 2',np.log((user.ZR - ministics.loc[i,'dh'])/ministics.loc[i,'hauteur']))
    #         print('3eme terme',ministics.loc[i,'hauteur'] / (2.5 * (ministics.loc[i,'hauteur'] - ministics.loc[i,'dh'])) * np.exp(2.5*(1-(ministics.loc[i,'dh'] + ministics.loc[i,'z0'])/ministics.loc[i,'hauteur'])-1))
    #         raise ValueError('raainf nan')
    # Différence entre tair et tmoy dans stics ??
    ministics.loc[i, "deltat"] = tvar(ministics.loc[i, "Temp"] + 0.5) - tvar(
        ministics.loc[i, "Temp"] - 0.5
    )

    # Calcul de la chaleur latente de vaporisation (MJ.kg-1) - eq 9.37
    ministics.loc[i, "L"] = (2500840 - 2358.6 * ministics.loc[i, "Temp"]) * 10e-6

    # Calcul du déficit de saturation de l'air
    ministics.loc[i, "dsat"] = tvar(ministics.loc[i, "Temp"]) - ministics.loc[i, "TPM"]

    ### Résistance de la surface ###

    if i == 0:
        pass
    elif (
        ministics.loc[i - 1, "LAI"] != 0
    ):  # Quand LAI est nul (notamment 1er jour) = pas de calcul, et pas de transpi
        # Calcul de la résistance de la couche limite de la canopée
        # Seuil min : rac = 12.5 sm-1
        ministics.loc[i, "rac"] = min(12.5, 50 / (2 * ministics.loc[i - 1, "LAI"]))

        # Calcul de la résistance de surface de la canopée
        ministics.loc[i, "rc"] = (
            CropParams.RSMIN
            * (0.5 * ministics.loc[i - 1, "LAI"] + 1)
            / ministics.loc[i - 1, "LAI"]
            * (0.039 * ministics.loc[i, "dsat"] + 0.45)
            * 28
            / (2.5 + ministics.loc[i, "trg"])
            * CropParams.FCO2S
        )

    return ministics


def intermed_sat_deficit(ministics, i, meteo_gamma):
    '''
    Calcul de l'évapotranspiration potentielle totale et du déficit de saturation.
    '''

    # Calcul de l'évapotranspiration potentielle totale (sol + canopée). Nécessaire pour calculer dos
    # Condition = 'soil conditions kept moist' ==> pas en cas de séchereesse ? mais du coup on pourrait pas calculer dos
    # Formule basée sur rnetS pourquoi ?
    ministics.loc[i, "ePT"] = (
        1.32
        * ministics.loc[i, "rnetS"]
        * ministics.loc[i, "deltat"]
        / (ministics.loc[i, "deltat"] + meteo_gamma)
    )

    # Calcul du déficit de sat(meteo_gamma*ministics.loc[i,'rnet'] - (ministics.loc[i,'deltat'] - meteo_gamma)* ministics.loc[i,'ePT'])uration dans la canopée
    # Définition = différence entre épaisseur de saturation et épaisseur dsatréelle des précipiations. Epaisseur de satiration = épaisseur à laquelle les précip devraient commencer pr 1 qtté donnée d'humidité dans l'atmosphère.
    # Basé sur eq. 11.5 mais y'a d'autre * ministics.loc[i,'raa'] /105.03s formules eq 9.25 et 9.37 où ça dépend de EVAPO qu'on ne peut pas calculer sans dos... --> iterative process ?
    ministics.loc[i, "dos"] = (
        ministics.loc[i, "dsat"]
        + (
            meteo_gamma * ministics.loc[i, "rnet"]
            - (ministics.loc[i, "deltat"] + meteo_gamma) * ministics.loc[i, "ePT"]
        )
        * ministics.loc[i, "raa"]
        / 105.03
    )

    return ministics


def foliage_water_evap_energy_budget(ministics, i, meteo_gamma, CropParams):
    """
    Calcul de l'évap de l'eau retenue par les feuilles, limitée par mouillePT
    """
    if (i == 0) or (CropParams.CODEINTERCEPT == 2):
        ministics.loc[i, "Emd"] = 0
    elif ministics.loc[i - 1, "LAI"] == 0:
        ministics.loc[i, "Emd"] = 0
    else:
        ministics.loc[i, "Emd"] = max(
            ministics.loc[i, "mouill"],
            (
                ministics.loc[i, "deltat"] * ministics.loc[i, "rnetP1"]
                + 105.03 * ministics.loc[i, "dos"] / ministics.loc[i, "rac"]
            )
            / (ministics.loc[i, "L"] * (ministics.loc[i, "deltat"] + meteo_gamma)),
        )

    return ministics


def potential_transpi_energy_budget(ministics, i, meteo_gamma):
    '''
    Calcul de la transpiration potentielle avec l'approche bilan d'énergie
    '''
    if i == 0:
        ministics.loc[i, "EOP"] = 0
    elif ministics.loc[i - 1, "LAI"] == 0:
        ministics.loc[
            i, "EOP"
        ] = 0  # C'est pas qu'il n'y a pas d'énergie dispo pr transpi, mais que rac n'est pas calculabele
    else:
        # Calcul de l'énergie dispo pour la transpi de la plante (= énergie dispo pour la plante - ce qui est utilisé pr évap eau sur plante)
        ministics.loc[i, "rnetP2"] = (
            ministics.loc[i, "rnetP1"] - ministics.loc[i, "Emd"]
        )

        # Calcul de la transpi plante maximale (potentiel ce flux --> réel va dépendre de eau dispo dans sol)
        # on l'appelle eop alors qu'on parle encore du potentiel. Dans Benj / miro j'avais mis eo = potentielle et eop = réel. Stics : eop = transpi max
        ministics.loc[i, "EOP"] = (
            ministics.loc[i, "deltat"] * ministics.loc[i, "rnetP2"]
            + 105.03 * ministics.loc[i, "dos"] / ministics.loc[i, "rac"]
        ) / (
            ministics.loc[i, "L"]
            * (
                ministics.loc[i, "deltat"]
                + meteo_gamma * (1 + ministics.loc[i, "rc"] / ministics.loc[i, "rac"])
            )
        )
        ministics.loc[i, "EOP"] = max(ministics.loc[i, "EOP"], 0) # Max evaporation can only be positive

    return ministics


def potential_soil_evap_energy_budget(ministics, i, meteo_gamma):
    '''
    Calcul de l'évaporation du sol potentielle avec l'approche bilan d'énergie (équation 11.3)
    '''

    ministics.loc[i, "EOS"] = (
        ministics.loc[i, "deltat"] * ministics.loc[i, "rnetS"]
        + 105.03 * ministics.loc[i, "dos"] / ministics.loc[i, "ras"]
    ) / (ministics.loc[i, "L"] * (ministics.loc[i, "deltat"] + meteo_gamma))
    return ministics


def potential_soil_evap_beer_law(ministics, CropParams, i):
    '''
    Calcul de l'évaporation du sol potentielle avec l'approche loi de Beer
    '''
    if i == 0:
        ministics.loc[i, "EOS"] = 0
    else:
        ministics.loc[i, "EOS"] = ministics.loc[i, "ETP"] * np.exp(
            -(CropParams.EXTIN - 0.2) * ministics.loc[i - 1, "LAI"]
        )  # eq 11.1

    return ministics


def foliage_water_evap_lai(ministics, i, CropParams):
    """
    Calcul de l'évaporation de l'eau retenue par les feuilles
    Eq 11.10 : jsp ce qu'est DELTA je l'ai pas misEO
    """

    if CropParams.CODEINTERCEPT == 2:
        ministics.loc[i, "Emd"] = 0
    elif CropParams.CODEINTERCEPT == 1:
        ministics.loc[i, "Emd"] = min(
            ministics.loc[i, "mouill"],
            ministics.loc[i, "EO"]
            - ministics.loc[i, "ETP"]
            * np.exp(-(CropParams.EXTIN - 0.2) * ministics.loc[i - 1, "LAI"]),
        )

    return ministics


def potential_transpi_crop_coef(ministics, i, CropParams):
    """
    Calcul de la transpiration potentielle de la plante avec l'approche Coefficient de culture.
    """
    if i == 0:
        ministics.loc[i, "EO"] = 0
        ministics.loc[i, "EOP"] = 0
    else:
        # Calcul de transpiration potentielle (ou besoin en eau) sans limite en eau dans sol ou plante
        # EO = ETP * sigmoïde du LAI. Palier  = KMAX = conductance stomatique max atteint quand LAI = 5, valeur spécifique à chaque plante.
        ministics.loc[i, "EO"] = ministics.loc[i, "ETP"] * (
            1
            + (CropParams.KMAX - 1)
            / (1 + np.exp(-1.5 * (ministics.loc[i - 1, "LAI"] - 3)))
        )  # eq 11.8

        # Calcul de l'évap de l'eau retenue par le feuillage
        ministics = foliage_water_evap_lai(ministics, i, CropParams)

        # Calcul de la transpi max de la plante en soustrayant composante liée au sol et prenant cpte manque d'eau
        # Eq 11.12 : Edirectm = EOS, et Edirect = esol_a d'après moi. Et parenthèse manquante d'après Benj je fais comme lui.
        ministics.loc[i, "EdirectM"] = ministics.loc[i, "EOS"] + ministics.loc[i, "Emd"]
        ministics.loc[i, "Edirect"] = ministics.loc[i, "esol"]
        ministics.loc[i, "EOP"] = 0 if ministics.loc[i, "ETP"] == 0 else (
            ministics.loc[i, "EO"] - ministics.loc[i, "EdirectM"]
        ) * (1.4 - 0.4 * (ministics.loc[i, "Edirect"] / ministics.loc[i, "EdirectM"]))
        # je rajoute condition avec etp sinon dénominateur nul

        ministics.loc[i, "EOP"] = max(ministics.loc[i, "EOP"], 0) # Max evaporation can only be positive

    return ministics


def actual_transpi(ministics, i, CropParams, SoilParams, LRAC_MATRIX):
    """'
    Calcul de la transpiration réelle.

    Explications :
    Principe = la plante peut transpirer (=prendre eau sol) à un taux max jusqu'à ce que contenu eau sol passe sous un seuil.
    Seuil dépend root density, paramètre stomate, demande évaporation.
    Rapport réel/max = ep/eop correspond à un stress hydrique stomatique SWFAC
    Cette fonction = exactement ce qu'a fait Benj car pas détaillé dans STICS.
    """

    # 1er jour, avant émergence, et après maturation. Et aussi si EOP = 0 = LAI nul donc pas de transpi
    if (i == 0) or (ministics.loc[i, "EOP"] == 0):
        ministics.loc[i, "EP"] = 0
        ministics.loc[i, "TETA"] = 0
        ministics.loc[i, "TETA_moy"] = 0
        ministics.loc[i, "EP_1_th"] = 0
        ministics.loc[i, "EP_2_th"] = 0
        ministics.loc[i, "EP_1_a"] = 0
        ministics.loc[i, "EP_2_a"] = 0
        ministics.loc[i, "SWFAC"] = 1

    elif (ministics.loc[i - 1, "DOY_2EMERG"] == 1) & (
        ministics.loc[i - 1, "DOY_2Z89"] == 0
    ):
        # Calcul du seuil de contenu en eau en-dessous duquel la transpi réelle < transpi max
        ministics.loc[i, "TETSTOMATE"] = (
            1
            / 80
            * np.log(
                ministics.loc[i, "EOP"]
                / (
                    2
                    * np.pi
                    * ministics.loc[i - 1, "CUMLRACZ"]
                    * CropParams.PSISTO
                    * 0.0001
                )
                * np.log(
                    1
                    / (
                        CropParams.RAYON
                        * (
                            ministics.loc[i - 1, "CUMLRACZ"]
                            / ministics.loc[i - 1, "ZRac"]
                        )
                        ** (1 / 2)
                    )
                )
            )
        )

        # /!\ At this point the concept is clearly different form miniSTICS/STICS concepts. coment ça ??

        # Calcul de l'eau dispo dans les cm de sol où il y a bien des racines, pour chaque xo
        # max(eau dispo * nb de cm visités par racines , 0)
        ministics.loc[i, "TETA1"] = max(
            (ministics.loc[i, "HUR1"] - SoilParams.HMINF_1)
            * sum(np.isnan(LRAC_MATRIX[i - 1, 0 : SoilParams.EPC_1]) == 0),
            0,
        )
        ministics.loc[i, "TETA2"] = (
            max(
                (SoilParams.HUR2_Init - SoilParams.HMINF_2)
                * sum(
                    np.isnan(
                        LRAC_MATRIX[
                            i - 1,
                            SoilParams.EPC_1 : SoilParams.EPC_1
                            + SoilParams.EPC_2,
                        ]
                    )
                    == 0
                ),
                0,
            )
            if i == 0
            else max(
                (ministics.loc[i - 1, "HUR2"] - SoilParams.HMINF_2)
                * sum(
                    np.isnan(
                        LRAC_MATRIX[
                            i - 1,
                            SoilParams.EPC_1 : SoilParams.EPC_1
                            + SoilParams.EPC_2,
                        ]
                    )
                    == 0
                ),
                0,
            )
        )

        # Calcul de somme totale des contenus en eau dans les cm de cols visités par les racines
        ministics.loc[i, "TETA"] = ministics.loc[i, "TETA1"] + ministics.loc[i, "TETA2"]

        # Calcul de moyenne / cm de sol visité par racine des contenus en eau sur les 2 horizons
        if (
            sum(
                pd.isna(
                    LRAC_MATRIX[
                        i - 1, 0 : (SoilParams.EPC_1 + SoilParams.EPC_2)
                    ]
                )
                == 0
            )
            == 0
        ):
            ministics.loc[i, "TETA_moy"] = 0
        else:
            ministics.loc[i, "TETA_moy"] = (
                ministics.loc[i, "TETA1"] + ministics.loc[i, "TETA2"]
            ) / sum(
                pd.isna(
                    LRAC_MATRIX[
                        i - 1, 0 : (SoilParams.EPC_1 + SoilParams.EPC_2)
                    ]
                )
                == 0
            )

        # Calcul de l'indice de stress hydrique stomatique = basé sur le contenu en eau moyen du sol
        # 1 si le contenu en eau moyen est au-dessus du seuil de fermeture des stomates
        # sinon : contenu en eau moyen / seuil de fermeture des stomates
        ministics.loc[i, "SWFAC"] = (
            1
            if ministics.loc[i, "TETA_moy"] > ministics.loc[i, "TETSTOMATE"]
            else ministics.loc[i, "TETA_moy"] / ministics.loc[i, "TETSTOMATE"]
        )

        # Calcul de la transpi en prenant cpte stress stomatique
        # Intuition : le stress hydrique stomatique correspond au rapport transpiration réelle / transpiration potentielle/max
        ministics.loc[i, "EP"] = ministics.loc[i, "SWFAC"] * ministics.loc[i, "EOP"]

        # Calcul de la répartition de la transpi par horizon 1 et 2 = proportionnel au nombre de racines dans chaque horizon
        # = transpi réelle * somme densités racinaires horiz X / somme des densités racinaires sur les 2 horiz
        ministics.loc[i, "EP_1_th"] = (
            ministics.loc[i, "EP"]
            * np.nansum(LRAC_MATRIX[i - 1, 0 : SoilParams.EPC_1])
            / ministics.loc[i, "CUMLRACZ"]
        )
        ministics.loc[i, "EP_2_th"] = (
            ministics.loc[i, "EP"]
            * np.nansum(
                LRAC_MATRIX[
                    i - 1,
                    (SoilParams.EPC_1 + 1) : (
                        SoilParams.EPC_1 + SoilParams.EPC_2
                    ),
                ]
            )
            / ministics.loc[i, "CUMLRACZ"]
        )

        # Calcul de la transpi réelle par horizon prenant en compte limite eau sol
        # contenu en eau < pt flétriss -> EP_a = 0
        # contenu en eau > pt flétriss -> EP_a = min(eau dispo horiz x, transpi horiz x)
        # Question : pq décalage i-1 que pour hoziron 2 ?
        ministics.loc[i, "EP_1_a"] = (
            0
            if ministics.loc[i, "HUR1"] < SoilParams.HMINF_1
            else min(
                (ministics.loc[i, "HUR1"] - SoilParams.HMINF_1)
                * 1000
                * SoilParams.EPC_1
                / 100,
                ministics.loc[i, "EP_1_th"],
            )
        )
        if i == 0:
            ministics.loc[i, "EP_2_a"] = (
                0
                if ministics.loc[i - 1, "HUR2"] < SoilParams.HMINF_2
                else min(
                    (SoilParams.HUR2_Init - SoilParams.HMINF_2)
                    * 1000
                    * SoilParams.EPC_2
                    / 100,
                    ministics.loc[i, "EP"] - ministics.loc[i, "EP_1_a"],
                )
            )
        else:
            ministics.loc[i, "EP_2_a"] = (
                0
                if ministics.loc[i - 1, "HUR2"] < SoilParams.HMINF_2
                else min(
                    (ministics.loc[i - 1, "HUR2"] - SoilParams.HMINF_2)
                    * 1000
                    * SoilParams.EPC_2
                    / 100,
                    ministics.loc[i, "EP"] - ministics.loc[i, "EP_1_a"],
                )
            )

    # Calcul de la transpi réelle en prenant en compte limite eau sol. En gros si eau a effectivement limité transpi, on réduit EP, sinon rien ne change.
    ministics.loc[i, "EP"] = (
        ministics.loc[i, "EP_1_a"] + ministics.loc[i, "EP_2_a"]
        if (ministics.loc[i, "EP_1_a"] + ministics.loc[i, "EP_2_a"])
        < ministics.loc[i, "EP"]
        else ministics.loc[i, "EP"]
    )

    # Calcul transpi réelle (limite eau prise en compte) cumulée entre jour 1 et jour i
    ministics.loc[i, "EP_sum"] = sum(ministics["EP"][0:i])

    # Calcul de teneur en eau de horiz 1 : teneur en eau horiz 1 - eau transpirée par cm de sol
    ministics.loc[i, "HUR1"] = ministics.loc[i, "HUR1"] - (
        (ministics.loc[i, "EP_1_a"] * 0.001) / (SoilParams.EPC_1 / 100)
    )
    # Calcul de teneur en eau de horiz 2 : teneur en eau horiz 2 veille - eau transpirée par cm de sol
    # pq i-1 ?
    ministics.loc[i, "HUR2"] = (
        SoilParams.HUR2_Init
        - ((ministics.loc[i, "EP_2_a"] * 0.001) / (SoilParams.EPC_2 / 100))
        if i == 0
        else ministics.loc[i - 1, "HUR2"]
        - ((ministics.loc[i, "EP_2_a"] * 0.001) / (SoilParams.EPC_2 / 100))
    )

    return ministics


def actual_soil_evap(ministics, i, SoilParams, user, station):
    """
    Calcul de l'évaporation du sol réelle

    Explications :
    Principe du calcul de l'évaporation réelle :
       - 1ère phase : evap réelle = evap maximum
       - 2ème phase = pas de pluie depuis plsrs jours et evap maximum cumulée dépasse seuil Q0 : evap sol ralentit
    Q0 = quantité d'eau évaporée du sol (sans nouvel évènement pluie) tq l'évaporation va ralentir (moins d'eau à évaporer)
    C'est bien base sur un 'évènement pluie' (Benj=STICS) et pas une quantité tombée, je vois pas trop pq, on pourrait calculer l'eau dans le sol suite à la pluie
    """

    # Calcul de l'évap du sol maximum 'cumulée' SEOS : s'il pleut un jour, SEOS = EOS, mais tant qu'il ne pleut pas SEOS = accumulation des EOS des jours d'avant et jour i
    # ATTENTION 'cumulé' = pas cumulé depuis début : c'est cumulé depuis le dernier jour de pluie donc pas cumulé pr 1 jour de pluie
    ministics.loc[i, "FlagRain"] = (
        1 if ministics.loc[i, "Rain"] > 0 else 0
    )  # 1 = évènement pluie, 0 = pas de pluie

    if i == 0:  # 1er jour : SEOS = EOS
        ministics.loc[i, "SEOS"] = ministics.loc[i, "EOS"]
    else:
        if ministics.loc[i, "FlagRain"] == 1:
            ministics.loc[i, "SEOS"] = ministics.loc[
                i, "EOS"
            ]  # si évènement pluie -> SEOS = EOS
        else:
            ministics.loc[i, "SEOS"] = (
                ministics.loc[i, "EOS"] + ministics.loc[i - 1, "SEOS"]
            )  # s'il n'a pas plu jour i -> SEOS = EOS + SEOS jour i-1

    # Calcul passage phase 1->2 évap = quand SEOS cumulé dépasse seuil Q0.
    ministics.loc[i, "FlagPhase"] = (
        1 if ministics.loc[i, "SEOS"] <= SoilParams.Q0 else 2
    )

    # Calcul de AEVAP
    SoilParams.AEVAP = ( 0.5 * station.ACLIM * ((0.63 - SoilParams.HA) ** (5 / 3)) * (SoilParams.HCCF_1 - SoilParams.HA))  # (aevap) param d'évap du sol combinant aspects climatiques et sol (mm). Vent * -humidité * (capacité au champ - humidité)

    # Calcul de l'évap du sol 'cumulée' - prend en cpte selon si on est phase 1 ou 2 = réduction de l'évap
    # eq 11.6 : pas de somme car SEOS déjà cumulé ici pour les jours sans pluie.
    # ATTENTION 'cumulé' : c'est cumulé depuis le dernier jour de pluie, avec réduction prise en compte
    if ministics.loc[i, "FlagPhase"] == 1:  # phase 1 --> evap à son max = SEOS
        ministics.loc[i, "SES"] = ministics.loc[i, "SEOS"]
    else:  # phase 2 --> cumulé = total potentiel phase 1 (Q0) + potentiel évaporé pdt phase 2 (SEOS - Q0) réduit par formule
        ministics.loc[i, "SES"] = (
            SoilParams.Q0
            + (
                2 * SoilParams.AEVAP * (ministics.loc[i, "SEOS"] - SoilParams.Q0)
                + (SoilParams.AEVAP**2)
            )
            ** 0.5
            - SoilParams.AEVAP
        )

    # Calcul de l'évap selon pluie/pas pluie et phase1/2
    # j'ai changé les conditions de Benj que je ne comprenais pas, voir avec Roxane. J'ai remis ancien ca rmachait pas
    if i == 0:
        ministics.loc[i, "esol"] = ministics.loc[i, "SES"]

    else:
        if (
            ministics.loc[i, "FlagRain"] == 1 and ministics.loc[i, "FlagPhase"] == 1
        ):  # s'il a plu et phase 1 -> ES = SEOS = EOS (SEOS est journalier déjà)
            ministics.loc[i, "esol"] = ministics.loc[i, "SES"]
        elif (
            ministics.loc[i, "FlagRain"] == 0 and ministics.loc[i, "FlagPhase"] == 1
        ):  # s'il n'a pas plu et phase 1 -> ES = SES - SES j-1 = on passe de cumulé à journalier
            ministics.loc[i, "esol"] = (
                ministics.loc[i, "SES"] - ministics.loc[i - 1, "SES"]
            )
        elif (
            ministics.loc[i - 1, "FlagPhase"] == 1
            and ministics.loc[i, "FlagPhase"] == 2
        ):  # si hier on avait pas atteint Q0 mais auj si, -> ES = SES - Q0 : c'est rare ça veut dire qu'en 1 jour EOS > Q0 -> evap jour i > Q0
            ministics.loc[i, "esol"] = (
                ministics.loc[i, "SES"] - SoilParams.Q0
            )  # pq on soustrait q0 ?
        else:  # s'il n'a pas plu et phase 2 -> ES = SES - SES j-1 = on passe de cumulé à journalier
            ministics.loc[i, "esol"] = (
                ministics.loc[i, "SES"] - ministics.loc[i - 1, "SES"]
            )

    # Calcul de l'évap réelle journalière prenant en cpte la limite en eau du sol
    # Intuition HUR1 - HA = eau vraiment disponible pr la plante, parce que HA = humidité qu'il y a dans un sol qu'on considère 'sec'
    # Méthode 1 (flag_HA true) = on utilise l'humid relative pr limiter evap : min (contenu eau - humid résid, évap journalière)
    # Méthode 2 (flag_HA false) = on utilise le pt de flétriss donné dans biblio : min(contenu eau - pt flétriss, évap journalière)
    if (
        i == 0
    ):  # 1er jour simulation = hyp pas de valeur négative positive pour HUR1 -HA ou HUMIN1
        ministics.loc[i, "esol_a"] = (
            min(
                (SoilParams.HUR1_Init - SoilParams.HA) * SoilParams.EPC_1 * 10,
                ministics.loc[i, "esol"],
            )
            if user.FLAG_HA
            else min(
                (SoilParams.HUR1_Init - SoilParams.HMINF_1)
                * SoilParams.EPC_1
                * 10,
                ministics.loc[i, "esol"],
            )
        )
    else:
        if user.FLAG_HA:  # Méthode 1
            if (ministics.loc[i - 1, "HUR1"] - SoilParams.HA) > 0:
                ministics.loc[i, "esol_a"] = min(
                    (ministics.loc[i - 1, "HUR1"] - SoilParams.HA)
                    * SoilParams.EPC_1
                    * 10,
                    ministics.loc[i, "esol"],
                )
            else:
                ministics.loc[i, "esol_a"] = 0
        else:  # Méthode 2
            if (ministics.loc[i - 1, "HUR1"] - SoilParams.HMINF_1) > 0:
                ministics.loc[i, "esol_a"] = min(
                    (ministics.loc[i - 1, "HUR1"] - SoilParams.HMINF_1)
                    * SoilParams.EPC_1
                    * 10,
                    ministics.loc[i, "esol"],
                )
            else:
                ministics.loc[i, "esol_a"] = 0

    # Calcul de l'eau restant dans le sol après évpaoration
    if i == 0:
        ministics.loc[i, "HUR1"] = SoilParams.HUR1_Init - (
            (ministics.loc[i, "esol_a"] * 0.001) / (SoilParams.EPC_1 / 100)
        )
    else:
        ministics.loc[i, "HUR1"] = ministics.loc[i - 1, "HUR1"] - (
            (ministics.loc[i, "esol_a"] * 0.001) / (SoilParams.EPC_1 / 100)
        )

    return ministics


def radiation_intercep(ministics, CropParams, Constants, i):
    '''
    Calcul de la radiation interceptée par la plante.
    
    Améliorations :
        - EXTIN doit varier journalièrement d'après STICS, reste à modéliser ça (dépend épaisseur feuille par exemple)
        - FSPM comme Caribu ou HELIOS peut faire mieux ?
    '''

    # Calcul de la radiation interceptée
    # Intuition : fraction de PAR dans radiation * radiation * (1 - exp( -coeff extinction * LAI))
    if i == 0:
        ministics.loc[i, "raint"] = 0
        ministics.loc[i, "raint_unstressed"] = 0
    else:
        ministics.loc[i, "raint"] = (
            0.95
            * Constants.PARSURRG
            * ministics.loc[i, "trg"]
            * (1 - np.exp(-CropParams.EXTIN * ministics.loc[i - 1, "LAI"]))
        )  # eq 9.1
        ministics.loc[i, "raint_unstressed"] = (
            0.95
            * Constants.PARSURRG
            * ministics.loc[i, "trg"]
            * (1 - np.exp(-CropParams.EXTIN * ministics.loc[i - 1, "LAI_unstressed"]))
        )  # eq 9.1

    return ministics


def total_biomass_growth(ministics, CropParams, Constants, i, ind_Z69):
    '''
    Calcul de la production journalière de biomasse

    Explications : 
    Principe : basé sur PAR interceptée, RUEmax et LAI
    Paramètre central = RUE (synthéthise PS et respiration) : dépend des stress, Temp, phénologie et coeff d'allocation carbone entre organes aériens et souterrrains (??)
    --> RUE considéré comme une fonction des stress : stress thermique (ftemp), stress hydrique stomatique (swfac), inns et exobiom
    '''

    if (
        ministics.loc[i, "DOY_2EMERG"] == 1 and ministics.loc[i, "DOY_2Z89"] == 0
    ):  # entre émergence et maturité
        # Calcul du RUE max (dépend du stade)
        # ATTENTION : doit être fait après LAI pr que ULAI soit calculé au jour i
        if ministics.loc[i, "ULAI"] < CropParams.VLAIMAX:  # si avant fin stade juvénile
            ministics.loc[i, "ebmax"] = CropParams.EFCROIJUV / 100
        elif i < ind_Z69:  # si avant stade idrp
            ministics.loc[i, "ebmax"] = CropParams.EFCROIVEG / 100
        else:  # si avant stade maturité (déjà cette condition dans grd if au-dessus)
            ministics.loc[i, "ebmax"] = CropParams.EFCROIREPRO / 100

        ministics.loc[i, "DELTAMS"] = shoot_biomass_production(
            ministics.loc[i, "raint"],
            ministics.loc[i, "ebmax"],
            Constants.COEFB,
            ministics.loc[i, "ftemp"],
            ministics.loc[i - 1, "SWFAC"],
            CropParams.FCO2,
        )
        #    + ministics.loc[i,'dltaremobil'] --> c'est fait dans biom partitioning
        ministics.loc[i, "DELTAMS_unstressed"] = shoot_biomass_production(
            ministics.loc[i, "raint_unstressed"],
            ministics.loc[i, "ebmax"],
            Constants.COEFB,
            ministics.loc[i, "ftemp"],
            ministics.loc[i - 1, "SWFAC"],
            CropParams.FCO2,
        )

    else:  # avant émerg ou après maturité = pas de production
        ministics.loc[i, "DELTAMS"] = 0
        ministics.loc[i, "DELTAMS_unstressed"] = 0

    # Calcul de la biomasse produite entre jour 1 et i, stressé ou non
    ministics.loc[i, "MASEC"] = ministics.DELTAMS[0:i].sum()
    ministics.loc[i, "MASEC_unstressed"] = ministics.DELTAMS_unstressed[0:i].sum()

    return ministics


def water_drainage(ministics, SoilParams, i):
    '''
    Calcul de l'eau drainée dans le sol en fin de journée.

    Explications :
    Principe : toute l'eau restante dans H1 et H2 est drainée en-dessous
    Question : au départ d'une journée on a donc tjrs HURx = pt flétriss ?
    '''

    # Calcul de l'eau restante dans H1 en fin de boucle qui va être drainée dans H2
    # S'il reste de l'eau dispo dans H1 -> DRAIN1 = eau H1 + précip - pt flétriss
    ministics.loc[i, "DRAIN1"] = (
        ministics.loc[i, "HUR1"] * SoilParams.EPC_1 * 10
        + ministics.loc[i, "Rain"]
        - SoilParams.HCCF_1 * SoilParams.EPC_1 * 10
        if (
            (
                ministics.loc[i, "HUR1"] * SoilParams.EPC_1 * 10
                + ministics.loc[i, "Rain"]
            )
            > (SoilParams.HCCF_1 * SoilParams.EPC_1 * 10)
        )
        else 0
    )

    # contenu eau H1 = contenu eau H1 + précip - DRAIN1 --> tout ce qui n'a pas été absorbé s'en va dans l'horizon 2 ?
    ministics.loc[i, "HUR1"] = (
        (
            ministics.loc[i, "HUR1"] * SoilParams.EPC_1 * 10
            + ministics.loc[i, "Rain"]
            - ministics.loc[i, "DRAIN1"]
        )
        / SoilParams.EPC_1
        / 10
    )

    # Calcul de l'eau restante dans H2 en fin de boucle qui va être drainée en-dessous
    # S'il reste de l'eau dispo dans H2 (eau H2 + drain1) -> DRAIN2 = eau H2 + DRAIN1 - pt flétriss
    ministics.loc[i, "DRAIN2"] = (
        ministics.loc[i, "HUR2"] * SoilParams.EPC_2 * 10
        + ministics.loc[i, "DRAIN1"]
        - SoilParams.HCCF_2 * SoilParams.EPC_2 * 10
        if (
            (
                ministics.loc[i, "HUR2"] * SoilParams.EPC_2 * 10
                + ministics.loc[i, "DRAIN1"]
            )
            > (SoilParams.HCCF_2 * SoilParams.EPC_2 * 10)
        )
        else 0
    )

    # conteneau eau H2 = contenu eau H2 + DRAIN1 - DRAIN2
    ministics.loc[i, "HUR2"] = (
        (
            ministics.loc[i, "HUR2"] * SoilParams.EPC_2 * 10
            + ministics.loc[i, "DRAIN1"]
            - ministics.loc[i, "DRAIN2"]
        )
        / SoilParams.EPC_2
        / 10
    )

    return ministics


def senescence_stress(ministics, CropParams, i, ind_EMERG):
    '''
    Calcul des indices de stress accélérant la sénescence.
    '''
    # A noter : le seuil du stress qui ralentit croiss feuille (TETURG) est plus haut que celui qui fait se fermer les stomates (TETSTOMATE). Donc la croiss des feuilles peut être inhinbée même pdt que la transpiration est au max de son potentiel

    if i >= ind_EMERG:
        ### STRESS ACCELERANT LA SENESCENCE ###

        # Calcul du seuil de teneur en eau tel que croiss feuilles ralentit
        ministics.loc[i, "TETURG"] = (
            1
            / 80
            * np.log(
                ministics.loc[i, "EOP"]
                / (
                    2
                    * np.pi
                    * ministics.loc[i, "CUMLRACZ"]
                    * CropParams.PSITURG
                    * 0.0001
                )
                * np.log(
                    1
                    / (
                        CropParams.RAYON
                        * (ministics.loc[i, "CUMLRACZ"] / ministics.loc[i, "ZRac"])
                        ** (1 / 2)
                    )
                )
            )
        )

        # Calcul du seuil de contenu en eau tel que sénescence accélérée
        ministics.loc[i, "tetsen"] = (
            CropParams.RAPSENTURG * ministics.loc[i, "TETURG"]
        )  # teturg = seuil de contenu en eau limitant la croissance des feuilles (en m2)

        # Calcul du stress hydrique -> ce facteur augmente taux de sénescence
        if ministics.loc[i, "TETA"] < ministics.loc[i, "tetsen"]:
            ministics.loc[i, "senfac"] = max(
                CropParams.SWFACMIN,
                ministics.loc[i, "TETA"] / ministics.loc[i, "tetsen"],
            )
        else:
            ministics.loc[i, "senfac"] = 1

        # Calcul de l'indice de stress de turgescence, j'ai mis TETA_moy parce que TETA c'est la teneur en eau totale, ça varie pas entre 0 et 1
        if ministics.loc[i, "TETA_moy"] < ministics.loc[i, "TETURG"]:
            ministics.loc[i, "turfac"] = max(
                CropParams.SWFACMIN,
                ministics.loc[i, "TETA_moy"] / ministics.loc[i, "TETURG"],
            )
        else:
            ministics.loc[i, "turfac"] = 1

        # Calcul du stress lié au gel qui impacte la scénescene des feuilles (avant ou après AMF)
        # Principe (FIGURE 4.16) : le stress lié se calcule à partir de 4 param : 2 param indép du stade phéno = Tmin de début d'action du gel où stress = 1 et Tlethal où stress = 0. Et 2 param dépendant du dommage : temp associé à 10% et 90% de dommage sur les feuilles.
        # Et entre chacun des 4 points la variation du stress est linéaire d'où mon calcul = on regarde où se situe le Tmin de la journée par rapport aux 4 params et on calcule le stress associé.
        # Impact différent selon si on est avant ou après le stade de croiss max des feuilles (= fin phase juvénile)

        if (
            ministics.loc[i, "ULAI"] < CropParams.VLAIMAX
        ):  # si on est avant stade iamf --> fgeljuv
            if ministics.loc[i, "Temp_min"] < CropParams.TGELJUV90:
                a = 0.1 / (CropParams.TGELJUV90 - CropParams.TLETALE)
                b = (
                    0.1
                    * np.abs(CropParams.TLETALE)
                    / (CropParams.TGELJUV90 - CropParams.TLETALE)
                )
                ministics.loc[i, "fstressgel"] = max(
                    0, ministics.loc[i, "Temp_min"] * a + b
                )
            elif ministics.loc[i, "Temp_min"] < CropParams.TGELJUV10:
                a = 0.8 / (CropParams.TGELJUV10 - CropParams.TGELJUV90)
                b = (
                    0.8
                    * np.abs(CropParams.TGELJUV90)
                    / (CropParams.TGELJUV10 - CropParams.TGELJUV90)
                    + 0.1
                )
                ministics.loc[i, "fstressgel"] = ministics.loc[i, "Temp_min"] * a + b
            else:
                a = 0.1 / (CropParams.TDEBGEL - CropParams.TGELJUV10)
                b = (
                    0.1
                    * np.abs(CropParams.TGELJUV10)
                    / (CropParams.TDEBGEL - CropParams.TGELJUV10)
                    + 0.9
                )
                ministics.loc[i, "fstressgel"] = min(
                    1, ministics.loc[i, "Temp_min"] * a + b
                )

        else:  # après stade iamf -> fgelveg
            if ministics.loc[i, "Temp_min"] < CropParams.TGELVEG90:
                a = 0.1 / (CropParams.TGELVEG90 - CropParams.TLETALE)
                b = (
                    0.1
                    * np.abs(CropParams.TLETALE)
                    / (CropParams.TGELJUV90 - CropParams.TLETALE)
                )
                ministics.loc[i, "fstressgel"] = max(
                    0, ministics.loc[i, "Temp_min"] * a + b
                )
            elif ministics.loc[i, "Temp_min"] < CropParams.TGELVEG10:
                a = 0.8 / (CropParams.TGELVEG10 - CropParams.TGELVEG90)
                b = (
                    0.8
                    * np.abs(CropParams.TGELVEG90)
                    / (CropParams.TGELVEG10 - CropParams.TGELVEG90)
                    + 0.1
                )
                ministics.loc[i, "fstressgel"] = ministics.loc[i, "Temp_min"] * a + b
            else:
                a = 0.1 / (CropParams.TDEBGEL - CropParams.TGELVEG10)
                b = (
                    0.1
                    * np.abs(CropParams.TGELVEG10)
                    / (CropParams.TDEBGEL - CropParams.TGELVEG10)
                    + 0.9
                )
                ministics.loc[i, "fstressgel"] = min(
                    1, ministics.loc[i, "Temp_min"] * a + b
                )

    return ministics


def senescence(ministics, CropParams, i, ind_EMERG):
    '''
    Module de sénescence EN DEVELOPPEMENT 
    
    '''
    if i >= ind_EMERG:
        # Calcul de la durée de vie max des feuilles = ce qu'on ajoute chaque jour diminue jusqu'à stade ILAX, et ensuite on ajoute plus rien
        if (ministics.loc[i, "ULAI"] > 0) and (
            ministics.loc[i, "ULAI"] < CropParams.VLAIMAX
        ):  # Avant le max de croiss LAI, la durée de vie des feuilles = durée de vie des feuilles précoces
            ministics.loc[
                i, "durage"
            ] = (
                CropParams.DURVIEI
            )  # 'calculé quand leaf emitted' ça veut dire qu'avant = 0. C quoi leaf emitted ?
        elif (ministics.loc[i, "ULAI"] >= CropParams.VLAIMAX) and (
            ministics.loc[i, "ULAI"] < 3
        ):  # Après le max de croiss LAI, durée vie feuilles = durée vie feuilles précoces + avancement relatif entre max croiss et max LAI * (durée vie dernières feuilles - durée vie feuilles précoce)
            # durée se rapproche de DURVIEF entre le stade iamf et ilax
            ministics.loc[i, "durage"] = CropParams.DURVIEI + (
                CropParams.DURVIEF - CropParams.DURVIEI
            ) * (ministics.loc[i, "ULAI"] - CropParams.VLAIMAX) / (
                3 - CropParams.VLAIMAX
            )

        # Calcul de l'indice de stress hydrique qui augmente le taux de sénescence
        ministics.loc[i, "senstress"] = min(
            ministics.loc[i, "senfac"], ministics.loc[i, "fstressgel"]
        )  # il y a N normalement aussi

        # Calcul de la durée de vie (avec stress pris en compte) au bout de laquelle la matière produite est perdue par sénescence
        #  t0 porte à confusion ds STICS Eq 4.10
        # En cas de grde dispo d'azote, STICS ajotue qqch, ici je ne fais rien, à voir.
        # Intuition = durée de vie restante ou déjà vécue ? car on la fait diminuer là, donc ça doit être la restante mais on * par durage
        ministics.loc[i, "durvie"] = (
            ministics.loc[i, "durage"] * ministics.loc[i, "senstress"]
        )
        # Rq : le durage t ne dépend pas de t-1. Donc si y'a eu plein de stress avant ça change rien pr la suite ?

        # Calcul du temps thermique cumulé (adapté à sénescence, pas gdd) au jour i (ensuite on cumule)
        ministics.loc[i, "somsen"] = (
            2
            * ministics.loc[ind_EMERG:i, "UDEVCULT"]
            * (
                CropParams.STRESSDEV * ministics.loc[ind_EMERG:i, "turfac"]
                + 1
                - CropParams.STRESSDEV
            )
        ).sum()

        # Calcul de la biomasse et lai perdus par sénescence quand durvie est dépassée
        # C'est pas progressif, tout devient sénescent d'un coup quand durée de vie globale dépassée ?
        if ministics.loc[i, "somsen"] > ministics.loc[ind_EMERG:i, "durvie"].sum():
            ministics.loc[i, "dltaisen"] = ministics.loc[ind_EMERG:i, "DELTAI"].sum()
            ministics.loc[i, "dltamsen"] = (
                CropParams.RATIOSEN * ministics.loc[ind_EMERG:i, "dltafv"].sum()
            )  # calcul de pfeuilverte improvisé : dltafv / deltams, les deltams s'annulent

            # diff entre dltaisen et laisen ? parce que là dltaisen est déjà cumulé
            # ATTENTION : j'ai peut être fait une erreur ici, les t0 il y en a peut être plein

            # laisen = LAI des feuilles sénescentes = en m2/m2 = foliage sénescent cumulé
            # dltaisen = chgt journalier du LAI sénescent = en m2/m2/jour

            # je pense que chaque jour y'a x matière produite et on calcul sa durée de vie. Quand sa durée de vie est dépassée, on l'ajoute à dltamsen et dltaisen.

            # Questions PAC :
            #   - 'calculated for the day when leaves are emitted (t0)' = C'est quoi t0 car c important ? ça veut dire qu'on calcule DURAGE une seule fois ou chaque jour à partir du LAI journalier ?
            #               --> Et du coup le somsen depuis t0 ce serait à calculé depuis t01, t02 etc ou y'a 1 seul t0 qui existe
            #   - durage augmente au cours du temps, est-ce logique avec ce que tu m'as dis ? montrer dessin. Et durage est d'autant plus grd qu'on se rapproche de LAI max, donc ça va un peu à l'encontre de la sensibilité dont tu parles
            #   - quand somsen dépasse somme de durage --> sénescence. Déjà pq c'est somme de durage et pas juste durage du jour t0 (comme écrit texte jusqute au-dessus).
            #      --> Autre conséq : la biomasse produite le jour j n'est pas perdue au cours d'un processus, mais tout est pérdue (enfin RATIO*biom) le jour où somsen dépasse durage
            #   - et pq y'a dltaisen et laisen (cumulé) alors que dltaisen il est déjà cumulé. dltaisen serait calculé chaque jour et on le recumule pour voir le sénescent total?
            #

        # On fait ça le temps d'avoir une réponse sur le forum stics
        ministics.loc[i, "dltaisen"] = 0
        ministics.loc[i, "dltamsen"] = 0

    return ministics
