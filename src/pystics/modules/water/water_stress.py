import numpy as np


def water_stress(eop, swfacmin, rapsenturg, cumlracz, psiturg, psisto, rayon, zrac, lev_i_prev, hmin, hur_i_prev, profsem, depth):
    """'
    This modules computes water content in the root zone (teta, resrac), water content thresholds (tetstomate, teturg, tetsen) and water stress indexes (swfac, turfac, senfac)
    """

    teta, tetstomate, teturg, tetsen, resrac = 0, 0, 0, 0, 0
    swfac, turfac, senfac = 1, 1, 1

    if zrac > 0:
        
        # Available water in the root zone
        root_range = range(max(0,int(profsem)-1), min(depth, int(zrac)+1))
        resrac = np.maximum(hur_i_prev[root_range] - hmin[root_range], 0.).sum()

        teta = resrac / (10 * zrac)


    if (lev_i_prev == 1):

        if cumlracz == 0:
            swfac = swfacmin
            if rapsenturg != 0:
                senfac = swfacmin
            turfac = swfacmin
        
        elif eop > 0:

            # Temporary wilting point for swfac
            tetstomate = temporary_wilting_point(psisto, eop, cumlracz,rayon, zrac)

           # Water content threshold for turfac
            teturg = temporary_wilting_point(psiturg, eop, cumlracz,rayon, zrac)

            # Water content threshold for senfac
            tetsen = rapsenturg * teturg

            # Turgescence water stress index
            turfac = (
                1
                if teta  >= teturg
                else max(swfacmin, teta / teturg)
            )
            
            # Stomatic water stress index
            swfac = (
                1
                if teta  > tetstomate
                else max(swfacmin, teta / tetstomate)
            )

            # Water stress index increasing senescence rate
            if rapsenturg > 0:
                senfac = (
                    1
                    if teta  >= tetsen
                    else max(swfacmin, teta / tetsen)
                )


    return teta, swfac, tetstomate, turfac, teturg, senfac, tetsen, resrac


def temporary_wilting_point(psisto_or_psiturg, eop, cumlracz,rayon, zrac):
    return (1
                / 80
                * np.log(
                    eop
                    / (
                        2
                        * np.pi
                        * cumlracz
                        * psisto_or_psiturg
                        * 0.0001
                    )
                    * np.log(
                        1
                        / (
                            rayon
                            * (np.pi * cumlracz / zrac)
                            ** (1 / 2)
                        )
                    )
                )
            )


def water_stress_on_root_growth(hur, hmin, hcc, sensrsec, code_humirac):
    '''
    This module computes the effect of soil water content on germination and root growth.
    '''

    if code_humirac == 0:
        return 1

    if code_humirac == 2:
        if hur > hmin:
            x = (hur - hmin) / (hcc - hmin)
            humirac = sensrsec + (1-sensrsec) * hur
        else:
            humirac = sensrsec * hur / hmin
    else:
        if hur >= hmin:
            humirac = 1
        else:
            humirac = sensrsec * hur / hmin

    humirac = np.clip(humirac, 0, 1)

    return humirac