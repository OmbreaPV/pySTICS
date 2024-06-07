import numpy as np


def water_stress(eop, swfacmin, rapsenturg, cumlracz, psiturg, psisto, rayon, zrac, lev_i_prev, hmin, hur_i_prev, profsem, depth):
    """'
    This modules computes water content in the root zone (teta, resrac), water content thresholds (tetstomate, teturg, tetsen) and water stress indexes (swfac, turfac, senfac)
    """

    teta, tetstomate, teturg, tetsen, resrac = 0, 0, 0, 0, 0
    swfac, turfac, senfac = 1, 1, 1

    if zrac > 0:
        
        # Available water in the root zone
        resrac = 0.
        for z_index in range(max(0,int(profsem)-1), min(depth, int(zrac)+1)):
            resrac = resrac + max(hur_i_prev[z_index] - hmin[z_index], 0.)

        teta = resrac / (10 * zrac)


    if (lev_i_prev == 1):

        if cumlracz == 0:
            swfac = swfacmin
            senfac = swfacmin
            turfac = swfacmin
        
        elif eop > 0:

            # Temporary wilting point for swfac
            tetstomate = (
                1
                / 80
                * np.log(
                    eop
                    / (
                        2
                        * np.pi
                        * cumlracz
                        * psisto
                        * 0.0001
                    )
                    * np.log(
                        1
                        / (
                            rayon
                            * (np.pi * 
                                cumlracz
                                / zrac
                            )
                            ** (1 / 2)
                        )
                    )
                )
            )
            
           # Water content threshold for turfac
            teturg = (
                1
                / 80
                * np.log(
                    eop
                    / (
                        2
                        * np.pi
                        * cumlracz
                        * psiturg
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

            # Water content threshold for senfac
            tetsen = (
                rapsenturg * teturg
            ) 

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
            senfac = (
                1
                if teta  >= tetsen
                else max(swfacmin, teta / tetsen)
            )

    return teta, swfac, tetstomate, turfac, teturg, senfac, tetsen, resrac