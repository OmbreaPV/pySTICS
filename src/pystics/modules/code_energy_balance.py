import numpy as np
from pystics.modules.water.potential_evapotranspiration import tvar

def energy_budget(raint, parsurrg, trg, rnet, zr, z0solnu, z0, dh, wind, hauteur, lai_prev, temp, tpm, rsmin, fco2s):


    ### Available energy and soil/plant distribution ###
    # pas calculée encore : Emd

    fapar = raint / (
        parsurrg * trg
    )

    # Available energy for the plant
    rnetp1 = (
        0.83 * fapar * rnet
    )  # eq 9.26

    # Available energy for the soil
    rnets = rnet - rnetp1

    ### Diffusion resistances ###

    ras0 = (
        np.log(zr / z0solnu)
        * np.log((z0 + dh) / z0solnu)
        / (0.41**2 * wind)
    )
    raa0 = (
        np.log(zr / z0solnu) ** 2
        / (0.41**2 * wind)
        - ras0
    )

    rasinf = np.log(
        (zr - dh) / z0
    ) / (0.41**2 * wind) - hauteur / (
        2.5 * (hauteur - dh)
    ) * (
        12.18
        - np.exp(
            2.5
            * (
                1
                - (dh + z0)
                / hauteur
            )
        )
    )
    raainf = (
        np.log((zr - dh) / z0)
        / (0.41**2 * wind)
        * (
            np.log((zr - dh) / hauteur)
            + hauteur
            / (2.5 * (hauteur - dh))
            * np.exp(
                2.5
                * (
                    1
                    - (dh + z0)
                    / hauteur
                )
                - 1
            )
        )
    )

    ras = 0.25 * (
        lai_prev * rasinf
        + (4 - lai_prev) * ras0
    )
    raa = 0.25 * (
        lai_prev * raainf
        + (4 - lai_prev) * raa0
    )


    deltat = tvar(temp + 0.5) - tvar(
        temp - 0.5
    )

    L = (2500840 - 2358.6 * temp) * 1e-6

    dsat = tvar(temp) - tpm

    ### Surface resistances ###

    if (
        lai_prev != 0
    ):
        rac = min(12.5, 50 / (2 * lai_prev))

        rc = (
            rsmin
            * (0.5 * lai_prev + 1)
            / lai_prev
            * (0.039 * dsat + 0.45)
            * 28
            / (2.5 + trg)
            * fco2s
        )

    return rnetp1, rnets, ras0, rasinf, raa0, raainf, ras, raa, deltat, L, dsat, rac, rc


def intermed_sat_deficit(rnets, deltat, dsat, rnet, raa, gamma):

    ept = (
        1.32
        * rnets
        * deltat
        / (deltat + gamma)
    )

    dos = (
        dsat
        + (
            gamma * rnet
            - (deltat + gamma) * ept
        )
        * raa
        / 105.03
    )

    return ept, dos

def foliage_water_evap_energy_budget(codeintercept, lai_prev, gamma, mouill, deltat, rnetp1, dos, rac, L):

    if codeintercept == 2:
        emd = 0
    elif lai_prev == 0:
        emd = 0
    else:
        emd = max(
            mouill,
            (
                deltat * rnetp1
                + 105.03 * dos / rac
            )
            / (L * (deltat + gamma)),
        )

    return emd