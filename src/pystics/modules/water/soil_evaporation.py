import numpy as np


def potential_soil_evap_beer_law(etp, extin, lai):
    '''
    This module computes the potential soil evaporation (eos) with Beer law approach.
    '''
    eos = etp * np.exp(
        -(extin- 0.2) * lai
    )

    return eos


def potential_soil_evap_energy_budget(deltat, rnets, dos, ras, L, gamma):
    '''
    This module computes the potential soil evaporation (eos) with the resistance approach.
    '''

    eos = (
        deltat * rnets
        + 105.03 * dos/ ras
    ) / (L * (deltat + gamma))
    
    return eos


def fper(x, a):
    return (a**2 + 2 * a*x)**(1/2) - a

def finv(x, a):
    return ((x+a)**2 - a**2) / (2 * a)

def actual_soil_evaporation(trr, airg, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc):
    '''
    This module computes the actual soil evaporation (esol).
    See section 11.2.2 of STICS book.
    '''

    precipsol = trr + airg

    sum2 = sumes2 + ses2j0

    esol = 0

    if (sumes1 >= q0) & (precipsol >= sum2):
        return phase1(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc, precipsol, sum2)

    if (sumes1 >= q0) & (precipsol < sum2):
        return phase2(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc, precipsol, sum2)

    if precipsol >= sumes1:
        return phase1_end(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc)
    
    sumes1 = sumes1 - precipsol

    return phase2_test(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc)


def phase1(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc, precipsol, sum2):

    if precipsol < sumes2:
        return phase2(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc, precipsol, sum2)
    sumes1 = q0 - (precipsol - sumes2)
    sumes0 = 0.
    sumes2 = 0.
    if precipsol > q0:
        return phase1_end(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc)
    return phase2_test(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc)



def phase2(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc, precipsol, sum2):


    if (precipsol <= 0) & (nstoc == 0):
        sumes0 = sumes0 + eos
        esol = fper(sumes0, aevap) - sumes2
        sumes2 = fper(sumes0, aevap)

        return esol, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc

    if (precipsol > 0) | (nstoc > 0):
        stoc = stoc + precipsol
        if nstoc == 0:
            sesj0 = sumes0
            ses2j0 = sumes2
        if precipsol > 0:
            sumes2 = 0
            smes02 = 0
        if (nstoc > 0) & (stoc >= ses2j0):
            precipsol = stoc - ses2j0
            ses2j0 = 0
            sumes2 = 0
            nstoc = 0
            stoc = 0
            smes02 = 0
            return phase1(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc, precipsol, sum2)
        smes02 = smes02 + eos
        esol = fper(smes02, aevap) - sumes2
        sumes2 = fper(smes02, aevap)
        if sumes2 < stoc:
            nstoc = nstoc + 1
        else:
            nstoc = 0
            stoc = 0
            smes02 = 0
            sumes0 = sesj0 + finv(sumes2 - stoc, aevap)
            sumes2 = fper(sumes0, aevap)
            ses2j0 = 0.

    return esol, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc


def phase1_end(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc):

    sumes1 = 0.
    return phase2_test(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc)

def phase2_test(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc):

    sumes1 = sumes1 + eos

    if sumes1 > q0:
        return phase2_begin(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc)
    esol = eos

    return esol, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc


def phase2_begin(esol, eos, aevap, q0, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc):

    sumes0 = sumes1 - q0
    sumes2 = fper(sumes0, aevap)
    esol = eos - (sumes0 - sumes2)

    return esol, sumes0, sumes1, sumes2, sesj0, ses2j0, smes02, nstoc, stoc


def soil_layers_contribution_to_evaporation(esol, zesx, hurlim, cfes, hcc, wi_i, hur_i, esz_i):
    '''
    This modules computes the soil layers contributions to evaporation (esz).
    See section 11.2.3 of STICS book.
    '''
    
    sum = 0

    # Soil water availability and soil contribution to evaporation
    if cfes > 0:
        wi_i[range(zesx)] = (hur_i[range(zesx)] - hurlim) / (hcc[range(zesx)] - hurlim)
        esz_i[range(zesx)] = wi_i[range(zesx)] * (1.-np.array([z_index for z_index in range(zesx)])/zesx)**(abs(cfes))
        sum = esz_i[range(zesx)].sum()


    if sum == 0:
        sum = 1
        esol = 0
    else:
        esz_i = esz_i * esol / sum
    
    esz_i[zesx:] = 0 

    # Check if soil layer contribution does not exceed available water
    esz_i[0:zesx] = np.minimum(esz_i[0:zesx], np.maximum(hur_i[0:zesx] - hurlim, np.zeros(zesx)))
    esol = esz_i[0:zesx].sum()

    return esol, esz_i, wi_i
