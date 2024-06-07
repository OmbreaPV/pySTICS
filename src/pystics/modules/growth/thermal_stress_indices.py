import numpy as np

def thermal_stress_on_biomass_broduction(tcult_prev, temin, temax, teopt, teoptbis):
    '''
    This function computes the thermal stress index affecting biomass production (through RUE).
    '''

    ftemp = np.array(1., dtype = np.float64) if (type(tcult_prev) == np.float64) else np.zeros_like(tcult_prev)+1 
    mask1 = tcult_prev <= teopt
    mask2 = tcult_prev >= teoptbis

    ftemp[mask1] = 1 - ((tcult_prev[mask1] - teopt) / (temin - teopt))**2
    ftemp[mask2] = 1 - ((tcult_prev[mask2] - teoptbis) / (temax - teoptbis))**2

    ftemp = np.clip(ftemp, 0, 1)

    return ftemp

def thermal_stress_on_grain_filling(codetremp, temp_min_prev, tcultmax_prev, tminremp, tmaxremp):
    '''
    This function computes the thermal stress index affecting grains/fruits filling.
    '''

    if codetremp == 1:
        ftempremp = np.where((temp_min_prev > tminremp) & (tcultmax_prev < tmaxremp), 1, 0)

    elif codetremp == 2:
        ftempremp = 1

    return ftempremp

def frost_stress_on_fruit_number(temp_min, tgelflo90, tgelflo10, tletale, tdebgel,codgelflo):
    '''
    This function computes the frost stress index affecting grain/fruit number.
    '''

    if codgelflo == 2:
        fgelflo = frost_stress(temp_min,
                                tgelflo90,
                                tgelflo10,
                                tletale,
                                tdebgel)
    elif codgelflo == 1:
        fgelflo = 1

    return fgelflo


def frost_stress(t, tgel90, tgel10, tletale, tdebgel):
    '''
    This function computes the frost stress index.
    '''

    stress = np.zeros_like(t)
    a = np.zeros_like(t)
    b = np.zeros_like(t)

    mask1 = t >= tdebgel
    mask2 = (t >= tgel10) & (tgel10 < tdebgel)
    mask3 = (t < tgel10) & (t >= tgel90)
    mask4 = t < tgel90
    
    try:
        a[mask2] = (0.9 - 1.0) / (tgel10 - tdebgel)
    except ZeroDivisionError:
        pass
            
    a[mask3] = (0.9 - 0.1) / (tgel10 - tgel90)
    try:
        a[mask4] = (0. - 0.1) / (tletale - tgel90)
    except ZeroDivisionError:
        pass
    
    b[mask2] = 1.0 - (a[mask2] * tdebgel)
    b[mask3] = 0.9 - (a[mask3] * tgel10)
    b[mask4] = 0. - (a[mask4] * tletale)
    
    stress[mask1] = 1

    stress[mask2 | mask3 | mask4] = (a[mask2 | mask3 | mask4] * t[mask2 | mask3 | mask4]) + b[mask2 | mask3 | mask4]
    
    stress = np.clip(stress, 0, 1)

    return stress
