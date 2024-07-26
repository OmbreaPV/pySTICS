


def water_balance(trr, airg, mouill, depth, hur_i, epz_i_prev, esz_i, hcc, hurlim):
    '''
    This module :
        - allocates daily rain (trr) and irrigation (airg) cm per cm (top-down) in hur.
        - substracts esz (evaporation) and epz (transpiration) cm per cm in hur.
        - computes drained water.
    '''
    
    dispo = trr + airg - mouill
    for z_index in range(depth):
        hur_i[z_index] = hur_i[z_index] - epz_i_prev[z_index] - esz_i[z_index]
        distr = min(hcc[z_index] - hur_i[z_index], dispo)
        hur_i[z_index] = hur_i[ z_index] + distr
        if hur_i[z_index] < hurlim:
            hur_i[z_index] = hurlim
        dispo = max(0,dispo - distr)
    drain = dispo

    return hur_i, drain