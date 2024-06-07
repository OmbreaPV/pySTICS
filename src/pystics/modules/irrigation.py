import numpy as np
from ..exceptions import pysticsException

def irrigation_dates_and_amounts(weather,manage):
    ''' 
    This module adds irrigation amounts to associated julian days given in inputs (codecalirrig = 2) in output dataframe.
    '''
    
    weather.insert(0,'airg',[0 for i in range(len(weather))])
    
    for day,amount in manage.IRRIGATION_INTERVENTIONS.items():
        if int(day) not in weather.index:
            raise pysticsException(
                detail="Irrigation intervention day is not in weather index",
            )
        weather.loc[int(day),'airg'] = amount

    return weather

def dates_irrig(daily_outputs, datedeb_irrigauto, datefin_irrigauto):
    ''''
    This module is called when codedate_irrigauto = 1 --> dates of irrigation activation (automatic computation).
    It returns a list of booleans : True = active irrigation, False = no irrigation.
    '''

    compute_airg_list = False
    index = np.array([i for i in range(len(daily_outputs))])
    compute_airg_list[(index > datedeb_irrigauto) & (index <= datefin_irrigauto)] = True

    return compute_airg_list

def automated_irrig(codedate_irrigauto, compute_airg, ratiol, swfac_prev, effirr, doseirrigmin, dosimx, hcc, hur_i_prev, zrac_prev, doy, iplt0, irrlev):
    '''
    This module is called if codecalirrig = 1 and automatically computes irrigation each day.

    Simplifications compared to STICS :
        - Undergound irrigation water can be evaporated
    '''
    if doy == iplt0:
        airg = irrlev
    else:
        airg = 0
        if (codedate_irrigauto == 3) | ((codedate_irrigauto  == 1) & compute_airg):
            if swfac_prev < ratiol:
                airg = effirr * max(0, - doseirrigmin + min(dosimx, (hcc[0:int(zrac_prev)] - hur_i_prev[0:int(zrac_prev)]).sum()))
            
    return airg