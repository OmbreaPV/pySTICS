import math
import pandas as pd
import numpy as np
from datetime import datetime
from .et0 import fao
import ast
from ..exceptions import pysticsException
from pystics.modules.water.potential_evapotranspiration import potential_etp

def check_weather_variables(weather, codeetp):

    # Check mandatory columns
    if not all([var in weather.columns for var in ['doy','temp_min','temp_max','radiation','rain','co2']]):
        raise pysticsException(
                detail="Not all mandatory weather inputs are in weather file",
            )
    
    # Rename Radiation and Rain
    weather = weather.rename(columns={'radiation':'trg', 'rain':'trr'})
    
    # Check for wind if codeetp = 2
    if (codeetp == 2) & (not 'wind' in weather):
        raise pysticsException(
            detail="Penman evapotranspiration (codeetp = 2) but wind variable is not in weather file",
        )
    
    # Check for etp if codeetp = 1
    if (codeetp == 1) & (not 'etp' in weather):
        raise pysticsException(
            detail="Forced etp (codeetp = 1) but etp variable is not in weather file",
        )
    
    # Add variables if not in weather file
    if 'wind' not in weather.columns:
        weather.insert(0,'wind',[np.nan for i in range(len(weather))])

    if 'etp' not in weather.columns:
        weather.insert(0,'etp',[0 for i in range(len(weather))])
    elif weather.etp[0] == -999.9:
        weather['etp'] = [0 for i in range(len(weather))]

    if 'tpm' not in weather.columns:
        weather.insert(0,'tpm',[0 for i in range(len(weather))])

    return weather

def compute_weather_variables(df, station, gamma):
    """
    This function computes climatic variables from weather variables.
    List of columns of weather: year / month / day / doy / temp_min / temp_max / radiation / etp / rain / wind / tpm / co2
    etp / tpm are not mandatory (computed here if not in weather)
    """

    # Mean temperature
    df['temp'] = (df.temp_min + df.temp_max) / 2

    # STICS method for extraterrestrial radiations. (FAO method is defined but not used).
    df['rgex'] = [rgex_stics(station.LATITUDE, df.doy[i]) for i in df.index]
    # df['rgex'] = [rgex_fao(station.LATITUDE, df.doy[i]) for i in df.index]

    df["esat"] = 0.6108 * np.exp(17.27 * (df["temp_max"] + df["temp_min"]) / 2 / ((df["temp_max"] + df["temp_min"]) / 2 + 237.3))
    df["delta"] = 4098 * df["esat"] / ((df["temp"] + 237.3) ** 2)

    # Water vapour pressure in air
    df['tpm'] = (df['temp_min'] - station.CORECTROSEE).apply(tvar) # no option to give it as input in weather file

    
    # Photoperiod
    df['phoi'] = df.doy.apply(lambda x:photoperiod(station.LATITUDE, x))

    # Hourly temperature
    if 'hourly_temp' not in df.columns:
        df['hourly_temp'] = 0
    if type(df.hourly_temp[0]) is str:
        df["hourly_temp"] = [ast.literal_eval(i) for i in df.hourly_temp]


    # Insolation fraction : STICS methid (FAO coefficients are different, see below)
    df['fracinsol'] = df.loc[:,['trg','rgex']].apply(lambda x : min(max(((x['trg'] / x['rgex']) - station.AANGST) / station.BANGST, 0), 1),
                                   axis=1)
    # fracinsol = (1.35 * (trg / rgex) - 0.35) # FAO method see net_out_lw_rad function https://github.com/woodcrafty/PyETo/blob/master/pyeto/fao.py

    # Potential evapotranspiration
    if station.CODEETP != 1:
        df['etp'] = df.apply(lambda x : potential_etp(x['trg'], x['temp'], x['wind'], x['tpm'], gamma, station.CODEETP, station.ALPHAPT, x['fracinsol']),
                                                   axis=1)

    return df

def tvar(x):
    '''
    This function computes water vapour pressure in air.
    '''
    return 6.1070 * (1 + 2 ** (1 / 2) * np.sin(0.017453293 * x / 3)) ** 8.827


def photoperiod(zlat,jday):
    ''' 
    This function computes the photoperiod, see STICS implementation.
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

    photp = [0.0, 0.0, 0.0]

    for i in range(3):
        theta1 = 2.0 * pi * (days[i] - 80.0) / 365.0
        theta2 = 0.034 * (math.sin(2.0 * pi * days[i] / 365.0) - math.sin(2.0 * pi * 80.0 / 365.0))
        theta = theta1 + theta2
        dec = math.asin(0.3978 * math.sin(theta))
        d = -1.0 * 0.10453 / (math.cos(alat) * math.cos(dec))
        p = d - (math.tan(alat) * math.tan(dec))
        p = math.acos(p)
        photp[i] = 24.0 * p / pi
        if photp[i] > 24.0:
            photp[i] = 24.0

    photp = photp[1]
    return photp


def rgex_fao(latitude, doy):
    '''
    This function computes extraterrestrial radiations with FAO method.
    '''
    latitude_rad = latitude * np.pi / 180 
    sol_dec = fao.sol_dec(doy)
    sunset_hour_angle = fao.sunset_hour_angle(latitude_rad, sol_dec)
    inv_dist_earth_sun = fao.inv_rel_dist_earth_sun(doy)

    return fao.et_rad(latitude_rad, sol_dec, sunset_hour_angle, inv_dist_earth_sun) * 0.75 # *0.75 = from extraterrestrial radiations to clear sky radiation (fao cs_rad)

def rgex_stics(lat, jul):
    '''
    This function computes extraterrestrial radiations with STICS method.
    '''
    lat = lat / 180 * 3.14

    pi = 4 * np.arctan(1.0)
    theta1 = 2 * pi * (jul - 80) / 365
    theta2 = 0.034 * (np.sin(2 * pi * jul / 365) - np.sin(2 * pi * 80 / 365))
    theta = theta1 - theta2
    z = np.arcsin(0.3978 * np.sin(theta))

    x = np.sin(z)
    y = np.sin(lat)

    solar = 1370 * 3600 * 24 / 3.14159
    a = -x * y / ((1 - x**2) * (1 - y**2))**(1/2)
    rgex = 0. * x
    if (a < 0.0):
        rgex = 3.14159
        if (a + 1.0 < 0.0):
            a = -1.0
        u = (1 - a**2)**(1/2) / a
        rgex = x * y * (rgex + np.arctan(u) - u)

    if (abs(a) < 1.0e-8):
        if (abs(x) < 1.0e-8):
            rgex = (1 - y**2)**(1/2)
        else:
            rgex = (1 - x**2)**(1/2)    

    if (a > 0.0):
        u = (1 - a**2)**(1/2) / a
        rgex = x * y * (rgex + np.arctan(u) - u)

    rgex = solar * (1 + (0.033 * np.cos(0.0172 * jul))) * rgex * 1e-6

    return rgex






