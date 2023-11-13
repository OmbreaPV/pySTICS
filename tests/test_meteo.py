from pystics.meteo import generate_meteo
import numpy as np

def test_meteo():
    lat = 39.742043
    lon = -104.991531
    meteo = generate_meteo(lat, lon, start_year=2012, end_year = 2015, save=False)
    assert np.isclose(meteo.Temp_moy.mean(), 10.314397102441252, rtol=0.01)
    assert np.isclose(meteo.Radiation.mean(), 17.26999175359343, rtol=0.01)
    assert np.isclose(meteo.Rain.mean(), 1.482014, rtol=0.01)
    assert np.isclose(meteo.Wind.mean(), 2.719121, rtol=0.01)
    assert np.isclose(meteo.Rel_Hum.mean(), 42.396193, rtol=0.01)
    assert np.isclose(meteo.Photoperiod.mean(), 12.232950, rtol=0.01)