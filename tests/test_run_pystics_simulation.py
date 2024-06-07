from pystics.params import parametrization_from_stics_example_files
from pystics.simulation import run_pystics_simulation
import numpy as np
import pandas as pd
import os

def test_tun_pystics_simulation():

    species = 'wheat'
    variety = 'Talent'

    mocked_dir = os.path.dirname(os.path.abspath(__file__)) + "/mocked_data"
    
    weather_pystics, crop, manage, soil, station, constants, initial = parametrization_from_stics_example_files(species=species, variety=variety, xml_folder_path = mocked_dir + '/mocked_param_files')
    pystics_ble_test, mat = run_pystics_simulation(weather_pystics, crop, soil, constants, manage, station, initial)

    pystics_ble_mocked = pd.read_csv(mocked_dir + '/pystics_simu_results.csv')

    assert np.allclose(pystics_ble_test.mafruit.values, pystics_ble_mocked.mafruit.values)
    assert  np.allclose(pystics_ble_test.resrac.values, pystics_ble_mocked.resrac.values)
    assert np.allclose(pystics_ble_test.swfac.values, pystics_ble_mocked.swfac.values)
    assert np.allclose(pystics_ble_test.dltams.values, pystics_ble_mocked.dltams.values)
    assert np.allclose(pystics_ble_test.ftemp.values, pystics_ble_mocked.ftemp.values)