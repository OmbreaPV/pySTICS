import pandas as pd
import numpy as np
from pystics.get_params import parametrization_from_usm_example
from pystics.tasks import run_stics_simulation

from pathlib import Path

def test_wheat_biomass():

    # Lecture des fichiers exemples de STICS de l'USM associé aux espèce et variété choisies
    meteo, crop, soil, manage, station, constants, user = parametrization_from_usm_example('wheat', 'Thesee', meteo_source='stics')

    # Lancement de la simulation
    stics, mat = run_stics_simulation(crop, soil, constants, meteo, user, manage, station)

    # assert np.isclose(stics.MASEC.sum(), 904.7551820694706, rtol=0.01)
    assert stics.mafruit.sum() > 5