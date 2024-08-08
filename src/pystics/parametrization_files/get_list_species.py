# %%
import pandas as pd
import xmltodict
import os
import sys
from pystics.params import gen_dict_extract
# %%
dict_list = []
for i in os.listdir('plant/'):
    with open('plant/' + i, 'r') as f:
        data = f.read()
        dict_list.append(xmltodict.parse(data))
# %%
results = list(set(gen_dict_extract('#text','@nom',dict_list[0])))
results
# %%
with open('example/usms.xml', 'r') as f:
    data = f.read()
dico = xmltodict.parse(data)
species_list = [dic['@nom'] for dic in dico['usms']['usm']]
# %%
from pystics.params import parametrization_from_stics_example_files
weather, crop, manage, soil, station, constants, initial = parametrization_from_stics_example_files(species_list[0], 'Talent')
# %%
