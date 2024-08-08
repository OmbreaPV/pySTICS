# %%
import pystics
import pandas as pd
import numpy as np
from pystics.simulation import run_pystics_simulation
import seaborn as sns
import matplotlib.pyplot as plt
# %%
##### Choice of species and variety --> the USM associated to the species 
species = 'wheat'
variety = 'Talent'
# %%
# Read example input files from STICS for the USM associated to chosen species and variety
from pystics.params import parametrization_from_stics_example_files
weather, crop, manage, soil, station, constants, initial = parametrization_from_stics_example_files(species, variety)
station.CODECALTEMP = 1
weather['wind'] = np.abs(np.round(np.random.normal(2.5, 1, weather.shape[0]), 1))
sns.histplot(data=weather, x='wind')
# %%
# Run the simulation
pystics_df1, pystics_mat_list = run_pystics_simulation(weather, crop, soil, constants, manage, station, initial)
# %%
station.CODECALTEMP = 2
pystics_df2, pystics_mat_list = run_pystics_simulation(weather, crop, soil, constants, manage, station, initial)
# %%
# pystics_df['']
# plot_df = pd.concat([pystics_df1.tcult,pystics_df2.tcult])
sns.kdeplot(data=pystics_df1, x='tcult', label='codecaltemp = 1')
sns.kdeplot(data=pystics_df2, x='tcult', label='codecaltemp = 2')
plt.legend()
# %%
pystics_df1['Date'] = pd.to_datetime(pystics_df1.year.astype(int) * 1000 + pystics_df1.doy, format='%Y%j')
pystics_df2['Date'] = pd.to_datetime(pystics_df2.year.astype(int) * 1000 + pystics_df2.doy, format='%Y%j')
sns.lineplot(data=pystics_df1, x='Date', y='tcult')
sns.lineplot(data=pystics_df2, x='Date', y='tcult')
# %%
###################################
############### old ###############
###################################
# Yield (t.ha-1)
pystics_df.mafruit.max()
# %%
pystics_df.mafruit
# %%
import os
import pandas as pd
parametrization_files_path = os.path.abspath('../src/pystics/parametrization_files')
# %%
obs = pd.read_csv(parametrization_files_path + f"/example/wheat.obs", sep=';', na_values='-999.99')
obs['Date'] = pd.to_datetime(obs.ian * 1000 + obs.jul, format='%Y%j')
obs_dm = obs[['Date','masec(n)']]
# %%
pystics_df['Date'] = pd.to_datetime(pystics_df.year.astype(int) * 1000 + pystics_df.doy, format='%Y%j')
sim_dm = pystics_df[['Date','masec']]
# %%
df = pd.merge(sim_dm,obs_dm,on='Date',how='left')
df.drop(df.index[-1],axis=0, inplace=True)
# df_plot = df.dropna().reset_index(drop=True)
# %%
import seaborn as sns
sns.scatterplot(data=df, x='masec',y='masec(n)')
# %%
sns.lineplot(data=df,x='Date',y='masec')
sns.scatterplot(data=df,x='Date',y='masec(n)')
# %%
