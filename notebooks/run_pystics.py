# %%
import pystics  
# %%
##### Choice of species and variety --> the USM associated to the species 
species = 'wheat'
variety = 'Talent'
# %%
##### Run a simulation

# Read example input files from STICS for the USM associated to chosen species and variety
from pystics.params import parametrization_from_stics_example_files
weather, crop, manage, soil, station, constants, initial = parametrization_from_stics_example_files(species, variety)
station
# %%
# Run the simulation
from pystics.simulation import run_pystics_simulation
pystics_df, pystics_mat_list = run_pystics_simulation(weather, crop, soil, constants, manage, station, initial)
# %%
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
