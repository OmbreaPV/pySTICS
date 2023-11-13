import pandas as pd
import numpy as np
import dask
import xarray as xr
import s3fs
from pathlib import Path

from pvlib.iotools.pvgis import get_pvgis_tmy, get_pvgis_hourly
from pvlib.solarposition import sun_rise_set_transit_ephem

from sklearn.ensemble import RandomForestRegressor

@dask.delayed
def s3open(path):
    fs = s3fs.S3FileSystem(anon=True, default_fill_cache=False, 
                           config_kwargs = {'max_pool_connections': 20})
    return s3fs.S3Map(path, s3=fs)

def open_era5_range(start_year, end_year, variables):
    ''' Opens ERA5 monthly Zarr files in S3, given a start and end year (all months loaded) and a list of variables'''
    
    
    file_pattern = 'era5-pds/zarr/{year}/{month}/data/{var}.zarr/' 
    years = list(np.arange(start_year, end_year+1, 1))
    months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
    
    l = []
    for var in variables:
        # print(var)        
        # Get files
        files_mapper = [s3open(file_pattern.format(year=year, month=month, var=var)) for year in years for month in months]
        # Look up correct time dimension by variable name
        if var in ['precipitation_amount_1hour_Accumulation']:
            concat_dim='time1'
        else:
            concat_dim='time0'
        # Lazy load
        ds = xr.open_mfdataset(files_mapper, engine='zarr', 
                               concat_dim=concat_dim, combine='nested', 
                               coords='minimal', compat='override', parallel=True)
        
        # Fix dimension names
        try:
            ds = ds.rename({'time0':'time'})
        except ValueError:
            ds = ds.rename({'time1':'time'})
        #ds = fix_accum_var_dims(ds, var)
        l.append(ds)
        
    ds_out = xr.merge(l)
    
    return ds_out


def generate_meteo(latitude,longitude,start_year=2005,end_year = 2015, save=True):

    meteo = pd.DataFrame()

    # get PVGIS hourly Data

    # print('hourly tmy')
    hourly_meteo = get_pvgis_hourly(latitude,longitude,start_year,end_year,timeout=300,map_variables=False)
    df = hourly_meteo[0]
    df['G(h)'] = df['Gb(i)']+ df['Gd(i)']
    # print('hourly tmy - OK')

    # Simple ML to get RH

    # print ('ML for RH')
    tmy = get_pvgis_tmy(latitude,longitude,map_variables=False)[0]
    features = ['T2m','G(h)', 'WS10m']
    X = tmy[features]
    y = tmy.RH

    clf = RandomForestRegressor(n_estimators=20,min_samples_split=10)
    clf.fit(X,y)
    df['RH'] = clf.predict(df[features])
    # print ('ML for RH - OK')

    # Photoperiod 
    # print ('compute Photoperiod')
    photoperiod = sun_rise_set_transit_ephem(pd.DatetimeIndex(np.unique(df.index.date),tz='UTC'), 44, 4)

    photoperiod['photoperiod'] = (photoperiod.sunset - photoperiod.sunrise).astype("timedelta64[ms]").astype(int) / (1000*3600)
    # print ('compute Photoperiod - OK')

    # Total Precipitation ERA5

    # print('Total precipitation ERA5')
    ds_out = open_era5_range(start_year, end_year, ['precipitation_amount_1hour_Accumulation'])
    ds_out = ds_out.assign_coords(lon=(((ds_out.lon + 180) % 360) - 180)).sortby(['time','lat','lon'])
    total_precipitation = ds_out.sel(lat=latitude, lon=longitude, method='nearest').precipitation_amount_1hour_Accumulation.resample({'time':'1d'}).sum().values
    # print('Total precipitation ERA5 - OK')

    # Daily Agregation for STICS 

    # print ('Daily Agregation')

    meteo['Temp_min'] = df['T2m'].resample('d').min()
    meteo['Temp_moy'] = df['T2m'].resample('d').mean()
    meteo['Temp_max'] = df['T2m'].resample('d').max()
    meteo['Daily_Temp'] = df['T2m'].resample('d').apply(list)
    meteo['Radiation'] = df['G(h)'].resample('d').sum()*36/10000 # convert W/M-2 to MJ/m-2
    meteo['Rel_Hum'] = df['RH'].resample('d').mean()
    meteo['Rel_Hum_min'] = df['RH'].resample('d').min()
    meteo['Rel_Hum_max'] = df['RH'].resample('d').max()
    meteo['Rain'] = total_precipitation * 1000 # convert m.day-1 to mm.day-1
    meteo['Photoperiod'] = photoperiod.photoperiod
    meteo['ANNEE'] = meteo.index.year
    meteo['DOY'] = meteo.index.dayofyear
    meteo['Wind'] = df['WS10m'].resample('d').mean()
    meteo["Date"] = meteo.index.strftime("%d/%m/%Y")
    meteo["latRad"] = (latitude*np.pi)/180 

    # print ('Daily Agregation - OK')

    meteo.reset_index(inplace=True)
    meteo = meteo[["Date","ANNEE", "DOY", "Temp_min", "Temp_moy", "Temp_max", "Daily_Temp", "Radiation", "Rain", "Wind", "Rel_Hum", "Rel_Hum_min","Rel_Hum_max","Photoperiod","latRad"]]

    if save:
        folder = Path('../meteo/tmy')
        folder.mkdir(exist_ok=True)
        file = f"meteo_{latitude}_{longitude}_{start_year}_{end_year}.csv"
        file_path = folder / file
        meteo.to_csv(file_path ,header=True)

    return meteo

def get_cached_meteo(latitude,longitude,start_year=2005,end_year = 2015):
    
    cache_folder = Path('../meteo')
    cache_folder.mkdir(exist_ok=True)
    file = f"meteo_{latitude}_{longitude}_{start_year}_{end_year}.csv"
    file_path = cache_folder / file

    if file_path.exists():
        with file_path.open() as f:
            df = pd.read_csv(f, index_col=0,header= 0)
    else:
        df = generate_meteo(latitude, longitude,start_year,end_year)
        with file_path.open("w") as f:
            df.to_csv(f)
    
    return df
