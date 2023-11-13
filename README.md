# pySTICS beta version
<img width="740" alt="OMBREA(16)" src="https://github.com/OmbreaPV/pySTICS/assets/105670904/d68a7c73-4bb7-4a15-8385-dd82508ce496">


Python Implementation of the STICS crop model (https://eng-stics.paca.hub.inrae.fr).

# Notebook examples
See Notebooks/run_simulation.ipynb to see a simulation example on wheat.

Two simulation types are possible :
- Simulation from STICS USM examples associated to a species and variety chosen by the user.
- Simulation with plant/soil parameters from USM examples and weather data from a location and a year chosen by the user. Weather data are then requested from ERA5 API (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels).


# How to use pySTICS
1. Create and activate a conda environment:```console
conda create --name envi_name  python=3.11 pip
conda activate envi_name
```
2. Clone repository and install pystics with pip:
```console
git clone git@github.com:OmbreaPV/pySTICS.git
cd pySTICS
pip install -e .
```
