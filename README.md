# pySTICS beta version
<img width="740" alt="OMBREA(16)" src="https://github.com/OmbreaPV/pySTICS/assets/105670904/d68a7c73-4bb7-4a15-8385-dd82508ce496">


Python Implementation of the STICS crop model (https://eng-stics.paca.hub.inrae.fr).

# Notebook examples
See Notebooks/run_simulation.ipynb to see a simulation example on wheat.

Two simulation types are possible :
- Simulation from STICS USM examples associated to a species and veriety chosen by the user.
- Simulation with plant/soil parameters from USM examples and weather data from a location and a year chosen by the user. Weather data are then requested from ERA5 API (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels).


# How to use pySTICS
## Clone gitlab repository and install pystics in a conda environment

Complete procedure to install pystics with conda environment `envi_name`:
```console
conda create --name envi_name  python=3.11 pip
conda activate envi_name
python -m pip install -U pip setuptools wheel
```
Clone gitlab repository in local and be in this repository to install pystics with pip:
```console
git clone git@gitlab.com:ombrea/ministics_uliege.git
cd ministics_uliege
pip install -e .
```
