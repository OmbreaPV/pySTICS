# pySTICS beta version
<img width="740" alt="OMBREA(16)" src="https://github.com/OmbreaPV/pySTICS/assets/105670904/d68a7c73-4bb7-4a15-8385-dd82508ce496">


PySTICS is a Python implementation of a simplified version of [STICS v10.0 crop model](https://eng-stics.paca.hub.inrae.fr).
Full-version or simplified modules are currently being implemented with a focus on determinate crops and water deficit effect.
![](https://github.com/OmbreaPV/pySTICS/blob/f29772c30cadb4847ad3e43852a780007e2bd2e1/docs/source/_static/table_progress.png)

# Notebook examples
See [run_simulation.ipynb](Notebooks/run_simulation.ipynb) to see a simulation example on wheat.

Two simulation types are possible :
- Simulation from STICS USM examples associated to a species and variety chosen by the user.
- Simulation with plant/soil parameters from USM examples and weather data from a location and a year chosen by the user. Weather data are then requested from ERA5 API (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels).


# How to use pySTICS
1. Create and activate a conda environment (or use your prefered environment manager):
```
conda create --name pystics python=3.11 pip
conda activate pystics
```
2. Clone repository and install pystics with pip (here in development mode):
```
git clone git@github.com:OmbreaPV/pySTICS.git
cd pySTICS
pip install -e .
```
