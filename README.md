# pySTICS v1.0
<img width="740" alt="OMBREA(16)" src="https://github.com/OmbreaPV/pySTICS/assets/105670904/d68a7c73-4bb7-4a15-8385-dd82508ce496">


pySTICS is an open-source collaborative Python implementation of a simplified version of [STICS v10.0 crop model](https://stics.inrae.fr/). This project is supported by Ombrea (https://www.ombrea.fr/), the international agri-energy centre of expertise of Total Energies, and the Digital Energy and Agriculture Lab, from the University of Liege.


pySTICS is being developed and tested in the context of agronomic impact modelling of agrivoltaic projects, but is even more valid for openfield plots as STICS model was designed and calibrated in such conditions.

For any inquiries or collaboration opportunity, feel free to reach:
- Ombrea : Etienne Perez - eperez@ombrea.fr
- Digital Energy and Agricultural Lab : Roxane Bruhwyler - roxane.bruhwyler@uliege.be


# version 1.0
pySTICS v1.0 includes STICS v10.0 main modules for determinate growth plants and has been validated by comparing its outputs to STICS (with the same options) outputs on many soil and weather data for wheat crop. Other annual plants and forage crops available in STICS are currently being tested with pySTICS, updates about validated species will be released here as soon as the tests are finalized.

The figure below shows the different modules implemented and simplifications made. For the modules description, please refer to STICS documentation ([Beaudoin N. et al.2022. STICS soil-crop model. Conceptual framework, equations and use. Editions Qae.](docs/source/_static/STICS%20soil-crop%20mode,%20Conceptual%20framework,%20equations%20and%20uses.html)). 'Full version' means that the module has been implemented as described in STICS documentation without major simplifying approach. 

![](https://github.com/OmbreaPV/pySTICS/blob/2b419b0c92a3789dcfa86b6fd7cb018d5160caf2/docs/source/_static/table_modules.png)


# Notebook examples
Two notebook examples are available to understand how to run your simulations:
- [parametrization_methods.ipynb](notebooks/parametrization_methods.ipynb) shows how to parametrize pySTICS inputs : weather data, crop parameters, soil parameters, technical parameters, station parameters and constants.
- [run_pystics_simulation_stics_example.ipynb](notebooks/run_pystics_simulation_stics_example.ipynb) shows how to run a simulation on parametrization files provided by STICS project team.


# Installation
1. Create and activate a conda environment (or use your prefered environment manager):
```
conda create --name pystics python=3.11 pip
conda activate pystics
```
2. Install from PyPI
```
pip install pySTICS
```
2bis. Or clone the repository to install it in development mode:
```
git clone git@github.com:OmbreaPV/pySTICS.git
cd pySTICS
pip install -e .
```
# Perspectives
pySTICS is an open-source collaborative project born from a statement : sharing efforts and expertise is the more efficient way to build a reliable crop model capable to assess the agronomic impact of agrivoltaic projects, and thus build virtuous projects. With this in mind, some developments are still needed to model the complex interactions between the plant and its environment under panels shading, particularly regarding intraday physiological responses to fluctuating microclimate.

The very next developpements will focus on including other annual and forage crops, multi-year modelling and model validation with field data. Feel free to reach out if you want to take part of the project !
