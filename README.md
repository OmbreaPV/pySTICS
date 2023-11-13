# pySTICS
<img width="740" alt="OMBREA(16)" src="https://github.com/OmbreaPV/pySTICS/assets/105670904/d68a7c73-4bb7-4a15-8385-dd82508ce496">


Python Implementation of the STICS crop model 


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
