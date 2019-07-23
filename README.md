# LRT

Code for [Higher-Order Income Dynamics with Linked Regression Trees](http://web.econ.ku.dk/druedahl/papers/2019_LRT.pdf), [Druedahl](http://web.econ.ku.dk/druedahl) and [Munk-Nielsen](http://web.econ.ku.dk/munk-nielsen), 2019).

## Requirements

1. SAS (Used: SAS 9.4)
2. Stata (Used: Stata 15)
3. MATLAB + C++ compiler (Used: MATLAB 2019a with Microsof Visual Studio 2017 Community Edition)
4. Python 3.7 (Used: Anaconda 2019.03)

## Data inputs

1. Danish register data (INDUPD, BEF, FAIN, IDAP, with views located in H:/Rawdata/706248/views/)
2. Consumer price data: data/PRIS61.xlsx (from Statistikbanken.dk)
3. House price data: data/EJEN6.xlsx (from Statistikbanken.dk)

## ReadMe

Everything can be run from **main.ipynb.**.

It calls:

1. SAS to fetch the data (0_*.sas)
2. Stata to structure the data (1_*.do)
3. MATLAB to run ABB estimator (ABB/run.m)
4. Python to run all other estimations and plot the results (2_*.ipynb)
5. MATLAB to solve and simulate the consumption-saving model (ConSav/run.m)

All results are saved in **output/**.

**Paths:**

It is necessary to adjust the following paths:

1. libname raw in 0_1_import_datasets.sas should point to the views for your project.
2. libname data in the 0_*.sas files should point to data/ in this directory.
3. The cd's in the 1_*.do files should be to this directory.