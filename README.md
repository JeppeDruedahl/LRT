# LRT

Code for [Higher-Order Income Dynamics with Linked Regression Trees](https://drive.google.com/open?id=18j-uGjIBF3YzQolRQgKUmkSg3JsPiqjp), [Druedahl](http://web.econ.ku.dk/druedahl) and [Munk-Nielsen](http://web.econ.ku.dk/munk-nielsen), 2020).

## Requirements

1. SAS (Used: SAS 9.4)
2. Stata (Used: Stata 15)
3. MATLAB + C++ compiler (Used: MATLAB 2019a with Microsof Visual Studio 2017 Community Edition)
4. Python 3.7 (Used: Anaconda 2019.10)

## ReadMe

Everything can be run from **main.ipynb.**.

It calls:

1. SAS to fetch the data (0_*.sas)
2. Stata to structure the data (1_*.do)
3. MATLAB to run ABB estimator (ABB/run.m)
4. Python to run all other estimations and plot the results (2_*.ipynb)
5. MATLAB to solve and simulate the consumption-saving model (ConSav/run.m)

All results are saved in **output/**.

The code can run in either **online-mode** or **offline-mode** (default). The mode is changed in the top of **main.ipynb**.

**Online-mode:** Produces the exact same results as shown in the paper, but can only be run on the Danmark Statistics severs in a project with the relvant data access (INDUPD, BEF, FAIN, IDAP). When running in online-mode the following paths should be adjusted:

1. libname raw in 0_1_import_datasets.sas should point to the views for your project (default: H:/Rawdata/706248/views/).
2. libname data in the 0_*.sas files should point to data/ in this directory (default:K:/workdata/706248/LRT/).
3. The cd's in the 1_*.do files should be to this directory (default: K:/workdata/706248\LRT/).

**Offline-mode:** This option is used to run the code when not having access to the register data. This implies:

1. ABB results are only simulated using the ABB process estimated on register data (saved in ABB/read_par_estimates.m)
2. Data is simulated from the LRT, depth = 6, process estimated on the register data + classical measurement error (see details replication.ibynb)
3. All other estimators are applied to the simulated data
4. The same tables and figures are produced as in online-mode.

## Data inputs

1. Danish register data (INDUPD, BEF, FAIN, IDAP)
2. Consumer price data: data/PRIS61.xlsx (from Statistikbanken.dk)
3. House price data: data/EJEN6.xlsx (from Statistikbanken.dk)

## Misc.

1. **ABB/:** Contains code for the ABB estimator
2. **censored_estimates/:** Contains censored estimates
3. **ConSav/:** Contains code for the consumption-saving model
4. **LRT/:** Contains code for the LRT estimator
5. **data/:** Contains raw data and some intermediate results
6. **dst/:** Contains code for taking home code and output from the DST servers and unpacking them