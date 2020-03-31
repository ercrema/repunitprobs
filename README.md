# Modifiable reporting unit problems and time series of long-term human activity: source code, data, and scripts

This repository contains all data and scripts required to fully reproduce figures and simulations in the manuscript "Modifiable reporting unit problems and time series of long-term human activity" authored by Crema, E.R., and Bevan, A.

## Simulation Scripts
R scripts for running the simulations are stored in the file [simulation_log.R](https://github.com/ercrema/repunitprobs/blob/master/simulation_log.R). This requires some additional functions stored in [/src/utility.R](https://github.com/ercrema/repunitprobs/blob/master/src/utility.R). Simulation outputs required for the figures are stored as R image files in the folder [R_Images](https://github.com/ercrema/repunitprobs/tree/master/R_Images).

## Figures and Data
R scripts for generating the figures in the manuscript are contained in the file [figures/figurelog.R](https://github.com/ercrema/repunitprobs/blob/master/figures/figurelog.R). Raw data required for figure 1 are stored in the folder [data](https://github.com/ercrema/repunitprobs/tree/master/data). 

## ESM
The electronic supplementary material ([ESM.pdf](https://github.com/ercrema/repunitprobs/blob/master/ESM.pfg)) contains additional details about the three simulations. This is just a rendered version of the Rmardkdown file [ESM.Rmd](https://github.com/ercrema/repunitprobs/blob/master/ESM.Rmd).

## Required R Packages & Software

### R Packages
* [rcarbon 1.3.2](https://github.com/ahb108/rcarbon) (this is the development version: To install: `devtools::install_github('ahb108/rcarbon@5bb28cf')`)
* [trapezoid 2.0.0](https://cran.r-project.org/web/packages/trapezoid/index.html)
* [oxcAAR 1.0.0](https://cran.r-project.org/web/packages/oxcAAR/index.html)

### Software 
* [OxCal 4.3](https://c14.arch.ox.ac.uk/oxcal.html) 
