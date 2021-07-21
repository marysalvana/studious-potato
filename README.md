# studious-potato
consolidated thesis files

## Table of Contents in Dissertation
 Chapter 1 - Introduction
<br> Chapter 2 - Lagrangian Spatio-Temporal Nonstationary Covariance Function
<br> Chapter 3 - Lagrangian Spatio-Temporal Nonstationary Cross-Covariance Function
<br> Chapter 4 - Lagrangian Spatio-Temporal Cross-Covariance Function with Multiple Advections
<br> Chapter 5 - Nonstationary Taylor's Hypothesis, Multivariate ExaGeoStat, and Spatio-Temporal ExaGeoStat

## Table of Contents in Codes
chapter2 - Lagrangian Spatio-Temporal Nonstationary Covariance Function
<br> chapter3 - Lagrangian Spatio-Temporal Nonstationary Cross-Covariance Function
<br> chapter4 - Lagrangian Spatio-Temporal Cross-Covariance Function with Multiple Advections
<br> chapter5 - Nonstationary Taylor's Hypothesis
<br> chapter6 - Multivariate ExaGeoStat
<br> chapter7 - Spatio-Temporal ExaGeoStat

## Codes

- chapter0-thesis-presentation-plots.R -- codes to plot figures for thesis defense presentation beamer
- chapter0-thesis-presentation.R -- codes to produce covariance and realizations for thesis defense presentation
<br> <br>
- chapter2-3.R -- codes for thin plate splines figure
- chapter2-4.R -- codes for simulation study
- chapter2-5-plots.R -- codes to plot space-time images in manuscript
- chapter2-5-estimation.R -- codes to fit models on real data 
<br> <br>
- chapter5-2.R -- codes for simulating theoretical and empirical Taylor's hypothesis test functions
- chapter5-3.R -- codes for simulating Lagrangian nonstationary covariance (spatially varying and deformation models) and nonparametric estimation of the advection velocity parameters
- chapter5-4.R -- codes for real data application
<br> <br>
- chapter7-1-synthetic-data.R -- codes for plotting synthetic data for space-time ExaGeoStat
- chapter7-2-extract-data-from-ncdf.R -- codes for extracting data from ncdf
- chapter7-3-real-data.R -- codes for plotting real data and processing results from Shaheen 
<br> <br>
- chapter8-1-synthetic-data.R -- codes for plotting synthetic data for ExaGeoStat Lagrangian

## Packages

library(mvnfast)
<br> library(future.apply)


## Offline Data Links

SC21 Netcdf SAUDI Data / Taylor's Hypothesis
<br> /yourlainess/phd/data/sc21/SAUDI

SC21 Netcdf US Data
<br> /yourlainess/phd/data/sc21/US

## Download Links:

pm2.5: https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary
<br> wind: https://disc.gsfc.nasa.gov/datasets/M2I1NXASM_5.12.4/summary
<br> pm2.5 3D: https://disc.gsfc.nasa.gov/datasets/M2I3NVAER_5.12.4/summary
<br> wind 3D: https://disc.gsfc.nasa.gov/datasets/M2I3NPASM_5.12.4/summary
<br> covariates: https://disc.gsfc.nasa.gov/datasets/M2I3NVASM_5.12.4/summary

subregion for saudi land and ocean region

Spatial subset: 26.719,5.625,85.078,42.188



## Note:

Run the Rscript inside the folder R_codes


