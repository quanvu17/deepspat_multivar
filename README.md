# Modeling Nonstationary and Asymmetric Multivariate Spatial Covariances via Deformations

<img align="right" src="https://github.com/quanvu17/deepspat_multivar/blob/master/figures/simulation_asymm_comparison_plot.png" alt="drawing" width="400"/>

This repository provides reproducible code for the manuscript titled *Modeling Nonstationary and Asymmetric Multivariate Spatial Covariances via Deformations* by Quan Vu, Andrew Zammit-Mangion, and Noel Cressie. The manuscript describes a new class of nonstationary and asymmetric multivariate spatial covariance models that are constructed via deformations. Specifically, nonstationary and asymmetric covariances on a geographic domain are modeled as simpler, more familiar, stationary and symmetric covariances on a warped domain through deep injective warping functions.

## Instructions

To reproduce the results and the plots in the manuscript, please download this repository.

Ensure the required packages are installed. The required packages for these scripts are `data.table`, `deepspat`, `dplyr`, `fields`, `ggplot2`, `gridExtra`, `ncdf4`, `sp`, and `verification`. The `deepspat` package can be installed from [this github repository](https://github.com/andrewzm/deepspat), while other packages can be installed using the function `install.packages` in R.

Once the repository is downloaded, run through the script files in the folder `scripts/` to reproduce either numerical results or figures in the manuscript. The data and results used in the manuscript are saved in the folders `data/` and `results/` while the figures used in the manuscript are saved in the folder `figures/`.

## Abstract

Multivariate spatial-statistical models are useful for modeling environmental and socio-demographic processes. The most commonly used models for multivariate spatial covariances assume both stationarity and symmetry for the cross-covariances, but these assumptions are rarely tenable in practice. In this article we introduce a new and highly flexible class of nonstationary and asymmetric multivariate spatial covariance models that are constructed by modeling the simpler and more familiar stationary and symmetric multivariate covariances on a warped domain. Inspired by recent developments in the univariate case, we propose modeling these deformation functions as the composition of a number of simple injective warping functions in a deep-learning framework. Importantly, covariance model validity is guaranteed by construction. We establish the types of warpings that allow for symmetry and asymmetry, and we use likelihood-based methods for inference that are computationally efficient. The utility of this new class of models is shown through various data illustrations, including a simulation study on nonstationary data and an application on ocean temperatures at two different depths.
