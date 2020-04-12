# Reproducible code for "Modeling Nonstationary and Asymmetric
# Multivariate Spatial Covariances via Deformations"
# Copyright (c) 2020 Quan Vu
# Author: Quan Vu, quanv (at) uow.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# Load packages
library("data.table")
library("deepspat")
library("dplyr")
library("fields")
library("ggplot2")
library("gridExtra")
library("ncdf4")
library("sp")
library("verification")

# RMSPE
RMSPE <- function(true, pred){
  sqrt(mean((true - pred)^2))
}

# CRPS
CRPS <- function(true, pred, pred_var){
  crps(true, cbind(pred, sqrt(pred_var)))$CRPS
}

# Homogenizing function
b1 <- function(sw, k, l, m){
  sw <- sw / sqrt((sw[k,1] - sw[l,1])^2 + (sw[k,2] - sw[l,2])^2)
  sw <- sw - matrix(rep(sw[k,], nrow(sw)), nrow=nrow(sw), byrow=T)
  sw
}

b2 <- function(sw, k, l, m){
  rot <- matrix(c(sw[l,1], sw[l,2], -sw[l,2], sw[l,1])/sqrt(sw[l,1]^2 + sw[l,2]^2), nrow=2, byrow=T)
  sw <- t(rot %*% t(sw))
  sw
}

b3 <- function(sw, k, l, m){
  refl <- matrix(c(1, 0, 0, sign(sw[m,2])), nrow=2)
  sw <- t(refl %*% t(sw))
  sw
}

b0 <- function(sw, k, l, m){
  sw <- b1(sw, k, l, m)
  sw <- b2(sw, k, l, m)
  sw <- b3(sw, k, l, m)
  sw
}