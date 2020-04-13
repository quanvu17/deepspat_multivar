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

# Load source
source("scripts/utils.R")

# Simulate data for the experiment
# Simulate warped location
set.seed(1)
r1 <- 50
a1 <- 1
a2 <- 0
a3 <- 1
a4 <- 1 + 1*1i
layers <- c(AWU(r = r1, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = r1, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(res = 1L),
            LFT(a=c(a1,a2,a3,a4)))
layers_asym <- c(AFF_2D(a=c(0.1, 1, 0, -0.1, 0, 1)))

nlayers <- length(layers)
nlayers_asym <- length(layers_asym)

eta <- list()
eta[[1]] <- sin(seq(0, pi, length.out = r1))
eta[[2]] <- c(1, rep(0, r1-1))
for(j in 3:(nlayers-1)) eta[[j]] <- runif(n = 1, min = -1, max = exp(3/2)/2)
eta[[nlayers]] <- 1

length.out <- 101
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped1 <- swarped2 <- s

for(j in 1: (nlayers_asym)){
  swarped1 <- swarped1
  swarped2 <- layers_asym[[j]]$fR(swarped2, eta[[j]])
  swarped1[,1] <- (swarped1[,1] - min(c(swarped1[,1], swarped2[,1]))) / (max(c(swarped1[,1], swarped2[,1])) - min(c(swarped1[,1], swarped2[,1]))) - 0.5
  swarped1[,2] <- (swarped1[,2] - min(c(swarped1[,2], swarped2[,2]))) / (max(c(swarped1[,2], swarped2[,2])) - min(c(swarped1[,2], swarped2[,2]))) - 0.5
  swarped2[,1] <- (swarped2[,1] - min(c(swarped1[,1], swarped2[,1]))) / (max(c(swarped1[,1], swarped2[,1])) - min(c(swarped1[,1], swarped2[,1]))) - 0.5
  swarped2[,2] <- (swarped2[,2] - min(c(swarped1[,2], swarped2[,2]))) / (max(c(swarped1[,2], swarped2[,2])) - min(c(swarped1[,2], swarped2[,2]))) - 0.5
}

for(j in 1: (nlayers)){
  swarped1 <- layers[[j]]$fR(swarped1, eta[[j]]) %>% scal_0_5_mat()
  swarped2 <- layers[[j]]$fR(swarped2, eta[[j]]) %>% scal_0_5_mat()
}

# Simulate symmetric, stationary covariance on the warped domain
D11 <- fields::rdist(swarped1)
D22 <- fields::rdist(swarped2)
D12 <- fields::rdist(swarped1, swarped2)

C11 <- fields::Matern(D11, range=0.2, nu=0.5, phi=1)
C12 <- fields::Matern(D12, range=0.2, nu=0.75, phi=0.8*0.9*1)
C22 <- fields::Matern(D22, range=0.2, nu=1, phi=0.9^2)
C <- rbind(cbind(C11, C12), cbind(t(C12), C22))
K <- t(chol(C))
i <- 1
set.seed(i)
y <- K %*% rnorm(nrow(s)*2)
y1 <- y[1: nrow(s),]
y2 <- y[(nrow(s)+1) : (nrow(s)*2),]
z1 <- y1 + 0.2*rnorm(length(y1))
z2 <- y2 + 0.1*rnorm(length(y2))
z <- c(z1, z2)
df <- data.frame(s, y1, y2, z1, z2)
names(df) <- c("s1", "s2", "y1", "y2", "z1", "z2")

# Save data
save(df, file="results/simulation_asymm_dataset.rda")

# Sample a subset of data for the experiment
RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
df2 <- df[sam2,]

# Five fold Cross validation experiment
RNGkind(sample.kind = "Rounding")
set.seed(1)
all_data <- df2
sam <- sample(1:nrow(all_data), nrow(all_data))
all_data <- all_data[sam,]
groups_data <- split(all_data, (seq(nrow(all_data))-1) %/% (nrow(all_data)/5))

rmse1_model1 <- rmse2_model1 <- crps1_model1 <- crps2_model1 <- rmse1_model2 <- rmse2_model2 <- crps1_model2 <- crps2_model2 <-
  rmse1_model3 <- rmse2_model3 <- crps1_model3 <- crps2_model3 <- rmse1_model4 <- rmse2_model4 <- crps1_model4 <- crps2_model4 <-
  rep(0,5)

layers <- c(AWU(r = 50, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 50, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())

layers_asym <- c(AFF_2D())

for (i in 1:5){
  test_data <- groups_data[[i]]
  train_data <- setdiff(all_data, test_data)
  df <- dplyr::select(train_data, s1, s2, z1, z2)
  newdata <- test_data
  
  # Fit Model 1 (stationary, symmetric)
  d1 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          family = "matern_stat_symm",
                          method = "REML", nsteps = 150L
  )
  predML1 <- predict(d1, newdata = newdata)
  
  # Fit Model 2 (stationary, asymmetric)
  d2 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          layers_asym = layers_asym,
                          family = "matern_stat_asymm",
                          method = "REML", nsteps = 150L
  )
  predML2 <- predict(d2, newdata = newdata)
  
  # Fit Model 3 (nonstationary, symmetric)
  d3 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          layers = layers,
                          family = "matern_nonstat_symm",
                          method = "REML", nsteps = 150L
  )
  predML3 <- predict(d3, newdata = newdata)
  
  # Fit Model 4 (nonstationary, asymmetric)
  d4 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          layers = layers, layers_asym = layers_asym,
                          family = "matern_nonstat_asymm",
                          method = "REML", nsteps = 150L
  )
  predML4 <- predict(d4, newdata = newdata)
  
  df_pred1 <- data.frame(predML1$df_pred, y1=test_data$y1, y2=test_data$y2)
  df_pred2 <- data.frame(predML2$df_pred, y1=test_data$y1, y2=test_data$y2)
  df_pred3 <- data.frame(predML3$df_pred, y1=test_data$y1, y2=test_data$y2)
  df_pred4 <- data.frame(predML4$df_pred, y1=test_data$y1, y2=test_data$y2)
  
  # Save RMSPE and CRPS
  rmse1_model1[i] <- RMSPE(df_pred1$y1, df_pred1$pred_mean_1)
  rmse2_model1[i] <- RMSPE(df_pred1$y2, df_pred1$pred_mean_2)
  crps1_model1[i] <- CRPS(df_pred1$y1, df_pred1$pred_mean_1, df_pred1$pred_var_1)
  crps2_model1[i] <- CRPS(df_pred1$y2, df_pred1$pred_mean_2, df_pred1$pred_var_2)
  
  rmse1_model2[i] <- RMSPE(df_pred2$y1, df_pred2$pred_mean_1)
  rmse2_model2[i] <- RMSPE(df_pred2$y2, df_pred2$pred_mean_2)
  crps1_model2[i] <- CRPS(df_pred2$y1, df_pred2$pred_mean_1, df_pred2$pred_var_1)
  crps2_model2[i] <- CRPS(df_pred2$y2, df_pred2$pred_mean_2, df_pred2$pred_var_2)
  
  rmse1_model3[i] <- RMSPE(df_pred3$y1, df_pred3$pred_mean_1)
  rmse2_model3[i] <- RMSPE(df_pred3$y2, df_pred3$pred_mean_2)
  crps1_model3[i] <- CRPS(df_pred3$y1, df_pred3$pred_mean_1, df_pred3$pred_var_1)
  crps2_model3[i] <- CRPS(df_pred3$y2, df_pred3$pred_mean_2, df_pred3$pred_var_2)
  
  rmse1_model4[i] <- RMSPE(df_pred4$y1, df_pred4$pred_mean_1)
  rmse2_model4[i] <- RMSPE(df_pred4$y2, df_pred4$pred_mean_2)
  crps1_model4[i] <- CRPS(df_pred4$y1, df_pred4$pred_mean_1, df_pred4$pred_var_1)
  crps2_model4[i] <- CRPS(df_pred4$y2, df_pred4$pred_mean_2, df_pred4$pred_var_2)
  
}

# Save cross validation results
crossvalid <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1,
                   rmse1_model2, crps1_model2, rmse2_model2, crps2_model2,
                   rmse1_model3, crps1_model3, rmse2_model3, crps2_model3,
                   rmse1_model4, crps1_model4, rmse2_model4, crps2_model4
)

save(crossvalid, file="results/simulation_asymm_crossvalidation.rda")

# Cross validation results
matrix(unlist(lapply(crossvalid, mean)), nrow = 4, byrow=T)

# Fit models for plotting the comparison
load("results/simulation_asymm_dataset.rda")
newdata <- df

# Fit Model 2 (stationary, asymmetric)
d2 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df2, g = ~ 1,
                        layers_asym = layers_asym,
                        family = "matern_stat_asymm",
                        method = "REML", nsteps = 150L
)
predML2 <- predict(d2, newdata = newdata)

# Fit Model 3 (nonstationary, symmetric)
d3 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df2, g = ~ 1,
                        layers = layers,
                        family = "matern_nonstat_symm",
                        method = "REML", nsteps = 150L
)
predML3 <- predict(d3, newdata = newdata)

# Fit Model 4 (nonstationary, asymmetric)
d4 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df2, g = ~ 1,
                        layers = layers, layers_asym = layers_asym,
                        family = "matern_nonstat_asymm",
                        method = "REML", nsteps = 150L
)
predML4 <- predict(d4, newdata = newdata)

predML <- list(predML2, predML3, predML4)
save(predML, file = "results/simulation_asymm_predictions.rda")
