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

# Simulate dataset
length.out <- 101
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped <- s

D <- fields::rdist(swarped)
C11 <- fields::Matern(D, range=0.2, nu=0.5, phi=1)
C12 <- fields::Matern(D, range=0.2, nu=0.75, phi=0.8*0.9*1)
C22 <- fields::Matern(D, range=0.2, nu=1, phi=0.9^2)
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
save(df, file="results/simulation_misspecified1_dataset.rda")

# Sample a subset of data for the experiment
RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
train_data <- df[sam2,]
test_data <- dplyr::setdiff(df, train_data)

df <- dplyr::select(train_data, s1, s2, z1, z2)
newdata <- dplyr::select(test_data, s1, s2)
groups_data <- split(newdata, (seq(nrow(newdata))-1) %/% (nrow(newdata)/5))

layers <- c(AWU(r = 50, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 50, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())

layers_asym <- c(AFF_2D())

# Fit Model 1 (stationary, symmetric)
t1 <- proc.time()
d1 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                        family = "matern_stat_symm",
                        method = "REML", nsteps = 150L
)
t2 <- proc.time()
time1 <- (t2 - t1)[3]
predML1 <- predict(d1, newdata = newdata)

# Fit Model 2 (nonstationary, asymmetric)
t1 <- proc.time()
d2 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                        layers = layers, layers_asym = layers_asym,
                        family = "matern_nonstat_asymm",
                        method = "REML", nsteps = 150L
)
t2 <- proc.time()
time2 <- (t2 - t1)[3]
predML2 <- predict(d2, newdata = newdata)

# Save parameter estimates
eta <- d2$run(d2$eta_tf)
AFF_pars <- d2$run(d2$layers_asym[[1]]$pars)
LFT_pars <- d2$run(d2$layers[[12]]$pars)
scalings <- d2$run(d2$scalings)
scalings_asym <- d2$run(d2$scalings_asym)
nu_1 <- d2$run(d2$nu_tf_1)
nu_2 <- d2$run(d2$nu_tf_2)
sigma2_1 <- d2$run(d2$sigma2_tf_1)
sigma2_2 <- d2$run(d2$sigma2_tf_2)
sigma2_12 <- d2$run(d2$sigma2_tf_12)
l <- as.numeric(d2$run(d2$l_tf_1))
precy_1 <- d2$run(d2$precy_tf_1)
precy_2 <- d2$run(d2$precy_tf_2)
s_warped1 <- d2$run(d2$swarped_tf1)
s_warped2 <- d2$run(d2$swarped_tf2)
beta <- d2$run(d2$beta)
parameter_est <- list(eta, AFF_pars, LFT_pars, scalings, scalings_asym,
                      nu_1, nu_2, sigma2_1, sigma2_2, sigma2_12, l, precy_1, precy_2,
                      s_warped1, s_warped2, beta)
save(parameter_est, file = paste0("results/simulation_misspecified1_parameter_est.rda"))

df_pred1 <- data.frame(predML1$df_pred, y1=test_data$y1, y2=test_data$y2)
df_pred2 <- data.frame(predML2$df_pred, y1=test_data$y1, y2=test_data$y2)

# Save RMSPE and CRPS
rmse1_model1 <- RMSPE(df_pred1$y1, df_pred1$pred_mean_1)
rmse2_model1 <- RMSPE(df_pred1$y2, df_pred1$pred_mean_2)
crps1_model1 <- CRPS(df_pred1$y1, df_pred1$pred_mean_1, df_pred1$pred_var_1)
crps2_model1 <- CRPS(df_pred1$y2, df_pred1$pred_mean_2, df_pred1$pred_var_2)

rmse1_model2 <- RMSPE(df_pred2$y1, df_pred2$pred_mean_1)
rmse2_model2 <- RMSPE(df_pred2$y2, df_pred2$pred_mean_2)
crps1_model2 <- CRPS(df_pred2$y1, df_pred2$pred_mean_1, df_pred2$pred_var_1)
crps2_model2 <- CRPS(df_pred2$y2, df_pred2$pred_mean_2, df_pred2$pred_var_2)

Cost1 <- d1$run(d1$Cost)
Cost2 <- d2$run(d2$Cost)

# Save cross validation results
validation <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1, Cost1, time1,
                   rmse1_model2, crps1_model2, rmse2_model2, crps2_model2, Cost2, time2
)

save(validation, file=paste0("results/simulation_misspecified1_validation.rda"))

# Cross validation results
matrix(unlist(lapply(validation, mean)), nrow = 2, byrow=T)
