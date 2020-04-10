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
nlayers <- length(layers)
eta <- list()
eta[[1]] <- sin(seq(0, pi, length.out = r1))
eta[[2]] <- c(1, rep(0, r1-1))
for(j in 3:(nlayers-1)) eta[[j]] <- runif(n = 1, min = -1, max = exp(3/2)/2)
eta[[nlayers]] <- 1

length.out <- 101
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped <- s
for(j in 1: (nlayers))
  swarped <- layers[[j]]$fR(swarped, eta[[j]]) %>% scal_0_5_mat()

# Simulate symmetric, stationary covariance on the warped domain
D <- fields::rdist(swarped)
C11 <- fields::Matern(D, range=0.1, nu=0.5, phi=1)
C12 <- fields::Matern(D, range=0.1, nu=1, phi=0.45*0.9*1)
C22 <- fields::Matern(D, range=0.1, nu=1.5, phi=0.9^2)
C <- rbind(cbind(C11, C12), cbind(t(C12), C22))
K <- t(chol(C))
i <- 17
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
save(df, file="results/simulation_symm_dataset.rda")

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

rmse1_model1 <- rmse2_model1 <- crps1_model1 <- crps2_model1 <-
  rmse1_model2 <- rmse2_model2 <- crps1_model2 <- crps2_model2 <- rep(0,5)

layers <- c(AWU(r = 50L, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 50L, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())

for (i in 1:5){
  test_data <- groups_data[[i]]
  train_data <- setdiff(all_data, test_data)
  
  df <- dplyr::select(train_data, s1, s2, z1, z2)
  newdata <- dplyr::select(test_data, s1, s2)
  
  # Fit Model 1 (stationary, symmetric)
  d1 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          family = "matern_stat_symm",
                          method = "REML", nsteps = 150L
  )
  predML1 <- predict.deepspat_multivar(d1, newdata = newdata)
  
  # Fit Model 2 (nonstationary, symmetric)
  d2 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          layers = layers,
                          family = "matern_nonstat_symm",
                          method = "REML", nsteps = 150L
  )
  predML2 <- predict.deepspat_multivar(d2, newdata = newdata)
  
  df_pred1 <- data.frame(predML1$df_pred, y1=test_data$y1, y2=test_data$y2)
  df_pred2 <- data.frame(predML2$df_pred, y1=test_data$y1, y2=test_data$y2)
  
  # Save RMSPE and CRPS
  rmse1_model1[i] <- RMSPE(df_pred1$y1, df_pred1$pred_mean_1)
  rmse2_model1[i] <- RMSPE(df_pred1$y2, df_pred1$pred_mean_2)
  crps1_model1[i] <- CRPS(df_pred1$y1, df_pred1$pred_mean_1, df_pred1$pred_var_1)
  crps2_model1[i] <- CRPS(df_pred1$y2, df_pred1$pred_mean_2, df_pred1$pred_var_2)
  
  rmse1_model2[i] <- RMSPE(df_pred2$y1, df_pred2$pred_mean_1)
  rmse2_model2[i] <- RMSPE(df_pred2$y2, df_pred2$pred_mean_2)
  crps1_model2[i] <- CRPS(df_pred2$y1, df_pred2$pred_mean_1, df_pred2$pred_var_1)
  crps2_model2[i] <- CRPS(df_pred2$y2, df_pred2$pred_mean_2, df_pred2$pred_var_2)
  
}

# Save cross validation results
crossvalid <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1,
                   rmse1_model2, crps1_model2, rmse2_model2, crps2_model2
)

save(crossvalid, file="results/simulation_symm_crossvalidation.rda")

# Cross validation results
matrix(unlist(lapply(crossvalid, mean)), nrow = 2, byrow=T)

# Plotting the comparison
load("results/simulation_symm_dataset.rda")

newdata <- df

# Fit Models

# Fit Model 1 (stationary, symmetric)
d1 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df2, g = ~ 1,
                        family = "matern_stat_symm",
                        method = "REML", nsteps = 150L
)
predML1 <- predict.deepspat_multivar(d1, newdata = newdata)

# Fit Model 2 (nonstationary, symmetric)
d2 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df2, g = ~ 1,
                        layers = layers,
                        family = "matern_nonstat_symm",
                        method = "REML", nsteps = 150L
)
# Save parameter estimates
eta <- d2$run(d2$eta_tf)
LFT_pars <- d2$run(d2$layers[[12]]$pars)
scalings <- d2$run(d2$scalings)
nu_1 <- d2$run(d2$nu_tf_1)
nu_2 <- d2$run(d2$nu_tf_2)
sigma2_1 <- d2$run(d2$sigma2_tf_1)
sigma2_2 <- d2$run(d2$sigma2_tf_2)
sigma2_12 <- d2$run(d2$sigma2_tf_12)
l <- as.numeric(d2$run(d2$l_tf_1))
precy_1 <- d2$run(d2$precy_tf_1)
precy_2 <- d2$run(d2$precy_tf_2)
s_warped <- d2$run(d2$swarped_tf1)
beta <- d2$run(d2$beta)
negcost <- d2$negcost
parameter_est <- list(eta,LFT_pars,scalings,nu_1,nu_2,sigma2_1,sigma2_2,sigma2_12,l,precy_1,precy_2,s_warped,beta,negcost)
save(parameter_est, file = paste0("results/simulation_symm_parameter_estimate.rda"))

predML2 <- predict.deepspat_multivar(d2, newdata = newdata)

df_pred1 <- data.frame(predML1$df_pred, y1=df$y1, y2=df$y2)
df_pred2 <- data.frame(predML2$df_pred, y1=df$y1, y2=df$y2)

# Gap experiment
df1 <- df[df$s1 > -0.28 & df$s1 < -0.08 & df$s2 > -0.48 & df$s2 < -0.28,]
df0 <- setdiff(df,df1)
RNGkind(sample.kind = "Rounding")
set.seed(2)
sam <- sample(1:nrow(df0), 1000)
df3 <- df0[sam,]

# Fit Models
newdata <- df

# Fit Model 1 (stationary, symmetric)
d1 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df3, g = ~ 1,
                        family = "matern_stat_symm",
                        method = "REML", nsteps = 150L
)
predML3 <- predict.deepspat_multivar(d1, newdata = newdata)

# Fit Model 2 (nonstationary, symmetric)
d2 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df3, g = ~ 1,
                        layers = layers,
                        family = "matern_nonstat_symm",
                        method = "REML", nsteps = 150L
)
predML4 <- predict.deepspat_multivar(d2, newdata = newdata)

df_pred3 <- data.frame(predML3$df_pred, y1=df$y1, y2=df$y2)
df_pred4 <- data.frame(predML4$df_pred, y1=df$y1, y2=df$y2)

predML3.0 <- predML3$df_pred[predML3$df_pred$s1 > -0.28 & predML3$df_pred$s1 < -0.08 & predML3$df_pred$s2 > -0.48 & predML3$df_pred$s2 < -0.28,]
predML4.0 <- predML4$df_pred[predML4$df_pred$s1 > -0.28 & predML4$df_pred$s1 < -0.08 & predML4$df_pred$s2 > -0.48 & predML4$df_pred$s2 < -0.28,]

# Save RMSPE and CRPS for the gap experiment
rmse1_model1 <- RMSPE(df1$y1, predML3.0$pred_mean_1)
rmse2_model1 <- RMSPE(df1$y2, predML3.0$pred_mean_2)
crps1_model1 <- CRPS(df1$y1, predML3.0$pred_mean_1, predML3.0$pred_var_1)
crps2_model1 <- CRPS(df1$y2, predML3.0$pred_mean_2, predML3.0$pred_var_2)

rmse1_model2 <- RMSPE(df1$y1, predML4.0$pred_mean_1)
rmse2_model2 <- RMSPE(df1$y2, predML4.0$pred_mean_2)
crps1_model2 <- CRPS(df1$y1, predML4.0$pred_mean_1, predML4.0$pred_var_1)
crps2_model2 <- CRPS(df1$y2, predML4.0$pred_mean_2, predML4.0$pred_var_2)

block_holdout <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1,
                      rmse1_model2, crps1_model2, rmse2_model2, crps2_model2
)

save(block_holdout, file="results/simulation_symm_block_holdout.rda")

# Save prediction results for plotting
predML <- list(predML1, predML2, predML3, predML4)
save(predML, file = "results/simulation_symm_predictions.rda")

# Bootstrap
# Estimated covariance matrix 
D <- fields::rdist(s_warped)
C11 <- fields::Matern(D, range=l, nu=nu_1, phi=sigma2_1)
C12 <- fields::Matern(D, range=l, nu=(nu_1+nu_2)/2, phi=sigma2_12)
C22 <- fields::Matern(D, range=l, nu=nu_2, phi=sigma2_2)
C <- rbind(cbind(C11, C12), cbind(t(C12), C22))
C0 <- rbind(cbind(1/precy_1 * diag(nrow(df2)), matrix(rep(0, nrow(df2)*nrow(df2)), nrow=nrow(df2))),
            cbind(matrix(rep(0, nrow(df2)*nrow(df2)), nrow=nrow(df2)), 1/precy_2 * diag(nrow(df2))))
K <- t(chol(C + C0))
K_inv <- solve(K)
z0 <- K_inv %*% c(df2$z1 - beta[1], df2$z2 - beta[2])
s <- s[sam2,]
# Bootstrap sample and find bootstrap estimates
for (i in 1:1000){  # this can be done in parallel
RNGkind(sample.kind = "Rounding")
set.seed(i)
resample_z0 <- sample(z0, length(z0), replace=T)
resample_z <- K %*% resample_z0
z1 <- resample_z[1: nrow(s)] + beta[1]
z2 <- resample_z[(nrow(s)+1) : (nrow(s)*2)] + beta[2]
df <- data.frame(s1 = s[,1], s2 = s[,2], z1, z2)
layers <- c(AWU(r = r1, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = r1, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())
d <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1,  data = df, g = ~ 1,
                       layers = layers,
                       family = "matern_nonstat_symm",
                       method = "REML", nsteps = 150L
)
eta <- d$run(d$eta_tf)
LFT_pars <- d$run(d$layers[[12]]$pars)
scalings <- d$run(d$scalings)
nu_1 <- d$run(d$nu_tf_1)
nu_2 <- d$run(d$nu_tf_2)
sigma2_1 <- d$run(d$sigma2_tf_1)
sigma2_2 <- d$run(d$sigma2_tf_2)
sigma2_12 <- d$run(d$sigma2_tf_12)
l <- as.numeric(d$run(d$l_tf_1))
precy_1 <- d$run(d$precy_tf_1)
precy_2 <- d$run(d$precy_tf_2)
swarped <- d$run(d$swarped_tf1)
beta <- d$run(d$beta)
negcost <- d$negcost
parameter_est <- list(eta,LFT_pars,scalings,nu_1,nu_2,sigma2_1,sigma2_2,sigma2_12,l,precy_1,precy_2,swarped,beta,negcost)
save(parameter_est, file = paste0("results/simulation_symm_bootstrap_parameter_estimate/bootstrap_par_est_sample", i, ".rda"))
}

# Summarize bootstrap results
nu_1 <- list()
nu_2 <- list()
sigma_1 <- list()
sigma_2 <- list()
rho <- list()
tau_1 <- list()
tau_2 <- list()
beta_1 <- list()
beta_2 <- list()
l <- list()
for (i in 1:1000){
  load(paste0("results/simulation_symm_bootstrap_parameter_estimate/bootstrap_par_est_sample", i, ".rda"))
  nu_1[[i]] <- parameter_est[[4]]
  nu_2[[i]] <- parameter_est[[5]]
  sigma_1[[i]] <- sqrt(parameter_est[[6]])
  sigma_2[[i]] <- sqrt(parameter_est[[7]])
  rho[[i]] <- parameter_est[[8]] / sqrt(parameter_est[[6]]) / sqrt(parameter_est[[7]])
  tau_1[[i]] <- 1/sqrt(parameter_est[[10]])
  tau_2[[i]] <- 1/sqrt(parameter_est[[11]])
  beta_1[[i]] <- parameter_est[[13]][1]
  beta_2[[i]] <- parameter_est[[13]][2]
  sw <- parameter_est[[12]]
  l[[i]] <- parameter_est[[9]] / sqrt((sw[127,1] - sw[822,1])^2 + (sw[127,2] - sw[822,2])^2)
}
load("results/simulation_symm_parameter_estimate.rda")
# Boostrap confidence intervals
c(0.5, parameter_est[[4]], 2*parameter_est[[4]]- quantile(unlist(nu_1), c(0.975,0.025)))
c(1.5, parameter_est[[5]], 2*parameter_est[[5]]- quantile(unlist(nu_2), c(0.975,0.025)))
c(1, sqrt(parameter_est[[6]]), 2*sqrt(parameter_est[[6]])- quantile(unlist(sigma_1), c(0.975,0.025)))
c(0.9, sqrt(parameter_est[[7]]), 2*sqrt(parameter_est[[7]])- quantile(unlist(sigma_2), c(0.975,0.025)))
c(0.45, parameter_est[[8]] / sqrt(parameter_est[[6]]) / sqrt(parameter_est[[7]]), 2*parameter_est[[8]] / sqrt(parameter_est[[6]]) / sqrt(parameter_est[[7]])- quantile(unlist(rho), c(0.975,0.025)))
c(0.2, 1/sqrt(parameter_est[[10]]), 2*1/sqrt(parameter_est[[10]])- quantile(unlist(tau_1), c(0.975,0.025)))
c(0.1, 1/sqrt(parameter_est[[11]]), 2*1/sqrt(parameter_est[[11]])- quantile(unlist(tau_2), c(0.975,0.025)))
c(0, parameter_est[[13]][1], 2*parameter_est[[13]][1]- quantile(unlist(beta_1), c(0.975,0.025)))
c(0, parameter_est[[13]][2], 2*parameter_est[[13]][2]- quantile(unlist(beta_2), c(0.975,0.025)))
sw <- parameter_est[[12]]
l_star <- parameter_est[[9]] / sqrt((sw[127,1] - sw[822,1])^2 + (sw[127,2] - sw[822,2])^2)
c(0.3286223, l_star, 2*l_star - quantile(unlist(l), c(0.975,0.025)))
