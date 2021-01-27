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

length.out <- 51
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped <- s

for(j in 1: (nlayers)){
  swarped <- layers[[j]]$fR(swarped, eta[[j]]) %>% scal_0_5_mat()
}

# Simulate symmetric, stationary covariance on the warped domain
D11 <- fields::rdist(swarped)
D22 <- fields::rdist(swarped)
D33 <- fields::rdist(swarped)

D12 <- fields::rdist(swarped, swarped)
D13 <- fields::rdist(swarped, swarped)
D23 <- fields::rdist(swarped, swarped)

C11 <- fields::Matern(D11, range=0.2, nu=0.5, phi=1)
C22 <- fields::Matern(D22, range=0.2, nu=1, phi=0.9^2)
C33 <- fields::Matern(D33, range=0.2, nu=1.5, phi=0.8^2)

C12 <- fields::Matern(D12, range=0.2, nu=0.75, phi=0.3)
C13 <- fields::Matern(D13, range=0.2, nu=1, phi=0.2)
C23 <- fields::Matern(D23, range=0.2, nu=1.25, phi=0.2)

C <- rbind(cbind(C11, C12, C13), cbind(t(C12), C22, C23), cbind(t(C13), t(C23), C33))
K <- t(chol(C))
i <- 11
set.seed(i)
y <- K %*% rnorm(nrow(s)*3)
y1 <- y[1: nrow(s),]
y2 <- y[(nrow(s)+1) : (nrow(s)*2),]
y3 <- y[(2*nrow(s)+1) : (nrow(s)*3),]

z1 <- y1 + 0.2*rnorm(length(y1))
z2 <- y2 + 0.1*rnorm(length(y2))
z3 <- y3 + 0.05*rnorm(length(y3))

z <- c(z1, z2, z3)
df <- data.frame(s, y1, y2, y3, z1, z2, z3)
names(df) <- c("s1", "s2", "y1", "y2", "y3", "z1", "z2", "z3")

# Save data
save(df, file="results/simulation_trivar_dataset.rda")

# Cross validation experiment
layers <- c(AWU(r = 50, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 50, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())

# Sample a subset of data for the experiment
RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
all_data <- df
train_data <- df[sam2,]
test_data <- dplyr::setdiff(df, train_data)

df <- dplyr::select(train_data, s1, s2, z1, z2, z3)
newdata <- dplyr::select(test_data, s1, s2)

# Fit Model 1 (stationary, symmetric)
t1 <- proc.time()
d1 <- deepspat_trivar_GP(f = z1 + z2 + z3 ~ s1 + s2 - 1, data = df, g = ~ 1,
                         family = "matern_stat_symm",
                         method = "REML", nsteps = 150L
)
t2 <- proc.time()
time1 <- (t2 - t1)[3]
predML1 <- predict(d1, newdata = newdata)

# Fit Model 2 (nonstationary, symmetric)
t1 <- proc.time()
d2 <- deepspat_trivar_GP(f = z1 + z2 + z3 ~ s1 + s2 - 1, data = df, g = ~ 1,
                         layers = layers,
                         family = "matern_nonstat_symm",
                         method = "REML", nsteps = 500L
)
t2 <- proc.time()
time2 <- (t2 - t1)[3]
predML2 <- predict(d2, newdata = newdata)

df_pred1 <- data.frame(predML1$df_pred, y1=test_data$y1, y2=test_data$y2, y3=test_data$y3)
df_pred2 <- data.frame(predML2$df_pred, y1=test_data$y1, y2=test_data$y2, y3=test_data$y3)

# Save RMSPE and CRPS
rmse1_model1 <- RMSPE(df_pred1$y1, df_pred1$pred_mean_1)
rmse2_model1 <- RMSPE(df_pred1$y2, df_pred1$pred_mean_2)
rmse3_model1 <- RMSPE(df_pred1$y3, df_pred1$pred_mean_3)
crps1_model1 <- CRPS(df_pred1$y1, df_pred1$pred_mean_1, df_pred1$pred_var_1)
crps2_model1 <- CRPS(df_pred1$y2, df_pred1$pred_mean_2, df_pred1$pred_var_2)
crps3_model1 <- CRPS(df_pred1$y3, df_pred1$pred_mean_3, df_pred1$pred_var_3)

rmse1_model2 <- RMSPE(df_pred2$y1, df_pred2$pred_mean_1)
rmse2_model2 <- RMSPE(df_pred2$y2, df_pred2$pred_mean_2)
rmse3_model2 <- RMSPE(df_pred2$y3, df_pred2$pred_mean_3)
crps1_model2 <- CRPS(df_pred2$y1, df_pred2$pred_mean_1, df_pred2$pred_var_1)
crps2_model2 <- CRPS(df_pred2$y2, df_pred2$pred_mean_2, df_pred2$pred_var_2)
crps3_model2 <- CRPS(df_pred2$y3, df_pred2$pred_mean_3, df_pred2$pred_var_3)

Cost1 <- d1$run(d1$Cost)
Cost2 <- d2$run(d2$Cost)

# Save parameter estimates
eta <- d2$run(d2$eta_tf)
LFT_pars <- d2$run(d2$layers[[12]]$pars)
scalings <- d2$run(d2$scalings)
nu_1 <- d2$run(d2$nu_tf_1)
nu_2 <- d2$run(d2$nu_tf_2)
nu_3 <- d2$run(d2$nu_tf_3)
sigma2_1 <- d2$run(d2$sigma2_tf_1)
sigma2_2 <- d2$run(d2$sigma2_tf_2)
sigma2_3 <- d2$run(d2$sigma2_tf_3)
sigma2_12 <- d2$run(d2$sigma2_tf_12)
sigma2_13 <- d2$run(d2$sigma2_tf_13)
sigma2_23 <- d2$run(d2$sigma2_tf_23)
l <- as.numeric(d2$run(d2$l_tf_1))
precy_1 <- d2$run(d2$precy_tf_1)
precy_2 <- d2$run(d2$precy_tf_2)
precy_3 <- d2$run(d2$precy_tf_3)
s_warped <- d2$run(d2$swarped_tf1)
beta <- d2$run(d2$beta)
parameter_est <- list(eta, LFT_pars, scalings,
                      nu_1, nu_2, nu_3, sigma2_1, sigma2_2, sigma2_3, sigma2_12, sigma2_13, sigma2_23,
                      l, precy_1, precy_2, precy_3,
                      s_warped, beta)
save(parameter_est, file = paste0("results/simulation_trivar_parameter_estimate.rda"))

# Save cross validation results
crossvalid <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1, rmse3_model1, crps3_model1, Cost1, time1,
                   rmse1_model2, crps1_model2, rmse2_model2, crps2_model2, rmse3_model2, crps3_model2, Cost2, time2
)

save(crossvalid, file="results/simulation_trivar_crossvalidation.rda")

# Cross validation results
matrix(unlist(lapply(crossvalid, mean)), nrow = 2, byrow=T)
