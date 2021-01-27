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
  minsw1 <- min(c(swarped1[,1], swarped2[,1]))
  maxsw1 <- max(c(swarped1[,1], swarped2[,1]))
  minsw2 <- min(c(swarped1[,2], swarped2[,2]))
  maxsw2 <- max(c(swarped1[,2], swarped2[,2]))
  swarped1[,1] <- (swarped1[,1] - minsw1) / (maxsw1 - minsw1) - 0.5
  swarped1[,2] <- (swarped1[,2] - minsw2) / (maxsw2 - minsw2) - 0.5
  swarped2[,1] <- (swarped2[,1] - minsw1) / (maxsw1 - minsw1) - 0.5
  swarped2[,2] <- (swarped2[,2] - minsw2) / (maxsw2 - minsw2) - 0.5
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

# Cross validation experiment
all_data <- df

rmse1_model1 <- rmse2_model1 <- crps1_model1 <- crps2_model1 <- Cost1 <- time1 <-
  rmse1_model2 <- rmse2_model2 <- crps1_model2 <- crps2_model2 <- Cost2 <- time2 <-
  rmse1_model3 <- rmse2_model3 <- crps1_model3 <- crps2_model3 <- Cost3 <- time3 <-
  rmse1_model4 <- rmse2_model4 <- crps1_model4 <- crps2_model4 <- Cost4 <- time4 <-
  rmse1_model5 <- rmse2_model5 <- crps1_model5 <- crps2_model5 <- Cost5 <- time5 <-
  rep(0,30)

layers <- c(AWU(r = 50, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 50, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())

layers_asym <- c(AFF_2D())

for (i in 1:30){
# Sample a subset of data for the experiment
RNGkind(sample.kind = "Rounding")
set.seed(i)
df <- all_data
sam2 <- sample(1:nrow(df), 1000)
train_data <- df[sam2,]
test_data <- dplyr::setdiff(df, train_data)

df <- dplyr::select(train_data, s1, s2, z1, z2)
newdata <- dplyr::select(test_data, s1, s2)

  # Fit Model 1 (stationary, symmetric)
  t1 <- proc.time()
  d1 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          family = "matern_stat_symm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time1[i] <- (t2 - t1)[3]
  predML1 <- predict(d1, newdata = newdata)
  
  # Fit Model 2 (stationary, asymmetric)
  t1 <- proc.time()
  d2 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          layers_asym = layers_asym,
                          family = "matern_stat_asymm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time2[i] <- (t2 - t1)[3]
  predML2 <- predict(d2, newdata = newdata)
  
  # Fit Model 3 (univariate nonstationary)
  t1 <- proc.time()
  d3.1 <- deepspat_GP(f = z1 ~ s1 + s2 - 1, data = df, g = ~ 1,
                      layers = layers,
                      family = "matern_nonstat",
                      method = "REML", nsteps = 150L
  )
  
  d3.2 <- deepspat_GP(f = z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                      layers = layers,
                      family = "matern_nonstat",
                      method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time3[i] <- (t2 - t1)[3]
  predML3.1 <- predict(d3.1, newdata = newdata)
  predML3.2 <- predict(d3.2, newdata = newdata)
  
  # Fit Model 4 (nonstationary, symmetric)
  t1 <- proc.time()
  d4 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          layers = layers,
                          family = "matern_nonstat_symm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time4[i] <- (t2 - t1)[3]
  predML4 <- predict(d4, newdata = newdata)
  
  # Fit Model 5 (nonstationary, asymmetric)
  t1 <- proc.time()
  d5 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          layers = layers, layers_asym = layers_asym,
                          family = "matern_nonstat_asymm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time5[i] <- (t2 - t1)[3]
  predML5 <- predict(d5, newdata = newdata)
  
  df_pred1 <- data.frame(predML1$df_pred, y1=test_data$y1, y2=test_data$y2)
  df_pred2 <- data.frame(predML2$df_pred, y1=test_data$y1, y2=test_data$y2)
  df_pred3.1 <- data.frame(predML3.1$df_pred, y1=test_data$y1)
  df_pred3.2 <- data.frame(predML3.2$df_pred, y2=test_data$y2)
  df_pred4 <- data.frame(predML4$df_pred, y1=test_data$y1, y2=test_data$y2)
  df_pred5 <- data.frame(predML5$df_pred, y1=test_data$y1, y2=test_data$y2)
  
  # Save RMSPE and CRPS
  rmse1_model1[i] <- RMSPE(df_pred1$y1, df_pred1$pred_mean_1)
  rmse2_model1[i] <- RMSPE(df_pred1$y2, df_pred1$pred_mean_2)
  crps1_model1[i] <- CRPS(df_pred1$y1, df_pred1$pred_mean_1, df_pred1$pred_var_1)
  crps2_model1[i] <- CRPS(df_pred1$y2, df_pred1$pred_mean_2, df_pred1$pred_var_2)
  
  rmse1_model2[i] <- RMSPE(df_pred2$y1, df_pred2$pred_mean_1)
  rmse2_model2[i] <- RMSPE(df_pred2$y2, df_pred2$pred_mean_2)
  crps1_model2[i] <- CRPS(df_pred2$y1, df_pred2$pred_mean_1, df_pred2$pred_var_1)
  crps2_model2[i] <- CRPS(df_pred2$y2, df_pred2$pred_mean_2, df_pred2$pred_var_2)
  
  rmse1_model3[i] <- RMSPE(df_pred3.1$y1, df_pred3.1$pred_mean)
  rmse2_model3[i] <- RMSPE(df_pred3.2$y2, df_pred3.2$pred_mean)
  crps1_model3[i] <- CRPS(df_pred3.1$y1, df_pred3.1$pred_mean, df_pred3.1$pred_var)
  crps2_model3[i] <- CRPS(df_pred3.2$y2, df_pred3.2$pred_mean, df_pred3.2$pred_var)
  
  rmse1_model4[i] <- RMSPE(df_pred4$y1, df_pred4$pred_mean_1)
  rmse2_model4[i] <- RMSPE(df_pred4$y2, df_pred4$pred_mean_2)
  crps1_model4[i] <- CRPS(df_pred4$y1, df_pred4$pred_mean_1, df_pred4$pred_var_1)
  crps2_model4[i] <- CRPS(df_pred4$y2, df_pred4$pred_mean_2, df_pred4$pred_var_2)
  
  rmse1_model5[i] <- RMSPE(df_pred5$y1, df_pred5$pred_mean_1)
  rmse2_model5[i] <- RMSPE(df_pred5$y2, df_pred5$pred_mean_2)
  crps1_model5[i] <- CRPS(df_pred5$y1, df_pred5$pred_mean_1, df_pred5$pred_var_1)
  crps2_model5[i] <- CRPS(df_pred5$y2, df_pred5$pred_mean_2, df_pred5$pred_var_2)
  
  Cost1[i] <- d1$run(d1$Cost)
  Cost2[i] <- d2$run(d2$Cost)
  Cost3[i] <- d3.1$run(d3.1$Cost) + d3.2$run(d3.2$Cost)
  Cost4[i] <- d4$run(d4$Cost)
  Cost5[i] <- d5$run(d5$Cost)
  
  if (i == 1){
    
    # Save parameter estimates
    eta <- d5$run(d5$eta_tf)
    AFF_pars <- d5$run(d5$layers_asym[[1]]$pars)
    LFT_pars <- d5$run(d5$layers[[12]]$pars)
    scalings <- d5$run(d5$scalings)
    scalings_asym <- d5$run(d5$scalings_asym)
    nu_1 <- d5$run(d5$nu_tf_1)
    nu_2 <- d5$run(d5$nu_tf_2)
    sigma2_1 <- d5$run(d5$sigma2_tf_1)
    sigma2_2 <- d5$run(d5$sigma2_tf_2)
    sigma2_12 <- d5$run(d5$sigma2_tf_12)
    l <- as.numeric(d5$run(d5$l_tf_1))
    precy_1 <- d5$run(d5$precy_tf_1)
    precy_2 <- d5$run(d5$precy_tf_2)
    s_warped1 <- d5$run(d5$swarped_tf1)
    s_warped2 <- d5$run(d5$swarped_tf2)
    beta <- d5$run(d5$beta)
    parameter_est <- list(eta, AFF_pars, LFT_pars, scalings, scalings_asym,
                          nu_1, nu_2, sigma2_1, sigma2_2, sigma2_12, l, precy_1, precy_2,
                          s_warped1, s_warped2, beta)
    save(parameter_est, file = paste0("results/simulation_asymm_parameter_estimate.rda"))
    
    newdata <- all_data
    
    # Save predictions for plotting the comparison
    predML2 <- predict(d2, newdata = newdata)
    predML4 <- predict(d4, newdata = newdata)
    predML5 <- predict(d5, newdata = newdata)
    
    predML <- list(predML2$df_pred, predML4$df_pred, predML5$df_pred)
    save(predML, file = "results/simulation_asymm_predictions.rda")
  }
}

# Save cross validation results
crossvalid <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1, Cost1, time1,
                   rmse1_model2, crps1_model2, rmse2_model2, crps2_model2, Cost2, time2,
                   rmse1_model3, crps1_model3, rmse2_model3, crps2_model3, Cost3, time3,
                   rmse1_model4, crps1_model4, rmse2_model4, crps2_model4, Cost4, time4,
                   rmse1_model5, crps1_model5, rmse2_model5, crps2_model5, Cost5, time5
)

save(crossvalid, file="results/simulation_asymm_crossvalidation.rda")

# Cross validation results
matrix(unlist(lapply(crossvalid, mean)), nrow = 5, byrow=T)
