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

# Load data
sim <- read.csv("data/1796013.csv")
sim <- sim[!is.na(sim$TMAX),]
sim <- sim[!is.na(sim$TMIN),]

# Five fold Cross validation experiment
RNGkind(sample.kind = "Rounding")
set.seed(1)
sam <- sample(1:nrow(sim), floor(nrow(sim)/5)*5)
all_data <- sim[sam,]
groups_data <- split(all_data, (seq(nrow(all_data))-1) %/% (nrow(all_data)/5))

rmse1_model1 <- rmse2_model1 <- crps1_model1 <- crps2_model1 <-
  rmse1_model2 <- rmse2_model2 <- crps1_model2 <- crps2_model2 <-
  rmse1_model3 <- rmse2_model3 <- crps1_model3 <- crps2_model3 <-
  rmse1_model4 <- rmse2_model4 <- crps1_model4 <- crps2_model4 <-
  Cost1 <- Cost2 <- Cost3 <- Cost4 <- time1 <- time2 <- time3 <- time4 <-
  rep(0,5)

layers <- c(AWU(r = 10L, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 10L, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(1L),
            LFT())

for (i in 1:5){
  test_data <- groups_data[[i]]
  train_data <- setdiff(all_data, test_data)
  
  df <- dplyr::select(train_data, LONGITUDE, LATITUDE, ELEVATION, TMAX, TMIN) %>%
    rename(s1 = LONGITUDE, s2 = LATITUDE, s3 = ELEVATION, z1 = TMAX, z2 = TMIN)
  newdata <- dplyr::select(test_data, LONGITUDE, LATITUDE, ELEVATION) %>%
    rename(s1 = LONGITUDE, s2 = LATITUDE, s3 = ELEVATION)
  
  # Fit Model 1 (stationary, symmetric)
  t1 <- proc.time()
  d1 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          family = "matern_stat_symm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time1[i] <- (t2 - t1)[3]
  predML1 <- predict(d1, newdata = newdata)
  
  # Fit Model 2 (nonstationary, symmetric)
  t1 <- proc.time()
  d2 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                          layers = layers,
                          family = "matern_nonstat_symm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time2[i] <- (t2 - t1)[3]
  predML2 <- predict(d2, newdata = newdata)
  
  # Fit Model 3 (stationary, symmetric, with covariate in trend)
  t1 <- proc.time()
  d3 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ s3 + 1,
                          family = "matern_stat_symm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time3[i] <- (t2 - t1)[3]
  predML3 <- predict(d3, newdata = newdata)
  
  # Fit Model 4 (nonstationary, symmetric, with covariate in trend)
  t1 <- proc.time()
  d4 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ s3 + 1,
                          layers = layers,
                          family = "matern_nonstat_symm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time4[i] <- (t2 - t1)[3]
  predML4 <- predict(d4, newdata = newdata)
  
  df_pred1 <- data.frame(predML1$df_pred, z1=test_data$TMAX, z2=test_data$TMIN)
  df_pred2 <- data.frame(predML2$df_pred, z1=test_data$TMAX, z2=test_data$TMIN)
  df_pred3 <- data.frame(predML3$df_pred, z1=test_data$TMAX, z2=test_data$TMIN)
  df_pred4 <- data.frame(predML4$df_pred, z1=test_data$TMAX, z2=test_data$TMIN)
  
  pred_var_1_model1 <- df_pred1$pred_var_1 + 1/d1$run(d1$precy_tf_1)
  pred_var_2_model1 <- df_pred1$pred_var_2 + 1/d1$run(d1$precy_tf_2)
  pred_var_1_model2 <- df_pred2$pred_var_1 + 1/d2$run(d2$precy_tf_1)
  pred_var_2_model2 <- df_pred2$pred_var_2 + 1/d2$run(d2$precy_tf_2)
  pred_var_1_model3 <- df_pred3$pred_var_1 + 1/d3$run(d3$precy_tf_1)
  pred_var_2_model3 <- df_pred3$pred_var_2 + 1/d3$run(d3$precy_tf_2)
  pred_var_1_model4 <- df_pred4$pred_var_1 + 1/d4$run(d4$precy_tf_1)
  pred_var_2_model4 <- df_pred4$pred_var_2 + 1/d4$run(d4$precy_tf_2)
  
  # Save RMSPE and CRPS
  rmse1_model1[i] <- RMSPE(df_pred1$z1, df_pred1$pred_mean_1)
  rmse2_model1[i] <- RMSPE(df_pred1$z2, df_pred1$pred_mean_2)
  crps1_model1[i] <- CRPS(df_pred1$z1, df_pred1$pred_mean_1, pred_var_1_model1)
  crps2_model1[i] <- CRPS(df_pred1$z2, df_pred1$pred_mean_2, pred_var_2_model1)

  rmse1_model2[i] <- RMSPE(df_pred2$z1, df_pred2$pred_mean_1)
  rmse2_model2[i] <- RMSPE(df_pred2$z2, df_pred2$pred_mean_2)
  crps1_model2[i] <- CRPS(df_pred2$z1, df_pred2$pred_mean_1, pred_var_1_model2)
  crps2_model2[i] <- CRPS(df_pred2$z2, df_pred2$pred_mean_2, pred_var_2_model2)

  rmse1_model3[i] <- RMSPE(df_pred3$z1, df_pred3$pred_mean_1)
  rmse2_model3[i] <- RMSPE(df_pred3$z2, df_pred3$pred_mean_2)
  crps1_model3[i] <- CRPS(df_pred3$z1, df_pred3$pred_mean_1, pred_var_1_model3)
  crps2_model3[i] <- CRPS(df_pred3$z2, df_pred3$pred_mean_2, pred_var_2_model3)

  rmse1_model4[i] <- RMSPE(df_pred4$z1, df_pred4$pred_mean_1)
  rmse2_model4[i] <- RMSPE(df_pred4$z2, df_pred4$pred_mean_2)
  crps1_model4[i] <- CRPS(df_pred4$z1, df_pred4$pred_mean_1, pred_var_1_model4)
  crps2_model4[i] <- CRPS(df_pred4$z2, df_pred4$pred_mean_2, pred_var_2_model4)
  
  Cost1[i] <- d1$run(d1$Cost)
  Cost2[i] <- d2$run(d2$Cost)
  Cost3[i] <- d3$run(d3$Cost)
  Cost4[i] <- d4$run(d4$Cost)
  
}

# Save cross validation results
crossvalid <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1, Cost1, time1,
                   rmse1_model2, crps1_model2, rmse2_model2, crps2_model2, Cost2, time2,
                   rmse1_model3, crps1_model3, rmse2_model3, crps2_model3, Cost3, time3,
                   rmse1_model4, crps1_model4, rmse2_model4, crps2_model4, Cost4, time4
)

save(crossvalid, file="results/application_maxmin_crossvalidation.rda")

# Cross validation results
matrix(unlist(lapply(crossvalid, mean)), nrow = 4, byrow=T)

# Gap experiment
all_data <- sim
test_data <- sim[(sim$LATITUDE > 36) & (sim$LATITUDE < 39) & (sim$LONGITUDE > -108) & (sim$LONGITUDE < -104),]
train_data <- setdiff(all_data, test_data)
df <- dplyr::select(train_data, LONGITUDE, LATITUDE, ELEVATION, TMAX, TMIN) %>%
  rename(s1 = LONGITUDE, s2 = LATITUDE, s3 = ELEVATION, z1 = TMAX, z2 = TMIN)
newdata <- dplyr::select(test_data, LONGITUDE, LATITUDE, ELEVATION) %>%
  rename(s1 = LONGITUDE, s2 = LATITUDE, s3 = ELEVATION)


layers <- c(AWU(r = 10L, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 10L, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(1L),
            LFT())

# Fit Model 1 (stationary, symmetric)
t1 <- proc.time()
d1 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                        family = "matern_stat_symm",
                        method = "REML", nsteps = 150L
)
t2 <- proc.time()
time1 <- (t2 - t1)[3]
predML1 <- predict(d1, newdata = newdata)

# Fit Model 2 (nonstationary, symmetric)
t1 <- proc.time()
d2 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                        layers = layers,
                        family = "matern_nonstat_symm",
                        method = "REML", nsteps = 150L
)
t2 <- proc.time()
time2 <- (t2 - t1)[3]
predML2 <- predict(d2, newdata = newdata)

# Fit Model 3 (stationary, symmetric, with covariate in trend)
t1 <- proc.time()
d3 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ s3 + 1,
                        family = "matern_stat_symm",
                        method = "REML", nsteps = 150L
)
t2 <- proc.time()
time3 <- (t2 - t1)[3]
predML3 <- predict(d3, newdata = newdata)

# Fit Model 4 (nonstationary, symmetric, with covariate in trend)
t1 <- proc.time()
d4 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ s3 + 1,
                        layers = layers,
                        family = "matern_nonstat_symm",
                        method = "REML", nsteps = 150L
)
t2 <- proc.time()
time4 <- (t2 - t1)[3]
predML4 <- predict(d4, newdata = newdata)

Cost1 <- d1$run(d1$Cost)
Cost2 <- d2$run(d2$Cost)
Cost3 <- d3$run(d3$Cost)
Cost4 <- d4$run(d4$Cost)

predML1.0 <- predML1$df_pred[predML1$df_pred$s1 > -108 & predML1$df_pred$s1 < -104 & predML1$df_pred$s2 > 36 & predML1$df_pred$s2 < 39,]
predML2.0 <- predML2$df_pred[predML2$df_pred$s1 > -108 & predML2$df_pred$s1 < -104 & predML2$df_pred$s2 > 36 & predML2$df_pred$s2 < 39,]
predML3.0 <- predML3$df_pred[predML3$df_pred$s1 > -108 & predML3$df_pred$s1 < -104 & predML3$df_pred$s2 > 36 & predML3$df_pred$s2 < 39,]
predML4.0 <- predML4$df_pred[predML4$df_pred$s1 > -108 & predML4$df_pred$s1 < -104 & predML4$df_pred$s2 > 36 & predML4$df_pred$s2 < 39,]

pred_var_1_model1 <- predML1.0$pred_var_1 + 1/d1$run(d1$precy_tf_1)
pred_var_2_model1 <- predML1.0$pred_var_2 + 1/d1$run(d1$precy_tf_2)
pred_var_1_model2 <- predML2.0$pred_var_1 + 1/d2$run(d2$precy_tf_1)
pred_var_2_model2 <- predML2.0$pred_var_2 + 1/d2$run(d2$precy_tf_2)
pred_var_1_model3 <- predML3.0$pred_var_1 + 1/d3$run(d3$precy_tf_1)
pred_var_2_model3 <- predML3.0$pred_var_2 + 1/d3$run(d3$precy_tf_2)
pred_var_1_model4 <- predML4.0$pred_var_1 + 1/d4$run(d4$precy_tf_1)
pred_var_2_model4 <- predML4.0$pred_var_2 + 1/d4$run(d4$precy_tf_2)

# Save RMSPE and CRPS for the gap experiment
rmse1_model1 <- RMSPE(test_data$TMAX, predML1.0$pred_mean_1)
rmse2_model1 <- RMSPE(test_data$TMIN, predML1.0$pred_mean_2)
crps1_model1 <- CRPS(test_data$TMAX, predML1.0$pred_mean_1, pred_var_1_model1)
crps2_model1 <- CRPS(test_data$TMIN, predML1.0$pred_mean_2, pred_var_2_model1)

rmse1_model2 <- RMSPE(test_data$TMAX, predML2.0$pred_mean_1)
rmse2_model2 <- RMSPE(test_data$TMIN, predML2.0$pred_mean_2)
crps1_model2 <- CRPS(test_data$TMAX, predML2.0$pred_mean_1, pred_var_1_model2)
crps2_model2 <- CRPS(test_data$TMIN, predML2.0$pred_mean_2, pred_var_2_model2)

rmse1_model3 <- RMSPE(test_data$TMAX, predML3.0$pred_mean_1)
rmse2_model3 <- RMSPE(test_data$TMIN, predML3.0$pred_mean_2)
crps1_model3 <- CRPS(test_data$TMAX, predML3.0$pred_mean_1, pred_var_1_model3)
crps2_model3 <- CRPS(test_data$TMIN, predML3.0$pred_mean_2, pred_var_2_model3)

rmse1_model4 <- RMSPE(test_data$TMAX, predML4.0$pred_mean_1)
rmse2_model4 <- RMSPE(test_data$TMIN, predML4.0$pred_mean_2)
crps1_model4 <- CRPS(test_data$TMAX, predML4.0$pred_mean_1, pred_var_1_model4)
crps2_model4 <- CRPS(test_data$TMIN, predML4.0$pred_mean_2, pred_var_2_model4)

block_holdout <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1, Cost1, time1,
                      rmse1_model2, crps1_model2, rmse2_model2, crps2_model2, Cost2, time2,
                      rmse1_model3, crps1_model3, rmse2_model3, crps2_model3, Cost3, time3,
                      rmse1_model4, crps1_model4, rmse2_model4, crps2_model4, Cost4, time4
)

save(block_holdout, file="results/application_maxmin_block_holdout.rda")

# Block holdout results
matrix(unlist(lapply(block_holdout, mean)), nrow = 4, byrow=T)

# Predictions for the whole region (plotting of the gap experiment)
s <- data.frame(s1 = sim$LONGITUDE,
                s2 = sim$LATITUDE)
chx <- chull(s)
chx <- rbind(s = s[chx, ], s[chx[1], ])

s <- expand.grid(seq(min(sim$LONGITUDE), max(sim$LONGITUDE), length.out = 100),
                 seq(min(sim$LATITUDE), max(sim$LATITUDE), length.out = 100)) %>% as.matrix()

inside <- point.in.polygon(s[,1], s[,2], chx[,1], chx[,2], mode.checked=FALSE)
s <- s[inside != 0,]
newdata <- data.frame(s1=s[,1], s2=s[,2])

predML1 <- predict(d1, newdata = newdata)
predML2 <- predict(d2, newdata = newdata)

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
save(parameter_est, file = paste0("results/application_maxmin_parameter_estimate.rda"))

noise_var_1_model1 <- 1/d1$run(d1$precy_tf_1)
noise_var_2_model1 <- 1/d1$run(d1$precy_tf_2)
noise_var_1_model2 <- 1/d2$run(d2$precy_tf_1)
noise_var_2_model2 <- 1/d2$run(d2$precy_tf_2)

data_for_plot <- list(sim, predML1, predML2,
                      noise_var_1_model1, noise_var_2_model1,
                      noise_var_1_model2, noise_var_2_model2)
save(data_for_plot, file = "results/application_maxmin_data_for_plot.rda")
