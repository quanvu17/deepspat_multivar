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
load("results/simulation_symm_dataset.rda")

# Gap experiment
rmse1_model1 <- rmse2_model1 <- crps1_model1 <- crps2_model1 <-
  rmse1_model2 <- rmse2_model2 <- crps1_model2 <- crps2_model2 <-
  Cost1 <- time1 <- Cost2 <- time2 <- rep(0,30)

layers <- c(AWU(r = 50L, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 50L, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())

for (i in 1:30){
  df1 <- df[df$s1 > -0.28 & df$s1 < -0.08 & df$s2 > -0.48 & df$s2 < -0.28,]
  df0 <- setdiff(df,df1)
  RNGkind(sample.kind = "Rounding")
  set.seed(i)
  sam <- sample(1:nrow(df0), 1000)
  df3 <- df0[sam,]
  
  # Fit Models
  newdata <- df
  
  # Fit Model 1 (stationary, symmetric)
  t1 <- proc.time()
  d1 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df3, g = ~ 1,
                          family = "matern_stat_symm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time1[i] <- (t2 - t1)[3]
  predML3 <- predict(d1, newdata = newdata)
  
  # Fit Model 2 (nonstationary, symmetric)
  t1 <- proc.time()
  d2 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df3, g = ~ 1,
                          layers = layers,
                          family = "matern_nonstat_symm",
                          method = "REML", nsteps = 150L
  )
  t2 <- proc.time()
  time2[i] <- (t2 - t1)[3]
  predML4 <- predict(d2, newdata = newdata)
  
  predML3.0 <- predML3$df_pred[predML3$df_pred$s1 > -0.28 & predML3$df_pred$s1 < -0.08 & predML3$df_pred$s2 > -0.48 & predML3$df_pred$s2 < -0.28,]
  predML4.0 <- predML4$df_pred[predML4$df_pred$s1 > -0.28 & predML4$df_pred$s1 < -0.08 & predML4$df_pred$s2 > -0.48 & predML4$df_pred$s2 < -0.28,]
  
  # Save RMSPE and CRPS for the gap experiment
  rmse1_model1[i] <- RMSPE(df1$y1, predML3.0$pred_mean_1)
  rmse2_model1[i] <- RMSPE(df1$y2, predML3.0$pred_mean_2)
  crps1_model1[i] <- CRPS(df1$y1, predML3.0$pred_mean_1, predML3.0$pred_var_1)
  crps2_model1[i] <- CRPS(df1$y2, predML3.0$pred_mean_2, predML3.0$pred_var_2)
  
  rmse1_model2[i] <- RMSPE(df1$y1, predML4.0$pred_mean_1)
  rmse2_model2[i] <- RMSPE(df1$y2, predML4.0$pred_mean_2)
  crps1_model2[i] <- CRPS(df1$y1, predML4.0$pred_mean_1, predML4.0$pred_var_1)
  crps2_model2[i] <- CRPS(df1$y2, predML4.0$pred_mean_2, predML4.0$pred_var_2)
  
  Cost1[i] <- d1$run(d1$Cost)
  Cost2[i] <- d2$run(d2$Cost)
  
  if (i == 2){
    # Save predictions for plotting the comparison
    predML <- list(predML3$df_pred, predML4$df_pred)
    save(predML, file = "results/simulation_symm_predictions_gap.rda")
  }
  
}

block_holdout <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1, Cost1, time1,
                      rmse1_model2, crps1_model2, rmse2_model2, crps2_model2, Cost2, time2
)

save(block_holdout, file="results/simulation_symm_block_holdout.rda")

# Gap results
matrix(unlist(lapply(block_holdout, mean)), nrow = 2, byrow=T)
