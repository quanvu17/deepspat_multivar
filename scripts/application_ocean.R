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
ocean_temp <- nc_open("data/global-analysis-forecast-phy-001-024_1575868462143.nc")
temp_1m <- ncvar_get(ocean_temp, "thetao")[, , 1]
temp_318m <- ncvar_get(ocean_temp, "thetao")[, , 29]
lat <- ncvar_get(ocean_temp, "latitude")
lon <- ncvar_get(ocean_temp, "longitude")
sim <- data.frame(temp_1m = as.vector(temp_1m), temp_318m = as.vector(temp_318m),
                  lon = rep(lon,length(lat)),
                  lat = as.vector(t(matrix(rep(lat,length(lon)), ncol=length(lon))))
)
sim <- sim[!is.na(sim$temp_1m),]
sim <- sim[!is.na(sim$temp_318m),]

# Five fold Cross validation experiment
RNGkind(sample.kind = "Rounding")
set.seed(1)
sam <- sample(1:nrow(sim), floor(nrow(sim)/5)*5)
all_data <- sim[sam, ]
groups_data <- split(all_data, (seq(nrow(all_data))-1) %/% (nrow(all_data)/5))

rmse1_model1 <- rmse2_model1 <- crps1_model1 <- crps2_model1 <- rmse1_model2 <- rmse2_model2 <- crps1_model2 <- crps2_model2 <-
  rmse1_model3 <- rmse2_model3 <- crps1_model3 <- crps2_model3 <- rep(0,5)

layers <- c(AWU(r = 10, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 10, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())

layers_asym <- c(AFF_2D())

for (i in 1 : 5) {
  test_data <- groups_data[[i]]
  train_data <- setdiff(all_data, test_data) 
  df <- cbind(train_data$lon, train_data$lat, train_data$temp_1m, train_data$temp_318m) %>% as.data.frame()
  names(df) <- c("s1", "s2", "z1", "z2")
  newdata <- data.frame(s1 = test_data$lon,
                        s2 = test_data$lat,
                        y1 = test_data$temp_1m,
                        y2 = test_data$temp_318m)

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
    
    # Fit Model 3 (nonstationary, asymmetric)
    d3 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                            layers = layers, layers_asym = layers_asym,
                            family = "matern_nonstat_asymm",
                            method = "REML", nsteps = 150L
    )
    predML3 <- predict.deepspat_multivar(d3, newdata = newdata)
    
    df_pred1 <- data.frame(predML1$df_pred, z1=test_data$temp_1m, z2=test_data$temp_318m)
    df_pred2 <- data.frame(predML2$df_pred, z1=test_data$temp_1m, z2=test_data$temp_318m)
    df_pred3 <- data.frame(predML3$df_pred, z1=test_data$temp_1m, z2=test_data$temp_318m)

    pred_var_1_model1 <- df_pred1$pred_var_1 + 1/d1$run(d1$precy_tf_1)
    pred_var_2_model1 <- df_pred1$pred_var_2 + 1/d1$run(d1$precy_tf_2)
    pred_var_1_model2 <- df_pred2$pred_var_1 + 1/d2$run(d2$precy_tf_1)
    pred_var_2_model2 <- df_pred2$pred_var_2 + 1/d2$run(d2$precy_tf_2)
    pred_var_1_model3 <- df_pred3$pred_var_1 + 1/d3$run(d3$precy_tf_1)
    pred_var_2_model3 <- df_pred3$pred_var_2 + 1/d3$run(d3$precy_tf_2)
    
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
    
}
  
# Save cross validation results
crossvalid <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1,
                   rmse1_model2, crps1_model2, rmse2_model2, crps2_model2,
                   rmse1_model3, crps1_model3, rmse2_model3, crps2_model3
)
  
save(crossvalid, file="results/application_ocean_crossvalidation.rda")
  
# Cross validation results
matrix(unlist(lapply(crossvalid, mean)), nrow = 3, byrow=T)

# Gap experiment
all_data <- sim
test_data <- sim[sim$lat > 37.5 & sim$lat < 38.2,]
train_data <- setdiff(sim, test_data) 

df <- dplyr::select(train_data, lon, lat, temp_1m, temp_318m) %>% 
      rename(s1 = lon, s2 = lat, z1 = temp_1m, z2 = temp_318m)
newdata <- dplyr::select(sim, lon, lat) %>% rename(s1 = lon, s2 = lat)

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

# Fit Model 3 (nonstationary, asymmetric)
d3 <- deepspat_multivar(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                        layers = layers, layers_asym = layers_asym,
                        family = "matern_nonstat_asymm",
                        method = "REML", nsteps = 150L
)
predML3 <- predict.deepspat_multivar(d3, newdata = newdata)

predML1.0 <- predML1$df_pred[predML1$df_pred$s2 > 37.5 & predML1$df_pred$s2 < 38.2,]
predML2.0 <- predML2$df_pred[predML2$df_pred$s2 > 37.5 & predML2$df_pred$s2 < 38.2,]
predML3.0 <- predML3$df_pred[predML3$df_pred$s2 > 37.5 & predML3$df_pred$s2 < 38.2,]

pred_var_1_model1 <- predML1.0$pred_var_1 + 1/d1$run(d1$precy_tf_1)
pred_var_2_model1 <- predML1.0$pred_var_2 + 1/d1$run(d1$precy_tf_2)
pred_var_1_model2 <- predML2.0$pred_var_1 + 1/d2$run(d2$precy_tf_1)
pred_var_2_model2 <- predML2.0$pred_var_2 + 1/d2$run(d2$precy_tf_2)
pred_var_1_model3 <- predML3.0$pred_var_1 + 1/d3$run(d3$precy_tf_1)
pred_var_2_model3 <- predML3.0$pred_var_2 + 1/d3$run(d3$precy_tf_2)

# Save RMSPE and CRPS for the gap experiment
rmse1_model1 <- RMSPE(test_data$temp_1m, predML1.0$pred_mean_1)
rmse2_model1 <- RMSPE(test_data$temp_318m, predML1.0$pred_mean_2)
crps1_model1 <- CRPS(test_data$temp_1m, predML1.0$pred_mean_1, pred_var_1_model1)
crps2_model1 <- CRPS(test_data$temp_318m, predML1.0$pred_mean_2, pred_var_2_model1)

rmse1_model2 <- RMSPE(test_data$temp_1m, predML2.0$pred_mean_1)
rmse2_model2 <- RMSPE(test_data$temp_318m, predML2.0$pred_mean_2)
crps1_model2 <- CRPS(test_data$temp_1m, predML2.0$pred_mean_1, pred_var_1_model2)
crps2_model2 <- CRPS(test_data$temp_318m, predML2.0$pred_mean_2, pred_var_2_model2)

rmse1_model3 <- RMSPE(test_data$temp_1m, predML3.0$pred_mean_1)
rmse2_model3 <- RMSPE(test_data$temp_318m, predML3.0$pred_mean_2)
crps1_model3 <- CRPS(test_data$temp_1m, predML3.0$pred_mean_1, pred_var_1_model3)
crps2_model3 <- CRPS(test_data$temp_318m, predML3.0$pred_mean_2, pred_var_2_model3)

block_holdout <- list(rmse1_model1, crps1_model1, rmse2_model1, crps2_model1,
                      rmse1_model2, crps1_model2, rmse2_model2, crps2_model2,
                      rmse1_model3, crps1_model3, rmse2_model3, crps2_model3)

save(block_holdout, file="results/application_ocean_block_holdout.rda")

# Save data for plotting
data_for_plot <- list(sim, predML1, predML2, predML3)
save(data_for_plot, file = "results/application_ocean_data_for_plot.rda")
