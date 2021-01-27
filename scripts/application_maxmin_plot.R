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

# Load data for plotting
load("results/application_maxmin_data_for_plot.rda")
sim <- data_for_plot[[1]]
predML1 <- data_for_plot[[2]]
predML2 <- data_for_plot[[3]]
noise_var_1_model1 <- data_for_plot[[4]]
noise_var_2_model1 <- data_for_plot[[5]]
noise_var_1_model2 <- data_for_plot[[6]]
noise_var_2_model2 <- data_for_plot[[7]]
test_data <- sim[(sim$LATITUDE > 36) & (sim$LATITUDE < 39) & (sim$LONGITUDE > -108) & (sim$LONGITUDE < -104),]
train_data <- setdiff(sim, test_data)

# Plotting
b1min <- min(c(predML1$df_pred$pred_mean_1, predML2$df_pred$pred_mean_1, sim$TMAX))
b1max <- max(c(predML1$df_pred$pred_mean_1, predML2$df_pred$pred_mean_1, sim$TMAX))
b2min <- min(c(predML1$df_pred$pred_mean_2, predML2$df_pred$pred_mean_2, sim$TMIN))
b2max <- max(c(predML1$df_pred$pred_mean_2, predML2$df_pred$pred_mean_2, sim$TMIN))
c1min <- min(c(sqrt(predML1$df_pred$pred_var_1 + noise_var_1_model1), sqrt(predML2$df_pred$pred_var_1 + noise_var_1_model2)))
c1max <- max(c(sqrt(predML1$df_pred$pred_var_1 + noise_var_1_model1), sqrt(predML2$df_pred$pred_var_1 + noise_var_1_model2)))
c2min <- min(c(sqrt(predML1$df_pred$pred_var_2 + noise_var_2_model1), sqrt(predML2$df_pred$pred_var_2 + noise_var_2_model2)))
c2max <- max(c(sqrt(predML1$df_pred$pred_var_2 + noise_var_2_model1), sqrt(predML2$df_pred$pred_var_2 + noise_var_2_model2)))

dt <- data.frame(xmin=-108, xmax=-104, ymin=36, ymax=39, b=0) 
a <-  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=b), inherit.aes=FALSE, color="black", alpha=0.5)

g2a1 <- ggplot(sim) +
  geom_point(aes(LONGITUDE, LATITUDE, color = TMAX)) +
  scale_color_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "Z1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0) +
  labs(x = "s1", y = "s2")

g2a2 <- ggplot(sim) +
  geom_point(aes(LONGITUDE, LATITUDE, color = TMIN)) +
  scale_color_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "Z2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0) +
  labs(x = "s1", y = "s2")

g2b1.1 <- ggplot(predML1$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2b2.1 <- ggplot(predML1$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b1.1 <- ggplot(predML1$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1 + noise_var_1_model1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b2.1 <- ggplot(predML1$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2 + noise_var_2_model1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2b1.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2b2.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b1.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1 + noise_var_1_model2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b2.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2 + noise_var_2_model2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2a <- ggplot(train_data) +
  geom_point(aes(LONGITUDE, LATITUDE, color = "black")) +
  scale_color_identity(name = "", breaks = c("black"), labels = c("locs"), guide = "legend") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  labs(x = "s1", y = "s2") +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

gplot <- grid.arrange(g2a1, g2b1.1, g2b1.2,
                      g2a, g3b1.1, g3b1.2,
                      g2a2, g2b2.1, g2b2.2,
                      g2a, g3b2.1, g3b2.2,
                      nrow = 4
)

ggsave(filename="figures/application_maxmin_comparison_block_holdout_plot.png", plot=gplot, device="png", width=33, height=44, scale=1, units="cm", dpi=300)

# Plot the true and estimated warping function
# Pick 3 points
k <- 714
l <- 190
m <- 508

# True warping
sim <- read.csv("data/1796013.csv")
sim <- sim[!is.na(sim$TMAX),]
sim <- sim[!is.na(sim$TMIN),]
all_data <- sim
test_data <- sim[(sim$LATITUDE > 36) & (sim$LATITUDE < 39) & (sim$LONGITUDE > -108) & (sim$LONGITUDE < -104),]
train_data <- setdiff(all_data, test_data)

sw <- matrix(c(train_data$LONGITUDE, train_data$LATITUDE), ncol=2)
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])

warp.plot1 <- ggplot() + geom_point(data=sw.df, aes(x, y), colour = "black", size = 0.5) +
  theme_bw() + coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.5,1)) +
  theme(text = element_text(size=15)) +
  labs(x=expression(paste('b'['0,1'], '(', bold(s), ')')), y=expression(paste('b'['0,2'], '(', bold(s), ')')))

# Estimated warping
load("results/application_maxmin_parameter_estimate.rda")

sw <- parameter_est[[12]]
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
warp.plot2 <- ggplot() + geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=0.5, alpha=1) +
  theme_bw() + coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.5,1)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

gplot <- grid.arrange(warp.plot1, warp.plot2, nrow = 1)

ggsave(filename="figures/application_maxmin_warping_plot.png", plot=gplot, device="png", width=20, height=10, scale=1, units="cm", dpi=300)