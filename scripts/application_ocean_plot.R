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
load("results/application_ocean_data_for_plot.rda")
sim <- data_for_plot[[1]]
predML1 <- data_for_plot[[2]]
predML2 <- data_for_plot[[3]]
predML3 <- data_for_plot[[4]]

# Plot for gap experiment
b1min <- min(c(predML1$df_pred$pred_mean_1, predML2$df_pred$pred_mean_1, predML3$df_pred$pred_mean_1, sim$temp_1m))
b1max <- max(c(predML1$df_pred$pred_mean_1, predML2$df_pred$pred_mean_1, predML3$df_pred$pred_mean_1, sim$temp_1m))
b2min <- min(c(predML1$df_pred$pred_mean_2, predML2$df_pred$pred_mean_2, predML3$df_pred$pred_mean_2, sim$temp_318m))
b2max <- max(c(predML1$df_pred$pred_mean_2, predML2$df_pred$pred_mean_2, predML3$df_pred$pred_mean_2, sim$temp_318m))

dt <- data.frame(xmin=-63.3, xmax=-60, ymin=37.5, ymax=38.2, b=0) 
a <-  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=b), inherit.aes=FALSE, color="black", alpha=0.5)

g2a1 <- ggplot(sim) +
  geom_tile(aes(lon, lat, fill = temp_1m)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "Z1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0) +
  labs(x = "s1", y = "s2")

g2a2 <- ggplot(sim) +
  geom_tile(aes(lon, lat, fill = temp_318m)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "Z2") +
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

g2b1.3 <- ggplot(predML3$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2b2.3 <- ggplot(predML3$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

gplot <- grid.arrange(g2a1, g2b1.1, g2b1.2, g2b1.3,
                      g2a2, g2b2.1, g2b2.2, g2b2.3,
                      nrow = 2)

ggsave(filename="figures/application_ocean_comparison_block_holdout_plot.png", plot=gplot, device="png", width=40, height=20, scale=1, units="cm", dpi=300)