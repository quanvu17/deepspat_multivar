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

# Load packages
library("ggplot2")
library("gridExtra")

# Load data and predictions
load("results/simulation_asymm_dataset.rda")
load("results/simulation_asymm_predictions.rda")
predML2 <- predML[[1]]
predML3 <- predML[[2]]
predML4 <- predML[[3]]

RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
df2 <- df[sam2,]

df_pred2 <- data.frame(predML2$df_pred, y1=df$y1, y2=df$y2)
df_pred3 <- data.frame(predML3$df_pred, y1=df$y1, y2=df$y2)
df_pred4 <- data.frame(predML4$df_pred, y1=df$y1, y2=df$y2)

# Plot for the comparison
b1min <- min(c(predML2$df_pred$pred_mean_1, predML3$df_pred$pred_mean_1, predML4$df_pred$pred_mean_1, df$y1))
b1max <- max(c(predML2$df_pred$pred_mean_1, predML3$df_pred$pred_mean_1, predML4$df_pred$pred_mean_1, df$y1))
b2min <- min(c(predML2$df_pred$pred_mean_2, predML3$df_pred$pred_mean_2, predML4$df_pred$pred_mean_2, df$y2))
b2max <- max(c(predML2$df_pred$pred_mean_2, predML3$df_pred$pred_mean_2, predML4$df_pred$pred_mean_2, df$y2))
c1min <- min(c(sqrt(predML2$df_pred$pred_var_1), sqrt(predML3$df_pred$pred_var_1), sqrt(predML4$df_pred$pred_var_1)))
c1max <- max(c(sqrt(predML2$df_pred$pred_var_1), sqrt(predML3$df_pred$pred_var_1), sqrt(predML4$df_pred$pred_var_1)))
c2min <- min(c(sqrt(predML2$df_pred$pred_var_2), sqrt(predML3$df_pred$pred_var_2), sqrt(predML4$df_pred$pred_var_2)))
c2max <- max(c(sqrt(predML2$df_pred$pred_var_2), sqrt(predML3$df_pred$pred_var_2), sqrt(predML4$df_pred$pred_var_2)))

g2a1 <- ggplot(df) +
  geom_tile(aes(s1, s2, fill = y1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "Y1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2a2 <- ggplot(df) +
  geom_tile(aes(s1, s2, fill = y2))+
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "Y2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b1.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b1.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b2.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b2.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b1.3 <- ggplot(predML3$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b1.3 <- ggplot(predML3$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b2.3 <- ggplot(predML3$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b2.3 <- ggplot(predML3$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b1.4 <- ggplot(predML4$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b1.4 <- ggplot(predML4$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b2.4 <- ggplot(predML4$df_pred) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b2.4 <- ggplot(predML4$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2a <- ggplot(df2) +
  geom_point(aes(s1, s2, color = "black")) +
  scale_color_identity(name = "", breaks = c("black"), labels = c("locs"), guide = "legend") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

gplot <- grid.arrange(g2a1, g2b1.2, g2b1.3, g2b1.4,
                      g2a, g3b1.2, g3b1.3, g3b1.4,
                      g2a2, g2b2.2, g2b2.3, g2b2.4,
                      g2a, g3b2.2, g3b2.3, g3b2.4,
                      nrow = 4
)

ggsave(filename="figures/simulation_asymm_comparison_plot.png", plot=gplot, device="png", width=40, height=40, scale=1, units="cm", dpi=300)
