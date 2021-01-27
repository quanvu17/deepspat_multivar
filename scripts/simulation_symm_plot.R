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

# Load data and predictions
load("results/simulation_symm_crossvalidation.rda")
load("results/simulation_symm_dataset.rda")
load("results/simulation_symm_predictions.rda")
predML1 <- predML[[1]]
predML2 <- predML[[2]]
load("results/simulation_symm_predictions_gap.rda")
predML3 <- predML[[1]]
predML4 <- predML[[2]]

## Boxplot

results_summary <- data.frame(rmse1 = c(crossvalid[[1]], crossvalid[[7]]),
                              crps1 = c(crossvalid[[2]], crossvalid[[8]]),
                              rmse2 = c(crossvalid[[3]], crossvalid[[9]]),
                              crps2 = c(crossvalid[[4]], crossvalid[[10]]),
                              model = c(rep(1, 30), rep(2, 30)))

boxplot1 <- ggplot() + geom_boxplot(data = results_summary, aes(x = as.factor(model), y = rmse1)) +
  theme_bw() + labs(x='Model (4.x)', y='RMSPE1') +
  theme(text = element_text(size=15))

boxplot2 <- ggplot() + geom_boxplot(data = results_summary, aes(x = as.factor(model), y = crps1)) +
  theme_bw() + labs(x='Model (4.x)', y='CRPS1') +
  theme(text = element_text(size=15))

boxplot3 <- ggplot() + geom_boxplot(data = results_summary, aes(x = as.factor(model), y = rmse2)) +
  theme_bw() + labs(x='Model (4.x)', y='RMSPE2') +
  theme(text = element_text(size=15))

boxplot4 <- ggplot() + geom_boxplot(data = results_summary, aes(x = as.factor(model), y = crps2)) +
  theme_bw() + labs(x='Model (4.x)', y='CRPS2') +
  theme(text = element_text(size=15))

gplot <- grid.arrange(boxplot1, boxplot2, boxplot3, boxplot4, nrow = 1)

ggsave(filename="figures/simulation_symm_boxplot.png", plot=gplot, device="png", width=40, height=10, scale=1, units="cm", dpi=300)

###################

RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
df2 <- df[sam2,]

df1 <- df[df$s1 > -0.28 & df$s1 < -0.08 & df$s2 > -0.48 & df$s2 < -0.28,]
df0 <- setdiff(df,df1)
RNGkind(sample.kind = "Rounding")
set.seed(2)
sam <- sample(1:nrow(df0), 1000)
df3 <- df0[sam,]

# Plot for the comparision
b1min <- min(c(predML1$pred_mean_1, predML2$pred_mean_1, predML3$pred_mean_1, predML4$pred_mean_1, df$y1))
b1max <- max(c(predML1$pred_mean_1, predML2$pred_mean_1, predML3$pred_mean_1, predML4$pred_mean_1, df$y1))
b2min <- min(c(predML1$pred_mean_2, predML2$pred_mean_2, predML3$pred_mean_2, predML4$pred_mean_2, df$y2))
b2max <- max(c(predML1$pred_mean_2, predML2$pred_mean_2, predML3$pred_mean_2, predML4$pred_mean_2, df$y2))
c1min <- min(c(sqrt(predML1$pred_var_1), sqrt(predML2$pred_var_1), sqrt(predML3$pred_var_1), sqrt(predML4$pred_var_1)))
c1max <- max(c(sqrt(predML1$pred_var_1), sqrt(predML2$pred_var_1), sqrt(predML3$pred_var_1), sqrt(predML4$pred_var_1)))
c2min <- min(c(sqrt(predML1$pred_var_2), sqrt(predML2$pred_var_2), sqrt(predML3$pred_var_2), sqrt(predML4$pred_var_2)))
c2max <- max(c(sqrt(predML1$pred_var_2), sqrt(predML2$pred_var_2), sqrt(predML3$pred_var_2), sqrt(predML4$pred_var_2)))

g2a1 <- ggplot(df) +
  geom_tile(aes(s1, s2, fill = y1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "Y1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2a2 <- ggplot(df) +
  geom_tile(aes(s1, s2, fill = y2))+
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "Y2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b1.1 <- ggplot(predML1) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b1.1 <- ggplot(predML1) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b2.1 <- ggplot(predML1) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b2.1 <- ggplot(predML1) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b1.2 <- ggplot(predML2) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b1.2 <- ggplot(predML2) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2b2.2 <- ggplot(predML2) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g3b2.2 <- ggplot(predML2) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

g2a <- ggplot(df2) +
  geom_point(aes(s1, s2, color = "black")) +
  scale_color_identity(name = "", breaks = c("black"), labels = c("locs"), guide = "legend") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

gplot <- grid.arrange(g2a1, g2b1.1, g2b1.2,
                      g2a, g3b1.1, g3b1.2,
                      g2a2, g2b2.1, g2b2.2,
                      g2a, g3b2.1, g3b2.2,
                      nrow = 4
)

ggsave(filename="figures/simulation_symm_comparison_plot.png", plot=gplot, device="png", width=30, height=40, scale=1, units="cm", dpi=300)

# Plot for the comparision when there is a missing block
dt <- data.frame(xmin=-0.28, xmax=-0.08, ymin=-0.48, ymax=-0.28, b=0) 
a <-  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=b), inherit.aes=FALSE, color="black", alpha=0.5)

g2a1 <- ggplot(df) +
  geom_tile(aes(s1, s2, fill = y1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "Y1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2a2 <- ggplot(df) +
  geom_tile(aes(s1, s2, fill = y2))+
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "Y2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2b1.1 <- ggplot(predML3) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b1.1 <- ggplot(predML3) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2b2.1 <- ggplot(predML3) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b2.1 <- ggplot(predML3) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2b1.2 <- ggplot(predML4) +
  geom_tile(aes(s1, s2, fill = pred_mean_1)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b1min, b1max), name = "pred1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b1.2 <- ggplot(predML4) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2b2.2 <- ggplot(predML4) +
  geom_tile(aes(s1, s2, fill = pred_mean_2)) +
  scale_fill_distiller(palette = "Spectral", limits=c(b2min, b2max), name = "pred2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b2.2 <- ggplot(predML4) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2a <- ggplot(df3) +
  geom_point(aes(s1, s2, color = "black")) +
  scale_color_identity(name = "", breaks = c("black"), labels = c("locs"), guide = "legend") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

gplot <- grid.arrange(g2a1, g2b1.1, g2b1.2,
                      g2a, g3b1.1, g3b1.2,
                      g2a2, g2b2.1, g2b2.2,
                      g2a, g3b2.1, g3b2.2,
                      nrow = 4
)

ggsave(filename="figures/simulation_symm_comparison_block_holdout_plot.png", plot=gplot, device="png", width=30, height=40, scale=1, units="cm", dpi=300)


# Plot the true, estimated and boostrap warping function
# Pick 3 points
k <- 822
l <- 127
m <- 121

# True warping
load("results/simulation_symm_dataset.rda")
set.seed(1)
r1 <- 50
a1 <- 1
a2 <- 0
a3 <- 1
a4 <- 1 + 1*1i
layers <- c(AWU(r = r1, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = r1, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(res = 1L),
            LFT(a=c(a1, a2, a3, a4)))
nlayers <- length(layers)
eta <- list()
eta[[1]] <- sin(seq(0, pi, length.out = r1))
eta[[2]] <- c(1, rep(0, r1-1))
for(j in 3:(nlayers-1)) eta[[j]] <- runif(n = 1, min = -1, max = exp(3/2)/2)
eta[[nlayers]] <- 1

RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
s <- dplyr::select(df[sam2,], s1, s2) %>% as.matrix

swarped <- s
for(j in 1: (nlayers))
  swarped <- layers[[j]]$fR(swarped, eta[[j]]) %>% scal_0_5_mat()

sw <- swarped
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
warp.plot1 <- ggplot() + geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=0.5, alpha=1) +
  theme_bw() + coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

# Estimated warping
load("results/simulation_symm_parameter_estimate.rda")
sw <- parameter_est[[12]]
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
warp.plot2 <- ggplot() + geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=0.5, alpha=1) +
  theme_bw() + coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

# Boostrap warpings
g <- seq(0.2, 0.7, length.out = 100)
i = 1
load(paste0("results/simulation_symm_bootstrap_parameter_estimate/bootstrap_par_est_sample", i, ".rda"))
sw <- parameter_est[[12]]
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
warp.plot3 <- ggplot(sw.df) + geom_point(aes(x,y), col=rgb(0,0,1), size=0.05, alpha=0.05) +
  theme_bw() + coord_cartesian(xlim = c(-2,2), ylim = c(-2,2)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

for (i in 2:100){
  load(paste0("results/simulation_symm_bootstrap_parameter_estimate/bootstrap_par_est_sample", i, ".rda"))
  sw <- parameter_est[[12]]
  sw <- b0(sw, k, l, m)
  sw.df <- data.frame(x=sw[,1], y=sw[,2])
  warp.plot3 <- warp.plot3 + geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=0.05, alpha=0.05)
}

original.plot <- ggplot(df2) +
  geom_point(aes(s1, s2), size=0.5) +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  labs(x=expression(paste('s'['1'])), y=expression(paste('s'['2'])))

gplot <- grid.arrange(original.plot, warp.plot1, warp.plot2, warp.plot3, nrow = 2)

ggsave(filename="figures/simulation_symm_warping_plot.png", plot=gplot, device="png", width=20, height=20, scale=1, units="cm", dpi=300)