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

test_data <- sim[sim$lat > 37.5 & sim$lat < 38.2,]
train_data <- setdiff(sim, test_data) 

# Plot for gap experiment
b1min <- min(c(predML1$df_pred$pred_mean_1, predML2$df_pred$pred_mean_1, predML3$df_pred$pred_mean_1, sim$temp_1m))
b1max <- max(c(predML1$df_pred$pred_mean_1, predML2$df_pred$pred_mean_1, predML3$df_pred$pred_mean_1, sim$temp_1m))
b2min <- min(c(predML1$df_pred$pred_mean_2, predML2$df_pred$pred_mean_2, predML3$df_pred$pred_mean_2, sim$temp_318m))
b2max <- max(c(predML1$df_pred$pred_mean_2, predML2$df_pred$pred_mean_2, predML3$df_pred$pred_mean_2, sim$temp_318m))
c1min <- min(c(sqrt(predML1$df_pred$pred_var_1), sqrt(predML2$df_pred$pred_var_1), sqrt(predML3$df_pred$pred_var_1)))
c1max <- max(c(sqrt(predML1$df_pred$pred_var_1), sqrt(predML2$df_pred$pred_var_1), sqrt(predML3$df_pred$pred_var_1)))
c2min <- min(c(sqrt(predML1$df_pred$pred_var_2), sqrt(predML2$df_pred$pred_var_2), sqrt(predML3$df_pred$pred_var_2)))
c2max <- max(c(sqrt(predML1$df_pred$pred_var_2), sqrt(predML2$df_pred$pred_var_2), sqrt(predML3$df_pred$pred_var_2)))

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

g3b1.1 <- ggplot(predML1$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b2.1 <- ggplot(predML1$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b1.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b2.2 <- ggplot(predML2$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b1.3 <- ggplot(predML3$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_1) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c1min, c1max), name = "s.e.1") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g3b2.3 <- ggplot(predML3$df_pred) +
  geom_tile(aes(s1, s2, fill = sqrt(pred_var_2) ))+
  scale_fill_distiller(palette = "BrBG", limits=c(c2min, c2max), name = "s.e.2") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  geom_rect(data=dt, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), inherit.aes=FALSE, color="black", alpha=0)

g2a <- ggplot(train_data) +
  geom_point(aes(lon, lat, color = "black")) +
  scale_color_identity(name = "", breaks = c("black"), labels = c("locs"), guide = "legend") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15))

gplot <- grid.arrange(g2a1, g2b1.1, g2b1.2, g2b1.3,
                      g2a, g3b2.1, g3b2.2, g3b2.3,
                      g2a2, g2b2.1, g2b2.2, g2b2.3,
                      g2a, g3b1.1, g3b1.2, g3b1.3,
                      nrow = 4)

ggsave(filename="figures/application_ocean_comparison_block_holdout_plot.png", plot=gplot, device="png", width=40, height=40, scale=1, units="cm", dpi=300)

# Plot the true and estimated warping function
# Pick 3 points
k <- 628
l <- 640
m <- 1280

# Estimated warping
load("results/application_ocean_parameter_estimate.rda")

layers_asym <- c(AFF_2D(a=c(unlist(parameter_est[[2]]))))
nlayers_asym <- length(layers_asym)

all_data <- sim
test_data <- sim[sim$lat > 37.5 & sim$lat < 38.2,]
train_data <- setdiff(sim, test_data) 

swarped1 <- swarped2 <- matrix(c(train_data$lon, train_data$lat), ncol=2)

minsw1 <- min(c(swarped1[,1], swarped2[,1]))
maxsw1 <- max(c(swarped1[,1], swarped2[,1]))
minsw2 <- min(c(swarped1[,2], swarped2[,2]))
maxsw2 <- max(c(swarped1[,2], swarped2[,2]))
swarped1[,1] <- (swarped1[,1] - minsw1) / (maxsw1 - minsw1) - 0.5
swarped1[,2] <- (swarped1[,2] - minsw2) / (maxsw2 - minsw2) - 0.5
swarped2[,1] <- (swarped2[,1] - minsw1) / (maxsw1 - minsw1) - 0.5
swarped2[,2] <- (swarped2[,2] - minsw2) / (maxsw2 - minsw2) - 0.5

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

swarped_asym1 <- data.frame(s1 = swarped1[,1], s2 = swarped1[,2])
swarped_asym2 <- data.frame(s1 = swarped2[,1], s2 = swarped2[,2])

warp.plot1 <- ggplot() + geom_point(data=swarped_asym1, aes(x=s1, y=s2), colour = "black", size = 0.5) +
  geom_point(data=swarped_asym2, aes(x=s1, y=s2), colour = "red", size = 0.5) +
  theme_bw() + coord_cartesian(xlim = c(-0.5,0.5), ylim = c(-0.5,0.5)) +
  labs(x=expression(paste('g'['1'],'(', bold(s), ')')), y=expression(paste('g'['2'],'(', bold(s), ')')))

sw <- parameter_est[[14]]
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
warp.plot2 <- ggplot() + geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=0.5, alpha=1) +
  theme_bw() + coord_cartesian(xlim = c(-3,1.5), ylim = c(-2.5,2.5)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

gplot <- grid.arrange(warp.plot1, warp.plot2, nrow = 1)

ggsave(filename="figures/application_ocean_warping_plot.png", plot=gplot, device="png", width=20, height=10, scale=1, units="cm", dpi=300)

################

## Contour plot

load("results/application_ocean_parameter_estimate.rda")

layers_asym <- c(AFF_2D(a=c(unlist(parameter_est[[2]]))))
nlayers_asym <- length(layers_asym)

length.out <- 36
s <- expand.grid(seq(-62.3, -60, length.out = length.out),
                 seq(36.4, 39.5, length.out = length.out)) %>% as.matrix()
swarped1 <- swarped2 <- s

minsw1 <- min(c(swarped1[,1], swarped2[,1]))
maxsw1 <- max(c(swarped1[,1], swarped2[,1]))
minsw2 <- min(c(swarped1[,2], swarped2[,2]))
maxsw2 <- max(c(swarped1[,2], swarped2[,2]))
swarped1[,1] <- (swarped1[,1] - minsw1) / (maxsw1 - minsw1) - 0.5
swarped1[,2] <- (swarped1[,2] - minsw2) / (maxsw2 - minsw2) - 0.5
swarped2[,1] <- (swarped2[,1] - minsw1) / (maxsw1 - minsw1) - 0.5
swarped2[,2] <- (swarped2[,2] - minsw2) / (maxsw2 - minsw2) - 0.5

for(j in 1: nlayers_asym){
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

layers <- c(AWU(r = 10L, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 10L, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT(a=c(as.numeric(parameter_est[[3]][1]) + 1i * as.numeric(parameter_est[[3]][5]),
                    as.numeric(parameter_est[[3]][2]) + 1i * as.numeric(parameter_est[[3]][6]),
                    as.numeric(parameter_est[[3]][3]) + 1i * as.numeric(parameter_est[[3]][7]),
                    as.numeric(parameter_est[[3]][4]) + 1i * as.numeric(parameter_est[[3]][8])))
)
nlayers <- length(layers)

eta <- parameter_est[[1]]

for(j in 1: (nlayers)){
  swarped1 <- layers[[j]]$fR(swarped1, as.numeric(eta[[j]])) %>% scal_0_5_mat()
  swarped2 <- layers[[j]]$fR(swarped2, as.numeric(eta[[j]])) %>% scal_0_5_mat()
}

center <- c((0:3)*9 + 4*36 + 5, (0:3)*9 + (9+4)*36 + 5, (0:3)*9 + (9*2+4)*36 + 5, (0:3)*9 + (9*3+4)*36 + 5)

idx = c(rep(c(rep(1,9), rep(2,9), rep(3,9), rep(4,9)), 9), 
        rep(c(rep(5,9), rep(6,9), rep(7,9), rep(8,9)), 9),
        rep(c(rep(9,9), rep(10,9), rep(11,9), rep(12,9)), 9),
        rep(c(rep(13,9), rep(14,9), rep(15,9), rep(16,9)), 9)
)

nu_1 <- parameter_est[[6]]
nu_2 <- parameter_est[[7]]
sigma2_1 <- parameter_est[[8]]
sigma2_2 <- parameter_est[[9]]
sigma2_12 <- parameter_est[[10]]
l <- parameter_est[[11]]

df.new <- data.frame(s1 = NULL, s2 = NULL,
                     C11 = NULL, C12 = NULL, C21 = NULL, C22 = NULL)

for (i in 1:16){
  D11 <- fields::rdist(swarped1[idx == i,], matrix(swarped1[center[i],], nrow=1))
  D22 <- fields::rdist(swarped2[idx == i,], matrix(swarped2[center[i],], nrow=1))
  D12 <- fields::rdist(swarped2[idx == i,], matrix(swarped1[center[i],], nrow=1))
  D21 <- fields::rdist(swarped1[idx == i,], matrix(swarped2[center[i],], nrow=1))
  
  C11 <- fields::Matern(D11, range=l, nu=nu_1, phi=sigma2_1)
  C12 <- fields::Matern(D12, range=l, nu=(nu_1 + nu_2)/2, phi=sigma2_12)
  C21 <- fields::Matern(D21, range=l, nu=(nu_1 + nu_2)/2, phi=sigma2_12)
  C22 <- fields::Matern(D22, range=l, nu=nu_2, phi=sigma2_2)
  
  df.new.i <- data.frame(s1 = s[idx == i, 1], s2 = s[idx == i, 2],
                         C11, C12, C21, C22)
  df.new <- rbind(df.new, df.new.i)
  
}

g1 <- ggplot() +
  geom_point(aes(s[center,1], s[center,2]), col="black", size=0.5, alpha=1) +
  geom_contour(data = df.new, aes(s1, s2, z = C11), binwidth = sigma2_1 * 0.4, col = "blue") +
  geom_contour(data = df.new, aes(s1, s2, z = C11), binwidth = sigma2_1 * 0.8, col = "red") +
  theme_bw() + coord_cartesian(xlim = c(-62.3, -60), ylim = c(36.4, 39.5)) +
  labs(x='lon', y='lat') +
  theme(text = element_text(size=15))

g2 <- ggplot() +
  geom_point(aes(s[center,1], s[center,2]), col="black", size=0.5, alpha=1) +
  geom_contour(data = df.new, aes(s1, s2, z = C22), binwidth = sigma2_2 * 0.4, col = "blue") +
  geom_contour(data = df.new, aes(s1, s2, z = C22), binwidth = sigma2_2 * 0.8, col = "red") +
  theme_bw() + coord_cartesian(xlim = c(-62.3, -60), ylim = c(36.4, 39.5)) +
  labs(x='lon', y='lat') +
  theme(text = element_text(size=15))

g12 <- ggplot() +
  geom_point(aes(s[center,1], s[center,2]), col="black", size=0.5, alpha=1) +
  geom_contour(data = df.new, aes(s1, s2, z = C12), binwidth = sigma2_12 * 0.4, col = "blue") +
  geom_contour(data = df.new, aes(s1, s2, z = C12), binwidth = sigma2_12 * 0.8, col = "red") +
  theme_bw() + coord_cartesian(xlim = c(-62.3, -60), ylim = c(36.4, 39.5)) +
  labs(x='lon', y='lat') +
  theme(text = element_text(size=15))

g21 <- ggplot() +
  geom_point(aes(s[center,1], s[center,2]), col="black", size=0.5, alpha=1) +
  geom_contour(data = df.new, aes(s1, s2, z = C21), binwidth = sigma2_12 * 0.4, col = "blue") +
  geom_contour(data = df.new, aes(s1, s2, z = C21), binwidth = sigma2_12 * 0.8, col = "red") +
  theme_bw() + coord_cartesian(xlim = c(-62.3, -60), ylim = c(36.4, 39.5)) +
  labs(x='lon', y='lat') +
  theme(text = element_text(size=15))

gplot <- grid.arrange(g1, g2, g12, g21, nrow = 2)

ggsave(filename="figures/application_ocean_heatmap.png", plot=gplot, device="png", width=20, height=20, scale=1, units="cm", dpi=300)