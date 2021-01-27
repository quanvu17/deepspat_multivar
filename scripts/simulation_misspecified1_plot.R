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

# Plot the true and estimated warping function
# Pick 3 points
k <- 822
l <- 127
m <- 121

# Misspecification 1 (stationary, symmetric data)

# True warping
load("results/simulation_misspecified1_dataset.rda")

RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
swarped1 <- swarped2 <- dplyr::select(df[sam2,], s1, s2) %>% as.matrix

swarped_asym1 <- data.frame(s1 = swarped1[,1], s2 = swarped1[,2])
swarped_asym2 <- data.frame(s1 = swarped2[,1], s2 = swarped2[,2])

warp.plot1 <- ggplot() + geom_point(data=swarped_asym1, aes(x=s1, y=s2), colour = "black", size = 1) +
  geom_point(data=swarped_asym2, aes(x=s1, y=s2), colour = "red", size = 0.25) +
  theme_bw() + coord_cartesian(xlim = c(-0.5,0.5), ylim = c(-0.5,0.5)) +
  labs(x=expression(paste('g'['1'],'(', bold(s), ')')), y=expression(paste('g'['2'],'(', bold(s), ')')))

sw <- swarped1
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
warp.plot2 <- ggplot() + geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=0.5, alpha=1) +
  theme_bw() + coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.2,1.2)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

# Estimated warping
load("results/simulation_misspecified1_parameter_est.rda")

layers_asym <- c(AFF_2D(a=c(unlist(parameter_est[[2]]))))

RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
swarped1 <- swarped2 <- dplyr::select(df[sam2,], s1, s2) %>% as.matrix

for(j in 1: (nlayers_asym)){
  swarped1 <- swarped1
  swarped2 <- layers_asym[[j]]$fR(swarped2, eta[[j]])
  swarped1[,1] <- (swarped1[,1] - min(c(swarped1[,1], swarped2[,1]))) / (max(c(swarped1[,1], swarped2[,1])) - min(c(swarped1[,1], swarped2[,1]))) - 0.5
  swarped1[,2] <- (swarped1[,2] - min(c(swarped1[,2], swarped2[,2]))) / (max(c(swarped1[,2], swarped2[,2])) - min(c(swarped1[,2], swarped2[,2]))) - 0.5
  swarped2[,1] <- (swarped2[,1] - min(c(swarped1[,1], swarped2[,1]))) / (max(c(swarped1[,1], swarped2[,1])) - min(c(swarped1[,1], swarped2[,1]))) - 0.5
  swarped2[,2] <- (swarped2[,2] - min(c(swarped1[,2], swarped2[,2]))) / (max(c(swarped1[,2], swarped2[,2])) - min(c(swarped1[,2], swarped2[,2]))) - 0.5
}

swarped_asym1 <- data.frame(s1 = swarped1[,1], s2 = swarped1[,2])
swarped_asym2 <- data.frame(s1 = swarped2[,1], s2 = swarped2[,2])

warp.plot3 <- ggplot() + geom_point(data=swarped_asym1, aes(x=s1, y=s2), colour = "black", size = 1) +
  geom_point(data=swarped_asym2, aes(x=s1, y=s2), colour = "red", size = 0.25) +
  theme_bw() + coord_cartesian(xlim = c(-0.5,0.5), ylim = c(-0.5,0.5)) +
  labs(x=expression(paste('g'['1'],'(', bold(s), ')')), y=expression(paste('g'['2'],'(', bold(s), ')')))
  
sw <- parameter_est[[14]]
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
warp.plot4 <- ggplot() + geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=0.5, alpha=1) +
  theme_bw() + coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.2,1.2)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

gplot <- grid.arrange(warp.plot1, warp.plot2, warp.plot3, warp.plot4, nrow = 2)

ggsave(filename="figures/simulation_misspecified1_warping_plot.png", plot=gplot, device="png", width=20, height=20, scale=1, units="cm", dpi=300)