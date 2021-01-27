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

# Plot the true, estimated and boostrap warping function
# Pick 3 points
k <- 822
l <- 127
m <- 121

# True warping
load("results/simulation_trivar_dataset.rda")
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
  theme_bw() + coord_cartesian(xlim = c(-2.5,2.5), ylim = c(-2.5,2.5)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

# Estimated warping
load("results/simulation_trivar_parameter_estimate.rda")
sw <- parameter_est[[17]]
sw <- b0(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
warp.plot2 <- ggplot() + geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=0.5, alpha=1) +
  theme_bw() + coord_cartesian(xlim = c(-2.5,2.5), ylim = c(-2.5,2.5)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15))

gplot <- grid.arrange(warp.plot1, warp.plot2, nrow = 1)

ggsave(filename="figures/simulation_trivar_warping_plot.png", plot=gplot, device="png", width=20, height=10, scale=1, units="cm", dpi=300)