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

# Simulate data for the experiment
# Simulate warped location
set.seed(11)
r1 <- 10
a1 <- 1
a2 <- 0
a3 <- 1
a4 <- 1 + 1*1i
layers <- c(AWU(r = r1, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = r1, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(res = 1L),
            LFT(a=c(a1,a2,a3,a4)))
nlayers <- length(layers)
eta <- list()
eta[[1]] <- sin(seq(0, pi, length.out = r1))
eta[[2]] <- c(1, rep(0, r1-1))
for(j in 3:(nlayers-1)) eta[[j]] <- runif(n = 1, min = -1, max = exp(3/2)/2)
eta[[nlayers]] <- 1

length.out <- 8
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped <- s
for(j in 1: (nlayers))
  swarped <- layers[[j]]$fR(swarped, eta[[j]]) %>% scal_0_5_mat()

# Plot the visualization of the homogenizing function
sw <- swarped
sw.df <- data.frame(x=sw[,1], y=sw[,2])

# Choose 3 points (randomly)
set.seed(1)
RNGkind(sample.kind = "Rounding")
sam <- sample(1:nrow(sw), 3)
l <- sam[1]
k <- sam[2]
m <- sam[3]

plot1 <- ggplot() +
  geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=2, alpha=1) + theme_bw() +
  coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.2,1.2)) +
  labs(x=expression(paste('f'[1], '(', bold(s), ')')), y=expression(paste('f'[2], '(', bold(s), ')'))) +
  theme(text = element_text(size=15)) +
  geom_point(data=dplyr::slice(sw.df, k), aes(x,y), col=rgb(1,0,0), size=2, alpha=1) +
  geom_point(data=dplyr::slice(sw.df, l), aes(x,y), col=rgb(1,0.8,0.2), size=2, alpha=1) +
  geom_point(data=dplyr::slice(sw.df, m), aes(x,y), col=rgb(0,1,0), size=2, alpha=1)

sw <- b1(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
plot2 <- ggplot() +
  geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=2, alpha=1) + theme_bw() +
  coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.2,1.2)) +
  labs(x=expression(paste('b'['1,1'],'(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['1,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15)) +
  geom_point(data=dplyr::slice(sw.df, k), aes(x,y), col=rgb(1,0,0), size=2, alpha=1) +
  geom_point(data=dplyr::slice(sw.df, l), aes(x,y), col=rgb(1,0.8,0.2), size=2, alpha=1) +
  geom_point(data=dplyr::slice(sw.df, m), aes(x,y), col=rgb(0,1,0), size=2, alpha=1)

sw <- b2(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
plot3 <- ggplot() +
  geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=2, alpha=1) + theme_bw() +
  coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.2,1.2)) +
  labs(x=expression(paste('b'['2,1'],'(', bold(b[1]), '(', bold(f), '(', bold(s), '))')), y=expression(paste('b'['2,2'],'(', bold(b[1]), '(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15)) +
  geom_point(data=dplyr::slice(sw.df, k), aes(x,y), col=rgb(1,0,0), size=2, alpha=1) +
  geom_point(data=dplyr::slice(sw.df, l), aes(x,y), col=rgb(1,0.8,0.2), size=2, alpha=1) +
  geom_point(data=dplyr::slice(sw.df, m), aes(x,y), col=rgb(0,1,0), size=2, alpha=1)

sw <- b3(sw, k, l, m)
sw.df <- data.frame(x=sw[,1], y=sw[,2])
plot4 <- ggplot() +
  geom_point(data=sw.df, aes(x,y), col=rgb(0,0,1), size=2, alpha=1) +
  theme_bw() + coord_cartesian(xlim = c(-1.2,1.2), ylim = c(-1.2,1.2)) +
  labs(x=expression(paste('b'['0,1'],'(', bold(f),, '(', bold(s), '))')), y=expression(paste('b'['0,2'],'(', bold(f), '(', bold(s), '))'))) +
  theme(text = element_text(size=15)) +
  geom_point(data=dplyr::slice(sw.df, k), aes(x,y), col=rgb(1,0,0), size=2, alpha=1) +
  geom_point(data=dplyr::slice(sw.df, l), aes(x,y), col=rgb(1,0.8,0.2), size=2, alpha=1) +
  geom_point(data=dplyr::slice(sw.df, m), aes(x,y), col=rgb(0,1,0), size=2, alpha=1)

gplot <- grid.arrange(plot1, plot2, plot3, plot4, nrow = 2)

ggsave(filename="figures/frame_homogenizing_function.png", plot=gplot, device="png", width=25, height=25, scale=1, units="cm", dpi=300)
