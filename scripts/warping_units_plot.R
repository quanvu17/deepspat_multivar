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

# Affine transformation

layers_asym <- c(AFF_2D(a=c(0.15, 0.95, 0.05, -0.15, 0, 1.05)))
nlayers_asym <- length(layers_asym)
length.out <- 21
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped1 <- swarped2 <- s

for(j in 1: (nlayers_asym)){
  swarped1 <- swarped1
  swarped2 <- layers_asym[[j]]$fR(swarped2, eta[[j]])
  swarped1[,1] <- (swarped1[,1] - min(c(swarped1[,1], swarped2[,1]))) / (max(c(swarped1[,1], swarped2[,1])) - min(c(swarped1[,1], swarped2[,1]))) - 0.5
  swarped1[,2] <- (swarped1[,2] - min(c(swarped1[,2], swarped2[,2]))) / (max(c(swarped1[,2], swarped2[,2])) - min(c(swarped1[,2], swarped2[,2]))) - 0.5
  swarped2[,1] <- (swarped2[,1] - min(c(swarped1[,1], swarped2[,1]))) / (max(c(swarped1[,1], swarped2[,1])) - min(c(swarped1[,1], swarped2[,1]))) - 0.5
  swarped2[,2] <- (swarped2[,2] - min(c(swarped1[,2], swarped2[,2]))) / (max(c(swarped1[,2], swarped2[,2])) - min(c(swarped1[,2], swarped2[,2]))) - 0.5
}

df = data.frame(sw1.1 = swarped1[,1], sw1.2 = swarped1[,2],
                sw2.1 = swarped2[,1], sw2.2 = swarped2[,2])

warpings <- ggplot() + geom_point(data=df, aes(x=sw1.1, y=sw1.2), colour = "black") +
  geom_point(data=df, aes(x=sw2.1, y=sw2.2), colour = "red") +
  theme_bw() +
  xlab(expression(paste("g"[1], "(s)"))) + ylab(expression(paste("g"[2], "(s)")))

ggsave("figures/warp_11.png", plot=warpings, device="png", width=10, height=10, scale=1, units="cm", dpi=300)

###########

# Axial warping

set.seed(1)
layers <- c(AWU(r = 50L, dim = 1L, grad = 200, lims = c(-0.5, 0.5))
)

nlayers <- length(layers)

eta <- list()
eta[[1]] <- rnorm(50)

length.out <- 21
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped1 <- s

for(j in 1: (nlayers)){
  swarped1 <- layers[[j]]$fR(swarped1, eta[[j]]) %>% scal_0_5_mat()
}

swarped2 <- swarped1

df = data.frame(sw1.1 = swarped1[,1], sw1.2 = swarped1[,2],
                sw2.1 = swarped2[,1], sw2.2 = swarped2[,2])

warpings <- ggplot() + geom_point(data=df, aes(x=sw1.1, y=sw1.2), colour = "blue") +
  theme_bw() +
  xlab(expression(paste("f"[1], "(s)"))) + ylab(expression(paste("f"[2], "(s)")))

ggsave("figures/warp_21.png", plot=warpings, device="png", width=10, height=10, scale=1, units="cm", dpi=300)

###########

# Radial basis function

length.out <- 21
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped1 <- s

swarped1 <- RBF(s, c(0, 0, -2))

swarped2 <- swarped1

df = data.frame(sw1.1 = swarped1[,1], sw1.2 = swarped1[,2],
                sw2.1 = swarped2[,1], sw2.2 = swarped2[,2])

warpings <- ggplot() + geom_point(data=df, aes(x=sw1.1, y=sw1.2), colour = "blue") +
  theme_bw() +
  xlab(expression(paste("f"[1], "(s)"))) + ylab(expression(paste("f"[2], "(s)")))

ggsave("figures/warp_31.png", plot=warpings, device="png", width=10, height=10, scale=1, units="cm", dpi=300)

###########

# Mobius transformation

set.seed(1)
r1 <- 50
a1 <- 2
a2 <- 0
a3 <- 1
a4 <- 1 + 1*1i
layers <- c(LFT(a=c(a1,a2,a3,a4)))

nlayers <- length(layers)

eta[[nlayers]] <- 1

length.out <- 21
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped1 <- s

for(j in 1: (nlayers)){
  swarped1 <- layers[[j]]$fR(swarped1, eta[[j]]) %>% scal_0_5_mat()
}

swarped2 <- swarped1

df = data.frame(sw1.1 = swarped1[,1], sw1.2 = swarped1[,2],
                sw2.1 = swarped2[,1], sw2.2 = swarped2[,2])

warpings <- ggplot() + geom_point(data=df, aes(x=sw1.1, y=sw1.2), colour = "blue") +
  theme_bw() +
  xlab(expression(paste("f"[1], "(s)"))) + ylab(expression(paste("f"[2], "(s)")))

ggsave("figures/warp_41.png", plot=warpings, device="png", width=10, height=10, scale=1, units="cm", dpi=300)

###########
