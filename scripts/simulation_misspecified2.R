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

# Simulate dataset
layers_asym <- c(AFF_2D(a=c(0.1, 1, 0, -0.1, 0, 1)))

RBF_Fouedjio <- function(x, theta) {
  
  theta11 <- theta[1]
  theta12 <- theta[2]
  
  sep1 <- (x[, 1, drop = FALSE] - theta11)
  sep2 <- (x[, 2, drop = FALSE] - theta12)
  sepsq <- sep1^2 + sep2^2
  
  theta + (x - theta) * matrix(rep(sqrt(sepsq), 2), ncol = 2)
}

length.out <- 101
s <- expand.grid(seq(-0.5, 0.5, length.out = length.out),
                 seq(-0.5, 0.5, length.out = length.out)) %>% as.matrix()
swarped1 <- swarped2 <- s

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

swarped1 <- RBF_Fouedjio(swarped1, theta = c(0,0))
swarped2 <- RBF_Fouedjio(swarped2, theta = c(0,0))

# Simulate symmetric, stationary covariance on the warped domain
D11 <- fields::rdist(swarped1)
D22 <- fields::rdist(swarped2)
D12 <- fields::rdist(swarped1, swarped2)

C11 <- fields::Matern(D11, range=0.2, nu=0.5, phi=1)
C12 <- fields::Matern(D12, range=0.2, nu=0.75, phi=0.8*0.9*1)
C22 <- fields::Matern(D22, range=0.2, nu=1, phi=0.9^2)
C <- rbind(cbind(C11, C12), cbind(t(C12), C22))
K <- t(chol(C))
i <- 2
set.seed(i)
y <- K %*% rnorm(nrow(s)*2)
y1 <- y[1: nrow(s),]
y2 <- y[(nrow(s)+1) : (nrow(s)*2),]
z1 <- y1 + 0.2*rnorm(length(y1))
z2 <- y2 + 0.1*rnorm(length(y2))
z <- c(z1, z2)
df <- data.frame(s, y1, y2, z1, z2)
names(df) <- c("s1", "s2", "y1", "y2", "z1", "z2")

# Save data
save(df, file="results/simulation_misspecified2_dataset.rda")

# Sample a subset of data for the experiment
RNGkind(sample.kind = "Rounding")
set.seed(1)
sam2 <- sample(1:nrow(df), 1000)
train_data <- df[sam2,]
test_data <- dplyr::setdiff(df, train_data)

df <- dplyr::select(train_data, s1, s2, z1, z2)
newdata <- dplyr::select(test_data, s1, s2)
groups_data <- split(newdata, (seq(nrow(newdata))-1) %/% (nrow(newdata)/5))

layers <- c(AWU(r = 50, dim = 1L, grad = 200, lims = c(-0.5, 0.5)),
            AWU(r = 50, dim = 2L, grad = 200, lims = c(-0.5, 0.5)),
            RBF_block(),
            LFT())

layers_asym <- c(AFF_2D())

# Fit Model (nonstationary, asymmetric)
t1 <- proc.time()
d2 <- deepspat_bivar_GP(f = z1 + z2 ~ s1 + s2 - 1, data = df, g = ~ 1,
                        layers = layers, layers_asym = layers_asym,
                        family = "matern_nonstat_asymm",
                        method = "REML", nsteps = 150L
)
t2 <- proc.time()
time2 <- (t2 - t1)[3]

# Save parameter estimates
eta <- d2$run(d2$eta_tf)
AFF_pars <- d2$run(d2$layers_asym[[1]]$pars)
LFT_pars <- d2$run(d2$layers[[12]]$pars)
scalings <- d2$run(d2$scalings)
scalings_asym <- d2$run(d2$scalings_asym)
nu_1 <- d2$run(d2$nu_tf_1)
nu_2 <- d2$run(d2$nu_tf_2)
sigma2_1 <- d2$run(d2$sigma2_tf_1)
sigma2_2 <- d2$run(d2$sigma2_tf_2)
sigma2_12 <- d2$run(d2$sigma2_tf_12)
l <- as.numeric(d2$run(d2$l_tf_1))
precy_1 <- d2$run(d2$precy_tf_1)
precy_2 <- d2$run(d2$precy_tf_2)
s_warped1 <- d2$run(d2$swarped_tf1)
s_warped2 <- d2$run(d2$swarped_tf2)
beta <- d2$run(d2$beta)
parameter_est <- list(eta, AFF_pars, LFT_pars, scalings, scalings_asym,
                      nu_1, nu_2, sigma2_1, sigma2_2, sigma2_12, l, precy_1, precy_2,
                      s_warped1, s_warped2, beta)
save(parameter_est, file = paste0("results/simulation_misspecified2_parameter_est.rda"))
