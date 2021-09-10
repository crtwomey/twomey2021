#
# rd.R
# Copyright (c) 2016-2021 Colin Twomey
# All rights reserved.
#
# This file is part of Twomey2021.
#
# Twomey2021 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Twomey2021 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Twomey2020  If not, see <https://www.gnu.org/licenses/>.
#

#?
#  Minimal demonstration of the inverse inference algorithm, based on
#  the example given in Twomey et al. (2021) SI Fig. B1.
#

source('rd.R')
source('inv.R')

# setup example using a 100 x 100 grid
nx <- 100
x  <- cbind(rep(1:nx, each=nx), rep(1:nx, nx)) / nx
nx <- nrow(x)

# create (arbitrary) ground-truth p(x)
px <- matrix((1.5*(x[,1] - 0.5))^3 + x[,2]^2, nrow(x), 1)
px <- px - ifelse(min(px) < 0, min(px), 0) + 0.01
px <- exp(2*log(px / sum(px)))
px <- px / sum(px)

# rate-distortion parameters
beta  <- 100
tmax  <- 500
nxhat <- 8

# get rate-distortion optimal xhat
rd <- init_rd(x, nxhat, beta, px)
rd <- run_rd(rd, tmax)

# inverse inference (true p(x) unknown)
inv <- init_inv(x, rd$xhat, rd$pxhat)
inv <- run_inv(inv, tmax)

# use a standard perceptually uniform color scale
library(viridis)
color_scale <- magma(500)

# plot the results side-by-side
op <- par(mar=c(3,3,3,1), mfrow=c(1,2))
image(matrix(px, sqrt(nx), byrow=TRUE),
	useRaster = TRUE,
	col       = color_scale,
	main      = 'ground truth'
)
image(matrix(inv$qx, sqrt(nx), byrow=TRUE),
	useRaster = TRUE,
	col       = color_scale,
	main      = 'inferred'
)
par(op)

# compute the KL-divergence between ground truth and inferred
message('KL-divergence = ', sum(px*log(px/inv$qx)))

