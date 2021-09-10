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
# along with Twomey2021  If not, see <https://www.gnu.org/licenses/>.
#

#?
#  Reference implementation of the Blahut-Arimoto algorithm
#  for rate-distortion quantization.
#
#  Blahut-Arimoto with an Expectation-Maximization-like step
#  for computing optimal quantization centers. Achieves a
#  local optimum for Bregman divergences (e.g. squared
#  Euclidean distance and KL-divergence).
#

# all-pairs distances
squared_euclidean_distance <- function(x, y) {
	s <- matrix(0, nrow(x), nrow(y))
	for (i in 1:ncol(x)) {
		s <- s + outer(x[,i], y[,i], "-")^2
	}
	return(s)
}

fqxhat_x <- function(beta, pxhat, x, xhat, d) {
	pxhat * exp(-beta * t(d(x,xhat)))
}

fpxhat_x <- function(beta, pxhat, x, xhat, d) {
	qxhat_x <- fqxhat_x(beta, pxhat, x, xhat, d)
	qxhat_x <- qxhat_x + 1E-16
	Zxhat   <- colSums(qxhat_x)
	sweep(qxhat_x, 2, Zxhat, "/")
}

# initialize solution
init_rd <- function(x,
	nxhat = 5,
	beta  = 50,
	px    = matrix(1/nrow(x), nrow(x), 1),
	d     = squared_euclidean_distance)
{
	# random initialization of RDBC centroids
	s       <- sample(1:nrow(x), nxhat)
	xhat    <- as.matrix(x[s,], nxhat, ncol(x))
	nx      <- length(px)
	pxhat   <- rep(1/nxhat, nxhat)
	pxhat_x <- matrix(1/nxhat, nxhat, nx)
	return(list(
		nxhat   = nxhat,
		beta    = beta,
		d       = d,
		pxhat_x = pxhat_x,
		pxhat   = pxhat,
		xhat    = xhat,
		x       = x,
		px      = px
	))
}

# iterative update step
update_rd <- function(rd) {
	rd$pxhat_x <- fpxhat_x(rd$beta, rd$pxhat, rd$x, rd$xhat, rd$d)
	rd$pxhat   <- as.vector(rd$pxhat_x %*% rd$px)
	px_xhat    <- sweep(rd$pxhat_x, 2, rd$px, "*")
	px_xhat    <- sweep(px_xhat, 1, rowSums(px_xhat), "/")
	rd$xhat    <- px_xhat %*% rd$x
	return(rd)
}

# run rate-distortion updates for a given number of iterations (tmax)
run_rd <- function(rd, tmax) {
	for (t in 1:tmax) {
		cat('\r', 't = ', t, ' / ', tmax)
		flush.console()
		rd <- update_rd(rd)
	}
	cat('\n')
	flush.console()
	return(rd)
}

