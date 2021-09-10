#
# inv.R
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
# along with Twomey2021.  If not, see <https://www.gnu.org/licenses/>.
#

#?
#  Implementation of the inverse inference algorithm derived in
#  Twomey et al. (2021) SI Sec. B.
#

# initialize solution
init_inv <- function(
	x,         # locations of the points corresponding to q(x)
	xhat,      # locations of centroids corresponding to p(xhat)
	pxhat,     # compressed representation frequencies
	tol = Inf) # maximum allowed magnitude for vhats
{
	nx      <- nrow(x)
	nd      <- ncol(x)
	nxhat   <- nrow(xhat)
	qx_xhat <- matrix(1/nx, nx, nxhat)
	qxhat_x <- matrix(1/nxhat, nx, nxhat)
	qx      <- rep(1/nx, nx)
	xtilde  <- xhat
	vhat    <- matrix(0, nxhat, nd)
	return(list(
		x       = x,
		qx_xhat = qx_xhat,
		qxhat_x = qxhat_x,
		qx      = qx,
		xtilde  = xtilde,
		vhat    = vhat,
		xhat    = xhat,
		pxhat   = pxhat,
		tol     = tol
	))
}

# iterative update step
update_inv <- function(inv) {
	x    <- inv$x
	xhat <- inv$xhat

	# log sum exp trick to avoid underflow
	lse <- function(u) {
		i <- head(which(u == max(u)),1)
		u[i] + log1p(sum(exp(u[-i] - u[i])))
	}

	# exp normalize trick to avoid overflow
	expnorm <- function(u) {
		umax <- apply(u, 2, max)
		expu <- exp(sweep(u, 2, umax, '-'))
		sweep(expu, 2, colSums(expu), '/')
	}

	# find optimal vhats
	inv$vhat <- t(sapply(1:nrow(xhat), function(i) {
		xtilde <- function(v) {
			qx_xhat <- as.vector(inv$qxhat_x[,i]*exp(inv$x %*% v))
			qx_xhat <- qx_xhat / sum(qx_xhat)
			colSums(inv$x * qx_xhat)
		}
		grad <- function(v) xtilde(v) - inv$xhat[i,]
		G    <- function(v) lse(inv$x %*% v + log(inv$qxhat_x[,i]))
		phi  <- function(v) G(v) - sum(inv$xhat[i,]*v)
		opt  <- optim(inv$vhat[i,], phi, grad,
			method  = 'BFGS',
			hessian = TRUE
		)
		return(opt$par)
	}))

	# restrict maximum magnitude (when inv$tol < Inf)
	m        <- sqrt(rowSums(inv$vhat^2))
	inv$vhat <- pmin(m, inv$tol) * inv$vhat / m

	# update qx_xhat
	inv$qx_xhat <- expnorm(inv$x %*% t(inv$vhat) + log(inv$qxhat_x))

	# update xtilde (depends on qx_xhat)
	inv$xtilde <- t(inv$qx_xhat) %*% x

	# update qx (depends on qx_xhat)
	inv$qx <- colSums(t(inv$qx_xhat) * inv$pxhat)

	# update qxhat_x (depends on qx_xhat and qx)
	inv$qxhat_x <- sweep(inv$qx_xhat, 2, inv$pxhat, '*') / inv$qx

	return(inv)
}

# run inverse inference for a given number of iterations (tmax)
run_inv <- function(inv, tmax) {
	for (t in 1:tmax) {
		cat('\r', 't = ', t, ' / ', tmax)
		flush.console()
		inv <- update_inv(inv)
	}
	cat('\n')
	flush.console()
	return(inv)
}


