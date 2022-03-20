## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2022 Stefan Widgren
##
## SimInf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SimInf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

KLIEP_kernel <- function(x, y, sigma) {
    i <- rep(seq_len(nrow(x)), nrow(y))
    j <- rep(seq_len(nrow(y)), each = nrow(x))
    matrix(exp(-rowSums((x[i, , drop = FALSE] - y[j, , drop = FALSE])^2) /
               (2*sigma^2)), nrow = nrow(x), ncol = nrow(y))
}

##' @importFrom MASS ginv
##' @noRd
KLIEP_alpha <- function(alpha, mean_xde) {
    alpha <- alpha + mean_xde * (1 - (t(mean_xde) %*% alpha)[1, 1]) *
        ginv(t(mean_xde) %*% mean_xde)[1, 1]
    alpha <- pmax(alpha, 0)
    alpha * ginv(t(mean_xde) %*% alpha)[1, 1]
}

KLIEP_learning <- function(mean_xde, xnu) {
    alpha <- matrix(1, ncol(xnu))
    alpha <- KLIEP_alpha(alpha, mean_xde)
    xnu_alpha <- xnu %*% alpha
    score <- mean(log(xnu_alpha))

    for (epsilon in 10^seq(from = 3, to = -3, by = -1)) {
        for (i in seq_len(100)) {
            alpha_new <- alpha + (epsilon * t(xnu)) %*% (1 / (xnu_alpha))
            alpha_new <- KLIEP_alpha(alpha_new, mean_xde)
            xnu_alpha_new <- xnu %*% alpha_new
            score_new <- mean(log(xnu_alpha_new))

            if(score_new <= score)
                break

            alpha <- alpha_new
            xnu_alpha <- xnu_alpha_new
            score <- score_new
        }
    }

    alpha
}

##' Kullback-Leibler Importance Estimation Procedure
##'
##' Kullback-Leiblar importance estimation procedure (with cross
##' validation)
##'
##' Estimating ratio of probability densities
##' \deqn{\frac{p_{nu}(x)}{p_{de}(x)}} from samples \eqn{{xde_i |
##' xde_i\in R^{d}}_{i=1}^{n_{de}}} drawn independently from
##' \eqn{p_{de}(x)} and samples \eqn{{xnu_i | xnu_i\in
##' R^{d}}_{i=1}^{n_{nu}}} drawn independently from \eqn{p_{nu}(x)}.
##' @param xnu numeric matrix with data for the numerator
##'     distribution.
##' @param xde numeric matrix with data for the denominator
##'     distribution.
##' @param fold number of folds for the cross validation. Default is
##'     5.
##' @param kernels number of kernels. Default is 100.
##' @return a list with centers, sigma, and weights.
##' @references
##'
##' \Sugiyama2008
##' @noRd
KLIEP <- function(xnu, xde, fold = 5L, kernels = 100L) {
    if (!is.matrix(xnu) || !is.numeric(xnu))
        stop("'xnu' must be a numeric matrix.", call. = FALSE)
    if (!is.matrix(xde) || !is.numeric(xde))
        stop("'xde' must be a numeric matrix.", call. = FALSE)
    if (ncol(xnu) != ncol(xde))
        stop("'xnu' and 'xde' must have the same dimension.", call. = FALSE)

    kernels <- min(kernels, nrow(xnu))
    centers <- xnu[sample.int(nrow(xnu), size = kernels), , drop = FALSE]
    sigma <- 10
    score <- -Inf

    ## Searching optimal sigma.
    for (epsilon in seq(from = log10(sigma) - 1, to = -5, by = -1)) {
        for (iteration in seq_len(9)) {
            sigma_new <- sigma - 10 ^ epsilon
            score_new <- 0

            Xnu <- KLIEP_kernel(xnu, centers, sigma_new)
            mean_Xde <- matrix(colMeans(KLIEP_kernel(xde, centers, sigma_new)))

            cv_i <- sample.int(nrow(xnu)) %% fold + 1L
            for (i in seq_len(fold)) {
                alpha_cv <- KLIEP_learning(
                    mean_Xde, Xnu[cv_i != i, , drop = FALSE])
                wh_cv <- Xnu[cv_i == i, , drop = FALSE] %*% alpha_cv
                score_new <- score_new + mean(log(wh_cv)) / fold
            }

            if(score_new <= score)
                break

            sigma <- sigma_new
            score <- score_new
        }
    }

    ## Optimizing kernel weights.
    Xnu <- KLIEP_kernel(xnu, centers, sigma)
    mean_Xde <- matrix(colMeans(KLIEP_kernel(xde, centers, sigma)))
    weights <- KLIEP_learning(mean_Xde, Xnu)

    list(centers = centers, sigma = sigma, weights = weights)
}

KLIEP_density_ratio <- function(x, centers, sigma, weights) {
    as.numeric(KLIEP_kernel(x, centers, sigma) %*% weights)
}
