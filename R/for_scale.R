#' Internal validation for zero-based range normalization
#'
#' @param quants Reported quantiles used in the detected scenario.
#' @param range Theoretical range of a zero-based bounded outcome. The lower
#'   bound is assumed to be 0.
#' @keywords internal
#' @noRd

validate.range <- function(quants, range) {
  if (missing(range)) {
    stop("range must be supplied when zero-based range normalization is requested.")
  }
  if (!is.numeric(range)) {
    stop("range must be numeric.")
  }
  if (length(range) != 1) {
    stop("range must be a single number.")
  }
  if (is.na(range) || !is.finite(range)) {
    stop("range cannot be NA, NaN, or infinite.")
  }
  if (range <= 0) {
    stop("range must be positive.")
  }
  if (any(is.na(quants)) || any(!is.finite(quants))) {
    stop("Summary statistics used for fitting must be finite numeric values.")
  }
  if (any(quants < 0 | quants > range)) {
    stop("Summary statistics must fall within [0, range] when range is supplied.")
  }
}


#' Internal zero-based normalization by theoretical range
#'
#' Normalizes reported summary statistics by assuming the theoretical lower bound is 0.
#' @param x Numeric vector.
#' @param range Theoretical range of a zero-based bounded outcome.
#' @return Numeric vector divided by \code{range}.
#' @keywords internal
#' @noRd

scale_by_range <- function(x, range) {
  x / range
}

#' Internal check for beta-distribution applicability
#'
#' The beta distribution is only fitted when the quantiles used for fitting
#' are on the [0, 1] scale. If zero-based range normalization is used, this check is applied
#' after normalization.
#' @param beta.quants Numeric vector of quantiles used for beta fitting.
#' @return Logical value indicating whether beta fitting is applicable.
#' @keywords internal
#' @noRd

is.beta.applicable <- function(beta.quants) {
  all(is.finite(beta.quants)) &&
    all(beta.quants >= 0 & beta.quants <= 1)
}



#' Internal mean back-transformation after zero-based range normalization
#'
#' @param x Mean estimate on the scaled scale.
#' @param range Theoretical range of a zero-based bounded outcome.
#' @return Mean estimate on the original scale.
#' @keywords internal
#' @noRd

restore.mean.from.range <- function(x, range) {
  x * range
}


#' Internal standard deviation back-transformation after zero-based range normalization
#'
#' @param x Standard deviation estimate on the scaled scale.
#' @param range Theoretical range of a zero-based bounded outcome.
#' @return Standard deviation estimate on the original scale.
#' @keywords internal
#' @noRd

restore.sd.from.range <- function(x, range) {
  x * range
}


#' Internal starting values for beta-distribution fitting
#'
#' @param beta.quants Reported quantiles scaled to approximately \code{[0, 1]}.
#' @param n Numeric value giving the sample size.
#' @param scenario Summary-data scenario.
#' @return A list with elements \code{shape1} and \code{shape2}.
#' @keywords internal
#' @noRd

get.beta.start <- function(beta.quants, n, scenario) {
  eps <- 1e-6

  mean.hat <- metaBLUE::Luo.mean(beta.quants, n, scenario)$muhat
  sd.hat <- metaBLUE::Wan.std(beta.quants, n, scenario)$sigmahat

  mean.hat <- min(max(mean.hat, eps), 1 - eps)

  max.var <- mean.hat * (1 - mean.hat)
  max.sd <- sqrt(max(max.var - eps, eps))

  if (!is.finite(sd.hat) || sd.hat <= 0) {
    sd.hat <- max.sd / 2
  }
  if (sd.hat >= max.sd) {
    sd.hat <- 0.95 * max.sd
  }

  common.term <- mean.hat * (1 - mean.hat) / (sd.hat^2) - 1
  shape1 <- max(mean.hat * common.term, 1e-3)
  shape2 <- max((1 - mean.hat) * common.term, 1e-3)

  list(shape1 = shape1, shape2 = shape2)
}
