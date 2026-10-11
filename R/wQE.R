#' Weighted quantile estimation method for estimating the sample mean and standard deviation
#'
#' This function applies the weighted quantile estimation (wQE) method to
#' estimate the sample mean and standard deviation from a study reporting one
#' of the following sets of summary statistics:
#' \itemize{
#'   \item \strong{S1}: minimum, median, maximum, and sample size;
#'   \item \strong{S2}: first quartile, median, third quartile, and sample size;
#'   \item \strong{S3}: minimum, first quartile, median, third quartile,
#'   maximum, and sample size.
#' }
#'
#' The wQE method fits candidate parametric distributions to the reported
#' quantiles using \code{\link{wQE.fit}}. Candidate distributions are fitted
#' using the method-specific objective. The selected candidate distribution is
#' then chosen using the smallest unweighted quantile-scale squared discrepancy
#' between the reported quantiles and fitted quantiles.
#'
#' Because different reporting scenarios require different summary statistics,
#' users should supply arguments by name.
#'
#' Candidate distributions include the normal, log-normal, gamma, Weibull, and
#' beta families. When \code{range} is supplied, the current implementation
#' assumes a theoretical lower bound of 0, so \code{range} represents the upper
#' endpoint of the theoretical scale and all reported summary statistics must
#' fall within \code{[0, range]}. All candidate distributions are fitted on the
#' normalized scale, and the final mean and SD estimates are transformed back to
#' the original scale.
#'
#' @param min Numeric value giving the sample minimum.
#' @param q1 Numeric value giving the sample first quartile.
#' @param med Numeric value giving the sample median.
#' @param q3 Numeric value giving the sample third quartile.
#' @param max Numeric value giving the sample maximum.
#' @param n Numeric value giving the sample size.
#' @param range Optional numeric scalar specifying the theoretical upper
#'   endpoint of a zero-based bounded outcome. The current implementation
#'   assumes a lower bound of 0. When supplied, all reported summary statistics
#'   must fall within \code{[0, range]}. The summaries are divided by
#'   \code{range} before fitting, and the estimated mean and SD are transformed
#'   back to the original scale. This value should reflect the theoretical
#'   scale of the outcome, not the observed sample maximum or sample range.
#'
#' @return
#' An object of class \code{"wQE"}, which is a list with the following
#' components:
#' \describe{
#'   \item{\code{est.mean}}{Estimated sample mean.}
#'   \item{\code{est.sd}}{Estimated sample standard deviation.}
#'   \item{\code{selected.dist}}{Selected outcome distribution.}
#'   \item{\code{valuesls}}{Values of the selection statistic used to compare
#'   candidate distributions.}
#'   \item{\code{args}}{List of arguments supplied to the function.}
#'   \item{\code{scenario}}{Detected summary-data scenario: \code{"S1"},
#'   \code{"S2"}, or \code{"S3"}.}
#'   \item{\code{fitted.dists}}{Object returned by \code{\link{wQE.fit}}
#'   containing distribution-specific fitting results.}
#'   \item{\code{normal.est.mean}, \code{log-normal.est.mean},
#'   \code{gamma.est.mean}, \code{weibull.est.mean},
#'   \code{beta.est.mean}}{Estimated sample means.}
#'   \item{\code{normal.est.sd}, \code{log-normal.est.sd},
#'   \code{gamma.est.sd}, \code{weibull.est.sd},
#'   \code{beta.est.sd}}{Estimated sample standard deviations.}
#' }
#'
#' @seealso \code{\link{wQE.fit}}
#'
#' @examples
#' ## Example 1: full five-number summary (scenario S3)
#' wQE(
#'   min = 1,
#'   q1 = 3,
#'   med = 5,
#'   q3 = 7,
#'   max = 10,
#'   n = 50
#' )
#'
#' ## Example 2: zero-based bounded outcome with known theoretical range
#' wQE(
#'   min = 5,
#'   q1 = 20,
#'   med = 35,
#'   q3 = 50,
#'   max = 80,
#'   n = 100,
#'   range = 100
#' )
#'
#' @export

#weighted QE
wQE <- function(min, q1, med, q3, max, n, range) {
  args <- list()
  if (!missing(min)) args$min <- min
  if (!missing(q1)) args$q1 <- q1
  if (!missing(med)) args$med <- med
  if (!missing(q3)) args$q3 <- q3
  if (!missing(max)) args$max <- max
  if (!missing(n)) args$n <- n
  if (!missing(range)) args$range <- range

  x <- do.call(wQE.fit, args)

  checked <- check.fit.estimates(x)
  x <- checked$fit
  all.ests <- checked$estimates

  if (all(!is.finite(x$valuesls))) {
    selected.dist <- NA_character_
    ests <- list(est.mean = NA_real_, est.sd = NA_real_)
  } else {
    selected.dist <- names(which.min(x$valuesls))
    ests <- all.ests[[selected.dist]]
  }


  output <- list(
    est.mean = ests$est.mean,
    est.sd = ests$est.sd,
    selected.dist = selected.dist,
    valuesls = x$valuesls,
    args = args,
    scenario = x$scenario,
    fitted.dists = x
  )

  for (dist.name in names(x$valuesls)) {
    ests.all <- all.ests[[dist.name]]
    output[[paste0(dist.name, ".est.mean")]] <- ests.all$est.mean
    output[[paste0(dist.name, ".est.sd")]] <- ests.all$est.sd
  }

  class(output) <- "wQE"
  return(output)
}

#' Print method for objects of class \code{"wQE"}
#'
#' Displays the main estimation results from an object of class
#' \code{"wQE"}, including the detected scenario, selected distribution,
#' estimated mean, and estimated standard deviation.
#'
#' @param x An object of class \code{"wQE"}.
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The input object, invisibly.
#'
#' @seealso \code{\link{wQE}}
#'
#' @examples
#' x <- list(
#'   scenario = "S2",
#'   selected.dist = "normal",
#'   est.mean = 10.12345,
#'   est.sd = 2.34567
#' )
#' class(x) <- "wQE"
#' print(x)
#'
#' @method print wQE
#' @export


## S3 method for class 'wQE'
print.wQE <- function(x, digits = 4, ...) {
  cat("wQE mean and SD estimation\n")
  cat("------------------------\n")
  cat("Scenario:", x$scenario, "\n")
  cat("Selected distribution:", x$selected.dist, "\n\n")

  cat("Estimated mean: ", round(x$est.mean, digits), "\n")
  cat("Estimated SD: ", round(x$est.sd, digits), "\n")

  invisible(x)
}
