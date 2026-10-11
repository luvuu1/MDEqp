#' Minimum Distance Estimation Fitting
#'
#' Fits candidate parametric distributions to reported quantiles using the
#' minimum distance estimation (MDE) method.
#'
#' This lower-level function is mainly intended for inspecting
#' distribution-specific fitting results. Most users should call
#' \code{\link{MDE}} directly.
#'
#' Three summary-data scenarios are supported:
#' \itemize{
#'   \item \strong{S1}: minimum, median, maximum, and sample size;
#'   \item \strong{S2}: first quartile, median, third quartile, and sample size;
#'   \item \strong{S3}: minimum, first quartile, median, third quartile,
#'   maximum, and sample size.
#' }
#'
#' @param min Numeric value giving the sample minimum. Required for scenarios
#'   \code{"S1"} and \code{"S3"}.
#' @param q1 Numeric value giving the sample first quartile (25th percentile).
#'   Required for scenarios \code{"S2"} and \code{"S3"}.
#' @param med Numeric value giving the sample median (50th percentile).
#'   Required in all scenarios.
#' @param q3 Numeric value giving the sample third quartile (75th percentile).
#'   Required for scenarios \code{"S2"} and \code{"S3"}.
#' @param max Numeric value giving the sample maximum. Required for scenarios
#'   \code{"S1"} and \code{"S3"}.
#' @param n Numeric value giving the sample size.
#' @param range Optional numeric scalar specifying the theoretical upper
#'   endpoint of a zero-based bounded outcome. The current implementation
#'   assumes a lower bound of 0. When supplied, all reported summary statistics
#'   must fall within \code{[0, range]}. The summaries are divided by
#'   \code{range} before fitting, and the estimated mean and SD are transformed
#'   back to the original scale. This value should reflect the theoretical
#'   scale of the outcome, not the observed sample maximum or sample range.
#'
#' @details
#' Candidate distributions include the normal, log-normal, gamma, Weibull, and
#' beta families.
#'
#' Optimization is performed using \code{"L-BFGS-B"} in
#' \code{\link[stats]{optim}} with internally generated starting values and
#' parameter bounds.
#'
#' If \code{range} is supplied, zero-based range normalization is applied after
#' the summary-data scenario has been detected and before candidate
#' distributions are evaluated. This means S2 calls do not need \code{min} or
#' \code{max}, but fitted distributions are compared on the same normalized
#' outcome scale.
#'
#' @return
#' An object of class \code{"MDE.fit"}, which is a list containing raw
#' distribution-specific fitting results. The list includes the following
#' components:
#' \describe{
#'   \item{\code{norm.par}}{Estimated parameters for the normal distribution.}
#'   \item{\code{lnorm.par}}{Estimated parameters for the log-normal distribution.}
#'   \item{\code{gamma.par}}{Estimated parameters for the gamma distribution.}
#'   \item{\code{weibull.par}}{Estimated parameters for the Weibull distribution.}
#'   \item{\code{beta.par}}{Estimated parameters for the beta distribution.}
#'   \item{\code{valuesls}}{Values of the selection statistic used to compare
#'   candidate distributions.}
#'   \item{\code{num.input}}{Input data summary statistics used in fitting.}
#'   \item{\code{scenario}}{Detected summary-data scenario: \code{"S1"},
#'   \code{"S2"}, or \code{"S3"}.}
#'   \item{\code{scaled}}{Whether zero-based range normalization was used.}
#'   \item{\code{range}}{Theoretical range used for zero-based normalization, or \code{NULL} if not provided.}
#' }
#'
#' @seealso \code{\link{MDE}}
#'
#' @examples
#' ## Example 1: full five-number summary (scenario S3)
#' fit1 <- MDE.fit(
#'   min = 1,
#'   q1 = 3,
#'   med = 5,
#'   q3 = 7,
#'   max = 10,
#'   n = 50
#' )
#'
#' ## Example 2: zero-based bounded outcome with known theoretical range
#' fit2 <- MDE.fit(
#'   q1 = 3,
#'   med = 5,
#'   q3 = 7,
#'   n = 100,
#'   range = 27
#' )
#'
#' @export
#' @importFrom stats pnorm plnorm pweibull pbeta pgamma
#' @importFrom stats dnorm dlnorm dweibull dbeta dgamma
#' @importFrom stats qnorm qlnorm qweibull qbeta qgamma optim



#MDE
MDE.fit <- function(min, q1, med, q3, max, n, range){
  scenario <- get.scenario(min, q1, med, q3, max)
  check_errors(min, q1, med, q3, max, n, scenario)

  if (scenario == "S1") {
    probs <- c(0.625 / n, 0.5, 1 - 0.625 / n)
    quants <- c(min, med, max)
  } else if (scenario == "S2") {
    probs <- c(0.25, 0.5, 0.75)
    quants <- c(q1, med, q3)
  } else if (scenario == "S3") {
    probs <- c(0.625 / n, 0.25, 0.5, 0.75, 1 - 0.625 / n)
    quants <- c(min, q1, med, q3, max)
  }

  original.quants <- quants
  scaled <- FALSE
  scale.range <- NULL

  ## Zero-based normalization must happen after scenario detection. This keeps S2 valid
  ## without min/max, while still evaluating candidate distributions on the
  ## same normalized outcome scale when range is supplied.
  if (!missing(range)) {
    validate.range(original.quants, range)
    quants <- scale_by_range(original.quants, range)
    scaled <- TRUE
    scale.range <- range
  }

  con <- set.qe.fit.control(quants, n, scenario)
  no.fit <- function(e) {
    return(list(par = NA, value = NA))
  }
  get.qfit.value <- function(pre) {
    # Failed fits are kept at Inf so that they cannot be selected.
    if (any(!is.finite(pre))) {
      return(Inf)
    }
    sum((quants - pre)^2)
  }

  if (scaled) {
    beta.start <- get.beta.start(quants, n, scenario)
    con$beta.shape1.start <- beta.start$shape1
    con$beta.shape2.start <- beta.start$shape2
  }
  beta.applicable <- is.beta.applicable(quants)

  fit.norm <- tryCatch({
    N <- n
    start_params <- c(mu = con$norm.mu.start, sigma = con$norm.sigma.start)
    p <- probs
    q <- quants

    objective_function_w <- function(params, q, p, N) {
      mu <- params[1]
      sigma <- params[2]

      if (sigma <= 0) return(Inf)

      q_pred <- qnorm(p, mean = mu, sd = sigma)
      x <- qnorm(p, mu, sigma)
      f <- dnorm(x, mean = mu, sd = sigma)

      cov <- matrix(NA, nrow = length(p), ncol = length(p))
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i, j] <- p[i] * (1 - p[j]) / (f[i] * f[j])
          cov[j, i] <- cov[i, j]
        }
      }
      weights <- qr.solve(cov)
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }

    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$norm.mu.bounds[1], con$norm.sigma.bounds[1]),
      upper = c(con$norm.mu.bounds[2], con$norm.sigma.bounds[2]),
      hessian = FALSE)
    result_w
  }, error = no.fit)

  pre.n <- qnorm(probs, mean = as.numeric(fit.norm$par[1]), sd = as.numeric(fit.norm$par[2]))
  value.n <- get.qfit.value(pre.n)

  fit.lnorm <- tryCatch({
    N <- n
    start_params <- c(mu = con$lnorm.mu.start, sigma = con$lnorm.sigma.start)
    p <- probs
    q <- quants

    objective_function_w <- function(params, q, p, N) {
      mu <- params[1]
      sigma <- params[2]

      q_pred <- qlnorm(p, meanlog = mu, sdlog = sigma)
      x <- qlnorm(p, meanlog = mu, sdlog = sigma)
      f <- dlnorm(x, meanlog = mu, sdlog = sigma)

      cov <- matrix(NA, nrow = length(p), ncol = length(p))
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i, j] <- p[i] * (1 - p[j]) / (f[i] * f[j])
          cov[j, i] <- cov[i, j]
        }
      }
      weights <- qr.solve(cov)
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }

    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$lnorm.mu.bounds[1], con$lnorm.sigma.bounds[1]),
      upper = c(con$lnorm.mu.bounds[2], con$lnorm.sigma.bounds[2]),
      hessian = FALSE)
    result_w
  }, error = no.fit)

  pre.l <- qlnorm(probs, meanlog = as.numeric(fit.lnorm$par[1]), sdlog = as.numeric(fit.lnorm$par[2]))
  value.l <- get.qfit.value(pre.l)

  fit.weibull <- tryCatch({
    N <- n
    start_params <- c(bsa = con$weibull.shape.start, bsb = con$weibull.scale.start)
    p <- probs
    q <- quants

    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]

      q_pred <- qweibull(p, shape = bsa, scale = bsb)
      x <- qweibull(p, shape = bsa, scale = bsb)
      f <- dweibull(x, shape = bsa, scale = bsb)

      cov <- matrix(NA, nrow = length(p), ncol = length(p))
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i, j] <- p[i] * (1 - p[j]) / (f[i] * f[j])
          cov[j, i] <- cov[i, j]
        }
      }
      weights <- qr.solve(cov)
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }

    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$weibull.shape.bounds[1], con$weibull.scale.bounds[1]),
      upper = c(con$weibull.shape.bounds[2], con$weibull.scale.bounds[2]),
      hessian = FALSE)
    result_w
  }, error = no.fit)

  pre.w <- qweibull(probs, shape = as.numeric(fit.weibull$par[1]), scale = as.numeric(fit.weibull$par[2]))
  value.w <- get.qfit.value(pre.w)

  if (!beta.applicable) {
    fit.beta <- list(par = NA_real_, value = Inf)
    value.b <- Inf
  } else {
    fit.beta <- tryCatch({
      N <- n
      start_params <- c(con$beta.shape1.start, con$beta.shape2.start)
      p <- probs
      q <- quants

      objective_function_w <- function(params, q, p, N) {
        bsa <- params[1]
        bsb <- params[2]

        q_pred <- qbeta(p, shape1 = bsa, shape2 = bsb)
        x <- qbeta(p, shape1 = bsa, shape2 = bsb)
        f <- dbeta(x, shape1 = bsa, shape2 = bsb)

        cov <- matrix(NA, nrow = length(p), ncol = length(p))
        for (i in 1:length(p)){
          for (j in i:length(p)){
            cov[i, j] <- p[i] * (1 - p[j]) / (f[i] * f[j])
            cov[j, i] <- cov[i, j]
          }
        }
        weights <- qr.solve(cov)
        residuals <- q - q_pred
        weighted_residuals <- t(residuals) %*% weights %*% residuals
        weighted_residuals
      }

      result_w <- optim(
        par = start_params,
        fn = objective_function_w,
        q = q,
        p = p,
        N = N,
        method = "L-BFGS-B",
        lower = c(con$beta.shape1.bounds[1], con$beta.shape2.bounds[1]),
        upper = c(con$beta.shape1.bounds[2], con$beta.shape2.bounds[2]),
        hessian = FALSE)
      result_w
    }, error = no.fit)

    pre.b <- qbeta(probs, shape1 = as.numeric(fit.beta$par[1]), shape2 = as.numeric(fit.beta$par[2]))
    value.b <- get.qfit.value(pre.b)
  }

  fit.gamma <- tryCatch({
    N <- n
    start_params <- c(con$gamma.shape.start, con$gamma.rate.start)
    p <- probs
    q <- quants

    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]

      q_pred <- qgamma(p, shape = bsa, rate = bsb)
      x <- qgamma(p, shape = bsa, rate = bsb)
      f <- dgamma(x, shape = bsa, rate = bsb)

      cov <- matrix(NA, nrow = length(p), ncol = length(p))
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i, j] <- p[i] * (1 - p[j]) / (f[i] * f[j])
          cov[j, i] <- cov[i, j]
        }
      }
      weights <- qr.solve(cov)
      residuals <- q - q_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }

    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$gamma.shape.bounds[1], con$gamma.rate.bounds[1]),
      upper = c(con$gamma.shape.bounds[2], con$gamma.rate.bounds[2]),
      hessian = FALSE)
    result_w
  }, error = no.fit)

  pre.g <- qgamma(probs, shape = as.numeric(fit.gamma$par[1]), rate = as.numeric(fit.gamma$par[2]))
  value.g <- get.qfit.value(pre.g)

  valuesls <- c(value.n, value.l, value.g, value.w, value.b)
  names(valuesls) <- c("normal", "log-normal", "gamma", "weibull", "beta")

  norm.par <- fit.norm$par
  lnorm.par <- fit.lnorm$par
  gamma.par <- fit.gamma$par
  weibull.par <- fit.weibull$par
  beta.par <- fit.beta$par

  num.input <- get.num.input(min, q1, med, q3, max, n, range)

  output <- list(
    norm.par = norm.par,
    lnorm.par = lnorm.par,
    gamma.par = gamma.par,
    weibull.par = weibull.par,
    beta.par = beta.par,
    valuesls = valuesls,
    num.input = num.input,
    scenario = scenario,
    original.quants = original.quants,
    fitted.quants = quants,
    scaled = scaled,
    range = scale.range
  )

  class(output) <- "MDE.fit"
  return(output)
}
