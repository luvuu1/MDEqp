#' Internal helper functions
#'
#' Utility functions used to identify the scenario, generate starting values
#' and parameter bounds, recover distribution-based mean and standard deviation
#' estimates, and collect supplied inputs.
#'
#' These functions are not intended to be called directly by users.
#'
#' @param min,q1,med,q3,max Reported summary statistics.
#' @param n Sample size.
#' @param scenario Summary-data scenario.
#' @param quants Numeric vector of reported quantiles.
#' @param x Numeric vector.
#' @param family Candidate distribution family.
#'
#' @return Internal objects used in model fitting and post-processing.
#'
#' @name helper_functions
#' @keywords internal
#' @noRd
get.scenario <- function(min, q1, med, q3, max) {
  if (!(missing(min) | missing(q1) | missing(med) | missing(q3)
        | missing(max))) {
    if (!any(is.na(c(min, q1, med, q3, max)))) {
      return("S3")
    }
  }
  if (!(missing(q1) | missing(med) | missing(q3))) {
    if (!any(is.na(c(q1, med, q3)))) {
      return("S2")
    }
  }
  if (!(missing(min) | missing(med) | missing(max))) {
    if (!any(is.na(c(min, med, max)))) {
      return("S1")
    }
  }
  stop("Summary measures not in appropriate form. See documentation for appropriate forms.")
}


#' @rdname helper_functions
#' @keywords internal
#' @noRd
set.qe.fit.control <- function(quants, n, scenario){
  con <- list()
  if (scenario == "S1" | scenario == "S2") {
    con$norm.mu.bounds <- c(quants[1], quants[3])
    med <- quants[2]
    if (min(quants) > 0) {
      con$lnorm.mu.start <- log(med)
      con$lnorm.mu.bounds <- c(log(quants[1]), log(quants[3]))
    }
  }
  if (scenario == "S3") {
    con$norm.mu.bounds <- c(quants[2], quants[4])
    med <- quants[3]
    if (min(quants) > 0) {
      con$lnorm.mu.start <- log(med)
      con$lnorm.mu.bounds <- c(log(quants[2]), log(quants[4]))
    }
  }

  con$norm.sigma.bounds <- c(1e-3, 50)
  con$lnorm.sigma.bounds <- c(1e-3, 10)

  mean.hat <- metaBLUE::Luo.mean(quants, n, scenario)$muhat
  sd.hat <- metaBLUE::Wan.std(quants, n, scenario)$sigmahat

  con$norm.mu.start <- mean.hat
  con$norm.sigma.start <- sd.hat

  con$lnorm.mu.start <- log(mean.hat / sqrt(1 + (sd.hat / mean.hat)^2))
  con$lnorm.sigma.start <- sqrt(log(1 + (sd.hat / mean.hat)^2))

  con$gamma.shape.start <- mean.hat^2 / sd.hat^2
  con$gamma.rate.start <- mean.hat / sd.hat^2

  con$beta.shape1.start <- mean.hat *
    (((mean.hat * (1 - mean.hat)) / (sd.hat^2)) - 1)
  con$beta.shape2.start <- con$beta.shape1.start * (1 - mean.hat) / mean.hat

  start.val <- (mean.hat / sd.hat)^1.086
  obj.fun <- function(k, mean.hat, sd.hat) {
    temp1 <- gamma((k + 2) / k)
    temp2 <- gamma((k + 1) / k)
    temp3 <- sqrt((temp1 / (temp2^2)) - 1)
    ((sd.hat / mean.hat) - temp3)^2
  }
  no.fit <- function(e) {
    return(list(par = NA))
  }
  get.weibull.start <- tryCatch({
    stats::nlminb(start = start.val, objective = obj.fun, lower = 1e-3,
                  mean.hat = mean.hat, sd.hat = sd.hat)
  },
  error = no.fit,
  warning = no.fit
  )
  if (!is.na(get.weibull.start$par)){
    con$weibull.shape.start <- get.weibull.start$par
  } else {
    con$weibull.shape.start <- start.val
  }
  con$weibull.scale.start <- mean.hat / gamma((con$weibull.shape.start + 1)
                                              / con$weibull.shape.start)

  con$gamma.shape.bounds <- c(1e-3, 1e2)
  con$gamma.rate.bounds <- c(1e-3, 1e2)
  con$weibull.shape.bounds <- c(1e-3, 1e2)
  con$weibull.scale.bounds <- c(1e-3, 1e2)
  con$beta.shape1.bounds <- c(1e-3, 40)
  con$beta.shape2.bounds <- c(1e-3, 40)
  return(con)
}


#' @rdname helper_functions
#' @keywords internal
#' @noRd
get.mean.sd <- function(x, family) {
  if (!(family %in% c("normal", "log-normal", "gamma", "weibull", "beta"))) {
    stop("family must be either normal, log-normal, gamma, Weibull, or beta.")
  }

  if (family == "normal") {
    par <- unname(x$norm.par)
    est.mean <- par[1]
    est.sd <- par[2]
  }
  else if (family == "log-normal") {
    par <- unname(x$lnorm.par)
    est.mean <- exp(par[1] + par[2]^2 / 2)
    est.sd <- sqrt((exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2))
  }
  else if (family == "gamma") {
    par <- unname(x$gamma.par)
    est.mean <- par[1] / par[2]
    est.sd <- sqrt(par[1] / (par[2]^2))
  }
  else if (family == "weibull") {
    par <- unname(x$weibull.par)
    est.mean <- par[2] * gamma(1 + 1 / par[1])
    est.sd <- sqrt(par[2]^2 * (gamma(1 + 2 / par[1]) -
                                 (gamma(1 + 1 / par[1]))^2))
  }
  else if (family == "beta") {
    par <- unname(x$beta.par)
    est.mean <- par[1] / (par[1] + par[2])
    est.sd <- sqrt(par[1] * par[2] / ((par[1] + par[2])^2 *
                                        (par[1] + par[2] + 1)))
  }

  ## When zero-based range normalization is used, fitted candidate distributions
  ## are evaluated on the normalized data. Therefore distribution-specific mean
  ## and SD estimates must be restored to the original scale, not only beta.
  if (!is.null(x$scaled) && isTRUE(x$scaled)) {
    est.mean <- restore.mean.from.range(est.mean, x$range)
    est.sd <- restore.sd.from.range(est.sd, x$range)
  }

  return(list(est.mean = est.mean, est.sd = est.sd))
}


#' @rdname helper_functions
#' @keywords internal
#' @noRd
get.num.input <- function(min, q1, med, q3, max, n, range){
  res <- list()
  if (!missing(min)){
    res$min <- min
  }
  if (!missing(q1)){
    res$q1 <- q1
  }
  if (!missing(med)){
    res$med <- med
  }
  if (!missing(q3)){
    res$q3 <- q3
  }
  if (!missing(max)){
    res$max <- max
  }
  if (!missing(n)){
    res$n <- n
  }
  if (!missing(range)){
    res$range <- range
  }
  return(res)
}


#' Get theoretical range for zero-based normalization
#' @rdname helper_functions
#' @keywords internal
#' @noRd
get.estimate.bounds <- function(x) {
  if (!is.null(x$range)) {
    return(c(0, x$range))
  }

  input <- x$num.input

  if (!is.null(x$scenario) &&
      x$scenario %in% c("S1", "S3") &&
      !is.null(input$min) &&
      !is.null(input$max)) {
    return(c(input$min, input$max))
  }

  NULL
}

#' Check whether distribution-specific mean and SD estimates are implausible
#' @rdname helper_functions
#' @keywords internal
#' @noRd
is.implausible.estimate <- function(est.mean, est.sd, x) {
  if (length(est.mean) != 1 || length(est.sd) != 1) {
    return(TRUE)
  }

  if (is.na(est.mean) || is.na(est.sd) ||
      !is.finite(est.mean) || !is.finite(est.sd)) {
    return(TRUE)
  }

  if (est.sd < 0) {
    return(TRUE)
  }

  bounds <- get.estimate.bounds(x)

  if (is.null(bounds)) {
    return(FALSE)
  }

  lower <- bounds[1]
  upper <- bounds[2]

  if (upper < lower) {
    return(TRUE)
  }

  if (est.mean < lower || est.mean > upper) {
    return(TRUE)
  }

  FALSE
}


#' Review invalid fit statistics and implausible estimates
#' @rdname helper_functions
#' @keywords internal
#' @noRd
check.fit.estimates <- function(x) {
  all.ests <- list()

  for (dist.name in names(x$valuesls)) {
    ests <- get.mean.sd(x, dist.name)

    bad.fit <- is.na(x$valuesls[[dist.name]]) ||
      !is.finite(x$valuesls[[dist.name]])
    bad.est <- is.implausible.estimate(ests$est.mean, ests$est.sd, x)

    if (bad.fit || bad.est) {
      x$valuesls[dist.name] <- Inf
      ests <- list(est.mean = NA_real_, est.sd = NA_real_)
    }

    all.ests[[dist.name]] <- ests
  }

  list(fit = x, estimates = all.ests)
}
