check_errors <- function(min, q1, med, q3, max, n, scenario){
  if (missing(n)) {
    stop("n must be supplied.")
  }
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || !is.finite(n) ||
      n != as.integer(n) || n < 3 || n > 1e6) {
    stop("n must be a single finite integer between 3 and 1,000,000.")
  }
  if (scenario == 'S1'){
    quants <- c(min, med, max)
    if (!is.numeric(quants) || length(quants) != 3 ||
        any(is.na(quants)) || any(!is.finite(quants))) {
      stop("min, med, and max must be single finite numeric values.")
    }
    if (is.unsorted(quants)){
      stop(paste('The quantiles are not in increasing order. The values of min, med, and max are:', min, med, max))
    }
    if (length(unique(c(quants))) < 3){
      warning(paste('Some of the quantiles are equal to each other. This may lead to unexpected results. The values of min, med, and max are:', min, med, max))
    }
  } else if (scenario == 'S2'){
    quants <- c(q1, med, q3)
    if (!is.numeric(quants) || length(quants) != 3 ||
        any(is.na(quants)) || any(!is.finite(quants))) {
      stop("q1, med, and q3 must be single finite numeric values.")
    }
    if (is.unsorted(quants)){
      stop(paste('The quantiles are not in increasing order. The values of q1, med, and q3 are:', q1, med, q3))
    }
    if (length(unique(quants)) < 3){
      warning(paste('Some of the quantiles are equal to each other. This may lead to unexpected results. The values of q1, med, and q3 are:', q1, med, q3))
    }
  } else if (scenario == 'S3'){
    quants <- c(min, q1, med, q3, max)
    if (!is.numeric(quants) || length(quants) != 5 ||
        any(is.na(quants)) || any(!is.finite(quants))) {
      stop("min, q1, med, q3, and max must be single finite numeric values.")
    }
    if (is.unsorted(quants)){
      stop(paste('The quantiles are not in increasing order. The values of min, q1, med, q3, and max are:', min, q1, med, q3, max))
    }
    if (length(unique(quants)) < 5){
      warning(paste('Some of the quantiles are equal to each other. This may lead to unexpected results. The values of min, q1, med, q3, and max are:', min, q1, med, q3, max))
    }
  }
}
