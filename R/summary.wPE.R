#' Summary method for objects of class \code{"wPE"}
#'
#' @param object An object of class \code{"wPE"}.
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{"summary.wPE"} and \code{"data.frame"}.
#' Rows correspond to candidate distributions and columns contain the
#' estimated mean and standard deviation. The result is ordered from best
#' to worst fit and includes attributes \code{scenario},
#' \code{selected.dist}, and \code{args}.
#'
#' @seealso \code{\link{wPE}}
#'
#' @examples
#' res <- wPE(min = 1, q1 = 3, med = 5, q3 = 7, max = 10, n = 50)
#' summary(res)
#'
#' @method summary wPE
#' @export
summary.wPE <- function(object, digits = 3, ...) {
  dists <- names(object$valuesls)
  objective <- as.numeric(object$valuesls)

  out <- data.frame(
    distribution = dists,
    est.mean = NA_real_,
    est.sd = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(dists)) {
    d <- dists[i]
    out$est.mean[i] <- object[[paste0(d, ".est.mean")]]
    out$est.sd[i]   <- object[[paste0(d, ".est.sd")]]
  }

  out <- out[order(objective), ]
  rownames(out) <- out$distribution
  out$distribution <- NULL

  out$est.mean <- round(out$est.mean, digits)
  out$est.sd <- round(out$est.sd, digits)

  colnames(out) <- c("Mean", "SD")

  attr(out, "selected.dist") <- object$selected.dist
  attr(out, "scenario") <- object$scenario
  attr(out, "args") <- object$args

  class(out) <- c("summary.wPE", "data.frame")

  out
}

#' @export
print.summary.wPE <- function(x, ...) {
  cat("Summary of wPE method:\n")
  cat("Scenario:", attr(x, "scenario"), "\n")
  cat("Selected distribution:", attr(x, "selected.dist"), "\n\n")

  print.data.frame(x)

  invisible(x)
}
