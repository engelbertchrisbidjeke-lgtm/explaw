#' Total exposure time
#'
#' Computes the total time units have been under observation. Because
#' censored units contribute their observed lifetime (the time until they
#' left the study), simply summing the `time` column gives the total
#' exposure, which is the denominator of the maximum-likelihood estimator
#' of the exponential rate.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A single positive `numeric`, in the same time unit as `d$time`
#'   (e.g. hours).
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' exposure_time(d)
#'
#' @seealso [estimate_lambda()], [confidence_interval()]
#' @export
exposure_time <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  sum(d$time)
}

#' Estimate the exponential failure rate
#'
#' Computes the maximum-likelihood estimate of the failure rate \eqn{\lambda}
#' under an exponential model: \eqn{\hat{\lambda} = D / E}, where \eqn{D}
#' is the number of observed failures and \eqn{E} is the total exposure
#' time (see [exposure_time()]).
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A `list` with three elements:
#'   \describe{
#'     \item{`lambda`}{Numeric. Estimated failure rate (failures per unit
#'       of time).}
#'     \item{`n_events`}{Integer. Number of observed failures.}
#'     \item{`exposure`}{Numeric. Total exposure time.}
#'   }
#'   If no failure is observed, `lambda` is `0` and a warning is issued.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' estimate_lambda(d)
#'
#' @seealso [exposure_time()], [confidence_interval()], [indicators_table()]
#' @export
estimate_lambda <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  n_events <- sum(d$event == 1L)
  exposure <- sum(d$time)
  if (exposure <= 0) {
    stop("Total exposure time is zero; cannot estimate lambda.")
  }
  lambda <- n_events / exposure
  if (n_events == 0) {
    warning("No events observed: lambda estimated at 0. ",
            "Only a lower bound on lifetime can be derived.")
  }
  list(
    lambda   = lambda,
    n_events = n_events,
    exposure = exposure
  )
}

#' Confidence interval for the exponential rate and mean lifetime
#'
#' Builds an exact confidence interval for \eqn{\lambda} under an
#' exponential model, using the pivotal quantity
#' \eqn{2 D \lambda \sim \chi^2_{2D}}, where \eqn{D} is the number of
#' observed failures. The interval on the mean lifetime \eqn{1/\lambda}
#' is derived by inverting the bounds.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#' @param level Numeric. Confidence level, strictly between 0 and 1.
#'   Default is `0.95`.
#'
#' @return A one-row `data.frame` with four columns:
#'   \describe{
#'     \item{`lambda_lower`, `lambda_upper`}{Bounds of the confidence
#'       interval for the failure rate.}
#'     \item{`mean_lower`, `mean_upper`}{Bounds of the confidence
#'       interval for the mean lifetime, in the same time unit as `d$time`.}
#'   }
#'   If no failure is observed, all bounds are returned as `0` or `Inf`
#'   and a warning is issued.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' confidence_interval(d)
#' confidence_interval(d, level = 0.90)
#'
#' @seealso [estimate_lambda()]
#' @export
confidence_interval <- function(d, level = 0.95) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  if (!is.numeric(level) || length(level) != 1 ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.")
  }
  n_events <- sum(d$event == 1L)
  exposure <- sum(d$time)
  if (n_events == 0) {
    warning("No events observed: confidence interval on lambda ",
            "is [0, +Inf].")
    return(data.frame(
      lambda_lower = 0,
      lambda_upper = Inf,
      mean_lower   = 0,
      mean_upper   = Inf
    ))
  }
  alpha <- (1 - level) / 2
  lambda_lower <- stats::qchisq(alpha,     2 * n_events) / (2 * exposure)
  lambda_upper <- stats::qchisq(1 - alpha, 2 * n_events) / (2 * exposure)
  data.frame(
    lambda_lower = lambda_lower,
    lambda_upper = lambda_upper,
    mean_lower   = 1 / lambda_upper,
    mean_upper   = 1 / lambda_lower
  )
}
