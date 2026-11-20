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

# Internal: checks that `d` is a usable prepared data frame.
# Stops with an explicit message if `d` lacks the `time` / `event` columns,
# if `time` is not numeric, finite and non-negative, or if `event` contains
# anything other than 0 and 1.
#' @noRd
.validate_d <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  if (!is.numeric(d$time) || any(!is.finite(d$time)) || any(d$time < 0)) {
    stop("'d$time' must be numeric, finite and non-negative. ",
         "Use prepare_data() first.")
  }
  if (!all(d$event %in% c(0, 1))) {
    stop("'d$event' must contain only 0 and 1 (no NA). ",
         "Use prepare_data() first.")
  }
  invisible(d)
}

# Internal: maximum-likelihood fit of the exponential rate, WITHOUT any
# warning. Public functions decide themselves whether and how to warn when
# no event is observed, so that a single call emits a single warning.
# Returns list(lambda, n_events, exposure).
#' @noRd
.fit_lambda <- function(d) {
  .validate_d(d)
  n_events <- sum(d$event == 1L)
  exposure <- sum(d$time)
  if (exposure <= 0) {
    stop("Total exposure time is zero; cannot estimate lambda.")
  }
  list(
    lambda   = n_events / exposure,
    n_events = n_events,
    exposure = exposure
  )
}

#' Estimate the exponential failure rate
#'
#' Computes the maximum-likelihood estimate of the failure rate \eqn{\lambda}
#' under an exponential model: \eqn{\hat{\lambda} = D / E}, where \eqn{D}
#' is the number of observed failures and \eqn{E} is the total exposure
#' time (see [exposure_time()]).
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()]. It
#'   must contain a `time` column (numeric, finite, non-negative) and an
#'   `event` column (only `0` and `1`, no `NA`). Otherwise the function
#'   stops with an error.
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
#' @details
#' The function stops with an error when `d` is not a valid prepared data
#' frame (missing columns, negative or non-finite times, `event` outside
#' `{0, 1}`) or when the total exposure time is zero, since \eqn{\lambda}
#' cannot be estimated in that case.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' estimate_lambda(d)
#'
#' @seealso [exposure_time()], [confidence_interval()], [indicators_table()]
#' @export
estimate_lambda <- function(d) {
  est <- .fit_lambda(d)
  if (est$n_events == 0) {
    warning("No events observed: lambda estimated at 0. ",
            "Only a lower bound on lifetime can be derived.")
  }
  est
}

#' Confidence interval for the exponential rate and mean lifetime
#'
#' Builds an exact confidence interval for \eqn{\lambda} under an
#' exponential model, using the pivotal quantity
#' \eqn{2 D \lambda \sim \chi^2_{2D}}, where \eqn{D} is the number of
#' observed failures. The interval on the mean lifetime \eqn{1/\lambda}
#' is derived by inverting the bounds.
#'
#' When no failure is observed (\eqn{D = 0}), the lower bound of
#' \eqn{\lambda} is `0` and the upper bound is the finite value
#' \eqn{-\ln(\alpha)/E}, where \eqn{\alpha = (1 - \text{level}) / 2} and
#' \eqn{E} is the total exposure time. Accordingly, the mean lifetime
#' interval is \eqn{[1/\lambda_{upper}, +\infty)}.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()]. It
#'   must contain a `time` column (numeric, finite, non-negative) and an
#'   `event` column (only `0` and `1`, no `NA`). Otherwise the function
#'   stops with an error.
#' @param level Numeric. Confidence level, a single number strictly
#'   between 0 and 1 (not `NA`). Default is `0.95`.
#'
#' @return A one-row `data.frame` with four columns:
#'   \describe{
#'     \item{`lambda_lower`, `lambda_upper`}{Bounds of the confidence
#'       interval for the failure rate.}
#'     \item{`mean_lower`, `mean_upper`}{Bounds of the confidence
#'       interval for the mean lifetime, in the same time unit as `d$time`.}
#'   }
#'   If no failure is observed, `lambda_lower` is `0`, `lambda_upper` is
#'   finite, `mean_upper` is `Inf`, and a warning is issued.
#'
#' @details
#' The function stops with an error when `level` is invalid, when `d` is
#' not a valid prepared data frame, or when the total exposure time is
#' zero (no finite bound can be derived in that case).
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
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.")
  }
  est      <- .fit_lambda(d)
  alpha    <- (1 - level) / 2
  exposure <- est$exposure
  n_events <- est$n_events

  if (n_events == 0) {
    warning("No events observed: the confidence interval on lambda ",
            "has a lower bound of 0 and an upper bound of -ln(alpha)/E.")
    lambda_upper <- -log(alpha) / exposure
    return(data.frame(
      lambda_lower = 0,
      lambda_upper = lambda_upper,
      mean_lower   = 1 / lambda_upper,
      mean_upper   = Inf
    ))
  }

  lambda_lower <- stats::qchisq(alpha,     2 * n_events) / (2 * exposure)
  lambda_upper <- stats::qchisq(1 - alpha, 2 * n_events) / (2 * exposure)
  data.frame(
    lambda_lower = lambda_lower,
    lambda_upper = lambda_upper,
    mean_lower   = 1 / lambda_upper,
    mean_upper   = 1 / lambda_lower
  )
}
