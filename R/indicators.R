#' Mean lifetime under the exponential model
#'
#' Returns the expected lifetime \eqn{E(T) = 1/\lambda} under an exponential
#' model with rate estimated from the data.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A single positive `numeric` (mean lifetime, in the same time
#'   unit as `d$time`). Returns `Inf` with a single warning if no event was
#'   observed.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' mean_lifetime(d)
#'
#' @seealso [median_lifetime()], [quantile_lifetime()], [indicators_table()]
#' @export
mean_lifetime <- function(d) {
  .validate_d(d)
  est <- suppressWarnings(estimate_lambda(d))
  if (est$n_events == 0) {
    warning("lambda is 0 (no events observed): mean lifetime is infinite.")
    return(Inf)
  }
  1 / est$lambda
}

#' Median lifetime under the exponential model
#'
#' Returns the median lifetime \eqn{\ln(2)/\lambda}: the time at which half
#' of the units are expected to have failed. Always smaller than the mean
#' lifetime, because the exponential distribution is right-skewed.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A single positive `numeric` (median lifetime). Returns `Inf`
#'   with a single warning if no event was observed.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' median_lifetime(d)
#'
#' @seealso [mean_lifetime()], [quantile_lifetime()]
#' @export
median_lifetime <- function(d) {
  .validate_d(d)
  est <- suppressWarnings(estimate_lambda(d))
  if (est$n_events == 0) {
    warning("lambda is 0 (no events observed): median lifetime is infinite.")
    return(Inf)
  }
  log(2) / est$lambda
}

#' Lifetime quantiles under the exponential model
#'
#' Returns the time at which a given proportion `p` of units is expected to
#' have failed. The special case `p = 0.10` gives the B10 life, a common
#' reliability indicator.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#' @param p Numeric vector of probabilities, each strictly between 0 and 1
#'   (no `NA`). Default is `c(0.10, 0.50, 0.90)`.
#'
#' @return A `data.frame` with columns `p` (input probabilities) and
#'   `time` (corresponding lifetimes). If no event was observed, all
#'   lifetimes are `Inf` and a single warning is issued.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' quantile_lifetime(d)
#' quantile_lifetime(d, p = 0.10)   # B10 life
#'
#' @seealso [median_lifetime()], [mean_lifetime()]
#' @export
quantile_lifetime <- function(d, p = c(0.10, 0.50, 0.90)) {
  .validate_d(d)
  if (!is.numeric(p) || length(p) == 0 || anyNA(p) ||
      any(p <= 0) || any(p >= 1)) {
    stop("'p' must be a numeric vector with values strictly between 0 and 1.")
  }
  est <- suppressWarnings(estimate_lambda(d))
  if (est$n_events == 0) {
    warning("lambda is 0 (no events observed): quantiles are infinite.")
    return(data.frame(p = p, time = rep(Inf, length(p))))
  }
  data.frame(
    p    = p,
    time = -log(1 - p) / est$lambda
  )
}

#' Survival and failure probabilities at given times
#'
#' For each time `t`, returns the probability that a unit is still running
#' (survival) and the probability that it has failed (failure) under the
#' fitted exponential model.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#' @param t Numeric vector of non-negative times (no `NA`), in the same
#'   time unit as `d$time`. An empty vector is accepted and gives a
#'   zero-row result.
#'
#' @return A `data.frame` with three columns: `t` (input times), `surv`
#'   (survival probability, between 0 and 1) and `failure` (failure
#'   probability, equal to `1 - surv`). If no event was observed, a
#'   warning is issued (the fitted rate is `0`).
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' survival_prob(d, c(1000, 4000, 8000))
#'
#' @seealso [expected_failures()]
#' @export
survival_prob <- function(d, t) {
  .validate_d(d)
  if (!is.numeric(t) || anyNA(t) || any(t < 0)) {
    stop("'t' must be a numeric vector of non-negative values.")
  }
  lambda <- estimate_lambda(d)$lambda
  data.frame(
    t       = t,
    surv    = exp(-lambda * t),
    failure = 1 - exp(-lambda * t)
  )
}

#' Expected number of failures in a fleet
#'
#' For a fleet of identical units, returns the expected number of failures
#' accumulated by time `t`. Useful to size spare-parts stock or plan
#' maintenance campaigns.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#' @param fleet_size Numeric. Size of the fleet (number of units): a single
#'   finite positive number.
#' @param t Numeric vector of non-negative times (no `NA`), in the same
#'   time unit as `d$time`. An empty vector is accepted and gives a
#'   zero-row result.
#'
#' @return A `data.frame` with columns `fleet_size`, `t` and
#'   `expected_fail`, giving the expected number of failures at each time.
#'   If no event was observed, a warning is issued (the fitted rate is
#'   `0`).
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' expected_failures(d, fleet_size = 500, t = c(1000, 4000, 8000))
#'
#' @seealso [survival_prob()]
#' @export
expected_failures <- function(d, fleet_size, t) {
  .validate_d(d)
  if (!is.numeric(fleet_size) || length(fleet_size) != 1 ||
      !is.finite(fleet_size) || fleet_size <= 0) {
    stop("'fleet_size' must be a single positive number.")
  }
  if (!is.numeric(t) || anyNA(t) || any(t < 0)) {
    stop("'t' must be a numeric vector of non-negative values.")
  }
  lambda <- estimate_lambda(d)$lambda
  data.frame(
    fleet_size    = rep(fleet_size, length(t)),
    t             = t,
    expected_fail = fleet_size * (1 - exp(-lambda * t))
  )
}

#' Compare reliability indicators between groups
#'
#' Splits the prepared data frame by the `group` column and, for each
#' group, computes the number of events, exposure time, failure rate,
#' mean lifetime, median and B10 life. Useful to spot which lot or test
#' condition is the least reliable.
#'
#' Rows whose group is `NA` are ignored. Groups are returned in sorted
#' order (not in order of appearance). A group with no observed failure
#' has a rate of `0` and infinite lifetimes, and triggers one warning.
#'
#' @param d A prepared `data.frame` that also contains a `group` column,
#'   as produced by `prepare_data(..., group_col = ...)`.
#'
#' @return A `data.frame` with one row per group and columns `group`,
#'   `n_total`, `n_events`, `exposure`, `lambda`, `mean_life`, `median`
#'   and `b10`.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' compare_groups(d)
#'
#' @seealso [indicators_table()], [plot_groups()]
#' @export
compare_groups <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event", "group") %in% names(d))) {
    stop("'d' must contain 'time', 'event' and 'group' columns. ",
         "Use prepare_data(..., group_col = ...) first.")
  }
  .validate_d(d)
  groups <- sort(unique(d$group[!is.na(d$group)]))
  if (length(groups) == 0) {
    stop("'group' column contains only NA.")
  }
  rows <- lapply(seq_along(groups), function(i) {
    g   <- groups[i]
    sub <- d[which(d$group == g), , drop = FALSE]
    est <- suppressWarnings(estimate_lambda(sub))
    if (est$n_events == 0) {
      warning("Group '", as.character(g), "': no events observed; ",
              "lambda is 0 and lifetimes are infinite.", call. = FALSE)
    }
    positive <- est$lambda > 0
    data.frame(
      group     = as.character(g),
      n_total   = nrow(sub),
      n_events  = est$n_events,
      exposure  = est$exposure,
      lambda    = est$lambda,
      mean_life = if (positive) 1 / est$lambda else Inf,
      median    = if (positive) log(2) / est$lambda else Inf,
      b10       = if (positive) -log(0.90) / est$lambda else Inf
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Full reliability indicator table
#'
#' Assembles in a single-row `data.frame` all the main indicators of the
#' fitted exponential model: sample size, number of events and censored
#' units, exposure time, failure rate with confidence interval, mean
#' lifetime with confidence interval, median and B10 life. This is the
#' function to use for a quick report.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#' @param level Numeric. Confidence level of the intervals, a single number
#'   strictly between 0 and 1. Default is `0.95`. See
#'   [confidence_interval()].
#'
#' @return A one-row `data.frame` with columns `n_total`, `n_events`,
#'   `n_censored`, `exposure`, `lambda`, `lambda_lower`, `lambda_upper`,
#'   `mean_life`, `mean_lower`, `mean_upper`, `median` and `b10`. If no
#'   event was observed, `mean_life`, `median` and `b10` are `Inf`, the
#'   upper bound of `lambda` is finite, and a single warning is issued.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' indicators_table(d)
#' indicators_table(d, level = 0.90)
#'
#' @seealso [estimate_lambda()], [confidence_interval()],
#'   [compare_groups()], [quantile_lifetime()]
#' @export
indicators_table <- function(d, level = 0.95) {
  .validate_d(d)
  est <- suppressWarnings(estimate_lambda(d))
  ci  <- suppressWarnings(confidence_interval(d, level = level))
  if (est$n_events == 0) {
    warning("No events observed: lambda is 0 and lifetimes are infinite.")
  }
  positive <- est$lambda > 0
  data.frame(
    n_total      = nrow(d),
    n_events     = est$n_events,
    n_censored   = sum(d$event == 0L),
    exposure     = est$exposure,
    lambda       = est$lambda,
    lambda_lower = ci$lambda_lower,
    lambda_upper = ci$lambda_upper,
    mean_life    = if (positive) 1 / est$lambda else Inf,
    mean_lower   = ci$mean_lower,
    mean_upper   = ci$mean_upper,
    median       = if (positive) log(2) / est$lambda else Inf,
    b10          = if (positive) -log(0.90) / est$lambda else Inf
  )
}
