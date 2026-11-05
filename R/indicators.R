#' Mean lifetime under the exponential model
#'
#' Returns the expected lifetime \eqn{E(T) = 1/\lambda} under an exponential
#' model with rate estimated from the data.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A single positive `numeric` (mean lifetime, in the same time
#'   unit as `d$time`). Returns `Inf` with a warning if no event was
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
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  est <- estimate_lambda(d)
  if (est$lambda <= 0) {
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
#'   with a warning if no event was observed.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' median_lifetime(d)
#'
#' @seealso [mean_lifetime()], [quantile_lifetime()]
#' @export
median_lifetime <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  est <- estimate_lambda(d)
  if (est$lambda <= 0) {
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
#' @param p Numeric vector of probabilities, each strictly between 0 and 1.
#'   Default is `c(0.10, 0.50, 0.90)`.
#'
#' @return A `data.frame` with columns `p` (input probabilities) and
#'   `time` (corresponding lifetimes). If no event was observed, all
#'   lifetimes are `Inf` and a warning is issued.
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
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  if (!is.numeric(p) || any(p <= 0) || any(p >= 1)) {
    stop("'p' must be a numeric vector with values strictly between 0 and 1.")
  }
  est <- estimate_lambda(d)
  if (est$lambda <= 0) {
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
#' @param t Numeric vector of non-negative times, in the same time unit as
#'   `d$time`.
#'
#' @return A `data.frame` with three columns: `t` (input times), `surv`
#'   (survival probability, between 0 and 1) and `failure` (failure
#'   probability, equal to `1 - surv`).
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' survival_prob(d, c(1000, 4000, 8000))
#'
#' @seealso [expected_failures()]
#' @export
survival_prob <- function(d, t) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  if (!is.numeric(t) || any(t < 0)) {
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
#' @param fleet_size Numeric. Size of the fleet (number of units).
#' @param t Numeric vector of non-negative times, in the same time unit as
#'   `d$time`.
#'
#' @return A `data.frame` with columns `fleet_size`, `t` and
#'   `expected_fail`, giving the expected number of failures at each time.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' expected_failures(d, fleet_size = 500, t = c(1000, 4000, 8000))
#'
#' @seealso [survival_prob()]
#' @export
expected_failures <- function(d, fleet_size, t) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  if (!is.numeric(fleet_size) || length(fleet_size) != 1 ||
      fleet_size <= 0) {
    stop("'fleet_size' must be a single positive number.")
  }
  if (!is.numeric(t) || any(t < 0)) {
    stop("'t' must be a numeric vector of non-negative values.")
  }
  lambda <- estimate_lambda(d)$lambda
  data.frame(
    fleet_size    = fleet_size,
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
  groups <- unique(d$group)
  groups <- groups[!is.na(groups)]
  if (length(groups) == 0) {
    stop("'group' column contains only NA.")
  }
  rows <- lapply(groups, function(g) {
    sub <- d[d$group == g, , drop = FALSE]
    est <- estimate_lambda(sub)
    med <- if (est$lambda > 0) log(2) / est$lambda else Inf
    b10 <- if (est$lambda > 0) -log(0.90) / est$lambda else Inf
    data.frame(
      group     = as.character(g),
      n_total   = nrow(sub),
      n_events  = est$n_events,
      exposure  = est$exposure,
      lambda    = est$lambda,
      mean_life = if (est$lambda > 0) 1 / est$lambda else Inf,
      median    = med,
      b10       = b10
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
#'
#' @return A one-row `data.frame` with columns `n_total`, `n_events`,
#'   `n_censored`, `exposure`, `lambda`, `lambda_lower`, `lambda_upper`,
#'   `mean_life`, `mean_lower`, `mean_upper`, `median` and `b10`.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' indicators_table(d)
#'
#' @seealso [estimate_lambda()], [confidence_interval()],
#'   [compare_groups()], [quantile_lifetime()]
#' @export
indicators_table <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  est <- estimate_lambda(d)
  ci  <- confidence_interval(d)
  med <- if (est$lambda > 0) log(2) / est$lambda else Inf
  b10 <- if (est$lambda > 0) -log(0.90) / est$lambda else Inf
  data.frame(
    n_total      = nrow(d),
    n_events     = est$n_events,
    n_censored   = sum(d$event == 0L),
    exposure     = est$exposure,
    lambda       = est$lambda,
    lambda_lower = ci$lambda_lower,
    lambda_upper = ci$lambda_upper,
    mean_life    = if (est$lambda > 0) 1 / est$lambda else Inf,
    mean_lower   = ci$mean_lower,
    mean_upper   = ci$mean_upper,
    median       = med,
    b10          = b10
  )
}
