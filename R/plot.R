#' Plot survival probability over time
#'
#' Draws the exponential survival curve \eqn{S(t) = e^{-\lambda t}} for the
#' fitted model: the probability that a unit is still running at time `t`.
#' The curve starts at 1 and decreases smoothly toward 0.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#' @param t_max Numeric or `NULL`. Upper bound of the time axis. Defaults
#'   to the maximum observed time in `d`.
#'
#' @return A `ggplot` object. It can be printed, saved with
#'   [ggplot2::ggsave()], or further customized with `+`.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' p <- plot_survival(d)
#' p
#'
#' @seealso [plot_failures()], [plot_groups()], [survival_prob()]
#' @export
plot_survival <- function(d, t_max = NULL) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  if (is.null(t_max)) t_max <- max(d$time)
  if (!is.numeric(t_max) || length(t_max) != 1 || t_max <= 0) {
    stop("'t_max' must be a single positive number.")
  }
  lambda <- estimate_lambda(d)$lambda
  t_seq <- seq(0, t_max, length.out = 200)
  df <- data.frame(t = t_seq, surv = exp(-lambda * t_seq))
  ggplot2::ggplot(df, ggplot2::aes(x = t, y = surv)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::labs(
      title = "Survival probability",
      subtitle = sprintf("lambda = %.5f per hour", lambda),
      x = "Time (hours)",
      y = "P(still running)"
    ) +
    ggplot2::theme_minimal()
}

#' Plot expected cumulative failures in a fleet
#'
#' Draws the expected number of failures accumulated over time in a fleet
#' of identical units: \eqn{N(t) = \text{fleet\_size} \cdot (1 - e^{-\lambda t})}.
#' Useful to justify spare-parts stock or plan maintenance campaigns.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#' @param fleet_size Numeric. Size of the fleet. Default is `1000`.
#' @param t_max Numeric or `NULL`. Upper bound of the time axis. Defaults
#'   to the maximum observed time in `d`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' p <- plot_failures(d, fleet_size = 500)
#' p
#'
#' @seealso [plot_survival()], [expected_failures()]
#' @export
plot_failures <- function(d, fleet_size = 1000, t_max = NULL) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  if (!is.numeric(fleet_size) || length(fleet_size) != 1 ||
      fleet_size <= 0) {
    stop("'fleet_size' must be a single positive number.")
  }
  if (is.null(t_max)) t_max <- max(d$time)
  if (!is.numeric(t_max) || length(t_max) != 1 || t_max <= 0) {
    stop("'t_max' must be a single positive number.")
  }
  lambda <- estimate_lambda(d)$lambda
  t_seq <- seq(0, t_max, length.out = 200)
  df <- data.frame(
    t     = t_seq,
    fails = fleet_size * (1 - exp(-lambda * t_seq))
  )
  ggplot2::ggplot(df, ggplot2::aes(x = t, y = fails)) +
    ggplot2::geom_line(color = "darkorange", linewidth = 1) +
    ggplot2::labs(
      title = sprintf("Expected failures in a fleet of %d motors", fleet_size),
      x = "Time (hours)",
      y = "Expected number of failures"
    ) +
    ggplot2::theme_minimal()
}

#' Plot survival curves by group
#'
#' Superimposes the exponential survival curve of each group on a single
#' plot. The group whose curve drops fastest is the least reliable one.
#'
#' @param d A prepared `data.frame` that also contains a `group` column,
#'   as produced by `prepare_data(..., group_col = ...)`.
#' @param t_max Numeric or `NULL`. Upper bound of the time axis. Defaults
#'   to the maximum observed time in `d`.
#'
#' @return A `ggplot` object with one colored line per group.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' p <- plot_groups(d)
#' p
#'
#' @seealso [plot_survival()], [compare_groups()]
#' @export
plot_groups <- function(d, t_max = NULL) {
  if (!is.data.frame(d) ||
      !all(c("time", "event", "group") %in% names(d))) {
    stop("'d' must contain 'time', 'event' and 'group' columns. ",
         "Use prepare_data(..., group_col = ...) first.")
  }
  if (is.null(t_max)) t_max <- max(d$time)
  if (!is.numeric(t_max) || length(t_max) != 1 || t_max <= 0) {
    stop("'t_max' must be a single positive number.")
  }
  groups <- unique(d$group)
  groups <- groups[!is.na(groups)]
  if (length(groups) == 0) {
    stop("'group' column contains only NA.")
  }
  t_seq <- seq(0, t_max, length.out = 200)
  rows <- lapply(groups, function(g) {
    sub <- d[d$group == g, , drop = FALSE]
    lambda <- estimate_lambda(sub)$lambda
    data.frame(
      t     = t_seq,
      surv  = exp(-lambda * t_seq),
      group = as.character(g)
    )
  })
  df <- do.call(rbind, rows)
  ggplot2::ggplot(df, ggplot2::aes(x = t, y = surv, color = group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = "Survival probability by group",
      x = "Time (hours)",
      y = "P(still running)",
      color = "Group"
    ) +
    ggplot2::theme_minimal()
}

#' Q-Q plot against the exponential distribution
#'
#' Compares the observed lifetimes (only uncensored units) to the quantiles
#' of the fitted exponential distribution. If the points lie close to the
#' dashed line, the exponential model is a plausible fit. Curvature
#' suggests wear-out or infant-mortality patterns.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A `ggplot` object.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' p <- plot_qq(d)
#' p
#'
#' @seealso [check_exponential()]
#' @export
plot_qq <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  obs <- d$time[d$event == 1L]
  if (length(obs) < 3) {
    stop("Not enough observed failures (need at least 3) for a Q-Q plot.")
  }
  lambda <- 1 / mean(obs)
  n <- length(obs)
  probs <- (seq_len(n) - 0.5) / n
  theo <- -log(1 - probs) / lambda
  df <- data.frame(
    theoretical = sort(theo),
    observed    = sort(obs)
  )
  ggplot2::ggplot(df, ggplot2::aes(x = theoretical, y = observed)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         color = "grey50", linetype = "dashed") +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::labs(
      title = "Q-Q plot vs exponential distribution",
      subtitle = "Points close to the line suggest an exponential fit",
      x = "Theoretical quantiles (hours)",
      y = "Observed lifetimes (hours)"
    ) +
    ggplot2::theme_minimal()
}
