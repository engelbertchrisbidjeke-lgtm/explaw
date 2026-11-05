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
    t      = t_seq,
    fails  = fleet_size * (1 - exp(-lambda * t_seq))
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
