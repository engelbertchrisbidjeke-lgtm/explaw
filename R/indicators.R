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
