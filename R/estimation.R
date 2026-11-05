exposure_time <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  sum(d$time)
}

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
      lambda_lower  = 0,
      lambda_upper  = Inf,
      mean_lower    = 0,
      mean_upper    = Inf
    ))
  }
  alpha <- (1 - level) / 2
  lambda_lower <- qchisq(alpha,       2 * n_events) / (2 * exposure)
  lambda_upper <- qchisq(1 - alpha,   2 * n_events) / (2 * exposure)
  data.frame(
    lambda_lower = lambda_lower,
    lambda_upper = lambda_upper,
    mean_lower   = 1 / lambda_upper,
    mean_upper   = 1 / lambda_lower
  )
}
