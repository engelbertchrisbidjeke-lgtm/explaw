check_exponential <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  obs <- d$time[d$event == 1L]
  if (length(obs) < 5) {
    return(list(
      statistic  = NA_real_,
      p_value    = NA_real_,
      n_observed = length(obs),
      conclusion = paste0(
        "Too few observed failures (n = ", length(obs),
        "). At least 5 are needed to test the exponential fit."
      )
    ))
  }
  lambda <- 1 / mean(obs)
  test <- suppressWarnings(ks.test(obs, "pexp", rate = lambda))
  conclusion <- if (test$p.value < 0.05) {
    paste0(
      "Exponential model rejected at 5% level (p = ",
      format.pval(test$p.value, digits = 3),
      "). Consider a Weibull or gamma model."
    )
  } else {
    paste0(
      "No incompatibility detected with the exponential model (p = ",
      format.pval(test$p.value, digits = 3),
      ")."
    )
  }
  list(
    statistic  = unname(test$statistic),
    p_value    = unname(test$p.value),
    n_observed = length(obs),
    conclusion = conclusion
  )
}
