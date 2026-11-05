#' Check compatibility with the exponential model
#'
#' Tests whether the observed lifetimes (uncensored units only) are
#' compatible with an exponential distribution, using a
#' Kolmogorov-Smirnov test against the fitted exponential distribution.
#' A small p-value (typically below 0.05) means the exponential model
#' is rejected, and a more flexible model (Weibull, gamma, log-normal)
#' should be considered. This is a diagnostic check, not a formal proof:
#' failing to reject does not prove the model is correct.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A `list` with four elements:
#'   \describe{
#'     \item{`statistic`}{Numeric. Kolmogorov-Smirnov test statistic, or
#'       `NA` if too few observed failures.}
#'     \item{`p_value`}{Numeric. Test p-value, or `NA` if too few
#'       observed failures.}
#'     \item{`n_observed`}{Integer. Number of uncensored observations
#'       used in the test.}
#'     \item{`conclusion`}{Character. Ready-to-read sentence stating
#'       whether the exponential model is rejected or not.}
#'   }
#'   If fewer than 5 failures are observed, `statistic` and `p_value` are
#'   `NA`, and `conclusion` explains why.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' result <- check_exponential(d)
#' result$conclusion
#'
#' @seealso [plot_qq()]
#' @export
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
  test <- suppressWarnings(stats::ks.test(obs, "pexp", rate = lambda))
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
