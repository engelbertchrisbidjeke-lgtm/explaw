# Internal helper (not exported, no documentation page).
# Kolmogorov-Smirnov distance between a sample and an exponential
# distribution with the given rate (default: 1 / mean(x)).
.exp_ks_stat <- function(x, rate = 1 / mean(x)) {
  n <- length(x)
  u <- stats::pexp(sort(x), rate = rate)
  i <- seq_len(n)
  max(i / n - u, u - (i - 1) / n)
}

# Internal helper (not exported, no documentation page).
# Cache of simulated null distributions, one entry per sample size.
.lillie_env <- new.env(parent = emptyenv())

# Internal helper (not exported, no documentation page).
# Simulated null distribution of the KS statistic when the exponential rate
# is estimated from the sample itself (Lilliefors-type correction).
# Uses a fixed internal seed and restores the caller's RNG state, so the
# result is reproducible and the user's random numbers are not disturbed.
# Computed once per sample size, then cached.
.lillie_null <- function(n, B = 5000L) {
  key <- paste(n, B, sep = "_")
  if (!is.null(.lillie_env[[key]])) return(.lillie_env[[key]])

  glob <- globalenv()
  had_seed <- exists(".Random.seed", envir = glob, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = glob, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = glob)
    } else if (exists(".Random.seed", envir = glob, inherits = FALSE)) {
      rm(".Random.seed", envir = glob)
    }
  }, add = TRUE)

  set.seed(2024L)
  sim   <- numeric(B)
  chunk <- max(1L, floor(1e6 / n))
  done  <- 0L
  while (done < B) {
    k <- min(chunk, B - done)
    x <- matrix(stats::rexp(n * k), nrow = n, ncol = k)
    sim[done + seq_len(k)] <- apply(x, 2L, .exp_ks_stat)
    done <- done + k
  }
  .lillie_env[[key]] <- sim
  sim
}

#' Check compatibility with the exponential model
#'
#' Tests whether the observed lifetimes (uncensored units only) are
#' compatible with an exponential distribution, using the
#' Kolmogorov-Smirnov distance to the fitted exponential distribution.
#' The rate used is the same as in [estimate_lambda()] (number of failures
#' divided by total exposure time).
#' A small p-value (typically below 0.05) means the exponential model
#' is rejected, and a more flexible model (Weibull, gamma, log-normal)
#' should be considered. This is a diagnostic check, not a formal proof:
#' failing to reject does not prove the model is correct.
#'
#' @details
#' The way the p-value is computed depends on the data:
#' \itemize{
#'   \item **No censored unit**: the rate is estimated from the very
#'     sample being tested, which makes the classical Kolmogorov-Smirnov
#'     p-value too conservative. A Lilliefors-type correction is applied:
#'     the p-value is obtained from a simulated null distribution (5000
#'     samples of the same size, fixed internal seed, computed once per
#'     sample size). Your own random-number generator state is left
#'     unchanged.
#'   \item **Censored units present**: the classical
#'     [stats::ks.test()] p-value is returned, computed on the uncensored
#'     lifetimes against the exponential distribution with the rate
#'     \eqn{D / E}. It is only an approximation, because censoring is not
#'     taken into account in the test.
#' }
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A `list` with four elements:
#'   \describe{
#'     \item{`statistic`}{Numeric. Kolmogorov-Smirnov test statistic, or
#'       `NA` if too few observed failures.}
#'     \item{`p_value`}{Numeric. Test p-value (see Details), or `NA` if
#'       too few observed failures.}
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
#' @seealso [plot_qq()], [estimate_lambda()]
#' @export
check_exponential <- function(d) {
  .validate_d(d)
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

  lambda <- estimate_lambda(d)$lambda

  if (all(d$event == 1L)) {
    statistic  <- .exp_ks_stat(obs, rate = lambda)
    null_stats <- .lillie_null(length(obs))
    p_value    <- (1 + sum(null_stats >= statistic)) / (length(null_stats) + 1)
  } else {
    test      <- suppressWarnings(stats::ks.test(obs, "pexp", rate = lambda))
    statistic <- unname(test$statistic)
    p_value   <- unname(test$p.value)
  }

  conclusion <- if (p_value < 0.05) {
    paste0(
      "Exponential model rejected at 5% level (p = ",
      format.pval(p_value, digits = 3),
      "). Consider a Weibull or gamma model."
    )
  } else {
    paste0(
      "No incompatibility detected with the exponential model (p = ",
      format.pval(p_value, digits = 3),
      ")."
    )
  }
  list(
    statistic  = statistic,
    p_value    = p_value,
    n_observed = length(obs),
    conclusion = conclusion
  )
}
