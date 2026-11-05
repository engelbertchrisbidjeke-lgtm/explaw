test_that("prepare_data renames columns", {
  d <- prepare_data(MASS::motors, "time", "cens", "temp")
  expect_true(all(c("time", "event", "group") %in% names(d)))
  expect_false("cens" %in% names(d))
  expect_false("temp" %in% names(d))
})

test_that("prepare_data rejects non-binary event", {
  bad <- MASS::motors
  bad$cens[1] <- 2
  expect_error(prepare_data(bad, "time", "cens"))
})

test_that("estimate_lambda computes events over exposure", {
  d <- prepare_data(MASS::motors, "time", "cens")
  est <- estimate_lambda(d)
  expect_equal(est$lambda, est$n_events / est$exposure)
})

test_that("estimate_lambda warns when no events", {
  no_ev <- MASS::motors
  no_ev$cens <- 0
  d <- prepare_data(no_ev, "time", "cens")
  expect_warning(estimate_lambda(d))
})

test_that("confidence_interval contains point estimate", {
  d <- prepare_data(MASS::motors, "time", "cens")
  ci <- confidence_interval(d)
  est <- estimate_lambda(d)
  expect_true(ci$lambda_lower < est$lambda)
  expect_true(ci$lambda_upper > est$lambda)
})

test_that("mean_lifetime equals 1 over lambda", {
  d <- prepare_data(MASS::motors, "time", "cens")
  est <- estimate_lambda(d)
  expect_equal(mean_lifetime(d), 1 / est$lambda)
})

test_that("survival_prob equals 1 at time zero", {
  d <- prepare_data(MASS::motors, "time", "cens")
  sp <- survival_prob(d, 0)
  expect_equal(sp$surv, 1)
  expect_equal(sp$failure, 0)
})

test_that("indicators_table has one row and 12 columns", {
  d <- prepare_data(MASS::motors, "time", "cens")
  tab <- indicators_table(d)
  expect_equal(nrow(tab), 1)
  expect_equal(ncol(tab), 12)
})

test_that("compare_groups returns one row per group", {
  d <- prepare_data(MASS::motors, "time", "cens", "temp")
  comp <- compare_groups(d)
  expect_equal(nrow(comp), 4)
  expect_true("lambda" %in% names(comp))
})

test_that("plot_survival returns a ggplot", {
  d <- prepare_data(MASS::motors, "time", "cens")
  p <- plot_survival(d)
  expect_s3_class(p, "ggplot")
})

test_that("check_exponential returns a 4-element list", {
  d <- prepare_data(MASS::motors, "time", "cens")
  res <- check_exponential(d)
  expect_type(res, "list")
  expect_length(res, 4)
  expect_true("conclusion" %in% names(res))
})

test_that("full pipeline runs without error", {
  d <- prepare_data(MASS::motors, "time", "cens", "temp")
  tab <- indicators_table(d)
  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 1)
})
