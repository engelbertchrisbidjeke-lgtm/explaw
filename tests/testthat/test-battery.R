# Batterie de tests de découverte — explaw (loi exponentielle uniquement)
# Placer dans tests/testthat/ puis lancer : devtools::test()
# Étiquettes : [REF] valeur de référence | [BUG] erreur probable | [STAT] validité statistique
#              [ROBUSTESSE] entrées invalides | [DESIGN] choix à trancher | [HYGIENE] paquet/CRAN
# Un test qui échoue = une lacune confirmée. Valeurs de référence calculées indépendamment (scipy).

mk <- function(time, event, group = NULL) {
  d <- data.frame(time = as.numeric(time), event = as.integer(event))
  if (!is.null(group)) d$group <- group
  d
}
n_warn <- function(expr) {
  n <- 0L
  withCallingHandlers(expr, warning = function(w) { n <<- n + 1L; invokeRestart("muffleWarning") })
  n
}
plot_label <- function(p, which) {
  labs <- if (exists("get_labs", envir = asNamespace("ggplot2"))) ggplot2::get_labs(p) else p$labels
  labs[[which]]
}
A <- mk(c(120, 340, 560, 780, 1000, 1500), c(1, 1, 1, 1, 1, 0))   # D = 5, E = 4300

# ---------------------------------------------------------------- [REF] valeurs de référence
test_that("[REF] lambda, exposition et IC exact sur un jeu calculé à la main", {
  expect_equal(exposure_time(A), 4300)
  est <- estimate_lambda(A)
  expect_equal(est$lambda, 5 / 4300, tolerance = 1e-12)
  expect_equal(est$n_events, 5)
  ci <- confidence_interval(A)
  expect_equal(ci$lambda_lower, 0.0003775549744, tolerance = 1e-7)
  expect_equal(ci$lambda_upper, 0.002381764808, tolerance = 1e-7)
  expect_equal(ci$mean_lower, 1 / 0.002381764808, tolerance = 1e-7)
  expect_equal(ci$mean_upper, 1 / 0.0003775549744, tolerance = 1e-7)
})

test_that("[REF] indicateurs dérivés", {
  expect_equal(mean_lifetime(A), 860)
  expect_equal(median_lifetime(A), 860 * log(2))
  q <- quantile_lifetime(A, p = c(0.1, 0.5))
  expect_equal(q$time, 860 * -log(1 - c(0.1, 0.5)))
  expect_equal(survival_prob(A, 1000)$surv, 0.3126125553, tolerance = 1e-8)
  expect_equal(expected_failures(A, 500, 1000)$expected_fail, 343.6937, tolerance = 1e-6)
  tab <- indicators_table(A)
  expect_named(tab, c("n_total", "n_events", "n_censored", "exposure", "lambda", "lambda_lower",
                      "lambda_upper", "mean_life", "mean_lower", "mean_upper", "median", "b10"))
  expect_equal(tab$b10, 860 * -log(0.9))
  expect_true(tab$mean_lower < tab$mean_life && tab$mean_life < tab$mean_upper)
})

test_that("[REF] propriétés : survie décroissante, S + F = 1, IC plus large si niveau plus haut", {
  s <- survival_prob(A, c(0, 10, 100, 1000, Inf))
  expect_true(all(diff(s$surv) <= 0))
  expect_equal(s$surv + s$failure, rep(1, 5))
  w <- function(l) { ci <- confidence_interval(A, l); ci$lambda_upper - ci$lambda_lower }
  expect_gt(w(0.99), w(0.95))
  expect_gt(w(0.95), w(0.80))
})

test_that("[REF] compare_groups : cohérence avec les données poolées", {
  d <- mk(c(10, 20, 30, 40, 50, 60), c(1, 1, 0, 1, 1, 1), c("a", "a", "b", "b", "c", "c"))
  cg <- compare_groups(d)
  expect_equal(sum(cg$n_total), nrow(d))
  expect_equal(sum(cg$exposure), exposure_time(d))
  expect_equal(sum(cg$n_events), sum(d$event))
  expect_true(estimate_lambda(d)$lambda >= min(cg$lambda) && estimate_lambda(d)$lambda <= max(cg$lambda))
})

# ---------------------------------------------------------------- prepare_data
test_that("[BUG] prepare_data : noms de colonnes croisés (time_col = 'event', event_col = 'time')", {
  raw <- data.frame(time = c(1, 0, 1), event = c(10, 20, 30))
  d <- prepare_data(raw, time_col = "event", event_col = "time")
  expect_equal(d$time, c(10, 20, 30))
  expect_equal(d$event, c(1L, 0L, 1L))
})

test_that("[BUG] prepare_data : événement en facteur '0'/'1' ne doit pas être recodé en 1/2", {
  raw <- data.frame(t = c(5, 6, 7), e = factor(c(0, 1, 1)))
  d <- prepare_data(raw, "t", "e")
  expect_equal(d$event, c(0L, 1L, 1L))
})

test_that("[ROBUSTESSE] prepare_data : temps infini rejeté ou retiré", {
  raw <- data.frame(t = c(5, Inf, 7), e = c(1, 1, 0))
  d <- tryCatch(prepare_data(raw, "t", "e"), error = function(e) NULL)
  expect_true(is.null(d) || all(is.finite(d$time)))
})

test_that("[ROBUSTESSE] prepare_data : aucune ligne valide -> erreur claire", {
  raw <- data.frame(t = c(NA, -1), e = c(1, 1))
  expect_error(prepare_data(raw, "t", "e"))
})

test_that("[DESIGN] prepare_data : une colonne 'time' non utilisée n'est pas écrasée en silence", {
  raw <- data.frame(time = c(1, 2, 3), dur = c(10, 20, 30), ev = c(1, 0, 1))
  expect_condition(prepare_data(raw, "dur", "ev"))
})

test_that("prepare_data : entrée non modifiée, logiques acceptés, idempotent", {
  raw <- data.frame(t = c(5, 6, 7), e = c(TRUE, FALSE, TRUE), g = c("a", "b", "a"))
  copie <- raw
  d <- prepare_data(raw, "t", "e", "g")
  expect_identical(raw, copie)
  expect_equal(d$event, c(1L, 0L, 1L))
  expect_equal(prepare_data(d, "time", "event", "group"), d)
})

test_that("prepare_data : NA et temps négatifs retirés, 0 et 2 rejetés/acceptés comme documenté", {
  raw <- data.frame(t = c(5, NA, -1, 8, 0), e = c(1, 1, 1, NA, 0))
  d <- prepare_data(raw, "t", "e")
  expect_equal(nrow(d), 2)
  expect_error(prepare_data(data.frame(t = 1:3, e = c(0, 1, 2)), "t", "e"))
  expect_error(prepare_data(data.frame(t = 1:3, e = c(0, 1, 0.5)), "t", "e"))
})

# ---------------------------------------------------------------- données non préparées
test_that("[ROBUSTESSE] les fonctions refusent un 'd' non conforme (event hors {0,1}, temps < 0)", {
  expect_error(estimate_lambda(data.frame(time = c(1, 2, 3), event = c(0, 1, 2))))
  expect_error(estimate_lambda(data.frame(time = c(-1, 2), event = c(1, 1))))
})

test_that("[ROBUSTESSE] exposition nulle avec pannes : erreur, pas de bornes infinies silencieuses", {
  z <- mk(c(0, 0, 0), c(1, 1, 0))
  expect_error(estimate_lambda(z))
  expect_error(confidence_interval(z))
})

# ---------------------------------------------------------------- aucune panne (D = 0)
test_that("D = 0 : lambda 0, une seule alerte par appel, bornes documentées", {
  z <- mk(c(100, 200, 300), c(0, 0, 0))
  expect_warning(estimate_lambda(z))
  expect_equal(suppressWarnings(mean_lifetime(z)), Inf)
  expect_equal(suppressWarnings(confidence_interval(z))$lambda_lower, 0)
})
test_that("[DESIGN] D = 0 : un seul avertissement par appel (pas de doublon)", {
  z <- mk(c(100, 200, 300), c(0, 0, 0))
  expect_equal(n_warn(mean_lifetime(z)), 1)
  expect_equal(n_warn(median_lifetime(z)), 1)
  expect_equal(n_warn(indicators_table(z)), 1)
})
test_that("[STAT] D = 0 : borne supérieure finie de lambda (-ln(alpha)/E)", {
  z <- mk(c(1000, 2000, 2000), c(0, 0, 0))
  ci <- suppressWarnings(confidence_interval(z, level = 0.90))
  expect_true(is.finite(ci$lambda_upper))
})

# ---------------------------------------------------------------- validation des arguments
test_that("[ROBUSTESSE] arguments NA / hors bornes -> message explicite", {
  expect_error(confidence_interval(A, level = NA_real_), "level")
  expect_error(confidence_interval(A, level = c(0.9, 0.95)), "level")
  expect_error(confidence_interval(A, level = 1), "level")
  expect_error(quantile_lifetime(A, p = NA_real_), "'p'")
  expect_error(quantile_lifetime(A, p = 1), "'p'")
  expect_error(survival_prob(A, t = NA_real_), "'t'")
  expect_error(survival_prob(A, t = -1), "'t'")
  expect_error(expected_failures(A, fleet_size = 0, t = 10), "fleet_size")
  expect_error(expected_failures(A, fleet_size = NA_real_, t = 10), "fleet_size")
})

test_that("survival_prob / expected_failures acceptent un vecteur vide", {
  expect_equal(nrow(survival_prob(A, numeric(0))), 0)
  expect_equal(nrow(expected_failures(A, 100, numeric(0))), 0)
})

# ---------------------------------------------------------------- groupes
test_that("[BUG] compare_groups et plot_groups : groupe avec NA ne fait pas planter", {
  raw <- data.frame(t = c(10, 20, 30, 40, 50, 60), e = c(1, 1, 0, 1, 1, 1),
                    g = c("a", "a", "b", "b", NA, "a"))
  d <- prepare_data(raw, "t", "e", "g")
  expect_error(cg <- compare_groups(d), NA)
  expect_equal(sum(cg$n_total), 5)
  expect_error(plot_groups(d), NA)
})

test_that("compare_groups : groupe sans panne -> lambda 0 et Inf, avertissement", {
  d <- mk(c(10, 20, 30, 40), c(1, 1, 0, 0), c("a", "a", "b", "b"))
  cg <- suppressWarnings(compare_groups(d))
  expect_equal(cg$lambda[cg$group == "b"], 0)
  expect_equal(cg$mean_life[cg$group == "b"], Inf)
  expect_warning(compare_groups(d))
})

test_that("[DESIGN] compare_groups : groupes triés (et non dans l'ordre d'apparition)", {
  d <- mk(c(10, 20, 30, 40), c(1, 1, 1, 1), c("B", "A", "B", "A"))
  expect_equal(compare_groups(d)$group, c("A", "B"))
})

test_that("[DESIGN] indicators_table accepte un niveau de confiance", {
  expect_true("level" %in% names(formals(indicators_table)))
})

# ---------------------------------------------------------------- vérification du modèle
test_that("check_exponential : moins de 5 pannes -> NA et message", {
  r <- check_exponential(mk(c(1, 2, 3, 4), c(1, 1, 1, 1)))
  expect_true(is.na(r$p_value))
  expect_equal(r$n_observed, 4)
})

test_that("check_exponential : statistique et p-value dans [0, 1], champs attendus", {
  set.seed(1)
  x <- rexp(30)
  r <- check_exponential(mk(x, rep(1, 30)))
  expect_named(r, c("statistic", "p_value", "n_observed", "conclusion"))
  expect_true(r$statistic >= 0 && r$statistic <= 1 && r$p_value >= 0 && r$p_value <= 1)
})

test_that("[BUG] check_exponential et plot_qq utilisent le même lambda que estimate_lambda (D/E)", {
  set.seed(10)
  x <- rexp(200); cc <- 0.7                                   # censure à date fixe
  d <- mk(pmin(x, cc), as.integer(x <= cc))
  lam <- estimate_lambda(d)$lambda
  obs <- d$time[d$event == 1L]
  ref <- suppressWarnings(stats::ks.test(obs, "pexp", rate = lam))$p.value
  expect_equal(check_exponential(d)$p_value, ref)
  n <- length(obs); probs <- (seq_len(n) - 0.5) / n
  expect_equal(plot_qq(d)$data$theoretical, sort(-log(1 - probs) / lam))
})

test_that("[STAT] check_exponential : taux de rejet proche de 5 % sous H0 (données complètes)", {
  skip_on_cran(); set.seed(2)
  rej <- mean(vapply(1:600, function(i) check_exponential(mk(rexp(30), rep(1, 30)))$p_value < 0.05, logical(1)))
  expect_gt(rej, 0.02)      # mesuré hors R : ~0,5 % (test trop conservateur car lambda estimé)
  expect_lt(rej, 0.09)
})

test_that("[STAT] check_exponential : détecte une usure marquée (Weibull, forme 3)", {
  skip_on_cran(); set.seed(3)
  pw <- mean(vapply(1:300, function(i) check_exponential(mk(rweibull(40, 3, 1), rep(1, 40)))$p_value < 0.05, logical(1)))
  expect_gt(pw, 0.8)
})

test_that("[STAT] IC de lambda : couverture ~95 % sur données complètes", {
  skip_on_cran(); set.seed(4)
  cov <- mean(vapply(1:3000, function(i) {
    ci <- confidence_interval(mk(rexp(30, 0.01), rep(1, 30)))
    ci$lambda_lower <= 0.01 && 0.01 <= ci$lambda_upper
  }, logical(1)))
  expect_gt(cov, 0.935); expect_lt(cov, 0.965)
})

# ---------------------------------------------------------------- graphiques
test_that("graphiques : objets ggplot constructibles et arguments validés", {
  expect_s3_class(ggplot2::ggplot_build(plot_survival(A)), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(plot_failures(A, 500)), "ggplot_built")
  expect_s3_class(plot_qq(mk(c(1, 2, 3, 4, 5), rep(1, 5))), "ggplot")
  expect_error(plot_survival(A, t_max = 0), "t_max")
  expect_error(plot_failures(A, fleet_size = -1), "fleet_size")
  expect_error(plot_qq(mk(c(1, 2), c(1, 1))))
})

test_that("[BUG] plot_failures : fleet_size non entier ou très grand ne plante pas", {
  expect_error(plot_failures(A, fleet_size = 500.5), NA)
  expect_error(plot_failures(A, fleet_size = 1e10), NA)
})

test_that("[DESIGN] graphiques : ni 'hours' ni 'motors' codés en dur", {
  p <- plot_failures(A, 500)
  expect_false(grepl("motors", plot_label(p, "title")))
  expect_false(grepl("hours", plot_label(p, "x"), ignore.case = TRUE) &&
                 !("time_unit" %in% names(formals(plot_failures))))
})

# ---------------------------------------------------------------- jeu motors (cohérence de l'exemple du README)
test_that("[BUG?] motors : les 10 moteurs à 150 °C n'ont connu aucune panne (sens de 'cens')", {
  skip_if_not_installed("MASS")
  d <- prepare_data(MASS::motors, "time", "cens", "temp")
  cg <- suppressWarnings(compare_groups(d))
  expect_equal(cg$n_events[cg$group == "150"], 0)
})

# ---------------------------------------------------------------- hygiène du paquet
test_that("[HYGIENE] DESCRIPTION : stats déclaré, version de R, pas de LazyData sans data/", {
  skip_on_cran()
  desc <- utils::packageDescription("explaw")
  expect_match(desc$Imports, "stats")
  expect_false(is.null(desc$Depends))
  expect_false(!is.null(desc$LazyData) && !dir.exists(system.file("data", package = "explaw")))
})
test_that("[HYGIENE] 17 fonctions exportées", {
  expect_length(getNamespaceExports("explaw"), 17)
})
