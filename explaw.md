# explaw — Cahier des charges

Package R d'analyse de fiabilité fondé sur la **loi exponentielle** (taux de panne constant λ).
Ce document décrit l'arborescence, les fonctions et leur répartition dans les fichiers.

---

## 0. Décisions de conception validées

| Sujet | Décision | Justification |
|---|---|---|
| Censure | À droite uniquement | Cas le plus courant en fiabilité ; gauche et par intervalle hors périmètre v0.1 |
| Graphiques | `plot.explaw_fit()` public ; `.plot_reliability()`, `.plot_cumhaz()`, `.plot_qq()` internes | Interface unique : `plot(fit, which = ...)` |
| `failure_rate()` | Conservée dans `metrics.R` | Estimateur λ̂ = r/T utile seul (données agrégées « r pannes sur T heures ») ; centralise la formule, appelée par `fit_exp()` |
| `grDevices` | Déclarée dans `Imports` | Fournit `adjustcolor()` pour la bande de confiance semi-transparente ; livré avec R, aucune installation |
| Kaplan-Meier | Via le package `survival` | Fiable, livré avec R, moins de code à maintenir |
| Jeu d'exemple | **Simulé uniquement** (`explaw_example`) | Vrai λ connu → validation des estimations et des IC ; censure contrôlée ; aucune dépendance ni licence tierce |

---

## 1. Conventions générales

- **Unité de temps** : libre (heures, jours, cycles…). Le package ne convertit rien ; l'utilisateur reste cohérent.
- **Statut** : `1` = panne observée, `0` = censuré à droite.
- **Notations** : `r` = nombre de pannes, `T` = temps cumulé de fonctionnement (pannes + censures), `λ` = `rate`.
- **Classes S3** : `explaw_fit`, `summary.explaw_fit`, `explaw_check`, `explaw_plan`, `explaw_system`.
- **Fonctions internes** : préfixe `.` (non exportées, non documentées dans le manuel public).
- **Erreurs** : `stop(..., call. = FALSE)` avec message clair ; avertissements via `warning(..., call. = FALSE)`.
- **Vectorisation** : toutes les métriques sont vectorisées (`t`, `p`, `rate` peuvent être des vecteurs).
- **Documentation** : roxygen2 (`#' @export`, `@param`, `@return`, `@examples`) ; `NAMESPACE` et `man/` sont **générés**, jamais édités à la main.
- **Tests** : testthat 3e édition, un fichier de tests par fichier de code.

### Formules de référence

| Grandeur | Formule |
|---|---|
| Fiabilité | R(t) = exp(−λt) |
| MTTF / MTBF | 1/λ |
| Vie B_p | t_p = −ln(1−p)/λ (= `qexp(p, λ)`) |
| Médiane | t_0.5 = ln(2)/λ |
| Seuil → t\* | t\* = −ln(seuil)/λ = B_(1−seuil) |
| MLE | λ̂ = r/T |
| IC de λ (bilatéral) | borne basse = χ²(α/2 ; 2r)/(2T) ; borne haute = χ²(1−α/2 ; df)/(2T) |
| Degrés de liberté (borne haute) | df = 2r+2 si essai arrêté à date fixe (type I) ; df = 2r si arrêté à r pannes (type II) |
| IC de r = 0 | λ ≤ −ln(1−conf)/T (unilatéral) |
| Disponibilité | A = MTBF / (MTBF + MTTR) |
| Système série | λ_sys = Σ λᵢ |
| Système parallèle | R(t) = 1 − Π (1 − exp(−λᵢt)) |
| k parmi n (identiques) | R(t) = P(X ≥ k), X ~ Binomiale(n, exp(−λt)) |
| Plan de test | T_requis = MTBF_cible × χ²(conf ; 2c+2) / 2 (c = pannes tolérées) |

---

## 2. Arborescence complète

```
Explaw/
├── .Rbuildignore
├── .gitignore
├── DESCRIPTION
├── LICENSE
├── NAMESPACE                    (généré par roxygen2)
├── NEWS.md
├── README.md
├── explaw.md                    (ce document)
├── explaw.Rproj
├── R/
│   ├── explaw-package.R
│   ├── utils.R
│   ├── metrics.R
│   ├── fit_exp.R
│   ├── methods.R
│   ├── plot.R
│   ├── check_exp.R
│   ├── test_plan.R
│   ├── systems.R
│   └── data.R
├── data/
│   └── explaw_example.rda       (généré depuis data-raw/)
├── data-raw/
│   └── make_explaw_data.R
├── man/                         (généré par roxygen2)
├── tests/
│   ├── testthat.R
│   └── testthat/
│       ├── test-utils.R
│       ├── test-metrics.R
│       ├── test-fit_exp.R
│       ├── test-methods.R
│       ├── test-plot.R
│       ├── test-check_exp.R
│       ├── test-test_plan.R
│       └── test-systems.R
└── vignettes/
    └── intro-explaw.Rmd
```

### Fichiers de configuration à prévoir

**DESCRIPTION** (champs clés)

```
Package: explaw
Title: Reliability Analysis Under the Exponential Distribution
Depends: R (>= 4.1.0)
Imports: stats, graphics, grDevices, survival
Suggests: testthat (>= 3.0.0), knitr, rmarkdown
Config/testthat/edition: 3
Encoding: UTF-8
LazyData: true
Roxygen: list(markdown = TRUE)
VignetteBuilder: knitr
```

**.Rbuildignore** : `^explaw\.Rproj$`, `^\.Rproj\.user$`, `^data-raw$`, `^explaw\.md$`, `^\.git`.

---

## 3. Liste des fonctions

Légende : 🔓 exportée · 🔒 interne.

### 3.1 Fonctions internes (`utils.R`)

#### 🔒 `.check_positive(x, name, allow_zero = FALSE)`
Vérifie qu'un argument numérique est valide.
- **Entrées** : `x` (numérique), `name` (chaîne, pour le message d'erreur), `allow_zero` (logique).
- **Traitement** : refuse les non-numériques, `NA`, valeurs négatives (ou nulles si `allow_zero = FALSE`).
- **Sortie** : `invisible(x)` ou erreur.

#### 🔒 `.check_prob(p, name = "p")`
Vérifie que `p` est une probabilité dans l'intervalle ouvert ]0, 1[.
- **Sortie** : `invisible(p)` ou erreur.

#### 🔒 `.check_data(time, status = NULL)`
Valide et normalise les données de durée de vie.
- **Entrées** : `time` (numérique > 0), `status` (0/1, ou `NULL` = tout en panne).
- **Traitement** : vérifie longueurs identiques, absence de `NA`, `time > 0`, `status ∈ {0,1}` ; convertit les logiques en entiers.
- **Sortie** : liste `list(time, status)` ; erreur si aucune donnée.

#### 🔒 `.exp_ci(r, total_time, conf_level = 0.95, test_type = "time")`
Calcule les IC exacts de λ par le χ².
- **Entrées** : `r` (entier ≥ 0), `total_time` (T > 0), `conf_level` (]0,1[), `test_type` (`"time"` ou `"failure"`).
- **Traitement** :
  - `r > 0` : bornes bilatérales selon la table de formules (§1).
  - `r = 0` : borne basse = 0, borne haute unilatérale = −ln(1−conf)/T ; avertissement.
  - `test_type = "failure"` avec `r = 0` : erreur (essai de type II impossible sans panne).
- **Sortie** : `list(rate = c(lower, upper), mttf = c(lower, upper))` (le MTTF est l'inverse, bornes échangées ; `Inf` si borne basse de λ = 0).

---

### 3.2 Métriques (`metrics.R`)

Fonctions pures, vectorisées, sans dépendance entre elles (sauf `stats::qexp`/`pexp`).

#### 🔓 `failure_rate(r, total_time)`
Estime λ à partir d'un comptage (estimateur du maximum de vraisemblance λ̂ = r/T).
- **Usage** : analyse rapide quand on ne dispose que de données agrégées (ex. « 5 pannes sur 45 000 h »), et brique utilisée par `fit_exp()`.
- **Paramètres** : `r` (nombre de pannes, ≥ 0), `total_time` (> 0).
- **Dépend de** : `.check_positive()`.
- **Sortie** : numérique, `r / total_time`.

#### 🔓 `reliability(t, rate)`
Fiabilité R(t).
- **Paramètres** : `t` (≥ 0), `rate` (> 0).
- **Traitement** : `pexp(t, rate, lower.tail = FALSE)`.
- **Sortie** : numérique dans [0, 1], longueur = recyclage de `t` et `rate`.

#### 🔓 `mttf(rate)`
Temps moyen avant panne (identique au MTBF sous l'exponentielle).
- **Paramètres** : `rate` (> 0).
- **Sortie** : `1 / rate`.

#### 🔓 `bx_life(p, rate)`
Temps où une proportion `p` des appareils est en panne. Sert pour t₁₀ (p = 0,1), la médiane (p = 0,5) et t\* (p = 1 − seuil).
- **Paramètres** : `p` (]0,1[), `rate` (> 0).
- **Traitement** : `qexp(p, rate)`.
- **Sortie** : numérique (temps).

#### 🔓 `availability(mtbf, mttr)`
Disponibilité en régime permanent.
- **Paramètres** : `mtbf` (> 0), `mttr` (≥ 0).
- **Sortie** : `mtbf / (mtbf + mttr)`, dans ]0, 1].

---

### 3.3 Estimation (`fit_exp.R`)

#### 🔓 `fit_exp(time, status = NULL, data = NULL, conf_level = 0.95, test_type = c("time", "failure"))`
Fonction centrale : ajuste une exponentielle par maximum de vraisemblance, avec censure à droite.

- **Paramètres** :
  - `time` : vecteur de temps, ou nom de colonne (chaîne) si `data` est fourni.
  - `status` : vecteur 0/1, ou nom de colonne ; `NULL` = toutes les observations sont des pannes.
  - `data` : data.frame optionnel.
  - `conf_level` : niveau de confiance (défaut 0,95).
  - `test_type` : `"time"` (type I, arrêt à date fixe) ou `"failure"` (type II, arrêt à r pannes).
- **Dépend de** : `.check_data()`, `.check_prob()`, `.exp_ci()`, `failure_rate()`, `mttf()`.
- **Traitement** :
  1. Extraction des colonnes si `data` fourni, puis validation.
  2. `r = sum(status)`, `T = sum(time)`, `λ̂ = failure_rate(r, T)`.
  3. IC via `.exp_ci()`.
  4. Assemblage de l'objet.
- **Sortie** : objet de classe `explaw_fit` (liste) :

| Élément | Contenu |
|---|---|
| `rate` | λ̂ |
| `rate_ci` | vecteur c(lower, upper) |
| `mttf` | 1/λ̂ (`Inf` si r = 0) |
| `mttf_ci` | vecteur c(lower, upper) |
| `n` | nombre d'appareils |
| `r` | nombre de pannes |
| `total_time` | T |
| `conf_level` | niveau utilisé |
| `test_type` | type d'essai |
| `data` | data.frame(time, status) nettoyé |
| `call` | appel d'origine |

---

### 3.4 Méthodes S3 (`methods.R`)

Toutes prennent un objet `explaw_fit`.

#### 🔓 `print.explaw_fit(x, digits = 4, ...)`
Affichage court : n, r, λ̂ et MTTF avec IC. **Sortie** : `invisible(x)`.

#### 🔓 `summary.explaw_fit(object, ...)`
Résumé enrichi : calcule en plus t₁₀ et la médiane (avec IC), nombre de censures, taux de censure.
- **Dépend de** : `bx_life()`.
- **Sortie** : objet `summary.explaw_fit` (liste).

#### 🔓 `print.summary.explaw_fit(x, digits = 4, ...)`
Affichage du résumé. **Sortie** : `invisible(x)`.

#### 🔓 `coef.explaw_fit(object, ...)`
**Sortie** : vecteur nommé `c(rate = λ̂)`.

#### 🔓 `confint.explaw_fit(object, parm = c("rate", "mttf"), level = NULL, ...)`
IC du paramètre demandé.
- **Paramètres** : `level` (si `NULL`, utilise `conf_level` du fit ; sinon recalcule via `.exp_ci()`).
- **Sortie** : matrice 1×2 avec colonnes `lower` / `upper`.

#### 🔓 `predict.explaw_fit(object, newdata, type = c("reliability", "quantile"), level = NULL, ...)`
Prédictions avec IC.
- **Paramètres** : `newdata` = temps (si `reliability`) ou probabilités p (si `quantile`).
- **Traitement** : l'IC se déduit de celui de λ par transformation monotone :
  - R(t) : `[exp(−λ_haut·t), exp(−λ_bas·t)]`.
  - t_p : `[bx_life(p, λ_haut), bx_life(p, λ_bas)]`.
- **Dépend de** : `reliability()`, `bx_life()`, `.exp_ci()`.
- **Sortie** : data.frame (`x`, `fit`, `lower`, `upper`).

---

### 3.5 Graphiques (`plot.R`)

#### 🔓 `plot.explaw_fit(x, which = c("reliability", "cumhaz", "qq"), ...)`
Point d'entrée : aiguille vers la fonction interne choisie. **Sortie** : `invisible(x)`.

#### 🔒 `.plot_reliability(fit, threshold = NULL, km = TRUE, conf_band = TRUE, ...)`
R(t) ajusté, avec :
- courbe de Kaplan-Meier superposée (`survival::survfit(Surv(time, status) ~ 1)`) ;
- bande de confiance issue de l'IC de λ ;
- si `threshold` fourni : ligne horizontale au seuil et ligne verticale à t\* = `bx_life(1 − threshold, rate)`, valeur annotée.

#### 🔒 `.plot_cumhaz(fit, ...)`
Graphe de log-survie : points (t, −ln R̂_KM(t)) et droite théorique de pente λ̂ passant par l'origine. Diagnostic visuel principal, valable avec censure.

#### 🔒 `.plot_qq(fit, ...)`
QQ-plot : temps ordonnés contre `qexp(ppoints(n), rate = λ̂)` avec droite y = x. Réservé aux données **sans censure** ; si censures présentes → avertissement et arrêt (ou usage des seules pannes si l'utilisateur force).

- **Dépendances** : `graphics` (plot, lines, abline, legend, text), `grDevices::adjustcolor` (bande semi-transparente ; `R CMD check` impose de déclarer `grDevices` dans `Imports`), `survival` (KM).
- **Repli possible** : si `grDevices` est retiré un jour, la bande peut être remplacée par deux courbes en pointillés.

---

### 3.6 Diagnostic (`check_exp.R`)

#### 🔓 `check_exp(fit, alpha = 0.05)`
Teste si l'hypothèse exponentielle est plausible, en testant β = 1 dans une Weibull.
- **Paramètres** : `fit` (`explaw_fit`), `alpha` (niveau du test).
- **Dépend de** : `survival::survreg`, `survival::Surv`.
- **Traitement** :
  1. Ajuste `survreg(Surv(time, status) ~ 1, dist = "weibull")` et `dist = "exponential"` sur `fit$data`.
  2. β̂ = 1/scale ; IC de Wald sur log(scale), puis transformation.
  3. Test du rapport de vraisemblance : `2·(loglik_weibull − loglik_exp)` ~ χ²(1) ; p-value via `pchisq(lower.tail = FALSE)`.
  4. Verdict textuel : « exponentielle plausible » / « rejetée » ; β̂ > 1 → usure, β̂ < 1 → rodage.
  5. Avertissement si r < 5 (puissance faible).
- **Sortie** : objet `explaw_check` (`beta`, `beta_ci`, `lrt_stat`, `p_value`, `alpha`, `verdict`, `n_failures`).

#### 🔓 `print.explaw_check(x, ...)`
Affichage du verdict. **Sortie** : `invisible(x)`.

---

### 3.7 Plans de test (`test_plan.R`)

#### 🔓 `test_plan(mtbf_target, conf_level = 0.95, failures_allowed = 0, n = NULL, duration = NULL)`
Dimensionne un essai de démonstration de MTBF.
- **Paramètres** : `mtbf_target` (> 0), `conf_level`, `failures_allowed` (c, entier ≥ 0), `n` (nombre d'appareils, optionnel), `duration` (durée par appareil, optionnel).
- **Traitement** :
  1. `T_requis = mtbf_target · qchisq(conf_level, 2c+2) / 2`.
  2. Si `n` donné : `duration = T_requis / n`.
  3. Si `duration` donnée : `n = ceiling(T_requis / duration)`.
  4. Si `n` et `duration` donnés : confiance atteinte = `pchisq(2·n·duration / mtbf_target, 2c+2)`.
  5. Si rien : renvoie seulement `T_requis`.
- **Hypothèse** : temps cumulé = `n × duration` (pire cas, sans panne anticipée).
- **Sortie** : objet `explaw_plan` (`mtbf_target`, `conf_level`, `failures_allowed`, `total_time`, `n`, `duration`, `achieved_conf`).

#### 🔓 `print.explaw_plan(x, ...)`
Résumé lisible du plan. **Sortie** : `invisible(x)`.

---

### 3.8 Systèmes (`systems.R`) — *optionnel pour la v0.1*

Chaque constructeur renvoie un objet `explaw_system` : `list(type, reliability = function(t), mttf, rate, call)` (`rate` = `NA` si le système n'est pas exponentiel).

#### 🔓 `system_series(rates)`
- **Paramètres** : `rates` (vecteur de λᵢ > 0).
- **Traitement** : λ_sys = Σ λᵢ ; le système reste exponentiel.
- **Sortie** : `explaw_system` avec `rate = λ_sys`, `mttf = 1/λ_sys`.

#### 🔓 `system_parallel(rates)`
- **Paramètres** : `rates` (vecteur de λᵢ > 0, n ≤ 15).
- **Traitement** : R(t) = 1 − Π(1 − exp(−λᵢt)) ; MTTF par inclusion-exclusion : Σ(−1)^(|S|+1) / Σ_{i∈S} λᵢ sur les sous-ensembles non vides S.
- **Sortie** : `explaw_system` (non exponentiel, `rate = NA`).

#### 🔓 `system_kofn(k, n, rate)`
- **Paramètres** : `k` (survivants minimum), `n` (composants), `rate` (λ commun).
- **Traitement** : R(t) = `pbinom(k − 1, n, exp(−rate·t), lower.tail = FALSE)` ; MTTF = (1/λ)·Σ_{i=k}^{n} 1/i.
- **Sortie** : `explaw_system`.

#### 🔓 `print.explaw_system(x, ...)`
Affiche type, MTTF et quelques valeurs de R(t). **Sortie** : `invisible(x)`.

---

### 3.9 Données (`data.R` et `data-raw/`)

#### Jeu `explaw_example`
**Décision : jeu simulé uniquement.** Aucun jeu de données externe (Kaggle, `boot`, etc.) n'est embarqué ni requis : pas de dépendance supplémentaire, pas de question de licence. Le vrai paramètre est connu (λ = 1/2000 h⁻¹), ce qui permet de vérifier que l'estimation et les IC sont corrects.

Data.frame de 30 appareils : `id` (entier), `time` (numérique, heures), `status` (0/1).
- **Génération** (`data-raw/make_explaw_data.R`) : `set.seed(2026)`, `rexp(30, rate = 1/2000)`, censure à 2500 h (`status = 0`, `time = 2500` pour les appareils encore en vie).
- **Sauvegarde** : `usethis::use_data(explaw_example, overwrite = TRUE)`.
- **Documentation** : bloc roxygen dans `data.R` (`@format`, `@source`, `"explaw_example"`).

---

## 4. Liste des fichiers et de leur contenu

### Code (`R/`)

| Fichier | Contenu | Dépend de |
|---|---|---|
| `explaw-package.R` | Documentation du package, `@importFrom` (stats, graphics, grDevices, survival), `utils::globalVariables` si besoin | — |
| `utils.R` | 🔒 `.check_positive()`, `.check_prob()`, `.check_data()`, `.exp_ci()` | stats |
| `metrics.R` | 🔓 `failure_rate()`, `reliability()`, `mttf()`, `bx_life()`, `availability()` | utils.R, stats |
| `fit_exp.R` | 🔓 `fit_exp()` | utils.R, metrics.R (`failure_rate()`, `mttf()`) |
| `methods.R` | 🔓 `print.explaw_fit()`, `summary.explaw_fit()`, `print.summary.explaw_fit()`, `coef.explaw_fit()`, `confint.explaw_fit()`, `predict.explaw_fit()` | metrics.R, utils.R |
| `plot.R` | 🔓 `plot.explaw_fit()` ; 🔒 `.plot_reliability()`, `.plot_cumhaz()`, `.plot_qq()` | metrics.R, survival, graphics, grDevices |
| `check_exp.R` | 🔓 `check_exp()`, `print.explaw_check()` | survival, stats |
| `test_plan.R` | 🔓 `test_plan()`, `print.explaw_plan()` | utils.R, stats |
| `systems.R` | 🔓 `system_series()`, `system_parallel()`, `system_kofn()`, `print.explaw_system()` | utils.R, stats |
| `data.R` | Documentation roxygen du jeu `explaw_example` | — |

### Hors `R/`

| Fichier | Rôle |
|---|---|
| `data-raw/make_explaw_data.R` | Script de génération du jeu d'exemple |
| `data/explaw_example.rda` | Jeu de données binaire généré |
| `tests/testthat.R` | Lanceur testthat |
| `tests/testthat/test-*.R` | Un fichier de tests par fichier de code (valeurs connues, cas limites, erreurs) |
| `vignettes/intro-explaw.Rmd` | Tutoriel complet sur `explaw_example` (jeu simulé) : données → `fit_exp()` → `summary` → graphes → `check_exp()` → `test_plan()` |
| `NEWS.md` | Journal des versions |
| `README.md` | Présentation, installation, exemple minimal |

### Tests à prévoir par fichier

| Fichier de test | Cas clés |
|---|---|
| `test-utils.R` | Entrées invalides (NA, négatifs, statut hors 0/1), longueurs différentes |
| `test-metrics.R` | `failure_rate(5, 45000) = 5/45000`, `reliability(0) = 1`, `mttf(0.5) = 2`, `bx_life(0.5) = log(2)/rate`, vectorisation |
| `test-fit_exp.R` | λ̂ = r/T sur petit exemple calculé à la main ; r = 0 ; type I vs II ; data.frame vs vecteurs ; sur `explaw_example`, le vrai MTTF (2000 h) tombe dans l'IC |
| `test-methods.R` | `confint` cohérent avec `fit$rate_ci` ; `predict` monotone ; IC encadrent l'estimation |
| `test-plot.R` | Les graphes s'exécutent sans erreur ; avertissement du QQ-plot avec censure |
| `test-check_exp.R` | Données simulées exponentielles → non rejet ; Weibull β = 3 → rejet |
| `test-test_plan.R` | c = 0, conf = 0,9 → T = MTBF × 2,303 ; cohérence n/durée |
| `test-systems.R` | Série de 2 composants identiques → λ doublé ; k = 1 parmi n ≡ parallèle |

---

## 5. Ordre de développement conseillé

1. **Socle** : `explaw-package.R`, `utils.R`, `DESCRIPTION`, `.Rbuildignore`.
2. **Métriques** : `metrics.R` + tests.
3. **Estimation** : `fit_exp.R` + tests.
4. **Méthodes S3** : `methods.R` + tests.
5. **Graphiques** : `plot.R`.
6. **Diagnostic** : `check_exp.R`.
7. **Plans de test** : `test_plan.R`.
8. **Données** : `data-raw/`, `data.R`.
9. **Systèmes** : `systems.R` (si retenu pour la v0.1).
10. **Documentation** : vignette, README, `devtools::document()`, `devtools::check()`.

Commandes usuelles : `devtools::document()`, `devtools::test()`, `devtools::check()`.
