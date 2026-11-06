# explaw

**Reliability Analysis Package Based on Exponential Law Principles**  
**Package d'Analyse Fiabilité via les Principes de Loi Exponentielles**

`explaw` fournit une boîte à outils complète pour l'analyse de fiabilité sous hypothèse de loi exponentielle : préparation des données, estimation, intervalles de confiance, indicateurs (moyenne, médiane, quantiles, probabilités de survie), diagnostics, comparaisons de groupes et graphiques prêts à publier.

Le modèle exponentiel (taux de défaillance constant, propriété d'absence de mémoire) est la référence de base en fiabilité. Ce package permet de l'ajuster, le tester et le rapporter avec des IC exacts (loi χ²) et des diagnostics visuels clairs.

---

## Installation

```r
devtools::install_local("path/to/explaw")
# ou après build :
install.packages("explaw_0.1.0.tar.gz", repos = NULL, type = "source")
```

**Dépendances** : `ggplot2` (graphiques), `MASS` (jeu de données `motors` pour exemples/tests).

---

## Démarrage rapide / Quick Start

```r
library(explaw)
data(motors, package = "MASS")   # 40 moteurs, temps en heures, censure (1=panne, 0=censuré), groupe temp

# 1. Préparer les données
d <- prepare_data(motors, time_col = "time", event_col = "cens", group_col = "temp")

# 2. Résumé rapide
summarize_data(d)

# 3. Tableau complet d'indicateurs (rapport one-liner)
indicators_table(d)

# 4. Vérifier l'hypothèse exponentielle
check_exponential(d)$conclusion

# 5. Diagnostics visuels
plot_qq(d)
plot_survival(d)
plot_groups(d)   # nécessite group_col dans prepare_data()
```

---

## Référence des fonctions par module

---

### 1. Préparation des données

Standardisent les données brutes au format canonique (`time`, `event`, `group` optionnel) attendu par toutes les autres fonctions.

| Fonction | Rôle |
|----------|------|
| `prepare_data(data, time_col, event_col, group_col = NULL)` | Renomme/valide les colonnes, supprime lignes incomplètes ou temps négatifs, assure `event` ∈ {0,1}, place les colonnes clés en premier. |
| `summarize_data(d)` | Aperçu une ligne : `n_total`, `n_events`, `n_censored`, `pct_censored`, `time_min`, `time_max`. |

**Points clés**  
- Point d'entrée unique : appeler `prepare_data()` une fois, transmettre son résultat partout.  
- `group_col` active les fonctions de groupe (`compare_groups()`, `plot_groups()`).  
- Lignes avec `NA` ou temps < 0 supprimées silencieusement.

---

### 2. Estimation de base

Quantités fondamentales du modèle exponentiel : temps d'exposition, taux de défaillance (λ), IC exacts.

| Fonction | Retour | Note |
|----------|--------|------|
| `exposure_time(d)` | `numeric` (temps total à risque) | Somme de `time` ; unités censurées contribuent leur temps observé. |
| `estimate_lambda(d)` | `list(lambda, n_events, exposure)` | MLE : `λ = D / E`. Avertissement si `D = 0` (retourne `lambda = 0`). |
| `confidence_interval(d, level = 0.95)` | `data.frame(lambda_lower, lambda_upper, mean_lower, mean_upper)` | IC exact via `2 D λ ~ χ²_{2D}`. Si `D = 0` → `[0, Inf]`. |

**Base statistique**  
- Exposition `E = Σ time_i` (toutes unités).  
- Événements `D = Σ event_i`.  
- MLE `λ̂ = D / E`.  
- Quantité pivote : `2 D λ̂ / λ ~ χ²_{2D}` → IC exact pour λ ; inversion pour moyenne `1/λ`.

---

### 3. Indicateurs de fiabilité

Grandeurs dérivées pour rapport, décision, gestion de flotte.

| Fonction | Sortie | Usage typique |
|----------|--------|---------------|
| `mean_lifetime(d)` | `numeric` (E[T] = 1/λ) | MTBF / MTTF. |
| `median_lifetime(d)` | `numeric` (ln(2)/λ) | Temps où 50 % de la population a défailli. |
| `quantile_lifetime(d, p = c(0.10, 0.50, 0.90))` | `data.frame(p, time)` | B10 (`p=0.10`), médiane, B90. |
| `survival_prob(d, t)` | `data.frame(t, surv, failure)` | Probabilité de fonctionnement aux instants `t`. |
| `expected_failures(d, fleet_size, t)` | `data.frame(fleet_size, t, expected_fail)` | Dimensionnement pièces de rechange, planification maintenance. |
| `compare_groups(d)` | `data.frame` (1 ligne/groupe) | Comparaison lots / conditions d'essai. |
| `indicators_table(d)` | `data.frame` (1 ligne, 12 colonnes) | **Tableau de bord exécutif** — tout en un. |

**Colonnes de `indicators_table()`** : `n_total`, `n_events`, `n_censored`, `exposure`, `lambda`, `lambda_lower`, `lambda_upper`, `mean_life`, `mean_lower`, `mean_upper`, `median`, `b10`.

**Notes**  
- Toutes retournent `Inf` (avec avertissement) si `D = 0`.  
- `compare_groups()` requiert colonne `group` issue de `prepare_data(..., group_col = ...)`.  
- `survival_prob()` et `expected_failures()` acceptent `t` vectoriel.

---

### 4. Vérification du modèle

Diagnostics pour juger la tenue de l'hypothèse exponentielle.

| Fonction | Sortie | Interprétation |
|----------|--------|----------------|
| `check_exponential(d)` | `list(statistic, p_value, n_observed, conclusion)` | Test KS sur durées non censurées vs exponentielle ajustée. `p < 0.05` → rejeter ; envisager Weibull, gamma, log-normale. |
| `plot_qq(d)` | `ggplot` | Q-Q plot observé vs théorique exponentiel. Points sur la ligne = bon ajustement ; courbure systématique = misspécification. |

**Mises en garde**  
- KS utilise λ *estimé* (pas fixé) → p-values approchées.  
- Exige ≥ 5 pannes observées ; sinon `NA` + message explicatif.  
- Ne pas rejeter ≠ preuve d'exponentiel ; seulement absence de preuve forte contre.  
- Coupler test et `plot_qq()` visuel.

---

### 5. Visualisation

Tous les graphiques retournent des objets `ggplot` — affichables, sauvables (`ggsave()`), extensibles (`+ theme_*()`, `+ labs()`, etc.).

| Fonction | Description | Arguments clés |
|----------|-------------|----------------|
| `plot_survival(d, t_max = NULL)` | Courbe de survie `S(t) = exp(-λt)` | `t_max` = max observé par défaut. |
| `plot_failures(d, fleet_size = 1000, t_max = NULL)` | Pannes cumulées attendues `N(t) = fleet_size × (1 - exp(-λt))` | `fleet_size` échelle la courbe ; justification stocks. |
| `plot_groups(d, t_max = NULL)` | Courbes de survie superposées par groupe | Colonne `group` requise ; courbe la plus pentue = groupe le moins fiable. |
| `plot_qq(d)` | Q-Q plot vs exponentiel | Voir section 4. |

**Style** : `theme_minimal()`, couleurs steelblue (survie), darkorange (pannes), palette qualitative (groupes), axes étiquetés avec unités (heures par défaut).

---

## Exemple de workflow

```r
library(explaw)
data(motors, package = "MASS")

# --- Préparation -----------------------------------------------------
d <- prepare_data(motors, "time", "cens", "temp")

# --- Exploration -----------------------------------------------------
summarize_data(d)
#   n_total n_events n_censored pct_censored time_min time_max
# 1      40       39          1          2.5      150    8933

# --- Estimation & Rapport -------------------------------------------
indicators_table(d)
#   n_total n_events n_censored exposure    lambda lambda_lower lambda_upper
# 1      40       39          1   174793 0.0002231   0.0001586   0.0003072
#   mean_life mean_lower mean_upper  median     b10
# 1      4482       3255       6303   3104   473.1

# --- Vérification hypothèse -----------------------------------------
check_exponential(d)$conclusion
# "No incompatibility detected with the exponential model (p = 0.487)."

plot_qq(d)            # confirmation visuelle
plot_survival(d)      # courbe S(t)
plot_groups(d)        # comparaison groupes température
plot_failures(d, fleet_size = 500)  # flotte de 500 moteurs
```

---

## Formules de référence

| Grandeur | Formule | Estimation |
|----------|---------|------------|
| Taux de défaillance | `λ` | `λ̂ = D / E` (MLE) |
| Durée de vie moyenne | `1/λ` | `1/λ̂` |
| Médiane | `ln(2)/λ` | `ln(2)/λ̂` |
| Quantile `p` | `-ln(1-p)/λ` | `-ln(1-p)/λ̂` |
| Survie `S(t)` | `exp(-λt)` | `exp(-λ̂ t)` |
| Pannes attendues | `N × (1 - exp(-λt))` | `N × (1 - exp(-λ̂ t))` |
| IC pour `λ` (niveau `1-α`) | `[ χ²_{α/2, 2D} / 2E , χ²_{1-α/2, 2D} / 2E ]` | Exact (pivote χ²) |
| IC pour moyenne `1/λ` | Inversion bornes λ | Exact |

---

## Cas fréquents

**Aucun événement observé (`D = 0`)**  
- `estimate_lambda()` → `lambda = 0` + avertissement.  
- `confidence_interval()` → `[0, Inf]` pour λ et pour la moyenne.  
- `mean_lifetime()`, `median_lifetime()`, `quantile_lifetime()` → `Inf`.  
- Seules des *bornes inférieures* sur la durée de vie ont du sens (ex. `exposure_time()` = temps total d'essai sans panne).

**Données censurées**  
- Censure à droite gérée naturellement : unités censurées apportent leur temps à l'exposition `E` mais pas aux événements `D`.  
- `prepare_data()` valide que `event` ne contient que 0/1.

**Comparaisons de groupes**  
- Passer `group_col` à `prepare_data()`.  
- `compare_groups()` ajuste un λ par groupe.  
- `plot_groups()` affiche les courbes de survie côte à côte.

---

## Tests

```r
devtools::test()
# ou
testthat::test_dir("tests/testthat")
```

Couverture : préparation (cas limites), estimation avec/sans événements, IC, indicateurs, groupes, diagnostics, structure des objets graphiques.

---

## Licence

Apache License 2.0 — voir [LICENSE.md](LICENSE.md).

---

*Développé avec R ≥ 4.0, `ggplot2`, `roxygen2`.*
