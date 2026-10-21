# miMeta 2.1.0

## New features

* Added **SMESH**, a structure-learning meta-analysis that learns latent clusters of contexts and
  estimates all-cluster-shared and context-dependent microbial signatures. Unlike a fixed-effect
  meta-analysis, which summarizes each feature through one overall association, SMESH separates
  signatures whose association is selected as nonzero in every cluster from those selected in only
  a subset of clusters.

* `smesh()` runs the whole procedure from sequencing read counts: it generates per-study
  association estimates and standard errors with PALM, combines the studies of each context by
  fixed-effect meta-analysis into one context-level summary, and fits the SMESH model. The `context` argument maps studies to contexts;
  `context = NULL` treats every study as its own context.

* `smesh.summary()` fits the mixture directly on context-level summary statistics, for summaries
  you already have or alternative summary statistics. By default the number of clusters `G` is
  selected by the minimum generalized information criterion (GIC) over the candidate values, and
  `support.ic` sets the information criterion for the sparse association model at each candidate
  `G`. The clustering is fitted from `nstart = 20` random starts by default.

* Both return, per covariate of interest, `effects` (the final cluster-level effect estimates from
  the unpenalized post-selection refit), `context_assignments`, `posterior_membership`,
  `mixing_proportions`, `gic`, and `diagnostics`. The diagnostics hold the GIC reached by each
  random start (`restart_gic`), whether it converged (`restart_convergence`), the start retained
  for the reported fit (`selected_restart`) and, when `G` is selected, the GIC of every candidate
  `G` (`gic_path`).

* With more than two clusters, the second stage builds a dense design over all `2^G - 1` cluster
  subsets. `smesh.summary()` estimates its working memory first and stops before allocation when
  the estimate exceeds 8 GB; the limit is set with `options(miMeta.max.stage2.gb = )`.

* `signature_table()` lists the selected signatures with their sharing category, active
  clusters, final effect and estimability in each cluster; `plot_signatures()` draws them as a
  heat map. A signature is all-cluster shared when it is selected in every cluster in which it is
  estimable, with at least two such clusters.

* Added the `sim_data` example data: 10 simulated contexts in two latent clusters, with 100
  features and 30 true signatures, used in the examples of the SMESH functions.

* Added a SMESH tutorial vignette, run on `sim_data` so that it builds in a few minutes.

## Breaking changes

* The count-matrix argument of `melody()` and `melody.null.model()` is renamed from `rel.abd` to
  `feature.counts`, and is used under that name by `smesh()` as well. It takes nonnegative
  sequencing read counts, which should not be converted to proportions or normalized. Code that
  passes `rel.abd =` must be updated.

## Dependencies

* `PALM` (>= 0.1.1) is now imported, for the per-study summary statistics `smesh()` generates.
  It is installed from <https://github.com/ZjpWei/PALM_package>;
  `devtools::install_github("ZjpWei/miMeta")` pulls it in automatically. PALM 0.1.1 names its
  count-matrix argument `feature.counts`, as miMeta does.

* Added `Matrix` to Imports.

# miMeta 1.1.0

* Melody: meta-analysis of microbiome association studies for discovering generalizable microbial
  signatures. `melody()`, `melody.null.model()`, `melody.get.summary()` and
  `melody.meta.summary()`.
