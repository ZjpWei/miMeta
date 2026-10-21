#' @title Signature table of a fitted SMESH model
#'
#' @description Summarizes the microbial signatures selected by SMESH, one row per feature: its
#' sharing category, the clusters in which it is selected, its final effect in each cluster, and
#' whether each cluster could estimate it.
#'
#' @author Zhoujingpeng Wei, Guanhua Chen, Zheng-Zheng Tang
#'
#' @param fit A fitted object returned by \code{\link{smesh}} or \code{\link{smesh.summary}}:
#' either the full list, with one element per covariate of interest, or one of its elements.
#' @param covariate The covariate of interest to summarize when \code{fit} holds several. Default
#' is NULL, which uses the first.
#' @param selected.only Whether to list only the selected signatures, that is, features with a
#' nonzero effect in at least one cluster. Default is TRUE; set it to FALSE to list every feature
#' in the fit.
#'
#' @return A data frame with one row per feature, in the order drawn by
#' \code{\link{plot_signatures}}: all-cluster-shared signatures first, by the sum of the cluster
#' effects; then subset-shared and then cluster-specific signatures, grouped by the pattern of
#' active clusters and, within one pattern, by the sum of the cluster effects.
#' The columns are:
#' \describe{
#'   \item{\code{feature}}{The microbial feature identifier.}
#'   \item{\code{sharing}}{"all-cluster shared" (selected in every cluster in which it is
#'   estimable), "subset-shared" (selected in more than one cluster but not in every cluster in
#'   which it is estimable), or "cluster-specific" (selected in exactly one cluster). With
#'   \code{selected.only = FALSE}, unselected features are "not selected".}
#'   \item{\code{active_clusters}}{The clusters in which the feature is selected, separated by
#'   "|", for example "C1|C3".}
#'   \item{\code{effect_C1}, \code{effect_C2}, ...}{The final effect estimate in each cluster, as
#'   in the \code{effects} component of the fit: nonzero where the feature is selected, zero where
#'   it is estimable but not selected, and \code{NA} where no context assigned to the cluster
#'   provided an estimable association summary for it.}
#'   \item{\code{estimable_C1}, \code{estimable_C2}, ...}{Whether at least one context assigned
#'   to the cluster provided an estimable association summary for the feature.}
#' }
#'
#' @details
#' The sharing category is decided by the clusters in which the effect is nonzero, among the
#' clusters in which the feature is estimable. A signature is all-cluster shared when it is
#' selected in every cluster in which it is estimable and at least two clusters can estimate it
#' (with a single cluster, every selected feature is all-cluster shared); its effect magnitudes
#' and directions may differ across clusters. A cluster in which no assigned context could
#' estimate the feature has an \code{NA} effect and is not counted; the \code{estimable_} columns
#' show when that is the case, and \code{active_clusters} then lists only the clusters that
#' selected the feature. A signature selected in exactly one cluster is cluster-specific.
#'
#' @seealso \code{\link{plot_signatures}},
#' \code{\link{smesh}},
#' \code{\link{smesh.summary}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' library("miMeta")
#' data("sim_data")
#' ## 10 simulated contexts in 2 latent clusters, 100 features, 30 true signatures
#' feature.counts <- sim_data$sim_abd
#' sim_meta <- sim_data$sim_meta
#'
#' ## Disease status of each sample, one single-column matrix per context
#' covariate.interest <- lapply(names(feature.counts), function(d){
#'   m <- sim_meta[sim_meta$Context == d, ]
#'   matrix(m$Group, ncol = 1, dimnames = list(m$Sample_ID, "disease"))
#' })
#' names(covariate.interest) <- names(feature.counts)
#'
#' ## Each context is one study, so `context` is left at NULL. G = 2, the true
#' ## number of clusters, keeps the example quick; with the default G = NULL,
#' ## G is selected by the minimum GIC.
#' fit <- smesh(feature.counts = feature.counts,
#'              covariate.interest = covariate.interest, G = 2)
#' sig <- signature_table(fit)
#' head(sig)
#' table(sig$sharing)
#' }
#'
signature_table <- function(fit, covariate = NULL, selected.only = TRUE) {

  obj <- get_smesh_component(fit, covariate)
  E <- obj$effects
  clusters <- colnames(E)
  G <- length(clusters)

  ## Sharing is judged among the clusters that can estimate the feature: a
  ## cluster with no data for it (NA) neither selects it nor counts against it.
  ## All-cluster shared needs at least two estimable clusters (one when G = 1),
  ## so a feature estimable in a single cluster is cluster-specific.
  active <- !is.na(E) & E != 0
  active <- matrix(active, nrow = nrow(E), dimnames = dimnames(E))
  n.active <- rowSums(active)
  n.estimable <- rowSums(!is.na(E))

  sharing <- ifelse(n.active == 0, "not selected",
             ifelse(n.active == n.estimable & n.estimable >= min(2, G), "all-cluster shared",
             ifelse(n.active == 1, "cluster-specific", "subset-shared")))
  pattern <- apply(active, 1, function(a) {
    if (any(a)) paste(clusters[a], collapse = "|") else NA_character_
  })

  ## Order: all-cluster shared, subset-shared, cluster-specific, not selected.
  ## The all-cluster-shared block is one group ordered by the sum of the effects,
  ## since its patterns differ only in which cluster had no data.  The others go
  ## by more active clusters first; among patterns of the same size, in the order
  ## of the cluster indices; within a pattern, by the sum of the effects.
  shared <- sharing == "all-cluster shared"
  rank <- match(sharing, c("all-cluster shared", "subset-shared", "cluster-specific",
                           "not selected"))
  key <- apply(active, 1, function(a) paste(sprintf("%04d", which(a)), collapse = "-"))
  key[shared] <- ""
  ord <- order(rank, -ifelse(shared, G, n.active), key, rowSums(E, na.rm = TRUE))

  effect <- as.data.frame(E)
  names(effect) <- paste0("effect_", clusters)
  estimable <- as.data.frame(matrix(!is.na(E), nrow = nrow(E)))
  names(estimable) <- paste0("estimable_", clusters)

  tab <- data.frame(feature = rownames(E), sharing = sharing, active_clusters = pattern,
                    stringsAsFactors = FALSE)
  tab <- cbind(tab, effect, estimable)[ord, , drop = FALSE]
  if (selected.only) {
    tab <- tab[tab$sharing != "not selected", , drop = FALSE]
  }
  rownames(tab) <- NULL
  tab
}


#' @title Heat map of the signatures selected by SMESH
#'
#' @description Draws the final cluster-level effect of every selected signature, one row per
#' signature and one column per cluster, grouped by sharing category.
#'
#' @author Zhoujingpeng Wei, Guanhua Chen, Zheng-Zheng Tang
#'
#' @param fit A fitted object returned by \code{\link{smesh}} or \code{\link{smesh.summary}}:
#' either the full list, with one element per covariate of interest, or one of its elements.
#' @param covariate The covariate of interest to plot when \code{fit} holds several. Default is
#' NULL, which plots the first.
#' @param title The plot title. Default is "Cluster-level association effects".
#'
#' @return A \code{ggplot} object. Rows follow the order of \code{\link{signature_table}}. Cells
#' where the feature is estimable but not selected (effect zero) are drawn in the neutral colour
#' of the scale; cells where no context assigned to the cluster could estimate the feature
#' (\code{NA}) are grey.
#'
#' @seealso \code{\link{signature_table}},
#' \code{\link{smesh}},
#' \code{\link{smesh.summary}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' library("miMeta")
#' data("sim_data")
#' ## 10 simulated contexts in 2 latent clusters, 100 features, 30 true signatures
#' feature.counts <- sim_data$sim_abd
#' sim_meta <- sim_data$sim_meta
#'
#' ## Disease status of each sample, one single-column matrix per context
#' covariate.interest <- lapply(names(feature.counts), function(d){
#'   m <- sim_meta[sim_meta$Context == d, ]
#'   matrix(m$Group, ncol = 1, dimnames = list(m$Sample_ID, "disease"))
#' })
#' names(covariate.interest) <- names(feature.counts)
#'
#' ## Each context is one study, so `context` is left at NULL. G = 2, the true
#' ## number of clusters, keeps the example quick; with the default G = NULL,
#' ## G is selected by the minimum GIC.
#' fit <- smesh(feature.counts = feature.counts,
#'              covariate.interest = covariate.interest, G = 2)
#' plot_signatures(fit)
#' }
#'
plot_signatures <- function(fit, covariate = NULL, title = "Cluster-level association effects") {

  tab <- signature_table(fit, covariate = covariate, selected.only = TRUE)
  if (nrow(tab) == 0) {
    stop("the fit selected no signatures.")
  }
  eff.cols <- grep("^effect_", names(tab), value = TRUE)
  clusters <- sub("^effect_", "", eff.cols)

  df <- data.frame(feature = rep(tab$feature, times = length(clusters)),
                   sharing = rep(tab$sharing, times = length(clusters)),
                   cluster = rep(clusters, each = nrow(tab)),
                   effect  = unlist(tab[eff.cols], use.names = FALSE),
                   stringsAsFactors = FALSE)
  df$feature <- factor(df$feature, levels = rev(tab$feature))     # table order, top down
  df$cluster <- factor(df$cluster, levels = clusters)
  df$sharing <- factor(df$sharing, levels = intersect(
    c("all-cluster shared", "subset-shared", "cluster-specific"), unique(tab$sharing)))
  mx <- max(abs(df$effect), na.rm = TRUE)

  ggplot(df, aes(x = cluster, y = feature, fill = effect)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_grid(rows = vars(sharing), scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#3B6FB6", mid = "#F7F7F7", high = "#B35836", midpoint = 0,
                         limits = c(-mx, mx), na.value = "grey85", name = "Effect") +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0))
}


## One covariate's component of a fitted SMESH object.  Accepts the full fit
## (one element per covariate of interest) or a single element of it.
get_smesh_component <- function(fit, covariate = NULL) {
  if (!is.null(fit$effects) && !is.null(fit$diagnostics)) {
    return(fit)
  }
  if (is.null(covariate)) covariate <- names(fit)[1]
  if (is.null(covariate) || !covariate %in% names(fit)) {
    stop("`fit` has no covariate of interest named '", covariate, "'.")
  }
  obj <- fit[[covariate]]
  if (is.null(obj$effects) || is.null(obj$diagnostics)) {
    stop("`fit` is not a fitted SMESH object returned by smesh() or smesh.summary().")
  }
  obj
}

## Columns referred to by name inside the ggplot2 calls above.
utils::globalVariables(c("cluster", "feature", "effect", "sharing"))
