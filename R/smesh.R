#' @title Structure-learning meta-analysis for discovering shared and context-dependent microbial signatures
#'
#' @description Runs the end-to-end SMESH pipeline: generates study-level association estimates and
#' standard errors using PALM, combines studies representing the same context, learns latent context
#' clusters, and estimates all-cluster-shared and context-dependent microbial signatures.
#'
#' @author Zhoujingpeng Wei, Guanhua Chen, Zheng-Zheng Tang
#'
#' @param feature.counts A named list of sample-by-feature matrices of nonnegative sequencing read
#' counts, with one matrix per study. The list names must be the study identifiers; matrix row names
#' and column names must be sample and microbial-feature identifiers, respectively. Row sums are
#' used as sequencing-depth measures. The counts should not be converted to proportions or
#' normalized. Missing values are not allowed.
#' @param covariate.interest A named list of matrices for single or multiple covariates of interest,
#' with one element per study named by the study identifier. In a matrix of a given study, each row
#' is a sample and each column is a covariate of interest and can only be a numeric vector. In a
#' given study, the order of samples should be matched with the order in "feature.counts". No
#' missing value is allowed.
#' @param context A vector that maps studies to contexts, with study identifiers as names and
#' context labels as values. Studies sharing a context label are combined by fixed-effect
#' meta-analysis into one context-level summary before the clusters are learned, which is how
#' several studies of the same disease or tumor type are combined. Default is NULL, in which case
#' every study is its own context and no combining is needed.
#' @param covariate.adjust A named list of data frames for covariates that will be adjusted for in
#' the model, with one element per study named by the study identifier. In a data frame of a given
#' study, each row is a sample and each column is a covariate to be adjusted and can only be a
#' factor or numeric vector. In a given study, the order of samples should be matched with the order
#' in "feature.counts". No missing value is allowed. Default is NULL (not adjusting for any
#' covariates).
#' @param cluster A named list of vectors that define the sample clusters for studies with
#' correlated samples, with one element per study that has correlated samples, named by the study
#' identifier. The order of samples should be matched with the order in "feature.counts" in the
#' corresponding study. For example, the values of this variable are subject IDs if each subject has
#' multiple correlated samples. Default is NULL (all samples in all studies are independent).
#' @param depth.filter A cutoff value to remove samples with sequencing depth less than or equal to
#' the cutoff. Default is 0.
#' @param prev.filter A cutoff value to remove microbial features with prevalence (proportion of
#' nonzero observations) less than or equal to the cutoff. This cutoff is applied to each study, so
#' a feature could be removed in a subset of the studies. The cutoff value must be in the range of
#' 0-1. Default is 0.1.
#' @param correct The bias correction used by PALM when generating the per-study summary
#' statistics. Available options are "tune" and "median". Default is "tune"; "median" is faster and
#' can be used for exploratory runs.
#' @param min.context The minimum number of contexts that must be able to estimate a feature for that
#' feature to enter the fit. Default is 1, which keeps every feature that at least one context can
#' estimate. Raising it restricts the analysis to features observed more widely.
#' @param G The number of context clusters. Default is NULL, in which case G is selected by the
#' minimum generalized information criterion (GIC): for each candidate G from 1 to
#' \eqn{\min(L - 1, \mathrm{round}(2L/3))}, where \eqn{L} is the number of contexts, the sparse
#' association model is first selected with \code{support.ic}, and the G whose selected model has
#' the minimum GIC is kept, with ties resolved in favor of the smaller G. If a number is supplied,
#' that G is fitted directly; it must be smaller than the number of contexts.
#' @param support.ic The information criterion used to select the sparse association model, that
#' is, which cluster-level coefficients are nonzero, conditional on a candidate G. Available options
#' are "BIC", "HBIC", "KBIC" and "EBIC". Default is "BIC". Selection of G is performed by the outer
#' GIC procedure described under \code{G}; the penalty in that GIC is also computed with
#' \code{support.ic}.
#' @param tol Convergence tolerance for the mixture fit. Default is 1e-3.
#' @param NMAX The maximum number of iterations in the mixture fit. Default is 20.
#' @param nstart The number of random initializations of the context clustering. Default is 20.
#' SMESH uses multiple random starts to reduce sensitivity to local optima; see the Details of
#' \code{\link{smesh.summary}}.
#' @param verbose Whether to print the progress of the analysis. Default is FALSE.
#'
#' @return
#' A named list with one element for each covariate of interest. Each element
#' contains the following components:
#'
#' \describe{
#'
#'   \item{\code{effects}}{
#'   A numeric matrix containing the final cluster-level association effect
#'   estimates. Rows correspond to microbial features and columns correspond
#'   to the learned context clusters. The estimates are obtained from the
#'   post-selection refit. A nonzero entry indicates that the
#'   feature was selected as associated in that cluster, whereas zero indicates
#'   that it was not selected in a cluster for which the feature was estimable.
#'   \code{NA} indicates that no context assigned to the cluster provided an
#'   estimable association summary for that feature. A feature with at least
#'   one nonzero entry is a selected microbial signature, and its pattern of
#'   nonzero entries determines whether it is all-cluster shared,
#'   subset-shared, or cluster-specific. A feature is all-cluster shared when its
#'   association is selected as nonzero in every cluster in which it is estimable;
#'   its effects need not be equal and may differ in magnitude and direction across
#'   clusters.
#'   \code{\link{signature_table}} lists the selected signatures with their sharing
#'   categories.
#'   }
#'
#'   \item{\code{context_assignments}}{
#'   A data frame with one row per context and the columns \code{context},
#'   \code{cluster}, and \code{posterior_probability}. The \code{cluster}
#'   column gives the cluster each context is assigned to in the reported fit,
#'   and \code{posterior_probability} gives the posterior membership probability
#'   of that cluster under the fitted mixture model. Cluster labels are
#'   consistent across all components of the fitted object but are otherwise
#'   arbitrary up to permutation.
#'   }
#'
#'   \item{\code{posterior_membership}}{
#'   A numeric context-by-cluster matrix of posterior membership probabilities
#'   from the finite-mixture model, evaluated at the reported fit. Rows
#'   correspond to contexts and columns to clusters, with each row summing to
#'   one. This matrix provides the complete posterior membership distribution,
#'   whereas \code{context_assignments} reports only the assigned cluster and
#'   its probability.
#'   }
#'
#'   \item{\code{mixing_proportions}}{
#'   A named numeric vector containing the estimated mixture proportion for
#'   each cluster. The values are nonnegative and sum to one. They represent
#'   the estimated prevalence of the latent clusters in the collection of
#'   contexts.
#'   }
#'
#'   \item{\code{gic}}{
#'   The generalized information criterion for the reported SMESH fit. When
#'   \code{G} is supplied by the user, this is the GIC for that fitted cluster
#'   resolution. When \code{G = NULL}, it is the GIC associated with the
#'   selected number of clusters.
#'   }
#'
#'   \item{\code{diagnostics}}{
#'   A named list of computational diagnostics for model fitting and
#'   cluster-number selection. The restart entries refer to the selected
#'   number of clusters and support size:
#'   \describe{
#'     \item{\code{restart_gic}}{
#'     The generalized information criterion (GIC) reached by each random
#'     start at the end of its EM run.
#'     }
#'     \item{\code{restart_convergence}}{
#'     A logical vector indicating, for each random start, whether the
#'     mixture fit met the convergence tolerance within \code{NMAX}
#'     iterations.
#'     }
#'     \item{\code{selected_restart}}{
#'     The index of the random start retained for the reported fit: the
#'     converged start with the smallest GIC. If no start converged, the start
#'     with the smallest GIC is retained, with a warning.
#'     }
#'     \item{\code{gic_path}}{
#'     When \code{G} is selected automatically, a data frame containing the
#'     candidate values of \code{G} that were fitted, their GIC values, and an
#'     indicator of the selected model. \code{NULL} when \code{G} is supplied.
#'     }
#'   }
#'   These quantities assess computational behavior across initializations of
#'   the same data and should not be interpreted as robustness to changes in
#'   the observed data.
#'   }
#' }
#'
#' @details
#' The function runs three steps. First it fits a PALM null model on the read counts and generates
#' association summary statistics for every study. Second, if "context" is supplied, the studies
#' belonging to one context are combined by a fixed-effect meta-analysis into a single
#' context-level estimate and standard error. Each context therefore contributes one association
#' profile to the clustering, while evidence from multiple studies is reflected through increased
#' precision. Third, the context-level
#' summaries are passed to \code{\link{smesh.summary}}, which learns the clusters and the
#' signatures.
#'
#' To reuse summary statistics you already have, or to produce them with a different association
#' model, call \code{\link{smesh.summary}} directly.
#'
#' @seealso \code{\link{smesh.summary}},
#' \code{\link{signature_table}},
#' \code{\link{plot_signatures}},
#' \code{\link{melody}}
#'
#' @import PALM
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
#' ########## Structure-learning meta-analysis ##########
#' ## Each context is one study, so `context` is left at NULL. G = 2, the true
#' ## number of clusters, keeps the example quick; with the default G = NULL,
#' ## G is selected by the minimum GIC.
#' fit <- smesh(feature.counts = feature.counts,
#'              covariate.interest = covariate.interest, G = 2)
#' fit$disease$context_assignments
#' head(fit$disease$effects)
#'
#' ## The truth the data were simulated from
#' sim_data$sim_truth$cluster
#' table(sim_data$sim_truth$signatures$sharing)
#' }
#'

smesh <- function(feature.counts,
                  covariate.interest,
                  context = NULL,
                  covariate.adjust = NULL,
                  cluster = NULL,
                  depth.filter = 0,
                  prev.filter = 0.1,
                  correct = c("tune", "median"),
                  min.context = 1,
                  G = NULL,
                  support.ic = c("BIC", "HBIC", "KBIC", "EBIC"),
                  tol = 1e-3,
                  NMAX = 20,
                  nstart = 20,
                  verbose = FALSE
) {

  correct   <- match.arg(correct)
  support.ic <- match.arg(support.ic)

  #=== Check the study-to-context map ===#
  study.ID <- names(feature.counts)
  if (is.null(context)) {
    ## Every study is its own context; nothing to combine.
    context <- stats::setNames(study.ID, study.ID)
  } else {
    if (is.null(names(context))) {
      stop("`context` must be named by study identifier.")
    }
    missing.ID <- setdiff(study.ID, names(context))
    if (length(missing.ID) > 0) {
      stop("`context` has no entry for the following studies in `feature.counts`: ",
           paste(missing.ID, collapse = ", "))
    }
    context <- context[study.ID]
    if (anyNA(context)) {
      stop("`context` has missing values for: ",
           paste(study.ID[is.na(context)], collapse = ", "))
    }
  }
  context.ID <- unique(as.character(context))

  ## Checked here as well as in smesh.summary(), so a G that cannot be fitted is
  ## caught before the summary statistics are generated rather than after.
  if (!is.null(G) && G >= length(context.ID)) {
    stop("`G` must be smaller than the number of contexts (", length(context.ID),
         "); every cluster would hold a single context.")
  }

  #=== Generate per-study summary statistics ===#
  if (verbose) {
    message("++ Generating summary statistics for ", length(study.ID), " studies. ++")
  }

  null.obj <- PALM::palm.null.model(feature.counts = feature.counts,
                                    covariate.adjust = covariate.adjust,
                                    depth.filter = depth.filter,
                                    prev.filter = prev.filter)

  summary.stats <- PALM::palm.get.summary(null.obj = null.obj,
                                          covariate.interest = covariate.interest,
                                          cluster = cluster,
                                          correct = correct)

  #=== Combine the studies of each context by fixed-effect meta-analysis ===#
  if (length(context.ID) == length(study.ID) &&
      all(sort(context.ID) == sort(study.ID))) {
    context.stats <- summary.stats
  } else {
    if (verbose) {
      message("++ Combining ", length(study.ID), " studies into ",
              length(context.ID), " contexts by fixed-effect meta-analysis. ++")
    }
    context.stats <- list()
    for (ctx in context.ID) {
      keep <- names(summary.stats)[context[names(summary.stats)] == ctx]
      res  <- PALM::palm.meta.summary(summary.stats = summary.stats[keep])

      ## `res` is one data frame per covariate of interest, with columns
      ## feature / coef / stderr.  Lay them out as the feature by covariate
      ## matrices smesh.summary() expects, NA where a covariate has no estimate.
      cov.int.ID <- names(res)
      feature.ID <- Reduce(union, lapply(res, function(x) x$feature))
      est <- stderr <- matrix(
        NA_real_, length(feature.ID), length(cov.int.ID),
        dimnames = list(feature.ID, cov.int.ID))
      for (cov.int in cov.int.ID) {
        est[res[[cov.int]]$feature, cov.int]    <- res[[cov.int]]$coef
        stderr[res[[cov.int]]$feature, cov.int] <- res[[cov.int]]$stderr
      }

      context.stats[[ctx]] <- list(
        est    = est,
        stderr = stderr,
        n      = sum(vapply(summary.stats[keep], function(d) d$n, numeric(1))))
    }
  }

  #=== Keep features estimable in at least `min.context` contexts ===#
  if (min.context > 1) {
    feature.ID <- Reduce(union, lapply(context.stats, function(d) rownames(d$est)))
    n.estimable <- vapply(feature.ID, function(k) {
      sum(vapply(context.stats, function(d)
        k %in% rownames(d$est) && any(!is.na(d$est[k, ])), logical(1)))
    }, numeric(1))
    keep.ID <- feature.ID[n.estimable >= min.context]
    if (length(keep.ID) == 0) {
      stop("no microbial feature is estimable in at least ", min.context,
           " contexts; lower `min.context`.")
    }
    if (verbose) {
      message("++ ", length(keep.ID), " of ", length(feature.ID),
              " features are estimable in at least ", min.context, " contexts. ++")
    }
    context.stats <- lapply(context.stats, function(d) {
      id <- intersect(keep.ID, rownames(d$est))
      d$est    <- d$est[id, , drop = FALSE]
      d$stderr <- d$stderr[id, , drop = FALSE]
      d
    })
  }

  #=== Structure-learning meta-analysis ===#
  SMESH.model <- smesh.summary(summary.stats = context.stats,
                               G = G,
                               support.ic = support.ic,
                               tol = tol,
                               NMAX = NMAX,
                               nstart = nstart,
                               verbose = verbose)

  return(SMESH.model)
}
