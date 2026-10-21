#' @title Structure-learning meta-analysis of context-level summary statistics
#'
#' @description Fits the SMESH structure-learning model to user-supplied context-level summary
#' statistics.
#'
#' @author Zhoujingpeng Wei, Guanhua Chen, Zheng-Zheng Tang
#'
#' @param summary.stats A named list of context-level summary statistics with one element per
#' context; the list names are the context identifiers. Each element is a list with three
#' components: \code{est}, a feature-by-covariate matrix of association coefficient estimates;
#' \code{stderr}, the matching matrix of their standard errors; and \code{n}, the number of samples
#' contributing to the context. Row names of \code{est} and \code{stderr} are microbial feature
#' identifiers and column names are the covariates of interest. Feature identifiers must be
#' harmonized so that the same identifier denotes the same microbial feature in every context,
#' although a feature need not be estimable in every context: a feature that a context cannot
#' estimate is recorded as \code{NA}.
#' @param G,support.ic,tol,NMAX,nstart,verbose See function \code{\link{smesh}}.
#'
#' @return Same output as the function \code{\link{smesh}}. In particular, a feature is all-cluster
#' shared when its association is selected as nonzero in every cluster in which it is estimable; its
#' effects need not be equal and may differ in magnitude and direction across clusters.
#'
#' @details
#' SMESH first learns latent context clusters jointly with sparse association effects. Conditional
#' on the learned clustering, it estimates the final signature-sharing patterns and refits the
#' selected cluster-level effects without penalization. A signature is all-cluster shared when its
#' association is selected as nonzero in every cluster in which it is estimable; the estimated
#' effects need not be equal in magnitude across clusters.
#'
#' With one or two clusters, the Stage 1 basis already represents every nonempty subset of
#' clusters. With \eqn{G > 2}, Stage 2 holds the learned clustering fixed and considers all
#' \eqn{2^G - 1} nonempty cluster subsets. Because the number of candidate sharing patterns grows
#' exponentially with \eqn{G}, the function estimates the required memory before constructing the
#' Stage 2 design matrix. By default, it stops if the estimated working memory exceeds 8 GB. The
#' limit can be changed, for example, with \code{options(miMeta.max.stage2.gb = 16)}.
#'
#' Conditional on each candidate \eqn{G}, the sparse support size is selected using the
#' information criterion specified by \code{support.ic}. When \code{G = NULL}, the candidate
#' \eqn{G} with the smallest Stage 1 GIC is retained. SMESH uses multiple random starts and
#' reports the converged fit with the smallest GIC. The restart information stored in
#' \code{diagnostics} describes computational behavior across initializations.
#'
#' @seealso \code{\link{smesh}},
#' \code{\link{signature_table}},
#' \code{\link{plot_signatures}},
#' \code{\link{melody.meta.summary}}
#'
#' @importFrom MASS ginv
#' @importFrom stats runif setNames
#' @importFrom utils combn
#' @importFrom Matrix Diagonal
#' @import abess
#' @import dplyr
#' @import ggplot2
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
#' ########## Generate context-level summary statistics ##########
#' ## One summary per context, from PALM. Alternative summary statistics can be
#' ## used instead, as long as they give an estimate and a standard error per
#' ## feature. A context covered by several studies would first combine them by
#' ## fixed-effect meta-analysis, for example with PALM::palm.meta.summary().
#' null.obj <- PALM::palm.null.model(feature.counts = feature.counts)
#' summary.stats <- PALM::palm.get.summary(null.obj = null.obj,
#'                                         covariate.interest = covariate.interest,
#'                                         correct = "tune")
#'
#' ########## Structure-learning meta-analysis ##########
#' ## G = 2 keeps the example quick; with the default G = NULL, G is selected by
#' ## the minimum GIC.
#' fit <- smesh.summary(summary.stats = summary.stats, G = 2)
#' fit$disease$context_assignments
#' head(fit$disease$effects)
#' }
#'
smesh.summary <- function(summary.stats,
                          G = NULL,
                          support.ic = c("BIC", "HBIC", "KBIC", "EBIC"),
                          tol = 1e-3,
                          NMAX = 20,
                          nstart = 20,
                          verbose = FALSE) {

  support.ic <- match.arg(support.ic)

  if (length(summary.stats) < 2) {
    stop("`summary.stats` must hold at least two contexts.")
  }
  ## At least one cluster must hold two or more contexts; a partition of
  ## singletons carries no sharing structure to estimate.
  if (!is.null(G)) {
    if (G >= length(summary.stats)) {
      stop("`G` must be smaller than the number of contexts (", length(summary.stats),
           "); every cluster would hold a single context.")
    }
    if (G < 1) {
      stop("`G` must be at least 1.")
    }
  }
  if (length(nstart) != 1 || is.na(nstart) || nstart < 1 || nstart != round(nstart)) {
    stop("`nstart` must be a single positive integer.")
  }

  study.ID <- names(summary.stats)
  cov.int.ID <- sort(unique(unlist(
    lapply(summary.stats, function(d) colnames(d$est))
  )))

  feature.ID <- sort(unique(unlist(
    lapply(summary.stats, function(d) {
      rownames(d$est)[rowSums(!is.na(d$est)) > 0]
    })
  )))

  output.result <- vector("list", length(cov.int.ID))
  names(output.result) <- cov.int.ID

  ## Hard (maximum-probability) membership, as the fit uses it.  Built without
  ## apply(): at G = 1 apply() returns a bare vector and t() of it loses the
  ## context names.
  hard <- function(W) {
    H <- (W == apply(W, 1, max)) * 1
    dimnames(H) <- list(rownames(W), as.character(seq_len(ncol(W))))
    H
  }

  for (cov.name in cov.int.ID) {
    if (verbose) {
      message("++ Search for the best model for covariate of interest: ", cov.name, " ++")
    }
    summary.stat.study <- list()
    for (d in study.ID) {
      est.mat <- summary.stats[[d]]$est
      se.mat  <- summary.stats[[d]]$stderr

      if (!cov.name %in% colnames(est.mat)) {
        next
      }
      nonna.id <- intersect(feature.ID, rownames(est.mat)[!is.na(est.mat[, cov.name])])
      if (length(nonna.id) == 0) {
        next
      }

      summary.stat.study[[d]] <- list(
        est = est.mat[nonna.id, cov.name],
        cov = diag(se.mat[nonna.id, cov.name]^2),
        n   = summary.stats[[d]]$n
      )
      colnames(summary.stat.study[[d]]$cov) = rownames(summary.stat.study[[d]]$cov) <- nonna.id
    }

    if (length(summary.stat.study) == 0) {
      warning("No valid study found for covariate: ", cov.name)
      next
    }

    ## ---- how many clusters to try ------------------------------------------
    L <- length(summary.stat.study)
    if (is.null(G)) {
      ## Candidate G runs from 1 to round(2L/3); a partition always keeps fewer
      ## clusters than contexts.
      G.grid <- seq_len(max(1, min(L - 1, round(2 * L / 3))))
    } else {
      if (length(G) != 1 || is.na(G) || G != round(G) || G < 1) {
        stop("G must be NULL or a single positive integer.")
      }
      if (G >= L) {
        stop("G = ", G, " is not smaller than the number of contexts available for '",
             cov.name, "' (", L, ").")
      }
      G.grid <- G
      ## A supplied G fixes the second-stage design now, so its memory is checked
      ## before any time is spent on the clustering.
      if (G > 2) check_stage2_memory(summary.stat.study, G)
    }

    ## ---- cluster learning, over the grid ------------------------------------
    ## Every candidate G is fitted, and the one with the minimum GIC is kept,
    ## with ties resolved in favor of the smaller G.  With an explicit G the
    ## grid has one element and the loop simply fits it.
    fits <- list()
    gic.path <- NULL
    for (G.try in G.grid) {
      if(verbose){
        message("++ Learning context clusters for G = ", G.try, ". ++")
      }

      ## The dictionary must match the one meta_stage1() builds: a single atom
      ## at G = 1, otherwise one atom per cluster plus a shared one.
      lasso.mat <- get_lasso_pre(summary.stat.study = summary.stat.study,
                                 G = G.try,
                                 D_mat = if (G.try == 1) matrix(1, 1, 1)
                                         else rbind(diag(G.try), rep(1, G.try)),
                                 feature.ID = feature.ID,
                                 feature.set = lapply(summary.stat.study, function(d){
                                   taxa.vec.tmp <- rep(FALSE, length(feature.ID))
                                   names(taxa.vec.tmp) <- feature.ID
                                   taxa.vec.tmp[names(d$est)[!is.na(d$est)]] <- TRUE
                                   return(taxa.vec.tmp)
                                 }),
                                 study.ID = study.ID)

      model_stage1 <- meta_stage1(
        summary.stat.study = summary.stat.study,
        lasso.mat = lasso.mat,
        G = G.try,
        feature.ID = feature.ID,
        tune.type = support.ic,
        NMAX = NMAX,
        tol = tol,
        nstart = nstart,
        verbose = verbose
      )
      gic.path <- rbind(gic.path, data.frame(G = G.try, gic = model_stage1$gic))
      fits[[as.character(G.try)]] <- model_stage1
    }

    G0 <- gic.path$G[order(gic.path$gic, gic.path$G)[1]]
    model0 <- fits[[as.character(G0)]]
    rm(fits)
    if (verbose && length(G.grid) > 1) {
      message("++ G = ", G0, " selected: the minimum GIC over G = ",
              min(G.grid), ", ..., ", max(G.grid), ". ++")
    }

    ## ---- final sharing patterns ---------------------------------------------
    ## With one or two clusters the dictionary above already holds every subset
    ## of clusters, so its selection is final.  With more, the patterns are
    ## re-selected against every subset, holding the partition fixed.
    if(G0 > 2){
      check_stage2_memory(summary.stat.study, G0)
      if(verbose){
        message("++ Estimating signature-sharing patterns across all cluster subsets. ++")
      }

      generate_D_mat <- function(G) {
        D_mat <- do.call(rbind, lapply(seq_len(G), function(k) {
          combs <- combn(G, k)
          t(apply(combs, 2, function(idx) {
            x <- rep(0, G)
            x[idx] <- 1
            x
          }))
        }))
        rownames(D_mat) <- NULL
        colnames(D_mat) <- paste0("g", seq_len(G))
        return(D_mat)
      }

      lasso.mat <- get_lasso_pre(summary.stat.study = summary.stat.study,
                                 G = G0,
                                 D_mat = generate_D_mat(G0),
                                 feature.ID = feature.ID,
                                 feature.set = lapply(summary.stat.study, function(d){
                                   taxa.vec.tmp <- rep(FALSE, length(feature.ID))
                                   names(taxa.vec.tmp) <- feature.ID
                                   taxa.vec.tmp[names(d$est)[!is.na(d$est)]] <- TRUE
                                   return(taxa.vec.tmp)
                                 }),
                                 study.ID = study.ID)

      model_stage2 <- meta_stage2(
        summary.stat.study = summary.stat.study,
        lasso.mat = lasso.mat,
        tune.type = support.ic,
        D_mat = generate_D_mat(G0),
        w_lg = hard(model0$W),
        s_start = model0$s_start
      )

      support <- model_stage2$B
      D.final <- generate_D_mat(G0)
      W.final <- hard(model_stage2$W)
    } else {
      support <- model0$B
      D.final <- if (G0 == 1) matrix(1, 1, 1) else rbind(diag(G0), rep(1, G0))
      W.final <- hard(model0$W)
    }

    ## ---- unpenalized refit on the selected patterns ---------------------------
    if(verbose){
      message("++ Refitting the selected effects without penalization. ++")
    }
    W.final <- W.final[names(summary.stat.study), , drop = FALSE]
    clusters <- paste0("C", seq_len(G0))
    colnames(W.final) <- clusters

    effects <- refit_effects(summary.stat.study = summary.stat.study,
                             W = W.final, B = support, D_mat = D.final)
    colnames(effects) <- clusters

    mixing <- setNames(colMeans(W.final), clusters)
    posterior <- posterior_membership(summary.stat.study = summary.stat.study,
                                      effects = effects, mixing = mixing)

    ## NA where no context assigned to the cluster could estimate the feature.
    estimable <- sapply(clusters, function(g) {
      members <- rownames(W.final)[W.final[, g] == 1]
      rownames(effects) %in% unique(unlist(lapply(summary.stat.study[members],
                                                  function(d) names(d$est))))
    })
    effects[!matrix(estimable, nrow = nrow(effects))] <- NA_real_

    assigned <- apply(W.final, 1, which.max)
    assignments <- data.frame(
      context = rownames(W.final),
      cluster = clusters[assigned],
      posterior_probability = posterior[cbind(seq_along(assigned), assigned)],
      row.names = NULL, stringsAsFactors = FALSE)

    output.result[[cov.name]] <- list(
      effects              = effects,
      context_assignments  = assignments,
      posterior_membership = posterior,
      mixing_proportions   = mixing,
      gic                  = model0$gic,
      diagnostics = list(
        restart_gic         = model0$restart_gic,
        restart_convergence = model0$restart_converged,
        selected_restart    = model0$selected_restart,
        gic_path = if (is.null(G)) transform(gic.path, selected = G == G0) else NULL
      )
    )
  }

  return(output.result)
}
