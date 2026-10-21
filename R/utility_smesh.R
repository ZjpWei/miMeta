get_lasso_pre <- function(summary.stat.study,
                          G,
                          D_mat,
                          feature.ID,
                          feature.set,
                          study.ID){

  L <- length(study.ID)
  if(L == 1){
    tmp.taxa.mat <- feature.set[[study.ID]]
    taxa.mat <- matrix(tmp.taxa.mat, nrow = 1)
    colnames(taxa.mat) <- names(tmp.taxa.mat)
  }else{
    taxa.mat <- NULL
    for(d in study.ID){
      tmp.taxa.mat <- feature.set[[d]]
      taxa.mat <- rbind(taxa.mat, tmp.taxa.mat)
    }
  }
  rownames(taxa.mat) <- study.ID
  nonempty.id <- names(which(colSums(taxa.mat) > 0))
  empty.id <- names(which(colSums(taxa.mat) == 0))
  mu.len <- length(nonempty.id)
  taxa.mat <- taxa.mat[,nonempty.id]
  if(L == 1){
    taxa.mat <- matrix(taxa.mat, nrow = 1)
    colnames(taxa.mat) <- nonempty.id
    rownames(taxa.mat) <- study.ID
  }
  N <- sum(sapply(summary.stat.study, function(d){d$n}))
  k.l <- apply(taxa.mat, 1, sum)

  if(G == 1){
    D_sparse <- Matrix::Matrix(1, sparse = TRUE)
  }else{
    D_sparse <- Matrix::Matrix(t(D_mat), sparse = TRUE)
  }

  Y_list <- list()
  R_mat <- matrix(0, nrow = L, ncol = length(feature.ID), dimnames = list(study.ID, feature.ID))
  for (d in study.ID) {
    loc.k <- which(d == study.ID)
    Sigma.chol <- Matrix::Diagonal(
      x = sqrt(1 / diag(summary.stat.study[[d]]$cov))
    )

    ## Record R matrix
    R_mat[d, names(summary.stat.study[[d]]$est)] <- 1

    ## New Y_l = 1_G %x% weighted beta_hat_l
    y_l <- Sigma.chol %*% summary.stat.study[[d]]$est
    Y_list[[loc.k]] <- as.vector(rep(y_l, times = G))
  }
  Y.enlarge <- unlist(Y_list)

  ## Number of coefficients in the enlarged design, i.e. ncol(D %x% A_l).
  ## byabess_cluster rebuilds that design itself (with weight-scaled atoms), so only
  ## the count is needed here - it sets the upper end of the sparsity search.
  n.coef <- ncol(D_sparse) * mu.len

  lasso.mat <- list(Y.enlarge = Y.enlarge,
                    n.coef = n.coef,
                    R_mat = R_mat,
                    taxa.mat = taxa.mat,
                    k.l = k.l,
                    N = N,
                    mu.len = mu.len,
                    nonempty.id = nonempty.id,
                    empty.id = empty.id,
                    feature.ID = feature.ID,
                    study.ID = study.ID)

  return(lasso.mat)
}

# Optimization function
byabess_cluster <- function(summary.stat.study,
                       D_mat,
                       w_lg,
                       lasso.mat,
                       support.size,
                       tune.type = c("BIC", "HBIC", "KBIC", "EBIC"),
                       stage = c("1", "2")){

  G <- ncol(D_mat)
  M <- nrow(D_mat)

  ## customize the X matrix
  W_vec <- list()
  for (d in lasso.mat$study.ID) {
    loc.k <- which(d == lasso.mat$study.ID)

    ## W_lst
    W_vec[[loc.k]] <- rep(w_lg[d,], each = length(summary.stat.study[[d]]$est))
  }
  W_hat <- unlist(W_vec)

  X_list <- list()
  ## Each atom is scaled by the sqrt of the posterior mass of the clusters it
  ## spans.  DO NOT drop the 1e-5: if a cluster empties during EM its mass is
  ## exactly 0, and an unguarded division puts Inf into the design.  glmtlp
  ## accepted that silently; abess validates its input and stops with
  ## "x has missing value or infinite value!".  This mirrors the guard in
  ## SMESH_v4.R, which divided by (sqrt(pi_h[g]) + 1e-5).
  atom.mass <- sqrt(sapply(1:nrow(D_mat), function(d){
    sum(colMeans(w_lg[, D_mat[d, ] == 1, drop = FALSE]))
  }))
  D_sparse <- Matrix::Matrix(t(D_mat / (atom.mass + 1e-5)), sparse = TRUE)
  for (d in lasso.mat$study.ID) {
    loc.k <- which(d == lasso.mat$study.ID)
    Sigma.chol <- Matrix::Diagonal(
      x = sqrt(1 / diag(summary.stat.study[[d]]$cov))
    )
    mu.mat <- Matrix::Diagonal(lasso.mat$mu.len)
    id.mu <- lasso.mat$taxa.mat[d, ]
    mu.mat <- mu.mat[id.mu, , drop = FALSE]
    A_l <- Sigma.chol %*% mu.mat

    ## New X_l = D %x% A_l
    X_list[[loc.k]] <- kronecker(D_sparse, A_l)

  }
  X.enlarge <- as.matrix(do.call(rbind, X_list))

  ## optim model
  suppressMessages(
    result <- abess::abess(x = X.enlarge,
                           y = as.vector(lasso.mat$Y.enlarge),
                           weight = W_hat,
                           family = "gaussian",
                           tune.type = "bic",
                           tune.path = "sequence",
                           normalize = 0,
                           support.size = support.size,
                           fit.intercept = FALSE,
                           standardize = FALSE)

    # result <- glmtlp::glmtlp(X = X.enlarge,
    #                          y = as.vector(lasso.mat$Y.enlarge),
    #                          family = "gaussian",
    #                          penalty = "l0",
    #                          weights = W_hat,
    #                          kappa = support.size,
    #                          tol = 1e-3,      ## Default: 1e-4
    #                          dc.maxit = 10,   ## Default: 20
    #                          cd.maxit = 1000, ## Default: 10000
    #                          nr.maxit = 10,   ## Default: 20
    #                          standardize = FALSE)
  )

  ## Compute residuals and q_loss
  res_beta <- result$beta[,1]
  resid <- as.vector(lasso.mat$Y.enlarge) - X.enlarge %*% res_beta
  q_loss <- sum(W_hat * resid^2) / lasso.mat$N

  ## Re-construct estimates
  B_mat <- matrix(res_beta, nrow = lasso.mat$mu.len, ncol = M, dimnames = list(lasso.mat$nonempty.id, as.character(1:M)))
  theta_mat <- B_mat %*% D_mat

  ## Calculate ic
  if(stage == "1"){
    df <- sum(res_beta!=0) + (G - 1)
  }else{
    df <- sum(res_beta!=0)
  }

  if(tune.type == "BIC"){
    ic <- df * log(lasso.mat$N)/lasso.mat$N
  }else if(tune.type == "KBIC"){
    ic <- df * log(max(exp(1), log(length(res_beta)))) * log(lasso.mat$N) / lasso.mat$N
  }else if(tune.type == "HBIC"){
    ic <- df * log(length(res_beta)) * log(log(lasso.mat$N))/lasso.mat$N
  }else if (tune.type == "EBIC"){
    ic <- (df * log(lasso.mat$N) + (lfactorial(length(res_beta)) - lfactorial(df) - lfactorial(length(res_beta) - df))) / lasso.mat$N
  }
  return(list(theta_mat = theta_mat, mu = B_mat, res_beta = res_beta, q_loss = q_loss, ic = ic, GIC = q_loss + ic))
}

# Meta-analysis
meta_stage1 <- function(summary.stat.study,
                        lasso.mat,
                        tune.type,
                        G,
                        feature.ID,
                        tol,
                        NMAX,
                        nstart,
                        verbose){

  study.ID <- names(summary.stat.study)
  L <- length(study.ID)

  ## Stage-1 atom dictionary: one atom per cluster plus one all-cluster-shared
  ## atom.  This must match the D_mat that smesh.meta.summary() passed to
  ## get_lasso_pre() when it built `lasso.mat`, or byabess_cluster() will index the
  ## design matrix with the wrong number of atoms.
  ## At G > 1 the dictionary is one atom per cluster plus one shared atom.  At
  ## G = 1 those two collapse to the same column, so a single atom is used
  ## instead - which is also what get_lasso_pre() builds in its G == 1 branch.
  D_mat <- if (G == 1) matrix(1, 1, 1) else rbind(diag(G), rep(1, G))

  ## Progress reporting.  Each sparsity level runs nstart random starts of an EM loop,
  ## which is far too many lines to print one by one.  The counter is therefore
  ## rewritten in place with a carriage return, so a whole sparsity level
  ## occupies a single live line; when the level finishes that line is cleared
  ## and one permanent line is committed.
  ##
  ## Padded to a fixed width and cleared with spaces rather than an ANSI escape,
  ## so a redirected log picks up no control codes beyond the "\r" itself.
  PROG_W    <- 56
  prog_line <- function(s, sd, nstart, shift) {
    formatC(sprintf("      [s = %s | restart %d/%d] pi shift %s",
                    s, sd, nstart, signif(shift, 3)),
            width = -PROG_W)
  }
  prog_clear <- function() {
    cat("\r", strrep(" ", PROG_W), "\r", sep = "", file = stderr())
  }

  ## Start to tune
  loop.search <- TRUE
  if(G == 1){
    search.vec <- round( lasso.mat$n.coef * c(1, 2 / (1 + sqrt(5)), 1 - 2 / (1 + sqrt(5)), 0) / 2 )
  }else{
    search.vec <- round( lasso.mat$n.coef * c(1, 2 / (1 + sqrt(5)), 1 - 2 / (1 + sqrt(5)), 0) / (G+2) )
  }
  search.vec.loc <- rep(TRUE, 4)
  initial_sort <- list()
  gic_sort <- rep(NA, 4)

  while(loop.search){
    for(k in 1:sum(search.vec.loc)){
      s.lambda <- (search.vec[search.vec.loc])[k]

      if(L == 1){
        tune.deltas <- 1
        names(tune.deltas) <- study.ID
      }else{

        ## One entry per random start: its cluster assignment, the GIC at the end
        ## of its EM run, and whether EM met `tol`.
        cluster_list <- list()
        restart_gic <- rep(NA_real_, nstart)
        restart_converged <- rep(NA, nstart)
        for(sd in 1:nstart){
          s.lambda.local <- s.lambda

          ## Random initial cluster
          set.seed(sd)
          tune.deltas <- sample(1:G, size = L, replace = TRUE)
          names(tune.deltas) <- study.ID

          ## initial matrix
          beta_mat <- matrix(
            NA,
            nrow = length(feature.ID),
            ncol = L,
            dimnames = list(feature.ID, study.ID)
          )
          for (d in study.ID) {
            beta_mat[names(summary.stat.study[[d]]$est), d] <- summary.stat.study[[d]]$est
          }

          ## initialize mu_g and pi_g
          mu_g <- pi_g <- NULL
          for (l in 1:G) {
            mu_g <- cbind(
              mu_g,
              rowMeans(
                beta_mat[, names(tune.deltas)[tune.deltas == l], drop = FALSE],
                na.rm = TRUE
              )
            )
            pi_g <- c(pi_g, mean(tune.deltas == l))
          }
          colnames(mu_g) <- names(pi_g) <- as.character(1:G)

          loops <- TRUE
          loop.count <- 0
          converged <- FALSE

          while (loops) {
            ## E step
            log_dens <- matrix(
              NA, nrow = L, ncol = G,
              dimnames = list(colnames(beta_mat), as.character(1:G))
            )

            for (g in 1:G) {
              for (d in colnames(beta_mat)) {
                mu.tmp <- mu_g[!is.na(beta_mat[, d]), g]
                mu.tmp[is.na(mu.tmp)] <- 0
                log_dens[d, g] <- mvn_density(
                  x = beta_mat[!is.na(beta_mat[, d]), d],
                  mu = mu.tmp,
                  Sigma = summary.stat.study[[d]]$cov
                )
              }
            }

            ## apply() returns a bare vector when G == 1, so t() would give a
            ## 1 x L matrix instead of L x G.  Rebuild the shape explicitly.
            w_raw <- apply(log_dens, 1, function(d) {
              log_weighted_normalize(log_vals = d, weights = pi_g)
            })
            w_lg <- matrix(w_raw, nrow = L, ncol = G, byrow = TRUE,
                           dimnames = list(rownames(log_dens), as.character(1:G)))

            ## M step
            pi_h <- colMeans(w_lg)
            s.result <- byabess_cluster(
              summary.stat.study = summary.stat.study,
              D_mat = D_mat,
              w_lg = w_lg,
              lasso.mat = lasso.mat,
              support.size = s.lambda.local,
              tune.type = tune.type,
              stage = "1"
            )

            ## theta_mat is B_mat %*% D_mat: the cluster-level effects
            ## (features x G).  `mu` from byabess_cluster() is the basis coefficient
            ## matrix (features x M) and is not what the EM step needs.
            mu_h <- s.result$theta_mat

            ## `sd` indexes the random start (a different random initial cluster
            ## assignment); `loop.count` counts EM iterations within it.
            pi.shift <- sqrt(sum((pi_h - pi_g)^2, na.rm = TRUE))
            if (verbose) {
              cat("\r", prog_line(s.lambda.local, sd, nstart, pi.shift),
                  sep = "", file = stderr())
            }
            if (pi.shift <= tol) {
              #if (sum(s.result$res_beta != 0) == s.lambda.local) {
              loops <- FALSE
              converged <- TRUE
              pi_g <- pi_h
              mu_g <- mu_h
            } else if (loop.count > NMAX) {
              warning(sprintf(
                "EM algorithm doesn't converge with given tolerance in s = %s.",
                s.lambda.local
              ))
              loops <- FALSE
              pi_g <- pi_h
              mu_g <- mu_h
            } else {
              pi_g <- pi_h
              mu_g <- mu_h
              loop.count <- loop.count + 1
            }
          }

          cluster_list[[sd]] <-   apply(w_lg, 1, which.max)
          restart_gic[sd] <- s.result$GIC
          restart_converged[sd] <- converged
        }

        ## All restarts for this sparsity level are done: clear the live line
        ## and commit one permanent line showing where it ended.
        if (verbose) {
          prog_clear()
          message(prog_line(s.lambda.local, nstart, nstart, pi.shift))
        }

        ## The reported partition is the best converged start: the one whose EM
        ## reached the smallest GIC at this (G, s).  Random starts guard against
        ## local optima but are NOT averaged into the estimator.
        eligible <- which(restart_converged %in% TRUE)
        if (!length(eligible)) {
          ## No start met `tol` at this sparsity level; fall back to the best of
          ## the unconverged ones rather than failing, and say so.
          eligible <- seq_len(nstart)
          warning(sprintf(
            "No random start converged at s = %s; the smallest GIC among the unconverged starts is used.",
            s.lambda))
        }
        best_start <- eligible[which.min(restart_gic[eligible])]

        #=== Search best model by sequence ===#
        tune.deltas <- cluster_list[[best_start]]

        ## Same shape guard as above.  Width is G, not max(tune.deltas): cutree()
        ## returns labels 1..G so they agree, but G cannot silently under-count.
        w_raw <- sapply(tune.deltas, function(d){
          dd <- rep(0, G)
          dd[d] <- 1
          dd
        })
        w_lg <- matrix(w_raw, nrow = L, ncol = G, byrow = TRUE,
                       dimnames = list(names(tune.deltas), as.character(1:G)))

        ## Estimate:
        pi_h <- colMeans(w_lg)
        s.result <- byabess_cluster(summary.stat.study = summary.stat.study,
                               D_mat = D_mat,
                               w_lg = w_lg,
                               lasso.mat = lasso.mat,
                               support.size = s.lambda,
                               tune.type = tune.type,
                               stage = "1")

        mu_h <- s.result$theta_mat
      }

      ####
      ct <- which(search.vec == s.lambda)
      initial_sort[[ct]] <- list(pi_g = pi_h, mu_g = mu_h, w_lg = w_lg,
                                 B = s.result$mu,
                                 gic = s.result$GIC, ic = s.result$ic,
                                 q_loss = s.result$q_loss,
                                 cluster_list = cluster_list,
                                 restart_gic = restart_gic,
                                 restart_converged = restart_converged,
                                 selected_restart = best_start,
                                 s.lambda = s.lambda)

      gic_sort[ct] <- initial_sort[[ct]]$gic
    }

    ## gsection
    min.gic.id <- which.min(gic_sort)[1]
    if(min.gic.id == 1 | min.gic.id == 2){
      if(abs(search.vec[1] - search.vec[2]) <= 2){
        loop.search <- FALSE
        results.all <- list(mu = initial_sort[[2]]$mu_g,
                            W = initial_sort[[2]]$w_lg,
                            Pi = initial_sort[[2]]$pi_g,
                            gic = initial_sort[[2]]$gic,
                            ic = initial_sort[[2]]$ic,
                            q_loss = initial_sort[[2]]$q_loss,
                            cluster_list = initial_sort[[2]]$cluster_list,
                            B = initial_sort[[2]]$B,
                            restart_gic = initial_sort[[2]]$restart_gic,
                            restart_converged = initial_sort[[2]]$restart_converged,
                            selected_restart = initial_sort[[2]]$selected_restart,
                            ## Support size stage 1 settled on.  meta_stage2()
                            ## takes this as `s_start`, the top of its own
                            ## sparsity search - the role `s_tune` played in the
                            ## original cluster_tune() implementation.
                            s_start = initial_sort[[2]]$s.lambda)
      }else{
        search.vec[3:4] <- search.vec[2:3]
        search.vec[2] <- search.vec[1] + round((1 - 2 / (1 + sqrt(5))) * (search.vec[4] - search.vec[1]))
        initial_sort[[4]] <- initial_sort[[3]]
        initial_sort[[3]] <- initial_sort[[2]]
        initial_sort[[2]] <- NA
        gic_sort[3:4] <- gic_sort[2:3]
        gic_sort[2] <- NA
        search.vec.loc <- c(FALSE, TRUE, FALSE, FALSE)
      }
    }else{
      if(abs(search.vec[3] - search.vec[4]) <= 2){
        loop.search <- FALSE
        results.all <- list(mu = initial_sort[[3]]$mu_g,
                            W = initial_sort[[3]]$w_lg,
                            Pi = initial_sort[[3]]$pi_g,
                            gic = initial_sort[[3]]$gic,
                            ic = initial_sort[[3]]$ic,
                            q_loss = initial_sort[[3]]$q_loss,
                            cluster_list = initial_sort[[3]]$cluster_list,
                            B = initial_sort[[3]]$B,
                            restart_gic = initial_sort[[3]]$restart_gic,
                            restart_converged = initial_sort[[3]]$restart_converged,
                            selected_restart = initial_sort[[3]]$selected_restart,
                            ## Support size stage 1 settled on.  meta_stage2()
                            ## takes this as `s_start`, the top of its own
                            ## sparsity search - the role `s_tune` played in the
                            ## original cluster_tune() implementation.
                            s_start = initial_sort[[3]]$s.lambda)
      }else{
        search.vec[1:2] <- search.vec[2:3]
        search.vec[3] <- search.vec[1] + round((2 / (1 + sqrt(5))) * (search.vec[4] - search.vec[1]))
        initial_sort[[1]] <- initial_sort[[2]]
        initial_sort[[2]] <- initial_sort[[3]]
        initial_sort[[3]] <- NA
        gic_sort[1:2] <- gic_sort[2:3]
        gic_sort[3] <- NA
        search.vec.loc <- c(FALSE, FALSE, TRUE, FALSE)
      }
    }
  }
  return(results.all)
}

## Stage 2: the cluster assignment from stage 1 is held fixed (w_lg is a hard
## 0/1 membership matrix), so there is no EM loop here - only a golden-section
## search over the sparsity level on the 2^G - 1 subset basis.  That is why this
## function takes no tol/NMAX.
##
## `s_start` sets the upper end of that search: the support size stage 1
## selected, carried over in model_stage1$s_start.  This is the same quantity
## the original cluster_tune() returned as `s_tune`.
meta_stage2 <- function(summary.stat.study,
                        lasso.mat,
                        tune.type,
                        D_mat,
                        w_lg,
                        s_start){


  # Start to tune
  search.vec <- round(s_start * c(1, 2 / (1 + sqrt(5)), 1 - 2 / (1 + sqrt(5)), 0))
  search.vec.loc <- rep(TRUE, 4)

  ## Step 2: 2^G -1 model fit
  loop.search <- TRUE
  initial_sort <- list()
  gic_sort <- rep(NA, 4)
  while(loop.search){
    for(k in 1:sum(search.vec.loc)){
      s.lambda <- (search.vec[search.vec.loc])[k]

      ## M step:
      pi_h <- colMeans(w_lg)
      s.result <- byabess_cluster(summary.stat.study = summary.stat.study,
                             D_mat = D_mat,
                             w_lg = w_lg,
                             lasso.mat = lasso.mat,
                             support.size = s.lambda,
                             tune.type = tune.type,
                             stage = "2")

      theta_h <- s.result$theta_mat
      mu_h <- s.result$mu

      ## Record model:
      ct <- which(search.vec == s.lambda)
      initial_sort[[ct]] <- list(pi_g = pi_h,
                                 mu_g = mu_h,
                                 theta_g = theta_h,
                                 w_lg = w_lg,
                                 gic = s.result$GIC,
                                 ic = s.result$ic,
                                 q_loss = s.result$q_loss,
                                 s.lambda = s.lambda)

      gic_sort[ct] <- initial_sort[[ct]]$gic
    }

    ## Golden-section search step
    min.gic.id <- which.min(gic_sort)
    extract_result <- function(fit) {
      list(
        mu = fit$theta_g,
        B = fit$mu_g,
        W = fit$w_lg,
        Pi = fit$pi_g,
        gic = fit$gic,
        ic = fit$ic,
        q_loss = fit$q_loss,
        s.lambda = fit$s.lambda
      )
    }
    left.side <- min.gic.id %in% c(1, 2)

    if (left.side) {
      ## The minimum is on the left side of the current search interval.
      ## Shrink the interval from the right.
      if (abs(search.vec[1] - search.vec[2]) <= 2) {
        loop.search <- FALSE
        best.id <- which.min(gic_sort)
        results.all <- extract_result(initial_sort[[best.id]])
      } else {
        search.vec[3:4] <- search.vec[2:3]
        search.vec[2] <- search.vec[1] + round((1 - 2 / (1 + sqrt(5))) * (search.vec[4] - search.vec[1]))
        initial_sort[[4]] <- initial_sort[[3]]
        initial_sort[[3]] <- initial_sort[[2]]
        initial_sort[[2]] <- NA
        gic_sort[3:4] <- gic_sort[2:3]
        gic_sort[2] <- NA
        search.vec.loc <- c(FALSE, TRUE, FALSE, FALSE)
      }
    } else {
      ## The minimum is on the right side of the current search interval.
      ## Shrink the interval from the left.
      if (abs(search.vec[3] - search.vec[4]) <= 2) {
        loop.search <- FALSE
        best.id <- which.min(gic_sort)
        results.all <- extract_result(initial_sort[[best.id]])
      } else {
        search.vec[1:2] <- search.vec[2:3]
        search.vec[3] <- search.vec[1] + round((2 / (1 + sqrt(5))) * (search.vec[4] - search.vec[1]))
        initial_sort[[1]] <- initial_sort[[2]]
        initial_sort[[2]] <- initial_sort[[3]]
        initial_sort[[3]] <- NA
        gic_sort[1:2] <- gic_sort[2:3]
        gic_sort[3] <- NA
        search.vec.loc <- c(FALSE, FALSE, TRUE, FALSE)
      }
    }
  }

  return(results.all)
}
mvn_density <- function(x, mu, Sigma) {
  n <- length(mu)
  centered <- x - mu
  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus
  mahalanobis_sq <- sum(centered^2 / diag(Sigma))
  log_density <- -0.5 * (n * log(2 * pi) + log_det_Sigma + mahalanobis_sq)
  return(as.numeric(log_density))
}

log_weighted_normalize <- function(log_vals, weights) {
  if (length(log_vals) != length(weights)) {
    stop("log_vals and weights must have the same length")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    warning("Weights do not sum to 1; normalizing them.")
    weights <- weights / sum(weights)
  }
  log_weighted <- log(weights) + log_vals
  log_max <- max(log_weighted)
  stabilized <- log_weighted - log_max
  unnormalized <- exp(stabilized)
  normalized_weights <- unnormalized / sum(unnormalized)
  return(normalized_weights)
}

## Memory guard for the second stage.  Its design is dense, with G rows per
## context-feature estimate and one column per (nonempty cluster subset, feature),
## so it grows as G 2^G.  The working memory is taken as three times the design
## (the design plus the copies made while it is built and fitted); the fit stops
## before allocation when that exceeds getOption("miMeta.max.stage2.gb", 8) GB.
check_stage2_memory <- function(summary.stat.study, G) {
  max.gb <- getOption("miMeta.max.stage2.gb", 8)
  if (!is.numeric(max.gb) || length(max.gb) != 1 || is.na(max.gb) || max.gb <= 0) {
    stop("`options(miMeta.max.stage2.gb)` must be a single positive number of GB.", call. = FALSE)
  }
  n.rows <- G * sum(vapply(summary.stat.study, function(d) length(d$est), numeric(1)))
  n.cols <- (2^G - 1) * length(unique(unlist(lapply(summary.stat.study, function(d) names(d$est)))))
  design.gb <- 8 * n.rows * n.cols / 1024^3
  working.gb <- 3 * design.gb
  if (working.gb > max.gb) {
    stop(sprintf(paste0(
      "The full Stage 2 subset basis is estimated to require %.1f GB of working memory. ",
      "Fit a smaller G, reduce the feature set, or increase the configured memory limit ",
      "(options(miMeta.max.stage2.gb = ), currently %s GB)."), working.gb, format(max.gb)),
      call. = FALSE)
  }
  invisible(working.gb)
}

## Unpenalized refit of the cluster-level effects on the selected patterns.
##
## For feature k with selected atoms S (the rows of D_mat where B[k, ] != 0), the
## cluster effects are the weighted least-squares fit constrained to the span of
## those atoms:
##   theta_k = D_S' (D_S A_k D_S')^+ D_S c_k,
##   A_k = diag(sum_l w_lg / v_lk),   c_k = sum_l w_lg beta_lk / v_lk,
## with W the hard context-by-cluster membership.  This is the fit the selection
## step makes on its own support once the atom scaling byabess_cluster() uses
## is undone, so it carries no shrinkage.
refit_effects <- function(summary.stat.study, W, B, D_mat) {
  ctx <- rownames(W)
  feat <- rownames(B)
  G <- ncol(W)
  est <- sapply(ctx, function(d) summary.stat.study[[d]]$est[feat])
  var <- sapply(ctx, function(d) {
    v <- diag(summary.stat.study[[d]]$cov)
    names(v) <- rownames(summary.stat.study[[d]]$cov)
    v[feat]
  })
  est <- matrix(est, nrow = length(feat), dimnames = list(feat, ctx))
  var <- matrix(var, nrow = length(feat), dimnames = list(feat, ctx))
  prec <- ifelse(is.na(var), 0, 1 / var)
  y    <- ifelse(is.na(est), 0, est)
  A <- prec %*% W            # feature x cluster: sum_l w_lg / v_lk
  C <- (prec * y) %*% W      # feature x cluster: sum_l w_lg beta_lk / v_lk

  theta <- matrix(0, length(feat), G, dimnames = list(feat, colnames(W)))
  for (k in seq_along(feat)) {
    sel <- B[k, ] != 0
    if (!any(sel)) next
    Ds <- D_mat[sel, , drop = FALSE]
    theta[k, ] <- t(Ds) %*% MASS::ginv(Ds %*% diag(A[k, ], nrow = G) %*% t(Ds)) %*%
      Ds %*% C[k, ]
  }
  theta
}

## Posterior membership of each context under the fitted mixture: the E-step of
## the clustering, evaluated once at the reported effects and mixing proportions.
posterior_membership <- function(summary.stat.study, effects, mixing) {
  ctx <- names(summary.stat.study)
  G <- ncol(effects)
  post <- t(sapply(ctx, function(d) {
    x  <- summary.stat.study[[d]]$est
    ld <- sapply(seq_len(G), function(g) {
      mu <- effects[names(x), g]
      mu[is.na(mu)] <- 0
      mvn_density(x = x, mu = mu, Sigma = summary.stat.study[[d]]$cov)
    })
    log_weighted_normalize(log_vals = ld, weights = mixing)
  }))
  post <- matrix(post, nrow = length(ctx), dimnames = list(ctx, colnames(effects)))
  post
}
