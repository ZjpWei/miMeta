#' @import abess

######## Utility file for meta-analysis ########
Get_lasso_pre = function(summary.stat.study,
                         feature.ID,
                         feature.set,
                         study.ID,
                         ref,
                         K){

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
  N <- 0
  k.l <- NULL
  for(d in study.ID){
    N <- N + summary.stat.study[[d]]$n
    k.l <- c(k.l, sum(taxa.mat[d,]))
  }
  for(d in study.ID){
    # Generate design matrix X and Y
    Sigma.chol <- chol(solve(summary.stat.study[[d]]$cov)/N)
    F.mat <- diag(mu.len)
    id.F <- taxa.mat[d,]
    F.mat <- F.mat[id.F,]
    if(d == study.ID[1]){
      X.col.1 <- Sigma.chol %*% F.mat
      Y.enlarge.1 <- Sigma.chol %*% (summary.stat.study[[d]]$est)
      Y.enlarge.2 <- Sigma.chol %*% rep(1, sum(taxa.mat[d,]))
    }else{
      X.col.1 <- rbind(X.col.1, Sigma.chol %*% F.mat)
      Y.enlarge.1 <- c(Y.enlarge.1, Sigma.chol %*% (summary.stat.study[[d]]$est) )
      Y.enlarge.2 <- c(Y.enlarge.2, Sigma.chol %*% rep(1, sum(taxa.mat[d,])) )
    }
  }
  colnames(X.col.1) <- nonempty.id

  lasso.mat <- list(X.enlarge = X.col.1, Y.enlarge.1 = Y.enlarge.1, Y.enlarge.2 = Y.enlarge.2,
                    taxa.mat = taxa.mat, k.l = k.l, N = N, mu.len = mu.len,
                    nonempty.id = nonempty.id, empty.id = empty.id,
                    feature.ID = feature.ID, feature.set = feature.set, study.ID = study.ID, ref = ref)

  return(lasso.mat)
}

# Optimization function
byabess_AA <- function(summary.stat.study,
                       lasso.mat,
                       delta,
                       support.size){
  #=== Load data ===#
  feature.set <- lasso.mat$feature.set
  feature.ID <- lasso.mat$feature.ID
  study.ID <- lasso.mat$study.ID
  X.enlarge <- lasso.mat$X.enlarge
  Y.enlarge.1 <- lasso.mat$Y.enlarge.1
  Y.enlarge.2 <- lasso.mat$Y.enlarge.2
  ref <- lasso.mat$ref
  N <- lasso.mat$N
  idx.lst <- lasso.mat$idx.lst
  taxa.mat <- lasso.mat$taxa.mat
  k.l <- lasso.mat$k.l
  mu.len <- lasso.mat$mu.len
  nonempty.id <- lasso.mat$nonempty.id
  empty.id <- lasso.mat$empty.id
  L <- length(study.ID)

  ref.shift <- c()
  if(length(delta) == 1){
    for(l in 1:L){
      ref.shift <- c(ref.shift, rep(delta, sum(feature.set[[l]])))
    }
  }else{
    for(l in 1:L){
      ref.shift <- c(ref.shift, rep(delta[l], sum(feature.set[[l]])))
    }
  }
  Y.enlarge <- Y.enlarge.1 + Y.enlarge.2 * ref.shift

  #=== Abess model ===#
  result <- abess(x = X.enlarge,
                  y = Y.enlarge,
                  family = "gaussian",
                  tune.type = "bic",
                  tune.path = "sequence",
                  normalize = 0,
                  fit.intercept = FALSE,
                  support.size = support.size)

  stopifnot(mu.len == length(result$beta))
  fit.coef <- result$beta[1:mu.len]
  q_loss <- sum( (Y.enlarge - X.enlarge %*% result$beta)^2 )
  zero.fit <- rep(0, length(empty.id))
  names(zero.fit) <- empty.id
  nonzero.fit <- fit.coef
  names(nonzero.fit) <- nonempty.id
  mu.fit <- c(nonzero.fit, zero.fit)[feature.ID]
  names(delta) <- ref
  if(length(unique(ref)) == 1){
    mu.fit <- mu.fit[setdiff(feature.ID, ref[1])]
  }
  return(list(mu.fit = mu.fit, delta = delta, q_loss = q_loss, N = N))
}

# Calculate GIC
GIC.cal <- function(result,
                    tune.type = c("BIC", "HBIC", "KBIC", "EBIC")
){
  tune.type <- match.arg(tune.type)
  df <- sum(result$mu.fit!=0)
  K <- length(result$mu.fit)
  L <- length(result$delta)
  P <- K
  if(tune.type == "BIC"){
    return(result$q_loss + df * log(result$N)/result$N )
  }else if(tune.type == "KBIC"){
    return(result$q_loss + df * log(max(exp(1), log(length(result$mu.fit)))) * log(result$N)/result$N )
  }else if(tune.type == "HBIC"){
    return(result$q_loss + df * log(length(result$mu.fit)) * log(log(result$N))/result$N)
  }else if (tune.type == "EBIC"){
    return(result$q_loss + (df * log(result$N) + (lfactorial(P) - lfactorial(df) - lfactorial(P-df))) / result$N)
  }
}

# Search delta
# search.ref.loc.gsec <- function(summary.stat.study,
#                                 lasso.mat,
#                                 study.l,
#                                 quantile_sm,
#                                 search.prop,
#                                 lambda,
#                                 result,
#                                 initial.range = 0.1,
#                                 tune.type,
#                                 tol = 1e-3,
#                                 verbose = FALSE){
#
#   feature.set <- lasso.mat$feature.set
#   study.ID <- lasso.mat$study.ID
#   feature.ID <- lasso.mat$feature.ID
#   ref <- lasso.mat$ref
#   L <- length(study.ID)
#   p.pi <- 2 / (sqrt(5) + 1)
#   if(initial.range > 1){
#     stop("Searching range should be smaller than 1. \n")
#   }
#   if(is.null(result)){
#     s.1.prop <- rep(1/2, L)
#     s.2.prop <- rep(1/2, L)
#     s.left.prop <- rep(1/2, L)
#     s.right.prop <- rep(1/2, L)
#     ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.1.prop)
#     initial.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                  delta = - ref.quant, support.size = lambda)
#     GIC.result <- GIC.cal(result = initial.result, tune.type = tune.type)
#   }else{
#     s.1.prop <- search.prop
#     s.2.prop <- search.prop
#     s.left.prop <- search.prop
#     s.right.prop <- search.prop
#     GIC.result <- result$GIC.result
#   }
#   loop.range <- TRUE
#   s.1.GIC <- -Inf
#   s.2.GIC <- -Inf
#   while(loop.range){
#     if(s.1.GIC < GIC.result){
#       s.1.prop[study.l] <- max(s.1.prop[study.l] - initial.range / 2, 0)
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.1.prop)
#       s.1.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                delta = - ref.quant, support.size = lambda)
#       s.1.GIC <- GIC.cal(result = s.1.result, tune.type = tune.type)
#     }
#     if(s.2.GIC < GIC.result){
#       s.2.prop[study.l] <- min(s.2.prop[study.l] + initial.range / 2, 1)
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.2.prop)
#       s.2.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                delta = - ref.quant, support.size = lambda)
#       s.2.GIC <- GIC.cal(result = s.2.result, tune.type = tune.type)
#     }
#     if(s.1.GIC > GIC.result & s.2.GIC > GIC.result){
#       loop.range <- FALSE
#     }else if(s.1.prop[study.l] == 0 | s.2.prop[study.l] == 1){
#       loop.range <- FALSE
#     }
#   }
#   #===============================================#
#   s.left.prop[study.l] <- s.1.prop[study.l] * p.pi + s.2.prop[study.l] * (1 - p.pi)
#   ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.left.prop)
#   s.left.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                               delta = - ref.quant, support.size = lambda)
#   s.left.GIC <- GIC.cal(result = s.left.result, tune.type = tune.type)
#   #===============================================#
#   s.right.prop[study.l] <- s.1.prop[study.l] * (1 - p.pi) + s.2.prop[study.l] * p.pi
#   ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.right.prop)
#   s.right.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                delta = - ref.quant, support.size = lambda)
#   s.right.GIC <- GIC.cal(result = s.right.result, tune.type = tune.type)
#   #===============================================#
#   # we have:
#   # s.1.result; s.1.GIC
#   # s.2.result; s.2.GIC
#   # s.left.result; s.left.GIC
#   # s.right.result; s.right.GIC
#   #===============================================#
#   loop.study <- TRUE
#   loop.times <- 0
#   while(loop.study){
#     #=== check BIC; or prop ===#
#     s.all.prop <- c(s.1.prop[study.l], s.left.prop[study.l], s.right.prop[study.l], s.2.prop[study.l])
#     s.all.GIC <-  c(s.1.GIC, s.left.GIC, s.right.GIC, s.2.GIC)
#     if(s.all.GIC[1] == min(s.all.GIC) & abs(s.all.prop[1] - s.all.prop[2]) <= tol){
#       loop.study <- FALSE
#       s.new.prop <- s.1.prop
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
#       s.new.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                  delta = - ref.quant, support.size = lambda)
#       s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
#       s.new.result$GIC.result <- s.new.GIC
#     }else if(s.all.GIC[4] == min(s.all.GIC) & abs(s.all.prop[4] - s.all.prop[3]) <= tol){
#       loop.study <- FALSE
#       s.new.prop <- s.2.prop
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
#       s.new.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                  delta = - ref.quant, support.size = lambda)
#       s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
#       s.new.result$GIC.result <- s.new.GIC
#     }else if(s.all.GIC[2] == min(s.all.GIC) & abs(s.all.prop[1] - s.all.prop[3]) <= tol){
#       loop.study <- FALSE
#       s.new.prop <- (s.1.prop + s.right.prop) / 2
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
#       s.new.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                  delta = - ref.quant, support.size = lambda)
#       s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
#       s.new.result$GIC.result <- s.new.GIC
#     }else if(s.all.GIC[3] == min(s.all.GIC) & abs(s.all.prop[2] - s.all.prop[4]) <= tol){
#       loop.study <- FALSE
#       s.new.prop <- (s.left.prop + s.2.prop) / 2
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
#       s.new.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                  delta = - ref.quant, support.size = lambda)
#       s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
#       s.new.result$GIC.result <- s.new.GIC
#     }else if(loop.times >= 100){
#       loop.study <- FALSE
#       warning("Search delta over 100 loops. Not converge. \n")
#       s.new.prop <- (s.1.prop + s.2.prop) / 2
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
#       s.new.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                  delta = - ref.quant, support.size = lambda)
#       s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
#       s.new.result$GIC.result <- s.new.GIC
#     }
#     if(min(s.all.GIC[1:2]) == min(s.all.GIC)){
#       #=== minimize GIC at s.left or s.1 ===#
#       s.2.prop[study.l] <- s.right.prop[study.l]
#       s.2.GIC <- s.right.GIC
#       s.right.prop[study.l] <- s.left.prop[study.l]
#       s.right.GIC <- s.left.GIC
#       s.left.prop[study.l] <- s.1.prop[study.l] * p.pi + s.2.prop[study.l] * (1 - p.pi)
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.left.prop)
#       s.left.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                   delta = - ref.quant, support.size = lambda)
#       s.left.GIC <- GIC.cal(result = s.left.result, tune.type = tune.type)
#     }else if(min(s.all.GIC[3:4]) == min(s.all.GIC)){
#       #=== minimize GIC at s.right or s.2 ===#
#       s.1.prop[study.l] <- s.left.prop[study.l]
#       s.1.GIC <- s.left.GIC
#       s.left.prop[study.l] <- s.right.prop[study.l]
#       s.left.GIC <- s.right.GIC
#       s.right.prop[study.l] <- s.1.prop[study.l] * (1 - p.pi) + s.2.prop[study.l] * p.pi
#       ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.right.prop)
#       s.right.result <- byabess_AA(summary.stat.study = summary.stat.study, lasso.mat = lasso.mat,
#                                    delta = - ref.quant, support.size = lambda)
#       s.right.GIC <- GIC.cal(result = s.right.result, tune.type = tune.type)
#     }
#     #=== Calculate loop time ===#
#     loop.times <- loop.times + 1
#   }
#   return(list(study.min = s.new.prop, study.min.GIC = s.new.GIC, result = s.new.result))
# }

# Utility function
# search.subset.s <- function(summary.stat.study,
#                             lasso.mat,
#                             L,
#                             quantile_sm,
#                             s.size,
#                             tmp.result,
#                             search.loc,
#                             GIC.result,
#                             initial.range,
#                             tune.type,
#                             tol,
#                             verbose){
#   loop <- TRUE
#   loop.time <- 0
#   if(is.null(tmp.result)){
#     for(l in 1:L){
#       tmp.search <- search.ref.loc(summary.stat.study = summary.stat.study,
#                                    lasso.mat = lasso.mat,
#                                    study.l = l,
#                                    quantile_sm = quantile_sm,
#                                    search.prop = search.loc,
#                                    lambda = s.size,
#                                    result = tmp.result,
#                                    initial.range = initial.range,
#                                    tune.type = tune.type,
#                                    tol = tol,
#                                    verbose = verbose)
#
#       GIC.result <- tmp.search$study.min.GIC
#       tmp.result <- tmp.search$result
#       search.loc <- tmp.search$study.min
#     }
#   }
#   while(loop){
#     search.loc.tmp <- search.loc
#     GIC.result.tmp <- GIC.result
#     for(l in 1:L){
#       tmp.search <- search.ref.loc(summary.stat.study = summary.stat.study,
#                                    lasso.mat = lasso.mat,
#                                    study.l = l,
#                                    quantile_sm = quantile_sm,
#                                    search.prop = search.loc.tmp,
#                                    lambda = s.size,
#                                    result = tmp.result,
#                                    initial.range = initial.range,
#                                    tune.type = tune.type,
#                                    tol = tol,
#                                    verbose = verbose)
#       if(tmp.search$study.min.GIC < GIC.result.tmp){
#         search.loc.tmp <- tmp.search$study.min
#         GIC.result.tmp <- tmp.search$study.min.GIC
#         tmp.result <- tmp.search$result
#       }
#     }
#     loop.time <- loop.time + 1
#     # cat(paste0("loop time ", loop.time, ":", paste0(search.loc - search.loc.tmp, collapse = ","), ".\n"))
#     search.loc <- search.loc.tmp
#     if(abs(GIC.result - GIC.result.tmp) <= tol){
#       loop <- FALSE
#     }else if(loop.time >= 100){
#       loop <- FALSE
#       warning(paste0("Loop over 100 times for searching best subset of ", s, ".\n"))
#     }else{
#       GIC.result <- GIC.result.tmp
#     }
#   }
#   return(list(tmp.result = tmp.result, GIC.result = GIC.result, search.loc = search.loc))
# }

# Utility function
cal.quantile <- function(quantile_sm, prop){
  tmp_quant <- c()
  for(l in 1:length(prop)){
    tmp_quant <- c(tmp_quant, prop[l] * quantile_sm[l,2] + (1 - prop[l]) * quantile_sm[l,1])
  }
  return(tmp_quant)
}

# Meta-analysis
meta.analysis <- function(summary.stat.study = summary.stat.study,
                          lasso.mat = lasso.mat,
                          quantile_sm = quantile_sm,
                          cov.name,
                          ref = ref,
                          tune.path = tune.path,
                          tune.size.sequence = tune.size.range,
                          tune.size.range = tune.size.range,
                          tune.type = tune.type,
                          output.best.one = output.best.one,
                          tol = tol,
                          NMAX = NMAX,
                          verbose = verbose){

  study.ID <- names(summary.stat.study)

  L <- length(study.ID)
  tune.set <- list()
  GIC.tau <- c()
  #===========================================#
  search.loc <- rep(-Inf, L)
  tmp.result <- NULL
  GIC.result <- Inf
  if(tune.path == "sequence"){
    if(is.null(tune.size.sequence) & is.null(tune.size.range)){
      dims <- dim(lasso.mat$X.enlarge)
      support.sizes <- 1:min(dims[1], round(dims / 2))
    }else if(!is.null(tune.size.sequence) & is.null(tune.size.range)){
      support.sizes <- tune.size.sequence
    }else if(is.null(tune.size.sequence) & !is.null(tune.size.range)){
      if(!all(tune.size.range == round(tune.size.range))){
        stop("The tune.size.range include non-interger.\n")
      }
      if(tune.size.range[1] < 0 | (tune.size.range[2] - tune.size.range[1]) < 2){
        stop("The range for tune.size.range include negative value or the second value is not large enough than the first value.\n")
      }
      support.sizes <- tune.size.range[1]:tune.size.range[2]
    }else if(!is.null(tune.size.sequence) & !is.null(tune.size.range)){
      stop("Please only provide `tune.size.range` or  `tune.size.sequence`.\n")
    }
    if(verbose == TRUE){
      pb <- txtProgressBar(max=length(support.sizes), style=3)
    }
    #=== Search best model by sequence ===#

    for(s.lambda in support.sizes){
      tmp.subset <- search.subset.s(summary.stat.study = summary.stat.study,
                                    lasso.mat = lasso.mat,
                                    L = L,
                                    quantile_sm = quantile_sm,
                                    s.size = s.lambda,
                                    tmp.result = tmp.result,
                                    search.loc = search.loc,
                                    GIC.result = GIC.result,
                                    initial.range = 0.1,
                                    tune.type = tune.type,
                                    tol = tol,
                                    NMAX = NMAX,
                                    verbose = verbose)
      if(s.lambda != sum(tmp.subset$tmp.result$mu.fit!=0)){
        tmp.subset <- search.subset.s(summary.stat.study = summary.stat.study,
                                      lasso.mat = lasso.mat,
                                      L = L,
                                      quantile_sm = quantile_sm,
                                      s.size = s.lambda,
                                      tmp.result = NULL,
                                      search.loc = rep(-Inf, L),
                                      GIC.result = Inf,
                                      initial.range = 0.1,
                                      tune.type = tune.type,
                                      tol = tol,
                                      NMAX = NMAX,
                                      verbose = verbose)
      }
      tmp.result <- tmp.subset$tmp.result
      GIC.result <- tmp.subset$GIC.result
      search.loc <- tmp.subset$search.loc
      tune.set[[as.character(s.lambda)]] <- tmp.result
      GIC.result <- GIC.cal(result = tmp.result, tune.type = tune.type)
      GIC.tau <- c(GIC.tau, GIC.result)
      if(min(GIC.tau) == GIC.result){
        min.result <- tmp.result
      }
      if(verbose){
        setTxtProgressBar(pb, (pb$getVal()+1))
      }
    }
    if(verbose == TRUE){
      cat("\n")
    }
    if(output.best.one){
      coef <- min.result$mu.fit
      delta <- min.result$delta
      dev <- min.result$q_loss
      ic <- min.result$GIC.result - dev
      names(dev) <- as.character(sum(min.result$mu.fit))
      names(ic) <- names(dev)
      names(delta) <- paste0("<",study.ID,">_<" ,ref,">")
      results.all <- list(coef = coef,
                          delta = delta,
                          dev = dev,
                          ic = ic)
      if(sum(coef!=0) == support.sizes[length(support.sizes)]){
        warning("The searching hits the boundary for covariate of interest ", cov.name, ". ++")
        boud.hit <- TRUE
      }else{
        boud.hit <- FALSE
      }
    }else{
      coef <- NULL
      delta <- NULL
      dev <- NULL
      ic <- NULL
      for(l in 1:length(tune.set)){
        coef <- cbind(coef, tune.set[[l]]$mu.fit)
        delta <- cbind(delta, tune.set[[l]]$delta)
        dev <- c(dev, tune.set[[l]]$q_loss)
      }
      ic <- GIC.tau - dev
      names(dev) <- as.character(support.sizes)
      names(ic) <- names(dev)
      rownames(delta) <- paste0("<",study.ID,">_<" ,ref,">")
      results.all <- list(coef = coef,
                          delta = delta,
                          dev = dev,
                          ic = ic)
      if(sum(coef[,which.min(GIC.tau)]!=0) == support.sizes[length(support.sizes)]){
        warning("The searching hits the boundary for covariate of interest ", cov.name, ". ++")
        boud.hit <- TRUE
      }else{
        boud.hit <- FALSE
      }
    }
  }else if(tune.path == "gsection"){
    #=== Search best model by g-section ===#
    if(is.null(tune.size.range)){
      dims <- dim(lasso.mat$X.enlarge)
      s.1.lambda <- 0
      s.2.lambda <- min(dims[1], round(dims / 2))
      support.sizes <- s.1.lambda:s.2.lambda
    }else{
      if(!all(tune.size.range == round(tune.size.range))){
        stop("The tune.size.range include non-interger.\n")
      }
      if(tune.size.range[1] < 0 | (tune.size.range[2] - tune.size.range[1]) < 2){
        stop("The range for tune.size.range include negative value or the second value is not large enough than the first value.\n")
      }
      s.1.lambda <- tune.size.range[1]
      s.2.lambda <- tune.size.range[2]
      support.sizes <- s.1.lambda : s.2.lambda
    }
    loops <- TRUE
    g.section <- 2 / (sqrt(5) + 1)
    g.sequence <- c()
    g.id <- 1
    s.left.lambda <- round(s.1.lambda * g.section + s.2.lambda * (1 - g.section))
    s.right.lambda <- round(s.1.lambda * (1 - g.section) + s.2.lambda * g.section)
    for(s.lambda in c(s.1.lambda, s.left.lambda, s.right.lambda, s.2.lambda)){
      tmp.subset <- search.subset.s(summary.stat.study = summary.stat.study,
                                    lasso.mat = lasso.mat,
                                    L = L,
                                    quantile_sm = quantile_sm,
                                    s.size = s.lambda,
                                    tmp.result = NULL,
                                    search.loc = rep(-Inf, L),
                                    GIC.result = Inf,
                                    initial.range = 0.2,
                                    tune.type = tune.type,
                                    tol = tol,
                                    NMAX = NMAX,
                                    verbose = verbose)
      tune.set[[as.character(s.lambda)]] <- tmp.subset$tmp.result
      GIC.result <- GIC.cal(result = tmp.subset$tmp.result, tune.type = tune.type)
      GIC.tau <- c(GIC.tau, GIC.result)
      g.sequence <- c(g.sequence, s.lambda)
      g.id <- g.id + 1
    }
    s.1.GIC <- GIC.tau[1]
    s.left.GIC <- GIC.tau[2]
    s.right.GIC <- GIC.tau[3]
    s.2.GIC <- GIC.tau[4]
    while (loops) {
      s.all.GIC <- c(s.1.GIC, s.left.GIC, s.right.GIC, s.2.GIC)
      if(s.1.GIC == min(s.all.GIC) & s.left.lambda - s.1.lambda <= 1){
        min.result <- tune.set[[as.character(s.1.lambda)]]
        loops <- FALSE
      }else if(s.2.GIC == min(s.all.GIC) & s.2.lambda - s.right.lambda <= 1){
        min.result <- tune.set[[as.character(s.2.lambda)]]
        loops <- FALSE
      }else if(s.left.GIC == min(s.all.GIC) & s.left.lambda - s.1.lambda <= 1 & s.right.lambda - s.left.lambda <= 1){
        min.result <- tune.set[[as.character(s.left.lambda)]]
        loops <- FALSE
      }else if(s.right.GIC == min(s.all.GIC) & s.right.lambda - s.left.lambda <= 1 & s.2.lambda - s.right.lambda <= 1){
        min.result <- tune.set[[as.character(s.right.lambda)]]
        loops <- FALSE
      }else{
        if(min(s.all.GIC[1:2]) == min(s.all.GIC)){
          #=== minimize GIC at s.1 or s.left ===#
          s.2.lambda <- s.right.lambda
          s.2.GIC <- s.right.GIC
          s.right.lambda <- s.left.lambda
          s.right.GIC <- s.left.GIC
          s.left.lambda <- round(s.1.lambda * g.section + s.2.lambda * (1 - g.section))
          tmp.subset <- search.subset.s(summary.stat.study = summary.stat.study,
                                        lasso.mat = lasso.mat,
                                        L = L,
                                        quantile_sm = quantile_sm,
                                        s.size = s.left.lambda,
                                        tmp.result = NULL,
                                        search.loc = rep(-Inf, L),
                                        GIC.result = Inf,
                                        initial.range = 0.2,
                                        tune.type = tune.type,
                                        tol = tol,
                                        NMAX = NMAX,
                                        verbose = verbose)
          tmp.result <- tmp.subset$tmp.result
          GIC.result <- tmp.subset$GIC.result
          search.loc <- tmp.subset$search.loc
          if(!as.character(s.left.lambda) %in% names(tune.set)){
            tune.set[[as.character(s.left.lambda)]] <- tmp.result
            s.left.GIC <- GIC.cal(result = tmp.result, tune.type = tune.type)
            GIC.tau <- c(GIC.tau, s.left.GIC)
            g.sequence <- c(g.sequence, s.left.lambda)
          }
        }else if(min(s.all.GIC[3:4]) == min(s.all.GIC)){
          #=== minimize GIC at s.right or s.2 ===#
          s.1.lambda <- s.left.lambda
          s.1.GIC <- s.left.GIC
          s.left.lambda <- s.right.lambda
          s.left.GIC <- s.right.GIC
          s.right.lambda <- round(s.1.lambda * (1 - g.section) + s.2.lambda * g.section)
          tmp.subset <- search.subset.s(summary.stat.study = summary.stat.study,
                                        lasso.mat = lasso.mat,
                                        L = L,
                                        quantile_sm = quantile_sm,
                                        s.size = s.right.lambda,
                                        tmp.result = NULL,
                                        search.loc = rep(-Inf, L),
                                        GIC.result = Inf,
                                        initial.range = 0.2,
                                        tune.type = tune.type,
                                        tol = tol,
                                        NMAX = NMAX,
                                        verbose = verbose)
          tmp.result <- tmp.subset$tmp.result
          GIC.result <- tmp.subset$GIC.result
          search.loc <- tmp.subset$search.loc
          if(!as.character(s.right.lambda) %in% names(tune.set)){
            tune.set[[as.character(s.right.lambda)]] <- tmp.result
            s.right.GIC <- GIC.cal(result = tmp.result, tune.type = tune.type)
            GIC.tau <- c(GIC.tau, s.right.GIC)
            g.sequence <- c(g.sequence, s.right.lambda)
          }
        }
        g.id <- g.id + 1
      }
    }
    if(output.best.one){
      coef <- min.result$mu.fit
      delta <- min.result$delta
      names(delta) <- paste0("<",study.ID,">_<" ,ref,">")
      dev <- min.result$q_loss
      ic <- min.result$GIC.result - dev
      names(dev) <- as.character(sum(min.result$mu.fit))
      names(ic) <- names(dev)
      results.all <- list(coef = coef,
                          delta = delta,
                          dev = dev,
                          ic = ic)
      if(sum(coef!=0) == support.sizes[length(support.sizes)]){
        warning("The searching hits the boundary for covariate of interest ", cov.name, ". ++")
        boud.hit <- TRUE
      }else{
        boud.hit <- FALSE
      }
    }else{
      coef <- NULL
      delta <- NULL
      dev <- NULL
      ic <- NULL
      for(l in 1:length(tune.set)){
        coef <- cbind(coef, tune.set[[l]]$mu.fit)
        delta <- cbind(delta, tune.set[[l]]$delta)
        dev <- c(dev, tune.set[[l]]$q_loss)
      }
      ic <- GIC.tau - dev
      rownames(delta) <- paste0("<",study.ID,">_<" ,ref,">")
      names(dev) <- as.character(g.sequence)
      names(ic) <- names(dev)
      results.all <- list(coef = coef,
                          delta = delta,
                          dev = dev,
                          ic = ic)
      if(sum(coef[,which.min(GIC.tau)]!=0) == support.sizes[length(support.sizes)]){
        warning("The searching hits the boundary for covariate of interest ", cov.name, ". ++")
        boud.hit <- TRUE
      }else{
        boud.hit <- FALSE
      }
    }
  }
  return(list(results.all = results.all, boud.hit = boud.hit))
}

################################### Powell's method ###########################
f_xy <- function(x, y, GIC_x, GIC_y){
  return((GIC_y - GIC_x) / (y - x))
}

f_xyz <- function(x, y, z, GIC_x, GIC_y, GIC_z){
  return((f_xy(x = y, y = z, GIC_x = GIC_y, GIC_y = GIC_z) -
            f_xy(x = x, y = y, GIC_x = GIC_x, GIC_y = GIC_y)) / (z - x))
}

x_t <- function(x0, x1,x2, GIC_x0, GIC_x1, GIC_x2){
  return((f_xyz(x = x0, y = x1, z = x2, GIC_x = GIC_x0, GIC_y = GIC_x1, GIC_z = GIC_x2) * (x0 + x1) -
            f_xy(x = x0, y = x1, GIC_x = GIC_x0, GIC_y = GIC_x1))/(2 * f_xyz(x = x0, y = x1, z = x2, GIC_x = GIC_x0, GIC_y = GIC_x1, GIC_z = GIC_x2)))
}

nt <- function(seqs, xt){
  return(seqs[which.min(abs(seqs - xt))[1]])
}

ft <- function(seqs, xt){
  return(seqs[which.max(abs(seqs - xt))[1]])
}

######################## version 2: Powell + quadratic interpolation #############
search.subset.s <- function(summary.stat.study,
                            lasso.mat,
                            L,
                            quantile_sm,
                            s.size,
                            tmp.result,
                            search.loc,
                            GIC.result,
                            initial.range,
                            tune.type,
                            tol,
                            NMAX,
                            verbose){

  loop <- TRUE
  loop.time <- 0
  loop.mat <- diag(L)
  if(is.null(tmp.result)){
    for(l in 1:L){
      tmp.search <- search.ref.loc(summary.stat.study = summary.stat.study,
                                   lasso.mat = lasso.mat,
                                   study.l = loop.mat[l,],
                                   quantile_sm = quantile_sm,
                                   search.prop = search.loc,
                                   lambda = s.size,
                                   result = tmp.result,
                                   initial.range = initial.range,
                                   tune.type = tune.type,
                                   tol = tol,
                                   verbose = verbose)

      GIC.result <- tmp.search$study.min.GIC
      tmp.result <- tmp.search$result
      search.loc <- tmp.search$study.min
    }
  }
  while(loop){
    search.loc.tmp <- search.loc
    GIC.result.tmp <- GIC.result
    for(l in 1:L){
      tmp.search <- search.ref.loc(summary.stat.study = summary.stat.study,
                                   lasso.mat = lasso.mat,
                                   study.l = loop.mat[l,],
                                   quantile_sm = quantile_sm,
                                   search.prop = search.loc.tmp,
                                   lambda = s.size,
                                   result = tmp.result,
                                   initial.range = initial.range,
                                   tune.type = tune.type,
                                   tol = tol,
                                   verbose = verbose)

      if(tmp.search$study.min.GIC < GIC.result.tmp){
        search.loc.tmp <- tmp.search$study.min
        GIC.result.tmp <- tmp.search$study.min.GIC
        tmp.result <- tmp.search$result
      }
    }
    ## Search on the new direction
    if(max(abs(search.loc - search.loc.tmp)) > tol){
      loop.new.vec <- (search.loc - search.loc.tmp) / max(abs(search.loc - search.loc.tmp))
      tmp.search <- search.ref.loc(summary.stat.study = summary.stat.study,
                                   lasso.mat = lasso.mat,
                                   study.l = loop.new.vec,
                                   quantile_sm = quantile_sm,
                                   search.prop = search.loc.tmp,
                                   lambda = s.size,
                                   result = tmp.result,
                                   initial.range = initial.range,
                                   tune.type = tune.type,
                                   tol = tol,
                                   verbose = verbose)
      if(tmp.search$study.min.GIC < GIC.result.tmp){
        search.loc.tmp <- tmp.search$study.min
        GIC.result.tmp <- tmp.search$study.min.GIC
        tmp.result <- tmp.search$result
      }
    }
    loop.time <- loop.time + 1

    # cat(paste0("loop time ", loop.time, ":", paste0(search.loc - search.loc.tmp, collapse = ","), ".\n"))
    if(max(abs(search.loc - search.loc.tmp)) <= tol){
      loop <- FALSE
    }else if(loop.time > NMAX){
      loop <- FALSE
      warning(paste0("Loops exceed the given loop times for searching best subset of ", s, ".\n"))
    }else{
      ## check loop.time is whether a multiplier of L
      if(loop.time %% L == 0){
        loop.mat <- diag(L)
      }else{
        rm.loop.idx <- which.max(apply(loop.mat, 1, function(d){abs(sum(d * (search.loc - search.loc.tmp))/sqrt(sum(d^2)))}))[1]
        loop.mat[rm.loop.idx,] <- loop.new.vec
      }
      GIC.result <- GIC.result.tmp
      search.loc <- search.loc.tmp
    }
  }
  return(list(tmp.result = tmp.result, GIC.result = GIC.result, search.loc = search.loc))
}


search.ref.loc <- function(summary.stat.study,
                           lasso.mat,
                           study.l,
                           quantile_sm,
                           search.prop,
                           lambda,
                           result,
                           initial.range = 0.1,
                           tune.type,
                           tol = 1e-3,
                           verbose = FALSE){

  feature.set <- lasso.mat$feature.set
  study.ID <- lasso.mat$study.ID
  feature.ID <- lasso.mat$feature.ID
  ref <- lasso.mat$ref
  L <- length(study.ID)
  non.zero.dir <- which(study.l!=0)[1]
  if(initial.range > 1){
    initial.range <- 1
  }
  if(is.null(result)){
    search.prop <- rep(1/2, L)
    ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = search.prop)
    search.result <- byabess_AA(summary.stat.study = summary.stat.study,
                                lasso.mat = lasso.mat,
                                delta = - ref.quant,
                                support.size = lambda)
    GIC.x <- GIC.cal(result = search.result, tune.type = tune.type)
  }else{
    GIC.x <- GIC.cal(result = result, tune.type = tune.type)
  }
  x.prop <- search.prop
  s_step <- 0.1
  m.mean <- 0.1
  epsilon <- tol
  s.prop <- search.prop
  s.prop <- s.prop + s_step * study.l
  ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.prop)
  s.result <- byabess_AA(summary.stat.study = summary.stat.study,
                         lasso.mat = lasso.mat,
                         delta = - ref.quant,
                         support.size = lambda)
  GIC.s <- GIC.cal(result = s.result, tune.type = tune.type)

  if(GIC.x < GIC.s){
    x0.prop <- x.prop
    x1.prop <- x.prop
    x2.prop <- x.prop
    x0.prop <- x0.prop - s_step * study.l
    x2.prop <- x2.prop + s_step * study.l
    GIC_x1 <- GIC.x
    ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x0.prop)
    x0.result <- byabess_AA(summary.stat.study = summary.stat.study,
                            lasso.mat = lasso.mat,
                            delta = - ref.quant,
                            support.size = lambda)
    GIC_x0 <- GIC.cal(result = x0.result, tune.type = tune.type)
    ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x2.prop)
    x2.result <- byabess_AA(summary.stat.study = summary.stat.study,
                            lasso.mat = lasso.mat,
                            delta = - ref.quant,
                            support.size = lambda)
    GIC_x2 <- GIC.cal(result = x2.result, tune.type = tune.type)
  }else{
    x0.prop <- x.prop
    x1.prop <- x.prop
    x2.prop <- x.prop
    x1.prop <- x1.prop + s_step * study.l
    x2.prop <- x2.prop + 2 * s_step * study.l
    GIC_x0 <- GIC.x
    ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x1.prop)
    x1.result <- byabess_AA(summary.stat.study = summary.stat.study,
                            lasso.mat = lasso.mat,
                            delta = - ref.quant,
                            support.size = lambda)
    GIC_x1 <- GIC.cal(result = x1.result, tune.type = tune.type)
    ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x2.prop)
    x2.result <- byabess_AA(summary.stat.study = summary.stat.study,
                            lasso.mat = lasso.mat,
                            delta = - ref.quant,
                            support.size = lambda)
    GIC_x2 <- GIC.cal(result = x2.result, tune.type = tune.type)
  }

  # seqs <- c(sqrt(sum((x0.prop - x.prop)^2)), sqrt(sum((x1.prop - x.prop)^2)), sqrt(sum((x2.prop - x.prop)^2))) *
  #   sign(c((x0.prop - x.prop)[non.zero.dir], (x1.prop - x.prop)[non.zero.dir], (x2.prop - x.prop)[non.zero.dir]))

  seqs <- c((x0.prop - x.prop)[non.zero.dir], (x1.prop - x.prop)[non.zero.dir], (x2.prop - x.prop)[non.zero.dir])
  GIC.seqs <- c(GIC_x0, GIC_x1, GIC_x2)

  loop.range <- TRUE
  inter.times <- 0
  while(loop.range){
    tmp.seqs <- seqs
    inter.times <-  inter.times  + 1
    m <- m.mean * runif(1, 0.8, 1.2)
    # cat(paste0("loc: ", paste(seqs, collapse = ","), ".\n"))
    # cat(paste0("GIC: ", paste(GIC.seqs, collapse = ","), ".\n"))
    s.iter <- x_t(x0 = seqs[1], x1 = seqs[2], x2 = seqs[3],
                  GIC_x0 = GIC.seqs[1], GIC_x1 = GIC.seqs[2], GIC_x2 = GIC.seqs[3])
    s.n <- nt(seqs = seqs, xt = s.iter)
    s.f <- ft(seqs = seqs, xt = s.iter)
    f2 <- f_xyz(x =seqs[1], y =seqs[2], z = seqs[3], GIC_x = GIC.seqs[1],
                GIC_y = GIC.seqs[2], GIC_z = GIC.seqs[3])
    if(f2 > 0 & abs(s.f - s.n) > m){
      rm.idx <- which(seqs == s.f)
      seqs <- seqs[-rm.idx]
      GIC.seqs <- GIC.seqs[-rm.idx]
      if(GIC.seqs[2] > GIC.seqs[1]){
        s.new <- seqs[1] - m
        x.new.prop <- x.prop + study.l * s.new
        ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x.new.prop)
        x.new.result <- byabess_AA(summary.stat.study = summary.stat.study,
                                   lasso.mat = lasso.mat,
                                   delta = - ref.quant,
                                   support.size = lambda)
        GIC_x_new <- GIC.cal(result = x.new.result, tune.type = tune.type)
        seqs <- c(s.new, seqs)
        GIC.seqs <- c(GIC_x_new, GIC.seqs)
      }else{
        s.new <- seqs[2] + m
        x.new.prop <- x.prop + study.l * s.new
        ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x.new.prop)
        x.new.result <- byabess_AA(summary.stat.study = summary.stat.study,
                                   lasso.mat = lasso.mat,
                                   delta = - ref.quant,
                                   support.size = lambda)
        GIC_x_new <- GIC.cal(result = x.new.result, tune.type = tune.type)
        seqs <- c(seqs, s.new)
        GIC.seqs <- c(GIC.seqs, GIC_x_new)
      }
    }else if(f2 < 0){
      rm.idx <- which(seqs == s.n)
      seqs <- seqs[-rm.idx]
      GIC.seqs <- GIC.seqs[-rm.idx]
      if(GIC.seqs[2] > GIC.seqs[1]){
        s.new <- seqs[1] - m
        x.new.prop <- x.prop + study.l * s.new
        ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x.new.prop)
        x.new.result <- byabess_AA(summary.stat.study = summary.stat.study,
                                   lasso.mat = lasso.mat,
                                   delta = - ref.quant,
                                   support.size = lambda)
        GIC_x_new <- GIC.cal(result = x.new.result, tune.type = tune.type)
        seqs <- c(s.new, seqs)
        GIC.seqs <- c(GIC_x_new, GIC.seqs)
      }else{
        s.new <- seqs[2] + m
        x.new.prop <- x.prop + study.l * s.new
        ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x.new.prop)
        x.new.result <- byabess_AA(summary.stat.study = summary.stat.study,
                                   lasso.mat = lasso.mat,
                                   delta = - ref.quant,
                                   support.size = lambda)
        GIC_x_new <- GIC.cal(result = x.new.result, tune.type = tune.type)
        seqs <- c(seqs, s.new)
        GIC.seqs <- c(GIC.seqs, GIC_x_new)
      }
    }else{
      if(abs(s.iter - s.n) < epsilon | inter.times > 30){
        loop.range <- FALSE
        s.star <- (s.n + s.iter) / 2
        x.star.prop <- x.prop + study.l * s.star
        ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x.star.prop)
        x.star.result <- byabess_AA(summary.stat.study = summary.stat.study,
                                    lasso.mat = lasso.mat,
                                    delta = - ref.quant,
                                    support.size = lambda)
        GIC_x_star <- GIC.cal(result = x.star.result, tune.type = tune.type)
        x.star.result$GIC.result <- GIC_x_star
      }else{
        max.idx <- which.max(GIC.seqs)
        seqs[max.idx] <- s.iter
        x.iter.prop <- x.prop + study.l * s.iter
        ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = x.iter.prop)
        x.iter.result <- byabess_AA(summary.stat.study = summary.stat.study,
                                    lasso.mat = lasso.mat,
                                    delta = - ref.quant,
                                    support.size = lambda)
        GIC.seqs[max.idx] <- GIC.cal(result = x.iter.result, tune.type = tune.type)
      }
    }
    order.seqs <- order(seqs)
    seqs <- seqs[order.seqs]
    GIC.seqs <- GIC.seqs[order.seqs]
  }
  return(list(study.min = x.star.prop, study.min.GIC = GIC_x_star, result = x.star.result))
}
