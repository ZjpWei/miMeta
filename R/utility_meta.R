#' @import abess

######## Utility file for meta-analysis ########
Get_lasso_pre = function(Melody){
  summary.stat.study <- Melody$summary.stat.study
  taxa.set <- Melody$taxa.set
  L <- Melody$dat.inf$L
  K <- Melody$dat.inf$K
  ref <- c()
  for(l in 1:L){
    ref <- c(ref, match(summary.stat.study[[l]]$ref, Melody$dat.inf$taxa.names))
  }
  taxa.mat <- NULL
  if(L == 1){
    idx = c(setdiff(1:K, ref[l]), ref[l])
    taxa.mat <- matrix(c(taxa.set[[1]], TRUE)[order(idx)], nrow = 1)
  }else{
    for(l in 1:L){
      idx = c(setdiff(1:K, ref[l]), ref[l])
      taxa.mat <- rbind(taxa.mat, c(taxa.set[[l]], TRUE)[order(idx)])
    }
  }

  taxa.tmp <- taxa.mat
  ref.tmp <- ref
  mu.len <- sum( colSums(taxa.mat) > 0)
  nonempty.id <- which(colSums(taxa.mat) > 0)
  empty.id <- which(colSums(taxa.mat) == 0)
  idx.lst <- list()
  for(l in 1:L){
    ref[l] <- sum(colSums(taxa.mat)[1:ref[l]] > 0)
    idx = c(setdiff(1:mu.len, ref[l]), ref[l])
    idx.lst[[l]] = order(idx)
  }
  taxa.mat <- taxa.mat[,colSums(taxa.mat, na.rm = TRUE) > 0]
  if(L == 1){
    taxa.mat <- matrix(taxa.mat, nrow = 1)
  }

  N <- 0
  k.l <- NULL
  for(l in 1:L){
    N <- N + summary.stat.study[[l]]$n
    k.l <- c(k.l, sum(taxa.set[[l]]))
  }
  for(l in 1:L){
    ## Generate design matrix X and Y
    Sigma.chol <- chol(solve(summary.stat.study[[l]]$cov)/N)
    F.mat <- diag(mu.len)
    id.F <- taxa.mat[l,]
    id.F[ref[l]] <- FALSE
    F.mat <- F.mat[id.F,]
    if(l == 1){
      X.col.1 <- Sigma.chol %*% F.mat
      Y.enlarge.1 <- Sigma.chol %*% (summary.stat.study[[l]]$est)
      Y.enlarge.2 <- Sigma.chol %*% rep(1, sum(taxa.set[[l]]))
    }else{
      X.col.1 <- rbind(X.col.1, Sigma.chol %*% F.mat)
      Y.enlarge.1 <- c(Y.enlarge.1, Sigma.chol %*% (summary.stat.study[[l]]$est) )
      Y.enlarge.2 <- c(Y.enlarge.2, Sigma.chol %*% rep(1, sum(taxa.set[[l]])) )
    }
  }
  Melody$lasso.mat$X.enlarge <- X.col.1
  Melody$lasso.mat$Y.enlarge.1 <- Y.enlarge.1
  Melody$lasso.mat$Y.enlarge.2 <- Y.enlarge.2
  Melody$lasso.mat$taxa.mat <- taxa.mat
  Melody$lasso.mat$idx.lst <- idx.lst
  Melody$lasso.mat$k.l <- k.l
  Melody$lasso.mat$N <- N
  Melody$lasso.mat$mu.len <- mu.len
  Melody$lasso.mat$nonempty.id <- nonempty.id
  Melody$lasso.mat$empty.id <- empty.id
  return(Melody)
}

### Optimization function
byabess_AA <- function(Melody, ref_est, support.size){
  ### Load data
  summary.stat.study <- Melody$summary.stat.study
  taxa.set <- Melody$taxa.set
  X.enlarge <- Melody$lasso.mat$X.enlarge
  Y.enlarge.1 <- Melody$lasso.mat$Y.enlarge.1
  Y.enlarge.2 <- Melody$lasso.mat$Y.enlarge.2
  idx.lst <- Melody$lasso.mat$idx.lst
  taxa.mat <- Melody$lasso.mat$taxa.mat
  k.l <- Melody$lasso.mat$k.l
  mu.len <- Melody$lasso.mat$mu.len
  nonempty.id <- Melody$lasso.mat$nonempty.id
  empty.id <- Melody$lasso.mat$empty.id
  L <- Melody$dat.inf$L

  ref.shift <- c()
  if(length(ref_est) == 1){
    for(l in 1:L){
      ref.shift <- c(ref.shift, rep(ref_est, sum(taxa.set[[l]])))
    }
  }else{
    for(l in 1:L){
      ref.shift <- c(ref.shift, rep(ref_est[l], sum(taxa.set[[l]])))
    }
  }
  Y.enlarge <- Y.enlarge.1 + Y.enlarge.2 * ref.shift

  ### Abess model
  result <- abess(x = X.enlarge,
                  y = Y.enlarge,
                  family = "gaussian",
                  tune.type = "bic",
                  tune.path = "sequence",
                  normalize = 0,
                  support.size = support.size)

  fit.coef <- result$beta
  q_loss <- sum( (Y.enlarge - X.enlarge %*% result$beta)^2 )
  mu.fit <- fit.coef[1:mu.len]
  mu.fit <- c(mu.fit, rep(0, length(empty.id)))[order(c(nonempty.id ,empty.id))]
  names(mu.fit) <- Melody$dat.inf$taxa.names
  names(ref_est) <- Melody$dat.inf$ref
  if(length(unique(Melody$dat.inf$ref)) == 1){
    mu.fit <- mu.fit[setdiff(names(mu.fit), unique(Melody$dat.inf$ref))]
  }
  return(list(mu.fit = mu.fit,
              ref_est = ref_est,
              N = Melody$lasso.mat$N,
              q_loss = q_loss))
}

### Calculate GIC
GIC.cal <- function(result,
                    tune.type = c("HBIC", "BIC", "KBIC", "EBIC")
){
  tune.type <- match.arg(tune.type)
  df <- sum(result$mu.fit!=0)
  K <- length(result$mu.fit)
  L <- length(result$ref_est)
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

### Search delta
search.ref.loc <- function(Melody,
                           study.l,
                           quantile_sm,
                           search.prop,
                           lambda,
                           result,
                           initial.range = 0.1,
                           tune.type,
                           tol = 1e-3,
                           verbose = FALSE){
  summary.stat.study <- Melody$summary.stat.study
  taxa.set <- Melody$taxa.set
  L <- Melody$dat.inf$L
  K <- Melody$dat.inf$K
  p.pi <- 2 / (sqrt(5) + 1)
  if(initial.range > 1){
    stop("Searching range should be smaller than 1 \n")
  }
  if(is.null(result)){
    s.1.prop <- rep(1/2, L)
    s.2.prop <- rep(1/2, L)
    s.left.prop <- rep(1/2, L)
    s.right.prop <- rep(1/2, L)
    ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.1.prop)
    initial.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
    GIC.result <- GIC.cal(result = initial.result, tune.type = tune.type)
  }else{
    s.1.prop <- search.prop
    s.2.prop <- search.prop
    s.left.prop <- search.prop
    s.right.prop <- search.prop
    GIC.result <- result$GIC.result
  }
  loop.range <- TRUE
  s.1.GIC <- -Inf
  s.2.GIC <- -Inf
  while(loop.range){
    if(s.1.GIC < GIC.result){
      s.1.prop[study.l] <- max(s.1.prop[study.l] - initial.range / 2, 0)
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.1.prop)
      s.1.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.1.GIC <- GIC.cal(result = s.1.result, tune.type = tune.type)
    }
    if(s.2.GIC < GIC.result){
      s.2.prop[study.l] <- min(s.2.prop[study.l] + initial.range / 2, 1)
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.2.prop)
      s.2.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.2.GIC <- GIC.cal(result = s.2.result, tune.type = tune.type)
    }
    if(s.1.GIC > GIC.result & s.2.GIC > GIC.result){
      loop.range <- FALSE
    }else if(s.1.prop[study.l] == 0 | s.2.prop[study.l] == 1){
      loop.range <- FALSE
    }
  }
  #####
  s.left.prop[study.l] <- s.1.prop[study.l] * p.pi + s.2.prop[study.l] * (1 - p.pi)
  ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.left.prop)
  s.left.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
  s.left.GIC <- GIC.cal(result = s.left.result, tune.type = tune.type)
  #####
  s.right.prop[study.l] <- s.1.prop[study.l] * (1 - p.pi) + s.2.prop[study.l] * p.pi
  ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.right.prop)
  s.right.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
  s.right.GIC <- GIC.cal(result = s.right.result, tune.type = tune.type)
  #####
  # we have:
  # s.1.result; s.1.GIC
  # s.2.result; s.2.GIC
  # s.left.result; s.left.GIC
  # s.right.result; s.right.GIC
  #####
  loop.study <- TRUE
  loop.times <- 0
  while(loop.study){
    ### check BIC; or prop
    s.all.prop <- c(s.1.prop[study.l], s.left.prop[study.l], s.right.prop[study.l], s.2.prop[study.l])
    s.all.GIC <-  c(s.1.GIC, s.left.GIC, s.right.GIC, s.2.GIC)
    if(s.all.GIC[1] == min(s.all.GIC) & abs(s.all.prop[1] - s.all.prop[2]) <= tol){
      loop.study <- FALSE
      s.new.prop <- s.1.prop
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
      s.new.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
      s.new.result$GIC.result <- s.new.GIC
    }else if(s.all.GIC[4] == min(s.all.GIC) & abs(s.all.prop[4] - s.all.prop[3]) <= tol){
      loop.study <- FALSE
      s.new.prop <- s.2.prop
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
      s.new.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
      s.new.result$GIC.result <- s.new.GIC
    }else if(s.all.GIC[2] == min(s.all.GIC) & abs(s.all.prop[1] - s.all.prop[3]) <= tol){
      loop.study <- FALSE
      s.new.prop <- (s.1.prop + s.right.prop) / 2
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
      s.new.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
      s.new.result$GIC.result <- s.new.GIC
    }else if(s.all.GIC[3] == min(s.all.GIC) & abs(s.all.prop[2] - s.all.prop[4]) <= tol){
      loop.study <- FALSE
      s.new.prop <- (s.left.prop + s.2.prop) / 2
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
      s.new.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
      s.new.result$GIC.result <- s.new.GIC
    }else if(loop.times >= 100){
      loop.study <- FALSE
      warning("Search delta over 100 loops. Not converge. \n")
      s.new.prop <- (s.1.prop + s.2.prop) / 2
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.new.prop)
      s.new.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.new.GIC <- GIC.cal(result = s.new.result, tune.type = tune.type)
      s.new.result$GIC.result <- s.new.GIC
    }
    if(min(s.all.GIC[1:2]) == min(s.all.GIC)){
      ### minimize GIC at s.left or s.1
      s.2.prop[study.l] <- s.right.prop[study.l]
      s.2.GIC <- s.right.GIC
      s.right.prop[study.l] <- s.left.prop[study.l]
      s.right.GIC <- s.left.GIC
      s.left.prop[study.l] <- s.1.prop[study.l] * p.pi + s.2.prop[study.l] * (1 - p.pi)
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.left.prop)
      s.left.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.left.GIC <- GIC.cal(result = s.left.result, tune.type = tune.type)
    }else if(min(s.all.GIC[3:4]) == min(s.all.GIC)){
      ### minimize GIC at s.right or s.2
      s.1.prop[study.l] <- s.left.prop[study.l]
      s.1.GIC <- s.left.GIC
      s.left.prop[study.l] <- s.right.prop[study.l]
      s.left.GIC <- s.right.GIC
      s.right.prop[study.l] <- s.1.prop[study.l] * (1 - p.pi) + s.2.prop[study.l] * p.pi
      ref.quant <- cal.quantile(quantile_sm = quantile_sm, prop = s.right.prop)
      s.right.result <- byabess_AA(Melody = Melody, ref_est = - ref.quant, support.size = lambda)
      s.right.GIC <- GIC.cal(result = s.right.result, tune.type = tune.type)
    }
    ### Calculate loop time
    loop.times <- loop.times + 1
  }
  return(list(study.min = s.new.prop, study.min.GIC = s.new.GIC, result = s.new.result))
}

### Utility function
search.subset.s <- function(Melody, L, quantile_sm, s.size, tmp.result, search.loc, GIC.result, initial.range, tune.type, tol, verbose){
  loop <- TRUE
  loop.time <- 0
  # if(verbose){
  #   message(paste0("-- searching model with subset size:", s.size, " --"))
  # }
  while(loop){
    search.loc.tmp <- search.loc
    GIC.result.tmp <- GIC.result
    for(l in 1:L){
      tmp.search <- search.ref.loc(Melody = Melody,
                                   study.l = l,
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
    search.loc <- search.loc.tmp
    if(abs(GIC.result - GIC.result.tmp) <= tol){
      loop <- FALSE
    }else if(loop.time >= 100){
      loop <- FALSE
      warning(paste0("Loop over 100 times for searching best subset of ", s))
    }else{
      GIC.result <- GIC.result.tmp
    }
  }
  return(list(tmp.result = tmp.result, GIC.result = GIC.result, search.loc = search.loc))
}

### Utility function
cal.quantile <- function(quantile_sm, prop){
  tmp_quant <- c()
  for(l in 1:length(prop)){
    tmp_quant <- c(tmp_quant, prop[l] * quantile_sm[l,2] + (1 - prop[l]) * quantile_sm[l,1])
  }
  return(tmp_quant)
}
