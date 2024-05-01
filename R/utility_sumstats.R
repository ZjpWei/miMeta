#' @import parallel
#' @import foreach
#' @import doParallel
#' @import brglm2

######## Utility functions for generating summary statistics ########
## This function summarizes the summary statistics
Get_gamma <- function(summary.stat.null,
                      Sample.info,
                      SUB.id){

  study.ID <- names(summary.stat.null)
  inv_gamma <- list()
  for(d in study.ID){
    inv.gamma.lst <- list()
    sample.lst <- Sample.info[[d]]$rm.sample.idx
    for(l in 1:length(sample.lst)){
      rm.sample.id <- sample.lst[[l]]
      if(length(rm.sample.id) == 0){
        pp_mat <- summary.stat.null[[d]]$p
        s.i.mat <- summary.stat.null[[d]]$res
        X.sub <- summary.stat.null[[d]]$X
        N <- summary.stat.null[[d]]$N
        SUBid <- as.character(SUB.id[[d]])
      }else{
        pp_mat <- summary.stat.null[[d]]$p[-rm.sample.id,]
        s.i.mat <- summary.stat.null[[d]]$res[-rm.sample.id,]
        X.sub <- summary.stat.null[[d]]$X[-rm.sample.id,]
        N <- summary.stat.null[[d]]$N[-rm.sample.id]
        SUBid <- as.character(SUB.id[[d]])[-rm.sample.id]
      }
      K <- ncol(s.i.mat) + 1
      V.i.lst <- list()
      for(i in 1:length(N)){
        pp <- pp_mat[i,]
        mu_hat <- pp * N[i]
        if(K == 2){
          V.i.lst[[i]] <- N[i] * (pp - pp * t(pp))
        }else{
          V.i.lst[[i]] <- N[i] * (diag(pp) - pp %*% t(pp))
        }
      }
      X.name <- colnames(X.sub)
      name.cov <- paste0(colnames(X.sub), ":", as.character(sort(rep(1:ncol(s.i.mat),ncol(X.sub)))))
      uniq.SUBid <- unique(SUBid)
      nn <- length(uniq.SUBid)
      s.i.lst <- NULL
      A <- 0
      for(ll in 1:nn){
        s.i.SUB <- 0
        for(i in which(uniq.SUBid[ll] == SUBid)){
          A <- A + kronecker(V.i.lst[[i]], X.sub[i,] %*% t(X.sub[i,]))
        }
      }
      colnames(A) <- name.cov
      rownames(A) <- name.cov
      beta.name <- paste0(X.name[length(X.name)],":", as.character(1:(K-1)))
      gamma.name <- setdiff(colnames(A), beta.name)
      inv.gamma.lst[[l]] <- solve(A[gamma.name, gamma.name])
    }
    inv_gamma[[d]] <- inv.gamma.lst
  }
  return(inv_gamma)
}

Get_summary <- function(summary.stat.null,
                        covariate.interest,
                        SUB.id,
                        inv_gamma,
                        Sample.info,
                        cov.type,
                        G,
                        parallel.core,
                        verbose){

  #=== Cluster samples from same subject ===#
  study.ID <- names(summary.stat.null)

  cores <- detectCores()
  if(is.null(parallel.core)){
    parallel.core <- max(cores[1]-1, 1)
  }else{
    if(parallel.core >= cores){
      warning("The number of cores excceed the capacity.\n")
      parallel.core <- max(cores[1]-1, 1)
    }
    if(parallel.core <= 0 | as.integer(parallel.core) != parallel.core){
      stop("The number of cores must be a positive interger.\n")
    }
  }
  if(verbose){
    message(paste0("++ ",parallel.core[1], " cores are using for generating summary statistics. ++"))
  }

  cl <- makeCluster(parallel.core[1])
  registerDoParallel(cl)
  summary.stat.study <- foreach(d = study.ID, .packages = "brglm2") %dopar% {
    cov.int.lst <- Sample.info[[d]]$rm.sample.cov
    cov.int.id <- Sample.info[[d]]$rm.sample.idx
    cov.int.nm <- colnames(covariate.interest[[d]])
    feat.id <- colnames(summary.stat.null[[d]]$p)
    est.mat <- matrix(NA, nrow = length(feat.id), ncol = length(cov.int.nm),
                      dimnames = list(feat.id, cov.int.nm))
    cov.mat <- matrix(NA, nrow = length(feat.id), ncol = length(cov.int.nm),
                      dimnames = list(feat.id, cov.int.nm))
    # n.vec <- NULL
    for(cov.name in cov.int.nm){
      if(verbose){
        message("++ Construct summary statistics for study ", d, " and covariate of interest ", cov.name, ". ++")
      }
      l <- which(unlist(lapply(cov.int.lst, function(r){cov.name %in% r})))
      rm.sample.id <- cov.int.id[[l]]
      if(length(rm.sample.id) == 0){
        pp_mat <- summary.stat.null[[d]]$p
        s.i.mat <- summary.stat.null[[d]]$res
        X.sub <- summary.stat.null[[d]]$X
        N <- summary.stat.null[[d]]$N
        SUBid <- as.character(SUB.id[[d]])
      }else{
        pp_mat <- summary.stat.null[[d]]$p[-rm.sample.id,]
        s.i.mat <- summary.stat.null[[d]]$res[-rm.sample.id,]
        X.sub <- summary.stat.null[[d]]$X[-rm.sample.id,]
        N <- summary.stat.null[[d]]$N[-rm.sample.id]
        SUBid <- as.character(SUB.id[[d]])[-rm.sample.id]
      }
      K <- ncol(s.i.mat) + 1
      V.i.lst <- list()
      for(i in 1:length(N)){
        pp <- pp_mat[i,]
        # mu_hat <- pp * N[i]
        if(K == 2){
          V.i.lst[[i]] <- N[i] * (pp - pp * t(pp))
        }else{
          V.i.lst[[i]] <- N[i] * (diag(pp) - pp %*% t(pp))
        }
      }
      if(!is.numeric(covariate.interest[[d]][,cov.name])){
        stop("Covariate.interest should be numeric, please check your input.")
      }
      X.sub[,ncol(X.sub)] <- covariate.interest[[d]][rownames(X.sub),cov.name]
      K <- ncol(s.i.mat) + 1
      n <- nrow(X.sub)
      X.name <- colnames(X.sub)
      name.cov <- paste0(colnames(X.sub), ":", rep(colnames(s.i.mat), each = ncol(X.sub)))
      uniq.SUBid <- unique(SUBid)
      nn <- length(uniq.SUBid)
      s.i.lst <- NULL
      A <- 0
      for(ll in 1:nn){
        s.i.SUB <- 0
        for(i in which(uniq.SUBid[ll] == SUBid)){
          s.i.SUB <- s.i.SUB + t(kronecker(matrix(unlist(s.i.mat[i,]), ncol=1), matrix(X.sub[i,], ncol=1)))
          A <- A + kronecker(V.i.lst[[i]], X.sub[i,] %*% t(X.sub[i,]))
        }
        s.i.lst <- rbind(s.i.lst, s.i.SUB)
      }
      cov_R <- solve(A)
      colnames(cov_R) <- name.cov
      rownames(cov_R) <- name.cov
      colnames(A) <- name.cov
      rownames(A) <- name.cov
      colnames(s.i.lst) <- name.cov
      R_lst <- list(s.i.lst = s.i.lst, cov_R = cov_R)

      #=== generate covriate matrix ===#
      ests <- cov_R[c(1:(K-1))*length(X.name), c(1:(K-1))*length(X.name)] %*% colSums(s.i.lst[,c(1:(K-1))*length(X.name),drop=FALSE], na.rm = TRUE)
      rownames(ests) <- gsub(paste0(X.name[length(X.name)],":"), "", colnames(s.i.mat))

      #=== loop for ridge regularization parameter ===#
      if(cov.type == "ridge"){
        set.seed(2023)
        sample.id <- split(sample(nn), (1:nn)%%G)
        loop_mat <- matrix(0, 3, 3)
        rownames(loop_mat) <- c("lambda", "Pearson.R","signal")
        loop_mat["lambda",] <- c(0, 0.05, 0.8)
        loop_mat <- Calc.pearson(loop_mat, R_lst, sample.id, lambda, G, shrehold = 1e-2)
        lambda <- (loop_mat["lambda",loop_mat["Pearson.R",] == max(loop_mat["Pearson.R",])])[1]
        colnames(loop_mat) <- paste0("Tn", as.character(1:ncol(loop_mat)))
        if(verbose){
          cat("Tuning ridge parameter : ",lambda, ".\n")
        }
        summary$loop_mat <- loop_mat
      }else if(cov.type == "diag"){
        lambda <- 0
      }

      #=== solve GEE equation ===#
      beta.name <- paste0(X.name[length(X.name)],":", colnames(s.i.mat))
      gamma.name <- setdiff(colnames(cov_R), beta.name)
      core.U <- 0
      for(ll in 1:nrow(s.i.lst)){
        tmp.U <- s.i.lst[ll,beta.name] - A[beta.name, gamma.name] %*% inv_gamma[[d]][[l]] %*% s.i.lst[ll,gamma.name]
        core.U <- core.U + tmp.U %*% t(tmp.U)
      }
      Sigma <- cov_R[beta.name,beta.name] %*% (core.U) %*% cov_R[beta.name,beta.name]
      if(K == 2){
        Sigma_d <- sqrt(Sigma)
      }else{
        Sigma_d <- diag(sqrt(diag(Sigma)))
      }
      R <- Sigma / (sqrt(diag(Sigma)) %*% t(sqrt(diag(Sigma))))
      R_lambda <- lambda * R + (1-lambda) * diag(nrow(R))
      Sigma_lambda <- diag(Sigma_d %*% R_lambda %*% Sigma_d)

      est.mat[rownames(ests), cov.name] <- ests[,1]
      cov.mat[rownames(ests), cov.name] <- Sigma_lambda
    }
    summary.stat.study.one <- list(est = est.mat, var = cov.mat, n = n, para.id = d)

    #=== output ===#
    summary.stat.study.one
  }
  #=== stop cluster ===#
  stopCluster(cl)

  #=== reorder output ===#
  order.vec <- NULL
  for(ll in 1:length(summary.stat.study)){
    order.vec <- c(order.vec, summary.stat.study[[ll]]$para.id)
    summary.stat.study[[ll]]$para.id <- NULL
  }
  names(summary.stat.study) <- order.vec
  summary.stat.study <- summary.stat.study[study.ID]

  return(summary.stat.study)
}

reg.fit = function(dat,
                   filter.threshold = 0.1,
                   ref = NULL,
                   parallel.core = NULL,
                   verbose = FALSE){

  study.ID <- names(dat)
  #=== check maximum C.V. 1.5 output warning if too large ===#
  feature.ID <- NULL
  for(d in study.ID){
    feature.ID <- c(feature.ID, colnames(dat[[d]]$Y))
  }
  feature.ID <- sort(unique(feature.ID))
  for(d in study.ID){
    dat[[d]]$Y <- dat[[d]]$Y[,sort(intersect(feature.ID, colnames(dat[[d]]$Y)))]
  }

  K <- length(feature.ID)
  maximum.CV <- matrix(NA, nrow = length(study.ID), ncol = length(feature.ID),
                       dimnames = list(study.ID, feature.ID))
  for(d in study.ID){
    data.prop <- dat[[d]]$Y / rowSums(dat[[d]]$Y)
    data.avr.prop <- colMeans(data.prop)
    if(length(study.ID) == 1){
      maximum.CV <- apply(data.prop, 2, sd)/data.avr.prop
    }else{
      maximum.CV[d, colnames(dat[[d]]$Y)] <- apply(data.prop, 2, sd)/data.avr.prop
    }
  }
  maximum.CV[is.na(maximum.CV)] <- Inf
  if(is.null(ref)){
    if(length(study.ID) == 1){
      ref <- feature.ID[which.min(maximum.CV)]
    }else{
      ref.loc <- apply(maximum.CV, 1, which.min)
      ref <- feature.ID[ref.loc]
    }
    names(ref) <- study.ID
  }else{
    if(length(ref) == 1){
      if(ref %in% feature.ID){
        ref <- rep(ref, length(study.ID))
        names(ref) <- study.ID
      }else{
        stop(paste0(as.character(ref)," isn't in the data.\n") )
      }
    }else{
      if(length(ref) != length(study.ID)){
        stop("The reference taxa numbers don't match the study numbers.\n")
      }
      if(!all(ref %in% feature.ID)){
        stop("some taxa aren't in the data.\n")
      }
    }
    if(length(study.ID) == 1){
      feature.lst <- feature.ID[rank(maximum.CV, ties.method = "random") <= 10]
      if(!(ref %in% feature.lst)){
        warning(paste0("The specified reference ", ref, " in study ", study.ID,
                       " has a large variation. Consider choosing a reference from this list for more stable results: ",
                       paste(feature.lst, collapse = ","), "."))
      }
    }else{
      for(d in study.ID){
        feature.lst <- feature.ID[rank(maximum.CV[d,], ties.method = "random") <= 10]
        if(!(ref[d] %in% feature.lst)){
          warning(paste0("The specified reference ", ref[d], " in study ", d,
                         " has a large variation. Consider choosing a reference from this list for more stable results: ",
                         paste(feature.lst, collapse = ","), "."))
        }
      }
    }
  }

  #=== Switch the reference to the last column ===#
  data.relative <- list()
  for(d in study.ID){
    idx = c(setdiff(colnames(dat[[d]]$Y), ref[d]), ref[d])
    data.relative[[d]] <- list(Y = dat[[d]]$Y[,idx], X = dat[[d]]$X)
  }

  #=== Generate summary statistics ===#
  cores <- detectCores()
  if(is.null(parallel.core)){
    parallel.core <- max(cores[1]-1, 1)
  }else{
    if(parallel.core >= cores){
      warning("The number of cores excceed the capacity.\n")
      parallel.core <- max(cores[1]-1, 1)
    }
    if(parallel.core <= 0 | as.integer(parallel.core) != parallel.core){
      stop("The number of cores must be a positive interger.\n")
    }
  }
  if(verbose){
    message(paste0("++ ",parallel.core[1], " cores are using for generating summary statistics. ++"))
  }

  #=== Setup parallel jobs ===#
  cl <- makeCluster(parallel.core[1])
  registerDoParallel(cl)
  reg.fit <- foreach(d = study.ID, .packages = "brglm2") %dopar% {
  # reg.fit <- list()
  # for(d in study.ID){
    Y.sub <- data.relative[[d]]$Y
    X.sub <- cbind(1, data.relative[[d]]$X)
    colnames(X.sub) <- c("Intercept", paste0("V_", as.character(1:(ncol(X.sub)-1))))
    rownames(X.sub) <- rownames(Y.sub)
    feature.set.tmp <- colSums(Y.sub != 0) > filter.threshold
    feature.set.tmp[ncol(Y.sub)] <- TRUE
    Y.sub.tmp <- Y.sub[,feature.set.tmp]

    ## Check if any correlations are 1.
    cors <- cor(Y.sub.tmp)
    lo_tri <- lower.tri(cors, diag = TRUE)
    cors[lo_tri] <- 0
    if(!all(abs(cors) != 1)){
      row_col <- which(cors == 1, arr.ind = TRUE)
      row_id <- unique(row_col[,"row"])
      col_id <- unique(row_col[,"col"])
      tax.rm <- colnames(Y.sub)[row_id]
      tax.kp <- colnames(Y.sub)[col_id]
      feature.set.tmp[tax.rm] <- FALSE
      warning("Some features have high correlation in study ", d, ", Rmove features:\n",
              paste0("    ",tax.rm,"\n"))

    }
    Y.sub <- Y.sub[,feature.set.tmp]
    tmp <- GetGlm(data.beta = list(Y = Y.sub, X = X.sub), X.idx = ncol(X.sub))

    #=== summary null model ===#
    est <- t(tmp$est)
    N <- rowSums(Y.sub)
    s.i.mat <- NULL
    pp_mat <- NULL
    for(i in 1:length(N)){
      dd <- colSums(matrix(rep(X.sub[i,], ncol(Y.sub)-1), nrow = ncol(X.sub)) * t(est), na.rm = TRUE)
      pp <- c(exp(dd - max(dd)),1/exp(max(dd)))
      pp <- (pp/sum(pp))[1:(ncol(Y.sub)-1)]
      pp_mat <- rbind(pp_mat, pp)
      mu_hat <- pp * N[i]
      s.i.mat <- rbind(s.i.mat, Y.sub[i,-ncol(Y.sub)] - mu_hat)
    }
    rownames(s.i.mat) <- rownames(Y.sub)
    rownames(pp_mat) <- rownames(Y.sub)
    reg.fit.one <- list(ref = ref[d], p = pp_mat, res = s.i.mat, N = N, X = X.sub, para.id = d)

    #=== output ===#
    reg.fit.one
  }
  #=== stop cluster ===#
  stopCluster(cl)

  #=== reorder output ===#
  order.vec <- NULL
  for(ll in 1:length(reg.fit)){
    order.vec <- c(order.vec, reg.fit[[ll]]$para.id)
    reg.fit[[ll]]$para.id <- NULL
  }
  names(reg.fit) <- order.vec
  reg.fit <- reg.fit[study.ID]
  return(reg.fit)
}

## This function generates the regression coefficient estimates and covariance estimates.
GetGlm <- function(data.beta, X.idx){
  n = nrow(data.beta$Y)
  K = ncol(data.beta$Y)
  d = ncol(data.beta$X)

  if(d == 2){
    est.single <- matrix(NA, nrow = 2, ncol = K-1,
                         dimnames = list(c(":(Intercept)", ":X"),
                                           colnames(data.beta$Y)[1:(K-1)]))
  }else if(d == 3){
    est.single <- matrix(NA, nrow = d, ncol = K-1,
                         dimnames = list(c(":(Intercept)", ":X", ":XV"),
                                         colnames(data.beta$Y)[1:(K-1)]))
  }else{
    est.single <- matrix(NA, nrow = d, ncol = K-1,
                         dimnames = list(c(":(Intercept)", paste0(":XV_", 1:(d-1))),
                                         colnames(data.beta$Y)[1:(K-1)]))
  }
  dim.name <- rownames(est.single)

  #=== Score-test summary statistics ===#
  suppressWarnings(
    for(k in 1:(K-1)){
      tmp = cbind(data.beta$Y[,K], data.beta$Y[,k])
      N = rowSums(tmp)
      idx.subj = which(N>0)
      tmp.X <- data.beta$X[idx.subj,]
      X.idx.id <- colnames(tmp.X)[X.idx]
      # Test the singularity
      if(all(is.na(tmp.X[,X.idx]))){
        qrstr <- qr(tmp.X[,c(1, setdiff(1:d, c(1,X.idx)))])
        # Check pivot
        if(qrstr$pivot[1] != 1){
          stop("This reference cannot work, please select another one.")
        }
      }else{
        qrstr <- qr(tmp.X[,c(1, X.idx, setdiff(1:d, c(1,X.idx)))])
        # Check pivot
        if(!all(qrstr$pivot[c(1,2)] %in% c(1,2))){
          stop("This reference cannot work, please select another one.")
        }
      }
      # Get non-singular sub-columns
      keep.covariate <- setdiff(colnames(tmp.X)[colnames(tmp.X) %in% colnames(qrstr$qr)[1:qrstr$rank]], colnames(tmp.X)[c(1,d)])

      # Try brglmFit model
      if(d == 2){
        input.data.tmp = list(Y=tmp[idx.subj,])
        glm.out.tmp =  glm(Y ~ 1, data = input.data.tmp,
                           family = binomial(logit), method = brglm2::brglmFit, type = "AS_mean")
      }else{
        input.data.tmp = list(Y=tmp[idx.subj,],
                              X = data.beta$X[idx.subj, keep.covariate])
        glm.out.tmp =  glm(Y ~ X, data = input.data.tmp,
                           family = binomial(logit), method = brglm2::brglmFit, type = "AS_mean")
      }
      if(glm.out.tmp$converged & max(abs(glm.out.tmp$coefficients)) < 100){
        names(glm.out.tmp$coefficients) <- paste0(":", names(glm.out.tmp$coefficients))
        est.single[names(glm.out.tmp$coefficients),k] <- - glm.out.tmp$coefficients
      }else{
        if(d == 2){
          input.data.tmp = list(Y=tmp[idx.subj,])
          multinom.out.tmp = brglm2::brmultinom(Y ~ 1 , input.data.tmp,  type = "AS_mean")
        }else{
          input.data.tmp = list(Y=tmp[idx.subj,],
                                X = data.beta$X[idx.subj,-c(1, which(colnames(tmp.X) == X.idx.id))])
          multinom.out.tmp = brglm2::brmultinom(Y ~ X , input.data.tmp,  type = "AS_mean")
        }
        est.single[intersect(names(multinom.out.tmp$coefficients), dim.name),k] <-
          multinom.out.tmp$coefficients[intersect(names(multinom.out.tmp$coefficients), dim.name)]
      }
    }
  )
  summary = list(est=est.single, n=n)
  return(summary)
}

# Calculate pearson for lambda
Calc.pearson <- function(loop_mat, R_lst ,sample.id, lambda, G, shrehold = 1e-2){
  loop.end <- TRUE
  para.NLoPR <- para_NLoPR(R_lst ,sample.id, G)
  while (loop.end) {
    empt.id <- which(loop_mat["signal",] == 0)
    Q <- NULL
    for(lambda in loop_mat["lambda",empt.id]){
      Q <- c(Q, get_NLoPR(R_lst ,sample.id, lambda, G, para.lst = para.NLoPR))
    }
    loop_mat["Pearson.R",empt.id] <- Q
    loop_mat["signal",empt.id] <- rep(1, length(empt.id))

    lambda.max <- which.max(loop_mat["Pearson.R",])
    if(lambda.max == 1){
      if(abs(loop_mat["lambda", lambda.max] - loop_mat["lambda", lambda.max + 1]) <= shrehold){
        loop.end <- FALSE
      }else{
        loop_mat <- cbind(loop_mat[,1:lambda.max],
                          c(mean(loop_mat["lambda", c(lambda.max, lambda.max + 1)]), 0, 0),
                          loop_mat[,(lambda.max + 1): ncol(loop_mat)])
      }
    }else if(lambda.max == ncol(loop_mat)){
      if(abs(loop_mat["lambda", lambda.max] - loop_mat["lambda", lambda.max - 1]) <= shrehold){
        loop.end <- FALSE
      }else{
        loop_mat <- cbind(loop_mat[,1:(lambda.max-1)],
                          c(mean(loop_mat["lambda", c(lambda.max - 1, lambda.max)]), 0, 0),
                          loop_mat[,lambda.max: ncol(loop_mat)])
      }
    }else{
      if(abs(loop_mat["lambda", lambda.max] - loop_mat["lambda", lambda.max + 1]) > shrehold){
        loop_mat <- cbind(loop_mat[,1:lambda.max],
                          c(mean(loop_mat["lambda", c(lambda.max, lambda.max + 1)]), 0, 0),
                          loop_mat[,(lambda.max + 1): ncol(loop_mat)])
      }
      if(abs(loop_mat["lambda", lambda.max] - loop_mat["lambda", lambda.max - 1]) > shrehold){
        loop_mat <- cbind(loop_mat[,1:(lambda.max-1)],
                          c(mean(loop_mat["lambda", c(lambda.max - 1, lambda.max)]), 0, 0),
                          loop_mat[,lambda.max: ncol(loop_mat)])
      }
      if(sum(loop_mat["signal",] == 0) == 0){
        loop.end <- FALSE
      }
    }
  }
  return(loop_mat)
}

# Normal likelihood
para_NLoPR <- function(R_lst ,sample.id, G){
  para.lst <- list()
  for(g in 1:G){
    S.id <- sort(sample.id[[g]])
    R.id <- NULL
    for(gg in which(1:G != g)){
      R.id <- c(R.id, sample.id[[gg]])
    }
    R.id <- sort(R.id)
    sand_slash_g <- 0
    for(ii in R.id){
      sand_slash_g  <- sand_slash_g + (R_lst$cov_R %*% R_lst$s.i.lst[ii,]) %*% t(R_lst$cov_R %*% R_lst$s.i.lst[ii,])
    }
    Sigma_slash_g <- sand_slash_g/length(R.id)
    Sigma_d <- diag(sqrt(diag(Sigma_slash_g)))
    R_slash_g <- Sigma_slash_g / (sqrt(diag(Sigma_slash_g)) %*% t(sqrt(diag(Sigma_slash_g))))
    para.lst[[g]] <- list(R_slash_g = R_slash_g, Sigma_d = Sigma_d)
  }
  return(para.lst)
}

get_NLoPR <- function(R_lst ,sample.id, lambda, G, para.lst){
  Q <- 0
  for(g in 1:G){
    S.id <- sort(sample.id[[g]])
    R_slash_g <- para.lst[[g]]$R_slash_g
    Sigma_d <- para.lst[[g]]$Sigma_d
    R_lambda_slash_g <- lambda * R_slash_g + (1-lambda) * diag(nrow(R_slash_g))
    Sigma_lambda_slash_g <- Sigma_d %*% R_lambda_slash_g %*% Sigma_d
    S_g <- R_lst$cov_R %*% t(R_lst$s.i.lst[S.id,])
    Quard <- t(S_g) %*% solve(Sigma_lambda_slash_g)  %*% S_g
    Q <- Q - (length(S.id) * ( 2 *sum(log(diag(Sigma_d))) + log(det(R_lambda_slash_g)) ) + sum(diag(Quard)))/2
  }
  if(abs(Q) == Inf){
    return(-Inf)
  }else{
    return(Q)
  }
}

########################### following code is the utility functions for paired version ######################
# Get_summary <- function(summary.stat.null,
#                         covariate.interest,
#                         SUB.id,
#                         # inv_gamma,
#                         Sample.info,
#                         cov.type,
#                         G,
#                         shrehold,
#                         verbose){
#
#   #=== Cluster samples from same subject ===#
#   study.ID <- names(summary.stat.null)
#   summary.stat.study <- list()
#   for(d in study.ID){
#     cov.int.lst <- Sample.info[[d]]$rm.sample.cov
#     cov.int.id <- Sample.info[[d]]$rm.sample.idx
#     cov.int.nm <- colnames(covariate.interest[[d]])
#     feat.id <- colnames(summary.stat.null[[d]]$p)
#     est.mat <- matrix(NA, nrow = length(feat.id), ncol = length(cov.int.nm),
#                       dimnames = list(feat.id, cov.int.nm))
#     cov.mat <- matrix(NA, nrow = length(feat.id), ncol = length(cov.int.nm),
#                       dimnames = list(feat.id, cov.int.nm))
#
#     # n.vec <- NULL
#     for(cov.name in cov.int.nm){
#       if(!is.numeric(covariate.interest[[d]][,cov.name])){
#         stop("Covariate.interest should be numeric, please check your input.")
#       }
#       if(verbose){
#         message("++ Construct summary statistics for study ", d, " and covariate of interest ", cov.name, ". ++")
#       }
#       l <- which(unlist(lapply(cov.int.lst, function(r){cov.name %in% r})))
#       rm.sample.id <- cov.int.id[[l]]
#       if(length(rm.sample.id) == 0){
#         pp.mat.all <- summary.stat.null[[d]]$p
#         s.i.mat.all <- summary.stat.null[[d]]$res
#         X.sub.all <- summary.stat.null[[d]]$X
#         N.all <- summary.stat.null[[d]]$N
#         SUBid.all <- SUB.id[[d]]
#         n <- sum(rowSums(N.all)!=0)
#       }else{
#         pp.mat.all <- summary.stat.null[[d]]$p[-rm.sample.id,,drop=FALSE]
#         s.i.mat.all <- summary.stat.null[[d]]$res[-rm.sample.id,,drop=FALSE]
#         X.sub.all <- summary.stat.null[[d]]$X[-rm.sample.id,,drop=FALSE]
#         N.all <- summary.stat.null[[d]]$N[-rm.sample.id,,drop=FALSE]
#         SUBid.all <- SUB.id[[d]][-rm.sample.id]
#         n <- sum(rowSums(N.all)!=0)
#       }
#       for(k in colnames(s.i.mat.all)){
#         rm.zero.seq <- names(which(N.all[,k]!=0))
#         pp.mat <- pp.mat.all[rm.zero.seq,k]
#         s.i.mat <- s.i.mat.all[rm.zero.seq,k]
#         X.sub <- X.sub.all[rm.zero.seq,,drop=FALSE]
#         N <- N.all[rm.zero.seq,k]
#         SUBid <- SUBid.all[rm.zero.seq]
#         V.i.lst <- list()
#         for(i in 1:length(N)){
#           pp <- pp.mat[i]
#           V.i.lst[[i]] <- N[i] * pp * (1 - pp)
#         }
#
#         X.sub[,ncol(X.sub)] <- covariate.interest[[d]][rownames(X.sub),cov.name]
#         X.name <- colnames(X.sub)
#         name.cov <- paste0(colnames(X.sub), ":", rep(k, each = ncol(X.sub)))
#         uniq.SUBid <- unique(SUBid)
#         nn <- length(uniq.SUBid)
#         s.i.lst <- NULL
#         A <- 0
#         for(ll in 1:nn){
#           s.i.SUB <- 0
#           for(i in which(uniq.SUBid[ll] == SUBid)){
#             s.i.SUB <- s.i.SUB + t(kronecker(s.i.mat[i], matrix(X.sub[i,], ncol=1)))
#             A <- A + kronecker(V.i.lst[[i]], X.sub[i,] %*% t(X.sub[i,]))
#           }
#           s.i.lst <- rbind(s.i.lst, s.i.SUB)
#         }
#         cov_R <- solve(A)
#         colnames(cov_R) <- name.cov
#         rownames(cov_R) <- name.cov
#         colnames(A) <- name.cov
#         rownames(A) <- name.cov
#         colnames(s.i.lst) <- name.cov
#
#         # R_lst <- list(s.i.lst = s.i.lst, cov_R = cov_R)
#         #=== generate covriate matrix ===#
#         ests <- cov_R[length(X.name), length(X.name)] * sum(s.i.lst[,length(X.name)])
#         names(ests) <- gsub(paste0(X.name[length(X.name)],":"), "", k)
#         #=== loop for ridge regularization parameter ===#
#         # if(cov.type == "ridge"){
#         #   set.seed(2023)
#         #   sample.id <- split(sample(nn), (1:nn)%%G)
#         #   loop_mat <- matrix(0, 3, 3)
#         #   rownames(loop_mat) <- c("lambda", "Pearson.R","signal")
#         #   loop_mat["lambda",] <- c(0, 0.05, 0.8)
#         #   loop_mat <- Calc.pearson(loop_mat, R_lst, sample.id, lambda, G, shrehold)
#         #   lambda <- (loop_mat["lambda",loop_mat["Pearson.R",] == max(loop_mat["Pearson.R",])])[1]
#         #   colnames(loop_mat) <- paste0("Tn", as.character(1:ncol(loop_mat)))
#         #   if(verbose){
#         #     cat("Tuning ridge parameter : ",lambda, ".\n")
#         #   }
#         #   summary$loop_mat <- loop_mat
#         # }else if(cov.type == "diag"){
#         lambda <- 0
#         # }
#
#         #=== solve GEE equation ===#
#         beta.name <- paste0(X.name[length(X.name)],":", k)
#         gamma.name <- setdiff(colnames(cov_R), beta.name)
#         core.U <- 0
#         for(ll in 1:nrow(s.i.lst)){
#           tmp.U <- s.i.lst[ll,beta.name] - A[beta.name, gamma.name] %*% solve(A[gamma.name, gamma.name]) %*% s.i.lst[ll,gamma.name]
#           core.U <- core.U + tmp.U %*% t(tmp.U)
#         }
#         Sigma <- cov_R[beta.name,beta.name] %*% (core.U) %*% cov_R[beta.name,beta.name]
#         Sigma_d <- sqrt(Sigma)
#         R <- Sigma / (sqrt(diag(Sigma)) %*% t(sqrt(diag(Sigma))))
#         R_lambda <- lambda * R + (1-lambda) * diag(nrow(R))
#         Sigma_lambda <- diag(Sigma_d %*% R_lambda %*% Sigma_d)
#         est.mat[names(ests), cov.name] <- ests
#         cov.mat[names(ests), cov.name] <- Sigma_lambda
#       }
#     }
#     summary.stat.study[[d]] <- list(est = est.mat, var = cov.mat, n = n)
#   }
#   return(summary.stat.study)
# }

# reg.fit = function(dat,
#                    filter.threshold = 0.1,
#                    ref = NULL,
#                    parallel.core = NULL,
#                    verbose = FALSE){
#
#   study.ID <- names(dat)
#   #=== check maximum C.V. 1.5 output warning if too large ===#
#   feature.ID <- NULL
#   for(d in study.ID){
#     feature.ID <- c(feature.ID, colnames(dat[[d]]$Y))
#   }
#   feature.ID <- sort(unique(feature.ID))
#   for(d in study.ID){
#     dat[[d]]$Y <- dat[[d]]$Y[,sort(intersect(feature.ID, colnames(dat[[d]]$Y)))]
#   }
#
#   K <- length(feature.ID)
#   maximum.CV <- matrix(NA, nrow = length(study.ID), ncol = length(feature.ID),
#                        dimnames = list(study.ID, feature.ID))
#   for(d in study.ID){
#     data.prop <- dat[[d]]$Y / rowSums(dat[[d]]$Y)
#     data.avr.prop <- colMeans(data.prop)
#     if(length(study.ID) == 1){
#       maximum.CV <- apply(data.prop, 2, sd)/data.avr.prop
#     }else{
#       maximum.CV[d, colnames(dat[[d]]$Y)] <- apply(data.prop, 2, sd)/data.avr.prop
#     }
#   }
#   maximum.CV[is.na(maximum.CV)] <- Inf
#   if(is.null(ref)){
#     if(length(study.ID) == 1){
#       ref <- feature.ID[which.min(maximum.CV)]
#     }else{
#       ref.loc <- apply(maximum.CV, 1, which.min)
#       ref <- feature.ID[ref.loc]
#     }
#     names(ref) <- study.ID
#   }else{
#     if(length(ref) == 1){
#       if(ref %in% feature.ID){
#         ref <- rep(ref, length(study.ID))
#         names(ref) <- study.ID
#       }else{
#         stop(paste0(as.character(ref)," isn't in the data.\n") )
#       }
#       feature.lst <- feature.ID[rank(maximum.CV, ties.method = "random") <= 10]
#       if(!(ref %in% feature.lst)){
#         warning(paste0("The specified reference ", ref, " in study ", feature.ID,
#                        " has a large variation. Consider choosing a reference from this list for more stable results: ",
#                        paste(feature.lst, collapse = ","), "."))
#       }
#     }else{
#       if(length(ref) != length(study.ID)){
#         stop("The reference taxa numbers don't match the study numbers.\n")
#       }
#       if(!all(ref %in% feature.ID)){
#         stop("some taxa aren't in the data.\n")
#       }
#       for(d in study.ID){
#         feature.lst <- feature.ID[rank(maximum.CV[d,], ties.method = "random") <= 10]
#         if(!(ref[d] %in% feature.lst)){
#           warning(paste0("The specified reference ", ref[d], " in study ", d,
#                          " has a large variation. Consider choosing a reference from this list for more stable results: ",
#                          paste(feature.lst, collapse = ","), "."))
#         }
#       }
#     }
#   }
#
#   #=== Switch the reference to the last column ===#
#   data.relative <- list()
#   for(d in study.ID){
#     idx = c(setdiff(colnames(dat[[d]]$Y), ref[d]), ref[d])
#     data.relative[[d]] <- list(Y = dat[[d]]$Y[,idx], X = dat[[d]]$X)
#   }
#   #=== Generate summary statistics ===#
#
#   # cores <- detectCores()
#   # if(is.null(parallel.core)){
#   #   parallel.core <- cores[1]-1
#   # }else{
#   #   if(parallel.core >= cores){
#   #     warning("The number of cores excceed the capacity.\n")
#   #     parallel.core <- cores[1]-1
#   #   }
#   #   if(parallel.core <= 0 | as.integer(parallel.core) != parallel.core){
#   #     stop("The number of cores must be a positive interger.\n")
#   #   }
#   # }
#   # if(verbose){
#   #   message(paste0("++ ",parallel.core[1], " cores are using for generating summary statistics. ++"))
#   # }
#
#   #=== Setup parallel jobs ===#
#   # cl <- makeCluster(parallel.core[1])
#   # registerDoParallel(cl)
#   # reg.fit <- foreach(l = 1:L, .packages = "brglm2") %dopar% {
#   reg.fit <- list()
#   for(d in study.ID){
#     X.sub <- cbind(1, data.relative[[d]]$X)
#     colnames(X.sub) <- c("Intercept", paste0("V_", as.character(1:(ncol(X.sub)-1))))
#     rownames(X.sub) <- rownames(data.relative[[d]]$Y)
#     nonref.features <- colnames(data.relative[[d]]$Y)[1:(ncol(data.relative[[d]]$Y)-1)]
#     ref.refeature <- colnames(data.relative[[d]]$Y)[ncol(data.relative[[d]]$Y)]
#
#     s.i.mat <- NULL
#     pp_mat <- NULL
#     N_mat <- NULL
#     feature.filter <- NULL
#     for(k in nonref.features){
#       if(sum(data.relative[[d]]$Y[, k]!=0) > filter.threshold){
#         Y.sub <- data.relative[[d]]$Y[, c(k,ref.refeature)]
#         if(abs(cor(Y.sub)[1,2]) !=1){
#           tmp <- GetGlm(data.beta = list(Y = Y.sub, X = X.sub), X.idx = ncol(X.sub))
#
#           #=== summary null model ===#
#           est <- t(tmp$est)
#           N <- rowSums(Y.sub)
#           pp.tmp <- NULL
#           s.i.mat.tmp <- NULL
#           for(i in 1:length(N)){
#             dd <- sum(matrix(X.sub[i,], nrow = ncol(X.sub)) * t(est), na.rm = TRUE)
#             pp <- c(exp(dd - max(dd)),1/exp(max(dd)))
#             pp <- (pp/sum(pp))[1]
#             pp.tmp <- c(pp.tmp, pp)
#             s.i.mat.tmp <- c(s.i.mat.tmp, Y.sub[i,1] - pp * N[i])
#           }
#           N_mat <- cbind(N_mat, N)
#           pp_mat <- cbind(pp_mat, pp.tmp)
#           s.i.mat <- cbind(s.i.mat, s.i.mat.tmp)
#           feature.filter <- c(feature.filter, k)
#         }
#       }
#     }
#     rownames(s.i.mat) <- rownames(data.relative[[d]]$Y)
#     rownames(pp_mat) <- rownames(data.relative[[d]]$Y)
#     colnames(s.i.mat) <- feature.filter
#     colnames(pp_mat) <- feature.filter
#     colnames(N_mat) <- feature.filter
#     reg.fit.one <- list(ref = ref[d], p = pp_mat, res = s.i.mat, N = N_mat, X = X.sub, para.id = d)
#
#     #=== output ===#
#     reg.fit[[d]] <- reg.fit.one
#   }
#   #=== stop cluster ===#
#   # stopCluster(cl)
#
#   #=== reorder output ===#
#   #order.vec <- c()
#   for(d in study.ID){
#     #order.vec <- c(order.vec, reg.fit[[l]]$para.id)
#     reg.fit[[d]]$para.id <- NULL
#   }
#   #reg.fit <- reg.fit[order(order.vec)]
#   return(reg.fit)
# }

