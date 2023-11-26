#' @import parallel
#' @import foreach
#' @import doParallel
#' @import brglm2

######## Utility functions for generating summary statistics ########
## This function summarizes the summary statistics
Get_summary = function(Melody,
                       G = 5,
                       shrehold = 1e-2,
                       cov.type = c("diag", "ridge"),
                       parallel.core = NULL,
                       verbose = FALSE){

  cov.type <- match.arg(cov.type)
  L <- Melody$dat.inf$L
  K <- Melody$dat.inf$K

  #================ Setup parallel jobs =================#
  cores <- detectCores()
  if(is.null(parallel.core)){
    parallel.core <- cores[1]-1
  }else{
    if(parallel.core >= cores){
      warning("The number of cores excceed the capacity.\n")
      parallel.core <- cores[1]-1
    }
    if(parallel.core <= 0 | as.integer(parallel.core) != parallel.core){
      stop("The number of cores must be a positive interger.\n")
    }
  }
  cl <- makeCluster(parallel.core[1])
  registerDoParallel(cl)
  summary.stat.study <- foreach(l = 1:L, .packages = "brglm2") %dopar% {
    data.beta <- Melody$reg.fit[[l]]$data.beta
    tmp <- Melody$reg.fit[[l]]$tmp
    ref <- Melody$reg.fit[[l]]$ref
    taxa.vec <- Melody$reg.fit[[l]]$taxa.set
    summary.stats <- get_ridge_sumstat(data.beta = data.beta,
                                       summarys = tmp,
                                       G = G, shrehold = shrehold,
                                       cov.type = cov.type,
                                       verbose = verbose)

    tmp.l.taxa.name <- names(taxa.vec)[taxa.vec]
    names(summary.stats$est) <- tmp.l.taxa.name
    rownames(summary.stats$cov) <- tmp.l.taxa.name
    colnames(summary.stats$cov) <- tmp.l.taxa.name
    summary.stat.swd = list(est=summary.stats$est,
                            cov=summary.stats$cov,
                            n = summary.stats$n,
                            ref = ref,
                            idx.rev = Melody$reg.fit[[l]]$idx.rev,
                            para.id = l)
    #=== output ===#
    summary.stat.swd
  }
  #=== stop cluster ===#
  stopCluster(cl)

  #=== reorder output ===#
  taxa.set <- list()
  order.vec <- c()
  for(l in 1:L){
    order.vec <- c(order.vec, summary.stat.study[[l]]$para.id)
    summary.stat.study[[l]]$para.id <- NULL
    taxa.set[[l]] <- Melody$reg.fit[[l]]$taxa.set
  }
  summary.stat.study <- summary.stat.study[order(order.vec)]
  Melody$summary.stat.study <- summary.stat.study
  Melody$taxa.set <- taxa.set
  return(Melody)
}

reg.fit = function(Melody, SUB.id, filter.threshold = 0, ref = NULL, parallel.core = NULL, verbose = FALSE){

  dat <- Melody$dat
  L <- Melody$dat.inf$L
  K <- Melody$dat.inf$K
  study.names <- Melody$dat.inf$study.names
  #=== check maximum C.V. 1.5 output warning if too large ===#
  if(is.null(ref)){
    maximum.CV <- NULL
    for(l in 1:L){
      data.prop <- dat[[l]]$Y / rowSums(dat[[l]]$Y)
      data.avr.prop <- colMeans(data.prop)
      maximum.CV <- rbind(maximum.CV, apply(data.prop, 2, sd)/data.avr.prop)
    }
    maximum.CV[is.na(maximum.CV)] <- Inf
    colnames(maximum.CV) <- Melody$dat.inf$taxa.names
    if(L == 1){
      ref <- Melody$dat.inf$taxa.names[which.min(maximum.CV)]
    }else{
      ref.loc <- which.min(apply(maximum.CV, 2, max))
      ref <- rep(Melody$dat.inf$taxa.names[ref.loc] ,L)
    }
  }else{
    if(length(ref) == 1){
      if(ref %in% Melody$dat.inf$taxa.names){
        ref <- rep(ref, L)
      }else{
        stop(paste0(as.character(ref)," isn't in the data.\n") )
      }
    }else{
      if(length(ref) != L){
        stop("The reference taxa numbers don't match the study numbers.\n")
      }
      if(!all(ref %in% Melody$dat.inf$taxa.names)){
        stop("some taxa aren't in the data.\n")
      }
    }
  }

  #=== Switch the reference to the last column ===#
  ref.id <- match(ref, Melody$dat.inf$taxa.names)
  data.relative <- list()
  for(l in 1:L){
    idx = c(setdiff(1:K, ref.id[l]), ref.id[l])
    idx.rev = order(idx)
    data.relative[[l]] <- list(Y = dat[[l]]$Y[,idx], X = dat[[l]]$X)
  }

  #=== Generate summary statistics ===#
  if(verbose){
    message('++ Generating summary statistics. ++')
  }

  cores <- detectCores()
  if(is.null(parallel.core)){
    parallel.core <- cores[1]-1
  }else{
    if(parallel.core >= cores){
      warning("The number of cores excceed the capacity.\n")
      parallel.core <- cores[1]-1
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
  reg.fit <- foreach(l = 1:L, .packages = "brglm2") %dopar% {
    Y.sub <- data.relative[[l]]$Y
    X.sub <- cbind(1, data.relative[[l]]$X)
    colnames(X.sub) <- c("Intercept", paste0("V_", as.character(1:(ncol(X.sub)-1))))
    taxa.set.tmp <- colSums(Y.sub != 0)[1:(ncol(Y.sub)-1)] > filter.threshold
    Y.sub.tmp <- Y.sub[,c(taxa.set.tmp,TRUE)]
    ### check if any taxa's correlation are 1.
    cors <- cor(Y.sub.tmp)
    lo_tri <- lower.tri(cors, diag = TRUE)
    cors[lo_tri] <- 0
    if(!all(abs(cors) != 1)){
      row_col <- which(cors == 1, arr.ind = TRUE)
      row_id <- unique(row_col[,"row"])
      col_id <- unique(row_col[,"col"])
      tax.rm <- colnames(Y.sub)[row_id]
      tax.kp <- colnames(Y.sub)[col_id]
      taxa.set.tmp[tax.rm] <- FALSE
      warning("Some features have high correlation in study ", study.names[l], ", Rmove features:\n",
              paste0("    ",tax.rm,"\n"))

    }
    Y.sub <- Y.sub[,c(taxa.set.tmp,TRUE)]
    data.beta <- list(Y = Y.sub, X = X.sub, SUB.id = SUB.id[[l]])
    tmp <- GetGlm(data.beta = data.beta, X.idx = ncol(X.sub))
    reg.fit.one <- list(tmp = tmp, data.beta = data.beta, taxa.set = taxa.set.tmp,
                        ref = ref[l], idx.rev = idx.rev, para.id = l)

    #=== output ===#
    reg.fit.one
  }
  #=== stop cluster ===#
  stopCluster(cl)

  #=== reorder output ===#
  order.vec <- c()
  for(l in 1:L){
    order.vec <- c(order.vec, reg.fit[[l]]$para.id)
    reg.fit[[l]]$para.id <- NULL
  }
  reg.fit <- reg.fit[order(order.vec)]
  Melody$dat.inf$ref <- ref
  Melody$reg.fit <- reg.fit
  return(Melody)
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
    est.single[":X",] <- 0
  }else{
    est.single <- matrix(NA, nrow = d, ncol = K-1,
                         dimnames = list(c(":(Intercept)", paste0(":XV_", 1:(d-1))),
                                         colnames(data.beta$Y)[1:(K-1)]))
    est.single[paste0(":XV_", d-1),] <- 0
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
      qrstr <- qr(tmp.X[,c(1, X.idx, setdiff(1:d, c(1,X.idx)))])
      # Check pivot
      if(!all(qrstr$pivot[c(1,2)] %in% c(1,2))){
        stop("This reference cannot work, please select another one.\n")
      }
      # Get non-singular sub-columns
      tmp.X <- tmp.X[,colnames(tmp.X) %in% colnames(qrstr$qr)[1:qrstr$rank]]
      # Try brglmFit model
      if(d == 2){
        input.data.tmp = list(Y=tmp[idx.subj,])
        glm.out.tmp =  glm(Y ~ 1, data = input.data.tmp,
                           family = binomial(logit), method = brglm2::brglmFit, type = "AS_mean")
      }else{
        input.data.tmp = list(Y=tmp[idx.subj,],
                              X = data.beta$X[idx.subj,-c(1, which(colnames(tmp.X) == X.idx.id))])
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

get_ridge_sumstat <- function(data.beta, summarys, G, shrehold, cov.type = "diag", verbose = FALSE){

  Y.sub <- data.beta$Y
  X.sub <- data.beta$X
  K <- ncol(Y.sub)
  X.name <- colnames(X.sub)
  SUBid <- as.character(data.beta$SUB.id)
  uniq.SUBid <- unique(SUBid)
  nn <- length(uniq.SUBid)
  R_lst <- get_cov_sand(summary.stat = summarys,data.beta = data.beta, K = K)
  #=== generate covriate matrix ===#
  set.seed(2023)
  sample.id <- split(sample(nn), (1:nn)%%G)
  summary <- list()
  if(K == 2){
    summary$est <- R_lst$cov_R[c(1:(K-1))*length(X.name), c(1:(K-1))*length(X.name)] %*% sum(R_lst$s.i.lst[,c(1:(K-1))*length(X.name)])
  }else{
    summary$est <- R_lst$cov_R[c(1:(K-1))*length(X.name), c(1:(K-1))*length(X.name)] %*% colSums(R_lst$s.i.lst[,c(1:(K-1))*length(X.name)])
  }
  summary$est <- as.vector(summary$est)

  #=== loop for ridge regularization parameter ===#
  if(cov.type == "ridge"){
    loop_mat <- matrix(0, 3, 3)
    rownames(loop_mat) <- c("lambda", "Pearson.R","signal")
    loop_mat["lambda",] <- c(0, 0.05, 0.8)
    loop_mat <- Calc.pearson(loop_mat, R_lst, sample.id, lambda, G, shrehold)
    lambda <- (loop_mat["lambda",loop_mat["Pearson.R",] == max(loop_mat["Pearson.R",])])[1]
    colnames(loop_mat) <- paste0("Tn", as.character(1:ncol(loop_mat)))
    print(loop_mat[1:2,])
    if(verbose){
      cat("Tuning ridge parameter : ",lambda, ".\n")
    }
    summary$loop_mat <- loop_mat
  }else if(cov.type == "diag"){
    lambda <- 0
  }
  cov.ridge <- R_lst$cov_R
  A.ridge <- R_lst$A
  colnames(cov.ridge) <- colnames(R_lst$cov_sand)
  rownames(cov.ridge) <- rownames(R_lst$cov_sand)

  #=== solve GEE equation ===#
  beta.name <- paste0(X.name[length(X.name)],":", as.character(1:(K-1)))
  gamma.name <- setdiff(colnames(cov.ridge), beta.name)
  core.U <- 0
  for(l in 1:nrow(R_lst$s.i.lst)){
    tmp.U <- R_lst$s.i.lst[l,beta.name] - A.ridge[beta.name, gamma.name]%*%solve(A.ridge[gamma.name, gamma.name]) %*% R_lst$s.i.lst[l,gamma.name]
    core.U <- core.U + tmp.U %*% t(tmp.U)
  }
  Sigma <- cov.ridge[beta.name,beta.name] %*% (core.U) %*% cov.ridge[beta.name,beta.name]
  if(K == 2){
    Sigma_d <- sqrt(Sigma)
  }else{
    Sigma_d <- diag(sqrt(diag(Sigma)))
  }
  R <- Sigma / (sqrt(diag(Sigma)) %*% t(sqrt(diag(Sigma))))
  R_lambda <- lambda * R + (1-lambda) * diag(nrow(R))
  Sigma_lambda <- Sigma_d %*% R_lambda %*% Sigma_d
  summary$cov <- Sigma_lambda
  summary$n <- R_lst$n
  return(summary)
}

# This find output the sandwich covariance
get_cov_sand <- function(summary.stat, data.beta, K){

  n <- summary.stat$n
  est <- t(summary.stat$est)
  X.name <- colnames(data.beta$X)

  #=== Cluster samples from same subject ===#
  SUBid <- as.character(data.beta$SUB.id)
  uniq.SUBid <- unique(SUBid)
  nn <- length(uniq.SUBid)
  X <- data.beta$X
  N <- rowSums(data.beta$Y)
  Y <- data.beta$Y
  name.cov <- paste0(X.name, ":", as.character(sort(rep(1:(K-1),ncol(data.beta$X)))))
  pp_mat <- NULL
  for(i in 1:n){
    dd <- colSums(matrix(rep(X[i,], K-1), nrow = ncol(data.beta$X)) * t(est), na.rm = TRUE)
    pp <- c(exp(dd - max(dd)),1/exp(max(dd)))
    pp <- pp/sum(pp)
    pp_mat <- rbind(pp_mat, pp[1:(K-1)])
  }

  s.i.lst <- NULL
  B <- matrix(0, ncol(data.beta$X)*(K-1), ncol(data.beta$X)*(K-1))
  A <- 0
  for(l in 1:nn){
    s.i.SUB <- 0
    for(i in which(uniq.SUBid[l] == SUBid)){
      mu_hat <- pp_mat[i,] * N[i]
      s.i.SUB <- s.i.SUB + t(kronecker( matrix(Y[i,-K] - mu_hat, ncol=1), matrix(X[i,], ncol=1)))
      if(K == 2){
        V_i <- N[i] * ( pp_mat[i,] - pp_mat[i,] * t(pp_mat[i,]))
      }else{
        V_i <- N[i] * ( diag(pp_mat[i,]) - pp_mat[i,] %*% t(pp_mat[i,]) )
      }
      A <- A + kronecker(V_i, X[i,] %*% t(X[i,]))
    }
    s.i.lst <- rbind(s.i.lst, s.i.SUB)
    B <- B + s.i.lst[l,]  %*% t(s.i.lst[l,])
  }

  cov_R <- solve(A)
  A_d <- diag(sqrt(diag(A/n)))
  R_A <- A / (sqrt(diag(A)) %*% t(sqrt(diag(A))))
  colnames(cov_R) <- name.cov
  rownames(cov_R) <- name.cov
  colnames(A) <- name.cov
  rownames(A) <- name.cov
  colnames(s.i.lst) <- name.cov
  B_d <- diag(sqrt(diag(B/n)))
  R_B <- B / (sqrt(diag(B)) %*% t(sqrt(diag(B))))
  cov_sand <- cov_R %*% B %*% cov_R
  return(list(s.i.lst = s.i.lst,
              A = A,
              cov_R = cov_R,
              cov_sand = cov_sand,
              est = est,
              n = nn,
              ids = paste0(X.name[length(X.name)],":",
                           as.character(1:(K-1)))))

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
