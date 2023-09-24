######## Utility file for generating summary statistics ########
## this function sumarizes the summary statistics
Get_summary_wald = function(Melody,
                            G = 5,
                            shrehold = 1e-2,
                            cov.type = c("diag", "ridge"),
                            verbose = FALSE){

  cov.type <- match.arg(cov.type)
  #######################################
  L <- Melody$dat.inf$L
  K <- Melody$dat.inf$K
  summary.stat.study = list()
  taxa.set <- list()
  for(l in 1:L){
    data.beta <- Melody$reg.fit[[l]]$data.beta
    tmp <- Melody$reg.fit[[l]]$tmp
    taxa.set[[l]] <- Melody$reg.fit[[l]]$taxa.set
    ref <- Melody$reg.fit[[l]]$ref
    summary.stats <- get_ridge_sumstat_wald(data.beta = data.beta, summarys = tmp,
                                            G = G, shrehold = shrehold, cov.type = cov.type, verbose = verbose)
    # if(verbose){
    #   message("++ Summarizing summary statistics for study ", as.character(l), ". ++")
    # }
    tmp.l.taxa.name <- names(taxa.set[[l]])[taxa.set[[l]]]
    names(summary.stats$est) <- tmp.l.taxa.name
    rownames(summary.stats$cov) <- tmp.l.taxa.name
    colnames(summary.stats$cov) <- tmp.l.taxa.name
    rownames(summary.stats$sandwich) <- tmp.l.taxa.name
    colnames(summary.stats$sandwich) <- tmp.l.taxa.name
    summary.stat.swd = list(est=summary.stats$est,
                            cov=summary.stats$cov,
                            n = summary.stats$n,
                            sandwich.cov = summary.stats$sandwich,
                            ref = ref,
                            idx.rev = Melody$reg.fit[[l]]$idx.rev)
    summary.stat.study[[l]] = summary.stat.swd
  }

  Melody$summary.stat.study <- summary.stat.study
  Melody$taxa.set <- taxa.set
  return(Melody)
}

## this function generates models
reg.fit.wald = function(Melody, SUB.id, filter.threshold = 0, ref = NULL, verbose = FALSE){

  dat <- Melody$dat
  L <- Melody$dat.inf$L
  K <- Melody$dat.inf$K
  study.names <- Melody$dat.inf$study.names
  # waring part:
  # if(verbose){
  #   message("++ Picking reference features for each study. ++")
  # }
  ### check maximum C.V. 1.5 output warning if too large
  if(is.null(ref)){
    maximum.CV <- NULL
    for(l in 1:L){
      data.prop <- dat[[l]]$Y / rowSums(dat[[l]]$Y)
      data.avr.prop <- colMeans(data.prop)
      maximum.CV <- rbind(maximum.CV, apply(data.prop, 2, sd)/data.avr.prop)
    }
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

  ### Switch the reference to the last column
  ref.id <- match(ref, Melody$dat.inf$taxa.names)
  data.relative <- list()
  for(l in 1:L){
    idx = c(setdiff(1:K, ref.id[l]), ref.id[l])
    idx.rev = order(idx)
    data.relative[[l]] <- list(Y = dat[[l]]$Y[,idx], X = dat[[l]]$X)
  }

  ### Generate summary statistics
  if(verbose){
    message('++ Generating summary statistics. ++')
  }

  taxa.set <- list()
  reg.fit <- list()
  for(l in 1:L){
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
    taxa.set[[l]] <- taxa.set.tmp
    data.beta <- list(Y = Y.sub, X = X.sub, SUB.id = SUB.id[[l]])
    # if(verbose){
    #   message(paste0("++ Fitting model for study ",  study.names[l], ". ++"))
    # }
    tmp <- GetGlm.wald(data.beta = data.beta, X.idx = ncol(X.sub))
    reg.fit[[l]] <- list(tmp = tmp, data.beta = data.beta, taxa.set = taxa.set[[l]], ref = ref[l], idx.rev = idx.rev)
  }
  Melody$dat.inf$ref <- ref
  Melody$reg.fit <- reg.fit
  return(Melody)
}

## this function output the regression coefficient estimate and covariance estimate of the coefficient.
GetGlm.wald <- function(data.beta, X.idx){

  n = nrow(data.beta$Y)
  K = ncol(data.beta$Y)
  d = ncol(data.beta$X)

  n.beta = (K-1)*d
  par.interest.index.beta = kronecker(1,((0:(K-2))*d)) + X.idx
  est.single = NULL
  for(k in 1:(K-1)){
    tmp = cbind(data.beta$Y[,K], data.beta$Y[,k])
    N = rowSums(tmp)
    idx.subj = which(N>0)
    input.data.tmp = list(Y=tmp[idx.subj,], X = data.beta$X[idx.subj,-1])
    multinom.out.tmp = brglm2::brmultinom(Y ~ X , input.data.tmp, type = "AS_mean")
    est.single = c(est.single,  multinom.out.tmp$coefficients[-c(1:(length(multinom.out.tmp$coefficients) - d))] )

  }
  est.single <- matrix(est.single, nrow = d)
  summary = list(est=est.single, n=n)
  return(summary)
}

## this function output the sandwich covariate matrix.
get_cov_sand_wald <- function(summary.stat, data.beta, K){

  n <- summary.stat$n
  est <- t(summary.stat$est)
  X.name <- colnames(data.beta$X)
  colnames(est) <- X.name
  ## cluster samples from same subject
  SUBid <- as.character(data.beta$SUB.id)
  uniq.SUBid <- unique(SUBid)
  nn <- length(uniq.SUBid)

  X <- data.beta$X
  N <- rowSums(data.beta$Y)
  Y <- data.beta$Y

  name.cov <- paste0(X.name, ":", as.character(sort(rep(1:(K-1),ncol(data.beta$X)))))
  pp_mat <- NULL
  for(i in 1:n){
    dd <- colSums(matrix(rep(X[i,],K-1), nrow = ncol(data.beta$X)) * t(est))
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
      V_i <- N[i] * ( diag(pp_mat[i,]) - pp_mat[i,] %*% t(pp_mat[i,]) )
      A <- A + kronecker(V_i, X[i,] %*% t(X[i,]))
    }
    s.i.lst <- rbind(s.i.lst, s.i.SUB)
    B <- B + s.i.lst[l,]  %*% t(s.i.lst[l,])
  }
  #return(A)
  cov_R <- solve(A)
  A_d <- diag(sqrt(diag(A/n)))
  R_A <- A / (sqrt(diag(A)) %*% t(sqrt(diag(A))))
  colnames(cov_R) <- name.cov
  rownames(cov_R) <- name.cov

  B_d <- diag(sqrt(diag(B/n)))
  R_B <- B / (sqrt(diag(B)) %*% t(sqrt(diag(B))))

  cov_sand <- cov_R %*% B %*% cov_R
  return(list(s.i.lst = s.i.lst, A = A,
              R_A = R_A, R_B = R_B, B_d = B_d, A_d = A_d,
              cov_R = cov_R, cov_sand = cov_sand, B = B,
              est = est, n = nn,
              ids = paste0(X.name[length(X.name)],":", as.character(1:(K-1)))))

}

### utility function for ridge regularization
para_NLoPR_wald <- function(R_lst ,sample.id, G){
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

### utility function for ridge regularization
get_NLoPR_wald <- function(R_lst ,sample.id, lambda, G, para.lst){
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

### utility function for ridge regularization
Calc.pearson.wald <- function(loop_mat, R_lst ,sample.id, lambda, G, shrehold = 1e-2){
  loop.end <- TRUE
  para.NLoPR <- para_NLoPR_wald(R_lst ,sample.id, G)
  while (loop.end) {
    empt.id <- which(loop_mat["signal",] == 0)
    Q <- NULL
    for(lambda in loop_mat["lambda",empt.id]){
      Q <- c(Q, get_NLoPR_wald(R_lst ,sample.id, lambda, G, para.lst = para.NLoPR))
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

### utility function for generating summary statistics
get_ridge_sumstat_wald <- function(data.beta, summarys, G, shrehold, cov.type = "ridge", verbose = FALSE){

  Y.sub <- data.beta$Y
  X.sub <- data.beta$X
  K <- ncol(Y.sub)
  X.name <- colnames(X.sub)
  SUBid <- as.character(data.beta$SUB.id)
  uniq.SUBid <- unique(SUBid)
  nn <- length(uniq.SUBid)
  R_lst <- get_cov_sand_wald(summary.stat = summarys,data.beta = data.beta, K = K)
  set.seed(2022)
  sample.id <- split(sample(nn), (1:nn)%%G)
  summary <- list()
  summary$est <- c(R_lst$est[,X.name[length(X.name)]])
  summary$est <- as.vector(summary$est)
  ### search best lambda
  if(cov.type == "ridge"){
    loop_mat <- matrix(0, 3, 3)
    rownames(loop_mat) <- c("lambda", "Pearson.R","signal")
    loop_mat["lambda",] <- c(0, 0.05, 0.8)
    loop_mat <- Calc.pearson.wald(loop_mat, R_lst, sample.id, lambda, G)
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
  sand <- 0
  for(ii in 1:nn){
    sand  <- sand + (R_lst$cov_R %*% R_lst$s.i.lst[ii,]) %*% t(R_lst$cov_R %*% R_lst$s.i.lst[ii,])
  }
  Sigma <- sand/nn
  Sigma_d <- diag(sqrt(diag(Sigma)))
  R <- Sigma / (sqrt(diag(Sigma)) %*% t(sqrt(diag(Sigma))))
  R_lambda <- lambda * R + (1-lambda) * diag(nrow(R))
  Sigma_lambda <- Sigma_d %*% R_lambda %*% Sigma_d
  cov.ridge <- Sigma_lambda * nn
  colnames(cov.ridge) <- colnames(R_lst$cov_sand)
  rownames(cov.ridge) <- rownames(R_lst$cov_sand)
  summary$cov <- cov.ridge[paste0(X.name[length(X.name)],":", as.character(1:(K-1))),
                           paste0(X.name[length(X.name)],":", as.character(1:(K-1)))]
  if(cov.type == "diag"){
    summary$cov <- diag(diag(summary$cov))
  }
  summary$R_lst <- R_lst
  summary$n <- R_lst$n
  summary$sandwich <- R_lst$cov_sand[paste0(X.name[length(X.name)],":", as.character(1:(K-1))),
                                     paste0(X.name[length(X.name)],":", as.character(1:(K-1)))]
  summary$A <- R_lst$A
  summary$B <- R_lst$B
  return(summary)
}
