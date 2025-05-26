#' @title Generate Melody summary statistics for individual studies.
#'
#' @description This function directly takes the output object from the "melody.null.model" function and
#' constructs summary statistics for covariates of interest in each study.
#' The output can be directly used by the function "melody.meta.summary" to perform meta-analysis.
#'
#' @param null.obj The output of function "melody.null.model".
#' @param ... See function melody. In covariate.interest and cluster,
#' the order of samples should be matched with the order in rel.abd used in the "melody.null.model" function.
#'
#' @return Output a list with each component for a study. The component includes the following elements.
#' \item{ref}{Reference feature ID for the study.}
#' \item{est}{A matrix contains the relative-abundance association coefficient estimates for the study. The row names are microbial feature IDs
#' and the column names are the covariates of interest IDs.}
#' \item{var}{A matrix contains the variances of the relative-abundance association coefficient estimates for the study.
#' The row names are microbial feature IDs
#' and the column names are the covariates of interest IDs.}
#' \item{n}{Sample size for the study.}
#'
#' @seealso \code{\link{melody.null.model}},
#' \code{\link{melody.meta.summary}},
#' \code{\link{melody}}
#'
#' @import doParallel
#' @importFrom brglm2 brglmFit
#' @export
#'
#' @examples
#' \donttest{
#' library("miMeta")
#' data("CRC_data")
#' CRC_abd <- CRC_data$CRC_abd
#' CRC_meta <- CRC_data$CRC_meta
#'
#' ########## Generate summary statistics ##########
#' rel.abd <- list()
#' covariate.interest <- list()
#' for(d in unique(CRC_meta$Study)){
#'   rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
#'   disease <- as.numeric(CRC_meta$Group[CRC_meta$Study == d] == "CRC")
#'   names(disease) <- CRC_meta$Sample_ID[CRC_meta$Study == d]
#'   covariate.interest[[d]] <- matrix(disease, ncol = 1, dimnames = list(names(disease), "disease"))
#' }
#'
#'  null.obj <- melody.null.model(rel.abd = rel.abd, ref = "Coprococcus catus [ref_mOTU_v2_4874]")
#'
#'  summary.stats <- melody.get.summary(null.obj = null.obj, covariate.interest = covariate.interest)
#' }
#'

melody.get.summary <- function(null.obj,
                               covariate.interest,
                               cluster = NULL,
                               parallel.core = NULL,
                               verbose = FALSE) {

  #=== Check input data ===#
  study.ID <- names(null.obj)
  feature.ID <- NULL
  ref <- null.obj$ref
  for(d in names(null.obj)){
    feature.ID <- c(feature.ID, colnames(null.obj[[d]]$res), null.obj[[d]]$ref)
    ref <- c(ref, null.obj[[d]]$ref)
    if(is.null(colnames(covariate.interest[[d]]))){
      stop("covariate.interest is not a matrix with column names. Please check your data.")
    }
  }
  names(ref) <- study.ID
  feature.ID <- sort(unique(feature.ID))
  K <- length(feature.ID)
  cov.names <- NULL
  for(d in study.ID){
    cov.names <- c(cov.names, colnames(covariate.interest[[d]]))
  }
  cov.names <- unique(cov.names)
  if(verbose){
    message.num <- 0
  }

  ### match rel.data, covariate.interest and covariate.adjust
  SUB.id <- list()
  for(d in study.ID){
    if(is.null(covariate.interest[[d]])){
      stop("Study IDs in rel.data and covariate.interest don't match, please check the input data.")
    }
    if(!is.matrix(covariate.interest[[d]])){
      if(!is.data.frame(covariate.interest[[d]])){
        stop("covariate.interest is not a list of matrices.\n")
      }
    }
    if(nrow(covariate.interest[[d]]) != (length(null.obj[[d]]$N) + length(null.obj[[d]]$rm.sample.idx))){
      stop("The sample size of covariate.interest is not correct, please check the input data.")
    }else{
      if(length(null.obj[[d]]$rm.sample.idx) > 0){
        covariate.interest[[d]] <- covariate.interest[[d]] %>% dplyr::slice(-null.obj[[d]]$rm.sample.idx)
      }
      rownames(covariate.interest[[d]]) <- names(null.obj[[d]]$N)
    }

    if(is.null(cluster[[d]])){
      cluster.nm <- names(null.obj[[d]]$N)
      names(cluster.nm) <- names(null.obj[[d]]$N)
      SUB.id[[d]] <- cluster.nm
    }else{
      if(!is.vector(cluster[[d]])){
        stop("cluster is not a list of vectors. \n")
      }
      if(length(cluster[[d]]) != (length(null.obj[[d]]$N) + length(null.obj[[d]]$rm.sample.idx))){
        stop("The sample size of cluster is not correct, please check the input data. \n")
      }else{
        if(length(null.obj[[d]]$rm.sample.idx) > 0){
          cluster[[d]] <- (cluster[[d]])[-null.obj[[d]]$rm.sample.idx]
        }
      }
      cluster.nm <- cluster[[d]]
      names(cluster.nm) <- names(null.obj[[d]]$N)
      SUB.id[[d]] <- cluster.nm
    }
  }

  # Transfer to data frame
  for(d in study.ID){
    covariate.interest[[d]] <- as.data.frame(covariate.interest[[d]])
  }

  #=== Sample combination in each study corresponded to covariate.interest
  Sample.info <- lapply(covariate.interest, function(d){
    if(ncol(d) > 1){
      tmp.cov <- apply(d,2,function(r){which(is.na(r))})
    }else{
      tmp.cov <- list(which(is.na(d)))
      names(tmp.cov) <- colnames(d)
    }
    if(length(tmp.cov) != 0){
      uni.tmp.cov <- unique(tmp.cov)
      rm.sample.idx <- list()
      rm.sample.cov <- list()
      for(l in 1:length(uni.tmp.cov)){
        rm.sample.idx[[l]] <- uni.tmp.cov[[l]]
        rm.sample.cov[[l]] <- names(tmp.cov)[tmp.cov %in% uni.tmp.cov[l]]
      }
      return(list(rm.sample.cov = rm.sample.cov, rm.sample.idx = rm.sample.idx))
    }else{
      return(list(rm.sample.cov = list(colnames(d)), rm.sample.idx = list(integer())))
    }
  })

  #=== Create gamma matrix ===#
  study.ID <- names(null.obj)
  inv_gamma <- list()
  for(d in study.ID){
    inv.gamma.lst <- list()
    sample.lst <- Sample.info[[d]]$rm.sample.idx
    for(l in 1:length(sample.lst)){
      rm.sample.id <- sample.lst[[l]]
      if(length(rm.sample.id) == 0){
        pp_mat <- null.obj[[d]]$p
        s.i.mat <- null.obj[[d]]$res
        X.sub <- null.obj[[d]]$X
        N <- null.obj[[d]]$N
        SUBid <- as.character(SUB.id[[d]])
      }else{
        pp_mat <- null.obj[[d]]$p[-rm.sample.id,]
        s.i.mat <- null.obj[[d]]$res[-rm.sample.id,]
        X.sub <- null.obj[[d]]$X[-rm.sample.id,]
        N <- null.obj[[d]]$N[-rm.sample.id]
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

  #=== Generate summary statistics ===#
  if(verbose){
    message("++ Construct summary statistics. ++")
  }

  #=== Cluster samples from same subject ===#
  study.ID <- names(null.obj)

  # cores <- detectCores()
  # if(is.null(parallel.core)){
  #   parallel.core <- max(cores[1]-1, 1)
  # }else{
  #   if(parallel.core >= cores){
  #     warning("The number of cores excceed the capacity.\n")
  #     parallel.core <- max(cores[1]-1, 1)
  #   }
  #   if(parallel.core <= 0 | as.integer(parallel.core) != parallel.core){
  #     stop("The number of cores must be a positive interger.\n")
  #   }
  # }
  # if(verbose){
  #   message(paste0("++ ",parallel.core[1], " cores are using for generating summary statistics. ++"))
  # }

  # cl <- makeCluster(parallel.core[1])
  # registerDoParallel(cl)
  summary.stat.study <- list()
  for(d in study.ID){
    # summary.stat.study <- foreach(d = study.ID, .packages = "brglm2") %dopar% {
    cov.int.lst <- Sample.info[[d]]$rm.sample.cov
    cov.int.id <- Sample.info[[d]]$rm.sample.idx
    cov.int.nm <- colnames(covariate.interest[[d]])
    feat.id <- colnames(null.obj[[d]]$p)
    est.mat <- matrix(NA, nrow = length(feat.id), ncol = length(cov.int.nm),
                      dimnames = list(feat.id, cov.int.nm))
    cov.mat <- matrix(NA, nrow = length(feat.id), ncol = length(cov.int.nm),
                      dimnames = list(feat.id, cov.int.nm))

    for(cov.name in cov.int.nm){
      if(verbose){
        message("++ Construct summary statistics for study ", d, " and covariate of interest ", cov.name, ". ++")
      }
      l <- which(unlist(lapply(cov.int.lst, function(r){cov.name %in% r})))
      rm.sample.id <- cov.int.id[[l]]
      if(length(rm.sample.id) == 0){
        pp_mat <- null.obj[[d]]$p
        s.i.mat <- null.obj[[d]]$res
        X.sub <- null.obj[[d]]$X
        N <- null.obj[[d]]$N
        SUBid <- as.character(SUB.id[[d]])
      }else{
        pp_mat <- null.obj[[d]]$p[-rm.sample.id,]
        s.i.mat <- null.obj[[d]]$res[-rm.sample.id,]
        X.sub <- null.obj[[d]]$X[-rm.sample.id,]
        N <- null.obj[[d]]$N[-rm.sample.id]
        SUBid <- as.character(SUB.id[[d]])[-rm.sample.id]
      }
      K <- ncol(s.i.mat) + 1
      V.i.lst <- list()
      for(i in 1:length(N)){
        pp <- pp_mat[i,]
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
      R_lambda <- diag(nrow(R))
      Sigma_lambda <- diag(Sigma_d %*% R_lambda %*% Sigma_d)

      est.mat[rownames(ests), cov.name] <- ests[,1]
      cov.mat[rownames(ests), cov.name] <- Sigma_lambda
    }
    summary.stat.study.one <- list(est = est.mat, var = cov.mat, n = n, para.id = d)

    #=== output ===#
    summary.stat.study[[d]] <- summary.stat.study.one
  }

  #=== stop cluster ===#
  # stopCluster(cl)

  #=== reorder output ===#
  # order.vec <- NULL
  # for(ll in 1:length(summary.stat.study)){
  #   order.vec <- c(order.vec, summary.stat.study[[ll]]$para.id)
  #   summary.stat.study[[ll]]$para.id <- NULL
  # }
  # names(summary.stat.study) <- order.vec
  # summary.stat.study <- summary.stat.study[study.ID]

  #=== Get a Melody object for meta-analysis ===#
  for(d in study.ID){
    summary.stat.study[[d]]$ref <- ref[d]
  }

  return(summary.stat.study)
}
