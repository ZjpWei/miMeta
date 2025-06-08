#' @title Get information from the null model.
#'
#' @description The "melody.null.model" function computes the estimated mean microbial feature proportions and residuals  of the relative abundance under the null model of
#' no associations for each study (so this function does not need the covariate.interest information).
#' The output of this function will be fed into the "melody.get.summary" function to construct summary statistics
#' for covariates of interest in each study.
#'
#' @param ... See function melody.
#'
#' @return Output a list with each component for a study. The component includes the following elements.
#' \item{ref}{Reference feature ID for the study.}
#' \item{p}{A matrix of the estimated mean microbial feature proportions for the study.}
#' \item{res}{A matrix of the residuals of the relative abundance for the study.}
#' \item{N}{A vector of microbiome sequencing depth for the study.}
#' \item{X}{A data frame of the cleaned covariate.adjust for the study.}
#' \item{rm.sample.idx}{The index of the removed samples for the study.}
#'
#' @seealso \code{\link{melody.get.summary}},
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
#' for(d in unique(CRC_meta$Study)){
#'   rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
#' }
#'
#' null.obj <- melody.null.model(rel.abd = rel.abd, ref = "Coprococcus catus [ref_mOTU_v2_4874]")
#' }
#'

melody.null.model <- function(rel.abd,
                              covariate.adjust = NULL,
                              depth.filter = 0,
                              prev.filter = 0.1,
                              ref = NULL,
                              parallel.core = NULL){

  #=== Check input data ===#
  study.ID <- names(rel.abd)
  if(is.null(study.ID)){
    stop("Please check the study name in rel.data.")
  }else{
    if(length(study.ID) != length(unique(study.ID))){
      stop("Multiple studies has the same names in rel.data, please check the input data.")
    }
  }

  ### match rel.data and covariate.adjust
  for(d in study.ID){
    if(!all(!is.na(rel.abd[[d]]))){
      stop("Detect NA in relative abundant counts.\n")
    }
    if(min(rel.abd[[d]]) < 0){
      stop("Detect negative value in relative abundant counts.\n")
    }
    if(!is.null(covariate.adjust[[d]])){
      if(!is.data.frame(covariate.adjust[[d]])){
        stop("covariate.adjust is not a list of data frames.\n")
      }
      if(nrow(covariate.adjust[[d]]) != nrow(rel.abd[[d]])){
        stop("The sample size doesn't match between rel.abd and covariate.adjust, please check the input data.")
      }else{
        rownames(covariate.adjust[[d]]) <- rownames(rel.abd[[d]])
      }
      if(!all(!is.na(covariate.adjust[[d]]))){
        stop("Detect NA in covariate.adjust.\n")
      }
    }
  }

  #=== Filter the samples using depth.filter ===#
  rm.sample.idx <- list()
  for(d in study.ID){
    depth.kp <- which(rowSums(rel.abd[[d]]) > depth.filter)
    rm.sample.idx[[d]] <- which(rowSums(rel.abd[[d]]) <= depth.filter)
    rel.abd[[d]] <- rel.abd[[d]][depth.kp,]
    if(!is.null(covariate.adjust[[d]])){
      covariate.adjust[[d]] <- covariate.adjust[[d]] %>% dplyr::slice(depth.kp)
    }
  }

  #=== Match samples in relative abundant counts and sample data ===#
  dat <- list()
  for(d in study.ID){
    Y.pool <- rel.abd[[d]]
    if(length(Y.pool) == 0){
      warning(paste0("Less than 20 samples in study ", d, ", the summary statistics is not stable.",
                     " Remove study ", d, "\n"))
    }else if(nrow(Y.pool) < 20){
      warning(paste0("Less than 20 samples in study ", d, ", the summary statistics is not stable.",
                     " Remove study ", d, "\n"))
    }else{
      X.pool <- matrix(NA, nrow = nrow(Y.pool))
      if(is.null(covariate.adjust)){
        dat[[d]] <- list(Y = Y.pool, X = X.pool)
      }else{
        cov.adjust <- covariate.adjust[[d]]
        for(cov_name in colnames(cov.adjust)){
          if(all(!is.na(cov.adjust[[cov_name]]))){
            if(is.factor(cov.adjust[[cov_name]])){
              class( cov.adjust[cov_name])
              dummys <- as.data.frame(model.matrix(formula(paste("~", cov_name)), data = cov.adjust[cov_name]))
              X.pool <- as.matrix(cbind(dummys[,-1], X.pool))
            }else if(is.character(cov.adjust[[cov_name]])){
              dummys <- as.data.frame(model.matrix(formula(paste("~", cov_name)), data = cov.adjust[cov_name]))
              X.pool <- as.matrix(cbind(dummys[,-1], X.pool))
            }else{
              X.pool <- cbind(cov.adjust[[cov_name]], X.pool)
            }
          }else{
            warning(paste0("NA presents in `", cov_name, "`, please check sample.data."))
          }
        }
        dat[[d]] <- list(Y = Y.pool, X = X.pool)
      }
    }
  }
  #=== Generate summary statistics ===#
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
  maximum.CV <- matrix(NA, nrow = length(study.ID), ncol = length(feature.ID), dimnames = list(study.ID, feature.ID))
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

  #=== Setup parallel jobs ===#
  # cl <- makeCluster(parallel.core[1])
  # registerDoParallel(cl)
  # reg.fit <- foreach(d = study.ID, .packages = "brglm2") %dopar% {
  reg.fit <- list()
  for(d in study.ID){
    Y.sub <- data.relative[[d]]$Y
    X.sub <- cbind(1, data.relative[[d]]$X)
    colnames(X.sub) <- c("Intercept", paste0("V_", as.character(1:(ncol(X.sub)-1))))
    rownames(X.sub) <- rownames(Y.sub)
    feature.set.tmp <- colMeans(Y.sub != 0) > prev.filter
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

    n = nrow(Y.sub)
    K = ncol(Y.sub)
    X.idx <- ncol(X.sub)
    if(ncol(X.sub) == 2){
      est.single <- matrix(NA, nrow = ncol(X.sub), ncol = K-1, dimnames = list(c(":(Intercept)", ":X"), colnames(Y.sub)[1:(K-1)]))
    }else if(ncol(X.sub) == 3){
      est.single <- matrix(NA, nrow = ncol(X.sub), ncol = K-1, dimnames = list(c(":(Intercept)", ":X", ":XV"), colnames(Y.sub)[1:(K-1)]))
    }else{
      est.single <- matrix(NA, nrow = ncol(X.sub), ncol = K-1, dimnames = list(c(":(Intercept)", paste0(":XV_", 1:(ncol(X.sub)-1))), colnames(Y.sub)[1:(K-1)]))
    }
    dim.name <- rownames(est.single)

    #=== Score-test summary statistics ===#
    Y_b <- matrix(0, nrow = n, ncol = K-1, dimnames = list(rownames(Y.sub), colnames(Y.sub)[1:(K-1)]))
    suppressWarnings(
      for(k in 1:(K-1)){
        tmp = cbind(Y.sub[,K], Y.sub[,k])
        N = rowSums(tmp)
        idx.subj = which(N>0)
        tmp.X <- X.sub[idx.subj,]
        X.idx.id <- colnames(tmp.X)[X.idx]
        # Test the singularity
        if(all(is.na(tmp.X[,X.idx]))){
          qrstr <- qr(tmp.X[,c(1, setdiff(1:ncol(X.sub), c(1,X.idx)))])
          # Check pivot
          if(qrstr$pivot[1] != 1){
            stop("This reference cannot work, please select another one.")
          }
        }else{
          qrstr <- qr(tmp.X[,c(1, X.idx, setdiff(1:ncol(X.sub), c(1,X.idx)))])
          # Check pivot
          if(!all(qrstr$pivot[c(1,2)] %in% c(1,2))){
            stop("This reference cannot work, please select another one.")
          }
        }
        # Get non-singular sub-columns
        keep.covariate <- setdiff(colnames(tmp.X)[colnames(tmp.X) %in% colnames(qrstr$qr)[1:qrstr$rank]], colnames(tmp.X)[c(1,X.idx)])

        ## Try brglmFit model
        if(ncol(X.sub) == 2){
          input.data.tmp = list(Y=tmp[idx.subj,])
          glm.out.tmp =  glm(Y ~ 1, data = input.data.tmp, family = binomial(logit), method = brglm2::brglmFit, type = "AS_mean")
        }else{
          input.data.tmp = list(Y=tmp[idx.subj,], X = X.sub[idx.subj, keep.covariate])
          glm.out.tmp =  glm(Y ~ X, data = input.data.tmp, family = binomial(logit), method = brglm2::brglmFit, type = "AS_mean")
        }
        if(glm.out.tmp$converged & max(abs(glm.out.tmp$coefficients)) < 100){
          names(glm.out.tmp$coefficients) <- paste0(":", names(glm.out.tmp$coefficients))
          est.single[names(glm.out.tmp$coefficients),k] <- - glm.out.tmp$coefficients

          ## compute the bias part
          Qmat <- qr.Q(glm.out.tmp$qr)
          Y_b[idx.subj,k] <- 0.5 * rowSums(Qmat * Qmat)
        }else{
          ## If brglmFit cannot converge, run brmultinom
          if(ncol(X.sub) == 2){
            input.data.tmp = list(Y=tmp[idx.subj,])
            multinom.out.tmp = brglm2::brmultinom(Y ~ 1 , input.data.tmp,  type = "AS_mean")
          }else{
            input.data.tmp = list(Y=tmp[idx.subj,], X = X.sub[idx.subj,-c(1, which(colnames(tmp.X) == X.idx.id))])
            multinom.out.tmp = brglm2::brmultinom(Y ~ X , input.data.tmp,  type = "AS_mean")
          }
          est.single[intersect(names(multinom.out.tmp$coefficients), dim.name),k] <-
            multinom.out.tmp$coefficients[intersect(names(multinom.out.tmp$coefficients), dim.name)]
        }
      }
    )

    #=== summary null model ===#
    N <- rowSums(Y.sub)
    s.i.mat <- NULL
    pp_mat <- NULL
    for(i in 1:length(N)){
      dd <- colSums(matrix(rep(X.sub[i,], ncol(Y.sub)-1), nrow = ncol(X.sub)) * est.single, na.rm = TRUE)
      pp <- c(exp(dd - max(dd)),1/exp(max(dd)))
      pp <- (pp/sum(pp))[1:(ncol(Y.sub)-1)]
      pp_mat <- rbind(pp_mat, pp)
      mu_hat <- pp * N[i]
      s.i.mat <- rbind(s.i.mat, Y.sub[i,-ncol(Y.sub)] - mu_hat)
    }
    s.i.mat <- s.i.mat + Y_b
    rownames(s.i.mat) <- rownames(Y.sub)
    rownames(pp_mat) <- rownames(Y.sub)
    reg.fit.one <- list(ref = ref[d], p = pp_mat, res = s.i.mat, N = N, X = X.sub, para.id = d)

    # #=== output ===#
    # reg.fit.one
    reg.fit[[d]] <- reg.fit.one
  }

  # #=== stop cluster ===#
  # stopCluster(cl)
  # #=== reorder output ===#
  # order.vec <- NULL
  # for(ll in 1:length(reg.fit)){
  #   order.vec <- c(order.vec, reg.fit[[ll]]$para.id)
  #   reg.fit[[ll]]$para.id <- NULL
  # }
  # names(reg.fit) <- order.vec
  # reg.fit <- reg.fit[study.ID]
  # return(reg.fit)

  for(d in study.ID){
    reg.fit[[d]]$rm.sample.idx <- rm.sample.idx[[d]]
  }
  return(reg.fit)
}
