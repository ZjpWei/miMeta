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
#' null.obj <- melody.null.model(rel.abd = rel.abd)
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
        # sample.lst <- Sample.info[[d]]$rm.sample.idx
        # for(l in 1:length(sample.lst)){
        #   rm.sample.id <- sample.lst[[l]]
        #   if(length(rm.sample.id) == 0){
             dat[[d]] <- list(Y = Y.pool, X = X.pool)
        #  }else{
        #     dat[[d]] <- list(Y = Y.pool[-rm.sample.id,], X = X.pool[-rm.sample.id,])
        #  }
        # }
      }
    }
  }

  #=== Generate summary statistics ===#
  reg.fit.result <- reg.fit(dat = dat,
                            filter.threshold = prev.filter,
                            ref = ref,
                            parallel.core = parallel.core)

  for(d in study.ID){
    reg.fit.result[[d]]$rm.sample.idx <- rm.sample.idx[[d]]
  }
  return(reg.fit.result)
}
