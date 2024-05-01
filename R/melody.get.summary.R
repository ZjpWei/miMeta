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
#'   covariate.interest[[d]] <- data.frame(disease = disease)
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
    }else{
      covariate.interest[[d]] <- data.frame(covariate.interest[[d]])
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
    if(!is.data.frame(covariate.interest[[d]])){
      stop("covariate.interest is not a list of data frames.\n")
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
  inv_gamma <- Get_gamma(summary.stat.null = null.obj,
                         Sample.info = Sample.info,
                         SUB.id = SUB.id)

  # return(inv_gamma)
  #=== Generate summary statistics ===#
  if(verbose){
    message("++ Construct summary statistics. ++")
  }

  #=== Get summary statistics with given covariate.interest and ridge regularization on covariate matrix ===#
  summary.stat.study <- Get_summary(summary.stat.null = null.obj,
                                    covariate.interest,
                                    SUB.id = SUB.id,
                                    inv_gamma = inv_gamma,
                                    Sample.info = Sample.info,
                                    cov.type = "diag",
                                    G = 5,
                                    parallel.core = parallel.core,
                                    verbose = verbose)

  #=== Get a Melody object for meta-analysis ===#
  for(d in study.ID){
    summary.stat.study[[d]]$ref <- ref[d]
  }

  return(summary.stat.study)
}
