#' @title Generate Melody summary statistics for one or multiple studies.
#'
#' @description The melody.get.summary, melody.merge.summary, and melody.meta.summary are three functions that represent individual components in the melody pipeline.
#' The function melody.get.summary takes individual-level data from one or multiple studies and outputs summary statistics for the studies in a Melody object.
#' The output can be directly used by the functions melody.merge.summary and melody.meta.summary.
#'
#' @param ... Same arguments as the function melody
#'
#' @return A Melody class.
#' \item{summary.stat.study}{A list includes summary statistics for each study.}
#' \item{dat.inf}{A list includes study number; reference feature for each study; feature number and feature names.}
#' \item{taxa.set}{A list includes the which features are included or not in each study.}
#'
#' @seealso \code{\link{melody}},
#' \code{\link{melody.meta.summary}},
#' \code{\link{melody.merge.summary}}
#' @export
#'
#' @examples
#' \donttest{
#' library("miMeta")
#' data("CRC_abd")
#' data("meta")
#'
#'########## Generate summary statistics for study FR ##########
#' meta_FR <- meta[meta$Study == "FR-CRC",]
#' CRC_abd_FR <- CRC_abd[,meta_FR$Sample_ID]
#'
#' sumstats <- melody.get.summary(rel.abd = CRC_abd_FR,
#'                                sample.data = meta_FR,
#'                                sample.id = "Sample_ID",
#'                                study = "Study",
#'                                disease = "Group",
#'                                verbose = TRUE)
#' }

melody.get.summary <- function(rel.abd,
                               sample.data,
                               sample.id,
                               study,
                               disease,
                               covariates = NULL,
                               cluster = NULL,
                               depth.filter = 0,
                               prev.filter = 0,
                               ref = NULL,
                               cov.type = c("diag", "ridge"),
                               keep.raw.data = FALSE,
                               verbose = FALSE) {

  cov.type <- match.arg(cov.type)
  ## Check input
  if(!all(!is.na(rel.abd))){
    stop("The relative abundant data cannot include NA.\n")
  }
  if(min(rel.abd) < 0){
    stop("The relative abundant data cannot be negative.\n")
  }
  if(length(study) != 1){
    stop("study variable is not a character.\n")
  }
  if(!(study %in% colnames(sample.data))){
    stop("study variable doesn't match any variables in sample.data. Please check your sample.data.\n")
  }
  if(length(disease) != 1){
    stop("disease variable is not a character.\n")
  }
  if(!(disease %in% colnames(sample.data))){
    stop("disease variable doesn't match any variables in sample.data. Please check your sample.data.\n")
  }
  if(!(sample.id %in% colnames(sample.data))){
    stop("sample.id doesn't match any variables in sample.data. Please check your sample.data.\n")
  }
  if(!is.null(cluster)){
    if(length(cluster) != 1){
      stop("cluster variable is not a character.\n")
    }
    if(!(cluster %in% colnames(sample.data))){
      stop("cluster variable doesn't match any variables in sample.data. Please check your sample.data.\n")
    }
  }
  if(!is.null(covariates)){
    if(!all(covariates %in% colnames(sample.data))){
      stop("covariates don't match variables in sample.data. Please check your sample.data.\n")
    }
  }
  ### Match relative abundant data and sample data
  study.names <- as.character(unique(sample.data[[study]]))
  IDs <- sample.data[[sample.id]]
  L <- length(study.names)
  dat <- list()
  SUB.id <- list()
  study.id <- c()
  k <- 1
  for(l in 1:L){
    sampleID <- IDs[study.names[l] == sample.data[[study]]]
    if(!all(sampleID %in% colnames(rel.abd))){
      stop(paste0("Some samples in study ", study.names[l], " don't match the sample in rel.abd. Please check your input data.\n"))
    }
    Y.pool <- t(rel.abd[,sampleID])
    if(length(Y.pool) == 0){
      stop(paste0("No sample in study ", study.names[l], " are included in rel.abd. Please check your sample ID in rel.abd and sample.data.\n"))
    }else if(nrow(Y.pool) < 20){
      warning(paste0("Less than 20 samples in study ", study.names[l], ", the summary statistics for this study may not stable.",
                     " Remove study ", study.names[l], "\n"))
    }else{
      study.id <- c(study.id, study.names[l])
      X.pool <- (sample.data[[disease]])[study.names[l] == sample.data[[study]]]
      if(is.character(X.pool)){
        X.pool <- factor(X.pool)
      }
      if(is.factor(X.pool)){
        X.pool <- as.numeric(X.pool)
        X.pool <- X.pool - min(X.pool)
      }
      if(is.null(cluster)){
        SUB.id[[k]] <- rownames(Y.pool)
      }else{
        SUB.id[[k]] <- (sample.data[[cluster]])[study.names[l] == sample.data[[study]]]
      }
      if(is.null(covariates)){
        dat[[k]] <- list(Y = Y.pool, X = X.pool)
      }else{
        for(cov_name in covariates){
          if(is.factor(sample.data[[cov_name]])){
            dummys <- as.data.frame(model.matrix(formula(paste("~", cov_name)), data = sample.data[study.names[l] == sample.data[[study]], ]))
            X.pool <- cbind(dummys[,-1], X.pool)
          }else{
            X.pool <- cbind(sample.data[[cov_name]], X.pool)
          }
        }
        dat[[k]] <- list(Y = Y.pool, X = X.pool)
      }
      k <- k + 1
    }
  }
  ### test the taxa number in all studies
  L <- length(dat)
  dat.dim <- NULL
  for(l in 1:L){
    dat.dim <- c(dat.dim, ncol(dat[[l]]$Y))
  }
  if(length(unique(dat.dim)) == 1){
    K <- unique(dat.dim)
    if(K <= 2){
      stop("Less than 3 taxa in the data.\n")
    }
  }else{
    stop("The taxa number are not same in all studies.\n")
  }

  ### taxon name should be same in each study
  if(L != 1){
    for(l in 1:(L-1)){
      if(!all(colnames(dat[[l]]$Y) == colnames(dat[[l+1]]$Y))){
        stop("The taxa name don't match in all studies, please cheack the taxa name in your data.\n")
      }
    }
  }
  ### remove low sequence depth
  for(l in 1:L){
    depth.kp <- rowSums(dat[[l]]$Y) > depth.filter
    if(is.null(covariates)){
      dat[[l]]$X <- dat[[l]]$X[depth.kp]
    }else{
      dat[[l]]$X <- dat[[l]]$X[depth.kp,]
    }
    dat[[l]]$Y <- dat[[l]]$Y[depth.kp,]
  }

  ### Create melody object
  summary.stat.study <- Melody$new(dat = dat)

  ### Generate summary statistics
  summary.stat.study <- reg.fit.wald(Melody = summary.stat.study,
                                     SUB.id = SUB.id,
                                     filter.threshold = prev.filter,
                                     ref = ref,
                                     verbose = verbose)

  ### Covariate matrix ridge regularization
  summary.stat.study <- Get_summary_wald(Melody = summary.stat.study,
                                         cov.type = cov.type,
                                         verbose = verbose)

  ### Check whether to keep input data
  if(!keep.raw.data){
    summary.stat.study$dat <- NULL
    summary.stat.study$reg.fit <- NULL
  }
  summary.stat.study$dat.inf$study.names <- study.id
  return(summary.stat.study)
}
