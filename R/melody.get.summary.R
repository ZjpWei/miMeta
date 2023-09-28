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
#'
#' @import UpSetR
#' @import doParallel
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
                               parallel.core = NULL,
                               verbose = FALSE) {

  cov.type <- match.arg(cov.type)
  ## Check input data
  if(verbose){
    message("++ Checking input data. ++")
  }
  if(!all(!is.na(rel.abd))){
    stop("Detect NA in relative abundant counts.\n")
  }
  if(min(rel.abd) < 0){
    stop("Detect negative value in relative abundant counts.\n")
  }
  if(length(study) != 1){
    stop("Detect more than one elements in study variable.\n")
  }
  if(!(study %in% colnames(sample.data))){
    stop("study variable doesn't match any variables in sample.data.\n")
  }
  if(length(disease) != 1){
    stop("Detect more than one elemnets in disease variable.\n")
  }
  if(!(disease %in% colnames(sample.data))){
    stop("disease variable doesn't match any variables in sample.data.\n")
  }
  if(length(sample.id) != 1){
    stop("Detect more than one elemnets in sample.id variable.\n")
  }
  if(!(sample.id %in% colnames(sample.data))){
    stop("sample.id doesn't match any variables in sample.data.\n")
  }
  if(!is.null(cluster)){
    if(length(cluster) != 1){
      stop("Detect more than one elemnets in cluster variable.\n")
    }
    if(!(cluster %in% colnames(sample.data))){
      stop("cluster variable doesn't match any variables in sample.data.\n")
    }
  }
  if(!is.null(covariates)){
    if(!all(covariates %in% colnames(sample.data))){
      stop("covariates don't match variables in sample.data.\n")
    }
  }
  ### Match samples in relative abundant counts and sample data
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
      stop(paste0("Some samples in study ", study.names[l], " don't match the sample in rel.abd.\n"))
    }
    Y.pool <- t(rel.abd[,sampleID])
    if(length(Y.pool) == 0){
      warning(paste0("Less than 20 samples in study ", study.names[l], ", the summary statistics is not stable.",
                     " Remove study ", study.names[l], "\n"))
    }else if(nrow(Y.pool) < 20){
      warning(paste0("Less than 20 samples in study ", study.names[l], ", the summary statistics is not stable.",
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
  ### Filter the samples using depth.filter
  for(l in 1:L){
    depth.kp <- rowSums(dat[[l]]$Y) > depth.filter
    if(is.null(covariates)){
      dat[[l]]$X <- dat[[l]]$X[depth.kp]
    }else{
      dat[[l]]$X <- dat[[l]]$X[depth.kp,]
    }
    dat[[l]]$Y <- dat[[l]]$Y[depth.kp,]
  }

  ### Create Melody object
  summary.stat.study <- Melody$new(dat = dat)
  summary.stat.study$dat.inf$study.names <- study.id

  ### Generate summary statistics
  summary.stat.study <- reg.fit.wald(Melody = summary.stat.study,
                                     SUB.id = SUB.id,
                                     filter.threshold = prev.filter,
                                     ref = ref,
                                     parallel.core = parallel.core,
                                     verbose = verbose)

  ### Ridge regularization on covariate matrix
  summary.stat.study <- Get_summary_wald(Melody = summary.stat.study,
                                         cov.type = cov.type,
                                         parallel.core = parallel.core,
                                         verbose = verbose)

  ### Remove raw data
  summary.stat.study$dat <- NULL
  summary.stat.study$reg.fit <- NULL

  if(verbose){
    ### Generate Upset plot
    taxa.mat <- matrix(FALSE,
                       nrow = summary.stat.study$dat.inf$L,
                       ncol = summary.stat.study$dat.inf$K)
    colnames(taxa.mat) <- summary.stat.study$dat.inf$taxa.names
    rownames(taxa.mat) <- summary.stat.study$dat.inf$study.names
    for(l in 1:length(summary.stat.study$taxa.set)){
      taxa.mat[l,names(summary.stat.study$taxa.set[[l]])] <- summary.stat.study$taxa.set[[l]]
      taxa.mat[l,summary.stat.study$dat.inf$ref[l]] <- TRUE
    }

    input <- list()
    taxa.names <- summary.stat.study$dat.inf$taxa.names
    for(l in 1:ncol(taxa.mat)){
      if(!paste0(summary.stat.study$dat.inf$study.names[taxa.mat[,l]], collapse = "&") %in% names(input)){
        input[[paste0(summary.stat.study$dat.inf$study.names[taxa.mat[,l]], collapse = "&")]] <- 1
      }else{
        input[[paste0(summary.stat.study$dat.inf$study.names[taxa.mat[,l]], collapse = "&")]] <- input[[paste0(summary.stat.study$dat.inf$study.names[taxa.mat[,l]], collapse = "&")]] + 1
      }
    }

    # Plotting
    print(upset(fromExpression(input),
                keep.order=T,
                sets = summary.stat.study$dat.inf$study.names,
                nintersects = 40,
                nsets = length(input),
                order.by = "freq",
                decreasing = T,
                mb.ratio = c(0.6, 0.4),
                number.angles = 0,
                text.scale = 1.1,
                point.size = 2.8,
                line.size = 1,
                set_size.scale_max = max(rowSums(taxa.mat)) * 1.2,
                set_size.show = TRUE
    ))
  }
  return(summary.stat.study)
}
