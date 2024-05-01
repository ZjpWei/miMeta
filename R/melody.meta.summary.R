#' @title Meta-analyze summary statistics across studies
#'
#' @description This function directly takes the summary statistics output from the "melody.get.summary"
#' function and combines the summary statistics across studies for selecting microbial signatures associated with each covariate of interest.
#'
#' @param summary.stats The output of function "melody.get.summary".
#' @param ... See function melody.
#'
#' @return Same output as the function "melody".
#'
#' @seealso \code{\link{melody.null.model}},
#' \code{\link{melody.get.summary}},
#' \code{\link{melody}}
#'
#' @import dplyr
#' @import UpSetR
#' @import ggplot2
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
#' null.obj <- melody.null.model(rel.abd = rel.abd)
#'
#' summary.stats <- melody.get.summary(null.obj = null.obj, covariate.interest = covariate.interest)
#'
#' ########## Meta-analysis ##########
#' meta.result <- melody.meta.summary(summary.stats = summary.stats)
#' }
#'

melody.meta.summary <- function(summary.stats,
                                tune.path = c("gsection", "sequence"),
                                tune.size.sequence = NULL,
                                tune.size.range = NULL,
                                tune.type = c("BIC", "HBIC", "KBIC", "EBIC"),
                                output.best.one = TRUE,
                                tol = 1e-3,
                                NMAX = 20,
                                verbose = FALSE){

  tune.path <- match.arg(tune.path)
  tune.type <- match.arg(tune.type)
  study.ID <- names(summary.stats)

  feature.ID <- NULL
  cov.int.ID <- NULL
  ref <- NULL
  for(d in study.ID){
    feature.ID <- c(feature.ID, rownames(summary.stats[[d]]$est), summary.stats[[d]]$ref)
    cov.int.ID <- c(cov.int.ID, colnames(summary.stats[[d]]$est))
    ref <- c(ref, summary.stats[[d]]$ref)
  }
  cov.int.ID <- sort(unique(cov.int.ID))
  feature.ID <- sort(unique(feature.ID))
  K <- length(feature.ID)

  feature.set <- list()
  for(d in study.ID){
    taxa.vec.tmp <- rep(FALSE, K)
    names(taxa.vec.tmp) <- feature.ID
    taxa.vec.tmp[rownames(summary.stats[[d]]$est)] <- TRUE
    feature.set[[d]] <- taxa.vec.tmp
  }

  output.result <- list()
  for(cov.name in cov.int.ID){
    summary.stat.study <- list()
    study.id.tmp <- c()
    for(d in study.ID){
      if(cov.name %in% colnames(summary.stats[[d]]$est)){
        summary.stat.study[[d]] <- list(est = summary.stats[[d]]$est[,cov.name],
                                        cov = diag(summary.stats[[d]]$var[,cov.name]),
                                        n = summary.stats[[d]]$n)
        study.id.tmp <- c(study.id.tmp, d)
      }
    }

    ##======= all covariates of interest ======##
    lasso.mat <- Get_lasso_pre(summary.stat.study = summary.stat.study,
                               feature.ID = feature.ID,
                               feature.set = feature.set[study.id.tmp],
                               study.ID = study.id.tmp,
                               ref = ref[study.id.tmp],
                               K = K)

    quantile_sm <- NULL
    # delta.median <- NULL
    for(d in study.id.tmp){
      tmp_qt <- quantile(summary.stat.study[[d]]$est, probs = c(0, 1))
      quantile_sm <- rbind(quantile_sm, tmp_qt)
    }
    if(verbose){
      message("++ Search for the best model for covariate of interest ", cov.name, ". ++")
    }

    meta.analysis <- meta.analysis(summary.stat.study = summary.stat.study,
                                   lasso.mat = lasso.mat,
                                   quantile_sm = quantile_sm,
                                   cov.name = cov.name,
                                   ref = ref[study.id.tmp],
                                   tune.path = tune.path,
                                   tune.size.sequence = tune.size.sequence,
                                   tune.size.range = tune.size.range,
                                   tune.type = tune.type,
                                   output.best.one = output.best.one,
                                   tol = tol,
                                   NMAX = NMAX,
                                   verbose = verbose)

    # return(meta.analysis)
    if(meta.analysis$boud.hit & is.null(tune.size.sequence) & is.null(tune.size.range)){
      if(tune.path == "gsection"){
        meta.analysis <- meta.analysis(summary.stat.study = summary.stat.study,
                                       lasso.mat = lasso.mat,
                                       quantile_sm = quantile_sm,
                                       cov.name = cov.name,
                                       ref = ref[study.id.tmp],
                                       tune.path = tune.path,
                                       tune.size.sequence = tune.size.sequence,
                                       tune.size.range = c(round(K/2), (ncol(lasso.mat$X.enlarge)-1)),
                                       tune.type = tune.type,
                                       output.best.one = output.best.one,
                                       tol = tol,
                                       NMAX = NMAX,
                                       verbose = verbose)
      }else{
        meta.analysis <- meta.analysis(summary.stat.study = summary.stat.study,
                                       lasso.mat = lasso.mat,
                                       quantile_sm = quantile_sm,
                                       cov.name = cov.name,
                                       ref = ref[study.id.tmp],
                                       tune.path = tune.path,
                                       tune.size.sequence = c(round(K/2):(ncol(lasso.mat$X.enlarge)-1)),
                                       tune.size.range = tune.size.range,
                                       tune.type = tune.type,
                                       output.best.one = output.best.one,
                                       tol = tol,
                                       NMAX = NMAX,
                                       verbose = verbose)
      }
    }
    output.result[[cov.name]] <- meta.analysis$results.all
  }

  if(verbose){
    taxa.mat <- matrix(FALSE, nrow = length(study.ID), ncol = K,
                       dimnames = list(study.ID, feature.ID))
    for(d in study.ID){
      taxa.mat[d,rownames(summary.stats[[d]]$est)] <- TRUE
    }

    if(verbose & length(study.ID) > 1){
      # Generate Upset plot
      input <- list()
      for(l in 1:ncol(taxa.mat)){
        if(paste0(study.ID[taxa.mat[,l]], collapse = "&") != ""){
          if(!paste0(study.ID[taxa.mat[,l]], collapse = "&") %in% names(input)){
            input[[paste0(study.ID[taxa.mat[,l]], collapse = "&")]] <- 1
          }else{
            input[[paste0(study.ID[taxa.mat[,l]], collapse = "&")]] <- input[[paste0(study.ID[taxa.mat[,l]], collapse = "&")]] + 1
          }
        }
      }

      # Plotting
      print(upset(fromExpression(input),
                  keep.order=T,
                  sets = study.ID,
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

    selected.num <- sort(unlist(lapply(output.result, function(d){sum(d$coef!=0)})), decreasing = TRUE)
    top.cov.name <- names(selected.num)[1:min(4, length(selected.num))]

    for(cov.name in top.cov.name){
      if(sum(output.result[[cov.name]]$coef!=0) > 0){
        if(output.best.one){
          taxa_tab <- data.frame(taxa = names(which(output.result[[cov.name]]$coef!=0)),
                                 coef = as.numeric(output.result[[cov.name]]$coef[output.result[[cov.name]]$coef!=0]))
        }else{
          min.id <- which.min(output.result[[cov.name]]$dev + output.result[[cov.name]]$ic)
          taxa_tab <- data.frame(taxa = names(which(output.result[[cov.name]]$coef[,min.id]!=0)),
                                 coef = as.numeric(output.result[[cov.name]]$coef[output.result[[cov.name]]$coef[,min.id]!=0,min.id]))
        }

        ggp1 <- taxa_tab %>% arrange(coef) %>%
          mutate(taxa = factor(taxa, levels = taxa)) %>%
          ggplot() + geom_point(aes(x= factor(taxa), y= coef)) +
          theme_classic() + coord_flip() + ylab("coef") +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 8),
                axis.text.x = element_text(size = 10),
                axis.title.y = element_blank(),
                axis.title.x = element_text(size = 10),
                panel.border = element_rect(colour = "black", fill=NA),
                legend.position = "right",
                plot.title = element_text(hjust = 0.5)) +
          ggtitle(paste0("Absolute-abundance coefficient estimates of the selected microbial features for ", cov.name)) +
          scale_x_discrete(position='bottom') +
          scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide="none") +
          geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed")

        # plotting
        print(ggp1)
      }
    }
  }
  return(output.result)
}
