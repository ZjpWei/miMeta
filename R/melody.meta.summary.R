#' @title Meta-analyze summary statistics across studies
#'
#' @description This function takes a Melody object that contains summary statistics of all studies participating the
#' meta-analysis, then combines these summary statistics to select microbial signatures (see the Description of the function melody)
#'
#' @param Melody Melody object.
#' @param ... Same arguments as the function melody
#'
#' @return Same output as the function melody
#'
#' @seealso \code{\link{melody}},
#' \code{\link{melody.get.summary}},
#' \code{\link{melody.merge.summary}}
#'
#' @import dplyr
#' @import ggplot2
#' @export
#'
#' @examples
#' \donttest{
#'
#' data("CRC_abd")
#' data("meta")
#'
#' ########## Generate summary statistics for study FR ##########
#' meta_FR <- meta[meta$Study == "FR-CRC",]
#' CRC_abd_FR <- CRC_abd[,meta_FR$Sample_ID]
#'
#' sumstats_FR <- melody.get.summary(rel.abd = CRC_abd_FR,
#'                                  sample.data = meta_FR,
#'                                  sample.id = "Sample_ID",
#'                                  study = "Study",
#'                                  disease = "Group",
#'                                  verbose = TRUE)
#'
#' ########## Generate summary statistics for study DE ##########
#' meta_DE <- meta[meta$Study == "DE-CRC",]
#' CRC_abd_DE <- CRC_abd[,meta_DE$Sample_ID]
#'
#' sumstats_DE <- melody.get.summary(rel.abd = CRC_abd_DE,
#'                                   sample.data = meta_DE,
#'                                   sample.id = "Sample_ID",
#'                                   study = "Study",
#'                                   disease = "Group",
#'                                   verbose = TRUE)
#'
#' ########## Merge summary statistics ##########
#' sumstats_merge <- melody.merge.summary(list(sumstats_FR, sumstats_DE))
#'
#' ########## Meta-analysis ##########
#' Melody.model <- melody.meta.summary(Melody = sumstats_merge)
#' }
#'

melody.meta.summary <- function(Melody,
                                tune.path = c("gsection", "sequence"),
                                tune.size.sequence = NULL,
                                tune.size.range = NULL,
                                tune.type = c("HBIC", "BIC", "KBIC", "EBIC"),
                                ouput.best.one = TRUE,
                                tol = 1e-3,
                                verbose = FALSE){

  tune.path <- match.arg(tune.path)
  tune.type <- match.arg(tune.type)
  if(verbose){
    message("++ Checking summary statistics. ++")
  }

  Melody <- Get_lasso_pre(Melody = Melody)
  summary.stat.study <- Melody$summary.stat.study
  taxa.set <- Melody$taxa.set
  L <- Melody$dat.inf$L
  K <- Melody$dat.inf$K

  quantile_sm <- NULL
  for(l in 1:L){
    tmp_qt <- quantile(summary.stat.study[[l]]$est, probs = c(0.1, 0.9))
    quantile_sm <- rbind(quantile_sm, tmp_qt)
  }

  tune.set <- list()
  GIC.tau <- c()
  if(verbose){
    message("++ Searching best model. ++")
  }
  search.loc <- rep(-Inf, L)
  tmp.result <- NULL
  GIC.result <- Inf
  if(tune.path == "sequence"){
    if(is.null(tune.size.sequence) & is.null(tune.size.range)){
      dims <- dim(Melody$lasso.mat$X.enlarge)
      support.sizes <- 1:min(dims[1], round(dims / 2))
    }else if(!is.null(tune.size.sequence) & is.null(tune.size.range)){
      support.sizes <- tune.size.sequence
    }else if(is.null(tune.size.sequence) & !is.null(tune.size.range)){
      if(!all(tune.size.range == round(tune.size.range))){
        stop("The tune.size.range include non-interger.\n")
      }
      if(tune.size.range[1] < 0 | (tune.size.range[2] - tune.size.range[1]) < 2){
        stop("The range for tune.size.range include negative value or the second value is not large enough than the first value.\n")
      }
      support.sizes <- tune.size.range[1]:tune.size.range[2]
    }else if(!is.null(tune.size.sequence) & !is.null(tune.size.range)){
      stop("Please only provide `tune.size.range` or  `tune.size.sequence`.\n")
    }
    if(verbose == TRUE){
      message("++ Searching progress. ++")
      pb <- txtProgressBar(max=length(support.sizes), style=3)
    }
    ### Search best model by sequence
    for(s.lambda in support.sizes){
      tmp.subset <- search.subset.s(Melody = Melody,
                                    L = L,
                                    quantile_sm = quantile_sm,
                                    s.size = s.lambda,
                                    tmp.result = tmp.result,
                                    search.loc = search.loc,
                                    GIC.result = GIC.result,
                                    initial.range = 0.1,
                                    tune.type = tune.type,
                                    tol = tol,
                                    verbose = FALSE)
      if(s.lambda != sum(tmp.subset$tmp.result$mu.fit!=0)){
        tmp.subset <- search.subset.s(Melody = Melody,
                                      L = L,
                                      quantile_sm = quantile_sm,
                                      s.size = s.lambda,
                                      tmp.result = NULL,
                                      search.loc = rep(-Inf, L),
                                      GIC.result = Inf,
                                      initial.range = 0.1,
                                      tune.type = tune.type,
                                      tol = tol,
                                      verbose = FALSE)
      }
      tmp.result <- tmp.subset$tmp.result
      GIC.result <- tmp.subset$GIC.result
      search.loc <- tmp.subset$search.loc
      tune.set[[as.character(s.lambda)]] <- tmp.result
      GIC.result <- GIC.cal(result = tmp.result, tune.type = tune.type)
      GIC.tau <- c(GIC.tau, GIC.result)
      if(min(GIC.tau) == GIC.result){
        min.result <- tmp.result
      }
      if(verbose){
        setTxtProgressBar(pb, (pb$getVal()+1))
      }
    }
    if(verbose == TRUE){
      cat("\n")
    }
    if(ouput.best.one){
      coef <- min.result$mu.fit
      delta <- min.result$ref_est
      dev <- min.result$q_loss
      ic <- min.result$GIC.result - dev
      names(dev) <- as.character(sum(min.result$mu.fit))
      names(ic) <- names(dev)
      names(delta) <- paste0("<",Melody$dat.inf$study.names,">_<" ,Melody$dat.inf$ref,">")
      results.all <- list(coef = coef,
                          delta = delta,
                          dev = dev,
                          ic = ic,
                          tune.path = tune.path,
                          tune.type = tune.type)
    }else{
      coef <- NULL
      delta <- NULL
      dev <- NULL
      ic <- NULL
      for(l in 1:length(tune.set)){
        coef <- cbind(coef, tune.set[[l]]$mu.fit)
        delta <- cbind(delta, tune.set[[l]]$ref_est)
        dev <- c(dev, tune.set[[l]]$q_loss)
      }
      ic <- GIC.tau - dev
      names(dev) <- as.character(support.sizes)
      names(ic) <- names(dev)
      rownames(delta) <- paste0("<",Melody$dat.inf$study.names,">_<" ,Melody$dat.inf$ref,">")
      results.all <- list(coef = coef,
                          delta = delta,
                          dev = dev,
                          ic = ic,
                          tune.path = tune.path,
                          tune.type = tune.type)
    }
  }else if(tune.path == "gsection"){
    ### Search best model by gsection
    if(is.null(tune.size.range)){
      dims <- dim(Melody$lasso.mat$X.enlarge)
      s.1.lambda <- 1
      s.2.lambda <- min(dims[1], round(dims / 2))
      support.sizes <- 1:min(dims[1], round(dims / 2))
    }else{
      if(!all(tune.size.range == round(tune.size.range))){
        stop("The tune.size.range include non-interger.\n")
      }
      if(tune.size.range[1] < 0 | (tune.size.range[2] - tune.size.range[1]) < 2){
        stop("The range for tune.size.range include negative value or the second value is not large enough than the first value.\n")
      }
      s.1.lambda <- tune.size.range[1]
      s.2.lambda <- tune.size.range[2]
      support.sizes <- s.1.lambda : s.2.lambda
    }
    loops <- TRUE
    g.section <- 2 / (sqrt(5) + 1)
    g.sequence <- c()
    g.id <- 1
    s.left.lambda <- support.sizes[round(s.1.lambda * g.section + s.2.lambda * (1 - g.section))]
    s.right.lambda <- support.sizes[round(s.1.lambda * (1 - g.section) + s.2.lambda * g.section)]
    for(s.lambda in c(s.1.lambda, s.left.lambda, s.right.lambda, s.2.lambda)){
      tmp.subset <- search.subset.s(Melody = Melody,
                                    L = L,
                                    quantile_sm = quantile_sm,
                                    s.size = s.lambda,
                                    tmp.result = NULL,
                                    search.loc = rep(-Inf, L),
                                    GIC.result = Inf,
                                    initial.range = 0.2,
                                    tune.type = tune.type,
                                    tol = tol,
                                    verbose = verbose)
      tmp.result <- tmp.subset$tmp.result
      GIC.result <- tmp.subset$GIC.result
      search.loc <- tmp.subset$search.loc
      tune.set[[as.character(s.lambda)]] <- tmp.result
      GIC.result <- GIC.cal(result = tmp.result, tune.type = tune.type)
      GIC.tau <- c(GIC.tau, GIC.result)
      g.sequence <- c(g.sequence, s.lambda)
      g.id <- g.id + 1
    }
    s.1.GIC <- GIC.tau[1]
    s.left.GIC <- GIC.tau[2]
    s.right.GIC <- GIC.tau[3]
    s.2.GIC <- GIC.tau[4]
    while (loops) {
      s.all.GIC <- c(s.1.GIC, s.left.GIC, s.right.GIC, s.2.GIC)
      if(s.1.GIC == min(s.all.GIC) & s.left.lambda - s.1.lambda <= 1){
        min.result <- tune.set[[as.character(s.1.lambda)]]
        loops <- FALSE
      }else if(s.2.GIC == min(s.all.GIC) & s.2.lambda - s.right.lambda <= 1){
        min.result <- tune.set[[as.character(s.2.lambda)]]
        loops <- FALSE
      }else if(s.left.GIC == min(s.all.GIC) & s.left.lambda - s.1.lambda <= 1 & s.right.lambda - s.left.lambda <= 1){
        min.result <- tune.set[[as.character(s.left.lambda)]]
        loops <- FALSE
      }else if(s.right.GIC == min(s.all.GIC) & s.right.lambda - s.left.lambda <= 1 & s.2.lambda - s.right.lambda <= 1){
        min.result <- tune.set[[as.character(s.right.lambda)]]
        loops <- FALSE
      }
      if(min(s.all.GIC[1:2]) == min(s.all.GIC)){
        ### minimize GIC at s.1 or s.left
        s.2.lambda <- s.right.lambda
        s.2.GIC <- s.right.GIC
        s.right.lambda <- s.left.lambda
        s.right.GIC <- s.left.GIC
        s.left.lambda <- round(s.1.lambda * g.section + s.2.lambda * (1 - g.section))
        GIC.result <- Inf
        tmp.subset <- search.subset.s(Melody = Melody,
                                      L = L,
                                      quantile_sm = quantile_sm,
                                      s.size = s.left.lambda,
                                      tmp.result = NULL,
                                      search.loc = rep(-Inf, L),
                                      GIC.result = Inf,
                                      initial.range = 0.2,
                                      tune.type = tune.type,
                                      tol = tol,
                                      verbose = verbose)
        tmp.result <- tmp.subset$tmp.result
        GIC.result <- tmp.subset$GIC.result
        search.loc <- tmp.subset$search.loc
        tune.set[[as.character(s.left.lambda)]] <- tmp.result
        s.left.GIC <- GIC.cal(result = tmp.result, tune.type = tune.type)
        GIC.tau <- c(GIC.tau, s.left.GIC)
        g.sequence <- c(g.sequence, s.left.lambda)
      }else if(min(s.all.GIC[3:4]) == min(s.all.GIC)){
        ### minimize GIC at s.right or s.2
        s.1.lambda <- s.left.lambda
        s.1.GIC <- s.left.GIC
        s.left.lambda <- s.right.lambda
        s.left.GIC <- s.right.GIC
        s.right.lambda <- round(s.1.lambda * (1 - g.section) + s.2.lambda * g.section)
        tmp.subset <- search.subset.s(Melody = Melody,
                                      L = L,
                                      quantile_sm = quantile_sm,
                                      s.size = s.right.lambda,
                                      tmp.result = NULL,
                                      search.loc = rep(-Inf, L),
                                      GIC.result = Inf,
                                      initial.range = 0.2,
                                      tune.type = tune.type,
                                      tol = tol,
                                      verbose = verbose)
        tmp.result <- tmp.subset$tmp.result
        GIC.result <- tmp.subset$GIC.result
        search.loc <- tmp.subset$search.loc
        tune.set[[as.character(s.right.lambda)]] <- tmp.result
        s.right.GIC <- GIC.cal(result = tmp.result, tune.type = tune.type)
        GIC.tau <- c(GIC.tau, s.right.GIC)
        g.sequence <- c(g.sequence, s.right.lambda)
      }
      g.id <- g.id + 1
    }
    if(ouput.best.one){
      coef <- min.result$mu.fit
      delta <- min.result$ref_est
      names(delta) <- paste0("<",Melody$dat.inf$study.names,">_<" ,Melody$dat.inf$ref,">")
      dev <- min.result$q_loss
      ic <- min.result$GIC.result - dev
      names(dev) <- as.character(sum(min.result$mu.fit))
      names(ic) <- names(dev)
      results.all <- list(coef = coef,
                          delta = delta,
                          dev = dev,
                          ic = ic,
                          tune.path = tune.path,
                          tune.type = tune.type)
    }else{
      coef <- NULL
      delta <- NULL
      dev <- NULL
      ic <- NULL
      for(l in 1:length(tune.set)){
        coef <- cbind(coef, tune.set[[l]]$mu.fit)
        delta <- cbind(delta, tune.set[[l]]$ref_est)
        dev <- c(dev, tune.set[[l]]$q_loss)
      }
      ic <- GIC.tau - dev
      rownames(delta) <- paste0("<",Melody$dat.inf$study.names,">_<" ,Melody$dat.inf$ref,">")
      names(dev) <- as.character(g.sequence)
      names(ic) <- names(dev)
      results.all <- list(coef = coef,
                          delta = delta,
                          dev = dev,
                          ic = ic,
                          tune.path = tune.path,
                          tune.type = tune.type)
    }
  }
  if(verbose){
    if(ouput.best.one){
      taxa_tab <- data.frame(taxa = names(which(results.all$coef!=0)),
                             coef = as.numeric(results.all$coef[results.all$coef!=0]))
    }else{
      min.id <- which.min(results.all$dev + results.all$ic)
      taxa_tab <- data.frame(taxa = names(which(results.all$coef[,min.id]!=0)),
                             coef = as.numeric(results.all$coef[results.all$coef[,min.id]!=0]))
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
            legend.position = "right") +
      scale_x_discrete(position='bottom') +
      scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide="none") +
      geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed")

    # plotting

    ggsave(filename = paste0(getwd(), "/miMeta.pdf"),
           plot = ggp1,
           width = 8,
           height = 0.2 * nrow(taxa_tab))

  }

  return(results.all)
}
