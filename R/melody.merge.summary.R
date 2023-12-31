#' @title Merge summary statistics across studies.
#'
#' @description This function takes a list of Melody objects generated by the function melody.get.summary,
#' align summary statistics by feature IDs, and output a Melody object containing summary statistics of all particpating studies.
#' The output can be directly used by the function melody.meta.summary.
#'
#' @param melody.obj.lst A list of summary statistics. Each element is a Melody object from \code{melody.get.summary}.
#' @param verbose whether to generate a plot for microbial feature overlap among studies. Default is FALSE.
#'
#' @return A Melody object that contains the summary statistics over multiple studies participating the meta-analysis.
#' \item{summary.stat.study}{A list includes summary statistics for each study.}
#' If verbose=TRUE, generate a plot for microbial feature overlap among studies.
#'
#' @seealso \code{\link{melody}},
#' \code{\link{melody.get.summary}},
#' \code{\link{melody.meta.summary}}
#'
#' @import UpSetR
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
#'                                   sample.data = meta_FR,
#'                                   sample.id = "Sample_ID",
#'                                   study = "Study",
#'                                   disease = "Group",
#'                                   verbose = TRUE)
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
#' }
#'

melody.merge.summary <- function(melody.obj.lst, verbose = FALSE){
  num.lst <- length(melody.obj.lst)
  # Warning: User should check the taxa names before merging summary statistics
  warning("Please ensure that each taxon in different studies uses the identical name.\n")
  # Create a merged melody object
  summary.stat.study.merged <- Melody$new(dat = NULL)
  merged.taxa.names <- c()
  merged.ref <- c()
  merged.study.names <- c()
  merged.L <- 0
  for(l in 1:num.lst){
    merged.taxa.names <- c(merged.taxa.names, melody.obj.lst[[l]]$dat.inf$taxa.names)
    merged.L <- merged.L + melody.obj.lst[[l]]$dat.inf$L
    merged.ref <- c(merged.ref, melody.obj.lst[[l]]$dat.inf$ref)
    merged.study.names <- c(merged.study.names, melody.obj.lst[[l]]$dat.inf$study.names)
  }
  merged.taxa.names <- unique(merged.taxa.names)
  K <- length(merged.taxa.names)
  summary.stat.study.merged$dat.inf$L <- merged.L
  summary.stat.study.merged$dat.inf$ref <- merged.ref
  summary.stat.study.merged$dat.inf$K <- K
  summary.stat.study.merged$dat.inf$taxa.names <- merged.taxa.names
  summary.stat.study.merged$dat.inf$study.names <- merged.study.names

  #=== Merge summary statistics and check the taxa names ===#
  merged.taxa.set <- list()
  merged.summary.stats <- list()
  tmp.len <- 1
  for(l in 1:num.lst){
    for(ll in 1:melody.obj.lst[[l]]$dat.inf$L){
      ### align taxa.set
      tmp.taxa.set <- rep(FALSE, K-1)
      tmp.names <- setdiff(merged.taxa.names, merged.ref[tmp.len])
      names(tmp.taxa.set) <- tmp.names
      idx.taxa.set <- match(names(melody.obj.lst[[l]]$taxa.set[[ll]]), names(tmp.taxa.set))
      tmp.taxa.set[idx.taxa.set] <- melody.obj.lst[[l]]$taxa.set[[ll]]
      merged.taxa.set[[tmp.len]] <- tmp.taxa.set
      ### align summary statistics
      est <- rep(NA, K-1)
      names(est) <- tmp.names
      cov <- matrix(NA, K-1, K-1)
      colnames(cov) <- tmp.names
      rownames(cov) <- tmp.names
      est[names(melody.obj.lst[[l]]$summary.stat.study[[ll]]$est)] <- melody.obj.lst[[l]]$summary.stat.study[[ll]]$est
      est <- est[tmp.taxa.set]
      cov[names(melody.obj.lst[[l]]$summary.stat.study[[ll]]$est),
          names(melody.obj.lst[[l]]$summary.stat.study[[ll]]$est)] <- melody.obj.lst[[l]]$summary.stat.study[[ll]]$cov
      cov <- cov[tmp.taxa.set,tmp.taxa.set]
      n <- melody.obj.lst[[l]]$summary.stat.study[[ll]]$n
      ref <- melody.obj.lst[[l]]$summary.stat.study[[ll]]$ref
      idx <- c(setdiff(1:K, match(ref, merged.taxa.names)), match(ref, merged.taxa.names))
      idx.rev = order(idx)
      merged.summary.stats[[tmp.len]] <- list(est = est, cov = cov, n = n, # sandwich.cov = sandwich.cov,
                                              ref = ref, idx.rev = idx.rev)
      tmp.len <- tmp.len + 1
    }
  }
  summary.stat.study.merged$taxa.set <- merged.taxa.set
  summary.stat.study.merged$summary.stat.study <- merged.summary.stats
  if(verbose){
    taxa.mat <- matrix(FALSE,
                       nrow = summary.stat.study.merged$dat.inf$L,
                       ncol = summary.stat.study.merged$dat.inf$K)
    colnames(taxa.mat) <- summary.stat.study.merged$dat.inf$taxa.names
    rownames(taxa.mat) <- summary.stat.study.merged$dat.inf$study.names
    for(l in 1:length(summary.stat.study.merged$taxa.set)){
      taxa.mat[l,names(summary.stat.study.merged$taxa.set[[l]])] <- summary.stat.study.merged$taxa.set[[l]]
      taxa.mat[l,summary.stat.study.merged$dat.inf$ref[l]] <- TRUE
    }

    input <- list()
    taxa.names <- summary.stat.study.merged$dat.inf$taxa.names
    for(l in 1:ncol(taxa.mat)){
      if(!paste0(summary.stat.study.merged$dat.inf$study.names[taxa.mat[,l]], collapse = "&") %in% names(input)){
        input[[paste0(summary.stat.study.merged$dat.inf$study.names[taxa.mat[,l]], collapse = "&")]] <- 1
      }else{
        input[[paste0(summary.stat.study.merged$dat.inf$study.names[taxa.mat[,l]], collapse = "&")]] <- input[[paste0(summary.stat.study.merged$dat.inf$study.names[taxa.mat[,l]], collapse = "&")]] + 1
      }
    }

    #=== Generate Upset plot ===#
    print(upset(fromExpression(input),
                keep.order=T,
                sets = summary.stat.study.merged$dat.inf$study.names,
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
  return(summary.stat.study.merged)
}
