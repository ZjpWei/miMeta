#' @import R6

Melody <- R6Class("Melody", list(
  dat = NULL,
  dat.inf = list(),
  summary.stat.study = NULL,
  taxa.set = NULL,
  lasso.mat = list(),
  reg.fit = list(),

  initialize = function(dat) {
    self$dat <- dat
    if(is.null(dat)){
      self$dat.inf$L <- NULL
      self$dat.inf$K <- NULL
      self$dat.inf$taxa.names <- NULL
    }else{
      self$dat.inf$L <- length(dat)
      self$dat.inf$K <- ncol(dat[[1]]$Y)
      self$dat.inf$taxa.names <- colnames(dat[[1]]$Y)
    }
  }
))
