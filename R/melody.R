#' @title meta-analysis for selecting microbial signatures
#'
#' @description Melody is a meta-analysis method designed to account for the unique features of compositional microbiome data
#' in selecting microbial signatrues.
#'
#' @author Zhoujingpeng Wei, Guanhua Chen, Zheng-Zheng Tang
#'
#' @param rel.abd Microbial feature-by-sample matrix of relative abundance counts. Must have microbial feature IDs
#' as row names and sample IDs as column names.
#' @param sample.data A dataframe. data frame of sample-level variables. Columns must include variables for study, disease,
#' covariates (if specified), and cluster (if specified). The row names must agree with column names in \code{rel.abd} matrix.
#' @param sample.id Name of the variable that defines the sample's ID in \code{sample.data}.
#' @param study Name of the variable that defines different studies in \code{sample.data}.
#' @param disease Name of the variable that defines the disease to identify microbial signature.
#' @param covariates A vector of names of the variables to adjust for in the microbiome-disease association model. Default is NULL.
#' @param cluster Name of the variable that define the sample clusters. For example, the values of this variable are subject IDs
#' if each subject has multiple correlated samples (e.g., measured in a longitudinal study). For a study only contains independent samples,
#' the value of this variable is unique for each sample (can set to sample ID) or set to NA for that study. Default is NULL
#' (all samples in all studies are independent).
#' @param depth.filter A cutoff value (>= 0) to remove samples with sequencing depth less than or equal to the cutoff. Default is 0 (not removing any samples).
#' @param prev.filter A cutoff value to remove microbial features with prevalence (proportion of nonzero observations in rel.abd matrix)
#' less than or equal to the cutoff. The cutoff value must be in the range of 0-1. Default is 0 (only removing microbial features that
#' have zero in all samples).
#' @param ref A vector with length L (Study number L>=1) or a character or NULL. The name of reference taxa/features in each study when
#' generating summary statistics by multinomial logistic regression. If ref is a character, all studies will use ref as reference taxon;
#' if ref is a vector with length L, the study l will use \code{ref_l} as reference taxon; if ref is NULL, melody will pick reference
#' taxa automatically (Check argument ref.iden to see how melody picks reference taxa).Default is NULL.
#' @param cov.type The type of covariate matrix of summary statistics in each study. Options include "ridge", "diag". If cov.type is "ridge",
#' melody will do ridge regularization to decide full covriate matrix of summary statistics, if cov.type is "diag", melody will use the diagonal covariate
#' matrix in summary statistics. Default is "diag".
#' @param tune.path Method to choose the optimal subset size in the best subset selection. If tune.path="sequence",
#' we solve the best subset selection problem for each size in tune.size.sequence. If tune.path="gsection",
#' we solve the best subset selection problem with the size ranged in tune.size.range, where the considered size
#' sequence within the range is determined by golden section search. Default is "gsection".
#' @param tune.size.sequence An integer vector containing the subset sizes to be considered. Only used when tune.path="sequence".
#' Default is  1:round(K/2), where K is number of microbial features.
#' @param tune.size.range An integer vector with two elements that define the range of the size. Default is c(1, K/2)
#' @param tune.type Type of information criterion for choosing the optimal tunning parameters. Available options are "BIC", "kBIC", and "mBIC".
#' Default is "HBIC".
#' @param tol Converge tolerance for detecting the best model. Default is 1e-3.
#' @param ouput.best.one Whether only output the best model. Default is TRUE.
#' @param verbose: whether to print verbose information. Default is FALSE. (see details in Value)
#'
#' @return If output.best.one=TRUE (default), output a list with the following components about the single best model over the subset sizes considered.
#' \item{coef}{A coefficient vector that contains the absolute-abundance coefficient estimates for the microbial features under the best subset size. The vector names are microbial features IDs.}
#' \item{delta}{A vector that contains the best values of the delta tunning parameters under the best subset size. The vector names are in the format "<study ID>_<reference microbial feature ID>"}
#' \item{dev}{The value of the deviance for the best subset size.}
#' \item{ic}{The value of the information criterion (specified in tune.type) for the best subset size.}
#' \item{best.size}{The best subset size.}
#'
#' If ouput.best.one=FALSE, output a list object with the following components for multiple best subset models, each of which is under a specific subset size considered (depends on tune.path, see details in Arguments).
#' \item{coef}{A coefficient matrix with each column contains the absolute-abundance coefficient estimates for the microbial features under a specific subset size. The row names are microbial features IDs and the column names are the subset sizes considered.}
#' \item{delta}{A matrix with each column contains the best values of the delta tunning parameters under a specific subset size. The row names are in the format "<study ID>_<reference microbial feature ID>" and the column names are the subset sizes considered.}
#' \item{dev}{A vector contains the values of the deviance for the subset sizes considered (shown as vector names).}
#' \item{ic}{A vector contains the values of the information criterion (specified in tune.type) for the subset sizes considered (shown as vector names).}
#' \item{best.size}{The best subset size.}
#'
#' If verbose=TRUE, Generate two plots and print information about the progress of meta-analysis.
#' plot for microbial feature overlap among studies: this plot shows the number of features shared among studies.
#' a plot showing the absolute-abundance coefficient estimates of the selected microbial features in the best model under the best subset size.
#'
#' @details
#' Melody first generates summary statistics (microbiome-disease association coefficient estimates and their variances) for
#' individual studies. Melody then combines the summary statistics across studies to select disease-associated microbial
#' signatures based on the average absolute-abundance association coefficients inferred from the summary statistics.
#' In particular, the selection of signature is operated through a best-subset selection.
#' There are two sets of tunning parameters in the best subset selection including the subset size
#' (i.e. the number of microbial features selected) and delta's which are introduced to recover the
#' absolute-abundance coefficients from the relative-abundance coefficients.
#'
#' @references Wei Z, Chen G, Tang ZZ. Melody identifies generalizable microbial signatures in microbiome association meta-analysis. Submitted.
#'
#' @seealso \code{\link{melody.get.summary}},
#' \code{\link{melody.meta.summary}},
#' \code{\link{melody.merge.summary}},
#'
#' @export
#'
#' @examples
#' \donttest{
#' library("miMeta")
#' data("CRC_abd")
#' data("meta")
#'
#' meta_FR <- meta[meta$Study == "FR-CRC",]
#' CRC_abd_FR <- CRC_abd[,meta_FR$Sample_ID]
#'
#' Melody.model <- melody(rel.abd = CRC_abd_FR,
#'                        sample.data = meta_FR,
#'                        sample.id = "Sample_ID",
#'                        study = "Study",
#'                        disease = "Group")
#' }
#'


melody <- function(rel.abd,
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
                   tune.path = c("gsection", "sequence"),
                   tune.size.sequence = NULL,
                   tune.size.range = NULL,
                   tune.type = c("HBIC", "BIC", "KBIC", "EBIC"),
                   ouput.best.one = TRUE,
                   tol = 1e-3,
                   verbose = FALSE
                   ) {

  cov.type <- match.arg(cov.type)
  tune.path <- match.arg(tune.path)
  tune.type <- match.arg(tune.type)

  ### Generate summary statistics
  summary.stat.study <- melody.get.summary(rel.abd,
                                           sample.data,
                                           sample.id,
                                           study,
                                           disease,
                                           covariates,
                                           cluster,
                                           depth.filter,
                                           prev.filter,
                                           ref,
                                           cov.type,
                                           verbose)

  ### Meta-analysis
  Melody.model <- melody.meta.summary(Melody = summary.stat.study,
                                      tune.path,
                                      tune.size.sequence,
                                      tune.size.range,
                                      tune.type,
                                      ouput.best.one,
                                      tol,
                                      verbose)

  return(Melody.model)
}







