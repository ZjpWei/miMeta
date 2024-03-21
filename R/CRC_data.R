#' Datasets from two metagenomics studies of colorectal cancer (CRC)
#'
#' The "CRC_data" includes "CRC_abd" and "CRC_meta". "CRC_abd" is a sample-by-feature matrix of
#' relative abundance counts from the two studies including 267 species under order "Clostridiales".
#' The "CRC_meta" is a data frame including the sample-level variables from the two studies.
#'
#' @usage data(CRC_data)
#'
#' @references Wirbel, Jakob et al. Nat Med. 2019 Apr;25(4):679-689.
#'
#' @examples
#' \donttest{
#' library("miMeta")
#' data("CRC_data")
#' CRC_abd <- CRC_data$CRC_abd
#' CRC_meta <- CRC_data$CRC_meta
#' }
#'
#' @source <https://github.com/zellerlab/crc_meta>
"CRC_data"
