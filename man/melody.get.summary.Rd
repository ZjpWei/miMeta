% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/melody.get.summary.R
\name{melody.get.summary}
\alias{melody.get.summary}
\title{Generate Melody summary statistics for one or multiple studies.}
\usage{
melody.get.summary(
  rel.abd,
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
  verbose = FALSE
)
}
\arguments{
\item{...}{Same arguments as the function melody}
}
\value{
A Melody class.
\item{summary.stat.study}{A list includes summary statistics for each study.}
\item{dat.inf}{A list includes study number; reference feature for each study; feature number and feature names.}
\item{taxa.set}{A list includes the which features are included or not in each study.}
}
\description{
The melody.get.summary, melody.merge.summary, and melody.meta.summary are three functions that represent individual components in the melody pipeline.
The function melody.get.summary takes individual-level data from one or multiple studies and outputs summary statistics for the studies in a Melody object.
The output can be directly used by the functions melody.merge.summary and melody.meta.summary.
}
\examples{
\donttest{
library("miMeta")
data("CRC_abd")
data("meta")

########## Generate summary statistics for study FR ##########
meta_FR <- meta[meta$Study == "FR-CRC",]
CRC_abd_FR <- CRC_abd[,meta_FR$Sample_ID]

sumstats <- melody.get.summary(rel.abd = CRC_abd_FR,
                               sample.data = meta_FR,
                               sample.id = "Sample_ID",
                               study = "Study",
                               disease = "Group",
                               verbose = TRUE)
}
}
\seealso{
\code{\link{melody}},
\code{\link{melody.meta.summary}},
\code{\link{melody.merge.summary}}
}
