## ----getPackage, echo=TRUE----------------------------------------------------
if(!require("miMeta", quietly = TRUE)){
  devtools::install_github("ZjpWei/miMeta")
}

## ----load, echo=TRUE, message=FALSE, warning=FALSE----------------------------
library("miMeta")
library("tidyverse")

## ----echo=TRUE----------------------------------------------------------------
data("CRC_data")
CRC_abd <- CRC_data$CRC_abd
CRC_meta <- CRC_data$CRC_meta

## ----echo=TRUE, message=TRUE, warning=FALSE-----------------------------------
# Prepare input data
rel.abd <- list()
covariate.interest <- list()
for(d in c("FR-CRC", "DE-CRC")){
  rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
  disease <- as.numeric(CRC_meta$Group[CRC_meta$Study == d] == "CRC")
  names(disease) <- CRC_meta$Sample_ID[CRC_meta$Study == d]
  covariate.interest[[d]] <- data.frame(disease = disease)
}

meta.result <- melody(rel.abd = rel.abd, covariate.interest = covariate.interest, 
                      ref = "Coprococcus catus [ref_mOTU_v2_4874]", verbose = TRUE)

## -----------------------------------------------------------------------------
head(data.frame(coef = meta.result$disease$coef), n = 20)

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
# Generate summary statistics for each study
null.obj.FR <- melody.null.model(rel.abd = rel.abd["FR-CRC"], ref = "Coprococcus catus [ref_mOTU_v2_4874]")
summary.stats.FR <- melody.get.summary(null.obj = null.obj.FR,
                                       covariate.interest = covariate.interest["FR-CRC"])

null.obj.DE <- melody.null.model(rel.abd = rel.abd["DE-CRC"], ref = "Coprococcus catus [ref_mOTU_v2_4874]")
summary.stats.DE <- melody.get.summary(null.obj = null.obj.DE,
                                       covariate.interest = covariate.interest["DE-CRC"])

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
# Concatenate summary statistics
summary.stats.all <- c(summary.stats.FR, summary.stats.DE)

## ----echo=TRUE----------------------------------------------------------------
# Meta-analysis to harmonize and combine summary statistics across studies
meta.result.2 <- melody.meta.summary(summary.stats = summary.stats.all, verbose = TRUE)

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
# Load metabolite data: 
# You may see the following website on how to directly load data from 
# github into R https://github.com/ZjpWei/Melody/raw/main/Metabolite.rda
load(file = url("https://github.com/ZjpWei/Melody/raw/main/Metabolite.rda"))

# Change genura names
for(d in names(otu_data_lst)){
  colnames(otu_data_lst[[d]]) <-  gsub(".*;g__", "", colnames(otu_data_lst[[d]]))
}

# Get null model
null.obj <- melody.null.model(rel.abd = otu_data_lst, covariate.adjust = covariates_adjust_lst)

# Get summary statistics
summary.stats <- melody.get.summary(null.obj = null.obj, covariate.interest = cmpd_data_lst, cluster = cluster_data_lst)

# Meta-analysis
meta.scan.result <- melody.meta.summary(summary.stats = summary.stats, verbose = TRUE)

## ----echo=TRUE----------------------------------------------------------------
selected.num <- sort(unlist(lapply(meta.scan.result, function(d){sum(d$coef!=0)})), decreasing = TRUE)
top.cov.name <- names(selected.num)[1:min(5, length(selected.num))]
coef_mat <- sapply(meta.scan.result[top.cov.name], function(d){d$coef}, simplify = TRUE)
rownames(coef_mat) <- gsub(".*;g__","g__",rownames(coef_mat))
head(coef_mat, n = 20)

## -----------------------------------------------------------------------------
sessionInfo()

