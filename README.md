## miMeta

`miMeta` is a R package and it implements the new methods for meta-analysis of microbiome association studies that respect the unique features of microbiome data such as compositionality.

See following items for more details:

* [`miMeta` Manual](https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_1.1.0.pdf).

* [`miMeta` Vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

## Reference
* Wei, Z., Chen, G. & Tang, ZZ. *Melody: meta-analysis of microbiome association studies for discovering generalizable microbial signatures*. Genome Biol 26, 245 (2025). https://doi.org/10.1186/s13059-025-03721-4.
  
* Wei, Z., Chen, G. & Tang, ZZ. *Structure-learning meta-analysis resolves cross-context heterogeneity in microbiome associations*.
  
## Author

Zhoujingpeng Wei @[Tang Lab](https://tangzheng1.github.io/tanglab/)

Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison

## System requirements

The `miMeta` package (version 2.1.0) should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the package, users should have R version 4.3.0 or higher.

The package depends on the following R packages: `Mass`, `brglm2`, `abess`, `UpSetR`, `dplyr`, `ggplot2`, `doParallel`, `parallel` and `foreach`.

## Installation guide

Install package from github.
```{r}
devtools::install_github("ZjpWei/miMeta")
```

## Quick start guide

The following are minimal examples on functionalities of `miMeta`. For more detailes, please refer to its [vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

### Melody

* Load package and data
```{r}
library("miMeta")

data("CRC_data")
CRC_abd <- CRC_data$CRC_abd
CRC_meta <- CRC_data$CRC_meta

rel.abd <- list()
covariate.interest <- list()
for(d in c("FR-CRC", "DE-CRC")){
  rel.abd[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
  disease <- as.numeric(CRC_meta$Group[CRC_meta$Study == d] == "CRC")
  names(disease) <- CRC_meta$Sample_ID[CRC_meta$Study == d]
  covariate.interest[[d]] <- data.frame(disease = disease)
}
```

* Perform meta-analysis
```{r}
meta.result <- melody(rel.abd = rel.abd, covariate.interest = covariate.interest, 
                      ref = "Coprococcus catus [ref_mOTU_v2_4874]", verbose = TRUE)
```
### SMESH
* Load package and data
```{r}
library("miMeta")

########## Download the pan-disease data ##########
## About 3 MB from the SMESH reproducibility repository, cached for the
## session, so several examples in one session download it only once.
repo <- "https://raw.githubusercontent.com/ZjpWei/SMESH/main/GMrepo_analysis/Data"
dest <- file.path(tempdir(), "GMrepo")
dir.create(dest, showWarnings = FALSE)
for(f in c("rel.abd.rds", "covariate.interest.rds", "species_lst.Rdata", "GMrepo_project_summary.xlsx")){
  target <- file.path(dest, f)
  if(!file.exists(target))
    download.file(file.path(repo, f), target, mode = "wb", quiet = TRUE)
}
rel.abd <- readRDS(file.path(dest, "rel.abd.rds"))
cov.int <- readRDS(file.path(dest, "covariate.interest.rds"))
load(file.path(dest, "species_lst.Rdata")) # dense_feature: 244 species
psum <- openxlsx::read.xlsx(file.path(dest, "GMrepo_project_summary.xlsx"))

########## Prepare input data ##########
## Restrict to the 244 species used in the pan-disease analysis
feature.counts <- lapply(rel.abd, function(m)
m[, intersect(colnames(m), dense_feature), drop = FALSE])
covariate.interest <- lapply(names(feature.counts), function(d){
  keep <- rownames(feature.counts[[d]])
  matrix(as.integer(cov.int[[d]][keep, 1]), ncol = 1,
  dimnames = list(keep, "disease"))
})
names(covariate.interest) <- names(feature.counts)
```

* Perform meta-analysis
```{r}
mesh <- unique(psum[, c("disease_mesh", "disease_name")])
context <- setNames(
  mesh$disease_name[match(sub("^.*_", "", names(feature.counts)), mesh$disease_mesh)],
  names(feature.counts))

########## Structure-learning meta-analysis ##########
## The number of clusters G is selected from the data (G = NULL, the default).
## "tune" (the default) is recommended for a real analysis; "median" is faster.
fit <- smesh(feature.counts = feature.counts, covariate.interest = covariate.interest,context = context, correct = "median")
fit$disease$context_assignments
head(fit$disease$effects)
fit$disease$diagnostics$gic_path
```

## Issues tracker

Please use the [issues tracker](https://github.com/ZjpWei/miMeta/issues) to report any bugs or give any suggestions.

## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3
