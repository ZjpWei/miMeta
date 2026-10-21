## miMeta

`miMeta` is a R package and it implements new methods for meta-analysis of microbiome association studies that respect the unique features of microbiome data such as compositionality.

The package provides two methods:

* **Melody** combines summary statistics across studies to select microbial signatures that generalize.
* **SMESH** learns latent clusters of contexts and separates signatures whose association is selected as nonzero in every cluster from those selected in only a subset of clusters.

See following items for more details:

* [`miMeta` Manual](https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_2.1.0.pdf).

* [`miMeta` Vignette for Melody](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

* [`miMeta` Vignette for SMESH](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/smesh_vignette.html).

* Article (Melody): Wei Z, Chen G, Tang ZZ. Melody: meta-analysis of microbiome association studies for discovering generalizable microbial signatures. Genome Biology. 2025;26:245. doi:10.1186/s13059-025-03721-4.

* Article (SMESH): Wei Z, Chen G, Tang ZZ. Structure-learning meta-analysis resolves cross-context heterogeneity in microbiome associations.

## Author

Zhoujingpeng Wei @[Tang Lab](https://tangzheng1.github.io/tanglab/)

Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison

## System requirements

The `miMeta` package (version 2.1.0) should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the package, users should have R version 4.3.0 or higher.

The package depends on the following R packages: `MASS`, `Matrix`, `brglm2`, `abess`, `UpSetR`, `dplyr`, `ggplot2`, `doParallel`, `parallel`, `foreach` and `PALM`.

`PALM` is installed from [github.com/ZjpWei/PALM_package](https://github.com/ZjpWei/PALM_package); the install command below pulls it in automatically.

## Installation guide

Install package from github.
```{r}
devtools::install_github("ZjpWei/miMeta")
```

## Quick start guide

The following are minimal examples on functionalities of `miMeta`. For more details, please refer to the [Melody vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html) and the [SMESH vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/smesh_vignette.html).

* Load package and data
```{r}
library("miMeta")

data("CRC_data")
CRC_abd <- CRC_data$CRC_abd
CRC_meta <- CRC_data$CRC_meta

feature.counts <- list()
covariate.interest <- list()
for(d in c("FR-CRC", "DE-CRC")){
  feature.counts[[d]] <- CRC_abd[CRC_meta$Sample_ID[CRC_meta$Study == d],]
  disease <- as.numeric(CRC_meta$Group[CRC_meta$Study == d] == "CRC")
  names(disease) <- CRC_meta$Sample_ID[CRC_meta$Study == d]
  covariate.interest[[d]] <- data.frame(disease = disease)
}
```

* Perform meta-analysis with Melody
```{r}
meta.result <- melody(feature.counts = feature.counts, covariate.interest = covariate.interest, 
                      ref = "Coprococcus catus [ref_mOTU_v2_4874]", verbose = TRUE)
```

* Perform structure-learning meta-analysis with SMESH

SMESH learns clusters of contexts, so it needs more contexts than the two studies above. This
example uses `sim_data`, 10 simulated contexts in two latent clusters with 100 features, of which
30 are true signatures.
```{r}
data("sim_data")
feature.counts <- sim_data$sim_abd
sim_meta <- sim_data$sim_meta

covariate.interest <- lapply(names(feature.counts), function(d){
  m <- sim_meta[sim_meta$Context == d, ]
  matrix(m$Group, ncol = 1, dimnames = list(m$Sample_ID, "disease"))
})
names(covariate.interest) <- names(feature.counts)

smesh.result <- smesh(feature.counts = feature.counts, covariate.interest = covariate.interest,
                      verbose = TRUE)
smesh.result$disease$context_assignments
```

By default the number of clusters is selected by the minimum generalized information criterion
(GIC); pass `G` to fix it instead. `correct = "tune"` (the default) is recommended, and
`correct = "median"` runs faster. To fit the model on context-level summary statistics you
already have, use `smesh.summary()`. The
[SMESH vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/smesh_vignette.html)
runs both routes on `sim_data`, checks the fit against the simulation truth, and shows how to build
context-level summary statistics.

## Issues tracker

Please use the [issues tracker](https://github.com/ZjpWei/miMeta/issues) to report any bugs or give any suggestions.

## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3
