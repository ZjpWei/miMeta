## miMeta

`miMeta` is a R package and it implements the new methods for meta-analysis of microbiome association studies that respect the unique features of microbiome data such as compositionality.

See following items for more details:

* [`miMeta` Manual](https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_0.1.0.pdf).

* [`miMeta` vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

* Article: Wei Z, Chen G, Tang ZZ. Melody: Meta-analysis of Microbiome Association Studies for Discovering Generalizable Microbial Signatures.

## Author

Zhoujingpeng Wei @[Tang Lab](https://tangzheng1.github.io/tanglab/)

Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison

## Installation guide

Install package from github.
```{r}
devtools::install_github("ZjpWei/miMeta")
```

Install package by [source code](https://github.com/ZjpWei/miMeta/blob/main/miMeta_0.1.0.tar.gz)
```{r}
install.packages("./miMeta_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick start guide

The following are minimal examples on functionalities of `miMeta`. For more detailes, please refer to its [vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

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
                      ref = "Coprococcus catus [ref_mOTU_v2_4874]")
```

## Issues tracker

Please use the [issues tracker](https://github.com/ZjpWei/miMeta/issues) to report any bugs or give any suggestions.

## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3
