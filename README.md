## miMeta

`miMeta` is a R package and it implements the new methods for meta-analysis of microbiome association studies that respect the unique features of microbiome data such as compositionality.

See following items for more details:

* [`miMeta` Manual](https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_0.1.0.pdf).

* [`miMeta` vignette](https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

* Article: Wei Z, Chen G, Tang ZZ. Melody identifies generalizable microbial signatures in microbiome association meta-analysis.

## Author

Zhoujingpeng Wei @[Tang](https://tangzheng1.github.io/tanglab/)

Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison

## Installation

Install package from github.
```{r}
devtools::install_github("ZjpWei/miMeta")
```

Install package by [source code](https://github.com/ZjpWei/miMeta/blob/main/miMeta_0.1.0.tar.gz)
```{r}
install.packages("./miMeta_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick start guide

The following are minimal examples on functionalities of `miMeta`. For more detailes, please refer to its [vignette](https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

* Load package and data
```{r}
library("miMeta")

data("CRC_abd", "meta", package = "miMeta")
```

* Perform meta-analysis
```{r}
meta.result <- melody(rel.abd = CRC_abd,
                      sample.data = meta,
                      sample.id = "Sample_ID",
                      study = "Study",
                      disease = "Group")
```

## Issues tracker

Please use the [issues tracker](https://github.com/ZjpWei/miMeta/issues) to report any bugs or give any suggestions.
