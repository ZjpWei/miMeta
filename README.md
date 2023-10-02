## miMeta

`miMeta` is a R package and it implements the new methods for meta-analysis of microbiome association studies that respect the unique features of microbiome data such as compositionality.

See following items for more details:

* [`miMeta` Manual](https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_0.1.0.pdf).

* [`miMeta` vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

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

The following are minimal examples on functionalities of `miMeta`. For more detailes, please refer to its [vignette](https://htmlpreview.github.io/?https://github.com/ZjpWei/miMeta/blob/main/doc/miMeta_vignette.html).

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

## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3
