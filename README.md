# miMeta
R pcakage


## Introduction

`miMeta` is a R package and it implements the new methods for meta-analysis of microbiome association studies that respect the unique features of microbiome data such as compositionality. It first generate summary statistics (microbiome-disease association coefficient estimates and their variances) for individual studies. And then it combines the summary statistics across studies to select disease-associated microbial signatures based on the average absolute-abundance association coefficients inferred from the summary statistics. In particular, the selection of signature is operated through a best-subset selection.

## Installation

Download package from github.

```{r getPackage, eval=FALSE}
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("miMeta")
```

Load the package. 

```{r load, include=FALSE}
library(miMeta)
library(tidyverse)
library(DT)
library(ggplot2)
```


## Example Data

`miMeta` requires a properly formatted collection of microbial community studies, 
with both feature abundances and accompanying metadata. Here we use the five published 
colorectal cancer (CRC) stool metagenomic studies and the studies are already packaged.

Importantly, `miMeta` asks that feature abundances be provided as a 
feature-by-sample matrix, and the metadata be provided as a data frame. 
The two objects should agree on sample IDs, that is, `colname` of the feature 
abundance matrix and `sample.id` of the metadata data frame must agree. 

To minimize users' efforts in loading data to run the examples, we only include 100 most abundant species in the example data. The feature abundances and metadata can be loaded with the following code chunk. For the interested user, the full data with 849 species can be found in the github data file (https://github.com/ZjpWei/miMeta).

We consider the following covriates in meta data:

* Sample identity: "Sample_ID"

* Study name: "Study"

* Disease status: "Group"

```{r}
data("CRC_abd", "meta", package = "miMeta")
table(meta[,c("Study", "Group")])
```

## meta-analysis by one function `melody`

One of the most common meta-analysis goals is to combine association effects across batches/studies to identify consistent overall effects. `melody` provides a straightforward interface to this task by first generating summary statistics in individual batches/studies, and then select the best-subset model by meta-analysis across the summary statistics. 

```{r message=FALSE, warning=FALSE}
# melody can analyze single study or multiple studies. Following is an example on 
# single study.
Melody.model <- melody(rel.abd = CRC_abd,
                       sample.data = meta,
                       sample.id = "Sample_ID",
                       study = "Study",
                       disease = "Group",
                       tune.type = "HBIC")
```

The `melody` result can determine the microbial features by the best model that's selected by information criterion according to the covriate of interest. It returns a list of more than one components. Depending on your input option for `ouput.best.one`, coef provides a estimates vector for best model or estimates matrix for all searched models. See help(melody) for more details. We can visualize the selected (non-zero in coef) species associated with CRC in these studies/samples as following.

```{r}
taxa_tab <- data.frame(taxa = names(which(Melody.model$coef!=0)),
                       coef = as.numeric(Melody.model$coef[Melody.model$coef!=0]))

taxa_tab %>% arrange(coef) %>% 
  mutate(taxa = factor(taxa, levels = taxa)) %>% 
  ggplot() + geom_point(aes(x= factor(taxa), y= coef)) + 
  theme_classic() + coord_flip() + ylab("coef") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "right") +
  scale_x_discrete(position='bottom') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide="none") +
  geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed")

```

## Calculate summary statsitics by `melody.get.summary`

miMeta package performs meta-analysis by two steps. the first step is calculating the summary statistics; the second step is analyzing the summary statsitics across the batches/studies. Different from using one function in analysis, `melody.get.summary` returns the summary statistics (`Melody` object) in first step. 

Followings are the result for `Melody` object.

* dat.inf: some information for the data

    + L: study number
    
    + K: taxa number
    
    + taxa.names: taxa names across the studies.
    
    + ref: reference taxa for generating the summary statistics
    
* taxa.set: a list of vectors that indicate the presentation of each taxon in studies.

* summary.stat.study: a list inlcudes summary statistis estimates and covriate matrix for each study.

```{r message=FALSE, warning=FALSE}
meta_FR <- meta[meta$Study == "FR-CRC",]
CRC_abd_FR <- CRC_abd[,meta_FR$Sample_ID]
summary.FR <- melody.get.summary(rel.abd = CRC_abd_FR,
                                 sample.data = meta_FR,
                                 sample.id = "Sample_ID",
                                 study = "Study",
                                 disease = "Group")
```


## Merge summary statsitics by `melody.merge.summary`

`melody.merge.summary` function can merge multiple summary statistics (`Melody` objects) from different batches/studies to single `Melody` object. The result can be directly applied to `melody.meta.summary` function for meta-analysis cross these batches/studies.

```{r message=FALSE, warning=FALSE}
# Generate summary statistics inclding study AT and CN
meta_other <- meta[meta$Study %in% c("AT-CRC", "CN-CRC", "DE-CRC", "US-CRC"),]
CRC_abd_other <- CRC_abd[,meta_other$Sample_ID]
summary.other <- melody.get.summary(rel.abd = CRC_abd_other,
                                    sample.data = meta_other,
                                    sample.id = "Sample_ID",
                                    study = "Study",
                                    disease = "Group")

summary.merge <- melody.merge.summary(list(summary.FR, summary.other))
```

## Meta-analysis by `melody.meta.summary`

`melody.meta.summary` need the `Melody` object from `melody.get.summary` or `melody.merge.summary` for further meta-analysis. Following example shows the meta-analysis cross three different batches/studies by merged summary statsitics.

```{r}
# meta-analysis for merged summary statistics
Melody.model.merge <- melody.meta.summary(Melody = summary.merge, tune.type = "HBIC")
```

We can visualize the selected (non-zero in coef) species associated with CRC in these studies/samples as following.

```{r}
taxa_tab <- data.frame(taxa = names(which(Melody.model.merge$coef!=0)),
                       coef = as.numeric(Melody.model.merge$coef[Melody.model.merge$coef!=0]))

taxa_tab %>% arrange(coef) %>% 
  mutate(taxa = factor(taxa, levels = taxa)) %>% 
  ggplot() + geom_point(aes(x= factor(taxa), y= coef)) + 
  theme_classic() + coord_flip() + ylab("coef") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "right") +
  scale_x_discrete(position='bottom') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey'), guide="none") +
  geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed")

```

## References
Zhoujingpeng Wei, Guanhua Chen and Zheng-Zheng Tang
