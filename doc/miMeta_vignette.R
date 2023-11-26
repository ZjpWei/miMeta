## ----getPackage, echo=TRUE----------------------------------------------------
if(!require("miMeta", quietly = TRUE)){
  devtools::install_github("ZjpWei/miMeta")
}

## ----load, echo=TRUE----------------------------------------------------------
library("miMeta")
library("DT")
library("ggplot2")
library("tidyverse")

## ----echo=TRUE----------------------------------------------------------------
data("CRC_abd", "meta", package = "miMeta")
table(meta[,c("Study", "Group")])
CRC_abd <- CRC_abd[,meta$Sample_ID]

## ----echo=TRUE, message=FALSE, warning=FALSE, fig.cap="a plot for microbial feature overlap among studies: this plot shows the number of features shared among studies. "----
# Pick a reference taxa
ref <- "Coprococcus catus [ref_mOTU_v2_4874]"

# melody can analyze single study or multiple studies. 
meta.result <- melody(rel.abd = CRC_abd,
                      sample.data = meta,
                      sample.id = "Sample_ID",
                      study = "Study",
                      disease = "Group",
                      ref = ref,
                      parallel.core = 1,
                      verbose = TRUE)

## ----echo=FALSE, fig.cap="a plot showing the absolute-abundance coefficient estimates of the selected microbial features in the best model under the best subset size."----
taxa_tab <- data.frame(taxa = names(which(meta.result$coef!=0)),
                       coef = as.numeric(meta.result$coef[meta.result$coef!=0]))

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


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
# Geneeate summary statistisc for each study
meta_FR <- meta[meta$Study == "FR-CRC",]
summary.FR <- melody.get.summary(rel.abd = CRC_abd[,meta_FR$Sample_ID],
                                 sample.data = meta_FR,
                                 sample.id = "Sample_ID",
                                 study = "Study",
                                 disease = "Group",
                                 ref = ref)

meta_DE <- meta[meta$Study == "DE-CRC",]
summary.DE <- melody.get.summary(rel.abd = CRC_abd[,meta_DE$Sample_ID],
                                 sample.data = meta_DE,
                                 sample.id = "Sample_ID",
                                 study = "Study",
                                 disease = "Group",
                                 ref = ref)

meta_CN <- meta[meta$Study == "CN-CRC",]
summary.CN <- melody.get.summary(rel.abd = CRC_abd[,meta_CN$Sample_ID],
                                 sample.data = meta_CN,
                                 sample.id = "Sample_ID",
                                 study = "Study",
                                 disease = "Group",
                                 ref = ref)

meta_US <- meta[meta$Study == "US-CRC",]
summary.US <- melody.get.summary(rel.abd = CRC_abd[,meta_US$Sample_ID],
                                 sample.data = meta_US,
                                 sample.id = "Sample_ID",
                                 study = "Study",
                                 disease = "Group",
                                 ref = ref)

meta_AT <- meta[meta$Study == "AT-CRC",]
summary.AT <- melody.get.summary(rel.abd = CRC_abd[,meta_AT$Sample_ID],
                                 sample.data = meta_AT,
                                 sample.id = "Sample_ID",
                                 study = "Study",
                                 disease = "Group",
                                 ref = ref)

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
# Merge summary statistics
summary.merge <- melody.merge.summary(list(summary.FR, summary.DE, summary.CN, summary.US, summary.AT))

## ----echo=TRUE----------------------------------------------------------------
# Meta-analysis for merged summary statistics
meta.result.merge <- melody.meta.summary(Melody = summary.merge, verbose = TRUE)

## ----echo=FALSE, fig.cap="a plot showing the absolute-abundance coefficient estimates of the selected microbial features in the best model under the best subset size by using merged summary statistics."----
taxa_tab <- data.frame(taxa = names(which(meta.result$coef!=0)),
                       coef = as.numeric(meta.result$coef[meta.result$coef!=0]))

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


## -----------------------------------------------------------------------------
sessionInfo()

