# Data Visualization: Paper Figure 2

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

tcga <- read.csv(paste0(OUT_DIR,"results/pca_tcga.csv"))
gcn <- read.csv(paste0(OUT_DIR,"results/pca_gcn.csv"))

df <- rbind(tcga, gcn)
df$Class <- ifelse(df$Subtype=="WT", "0", ifelse(df$Subtype=="MUT", "1A", "1B"))

table(df$Class)

ggplot(df, aes(Dimension1, Dimension2, color=Class)) +
  geom_point(size=5, alpha=0.5) +
  scale_color_brewer(palette="Set1", direction=-1) +
  scale_x_continuous(limits=c(-3,2.5)) +
  facet_wrap(~ Dataset, ncol=2) +
  THEME_SCATTER
