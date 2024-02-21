# TCGA-LGG: PCA Execution

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

load(paste0(OUT_DIR,"results/tcgalgg_objects.RData"))

## Execute PCA & normalize PC1-2:
pca <- prcomp(t(lgg450), center=TRUE, scale.=TRUE)
summary(pca)$importance[c(2,3), 1:5]

## Extract PC1-2:
resPca <- as.data.frame(scale(pca$x[ , c(1,2)]))
resPca$sample <- rownames(resPca)
colnames(resPca) <- gsub("PC", "Dimension", colnames(resPca), ignore.case=FALSE)
resPca <- merge(resPca, patients[,c("sample","Subtype")], by="sample")
resPca$Dataset <- "American"

write.csv(resPca, paste0(OUT_DIR,"results/pca_tcga.csv"), row.names=FALSE, quote=FALSE)
