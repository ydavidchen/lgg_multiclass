# German Glioma Network: PCA Execution

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

load(paste0(OUT_DIR,"results/gcn_objects.RData"))

## Execute PCA:
pca <- prcomp(t(lgg450), center=TRUE, scale.=TRUE)
summary(pca)$importance[c(2,3), 1:5]

## Extract PC1-2:
resPca <- as.data.frame(scale(pca$x[ , c(1,2)]))
resPca$sample <- rownames(resPca)
colnames(resPca) <- gsub("PC", "Dimension", colnames(resPca), ignore.case=FALSE)
resPca <- merge(resPca, patients[,c("sample","Subtype")], by="sample")
resPca$Dataset <- "German"

write.csv(resPca, paste0(OUT_DIR,"results/pca_gcn.csv"), row.names=FALSE, quote=FALSE)
