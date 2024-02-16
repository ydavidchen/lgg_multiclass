# TCGA-LGG: Unsupervised EDA

rm(list=ls())
library(pheatmap)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

load(paste0(OUT_DIR,"results/tcgalgg_objects.RData"))

## Unsupervised Learning as evidence for potential:
hm_annot_samps <- data.frame(
  row.names = patients$sample,
  Subtype = patients$Subtype
)
hm_cols <- list(Subtype=c(MUTCODEL="black", MUT="gray40", WT="lightgray"))

lSub <- select_most_var(lgg450, 739)

pheatmap(
  t(lSub),
  cutree_rows = 3,
  annotation_row = hm_annot_samps,
  annotation_colors = hm_cols,
  clustering_method = CL_PARAMS[1],
  clustering_distance_rows = CL_PARAMS[2],
  clustering_distance_cols = CL_PARAMS[2],
  color = HEAT_COLS_BW,
  fontsize = 12,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE
)

## Extract cluster membership for testing:
cl_lgg <- custom_hier_clust(t(lSub), CL_PARAMS[1], CL_PARAMS[2], 3)
patients <- merge(patients, cl_lgg, by.x="sample", by.y="row.names")

table(patients$Cluster, patients$Subtype)

ctabIdh <- table(MUT=patients$IDH, Cluster1ab=patients$Cluster!=2)
ctabIdh <- ctabIdh[c(2,1), c(2,1)]
ctabIdh
fisher.test(ctabIdh, conf.level=0.99)
fisher.test(ctabIdh + 1, conf.level=0.99)

ctabCod <- table(CODEL=patients$codel[patients$IDH], Cluster2=patients$Cluster[patients$IDH]==3)
ctabCod <- ctabCod[c(2,1), c(2,1)]
ctabCod
fisher.test(ctabCod, conf.level=0.99)
