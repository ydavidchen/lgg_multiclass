# German Glioma Network: HCLUST

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)

load(paste0(OUT_DIR,"results/gcn_objects.RData"))

## Visualize as heat map:
hm_annot_samps <- data.frame(
  row.names = patients$Accession,
  Subtype = patients$Subtype
)
hm_cols <- list(Subtype=c(MUTCODEL="black", MUT="gray40", WT="lightgray"))

lSub <- select_most_var(lgg450, 739)

pheatmap(
  t(lSub),
  cutree_row = 3,
  annotation_row = hm_annot_samps,
  annotation_colors = hm_cols,
  clustering_method = CL_PARAMS[1],
  clustering_distance_rows = CL_PARAMS[2],
  clustering_distance_cols = CL_PARAMS[2],
  color = HEAT_COLS_BW,
  fontsize = 9,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE
)

## Extract cluster membership for testing:
cl_lgg <- custom_hier_clust(t(lSub), CL_PARAMS[1], CL_PARAMS[2], 3)
patients <- merge(patients, cl_lgg, by.x="Accession", by.y="row.names")

table(patients$Cluster, patients$Subtype)

ctabIdh <- table(MUT=patients$IDH, Cluster23=patients$Cluster!=1)
ctabIdh <- ctabIdh[c(2,1), c(2,1)]
ctabIdh
fisher.test(ctabIdh, conf.level=0.99)

ctabCod <- table(CODEL=patients$codel[patients$IDH], Cluster2=patients$Cluster[patients$IDH]==3)
ctabCod <- ctabCod[c(2,1), c(2,1)]
ctabCod
fisher.test(ctabCod, conf.level=0.99)
