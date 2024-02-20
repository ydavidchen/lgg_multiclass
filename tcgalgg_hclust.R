# TCGA-LGG: Unsupervised EDA

rm(list=ls())
library(pheatmap)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

load(paste0(OUT_DIR,"results/tcgalgg_objects.RData"))
patients$Class <- ifelse(patients$Subtype=="WT", "0", ifelse(patients$Subtype=="MUT", "1A", "1B"))

## Cohort summary:
t1 <- tableone::CreateTableOne(c("IDH","codel","Subtype"), data=patients)
print(t1, showAllLevels=TRUE)

## Unsupervised ML & association test: 
hm_annot_samps <- data.frame(
  row.names = patients$sample,
  Class = patients$Class,
  Codel = ifelse(patients$codel, "Yes", "No"),
  IDH = ifelse(patients$IDH, "Mutated", "Normal")
)

lSub <- select_most_var(lgg450, 739)

pheatmap(
  t(lSub),
  cutree_rows = 3,
  annotation_row = hm_annot_samps,
  annotation_colors = HM_COLS,
  clustering_method = CL_PARAMS[1],
  clustering_distance_rows = CL_PARAMS[2],
  clustering_distance_cols = CL_PARAMS[2],
  color = HEAT_COLS,
  fontsize = 12,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE
)

cl_lgg <- custom_hier_clust(t(lSub), CL_PARAMS[1], CL_PARAMS[2], 3)
cl_lgg <- merge(cl_lgg, patients[ , c("sample","Class")], by.x="row.names", by.y="sample")

ctab <- table(cl_lgg$Cluster, cl_lgg$Class)
ctab

ft <- fisher.test(ctab, workspace=1e6)
ft$p.value
