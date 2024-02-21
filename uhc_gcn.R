# German Glioma Network (GCN): Unsupervised Exploration

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)

NUM_CGI <- 750
NUM_CLS <- 3

load(paste0(OUT_DIR,"results/gcn_objects.RData"))
patients$Class <- ifelse(patients$Subtype=="WT", "0", ifelse(patients$Subtype=="MUT", "1A", "1B"))

## Cohort summaries:
t1 <- tableone::CreateTableOne(c("IDH","codel","Subtype"), data=patients)
print(t1, showAllLevels=TRUE)

## UHC heatmap:
hm_annot_samps <- data.frame(
  row.names = patients$sample,
  Class = patients$Class
)

lSub <- winsorize(select_most_var(lgg450, NUM_CGI), -6, 6)

pheatmap(
  lSub,
  cutree_cols = NUM_CLS,
  annotation_col = hm_annot_samps,
  annotation_colors = HM_COLS,
  clustering_method = CL_PARAMS[1],
  clustering_distance_rows = CL_PARAMS[2],
  clustering_distance_cols = CL_PARAMS[2],
  color = HEAT_COLS,
  legend_breaks = HM_SEQ,
  fontsize = 15,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE
)

## Extract cluster assignment for Fisher's test w/ class:
cl_lgg <- custom_hier_clust(t(lSub), CL_PARAMS[1], CL_PARAMS[2], NUM_CLS)
cl_lgg <- merge(cl_lgg, patients[ , c("sample","Class")], by.x="row.names", by.y="sample")

ctab <- table(Cluster=cl_lgg$Cluster, Class=cl_lgg$Class)
ctab
round(100 * prop.table(ctab, margin=2), 1)

ft <- fisher.test(ctab)
ft$p.value
