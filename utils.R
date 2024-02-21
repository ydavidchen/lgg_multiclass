# Utility Module for LGG ML

library(data.table)
library(ggplot2)

# Paths ** Mask upon GitPush **:
DEU_DIR <- "******** MASKED ********"
TCGA_DIR <- "******** MASKED ********"
OUT_DIR <- "******** MASKED ********"
ILMN_PATH <- "******** MASKED ********/humanmethylation450_15017482_v1-2.csv"

## Constants:
STRNAS <- c("", "NA", "n/a")
CL_PARAMS <- c("ward.D", "euclidean")

## Data loading methods:
load_450k_annot <- function(path=ILMN_PATH, simplify=TRUE) {
  annot <- fread(path, na.strings=STRNAS, skip=7)
  annot$chr <- paste0("chr", annot$CHR)
  annot$pos <- as.integer(annot$MAPINFO)
  annot$strandPlus <- annot$Strand == "F"
  
  annot$Strand <- annot$MAPINFO <- annot$CHR <- NULL
  if(simplify) annot <- annot[ , c("Name","chr","pos","strandPlus",
                                   "Relation_to_UCSC_CpG_Island","UCSC_CpG_Islands_Name",
                                   "UCSC_RefGene_Name","UCSC_RefGene_Group","UCSC_RefGene_Accession",
                                   "Methyl27_Loci", "Probe_SNPs","Probe_SNPs_10")]
  return(annot)
}

select_most_var <- function(mat, size) {
  #'@description Select most variable rows from a matrix
  #'@param size Either integer>0 or a proportion<1
  mVars <- matrixStats::rowVars(mat)
  names(mVars) <- rownames(mat)
  
  mVars <- sort(mVars, decreasing=TRUE)
  if(size < 1) size <- round(size * nrow(mat))
  mVars <- mVars[1:size]
  
  return(mat[rownames(mat) %in% names(mVars), ])
}


## Statistical methods:
custom_scale <- function(vec) {
  #'@usage apply(<dataframe>, <margin>, FUN=custom_scale)
  mu <- mean(vec, na.rm=TRUE)
  sig <- sd(vec, na.rm=TRUE)
  return((vec-mu)/sig)
}

winsorize <- function(mat, lower, upper) {
  #'@description Constrains extreme values in matrix `mat`
  #'@describeIn Chen 2023 JOMES & Chen 2024 JMCCPL
  mat[mat < lower] <- lower
  mat[mat > upper] <- upper
  return(mat)
}

custom_hier_clust <- function(tMat, method_hc, method_dist, num_cl) {
  #'@param tMat Expression matrix where row=Samples(observations), col=genes(features)
  #'@describeIn Chen 2022 MedResArch
  myDist <- dist(tMat, method=method_dist)
  myHcl <- hclust(myDist, method=method_hc)
  
  plot(myHcl, labels=FALSE)
  rect.hclust(myHcl, k=num_cl)
  
  membership <- data.frame(cutree(myHcl, k=num_cl))
  colnames(membership) <- "Cluster"
  return(membership)
}

## Data visualization objects & helpers:
BINARY_COLORS <- c(Yes="black", No="lightgray")
HEAT_COLS <- colorRampPalette(c("blue","lightgray","red"))(1024)
HM_COLS <- list(
  Class = c(`0`="lightgray",`1A`="gray40",`1B`="black"),
  Codel = BINARY_COLORS,
  IDH = c(Mutated="black", Normal="lightgray")
)
HM_SEQ <- seq(-6, 4, 2)
THEME_SCATTER <- theme_gray() + 
  theme(axis.text.x=element_text(size=10,color="black"), axis.title.x=element_text(size=15,color="black"),
        axis.text.y=element_text(size=10,color="black"), axis.title.y=element_text(size=15,color="black"),
        strip.text.x=element_text(size=15,color="black"),
        legend.position="top", legend.text=element_text(size=15,color="black"), legend.title=element_text(size=15))

