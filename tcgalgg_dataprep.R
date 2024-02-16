# TCGA-LGG Data Preparation

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

# --------------- Part I: Patient- & Sample-level Covariates ---------------
patients <- read.csv(paste0(TCGA_DIR,"nationwidechildrens.org_LGG_bio.patient.tsv"), sep="\t", na.strings=c("NA",""))
colnames(patients)[1] <- "patient"
patients$age <- as.integer(patients$age_at_initial_pathologic_diagnosis)
patients$raceWhite <- patients$race == "WHITE"
patients$ethHisp <- patients$ethnicity == "HISPANIC OR LATINO"
patients <- patients[ , c("patient","gender","histologic_diagnosis","tumor_grade","age","raceWhite","ethHisp")]

## Samples:
samples <- read.csv(paste0(TCGA_DIR,"nationwidechildrens.org_LGG_bio.sample.tsv"), sep="\t")
samples <- subset(samples, sample_type == "Primary Tumor")
samples <- samples[ , "sample", drop=FALSE]
samples$patient <- substr(samples$sample, 1, 12)

patients <- merge(samples, patients, by="patient")
patients$SAMPLE_ID <- substr(patients$sample, 1, 15)

## Additional covariates from cBioPortal:
cbio <- read.csv(paste0(TCGA_DIR,"cbio_queries/lgg_tcga_pan_can_atlas_2018_clinical_data.tsv"), sep="\t", check.names=FALSE)
cbio$Subtype <- toupper(gsub("LGG_|-|IDH|non-codel", "", cbio$Subtype))
cbio <- cbio[ , c("Patient ID","Sample ID","Subtype","Radiation Therapy","Fraction Genome Altered")]

patients <- merge(patients, cbio, by.x="SAMPLE_ID", by.y="Sample ID", all.x=TRUE)
patients$Subtype <- factor(patients$Subtype, c("WT","MUT","MUTCODEL"))
patients$SAMPLE_ID <- NULL

# --------------- Part II. 450K Data ---------------
cpg_list <- read.csv(paste0(OUT_DIR,"results/cpg_list.csv"))
cgi_stats <- read.table(paste0(OUT_DIR,"results/cgi_sele.txt"))

## 450k Data Loading & Subsetting:
lgg450 <- fread(paste0(TCGA_DIR,"jhu-usc.edu_LGG_HumanMethylation450.betaValue.tsv"), data.table=FALSE)
rownames(lgg450) <- lgg450$V1
lgg450$V1 <- NULL
colnames(lgg450) <- substr(colnames(lgg450), 1, 16)
lgg450 <- lgg450[ , colnames(lgg450) %in% patients$sample]

## Aggregation:
lgg450 <- subset(lgg450, rownames(lgg450) %in% cpg_list$Name) #optional: speed things up
lgg450 <- merge(lgg450, cpg_list, by.x="row.names", by.y="Name")
lgg450$Row.names <- NULL

## Check to see if all CGIs available in TCGA data:
all(lgg450$UCSC_CpG_Islands_Name %in% cgi_stats$UCSC_CpG_Islands_Name)

lgg450 <- aggregate(. ~ UCSC_CpG_Islands_Name, data=lgg450, FUN=mean)

rownames(lgg450) <- lgg450$UCSC_CpG_Islands_Name
lgg450$UCSC_CpG_Islands_Name <- NULL
lgg450 <- data.matrix(lgg450)
dim(lgg450)


## Export: 
write.csv(patients, paste0(OUT_DIR,"results/tcgalgg_patients.csv"), row.names=FALSE, quote=FALSE)
write.csv(lgg450, paste0(OUT_DIR,"results/tcgalgg_cgi_avg.csv"), row.names=TRUE, quote=FALSE)
save(
  list = c("lgg450", "patients", "cpg_list", "cgi_stats"),
  file = paste0(OUT_DIR, "results/tcgalgg_objects.RData"),
  compress = TRUE
)
