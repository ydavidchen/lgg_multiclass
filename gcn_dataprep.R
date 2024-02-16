# German Glioma Network Data Preparation

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

## Clinical metadata:
patients <- readxl::read_excel(paste0(DEU_DIR,"SeriesMatrices.xlsx"), sheet=2)
patients <- as.data.frame(patients)
patients$SampleTitle <- as.integer(patients$SampleTitle)

patients$tumor_grade <- NA
patients$tumor_grade[grepl("2", patients$diagnosis)] <- "G2"
patients$tumor_grade[grepl("3", patients$diagnosis)] <- "G3"

patients$IDH <- as.logical(patients$`idh mutation`)
patients$codel <- patients$chr1_loss==1 & patients$chr19_loss==1

patients$Subtype <- ifelse(patients$IDH, "MUT", "WT")
patients$Subtype[patients$IDH & patients$codel] <- "MUTCODEL"

patients <- patients[ , c("Accession","tumor_grade","IDH","codel","Subtype")]
patients$dummy <- ifelse(patients$Subtype=="WT", 2, ifelse(patients$Subtype=="MUT", 0, 1))

## CpG & CGI selected:
cpg_list <- read.csv(paste0(OUT_DIR,"results/cpg_list.csv"))
cgi_stats <- read.table(paste0(OUT_DIR,"results/cgi_sele.txt"))

## 450K data:
lgg450 <- fread(paste0(DEU_DIR,"GSE129477.txt"), data.table=FALSE, header=TRUE)
rownames(lgg450) <- lgg450$V1
lgg450$V1 <- NULL

stopifnot(identical(colnames(lgg450), patients$Accession)) #checkpoint; if not: match

## Aggregation:
lgg450 <- subset(lgg450, rownames(lgg450) %in% cpg_list$Name) #optional: speed things up
lgg450 <- merge(lgg450, cpg_list, by.x="row.names", by.y="Name")
lgg450$Row.names <- NULL
all(lgg450$UCSC_CpG_Islands_Name %in% cgi_stats$UCSC_CpG_Islands_Name)

lgg450 <- aggregate(. ~ UCSC_CpG_Islands_Name, data=lgg450, FUN=mean)
rownames(lgg450) <- lgg450$UCSC_CpG_Islands_Name
lgg450$UCSC_CpG_Islands_Name <- NULL
lgg450 <- data.matrix(lgg450)
dim(lgg450)

## Export:
write.csv(patients, paste0(OUT_DIR,"results/gcn_samples.csv"), row.names=FALSE, quote=FALSE)
write.csv(lgg450, paste0(OUT_DIR,"results/gcn_cgi.csv"), row.names=TRUE, quote=FALSE)
save(
 list = c("patients","lgg450","cpg_list","cgi_stats"),
 file = paste0(OUT_DIR, "results/gcn_objects.RData"),
 compress = TRUE
)
