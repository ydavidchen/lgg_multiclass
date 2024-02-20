# German Glioma Network (GCN) Data Preparation

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

# --------------- Part I: Sample-level Clinical Metadata ---------------
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
patients$Subtype <- factor(patients$Subtype, c("WT","MUT","MUTCODEL"))
patients$dummy <- ifelse(patients$Subtype=="WT", 2, ifelse(patients$Subtype=="MUT", 0, 1))

# --------------- Part II. CGI Aggregation ---------------
cpg_list <- read.csv(paste0(OUT_DIR,"results/cpg_list.csv"))
cgi_stats <- read.table(paste0(OUT_DIR,"results/cgi_sele.txt"))

lgg450 <- fread(paste0(DEU_DIR,"GSE129477.txt"), data.table=FALSE, header=TRUE)
rownames(lgg450) <- lgg450$V1
lgg450$V1 <- NULL

lgg450 <- subset(lgg450, rownames(lgg450) %in% cpg_list$Name) #optional: speed things up
lgg450 <- merge(lgg450, cpg_list, by.x="row.names", by.y="Name")
lgg450$Row.names <- NULL
all(lgg450$UCSC_CpG_Islands_Name %in% cgi_stats$UCSC_CpG_Islands_Name)

lgg450 <- aggregate(. ~ UCSC_CpG_Islands_Name, data=lgg450, FUN=mean)
rownames(lgg450) <- lgg450$UCSC_CpG_Islands_Name
lgg450$UCSC_CpG_Islands_Name <- NULL
lgg450 <- data.matrix(lgg450)
lgg450 <- minfi::logit2(winsorize(lgg450, 0.0001, 0.9999)) #M-value
dim(lgg450)

# --------------- Part III. Proc. Data Export ---------------
save(
 list = c("patients","lgg450","cpg_list","cgi_stats"),
 file = paste0(OUT_DIR, "results/gcn_objects.RData"),
 compress = TRUE
)

lgg450 <- as.data.frame(t(lgg450))
colnames(patients)[colnames(patients)=="Accession"] <- "sample"
lgg450 <- merge(patients, lgg450, by.x="sample", by.y="row.names")
write.csv(lgg450, paste0(OUT_DIR, "results/gcnlgg.csv"), row.names=FALSE, quote=FALSE)
