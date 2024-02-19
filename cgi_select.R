# Autosomal CGI & CpG Selection from Illumina Annotation File

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(dplyr)

MIN_SIZE <- 800
MIN_CPGS <- 5
# CGI_EXCL <- "chr8:144635308-144636026"

## Data Filtering:
cpg_list <- load_450k_annot()
cpg_list <- subset(cpg_list, ! chr %in% c("chrX","chrY"))
cpg_list <- subset(cpg_list, ! is.na(UCSC_CpG_Islands_Name))
cpg_list <- subset(cpg_list, Relation_to_UCSC_CpG_Island == "Island")
cpg_list <- subset(cpg_list, ! is.na(UCSC_RefGene_Name))
cpg_list <- subset(cpg_list, grepl("TSS",UCSC_RefGene_Group))

## CpG statistics:
cgi_stats <- cpg_list %>% 
  select(Name, UCSC_CpG_Islands_Name) %>%
  group_by(UCSC_CpG_Islands_Name) %>%
  count() %>%
  arrange(desc(n))

cgi_stats$chr <- gsub("\\:.*", "", cgi_stats$UCSC_CpG_Islands_Name)
spl <- strsplit(gsub(".*:", "", cgi_stats$UCSC_CpG_Islands_Name), "-")
cgi_stats$hg19_start <- as.integer(sapply(spl, "[", 1))
cgi_stats$hg19_end <- as.integer(sapply(spl, "[", 2))
cgi_stats$size_bp <- abs(cgi_stats$hg19_end - cgi_stats$hg19_start)

barplot(cgi_stats$n, main="CpGs per Island")
summary(cgi_stats$n) #avg1015, med825
summary(cgi_stats$size_bp) #avg5.59, median5
ggplot(cgi_stats, aes(n, size_bp)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Number of CpGs", y="Island Size (bp)") +
  THEME_SCATTER

quantile(cgi_stats$n, c(0.75, 0.9, 0.99))

## Select Islands:
cgi_stats <- subset(cgi_stats, n >= MIN_CPGS & size_bp >= MIN_SIZE)

## Subset CGIs: 
## Downstream, retrospectively revisit the entire list of loci
cpg_list <- subset(cpg_list, UCSC_CpG_Islands_Name %in% cgi_stats$UCSC_CpG_Islands_Name)
cpg_list <- cpg_list[ , c("Name","UCSC_CpG_Islands_Name")] #columns for merging

## Export: 
write.csv(cpg_list, paste0(OUT_DIR,"results/cpg_list.csv"), row.names=FALSE, quote=FALSE)
write.table(cgi_stats, paste0(OUT_DIR,"results/cgi_sele.txt"))
