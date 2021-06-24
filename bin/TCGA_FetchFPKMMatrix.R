# Libraries
require(tidyverse)
require(readr)
require(data.table)
require(TCGAbiolinks)
library(biomaRt)

# Cancer Type
typeCancer = "COAD"
typeTCGA = paste("TCGA", typeCancer, sep ="-")

# Ancestry
ancestryTCGA = read.delim(file = "./data/anc/TCGA_consensus_ancestry.tsv",
                          sep = "\t",
                          header = T) %>% as.data.table()


# Query
query.exp.TCGA.anc.fpkm <- GDCquery(project = typeTCGA,
                                    legacy = F,
                                    data.category = "Transcriptome Profiling",
                                    data.type = "Gene Expression Quantification",
                                    workflow.type = "HTSeq - FPKM",                  
                                    sample.type = c("Solid Tissue Normal", "Primary Tumor"))

GDCdownload(query.exp.TCGA.anc.fpkm)
exp.TCGA.anc.fpkm <- GDCprepare(query.exp.TCGA.anc.fpkm)


# ancestry all
AFR = ancestryTCGA[tumor_type == typeCancer & consensus_ancestry %in% c("afr_admix","afr")]$patient
EUR = ancestryTCGA[tumor_type == typeCancer & consensus_ancestry %in% c("eur_admix","eur")]$patient

# why is the sample count the same? they aren't the same samples
colnames(exp.TCGA.anc.fpkm) %like% "TCGA-QG-A5Z1"

# ensembl to hugo
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
annot = getBM(attributes = "hgnc_symbol", "ensembl_gene_id", values = rownames(exp.TCGA.anc.fpkm), mart = mart)
rownames(exp.TCGA.anc.fpkm) = annot$hgnc_symbol

FPKM = as.data.table(exp.TCGA.anc.fpkm@assays@data)
FPKM[, group := NULL]
FPKM[, group_name := NULL]
colnames(FPKM) = colnames(exp.TCGA.anc.fpkm)
FPKM$gene = rownames(exp.TCGA.anc.fpkm)
FPKM = as.data.table(FPKM)
colstotake = c("gene",colnames(FPKM)[colnames(FPKM)!="gene"])
write.table(x = FPKM[, colstotake, with = F], file = "../scratch/FPKM/exp.TCGA.anc.ALLGROUPS.fpkm.txt", sep = "\t", row.names = F, col.names = T, quote = F)

AFRcols = c("gene",unlist(sapply(AFR, FUN = function(x) colnames(FPKM)[colnames(FPKM) %like% x])))
EURcols = c("gene",unlist(sapply(EUR, FUN = function(x) colnames(FPKM)[colnames(FPKM) %like% x])))

EURcols_samp67 = c("gene",sample(unlist(sapply(EUR, FUN = function(x) colnames(FPKM)[colnames(FPKM) %like% x])), size = 67, replace = FALSE))

write.table(x = FPKM[, AFRcols, with = F], file = "../scratch/FPKM/exp.TCGA.anc.AFR.fpkm.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = FPKM[, EURcols, with = F], file = "../scratch/FPKM/exp.TCGA.anc.EUR.fpkm.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(x = FPKM[, EURcols_samp67, with = F], file = "../scratch/FPKM/exp.TCGA.anc.EUR.random.fpkm.txt", sep = "\t", row.names = F, col.names = T, quote = F)
