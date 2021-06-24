# Workflow for differential gene expression stratified by ancestry
# Libraries
require(tidyverse)
require(readr)
require(data.table)
require(TCGAbiolinks)
#require(DESeq2)
require(EDAseq)

# Settings
# local: setwd("~/scratch/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/")
setwd("~/branch/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/")

# Cancer Type
typeCancer = "COAD"
typeTCGA = paste("TCGA", typeCancer, sep ="-")

# Import Data
# import clinical test gene sets:
clinicalTests = read.delim(file = "./data/gene_sets/OncotypeDX_Colon_RS.gmx.tsv",
                           sep = "\t",
                           header = T) %>% as.data.table()

# import inferred ancestry in TCGA participants:
ancestryTCGA = read.delim(file = "./data/anc/TCGA_consensus_ancestry.tsv",
                          sep = "\t",
                          header = T) %>% as.data.table()

# import TCGA data
query <- GDCquery(
  project = typeTCGA,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  experimental.strategy = "RNA-Seq",
  sample.type = c("Primary Tumor"),
  legacy = TRUE)

GDCdownload(query)
data <- GDCprepare(query, 
                   summarizedExperiment = TRUE)

# bin participants by ancestry:
AFR = ancestryTCGA[tumor_type == typeCancer & consensus_ancestry %in% c("afr_admix","afr")]$patient
EUR = ancestryTCGA[tumor_type == typeCancer & consensus_ancestry %in% c("eur_admix","eur")]$patient

# identify samples by participants:
AFRcols = unlist(sapply(AFR, FUN = function(x) colnames(data)[colnames(data) %like% x]))
EURcols = unlist(sapply(EUR, FUN = function(x) colnames(data)[colnames(data) %like% x]))

# stratify data by ancestry:
queryAFR <- GDCquery(
  project = typeTCGA,
  barcode = AFRcols,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  experimental.strategy = "RNA-Seq",
  sample.type = c("Primary Tumor"),
  legacy = TRUE)
GDCdownload(queryAFR)
dataAFR <- GDCprepare(queryAFR, 
                   summarizedExperiment = TRUE)

queryEUR <- GDCquery(
  project = typeTCGA,
  barcode = EURcols,
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  experimental.strategy = "RNA-Seq",
  sample.type = c("Primary Tumor"),
  legacy = TRUE)
GDCdownload(queryEUR)
dataEUR <- GDCprepare(queryEUR, 
                   summarizedExperiment = TRUE)

# data prep:
dataPrepEUR <- TCGAanalyze_Preprocessing(
  object = dataEUR,
  cor.cut = 0.6,    
  datatype = "raw_count",
  filename = "coad_EUR.png")

dataPrepAFR <- TCGAanalyze_Preprocessing(
  object = dataAFR,
  cor.cut = 0.6,    
  datatype = "raw_count",
  filename = "coad_AFR.png")

# normalize gene expression:
dataNormEUR <- TCGAanalyze_Normalization(
  tabDF = dataPrepEUR,
  geneInfo = TCGAbiolinks::geneInfo,
  method = "gcContent")

dataNormAFR <- TCGAanalyze_Normalization(
  tabDF = dataPrepAFR,
  geneInfo = TCGAbiolinks::geneInfo,
  method = "gcContent")

# differential gene expression
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataNormEUR,
  mat2 = dataNormAFR,
  Cond1type = "white",
  Cond2type = "black",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT")

# number of differentially expressed genes (DEG)
nrow(dataDEGs)
# dataDEGs that showed logFC > 1 and FDR < 0.01 for example, can be considered significantly up-regulated genes in cancer tissues by ancestry 

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"EUR","AFR",
                                          dataNormEUR, dataNormAFR)
write.table(dataDEGsFiltLevel, file = "./data/DEGs_EUR_AFR.tsv", sep = "\t", row.names = F, quote = F)

# label genes used in clinical test:
#clinicalTestGenes = clinicalTests %>% pivot_longer(cols = colnames(clinicalTests), names_to = "pathway", values_to = "gene")
#clinicalTestGenesAll = data.table(test = "test", gene = unique(clinicalTestGenes$gene))
#dataDEGsFiltLevel = merge(dataDEGsFiltLevel,clinicalTestGenesAll, by.x = "mRNA", by.y = "gene", all.x=T)
