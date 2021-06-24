# Workflow for differential gene expression stratified by ancestry
# Libraries
require(tidyverse)
require(readr)
require(data.table)
require(TCGAbiolinks)
require(DESeq2)

# Settings
# local: setwd("~/scratch/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/")
setwd("~/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/")

# Cancer Type
typeCancer = "COAD"
typeTCGA = paste("TCGA", typeCancer, sep ="-")

# Import Data
# import clinical test gene sets:
clinicalTests = read.delim(file = ".",
)

# import inferred ancestry in TCGA participants
ancestryTCGA = read.delim(file = "./data/anc/TCGA_consensus_ancestry.tsv",
                          sep = "\t",
                          header = T) %>% as.data.table()

# import 
query <- GDCquery(
  project = typeTCGA,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)
