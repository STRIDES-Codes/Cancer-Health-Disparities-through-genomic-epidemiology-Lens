install.packages("readr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install(c("TCGAbiolinks","TCGAutils","pathview","SummarizedExperiment","clusterProfiler","org.Hs.eg.d"))


library(readr)
TCGA_ANS <- read_csv("scratch/Cancer-Health-Disparities-through-genomic-epidemiology-Lens/data/anc/TCGA_consensus_ancestry.csv")
TCGA_ANS <-as.data.frame(TCGA_ANS)
TCGA_ANS <- TCGA_ANS[TCGA_ANS$tumor_type =="COAD",]
TCGA_ANS  <- TCGA_ANS[,c(1,3)]
colnames(TCGA_ANS)<-c("patcode","Ances")


query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  sample.type = c("Primary Tumor"),
  legacy = TRUE)

avail <- query[[1]][[1]]$cases
avail <- as.data.frame(avail)
avail$patcode <- substring(avail$avail,1,12)
colnames(avail)<-c("barcode", "patcode")
All_barcode <- (merge(TCGA_ANS, avail, by = 'patcode'))
Black<- unique(All_barcode[which(All_barcode$Ances == "AA"),]$barcode)
Black<-unique(Black)
Black <- as.data.frame(Black)
colnames(Black)<-"barcode"
Black <- Black[Black$barcode !="TCGA-A6-5661-01B-05R-2302-07",]


White<- unique(All_barcode[which(All_barcode$Ances == "EA"),]$barcode)
White<-unique(White)
White <- as.data.frame(White)
colnames(White)<-"barcode"
White <- White[White$barcode !="TCGA-A6-2672-01B-03R-2302-07",]

queryB <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  sample.type = "Primary Tumor",
  barcode = Black,
  legacy = TRUE)


# We will use only 20 samples to make the example faster
queryB$results[[1]] <-  queryB$results[[1]][1:20,]  

GDCdownload(queryB)
queryB$results[[1]] <- unique(queryB$results[[1]])[1:20,]  
coadBB.exp <- GDCprepare(
  query = queryB, 
  summarizedExperiment = TRUE)

dataPrep_BB <- TCGAanalyze_Preprocessing(
  object = coadBB.exp,
  cor.cut = 0.6,    
  datatype = "raw_count",
  filename = "coad_BB.png")



###############################EE


queryEE <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "results", 
  sample.type = "Primary Tumor",
  barcode = White,
  legacy = TRUE)


# We will use only 20 samples to make the example faster
queryEE$results[[1]] <-  queryEE$results[[1]][1:20,]  

GDCdownload(queryEE)
coadEE.exp <- GDCprepare(
  query = queryEE, 
  save = TRUE, 
  summarizedExperiment = TRUE, 
  save.filename = "coadEE.rda")

dataPrep_EE <- TCGAanalyze_Preprocessing(
  object = coadEE.exp,
  cor.cut = 0.6,    
  datatype = "raw_count",
  filename = "coad_EE.png")

#####################################
dataNormEE <- TCGAanalyze_Normalization(
  tabDF = dataPrep_EE,
  geneInfo = TCGAbiolinks::geneInfo,
  method = "gcContent")

dataNormBB <- TCGAanalyze_Normalization(
  tabDF = dataPrep_BB,
  geneInfo = TCGAbiolinks::geneInfo,
  method = "gcContent")



dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataNormEE,
  mat2 = dataNormBB,
  Cond1type = "Whitw",
  Cond2type = "Black",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT")

# Number of differentially expressed genes (DEG)
nrow(dataDEGs)


#dataDEGs that showed logFC > 1 and FDR < 0.01 for example, can be considered significantly Up-regulated genes in cancer (EE) tissues. 


# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataNormEE, dataNormBB)

