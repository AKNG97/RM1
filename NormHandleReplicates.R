#dir.create("PseudogenesAnalysis")
#This was run on biotin4, /storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps
#### 1. Load dependencies ####

library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)
library(biomaRt)
#library(xCell)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(Hmisc)
library(tidyr)

# 1.1 Load functions

Get_raw_matrix <- function(x,y) {
  z <- cbind(x, y)
  print(dim(z))
  print(head(rownames(z)))
  # rownames(z) <- rowData(x)$gene_name
  # print(head(rownames(z)))
  z <- z[rownames(z) %in% annot_TCGA$gene_id,]
  rownames(z) <- annot_TCGA$HGNC_symbol[match(rownames(z), annot_TCGA$gene_id)]
  print(dim(z))
  print(head(rownames(z)))
  print("Filter matrix")
  dataFilt <- TCGAanalyze_Filtering(tabDF = z,
                                    method = "quantile",
                                    qnt.cut = 0.25)
  threshold <- round(dim(z)[2]/2)
  ridx <- rowSums(dataFilt == 0) <= threshold
  dataFilt <- dataFilt[ridx, ]
  ridx <- rowMeans(dataFilt) >= 10
  dataFilt <- dataFilt[ridx, ]
  z <- z[rownames(z) %in% rownames(dataFilt), ]
  print(dim(z))
  return(z)
}

Get_factors_objects <- function(x,y){
  factors <- x[,c("patient", "primary_diagnosis")] %>%
    bind_rows(y[,c("patient", "primary_diagnosis")])
  rownames(factors) <- factors$barcode
  return(factors)
}

norm <- function(x, y, z) {
  y <- y[y$HGNC_symbol %in% rownames(x),]
  y <- y[match(rownames(x), y$HGNC_symbol),]
  ln.data <- withinLaneNormalization(x, y$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , y$GC, which = "full")
  norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(z$primary_diagnosis))
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  rnas2 <- exprs(mydata2corr1)
  return(rnas2)
}

#### 2.1 Get Biomart annotation  #####
#httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")

features <- c("ensembl_gene_id", "chromosome_name", 
              "start_position", "end_position", "hgnc_symbol",	
              "percentage_gene_gc_content", "gene_biotype", "ensembl_gene_id_version", "hgnc_id")
chrs <- c(1:22, "X", "Y")

annot <- getBM(attributes = features,
               filters = "chromosome_name",
               values = chrs, 
               mart = ensembl)

colnames(annot)<-c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type", "Ensembl_ID_Version", "HGNC_ID")
annot$Length <- abs(annot$End - annot$Start)
annot <- annot[annot$HGNC_symbol != "",]
annot <- annot[!duplicated(annot$ensembl_gene_id),]
dim(annot)
#[1] 40897    10

#### 2.2 Get RNA-seq data from TCGA #####

# #BALL and TALL
# qry.ALL <- GDCquery(project = "TARGET-ALL-P2",
#                 data.category= "Transcriptome Profiling",
#                 data.type = "Gene Expression Quantification",
#                 workflow.type = "STAR - Counts")
# GDCdownload(qry.ALL)
# ALL <- GDCprepare(qry.ALL, summarizedExperiment = TRUE)
# 
# #MM
# qry.MM <- GDCquery(project = "MMRF-COMMPASS",
#                     data.category= "Transcriptome Profiling",
#                     data.type = "Gene Expression Quantification",
#                     workflow.type = "STAR - Counts")
# GDCdownload(qry.MM)
# MM <- GDCprepare(qry.MM, summarizedExperiment = TRUE)
# 
# #AML and Normal bone marrow
# 
# qry.AML_NBM <- GDCquery(project = "TARGET-AML",
#                    data.category= "Transcriptome Profiling",
#                    data.type = "Gene Expression Quantification",
#                    workflow.type = "STAR - Counts",
#                    sample.type = c("Primary Blood Derived Cancer - Bone Marrow",
#                                    "Bone Marrow Normal"))
# # Warning: There are more than one file for the same case. Please verify query results. 
# # You can use the command View(getResults(query)) in rstudio
# # Remove duplicated cases
# qry.results.AML <- getResults(qry.AML_NBM)
# qry.results.AML <- qry.results.AML[!duplicated(qry.results.AML$cases),]
# qry.results.AML <- qry.results.AML[!duplicated(qry.results.AML$cases.submitter_id),]
# # Warning: There are more than one file for the same case. Please verify query results. 
# # You can use the command View(getResults(query)) in rstudio
# 
# #Remove unusual names, then download files and get SummarizedExperiment object
# cases <- qry.results.AML$cases
# dash_count <- str_count(cases, "-")
# qry.results.AML <- qry.results.AML[!(dash_count > 4),]
# qry.AML_NBM[[1,1]] <- qry.results.AML
# 
# GDCdownload(qry.AML_NBM)
# AML_NBM <- GDCprepare(qry.AML_NBM, summarizedExperiment = TRUE)

dir.create("RawData")
# saveRDS(ALL, "RawData/ALL_download_061023.RDS")
# saveRDS(MM, "RawData/MM_download_061023.RDS")
# saveRDS(AML_NBM, "RawData/AML_NBM_download_061023.RDS")

ALL <- readRDS("../ALL_download_061023.RDS")
MM <- readRDS("../MM_download_061023.RDS")
AML_NBM <- readRDS("../AML_NBM_download_061023.RDS")

#### 2.3 Merge annotations from Biomart and TCGA, retain only conserved genes based on Ensembl Version and get PS annotation #####
annot_TCGA <- as.data.frame(rowData(ALL))
annot_TCGA$Ensembl <- gsub("\\.[0-9]*", "", annot_TCGA$gene_id)
dim(annot_TCGA)
annot_TCGA <- annot_TCGA %>% inner_join(annot[,c("ensembl_gene_id","Ensembl_ID_Version","HGNC_symbol", "Type", "Chr", "GC", "Start", "End","Length")], 
                                        by = c("Ensembl" = "ensembl_gene_id"))
dim(annot)
dim(annot_TCGA)

which(duplicated(annot_TCGA$gene_name))
#[1] 35274 37941 38428 38641 39746 39747 39793 39822 39857 39862 40575

which(duplicated(annot_TCGA$HGNC_symbol))
#[1] 35274 37941 38428 38641 39746 39747 39793 39822 39857 39862 40575

duplicated_GeneNames <- annot_TCGA$HGNC_symbol[which(duplicated(annot_TCGA$HGNC_symbol))]
dim(annot_TCGA)
#40704    19
annot_TCGA <- annot_TCGA[!(annot_TCGA$HGNC_symbol %in% duplicated_GeneNames),]

dim(annot_TCGA)
#[1] 40669    19


#Get Pseudogenes conserved on both Biomart and TCGA annotations

length(grep(".*pseudogene", annot_TCGA$gene_type))
#[1] 10304
length(grep(".*pseudogene", annot_TCGA$Type))
#[1] 10262

annot_PS <- annot_TCGA[intersect(grep(".*pseudogene", annot_TCGA$gene_type), 
                                 grep(".*pseudogene", annot_TCGA$Type)),]
dim(annot_PS)
#[1] 10225    19


#### 2.4 Prepare samples #####
#### 2.4.1 Get primary samples of BALL, outputs are BALL_raw and CI_BALL ####
#Check for NAs in relevant clinical variables for subsetting samples

#ALL
dim(ALL)
#[1] 60660   532
which(is.na(ALL$sample_type))
# [1]  13  34  56  92 134 168 194 205 257 264 299 335 419 424 449 452 466 510 516
# [20] 523 524
which(is.na(ALL$definition))
#integer(0)
table(as.factor(ALL$primary_diagnosis))
# Acute lymphocytic leukemia 
# 507

table(as.factor(ALL$definition))
# Primary Blood Derived Cancer - Bone Marrow 
# 387 
# Primary Blood Derived Cancer - Peripheral Blood 
# 76 
# Recurrent Blood Derived Cancer - Bone Marrow 
# 68 
# Recurrent Blood Derived Cancer - Peripheral Blood 
# 1 

#ALL lacks primary diagnosis of previuos version of TCGA,
#also the sample type column is incomplete, but definition seems adequeate to work with
#Let's load previous annotation in order to get BALL and TALL datasets

#Get primary BM samples
ALL_primaryBM <- ALL[,ALL$definition == "Primary Blood Derived Cancer - Bone Marrow"]
#Load old annotation
previous.ALL <- readRDS("../rnas_raw_TAP2.rds")
table(as.factor(previous.ALL$primary_diagnosis))
# Precursor B-cell lymphoblastic leukemia       T lymphoblastic leukemia/lymphoma 
# 267                                     265 
table(as.factor(previous.ALL$sample_type))
# Primary Blood Derived Cancer - Bone Marrow 
# 387 
# Primary Blood Derived Cancer - Peripheral Blood 
# 76 
# Recurrent Blood Derived Cancer - Bone Marrow 
# 68 
# Recurrent Blood Derived Cancer - Peripheral Blood 
# 1

previous.ALL_primaryBM <- previous.ALL[,previous.ALL$definition == "Primary Blood Derived Cancer - Bone Marrow"]
length(intersect(colnames(ALL_primaryBM), colnames(previous.ALL_primaryBM)))
#[1] 387, Primary samples are the same

previous.BALL <- previous.ALL_primaryBM[,previous.ALL_primaryBM$primary_diagnosis == "Precursor B-cell lymphoblastic leukemia"]
length(intersect(colnames(ALL_primaryBM), colnames(previous.BALL)))
#[1] 142, these samples were classified as B cell leukemia

BALL <- ALL_primaryBM[,colnames(ALL_primaryBM) %in% colnames(previous.BALL)]
dim(BALL)
#[1] 60660   142

#Let's check for replicates
which(duplicated(BALL$patient))
#[1]  25  62  69  84 106 111 120 136 141

#Prepare to add them with DESeq2
CI_BALL <- as.data.frame(colData(BALL)) %>% mutate(primary_diagnosis = "BALL")

1:ncol(BALL) == match(colnames(BALL), CI_BALL$barcode)

dds <- DESeqDataSetFromMatrix(countData = assay(BALL),
                              colData = CI_BALL,
                              design = ~ patient)

ddsColl <- collapseReplicates(dds, dds$patient, dds$patient)

#Check if counts were added
CI_BALL$barcode[25]
#[1] "TARGET-10-PARBVI-09A-02R"
#Collapsed data
head(counts(ddsColl[,"TARGET-10-PARBVI"]))
#                     TARGET-10-PARBVI
# ENSG00000000003.15               21
# ENSG00000000005.6                 0
# ENSG00000000419.13             2866
# ENSG00000000457.14              317
# ENSG00000000460.17              165
# ENSG00000000938.13              141
#Replicate 2
head(counts(dds[,dds$barcode[25]]))
#                     TARGET-10-PARBVI-09A-02R
# ENSG00000000003.15                       17
# ENSG00000000005.6                         0
# ENSG00000000419.13                       34
# ENSG00000000457.14                       23
# ENSG00000000460.17                        1
# ENSG00000000938.13                       13
#Replicate 1
head(counts(dds[,"TARGET-10-PARBVI-09A-01R"]))
#                     TARGET-10-PARBVI-09A-01R
# ENSG00000000003.15                        4
# ENSG00000000005.6                         0
# ENSG00000000419.13                     2832
# ENSG00000000457.14                      294
# ENSG00000000460.17                      164
# ENSG00000000938.13                      128

#Data collapsed correctly

#Arrange samples metadata

BALL_raw <- counts(ddsColl)
dim(BALL_raw)
#[1] 60660   133
dim(CI_BALL)
#[1] 142  56
CI_BALL <- CI_BALL[!duplicated(CI_BALL$patient),]
dim(CI_BALL)

1:ncol(BALL_raw) == match(colnames(BALL_raw), CI_BALL$patient)
# [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [13] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
# [25] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [37] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [49] FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
# [61] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [73] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [85] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [97] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [109] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [121] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [133] FALSE

CI_BALL <- CI_BALL[match(colnames(BALL_raw), CI_BALL$patient),]

1:ncol(BALL_raw) == match(colnames(BALL_raw), CI_BALL$patient)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [46] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [76] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [91] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [106] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# [121] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

#Check for NAs in survival features
length(which(is.na(CI_BALL$days_to_last_follow_up)))
#[1] 133
length(which(is.na(CI_BALL$days_to_death)))
#[1] 73

#Let's complete DLF with previuos metadata

previous.BALL.MD <- as.data.frame(colData(previous.BALL)) %>% dplyr::select(patient, days_to_last_follow_up, days_to_death, vital_status) %>%
                                                              dplyr::mutate(P_DLF = days_to_last_follow_up,
                                                                    P_VS = vital_status,
                                                                    P_DTD = days_to_death,)
previous.BALL.MD <- previous.BALL.MD[!duplicated(previous.BALL.MD$patient),]                                                                     
dim(previous.BALL.MD)
#[1] 133   7

CI_BALL <- CI_BALL %>% inner_join(previous.BALL.MD, by = c("patient" = "patient"))
colnames(BALL_raw) <- CI_BALL$patient
#BALL done

#### 2.4.2 Get primary samples of TALL, outputs are TALL_raw and CI_TALL ####
#Get previous TALL from previous metadata
previous.TALL <- previous.ALL_primaryBM[,previous.ALL_primaryBM$primary_diagnosis == "T lymphoblastic leukemia/lymphoma"]
length(intersect(colnames(ALL_primaryBM), colnames(previous.TALL)))
#[1] 245, these samples were classified as T cell leukemia

TALL <- ALL_primaryBM[,colnames(ALL_primaryBM) %in% colnames(previous.TALL)]
dim(TALL)
#[1] 60660   245

which(duplicated(TALL$patient))
#integer(0)

TALL_raw <- assay(TALL)

CI_TALL <- as.data.frame(colData(TALL)) %>% mutate(primary_diagnosis = "TALL")

previous.TALL <- previous.ALL_primaryBM[,previous.ALL_primaryBM$primary_diagnosis == "T lymphoblastic leukemia/lymphoma"]
previous.TALL.MD <- as.data.frame(colData(previous.TALL)) %>% dplyr::select(patient, days_to_last_follow_up, days_to_death, vital_status) %>%
  dplyr::mutate(P_DLF = days_to_last_follow_up,
                P_VS = vital_status,
                P_DTD = days_to_death,)

CI_TALL <- CI_TALL %>% inner_join(previous.TALL.MD, by = c("patient" = "patient"))
1:ncol(TALL_raw) == match(colnames(TALL_raw), CI_TALL$barcode)
colnames(TALL_raw) <- CI_TALL$patient

#TALL done

#### 2.4.3 Get primary samples of MM, outputs are MM_raw and CI_MM ####

#Check MM
table(as.factor(MM$sample_type))
# Primary Blood Derived Cancer - Bone Marrow 
# 764 
# Primary Blood Derived Cancer - Peripheral Blood 
# 12 
# Recurrent Blood Derived Cancer - Bone Marrow 
# 80 
# Recurrent Blood Derived Cancer - Peripheral Blood 
# 3 

which(is.na(MM$sample_type))
#integer(0)

MM <- MM[,MM$sample_type=="Primary Blood Derived Cancer - Bone Marrow"]
which(duplicated(MM$patient))
#integer(0)
table(as.factor(MM$vital_status))
# Alive  Dead 
# 610   154 

MM_raw <- assay(MM)
CI_MM <- as.data.frame(colData(MM)) %>% mutate(primary_diagnosis = "MM")
which(is.na(CI_MM$days_to_last_follow_up))
# integer(0)

1:ncol(MM_raw) == match(colnames(MM_raw), CI_MM$barcode)
colnames(MM_raw) <- CI_MM$patient


#### 2.4.4 Get primary samples of AML, outputs are AML_raw and CI_AML ####

#Check AML and Normal Bone Marrow
table(as.factor(AML_NBM$sample_type))
# Bone Marrow Normal 
# 215 
# Primary Blood Derived Cancer - Bone Marrow 
# 1602 
which(is.na(AML_NBM$sample_type))
#integer(0)

AML <- AML_NBM[,AML_NBM$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
table(as.factor(AML$primary_diagnosis))
# Acute myeloid leukemia, NOS 
# 1485
which(duplicated(AML$patient))
#integer(0)
table(as.factor(AML$vital_status))
# Alive         Dead Not Reported 
# 955          511           16 

previous.AML <- readRDS("../AML_Normal_BM.rds")
previous.AML_primaryBM <- previous.AML[,previous.AML$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
length(intersect(colnames(AML), colnames(previous.AML_primaryBM)))
#[1] 387, Primary samples are the same

AML_raw <- assay(AML)
CI_AML <- as.data.frame(colData(AML)) %>% mutate(primary_diagnosis = "AML")
CI_AML <- CI_AML %>% left_join(as.data.frame(colData(previous.AML)) %>% dplyr::select(days_to_last_follow_up, vital_status, barcode), 
                                by = c("barcode" = "barcode"))
1:ncol(AML_raw) == match(as.data.frame(colnames(AML_raw)))
colnames(AML_raw) <- CI_AML$patient

CI_AML_2 <- CI_AML %>% inner_join(as.data.frame(colData(previous.AML)) %>% dplyr::select(days_to_last_follow_up, vital_status, barcode), 
                               by = c("barcode" = "barcode"))

#### 2.4.5 Get samples of NBM, outputs are Remissioned_NBM_raw, NBM_raw, CI_RemissionedNBM and CI_NBM ####
NBM <- AML_NBM[,AML_NBM$sample_type == "Bone Marrow Normal"]
table(as.factor(NBM$primary_diagnosis))
# Acute myeloid leukemia, NOS 
# 145 
#Hence, there are "normal" samples from AML patients
which(duplicated(NBM$patient))
#integer(0)
table(as.factor(NBM$vital_status))
# Alive  Dead 
# 90    55 
#From these "recovered patients" 55 died later.
#¿Should I manage these samples separately from samples without a primary diagnosis?
#Let's subset them, check on that later

Remissioned_NBM <- NBM[, !is.na(NBM$primary_diagnosis)]
dim(Remissioned_NBM)
# [1] 60660   145
Remissioned_NBM_raw <- assay(Remissioned_NBM)
CI_RemissionedNBM <- as.data.frame(colData(Remissioned_NBM)) %>% mutate(primary_diagnosis = "Remission_Normal_Bone_Marrow_AML")

1:ncol(Remissioned_NBM_raw) == match(colnames(Remissioned_NBM_raw), CI_RemissionedNBM$barcode)
colnames(Remissioned_NBM_raw) <- CI_RemissionedNBM$patient

NBM <- NBM[, is.na(NBM$primary_diagnosis)]
dim(NBM)
#[1] 60660    70
NBM_raw <- assay(NBM)
CI_NBM <- as.data.frame(colData(NBM)) %>% mutate(primary_diagnosis = "Normal_Bone_Marrow")
  
1:ncol(NBM_raw) == match(colnames(NBM_raw), CI_NBM$barcode)
colnames(NBM_raw) <- CI_NBM$patient

#### 2.4.6 Save raw data ####

saveRDS(BALL_raw, "RawData/BALL_primaryBM_raw.RDS")
saveRDS(TALL_raw, "RawData/TALL_primaryBM_raw.RDS")
saveRDS(MM_raw, "RawData/MM_primaryBM_raw.RDS")
saveRDS(AML_raw, "RawData/AML_primaryBM_raw.RDS")
saveRDS(Remissioned_NBM_raw, "RawData/Remissioned_NBM_primaryBM_raw.RDS")
saveRDS(NBM_raw, "RawData/NBM_primaryBM_raw.RDS")

saveRDS(CI_NBM, "CI_primaryBM_NBM.RDS")
saveRDS(CI_RemissionedNBM, "CI_primaryBM_RemissionedNBM.RDS")
saveRDS(CI_TALL, "CI_primaryBM_TALL.RDS")
saveRDS(CI_BALL, "CI_primaryBM_BALL.RDS")
saveRDS(CI_AML, "CI_primaryBM_AML.RDS")
#saveRDS(CI_AML, "CI_primaryBM_AML_with_previousData.RDS")
#saveRDS(CI_AML_2, "CI_primaryBM_AML_with_previousData_match.RDS")
saveRDS(CI_MM, "CI_primaryBM_MM.RDS")

saveRDS(annot_TCGA, "annot_TCGA.RDS")
saveRDS(annot_PS, "annot_PS.RDS")

#### 2.5 Get raw matrices, merge with NBM, then filter low expressedd genes across most samples ####

BALL_F <- Get_raw_matrix(BALL_raw, NBM_raw)
TALL_F <- Get_raw_matrix(TALL_raw, NBM_raw)
AML_F <- Get_raw_matrix(AML_raw, NBM_raw)
MM_F <- Get_raw_matrix(MM_raw, NBM_raw)
R_NBM_F <- Get_raw_matrix(Remissioned_NBM_raw, NBM_raw)

#### 3 Get norm data ####
#Get factors objs
R_NBM_factors <- Get_factors_objects(CI_RemissionedNBM, CI_NBM)
TALL_factors <- Get_factors_objects(CI_TALL, CI_NBM)
BALL_factors <- Get_factors_objects(CI_BALL, CI_NBM)
AML_factors <- Get_factors_objects(CI_AML, CI_NBM)
MM_factors <- Get_factors_objects(CI_MM, CI_NBM)

#Norm data
BALL_norm <- norm(BALL_F, annot_TCGA, BALL_factors)
TALL_norm <- norm(TALL_F, annot_TCGA, TALL_factors)
R_NBM_norm <- norm(R_NBM_F, annot_TCGA, R_NBM_factors)
MM_norm <- norm(MM_F, annot_TCGA, MM_factors)
AML_norm <- norm(AML_F, annot_TCGA, AML_factors)

#Save data
dir.create("NormData/")
saveRDS(BALL_norm, "NormData/BALL_norm.RDS")
saveRDS(TALL_norm, "NormData/TALL_norm.RDS")
saveRDS(R_NBM_norm, "NormData/R_NBM_norm.RDS")
saveRDS(MM_norm, "NormData/MM_norm.RDS")
saveRDS(AML_norm, "NormData/AML_norm.RDS")

#### 3.1 Get PS expression ####
R_NBM_PS <- R_NBM_norm[rownames(R_NBM_norm) %in% annot_PS$HGNC_symbol,]
AML_PS <- AML_norm[rownames(AML_norm) %in% annot_PS$HGNC_symbol,]
TALL_PS <- TALL_norm[rownames(TALL_norm) %in% annot_PS$HGNC_symbol,]
BALL_PS <- BALL_norm[rownames(BALL_norm) %in% annot_PS$HGNC_symbol,]
MM_PS <- MM_norm[rownames(MM_norm) %in% annot_PS$HGNC_symbol,]

dim(R_NBM_PS)
dim(AML_PS)
dim(TALL_PS)
dim(BALL_PS)
dim(MM_PS)

saveRDS(R_NBM_PS, "NormData/R_NBM_PS.RDS")
saveRDS(AML_PS, "NormData/AML_PS.RDS")
saveRDS(TALL_PS, "NormData/TALL_PS.RDS")
saveRDS(BALL_PS, "NormData/BALL_PS.RDS")
saveRDS(MM_PS, "NormData/MM_PS.RDS")

#### 4 Get GCNs #####
#Falto correr esto, solo hice TALL y su control
SP_correlation <- function(x){
  
  net <- rcorr(t(x), type = "spearman")
  
  coexpr_val <- net$r
  coexpr_val[lower.tri(coexpr_val, diag = TRUE)] <- NA
  coexpr_val <- coexpr_val %>% as.data.frame() %>% mutate(Gene = rownames(coexpr_val)) %>%
    pivot_longer(!Gene, names_to = "Target", values_to = "Sp_corr") %>% 
    rename(Source = Gene)  %>% drop_na(Sp_corr) %>%
    mutate(ID = paste(pmin(Source,Target),pmax(Source,Target),sep="-"),
           Abs_corr = abs(Sp_corr))
  
  coexpr_P <- net$P
  coexpr_P[lower.tri(coexpr_P, diag = TRUE)] <- NA
  coexpr_P <- coexpr_P %>% as.data.frame() %>% mutate(Gene = rownames(coexpr_P)) %>%
    pivot_longer(!Gene, names_to = "Target", values_to = "p_value") %>% 
    rename(Source = Gene)  %>% drop_na(p_value) %>% filter(p_value < 1e-08) %>%
    mutate(ID = paste(pmin(Source,Target),pmax(Source,Target),sep="-")) 
  
  coexpr_full <- coexpr_val %>% inner_join(coexpr_P[,c("ID", "p_value")], by = c("ID" = "ID"))
  
  return(coexpr_full)
  
}

dir.create("GCNs/")
#Consider that the last 70 samples of each matrix are NBM
# R_NBM_PS_GCN <- SP_correlation(R_NBM_PS[,1:(length(colnames(R_NBM_PS)) - 70)])
# dim(R_NBM_PS_GCN)
# saveRDS(R_NBM_PS_GCN, "R_NBM_PS_GCN.RDS")

# BALL_PS_GCN <- SP_correlation(BALL_PS[,1:(length(colnames(BALL_PS)) - 70)])
# dim(BALL_PS_GCN)
# saveRDS(BALL_PS_GCN, "BALL_PS_GCN.RDS")

TALL_PS_GCN <- SP_correlation(TALL_PS[,1:(length(colnames(TALL_PS)) - 70)])
dim(TALL_PS_GCN)
saveRDS(TALL_PS_GCN, "GCNs/TALL_PS_GCN.RDS")

# MM_PS_GCN <- SP_correlation(MM_PS[,1:(length(colnames(MM_PS)) - 70)])
# dim(MM_PS_GCN)
# saveRDS(MM_PS_GCN, "MM_PS_GCN.RDS")
# 
# AML_PS_GCN <- SP_correlation(AML_PS[,1:(length(colnames(AML_PS)) - 70)])
# dim(AML_PS_GCN)
# saveRDS(AML_PS_GCN, "AML_PS_GCN.RDS")

#Normal
R_NBM_vs_NBM_PS_GCN <- SP_correlation(R_NBM_PS[,(length(colnames(R_NBM_PS)) - 69):length(colnames(R_NBM_PS))])
dim(R_NBM_vs_NBM_PS_GCN)
saveRDS(R_NBM_vs_NBM_PS_GCN, "R_NBM_vs_NBM_PS_GCN.RDS")

BALL_vs_NBM_PS_GCN <- SP_correlation(BALL_PS[,(length(colnames(BALL_PS)) - 69):length(colnames(BALL_PS))])
dim(BALL_vs_NBM_PS_GCN)
saveRDS(BALL_vs_NBM_PS_GCN, "BALL_vs_NBM_PS_GCN.RDS")

TALL_vs_NBM_PS_GCN <- SP_correlation(TALL_PS[,(length(colnames(TALL_PS)) - 69):length(colnames(TALL_PS))])
dim(TALL_vs_NBM_PS_GCN)
saveRDS(TALL_vs_NBM_PS_GCN, "GCNs/TALL_vs_NBM_PS_GCN.RDS")

MM_vs_NBM_PS_GCN <- SP_correlation(MM_PS[,(length(colnames(MM_PS)) - 69):length(colnames(MM_PS))])
dim(MM_vs_NBM_PS_GCN)
saveRDS(MM_vs_NBM_PS_GCN, "MM_vs_NBM_PS_GCN.RDS")

AML_vs_NBM_PS_GCN <- SP_correlation(AML_PS[,(length(colnames(AML_PS)) - 69):length(colnames(AML_PS))])
dim(AML_vs_NBM_PS_GCN)
saveRDS(AML_vs_NBM_PS_GCN, "AML_vs_NBM_PS_GCN.RDS")
