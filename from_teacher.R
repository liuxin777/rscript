library(SummarizedExperiment)
query.exp.hg38 <- GDCquery(project = "TCGA-STAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38,files.per.chunk = 50)
exp.hg38 <- GDCprepare(query = query.exp.hg38)
exp.hg38.values <- assay(exp.hg38)
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
write.csv(exp.hg38.values,file = "stad_exp_hg38_FPKM.csv")

query.exp.hg38 <- GDCquery(project = "TCGA-STAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query.exp.hg38,files.per.chunk = 50)
exp.hg38 <- GDCprepare(query = query.exp.hg38)
exp.hg38.values <- assay(exp.hg38)
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
write.csv(exp.hg38.values,file = "stad_exp_hg38_FPKM_UQ.csv")

query.exp.hg38 <- GDCquery(project = "TCGA-STAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
GDCdownload(query.exp.hg38,files.per.chunk = 50)
exp.hg38 <- GDCprepare(query = query.exp.hg38)
exp.hg38.values <- assay(exp.hg38)
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
write.csv(exp.hg38.values,file = "stad_exp_hg38_htseq_counts.csv")

query.exp <- GDCquery(project = "TCGA-STAD", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq")
GDCdownload(query.exp,files.per.chunk = 50)
stad.hg19 <- GDCprepare(query = query.exp)
stad.hg19.exp <- assay(stad.hg19)
rownames(stad.hg19.exp) <- values(stad.hg19)$gene_id
write.csv(stad.hg19.exp,file = "stad_exp_hg19_rsem.csv")

query.exp <- GDCquery(project = "TCGA-STAD", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expres


query.exp.hg38 <38 <- GDCquery(project = "TCGA-LUAD", 
                              data.category = " = "Transcriptome Profiling", 
                      data.type = " = "Gene Expression Quantification", 
                      workflow.type = " = "HTSeq - FPKM")
GDCdownload(ad(query.exp.hg38,f38,38,files.per.chunk = 1 = 1)1)
exp.hg38 <38 <- GDCprepare(query =  = query.exp.hg38)38)8)
exp.hg38.values <-  <- assay(ay(exphel.hg38)38)
rownames(es(exp.hg38.values) <-) <- values(es(exp.hg38)$38)$external_gene_name


library(TCGAbiolinks)
query.exp.hg19 <- GDCquery(project = "TCGA-LUAD", 
                           legacy=TRUE,
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform="Illumina HiSeq",
                           file.type="results"
                           barcode=c("")
                           )
GDCdownload(query.exp.hg19,method="api",files.per.chunk = 10)
exp.hg19 <- GDCprepare(query.exp.hg19)


# You can define a list of samples to query and download providing relative TCGA barcodes.
listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
                 "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
                 "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
                 "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
                 "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07")


# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)


# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)


# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)


BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")


# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseqSE)

installed.packages()[, c("Package", "LibPath")]

# Downstream analysis using gene expression data  
# TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
# save(dataBRCA, geneInfo , file = "dataGeneExpression.rda")
library(TCGAbiolinks)


# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo =  geneInfo)


# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)


# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))


# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))


# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")


# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])









#TCGAanalyze_SurvivalKM: Correlating gene expression and Survival Analysis
library(TCGAbiolinks)
# Survival Analysis SA


clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
dataBRCAcomplete <- log2(BRCA_rnaseqv2)


tokenStop<- 1


tabSurvKMcomplete <- NULL


for( i in 1: round(nrow(dataBRCAcomplete)/100)){
  message( paste( i, "of ", round(nrow(dataBRCAcomplete)/100)))
  tokenStart <- tokenStop
  tokenStop <-100*i
  tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
                                    dataBRCAcomplete,
                                    Genelist = rownames(dataBRCAcomplete)[tokenStart:tokenStop],
                                    Survresult = F,
                                    ThreshTop=0.67,
                                    ThreshDown=0.33)
  
  tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
}


tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.01,]
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]


tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
  rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
  ]








#代码
#整理转录组的患者的ID号
#同时提取  癌组织的ID号
sign =c()
patient_id = colnames(exp.hg38.values)
for (i in 1:length(colnames(exp.hg38.values))){
  tmp = strsplit2(as.character(patient_id[i]),split = "-")
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "_")
  patient_id[i] = tmp
}
colnames(exp.hg38.values) = patient_id

#gene_name_exp = exp.hg38.values[rownames(exp.hg38.values)%in%gene_name$gene_name,]
#rownames(gene_name_exp)
gene_name_exp_carcer = exp.hg38.values_targeted_gene[,sign == "01A"]

#整理用于计算生存分析的数据代码
clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")

clinical_LUAD_m1 = DataFrame(clinical_LUAD$submitter_id,
                             clinical_LUAD$tumor_stage,
                             clinical_LUAD$days_to_last_follow_up,
                             clinical_LUAD$days_to_death)

survial_day=c()
for (i in 1:nrow(clinical_LUAD_m1)){
  survial_day[i] = ifelse(is.na(clinical_LUAD_m1$clinical_LUAD.days_to_last_follow_up[i]),clinical_LUAD$days_to_death[i],clinical_LUAD_m1$clinical_LUAD.days_to_last_follow_up[i])
}
survial_state = c()
for (i in 1:nrow(clinical_LUAD_m1)){
  survial_state[i] = ifelse(is.na(clinical_LUAD_m1$clinical_LUAD.days_to_last_follow_up[i]),1,0)
}

