#https://zhuanlan.zhihu.com/p/32607443
#library("TCGAbiolinks")
samplebarcode=c("TCGA-05-5429-01A-01R-1628-07",
"TCGA-55-6984-01A-11R-1949-07","TCGA-35-4122-01A-01R-1107-07",
"TCGA-55-A48X-01A-11R-A24H-07","TCGA-64-1678-01A-01R-0946-07",
"TCGA-99-AA5R-01A-11R-A39D-07","TCGA-64-1680-01A-02R-0946-07",
"TCGA-55-8087-01A-11R-2241-07","TCGA-69-7979-01A-11R-2187-07",
"TCGA-91-6836-01A-21R-1858-07","TCGA-44-2661-11A-01R-1758-07",
"TCGA-91-6836-11A-01R-1858-07","TCGA-55-6981-11A-01R-1949-07",
"TCGA-50-6595-11A-01R-1858-07","TCGA-55-6983-11A-01R-1949-07",
"TCGA-44-5645-11A-01R-1628-07","TCGA-50-5932-11A-01R-1755-07",
"TCGA-49-4490-11A-01R-1858-07","TCGA-55-6970-11A-01R-1949-07")

query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM",
                  barcode = samplebarcode)
GDCdownload(query,files.per.chunk = 10)
#library(SummarizedExperiment)
exp.hg38 <- GDCprepare(query)
samples.information=colData(exp.hg38)
#geneInfo=geneInfoHT，default其实是geneInfo，但由于我们前面选择的是HTseq，所以要选择geneInfoHT
dataNorm <- TCGAanalyze_Normalization(tabDF = exp.hg38, geneInfo = geneInfoHT)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  
                                  method ="quantile",
                                  
                                  qnt.cut = 0.25)
#定义对照组（这里的对照组是Solid normal tissue）
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   
                                   typesample = c("NT"))

#定义肿瘤组
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   
                                   typesample = c("TP"))
#进行DEA分析，用到DEA()

dataDEGs <- TCGAanalyze_DEA(mat1 =dataFilt[,samplesNT],
                              
                              mat2 = dataFilt[,samplesTP],
                              
                              Cond1type = "Normal",
                              
                              Cond2type = "Tumor",
                            pipeline="edgeR"
                              
                              #fdr.cut = 0.01 ,
                              
                            #  logFC.cut = 1,
                              
                              #method = "glmLRT"
                              )
write.csv(dataDEGsFiltLevel,file="DEA_LUAD.csv")
