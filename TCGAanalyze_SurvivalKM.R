#TCGAanalyze_SurvivalKM: Correlating gene expression and Survival Analysis
library(TCGAbiolinks)
# Survival Analysis SA
clinical_patient_Cancer <- GDCquery_clinic("TCGA-LUAD","clinical", save.csv = FALSE)
dataLUADcomplete <- log2(exp.hg38.values)
tokenStop<- 1
tabSurvKMcomplete <- NULL
for( i in 1: round(nrow(dataLUADcomplete)/100)){
  message( paste( i, "of ", round(nrow(dataLUADcomplete)/100)))
  tokenStart <- tokenStop
  tokenStop <-100*i
  tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
                                    dataLUADcomplete,
                                    Genelist = rownames(dataLUADcomplete)[tokenStart:tokenStop],
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












# Selecting only 20 genes for example
dataBRCAcomplete <- log2(exp.hg38.values[1:20,] + 1)

# clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")
clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")


group1 <- TCGAquery_SampleTypes(colnames(dataLUADcomplete), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(dataLUADcomplete), typesample = c("TP"))
clinical_LUAD_m1 = DataFrame(clinical_LUAD$submitter_id,
                             clinical_LUAD$days_to_last_follow_up,
                             clinical_LUAD$days_to_death,
                             clinical_LUAD$vital_status)
TCGAanalyze_SurvivalKM(clinical_LUAD_m1, dataLUADcomplete,Genelist = row[2:500])
clinical_LUAD=clinical_patient_Cancer
