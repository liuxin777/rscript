#survival analysis
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(survival)
#library(survminer)
#仅选取几个sample试验
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
exp.hg38 <- GDCprepare(query)
exp.hg38.values <- assay(exp.hg38)
#选取1：50个基因做试验
exp.hg38.values=exp.hg38.values[1:50,] 
clinical_LUAD=GDCquery_clinic("TCGA-LUAD","clinical",save.csv = F)
#这里的clinical_LUAD就是完整的临床数据，不需要再download
#可以通过head查看对象的结构
#head(clinical_LUAD[1:5,1:5])


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
#仅保留癌组织的表达量的sample ?
exp.values=exp.hg38.values[,sign=="11A"] 
rm(exp.hg38,samplebarcode,query,i,tmp,exp.hg38.values) 
#现在只有clinical_LUAD,exp.values两个data; 和   patient_id,sign两个values
#生成新的病人号
patient_id=substr(colnames(exp.values),1,12) 
#把exp.values的列名改成submitter_id，类似于TCGA-44-2661
colnames(exp.values)=patient_id
##temp=substr(colnames(exp.hg38.values),14,15)
##gene_name_exp = exp.hg38.values[rownames(exp.hg38.values)%in%gene_name$gene_name,]
##rownames(gene_name_exp)
##gene_name_exp_carcer = exp.hg38.values_targeted_gene[,sign == "01A"]

#整理用于计算生存分析的数据代码
newclinical=data.frame(clinical_LUAD,row.names = clinical_LUAD$submitter_id)
clinical_data=newclinical[patient_id,]
clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")
clinical_data_m1 = DataFrame(clinical_data$submitter_id,
                             clinical_data$tumor_stage,
                             clinical_data$days_to_last_follow_up,
                             clinical_data$days_to_death)
#head(clinical_LUAD_m1[1:5,])
survival_day=c()
for (i in 1:nrow(clinical_data_m1)){
  survival_day[i] = ifelse(is.na(clinical_data_m1$clinical_data.days_to_last_follow_up[i]),clinical_data$days_to_death[i],clinical_data_m1$clinical_data.days_to_last_follow_up[i])
}
survival_state = c()
for (i in 1:nrow(clinical_data_m1)){
  survival_state[i] = ifelse(is.na(clinical_data_m1$clinical_data.days_to_last_follow_up[i]),1,0)
}
submitter_id=clinical_data$submitter_id
clinical_data=data.frame(survival_state,survival_day,submitter_id)

#head(survival_data)
##temp=substr(colnames(exp.hg38.values),14,15)
##exp.hg38.values=exp.hg38.values[,temp=="01"]
#http://chuansong.me/n/2066050952022
# 创造适合生存分析的数据框
chi=c()
p.val=c()
indices=c()
group=c()
#j对基因循环，i对基因表达量循环用于分组，即生成group
#ma是随便取的，ma不是矩阵，是一个数据框，这个数据框的每一行是一个基因
#其实我感觉90行跟95行可以合并
for(j in 1:nrow(exp.values)){
ma=DataFrame(clinical_data,value=exp.values[j,])
for (i in 1:nrow(ma)){
  group[i] <- ifelse(ma$value[i]<=median(ma$value),0,1) 
  #0 for low , 1 for high
}
ma=DataFrame(ma,group)
re=survdiff(Surv(survival_day,survival_state)~group,data=ma)
chi[j]=re$chisq
p.val[j] = 1 - pchisq(re$chisq, length(re$n) - 1)
indices[j]=j
}
#把循环产生的indices,chi,p.val加上基因名称生成数据框
result=DataFrame(indices,rownames(exp.values)[indices],chi,p.val)
#按照p.val排序，小的在上，大的在下
result=result[order(p.val),]
write.table(result,sep="\t",file="sa.txt",row.names = F)
##sfit <- survfit(Surv(survival_day,survival_state)~group,data=ma)
##ggsurvplot(sfit, conf.int=F, pval=TRUE)
