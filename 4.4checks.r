##some checks

#call out data and packages-----------------------------------------------------------
rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
setwd(working_path)

library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pROC)
library(ROCR)
library(caret)
library(corrplot)
library(DescTools)
library(flextable)

source("/gpfs3/well/emberson/users/hma817/codes/R_TablesFunction_27Sep2022_TL.r")
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")

data<-readRDS(paste(data_path,"analysis_data_16May2023.rds",sep="/"))

variants_included<-readRDS(paste(data_path,"variants_included.rds",sep="/"))
variant_t<-data.frame(t(variants_included))
colnames(variant_t)<-variant_t[1,]
variant_t<-variant_t[2,]
variant_t<-cbind(rownames(variant_t),variant_t)
colnames(variant_t)[1]<-NA
colnames(variant_t)[8]<-"custom_Oni_Orisan"
variant_t[1,1]<-"Number of SNPs"
PCs<-c(paste0("PC",1:7))
full_adjustments = c("WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline")

#check analysis results----------------------------------------------------

#correlation with AGE
age<-cor(data[,c(3,33:39)],method = "pearson",use="complete.obs")[1,2:8]

#check difference between gender
gender<-NA
for (i in c(33:39)){
  a<-t.test(data[,i]~data$SEX)
  gender[i-32]<-a$p.value
  names(gender[i-32])<-colnames(data)[i]
  print(colnames(data)[i])
  print(a$p.value)
}
#pvalue<0.05 for PGS000018 and PGS001780
t.test(data$PGS000018~data$SEX)
t.test(data$PGS001780~data$SEX)

#correlation with PC1
pc1<-cor(data[,c(24,33:39)],method = "pearson",use="complete.obs")[1,2:8]

#correlation with PC2
pc2<-cor(data[,c(25,33:39)],method = "pearson",use="complete.obs")[1,2:8]

#correlation with PC3
pc3<-cor(data[,c(26,33:39)],method = "pearson",use="complete.obs")[1,2:8]


#correlation with WHRATIO
whr<-cor(data[,c(14,33:39)],method = "pearson",use="complete.obs")[1,2:8]

#correlation with SBP
sbp<-cor(data[,c(16,33:39)],method = "pearson",use="complete.obs")[1,2:8]

#correlation with DBP
dbp<-cor(data[,c(17,33:39)],method = "pearson",use="complete.obs")[1,2:8]





#check dfference between diabetes at baseline
t2d<-NA
for (i in c(33:39)){
  a<-t.test(data[,i]~data$diabetes_at_baseline)
  t2d[i-32]<-a$p.value
  names(t2d[i-32])<-colnames(data)[i]
  print(colnames(data)[i])
  print(a$p.value)
}

check_table<-rbind(age,pc1,pc2,pc3,whr,sbp,dbp,gender,t2d)



#modelling results check---------------------------------
EPA_prevalent<-readRDS(paste(data_path,"/Prevalent_CHD_EPA_model_outcomes_09May_2023.rds",sep=""))
EPA_cont<-EPA_prevalent[[1]]
EPA_cat<-EPA_prevalent[[2]]


#self-creating function check------------------------
run_glm(model_type = "simple",data=data,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",pgs_name = "PGS000018")
glm(prevalent_CHD_EPA~custom_PGS000337,data=data,family=binomial(link='logit'))


run_glm(model_type = "partial",data=data,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",pgs_name = "PGS000018",partial_adjustment = c("AGE","SEX"))
glm(prevalent_CHD_EPA~+AGE+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7,data=data,family=binomial(link='logit'))



run_glm(model_type = "full",data=data,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",pgs_name = "PGS000018",partial_adjustment = c("AGE","SEX"))
glm(prevalent_CHD_EPA~PGS000018+AGE+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+WHRATIO+SBP+DBP+EDU_LEVEL+smokegp+diabetes_at_baseline,data=data,family=binomial(link='logit'))



# check model output-------------------

exp(coef(run_glm(model_type = "simple",data=data,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",pgs_name = "PGS001780")))["PGS001780"]

exp(coef(run_glm(model_type = "partial",data=data,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",pgs_name = "PGS001780",partial_adjustment = c("AGE","SEX"))))["PGS001780"]

exp(coef(run_glm(model_type = "full",data=data,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",pgs_name = "PGS001780",partial_adjustment = c("AGE","SEX"))))["PGS001780"]


exp(coef(glm(prevalent_CHD_EPA~PGS001780,data=data,family=binomial(link='logit'))))
exp(coef(glm(prevalent_CHD_EPA~PGS001780+AGE+SEX,data=data,family=binomial(link='logit'))))

exp(coef(glm(prevalent_CHD_EPA~PGS001780+AGE+PC1,data=data,family=binomial(link='logit'))))

exp(coef(glm(prevalent_CHD_EPA~PGS001780+SEX,data=data,family=binomial(link='logit'))))

exp(coef(glm(prevalent_CHD_EPA~PGS001780+SEX+PC1,data=data,family=binomial(link='logit'))))




exp(coef(glm(prevalent_CHD_EPA~PGS001780+AGE+SEX+PC1,data=data,family=binomial(link='logit'))))

cor(data$AGE,data$PC1,method = "pearson",use="complete.obs")
t.test(data$PC1~data$SEX)
