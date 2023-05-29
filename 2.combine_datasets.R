#this file is to match PRS with individual baseline info
rm(list=ls())

#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
results_path<-paste(working_path,'/graphs',sep='')
prs_path<-paste(working_path,'/PRS_calculated',sep='')
phenotype_path<-"/well/emberson/projects/mcps/data/phenotypes/"

setwd(working_path)

library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pROC)
library(ROCR)
library(caret)

#prs data readin-------------------------------
name<-NA
number<-NA


prs_reading<-data.frame(NA)
readin_prs<-function(prs){
  
  for ( i in c(1:length(prs))) {
  cat("\n read in PRS data for ",prs[i], ".......")
  specific_prs_path<-paste(prs_path,prs[i],"score",sep="/")
  cat("\n path is ", specific_prs_path)
  setwd(specific_prs_path)
  dat=read.table("aggregated_scores.txt.gz", 
                 sep ="", header = TRUE, dec =".",fill=T)
  
  if(i==1){
    prs_reading<-data.frame(dat$IID,dat[4])
    colnames(prs_reading)<-c("IID",prs[i])
  }else{
    prs_reading<-merge(prs_reading,dat[,c(2,4)],by="IID")
    colnames(prs_reading)[i+1]<-prs[i]
  }
  
  
  
  name[i]<-prs[i]
  if(mean(dat$DENOM)==dat$DENOM[1]){
    number[i]<-mean(dat$DENOM)/2
  }
  
  }
  
  
 
  variant_use<-data.frame(name,number)
  return(list(prs_reading, variant_use))
}

prs<-c("PGS000011","PGS000013","PGS000018","PGS001780","PGS003446","custom_PGS000337","custom_Oni-Orisan")

dataset_combine<-readin_prs(prs)
prs_reading<-dataset_combine[[1]]
variant_use<-dataset_combine[[2]]

# phenotype data readin -------------------------------



phenotype<-read.table(paste(data_path, '/extracted-phenotypes.txt',sep=''), 
                      sep ="", header=T,dec =".",fill=T)

baseline<-data.frame(fread(paste(phenotype_path,"v2.0_BASELINE.csv",sep = "")))
base_keep<-baseline%>%select(REGISTRO,ACTIVITY)
reg_link<-data.frame(fread(paste(phenotype_path,"RGN_LINK_IID.csv",sep = "")))

base_keep<-merge(base_keep,reg_link[,c(1,3)],by="REGISTRO")
base_keep<-base_keep[,c(2,3)]
phenotype<-left_join(phenotype,base_keep,by="IID")




saveRDS(phenotype,paste(data_path,'/phenotype_extracted_07May2023.rds',sep=''))

nrow(prs_reading)
nrow(phenotype)

#merge datasets--------------------------------------------


full_data_with_prs<-merge(phenotype,prs_reading,by='IID')
saveRDS(full_data_with_prs,paste(data_path,'/FullData_07May2023.rds',sep=''))

# SUMMARY stats -------------

ggplot(data=full_data_with_prs)+geom_density(aes(x=PGS000011,group=factor(BASE_CHD),fill=as.factor(BASE_CHD),
),alpha=0.5)

ggplot(data=full_data_with_prs)+geom_density(aes(x=PGS000013,group=factor(BASE_CHD),fill=as.factor(BASE_CHD),
),alpha=0.5)

ggplot(data=full_data_with_prs)+geom_density(aes(x=PGS000018,group=factor(BASE_CHD),fill=as.factor(BASE_CHD),
),alpha=0.5)

ggplot(data=full_data_with_prs)+geom_density(aes(x=PGS001780,group=factor(BASE_CHD),fill=as.factor(BASE_CHD),
),alpha=0.5)
ggplot(data=full_data_with_prs)+geom_density(aes(x=PGS003446,group=factor(BASE_CHD),fill=as.factor(BASE_CHD),
),alpha=0.5)
ggplot(data=full_data_with_prs)+geom_density(aes(x=custom_PGS000337,group=factor(BASE_CHD),fill=as.factor(BASE_CHD),
),alpha=0.5)
ggplot(data=full_data_with_prs)+geom_density(aes(x=`custom_Oni-Orisan`,group=factor(BASE_CHD),fill=as.factor(BASE_CHD),
),alpha=0.5)


##phenotype summary--------------------
apply(full_data_with_prs%>%select(SEX, SMOKE,BASE_CHD,EPA001,EPA001A,EPO001,EPO001A),2,table)
apply(full_data_with_prs%>%select(BMI,SBP,DBP,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PGS000011,PGS000013,PGS000018,PGS000018,
                                  PGS001780,PGS003446,custom_PGS000337,`custom_Oni-Orisan`),2,summary)

