#this file is to prepare dataset for analysis
rm(list=ls())

#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
phenotype_path<-"/well/emberson/projects/mcps/data/phenotypes/"
setwd(working_path)

library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pROC)
library(ROCR)
library(caret)


#readin data -----------------------------------------------
data<-readRDS(paste(data_path,"FullData_07Aug2023.rds",sep="/"))

prs<-data%>%select(IID,contains(c("PGS","custom")))
ncol_prs<-ncol(prs)

for (i in 2:ncol_prs){
  prs[,ncol_prs-1+i]<-(prs[,i]-mean(prs[,i]))/sd(prs[,i])
  colnames(prs)[ncol_prs-1+i]<-paste(colnames(prs)[i],"standardised",sep="_")
}


apply(prs[2:ncol(prs)],2,summary)

prs_standardised<-readRDS(paste(data_path,"FullData_standardisedPRS_31May2023.rds",sep="/"))
prs_standardised_new<-left_join(prs_standardised[,1:34],prs[,c(1,(ncol_prs+1):ncol(prs))])
data_table_1_standardised<-readRDS(paste(data_path,"Table1_data_standardisedPRS_31May2023.rds",sep="/"))
data_table_1_standardised_new<-left_join(data_table_1_standardised[,1:34],prs[,c(1,(ncol_prs+1):ncol(prs))])
prs_standardised_75<-readRDS(paste(data_path,"FullData_standardisedPRS_75_16Jun2023.rds",sep="/"))
prs_standardised_75_new<-left_join(prs_standardised_75[,1:34],prs[,c(1,(ncol_prs+1):ncol(prs))])

saveRDS(prs_standardised_new,paste(data_path,"FullData_standardisedPRS_07Aug2023.rds",sep="/"))
saveRDS(data_table_1_standardised_new,paste(data_path,"Table1_data_standardisedPRS_07Aug2023.rds",sep="/"))
saveRDS(prs_standardised_75_new,paste(data_path,"FullData_standardisedPRS_75_07Aug2023.rds",sep="/"))




# plot density plots --------------------------------------



dt_melt_standardised<-melt(prs[,c(1,10:ncol(prs))],"IID")
ggplot(data = dt_melt_standardised, aes(x=value)) + geom_density( alpha = 0.4,fill="grey")+facet_wrap( ~ variable)



data_ancestry<-data.frame(fread(paste(data_path,"/ancestry_gen.csv",sep="")))
#data<-readRDS(paste(data_path,"/analysis_data_31May2023.rds",sep="/"))

#data_restricted_IID<-data.frame(fread("/well/emberson/users/bjk420/projects/popgen_v2/07_pca/mcps_only/no_references/king/ivs_unrelated-set_oxford.txt", header=F))

#colnames(data_ancestry)[1]<-"IID"
#data_ancestry_join<-left_join(data,data_ancestry,by="IID")

#data_restricted<-data[data$IID%in%data_restricted_IID$V2,]
#saveRDS(data_ancestry_join,paste(data_path,"/FullData_ancestry_22Jun2023.rds",sep=""))
#saveRDS(data_restricted,paste(data_path,"/RestrictedData_22Jun2023.rds",sep=""))
