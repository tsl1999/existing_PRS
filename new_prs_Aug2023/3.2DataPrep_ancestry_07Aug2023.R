#this file is to prepare dataset for analysis
rm(list=ls())

#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
phenotype_path<-"/well/emberson/projects/mcps/data/phenotypes/"
setwd(working_path)

library(data.table)

#readin data -----------------------------------------------
data_ancestry<-data.frame(fread(paste(data_path,"/ancestry_gen.csv",sep="")))
data<-readRDS(paste(data_path,"/analysis_data_07Aug2023.rds",sep="/"))

data_restricted_IID<-data.frame(fread("/well/emberson/users/bjk420/projects/popgen_v2/07_pca/mcps_only/no_references/king/ivs_unrelated-set_oxford.txt", header=F))

colnames(data_ancestry)[1]<-"IID"
data_ancestry_join<-left_join(data,data_ancestry,by="IID")

data_restricted<-data[data$IID%in%data_restricted_IID$V2,]
saveRDS(data_ancestry_join,paste(data_path,"/FullData_ancestry_07Aug2023.rds",sep=""))
saveRDS(data_restricted,paste(data_path,"/RestrictedData_07Aug2023.rds",sep=""))
