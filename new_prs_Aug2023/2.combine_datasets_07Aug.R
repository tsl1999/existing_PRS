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
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")
#prs data readin-------------------------------
name<-NA
number<-NA


prs_reading<-data.frame(NA)

prs<-c("PGS000011","PGS000013","PGS000018","PGS001780","PGS003446","custom_PGS000337","custom_Oni-Orisan","PGS003725")

dataset_combine<-readin_prs(prs)
prs_reading<-dataset_combine[[1]]
variant_use<-dataset_combine[[2]]
saveRDS(variant_use,paste(data_path,"variants_included_07Aug2023.rds",sep="/"))

# phenotype data readin -------------------------------
phenotype<-readRDS(paste(data_path,'/phenotype_extracted_07May2023.rds',sep=''))
nrow(prs_reading)
nrow(phenotype)

#merge datasets--------------------------------------------


full_data_with_prs<-merge(phenotype,prs_reading,by='IID')
saveRDS(full_data_with_prs,paste(data_path,'/FullData_07Aug2023.rds',sep=''))

