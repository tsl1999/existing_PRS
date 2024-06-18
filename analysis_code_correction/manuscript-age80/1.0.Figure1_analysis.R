##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# Figure 1
# Author: Tianshu Liu
# Date: 30 May 2024
#################################
#call packages and data
##package and self-create function

rm(list=ls())
library(survival)
library(survminer)
library(dynpred)
library(Epi)
library(forestploter)
library(ggfortify)
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
data_path<-paste(working_path,'/data',sep='')
data_save_path<-paste(working_path,'/data/manuscript-age80',sep='')
graphs_path<-paste(working_path,'/graphs/manuscript-age80',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  
#data---------------------
data<-readRDS(paste(data_path,"analysis_data_80_17Apr2024.rds",sep="/"))
remove_data<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data<-data[,-remove_data]
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1

# For primary analysis EPA+Baseline, logistic regression---------------
partial_adjustments=c("AGE","SEX",paste(rep("PC",7),seq(1,7),sep=""))#paste(rep("PC",7),seq(1,7),sep="")
full_adjustments = c(partial_adjustments,"WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline")

model_output_partial<-create_output_table(
  trainsplit = F,data,adjustments=partial_adjustments,outcome="prevalent_CHD_EPA",namew="partial",roc=F)
write.csv(model_output_partial,paste(data_save_path,"/EPA_mortality_partial.csv",sep=""))
saveRDS(model_output_partial,paste(data_save_path,"/EPA_mortality_partial.rds",sep=""))

model_output_full<-create_output_table(
  trainsplit = F,data,adjustments=full_adjustments,outcome="prevalent_CHD_EPA",namew="full",roc=F)
write.csv(model_output_full,paste(data_save_path,"/EPA_mortality_full.csv",sep=""))
saveRDS(model_output_full,paste(data_save_path,"/EPA_mortality_full.rds",sep=""))


### Figure 1 -------------
forest_table<-generate_forest_table(
  tables=list(model_output_partial,model_output_full),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.4),ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or = "OR",
               footnote_in=" ")
a<-get_wh(p,unit = "cm")+2
png(paste(graphs_path,"/Figure1_logistic_forestEPA_age80",".png",sep=""), res = 300, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()
# 
# pdf(paste(graphs_path,"/Figure1_logistic_forestEPA_age80_30May2024",".pdf",sep=""), width = a[1], height = a[2],paper="a4r")
# print(p)
# dev.off()
