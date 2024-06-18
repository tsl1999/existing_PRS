##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# Figure 2
# Quintile plots for primary analysis
# Author: Tianshu Liu
# Date: 03 June 2024
#################################
#call packages and data
##package and self-create function

rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
data_path<-paste(working_path,'/data',sep='')
data_save_path<-paste(working_path,'/data/manuscript-age80',sep='')
graphs_path<-paste(working_path,'/graphs/manuscript-age80',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  

#flextable of CAD-EPA--------------------------------------------------------
model_output_partial<-readRDS(paste(data_save_path,"/primary_partial.rds",sep=""))
model_output_full<-readRDS(paste(data_save_path,"/primary_full.rds",sep=""))
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1
prevalent_EPA_cont<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")
write.csv(prevalent_EPA_cont,paste(data_save_path,"/supplement/primary_flex.csv",sep=""))
saveRDS(prevalent_EPA_cont,paste(data_save_path,"/supplement/primary_flex.rds",sep=""))
flex_table<-prev_cont_flextable(prevalent_EPA_cont,
                                table_caption = "Association of PRSs with CAD 
                                before age 80(baseline CAD and CAD mortality 
                                anywhere mentioned on the death certificate)")

saveRDS(flex_table,paste(data_save_path,"/supplement/primary_results_table.rds",sep=""))


#quintile plot--------------------------------------------------------------
data<-readRDS(paste(data_path,"analysis_data_80_17Apr2024.rds",sep="/"))
remove_data<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data_up<-data[,c(1:32,remove_data)]

## EPA ----------------------
non_adjustment=NULL
partial_adjustments=c("AGE","SEX",paste(rep("PC",7),seq(1,7),sep=""))
full_adjustments = c(partial_adjustments,"WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline")


## analysis ------------------
model_output_partial<-create_output_table(
  trainsplit = F,data_up,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPA",namew="partial",type = "cat",roc = F)

saveRDS(model_output_partial,paste(data_save_path,"/primary_quintile_partial.rds",sep=""))

model_output_full<-create_output_table(
  trainsplit = F,data_up,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPA",namew="full",type = "cat",roc = F)

saveRDS(model_output_full,paste(data_save_path,"/primary_quintile_full.rds",sep=""))

## calculate floating absolute risk and plot -----------
model_output_list<-list(model_output_partial,model_output_full)
plot_FAR_gg(model_output_list,name=c("Partial","Full"),outcome="Figure2_primary",data = data,
            graphs_path=graphs_path)




