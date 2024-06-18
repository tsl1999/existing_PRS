##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# Supplementary part 1
# EPO+Baseline CAD
# Author: Tianshu Liu
# Date: 03 June 2024
#################################
#call packages and data

rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
data_path<-paste(working_path,'/data',sep='')
data_save_path<-paste(working_path,'/data/manuscript-age80/supplement',sep='')
graphs_path<-paste(working_path,'/graphs/manuscript-age80/supplement',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  
#data---------------------
data<-readRDS(paste(data_path,"analysis_data_80_17Apr2024.rds",sep="/"))
remove_data<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data<-data[,-remove_data]
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1

#EPO+Baseline, logistic regression---------------
partial_adjustments=c("AGE","SEX",paste(rep("PC",7),seq(1,7),sep=""))
full_adjustments = c(partial_adjustments,"WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline")

model_output_partial<-create_output_table(
  trainsplit = F,data,adjustments=partial_adjustments,outcome="prevalent_CHD_EPO",namew="partial",roc=F)

model_output_full<-create_output_table(
  trainsplit = F,data,adjustments=full_adjustments,outcome="prevalent_CHD_EPO",namew="full",roc=F)


write.csv(model_output_partial,paste(data_save_path,"/baselineEPO_partial.csv",sep=""))
saveRDS(model_output_partial,paste(data_save_path,"/baselineEPO_partial.rds",sep=""))
write.csv(model_output_full,paste(data_save_path,"/baselineEPO_full.csv",sep=""))
saveRDS(model_output_full,paste(data_save_path,"/baselineEPO_full.rds",sep=""))
## results table--------------------------------------------------------------
# prevalent_EPO_cont<-combine_model_tables(
#   tables=list(model_output_partial,model_output_full),
#   name = c("Partially adjusted","Fully adjusted"),log_cox = "log")
# write.csv(prevalent_EPO_cont,paste(data_save_path,"EPO_baseline_flex.csv",sep=""))
# saveRDS(prevalent_EPO_cont,paste(data_save_path,"EPO_baseline_flex.rds",sep=""))
# 
# flex_table<-prev_cont_flextable(prevalent_EPO_cont,table_caption = "Association of PRSs with CAD 
#                                 before age 80 (baseline CAD and CAD mortality as primary cause of death)")
# saveRDS(flex_table,paste(data_save_path,"/supplement/EPO_baseline_results_table.rds",sep=""))

##forest table-------------------------------------------------------------

forest_table_EPO<-generate_forest_table(
  tables=list(model_output_partial,model_output_full),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table_EPO,xlim=c(0.9,1.4),ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or="OR",
               footnote_in=" ")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/EPO_baseline_forest.png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

## quintile plots----------------------------------------------------------
data<-readRDS(paste(data_path,"analysis_data_80_17Apr2024.rds",sep="/"))
remove_data<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data_up<-data[,c(1:32,remove_data)]

model_output_partial<-create_output_table(
  trainsplit = F,data_up,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPO",namew="partial",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_up,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPO",namew="full",type = "cat",roc = F)


## calculate floating absolute risk and plot -----------
model_output_list<-list(model_output_partial,model_output_full)

plot_FAR_gg(model_output_list,name=c("Partial","Full"),outcome="EPO_baseline",data = data,
            graphs_path=graphs_path)
