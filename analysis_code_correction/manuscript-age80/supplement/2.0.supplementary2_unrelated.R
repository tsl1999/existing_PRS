##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# Supplementary part 3
# 3rd degree unrelated participants
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

#data readin-----------------------------------------------------------------
data<-readRDS(paste(data_path,"/RestrictedData_80_17Apr2024.rds",sep=""))
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1
remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data_cont<-data[,-remove_cat]
data_cat<-data[,c(1:32,remove_cat)]

#analysis EPA-------------------------------------------------------------------

## set adjustments -----------------------
partial_adjustments=c("AGE","SEX",paste(rep("PC",7),seq(1,7),sep=""))
full_adjustments = c(partial_adjustments,"WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline")

## analysis --------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data_cont,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPA",namew="partial3rd",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_cont,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPA",namew="full3rd",roc = F)
write.csv(model_output_partial,paste(data_save_path,"/thirdEPA_partial.csv",sep=""))
saveRDS(model_output_partial,paste(data_save_path,"/thirdEPA_partial.rds",sep=""))
write.csv(model_output_full,paste(data_save_path,"/thirdEPA_full.csv",sep=""))
saveRDS(model_output_full,paste(data_save_path,"/thirdEPA_full.rds",sep=""))

## results table ----------------------------------------------------------
prevalent_EPA_cont<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")

saveRDS(prevalent_EPA_cont,paste(data_save_path,"/primary_unrelated.rds",sep=""))
write.csv(prevalent_EPA_cont,paste(data_save_path,"/primary_unrelated.csv",sep=""))

flex_table<-prev_cont_flextable(prevalent_EPA_cont,table_caption = "Association of PRSs with CAD (baseline CAD and CAD mortality anywhere mentioned on the death certificate) for 3rd degree unrelated participants")
saveRDS(flex_table,paste(data_save_path,"/primary_unrelated_results_table.rds",sep=""))


##forest table-----------------------------------------------------------
forest_table<-generate_forest_table(
  tables=list(model_output_partial,model_output_full),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.45),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),hr_or = "OR",footnote_in=" ")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/primary_unrelated_forest",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

## quintile plots-------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data_cat,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPA",namew="partial3rd",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_cat,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPA",namew="full3rd",type = "cat",roc = F)


model_output_list<-list(model_output_partial,model_output_full)

plot_FAR_gg(model_output_list,name=c("Partial","Full"),
            outcome="unrelated",data=data,graphs_path = graphs_path)

