#call packages and data
##package and self-create function
rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/codes/R_TablesFunction_27Sep2022_TL.r")
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  
library(Epi)

#data readin-----------------------------------------------------------------
data<-readRDS(paste(data_path,"/RestrictedData_07Aug2023.rds",sep=""))
model_table_top<-readRDS(paste(data_path,"/model_table_top_07Aug2023.rds",sep=""))
npgs<-ncol(model_table_top)-1
remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data_cont<-data[,-remove_cat]
data_cat<-data[,c(1:32,remove_cat)]

#analysis EPA-------------------------------------------------------------------


## set adjustments -----------------------
non_adjustment=NULL
partial_adjustments=c("AGE","SEX")
full_adjustments = c("AGE","WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline")
## model without prs--------------------------------------------------------
model_withoutPRS_EPA<-discrimination_without_prs(
  train_data = data_cont,test_data = data_cont,
  full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",
  partial_adjustment = partial_adjustments)
print(model_withoutPRS_EPA)

saveRDS(model_withoutPRS_EPA,paste(data_path,"/Prevalent_CHD_EPA_withoutPRS_restricted3rd_07Aug2023.rds",sep=""))

## with prs--------------------------------------------------------------------


## analysis --------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data_cont,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPA",namew="partial3rd",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_cont,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPA",namew="full3rd",roc = F)

## results table ----------------------------------------------------------
prevalent_EPA_cont<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")

flex_table<-prev_cont_flextable(prevalent_EPA_cont,table_caption = "Association of PRSs with CAD (baseline CAD and CAD mortality anywhere mentioned on the death certificate) for 3rd degree unrelated participants")
saveRDS(flex_table,paste(data_path,"/Prevalent_CHD_EPA_cont_restricted3rd_flex_07Aug2023.rds",sep=""))

##forest table-----------------------------------------------------------
forest_table<-generate_forest_table(
  tables=list(model_output_partial,model_output_full),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.3),ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or = "OR")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/logistic/Logistic_forestplotEPA_cont_restricted3rd_07Aug2023",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
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
            outcome="Prevalent_CHD_EPA_restricted3rd_07Aug2023",data=data)

#analysis EPO-------------------------------------------------------------

model_withoutPRS_EPO<-discrimination_without_prs(
  train_data = data_cont,test_data = data_cont,
  full_adjustments = full_adjustments,outcome="prevalent_CHD_EPO",
  partial_adjustment = partial_adjustments)
print(model_withoutPRS_EPO)

saveRDS(model_withoutPRS_EPO,paste(data_path,"/Prevalent_CHD_EPO_withoutPRS_restricted3rd_07Aug2023.rds",sep=""))

## with prs--------------------------------------------------------------------


## analysis --------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data_cont,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPO",namew="partial3rd",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_cont,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPO",namew="full3rd",roc = F)

## results table ----------------------------------------------------------
prevalent_EPO_cont<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")

flex_table<-prev_cont_flextable(prevalent_EPO_cont,table_caption = "Association of PRSs with CAD (baseline CAD and CAD mortality as primary cause of death) for 3rd degree unrelated participants")
saveRDS(flex_table,paste(data_path,"/Prevalent_CHD_EPO_cont_restricted3rd_flex_07Aug2023.rds",sep=""))

##forest table-----------------------------------------------------------
forest_table<-generate_forest_table(
  tables=list(model_output_partial,model_output_full),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.3),ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or = "OR")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/logistic/Logistic_forestplotEPO_cont_restricted3rd_07Aug2023",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

## quintile plots-------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data_cat,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPO",namew="partial3rd",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_cat,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPO",namew="full3rd",type = "cat",roc = F)


model_output_list<-list(model_output_partial,model_output_full)

plot_FAR_gg(model_output_list,name=c("Partial","Full"),
            outcome="Prevalent_CHD_EPO_restricted3rd_07Aug2023",data=data)
