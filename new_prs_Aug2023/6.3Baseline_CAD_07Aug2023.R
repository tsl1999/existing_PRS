#call packages and data
##package and self-create function
rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
setwd(working_path)

library(dynpred)
library(Epi)
library(forestploter)
source("/gpfs3/well/emberson/users/hma817/codes/R_TablesFunction_27Sep2022_TL.r")
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  

# data readin ---------------------------------

data<-readRDS(paste(data_path,"FullData_standardisedPRS_07Aug2023.rds",sep="/"))
data<-data[,c(1:32,35:ncol(data))]
table1_data<-readRDS(paste(data_path,"Table1_data_standardisedPRS_07Aug2023.rds",sep="/"))
colnames(data)[33:ncol(data)]<-sub("_standardised","",colnames(data)[33:ncol(data)])#remove standardised after all PRS names to make column names shorter
colnames(data)[40]<-"custom_Oni_Orisan"
model_table_top<-readRDS(paste(data_path,"/model_table_top_07Aug2023.rds",sep=""))


# checks 

table(data$BASE_CHD)

#set adjustments---------------------------
non_adjustment=NULL
partial_adjustments=c("AGE","SEX")
full_adjustments = c("AGE","SEX","WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline")

# without prs-------------------------------
model_withoutPRS_baseCHD<-discrimination_without_prs(
  train_data = data,test_data = data,outcome="BASE_CHD",
  partial_adjustments = partial_adjustments,full_adjustments = full_adjustments)
print(model_withoutPRS_baseCHD)
#saveRDS(model_withoutPRS_baseCHD,paste(data_path,"/BASE_CHD_withoutPRS.rds",sep=""))


# analysis----------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data,adjustments=partial_adjustments,
  outcome="BASE_CHD",namew="partial",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data,adjustments=full_adjustments,
  outcome="BASE_CHD",namew="full",roc = F)

## results table---------------------
prevalent_base_cont<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")
#saveRDS(prevalent_base_cont,paste(data_path,"/Prevalent_CHD_EPA_cont_logistic_15Jun2023.rds",sep=""))

flex_table<-prev_cont_flextable(prevalent_base_cont,table_caption= "Association of PRSs with baseline CAD")
saveRDS(flex_table,paste(data_path,"/Prevalent_CHD_baseCAD_logistic_flex_07Aug2023.rds",sep=""))

## forest table-------------------------
forest_table<-generate_forest_table(tables=list(model_output_partial,model_output_full),model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.45),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),hr_or = "OR")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/logistic/Logistic_forestplotbaseCHD_cont_07Aug2023",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

## quintile table-------------------------
data<-readRDS(paste(data_path,"analysis_data_07Aug2023.rds",sep="/"))
remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))

data_up<-data[,c(1:32,remove_cat)]

model_output_partial<-create_output_table(
  trainsplit = F,data_up,adjustments=partial_adjustments,outcome="BASE_CHD",
  namew="partial",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_up,adjustments=full_adjustments,outcome="BASE_CHD",
  namew="full",type = "cat",roc = F)

model_output_list<-list(model_output_partial,model_output_full)
#saveRDS(model_output_list,paste(data_path,"/prevalent_CHD_EPA_cat_logistic_06Jun2023.rds",sep=""))

plot_FAR_gg(model_output_list,name=c("Partial","Full"),outcome="BASE_CHD_07Aug2023",data=data)












