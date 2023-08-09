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

# data readin-------------------------------------------
data<-readRDS(paste(data_path,"FullData_standardisedPRS_75_07Aug2023.rds",sep="/"))
data<-data[,c(1:32,35:ncol(data))]
colnames(data)[33:ncol(data)]<-sub("_standardised","",colnames(data)[33:ncol(data)])#remove standardised after all PRS names to make column names shorter
colnames(data)[40]<-"custom_Oni_Orisan"
model_table_top<-readRDS(paste(data_path,"/model_table_top_07Aug2023.rds",sep=""))
npgs<-ncol(model_table_top)-1
# checks---------------------------------------------
table(data$prevalent_CHD_EPA)
table(data$prevalent_CHD_EPO)
table(data$BASE_CHD)


# set adjustments-----------------------------------
non_adjustment=NULL
partial_adjustments=c("AGE","SEX")
full_adjustments = c("AGE","SEX","WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline")

## EPA--------------------------------------
outcome="prevalent_CHD_EPA"

### without prs------------------------------
model_without_75<-discrimination_without_prs(
  train_data = data,test_data = data,outcome=outcome,
  partial_adjustments = partial_adjustments,full_adjustments = full_adjustments)
print(model_without_75)
#saveRDS(model_without_75,paste(data_path,"/Prevalent_CHD_EPA_withoutPRS75.rds",sep=""))

### analysis ----------------------------------------
model_output_simple<-create_output_table(
  trainsplit = F,data,adjustments=non_adjustment,outcome=outcome,namew="simple_75",roc = F)


model_output_partial<-create_output_table(
  trainsplit = F,data,adjustments=partial_adjustments,outcome=outcome,namew="partial_75",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data,adjustments=full_adjustments,outcome=outcome,namew="full_75",roc = F)


### results table-----------------------------------
prevalent_base_cont_75<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")
#saveRDS(prevalent_base_cont_75,paste(data_path,"/Prevalent_CHD_EPA_cont75_logistic_15Jun2023.rds",sep=""))


flex_table<-prev_cont_flextable(prevalent_base_cont_75,table_caption = "Association of PRSs with CAD (baseline CAD and CAD mortality anywhere mentioned on the death certificate) among people aged over 75 at baseline")
saveRDS(flex_table,paste(data_path,"/Prevalent_CHD_EPA_cont75_logistic_flex_07Aug2023.rds",sep=""))

### forest table----------------------------------
forest_table<-generate_forest_table(tables=list(model_output_partial,model_output_full),model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.2),ticks_at = c(0.9,1,1.1,1.2),hr_or = "OR")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/logistic/Logistic_forestplotCHD_EPA_cont75_07Aug2023",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

### quintile table--------------------------------------

# quintile cut
ncol_data<-ncol(data)
for (i in 1:npgs){
  quintile_cutoff<-quantile(data[,32+i],probs=seq(0,1,0.2),na.rm=F)
  
  data[,ncol_data+i]<-cut(data[,32+i], breaks=quintile_cutoff,labels = c(1:5),include.lowest=TRUE)
  colnames(data)[ncol_data+i]<-paste(colnames(data)[32+i],"cat",sep="_")
}
data_up<-data[,c(1:32,(ncol_data+1):ncol(data))]


model_output_simple<-create_output_table(
  trainsplit = F,data_up,adjustments=non_adjustment,outcome="prevalent_CHD_EPA",
  namew=NULL,type = "cat",roc = F)

model_output_partial<-create_output_table(
  trainsplit = F,data_up,adjustments=partial_adjustments,outcome="prevalent_CHD_EPA",
  namew="partial",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_up,adjustments=full_adjustments,outcome="prevalent_CHD_EPA",
  namew="full",type = "cat",roc = F)


model_output_list<-list(model_output_partial,model_output_full)
#saveRDS(model_output_list,paste(data_path,"/prevalent_CHD_EPA_cat_logistic_06Jun2023.rds",sep=""))

plot_FAR_gg(model_output_list,name=c("Partial","Full"),outcome="Prevalent_CHD_EPA75_07Aug2023",data = data)


## EPO-------------------------------------
data_or<-data[,1:40]
outcome="prevalent_CHD_EPO"

### without prs----------------------------
model_without_75<-discrimination_without_prs(
  train_data = data_or,test_data = data_or,outcome=outcome,
  partial_adjustments = partial_adjustments,full_adjustments = full_adjustments)
print(model_without_75)
#saveRDS(model_without_75,paste(data_path,"/Prevalent_CHD_EPO_withoutPRS75.rds",sep=""))

### analysis ----------------------------------------
model_output_simple<-create_output_table(
  trainsplit = F,data_or,adjustments=non_adjustment,outcome=outcome,namew="simple_75",roc = F)


model_output_partial<-create_output_table(
  trainsplit = F,data_or,adjustments=partial_adjustments,outcome=outcome,namew="partial_75",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_or,adjustments=full_adjustments,outcome=outcome,namew="full_75",roc = F)


### results table-----------------------------------
prevalent_base_cont_75<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")
#saveRDS(prevalent_base_cont_75,paste(data_path,"/Prevalent_CHD_EPO_cont75_logistic_15Jun2023.rds",sep=""))


flex_table<-prev_cont_flextable(prevalent_base_cont_75,table_caption = "Association of PRSs with CAD (baseline CAD and CAD mortality as primary cause of death) among people aged over 75 at baseline")
saveRDS(flex_table,paste(data_path,"/Prevalent_CHD_EPO_cont75_logistic_flex_07Aug2023.rds",sep=""))

### forest table----------------------------------
forest_table<-generate_forest_table(tables=list(model_output_partial,model_output_full),model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.2),ticks_at = c(0.9,1,1.1,1.2),hr_or = "OR")
a<-get_wh(p,unit = "cm")

png(paste(graphs_path,"/logistic/Logistic_forestplotCHD_EPO_cont75_07Aug2023",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

### quintile table-------------------------------
model_output_simple<-create_output_table(
  trainsplit = F,data_up,adjustments=non_adjustment,outcome="prevalent_CHD_EPO",
  namew=NULL,type = "cat",roc = F)

model_output_partial<-create_output_table(
  trainsplit = F,data_up,adjustments=partial_adjustments,outcome="prevalent_CHD_EPO",
  namew="partial",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_up,adjustments=full_adjustments,outcome="prevalent_CHD_EPO",
  namew="full",type = "cat",roc = F)

model_output_list<-list(model_output_partial,model_output_full)
#saveRDS(model_output_list,paste(data_path,"/prevalent_CHD_EPO_cat_logistic_06Jun2023.rds",sep=""))

plot_FAR_gg(model_output_list,name=c("Partial","Full"),outcome="Prevalent_CHD_EPO75_07Aug2023",data = data)
