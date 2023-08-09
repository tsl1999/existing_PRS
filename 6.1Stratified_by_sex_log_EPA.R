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
data<-readRDS(paste(data_path,"analysis_data_13Jun2023.rds",sep="/"))
data_decile<-readRDS(paste(data_path,"/analysis_data_31May2023.rds",sep=""))
data_d<-data_decile[,c(1:32,40:46)]
model_table_top<-readRDS(paste(data_path,"/model_table_top.rds",sep=""))

data_female<-data[data$SEX=="Women",]
data_male<-data[data$SEX=="Men",]
data_female_d<-data_d[data_d$SEX=="Women",]
data_female_up<-data_decile[data_decile$SEX=="Women",]
data_male_d<-data_d[data_d$SEX=="Men",]
data_male_up<-data_decile[data_decile$SEX=="Men",]


# checks-----------------------------------
table(data$prevalent_CHD_EPA)
table(data$SEX)
prop.table(table(data_female$prevalent_CHD_EPA))
prop.table(table(data_male$prevalent_CHD_EPA))
prop.table(table(data_female_d$prevalent_CHD_EPA))
prop.table(table(data_male_d$prevalent_CHD_EPA))


# set adjustments -----------------------
non_adjustment=NULL
partial_adjustments=c("AGE")
full_adjustments = c("AGE","WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline")

# Female -----------------------------

## without prs
model_withoutPRS_EPA_female<-discrimination_without_prs(train_data = data_female,test_data = data_female,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",partial_adjustment = partial_adjustments)
print(model_withoutPRS_EPA_female)
#saveRDS(model_withoutPRS_EPA_female,paste(data_path,"/Prevalent_CHD_EPA_withoutPRSF.rds",sep=""))

## analysis --------------------------
model_output_partialf<-create_output_table(
  trainsplit = F,data_female,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPA",namew="partialF",roc = F)

model_output_fullf<-create_output_table(
  trainsplit = F,data_female,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPA",namew="fullF",roc = F)

## results table ------------------------
prevalent_EPA_contf<-combine_model_tables(
  tables=list(model_output_partialf,model_output_fullf),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")

#saveRDS(prevalent_EPA_contf,paste(data_path,"/Prevalent_CHD_EPA_contF_logistic_13Jun2023.rds",sep=""))

flex_table<-prev_cont_flextable(prevalent_EPA_contf,table_caption = "Association of PRSs with CAD (baseline CAD and CAD mortality anywhere mention on death certificate) in women")
saveRDS(flex_table,paste(data_path,"/Prevalent_CHD_EPA_contF_logistic_flex.rds",sep=""))

## forest table -----------------------------
forest_table<-generate_forest_table(
  tables=list(model_output_partialf,model_output_fullf),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.3),ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or = "OR")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/logistic/Logistic_forestplotEPA_contf",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

## quintile plots-------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data_female_d,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPA",namew="partialF",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_female_d,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPA",namew="fullF",type = "cat",roc = F)


model_output_list<-list(model_output_partial,model_output_full)

plot_FAR_gg(model_output_list,name=c("Partial","Full"),
            outcome="Prevalent_CHD_EPA_F",data=data_female_up)


# Male -----------------------------------------

## without prs
model_withoutPRS_EPA_male<-discrimination_without_prs(train_data = data_male,test_data = data_male,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA",partial_adjustment = partial_adjustments)
print(model_withoutPRS_EPA_male)
#saveRDS(model_withoutPRS_EPA_male,paste(data_path,"/Prevalent_CHD_EPA_withoutPRSM.rds",sep=""))

## analysis --------------------------
model_output_partialm<-create_output_table(
  trainsplit = F,data_male,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPA",namew="partialM",roc = F)

model_output_fullm<-create_output_table(
  trainsplit = F,data_male,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPA",namew="fullM",roc = F)

## results table ------------------------
prevalent_EPA_contm<-combine_model_tables(
  tables=list(model_output_partialm,model_output_fullm),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")

#saveRDS(prevalent_EPA_contm,paste(data_path,"/Prevalent_CHD_EPA_contM_logistic_13Jun2023.rds",sep=""))

flex_table<-prev_cont_flextable(prevalent_EPA_contm,table_caption  = "Association of PRSs with CAD (baseline CAD and CAD mortality anywhere mention on death certificate) in men")
saveRDS(flex_table,paste(data_path,"/Prevalent_CHD_EPA_contM_logistic_flex.rds",sep=""))

## forest table -----------------------------
forest_table<-generate_forest_table(
  tables=list(model_output_partialm,model_output_fullm),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.3),ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or = "OR")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/logistic/Logistic_forestplotEPA_contm",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

## quintile plots-------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data_male_d,adjustments=partial_adjustments,
  outcome="prevalent_CHD_EPA",namew="partialF",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_male_d,adjustments=full_adjustments,
  outcome="prevalent_CHD_EPA",namew="fullF",type = "cat",roc = F)


model_output_list<-list(model_output_partial,model_output_full)

plot_FAR_gg(model_output_list,name=c("Partial","Full"),
            outcome="Prevalent_CHD_EPA_M",data=data_male_up)


