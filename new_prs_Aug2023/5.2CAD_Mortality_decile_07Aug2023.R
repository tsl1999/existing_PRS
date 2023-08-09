#call packages and data
##package and self-create function
rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
setwd(working_path)
library(survival)
library(survminer)
library(dynpred)
library(Epi)
library(forestploter)
library(ggfortify)
source("/gpfs3/well/emberson/users/hma817/codes/R_TablesFunction_27Sep2022_TL.r")
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  

# data readin ------------------

data<-readRDS(paste(data_path,"FullData_timesplit_07Aug2023.rds",sep="/"))
data2<-readRDS(paste(data_path,"FullData_standardisedPRS_mortality_07Aug2023.rds",sep="/"))

model_table_top<-readRDS(paste(data_path,"/model_table_top_07Aug2023.rds",sep=""))
npgs<-ncol(model_table_top)-1
data$agegroup<-factor(data$agegroup,levels = c(1:5))

# checks----------------------
model_simple<-coxph(Surv(time_in,time_out,EPO001_up_75,type = "counting")~PGS000011_cat+SEX,data=data)
fl_a<-float(model_simple,factor="PGS000011_cat")

# EPO001------------------
## set adjustments --------------------
remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data_up<-data[,c(1:32,41:45,remove_cat[1]:ncol(data))]
outcome="EPO001_up_75"
non_adjustment=NULL
partial_adjustments=c("SEX")
strata="agegroup"
full_adjustments = c(partial_adjustments,c("WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline"))

## analysis ------------------
model_output_simple<-create_output_table_cox(
  trainsplit = F,data_up,adjustments=non_adjustment,outcome=outcome,namew=NULL,type = "cat",strata = strata,ph=F)

model_output_partial<-create_output_table_cox(
  trainsplit = F,data_up,adjustments=partial_adjustments,outcome=outcome,namew="partial",type = "cat",strata = strata,ph=F)

model_output_full<-create_output_table_cox(
  trainsplit = F,data_up,adjustments=full_adjustments,outcome=outcome,namew="full",type = "cat",strata = strata,ph=F)

## FAR plot ---------------------------
model_output_list<-list(model_output_partial,model_output_full)
#saveRDS(model_output_list,paste(data_path,"/EPO001_cat_cox_07Jun2023.rds",sep=""))
plot_FAR_gg(model_output_list,name=c("Partial","Full"),outcome="EPO001_07Aug2023",type="HR",data = data)

## results table ---------------------
model_output_partial_table<-model_output_partial[[1]]
model_output_full_table<-model_output_full[[1]]
EPO_cat<-combine_model_tables(tables=list(model_output_partial_table,model_output_full_table),name = c("Partially adjusted","Fully adjusted"),log_cox = "cox")
#saveRDS(EPO_cat,paste(data_path,"/EPO001_cat_AUC_table_06Jun2023.rds",sep=""))

flex_cat<-my_table(EPO_cat,y=c(1:3,5),x=1:(npgs+1))
flex_cat<-footnote(flex_cat,i=c(3,5),j=1,part="body",value = as_paragraph(c("Partially adjusted models adjust for Sex and baseline age. Fully adjusted models additionally adjust for waist-to-hip ratio, systolic and diastolic blood pressure, education attainmen level, smoking status,diabetes at baseline")),ref_symbols = c("*"))
flex_cat <- add_header_row(flex_cat, values =  c(" ", "European PRS","Non-European PRS"),
                           colwidths = c(1, 4,npgs-4), top = T)
flex_cat<-set_caption(
  flex_cat,
  caption = "Association of PRSs in quintile groups with CAD mortality(as primary cause of death)",align_with_table = FALSE,
  word_stylename = "Table Caption",
  fp_p = fpp
)

flex_cat<-my_theme(flex_cat,y=c(1,2,4),set_padding = 3,fontsize = 8,header_border = 7)
saveRDS(flex_cat,paste(data_path,"/EPO001_cat_AUC_table_flex_07Aug2023.rds",sep=""))


# EPA001------------------

outcome="EPA001_up_75"


## analysis ------------------
model_output_simple<-create_output_table_cox(
  trainsplit = F,data_up,adjustments=non_adjustment,outcome=outcome,namew=NULL,type = "cat",strata = strata,ph=F)

model_output_partial<-create_output_table_cox(
  trainsplit = F,data_up,adjustments=partial_adjustments,outcome=outcome,namew="partial",type = "cat",strata = strata,ph=F)

model_output_full<-create_output_table_cox(
  trainsplit = F,data_up,adjustments=full_adjustments,outcome=outcome,namew="full",type = "cat",strata = strata,ph=F)

## FAR plot ---------------------------
model_output_list<-list(model_output_partial,model_output_full)
#saveRDS(model_output_list,paste(data_path,"/EPO001_cat_cox_07Jun2023.rds",sep=""))
plot_FAR_gg(model_output_list,name=c("Partial","Full"),outcome="EPA001_07Aug2023",type="HR",data = data)

## results table ---------------------
model_output_partial_table<-model_output_partial[[1]]
model_output_full_table<-model_output_full[[1]]
EPA_cat<-combine_model_tables(
  tables=list(model_output_partial_table,model_output_full_table),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "cox")
#saveRDS(EPA_cat,paste(data_path,"/EPA001_cat_AUC_table_06Jun2023.rds",sep=""))

flex_cat<-my_table(EPA_cat,y=c(1:3,5),x=1:(npgs+1))
flex_cat<-footnote(flex_cat,i=c(3,5),j=1,part="body",value = as_paragraph(c("Partially adjusted models adjust for Sex and baseline age. Fully adjusted models additionally adjust for waist-to-hip ratio, systolic and diastolic blood pressure, education attainmen level, smoking status,diabetes at baseline")),ref_symbols = c("*"))
flex_cat <- add_header_row(flex_cat, values =  c(" ", "European PRS","Non-European PRS"),
                           colwidths = c(1, 4,npgs-4), top = T)
flex_cat<-set_caption(
  flex_cat,
  caption = "Association of PRSs in quintile groups with CAD mortality(anywhere mention on the death certificate)",align_with_table = FALSE,
  word_stylename = "Table Caption",
  fp_p = fpp
)

flex_cat<-my_theme(flex_cat,y=c(1,2,4),set_padding = 3,fontsize = 8,header_border = 7)
saveRDS(flex_cat,paste(data_path,"/EPA001_cat_AUC_table_flex_07Aug2023.rds",sep=""))


