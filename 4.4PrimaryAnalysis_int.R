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


# EPA ---------
# set adjustments first
non_adjustment=NULL
partial_adjustments=c("AGE","SEX")
full_adjustments = c("AGE","SEX","WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline")
int="SEX"
outcome=c("prevalent_CHD_EPA","prevalent_CHD_EPO")

## with prs-------------------

for(i in 1:2){
model_output_partial<-create_output_table(
  trainsplit = F,data,adjustments=partial_adjustments,outcome=outcome[i],namew="partial",roc =F,int=int)

model_output_full<-create_output_table(
  trainsplit = F,data,adjustments=full_adjustments,outcome=outcome[i],namew="full",roc=F,int=int)

### results table -------------
prevalent_EPA_cont_int<-combine_model_tables(
  tables=list(model_output_partial[,1:11],model_output_full[,1:11]),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "log")
#saveRDS(prevalent_EPA_cont,paste(data_path,"/Prevalent_CHD_EPA_cont_logistic_06Jun2023.rds",sep=""))

flex_table<-prev_cont_flextable(prevalent_EPA_cont_int,table_caption = "Association of PRSs with CAD, with PGS and Sex interaction term (baseline CAD and CAD mortality anywhere mentioned on the death certificate)")



### forest table -------------------
forest_table<-generate_forest_table(
  tables=list(model_output_partial,model_output_full),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T,int="Sex")

interaction_result<-forest_table[,c(1,7:10)]
colnames(interaction_result)[2]<-"OR per SD (95%CI)"

flex_table<-my_table(interaction_result,x=1:5,y=seq(1,21,by=3))
flex_table<-add_footer_lines(
  flex_table,value = as_paragraph(c("Partially adjusted models adjust for Sex and baseline age. Fully adjusted models additionally adjust for waist-to-hip ratio, systolic anddiastolic blood pressure, education attainmen level, smoking status,diabetes at baseline")))
flex_table<-set_caption(
  flex_table,
  caption = "Association of PRSs with CAD with an interaction term with Sex",align_with_table = FALSE,
  word_stylename = "Table Caption",
  fp_p = fpp
)

flex_table<-my_theme(flex_table,y=c(21),set_padding = 3,fontsize = 8)

saveRDS(flex_table,paste(data_path,"/",outcome[i],"_cont_int_logistic_flex.rds",sep=""))}

#p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.3),ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or = "OR",col_input = c(1:3,7,11,8:10))
#a<-get_wh(p,unit = "cm")+2

#png(paste(graphs_path,"/logistic/Logistic_forestplotEPA_int_cont",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
#print(p)
#dev.off()
