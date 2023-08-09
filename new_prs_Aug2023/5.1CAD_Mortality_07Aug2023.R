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

# data readin ------------------------
model_table_top<-readRDS(paste(data_path,"/model_table_top_07Aug2023.rds",sep=""))
npgs<-ncol(model_table_top)-1
data<-readRDS(paste(data_path,"FullData_timesplit_07Aug2023.rds",sep="/"))
remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data<-data[,-c(remove_cat)]
data2<-readRDS(paste(data_path,"FullData_standardisedPRS_mortality_07Aug2023.rds",sep="/"))



data$agegroup<-factor(data$agegroup,levels = c(1:5))

#checks --------------------------------


#check events in two data are matched up
table(data$EPA001_up_75)
table(data$EPO001_up_75)
table(data2$EPA001_up&data2$AGE_followup<75)
table(data2$EPO001_up&data2$AGE_followup<75)

#check outputs of two datasets are the same if don't account for time-dependent variables
fit<-coxph(Surv(time_in,time_out,EPO001_up_75,type = "counting")~SEX,data=data)
summary(fit)


fit2<-coxph(Surv(followup,data2$EPO001_up&data2$AGE_followup<75)~SEX,data=data2)
summary(fit2)

#no prs, adjust for age-bands---------------------------
model_no_prs_partial<-coxph(Surv(time_in,time_out,EPO001_up_75,type = "counting")~SEX+as.factor(agegroup),data=data)
c_index_partial<-round(ci_manual(x=concordance(model_no_prs_partial)$concordance,se=sqrt(concordance(model_no_prs_partial)$var)),3)
model_no_prs_full<-coxph(Surv(time_in,time_out,EPO001_up_75,type = "counting")~SEX+as.factor(agegroup)+EDU_LEVEL+as.factor(diabetes_at_baseline)+smokegp+ SBP + 
                           DBP + WHRATIO,data=data)
c_index_full<-round(ci_manual(x=concordance(model_no_prs_full)$concordance,se=sqrt(concordance(model_no_prs_full)$var)),3)

model_withoutPRS_EPO001<-data.frame(partial_cindex=paste(c_index_partial[1],"(",c_index_partial[2],"-",c_index_partial[3],")"),full_cindex=paste(c_index_full[1],"(",c_index_full[2],"-",c_index_full[3],")"))
print(model_withoutPRS_EPO001)
#saveRDS(model_withoutPRS_EPO001,paste(data_path,"/EPO001_withoutPRS.rds",sep=""))

# check if violated schofeld residuals
cox.zph(model_no_prs_partial)
cox.zph(model_no_prs_full)


# EPO001 -----------------
## set adjustments ----------------
outcome="EPO001_up_75"
partial_adjustments=c("SEX")
full_adjustments = c(partial_adjustments,c("WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline"))
strata="agegroup"

## analysis ----------------------
model_output_partial<-create_output_table_cox(
  data=data,adjustments=partial_adjustments,outcome=outcome,trainsplit = F,strata = strata)

model_output_full<-create_output_table_cox(
  data=data,adjustments=full_adjustments,outcome=outcome,trainsplit = F,namew = "_full",strata = strata)


### results table ------------------------
EPO001_cont<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "cox")
#saveRDS(EPO001_cont,paste(data_path,"/EPO001_cont_model_outcomes_05Jun2023.rds",sep=""))
flex_cont<-my_table(EPO001_cont,y=c(1:3,6),x=1:c(npgs+1))
flex_cont<-footnote(flex_cont,i=c(3,6),j=1,part="body",value = as_paragraph(c("Partial adjusted models adjust for Sex and age groups. Fully adjusted models additionally adjust for  waist-to-hip ratio, systolic and diastolic blood pressure, education attainmen level, smoking status,diabetes at baseline")),ref_symbols = c("*"))
flex_cont <- add_header_row(flex_cont, values =  c(" ", "European PRS","Non-European PRS"),
                            colwidths = c(1, 4,npgs-4), top = T)
flex_cont<-set_caption(
  flex_cont,
  caption = "Association of PRSs with CAD mortality(primary cause of death)",align_with_table = FALSE,
  word_stylename = "Table Caption",
  fp_p = fpp
)
flex_cont<-my_theme(flex_cont,y=c(1,2,5,8),set_padding = 3,fontsize = 8,header_border = 7)

saveRDS(flex_cont,paste(data_path,"/EPO001_cont_model_outcomes_flex_07Aug2023.rds",sep=""))

### forest plot ------------------------------
forest_table<-generate_forest_table(tables=list(model_output_partial,model_output_full),model_name = c("Partially adjusted","Fully adjusted"),or_hr = "HR",show_pgsname=T)
p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.45),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),hr_or = "HR")
a<-get_wh(p,unit = "cm")+2
png(paste(graphs_path,"/cox/cox_forestplotEPO001_cont_07Aug2023",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

# EPA001--------------------------
outcome="EPA001_up_75"

## analysis ----------------------
model_output_partial<-create_output_table_cox(
  data=data,adjustments=partial_adjustments,outcome=outcome,trainsplit = F,strata = strata)

model_output_full<-create_output_table_cox(
  data=data,adjustments=full_adjustments,outcome=outcome,trainsplit = F,namew = "_full",strata = strata)


### results table ------------------------
EPA001_cont<-combine_model_tables(
  tables=list(model_output_partial,model_output_full),
  name = c("Partially adjusted","Fully adjusted"),log_cox = "cox")
#saveRDS(EPA001_cont,paste(data_path,"/EPA001_cont_model_outcomes_05Jun2023.rds",sep=""))
flex_cont<-my_table(EPA001_cont,y=c(1:3,6),x=1:c(npgs+1))
flex_cont<-footnote(flex_cont,i=c(3,6),j=1,part="body",value = as_paragraph(c("Partial adjusted models adjust for Sex and age groups. Fully adjusted models additionally adjust for  waist-to-hip ratio, systolic and diastolic blood pressure, education attainmen level, smoking status,diabetes at baseline")),ref_symbols = c("*"))
flex_cont <- add_header_row(flex_cont, values =  c(" ", "European PRS","Non-European PRS"),
                            colwidths = c(1, 4,npgs-4), top = T)
flex_cont<-set_caption(
  flex_cont,
  caption = "Association of PRSs with CAD mortality(anywhere mention on the death certificate)",align_with_table = FALSE,
  word_stylename = "Table Caption",
  fp_p = fpp
)
flex_cont<-my_theme(flex_cont,y=c(1,2,5,8),set_padding = 3,fontsize = 8,header_border = 7)

saveRDS(flex_cont,paste(data_path,"/EPA001_cont_model_outcomes_flex_07Aug2023.rds",sep=""))

### forest plot ------------------------------
forest_table<-generate_forest_table(tables=list(model_output_partial,model_output_full),model_name = c("Partially adjusted","Fully adjusted"),or_hr = "HR",show_pgsname=T)
p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.45),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),hr_or = "HR")
a<-get_wh(p,unit = "cm")+2
png(paste(graphs_path,"/cox/cox_forestplotEPA001_cont_07Aug2023",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()








