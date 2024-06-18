##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# supplementary part 4 
# CAD mortality cox analysis
# Author: Tianshu Liu
# Date: 06 Jun 2024
#################################
#call packages and data
##package and self-create function
rm(list=ls())
library(survival)
library(dynpred)
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
data_path<-paste(working_path,'/data',sep='')
data_save_path<-paste(working_path,'/data/manuscript-age80/supplement',sep='')
graphs_path<-paste(working_path,'/graphs/manuscript-age80/supplement',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  
#data---------------------
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1

data<-readRDS(paste(data_path,"FullData_timesplit_80_17Apr2024.rds",sep="/"))
remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data<-data[,-c(remove_cat)]

# EPA001--------------------------
outcome="EPA001_up_80"
partial_adjustments=c("SEX",paste(rep("PC",7),seq(1,7),sep=""))
full_adjustments = c(partial_adjustments,c("WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline"))
strata="agegroup"

## analysis ----------------------
model_output_partial<-create_output_table_cox(
  data=data,adjustments=partial_adjustments,outcome=outcome,trainsplit = F,strata = strata,namew = "EPA_partial",ph=F)

model_output_full<-create_output_table_cox(
  data=data,adjustments=full_adjustments,outcome=outcome,trainsplit = F,namew = "EPA_full",strata = strata,ph=F)

write.csv(model_output_partial,paste(data_save_path,"/mortality_partial.csv",sep=""))
saveRDS(model_output_partial,paste(data_save_path,"/mortality_partial.rds",sep=""))
write.csv(model_output_full,paste(data_save_path,"/mortality_full.csv",sep=""))
saveRDS(model_output_full,paste(data_save_path,"/mortality_full.rds",sep=""))

### forest plot ------------------------------
forest_table<-generate_forest_table(tables=list(model_output_partial,model_output_full),
                                    model_name = c("Partially adjusted","Fully adjusted"),
                                    or_hr = "HR",show_pgsname=T)
p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.45),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),
               hr_or = "HR",footnote_in=" ")
a<-get_wh(p,unit = "cm")+2
png(paste(graphs_path,"/mortalityEPA_forest",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

##quintile plot----------------------------------------------------------

data<-readRDS(paste(data_path,"FullData_timesplit_80_17Apr2024.rds",sep="/"))
data$agegroup<-factor(data$agegroup,levels = c(1:5))
remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data_up<-data[,c(1:32,41:43,remove_cat[1]:ncol(data))]
partial_adjustments=c("SEX",paste(rep("PC",7),seq(1,7),sep=""))
strata="agegroup"
full_adjustments = c(partial_adjustments,c("WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline"))

outcome="EPA001_up_80"
## analysis ------------------
model_output_partial<-create_output_table_cox(
  trainsplit = F,data_up,adjustments=partial_adjustments,outcome=outcome,namew="partial",type = "cat",strata = strata,ph=F)

model_output_full<-create_output_table_cox(
  trainsplit = F,data_up,adjustments=full_adjustments,outcome=outcome,namew="full",type = "cat",strata = strata,ph=F)

## FAR plot ---------------------------
model_output_list<-list(model_output_partial,model_output_full)
#saveRDS(model_output_list,paste(data_path,"/EPO001_cat_cox_07Jun2023.rds",sep=""))
plot_FAR_gg(model_output_list,name=c("Partial","Full"),outcome="mortalityEPA",type="HR",
            data = data,graphs_path = graphs_path)







