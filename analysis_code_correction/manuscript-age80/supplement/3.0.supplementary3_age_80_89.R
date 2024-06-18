##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# Supplementary part 4
# Age 80-89
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
# data readin-------------------------------------------
data<-readRDS(paste(data_path,"FullData_standardisedPRS_after_80_17Apr2024.rds",sep="/"))
#data<-data[,c(1:32,35:ncol(data))]
colnames(data)[33:ncol(data)]<-sub("_standardised","",colnames(data)[33:ncol(data)])#remove standardised after all PRS names to make column names shorter
colnames(data)[40]<-"custom_Oni_Orisan"
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1
# checks---------------------------------------------
#80-90
table(data$BASE_CHD)
data<-data[data$AGE<90,]
data$EPA001_up<-ifelse(data$AGE_followup<90&data$EPA001==1,1,0)
data$prevalent_CHD_EPA<-ifelse(data$BASE_CHD==1,1,ifelse(data$EPA001_up==1&is.na(data$EPA001_up)==F,1,0))
# set adjustments-----------------------------------
partial_adjustments=c("AGE","SEX",paste(rep("PC",7),seq(1,7),sep=""))
full_adjustments = c(partial_adjustments,"WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline")

## EPA--------------------------------------
outcome="prevalent_CHD_EPA"

### without prs------------------------------
model_without_75<-discrimination_without_prs(
  train_data = data,test_data = data,outcome=outcome,
  partial_adjustments = partial_adjustments,full_adjustments = full_adjustments)
print(model_without_75)
#saveRDS(model_without_75,paste(data_path,"/Prevalent_CHD_EPA_withoutPRS80.rds",sep=""))

### analysis ----------------------------------------
model_output_partial<-create_output_table(
  trainsplit = F,data,adjustments=partial_adjustments,outcome=outcome,namew="partial_over80",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data,adjustments=full_adjustments,outcome=outcome,namew="full_over80",roc = F)

write.csv(model_output_partial,paste(data_save_path,"/age_80_89EPA_partial.csv",sep=""))
saveRDS(model_output_partial,paste(data_save_path,"/age_80_89EPA_partial.rds",sep=""))
write.csv(model_output_full,paste(data_save_path,"/age_80_89EPA_full.csv",sep=""))
saveRDS(model_output_full,paste(data_save_path,"/age_80_89EPA_full.rds",sep=""))

### forest table----------------------------------
forest_table<-generate_forest_table(tables=list(model_output_partial,model_output_full),
                                    model_name = c("Partially adjusted","Fully adjusted"),
                                    or_hr = "OR",show_pgsname = T)

p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.3),ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or = "OR",footnote_in=" ")
a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/age_80_89EPA_forest",".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
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
data_up<-data[,c(1:32,41:42,(ncol_data+1):ncol(data))]

model_output_partial<-create_output_table(
  trainsplit = F,data_up,adjustments=partial_adjustments,outcome="prevalent_CHD_EPA",
  namew="partial",type = "cat",roc = F)

model_output_full<-create_output_table(
  trainsplit = F,data_up,adjustments=full_adjustments,outcome="prevalent_CHD_EPA",
  namew="full",type = "cat",roc = F)


model_output_list<-list(model_output_partial,model_output_full)
#saveRDS(model_output_list,paste(data_path,"/prevalent_CHD_EPA_cat_logistic_06Jun2023.rds",sep=""))

plot_FAR_gg(model_output_list,name=c("Partial","Full"),
            outcome="age_80_89EPA",data = data,
            outcome_name = "CAD (baseline CAD and CAD mortality anywhere mentioned on the death certificate) for participants age 75-90",
            graphs_path = graphs_path)



