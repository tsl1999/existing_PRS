#call packages and data
##package and self-create function
rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
main_path<-paste(working_path,'/data/manuscript-age80/',sep='')
supp_path<-paste(working_path,'/data/manuscript-age80/supplement',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R") 

model_table_top<-readRDS(paste(working_path,"/data/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1
#compare partial analysis-----------------------------------------------
primary_partial<-readRDS(paste(main_path,"/primary_partial.rds",sep=""))
base_epo_partial<-readRDS(paste(supp_path,"/baselineEPO_partial.rds",sep=""))
third_partial<-readRDS(paste(supp_path,"/thirdEPA_partial.rds",sep=""))
base_partial<-readRDS(paste(supp_path,"/baseline_partial.rds",sep=""))
#age_80_89<-readRDS(paste(supp_path,"/age_80_89EPA_partial.rds",sep=""))

model_list<-list(primary_partial,base_epo_partial,base_partial,third_partial)
forest_table<-generate_forest_table(
  tables=model_list,
  model_name = c("Primary model","Baseline+primary cause of death","Baseline","Third degree unrelated"),or_hr = "OR",show_pgsname = T)
p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.4),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),hr_or = "OR",
               footnote_in=" ")

a<-get_wh(p,unit = "cm")+2
png(paste(working_path,"/graphs/manuscript-age80/supplement/all_model_compare_partial",".png",sep=""), res = 300, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

#compare full analysis-----------------------------------------------
primary_full<-readRDS(paste(main_path,"/primary_full.rds",sep=""))
base_epo_full<-readRDS(paste(supp_path,"/baselineEPO_full.rds",sep=""))
third_full<-readRDS(paste(supp_path,"/thirdEPA_full.rds",sep=""))
baseline_full<-readRDS(paste(supp_path,"/baseline_full.rds",sep=""))
#age_80_89_full<-readRDS(paste(supp_path,"/age_80_89EPA_full.rds",sep=""))

model_list<-list(primary_full,base_epo_full,third_full,baseline_full)
forest_table<-generate_forest_table(
  tables=model_list,
  model_name = c("Primary model","Baseline+primary cause of death","Baseline CAD","Third degree unrelated"),or_hr = "OR",show_pgsname = T)
p<-plot_forest(forest_table = forest_table,xlim=c(0.9,1.4),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),hr_or = "OR",
               footnote_in=" ")

a<-get_wh(p,unit = "cm")+2
png(paste(working_path,"/graphs/manuscript-age80/supplement/all_model_compare_full",".png",sep=""), res = 300, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()
