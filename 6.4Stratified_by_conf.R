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


# data readin and set factors----------------------------------------------
data<-readRDS(paste(data_path,"/FullData_ancestry_22Jun2023.rds",sep=""))
model_table_top<-readRDS(paste(data_path,"/model_table_top.rds",sep=""))
data$EDU_LEVEL<-factor(data$EDU_LEVEL,levels = c(1:4),c("Uni/College","High school","Elementary","Other"),labels = c("University or College","High school","Elementary","Other"))
data$smokegp<-factor(data$smokegp,levels = c(1:5),labels = c("Never","Ex-smoker","Current (<5/day)","Current (5-14/day)","Current (>14/day)"))


full_adjustment= c("AGE","SEX","WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline")
confounder_check<-c(full_adjustment,"BASE_HBA1C","BMI","eurscore","amrscore")

all_name<-c("Age","Gender","Waist-Hip Ratio","Systolic blood pressure","Diastolic blood pressure",
            "Education level","Smoking status group","Diabetes at baseline","Baseline HbA1c level","Body Mass Index","European ancestry proportion","Indigenous ancestry proportion")
data$diabetes_at_baseline<-factor(data$diabetes_at_baseline,levels = c(0,1),labels = c("No","Yes"))
data$AGE_cut<-cut(data$AGE, breaks=c(35,55,65,74),
                  labels = c("35-54","55-64","65-74"),include.lowest=TRUE,right=F)

data$SBP_cut<-cut(data$SBP, breaks=c(0,125,140,250),
                  labels = c("<125","125-139",">=140"),
                  include.lowest=TRUE,right = F)

data$DBP_cut<-cut(data$DBP, breaks=c(0,80,90,190),
                  labels = c("<80","80-89",">=90"),
                  include.lowest=TRUE,right = F)

data$WHRATIO_cut<-cut(data$WHRATIO, breaks=c(0,0.85,0.9,2.7),
                      labels = c("<0.85","0.85-0.9",">=0.9"),
                      include.lowest=TRUE,right = F)


data$eurscore_cut<-cut(data$eurscore, breaks=c(0,0.2,0.4,1),
                      labels = c("<20%","20%-40%",">=40%"),
                      include.lowest=TRUE,right = F)

data$amrscore_cut<-cut(data$amrscore, breaks=c(0,0.6,0.8,1),
                       labels = c("<60%","60%-80%",">=80%"),
                       include.lowest=TRUE,right = F)

data$BMI_cut<-cut(data$BMI, breaks=c(0,18.5,25,30,35,100),
                       labels = c("<18.5","18.5-25","25-30","30-35",">=35"),
                       include.lowest=TRUE,right = F)

data$BASE_HBA1C_cut<-cut(data$BASE_HBA1C, breaks=c(0,6,6.5,20),
                  labels = c("Normal (<6.0%)","Prediabetes (6.0%-6.5%)","Diabetes (>=6.5%)"),
                  include.lowest=TRUE,right = F)








strata_in<-confounder_check

data_analysis<-data[,c(1:39,47:ncol(data))]
results_full_confounder<-readRDS(paste(data_path,"/Prevalent_CHD_EPA_cont_logistic_06Jun2023.rds",sep=""))
#run analysis EPA-------------------------
all_strata_outcome<-run_stratified(data_analysis, strata_in = confounder_check,outcome="prevalent_CHD_EPA",adjustment = full_adjustment)

forest_table_list<-stratified_analysis_table(all_strata_outcome,all_name,confounder_check,data=data_analysis)


for (i in 1:7){
  forest_input<-forest_table_list[[i]]
  max_ur<-max(as.numeric(forest_input$UR),na.rm=T)
  min_ur<-min(as.numeric(forest_input$LR),na.rm=T)
  
  p<-plot_forest(forest_table = forest_input,xlim=c(0.9,1.4),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),hr_or = "OR",theme_in = "tm")
  p<- edit_plot(p, row = c(1,5,8,12,16,20,25,31,34,38,44,48,52), col=1,
                gp = gpar(fontface = "bold"))
  
 
  p <- edit_plot(p, col = 1:7, 
                 row = 52, 
                 which = "background", 
                 gp = gpar(fill = "#f6eff7"))
  a<-get_wh(p,unit = "cm")+2
  
  png(paste(graphs_path,"/logistic/Logistic_forestplotEPA_confounder",colnames(forest_input)[1],".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
  print(p)
  dev.off()
  
  
}

#run analysis EPO-------------------------
all_strata_outcome<-run_stratified(data_analysis, strata_in = confounder_check,outcome="prevalent_CHD_EPO",adjustment = full_adjustment)

forest_table_list<-stratified_analysis_table(all_strata_outcome,all_name,confounder_check,data=data_analysis)


for (i in 1:7){
  forest_input<-forest_table_list[[i]]
  max_ur<-max(as.numeric(forest_input$UR),na.rm=T)
  min_ur<-min(as.numeric(forest_input$LR),na.rm=T)
  
  p<-plot_forest(forest_table = forest_input,xlim=c(0.9,1.4),ticks_at = c(0.9,1,1.1,1.2,1.3,1.4),hr_or = "OR",theme_in = "tm")
  p<- edit_plot(p, row = c(1,5,8,12,16,20,25,31,34,38,44,48,52), col=1,
                gp = gpar(fontface = "bold"))
  
 
  p <- edit_plot(p, col = 1:7, 
                 row = 52, 
                 which = "background", 
                 gp = gpar(fill = "#f6eff7"))
  a<-get_wh(p,unit = "cm")+2
  
  png(paste(graphs_path,"/logistic/Logistic_forestplotEPO_confounder",colnames(forest_input)[1],".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
  print(p)
  dev.off()
  
  
}

