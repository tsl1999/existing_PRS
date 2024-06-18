##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# supplementary part 6
# confounder stratified analysis
# Author: Tianshu Liu
# Date: 10 June 2023
#################################
#call packages and data
##package and self-create function
#call packages and data
rm(list=ls())
library(xlsx)
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
data_path<-paste(working_path,'/data',sep='')
data_save_path<-paste(working_path,'/data/manuscript-age80/supplement',sep='')
graphs_path<-paste(working_path,'/graphs/manuscript-age80/supplement',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")

# data readin and set factors----------------------------------------------
data<-readRDS(paste(data_path,"/FullData_ancestry_80_17Apr2024.rds",sep=""))
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1

data$EDU_LEVEL<-factor(data$EDU_LEVEL,levels = c(1:4),c("Uni/College","High school","Elementary","Other"),labels = c("University or College","High school","Elementary","Other"))
data$smokegp2<-factor(data$smokegp2,levels = c(1:5),labels = c("Never","Ex-smoker","Current (<daily)","Current (<10/day)","Current (>=10/day)"))#setting this as category will change the OR of the main analysis,but only very slightly


full_adjustment= c("AGE","SEX","WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline")
confounder_check<-c(full_adjustment,"BASE_HBA1C","BMI","eurscore","amrscore")
full_adjustment= c("AGE","SEX",paste(rep("PC",7),seq(1,7),sep=""),"WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline")
all_name<-c("Age at last follow up","Gender","Waist-Hip Ratio","Systolic blood pressure","Diastolic blood pressure",
            "Education level","Smoking status group","Diabetes at baseline","Baseline HbA1c level","Body Mass Index","European ancestry proportion","Indigenous ancestry proportion")
data$diabetes_at_baseline<-factor(data$diabetes_at_baseline,levels = c(0,1),labels = c("No","Yes"))
data$AGE_cut<-cut(data$AGE_followup, breaks=c(35,55,65,79),
                  labels = c("35-54","55-64","65-79"),include.lowest=TRUE,right=F)

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

remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))
data_analysis<-data[,-remove_cat]

data_female<-data_analysis[data_analysis$SEX=="Women",]
data_male<-data_analysis[data_analysis$SEX=="Men",]


#results_full_confounder<-readRDS(paste(data_path,"/Prevalent_CHD_EPA_cont_logistic_flex_07Aug2023.rds",sep=""))
#run analysis EPA-------------------------
all_strata_outcome<-run_stratified(data_analysis, strata_in = confounder_check,outcome="prevalent_CHD_EPA",
                                   adjustment = full_adjustment,se=T)

for (i in 1:length(confounder_check)){
  data_save<-all_strata_outcome[[i]]
  data_save_name<-names(all_strata_outcome)[i]
  for(j in 1:length(data_save)){
    if(i==1&j==1){
      append_i=FALSE
    }else{append_i=TRUE}
    data_strata<-data_save[[j]]
    data_strata_name<-names(data_save)[j]
    data_strata_name<-sub("/","_",data_strata_name)
    write.xlsx(data_strata,file=paste(data_save_path,
                                      "/confounder_stratified_analysis.xlsx",sep=""),
               sheetName = paste(data_save_name,"_",data_strata_name,sep=""),
               append = append_i,row.names=T)
    
  }

}
write.xlsx(all_strata_outcome[[length(all_strata_outcome)]],file=paste(data_save_path,
                                  "/confounder_stratified_analysis.xlsx",sep=""),
           sheetName = paste("Overall"),
           append = append_i,row.names=T)


saveRDS(all_strata_outcome,paste(data_save_path,
                                 "/confounder_stratified_analysis.rds",sep="/"))
forest_table_list<-stratified_analysis_table(all_strata_outcome,all_name,confounder_check,data=data_analysis)
for (i in 1:length(forest_table_list)){
  if(i==1){
    append_i=FALSE
  }else{append_i=TRUE}
  write.xlsx(forest_table_list[[i]],file=paste(data_save_path,
                                    "/confounder_stratified_forest.xlsx",sep=""),
             sheetName = paste(names(forest_table_list)[i]),
             append = append_i,row.names=F)
}
saveRDS(forest_table_list,paste(data_save_path,
                                 "/confounder_forest_tablelist.rds",sep="/"))




for (i in 1:npgs){
  forest_input<-forest_table_list[[i]]
  max_ur<-max(as.numeric(forest_input$UR),na.rm=T)
  min_ur<-min(as.numeric(forest_input$LR),na.rm=T)
  
  p<-plot_forest(forest_table = forest_input,xlim=c(0.8,1.5),ticks_at = c(0.8,0.9,1,1.1,1.2,1.3,1.4,1.5),
                 hr_or = "OR",theme_in = "tm",footnote_in = " ")
  p<- edit_plot(p, row = c(1,5,8,12,16,20,25,31,34,38,44,48,52), col=1,
                gp = gpar(fontface = "bold"))
  
  
  p <- edit_plot(p, col = 1:npgs, 
                 row = 52, 
                 which = "background", 
                 gp = gpar(fill = "#f6eff7"))
  a<-get_wh(p,unit = "cm")+2
  
  png(paste(graphs_path,"/forestEPA_confounder_80_",names(forest_table_list)[i],".png",sep=""), res = 200, width = a[1], height = a[2], units = "cm")
  print(p)
  dev.off()
  
  
}

#each confounder ------------------------------------------------
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/analysis_code_correction/manuscript-age80/auto_forest.R")
overall<-all_strata_outcome[[length(all_strata_outcome)]]
for(i in 1:(length(all_strata_outcome)-1)){
  confounder<-all_strata_outcome[[i]]
  confounder[[length(confounder)+1]]<-overall
  names(confounder)[length(confounder)]<-"Overall"

  forest_table<-generate_forest_table(
    tables=confounder,
    model_name = c(names(confounder)),or_hr = "OR",show_pgsname = T)
  colnames(forest_table)[1]<-all_name[i]
  
  forest_table$se<-1
  if(i==10){
    ci_size=25
  }else{
    ci_size=40
  }
  summary_row_in<-c(which(grepl("\U{00A0}\U{00A0}",forest_table[,1])==F)-1,nrow(forest_table))[2:9]
  
  p<-plot_forestplot(forest_table = forest_table,col_input=c(1:3,10,7,9,8),ci_column=4,fontsize=7,
                     estimate_colname="estimate",ur_colname="UR",lr_colname = "LR",se_colname = "se",
                     summary_row_in= summary_row_in,ci_size = ci_size,summary_size = 20,
                     xjust = 0.35,log=F,xlim_min = 0.8,xlim_max = 1.5,
                     x_ticks=c(seq(0.8,1.5,0.1)))
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 1:7, gp = gpar(lwd = 1))
  p <- add_text(p, text = "OR per SD (95%CI)",
                part = "header", 
                col = 4:5,row=1,
                gp = gpar(fontface="bold",fontsize=7))
  p <- add_border(p, 
                  part = "body", 
                  col = c(2:3,5:7),row=1:nrow(forest_table),
                  gp = gpar(lwd = .5))
  
  
  a<-get_wh(p,unit = "cm")+2
  png(paste(graphs_path,"/",names(all_strata_outcome)[i],"_stratified_age80",".png",sep=""), res = 300, width = a[1], height = a[2], units = "cm")
  print(p)
  dev.off()
}


#35-89------------------------------------------------------
age_last_follow<-all_strata_outcome$AGE
age_80_89_full<-readRDS(paste(data_path,"/manuscript-age80/supplement/age_80_89EPA_full.rds",sep=""))

age_last_follow[[4]]<-age_80_89_full
names(age_last_follow)[4]<-"80-89"
forest_age<-generate_forest_table(tables=age_last_follow,
                                  model_name = names(age_last_follow),or_hr = "OR",
                                  show_pgsname = T)
colnames(forest_age)[1]<-"Age at last follow up"
p<-plot_forest(forest_table = forest_age,xlim=c(0.9,1.4),
               ticks_at = c(0.9,1,1.1,1.2,1.3),hr_or = "OR",
               footnote_in=" ")
a<-get_wh(p,unit = "cm")+2
png(paste(graphs_path,"/age_followup35-90_stratified",".png",sep=""), res = 300, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

