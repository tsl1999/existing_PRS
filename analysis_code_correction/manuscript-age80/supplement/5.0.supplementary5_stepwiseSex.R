##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# Supplementary part 5
# Stepwise by sex
# Author: Tianshu Liu
# Date: 10 June 2024
#################################
#call packages and data
library(xlsx)
rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
data_path<-paste(working_path,'/data',sep='')
data_save_path<-paste(working_path,'/data/manuscript-age80/supplement',sep='')
graphs_path<-paste(working_path,'/graphs/manuscript-age80/supplement',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  

# data readin and set factors----------------------------------------------
data<-readRDS(paste(data_path,"/analysis_data_80_17Apr2024.rds",sep=""))
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1
data<-data[is.na(data$smokegp2)==F&is.na(data$diabetes_at_baseline)==F&is.na(data$SBP)==F&is.na(data$DBP)==F&is.na(data$EDU_LEVEL)==F&is.na(data$WHRATIO)==F,]
#mediator analysis

outcome="prevalent_CHD_EPA"
adjustments_in=c("AGE","SEX",paste(rep("PC",7),seq(1,7),sep=""))
mediator=c("smokegp2","diabetes_at_baseline","SBP","DBP","WHRATIO","EDU_LEVEL")
mediator_name<-c("Smoking status","Diabetes at baseline",
                 "Systolic Blood Pressure","Diastolic Blood Pressure","Waist-to-hip ratio",
                 "Education attainment level")
mediator_results_male_female<-list()
mediator_results_sex<-list()
for (j in c("Women","Men")){
  adjustments_in=c("AGE")
  data_in<-data[data$SEX==j,]
  for(i in c(1:npgs)){
    pgs_name<-colnames(model_table_top)[i+1]
    forest_table<-run_mediator(data = data_in,adjustments = adjustments_in,mediator = mediator,pgs_name = pgs_name)
    colnames(forest_table)[1]<-"Model Adjustments"
    colnames(forest_table)[5]<-"X\u00b2 test p-value"
    mediator_results_male_female[[i]]<-forest_table
    if(i==1){
      append_i=FALSE
    }else{append_i=TRUE}
    
    write.xlsx(forest_table,file=paste(data_save_path,
                                                   "/stepwise_logistic_",j,".xlsx",sep=""),
               sheetName = paste(pgs_name),append = append_i,row.names=T)
  }
  
  mediator_results_sex[[j]]<-mediator_results_male_female
}

saveRDS(mediator_results_sex,paste(data_save_path,"/stepwise_logistic_allprs_women_men.rds",sep=""))


mediator_results_table_female<-data.frame()
mediator_results_table_male<-data.frame()
for (i in c(1:npgs)){
  mediator_results_table_female<-rbind(mediator_results_table_female,
                                       rep("",10),
                                       mediator_results_sex[[1]][[i]])
  mediator_results_table_male<-rbind(mediator_results_table_male,
                                     rep("",10),
                                     mediator_results_sex[[2]][[i]])
}

mediator_results_table_female<-mediator_results_table_female[-1,]
mediator_results_table_male<-mediator_results_table_male[-1,]
colnames(mediator_results_table_female)[c(2:7,9)]<-paste(colnames(mediator_results_table_female)[c(2:7,9)],sep=" ")

colnames(mediator_results_table_male)[c(2:7,9)]<-paste(" ",colnames(mediator_results_table_male)[c(2:7,9)],sep=" ")
colnames(mediator_results_table_male)[8]<-"    "
colnames(mediator_results_table_male)[10]<-paste(colnames(mediator_results_table_male)[10],"      ")
mediator_results_combined<-cbind(mediator_results_table_female,mediator_results_table_male[2:10])

colnames(mediator_results_combined)[1]<-"Model Adjustments"
rownames(mediator_results_combined)<-seq(1,nrow(mediator_results_combined))

write.csv(mediator_results_combined,paste(data_save_path,"/stepwise_sex_forest.csv",sep=""))
saveRDS(mediator_results_combined,paste(data_save_path,"/stepwise_sex_forest.rds",sep=""))

p <- forest(mediator_results_combined[c(1,9,10,8,5,18,19,17,14)],
            est =list(as.numeric(mediator_results_combined$estimate),as.numeric(mediator_results_combined$`  estimate`)),
            lower = list(as.numeric(mediator_results_combined$LR),as.numeric(mediator_results_combined$`  LR`)), 
            upper = list(as.numeric(mediator_results_combined$UR),as.numeric(mediator_results_combined$`  UR`)),
            ci_column = c(3,7),
            ref_line = 1,
            arrow_lab = c("Lower risk of CAD", "Higher risk of CAD"),
            xlim = c(0.9,1.5),
            ticks_at = seq(0.9,1.5,by=0.1),
            footnote = " ",
            theme=tm2)

p <- add_border(p, 
                part = "header", 
                row = 1,
                col = 1:9, gp = gpar(lwd = 1))

p <- add_border(p, 
                part = "body", 
                col = c(2,4:6,8,9),row=1:nrow(mediator_results_combined),
                gp = gpar(lwd = .5))
p <- add_border(p, 
                part = "body", 
                col = c(5),row=1:nrow(mediator_results_combined),where = c("right"),
                gp = gpar(lwd = .5))
p<- edit_plot(p, row = c(seq(1,nrow(mediator_results_combined),
                             by=(nrow(mediator_results_combined)+1)%/%(npgs))),col=1,part = "body", 
              gp = gpar(fontface = "bold"))

p <- add_text(p, text = paste("OR per SD (95% CI)",sep=" "),
              part = "header", 
              col = 3:4,row=1,
              gp = gpar(fontface="bold",fontsize=7))

p <- add_text(p, text = paste("Women",sep=" "),
              part = "header", 
              col = 2:6,row=0,
              gp = gpar(fontface="bold",fontsize=7))

p <- add_text(p, text = paste("OR per SD (95% CI)",sep=" "),
              part = "header", 
              col = 7:8,row=1,
              gp = gpar(fontface="bold",fontsize=7))
p <- add_text(p, text = paste("Men",sep=" "),
              part = "header", 
              col = 7:10,row=0,
              gp = gpar(fontface="bold",fontsize=7))

a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/mediator_plot_sex_80",".png",sep=""), res = 300, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()

