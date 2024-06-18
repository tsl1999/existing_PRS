##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# Figure 3
# Author: Tianshu Liu
# Date: 30 May 2024
#################################
#call packages and data
##package and self-create function
library(xlsx)
rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
phenotype_path<-"/well/emberson/projects/mcps/data/phenotypes/"
data_path<-paste(working_path,'/data',sep='')
data_save_path<-paste(working_path,'/data/manuscript-age80',sep='')
graphs_path<-paste(working_path,'/graphs/manuscript-age80',sep='')
setwd(working_path)
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")  
#data readin------------------------------------------------
data<-readRDS(paste(data_path,"/FullData_ancestry_80_17Apr2024.rds",sep=""))
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1


#Figure 2----------------------------------------------------------------
outcome="prevalent_CHD_EPA"
adjustments_in=c("AGE","SEX")
mediator=c("smokegp2","diabetes_at_baseline","SBP","DBP","WHRATIO","EDU_LEVEL")
mediator_name<-c("Smoking status","Diabetes at baseline",
                 "Systolic Blood Pressure","Diastolic Blood Pressure","Waist-to-hip ratio",
                 "Education attainment level")


data<-data[is.na(data$smokegp2)==F&is.na(data$diabetes_at_baseline)==F&is.na(data$SBP)==F&is.na(data$DBP)==F&is.na(data$EDU_LEVEL)==F&is.na(data$WHRATIO)==F,]
pgs=which(colnames(data)%in%colnames(data%>%select(contains(colnames(model_table_top)[2:ncol(model_table_top)]))))
mediator_results<-list()
for(i in 1:npgs){
  pgs_name<-colnames(model_table_top)[i+1]
  forest_table<-run_mediator(data = data,adjustments = adjustments_in,mediator = mediator,pgs_name = pgs_name)
  colnames(forest_table)[1]<-"Model Adjustments"
  colnames(forest_table)[5]<-
    "X\u00b2 test p-value"
  mediator_results[[i]]<-forest_table
  if(i==1){
    append_i=FALSE
  }else{append_i=TRUE}
  
  write.xlsx(forest_table,file=paste(data_save_path,
                "stepwise_pgs_logistic.xlsx",sep="/"),
             sheetName = paste(pgs_name),append = append_i,row.names=T)
#write.csv(forest_table,paste(data_save_path,"/stepwise_",pgs_name,".csv",sep=""))
}


mediator_results_table<-data.frame()
for (i in 1:npgs){
  mediator_results_table<-rbind(mediator_results_table,
                                rep(" ", 10),
                                mediator_results[[i]])
}


mediator_results_table<-mediator_results_table[-1,]
saveRDS(mediator_results_table,paste(data_save_path,"/stepwise_analysis_allprs.rds",sep=""))


p<-plot_forest(forest_table =mediator_results_table ,xlim=c(0.9,1.4),ticks_at=c(0.9,1,1.1,1.2,1.3,1.4),
               theme,hr_or="OR",col_input=c(1,10,8,5),ci_column=2,primary_analysis = F,
               footnote_in = " ",)
p <- add_border(p,
                part = "body",
                col = c(3,4),row=1:nrow(mediator_results_table),
                gp = gpar(lwd = .5))
p<- edit_plot(p, row = c(seq(1,nrow(mediator_results_table),
                             by=(nrow(mediator_results_table)+1)%/%(npgs))),col=1,part = "body", 
              gp = gpar(fontface = "bold"))

p <- add_text(p, text = paste("OR per SD (95% CI)",sep=" "),
              part = "header", 
              col = 2:3,row=1,
              gp = gpar(fontface="bold",fontsize=7))

a<-get_wh(p,unit = "cm")




png(paste(graphs_path,"/Figure3_mediator_plot.png",sep=""), res = 300, width = a[1], height = a[2], units = "cm")
#pdf(paste(graphs_path,"/Figure2_mediator_plot.pdf",sep=""), width = a[1]/2.54, height = a[2]/2.54,paper="a4")


print(p)
dev.off()

