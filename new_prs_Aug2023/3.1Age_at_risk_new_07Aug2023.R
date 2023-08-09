#check
rm(list=ls())

library(data.table)
library(tidyr)
library(tidyverse)
library(lubridate)
library(dplyr)
library(survival)
library(Epi)
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
setwd(working_path)

model_table_top<-readRDS(paste(data_path,"/model_table_top_07Aug2023.rds",sep=""))
npgs<-ncol(model_table_top)-1
data<-readRDS(paste(data_path,"FullData_standardisedPRS_07Aug2023.rds",sep="/"))
colnames(data)<-sub("_standardised","",colnames(data))
colnames(data)[42]<-"custom_Oni_Orisan"
data$EPA001_up<-ifelse(is.na(data$EPA001)==F&data$EPA001==1,1,0)
data$EPA001_up<-as.numeric(data$EPA001_up)
data$EPO001_up<-ifelse(is.na(data$EPO001)==F&data$EPO001==1,1,0)
data$EPO001_up<-as.numeric(data$EPO001_up)

data$DATE_OF_DEATH<-as.Date(data$DATE_OF_DEATH,"%d/%m/%Y")
data$DATE_RECRUITED<-as.Date(data$DATE_RECRUITED,"%d%b%Y")


data$age_at_entry<-data$AGE


#censor at age 75---------------------------------------------
data$EPO001_up_75<-ifelse(data$EPO001_up==1&data$AGE_followup<75,1,0)
data$EPA001_up_75<-ifelse(data$EPA001_up==1&data$AGE_followup<75,1,0)
data$AGE_followup<-ifelse(data$AGE_followup>=75,75,data$AGE_followup)
data$followup<-data$AGE_followup-data$age_at_entry
ncol_data<-ncol(data)
for (i in 1:npgs){
  
  quintile_cutoff<-quantile(data[,34+i],probs=seq(0,1,0.2),na.rm=F)
  
  data[,ncol_data+i]<-cut(data[,34+i], breaks=quintile_cutoff,labels = c(1:5),include.lowest=TRUE)
  colnames(data)[ncol_data+i]<-paste(colnames(data)[34+i],"cat",sep="_")
  
}

#saveRDS(data,paste(data_path,"/analysis_data_31May2023.rds",sep=""))

#cut data into age groups-----------------------------------------------------
timesplit <- survSplit(Surv(time=AGE,time2=AGE_followup, EPO001_up_75) ~., data,start="tstart",
                  cut=c(35,45,55,65,75), episode ="agegroup",end="tstop",zero=0)


#compute time at entry and time exit------------------------------------------
timesplit$time_in<-timesplit$tstart-timesplit$age_at_entry
timesplit$time_out<-timesplit$tstop-timesplit$age_at_entry

#timesplit<-timesplit[is.na(timesplit$tstart)==F,]


#include EPA outcome------------------------------------------------
timesplit$EPA001_up_75<-ifelse(timesplit$EPA001_up_75==1&timesplit$age_at_entry+timesplit$followup==timesplit$tstop,1,0)



table(data$EPO001_up_75)["1"]==table(timesplit$EPO001_up_75)["1"]

table(data$EPA001_up_75)["1"]==table(timesplit$EPA001_up_75)["1"]






#checks----------------


try1<-coxph(Surv(time=time_in,time2 = time_out,EPO001_up_75)~PGS000018,data=timesplit)
summary(try1)

try2<-coxph(Surv(followup,EPO001_up_75)~PGS000018,data=data)
summary(try2)

print(round(coef(try1)-coef(try2),2))


try3<-coxph(Surv(time=time_in,time2 = time_out,EPO001_up_75)~PGS000018_cat,data=timesplit)
summary(try3)

try4<-coxph(Surv(followup,EPO001_up_75)~PGS000018_cat,data=data)
summary(try4)

print(round(coef(try3)-coef(try4),2))


saveRDS(timesplit,paste(data_path,"/FullData_timesplit_07Aug2023.rds",sep=""))
saveRDS(data,paste(data_path,"/FullData_standardisedPRS_mortality_07Aug2023.rds",sep=""))
