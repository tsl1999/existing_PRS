
rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
setwd(working_path)

library(data.table)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(lubridate)
library(dplyr)
library(pROC)
library(ROCR)
library(caret)
library(corrplot)
library(DescTools)
library(flextable)
library(forestploter)

source("/gpfs3/well/emberson/users/hma817/codes/R_TablesFunction_27Sep2022_TL.r")
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")

# data readin-------------
data<-readRDS(paste(data_path,"FullData_standardisedPRS_29May2023.rds",sep="/"))
data$EPA001_up<-ifelse(is.na(data$EPA001)==F&data$EPA001==1,1,0)
data$EPA001_up<-as.numeric(data$EPA001_up)

data$DATE_OF_DEATH<-as.Date(data$DATE_OF_DEATH,"%d/%m/%Y")
data$DATE_RECRUITED<-as.Date(data$DATE_RECRUITED,"%d%b%Y")



data_age_risk<-data%>%select(IID,AGE,DATE_RECRUITED,DATE_OF_DEATH,AGE_followup)
data_age_risk$endpoint<-ifelse(is.na(data_age_risk$DATE_OF_DEATH)==T,0,1)

agegroup<-c(35,55,65,75,85,100)

# create row for each participant and each age group
data_age_risk_expand<-data_age_risk%>%crossing(tibble(XAgeGrp_start = agegroup[-length(agegroup)],
                                               XAgeGrp_end = agegroup[-1])) %>% 
  mutate(XAgeGrp = (XAgeGrp_start + XAgeGrp_end)/2)

#start each interval on the month and date of recruitment--------------------
#2/29
data_age_risk_expand$start_int<-(as.Date(data_age_risk_expand$DATE_RECRUITED)-(data_age_risk_expand$AGE)*365.25+(data_age_risk_expand$XAgeGrp_start)*365.25)
#remove intervals start after death date/censoring date and before recruitment date -----------------

data_age_risk_expand$end_date<-ifelse(data_age_risk_expand$endpoint==1,as.character(data_age_risk_expand$DATE_OF_DEATH),ifelse(
  data_age_risk_expand$endpoint==0,as.character("2021-01-01"),NA
))

data_age_risk_expand$end_date<-as.Date(data_age_risk_expand$end_date)

data_age_risk_expand<-data_age_risk_expand%>%filter(start_int <= as.Date(end_date))

data_age_risk_expand$end_int<-as.Date(data_age_risk_expand$DATE_RECRUITED)-(data_age_risk_expand$AGE)*365.25+(data_age_risk_expand$XAgeGrp_end)*365.25-1

data_age_risk_expand<-data_age_risk_expand%>%filter(end_int>= as.Date(DATE_RECRUITED))



data_age_risk_expand<-data_age_risk_expand%>%mutate(start_int = pmax(as.Date(DATE_RECRUITED), start_int),
                                                    end_int   = pmin(as.Date(end_date), end_int))

# calculate days / years from start of study---------------
data_age_risk_expand$endpoint<-ifelse(data_age_risk_expand$endpoint==1&data_age_risk_expand$end_int==data_age_risk_expand$DATE_OF_DEATH,1,0)

data_age_risk_expand<-data_age_risk_expand%>%mutate(
t_start_days   = -(as.Date(DATE_RECRUITED)-start_int) ,
t_end_days     = -(as.Date(DATE_RECRUITED)- end_int)  + 0.95,
time_in        = -round((as.Date(DATE_RECRUITED)-start_int) / 365.25, 4),
time_out       = -round((as.Date(DATE_RECRUITED)-end_int) / 365.25 + .95/365.25, 4)
)


data_age_risk_expand$years_in_int<-as.numeric(round(data_age_risk_expand$time_out-data_age_risk_expand$time_in,2))
data_age_risk_expand_final<-left_join(data_age_risk_expand,data)



#adjiust colnames ---------------
colnames(data_age_risk_expand_final)<-sub("_standardised","",colnames(data_age_risk_expand_final))
colnames(data_age_risk_expand_final)[53]<-"custom_Oni_Orisan"

#change endpoint from all mortality to CAD mortality--------------


data_age_risk_expand_final$endpoint_CAD<-ifelse(data_age_risk_expand_final$endpoint==1&data_age_risk_expand_final$EPA001_up==1,1,0)

#try<-coxph(Surv(years_in_int,endpoint)~PGS000011+strata(as.factor(XAgeGrp)),data=data_age_risk_expand_final)
#summary(try)



saveRDS(data_age_risk_expand_final,paste(data_path,"/FullDate_standardisedPRS_age_at_risk_29May2023.rds",sep=""))


data$IID[data$EPA001_up==1][!data$IID[data$EPA001_up==1]%in%data_age_risk_expand_final$IID[data_age_risk_expand_final$endpoint_CAD==1]==T]



