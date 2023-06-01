#### for age-at-risk analysis

rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')

setwd(working_path)

library(data.table)
library(tidyr)
library(tidyverse)
library(lubridate)
library(dplyr)
library(survival)
library(Epi)

# data readin-------------
data<-readRDS(paste(data_path,"FullData_standardisedPRS_31May2023.rds",sep="/"))
data$EPA001_up<-ifelse(is.na(data$EPA001)==F&data$EPA001==1,1,0)
data$EPA001_up<-as.numeric(data$EPA001_up)
data$EPO001_up<-ifelse(is.na(data$EPO001)==F&data$EPO001==1,1,0)
data$EPO001_up<-as.numeric(data$EPO001_up)
#data$death_01<-ifelse(is.na(data$DATE_OF_DEATH)==T,0,1)#19867 cases
#data$death_01_75<-ifelse(is.na(data$DATE_OF_DEATH)==F&data$AGE_followup<75,1,0)#12603 cases




data$DATE_OF_DEATH<-as.Date(data$DATE_OF_DEATH,"%d/%m/%Y")
data$DATE_RECRUITED<-as.Date(data$DATE_RECRUITED,"%d%b%Y")



#saveRDS(data,paste(data_path,"/FullData_standardisedPRS_mortality_31May2023.rds",sep=""))



#new codes available see 3.1 new--------------------------------------
data_age_risk<-data%>%select(IID,AGE,DATE_RECRUITED,DATE_OF_DEATH,AGE_followup)
data_age_risk$endpoint<-ifelse(is.na(data_age_risk$DATE_OF_DEATH)==F&data_age_risk$AGE_followup<75,1,0)
table(data_age_risk$endpoint)


agegroup<-c(35,45,55,65,75)#main analysis will use people die of CAD between 35-75, 85 and 100 age bands are set to ensure data completeness

# create row for each participant and each age group
data_age_risk_expand<-data_age_risk%>%crossing(tibble(XAgeGrp_start = agegroup[-length(agegroup)],
                                               XAgeGrp_end = agegroup[-1])) %>% 
  mutate(XAgeGrp = (XAgeGrp_start + XAgeGrp_end)/2)

#start each interval on the month and date of recruitment--------------------
##2/29 in leap years are solved using %m+% pt %m-%

data_age_risk_expand$start_int<-(as.Date(data_age_risk_expand$DATE_RECRUITED)%m-%years(data_age_risk_expand$AGE)%m+%years(data_age_risk_expand$XAgeGrp_start))
#remove intervals start after death date/censoring date and before recruitment date -----------------

data_age_risk_expand$end_date<-ifelse(data_age_risk_expand$endpoint==1,as.character(data_age_risk_expand$DATE_OF_DEATH),ifelse(
  data_age_risk_expand$endpoint==0,as.character("2021-01-01"),NA
))

data_age_risk_expand$end_date<-as.Date(data_age_risk_expand$end_date)

data_age_risk_expand<-data_age_risk_expand%>%filter(start_int <= as.Date(end_date))

data_age_risk_expand$end_int<-as.Date(data_age_risk_expand$DATE_RECRUITED)%m-%years(data_age_risk_expand$AGE)%m+%years(data_age_risk_expand$XAgeGrp_end)%m-%days(1)

data_age_risk_expand<-data_age_risk_expand%>%filter(end_int>= as.Date(DATE_RECRUITED))



data_age_risk_expand<-data_age_risk_expand%>%mutate(start_int = pmax(as.Date(DATE_RECRUITED), start_int),
                                                    end_int   = pmin(as.Date(end_date), end_int))

# calculate days / years from start of study---------------
data_age_risk_expand$endpoint<-ifelse(data_age_risk_expand$endpoint==1&data_age_risk_expand$end_int==data_age_risk_expand$DATE_OF_DEATH,1,0)





#calculate years in each age group
#data_age_risk_expand$years_int<-as.numeric((data_age_risk_expand$end_int-data_age_risk_expand$start_int)/365.25)
#data_age_risk_expand$years_int<-round(data_age_risk_expand$years_int,4)

data_age_risk_expand<-data_age_risk_expand%>%mutate(
  t_start_days   = interval(as.Date(DATE_RECRUITED), start_int) / days(1),
  t_end_days     = interval(as.Date(DATE_RECRUITED), end_int) / days(1) + 0.95,
  time_in        = round(interval(as.Date(DATE_RECRUITED), start_int) / years(1), 4),
  time_out       = round(interval(as.Date(DATE_RECRUITED), end_int) / years(1) + .95/365.25, 4))%>%mutate(
    years_int=time_out-time_in
  )



#change columns order
data_age_risk_expand<-data_age_risk_expand%>%select(IID,AGE,AGE_followup,DATE_RECRUITED,DATE_OF_DEATH,end_date,XAgeGrp,XAgeGrp_start,XAgeGrp_end,start_int,end_int,endpoint,time_in,time_out,years_int)


#merge with original data
data_age_risk_expand_final<-left_join(data_age_risk_expand,data)



#adjiust colnames ---------------
colnames(data_age_risk_expand_final)<-sub("_standardised","",colnames(data_age_risk_expand_final))
colnames(data_age_risk_expand_final)[51]<-"custom_Oni_Orisan"

#CAD mortality as endpoint--------------

data_age_risk_expand_final$endpoint_CAD_epa<-ifelse(data_age_risk_expand_final$endpoint==1&data_age_risk_expand_final$EPA001_up==1,1,0)
data_age_risk_expand_final$endpoint_CAD_epo<-ifelse(data_age_risk_expand_final$endpoint==1&data_age_risk_expand_final$EPO001_up==1,1,0)

#final checks--------------------------------------------
table(data_age_risk_expand_final$endpoint)==table(is.na(data_age_risk$DATE_OF_DEATH)==F&data_age_risk$AGE_followup<75)
table(data_age_risk_expand_final$endpoint_CAD_epa)==table(data$EPA001_up&data_age_risk$AGE_followup<75)
table(data_age_risk_expand_final$endpoint_CAD_epo)==table(data$EPO001_up&data_age_risk$AGE_followup<75)


try1<-coxph(Surv(time=time_in,time2 = time_out,endpoint_CAD_epo)~PGS000018,data=data_age_risk_expand_final)
summary(try1)



try2<-coxph(Surv(followup,data$EPO001_up&data_age_risk$AGE_followup<75)~PGS000018_standardised,data=data)
summary(try2)

print(round(coef(try1)-coef(try2),2))# expecting 0



#save dataset---------------------------
saveRDS(data_age_risk_expand_final,paste(data_path,"/FullData_standardisedPRS_age_at_risk_29May2023.rds",sep=""))


