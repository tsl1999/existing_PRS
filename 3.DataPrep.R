#this file is to prepare dataset for analysis
rm(list=ls())

#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
phenotype_path<-"/well/emberson/projects/mcps/data/phenotypes/"
setwd(working_path)

library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pROC)
library(ROCR)
library(caret)


#readin data -----------------------------------------------
data<-readRDS(paste(data_path,"FullData_07May2023.rds",sep="/"))
Mortality<-data.frame(fread(paste(phenotype_path,"v2.1_DEATHS.csv",sep = "")))
mortality<-Mortality[!Mortality$grade%in%c("D","E","F","U","Z"),]
mortality_remove<-Mortality[Mortality$grade%in%c("D","E","F","U","Z"),]
table(mortality$grade,exclude=NULL)
table(mortality_remove$grade)
reg_link<-data.frame(fread(paste(phenotype_path,"RGN_LINK_IID.csv",sep = "")))
mortality_birth<-merge(mortality[,c(1,4,13)],reg_link[,c(1,3)], by="REGISTRO")
mortality_remove<-left_join(mortality_remove[,c(1,4,13)],reg_link[,c(1,3)], by="REGISTRO")

data<-left_join(data,mortality_birth[,c(2:4)],by="IID")
data_remove<-data[data$IID%in%mortality_remove$IID,]
data<-data[!data$IID%in%mortality_remove$IID,]
data$HDL_C<-as.numeric(data$HDL_C)
data$LDL_C<-as.numeric(data$LDL_C)
data$EDU_LEVEL<-factor(data$EDU_LEVEL,levels = c(1:4))
data$BASE_CHD<-factor(data$BASE_CHD,levels=c(0:1))
data$EPA001<-factor(data$EPA001,levels=c(0:1))
data$EPO001<-factor(data$EPO001,levels=c(0:1))
data$COYOACAN<-factor(data$COYOACAN,levels=c(0:1))
data$BASE_DIABETES<-factor(data$BASE_DIABETES,levels=c(0:1))
data$BASE_CANCER<-factor(data$BASE_CANCER,levels=c(0:1))
data$BASE_EMPHYSEMA<-factor(data$BASE_EMPHYSEMA,levels=c(0:1))
data$BASE_CIRR<-factor(data$BASE_CIRR,levels=c(0:1))
data$BASE_PEP<-factor(data$BASE_PEP,levels=c(0:1))
data$BASE_CKD<-factor(data$BASE_CKD,levels=c(0:1))
data$BASE_PAD<-factor(data$BASE_PAD,levels=c(0:1))
data$ACTIVITY<-factor(data$ACTIVITY,levels = c("NONE","< ONCE A WEEK","1-2 TIMES A WEEK","AT LEAST 3 TIMES A WEEK"))


data$SEX<-factor(data$SEX,levels = c(1,2),labels = c("Men","Women"))

data$prevalent_CHD_EPA<-ifelse(data$BASE_CHD==1,1,ifelse(is.na(data$EPA001)==F&data$EPA001==1,1,0))
data$prevalent_CHD_EPA<-factor(data$prevalent_CHD_EPA,levels=c(0,1))
data$prevalent_CHD_EPO<-ifelse(data$BASE_CHD==1,1,ifelse(is.na(data$EPO001)==F&data$EPO001==1,1,0))
data$prevalent_CHD_EPO<-factor(data$prevalent_CHD_EPO,levels=c(0,1))

table(data$prevalent_CHD_EPA)#8059   7616
table(data$prevalent_CHD_EPO)#7222   6834


data$AGE_75<-ifelse(data$AGE>=75,NA,data$AGE)
sum(is.na(data$AGE_75))
data$date_since_recruitment<-as.Date("2020-12-31","%Y-%m-%d")-as.Date(data$DATE_RECRUITED,"%d%b%Y")
data$yr_since_recruitment<-as.numeric(data$date_since_recruitment)/365.25
data$yrs_died_recruitment<-(as.Date(data$DATE_OF_DEATH,"%d/%m/%Y")-as.Date(data$DATE_RECRUITED,"%d%b%Y"))/365.25

data$AGE_followup<-ifelse(is.na(data$yrs_died_recruitment)==F,round(data$AGE+data$yrs_died_recruitment,4),ifelse(
  is.na(data$yrs_died_recruitment)==T, round(data$AGE+data$yr_since_recruitment,4),NA
))

data$anti_diabetic<-rowSums(data%>%select(contains("DRUG")))
data$anti_diabetic_medication<-ifelse(data$anti_diabetic>=1,1,0)
data$anti_diabetic_medication<-factor(data$anti_diabetic_medication,levels=c(0,1))

#cross check age -----------------------
#data$date_of_birth_up<-as.Date(data$date_of_birth,"%d/%m/%Y")
#data$age_check<-ifelse(is.na(data$date_of_birth)==F, (as.Date("2023-04-01","%Y-%m-%d")-data$date_of_birth_up)/365.25,NA)
#data$dage_death_check<-ifelse(is.na(data$date_of_birth)==F, (as.Date(data$DATE_OF_DEATH,"%d/%m/%Y")-data$date_of_birth_up)/365.25,NA)

#age_check<-data%>%select(IID,grade,DATE_RECRUITED,DATE_OF_DEATH,date_of_birth_up,AGE,AGE_75,yr_since_recruitment,yrs_died_recruitment,age_check,AGE_now,dage_death_check)
#age_check$age_at_recruit<-as.numeric((as.Date(age_check$DATE_RECRUITED,"%d%b%Y")-as.Date(age_check$date_of_birth_up))/365.25)

#age_check<-age_check[is.na(data$date_of_birth_up)==F&is.na(data$AGE_75)==F,]

#age_check$discrep_death<-age_check$AGE_now-age_check$dage_death_check
#age_check$discrep_recruit<-age_check$AGE-age_check$age_at_recruit
#data_keep<-data[is.na(data$AGE_75)==F,]

#sum(abs(age_check$discrep_death)>=2)#2132
#sum(abs(age_check$discrep_recruit)>=2)#2027

#age_discrp<-age_check[abs(age_check$discrep_recruit)>=2,]

#age_discrp<-age_discrp%>%select(IID,grade,DATE_RECRUITED,date_of_birth_up,AGE,age_at_recruit,discrep_recruit)

#table(age_discrp$grade)
#write.csv(age_discrp,"age_discrep.csv") 

#--------------------------------------------------------------
data_keep<-data[is.na(data$AGE_75)==F,]

data_keep$diabetic_lab<-ifelse(is.na(data_keep$BASE_HBA1C)==F&data_keep$BASE_HBA1C>=6.5,1,0)

data_keep$diabetes_at_baseline<-ifelse(data_keep$BASE_DIABETES==1|data_keep$diabetic_lab==1|data_keep$anti_diabetic_medication==1,1,0)
data_keep$diabetic_lab<-factor(data_keep$diabetic_lab,levels = c(0,1))

data_keep$smokegp<-factor(data_keep$smokegp,levels=c(1:5))

data_keep_further<-data_keep%>%select(IID,SEX,AGE,AGE_followup,BMI,BASE_CHD,BASE_CVD,BASE_DIABETES,BASE_HBA1C,
                                   smokegp,EDU_UNI,EDU_LEVEL,COYOACAN,WHRATIO,anti_diabetic_medication,SBP,DBP,
                                   HDL_C,LDL_C,prevalent_CHD_EPO,prevalent_CHD_EPA,EPO001,EPA001,contains(c("PC","PGS","custom")),
                                   diabetic_lab,diabetes_at_baseline,DATE_OF_DEATH,DATE_RECRUITED)



data_table_1<-data_keep%>%select(IID,SEX,AGE,BMI,BASE_CHD,BASE_CVD,BASE_DIABETES,BASE_HBA1C,BASE_CANCER,BASE_EMPHYSEMA,
                                 BASE_CIRR,BASE_PEP,BASE_CKD,BASE_PAD,INCOME,WAISTC,HIPC,ACTIVITY,
                                 smokegp,EDU_UNI,EDU_LEVEL,COYOACAN,WHRATIO,anti_diabetic_medication,SBP,DBP,
                                 HDL_C,LDL_C,prevalent_CHD_EPO,prevalent_CHD_EPA,EPO001,EPA001,diabetic_lab,diabetes_at_baseline,contains(c("PGS","custom")))
#standardise prs-------------------------------------
prs<-data%>%select(IID,contains(c("PGS","custom")))

for (i in 2:8){
  prs[,7+i]<-(prs[,i]-mean(prs[,i]))/sd(prs[,i])
  colnames(prs)[7+i]<-paste(colnames(prs)[i],"standardised",sep="_")
}


apply(prs[2:15],2,summary)
# plot density plots --------------------------------------

dt_melt_unstandardised<-melt(prs[,1:8],"IID")
ggplot(data = dt_melt_unstandardised, aes(x=value)) + geom_density( alpha = 0.4,fill="grey")+facet_wrap( ~ variable,scales = "free")


dt_melt_standardised<-melt(prs[,c(1,9:15)],"IID")
ggplot(data = dt_melt_standardised, aes(x=value)) + geom_density( alpha = 0.4,fill="grey")+facet_wrap( ~ variable)



#check difference between PRS for case and control-------------------------

prs_standardised<-merge(data_keep_further[,c(1:23,25:31,39,40:42)],prs[,c(1,9:15)],by="IID")

data_table_1_standardised<-merge(data_table_1[,1:34],prs[,c(1,9:15)],by="IID")
dt_melt_standardised<-melt(prs_standardised%>%select(IID,SEX,prevalent_CHD_EPA,prevalent_CHD_EPO,contains(c("PGS","custom"))),c("IID","SEX","prevalent_CHD_EPA","prevalent_CHD_EPO"))

dt_melt_standardised$prevalent_CHD_EPA<-factor(dt_melt_standardised$prevalent_CHD_EPA,levels = c(0,1),labels = c("No CAD","CAD case"))
dt_melt_standardised$prevalent_CHD_EPO<-factor(dt_melt_standardised$prevalent_CHD_EPO,levels = c(0,1),labels = c("No CAD","CAD case"))

ggplot(data = dt_melt_standardised, aes(x=value)) + 
  geom_density(aes(group=prevalent_CHD_EPA,fill=prevalent_CHD_EPA), alpha = 0.4)+xlab("Polygenic risk score")+
  facet_wrap( ~ variable,scales = "free")+scale_fill_discrete(name = " ")+ggtitle("PRS distribution by case/non-case (EPA)")

#ggsave(paste(graphs_path,"/PRS_standardised_byCaseCtrl_EPA_31May2023.png",sep=""),width = 10, height = 8)



ggplot(data = dt_melt_standardised, aes(x=value)) + 
  geom_density(aes(group=SEX,fill=SEX), alpha = 0.4)+xlab("Polygenic risk score")+
  facet_wrap( ~ variable,scales = "free")+scale_fill_discrete(name = " ")+ggtitle("PRS distribution by gender")

  #ggsave(paste(graphs_path,"/PRS_standardised_bysex_31May2023.png",sep=""),width = 10, height = 8)



ggplot(data = dt_melt_standardised, aes(x=value)) + xlab("Polygenic risk score")+
  geom_density(aes(group=prevalent_CHD_EPO,fill=prevalent_CHD_EPO), alpha = 0.4)+
  facet_wrap( ~ variable,scales = "free")+scale_fill_discrete(name = " ")+ggtitle("PRS distribution by case/non-case (EPO)")
#ggsave(paste(graphs_path,"/PRS_standardised_byCaseCtrl_EPO_31May2023.png",sep=""),width = 10, height = 8)

data_keep$prevalent_CHD_EPA<-factor(data_keep$prevalent_CHD_EPA,levels = c(0,1),labels = c("No CAD","CAD case"))



ggplot(data = data_keep, aes(x=AGE)) + 
  geom_density(aes(group=prevalent_CHD_EPA,fill=prevalent_CHD_EPA), alpha = 0.4)+scale_fill_discrete(name = " ")+ggtitle("Age at baseline distribution by case/non-case (EPA)")
#ggsave(paste(graphs_path,"/Age_at_recruitment_byCaseCtrl_EPA_31May2023.png",sep=""),width = 10, height = 8)


ggplot(data = data_keep, aes(x=AGE_followup)) + xlab("Age at risk")+
  geom_density(aes(group=prevalent_CHD_EPA,fill=prevalent_CHD_EPA), alpha = 0.4)+scale_fill_discrete(name = " ")+ggtitle("Age at risk distribution by case/non-case (EPA)")

#ggsave(paste(graphs_path,"/Age_now_byCaseCtrl_EPA_31May2023.png",sep=""),width = 10, height = 8)


#saveRDS(prs,paste(data_path,"PRS_standardised_31Apr2023.rds",sep="/"))
#saveRDS(prs_standardised,paste(data_path,"FullData_standardisedPRS_31May2023.rds",sep="/"))
#saveRDS(data_table_1_standardised,paste(data_path,"Table1_data_standardisedPRS_31May2023.rds",sep="/"))
