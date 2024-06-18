rm(list=ls())
#set working directory load package-----------------
working_path<-'/well/emberson/users/hma817/projects/existing_PRS/correction'
data_path<-paste(working_path,'/data',sep='')
graphs_path<-paste(working_path,'/graphs',sep='')
setwd(working_path)

library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pROC)
library(ROCR)
library(caret)
library(corrplot)
library(DescTools)
library(flextable)
library(officedown)
library(scales)

source("/gpfs3/well/emberson/users/hma817/codes/R_TablesFunction_27Sep2022_TL.r")
source("/gpfs3/well/emberson/users/hma817/projects/existing_PRS/0.1.utils.R")

##data-----------------------------
data<-readRDS(paste(data_path,"FullData_standardisedPRS_80_17Apr2024.rds",sep="/"))
table1_data<-readRDS(paste(data_path,"Table1_data_standardisedPRS_80_17Apr2024.rds",sep="/"))
colnames(data)[35:ncol(data)]<-sub("_standardised","",colnames(data)[35:ncol(data)])#remove standardised after all PRS names to make column names shorter
colnames(data)[42]<-"custom_Oni_Orisan"# change custom_Oni-Orisan to Oni_Orisan so that formula() can read it
model_table_top<-readRDS(paste(data_path,"/model_table_top_17Apr2024.rds",sep=""))
npgs<-ncol(model_table_top)-1
model_table_top<-data.frame(model_table_top)

data_ancestry<-data.frame(fread(paste(data_path,"/ancestry_gen.csv",sep="")))
colnames(data_ancestry)[1]<-"IID"
table1_data<-left_join(table1_data,data_ancestry,by="IID")

##formatting-----------------------------------------------
table1_data$BASE_CHD<-factor(table1_data$BASE_CHD,levels = c(0,1),labels = c("No","Yes"))
table1_data$BASE_CVD<-factor(table1_data$BASE_CVD,levels = c(0,1),labels = c("No","Yes"))
table1_data$BASE_CANCER<-factor(table1_data$BASE_CANCER,levels = c(0,1),labels = c("No","Yes"))
table1_data$BASE_EMPHYSEMA<-factor(table1_data$BASE_EMPHYSEMA,levels = c(0,1),labels = c("No","Yes"))
table1_data$BASE_CIRR<-factor(table1_data$BASE_CIRR,levels = c(0,1),labels = c("No","Yes"))
table1_data$BASE_PEP<-factor(table1_data$BASE_PEP,levels = c(0,1),labels = c("No","Yes"))
table1_data$BASE_CKD<-factor(table1_data$BASE_CKD,levels = c(0,1),labels = c("No","Yes"))
table1_data$BASE_PAD<-factor(table1_data$BASE_PAD,levels = c(0,1),labels = c("No","Yes"))
table1_data$BASE_DIABETES<-factor(table1_data$BASE_DIABETES,levels = c(0,1),labels = c("No","Yes"))
table1_data$base_other<-ifelse(table1_data$BASE_EMPHYSEMA=="Yes"|table1_data$BASE_CIRR=="Yes"|table1_data$BASE_PEP=="Yes"|table1_data$BASE_PAD=="Yes"|table1_data$BASE_CKD=="Yes","Yes","No")
table1_data$base_other<-factor(table1_data$base_other,levels = c("No","Yes"))
table1_data$COYOACAN<-factor(table1_data$COYOACAN,levels = c(0,1),labels = c("No","Yes"))
table1_data$diabetes_at_baseline<-factor(table1_data$diabetes_at_baseline,levels = c(0,1),labels = c("No","Yes"))
table1_data$smokegp2<-factor(table1_data$smokegp2,levels = c(1:5),labels = c("Never","Ex-smoker","Current (<daily)","Current (<10/day)","Current (>=10/day)"))
table1_data$EDU_LEVEL<-factor(table1_data$EDU_LEVEL,levels = c(1:4),c("Uni/College","High school","Elementary","Other"),labels = c("University/College","High school","Elementary","Other"))
table1_data$ACTIVITY<-factor(table1_data$ACTIVITY,levels = c("NONE","< ONCE A WEEK","1-2 TIMES A WEEK","AT LEAST 3 TIMES A WEEK"),labels = c("None","< Once a week","1-2 times a week","At least 3 times a week"))
table1_data$ALCGP<-factor(table1_data$ALCGP,levels=c(1:5),c("Never", "Former","up to 3 times a month", "up to 2 times a week", "3+ times a week"))

##additional check for ancestry proportion-------------------------
table1_data$amrscore_cut<-cut(table1_data$amrscore, breaks=c(0,0.6,0.8,1),
                              labels = c("<60%","60%-80%",">=80%"),
                              include.lowest=TRUE,right = F)


variables_to_select<-c("AGE","SEX","COYOACAN","EDU_LEVEL",
                       "ACTIVITY","smokegp2","ALCGP","SBP","DBP","BMI","WHRATIO","HDL_C","LDL_C",
                       "BASE_HBA1C","BASE_CHD","BASE_CVD",
                       "BASE_CANCER","diabetes_at_baseline","base_other")
var_name<-c("Age,years","Sex","District (Coyoacan)","Education level",
            "Leisure-time physical activity","Smoking status","Alcohol intake",
            "SBP(mmHg)","DBP(mmHg)","BMI(kg/m\u00b2)",
            "Waist to Hip Ratio","NMR measured HDL-C (nm)","NMR measured LDL-C(nm)",
            "HbA1c (%)","Coronary Artery Disease",
            "Cardiovascular Disease","Cancer","Diabetes (HBA1C>6.5%, self-report or medication use)","Other chronic diseases")



baseline_ancestry<-create_table_1(data=table1_data,
                                  variables_to_select = variables_to_select,
                                  var_name = var_name,by="amrscore_cut",missing = F,compare=F)


baseline_ancestry<-rbind(baseline_ancestry[c(1:27),],
                         c("Physiological measures, Mean(SD)",rep(NA,4)),
                         baseline_ancestry[c(28:34),],
                         c("Reoprted medical history, N(%)",rep(NA,4)),
                         baseline_ancestry[c(35:39),])

rownames(baseline_ancestry)<-seq(1:nrow(baseline_ancestry))

for (i in c(29:35,37:41)){
  baseline_ancestry[i,1]<-paste('\U{00A0}\U{00A0}',gsub(",.*$", "",baseline_ancestry[i,1]),sep="")
}


flex_baseline<-my_table(baseline_ancestry, y=c(1:2,5,6,11,16,22,28,36),x=1:5)

flex_baseline <- add_header_row(flex_baseline, values =  c("Characteristic", "Indigenous American ancestry proportion "),
                                colwidths = c(1, 4), top = T)

flex_baseline<-add_footer_lines(flex_baseline,"Abbreviation \n SBP:systolic blood pressure, DBP: diastolic blood pressure, BMI: body mass index, NMR:nuclear magnetic resonance,HDL-C: high density lipoprotein cholesterol,LDL-C: low density lipoprotein cholesterol, HbA1c: glycosylated hemoglobin A1c")

flex_baseline<-footnote(flex_baseline,i=41,j=1,part="body",ref_symbols="§",
                        value = as_paragraph("Other diseases include self-reported emphysema, chronic kidney disease, peptic ulcer,liver cirrhosis,and peripheral arterial disease."))

flex_baseline<-set_caption(
  flex_baseline,
  caption = "Baseline Characteristics by ancestry proportion",align_with_table = FALSE,
  word_stylename = "Table Caption",
  fp_p = fpp
)

flex_baseline<-my_theme(flex_baseline,y=c(1:41))
saveRDS(flex_baseline,paste(data_path,"/manuscript-age80/supplement/baseline_char_by_ancestry_80.rds",sep=""))
flex_baseline
