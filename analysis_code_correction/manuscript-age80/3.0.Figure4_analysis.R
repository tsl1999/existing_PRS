##################################
# PRS prediction of CAD in MCPS
# Analysis script for manuscript
# Figure 3
# Author: Tianshu Liu
# Date: 03 June 2024
#################################
#call packages and data
##package and self-create function
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
data_decile=data
remove_cat<-which(colnames(data)%in%colnames(data%>%select(contains("_cat"))))

data_d<-data[,c(1:32,remove_cat)]
data<-data[,c(-remove_cat)]

# combine sex into one forest plot -----------------------
partial_adjustments=c("AGE","SEX",paste(rep("PC",7),seq(1,7),sep=""))#,paste(rep("PC",7),seq(1,7),sep="")
full_adjustments = c(partial_adjustments, "WHRATIO","SBP","DBP","EDU_LEVEL","smokegp2","diabetes_at_baseline")#,paste(rep("PC",7),seq(1,7),sep="")
confounder_check<-c("SEX")


## Partial ----------------

all_strata_outcome<-run_stratified(data, strata_in = confounder_check,
                                   outcome="prevalent_CHD_EPA",adjustment = partial_adjustments,
                                   se=T,show_int = T)
levels(data$SEX)
model_output_partialm<-all_strata_outcome[[1]][[1]]
model_output_partialf<-all_strata_outcome[[1]][[2]]
model_output_partial_int<-all_strata_outcome[[1]][[3]]
model_output_heter_partial<-model_output_partial_int$int_pval
names(model_output_heter_partial)<-rownames(model_output_partial_int)


##Full ----------------

all_strata_outcome<-run_stratified(data, strata_in = confounder_check,outcome="prevalent_CHD_EPA",
                                   adjustment = full_adjustments,se=T,show_int = T)
levels(data$SEX)
model_output_fullm<-all_strata_outcome[[1]][[1]]
model_output_fullf<-all_strata_outcome[[1]][[2]]
model_output_full_int<-all_strata_outcome[[1]][[3]]
model_output_heter_full<-model_output_full_int$int_pval
names(model_output_heter_full)<-rownames(model_output_full_int)

#plot forest plot-----------
tm <- forest_theme(base_size = 7,
                   core = list(bg_params=list(fill = c("white"))),
                   summary_col = "black",
                   refline_lty = "solid",
                   ci_pch = c(15, 16),
                   ci_col = c("skyblue", "tomato"),
                   footnote_col = "black",
                   footnote_cex = 0.9,
                   legend_name = "Sex",
                   legend_value = c("Men", "Women"),
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",xaxis_cex=1,
                   ci_lty =1 ,
                   ci_lwd = 1,
                   ci_Theight = 0)


forest_tablem<-generate_forest_table(
  tables=list(model_output_partialm,model_output_fullm),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)
colnames(forest_tablem)[4:6]<-paste(colnames(forest_tablem)[4:6],"_m",sep="")
colnames(forest_tablem)[7]<-"OR per SD (95%CI)"
forest_tablem$"CAD case/control"<-paste(forest_tablem$`CAD cases`,"/",forest_tablem$`CAD controls`,sep="")
forest_tablem$"CAD case/control"[seq(1,24,3)]<-""

forest_tablef<-generate_forest_table(
  tables=list(model_output_partialf,model_output_fullf),
  model_name = c("Partially adjusted","Fully adjusted"),or_hr = "OR",show_pgsname = T)
forest_tablef$"CAD case/control"<-paste(forest_tablef$`CAD cases`,"/",forest_tablef$`CAD controls`,sep="")
forest_tablef$"CAD case/control"[seq(1,24,3)]<-""
colnames(forest_tablef)[7]<-"OR per SD (95%CI) "
colnames(forest_tablef)[4:6]<-paste(colnames(forest_tablef)[4:6],"_f",sep="")
colnames(forest_tablef)[c(2,3,7:11)]<-paste(colnames(forest_tablef)[c(2,3,7:11)]," ",sep="")

forest_table_sex<-cbind(forest_tablem[,c(1,11,4:8,10)],forest_tablef[,c(11,4:8,10)])

heterogenity<-data.frame(rbind(model_output_heter_partial,model_output_heter_full))
rownames(heterogenity)<-c("Partial","Full")

n_table=2

for(i in 1:npgs){
  pgs<-colnames(model_table_top)[i+1]
  forest_table_sex$heterogenity_test[nrow(forest_table_sex)/npgs*i-2]<-" "
  for (j in 1:n_table){
    
    heterogenity_p<-heterogenity[j,colnames(model_table_top)[i+1]==colnames(heterogenity)]
    forest_table_sex$heterogenity_test[nrow(forest_table_sex)/npgs*i-2+j]<-heterogenity_p
  }
}

colnames(forest_table_sex)[ncol(forest_table_sex)]<-"Interaction p-value "

write.csv(forest_table_sex,paste(data_save_path,"/primary_sex_stratified.csv",sep=""))
saveRDS(forest_table_sex,paste(data_save_path,"/primary_sex_stratified.rds",sep=""))

p <- forest(forest_table_sex[,c(1,2,6,7,9,13,14,15,16)],
            est = list(as.numeric(forest_table_sex$estimate_m),
                       as.numeric(forest_table_sex$estimate_f)),
            lower = list(as.numeric(forest_table_sex$LR_m),
                         as.numeric(forest_table_sex$LR_f)), 
            upper = list(as.numeric(forest_table_sex$UR_m),
                         as.numeric(forest_table_sex$UR_f)),
            ci_column = c(8),
            ref_line = 1,
            arrow_lab =c("Lower risk of CAD", "Higher risk of CAD") ,
            xlim = c(0.9,1.5),
            ticks_at = c(0.9,1,1.1,1.2,1.3,1.4,1.5),
            footnote = "\n\n\n\n\n\n\n\n\n ",
            theme = tm,
            font.family="Calibri")



p <- add_text(p, text = c("Men"),
              part = "header", 
              col = c(2:4),row=0,
              gp = gpar(fontface="bold",fontsize=8))
p <- add_text(p, text = c("Women"),
              part = "header", 
              col = c(5:7),row=0,
              gp = gpar(fontface="bold",fontsize=8))

p <- add_text(p, text = c("OR comparison between gender"),
              part = "header", 
              col = c(8),row=1,
              gp = gpar(fontface="bold",fontsize=7))



p <- add_border(p, 
                part = "header", 
                row = 1,
                col = 1:9, gp = gpar(lwd = 1))


p <- add_border(p, 
                part = "body", 
                col = c(2:7),row=1:nrow(forest_table_sex),
                gp = gpar(lwd = .5))

p <- add_border(p, where=c("left"),
                part = "body",
                col = c(5),row=1:nrow(forest_table_sex),
                gp = gpar(lwd = .5))

p<- edit_plot(p, row = c(seq(1,nrow(forest_table_sex),by=nrow(forest_table_sex)/8)), col=1,
              gp = gpar(fontface = "bold"))

a<-get_wh(p,unit = "cm")+2

png(paste(graphs_path,"/Figure4_logistic_forestplotEPA_sex",".png",sep=""), res = 300, width = a[1], height = a[2], units = "cm")
print(p)
dev.off()
