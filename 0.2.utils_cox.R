library(gridExtra)

library(forestploter)
pgsname<-colnames(variant_t)[2:8]
data=data
outcome="EPO001_up_75"
adjustments<-"SEX"
#adjustments=c("WHRATIO","SBP","DBP","EDU_LEVEL","smokegp","diabetes_at_baseline")

#library(forestplot)

create_surv_formula<-function(adjustments,outcome,pgs_name){
  surv_formula<-paste("Surv(time_in,time_out,",outcome,")",sep="")
  output_formula=formula(paste(surv_formula,
                                 "~",paste(pgs_name,paste(adjustments,collapse="+"),
                                           sep="+"),sep=""))
output_formula
}

#create_surv_formula(adjustments=adjustments,outcome=outcome,pgs_name=pgs_name)

run_cox<-function(data,adjustments,outcome,pgs_name){
  cat("\n Modelling ",pgs_name)
  model_output<-coxph(create_surv_formula(adjustments,outcome,pgs_name),data=data)
  model_output
}

#run_cox(data=data,adjustments=adjustments,outcome=outcome,pgs_name=pgs_name)

create_output_table_cox<-function(trainsplit=F,data,train_data=NULL,test_data=NULL,adjustments,outcome,namew=NULL){
  
  if(trainsplit==F){
    train_data=data
    test_data=data
  }else if(trainsplit==T){
    train_data=train_data
    test_data=test_data
  }
    train<-na.omit(train_data%>%select(IID,time_in,time_out,all_of(c(outcome,adjustments)),contains(c("PGS","custom"))))
    test<-na.omit(test_data%>%select(IID,time_in,time_out,all_of(c(outcome,adjustments,pgsname)),contains(c("PGS","custom"))))

    start_column=which( colnames(train)%in%colnames(variant_t))
    model_output_table<-data.frame(HR=NA,LR=NA,UR=NA,AUC=NA,AUC_LR=NA,AUC_UR=NA,case_no=NA,control_no=NA,pval=NA)
    
    for (i in 1:7){
     
      model_output<-run_cox(data=train,adjustments,outcome,pgs_name=colnames(train)[start_column[i]])
      
      jpeg(paste(graphs_path,"/cox/",colnames(train)[start_column[i]],"/",colnames(train)[start_column[i]],"_Schoenfeld_residual",namew,".png",sep=""),width = 800, height = 600)
      cox.zph.fit <- cox.zph(model_output)
      p<-ggcoxzph(cox.zph.fit,font.main = 8)
      print(p)
      dev.off()
      
      hr_ci<-round(summary(model_output)$conf.int[1,c(1,3:4)],3)
      
      c_index<-round(ci_manual(x=concordance(model_output)$concordance,se=sqrt(concordance(model_output)$var)),3)
      case_no<-sum(train[,outcome]==1)
      case_id<-unique(train[train[,outcome]==1,"IID"])

      control<-train[!train$IID%in%case_id,]
      control_no<-length(unique(control[control[,outcome]==0,"IID"]))
     pval<-ifelse(summary(model_output)$coefficients[,"Pr(>|z|)"][colnames(train)[start_column[i]]]<0.001,
                  "<0.001",round(summary(model_output)$coefficients[,"Pr(>|z|)"][colnames(train)[start_column[i]]],3))
      
      model_output_table[i,]<-c(hr_ci[1],hr_ci[2],hr_ci[3],
                                c_index[1],c_index[2],c_index[3],case_no,control_no,pval)
      rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
      
    }
    
    model_output_table
    
  
    

    
}
  




model_to_table_cox<-function(model_output_table){
  colength<-ncol(model_output_table)#7,10
  
  if(colength==9){
    model_output_table[,1:8]<-apply(model_output_table[,1:8],2,as.numeric)
    model_output_table[,7]<-paste(round(model_output_table[,1],3),"(",round(model_output_table[,2],3),"-",round(model_output_table[,3],3),")",sep="")
    model_output_table[,8]<-paste(round(model_output_table[,4],3),"(",round(model_output_table[,5],3),"-",round(model_output_table[,6],3),")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(7,8),]
    rownames(model_output_table_t)[1:2]<-c("HR per SD","Concordance index")
    model_output_table_t<-cbind(rownames(model_output_table_t),model_output_table_t)
    colnames(model_output_table_t)[1]<-NA
  }else if (colength==12){
    model_output_table[,1:11]<-apply(model_output_table[,1:8],2,as.numeric)
    model_output_table[,10]<-paste(round(model_output_table[,1],3),"(",model_output_table[,2],"-",model_output_table[,3],")",sep="")
    model_output_table[,11]<-paste(round(model_output_table[,4],3),"(",model_output_table[,5],"-",model_output_table[,6],")",sep="")
    model_output_table[,12]<-paste(round(model_output_table[,7],3),"(",round(model_output_table[,8],3),"-",round(model_output_table[,9],3),")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(10:12),]
    colnames(model_output_table_t)<-sub("_cat","",colnames(model_output_table_t))
    rownames(model_output_table_t)[1:3]<-c("HR (intermediate vs Low)","HR (high vs Low)","Concordance index")
    model_output_table_t<-cbind(rownames(model_output_table_t),model_output_table_t)
    colnames(model_output_table_t)[1]<-NA
  }
  model_output_table_t
}











plot_martingale<-function(variable,outcome,include="all",save=NULL){

data_in<-data%>%select(IID,variable,time_in,time_out,outcome)
data_in<-na.omit(data_in)
if(include=="all"){
formula_in<-formula(paste("Surv(time_in,time_out, ",outcome,")~"
                       ,paste(variable,"+log(",variable,")+ I(",variable,"^2",")+sqrt(",variable,")",sep=""),sep=""))
}else if(include=="linear") {
  formula_in<-formula(paste("Surv(time_in,time_out, ",outcome,")~",variable,sep=""))
}
model<-coxph(formula_in, data = data_in)
jpeg(paste(graphs_path,save,"/Martingale_residual_",variable,".png",sep=""),width = 800, height = 600)
p<-ggcoxfunctional(model, data = data_in, point.col = "blue", point.alpha = 0.5,
                title = paste("Martingale residules for",variable,sep=" "))
print(p)
dev.off()}


library(grid)
tm2 <- forest_theme(base_size = 5,
                    core = list(bg_params=list(fill = c("white"))),
                    ci_pch = c(15, 16),
                    ci_col = c("black", "navy"),
                    footnote_col = "black",
                    footnote_cex = 0.9,
                    legend_name = "Model type",
                    legend_value = c("Partial", "Full"),
                    vertline_lwd = 1,
                    vertline_lty = "dashed",
                    vertline_col = "grey20",xaxis_cex=0.9,
                    ci_lty =1 ,
                    ci_lwd = 0.5,
                    ci_Theight = 0.1,)
