
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
library(PredictABEL)

ci<-function(x,number=95,roundup=3){
  actualci<-(100-(100-number)/2)/100
  ci_bound<-c(-qnorm(actualci,0,1),qnorm(actualci,0,1))
  coef<-summary(x)$coefficients[2,1]
  se<-summary(x)$coefficients[2,'Std. Error']
  ci<-exp(coef+ci_bound*se)
  c(round(exp(coef),roundup),round(ci,roundup))
}

ci_cat<-function(x,y,number=95,roundup=3){
  actualci<-(100-(100-number)/2)/100
  ci_bound<-c(-qnorm(actualci,0,1),qnorm(actualci,0,1))
  coef<-summary(x)$coefficients[y,1]
  se<-summary(x)$coefficients[y,'Std. Error']
  ci_lr<-exp(coef+ci_bound[1]*se)
  ci_ur<-exp(coef+ci_bound[2]*se)
  c(round(exp(coef),roundup),round(ci_lr,roundup),round(ci_ur,roundup))
}



create_glm_formula<-function(model_type,data,full_adjustments,outcome,pgs_name,partial_adjustment){
  
  
  if (model_type=="simple"){
    output_formula=formula(paste(outcome,"~",pgs_name,sep=""))
  }else if (model_type=="partial"){
    output_formula=formula(paste(outcome,"~",paste(pgs_name,paste(partial_adjustment,collapse = "+"),
                                                   paste(PCs,collapse="+"),sep="+"),sep=""))
  }else if (model_type=="full"){
    output_formula=formula(paste(outcome,"~",paste(pgs_name,paste(partial_adjustment,collapse = "+"),
                                                   paste(full_adjustments,collapse = "+"),
                                                   paste(PCs,collapse="+"),sep="+"),sep=""))
  }else (print("ERROR: NO SUCH MODEL"))
  
  output_formula
  
}

run_glm<-function(model_type,data,full_adjustments,outcome,pgs_name,partial_adjustment){
  cat("\n Modelling ",pgs_name)
  model_output<-glm(create_glm_formula(model_type,data,full_adjustments,outcome,pgs_name,partial_adjustment),data=data,family=binomial(link='logit'))
  model_output
}



create_output_table<-function(trainsplit=F,data,train_data=NULL,test_data=NULL,model_type,cat_or_cont="cont", full_adjustments,outcome,partial_adjustment,name){
  
  if(trainsplit==F){
    train_data=data
    test_data=data
  }else if(trainsplit==T){
    train_data=train_data
    test_data=test_data
  }
  if (model_type=="simple"){
    train=train_data
    test=test_data
  }else if(model_type=="partial"){
    train<-na.omit(train_data%>%select(IID,partial_adjustment,contains(c("PC","PGS","custom",outcome))))
    test<-na.omit(test_data%>%select(IID,partial_adjustment,contains(c("PC","PGS","custom",outcome))))
  }else if (model_type=="full"){
    train<-na.omit(train_data%>%select(IID,partial_adjustment,contains(c("PC","PGS","custom",outcome,full_adjustments))))
    test<-na.omit(train_data%>%select(IID,partial_adjustment,contains(c("PC","PGS","custom",outcome,full_adjustments))))
  }
  
  
  if (cat_or_cont=="cont"){
    start_column=which( colnames(train)%in%colnames(variant_t))
    model_output_table<-data.frame(OR=NA,LR=NA,UR=NA,AUC=NA,AUC_LR=NA,AUC_UR=NA,Nagelkerke_PseudoR2=NA,lee_PseudoR2=NA,case_no=NA,control_no=NA,pval=NA)
    jpeg(paste(graphs_path,"/AUCplot_",model_type,"_",outcome,"_",cat_or_cont,namew,".jpg",sep=""),width = 800, height = 600)
    
    for (i in 1:7){
      
      model_output<-run_glm(model_type,data=train,full_adjustments,outcome,pgs_name=colnames(train)[start_column[i]],partial_adjustment)
      pred<-predict(model_output,test,type='response')
      pred.obj<-prediction(pred,test[,outcome])
      perf.obj<-performance(pred.obj,"tpr","fpr")
      roc_prs<-roc(test[,outcome],pred)
      if(i==1){
        plot(perf.obj,col=i)}else{plot(perf.obj,col=i,add=T)}
      
      vr=runif(nrow(train),0,1)
      vsel=model_output$linear.predictors[model_output$y==0|vr<0.05]
      r2<-var(vsel)/(var(vsel)+pi^2/3)
      
      
      auc_point<-auc(roc_prs)
      auc_ci<-ci.auc(test[,outcome],pred)
      
      pval<-ifelse(summary(model_output)$coefficients[,"Pr(>|z|)"][2]<0.001,"<0.001",round(summary(model_output)$coefficients[,"Pr(>|z|)"][2],3))
      model_output_table[i,]<-c(exp(coef(model_output)[2]),ci(model_output)[2],
                                ci(model_output)[3],auc_point,auc_ci[1],auc_ci[3],
                                PseudoR2(model_output,c("Nagelkerke")),r2,
                                table(train[,outcome])["1"], table(train[,outcome])["0"],
                                pval)
      rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
      
    }
    abline(a=0,b=1,lty="dashed",col="gray")
    legend("bottomright",legend=c(paste(colnames(train)[start_column[1:length(start_column)]]," AUC=",round(as.numeric(model_output_table[,4]),3))),
           col=c(palette()[1:length(start_column)]),lty=1,cex=0.8)
    dev.off()
    model_output_table
    
  }else if (cat_or_cont=="cat"){
    jpeg(paste(graphs_path,"/AUCplot_",model_type,"_",outcome,"_",cat_or_cont,name,".jpg",sep=""),width = 800, height = 600)
    start_column=which( colnames(train)%in%paste(colnames(variant_t),"_cat",sep=""))
    model_output_table<-data.frame(intermediate_prs=NA,intermediate_lr=NA,intermediate_ur=NA,
                                   high_prs=NA,high_lr=NA,high_ur=NA,AUC=NA,AUC_LR=NA,AUC_UR=NA,
                                   Nagelkerke_PseudoR2=NA,lee_PseudoR2=NA,case_no=NA,control_no=NA,case_int=NA,case_high=NA,control_int=NA,control_high=NA,pval_int=NA,pval_high=NA)
    for (i in 1:7){
      model_output<-run_glm(model_type,data=train,full_adjustments,outcome,pgs_name=colnames(train)[start_column[i]],partial_adjustment)
      pred<-predict(model_output,test,type='response')
      pred.obj<-prediction(pred,test[,outcome])
      perf.obj<-performance(pred.obj,"tpr","fpr")
      if(i==1){
        plot(perf.obj,col=i)}else{plot(perf.obj,col=i,add=T)}
      roc_prs<-roc(test[,outcome],pred)
      auc_point<-auc(roc_prs)
      auc_ci<-ci.auc(test[,outcome],pred)
      vr=runif(nrow(train),0,1)
      vsel=model_output$linear.predictors[model_output$y==0|vr<0.05]
      r2<-var(vsel)/(var(vsel)+pi^2/3)
      
      pval_int<-ifelse(summary(model_output)$coefficients[,"Pr(>|z|)"][2]<0.001,"<0.001",round(summary(model_output)$coefficients[,"Pr(>|z|)"][2],3))
      pval_high<-ifelse(summary(model_output)$coefficients[,"Pr(>|z|)"][3]<0.001,"<0.001",round(summary(model_output)$coefficients[,"Pr(>|z|)"][3],3))
      pgs_name=colnames(train)[start_column[i]]
      model_output_table[i,]<-c(exp(coef(model_output)[2]),ci_cat(model_output,y=2)[2:3],
                                exp(coef(model_output)[3]),ci_cat(model_output,y=3)[2:3],auc_point,auc_ci[1],
                                auc_ci[3],PseudoR2(model_output,c("Nagelkerke")),r2,
                                table(train[,outcome])["1"], table(train[,outcome])["0"],
                                
                                table(train[train[,pgs_name]=="intermediate",outcome])["1"],
                                table(train[train[,pgs_name]=="high",outcome])["1"],
                                table(train[train[,pgs_name]=="intermediate",outcome])["0"],
                                table(train[train[,pgs_name]=="high",outcome])["0"],
                                pval_int,pval_high)
      rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
    }
    abline(a=0,b=1,lty="dashed",col="gray")
    legend("bottomright",legend=c(paste(colnames(train)[start_column[1:length(start_column)]]," AUC=",round(as.numeric(model_output_table[,7]),3))),
           col=c(palette()[1:length(start_column)]),lty=1,cex=0.8)
    dev.off()
    model_output_table
    
  }
  
  
}


model_to_table<-function(model_output_table){
  colength<-ncol(model_output_table)#7,10
  model_output_table_up<-apply(model_output_table,2,as.numeric)
  rownames(model_output_table_up)<-rownames(model_output_table)
  model_output_table<-data.frame(model_output_table_up)
  if(colength==8){
    model_output_table[,9]<-paste(round(model_output_table[,1],3),"(",round(model_output_table[,2],3),"-",round(model_output_table[,3],3),")",sep="")
    model_output_table[,10]<-paste(round(model_output_table[,4],3),"(",round(model_output_table[,5],3),"-",round(model_output_table[,6],3),")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(9,10,7,8),]
    rownames(model_output_table_t)[1:2]<-c("OR per SD","AUC")
    model_output_table_t<-cbind(rownames(model_output_table_t),model_output_table_t)
    colnames(model_output_table_t)[1]<-NA
  }else if (colength==11){
    model_output_table[,12]<-paste(round(model_output_table[,1],3),"(",model_output_table[,2],"-",model_output_table[,3],")",sep="")
    model_output_table[,13]<-paste(round(model_output_table[,4],3),"(",model_output_table[,5],"-",model_output_table[,6],")",sep="")
    model_output_table[,14]<-paste(round(model_output_table[,7],3),"(",round(model_output_table[,8],3),"-",round(model_output_table[,9],3),")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(12:14,10,11),]
    colnames(model_output_table_t)<-sub("_cat","",colnames(model_output_table_t))
    rownames(model_output_table_t)[1:3]<-c("OR (intermediate vs Low)","OR (high vs Low)","AUC")
    model_output_table_t<-cbind(rownames(model_output_table_t),model_output_table_t)
    colnames(model_output_table_t)[1]<-NA
  }
  model_output_table_t
}

output_table_summary<-function(outcome,data,full_adjustments,trainsplit=F,train_data=NULL,test_data=NULL,partial_adjustment=c("SEX","AGE"),namew=NULL,fp_layout=1){
  
  
  #continuous
  cat("\n Continuous PRS ")
  cat("\n Unadjusted model ")
  model_output_simple<-create_output_table(trainsplit ,data=data,train_data,test_data,model_type = "simple",cat_or_cont = "cont",full_adjustments = full_adjustments,outcome=outcome,partial_adjustment=partial_adjustment,namew)
  cat("\n Age and Sex Adjusted model (+7PCs) ")
  model_output_age_sex<-create_output_table(trainsplit, data,train_data,test_data,model_type = "partial",cat_or_cont = "cont",full_adjustments = full_adjustments,outcome=outcome,partial_adjustment=partial_adjustment,namew)
  cat("\n Fully Adjusted Model ")
  model_output_full<-create_output_table(trainsplit , data,train_data,test_data,model_type = "full",cat_or_cont = "cont",full_adjustments = full_adjustments,outcome=outcome,partial_adjustment=partial_adjustment,namew)
  
  
  
  model_list<-list(model_output_simple,model_output_age_sex,model_output_full)
  #create forest plot table
  forest_table<-data.frame(subgroup=NA,case_no=NA,control_no=NA,OR=NA,LR=NA,UR=NA,hr_ci=NA,pval=NA)
  
  if(fp_layout==1){
    hr_ci<-NA
    name<-c("\U{00A0}\U{00A0}Simple","\U{00A0}\U{00A0}Partially adjusted","\U{00A0}\U{00A0}Fully adjusted")
    for (i in 1:7){
      forest_table<-rbind(forest_table,c(rownames(model_output_simple)[i],rep(" ",7)))
      for (j in 1:3){
        hr_ci<-paste(round(as.numeric(model_list[[j]][i,1]),3),"(",model_list[[j]][i,2],"-",model_list[[j]][i,3],")",sep="")
        row_input<-unlist(c(name[j],model_list[[j]][i,c("case_no","control_no","OR","LR","UR")],
                            hr_ci,model_list[[j]][i,c("pval")]))
        forest_table<-rbind(forest_table,row_input,make.row.names=F)
      }
      
    }
    forest_table<-forest_table[2:29,]
    
  }else if(fp_layout==2){
    namem<-c("Simple","Partially adjusted","Fully adjusted")
    hr_ci<-NA
    for(i in 1:3){
      forest_table<-rbind(forest_table,c(namem[i],rep(" ",7)))
      hr_ci<-paste(round(as.numeric(model_list[[i]][,1]),3),"(",model_list[[i]][,2],"-",model_list[[i]][,3],")",sep="")
      for(j in 1:7){
        row_input<-unlist(c(paste("\U{00A0}\U{00A0}",rownames(model_list[[i]])[j]),model_list[[i]][j,c("case_no","control_no","OR","LR","UR")],
                            hr_ci[j],model_list[[i]][j,c("pval")]))
        forest_table<-rbind(forest_table,row_input,make.row.names=F)
      }
      
    }
    forest_table<-forest_table[2:25,]
    
  }
  
  
  
  forest_table$" "<- paste(rep(" ", 40), collapse = " ")
  colnames(forest_table)<-c("Models","CAD case","No CAD","OR","LR","UR","OR per SD(95%CI)","p-value","")
  library(forestploter)
  p <- forest(forest_table[,c(1:3, 7,9,8)],
              est = as.numeric(forest_table$OR),
              lower = as.numeric(forest_table$LR), 
              upper = as.numeric(forest_table$UR),
              ci_column = 5,
              ref_line = 1,
              arrow_lab = c("Lower risk of CAD", "Higher risk of CAD"),
              xlim = c(0.5, 1.3),
              ticks_at = c(0.5,0.8, 1, 1.3),
              footnote = "\n\n\n\n\n\n\n\n\n Simple model has no adjustments. Partially adjusted model adjusts for age, sex and 7 PCs. \n Fully adjusted model additionally adjusts for education level, smoking status, blood pressure and waist-hip ratio ",
              theme=tm(15))
  
  
  
  jpeg(paste(graphs_path,"/Forestplot_cont","_outcome",namew,".jpg",sep=""),width = 1000, height = 800)
  # Print plot
  plot(p)
  
  dev.off()
  
  model_output_simple<-model_output_simple[,1:8]
  model_output_age_sex<-model_output_age_sex[,1:8]
  model_output_full<-model_output_full[,1:8]
  model_output_simple_t<-model_to_table(model_output_simple)
  model_output_age_sex_t<-model_to_table(model_output_age_sex)
  model_output_full_t<-model_to_table(model_output_full)
  
  
  #categorical
  cat("\n Categorical PRS ")
  cat("\n Unadjusted model ")
  model_output_simple_cat<-create_output_table(trainsplit , data,train_data,test_data,model_type = "simple",cat_or_cont = "cat",full_adjustments = full_adjustments,outcome=outcome,partial_adjustment=partial_adjustment,namew)
  cat("\n Age and Sex Adjusted model (+7PCs) ")
  model_output_age_sex_cat<-create_output_table(trainsplit , data,train_data,test_data,model_type = "partial",cat_or_cont = "cat",full_adjustments = full_adjustments,outcome=outcome,partial_adjustment=partial_adjustment,namew)
  cat("\n Fully Adjusted Model")
  model_output_full_cat<-create_output_table(trainsplit , data,train_data,test_data,model_type = "full",cat_or_cont = "cat",full_adjustments = full_adjustments,outcome=outcome,partial_adjustment=partial_adjustment,namew)
  
  model_list<-list(model_output_simple_cat,model_output_age_sex_cat,model_output_full_cat)
  #create forest plot table
  forest_table<-data.frame(subgroup=NA,category=NA,case_no=NA,control_no=NA,OR=NA,LR=NA,UR=NA,hr_ci=NA,pval=NA)
  
  
  if (fp_layout==1){
    
    namem<-c("\U{00A0}\U{00A0}Simple","\U{00A0}\U{00A0}Partially adjusted","\U{00A0}\U{00A0}Fully adjusted")
    for (i in 1:7){
      forest_table<-rbind(forest_table,c(rownames(model_output_simple)[i],rep(" ",8)))
      for (j in 1:3){
        hr_ci_int<-paste(round(as.numeric(model_list[[j]][i,1]),3),"(",model_list[[j]][i,2],"-",model_list[[j]][i,3],")",sep="")
        hr_ci_high<-paste(round(as.numeric(model_list[[j]][i,4]),3),"(",model_list[[j]][i,5],"-",model_list[[j]][i,6],")",sep="")
        row_input_int<-unlist(c(" ","intermediate",model_list[[j]][i,c("case_int","control_int","intermediate_prs","intermediate_lr","intermediate_ur")],
                                hr_ci_int,model_list[[j]][i,c("pval_int")]))
        row_input_high<-unlist(c(" ","high",model_list[[j]][i,c("case_high","control_high","high_prs","high_lr","high_ur")],
                                 hr_ci_high,model_list[[j]][i,c("pval_high")]))
        prs_input<-unlist(c(namem[j]," ",model_list[[j]][i,c("case_no","control_no")],rep(NA,5)))
        forest_table<-rbind(forest_table,prs_input,row_input_int,row_input_high,make.row.names=F)
      }
      
    }
    forest_table<-forest_table[2:71,]
  }else if (fp_layout==2){
    
    namem<-c("Simple","Partially adjusted","Fully adjusted")
    
    for(i in 1:3){
      forest_table<-rbind(forest_table,c(namem[i],rep(" ",8)))
      hr_ci_int<-paste(round(as.numeric(model_list[[i]][,1]),3),"(",model_list[[i]][,2],"-",model_list[[i]][,3],")",sep="")
      hr_ci_high<-paste(round(as.numeric(model_list[[i]][,4]),3),"(",model_list[[i]][,5],"-",model_list[[i]][,6],")",sep="")
      
      for(j in 1:7){
        row_input_int<-unlist(c(" ","Intermediate",model_list[[i]][j,c("case_int","control_int","intermediate_prs","intermediate_lr","intermediate_ur")],
                                hr_ci_int[j],model_list[[i]][j,c("pval_int")]))
        row_input_high<-unlist(c(" ","High",model_list[[i]][j,c("case_high","control_high","high_prs","high_lr","high_ur")],
                                 hr_ci_high[j],model_list[[i]][j,c("pval_high")]))
        prs_input<-unlist(c(paste("\U{00A0}\U{00A0}",rownames(model_list[[i]])[j])," ",model_list[[i]][j,c("case_no","control_no")],rep(" ",5)))
        
        forest_table<-rbind(forest_table,prs_input,row_input_int,row_input_high,make.row.names=F)}
      
    }}
  forest_table<-forest_table[2:67,]
  
  
  forest_table$" "<- paste(rep(" ", 40), collapse = " ")
  colnames(forest_table)<-c("Models"," ","CAD case","No CAD","OR","LR","UR","OR per SD(95%CI)","p-value","")
  library(forestploter)
  p <- forest(forest_table[1:35,c(1:4,8,10,9)],
              est = as.numeric(forest_table$OR),
              lower = as.numeric(forest_table$LR), 
              upper = as.numeric(forest_table$UR),
              ci_column = 6,
              ref_line = 1,
              arrow_lab = c("Lower risk of CAD", "Higher risk of CAD"),
              xlim = c(0.5, 2),
              ticks_at = c(0.5,0.8, 1, 1.5,1.8,2),
              mar = unit(rep(1, times = 4), "mm"),
              footnote = "\n\n\n\n\n\n\n\n\n Simple model has no adjustments. Partially adjusted model adjusts for age, sex and 7 PCs. \n Fully adjusted model additionally adjusts for education level, smoking status, blood pressure and waist-hip ratio "
              ,theme=tm(10))
  
  p2 <- forest(forest_table[35:nrow(forest_table),c(1:4,8,10,9)],
               est = as.numeric(forest_table$OR),
               lower = as.numeric(forest_table$LR), 
               upper = as.numeric(forest_table$UR),
               ci_column = 6,
               ref_line = 1,
               arrow_lab = c("Lower risk of CAD", "Higher risk of CAD"),
               xlim = c(0.5, 2),
               ticks_at = c(0.5,0.8, 1, 1.5,1.8,2),
               mar = unit(rep(1, times = 4), "mm"),
               footnote = "\n\n\n\n\n\n\n\n\n Simple model has no adjustments. Partially adjusted model adjusts for age, sex and 7 PCs. \n Fully adjusted model additionally adjusts for education level, smoking status, blood pressure and waist-hip ratio "
               ,theme=tm(12))
  
  
  jpeg(paste(graphs_path,"/Forestplot_cat1","_outcome",namew,".jpg",sep=""),width = 1000, height = 800)
  # Print plot
  plot(p)
  dev.off()
  jpeg(paste(graphs_path,"/Forestplot_cat2","_outcome",namew,".jpg",sep=""),width = 1000, height = 800)
  # Print plot
  plot(p2)
  dev.off()
  
  
  
  
  model_output_simple_cat<-model_output_simple_cat[,1:11]
  model_output_age_sex_cat<-model_output_age_sex_cat[,1:11]
  model_output_full_cat<-model_output_full_cat[,1:11]
  
  model_output_simple_cat_t<-model_to_table(model_output_simple_cat)
  model_output_age_sex_cat_t<-model_to_table(model_output_age_sex_cat)
  model_output_full_cat_t<-model_to_table(model_output_full_cat)
  
  
  prs_table_cont<-rbind(variant_t,c("Simple Model", rep(NA,7)),
                        model_output_simple_t,
                        
                        c("Partially adjusted Model", rep(NA,7)),
                        model_output_age_sex_t,
                        
                        c("Fully Adjusted Model", rep(NA,7)),
                        model_output_full_t
  )
  rownames(prs_table_cont)<-seq(1:nrow(prs_table_cont))
  colnames(prs_table_cont)[1]<-"95% CI"
  prs_table_cat<-rbind(variant_t,c("Simple Model", rep(NA,7)),
                       
                       model_output_simple_cat_t,
                       c("Adjusted Model", rep(NA,7)),
                       
                       model_output_age_sex_cat_t,
                       c("Fully Adjusted Model", rep(NA,7)),
                       
                       model_output_full_cat_t)
  rownames(prs_table_cat)<-seq(1:nrow(prs_table_cat))
  colnames(prs_table_cat)[1]<-"95% CI"
  prs_table<-list(prs_table_cont,prs_table_cat)
  prs_table
  
  prs_table
}




discrimination_without_prs<-function(train_data,test_data,full_adjustments,outcome,partial_adjustment=c("SEX","AGE")){
  #partial#
  train_partial<-na.omit(train_data%>%select(IID,partial_adjustment,contains(c("PC","PGS","custom",outcome))))
  test_partial<-na.omit(test_data%>%select(IID,partial_adjustment,contains(c("PC","PGS","custom",outcome))))
  model_partial<-run_glm(model_type="partial",data=train_partial,full_adjustments,outcome,pgs_name=NULL,partial_adjustment)
  #full#
  train_full<-na.omit(train_data%>%select(IID,partial_adjustment,contains(c("PC","PGS","custom",outcome,full_adjustments))))
  test_full<-na.omit(train_data%>%select(IID,partial_adjustment,contains(c("PC","PGS","custom",outcome,full_adjustments))))
  model_full<-run_glm(model_type="full",data=train_full,full_adjustments,outcome,pgs_name=NULL,partial_adjustment)
  
  pred<-predict(model_partial,test_partial,type='response')
  roc_prs<-roc(test_partial[,outcome],pred)
  auc_point<-auc(roc_prs)
  auc_ci<-ci.auc(test_partial[,outcome],pred)
  model_output<-data.frame(AUC=NA,AUC_LR=NA,AUC_UR=NA,pseudoR=NA)
  model_output[1,]<-c(auc_point,auc_ci[1],auc_ci[3],PseudoR2(model_partial))
  rownames(model_output)[1]<-"partial"
  
  pred<-predict(model_full,test_full,type='response')
  roc_prs<-roc(test_full[,outcome],pred)
  auc_point<-auc(roc_prs)
  auc_ci<-ci.auc(test_full[,outcome],pred)
  model_output[2,]<-c(auc_point,auc_ci[1],auc_ci[3],PseudoR2(model_full))
  rownames(model_output)[2]<-"full"
  model_output
  model_output[,5]<-paste(round(model_output[,1],3),"(",round(model_output[,2],3),"-",round(model_output[,3],3),")",sep="")
  colnames(model_output)[5]<-"AUC(95%CI)"
  model_output
}


create_quntile_char<-function(data_prs, pgs,variables_to_select,var_name){
  data_keep<-data_prs[,c(variables_to_select,paste(pgs,"_quintile",sep=""))]
  data_cha<-NA
  for (i in 1:length(variables_to_select)){
    data_cha[i]<-class(data_keep[,i])
  }
  factor_var<-data_keep[,c(1,which(data_cha=="factor"),ncol(data_keep))]
  factor_var_name<-c(var_name[which(data_cha=="factor")])
  cat_table<-NA
  for (i in 2:(length(factor_var)-1)){
    cat_table<-rbind(cat_table,cate_table(x=factor_var[,i],data=factor_var[,paste(pgs,"_quintile",sep="")],name=factor_var_name[i-1],use_preset_level = T))
    
  }
  colnames(cat_table)<-c("name","Q1","Q2","Q3","Q4","Q5","Overall")
  cat_table<-cat_table[2:nrow(cat_table),]
  numeric_var<-data_keep[,c(1,which(data_cha=="numeric"|data_cha=="integer"),ncol(data_keep))]
  numeric_var_name<-var_name[which(data_cha=="numeric"|data_cha=="integer")]
  cont_table<-data.frame(name=NA,Q1=NA,Q2=NA,Q3=NA,Q4=NA,Q5=NA,Overall=NA)
  for (i in 2:(length(numeric_var)-1)){
    cont_table[i-1,]<-get_msd_multi(x=numeric_var[,i],y=numeric_var_name[i-1],data=numeric_var[,paste(pgs,"_quintile",sep="")])
    
  }
  
  
  char_table=rbind(cont_table,cat_table)
  pgs_quint<-paste(pgs,"_quintile",sep="")
  colnames(char_table)<-c(pgs,rep(NA,6))
  for (i in 2:6){
    colnames(char_table)[i]<-paste("Q",i-1,"(N=",sum(data_keep[,pgs_quint]==i-1),")",sep="")
  }
  colnames(char_table)[7]<-paste("Overall (N=",nrow(data_keep),")",sep="")
  
  char_table
  
}


create_table_1<-function(data,variables_to_select,var_name,by,compare=F,missing=T,self_define=NULL){
  data_keep<-data[,c(variables_to_select,by)]
  data_cha<-NA
  compare_level<-levels(data_keep[,by])
  level_length<-length(compare_level)
  for (i in 1:length(variables_to_select)){
    data_cha[i]<-class(data_keep[,i])
  }
  factor_var<-data_keep[,c(which(data_cha=="factor"))]
  factor_var_name<-c(var_name[which(data_cha=="factor")])
  cat_table<-NA
  numeric_var<-data_keep[,c(which(data_cha=="numeric"|data_cha=="integer"))]
  numeric_var_name<-var_name[which(data_cha=="numeric"|data_cha=="integer")]
  cont_table<-NA
  
  if (compare==F){
    for (i in 1:(length(factor_var))){
      if(missing==T){
        cat_table<-rbind(cat_table,cate_table(x=factor_var[,i],data=data_keep[,by],name=paste(factor_var_name[i],", N(%)",sep=""),use_preset_level = T))
      }else{
        cat_table<-rbind(cat_table,cate_table(x=factor_var[,i],data=data_keep[,by],name=paste(factor_var_name[i],", N(%)",sep=""),use_preset_level = T,missing = F))
      }
      
      
    }
    colnames(cat_table)<-c(" ",levels(data_keep[,by]),"Overall")
    for (i in 1:(length(numeric_var))){
      if(sum(is.na(numeric_var[,i]))==0|missing==F){
        cont_table<-rbind(cont_table,get_msd_multi(x=numeric_var[,i],y=paste(numeric_var_name[i],", Mean(SD)",sep=""),z=c(60,10),data=data_keep[,by],self_define=self_define))
      }else if (sum(is.na(numeric_var[,i]))!=0&missing==T){
        x_update<-ifelse(is.na(numeric_var[,i])==T,"Missing",numeric_var[,i])
        cont_table<-rbind(cont_table,get_msd_multi(x=numeric_var[,i],y=paste(numeric_var_name[i],", Mean(SD)",sep=""),z=c(60,10),data=data_keep[,by],self_define = self_define),get_prop_multi(x_update,"Missing",data_keep[,by]))
      }
    }
    
    colnames(cont_table)<-c(" ",levels(data_keep[,by]),"Overall")
    
  }else if (compare==T){
    for (i in 1:(length(factor_var))){
      chi_test<-chisq.test(factor_var[,i][data_keep[,by]%in%c(compare_level[1],compare_level[level_length])],data_keep[,by][data_keep[,by]%in%c(compare_level[1],compare_level[level_length])])
      cat_table_i<-cate_table(x=factor_var[,i],data=data_keep[,by],
                              name=factor_var_name[i],use_preset_level = T,missing = missing)
      if(is.null(nrow(cat_table_i))==T){
        cat_table_i<-c(cat_table_i,round(chi_test$p.value,2))
      }else{cat_table_i<-cbind(cat_table_i,c(round(chi_test$p.value,2),rep(NA,nrow(cat_table_i)-1)))
      colnames(cat_table_i)<-NA
      }
      
      
      cat_table<-rbind(cat_table,cat_table_i)
      
    }
    colnames(cat_table)<-c(" ",levels(data_keep[,by]),"Overall",paste(compare_level[level_length],"VS",compare_level[1],", p-value",sep=" "))
    
    for (i in 1:(length(numeric_var))){
      t_test<-t.test(numeric_var[,i][data_keep[,by]%in%c(compare_level[1],compare_level[level_length])]~data_keep[,by][data_keep[,by]%in%c(compare_level[1],compare_level[level_length])])
      #t_test_conf<-paste(round(t_test$estimate[1]-t_test$estimate[2],3),"(",round(t_test$conf.int,3)[1],",",round(t_test$conf.int,3)[2],")",sep="")
      
      if(sum(is.na(numeric_var[,i]))==0|missing==F){
        cont_table<-rbind(cont_table,c(get_msd_multi(x=numeric_var[,i],y=paste(numeric_var_name[i],", Mean(SD)",sep=""),z=c(60,10),data=data_keep[,by],self_define = self_define),round(t_test$p.value,2)))
      }else if (sum(is.na(numeric_var[,i]))!=0&missing==T){
        x_update<-ifelse(is.na(numeric_var[,i])==T,"Missing",numeric_var[,i])
        cont_table<-rbind(cont_table,c(get_msd_multi(x=numeric_var[,i],y=paste(numeric_var_name[i],", Mean(SD)",sep=""),z=c(60,10),data=data_keep[,by],self_define = self_define),round(t_test$p.value,2)),c(get_prop_multi(x_update,"Missing",data_keep[,by]),NA))
      }
      
      
    }
    
    colnames(cont_table)<-c(" ",levels(data_keep[,by]),"Overall",paste(compare_level[level_length],"VS",compare_level[1],", p-value",sep=" "))
    
    
    
  }
  
  
  cat_table<-cat_table[2:nrow(cat_table),]
  
  cont_table<-cont_table[2:nrow(cont_table),]
  
  char_table=rbind(c("Continuous (Mean(SD))",rep(NA,ncol(cont_table)-1)),cont_table,
                   c("Categorical(N(%))",rep(NA,ncol(cont_table)-1)),cat_table)
  length_N<-ifelse(compare==F,ncol(char_table),ncol(char_table)-1)
  
  for (i in 2:length_N){
    colnames(char_table)[i]<-paste(levels(data_keep[,by])[i-1]," (N=",prettyNum(sum(data_keep[,by]==levels(data_keep[,by])[i-1]),big.mark = ",", scientific = FALSE),")",sep="")
  }
  colnames(char_table)[length_N]<-paste("All (N=",prettyNum(nrow(data_keep),big.mark = ",", scientific = FALSE),")",sep="")
  rownames(char_table)<-seq(nrow(char_table))
  char_table
  
  
  
}

create_surv_formula<-function(model_type,data,full_adjustments,outcome,pgs_name){
  
  surv_formula<-paste("Surv(AGE_at_risk,",outcome,")",sep="")
  
  if (model_type=="simple"){
    output_formula=formula(paste(surv_formula,"~",pgs_name,sep=""))
  }else if (model_type=="partial"){
    output_formula=formula(paste(surv_formula,"~",paste(pgs_name,"SEX",
                                                        paste(PCs,collapse="+"),sep="+"),sep=""))
  }else if (model_type=="full"){
    output_formula=formula(paste(surv_formula,"~",paste(pgs_name,"SEX",
                                                        paste(full_adjustments,collapse = "+"),
                                                        paste(PCs,collapse="+"),sep="+"),sep=""))
  }else (print("ERROR: NO SUCH MODEL"))
  
  output_formula
  
}

run_cox<-function(model_type,data,full_adjustments,outcome,pgs_name){
  cat("\n Modelling ",pgs_name)
  model_output<-coxph(create_surv_formula(model_type,data,full_adjustments,outcome,pgs_name),data=data)
  model_output
}

create_output_table_cox<-function(trainsplit=F,data,train_data=NULL,test_data=NULL,model_type,cat_or_cont="cont", full_adjustments,outcome){
  
  if(trainsplit==F){
    train_data=data
    test_data=data
  }else if(trainsplit==T){
    train_data=train_data
    test_data=test_data
  }
  if (model_type=="simple"){
    train=train_data
    test=test_data
  }else if(model_type=="partial"){
    train<-na.omit(train_data%>%select(IID,AGE_at_risk,SEX,contains(c("PC","PGS","custom",outcome))))
    test<-na.omit(test_data%>%select(IID,AGE_at_risk,SEX,contains(c("PC","PGS","custom",outcome))))
  }else if (model_type=="full"){
    train<-na.omit(train_data%>%select(IID,AGE_at_risk,SEX,contains(c("PC","PGS","custom",outcome,full_adjustments))))
    test<-na.omit(train_data%>%select(IID,AGE_at_risk,SEX,contains(c("PC","PGS","custom",outcome,full_adjustments))))
  }
  
  
  if (cat_or_cont=="cont"){
    
    start_column=which( colnames(train)%in%colnames(variant_t))
    model_output_table<-data.frame(HR=NA,LR=NA,UR=NA,AUC=NA,AUC_LR=NA,AUC_UR=NA)
    
    for (i in 1:7){
      
      model_output<-run_cox(model_type,data=train,full_adjustments,outcome,pgs_name=colnames(train)[start_column[i]])
      
      
      hr_ci<-round(summary(model_output)$conf.int[1,3:4],3)
      
      c_index<-ci_manual(x=concordance(model_output)$concordance,se=sqrt(concordance(model_output)$var))
      model_output_table[i,]<-c(summary(model_output)$coefficient[1,2],hr_ci[1],hr_ci[2],
                                c_index[1],c_index[2],c_index[3])
      rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
      
    }
    
    model_output_table
    
  }else if (cat_or_cont=="cat"){
    
    start_column=which( colnames(train)%in%paste(colnames(variant_t),"_cat",sep=""))
    model_output_table<-data.frame(hr_int=NA,hr_intermediate_lr=NA,hr_intermediate_ur=NA,
                                   hr_high=NA,hr_high_lr=NA,hr_high_ur=NA,AUC=NA,AUC_LR=NA,AUC_UR=NA)
    for (i in 1:7){
      model_output<-run_cox(model_type,data=train,full_adjustments,outcome,pgs_name=colnames(train)[start_column[i]])
      
      hr_ci<-round(summary(model_output)$conf.int[1:2,3:4],3)
      c_index<-ci_manual(x=concordance(model_output)$concordance,se=sqrt(concordance(model_output)$var))
      model_output_table[i,]<-c(summary(model_output)$coefficient[1,2],hr_ci[1,1],
                                hr_ci[1,2],summary(model_output)$coefficient[2,2],hr_ci[2,1],hr_ci[2,2],
                                c_index[1],c_index[2],c_index[3])
      rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
    }
    
    model_output_table
    
  }
  
  
}



model_to_table_cox<-function(model_output_table){
  colength<-ncol(model_output_table)#7,10
  
  if(colength==6){
    model_output_table[,7]<-paste(round(model_output_table[,1],3),"(",round(model_output_table[,2],3),"-",round(model_output_table[,3],3),")",sep="")
    model_output_table[,8]<-paste(round(model_output_table[,4],3),"(",round(model_output_table[,5],3),"-",round(model_output_table[,6],3),")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(7,8),]
    rownames(model_output_table_t)[1:2]<-c("HR per SD","Concordance index")
    model_output_table_t<-cbind(rownames(model_output_table_t),model_output_table_t)
    colnames(model_output_table_t)[1]<-NA
  }else if (colength==9){
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







output_table_summary_cox<-function(outcome,data,full_adjustments,trainsplit=F,train_data=NULL,test_data=NULL){
  
  
  #continuous
  cat("\n Continuous PRS ")
  cat("\n Unadjusted model ")
  model_output_simple<-create_output_table_cox(trainsplit ,data=data,train_data,test_data,model_type = "simple",cat_or_cont = "cont",full_adjustments = full_adjustments,outcome=outcome)
  cat("\n Age and Sex Adjusted model (+7PCs) ")
  model_output_age_sex<-create_output_table_cox(trainsplit, data,train_data,test_data,model_type = "partial",cat_or_cont = "cont",full_adjustments = full_adjustments,outcome=outcome)
  cat("\n Fully Adjusted Model ")
  model_output_full<-create_output_table_cox(trainsplit , data,train_data,test_data,model_type = "full",cat_or_cont = "cont",full_adjustments = full_adjustments,outcome=outcome)
  
  model_output_simple_t<-model_to_table_cox(model_output_simple)
  model_output_age_sex_t<-model_to_table_cox(model_output_age_sex)
  model_output_full_t<-model_to_table_cox(model_output_full)
  
  
  #categorical
  cat("\n Categorical PRS ")
  cat("\n Unadjusted model ")
  model_output_simple_cat<-create_output_table_cox(trainsplit , data,train_data,test_data,model_type = "simple",cat_or_cont = "cat",full_adjustments = full_adjustments,outcome=outcome)
  cat("\n Age and Sex Adjusted model (+7PCs) ")
  model_output_age_sex_cat<-create_output_table_cox(trainsplit , data,train_data,test_data,model_type = "partial",cat_or_cont = "cat",full_adjustments = full_adjustments,outcome=outcome)
  cat("\n Fully Adjusted Model")
  model_output_full_cat<-create_output_table_cox(trainsplit , data,train_data,test_data,model_type = "full",cat_or_cont = "cat",full_adjustments = full_adjustments,outcome=outcome)
  
  model_output_simple_cat_t<-model_to_table_cox(model_output_simple_cat)
  model_output_age_sex_cat_t<-model_to_table_cox(model_output_age_sex_cat)
  model_output_full_cat_t<-model_to_table_cox(model_output_full_cat)
  
  
  prs_table_cont<-rbind(variant_t,c("Simple Model", rep(NA,7)),
                        model_output_simple_t,
                        
                        c("Partially adjusted Model", rep(NA,7)),
                        model_output_age_sex_t,
                        
                        c("Fully Adjusted Model", rep(NA,7)),
                        model_output_full_t
  )
  rownames(prs_table_cont)<-seq(1:nrow(prs_table_cont))
  colnames(prs_table_cont)[1]<-"95% CI"
  prs_table_cat<-rbind(variant_t,c("Simple Model", rep(NA,7)),
                       
                       model_output_simple_cat_t,
                       c("Adjusted Model", rep(NA,7)),
                       
                       model_output_age_sex_cat_t,
                       c("Fully Adjusted Model", rep(NA,7)),
                       
                       model_output_full_cat_t)
  rownames(prs_table_cat)<-seq(1:nrow(prs_table_cat))
  colnames(prs_table_cat)[1]<-"95% CI"
  prs_table<-list(prs_table_cont,prs_table_cat)
  prs_table
  
  prs_table
}








#create_output_table(trainsplit = F, model_type = "simple",cat_or_cont = "cont",full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA")

#output_table_summary(outcome="prevalent_CHD_EPA",data=data,full_adjustments,trainsplit=F,train_data=NULL,test_data=NULL)

#output_table_summary(outcome="prevalent_CHD_EPA",data=data,full_adjustments,trainsplit=F,train_data=NULL,test_data=NULL)

#discrimination_without_prs(train_data = train,test_data = test,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA")

#create_quntile_char(data_prs, pgs,variables_to_select)