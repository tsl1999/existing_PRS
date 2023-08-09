
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
library(forestploter)
library(grid)
library(officer)
# some theme-------------------------------
tm <- forest_theme(base_size = 5,
                   core = list(bg_params=list(fill = c("white"))),
                   summary_col = "black",
                   refline_lty = "solid",
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
                   ci_Theight = 0.1)

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
                    ci_Theight = 0.1)

fpp <- fp_par(text.align = "left", padding = 3)


readin_prs<-function(prs){
  
  for ( i in c(1:length(prs))) {
    cat("\n read in PRS data for ",prs[i], ".......")
    specific_prs_path<-paste(prs_path,prs[i],"score",sep="/")
    cat("\n path is ", specific_prs_path)
    setwd(specific_prs_path)
    dat=read.table("aggregated_scores.txt.gz", 
                   sep ="", header = TRUE, dec =".",fill=T)
    
    if(i==1){
      prs_reading<-data.frame(dat$IID,dat[4])
      colnames(prs_reading)<-c("IID",prs[i])
    }else{
      prs_reading<-merge(prs_reading,dat[,c(2,4)],by="IID")
      colnames(prs_reading)[i+1]<-prs[i]
    }
    
    
    
    name[i]<-prs[i]
    if(mean(dat$DENOM)==dat$DENOM[1]){
      number[i]<-mean(dat$DENOM)/2
    }
    
  }
  
  
  
  variant_use<-data.frame(name,number)
  return(list(prs_reading, variant_use))
}



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



create_glm_formula<-function(data,adjustments,outcome,pgs_name,int=NULL){
  if(is.null(adjustments)==T&is.null(int)==T){
    output_formula=formula(paste(outcome,"~",pgs_name))
  }else if (is.null(adjustments)==F&is.null(int)==T){
    output_formula=formula(paste(outcome,"~",paste(pgs_name,paste(adjustments,collapse = "+"),sep="+"),sep=""))
  } else {
    output_formula=formula(paste(outcome,"~",paste(pgs_name,paste(adjustments,collapse = "+"),paste(paste(int,pgs_name,sep="*"),collapse="+"),sep="+"),sep=""))
  }
  output_formula
  
}
#create_glm_formula(data,adjustments,outcome,"PGS000011")
#create_glm_formula(data,adjustments = full_adjustment ,outcome="prevalent_CHD_EPA",pgs_name = "PGS000011",int=c("SEX"))

run_glm<-function(data,adjustments,outcome,pgs_name,int=NULL){
  cat("\n Modelling ",pgs_name)
  model_output<-glm(create_glm_formula(data,adjustments,outcome,pgs_name,int),data=data,family=binomial(link='logit'))
  model_output
}

#run_glm(data,adjustments,outcome,"PGS000011")
#a<-run_glm(data,adjustments,outcome,"PGS000011",int="SEX")
create_output_table<-function(trainsplit=F,data,train_data=NULL,test_data=NULL,adjustments,outcome,namew,type="cont",roc=T,int=NULL){
  
  if(trainsplit==F){
    train_data=data
    test_data=data
  }else if(trainsplit==T){
    train_data=train_data
    test_data=test_data
  }
  train<-na.omit(train_data%>%select(IID,all_of(c(outcome,adjustments)),contains(c("PGS","custom"))))
  test<-na.omit(test_data%>%select(IID,all_of(c(outcome,adjustments)),contains(c("PGS","custom"))))
  
  
  
  start_column=which(colnames(train)%in%colnames(train%>%select(contains(colnames(model_table_top)[2:ncol(model_table_top)]))))
  
  if(type=="cont"){
    model_output_table<-data.frame(OR=NA,LR=NA,UR=NA,AUC=NA,AUC_LR=NA,AUC_UR=NA,Nagelkerke_PseudoR2=NA,lee_PseudoR2=NA,case_no=NA,control_no=NA,pval=NA,int_pval=NA)
    
    if(roc==T){
      png(paste(graphs_path,"/logistic/AUCplot_cont","_",outcome,"_",namew,".png",sep=""),res=150,width = 15, height = 15,units = "cm")}
    
    
    for (i in 1:length(start_column)){
      
      model_output<-run_glm(data=train,adjustments,outcome,colnames(train)[start_column[i]],int)
      pred<-predict(model_output,test,type='response')
      pred.obj<-prediction(pred,test[,outcome])
      perf.obj<-performance(pred.obj,"tpr","fpr")
      roc_prs<-roc(test[,outcome],pred)
      if(roc==T){
        if(i==1){
          plot(perf.obj,col=i)}else{plot(perf.obj,col=i,add=T)}}
      
      vr=runif(nrow(train),0,1)
      vsel=model_output$linear.predictors[model_output$y==0|vr<0.05]
      r2<-format(round(var(vsel)/(var(vsel)+pi^2/3),3),nsmall=3)
      
      
      auc_point<-auc(roc_prs)
      auc_ci<-format(round(ci.auc(test[,outcome],pred),3),nsmall=3)
      
      pval<-summary(model_output)$coefficients[,"Pr(>|z|)"][2]
      pval<-format(signif(pval,3),nsmall=3)
      #pval<-as.character(signif(pval,digits=2))
      #pval = sub("e"," 10^",pval) 
      if(is.null(int)==F){
        int_pval<-summary(model_output)$coefficients[,"Pr(>|z|)"][nrow(summary(model_output)$coefficients)]
        int_pval<-format(signif(int_pval,3),nsmall=3)}else{int_pval=NA}
      or_ci<-format(ci(model_output,roundup=3),nsmall=3)
      
      #pval<-ifelse(summary(model_output)$coefficients[,"Pr(>|z|)"][2]<0.001,"<0.001",round(summary(model_output)$coefficients[,"Pr(>|z|)"][2],3))
      model_output_table[i,]<-c(or_ci[1],or_ci[2],
                                or_ci[3],auc_ci[2],auc_ci[1],auc_ci[3],
                                format(round(PseudoR2(model_output,c("Nagelkerke")),3),nsmall=3),r2,
                                table(train[,outcome])["1"], table(train[,outcome])["0"],
                                pval,int_pval)
      
      
      
      rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
      
    }
    if(roc==T){
      abline(a=0,b=1,lty="dashed",col="gray")
      legend("bottomright",legend=c(paste(colnames(train)[start_column[1:length(start_column)]]," AUC=",round(as.numeric(model_output_table[,4]),3))),
             col=c(palette()[1:length(start_column)]),lty=1,cex=0.8)
      dev.off()}
    if(is.null(int)==F){model_output_table<-model_output_table}else{
      model_output_table<-model_output_table[,1:11]
    }
    
    output<-model_output_table
    
    
    output
  }else if(type=="cat"){
    model_output_table<-data.frame(AUC=NA,AUC_LR=NA,AUC_UR=NA,Nagelkerke_PseudoR2=NA,lee_PseudoR2=NA,case_no=NA,control_no=NA)   
    model_estimate_list<-list()
    if (roc==T){
      png(paste(graphs_path,"/logistic/AUCplot_cat","_",outcome,"_",namew,".png",sep=""),res=150,width = 15, height = 15,units = "cm")}
    for (i in 1:length(start_column)){
      model_output<-run_glm(data=train,adjustments,outcome,colnames(train)[start_column[i]],int=NULL)
      fl_a<-float(model_output,factor =colnames(train)[start_column[i]])
      model_estimate_output_table<-data.frame(estimate = round(exp(fl_a$coef), 3),
                                              LR = round(exp(fl_a$coef - 1.96 * sqrt(fl_a$var)), 3),
                                              UR = round(exp(fl_a$coef + 1.96 * sqrt(fl_a$var)), 3),
                                              SE=fl_a$var,
                                              case=table(train[,colnames(train)[start_column[i]]],train[,outcome])[,2],
                                              control=table(train[,colnames(train)[start_column[i]]],train[,outcome])[,1])
      model_estimate_output_table$category<-rownames(model_estimate_output_table)
      model_estimate_output_table<-model_estimate_output_table[,c(7,1:6)]
      colnames(model_estimate_output_table)[1]<-colnames(train)[start_column[i]]
      model_estimate_list[[i]]<-model_estimate_output_table
      pred<-predict(model_output,test,type='response')
      pred.obj<-prediction(pred,test[,outcome])
      perf.obj<-performance(pred.obj,"tpr","fpr")
      roc_prs<-roc(test[,outcome],pred)
      
      if(roc==T){
        if(i==1){
          plot(perf.obj,col=i)}else{plot(perf.obj,col=i,add=T)}}
      
      vr=runif(nrow(train),0,1)
      vsel=model_output$linear.predictors[model_output$y==0|vr<0.05]
      r2<-round(var(vsel)/(var(vsel)+pi^2/3),3)
      
      auc_ci<-round(ci.auc(test[,outcome],pred),3)
      
      model_output_table[i,]<-c(auc_ci[2],auc_ci[1],auc_ci[3],
                                round(PseudoR2(model_output,c("Nagelkerke")),3),r2,
                                table(train[,outcome])["1"], table(train[,outcome])["0"])
      rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
    }
    if(roc==T){
      abline(a=0,b=1,lty="dashed",col="gray")
      legend("bottomright",legend=c(paste(colnames(train)[start_column[1:length(start_column)]]," AUC=",round(as.numeric(model_output_table[,1]),3))),
             col=c(palette()[1:length(start_column)]),lty=1,cex=0.8)
      dev.off()  }
    output<-list(model_output_table,model_estimate_list)
  }
  
  
  
  
  output
  
}




plot_FAR<-function(model_output_list,name,outcome,type="OR"){
  or_hr<-ifelse(type=="OR","Odds Ratio","Hazard Ratio")
  saveor_hr<-ifelse(type=="OR","/logistic","/cox")
no_model<-length(model_output_list)
for (i in 1:7){
    model_name<-colnames(model_output_list[[1]][[2]][[i]])[1]
    png(paste(graphs_path,saveor_hr,"/FAR_plot_",outcome,"_",sub("_cat","",model_name),".png",sep=""),width = 15*no_model,height = 15,units = "cm",res=50*no_model)
    par(mfrow=c(1,no_model))
    ylim=data.frame(min=NA,max=NA)
      for(j in 1:no_model){
      ylim[j,1]<-min(model_output_list[[j]][[2]][[i]]$LR)
      ylim[j,2]<-max(model_output_list[[j]][[2]][[i]]$UR)
    }
    ylim_all=c(min(ylim$min),max(ylim$max))
    for(j in 1:no_model){
    model_estimate<-model_output_list[[j]][[2]][[i]]
   if(j==1){p<-plot(y=model_estimate$estimate,x=model_estimate[,1],pch=15,ylim=c(ylim_all[1]-0.1,ylim_all[2]+0.1),
                 xaxt="n",xlab="Polygenic risk score deciles" ,ylab=paste(or_hr, "floating absoluted risk (95%CI) of ",
                                                                          sub("_cat","",colnames(model_estimate)[1])),main=name[j],bty = "l")
   }else{plot(y=model_estimate$estimate,x=model_estimate[,1],pch=15,ylim=c(ylim_all[1]-0.1,ylim_all[2]+0.1),
              xaxt="n",xlab="Polygenic risk score deciles" ,ylab=" ",main=name[j],yaxt="n",bty="l")}
    
    axis(1,at=1:11,labels = T)
    lines(y=model_estimate$estimate,x=model_estimate[,1],pch=15,type="l",col="grey",lty="dashed")
    
    arrows(x0=as.numeric(model_estimate[,1]), y0=model_estimate$LR, x1=as.numeric(model_estimate[,1]), y1=model_estimate$UR, code=3, angle=90, length=0, col="black", lwd=1)
    text(x=as.numeric(model_estimate[,1]),y=(model_estimate$LR-0.05),labels=model_estimate$LR,cex=0.7)
    text(x=as.numeric(model_estimate[,1]),y=(model_estimate$UR+0.05),labels=model_estimate$UR,cex=0.7)
    }
    dev.off()
  
  }
par(mfrow=c(1,1))
}

plot_FAR_up<-function(model_output_list,name,outcome,type="OR"){
  or_hr<-ifelse(type=="OR","Odds Ratio","Hazard Ratio")
  saveor_hr<-ifelse(type=="OR","/logistic","/cox")
  no_model<-length(model_output_list)
  
  
  
  
  for(j in 1:no_model){
    png(paste(graphs_path,saveor_hr,"/FAR_plot_",outcome,"_",name[j],".png",sep=""),width = 30*no_model,height = 15,units = "cm",res=80*7)
    ylim=data.frame(min=NA,max=NA)
    par(mfrow=c(1,7),mar = c(4,4,2,2))
    for(i in 1:7){
      ylim[i,1]<-min(model_output_list[[j]][[2]][[i]]$LR)
      ylim[i,2]<-max(model_output_list[[j]][[2]][[i]]$UR)
    }
    
    ymax=max(ylim$max)
    ymin=min(ylim$min)
    for (i in 1:7){
      
      model_name<-paste(colnames(model_table_top)[i+1],"_cat",sep="")
      model_estimate<-model_output_list[[j]][[2]][[which(sapply(model_output_list[[j]][[2]], function(x) colnames(x)[1]==model_name))]]
      model_pgs_reading<-data[,c(sub("_cat","",colnames(model_estimate)[1]),model_name)]
      model_pgs_reading_decile<-data.frame(decile=NA,avg=NA)
      model_pgs_reading_decile[1:length(levels(model_pgs_reading[,2])),1]<-levels(model_pgs_reading[,2])
      for (k in 1:length(levels(model_pgs_reading[,2]))) {
        model_pgs_reading_decile[k,2]<-mean(model_pgs_reading[model_pgs_reading[,2]==model_pgs_reading_decile[k,1],1])
      }
      
      
      if(i==1){
        p<-plot(y=log(model_estimate$estimate),x=model_pgs_reading_decile$avg,pch=15,ylim=c(log(ymin-0.1),
                log(ymax+0.1)),xlim=c(-1.8,2),xaxt="n",yaxt="n",xlab="Polygenic risk score deciles" ,
                ylab=paste(or_hr, "on log-scale (95%CI)"),cex=0.0015/model_estimate$SE,
                main=paste(model_table_top[3,sub("_cat","",colnames(model_estimate)[1])],"\n Publication date:",model_table_top[2,sub("_cat","",colnames(model_estimate)[1])]),bty = "l",cex.main = 0.8)
        axis(2,at=log(c(seq(0.9,ymax,0.2))),labels = c(seq(0.9,ymax,0.2)))
      }else{
        p<- plot(y=log(model_estimate$estimate),x=model_pgs_reading_decile$avg,pch=15,
                 ylim=c(log(ymin-0.1),log(ymax+0.1)),xlim=c(-1.8,2),
                 xaxt="n",yaxt="n",xlab="Polygenic risk score deciles" ,ylab=" ",cex=0.0015/model_estimate$SE,
                 main=paste(model_table_top[3,sub("_cat","",colnames(model_estimate)[1])],
                            "\n Publication date:",model_table_top[2,sub("_cat","",colnames(model_estimate)[1])]),bty = "l",cex.main = 0.8)}
      axis(1,at=model_pgs_reading_decile$avg,labels = model_pgs_reading_decile$decile)
      
      lines(y=log(model_estimate$estimate),x=model_pgs_reading_decile$avg,pch=15,type="l",col="grey",lty="dashed")
      
      arrows(x0=model_pgs_reading_decile$avg, y0=log(model_estimate$LR), x1=model_pgs_reading_decile$avg, y1=log(model_estimate$UR), code=3, angle=90, length=0, col="black", lwd=1)
      text(x=model_pgs_reading_decile$avg,y=log(model_estimate$LR-0.05),labels=paste(model_estimate$case,"/",model_estimate$control),cex=0.6)
      text(x=model_pgs_reading_decile$avg,y=log(model_estimate$UR+0.05),labels=model_estimate$estimate,cex=0.7)
      
      
      
      
      
    }
    
    dev.off()
    
  }
  par(mfrow=c(1,1))
}

plot_FAR_gg<-function(model_output_list,name,outcome,type="OR",data=data){
  or_hr<-ifelse(type=="OR","Odds Ratio","Hazard Ratio")
  saveor_hr<-ifelse(type=="OR","/logistic","/cox")
  no_model<-length(model_output_list)
  
  
  library(grid)
  library(gridExtra)
  
  for(j in 1:no_model){
    png(paste(graphs_path,saveor_hr,"/FAR_ggplot_",outcome,"_",name[j],".png",sep=""),width = 20*no_model,height = 30,units = "cm",res=80*7)
    ylim=data.frame(min=NA,max=NA)
    
    for(i in 1:(ncol(model_table_top)-1)){
      ylim[i,1]<-min(model_output_list[[j]][[2]][[i]]$LR)
      ylim[i,2]<-max(model_output_list[[j]][[2]][[i]]$UR)
    }
    
    ymax=max(ylim$max)
    ymin=min(ylim$min)
    p<-list()
    for (i in 1:(ncol(model_table_top)-1)){
      
      model_name<-paste(colnames(model_table_top)[i+1],"_cat",sep="")
      model_estimate<-model_output_list[[j]][[2]][[which(sapply(model_output_list[[j]][[2]], function(x) colnames(x)[1]==model_name))]]
      model_pgs_reading<-data[,c(sub("_cat","",colnames(model_estimate)[1]),model_name)]
      model_pgs_reading_decile<-data.frame(decile=NA,avg=NA)
      model_pgs_reading_decile[1:length(levels(model_pgs_reading[,2])),1]<-levels(model_pgs_reading[,2])
      for (k in 1:length(levels(model_pgs_reading[,2]))) {
        model_pgs_reading_decile[k,2]<-mean(model_pgs_reading[model_pgs_reading[,2]==model_pgs_reading_decile[k,1],1])
      }
      
      
      if(i==1|(i-1)%%4==0){
        p[[i]]<-ggplot(model_estimate, aes(x=model_pgs_reading_decile$avg, y=estimate)) + 
          geom_line(linetype="dashed",color="gray") +
          geom_point(aes(size=0.0001/SE),shape=15)+
          annotate("text",x=model_pgs_reading_decile$avg, y=model_estimate$UR+0.05, label=model_estimate$estimate,size=3)+
          annotate("text",x=model_pgs_reading_decile$avg, y=model_estimate$LR-0.05, label=paste(model_estimate$case),size=3)+
          geom_pointrange(aes(ymin=LR, ymax=UR), 
                          position=position_dodge(0.05),shape=15)+theme_classic()+
          theme(legend.position = "none",plot.title = element_text(size = 7, face = "bold"),axis.title.y = element_text(face="bold"))+
          scale_fill_discrete(name = " ")+scale_y_continuous(trans = "log",limits = c(ymin-0.1,ymax+0.1), breaks = seq(0.9, ymax,0.2))+
          ggtitle(paste(model_table_top[3,sub("_cat","",colnames(model_estimate)[1])],"\n No. of SNPs:",model_table_top[1,sub("_cat","",colnames(model_estimate)[1])]))+
          xlab(" ")+ylab(paste(or_hr, "(95%CI)"))+
          scale_x_continuous(breaks=c(model_pgs_reading_decile$avg),labels=c(1:length(model_pgs_reading_decile$avg)),limits = c(-1.8,1.8))
        
        
        
        
      }else{
        p[[i]]<-ggplot(model_estimate, aes(x=model_pgs_reading_decile$avg, y=estimate)) + 
          geom_line(linetype="dashed",color="gray") +
          geom_point(aes(size=0.0001/SE),shape=15)+
          annotate("text",x=model_pgs_reading_decile$avg, y=model_estimate$UR+0.05, label=model_estimate$estimate,size=3)+
          annotate("text",x=model_pgs_reading_decile$avg, y=model_estimate$LR-0.05, label=paste(model_estimate$case),size=3)+
          geom_pointrange(aes(ymin=LR, ymax=UR),
                          position=position_dodge(0.05),shape=15)+
          theme_classic()+theme(legend.position = "none",plot.title = element_text(size = 7, face = "bold"),axis.title.y = element_text(face="bold"))+
          scale_fill_discrete(name = " ")+scale_y_continuous(trans = "log",limits = c(ymin-0.1,ymax+0.1), breaks = NULL)+
          ggtitle(paste(model_table_top[3,sub("_cat","",colnames(model_estimate)[1])],"\n No. of SNPs:",model_table_top[1,sub("_cat","",colnames(model_estimate)[1])]))+
          xlab(" ")+ylab(" ")+
          scale_x_continuous(breaks=c(model_pgs_reading_decile$avg),labels=c(1:length(model_pgs_reading_decile$avg)),limits = c(-1.8,1.8))
        
        
        
      }
      
      
    }
    
    
    grid.arrange(arrangeGrob(grobs=lapply(p, function(p) p + guides(scale="None")), ncol=4, 
                             bottom=textGrob("Polygenic risk score quintiles",
                                             gp=gpar(fontface="bold", col="Black", fontsize=10)),
                             top=textGrob(paste("Association of Polygenic risk scores and",outcome,"in",name[j],"model",sep=" "),gp=gpar(fontface="bold", col="Black", fontsize=14)),
                             sub = textGrob("Footnote", x = 2, hjust = 1, vjust=1, gp = gpar(fontface = "italic", fontsize = 10))))
    
    
    
    
    
    dev.off()
    
    
  }}
  
#create_output_table(trainsplit = F,data,adjustments=adjustments,outcome=outcome,namew=NULL)

create_surv_formula<-function(adjustments,strata=NULL,outcome,pgs_name){
  surv_formula<-paste("Surv(time_in,time_out,",outcome,")",sep="")
  if(is.null(adjustments)==T&is.null(strata)==T){
    output_formula=formula(paste(surv_formula,
                                 "~",pgs_name,sep=""))}else if (is.null(adjustments)==F&is.null(strata)==T){output_formula=formula(paste(surv_formula,
                               "~",paste(pgs_name,paste(adjustments,collapse = "+"),sep = "+"),sep=""))
}else{
  output_formula=formula(paste(surv_formula, "~",
                               paste(pgs_name,"+",paste(adjustments,collapse = "+"),"+ strata(",strata,")",sep = ""),sep=""))
}  
  output_formula
  }

#create_surv_formula(adjustments=adjustments,outcome=outcome,pgs_name=pgs_name)

run_cox<-function(data,adjustments,strata=NULL,outcome,pgs_name){
  cat("\n Modelling ",pgs_name)
  model_output<-coxph(create_surv_formula(adjustments,strata,outcome,pgs_name),data=data)
  model_output
}

#run_cox(data=data,adjustments=adjustments,outcome=outcome,pgs_name=pgs_name)


create_output_table_cox<-function(trainsplit=F,data,train_data=NULL,test_data=NULL,adjustments,strata=NULL,outcome,namew=NULL,type="cont",ph=T){
  
  if(trainsplit==F){
    train_data=data
    test_data=data
  }else if(trainsplit==T){
    train_data=train_data
    test_data=test_data
  }
  train<-na.omit(train_data%>%select(IID,time_in,time_out,all_of(c(outcome,adjustments,strata)),contains(c("PGS","custom"))))
  test<-na.omit(test_data%>%select(IID,time_in,time_out,all_of(c(outcome,adjustments,strata)),contains(c("PGS","custom"))))
  
  start_column=which(colnames(train)%in%colnames(train%>%select(contains(colnames(model_table_top)[2:ncol(model_table_top)]))))
  if(type=="cont"){
    model_output_table<-data.frame(HR=NA,LR=NA,UR=NA,AUC=NA,AUC_LR=NA,AUC_UR=NA,case_no=NA,control_no=NA,pval=NA)
    
    for (i in 1:length(start_column)){
      
      model_output<-run_cox(data=train,adjustments,strata,outcome,pgs_name=colnames(train)[start_column[i]])
      dir.create(paste(graphs_path,"/cox/",colnames(train)[start_column[i]],sep=""))
      if (ph==T){
      
      jpeg(paste(graphs_path,"/cox/",colnames(train)[start_column[i]],"/",colnames(train)[start_column[i]],"_Schoenfeld_residual",type,namew,".png",sep=""),width = 800, height = 600)
      cox.zph.fit <- cox.zph(model_output)
      p<-ggcoxzph(cox.zph.fit,font.main = 8)
      print(p)
      dev.off()}
      
      hr_ci<-format(round(summary(model_output)$conf.int[1,c(1,3:4)],3),nsmall=3)
      
      c_index<-format(round(ci_manual(x=concordance(model_output)$concordance,se=sqrt(concordance(model_output)$var)),3),nsmall=3)
      case_no<-sum(train[,outcome]==1)
      case_id<-unique(train[train[,outcome]==1,"IID"])
      
      control<-train[!train$IID%in%case_id,]
      control_no<-length(unique(control[control[,outcome]==0,"IID"]))
      pval<-summary(model_output)$coefficients[colnames(train)[start_column[i]],"Pr(>|z|)"]
      pval<-format(signif(pval,3),nsmall=3)
      
      model_output_table[i,]<-c(hr_ci[1],hr_ci[2],hr_ci[3],
                                c_index[1],c_index[2],c_index[3],case_no,control_no,pval)
      rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
      
    }
    output<-model_output_table
    }else if(type=="cat"){
      model_output_table<-data.frame(AUC=NA,AUC_LR=NA,AUC_UR=NA,case_no=NA,control_no=NA)
      model_estimate_list<-list()
      for (i in 1:length(start_column)){
        model_output<-run_cox(data = train,adjustments,strata,outcome,colnames(train)[start_column[i]])
        fl_a<-float(model_output,factor =colnames(train)[start_column[i]])

        model_estimate_output_table<-data.frame(estimate = round(exp(fl_a$coef), 3),
                                                LR = round(exp(fl_a$coef - 1.96 * sqrt(fl_a$var)), 3),
                                                UR = round(exp(fl_a$coef + 1.96 * sqrt(fl_a$var)), 3),
                                                SE=fl_a$var,
                                                case=NA,
                                                control=NA)
        for (j in 1: length(levels(train[,colnames(train)[start_column[i]]]))){
          cat_only_row<-train[train[,colnames(train)[start_column[i]]]==j,]
          case_id<-unique(cat_only_row[cat_only_row[,outcome]==1,"IID"])
          case_no<-sum(cat_only_row[,outcome]==1)
          control<-cat_only_row[!cat_only_row$IID%in%case_id,]
          control_no<-length(unique(control[control[,outcome]==0,"IID"]))
          model_estimate_output_table$case[j]<-case_no
          model_estimate_output_table$control[j]<-control_no
        }
        model_estimate_output_table$category<-rownames(model_estimate_output_table)
        model_estimate_output_table<-model_estimate_output_table[,c(7,1:6)]
        colnames(model_estimate_output_table)[1]<-colnames(train)[start_column[i]]
        model_estimate_list[[i]]<-model_estimate_output_table
        
        
        if (ph==T){
        jpeg(paste(graphs_path,"/cox/",sub("_cat","",colnames(train)[start_column[i]]),"/",colnames(train)[start_column[i]],"_Schoenfeld_residual",type,namew,".png",sep=""),width = 800, height = 600)
        cox.zph.fit <- cox.zph(model_output)
        p<-ggcoxzph(cox.zph.fit,font.main = 8)
        print(p)
        dev.off()}
        
        c_index<-format(round(ci_manual(x=concordance(model_output)$concordance,se=sqrt(concordance(model_output)$var)),3),nsmall=3)
        case_no<-sum(train[,outcome]==1)
        case_id<-unique(train[train[,outcome]==1,"IID"])
        
        control<-train[!train$IID%in%case_id,]
        control_no<-length(unique(control[control[,outcome]==0,"IID"]))
        
        model_output_table[i,]<-c(c_index[1],c_index[2],c_index[3],case_no,control_no)
        rownames(model_output_table)[i]<-colnames(train)[start_column[i]]
        
        
        
        
      }
      output<-list(model_output_table,model_estimate_list) 
    }  
  output
  
}



model_to_table_cox<-function(model_output_table){
  colength<-ncol(model_output_table)#7,10
  model_output_table[,1:(ncol(model_output_table)-1)]<-apply(model_output_table[,1:(ncol(model_output_table)-1)],2,as.numeric)
  
  if(colength==9){
    model_output_table[,7]<-paste(model_output_table[,1],"(",model_output_table[,2],"-",model_output_table[,3],")",sep="")
    model_output_table[,8]<-paste(model_output_table[,4],"(",model_output_table[,5],"-",model_output_table[,6],")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(7,8),]
    rownames(model_output_table_t)[1:2]<-c("HR per SD","Concordance index")
    model_output_table_t<-cbind(rownames(model_output_table_t),model_output_table_t)
    colnames(model_output_table_t)[1]<-NA
  }else if (colength==5){
    model_output_table[,6]<-paste(model_output_table[,1],"(",model_output_table[,2],"-",model_output_table[,3],")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(6),]
    colnames(model_output_table_t)<-sub("_cat","",colnames(model_output_table_t))
    rownames(model_output_table_t)[1]<-c("Concordance index")
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



model_to_table_log<-function(model_output_table){
  colength<-ncol(model_output_table)#7,11
  
  rownames(model_output_table)<-rownames(model_output_table)
  model_output_table<-data.frame(model_output_table)
  if(colength==11){
    model_output_table[,12]<-paste(model_output_table[,1],"(",model_output_table[,2],"-",model_output_table[,3],")",sep="")
    model_output_table[,13]<-paste(model_output_table[,4],"(",model_output_table[,5],"-",model_output_table[,6],")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(12,13,7,8),]
    rownames(model_output_table_t)[1:2]<-c("OR per SD","AUC")
    model_output_table_t<-cbind(rownames(model_output_table_t),model_output_table_t)
    colnames(model_output_table_t)[1]<-NA
  }else if (colength==7){
    model_output_table[,8]<-paste(model_output_table[,1],"(",model_output_table[,2],"-",model_output_table[,3],")",sep="")
    model_output_table_t<-data.frame(t(model_output_table))
    model_output_table_t<-model_output_table_t[c(8,4,5),]
    colnames(model_output_table_t)<-sub("_cat","",colnames(model_output_table_t))
    rownames(model_output_table_t)[1]<-c("AUC")
    model_output_table_t<-cbind(rownames(model_output_table_t),model_output_table_t)
    colnames(model_output_table_t)[1]<-NA
  }
  model_output_table_t
}





model_to_table<-function(model_output_table,log_cox){
  if(log_cox=="log"){
    model_output_table_t<-model_to_table_log(model_output_table)
  }else if( log_cox=="cox"){
    model_output_table_t<-model_to_table_cox(model_output_table)
  }
  model_output_table_t
}




generate_forest_table<-function(tables,model_name,or_hr="OR",show_pgsname=F,int=NULL){
  hr_ci<-NA
  if (is.null(int)==T){
  forest_table<-data.frame(subgroup=NA,case_no=NA,control_no=NA,OR=NA,LR=NA,UR=NA,hr_ci=NA,auc_ci=NA,pval=NA)}else{
    forest_table<-data.frame(subgroup=NA,case_no=NA,control_no=NA,OR=NA,LR=NA,UR=NA,hr_ci=NA,auc_ci=NA,pval=NA,int_pval=NA)
  }
  n_table<-length(tables)
  name<-c(rep(NA,n_table))
  for (i in 1:n_table){
    name[i]<-paste("\U{00A0}\U{00A0}",model_name[i],sep="")
  }
  
  model_list<-tables
  for (i in 1:(ncol(model_table_top)-1)){
    if (show_pgsname==T){
  forest_table<-rbind(forest_table,c(model_table_top[3,i+1],rep(" ",ncol(forest_table)-1)))
    }else{
    forest_table<-rbind(forest_table,c(colnames(model_table_top)[i+1],rep(" ",ncol(forest_table)-1)))}
    for (j in 1:n_table){
      row<-rownames(model_list[[j]])[colnames(model_table_top)[i+1]==rownames(model_list[[j]])]
      row_computing<- model_list[[j]][row,]
      hr_ci<-paste(format(row_computing[1],nsmall=3),"(",format(row_computing[2],nsmall=3),"-",format(row_computing[3],nsmall=3),")",sep="")
      auc_ci<-paste(format(row_computing[4],nsmall=3),"(",format(row_computing[5],nsmall=3),"-",format(row_computing[6],nsmall=3),")",sep="")
      
      if(is.null(int)==T){
      row_input<-unlist(c(name[j],model_list[[j]][row,c("case_no","control_no",or_hr,"LR","UR")],
                          hr_ci,auc_ci,model_list[[j]][i,c("pval")]))}else if (is.null(int)==F){
      row_input<-unlist(c(name[j],model_list[[j]][row,c("case_no","control_no",or_hr,"LR","UR")],
                            hr_ci,auc_ci,model_list[[j]][i,c("pval","int_pval")]))  
                          }
      forest_table<-rbind(forest_table,row_input,make.row.names=F)
    }
    
  }
  forest_table<-forest_table[2:nrow(forest_table),]
  
  forest_table$" "<- paste(rep(" ", 40), collapse = " ")
  if(is.null(int)==T){
    colnames(forest_table)<-c("PRS name","Number of CAD cases", "Number of CAD controls", "estimate","LR","UR"," ","AUC (95%CI)","p-value","   ")
  }else{
    colnames(forest_table)<-c("PRS name","Number of CAD cases", "Number of CAD controls", "estimate","LR","UR"," ","AUC (95%CI)","p-value",paste("Interaction with",int,"(p-value)" ,sep=" "),"   ")
  }
  
  
  
  
  forest_table
}


plot_forest<-function(forest_table,xlim,ticks_at,theme,hr_or="OR",theme_in="tm2",col_input=c(1:3,7,10,8,9)){
  if(theme_in=="tm"){
  p <- forest(forest_table[,col_input],
              est = as.numeric(forest_table$estimate),
              lower = as.numeric(forest_table$LR), 
              upper = as.numeric(forest_table$UR),
              ci_column = 5,
              ref_line = 1,
              vert_line = c(as.numeric(forest_table[forest_table[,1]=="Overall",]$estimate)),
              arrow_lab = c("Lower risk of CAD", "Higher risk of CAD"),
              xlim = xlim,
              ticks_at = ticks_at,
              is_summary = c(rep(FALSE, nrow(forest_table)-1), TRUE),
              footnote = "\n\n\n\n\n\n\n\n\n Odds Ratios estimated separately by levels of included variable, adjusted for sex, baseline age,education level, \n smoking status, blood pressure and waist-hip ratio if they are not included in the strata. \n For analysis on HbA1c strata, adjustments excluded baseline diabetes.",
              theme=tm)}else{
                p <- forest(forest_table[,col_input],
                            est = as.numeric(forest_table$estimate),
                            lower = as.numeric(forest_table$LR), 
                            upper = as.numeric(forest_table$UR),
                            ci_column = 5,
                            ref_line = 1,
                            arrow_lab = c("Lower risk of CAD", "Higher risk of CAD"),
                            xlim = xlim,
                            ticks_at = ticks_at,
                            footnote = ifelse(hr_or=="OR","\n\n\n\n\n\n\n\n\n Partially adjusted model adjusts for sex and baseline age. Fully adjusted model additionally adjusts for education level, smoking status, \n blood pressure and waist-hip ratio ",
                                              "\n\n\n\n\n\n\n\n\n Partially adjusted model adjusts for sex and age group. Fully adjusted model additionally adjusts for education level, smoking status, \n blood pressure and waist-hip ratio "),
                            theme=tm2)
              }
  
  p <- add_border(p, 
                  part = "header", 
                  row = 1,
                  col = 1:length(col_input), gp = gpar(lwd = 1))
  
  p <- add_border(p, 
                  part = "body", 
                  col = c(2:4,6:length(col_input)),row=1:nrow(forest_table),
                  gp = gpar(lwd = .5))
  
  
  if(theme_in=="tm2"){
  p<- edit_plot(p, row = c(seq(1,nrow(forest_table),by=nrow(forest_table)/npgs)), col=1,
                gp = gpar(fontface = "bold"))}
  p <- add_text(p, text = paste(hr_or, "per SD (95% CI)",sep=" "),
                part = "header", 
                col = 4:5,row=1,
                gp = gpar(fontface="bold",fontsize=5))
  p
  
}

combine_model_tables<-function(tables,name,log_cox){
  n_table<-length(tables)
  model_t<-list(rep(NA,n_table))
  for (i in 1:n_table){
    model_output_t<-model_to_table(tables[[i]],log_cox)
    model_output_t[,1]<-paste("\U{00A0}\U{00A0}",model_output_t[,1],sep="")
    model_t[[i]]<-model_output_t
  }
  
  model_output_combine<-rbind(model_table_top[1:2,])

  
  for (i in 1:n_table){
    
    model_output_combine<-rbind(model_output_combine,c(name[i], rep(NA,ncol(model_t[[i]])-1)),
                                model_t[[i]])
  }
  
  
  rownames(model_output_combine)<-seq(nrow(model_output_combine))
  
  model_output_combine<-data.frame(model_output_combine)
  colnames(model_output_combine)<-model_table_top[3,]
  model_output_combine
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




discrimination_without_prs<-function(train_data,test_data,partial_adjustments,full_adjustments,outcome){
  #partial#
  train_partial<-na.omit(train_data%>%select(IID,partial_adjustments,contains(c("PGS","custom",outcome))))
  test_partial<-na.omit(test_data%>%select(IID,partial_adjustments,contains(c("PC","PGS","custom",outcome))))
  model_partial<-run_glm(data=train_partial,partial_adjustments,outcome,pgs_name=NULL)
  #full#
  train_full<-na.omit(train_data%>%select(IID,full_adjustments,contains(c("PGS","custom",outcome))))
  test_full<-na.omit(train_data%>%select(IID,full_adjustments,contains(c("PGS","custom",outcome))))
  model_full<-run_glm(data=train_full,full_adjustments,outcome,pgs_name=NULL)
  
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






prev_cont_flextable<-function(table,table_caption=NULL){
  
  ncol_table<-ncol(table)
  flex_cont<-my_table(table,y=c(1:3,8),x=1:ncol_table)
  flex_cont<-footnote(flex_cont,i=c(3,8),j=1,part="body",value = as_paragraph(c("Partially adjusted models adjust for Sex and baseline age. Fully adjusted models additionally adjust for waist-to-hip ratio, systolic and diastolic blood pressure, education attainmen level, smoking status,diabetes at baseline")),ref_symbols = c("*"))
  
  flex_cont<-add_footer_lines(flex_cont,value = as_paragraph(c("The table is ordered by PRSs target ethnicity and publication date")))
  flex_cont <- add_header_row(flex_cont, values =  c(" ", "European PRS","Non-European PRS"),
                              colwidths = c(1, 4,ncol_table-5), top = T)
  
  flex_cont <- compose(flex_cont, i = c(6,11), j = 1, part = "body",
                       value = as_paragraph(
                         "\U{00A0}\U{00A0}Naglekerke pseudo R",
                         as_sup("2")   ) )
  
  flex_cont <- compose(flex_cont, i = c(7,12), j = 1, part = "body",
                       value = as_paragraph(
                         "\U{00A0}\U{00A0}Lee's pseudo R",
                         as_sup("2")   ) )

  flex_cont<-set_caption(
    flex_cont,
    caption = table_caption,align_with_table = FALSE,
    word_stylename = "Table Caption",
    fp_p = fpp
  )
  
  my_theme(flex_cont,y=c(1,2,7),set_padding = 3,fontsize = 8,header_border = 7)
  
  
  
}

run_stratified<-function(data,strata_in,outcome,adjustment){
  all_strata_outcome<-list()
  for(i in 1: length(strata_in)){
    strata_var<-strata_in[i]
    print(strata_var)
    var_in<-ifelse(class(data[,strata_var])%in%c("integer","numeric"),
                   paste(strata_var,"_cut",sep=""),strata_var)
    var<-data[,var_in]
    strata_outcome_table<-list()
    for(j in 1:length(levels(var))){
      print(levels(var)[j])
      data_stratified<-data[data[,var_in]==levels(var)[j],]
      print(table(data_stratified[,outcome]))
      if(strata_var=="BASE_HBA1C"){
        adjustments<-adjustment[!adjustment%in%c(strata_var,"diabetes_at_baseline")]
      }else{adjustments<-adjustment[!adjustment%in%c(strata_var)]}
      model_output_stratify<-create_output_table(
        trainsplit = F,data_stratified,adjustments=adjustments,
        outcome=outcome,namew=paste(strata_var,levels(var)[j],sep="_"),roc = F)
      
      strata_outcome_table[[j]]<-model_output_stratify
      
      
    }
    all_strata_outcome[[i]]<-strata_outcome_table
  }
  
  
  all_strata_outcome[[length(strata_in)+1]]<-create_output_table(
    trainsplit = F,data,adjustments=adjustment,
    outcome=outcome,namew="overall",roc = F)
  all_strata_outcome
}



stratified_analysis_table<-function(all_strata_outcome, all_name,strata_in,data){
  pgs=which(colnames(data)%in%colnames(data%>%select(contains(colnames(model_table_top)[2:ncol(model_table_top)]))))
  forest_table_list<-list()
  for (i in 1:length(pgs)){
    pgsname<-colnames(data)[pgs[i]]
    forest_table<-data.frame(Variable=NA,case_no=NA,control_no=NA,OR=NA,LR=NA,UR=NA,hr_ci=NA,auc_ci=NA,pval=NA)
    
    for(j in 1:(length(all_strata_outcome)-1)){
      forest_table<-rbind(forest_table,c(all_name[j],rep(" ",8)))
      strata_var<-strata_in[j]
      var_in<-ifelse(class(data[,strata_var])%in%c("integer","numeric"),
                     paste(strata_var,"_cut",sep=""),strata_var)
      var_level<-levels(data[,var_in])
      strata_outcome<-all_strata_outcome[[j]]
      for(k in 1:length(var_level)){
        stratified_var<-var_level[k]
        strata_outcome_row<-strata_outcome[[k]][pgsname,]
        hr_ci<-paste(format(strata_outcome_row[1],nsmall=3),"(",
                     format(strata_outcome_row[2],nsmall=3),"-",
                     format(strata_outcome_row[3],nsmall=3),")",sep="")
        auc_ci<-paste(format(strata_outcome_row[4],nsmall=3),"(",
                      format(strata_outcome_row[5],nsmall=3),"-",
                      format(strata_outcome_row[6],nsmall=3),")",sep="")
        row_input<-unlist(c(paste("\U{00A0}\U{00A0}",stratified_var),strata_outcome_row[c("case_no","control_no","OR","LR","UR")],
                            hr_ci,auc_ci,strata_outcome_row[c("pval")]))
        
        
        forest_table<-rbind(forest_table,row_input,make.row.names=F)
        
      }

      
    }
    strata_outcome<-all_strata_outcome[[length(all_strata_outcome)]]
    strata_outcome_row<-strata_outcome[pgsname,]
    hr_ci<-paste(format(strata_outcome_row[1],nsmall=3),"(",
                 format(strata_outcome_row[2],nsmall=3),"-",
                 format(strata_outcome_row[3],nsmall=3),")",sep="")
    auc_ci<-paste(format(strata_outcome_row[4],nsmall=3),"(",
                  format(strata_outcome_row[5],nsmall=3),"-",
                  format(strata_outcome_row[6],nsmall=3),")",sep="")
    row_input<-unlist(c("Overall",strata_outcome_row[c("case_no","control_no","OR","LR","UR")],
                        hr_ci,auc_ci,strata_outcome_row[c("pval")]))
    
    forest_table<-rbind(forest_table,row_input,make.row.names=F)
    forest_table<-forest_table[2:nrow(forest_table),]
    forest_table$" "<- paste(rep(" ", 40), collapse = " ")
    colnames(forest_table)<-c(pgsname,"CAD cases", "CAD controls", "estimate","LR","UR"," ","AUC (95%CI)","p-value","   ")
    forest_table_list[[i]]<-forest_table
    
    
  }
  forest_table_list
}





#create_output_table(trainsplit = F, model_type = "simple",cat_or_cont = "cont",full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA")

#output_table_summary(outcome="prevalent_CHD_EPA",data=data,full_adjustments,trainsplit=F,train_data=NULL,test_data=NULL)

#output_table_summary(outcome="prevalent_CHD_EPA",data=data,full_adjustments,trainsplit=F,train_data=NULL,test_data=NULL)

#discrimination_without_prs(train_data = train,test_data = test,full_adjustments = full_adjustments,outcome="prevalent_CHD_EPA")

#create_quntile_char(data_prs, pgs,variables_to_select)