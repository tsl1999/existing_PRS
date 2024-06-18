plot_forestplot<-function(forest_table,col_input=c(1,2,8,3),ci_column=3,title=NULL,
                          fontsize=10,font="",summary_row=T,summary_row_in=NULL,
                          summary_fill="white",log=F,
                          title_cex=1,background_col="white",
                          ci_col="black",ci_fill="black",refline_col="black",
                          estimate_colname="hriwv",ur_colname="hriwv_uci",
                          lr_colname="hriwv_lci",se_colname="seiwv",
                          below_xaxis=NULL,footnote=NULL,ci_size=NULL,
                          xlim_min=NULL,xlim_max=NULL,
                          summary_size=1,overall_row=F,summary_size_overall=1,
                          ci_pch=15,ci_alpha=1,linear=F,ref_line=1,xjust=0.4,
                          x_ticks=NULL,heterogeneity=F){
  
  
  parameter<-setup_parameters(forest_table=forest_table,log=log,se_colname=se_colname,
                              estimate_colname=estimate_colname,ur_colname=ur_colname,
                              lr_colname=lr_colname,xlim_min=xlim_min,xlim_max=xlim_max,
                              ci_size=ci_size,fontsize=fontsize,font=font,
                              summary_fill=summary_fill,
                              title_cex=title_cex,background_col=background_col,
                              ci_col=ci_col,ci_fill=ci_fill,refline_col=refline_col,
                              ci_pch=ci_pch,ci_alpha=ci_alpha,linear=linear,xjust=xjust,
                              x_ticks=x_ticks)
  
  
  
  ci_size_in=parameter[[1]]
  est=parameter[[2]]
  lower=parameter[[3]]
  upper=parameter[[4]]
  xlim=parameter[[5]]
  ticks_at=parameter[[6]]
  x_trans=parameter[[7]]
  xlim_vert=parameter[[8]]
  tm <- parameter[[9]]
  
  forest_table_in<-forest_table[,col_input]
  
  
  blank_row<-which(forest_table[,c(3)]=="")
  
  
  row_highlight<- which(sapply(forest_table[,1],function(x) grepl("\U{00A0}\U{00A0}\U{00A0}\U{00A0}",x))+
                          sapply(forest_table[,1],function(x) grepl("\U{00A0}\U{00A0}",x))==0)
  if(summary_row==T){
    
    vertline_list<-setup_vertline(summary_row_in=summary_row_in,forest_table=forest_table,
                                  estimate_colname=estimate_colname,log=log,
                                  overall_row=overall_row)
    
    vertline<-vertline_list[[1]]
    vertline_up<-vertline_list[[2]]
    vertline_down<-vertline_list[[3]]
    summary_in<-vertline_list[[4]]
    summary_true<-which(summary_in==T)
    
    
    for(j in 1:length(estimate_colname)){
      for(x in 1: length(summary_true)){
        ci_size_in[[j]][summary_true[x]]<-0.005*(summary_size/as.numeric(forest_table[summary_true[x],se_colname[j]])^1.5)
      }
      if(overall_row==T){
        ci_size_in[[j]][length(ci_size_in[[j]])]<-0.001*(summary_size_overall/as.numeric(forest_table[nrow(forest_table),se_colname[j]])^1.5)
      }
      
    }}else if(summary_row==F){
      
      summary_in<-rep(FALSE, nrow(forest_table))
      
      if(overall_row==T){
        summary_in[length(summary_in)]<-TRUE
        for(j in 1:length(estimate_colname)){
          ci_size_in[[j]][length(ci_size_in[[j]])]<-0.001*(summary_size_overall/as.numeric(forest_table[nrow(forest_table),se_colname[j]])^1.5)
          
          
        }
        vertline_up<-1
        vertline_down<-nrow(forest_table_in)-1
        vertline<-c()
        for(k in 1:length(se_colname)){
          if(log==T){
            vertline[k]<-c(log(as.numeric(forest_table[c(nrow(forest_table_in)),estimate_colname[k]])))
          }else{
            vertline[k]<-c(as.numeric(forest_table[c(nrow(forest_table_in)),estimate_colname[k]]))}
          
        }}else{
          
          vertline_up<-NULL
          vertline_down<-NULL
          vertline<-NULL
        }}
  
  
  if(length(se_colname)==1){
    ci_size_in<-unlist(ci_size_in)
    est<-unlist(est)
    lower<-unlist(lower)
    upper<-unlist(upper)
    xlim<-unlist(xlim)
    ticks_at<-unlist(ticks_at)
    xlim_vert<-unlist(xlim_vert)
    
  }
  if(length(unique(ci_pch))==1){
    tm<-tm[[1]]
    p<-setup_forestplot(forest_table_in,est,lower,upper,ci_column,ci_size_in,
                        summary_in,xlim,xlim_vert,ticks_at,footnote,tm,title,x_trans,
                        row_highlight,blank_row=blank_row,below_xaxis=below_xaxis,
                        se_colname=se_colname,fontsize,vertline=vertline,
                        summary_row=summary_row,vertline_up=vertline_up, 
                        vertline_down=vertline_down,ref_line=ref_line)
    
  }else if (length(unique(ci_pch))>=2){
    g<-list()
    front_row<-ci_column[1]-1
    group<-(length(col_input)-front_row)/length(ci_pch)#can be changed to more columns if neeeded
    for ( x in 1:length(ci_pch)){
      forest_table_in_multi<-forest_table_in[,c(1:front_row,group*x+front_row-1,group*x+front_row)]
      g[[x]]<-setup_forestplot(forest_table_in=forest_table_in_multi,est=est[[x]],lower=lower[[x]],
                               upper=upper[[x]],
                               ci_column=ci_column[1],ci_size_in=ci_size_in[[x]],
                               summary_in,xlim=xlim[[x]],xlim_vert=xlim_vert[[x]],
                               ticks_at=ticks_at[[x]], footnote[x],
                               tm=tm[[x]],title,x_trans=x_trans[x],
                               row_highlight,blank_row=blank_row,below_xaxis=below_xaxis[x],
                               se_colname=se_colname[x],fontsize,vertline=list(vertline[[x]]),
                               summary_row=summary_row,vertline_up=vertline_up, 
                               vertline_down=vertline_down,ref_line=ref_line)
      if(x==1){
        p<-g[[x]]
      }else{
        p<-cbind(p,g[[x]][,3:ncol(g[[x]])])
      }
      
    }
  }
  
  
  # p<-edit_plot(p,col=right_align_col,which="text",
  #              hjust = unit(1, "npc"),
  #              x = unit(0.9, "npc"))
  if(heterogeneity==T){
    p<-edit_plot(p,row=c(blank_row-1),which="text",
                 gp=gpar(fontface="italic",fontsize=fontsize-2))
  }
  
  p
  
}

setup_parameters<-function(forest_table,log=F,se_colname="seiwv",
                           estimate_colname="hriwv",ur_colname="hriwv_uci",
                           lr_colname="hriwv_lci",xlim_min=NULL,xlim_max=NULL,
                           ci_size=NULL,
                           fontsize=10,font="",
                           summary_fill="white",
                           title_cex=1,background_col="white",
                           ci_col="black",ci_fill="black",refline_col="black",
                           ci_pch=15,ci_alpha=1,linear=F,xjust=4,x_ticks=NULL){
  ci_size_in<-list()
  est<-list()
  lower<-list()
  upper<-list()
  xlim<-list()
  ticks_at<-list()
  x_trans=c()
  xlim_vert<-list()
  tm<-list()
  for(j in 1: length(se_colname)){
    if(log==F){
      x_trans[j]<-"none"
    }else if (log==T){
      x_trans[j]<-"log"
    }
  }
  
  # if(log==T){
  #   forest_table[,c(estimate_colname,ur_colname,lr_colname)]<-apply(
  #     forest_table[,c(estimate_colname,ur_colname,lr_colname)],2,function(x) log(as.numeric(x)))
  # }
  if(length(se_colname)>1&length(ci_alpha)==1){
    ci_alpha<-rep(ci_alpha,length(se_colname))
  }
  
  if(length(se_colname)>1&length(ci_pch)==1){
    ci_pch<-rep(ci_alpha,length(se_colname))
  }
  
  for(j in 1:length(se_colname)){
    est[[j]] = as.numeric(forest_table[,estimate_colname[j]])
    lower[[j]] = as.numeric(forest_table[,lr_colname[j]])
    upper[[j]] = as.numeric(forest_table[,ur_colname[j]])
    tm[[j]] <- setup_theme(fontsize=fontsize,font=font,
                           summary_fill=summary_fill,log=log,
                           title_cex=title_cex,background_col=background_col,
                           ci_col=ci_col,ci_fill=ci_fill,refline_col=refline_col,
                           ci_pch=ci_pch[j],ci_alpha=ci_alpha[j],xjust=xjust)
    
    min_x<-max(round(min(as.numeric(forest_table[,lr_colname[j]]),na.rm=T),1)-0.1,0.1)
    # min_x<-min(floor(min(as.numeric(forest_table[,lr_colname[j]]),na.rm=T)),0.5)
    max_x<-max(ceiling(max(as.numeric(forest_table[,ur_colname[j]]),na.rm=T)),1.5)
    if(linear==T){
      max_x<-ceiling(max(as.numeric(forest_table[,ur_colname[j]]),na.rm=T))
    }
    max_x_check<-ceiling(max(as.numeric(forest_table[,ur_colname[j]]),na.rm=T))
    if(is.null(xlim_min)==T&is.null(xlim_max)==T){
      xlim[[j]]<-c(min(min_x),max_x)
      
    }else if(is.null(xlim_min)==F&is.null(xlim_max)==T){
      xlim[[j]]<-c(xlim_min[j],max_x)
    }else if(is.null(xlim_max)==F&is.null(xlim_min)==T){
      xlim[[j]]<-c(min(min_x),xlim_max[j])
      
    }else{xlim[[j]]<-c(xlim_min[j],xlim_max[j])
    }
    
    
    if(log==T){
      if(is.null(xlim_min)==F&is.null(xlim_max)==T){
        xlim_vert[[j]]<-c(log(xlim_min[j]),log(max_x))
      }else if (is.null(xlim_min)==T&is.null(xlim_max)==F){
        xlim_vert[[j]]<-c(log(min_x),log(xlim_max[j]))
      }else{
        xlim[[j]]<-xlim[[j]]
        xlim_vert[[j]]<-c(log(xlim[[j]][1]),log(xlim[[j]][2]))
      }
    }else{
      xlim_vert<-xlim
    }
    
    if(is.null(x_ticks)==T){
      ticks_at[[j]]<-unique(c(xlim[[j]][1],round(seq(min(1,xlim[[j]][1]),xlim[[j]][2],length.out=4),2),
                              xlim[[j]][2],1))}else{
                                if(length(se_colname)==1){
                                  ticks_at[[j]]<-unique(c(xlim[[j]][1],x_ticks, xlim[[j]][2],1))
                                }else{
                                  ticks_at[[j]]<-unique(c(xlim[[j]][1],x_ticks[[j]], xlim[[j]][2],1))
                                }
                                
                              }
    
    
    
    control_size<-max_x_check-min_x
    if(linear==T){
      ci_size_in_j<-0.1/sqrt(as.numeric(forest_table[,se_colname[j]]))
    }else{ci_size_in_j<-0.01*control_size/(as.numeric(forest_table[,se_colname[j]])^1.5)}
    
    ci_size_in[[j]]<-ci_size_in_j
    if(is.null(ci_size)==F){
      ci_size_in[[j]]<-ci_size*ci_size_in_j}}
  
  
  result=list(ci_size_in=ci_size_in,
              est=est,
              lower=lower,
              upper=upper,
              xlim=xlim,
              ticks_at=ticks_at,
              x_trans=x_trans,
              xlim_vert=xlim_vert,
              tm=tm)
  return(result)
  
}


setup_vertline<-function(summary_row_in=NULL,forest_table,
                         estimate_colname="hriwv",log=F,
                         overall_row=F){
  
  blank_row<-which(forest_table[,1]=="")
  if(is.null(summary_row_in)==T){
    summary_in<-rep(FALSE, nrow(forest_table))
    summary_in[c(blank_row-1,nrow(forest_table))]<-TRUE
    summary_true<-which(summary_in==T)
    
    vertline<-list()
    vertline_up<-c()
    vertline_down<-c()
    for(k in 1:length(estimate_colname)){
      if (overall_row==T){
        vertline_up<-c(1,blank_row+1)[1:length(blank_row)]
        vertline_down<- c(blank_row-2,nrow(forest_table)-1)[1:length(blank_row)]
        summary_true<-summary_true[1:(length(summary_true)-1)]
        if(log==T){
          vertline[[k]]<-c(log(as.numeric(forest_table[c(blank_row-1),estimate_colname[k]])))
        }else{vertline[[k]]<-c(as.numeric(forest_table[c(blank_row-1),estimate_colname[k]]))}
      }else{
        vertline_up<-c(1,blank_row+1)
        vertline_down<- c(blank_row-2,nrow(forest_table)-1)
        if(log==T){
          vertline[[k]]<-c(log(as.numeric(forest_table[c(blank_row-1,nrow(forest_table)),estimate_colname[k]])))
        }else{
          vertline[[k]]<-c(as.numeric(forest_table[c(blank_row-1,nrow(forest_table)),estimate_colname[k]]))
        }
      }
      
    }
    
    
  }else{
    summary_in<-rep(FALSE, nrow(forest_table))
    summary_in[c(summary_row_in)]<-TRUE
    summary_true<-summary_row_in
    
    #vertline<-c(as.numeric(forest_table[c(summary_row_in),estimate_colname]))
    vertline<-list()
    vertline_up<-c()
    vertline_down<-c()
    for(k in 1:length(estimate_colname)){
      if (overall_row==T){
        vertline_up<-c(1,summary_row_in[1:(length(summary_row_in)-2)]+1)
        vertline_down<-c(summary_row_in)[1:(length(summary_row_in)-1)]-1
        summary_true<-summary_true[1:(length(summary_true)-1)]
        if(log==T){
          vertline[[k]]<-c(log(as.numeric(forest_table[c(summary_row_in)[1:c(length(summary_row_in)-1)],estimate_colname[k]])))
        }else{
          vertline[[k]]<-c(as.numeric(forest_table[c(summary_row_in)[1:c(length(summary_row_in)-1)],estimate_colname[k]]))}
      }else{
        vertline_up<-c(1,summary_row_in[1:(length(summary_row_in)-1)]+1)
        vertline_down<-c(summary_row_in)-1
        if(log==T){
          vertline[[k]]<-c(log(as.numeric(forest_table[c(summary_row_in),estimate_colname[k]])))
        }else{
          vertline[[k]]<-c(as.numeric(forest_table[c(summary_row_in),estimate_colname[k]]))}}
    }
    
    
    
    for(j in 1:length(vertline_up)){
      if(vertline_up[j]%in%blank_row){
        vertline_up[j]<-vertline_up[j]+1
      }
    }
    
    
    
  }
  result<-list(vertline=vertline,
               vertline_up=vertline_up,
               vertline_down=vertline_down,
               summary_in=summary_in) 
  return(result)
  
}


make_summary1 <- function(est, lower, upper, gp, xlim, sizes) {
  # Return NULL if the CI is outside
  if(upper < min(xlim) || lower > max(xlim))
    return(NULL)
  
  polygonGrob(x = unit(c(lower, est, upper, est), "native"),
              y = unit(0.5 + c(0, sqrt(sizes)+0.05 , 0, -sqrt(sizes)), "npc"),
              gp = gp,
              vp = viewport(xscale = xlim),
              name = "pooled.diamond")
}



setup_forestplot<-function(forest_table_in,est,lower,upper,ci_column,ci_size_in,
                           summary_in,xlim,xlim_vert,ticks_at,footnote,tm,title,x_trans,
                           row_highlight,blank_row=NULL,below_xaxis=NULL,se_colname="seiwv",
                           fontsize,vertline=NULL,summary_row=T,vertline_up=NULL, 
                           vertline_down=NULL,ref_line=1){
  
  
  suppressWarnings(  
    p<-forest(forest_table_in,
              est = est,
              lower = lower, 
              upper = upper,
              ci_column = ci_column,
              sizes=ci_size_in,
              ref_line = ref_line,
              #arrow_lab = below_arrow,
              is_summary = summary_in,
              xlim = xlim,
              ticks_at = ticks_at,
              footnote = paste("\n\n\n\n\n\n\n\n ",footnote,sep=""),
              theme=tm,
              title = title,
              x_trans = x_trans,
              fn_summary = make_summary1)
  )
  
  p<- edit_plot(p, row = row_highlight, col=1,
                gp = gpar(fontface = "bold"))
  if(is.null(blank_row)==F&summary_row==T){
    p<- edit_plot(p, row = c(blank_row-1), 
                  col=c(1:ncol(forest_table_in)),
                  gp = gpar(fontface = "bold")) 
    # p<- edit_plot(p, row = c(blank_row-1,nrow(forest_table_in)), 
    #               col=c(1:ncol(forest_table_in)),
    #               gp = gpar(fontface = "bold")) 
  }
  
  if(is.null(below_xaxis)==F){
    for (x in 1:length(se_colname)){
      p<-add_text(p, text = below_xaxis[x],col = ci_column[x],row = nrow(forest_table_in)+2,
                  just = "center",
                  gp = gpar(fontface = "bold",fontsize=fontsize))
    }
    
  }
  
  if(is.null(unlist(vertline))==F){
    for (i in 1:length(vertline)){
      if(length(vertline)==1){
        xlim_in<-xlim_vert
      }else{
        xlim_in<-xlim_vert[[i]]
      }
      for(j in 1:length(vertline_up)){
        p <- add_grob(p, 
                      row = c(vertline_up[j]:vertline_down[j]), 
                      col = ci_column[i],
                      gb_fn = segmentsGrob,
                      gp = gpar(col = "black", lty = "dotted"),
                      x0 = unit(c(vertline[[i]][j]),"native"), # This is the lines you want to put
                      x1 = unit(c(vertline[[i]][j]),"native"), # Repeat it here to get a vertical line
                      y0 = unit(0.01,"npc"),
                      y1 = unit(0.99,"npc"),
                      vp = viewport(xscale = xlim_in))
      }
    }}
  p
}

setup_theme<-function(fontsize=10,font="",
                      summary_fill="white",log=F,
                      title_cex=1,background_col="white",
                      ci_col="black",ci_fill="black",refline_col="black",
                      ci_pch=15,ci_alpha=1,xjust=0.4){
  
  suppressWarnings( 
    
    tm <- forest_theme(base_size = fontsize,
                       core = list(bg_params=list(fill = c(background_col)),
                                   list(fg_params=list(hjust = unit(1,"npc"), 
                                                       x = unit(1,"npc")))),
                       colhead=list(fg_params=list(hjust=0,x=0)),
                       summary_col = "black",
                       summary_fill=summary_fill,
                       refline_lty = "solid",
                       ci_pch = ci_pch,
                       ci_col =ci_col,
                       ci_fill=ci_fill,
                       footnote_col = "black",
                       footnote_cex = 0.9,
                       vertline_lwd = 1,
                       vertline_lty = "dashed",
                       vertline_col = "grey20",
                       xaxis_cex=1.1,
                       ci_lty =1 ,
                       ci_lwd = 1,
                       ci_Theight = 0,
                       title_cex =title_cex,
                       ref_lty="solid",
                       refline_col = refline_col,
                       base_family = font ,
                       title_fontfamily = font,
                       ci_alpha=ci_alpha,
                       xlab_fontface = "bold"))
  tm
}



