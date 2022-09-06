
## my_lodo.adj ----
#' @param add_group=F,re_scale,is_combat were same as my_adj()
#' @param feature.type = "filtered" standardization when modeling
#' @param imbalance=F,max_index=3 add the sample imbalance adjustment parameter
#' @param do_adj=T   adjust confounding factors
#' @param add_group=T keep Group differences
#' @param is_cross=T batch effect correction
#' @param bacth_group=T is added to the Group column as covariate

my_lodo.adj=function(feat,meta,top=NULL,models=NULL,method,
                     num.folds=5, num.resample=3,
                     check.con.before=NULL,check.con.after=NULL,
                     add_group=T,is_combat=T,re_scale=T,is_cross=T,
                     pca_plot=F,imbalance=F,max_index=3,batch_group=T,do_adj=T,
                     do.fs=F,nest_top=NULL){
  if(do_adj==T){
    feat_all_adj=list()
    if(!is.null(check.con.before)){
      check.con.before=paste0(check.con.before,names(feat),'.pdf')
    }
    if(!is.null(check.con.after)){
      check.con.after=paste0(check.con.after,names(feat),'.pdf')
    }
    
    for(i in names(feat)){
      feat_all_adj[[i]]=my_adj(feat[[i]],meta[[i]],check.con.before[i],check.con.after[i],add_group,is_combat,re_scale)
    }
    feat=lapply(feat_all_adj,function(x) x$feat_new)
    confounder=lapply(feat_all_adj,function(x) x$confounder)
    meta_list=lapply(feat_all_adj,function(x) x$meta_new)
    for (i in names(meta_list)) {
      cat('+++final model dataset:','\n')
      cat(i,'\n')
      cat(names(table(meta_list[[i]]$Group)),'\n')
      cat(table(meta_list[[i]]$Group),'\n')

    }
    

    feat_list=my_pair.table(feat)
    
  }else{
    feat_list=my_pair.table(feat)
    meta_list=meta
  }
  
  if(is_cross==F){
    
    meta_list=lapply(meta_list,function(x){
      if(ncol(x)==7){
        return(x[,c('Group','Bodysite','disease_stage','country','sex','host_age','BMI')])
      }else{
        x[,setdiff(c('Group','Bodysite','disease_stage','country','sex','host_age','BMI'),colnames(x))]=NA
        return(x[,c('Group','Bodysite','disease_stage','country','sex','host_age','BMI')])
      }
      
    })
  }else{
    meta_list=lapply(meta_list,function(x){
      if(ncol(x)==7){
        return(x[,c('Group','Bodysite','disease_stage','country','sex','host_age','BMI')])
      }else{
        x[,setdiff(c('Group','Bodysite','disease_stage','country','sex','host_age','BMI'),colnames(x))]=NA
        return(x[,c('Group','Bodysite','disease_stage','country','sex','host_age','BMI')])
      }
      
    })
    for (i in 1:length(meta_list)) {
      meta_list[[i]][,'batch']=paste0('pro',i)
      meta_list[[i]]$batch=as.factor(meta_list[[i]]$batch)
    }
    meta_all=do.call(rbind,meta_list)
    feat_all=do.call(rbind,feat_list)
    library(MMUPHin)
    confounder=table(unlist(confounder))
    confounder=names(which(confounder==length(meta_list)))
    if(length(confounder)==0){
      confounder=NULL
      message('####cross-model: no con-confounder in batch-adjust!')
    }else{message('####cross-model: con-confounder:',confounder)}
    if(batch_group==T){
      abd_adjust_batch <- adjust_batch(feature_abd = t(feat_all),
                                       batch = "batch",
                                       covariates = c("Group",confounder),
                                       data = meta_all,
                                       control = list(verbose = FALSE))$feature_abd_adj
      
    }else{
      abd_adjust_batch <- adjust_batch(feature_abd = t(feat_all),
                                       batch = "batch",
                                       covariates = c(confounder),
                                       data = meta_all,
                                       control = list(verbose = FALSE))$feature_abd_adj
      
    }
    
    if(pca_plot==T){
      library(vegan, quietly = TRUE)
      D_before <- vegdist(feat_all,method='euclidean')
      D_after <- vegdist(t(abd_adjust_batch),method='euclidean')
      
      set.seed(1)
      fit_adonis_before <- adonis(D_before ~ batch, data = meta_all)
      fit_adonis_after <- adonis(D_after ~ batch, data = meta_all)
      message('#1.fit_adonis_before: (euclidean)')
      print(fit_adonis_before)
      message('#2.fit_adonis_after: (euclidean)')
      print(fit_adonis_after)
      
      P1=prcomp(t(abd_adjust_batch),scale=F)
      P2=prcomp(feat_all,scale=F)
      
      plot_data=P1$x[,c('PC1','PC2')]
      plot_data=cbind(plot_data,meta_all)
      plot_data0=P2$x[,c('PC1','PC2')]
      plot_data0=cbind(plot_data0,meta_all)
      P1=ggplot(plot_data0, aes(x=PC1, y=PC2, color=Group, shape=batch)) +
        geom_point()+stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
      
      P2=ggplot(plot_data, aes(x=PC1, y=PC2, color=Group, shape=batch)) +
        geom_point()+stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group)) #+ geom_text_repel(aes(label=subjectID), show.legend = F)
      
      P=plot_grid(P1,P2,labels=c("A", "B"), ncol = 1)
      print(P)
    }
    
    
    d=abd_adjust_batch%>%t()%>%as.data.frame()
    for (i in names(feat_list)) {
      feat_list[[i]]=d[paste0(i,'.',rownames(feat_list[[i]])),]
      rownames(meta_list[[i]])=paste0(i,'.',rownames(meta_list[[i]]))
    }
  }
  if(imbalance==T){
    temp2=list(feat_list=feat_list,meta_list=meta_list)
    ba_data=my_SMOTE_balance(temp2,max_index=max_index)
    feat_list=ba_data$feat_list
    meta_list=ba_data$meta_list
  }
  my_lodo.adj<-my_lodo(feat_list,meta_list,top=top,models=models,
                       method=method,num.folds=num.folds,num.resample=num.resample,
                       do.fs=do.fs,nest_top=nest_top)
  return(my_lodo.adj)
}

## my_lodo_marker.adj ----
#' @param level genus or species
#' @param is_change=T whether to compare the NCBI ID with the bacteria name
#' @param do_adj=T adjust confounding factors
my_lodo_marker.adj <- function(feat,meta,level,do_adj=T,is_change=T){
  library(tidyverse)
  if(do_adj==T){
    feat_all_adj=list()
    
    for(i in names(feat)){
      feat_all_adj[[i]]=my_adj(feat[[i]],meta[[i]],NULL,NULL,add_group=T,is_combat=T,re_scale=T)
    }
    
    feat=lapply(feat_all_adj,function(x) x$feat_new)
    confounder=lapply(feat_all_adj,function(x) x$confounder)
    meta_list=lapply(feat_all_adj,function(x) x$meta_new)
  }else{
    feat=feat
    meta_list=meta
  }
  for (i in names(meta_list)) {
    cat('+++final model dataset:','\n')
    cat(i,'\n')
    cat(names(table(meta_list[[i]]$Group)),'\n')
    cat(table(meta_list[[i]]$Group),'\n')
  }
  
  feat_list=my_pair.table(feat)
  m0=my_marker.plot_new(markers_data=NULL,feat,meta_list,lda_cutoff = 2,
                    nproj_cutoff = 1,level = level,change_name = is_change)
  m0$class='adjust'
  
  project.list <- names(feat_list)
  all.feat <- feat_list
  lodo.feat <- list()
  for (project in project.list) {
    temp <- all.feat
    temp[[project]] <- NULL
    lodo.feat[[project]] <- temp
  }
  lodo.feat <- lapply(lodo.feat, function(data){reduce(data,rbind)})
  
  meta_list=lapply(meta_list,function(x){
    if(ncol(x)==7){
      return(x[,c('Group','Bodysite','disease_stage','country','sex','host_age','BMI')])
    }else{
      x[,setdiff(c('Group','Bodysite','disease_stage','country','sex','host_age','BMI'),colnames(x))]=NA
      return(x[,c('Group','Bodysite','disease_stage','country','sex','host_age','BMI')])
    }
    
  })
  
  lodo.meta <- list()
  for (project in project.list) {
    temp <- meta_list
    temp[[project]] <- NULL
    lodo.meta[[project]] <- temp
  }
  lodo.meta <- lapply(lodo.meta, function(data){reduce(data,rbind)})

  m1=my_marker.plot_new(markers_data=NULL,lodo.feat,lodo.meta,lda_cutoff = 2,
                    nproj_cutoff = 1,level = level,change_name = is_change)
  m1$class='lodo'
  m=rbind(m0,m1)
  
  m1$project_id=paste0(m1$project_id,'_lodo')
  x=my_marker.plot_new(markers_data=rbind(m0,m1),lda_cutoff = 2,
                   nproj_cutoff = 1,level = level)
  
  list(marker_data=m,confounder=confounder)
}
