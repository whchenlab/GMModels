## my_CCM.adj ----
my_CCM.adj=function(feat,meta,top=NULL,models=NULL,method,
                    num.folds=5, num.resample=3,
                    check.con.before=NULL,check.con.after=NULL,
                    add_group=T,is_combat=T,re_scale=T,is_cross=T,
                    pca_plot=F,imbalance=F,max_index=3,batch_group=T,do_adj=T,
                    do.fs=F,nest_top=NULL,save_model=F,model_file=NULL){
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
  
  
  
  feat_list0 <- feat_list
  meta_list0 <- meta_list
  
  # model
  library(tidyverse)
  
  #1.data
  #feat
  # feat_list <- my_pair.table(feat_list)
  l0 <- length(feat_list)
  
  #more
  for (i in 2:(l0-2)){
    pai <- combn(1:l0,i)
    d=dim(pai)[2]
    lodo.feat <- list()
    lodo.meta <- list()
    for (j in 1:d){
      p_train <- names(feat_list0)[pai[,j]]
      feat_list.train <- feat_list0[p_train]
      meta_list.train <- meta_list0[p_train]
      feat.train<- reduce(feat_list.train,rbind)
      meta.train<- reduce(meta_list.train,rbind)
      t_model_lodo<-my_siamcat(train_feat=feat.train,train_meta=meta.train,method="lasso",num.folds=num.folds,num.resample=num.resample,feature.type = "normalized")
      if(save_model==T){
        save(t_model_lodo,file = paste0(model_file,paste(p_train,collapse = '_'),'.RData'))
      }
      p_r <- rbind(p_r,data.frame(AUC=t_model_lodo[["aucinfo"]],method='train',num=i,project.test=0,project.train=paste(p_train,collapse = ';')))
      
      p_test <- names(feat_list0)[-pai[,j]]
      feat_list.test <- feat_list0[p_test]
      meta_list.test <- meta_list0[p_test]
      for(t_p in names(feat_list.test)){
        t_model_test <- my_external.vali_siamcat(siamcat.obj=t_model_lodo,test_feat=feat_list.test[[t_p]],
                                                 test_meta=meta_list.test[[t_p]],plot=F,feature.type = "normalized")
        p_r <- rbind(p_r,data.frame(AUC=t_model_test[["aucinfo"]],method='test',num=i,project.test=t_p,project.train=paste(p_train,collapse = ';')))
      }
    }
  }
  return(p_r)
}


## my_CCM.adj.data ----
my_CCM.adj.data=function(feat,meta,top=NULL,models=NULL,method='lasso',
                         check.con.before=NULL,check.con.after=NULL,
                         add_group=T,is_combat=T,re_scale=T,is_cross=T,
                         pca_plot=F,imbalance=T,max_index=3,batch_group=F,do_adj=T,
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
        geom_point()+stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group)) 
      
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
  
  library(tidyverse)
  
  
  return(list(feat_list=feat_list,meta_list=meta_list))
}

## my_SCM.adj ----
my_SCM.adj  <- function(feat0,meta0,method='lasso',
                             model_file = NULL,
                             save_model=F) {
  
  data<-my_CCM.adj.data(feat=feat0,meta=meta0,top=NULL,models=NULL,method=method,
                        check.con.before=NULL,check.con.after=NULL,
                        add_group=T,is_combat=T,re_scale=T,is_cross=T,
                        pca_plot=F,imbalance=F,max_index=3,batch_group=F,do_adj=T,
                        do.fs=F,nest_top=NULL)
  feat_list0 <- data$feat_list
  meta_list0 <- data$meta_list
  l0<-length(feat_list0)
  project0 <- names(feat_list0)
  
  p_r <- data.frame('test.project'=0,'num' = 0,'AUC.test'=0,'AUC.train'=0)
  
  for (i in 1:l0){
    p <- project0[i]
    project1 <- project0[-i]
    feat_list.train <- feat_list0
    meta_list.train <- meta_list0
    feat_list.train[[p]] <- NULL
    meta_list.train[[p]] <- NULL
    feat_list.test <- feat_list0[[p]]
    meta_list.test <- meta_list0[[p]]
    
    feat.train<- reduce(feat_list.train,rbind)
    meta.train<- reduce(meta_list.train,rbind)
    len_max <- 2*min(sum(meta.train$Group=='Case'),sum(meta.train$Group=='Control'))
    num_case_all <- subset(meta.train,Group=='Case')
    num_control_all <- subset(meta.train,Group=='Control')
    
    s <- c(seq(16,34,6),seq(from=40, to=len_max, by=20))
    
    for (k in 1:10){
      for (j in s) {
        #train
        s1 <- sample(1:nrow(num_case_all),j/2)
        s2 <- sample(1:nrow(num_control_all),j/2)
        meta.train_t <- rbind(num_case_all[s1,],num_control_all[s2,])
        feat.train_t <- feat.train[rownames(meta.train_t),]
        t_model_lodo<-my_siamcat(train_feat=feat.train_t,train_meta=meta.train_t,method="lasso",num.folds=5,num.resample=3,feature.type = "normalized")
        
        #test
        test_list <- list()
        num_case_test <- subset(meta_list.test,Group=='Case')
        num_control_test <- subset(meta_list.test,Group=='Control')
        for (m in 1:10) {
          j_min <- min(nrow(num_case_test),nrow(num_control_test))
          s3 <- sample(1:nrow(num_case_test),j_min)
          s4 <- sample(1:nrow(num_control_test),j_min)
          meta.test_t <- rbind(num_case_test[s3,],num_control_test[s4,])
          feat.test_t <- feat_list.test[rownames(meta.test_t),]
          test_list[[m]] <- rownames(meta.test_t)
          auc_t<-my_external.vali_siamcat(siamcat.obj=t_model_lodo,test_feat=feat.test_t,test_meta=meta.test_t,plot=F,feature.type = "normalized")
          p_r <- rbind(p_r,data.frame('test.project'=p,'num' = j,'AUC.test'=auc_t$aucinfo,'AUC.train'=t_model_lodo[["aucinfo"]]))
        }
        test_list[[11]] <- rownames(meta_list.test)
        auc_t<-my_external.vali_siamcat(siamcat.obj=t_model_lodo,test_feat=feat_list.test,test_meta=meta_list.test,plot=F,feature.type = "normalized")
        p_r <- rbind(p_r,data.frame('test.project'=p,'num' = j,'AUC.test'=auc_t$aucinfo,'AUC.train'=t_model_lodo[["aucinfo"]]))
        model=list(model=t_model_lodo,test =test_list)
        if(save_model==T){
          save(model,file = paste0(model_file,p,'_',j,'_',k,'.RData'))
        }
      }
    }
  }
  return(p_r)
}
## my_CCM.adj_MSI ----
my_CCM.adj_MSI=function(feat,meta,marker_data,p_r){
  
  feat_list0 <- feat
  meta_list0 <- meta
  
  # model
  library(tidyverse)
  
  
  l0 <- length(feat_list0)
  
  for (i in 2:(l0-2)){
    pai <- combn(1:l0,i)
    d=dim(pai)[2]
    lodo.feat <- list()
    lodo.meta <- list()
    for (j in 1:d){
      p_test <- names(feat_list0)[-pai[,j]]
      feat_list.test <- feat_list0[p_test]
      meta_list.test <- meta_list0[p_test]
      #marker
      m <- my_marker.plot_new(markers_data=NULL,feat_list=feat_list.test,meta_list=meta_list.test,lda_cutoff=2,nproj_cutoff=1,level='genus',change_name=F,cut=F)
      m$class <- 'adjust'
      
      p_train <- names(feat_list0)[pai[,j]]
      feat_list.train <- feat_list0[p_train]
      meta_list.train <- meta_list0[p_train]
      feat.train<- reduce(feat_list.train,rbind)
      meta.train<- reduce(meta_list.train,rbind)
      #marker
      a<- my_marker.plot_new(markers_data=NULL,feat_list=list(lodo=feat.train),meta_list=list(lodo=meta.train),lda_cutoff=2,nproj_cutoff=1,level='genus',change_name=F,cut=F)
      a$class <- 'lodo'
      a0 <- a[rep(1:nrow(a),length(p_test)),]
      a0$project_id <- rep(p_test,each=nrow(a))
      
      marker_all <- rbind(a0,m)
      
      c <- my_cor.auc_number(marker_data=marker_all,rm_TR=F,method.dis="euclidean",project0=p_test)
      c <- c$mat
      c$TR <- paste(p_train,collapse = ';')
      for (t_p in 1:length(p_test)){
        p_r[(p_r$project.train %in% c$TR[t_p])&(p_r$project.test==c$TE[t_p]),'dis'] <- 1/c[c$TE==c$TE[t_p],'dis_value']
      }
    }
  }
  return(p_r)
}
