
## my_adj ----

# After simple data correction, the low abundance was removed, and the Meta difference test was performed, and significant confounders were used for correction
#' @param check.con.before,check.con.after Whether to plot the data using siamcat check.confounder
#' @param add_group=F Whether to add grouping information during calibration
#' @param is_combat T:limma combat; F:MMUPHin
#' @param re_scale=F Whether to normalize the revised data again

my_adj=function(feat_list,meta_list,check.con.before=NULL,check.con.after=NULL,add_group=T,is_combat=T,re_scale=T){
  library(SIAMCAT)
  library(limma)
  library(MMUPHin)
  
  
  meta_list=meta_list[,apply(meta_list, 2, function(x) sum(is.na(x)))<0.1*nrow(meta_list)]
  rm_id=unique(as.numeric(unlist(apply(meta_list, 2, function(x) c(which(is.na(x)))))))#要去除的样本
  
  if(length(rm_id)>0){
    message(sprintf('#### 1.There are %g samples meta-info is NA!',length(rm_id)))
    meta_list=meta_list[-rm_id,]
    feat_list=feat_list[-rm_id,]
  }
  meta_list=meta_list[,apply(meta_list, 2, function(x) sum(x=='NA'))<0.1*nrow(meta_list)]
  rm_id=unique(as.numeric(unlist(apply(meta_list, 2, function(x) c(which(x=='NA'))))))
  
  if(length(rm_id)>0){
    message(sprintf('#### There are %g samples meta-info is char NA!',length(rm_id)))
    meta_list=meta_list[-rm_id,]
    feat_list=feat_list[-rm_id,]
  }
  
  
  if(sum(meta_list$BMI==0)>0.1*nrow(meta_list)){
    meta_list$BMI=NULL
  }else{
    rm_id=which(meta_list$BMI==0)
    
    if(length(rm_id)>0){
      message(sprintf('#### There are %g samples meta-info is 0 in BMI!',length(rm_id)))
      meta_list=meta_list[-rm_id,]
      feat_list=feat_list[-rm_id,]
    }
  }
  #host_age
  if(sum(meta_list$host_age==0)>0.1*nrow(meta_list)){
    meta_list$host_age=NULL
  }else{
    rm_id=which(meta_list$host_age==0)
    
    if(length(rm_id)>0){
      message(sprintf('#### There are %g samples meta-info is 0 in age!',length(rm_id)))
      meta_list=meta_list[-rm_id,]
      feat_list=feat_list[-rm_id,]
    }
  }
  siamcat.train <- siamcat(feat=t(feat_list), meta=as.data.frame(meta_list),
                           label='Group', case='Case')
  
  siamcat.train <- filter.features(
    siamcat.train,
    filter.method = 'abundance',
    cutoff = 0.001,
    rm.unmapped = TRUE,
    verbose=2
  )
  
  label <- SIAMCAT::label(siamcat.train)
  meta <- meta(siamcat.train)
  feat <- get.filt_feat.matrix(siamcat.train)
  cases <- which(label$label == max(label$info))#1
  controls <- which(label$label == min(label$info))#-1
  
  if(length(which(apply(meta,2,function(x){length(unique(x))})>1))>0){
    var.level.names <- meta[,which(apply(meta,2,function(x){length(unique(x))})>1)] #选出不同取值超过1的变量
    meta=meta[,colnames(var.level.names)]
    
    if(!is.null(check.con.before)){
      check.confounders(siamcat.train,fn.plot = check.con.before)
    }
    
    confounder=c()
    for (m in seq_along(meta)) {
      
      mname <- colnames(meta)[m]
      message(paste("+++ checking",mname,"as a potential confounder"))
      mvar <- meta[[m]]
      
      u.val <- sort(unique(mvar)[!is.na(unique(mvar))])
      
      if (is.character(u.val)) {
        Group=ifelse(label$label==1,'Case','Control')
        ct=table(Group,mvar)
        p.val <- fisher.test(ct)$p.value
        cat('1.this is a contingency table info:','\n')
        print(ct)
        cat('2.this is a Fisher Test:','\n')
        print(paste0(mname,": Fisher Test P Value:", format(p.val, digits = 4)))
        if(p.val<0.05){confounder=c(confounder,mname)}
      } else {
        #numeric
        p.val <- wilcox.test(mvar[controls], mvar[cases],
                             exact = FALSE)$p.value
        cat('1.this is a Wilcox Test:','\n')
        print(paste0(mname,": MWW Test P Value:",format(p.val, digits = 4)))
        if(p.val<0.05){confounder=c(confounder,mname)}
      }
    }
    message('final confounder: ',paste(confounder,seq=';'))
    cat('over:finish find counfounders!','\n')
    
    if(!is.null(confounder)){
      #library(limma)
      covariates=intersect(confounder,c('host_age','BMI'))
      batch=intersect(confounder,c('disease_stage','country','sex'))
      meta=data.frame(meta)
      if(is_combat==F){
        if(length(covariates)==0){
          covariates=NULL
        }
        if(length(batch)==0){
          batch=NULL
          stop('+++ the project without confounder adjust about batch!')
        }
        
        meta=data.frame(meta)
        for (i in batch) {
          meta[,i]=factor(meta[,i])
        }
        
        feat_new=MMUPHin::adjust_batch(feature_abd = feat,
                                       batch = batch,
                                       covariates = covariates,
                                       data = meta,
                                       control = list(verbose = FALSE))
        
      }else{
        if(length(covariates)==0){
          covariates=NULL
        }else{
          covariates=meta[,covariates]
        }
        if(length(batch)==0){
          batch1=NULL
          batch2=NULL
        }else{
          if(length(batch)==1){
            batch1=meta[,batch]
            batch2=NULL
          }else{
            batch1=meta[,batch[1]]
            batch2=meta[,batch[2]]
          }
          #batch=meta[,batch]
          #batch=data.frame(batch)
          # if(ncol(batch)==1){
          #   batch=batch[,1]
          #   
          # }
        }
        
        if(add_group==T){
          modcombat = model.matrix(~as.factor(meta_list$Group))#设计矩阵
          feat_new=removeBatchEffect(feat,covariates=covariates,batch=batch1,batch2=batch2,design = modcombat)
        }else{
          feat_new=removeBatchEffect(feat,covariates=covariates,batch=batch1,batch2=batch2)
        }
      }
      
      if(re_scale==T){
        
        back_transform_abd <- function(adj_data, feature_abd, type_feature_abd) {
          normalize_features <- function(features,
                                         normalization = "NONE",
                                         pseudo_count = 0) {
            TSS <- function(x) {
              if(all(x == 0)) return(x)
              return(x / sum(x))
            }
            if(any(features < 0))
              stop("Feature table must be non-negative for normalization!")
            features <- features + pseudo_count
            if (normalization=='TSS')
              features <- apply(features, 2, TSS)
            if (normalization=='NONE')
              features <- features
            return(features)
          }
          adj_data <- 2^adj_data
          adj_data[feature_abd == 0] <- 0
          adj_data <- normalize_features(adj_data, normalization = "TSS")
          adj_data <- t(t(adj_data) * apply(feature_abd, 2, sum))
          dimnames(adj_data) <- dimnames(feature_abd)
          
          return(adj_data)
        }
        feat_new=back_transform_abd(feat_new,feat,'proportions')
      }
      
      filt_feat(siamcat.train)$filt.feat=feat_new
      if(!is.null(check.con.after)){
        check.confounders(siamcat.train,fn.plot = check.con.after)
      }
    }else{
      if(!is.null(check.con.before)){
        check.confounders(siamcat.train,fn.plot = check.con.before)
      }
      #message('over: there are no significantly different confounders to adjuest it in this project!')
    }
    
  }else{
    confounder=NULL
    message('over: there are no meta-variates whose freqency over 1!')
  }
  feat_new=as.data.frame(t(get.filt_feat.matrix(siamcat.train)))#将其转化为行为样本，列为特征
  list(feat_new=feat_new,meta_new=meta_list,confounder=confounder)
}

## my_result_adj ----

#' @param add_group=F,re_scale,is_combat were same as my_adj()
#' @param feature.type = "filtered" standardization when modeling
#' @param imbalance=F,max_index=3 add the sample imbalance adjustment parameter
#' @param do.con=T confounding factors
#' @param is_cross=T batch effect correction
#' @param bacth_group=T is added to the Group column as covariate
#' @param do.con=T adjustment for confounding factors 
#' @param do.fs=F nested top
my_result_adj=function(feat,meta,method,label,top=NULL,models=NULL,
                       num.folds=5, num.resample=3,new=NULL,auc.matrix=NULL,feature.type = "filtered",
                       check.con.before=NULL,check.con.after=NULL,
                       add_group=T,is_combat=T,re_scale=T,is_cross=F,
                       pca_plot=F,imbalance=F,max_index=3,batch_group=T,do.con=T,do.fs=F,nest_top=NULL){
  
  if(do.con==T){
    feat_all_adj=list()
    if(!is.null(check.con.before)){
      check.con.before=paste0(check.con.before,names(feat),'.pdf')
    }
    if(!is.null(check.con.after)){
      check.con.after=paste0(check.con.after,names(feat),'.pdf')
    }
    
    for(i in names(feat)){
      feat_all_adj[[i]]=my_adj(feat_list=feat[[i]],meta_list=meta[[i]],check.con.before[i],check.con.after[i],add_group,is_combat,re_scale)
    }
    feat=lapply(feat_all_adj,function(x) x$feat_new)
    confounder=lapply(feat_all_adj,function(x) x$confounder)
    meta_list=lapply(feat_all_adj,function(x) x$meta_new)
  }else{
    meta_list=meta
    
    feat=lapply(feat, function(x){
      max_p=apply(x,2,max)
      x=x[,max_p>0.001]
      message('this step filter out ',sum(max_p<=0.001),' features!')
      return(x)
    })
    
  }
  
  for (i in names(meta_list)) {
    cat('+++final model dataset:','\n')
    cat(i,'\n')
    cat(names(table(meta_list[[i]]$Group)),'\n')
    cat(table(meta_list[[i]]$Group),'\n')
  }
  
 
  feat_list=my_pair.table(feat)
  
  if(is_cross==F){
    
    results=my_result(feat_list,meta_list,method,label,top,models,num.folds, num.resample,
                      new,auc.matrix,feature.type)
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
    
    if(do.con==T){
      confounder=table(unlist(confounder))
      confounder=names(which(confounder==length(meta_list)))
      if(length(confounder)==0){
        confounder=NULL
        message('####cross-model: no con-confounder in batch-adjust!')
      }else{message('####cross-model: con-confounder:',confounder)}
    }else{
      confounder=NULL
      message('####cross-model: you choose no confounder!')
    }
    
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
    if(imbalance==T){
     
      temp2=list(feat_list=feat_list,meta_list=meta_list)
      ba_data=my_SMOTE_balance(temp2,max_index=max_index)
      feat_list=ba_data$feat_list
      meta_list=ba_data$meta_list
    }
    results=my_result(feat_list,meta_list,method,label,top,models,num.folds,num.resample,
                      new,auc.matrix,feature.type,do.fs=do.fs,nest_top=nest_top)
  }
  return(results)
}


## my_marker_adj ----

my_marker_adj=function(feat_list,meta_list,
                       check.con.before=NULL,check.con.after=NULL,
                       add_group=T,is_combat=T,re_scale=T,
                       is_plot=F,lda_cutoff=2,nproj_cutoff=1,
                       level='genus',change_name=T,do.con=T){
  m0=my_marker.plot_new(markers_data=NULL,feat_list,meta_list,
                    lda_cutoff,nproj_cutoff,level,change_name)
  
  if(do.con==T){
    feat_all_adj=list()
    if(!is.null(check.con.before)){
      check.con.before=paste0(check.con.before,names(feat_list),'.pdf')
    }
    if(!is.null(check.con.after)){
      check.con.after=paste0(check.con.after,names(feat_list),'.pdf')
    }
    
    for(i in names(feat_list)){
      feat_all_adj[[i]]=my_adj(feat_list[[i]],meta_list[[i]],check.con.before[i],check.con.after[i],add_group,is_combat,re_scale)
    }
    feat=lapply(feat_all_adj,function(x) x$feat_new)
    confounder=lapply(feat_all_adj,function(x) x$confounder)
    meta_list=lapply(feat_all_adj,function(x) x$meta_new)
    
    m1=my_marker.plot_new(markers_data=NULL,feat,meta_list,
                      lda_cutoff,nproj_cutoff,level,change_name)
    
    m0$class='origin'
    m1$class='adjust'
    m=rbind(m0,m1)
  }else{
    m1=m0
    m0$class='origin'
    m1$class='adjust'
    m=rbind(m0,m1)
  }
  
  if(is_plot==T){
    m1$project_id=paste0(m1$project_id,'_adj')
    mm=rbind(m0,m1)
    mm$project_id=factor(mm$project_id,levels = c(unique(m0$project_id),paste0(unique(m0$project_id),'_adj')))
    mm=mm[order(mm$project_id),]
    mm=my_marker.plot_new(mm,lda_cutoff = lda_cutoff,
                      nproj_cutoff = nproj_cutoff,level = level)
  }
  
  list(marker_data=m,confounder=confounder)
}

## my_cor.auc ----

#' @param is_adj whether to correct the data
#' @param cor_adj whether to the correlation coefficient
#' @param rm_TR=T do not plot trainset's auc
#' @param method.dis="euclidean" 

my_cor.auc=function(result,marker_data,is_adj=T,cor_adj=T,method,method.test,
                    rm_TR=F,method.dis="euclidean",is_lodo=F){
  library(reshape2)
  library(vegan)
  if(is_lodo==F){
    a=melt(result[,-2],id.vars = c('TR') ,value.name = 'auc',variable.name = 'TE')
    if(is_adj==T){
      b=subset(marker_data,class=='adjust')
    }else{
      b=marker_data
    }
    
    #mm$marker_data[,c('project_id','scientific_name','LDA')]
    b=acast(b,project_id~scientific_name,value.var = 'LDA',fill=0)
    c=apply(b,2,as.numeric)
    rownames(c)=rownames(b)
    
    cor_mat=data.frame()
    for(i in rownames(c)) {
      stopifnot(sum(abs(c[i,])>0)>0)
      if(sum(abs(c[i,])>0)==1){
        
        temp=data.frame(marker=c[,abs(c[i,])>0])
      }else{
        temp=c[,abs(c[i,])>0]
      }
      
      if(cor_adj==T){
        
        k=1/max(apply(c,1,function(x) sum(abs(x)>0)))
        c_i=sum(abs(c[i,])>0)*k
        t_i=cor(t(temp),method=method)[,i]*c_i
      }else{
        
        t_i=cor(t(temp),method=method)[,which(rownames(c) %in% i)] 
        names(t_i)<- rownames(temp)
        
      }
      
      t_i[is.na(t_i)]=0
      
      
      if (method.dis=='jaccard'){
        temp <- temp-min(temp)
        d_i <- vegdist(temp, method=method.dis)%>%as.matrix()
        d_i=1-d_i[,i]
      }
      if (method.dis=='euclidean'){
        d_i=dist(temp,diag=T,upper=T,method = method.dis)%>%as.matrix()
        d_i=d_i[,i]/ncol(temp)
      }
      cor_mat=rbind(cor_mat,data.frame(TR=i,TE=names(t_i),cor_value=as.numeric(t_i),dis_value=as.numeric(d_i),stringsAsFactors = F))
    }
  }else{
    
    pro=colnames(result)
    result$class=rownames(result)
    a=melt(result,id.vars = c('class') ,value.name = 'auc',variable.name = 'TE')
    a$TE=as.vector(a$TE)
    a$TR=a$TE
    
   
    a$TR=paste0(a$TR,'_lodo')
    
    a$TE[a$class=='self']=a$TR[a$class=='self']
    a=a[,c('TR','TE','auc','class')]
    
    b=marker_data#对应adjust,lodo
    b$project_id=paste(b$project_id,b$class,sep='_')
    b$project_id=gsub('_adjust','',b$project_id)
    b=acast(b,project_id~scientific_name,value.var = 'LDA',fill=0)
    c=apply(b,2,as.numeric)
    rownames(c)=rownames(b)
   
    cor_mat=data.frame()
    for(i in paste0(pro,'_lodo')) {
      stopifnot(sum(abs(c[i,])>0)>0)
      if(sum(abs(c[i,])>0)==1){
       
        temp=data.frame(marker=c[,abs(c[i,])>0])
      }else{
        temp=c[,abs(c[i,])>0]
      }
      temp=temp[c(i,gsub('_lodo','',i)),]
      
   
      if(cor_adj==T){

        
        k=1/max(apply(c,1,function(x) sum(abs(x)>0)))
        
        c_i=sum(abs(c[i,])>0)*k
        t_i=cor(t(temp),method=method)[,i]*c_i
      }else{
       
        t_i=cor(t(temp),method=method)[,i] 
      }
      
      t_i[is.na(t_i)]=0
      
      
      if (method.dis=='jaccard'){
        temp <- temp-min(temp)
        stopifnot(temp>=0)
        d_i <- vegdist(temp, method=method.dis)%>%as.matrix()
        d_i=1-d_i[,i]
      }
      if (method.dis=='euclidean'){
        d_i=dist(temp,diag=T,upper=T,method = method.dis)%>%as.matrix()
        d_i=d_i[,i]/ncol(temp)
      }
      cor_mat=rbind(cor_mat,data.frame(TR=i,TE=names(t_i),cor_value=as.numeric(t_i),dis_value=as.numeric(d_i),stringsAsFactors = F))
    }
  }
  
  mat=merge(a,cor_mat,by=c('TR','TE'))
  if(rm_TR==T){
    
    mat=filter(mat,TR!=TE)
  }
  
  #plot(mat$cor_value,mat$auc)
  P1=ggplot(data=mat, mapping=aes(x=cor_value, y=auc, colour=TR))+geom_line(size=1)+geom_point(size=3) 
  
  
  
  r = paste("cor: ",round(cor(mat$cor_value,mat$auc,method = method.test),3), sep = "")
  p_v = paste("p: ",round(cor.test(mat$cor_value,mat$auc)$p.value,3), sep = "")
  title1=paste(r," ",p_v, sep = "")
  P2=ggplot(mat, aes(x=cor_value, y=auc))+geom_point() + stat_smooth(method=lm)+labs(title = title1)
  print(P2)
  
  list(P1=P1,P2=P2,mat=mat)
}

## my_marker.stat ----

my_marker.stat=function(result=NULL,marker_data,is_adj=T,is_lodo=F){
  library(reshape2)
  if(is_lodo==F){
    if(!is.null(result)){
      a=melt(result[,-2],id.vars = c('TR') ,value.name = 'auc',variable.name = 'TE')
    }else{
      n=length(unique(marker_data$project_id))
      pro=unique(marker_data$project_id)
      a=data.frame(TR=as.vector(gl(n,1,n*n,labels=pro)),TE=as.vector(gl(n,n,n*n,labels=pro)))
    }
    if(is_adj==T){
      b=subset(marker_data,class=='adjust')
    }else{
      b=marker_data
    }
  }else{  
   
    if(!is.null(result)){
      pro=colnames(result)
      result$class=rownames(result)
      a=melt(result,id.vars = c('class') ,value.name = 'auc',variable.name = 'TE')
  
      a$TE=as.vector(a$TE)
      a$TR=a$TE
      
      a$TR=paste0(a$TR,'_lodo')
      a$TE[a$class=='self']=a$TR[a$class=='self']
      a=a[,c('TR','TE','auc','class')]
    }else{
      n=length(unique(marker_data$project_id))
      pro=unique(marker_data$project_id)
      pro.lodo=paste0(pro,'_lodo')
      a=data.frame(TR=rep(pro.lodo,length=2*length(pro.lodo)),TE=c(pro,pro.lodo))
    }
    
    b=marker_data
    b$project_id=paste(b$project_id,b$class,sep='_')
    b$project_id=gsub('_adjust','',b$project_id)
  }
  

  temp=cbind(a,com.m=NA,spe.m=NA,rev.m=NA,TR.m=NA,TE.m=NA)

  for(i in 1:nrow(a)){
    m_TR=subset(b,project_id==a$TR[i])[,c('LDA','scientific_name')]
    rownames(m_TR)=m_TR$scientific_name
    m_TE=subset(b,project_id==a$TE[i])[,c('LDA','scientific_name')]
    rownames(m_TE)=m_TE$scientific_name
    m_TR$LDA=as.numeric(m_TR$LDA)
    m_TE$LDA=as.numeric(m_TE$LDA)
    temp[i,'com.m']=length(intersect(m_TR$scientific_name,m_TE$scientific_name))

    temp[i,'spe.m']=length(setdiff(m_TE$scientific_name,m_TR$scientific_name))
    if(temp[i,'com.m']==0){
      temp[i,'rev.m']=0
    }else{
      
      a1=m_TR[intersect(m_TR$scientific_name,m_TE$scientific_name),'LDA']
      a2=m_TE[intersect(m_TR$scientific_name,m_TE$scientific_name),'LDA']
      temp[i,'rev.m']=sum(a1*a2<0)
    }
    
    temp[i,'TR.m']=nrow(m_TR)
    temp[i,'TE.m']=nrow(m_TE)
  }
  
  return(temp)
}

## my_adj_auc.mar ----

#' @param dir the path of the models

my_adj_auc.mar=function(dir='',file,data_type,taxon,disease,
                        method='pearson',method.test='pearson',
                        rm_TR = T,method.dis="euclidean",is_lodo=F,label0=NULL){
  if(is_lodo==T){
    load(paste0(dir,'lodo_',file))
    result=model.adj$result
    load(paste0(dir,'lodo_m_',file))
    marker_data=marker.adj$marker_data
    #1.cor
    cor1=my_cor.auc(result,marker_data,is_adj =T,
                    cor_adj=F,method,method.test,
                    rm_TR,method.dis,is_lodo=is_lodo)
    cor1$mat$data_type=data_type
    cor1$mat$taxon=taxon
    cor1$mat$disease=disease
    
    mar1=my_marker.stat(result=NULL,marker_data,is_adj =T,is_lodo=is_lodo)
  }else{
    load(paste0(dir,file))
    result=model.adj$result
    load(paste0(dir,'m_',file))
    marker_data=marker.adj$marker_data
    #1.cor
    cor1=my_cor.auc(result,marker_data,is_adj =T,
                    cor_adj=F,method,method.test,rm_TR,method.dis)
    cor1$mat$data_type=data_type
    cor1$mat$taxon=taxon
    cor1$mat$disease=disease
    
    mar1=my_marker.stat(result=NULL,marker_data,is_adj =T)
    
    
    # mar1$data_type=data_type
    # mar1$taxon=taxon
  }
  
  all_stat=merge(cor1$mat,mar1,by=c('TR','TE'))

  if(rm_TR==T){
    all_stat=filter(all_stat,TR!=TE)
  }
  
  all_stat$method=label0
  return(all_stat)
}

my_cor.auc_number=function(marker_data,rm_TR=F,method.dis="euclidean",project0){
  library(reshape2)
  library(vegan)
  b=marker_data
  pro=project0
  b$project_id=paste(b$project_id,b$class,sep='_')
  b$project_id=gsub('_adjust','',b$project_id)
  b=acast(b,project_id~scientific_name,value.var = 'LDA',fill=0)
  c=apply(b,2,as.numeric)
  rownames(c)=rownames(b)
  
  
  
  cor_mat=data.frame()
  for(j in pro) {
    i=paste0(j,'_lodo')
    stopifnot(sum(abs(c[i,])>0)>0)
    if(sum(abs(c[i,])>0)==1){
      
      temp=data.frame(marker=c[,abs(c[i,])>0])
    }else{
      temp=c[,abs(c[i,])>0]
    }
    temp=temp[c(i,gsub('_lodo','',i)),]
   
    if (method.dis=='jaccard'){
      temp <- temp-min(temp)
      stopifnot(temp>=0)
      d_i <- vegdist(temp, method=method.dis)%>%as.matrix()
      d_i=1-d_i[,i]
    }
    if (method.dis=='euclidean'){
      d_i=dist(temp,diag=T,upper=T,method = method.dis)%>%as.matrix()
      d_i=d_i[j,i]/ncol(temp)
    }
    cor_mat=rbind(cor_mat,data.frame(TR=i,TE=j,dis_value=as.numeric(d_i),stringsAsFactors = F))
  }
  mat=cor_mat
  list(mat=mat)
}

## my_marker.plot_new ----
#' @param change_name Connect to the GMrepo database and change NCBI name to taxon name. Here we only provide F.
my_marker.plot_new <- function(markers_data=NULL,feat_list=NULL,meta_list=NULL,lda_cutoff=2,nproj_cutoff=1,level='genus',change_name=F,cut=F){
  marker.identify.abundance <- function(feat_list,meta_list,level,lda_cutoff){
    library(microbiomeMarker)
    library(phyloseq)
    library(ggplot2)
    library(RMySQL)
    lefse.obj <- list()
    for (i in 1:length(feat_list)) {
      otu_mat <- feat_list[[i]]
      samples <- meta_list[[i]]
      
      otu_mat <- otu_mat[,colSums(otu_mat>0)>2]
      samples <- samples[rownames(otu_mat),]

      otu_mat <- t(otu_mat)
      tax_mat <- rownames(otu_mat)
      rownames(otu_mat) <- 1:length(tax_mat)
      otu_mat <- as.matrix(otu_mat)
      tax_mat <- as.matrix(tax_mat)
      rownames(tax_mat) <- 1:nrow(tax_mat)
      if(level=='genus'){
        colnames(tax_mat) <- "Genus"
      }
      if(level=='species'){
        colnames(tax_mat) <- "Species"
      }
      ## combined to phyloseq format 
      OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
      TAX = tax_table(tax_mat)
      samples = sample_data(samples)
      phylodat <- phyloseq(OTU, TAX, samples)
      #lefse
      lefse.obj[[i]] <- run_lefse(
        ps=phylodat, 
        group = "Group", 
        multigrp_strat = TRUE,
        lda_cutoff = lda_cutoff
      )
    }
    if(!is.null(names(feat_list))){
      names(lefse.obj) <- names(feat_list)
    }
    return(lefse.obj)
  }
  if(is.null(markers_data)){
    lda.plot <- data.frame()
    lefse.obj <- marker.identify.abundance(feat_list,meta_list,level,lda_cutoff)
    for (i in 1:length(lefse.obj)) {
      lda <- lefse.obj[[i]]@marker_table
      if(!is.null(lda)){
        lda$ef_lda[lda$enrich_group == 'Control'] <- -lda$ef_lda[lda$enrich_group == 'Control']
        lda <- as.matrix(lda)
        lda <- as.data.frame(lda)
        lda <- dplyr::select(lda,'feature','ef_lda','pvalue')
        lda <- cbind( lda, rep(names(lefse.obj)[i],times=nrow(lda)))
        colnames(lda) <- c('scientific_name','LDA','p_value','project_id')
      }else{
        lda <- data.frame()
      }
      
      lda.plot <- rbind( lda.plot, lda)
     
    }
    rownames(lda.plot) <- NULL
    nrproj <- rep(NA,times=nrow(lda.plot))
    for (i in 1:nrow(lda.plot)) {
      nrproj[i] <- sum(lda.plot$scientific_name[i]==lda.plot$scientific_name)
    }
    lda.plot <- cbind(lda.plot,nrproj)
    lda.plot.standard <- lda.plot
    set.seed(1)
    lda.plot.standard <- lda.plot.standard[sample(1:nrow(lda.plot.standard),nrow(lda.plot.standard)),]
    
    if(change_name){
      #change_name: g_ncbi_taxon_id to taxon name
    }
  }else{
    lda.plot.standard=markers_data
    
    if(cut==T){
      lda.plot.standard <- merge(lda.plot.standard[,-5],data.frame(table(lda.plot.standard$scientific_name)),by.x = "scientific_name",by.y="Var1")
      colnames(lda.plot.standard)[5] <- 'nrproj'
    }
  }
 
  library(rjson)
  jsondata <- toJSON( list( data = unname(split(lda.plot.standard, 1:nrow(lda.plot.standard))),
                            axisTickLabelFontSize = 13,
                            lda_cutoff=lda_cutoff ,
                            nproj_cutoff=nproj_cutoff )
  )
  
  library(r2d3)
  print( r2d3( data = jsondata, script = "00_model_FUN/marker_comparison_2.js", d3_version = "5") )
  
  return(lda.plot.standard)
}