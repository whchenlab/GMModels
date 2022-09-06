library(tidyverse)
library(SIAMCAT)
library(randomForest)
library(dplyr)
library(tidyr)
library(reshape2)  
library(ggplot2)
library(cowplot)
library(smotefamily)
library(RMySQL)


## my_SMOTE_balance ----

#' @param my_data0 list(feat_list,meta_list)
my_SMOTE_balance <- function(my_data0,max_index=3){
  library(smotefamily)
  feat_list0 <- my_data0$feat_list
  meta_list0 <- my_data0$meta_list
  name_p <- names(feat_list0) 
  for (i in 1:length(feat_list0)) {
    #deal with imbalance data
    feat_list <- feat_list0[[i]]
    meta_list <- meta_list0[[i]]
    L_case <- sum(meta_list$Group=='Case')
    L_control <- sum(meta_list$Group=='Control')
    L_max <- max(L_case,L_control)
    L_min <- min(L_case,L_control)
    if(L_max>ceiling(L_min*max_index)){
      cat(name_p[i],' data will be balanced\n')
      Group <- meta_list[rownames(feat_list),]$Group
      feat_list<- apply(feat_list,2,as.numeric)
      feat_list <- as.data.frame(feat_list)
      ##at least 3 times
      newData <- SMOTE(feat_list,Group,dup_size=max(max_index,floor(max(table(Group))/min(table(Group))/2)))
      feat_list <- newData$data
      sname <- paste('s',1:nrow(feat_list),sep = '')
      Group <- feat_list[,'class']
      Bodysite <- rep('stool',nrow(feat_list))
      rownames(feat_list) <- sname
      meta_list <- data.frame(Bodysite,Group,row.names = sname)
      feat_list <- feat_list[,-which('class'==colnames(feat_list))]
      my_data0$feat_list[[i]] <- feat_list
      my_data0$meta_list[[i]] <- meta_list
    }
  }
  return(my_data0)
}


## my_pair.table ----

#' @param feat_list 

my_pair.table=function(feat_list){
  library(tidyverse)
  
  
  feat_name=lapply(feat_list,function(data){colnames(data)})
  feat_ID=feat_name%>%unlist()%>%as.vector()%>%unique()
  b=lapply(feat_list, function(data){
    add=setdiff(feat_ID,colnames(data))
  })
  
  my_add=function(b,feat){
    if(length(b)!=0){
      b=as.character(b)
      x_add=data.frame(matrix(0,dim(feat)[1],length(b)))
      colnames(x_add)=b
      xtest=cbind(feat,x_add)
    }else{
      xtest=feat
    }
    return(xtest)
  }
  
  #add
  data_add=list()
  for (i in 1:length(feat_list)) {
    data_add[[i]]=my_add(b[[i]],feat_list[[i]])
  }
  names(data_add)=names(feat_list)
  
  feat_seq=sort(feat_ID)
  data_add=lapply(data_add, function(data){new_data=data[,c(feat_seq)]})
  
  return(data_add)
}

## my_siamcat ----

#input：
#' @param method: default is lasso
#' @param num.folds=5,num.resample=3
#' @param feature.type = "normalized"(when adjusting could be "filtered")


my_siamcat=function(train_feat,train_meta,method="lasso",num.folds=5,num.resample=3,feature.type = "normalized",do.fs=F,nest_top=NULL){
  library(SIAMCAT)
  library('randomForest')
  siamcat.train <- siamcat(feat=t(train_feat), meta=as.data.frame(train_meta),
                           label='Group', case='Case')
  
  siamcat.train <- filter.features(
    siamcat.train,
    filter.method = 'abundance',
    cutoff = 0.001,
    rm.unmapped = TRUE,
    verbose=2
  )
  
  if(feature.type == "normalized"){
    siamcat.train <- normalize.features(
      siamcat.train,
      norm.method = "log.std",
      norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1),
      verbose = 2
    )
  }
  siamcat.train <-  create.data.split(
    siamcat.train,
    num.folds = num.folds,
    num.resample = num.resample
  )
  
  if(do.fs==T){
    siamcat.train <- check.associations(siamcat.train, detect.lim = 1e-05,
                                        fn.plot = 'assoc_plot.pdf')
    
    siamcat.train<- train.model(
      siamcat.train,
      method = method,
      feature.type = feature.type, perform.fs = TRUE,
      param.fs = list(thres.fs = nest_top,method.fs = "AUC",direction='absolute')
    )
  }else{
    
      siamcat.train<- train.model(
        siamcat.train,
        method = method,
        feature.type = feature.type
      )
  }
  
  
  siamcat.train <- make.predictions(siamcat.train)
  siamcat.train <-  evaluate.predictions(siamcat.train)
  
  # model.evaluation.plot(siamcat.train)
  
  auc=as.numeric(eval_data(siamcat.train)$auroc)
  cat(auc,'\n')
  
  list(model_back=siamcat.train,aucinfo=auc)
}

## my_external.vali_siamcat ----

#' @param plot=F do not plot
my_external.vali_siamcat=function(siamcat.obj,test_feat,test_meta,plot=F,feature.type = "normalized"){
  library(SIAMCAT)
  siamcat.train=siamcat.obj$model_back
  siamcat.test <- siamcat(feat=t(test_feat), meta=as.data.frame(test_meta),
                          label='Group', case='Case')
  if(feature.type == "normalized"){
    siamcat.test <- normalize.features(siamcat.test,
                                       norm.param=norm_params(siamcat.train),
                                       feature.type='original',
                                       verbose = 2)
  }
  if(feature.type == "filtered"){
    siamcat.test <- filter.features(
      siamcat.test,
      filter.method = 'abundance',
      cutoff = 0,
      rm.unmapped = TRUE,
      verbose=2
    )
  }
  
  siamcat.test <- make.predictions(
    siamcat = siamcat.train,
    siamcat.holdout = siamcat.test,
    normalize.holdout = FALSE)

  siamcat.test <- evaluate.predictions(siamcat.test)
  
  auc=as.numeric(eval_data(siamcat.test)$auroc)
  cat(auc,'\n')
  
  if(plot==T){
    model.evaluation.plot('train'=siamcat.train,
                          'test'=siamcat.test,
                          colours=c('dimgrey', 'orange'))
  }
  
  list(model_back=siamcat.test,aucinfo=auc)
}


## my_siamcat.sel do not use----
  
#' @param siamcat siamcat object
#' @param consens.thres the proportion of coefficients that are not zero
#' @param max.show the number of top features
my_siamcat.sel=function(siamcat,
                        consens.thres=0.5,
                        max.show=20,
                        verbose = 1,method='lasso',feature.type = "normalized"){
  library(SIAMCAT)
  library('randomForest')
  siamcat.select.features <-function(
    siamcat,
    consens.thres,
    max.show,
    verbose = 1) {#,feature.type
    
    library(SIAMCAT)
    library('randomForest')
    
    label <- label(siamcat)
    model.type <- model_type(siamcat)
    feature.type <- feature_type(siamcat)
    models <- models(siamcat)
    feature.weights <- feature_weights(siamcat)
    weight.matrix <- weight_matrix(siamcat)
    if (verbose > 2) message("+ model.interpretation.select.features")
    
    if (model.type != "randomForest") {
      sel.idx = which(feature.weights$percentage > consens.thres)
      median.sorted.features <-
        sort(feature.weights$median.rel.weight[sel.idx],
             decreasing = TRUE,
             index.return = TRUE)
      
      if (length(sel.idx) > max.show) {
        warning(paste0("WARNING: restricting amount of features",
                       " to be plotted to ", max.show))
        median.sorted.features.abs <- sort(
          abs(feature.weights$median.rel.weight),
          decreasing = TRUE,
          index.return = TRUE)
        idx <- head(median.sorted.features.abs$ix, n = max.show)
        median.sorted.features <- sort(
          feature.weights$mean.rel.weight[idx],
          decreasing = TRUE,
          index.return = TRUE)
        sel.idx <- idx[median.sorted.features$ix]
      } else {
        warning(paste0("OK:An amount of features",
                       " to be plotted to ", length(sel.idx)))
        sel.idx = sel.idx[median.sorted.features$ix]
      }
    } else {
      median.sorted.features <-
        sort(feature.weights$median.rel.weight,
             decreasing = FALSE,
             index.return = TRUE)
      consens.thres=0
      sel.idx <-
        median.sorted.features$ix[which(median.sorted.features$x >=
                                          consens.thres)]
      
      sel.idx <- tail(sel.idx, n = max.show)
      warning(paste0("WARNING:Rf restricting amount of features",
                     " to be plotted to ", max.show))
    }
    
    if (verbose > 2)
      message(paste(
        "+++ generating plot for a model with",
        length(sel.idx),
        "selected features"
      ))
    if (verbose > 2)
      message("+ finished model.interpretation.select.features")
    
    if(feature.type == "normalized"){
      feat.norm<-get.norm_feat.matrix(siamcat)[sel.idx,]
      norm_feat(siamcat)[[1]] <- feat.norm
    }
    if(feature.type == "filtered"){
      feat.norm<-get.filt_feat.matrix(siamcat)[sel.idx,]
      filt_feat(siamcat)$filt.feat <- feat.norm
    }
    
    return(siamcat)
  }
  siamcat.train <- siamcat.select.features(siamcat,consens.thres,max.show,verbose)
  
  siamcat.train<- train.model(
    siamcat.train,
    method = method,feature.type=feature.type
  )
  
  siamcat.train <- make.predictions(siamcat.train)
  
  siamcat.train <-  evaluate.predictions(siamcat.train)
  auc=as.numeric(eval_data(siamcat.train)$auroc)
  cat(auc,'\n')
  
  list(model_back=siamcat.train,aucinfo=auc)
}

## my_result ----

#' @param top select features
#' @param models the models (when new is not NULL, there must be NA)
#' @param label used as a factor for visualization
#' @param new  project_id
#' @param auc.matrix (when new is not NULL, it should be original auc matriax)

my_result=function(feat,meta,method,label,top=NULL,models=NULL,num.folds=5, num.resample=3,
                   new=NULL,auc.matrix=NULL,feature.type = "normalized",consens.thres=0.5,do.fs=F,nest_top=NULL){
  if(!is.null(new)){
    message('you add new datas to train new models!')
    
    if(is.null(models)){stop('models.obj can not be NULL!')}
    if(is.null(auc.matrix)){stop('auc.matrix.obj can not be NULL!')}
    for (i in new) {
      if(is.na(models[[i]])){
        models[[i]]=my_siamcat(train_feat = feat[[i]],train_meta =meta[[i]],method,num.folds,num.resample,
                               feature.type=feature.type,do.fs=do.fs,nest_top=nest_top)
        
      }
    }
    
    models_top=models
    
    if(!is.null(top)){
      for(i in new){
        models_top[[i]]=my_siamcat.sel(siamcat=models_top[[i]]$model_back,
                                       consens.thres=consens.thres,
                                       max.show=top,
                                       verbose = 1,method=method,feature.type = feature.type)
      }
      
    }else{
      message("you don't choose top features!")
    }
    
    task_id=names(feat)
    add.id=c()
    for (k in new) {
      add.id=c(add.id,which(task_id==k))
    }
    
    cross.siamcat <- matrix(NA,length(task_id),length(task_id))
    cross_models=list()
    for (i in 1:length(task_id)) {
      for (j in 1:length(task_id)) {
        if(i==j){
          cross.siamcat[i,i]=models_top[[i]]$aucinfo
        }else{
          if(i%in%add.id | j%in%add.id){
            
            cross_models[[paste(task_id[i],task_id[j],sep='-')]]=my_external.vali_siamcat(models_top[[i]],
                                                                                          test_feat = feat[[j]],
                                                                                          test_meta = meta[[j]],
                                                                                          feature.type=feature.type)
            cross.siamcat[i,j]=cross_models[[paste(task_id[i],task_id[j],sep='-')]]$aucinfo
          }else{
            old.i=which(rownames(auc.matrix)==task_id[i])
            old.j=which(colnames(auc.matrix)==task_id[j])
            cross.siamcat[i,j]=auc.matrix[old.i,old.j]
          }
        }
      }
    }
    
    
  }else{
    if(is.null(models)){
      models=list()
      for (i in 1:length(feat)) {
        models[[i]]=my_siamcat(train_feat = feat[[i]],train_meta =meta[[i]],
                               method,num.folds,num.resample,feature.type=feature.type,do.fs=do.fs,nest_top=nest_top)
      }
    }
    models_top=models
    
    if(!is.null(top)){
      if(is.null(models)){
        stop('you must input models this time!')
      }
      models_top=lapply(models, function(data){
        data_top=my_siamcat.sel(siamcat=data$model_back,
                                consens.thres=consens.thres,
                                max.show=top,
                                verbose = 1,method=method,feature.type=feature.type)
      })
    }else{
      message("you don't choose top features!")
    }
    
    task_id=names(feat)
    cross.siamcat <- matrix(NA,length(task_id),length(task_id))
    cross_models=list()
    for (i in 1:length(task_id)) {
      for (j in 1:length(task_id)) {
        if(i==j){
          cross.siamcat[i,i]=models_top[[i]]$aucinfo
        }else{
          cross_models[[paste(task_id[i],task_id[j],sep='-')]]=my_external.vali_siamcat(models_top[[i]],
                                                                                        test_feat = feat[[j]],
                                                                                        test_meta = meta[[j]],
                                                                                        feature.type=feature.type)
          cross.siamcat[i,j]=cross_models[[paste(task_id[i],task_id[j],sep='-')]]$aucinfo
        }
      }
    }
  }
  
  tb = round(cross.siamcat, 2)
  tb=as.data.frame(tb)
  colnames(tb) = task_id
  rownames(tb) = task_id
  tb=cbind(TR=rownames(tb),method=label,tb)

  names(models)=names(feat)
  names(models_top)=names(feat)
  list(result=tb,models=models,models_top=models_top,cross_models=cross_models)
}


## my_auc.heatmap ----

#' @param models my_result() output result
#' @param models my_result() label

my_auc.heatmap=function(datas,models,save.png=F){
  
 
  library(dplyr)#tibble
  library(tidyr)#gather
  library(reshape2)  #melt
  library("tidyverse")
  library(ggplot2)
  library("cowplot")
  
  #col.scheme.heatmap=c("black","darkgreen","forestgreen","chartreuse3","lawngreen","yellow")
  #col.scheme.heatmap= c('gray98','steelblue1','midnightblue')
  col.scheme.heatmap= c('white','pink','red')

  proj=unique(as.character(datas$TR))
  
  data_df=melt(datas,id=c('TR','method'),variable.name='study.test',value.name = 'AUC')
  data_df=tibble(data_df)
  
 
  g <- data_df %>% 
    mutate(study.train=TR) %>%
    filter(method == models) %>%
    mutate(study.test=factor(study.test, levels=proj)) %>% 
    mutate(study.train=factor(study.train, levels=rev(proj))) %>% 
    mutate(CV=study.train == study.test) %>%
    ggplot(aes(y=study.train, x=study.test, fill=AUC)) +
    geom_tile() + theme_bw() +
    geom_text(aes_string(label="format(AUC, digits=2)"), col='black', size=7)+
    # color scheme
    scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0.4, 1)) +
    # axis position/remove boxes/ticks/facet background/etc.
    scale_x_discrete(position='top') + 
    theme(axis.line=element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.x.top = element_text(angle=45, hjust=.1), 
          panel.grid=element_blank(), 
          panel.border=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    xlab('Test Set') + ylab('Training Set') + 
    scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
    scale_size_manual(values=c(0, 3), guide=FALSE)
  
  
  g2 <- data_df %>% 
    mutate(study.train=TR) %>%
    filter(method == models) %>%
    filter(study.test != study.train) %>% 
    group_by(study.train) %>% 
    summarise(AUROC=mean(AUC)) %>% 
    mutate(study.train=factor(study.train, levels=rev(proj))) %>% 
    ggplot(aes(y=study.train, x=1, fill=AUROC)) + 
    geom_tile() + theme_bw() +
    geom_text(aes_string(label="format(AUROC, digits=2)"), col='black', size=7)+
    scale_fill_gradientn(colours = col.scheme.heatmap, limits=c(0.4, 1), 
                         guide=FALSE) + 
    scale_x_discrete(position='top') + 
    theme(axis.line=element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid=element_blank(), 
          panel.border=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + xlab('Model Average') + ylab('')
  
  plot0=plot_grid(g, g2, rel_widths = c(5/6, 2/6), align = 'h')
  print(plot0)
  if(save.png!=FALSE){
   
    save_plot(save.png, plot0)
  }
  return(plot0)
}

## my_marker.plot ----

#' @param lda_cutoff
#' @param nproj_cutoff No. of same marker
#' @param change_name Connect to the GMrepo database and change NCBI name to taxon name. Here we only provide F.
my_marker.plot <- function(markers_data=NULL,feat_list=NULL,meta_list=NULL,lda_cutoff=2,nproj_cutoff=1,level='genus',change_name=F,cut=F,order=F){
  marker.identify.abundance <- function(feat_list,meta_list,level,lda_cutoff){
    library(microbiomeMarker)
    library(phyloseq)
    library(ggplot2)
    library(RMySQL)
    lefse.obj <- list()
    for (i in 1:length(feat_list)) {
      otu_mat <- feat_list[[i]]
      samples <- meta_list[[i]]
      #otu_mat <- otu_mat[rowSums(otu_mat>0)>2,]
      otu_mat <- otu_mat[,colSums(otu_mat>0)>1]
      samples <- samples[rownames(otu_mat),]
      #make phyloseq
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
    lda.plot <- NULL
    lefse.obj <- marker.identify.abundance(feat_list,meta_list,level,lda_cutoff)
    for (i in 1:length(lefse.obj)) {
      lda <- lefse.obj[[i]]@marker_table
      lda$lda[lda$enrich_group == 'Control'] <- -lda$lda[lda$enrich_group == 'Control']
      lda <- as.matrix(lda)
      lda <- as.data.frame(lda)
      lda <- dplyr::select(lda,'feature','ef_lda','p_value')
      lda <- cbind( lda, rep(names(lefse.obj)[i],times=nrow(lda)))
      colnames(lda) <- c('scientific_name','LDA','p_value','project_id')
      if (!is.null(lda.plot)){
        lda.plot <- rbind( lda.plot, lda)
      }else{
        lda.plot <- lda
      }
    }
    rownames(lda.plot) <- NULL

    nrproj <- rep(NA,times=nrow(lda.plot))
    for (i in 1:nrow(lda.plot)) {
      nrproj[i] <- sum(lda.plot$scientific_name[i]==lda.plot$scientific_name)
    }
    lda.plot <- cbind(lda.plot,nrproj)

    lda.plot.standard <- lda.plot
    
    if(order==F){
      set.seed(1)
      lda.plot.standard <- lda.plot.standard[sample(1:nrow(lda.plot.standard),nrow(lda.plot.standard)),]
    }
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
  
  ## plot
  library(r2d3)
  print( r2d3( data = jsondata, script = "00_model_FUN/marker_comparison_2.js", d3_version = "5") )
  return(lda.plot.standard)
}


## genus.read ----

genus.read <- function(genus.data0){
  genus.data <- genus.data0[-grep('_$|uncultured$|Incertae_Sedis$|(gut metagenome$)|(uncultured bacterium$)|(uncultured organism$)',rownames(genus.data0),value=FALSE),]
  # rownames(genus.data) <- sapply(strsplit(rownames(genus.data),';'),"[",6)
  # rownames(genus.data) <- gsub("g__","",rownames(genus.data))
  genus.data <- genus.data/matrix(colSums(genus.data), nrow = nrow(genus.data) , ncol = ncol(genus.data) , byrow=TRUE)
  genus.data <- t(genus.data)
  return(genus.data)
}