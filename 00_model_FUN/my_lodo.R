
## my_lodo ----
#' @param feat_list
#' @param meta_list
#' @param top=NULL how many features are selected ,NULL means the featurs are not selected
#' @param models shuold be offered when top is not NULL(which is the result of my.lodo(..,top=NULL,models=NULL))

my_lodo <- function(feat_list,meta_list,top=NULL,models=NULL,
                    method="lasso",num.folds=5,num.resample=3,do.fs=F,nest_top=NULL){
  library(tidyverse)
  
  #1.data
  #feat
  project.list <- names(feat_list)
  all.feat <- my_pair.table(feat_list)
  lodo.feat <- list()
  for (project in project.list) {
    temp <- all.feat
    temp[[project]] <- NULL
    lodo.feat[[project]] <- temp
  }
  lodo.feat <- lapply(lodo.feat, function(data){reduce(data,rbind)})
  #meta
  lodo.meta <- list()
  for (project in project.list) {
    temp <- meta_list
    temp[[project]] <- NULL
    lodo.meta[[project]] <- temp
  }
  lodo.meta <- lapply(lodo.meta, function(data){reduce(data,rbind)})
  
  #2.modeling
  lodo.model <- list()
  lodo.predict <- list()
  lodo.model.top <- list()
  lodo.predict.top <- list()
  if(!is.null(top)){
    if(!is.null(models)){
      lodo.model <- models$model
      for (project in project.list){
        lodo.model.top[[project]] <- my_siamcat.sel(lodo.model[[project]]$model_back,
                                                    consens.thres=0.5,
                                                    max.show=top,
                                                    verbose = 1,method=method)
        lodo.predict.top[[project]] <- my_external.vali_siamcat(siamcat.obj=lodo.model.top[[project]],test_feat=all.feat[[project]],test_meta=meta_list[[project]],plot=F)
      }
    }else{
      stop('you must input models this time!')
    }
  }else{
    for (project in project.list) {
      #train_feat
      lodo.model[[project]] <- my_siamcat(train_feat=lodo.feat[[project]],train_meta=lodo.meta[[project]],method=method,num.folds=num.folds,num.resample=num.resample,do.fs=do.fs,nest_top=nest_top)
      lodo.predict[[project]] <- my_external.vali_siamcat(siamcat.obj=lodo.model[[project]],test_feat=all.feat[[project]],test_meta=meta_list[[project]],plot=F)
    }
  }
  mylodo=list(model=lodo.model, predict=lodo.predict, model.top=lodo.model.top, predict.top=lodo.predict.top)
  if(!is.null(top)){
    result <- as.data.frame(rbind(self=sapply(mylodo$model.top, function(x){return(x[["aucinfo"]])}),
                                  lodo=sapply(mylodo$predict.top, function(x){return(x[["aucinfo"]])})))
  }else{
    result <- as.data.frame(rbind(self=sapply(mylodo$model, function(x){return(x[["aucinfo"]])}),
                                  lodo=sapply(mylodo$predict, function(x){return(x[["aucinfo"]])})))
  }
  mylodo$result <- result
  return(mylodo)
}


