### load model_FUN

project_file <- '/mnt/raid7/limin/project/01_2022_GMModel_lm&liu'
setwd(project_file)
{
  source('00_model_FUN/siamcat_models.R')
  source('00_model_FUN/siamcat_models_adj.R')
  source('00_model_FUN/my_lodo.R')
  source('00_model_FUN/my_lodo.adj.R')
}
load('./01_data/cohort_info.RData')

## intro-cohort & LODO: AUC and marker----

## modeling method

# method='enet'
# label='enet_adj_nc'
# method_dir <- '01enet'

method='lasso'
label='lasso_adj_nc'
method_dir <- '02lasso'

# method='ridge'
# label='ridge_adj_nc'
# method_dir <- '03ridge'

# method='randomForest'
# label='rf_adj_nc'
# method_dir <- '04rf'


pro1 <- unique(cross[,c('disease','num','data_type','taxon','group1','mesh')])
# 1:nrow(pro1)
for (i in 1:nrow(pro1)) {
  
  #input parameters
  project_stat0 <- subset(project_stat,disease==pro1$disease[i] & data_type==pro1$data_type[i])
  project_id0 <- project_stat0$project
  disease <- pro1[i,'disease'] 
  data_type <- pro1[i,'data_type'] 
  taxon <- pro1[i,'taxon'] 
  case_meshid <- pro1[i,'mesh'] #for example IBS:'D015212'
  control_meshid <- 'D006262'
  
  load(paste0('01_data/feat_meta/',disease,'_',data_type,'_',taxon,'.RData'))
  feat_list <- my_data0$feat_list
  meta_list <- my_data0$meta_list
  
  # single: AUC and marker
  model.adj <- my_result_adj(feat_list,meta_list,
                           method=method,label=label,top=NULL,models=NULL,
                           num.folds=5, num.resample=3,new=NULL,auc.matrix=NULL,
                           feature.type = "normalized",
                           add_group=T,is_combat=T,re_scale=T,is_cross=T,
                           pca_plot = F,imbalance = T,batch_group=F)
  save(model.adj,file = paste('04_model_raw_RData/single_LODO/',method_dir,'/',case_meshid,'_',data_type,'.',taxon,'_',length(project_id0),'.','RData',sep = ''))
  marker.adj=my_marker_adj(feat_list,meta_list,NULL,NULL,T,T,T,
                           is_plot=F,lda_cutoff=2,nproj_cutoff=1,level=taxon,change_name=F)
  save(marker.adj,file = paste('04_model_raw_RData/single_LODO/',method_dir,'/','m_',case_meshid,'_',data_type,
                               '.',taxon,'_',length(project_id0),'.','RData',sep = ''))
  #lodo: AUC and marker
  model.adj <- my_lodo.adj(feat=feat_list,meta=meta_list,top=NULL,models=NULL,method=method,
                           num.folds=10, num.resample=3,
                           check.con.before=NULL,check.con.after=NULL,
                           add_group=T,is_combat=T,re_scale=T,is_cross=T,
                           pca_plot=F,imbalance=F,max_index=3,batch_group=T,do_adj=T,
                           do.fs=F,nest_top=NULL)
  save(model.adj,file = paste('04_model_raw_RData/single_LODO/',method_dir,'/lodo_',case_meshid,'_',data_type,'.',taxon,'_',length(project_id0),'.','RData',sep = ''))
  marker.adj=my_lodo_marker.adj(feat_list,meta_list,level=taxon,is_change=F)
  save(marker.adj,file = paste('04_model_raw_RData/single_LODO/',method_dir,'/lodo_m_',case_meshid,'_',data_type,
                               '.',taxon,'_',length(project_id0),'.','RData',sep = ''))
}

# visualization
p <- my_auc.heatmap(datas=model.adj$result,models=label)
my_marker.plot(markers_data=marker.adj$marker_data,feat_list=NULL,meta_list=NULL,lda_cutoff=2,nproj_cutoff=2,level='species')


## nested feature select
#script33:重新跑嵌套特征选择[类比script31]
#2022.3.15 watson已成功外接VPN

## intro-cohort nested feature selection ----


## modeling method
method='lasso'
label='lasso_adj_nc'
method_dir <- '02lasso'

pro1 <- unique(cross[,c('disease','num','data_type','taxon','group1','mesh')])
pro1 <- subset(pro1,num>3)
# 1:nrow(pro1)
for (i in 1:nrow(pro1)) {
  model.adj <- list()
  for (num in c(11,15,20,25,30,35,40)){
  #input parameters
  project_stat0 <- subset(project_stat,disease==pro1$disease[i] & data_type==pro1$data_type[i])
  project_id0 <- project_stat0$project
  disease <- pro1[i,'disease'] 
  data_type <- pro1[i,'data_type'] 
  taxon <- pro1[i,'taxon'] 
  case_meshid <- pro1[i,'mesh'] #for example IBS:'D015212'
  control_meshid <- 'D006262'
  
  load(paste0('01_data/feat_meta/',disease,'_',data_type,'_',taxon,'.RData'))
  feat_list <- my_data0$feat_list
  meta_list <- my_data0$meta_list
  
  feat_list <- my_pair.table(feat_list)
  
  # all: 
  label1=paste0(method,'_all')
  results=my_result(feat=feat_list,meta=meta_list,method=method,label=label1,
                        top=NULL,models=NULL,
                        num.folds=5,num.resample=3,
                        new=NULL,auc.matrix=NULL)
  # top
  label2=paste0(method,'_top',top)
  results_top=my_result(feat=feat_list,meta=meta_list,method = method,
                        num.folds=5,num.resample=3,
                            label =label2,models = results$models,top = num,
                            new=NULL,auc.matrix=NULL)
  
  model.adj[as.character(num)] <- list(results_top$result)
  }
  save(model.adj,file = paste('04_model_raw_RData/single_LODO/02lasso/nest_top/',case_meshid,'_',data_type,'.',taxon,'_',length(project_id0),'.','RData',sep = ''))
}