### CCM SCM MSI

source('00_model_FUN/CCM & SCM Fun.R')

## CCM ----
{
  # ASD
  model <- 'D000067877_Amplicon.genus_6.RData'
  feat_meta_flie <- 'ASD_Amplicon_genus.RData'
  m_file <- '04_model_raw_RData/CCM&SCM/PD/model_CCM/'
  
  #PD
  model <- 'D010300_Amplicon.genus_6.RData'
  feat_meta_flie <- 'PD_Amplicon_genus.RData'
  m_file <- '04_model_raw_RData/CCM&SCM/PD/model_CCM/'
  
  # CRC
  model <- 'D015179_Metagenomics.species_7.RData'
  feat_meta_flie <- 'CRC_Metagenomics_species.RData'
  m_file <- '04_model_raw_RData/CCM&SCM/CRC/model_CCM/'
  
  # CD
  model <- 'D003424_Metagenomics.species_5.RData'
  feat_meta_flie <- 'CD_Metagenomics_species.RData'
  m_file <- '04_model_raw_RData/CCM&SCM/CD/model_CCM/'
} 


load(paste0('01_data/feat_meta/',feat_meta_flie))
feat_list <- my_data0$feat_list
meta_list <- my_data0$meta_list
feat_list=my_pair.table(feat_list)
project0 <- names(feat_list)

#one
load(paste('04_model_raw_RData/single_LODO/02lasso/',model,sep = ''))
for (i in 1:length(project0)){
  p=project0[i]
  if (i==1){
    p_r <- data.frame(AUC=model.adj[["result"]][,p][-which(p ==project0)],method='test',num=1,project.test=p,project.train=project0[-i])
  }else{
    p_r <- rbind(p_r,data.frame(AUC=model.adj[["result"]][,p][-which(p ==project0)],method='test',num=1,project.test=p,project.train=project0[-i]))
  }
}
r=model.adj[["result"]]
r=r[,3:ncol(r)]
r<-r[,rev(sequence(ncol(r)))]
p_r <- rbind(p_r,data.frame(AUC=rev(r[row(r) == NCOL(r) - col(r) + 1]),method='train',num=1,project.test=0,project.train=rownames(r)))

#most
l0 <- length(project0)
load(paste('04_model_raw_RData/single_LODO/02lasso/','lodo_',model,sep=''))
for (i in 1:length(project0)){
  p=project0[i]
  p_r <-  rbind(p_r,data.frame(AUC=model.adj[["result"]]['lodo',p],method='test',num=l0-1,project.test=project0[i],project.train=paste(project0[-i],collapse  = ';')))
  p_r <-  rbind(p_r,data.frame(AUC=model.adj[["result"]]['self',p],method='train',num=l0-1,project.test=0,project.train=paste(project0[-i],collapse  = ';')))
}

#more
# We suggest to run inside the function because the runtime was very long and the procession was easily broken.

p_r=my_CCM.adj(feat=feat_list,meta=meta_list,top=NULL,models=NULL,method='lasso',
               num.folds=10, num.resample=3,
               check.con.before=NULL,check.con.after=NULL,
               add_group=T,is_combat=T,re_scale=T,is_cross=T,
               imbalance=F,max_index=3,batch_group=F,p_r=p_r,save_model=F,model_file=m_file)


# save(p_r,file = '01_data/p_r_ASD.RData')
# save(p_r,file = '01_data/p_r_PD.RData')
# save(p_r,file = '01_data/p_r_CRC_species.RData')
# save(p_r,file = '01_data/p_r_CD_species.RData')



## SCM ----

#ASD
file1 <- 'D000067877_Amplicon.genus_6.RData'
load('01_data/feat_meta/ASD_Amplicon_genus.RData')
m_file <- '04_model_raw_RData/CCM&SCM/ASD/model_SCM/'

#PD
file1 <- 'D010300_Amplicon.genus_6.RData'
load('01_data/feat_meta/PD_Amplicon_genus.RData')
m_file <- '04_model_raw_RData/CCM&SCM/PD/model_SCM/'


#CRC
file1 <- 'D015179_Metagenomics.species_7.RData'
load('01_data/feat_meta/CRC_Metagenomics_species.RData')
m_file <- '04_model_raw_RData/CCM&SCM/CRC/model_SCM/'

#CD
file1 <- 'D003424_Metagenomics.species_5.RData'
load('01_data/feat_meta/CD_Metagenomics_species.RData')
m_file <- '04_model_raw_RData/CCM&SCM/CD/model_SCM/'


p_r <- my_SCM.adj(feat0=my_data0$feat_list,meta0=my_data0$meta_list,method='lasso',
                        model_file = m_file,
                        save_model=F) 
p_r <- get(load(('/mnt/raid7/limin/project/01_2022_GMModel_lm&liu/temp/p_r_PD_1_1.RData')))

p_r_temp <- data.frame()
for (i in names(my_data0$feat_list)){
  p_r_x <- subset(p_r,test.project ==i)
  print(table(subset(p_r_x,test.project==i)$num))
  for (j in unique(p_r_x$num)) {
    p_r_x_t <- subset(p_r_x,num==j)
    p_r_temp <- rbind(p_r_temp,p_r_x_t[cumsum(rep(11,10)),])
  }
}

table(p_r_temp$num)
p_r_temp <- subset(p_r_temp,num!=0)
p_r_temp <- p_r_temp[sort(p_r_temp$num,index.return =T)$ix,]
table(p_r_temp$num)
print(unique(p_r_temp$test.project))

for (p in unique(p_r_temp$test.project)) {
  print(table(subset(p_r_temp,test.project==p)$num))
}
p_r <- p_r_temp


# save(p_r,file = '01_data/p_r_ASD_lasso_genus_number_1_1.RData')
# save(p_r,file = '01_data/p_r_PD_lasso_genus_number_1_1.RData')
# save(p_r,file = '01_data/p_r_CRC_lasso_species_number_1_1.RData')
# save(p_r,file = '01_data/p_r_CD_lasso_species_number_1_1.RData')



## MSI ----

'D000067877_Amplicon.genus_6.RData'
'D010300_Amplicon.genus_6.RData'
'D015179_Metagenomics.species_7.RData'
'D003424_Metagenomics.species_5.RData'

#ASD


#PD


#CRC
cross_number_t <- cross[cross$num>4&cross$taxon=='species'&cross$disease=='CRC',]
file1 <- cross_number_t$models
marker_data_p <- get(load(paste0('04_model_raw_RData/single_LODO/02lasso/',cross_number_t$markers)))
load('01_data/p_r_CRC_species.RData')
load('01_data/feat_meta/CRC_Metagenomics_species.RData')


#CD
cross_number_t <- cross[cross$num>4&cross$taxon=='species'&cross$disease=='PD',]
file1 <- cross_number_t$models
marker_data_p <- get(load(paste0('04_model_raw_RData/single_LODO/02lasso/',cross_number_t$markers)))
load('01_data/p_r_CD_species.RData')
load('01_data/feat_meta/CD_Metagenomics_species.RData')




p_r$dis <- 0
data<-my_CCM.adj.data(feat=my_data0$feat_list,meta=my_data0$meta_list,top=NULL,models=NULL,method='lasso',
                      check.con.before=NULL,check.con.after=NULL,
                      add_group=T,is_combat=T,re_scale=T,is_cross=F,
                      pca_plot=F,imbalance=F,max_index=3,batch_group=F,do_adj=T,
                      do.fs=F,nest_top=NULL)
feat_list0 <- data$feat_list
meta_list0 <- data$meta_list
l0<-length(feat_list0)

project0 <- names(feat_list0)


#most
load(paste('04_model_raw_RData/single_LODO/02lasso/lodo_m_',file1,sep = ''))
c <- my_cor.auc_number(marker_data=marker.adj$marker_data,rm_TR=F,method.dis="euclidean",project0=project0)
c<-c[["mat"]]
# project.test


for (i in 1:length(project0)){
  p=project0[i]
  project.test_t=p
  project.train_t=paste(project0[-i],collapse  = ';')
  p_r[(p_r$project.test %in% project.test_t)&(p_r$project.train %in% project.train_t),'dis'] <- 1/c[c$TE==project.test_t,'dis_value']
}


#one
file <- cross_number_t$models
print(file)
disease <- cross_number_t$disease
data_type <- cross_number_t$data_type
taxon <- cross_number_t$taxon
## 01enet 02lasso 03ridge 04rf
c <- my_adj_auc.mar(dir='02lasso',file,data_type,taxon,disease,
                    method=method,method.test=method.test,
                    rm_TR = F,method.dis=method.dis,is_lodo=F)
c <- c[c$TR!=c$TE,]
if(method.dis=="euclidean"){
  c$dis_value <- 1/c$dis_value
}

for (i in 1:nrow(c)){
  p_r[(p_r$project.train==c$TR[i])&(p_r$project.test==c$TE[i]),'dis'] <- c$dis_value[i]
}


#more
p_r=my_CCM.adj_MSI(feat=feat_list0,meta=meta_list0,marker_data=marker_data_p$marker_data,p_r=p_r)

save(p_r,file='01_data/p_r_CD_marker_new_species.RData')


