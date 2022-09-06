### F0 disposal plot data

## cohort_info ----
project_stat <- read.csv('01_data/all_models_data.csv',header  = T)
colnames(project_stat) <- c('project',"mesh","disease","case","control","taxon","data_type","num")
project_stat$samples <- project_stat$case + project_stat$control
project_stat <- rbind(subset(project_stat,data_type=='Amplicon'&taxon=='genus'),
                      subset(project_stat,data_type=='Metagenomics'&taxon=='species'))
cross <- read.csv('01_data/disease_num.csv',header  = T)
colnames(cross)[4] <- 'group1'
group2 <- cross$group1
group2[group2!='Intestinal'] <- 'non-Intestinal'
markers <- paste('m_',cross$mesh,'_',cross$data_type,'.',cross$taxon,'_',cross$count,'.RData',sep='')
models <- paste(cross$mesh,'_',cross$data_type,'.',cross$taxon,'_',cross$count,'.RData',sep='')
cross <- data.frame(disease=cross$disease,data_type=cross$data_type,taxon=cross$taxon,num=cross$count,mesh=cross$mesh,markers=markers,models=models,group1=cross$group1,group2=group2) 
pro=unique(cross[,c('disease','num','data_type','group1')])
save(cross,pro,project_stat,file = './01_data/cohort_info.RData')
# load('./01_data/cohort_info.RData')

## single: AUC MSI ----
method.dis="euclidean"
method.m ='median'
method='spearman'
method.test='spearman'
auc_marker <- data.frame(auc=NULL,dis=NULL,level=NULL,disease=NULL)
auc0 <- data.frame(TR=NULL,TE=NULL,auc=NULL,disease=NULL,data_type=NULL,taxon=NULL,group1=NULL,group2=NULL)
MSI_all <- data.frame()

for (i in 1:nrow(cross)){
  if (cross$num[i]>1){
    file <- cross$models[i]
    print(file)
    disease <- cross$disease[i]
    data_type <- cross$data_type[i]
    taxon <- cross$taxon[i]
    ## 01enet 02lasso 03ridge 04rf
    c <- my_adj_auc.mar(dir='/mnt/raid2/liujinxin/datarepo2021/models_adj/02lasso/',file,data_type,taxon,disease,
                        method=method,method.test=method.test,
                        rm_TR = F,method.dis=method.dis,is_lodo=F)
    auc0 <- rbind(auc0,cbind(c[,c('TR','TE','auc')],method=rep('external',nrow(c)),cross[i,c('disease','data_type','taxon','group1','group2')][rep(1,times=nrow(c)),]))
    auc0$method[auc0$TR==auc0$TE] <- 'self'
    
    
    c <- c[c$TR!=c$TE,]
    c1 <- c
    
    if(method.dis=="euclidean"){
      c1$dis_value <- 1/c1$dis_value
    }
    c1$group1 <- cross$group1[i]
    c1$group2 <- cross$group2[i]
    MSI_all <- rbind(MSI_all,c1)
    
    
    if(method.dis=="euclidean"){
      c$dis_value <- 1/c$dis_value
    }
    if(method.m =='median'){
      auc_marker <- rbind(auc_marker,data.frame(auc=median(c$auc),dis = median(c$dis_value),level=paste(data_type,taxon,sep = '_'),disease=disease))
    }
    if(method.m =='mean'){
      auc_marker <- rbind(auc_marker,data.frame(auc=mean(c$auc),dis = mean(c$dis_value),level=paste(data_type,taxon,sep = '_'),disease=disease))
    }
    
  }
}
auc_marker <- cbind(auc_marker,disease_type1=cross$group1,disease_type2=cross$group2)
auc0$level <- paste(auc0$data_type,auc0$taxon,sep = '_')
auc_all <- auc0

setwd(project_file)

save(auc_marker,file = '01_data/AUC_MSI_euclidean_median.RData')    
save(auc_all,file = '01_data/AUC_all_lasso.RData')
save(MSI_all,file ='01_data/MSI_all.RData')

## LODO: AUC MSI ----
cross0 <- cross[cross$num>2,]
cross0$level <- paste(cross0$data_type,cross0$taxon,sep = '_')
a1=get(load('01_data/AUC_all_lasso.RData'))
auc_marker_lodo <- data.frame(auc=NULL,dis=NULL,level=NULL,disease=NULL)
for (i in 1:nrow(cross0)){
  if (cross0$num[i]>1){
    file <- cross0$models[i]
    print(file)
    disease0 <- cross0$disease[i]
    print(disease0)
    data_type <- cross0$data_type[i]
    taxon <- cross0$taxon[i]
    auc0 <- subset(a1,disease==disease0&level==cross0$level[i])
    c0 <- my_adj_auc.mar(dir='02lasso',file,data_type,taxon,disease0,
                         method=method,method.test=method.test,
                         rm_TR = F,method.dis=method.dis,is_lodo=T)
    c <- c0[c0$TR==c0$TE,]
    auc0 <- rbind(auc0,cbind(c[,c('TR','TE','auc')],method=rep('lodo_self',nrow(c)),cross0[i,c('disease','data_type','taxon','group1','group2','level')][rep(1,times=nrow(c)),]))
    c <- c0[c0$TR!=c0$TE,]
    auc0 <- rbind(auc0,cbind(c[,c('TR','TE','auc')],method=rep('lodo_external',nrow(c)),cross0[i,c('disease','data_type','taxon','group1','group2','level')][rep(1,times=nrow(c)),]))
    if(i==1){
      auc_lodo <- auc0
    }else{
      auc_lodo <- rbind(auc_lodo,auc0)
    }
    
    
    if(method.dis=="euclidean"){
      c$dis_value <- 1/c$dis_value
    }
    # msi was NA meaning that the dataset hadn't marker
    c$dis_value[is.na(c$dis_value)] <- 0
    
    
    if(method.m =='median'){
      auc_marker_lodo <- rbind(auc_marker_lodo,data.frame(auc=median(c$auc),dis = median(c$dis_value),level=paste(data_type,taxon,sep = '_'),disease=disease0))
    }
    if(method.m =='mean'){
      auc_marker_lodo <- rbind(auc_marker_lodo,data.frame(auc=mean(c$auc),dis = mean(c$dis_value),level=paste(data_type,taxon,sep = '_'),disease=disease0))
    }
    if(method.m =='entropy'){
      auc_marker_lodo <- rbind(auc_marker_lodo,data.frame(auc=mean(c$auc),dis = entropy(c$dis_value,method="ML",unit = "log2"),level=paste(data_type,taxon,sep = '_'),disease=disease0))
    }
  }
}

setwd(project_file)

save(auc_lodo,file='01_data/lodo_AUC_all_lasso.RData')
auc_marker_lodo <- cbind(auc_marker_lodo,disease_type1=cross0$group1,disease_type2=cross0$group2)
save(auc_marker_lodo,file='01_data/AUC_MSI_lodo_lasso.RData')


## LODO single (median) ----
auc_all <- get(load('./01_data/lodo_AUC_all_lasso.RData'))

self.e <- subset(auc_all,method %in%  c("external"))
p=unique(self.e[,c('disease','TE','level')])
self.lodo.e <- subset(auc_all,method %in%  c("lodo_external"))
for (i in 1:nrow(p)) {
  auc_temp  <- subset(self.e,disease==p[i,'disease']&TE==p[i,'TE']&level==p[i,'level'])
  auc_lodo_temp  <- subset(self.lodo.e,disease==p[i,'disease']&TE==p[i,'TE']&level==p[i,'level'])
  auc_temp1 <- auc_temp[1,]
  auc_temp1[,'TR'] <- paste(p[i,'TE'],'external_mean',sep = '_')
  # median
  auc_temp1$auc <- median(auc_temp$auc)
  auc_temp1 <- rbind(auc_temp1,auc_lodo_temp)
  if(i==1){
    auc_external_mean <- auc_temp1
  }else{
    auc_external_mean <- rbind(auc_external_mean,auc_temp1)
  }
} 
auc_external_mean_e <- auc_external_mean

self.e <- subset(auc_all,method %in%  c("self"))
p=unique(self.e[,c('disease','TE','level')])
self.lodo.e <- subset(auc_all,method %in%  c("lodo_self"))
for (i in 1:nrow(p)) {
  auc_temp  <- subset(self.e,disease==p[i,'disease']&TE==p[i,'TE']&level==p[i,'level'])
  auc_lodo_temp  <- subset(self.lodo.e,disease==p[i,'disease']&TE==paste(p[i,'TE'],'lodo',sep = '_')&level==p[i,'level'])
  auc_temp <- rbind(auc_temp,auc_lodo_temp)
  if(i==1){
    auc_external_mean <- auc_temp
  }else{
    auc_external_mean <- rbind(auc_external_mean,auc_temp)
  }
} 
auc_external_mean_s <- auc_external_mean

auc_external_median <- rbind(auc_external_mean_e,auc_external_mean_s)
auc_external_mean_e$method[auc_external_mean_e$method=='external'] <- 'Single cohort'
auc_external_mean_e$method[auc_external_mean_e$method=='lodo_external'] <- 'LODO'
auc_external_mean_e$method <- factor(auc_external_mean_e$method,levels=c('Single cohort','LODO'))
save(auc_external_mean_e,file = '01_data/AUC_ext_median.RData')

## AUC: single int- and ext- AUC ----
#disease category

a=get(load('./01_data/AUC_all_lasso.RData'))
a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
a$method[a$method == 'self'] <- 'internal'




# internal

auc_self <- subset(a,method == 'internal')
round(mean(auc_self$auc),3)
round(sd(auc_self$auc),3)
round(quantile(auc_self$auc, 1/4),3)
round(quantile(auc_self$auc, 3/4),3)



round(mean(subset(auc_self,group1=='Intestinal')$auc),3)
round(sd(subset(auc_self,group1=='Intestinal')$auc),3)
round(quantile(subset(auc_self,group1=='Intestinal')$auc, 1/4),3)
round(quantile(subset(auc_self,group1=='Intestinal')$auc, 3/4),3)


round(mean(subset(auc_self,group1=='Metabolic')$auc),3)
round(mean(subset(auc_self,group1=='Mental')$auc),3)
round(mean(subset(auc_self,group1=='Autoimmune')$auc),3)
round(mean(subset(auc_self,group1=='Liver')$auc),3)


round(mean(subset(auc_self,group1!='Intestinal')$auc),3)
round(sd(subset(auc_self,group1!='Intestinal')$auc),3)
round(quantile(subset(auc_self,group1!='Intestinal')$auc, 1/4),3)
round(quantile(subset(auc_self,group1!='Intestinal')$auc, 3/4),3)


# external
auc_external <- subset(a,method == 'external')

round(mean(auc_external$auc),3)
round(sd(auc_external$auc),3)
round(quantile(auc_external$auc, 1/4),3)
round(quantile(auc_external$auc, 3/4),3)



round(mean(subset(auc_external,group1=='Intestinal')$auc),3)
round(sd(subset(auc_external,group1=='Intestinal')$auc),3)
round(quantile(subset(auc_external,group1=='Intestinal')$auc, 1/4),3)
round(quantile(subset(auc_external,group1=='Intestinal')$auc, 3/4),3)

round(mean(subset(auc_external,group1=='Metabolic')$auc),3)
round(mean(subset(auc_external,group1=='Mental')$auc),3)
round(mean(subset(auc_external,group1=='Autoimmune')$auc),3)
round(mean(subset(auc_external,group1=='Liver')$auc),3)

round(mean(subset(auc_external,group1!='Intestinal')$auc),3)
round(sd(subset(auc_external,group1!='Intestinal')$auc),3)
round(quantile(subset(auc_external,group1!='Intestinal')$auc, 1/4),3)
round(quantile(subset(auc_external,group1!='Intestinal')$auc, 3/4),3)




#data_type
#int
auc_all <- get(load('./01_data/AUC_all_lasso.RData'))
auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","Amplicon_genus"))

auc0 <- auc_all
auc0 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
self.e <- subset(auc0,method=='self')

round(mean(subset(self.e,level=='mNGS_species')$auc),3)
round(mean(subset(self.e,level=='mNGS_genus')$auc),3)
round(mean(subset(self.e,level=='Amplicon_genus')$auc),3)
#ext
self.e <- subset(auc0,method=='external')

round(mean(subset(self.e,level=='mNGS_species')$auc),3)
round(mean(subset(self.e,level=='mNGS_genus')$auc),3)
round(mean(subset(self.e,level=='Amplicon_genus')$auc),3)


## AUC: intestinal disease AUC & reference-----------------------------------------
a=get(load('./01_data/AUC_all_lasso.RData'))
a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'

auc_external <- subset(a,method=='external')
{
  paste('CRC',round(mean(subset(auc_external,disease=='CRC')$auc),2))
  paste('CD',round(mean(subset(auc_external,disease=='CD')$auc),2))
  paste('IBD',round(mean(subset(auc_external,disease=='IBD')$auc),2))
  paste('CDI',round(mean(subset(auc_external,disease=='CDI')$auc),2))
  paste('UC',round(mean(subset(auc_external,disease=='UC')$auc),2))
  paste('Adenoma',round(mean(subset(auc_external,disease=='Adenoma')$auc),2))
  paste('IBS',round(mean(subset(auc_external,disease=='IBS')$auc),2))
}


paste('CD',round(mean(subset(auc_external,disease=='CD')$auc),3))
paste('CD',round(mean(subset(auc_external,disease=='CD'&level=='mNGS_genus')$auc),3))
paste('IBD',round(mean(subset(auc_external,disease=='IBD')$auc),3))
paste('CRC',round(mean(subset(auc_external,disease=='CRC')$auc),3))
paste('CDI',round(mean(subset(auc_external,disease=='CDI')$auc),3))
paste('UC',round(mean(subset(auc_external,disease=='UC')$auc),3))
paste('Adenoma',round(mean(subset(auc_external,disease=='Adenoma')$auc),3))
paste('IBS',round(mean(subset(auc_external,disease=='IBS')$auc),3))

paste('CD & IBD & CRC',round(mean(subset(auc_external,disease%in%c('CD','IBD','CRC'))$auc),2))

auc_external <- subset(a,method=='self')
paste('CD',round(mean(subset(auc_external,disease=='CD')$auc),3))
paste('IBD',round(mean(subset(auc_external,disease=='IBD')$auc),3))
paste('CRC',round(mean(subset(auc_external,disease=='CRC')$auc),3))
paste('CDI',round(mean(subset(auc_external,disease=='CDI')$auc),3))
paste('UC',round(mean(subset(auc_external,disease=='UC')$auc),3))
paste('Adenoma',round(mean(subset(auc_external,disease=='Adenoma')$auc),3))
paste('IBS',round(mean(subset(auc_external,disease=='IBS')$auc),3))


#Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer.
CRC_ex_1 <- c(0.71,0.79,0.72,0.67,0.60,0.50,
              0.83,0.77,0.86,0.66,0.66,0.68,
              0.76,0.66,0.81,0.71,0.59,0.55,
              0.78,0.80,0.82,0.62,0.73,0.56,
              0.63,0.78,0.65,0.72,0.75,0.55,
              0.77,0.79,0.75,0.83,0.80,0.59,
              0.54,0.62,0.54,0.63,0.59,0.55)
mean(CRC_ex_1)

CRC_ex_2 <- c(0.76,0.82,0.64,0.83,
              0.62,0.74,0.59,0.65,
              0.82,0.76,0.67,0.83,
              0.76,0.78,0.70,0.79,
              0.84,0.79,0.88,0.74)
mean(CRC_ex_2)




## AUC: LODO int- and ext- AUC improvement -----
#external
auc_all <- get(load('./01_data/lodo_AUC_all_lasso.RData'))
{
  self.e <- subset(auc_all,method %in%  c("external"))
  p=unique(self.e[,c('disease','TE','level')])
  self.lodo.e <- subset(auc_all,method %in%  c("lodo_external"))
  for (i in 1:nrow(p)) {
    auc_temp  <- subset(self.e,disease==p[i,'disease']&TE==p[i,'TE']&level==p[i,'level'])
    auc_lodo_temp  <- subset(self.lodo.e,disease==p[i,'disease']&TE==p[i,'TE']&level==p[i,'level'])
    auc_temp1 <- auc_temp[1,]
    auc_temp1[,'TR'] <- paste(p[i,'TE'],'external_mean',sep = '_')
    # median
    auc_temp1$auc <- median(auc_temp$auc)
    auc_temp1 <- rbind(auc_temp1,auc_lodo_temp)
    if(i==1){
      auc_external_mean <- auc_temp1
    }else{
      auc_external_mean <- rbind(auc_external_mean,auc_temp1)
    }
  } 
  auc_external_mean_e <- auc_external_mean
}

#internal
{
  self.e <- subset(auc_all,method %in%  c("self"))
  p=unique(self.e[,c('disease','TE','level')])
  self.lodo.e <- subset(auc_all,method %in%  c("lodo_self"))
  for (i in 1:nrow(p)) {
    auc_temp  <- subset(self.e,disease==p[i,'disease']&TE==p[i,'TE']&level==p[i,'level'])
    auc_lodo_temp  <- subset(self.lodo.e,disease==p[i,'disease']&TE==paste(p[i,'TE'],'lodo',sep = '_')&level==p[i,'level'])
    auc_temp <- rbind(auc_temp,auc_lodo_temp)
    if(i==1){
      auc_external_mean <- auc_temp
    }else{
      auc_external_mean <- rbind(auc_external_mean,auc_temp)
    }
  } 
  auc_external_mean_s <- auc_external_mean
}
auc_external_median <- rbind(auc_external_mean_e,auc_external_mean_s)
auc_external_mean_e$method[auc_external_mean_e$method=='external'] <- 'Single cohort'
auc_external_mean_e$method[auc_external_mean_e$method=='lodo_external'] <- 'LODO'
auc_external_mean_e$method <- factor(auc_external_mean_e$method,levels=c('Single cohort','LODO'))

#save(auc_external_median,file = 'auc_external_median.RData')

#statistic

for (i in c('T2D', 'AD', 'PD', 'ASD', 'NAFLD')){
  print(paste0(i,' : ','single:',round(median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc),2),', '
               ,'LODO:',round(median(subset(auc_external_mean_e,method=='LODO'&disease==i)$auc),2),', '
               ,round((median(subset(auc_external_mean_e,method=='LODO'&disease==i)$auc)-median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc))/median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc)*100,2),'%'))
  
}

for (i in c('CD', 'UC', 'IBS', 'CRC', 'Adenoma')){
  print(paste0(i,' : ','single:',round(median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc),2),', '
               ,'LODO:',round(median(subset(auc_external_mean_e,method=='LODO'&disease==i)$auc),2),', '
               ,round((median(subset(auc_external_mean_e,method=='LODO'&disease==i)$auc)-median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc))/median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc)*100,2),'%'))
  
}

print(paste0('non-Intestinal: single study ext-AUC median: ',round(median(subset(auc_external_mean_e,method=='Single cohort'&group2=='non-Intestinal')$auc),2)))
print(paste0('non-Intestinal: LODO ext-AUC median: ',round(median(subset(auc_external_mean_e,method=='LODO'&group2=='non-Intestinal')$auc),2)))

print(paste0('Intestinal: single study ext-AUC median: ',round(median(subset(auc_external_mean_e,method=='Single cohort'&group2=='Intestinal')$auc),2)))
print(paste0('Intestinal: LODO ext-AUC median: ',round(median(subset(auc_external_mean_e,method=='LODO'&group2=='Intestinal')$auc),2)))


# print(paste0('non-Intestinal: single study ext-AUC median: ',round(median(subset(auc_external_mean_s,method=='self'&group2=='non-Intestinal')$auc),2)))
# print(paste0('non-Intestinal: LODO ext-AUC median: ',round(median(subset(auc_external_mean_s,method=='lodo_self'&group2=='non-Intestinal')$auc),2)))
# 
# print(paste0('Intestinal: single study ext-AUC median: ',round(median(subset(auc_external_mean_s,method=='self'&group2=='Intestinal')$auc),2)))
# print(paste0('Intestinal: LODO ext-AUC median: ',round(median(subset(auc_external_mean_s,method=='lodo_self'&group2=='Intestinal')$auc),2)))

# ASD & PD 16S  LODO--

ASD <- subset(auc_external_mean_e,method=='LODO'&disease=="ASD"&level=='Amplicon_genus')
median(ASD$auc)
PD <- subset(auc_external_mean_e,method=='LODO'&disease=="PD"&level=='Amplicon_genus')
median(PD$auc)
CRC <- subset(auc_external_mean_e,method=='LODO'&disease=="CRC"&level=='Metagenomics_species')
median(CRC$auc)
CD <- subset(auc_external_mean_e,method=='LODO'&disease=="CD"&level=='Metagenomics_species')
median(CD$auc)





#
for (i in c('T2D', 'AD', 'PD', 'ASD', 'NAFLD')){
  print(paste0(i,' : ','single:',round(median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc),2),', '
               ,'LODO:',round(median(subset(auc_external_mean_e,method=='LODO'&disease==i)$auc),2),', '
               ,round((median(subset(auc_external_mean_e,method=='LODO'&disease==i)$auc)-median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc))/median(subset(auc_external_mean_e,method=='Single cohort'&disease==i)$auc)*100,2),'%'))
  
}
print(paste0('non-Intestinal: single study ext-AUC median: ',round(median(subset(auc_external_mean_e,method=='Single cohort'&group2=='non-Intestinal')$auc),2)))
print(paste0('non-Intestinal: LODO ext-AUC median: ',round(median(subset(auc_external_mean_e,method=='LODO'&group2=='non-Intestinal')$auc),2)))

print(paste0('Intestinal: single study ext-AUC median: ',round(median(subset(auc_external_mean_e,method=='Single cohort'&group2=='Intestinal')$auc),2)))
print(paste0('Intestinal: LODO ext-AUC median: ',round(median(subset(auc_external_mean_e,method=='LODO'&group2=='Intestinal')$auc),2)))



print(paste0('non-Intestinal: single study ext-AUC median: ',round(median(subset(auc_external_mean_s,method=='self'&group2=='non-Intestinal')$auc),2)))
print(paste0('non-Intestinal: LODO ext-AUC median: ',round(median(subset(auc_external_mean_s,method=='lodo_self'&group2=='non-Intestinal')$auc),2)))

print(paste0('Intestinal: single study ext-AUC median: ',round(median(subset(auc_external_mean_s,method=='self'&group2=='Intestinal')$auc),2)))
print(paste0('Intestinal: LODO ext-AUC median: ',round(median(subset(auc_external_mean_s,method=='lodo_self'&group2=='Intestinal')$auc),2)))





