### Figure S1

## FS1 A -----
a1=get(load('./01_data/AUC_all_enet.RData'))
a1$model_method <- 'enet'
a2=get(load('./01_data/AUC_all_lasso.RData'))
a2$model_method <- 'lasso'
a3=get(load('./01_data/AUC_all_ridge.RData'))
a3$model_method <- 'ridge'
a4=get(load('./01_data/AUC_all_rf.RData'))
a4$model_method <- 'rf'

library(dplyr)
auc_all_method <- bind_rows(list(a1, a2,a3,a4))
auc_all_method$method[auc_all_method$method=='self'] <- 'Internal'
auc_all_method$method[auc_all_method$method=='external'] <- 'External'
auc_self_method <- subset(auc_all_method,method=='Internal')
auc_external_method <- subset(auc_all_method,method=='External')
{
  l <- unique(auc_all_method$method)
  pailie <- t(combn(1:length(l),2))
  my_comparisons_f0 <- list()
  for (i in 1:nrow(pailie)) {
    my_comparisons_f0[[i]] <- c(l[pailie[i,]])
  }
}

save(auc_all_method,file='01_data/plot_data/FS1A.RData')

library(ggpubr)
p <- ggboxplot(auc_all_method, x = "model_method", y = "auc", fill="model_method",
               palette = "npj",width = 0.3)+ 
  facet_grid(~method)+
  stat_compare_means(label.y = 1.08)+
  theme(legend.position="none")+    
  xlab("Modeling method") + ylab("AUC")+
  labs(color = "method")+
  ylim(0.1,1.2)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),
        axis.text.x = element_blank(),axis.ticks.x=element_blank()) 
print(p)

pdf("./02_figure/FigS1A.pdf", height = 2.6, width = 4.5)
p
dev.off()


## FS1 B ----
load('./01_data/cohort_info.RData')
cross_top <- subset(cross,num>3)
cross_top <- subset(cross_top,(data_type == 'Metagenomics' & taxon=='species')|(data_type=='Amplicon' & taxon=='genus'))
cross_top <- subset(cross_top,disease!='T2D')
p_r_top_median <- data.frame()
for (d_i in 1:nrow(cross_top)){
  auc_top <- get(load(paste0('04_model_raw_RData/single_LODO/02lasso/nest_top/',cross_top[d_i,'models'])))
  p_top_r <- data.frame()
  for (i in names(auc_top)){
    top_r=auc_top[[i]]
    top_r=top_r[,3:ncol(top_r)]
    
    project0 <- colnames(top_r)
    for (p_i in 1:length(project0)){
      p=project0[p_i]
      p_top_r <- rbind(p_top_r,data.frame(AUC=top_r[,p][-which(p ==project0)],method='test',num=i,project.test=p,project.train=project0[-p_i]))
    }
    
    top_r<-top_r[,rev(sequence(ncol(top_r)))]
    p_top_r <- rbind(p_top_r,data.frame(AUC=rev(top_r[row(top_r) == NCOL(top_r) - col(top_r) + 1]),method='train',num=i,project.test=0,project.train=rownames(top_r)))
  }
  
  p_r_median <- data.frame()
  for (i in unique(p_top_r$num)){
    AUC_median <- median(subset(p_top_r,num==i&method=='train')$AUC)
    p_r_median <- rbind(p_r_median,data.frame(AUC=AUC_median,method='train',num=i))
    AUC_median <- median(subset(p_top_r,num==i&method=='test')$AUC)
    p_r_median <- rbind(p_r_median,data.frame(AUC=AUC_median,method='test',num=i))
  }
  p_r_median$disease <- cross_top[d_i,'disease']
  p_r_median$level <- paste(cross_top[d_i,'data_type'],cross_top[d_i,'taxon'],sep = '_')
  p_r_top_median <- rbind(p_r_top_median,p_r_median)
}


p_r_top_median$num[p_r_top_median$num=='all'] <- 0.5
p_r_top_median$num <- as.numeric(p_r_top_median$num)
for (i in unique(p_r_top_median$disease)){
  p_r_top_median$num[(p_r_top_median$num=='0.5')&(p_r_top_median$disease==i)] <- max(subset(p_r_top_median,disease==i)$num)+5
}


save(p_r_top_median,file='01_data/plot_data/FS1B.RData')

p_r_top_median_external <- subset(p_r_top_median,method=='test')
p_r_top_median_internal <- subset(p_r_top_median,method=='train')

p_r_top_median_external <- subset(p_r_top_median_external,num<=45)
p1<-ggplot(data = p_r_top_median_external,aes(x=num,y=AUC,color=disease,shape=level))+
  xlab('No. of top feature')+
  ylab('External AUC')+
  geom_line()+
  geom_point(size=3)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position="right")+
  geom_vline(xintercept =45,color='red')+
  # geom_vline(xintercept =20,color='blue')+
  ylim(0.42,0.95)
p1



p_r_top_median_internal <- subset(p_r_top_median_internal,num<=45)
p2<-ggplot(data = p_r_top_median_internal,aes(x=num,y=AUC,color=disease,shape=level))+
  xlab('No. of top feature')+
  ylab('internal AUC')+
  geom_line()+
  geom_point(size=3)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position="none")+
  geom_vline(xintercept =45,color='red')+
  # geom_vline(xintercept =20,color='blue')+
  ylim(0.42,0.95)
p2
library(patchwork)
p=p2+p1
p

pdf("./02_figure/FigS1B.pdf", height = 3.8, width = 10)
p
dev.off()






## FS1 C ----
p_r <- data.frame()
for (c_i in 1:nrow(cross_top)){
  file1 <- cross_top$models[c_i]
  message(paste0('04_model_raw_RData/single_LODO/02lasso/nest_top/',file1))
  load(paste0('04_model_raw_RData/single_LODO/02lasso/nest_top/',file1))
  project0 <- row.names(model.adj[["40"]])
  for (i in 1:length(project0)){
    p=project0[i]
    p_r <- rbind(p_r,data.frame(AUC=model.adj[["40"]][,p][-which(p ==project0)],method='test',num=1,project.test=p,project.train=project0[-i],top='40'))
  }
  r=model.adj[["40"]]
  r=r[,3:ncol(r)]
  r<-r[,rev(sequence(ncol(r)))]
  p_r <- rbind(p_r,data.frame(AUC=rev(r[row(r) == NCOL(r) - col(r) + 1]),method='train',num=1,project.test=0,project.train=rownames(r),top='40'))
}
for (c_i in 1:nrow(cross_top)){
  file1 <- cross_top$models[c_i]
  message(paste0('04_model_raw_RData/single_LODO/02lasso/nest_top/',file1))
  load(paste0('04_model_raw_RData/single_LODO/02lasso/nest_top/',file1))
  project0 <- row.names(model.adj[["all"]])
  for (i in 1:length(project0)){
    p=project0[i]
    p_r <- rbind(p_r,data.frame(AUC=model.adj[["all"]][,p][-which(p ==project0)],method='test',num=1,project.test=p,project.train=project0[-i],top='all'))
  }
  r=model.adj[["all"]]
  r=r[,3:ncol(r)]
  r<-r[,rev(sequence(ncol(r)))]
  p_r <- rbind(p_r,data.frame(AUC=rev(r[row(r) == NCOL(r) - col(r) + 1]),method='train',num=1,project.test=0,project.train=rownames(r),top='all'))
}
a <- p_r
a$method[a$method%in%'train'] <- 'internal'
a$method[a$method%in%'test'] <- 'external'

save(a,file='01_data/plot_data/FS1C.RData')

p1 <- ggboxplot(subset(a,method=='internal'), x = "top", y = "AUC", fill = "top",
                palette = "d3",width = 0.3)+ 
  stat_compare_means(label.y=1.05)+
  xlab("No. of top feature") + ylab("Internal AUC")+
  labs(color = "method")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")+
  ylim(0.2,1.1)
p2 <- ggboxplot(subset(a,method=='external'), x = "top", y = "AUC", fill = "top",
                palette = "d3",width = 0.3)+ 
  stat_compare_means(label.y=1.05)+
  xlab("No. of top feature") + ylab("External AUC")+
  labs(color = "method")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")+
  ylim(0.2,1.1)
p=p1+p2
print(p)


pdf("./02_figure/FigS1C.pdf", height = 3.8, width = 5.8)
p
dev.off()



## FS1 D ----
auroc_all_filtered <- get(load('./01_data/AUC_all_lasso_filtered.RData'))
colnames(auroc_all_filtered)[colnames(auroc_all_filtered)=='auprc'] <- 'auc' 
auroc_all_filtered$Group <- 'Rel.'
auroc_all_normalized <- get(load('./01_data/AUC_all_lasso.RData'))
auroc_all_normalized$Group <- 'Log.'

auroc_all <- rbind(auroc_all_filtered,auroc_all_normalized)
auroc_all$method[auroc_all$method=='self'] <- 'Internal'
auroc_all$method[auroc_all$method=='external'] <- 'External'
auroc_all$method <- factor(auroc_all$method,levels = c('Internal','External'))

save(auroc_all,file='01_data/plot_data/FS1D.RData')

p <- ggboxplot(auroc_all, x = 'Group', y = 'auc', fill = 'Group',palette = "aaas", 
               line.color = "gray", line.size = 0.4,
               short.panel.labs = FALSE,width = 0.4)+
  facet_grid(~method)+
  stat_compare_means(paired = T,label.y=1.1)+
  # theme_bw()+
  theme(legend.position = 'top')+
  labs(y='AUC',x='')+
  theme(text = element_text(family = '',size = 12,face = 'plain'))
# axis.text.x=element_text(angle=45, hjust=1,family = '',size = 13,face = 'plain'))

print(p)

pdf("./02_figure/FigS1D.pdf", height = 3.4, width = 4)
p
dev.off()


## FS1 E ----

kk=read.csv('01_data/auc_comparison.csv',header = T)

colnames(kk)
kk=kk[,c("projects","disease",'AUC','final_AUC',"AUC.final_AUC",'label')]
colnames(kk)=c("projects","disease",'paper_AUC','our_AUC','diff','label')

library(ggplot2)
library(ggrepel)
min(kk$paper_AUC)
min(kk$our_AUC)

p=kk%>%
  mutate(label=ifelse(label%in%c('type3','type4'),'type2','type1'))%>%#只要两类
  ggplot()+
  geom_point(aes(x=our_AUC,y=paper_AUC,color=label,size=abs(diff)))+theme_classic()+
  scale_y_continuous(limits = c(0.3,1),breaks = seq(0.3,1,0.1))+
  scale_x_continuous(limits = c(0.3,1),breaks = seq(0.3,1,0.1))+
  # scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+ 
  # scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+ 
  scale_color_manual(values = c('type1'='#98d09d','type2'='#e77381'),name='')+
  # labs(title='AUC comparison')+
  geom_abline(slope = 1,intercept = 0,colour='#8f888b',size=1)+
  geom_abline(slope = 1,intercept = 0.15,linetype="dashed",colour='#8f888b',size=1)+
  geom_abline(slope = 1,intercept = -0.15,linetype="dashed",colour='#8f888b',size=1)+
  geom_text_repel(aes(x=our_AUC,y=paper_AUC,label=disease),colour="black", size=2,max.overlaps = 30)+
  theme(plot.title = element_text(hjust = 0.5),legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", color = "black"),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+
  coord_fixed()+
  xlab('Our AUC')+
  ylab("Literature's AUC")
p


pdf("./02_figure/FigS1E.pdf", height = 3.4, width = 4)
p
dev.off()



