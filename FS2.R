### Figure S2

## FS2 A ----

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


save(p_r_top_median,file='01_data/plot_data/FS2A.RData')

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

# 
# p1<-ggboxplot(data = p_r_top_median_external,x='num',y='AUC')+
#   xlab('No. of top feature')+
#   ylab('External AUC')+
#   # geom_line()+
#   geom_point(size=3)+
#   theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
#   theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position="right")+
#   geom_vline(xintercept =45,color='red')+
#   geom_vline(xintercept =20,color='blue')
#   # ylim(0.42,0.95)
# p1
# for (i in unique(p_r_top_median_external$num)){
#   message(median(subset(p_r_top_median_external,num==i)$AUC))
#   }

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
p=p2+p1
p

pdf("./02_figure/FigS2A.pdf", height = 3.8, width = 10)
p
dev.off()






## FS2 B ----
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

save(a,file='01_data/plot_data/FS2B.RData')

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


pdf("./02_figure/FigS2B.pdf", height = 3.8, width = 5.8)
p
dev.off()


