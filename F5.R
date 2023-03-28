### Figure 5
{
  require(psych)
  require(ggcorrplot)
  library(showtext)
}

## F45 -----
{
  auc_marker <- get(load('01_data/AUC_MSI_euclidean_median.RData'))
  auc_marker$level[auc_marker$level=='Metagenomics_genus']	<- 'mNGS_genus'
  auc_marker$level[auc_marker$level=='Metagenomics_species']	<- 'mNGS_species'
  auc_marker$level[auc_marker$level=='Amplicon_genus']	<- '16S_genus'
  auc_marker$level <- factor(auc_marker$level,levels = c('mNGS_species','mNGS_genus','16S_genus'))
}
save(auc_marker,file = '01_data/plot_data/F5A.RData')

{
  cor <- corr.test(auc_marker[,1:2],method = 'spearman')
  cor$p[1,2]
  cor$r[1,2]
}
p1=ggscatter(data=auc_marker,x='dis',y='auc',shape='level',color='disease',palette = "igv",alpha=0.8,size=2)+
  geom_smooth(mapping=aes(x=dis,y=auc),data=auc_marker,formula =y~x, method=lm, inherit.aes = F,size=0.5)+
  geom_text_repel(data=auc_marker,aes(x=dis,y=auc,label = disease),inherit.aes = F,size=3)+
  xlab("MSI (median)") + ylab("External AUC (median)")+
  labs(color = "disease")+
  theme(legend.position = 'none')+
  annotate('text',x=0.7,y=0.84,label=paste0('r=',round(cor$r[1,2],2),'  ','p','=',format(cor$p[1,2],scientific=T,digits=3)),size=3)+
  border()

p2 <- ggdensity(auc_marker, 'dis', fill="disease_type2",palette = "npg")+
  clean_theme()+
  theme(legend.position = 'none')
p3 <- ggdensity(auc_marker, 'auc', fill="disease_type2", palette = "npg")+
  clean_theme()+
  theme(legend.position = 'none')+
  ggpubr::rotate()

p=ggarrange(p2, NULL, p1, p3, ncol = 2, nrow = 2, align = "hv", widths = c(3, 1), 
            heights = c(1, 3))

print(p)
pdf("./02_figure/Fig5A.pdf", height = 3.8, width = 5)
p
dev.off()


print(paste('ASD:',subset(auc_marker,disease=='ASD'& level=='16S_genus')$dis))
print(paste('CRC:',subset(auc_marker,disease=='CRC'& level=='mNGS_genus')$dis))

# suplement p value
wilcox.test(subset(auc_marker,disease_type2=='non-Intestinal')$auc,subset(auc_marker,disease_type2=='Intestinal')$auc)
wilcox.test(subset(auc_marker,disease_type2=='non-Intestinal')$dis,subset(auc_marker,disease_type2=='Intestinal')$dis)

## F5B ----

self.e <- get(load('./01_data/MSI_all.RData'))
self.e$data_type[self.e$data_type=='Metagenomics'] <- 'mNGS'
self.e$data_type[self.e$data_type=='Amplicon'] <- '16S'
self.e$level <- paste0(self.e$data_type,'_',self.e$taxon)
self.e$level <- factor(self.e$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))

save(self.e,file = '01_data/plot_data/F5B.RData')

stat.test <- compare_means(
  dis_value~group1,data = self.e, group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = c(seq(from=4, to=7.5,length.out=10),rep(seq(from=4, to=6,length.out=6), times=2)))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

# adj
p <- ggboxplot(self.e, x = "group1", y = "dis_value",fill = "group1", palette = "jco", 
               facet.by = "level",width = 0.5,outlier.shape = NA)+
  stat_compare_means(aes(label = paste0('p = ',..p.format..)),label.y = 6.9,)+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank(),legend.position="none")+
  xlab("Disease category") + ylab("MSI")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black')) 
p=p + stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p)



pdf("./02_figure/Fig5B.pdf", height = 5, width = 8)
p
dev.off()

## F5C ----
# data

self.e <- get(load('./01_data/MSI_all.RData'))
self.e$data_type[self.e$data_type=='Metagenomics'] <- 'mNGS'
self.e$data_type[self.e$data_type=='Amplicon'] <- '16S'
self.e$level <- paste0(self.e$data_type,'_',self.e$taxon)
self.e$level <- factor(self.e$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))

self.e <- subset(self.e,disease %in% names(which(table(unique(self.e[,c('disease','data_type')])$disease)>1)))
self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune"))
colnames(self.e)[colnames(self.e)=='group1']='subtype'

save(self.e,file = '01_data/plot_data/F5C.RData')

stat.test <- compare_means(
  dis_value~level,data = self.e, group.by = "subtype",
  method = "wilcox.test"
)%>%mutate(y.position = rep(seq(from=4.5, to=5.5,length.out=3), times=4))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

# adj
p <- ggboxplot(self.e, x = "level", y = "dis_value",fill = "level", 
               palette = c('#774ec7','#bd93cc','#a2c4b1'),
               facet.by = "subtype",nrow = 1,width = 0.5,outlier.shape = NA)+
  stat_compare_means(aes(label = paste0('p = ',..p.format..)),label.y = 6.9,)+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank(),legend.position="top")+ 
  xlab("Data type") + ylab("MSI")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position = 'none')
p=p + stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p)

pdf("./02_figure/Fig5C.pdf", height = 5, width = 8)
p
dev.off()


## F4D -----
{
  auc_marker_lodo <- get(load('./01_data/AUC_MSI_lodo_lasso.RData'))
  auc_marker_lodo$level[auc_marker_lodo$level %in% 'Metagenomics_genus'] <- 'mNGS_genus'
  auc_marker_lodo$level[auc_marker_lodo$level %in% 'Metagenomics_species'] <- 'mNGS_species'
  auc_marker_lodo$level[auc_marker_lodo$level %in% 'Amplicon_genus'] <- '16S_genus'
  auc_marker_lodo$method <- 'LODO'
  auc_marker <- get(load('./01_data/AUC_MSI_euclidean_median.RData'))
  auc_marker$method <- 'external'
  auc_marker$level[auc_marker$level=='Metagenomics_genus']	<- 'mNGS_genus'
  auc_marker$level[auc_marker$level=='Metagenomics_species']	<- 'mNGS_species'
  auc_marker$level[auc_marker$level=='Amplicon_genus']	<- '16S_genus'
  for (i in 1:nrow(auc_marker_lodo)){
    auc_marker_lodo <- rbind(auc_marker_lodo,subset(auc_marker,level==auc_marker_lodo$level[i] & disease==auc_marker_lodo$disease[i]))
  }
  auc_marker_lodo$method[auc_marker_lodo$method=='external'] <- 'single-study'
  auc_marker_lodo$method <- factor(auc_marker_lodo$method,levels = c('single-study','LODO'))
}

save(auc_marker_lodo,file = '01_data/plot_data/F5D.RData')

p <- ggpaired(auc_marker_lodo, x = 'method', y = 'dis',color = 'method', palette = "npg", 
              line.color = "gray", line.size = 0.4,
              short.panel.labs = FALSE,width = 0.4)+
  stat_compare_means(paired = T,label.y=4.6)+
  theme(legend.position = 'none')+
  labs(y='MSI (median)',x='')+
  theme(text = element_text(family = '',size = 13,face = 'plain'))

print(p)

pdf("./02_figure/Fig5D.pdf", height = 3.2, width = 1.9)
p
dev.off()


## F4E ------
load('01_data/p_r_ASD_marker_new.RData')
p_r_test <- subset(p_r,method=='test')
save(p_r_test,file = '01_data/plot_data/F5E_left.RData')
for (i in unique(p_r_test$num)){
  dis_median <- median(subset(p_r_test,num==i)$dis)
  if (i==1){
    p_r_test_median <- data.frame(dis=dis_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(dis=dis_median,method='test',num=i))
  } 
}


p1<-ggboxplot( p_r_test,x='num',y='dis',color=pal_npg("nrc")(10)[1],add = "jitter",width = 0.3)+
  xlab('No. of training cohorts')+
  ylab('MSI (median)')+
  labs(title = 'ASD')+
  geom_line(data = p_r_test_median,aes(x=num,y=dis),color=pal_npg("nrc")(10)[2])+
  geom_point(data = p_r_test_median,aes(x=num,y=dis),color=pal_npg("nrc")(10)[2])+
  # geom_hline(yintercept =0.5,color=pal_d3('category10')(10)[4])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
ylim(0,2.5)
p1


load('data/p_r_PD_marker_new.RData')
save(p_r_test,file = '01_data/plot_data/F5E_right.RData')

p_r_test <- subset(p_r,method=='test')
for (i in unique(p_r_test$num)){
  dis_median <- median(subset(p_r_test,num==i)$dis)
  if (i==1){
    p_r_test_median <- data.frame(dis=dis_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(dis=dis_median,method='test',num=i))
  } 
}
p2<-ggboxplot( p_r_test,x='num',y='dis',color=pal_npg("nrc")(10)[1],add = "jitter",width = 0.3)+
  xlab('No. of training cohorts')+
  ylab('MSI (median)')+
  labs(title = 'PD')+
  geom_line(data = p_r_test_median,aes(x=num,y=dis),color=pal_npg("nrc")(10)[2])+
  geom_point(data = p_r_test_median,aes(x=num,y=dis),color=pal_npg("nrc")(10)[2])+
  # geom_hline(yintercept =0.5,color=pal_d3('category10')(10)[4])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
# ylim(0.7,3)

p=ggarrange(p1,p2)
p
pdf("./02_figure/Fig5E.pdf", height = 3.8, width = 8)
p
dev.off()



















