### Figure S1

## FS1 A-----
a1=get(load('./01_data/AUC_all_enet.RData'))
a1$model_method <- 'enet'
a2=get(load('./01_data/AUC_all_lasso.RData'))
a2$model_method <- 'lasso'
a3=get(load('./01_data/AUC_all_ridge.RData'))
a3$model_method <- 'ridge'
a4=get(load('./01_data/AUC_all_rf.RData'))
a4$model_method <- 'rf'


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



## FS1 BC ----
file <- '04_model_raw_RData/single_LODO/02lasso/'
model.data <- 'D000067877_Amplicon.genus_6.RData'
model.adj1 <- get(load(paste0(file,model.data)))

save(model.adj1,file='01_data/plot_data/FS1B.RData')
model.adj1[["result"]]
p1 <- my_auc.heatmap(datas=model.adj1$result,models='lasso_adj_nc')


model.adj2 <- get(load(paste0(file,'D015179_Metagenomics.species_7.RData')))

save(model.adj2,file='01_data/plot_data/FS1C.RData')
p2 <- my_auc.heatmap(datas=model.adj2$result,models='lasso_adj_nc')

## FS1 D ----
a=get(load('./01_data/AUC_all_lasso.RData'))
a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
a$method[a$method == 'self'] <- 'internal'
a_d <- subset(a,group1=='Intestinal')
a_d$disease <- factor(a_d$disease,levels = c('CD','IBD','CRC','CDI','UC','Adenoma','IBS'))
save(a_d,file='01_data/plot_data/FS1D.RData')
p <- ggboxplot(a_d, x = "method", y = "auc", color = 'black',fill ="method" ,
               palette = c('#1fb8b4','#ff7f0e'),width = 0.4)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  geom_hline(yintercept =0.9,color='#e60020')+
  facet_wrap(~disease,nrow = 1)+
  stat_compare_means(label.y=1.05,label = "p.signif")+
  ylim(0.3,1.08)+
  theme(legend.position="right")+    
  xlab("") + ylab("AUC")+
  labs(fill = "Group")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position = 'none',axis.text.x = element_blank(),axis.ticks.x =element_blank())

print(p)

pdf("./02_figure/FigS1D.pdf", height = 3.4, width = 6)
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



## FS1 paired -----
a=get(load('./data/auc_all_lasso_COVID-19.RData'))
a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'

a$method[a$method == 'self'] <- 'internal'
a$group1[a$group1=='Liver disease'] <- 'Others'
auc_self <- subset(a,method == 'internal')
auc_external <- subset(a,method=='external')


all_p <- auc_self
all_p$TE <- NULL
all_p$auc <- 0
for (i in 1:nrow(all_p)){
  auc_t <- subset(auc_external,(TR==all_p[i,'TR'])&(disease==all_p[i,'disease'])&(level==all_p[i,'level']))
  all_p[i,'auc'] <- median(auc_t$auc)
}

# by disease
all_d <- all_p 
all_d$TR <- NULL
all_d$auc <- 0
all_d <- unique(all_d)
all_d_internal <- all_d
for (i in 1:nrow(all_d)){
  auc_t <- subset(all_p,(disease==all_d[i,'disease'])&(level==all_d[i,'level']))
  all_d[i,'auc'] <- median(auc_t$auc)
  auc_t <- subset(auc_self,(disease==all_d[i,'disease'])&(level==all_d[i,'level']))
  all_d_internal[i,'auc'] <- median(auc_t$auc)
}
all_d$model <- 'External'
all_d_internal$model <- 'Internal'
all_auc <- rbind(all_d,all_d_internal)

# plot
p <- ggpaired(all_auc, x = 'model', y = 'auc',color = 'model', palette = "lancet", 
              line.color = "gray", line.size = 0.4,
              short.panel.labs = FALSE,width = 0.4,facet.by = 'group2')+
  stat_compare_means(paired = T,label.y=1)+
  theme(legend.position = 'top')+
  labs(y='AUC',x='')+
  theme(text = element_text(family = '',size = 13,face = 'plain'))
all_d[which((all_d$auc-all_d_internal$auc)>0),]
all_d[which((all_d$auc-all_d_internal$auc)==0),]
print(p)
pdf("./02_figure/FigS1E.pdf", height = 5, width = 5)
p
dev.off()

