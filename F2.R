### Figure 2

## F2A ----
{
a=get(load('./01_data/AUC_all_lasso.RData'))
a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
a$method[a$method == 'self'] <- 'internal'
auc_self <- subset(a,method == 'internal')
auc_external <- subset(a,method=='external')
}
save(auc_external,file = '01_data/plot_data/F2A.RData')

stat.test <- compare_means(
  auc~group1,data = auc_external,
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.15, to=1.9,length.out=10))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')
p1 <- ggviolin(auc_external, x = "group1", y = "auc", fill = "group1",alpha = 0.3,
               palette = "jco",add = "boxplot",width = 0.6)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  # stat_compare_means()+
  ylim(0.05,1.9)+
  theme(legend.position="none")+   
  xlab("") + ylab("External AUC")+
  ggtitle("Disease category")+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p1)





## F2B ----
# only included diseases with three data types
{
auc_all <- get(load('01_data/AUC_all_lasso.RData'))
auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
auc0 <- auc_all
auc0 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
self.e <- subset(auc0,method=='external')
}
save(self.e,file = '01_data/plot_data/F2B.RData')
stat.test <- compare_means(
  auc~level,data = self.e, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.35, to=1.7,length.out=3))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p2 <- ggviolin(self.e, x = "level", y = "auc", fill = "level",alpha = 0.3,
           palette = c('#774ec7','#bd93cc','#a2c4b1'),add = "boxplot",width = 0.5)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  ylim(0.05,1.9)+
  # stat_compare_means()+
  theme(legend.position="none")+    
  ylab("External AUC")+xlab('')+
  ggtitle('Data type')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))

p2 <- p2+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p2)

#stat
{
t.test(subset(self.e,level=='16S_genus')$auc)
sd(subset(self.e,level=='16S_genus')$auc)

t.test(subset(self.e,level=='mNGS_genus')$auc)
sd(subset(self.e,level=='mNGS_genus')$auc)

t.test(subset(self.e,level=='mNGS_species')$auc)
sd(subset(self.e,level=='mNGS_species')$auc)
}
p=ggarrange(p1,p2,
            ncol = 2, nrow = 1,
            widths = c(7,4.5)
)
p
pdf("./02_figure/Fig2AB.pdf", height = 4.5, width = 6)
p
dev.off()



## F2C ----
# data
{
auc_all <- get(load('01_data/AUC_all_lasso.RData'))
auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
self.e <- subset(auc_all,method=='external')
self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune","Liver"))
}

save(self.e,file = '01_data/plot_data/F2C.RData')

stat.test <- compare_means(
  auc~group1,data = self.e, group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = c(seq(from=1, to=1.5,length.out=10),rep(seq(from=1, to=1.3,length.out=6), times=2)))
x <- stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p <- ggboxplot(self.e, x = "group1", y = "auc",color = "group1", palette = "jco", add = "jitter",
               width = 0.5,add.params = list(size=1))+
  facet_wrap(vars(level),scales = "free")+
  stat_compare_means(aes(label = paste0('p = ',..p.format..)),label.y = 1.5,hjust=-0.01)+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank(),legend.position="none")+ 
  xlab("Disease category") + ylab("External AUC")+
  labs(color = "disease type")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p=p + stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p)

pdf("./02_figure/Fig2C.pdf", height = 4.5, width = 6.7)
p
dev.off()


## F2D ----
#data
{
auc_all <- get(load('01_data/AUC_all_lasso.RData'))
auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
auc0 <- subset(auc_all,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
self.e <- subset(auc0,method=='external')
self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune"))
colnames(self.e)[8]='subtype'
}
save(self.e,file = '01_data/plot_data/F2D.RData')

stat.test <- compare_means(
  auc~level,data = self.e, group.by = "subtype",
  method = "wilcox.test"
)%>%mutate(y.position = rep(seq(from=1, to=1.15,length.out=3), times=4))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p <- ggboxplot(self.e, x = "level", y = "auc",color = "level", add = "jitter",
               palette = c('#774ec7','#bd93cc','#a2c4b1'),
               width = 0.5,add.params = list(size=1))+
  facet_wrap(vars(subtype),nrow = 1)+
  stat_compare_means(aes(label = paste0('p = ',..p.format..)),label.y = 1.25,hjust=-0.01)+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank(),legend.position="none")+
  xlab("Data type") + ylab("External AUC")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p=p + stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p)

pdf("./02_figure/Fig2D.pdf", height = 4.5, width = 6.7)
p
dev.off()

