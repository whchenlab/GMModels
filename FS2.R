### Figure S2

## FS2 A ----
load('./01_data/all_measure_lasso.RData')
{
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
}
{
  a=all_measure
  a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
  a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
  a$method[a$method == 'self'] <- 'internal'
  auc_self <- subset(a,method == 'internal')
}
save(auc_self,file = '01_data/plot_data/FS2A.RData')

stat.test <- compare_means(
  auprc~group1,data = auc_self, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.05, to=1.65,length.out=10))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p1 <- ggboxplot(auc_self, x = "group1", y = "auprc", fill = "group1",
                palette = "jco",width = 0.2)+ 
  # geom_hline(yintercept =0.5,color='#dbdcdc')+
  # geom_hline(yintercept =0.6,color='#ffd09a')+
  # geom_hline(yintercept =0.7,color='#ffcbd8')+
  # geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  theme_bw() +
  # stat_compare_means()+
  ylim(0.15,1.68)+
  theme(legend.position="none")+    
  ylab("Internal AUPRC")+xlab('')+
  ggtitle('Disease category')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p1)

#data_type
{
  auc_all <- all_measure
  auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
  auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
  auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
  auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
  auc0 <- auc_all
  auc0 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
  self.e <- subset(auc0,method=='self')
  self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune","Liver"))
}

# save(self.e,file='01_data/plot_data/F1E.RData')

stat.test <- compare_means(
  auprc~level,data = self.e, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.2, to=1.65,length.out=3))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p2 <- ggboxplot(self.e, x = "level", y = "auprc", fill = "level",
                width = 0.2,palette = c('#774ec7','#bd93cc','#a2c4b1'))+
  # geom_hline(yintercept =0.5,color='#dbdcdc')+
  # geom_hline(yintercept =0.6,color='#ffd09a')+
  # geom_hline(yintercept =0.7,color='#ffcbd8')+
  # geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  # ylim(0.05,1.68)+
  theme_bw() +
  ylim(0.15,1.68)+
  # stat_compare_means()+
  theme(legend.position="none")+    
  ylab("Internal AUPRC")+xlab('')+
  ggtitle('Data type')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p2 <- p2+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p2)

pA=ggarrange(p1,p2,
             ncol = 2, nrow = 1,
             widths = c(4,3)
)
pA
pdf("./02_figure/FigS2A.pdf", height = 4.5, width = 5)
pA
dev.off()


## FS2 B ---- 
load('./01_data/all_measure_lasso.RData')
{
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
}
{
  a=all_measure
  a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
  a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
  a$method[a$method == 'self'] <- 'internal'
  auc_self <- subset(a,method == 'external')
}
save(auc_self,file = '01_data/plot_data/FS2B.RData')

stat.test <- compare_means(
  auprc~group1,data = auc_self, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.05, to=1.65,length.out=10))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p1 <- ggboxplot(auc_self, x = "group1", y = "auprc", fill = "group1",
                palette = "jco",width = 0.2)+ 
  # geom_hline(yintercept =0.5,color='#dbdcdc')+
  # geom_hline(yintercept =0.6,color='#ffd09a')+
  # geom_hline(yintercept =0.7,color='#ffcbd8')+
  # geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  theme_bw() +
  # stat_compare_means()+
  ylim(0.15,1.68)+
  theme(legend.position="none")+    
  ylab("External AUPRC")+xlab('')+
  ggtitle('Disease category')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p1)

#data_type
{
  auc_all <- all_measure
  auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
  auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
  auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
  auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
  auc0 <- auc_all
  auc0 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
  self.e <- subset(auc0,method=='external')
  self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune","Liver"))
}

# save(self.e,file='01_data/plot_data/F1E.RData')

stat.test <- compare_means(
  auprc~level,data = self.e, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.2, to=1.65,length.out=3))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p2 <- ggboxplot(self.e, x = "level", y = "auprc", fill = "level",
                width = 0.2,palette = c('#774ec7','#bd93cc','#a2c4b1'))+
  # geom_hline(yintercept =0.5,color='#dbdcdc')+
  # geom_hline(yintercept =0.6,color='#ffd09a')+
  # geom_hline(yintercept =0.7,color='#ffcbd8')+
  # geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  theme_bw() +
  ylim(0.15,1.68)+
  # ylim(0.05,1.68)+
  # stat_compare_means()+
  theme(legend.position="none")+    
  ylab("External AUPRC")+xlab('')+
  ggtitle('Data type')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p2 <- p2+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p2)

pB=ggarrange(p1,p2,
             ncol = 2, nrow = 1,
             widths = c(4,3)
)
pB

pdf("./02_figure/FigS2B.pdf", height = 4.5, width = 5)
pB
dev.off()

## FS2 C ----
load('./01_data/all_measure_lasso.RData')
{
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
}

{
  a=all_measure
  a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
  a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
  a$method[a$method == 'self'] <- 'internal'
  auc_self <- subset(a,method == 'internal')
}
save(auc_self,file = '01_data/plot_data/FS2C.RData')

stat.test <- compare_means(
  mcc~group1,data = auc_self, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=0.9, to=1.65,length.out=10))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p1 <- ggboxplot(auc_self, x = "group1", y = "mcc", fill = "group1",
                palette = "jco",width = 0.2)+ 
  # geom_hline(yintercept =0.5,color='#dbdcdc')+
  # geom_hline(yintercept =0.6,color='#ffd09a')+
  # geom_hline(yintercept =0.7,color='#ffcbd8')+
  # geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  # stat_compare_means()+
  # ylim(0.05,1.68)+
  ylab("Internal MCC")+xlab('')+
  ggtitle('Disease category')+
  theme_bw() +
  theme(legend.position="none")+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p1)

#data_type
{
  auc_all <- all_measure
  auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
  auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
  auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
  auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
  auc0 <- auc_all
  auc0 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
  self.e <- subset(auc0,method=='self')
  self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune","Liver"))
}

# save(self.e,file='01_data/plot_data/F1E.RData')

stat.test <- compare_means(
  mcc~level,data = self.e, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.2, to=1.65,length.out=3))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p2 <- ggboxplot(self.e, x = "level", y = "mcc", fill = "level",
                width = 0.2,palette = c('#774ec7','#bd93cc','#a2c4b1'))+
  # geom_hline(yintercept =0.5,color='#dbdcdc')+
  # geom_hline(yintercept =0.6,color='#ffd09a')+
  # geom_hline(yintercept =0.7,color='#ffcbd8')+
  # geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  # ylim(0.05,1.68)+
  # stat_compare_means()+
  theme_bw() +
  theme(legend.position="none")+    
  ylab("Internal MCC")+xlab('')+
  ggtitle('Data type')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p2 <- p2+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p2)

pC=ggarrange(p1,p2,
             ncol = 2, nrow = 1,
             widths = c(4,3)
)
pC
pdf("./02_figure/FigS2C.pdf", height = 4.5, width = 5)
pC
dev.off()



## FS2 D ----
load('./01_data/all_measure_lasso.RData')
{
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
}
{
  a=all_measure
  a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
  a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
  a$method[a$method == 'self'] <- 'internal'
  auc_self <- subset(a,method == 'external')
}
save(auc_self,file = '01_data/plot_data/FS2D.RData')

stat.test <- compare_means(
  mcc~group1,data = auc_self, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=0.9, to=1.65,length.out=10))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p1 <- ggboxplot(auc_self, x = "group1", y = "mcc", fill = "group1",
                palette = "jco",width = 0.2)+ 
  # geom_hline(yintercept =0.5,color='#dbdcdc')+
  # geom_hline(yintercept =0.6,color='#ffd09a')+
  # geom_hline(yintercept =0.7,color='#ffcbd8')+
  # geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  # stat_compare_means()+
  # ylim(0.05,1.68)+
  ylab("External MCC")+xlab('')+
  ggtitle('Disease category')+
  theme_bw() +
  theme(legend.position="none")+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p1)

#data_type
{
  auc_all <- all_measure
  auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
  auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
  auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
  auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
  auc0 <- auc_all
  auc0 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
  self.e <- subset(auc0,method=='external')
  self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune","Liver"))
}

# save(self.e,file='01_data/plot_data/F1E.RData')

stat.test <- compare_means(
  mcc~level,data = self.e, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.2, to=1.65,length.out=3))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p2 <- ggboxplot(self.e, x = "level", y = "mcc", fill = "level",
                width = 0.2,palette = c('#774ec7','#bd93cc','#a2c4b1'))+
  # geom_hline(yintercept =0.5,color='#dbdcdc')+
  # geom_hline(yintercept =0.6,color='#ffd09a')+
  # geom_hline(yintercept =0.7,color='#ffcbd8')+
  # geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  # ylim(0.05,1.68)+
  # stat_compare_means()+
  theme_bw() +
  theme(legend.position="none")+    
  ylab("External MCC")+xlab('')+
  ggtitle('Data type')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p2 <- p2+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p2)

pD=ggarrange(p1,p2,
             ncol = 2, nrow = 1,
             widths = c(4,3)
)
pD
pdf("./02_figure/FigS2D.pdf", height = 4.5, width = 5)
pD
dev.off()


