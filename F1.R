### Figure 1 

{
  library(tidyr)
  library(tidyverse)
  library(dplyr)
}
## F1A ----
# statistic for NO. of samples
load('./01_data/cohort_info.RData')
project_stat0 <- project_stat
save(project_stat0,file='01_data/plot_data/F1A.RData')
{
  num_case=sum(project_stat0$case)
  num_control <- 0
  for (i in unique(project_stat0$project)){
    num_control=max(subset(project_stat0,project==i)$control)+num_control
  }
  message('n>=2',',',num_case,'cases',',',num_control,'control',',','total ',num_case+num_control,'samples',',',length(unique(project_stat0$project)),'projects',',',length(unique(project_stat0$disease)),'diseases') 
  
  project_stat0 <- subset(project_stat0,num>2)
  num_case=sum(project_stat0$case)
  num_control <- 0
  for (i in unique(project_stat0$project)){
    num_control=max(subset(project_stat0,project==i)$control)+num_control
  }
  message('n>=3',',',num_case,'cases',',',num_control,'control',',','total ',num_case+num_control,'samples',',',length(unique(project_stat0$project)),'projects',',',length(unique(project_stat0$disease)),'diseases') 
  
  project_stat0 <- subset(project_stat0,num>4)
  num_case=sum(project_stat0$case)
  num_control <- 0
  for (i in unique(project_stat0$project)){
    num_control=max(subset(project_stat0,project==i)$control)+num_control
  }
  message('n>=5',',',num_case,'cases',',',num_control,'control',',','total ',num_case+num_control,'samples',',',length(unique(project_stat0$project)),'projects',',',length(unique(project_stat0$disease)),'diseases') 
}



{
  message(length(unique(subset(project_stat,data_type=='Amplicon')$project)),
          ' 16S studies on ',nrow(subset(pro,data_type=='Amplicon')),
          ' deseases (sample-level: ',sum(subset(project_stat,data_type=='Amplicon')$case),' cases, ',
          sum(unique(subset(project_stat,data_type=='Amplicon')[,c('project','control')])$control),' controls)')
  message(length(unique(subset(project_stat,data_type=='Metagenomics')$project)),
          ' 16S studies on ',nrow(subset(pro,data_type=='Metagenomics')),
          ' deseases (sample-level: ',sum(subset(project_stat,data_type=='Metagenomics')$case),' cases, ',
          sum(unique(subset(project_stat,data_type=='Metagenomics')[,c('project','control')])$control),' controls)')
  
}

## F1B ----
# data
temp <- pro
colnames(temp)[which(colnames(temp)=='group1')] <- 'type'
colnames(temp)[which(colnames(temp)=='num')] <- 'count'
save(temp,file='01_data/plot_data/F1B.RData')

#plot
p=temp %>%
  mutate(type=factor(type, levels=c("Intestinal" ,"Metabolic" , "Mental"  ,  "Autoimmune", "Liver"))) %>%
  arrange(desc(count)) %>% 
  mutate(disease=factor(disease, levels=unique(disease))) %>%
  ggplot(aes(x=disease, y=count,group=data_type)) +

  geom_bar(stat="identity",position='stack', aes(fill=data_type)) +

  geom_text(aes(label=count),position=position_stack(vjust = 0.5),size=7)+
  facet_grid(~type, scales="free", space="free") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1,face = 'bold',size=12),
        axis.text.y=element_text(face = 'bold',size=12),
        plot.title=element_text(hjust=0.5)) +
  ylab("No. of project") +
  xlab('disease') +
  coord_cartesian(ylim=c(0,11),expand=FALSE) +
  scale_y_continuous(breaks=seq(0, 12, 2))+
  theme(panel.border = element_blank(), axis.line = element_line())+
  scale_fill_d3(alpha = 0.5)+
  theme(text = element_text(size=16,face = 'plain',family ='',colour = 'black'))
print(p)

pdf("./02_figure/Fig1B.pdf", height = 10, width = 12)
p
dev.off()


## F1C ----
project_stat0 <- project_stat
project_stat0$project <- 1:nrow(project_stat0)
save(project_stat0,file='01_data/plot_data/F1C.RData')
#plot
library(tidyr)
project_stat0 <-  gather(project_stat0,phenotype,num,c('case','control'))
project_stat0$phenotype <- factor(project_stat0$phenotype,levels = c('control','case'))
p<-ggdensity(project_stat0, 'num', color="phenotype",palette = "aaas",add = "median",alpha = 0.1,size=1,fill ="phenotype",rug = TRUE)+
  labs(x = 'No. of samples in each cohort',y='Density')+
  annotate("text", label = paste0("Median: ",median(subset(project_stat0,phenotype=='case')$num)), x = 150, y = 0.015, size = 4, colour = pal_aaas("default", alpha = 0.6)(10)[2])+
  annotate("text", label = paste0("Median: ",median(subset(project_stat0,phenotype=='control')$num)), x = 150, y = 0.013, size = 4, colour = pal_aaas("default", alpha = 0.6)(10)[1]) 
p
pdf("./02_figure/Fig1C.pdf", height = 5, width = 6)
p
dev.off()


## F1DE ----
{
a=get(load('./01_data/AUC_all_lasso.RData'))
a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
a$method[a$method == 'self'] <- 'internal'
auc_self <- subset(a,method == 'internal')
}
save(auc_self,file = '01_data/plot_data/F1D.RData')

stat.test <- compare_means(
  auc~group1,data = auc_self, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.05, to=1.65,length.out=10))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p1 <- ggboxplot(auc_self, x = "group1", y = "auc", fill = "group1",
               palette = "jco",width = 0.2)+ 
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  geom_hline(yintercept =0.9,color='#e60020')+
  # stat_compare_means()+
  ylim(0.05,1.68)+
  theme(legend.position="none")+    
  ylab("Internal AUC")+xlab('')+
  ggtitle('Disease category')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p1 <- p1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p1)

#data_type
{
auc_all <- get(load('01_data/AUC_all_lasso.RData'))
auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
auc0 <- auc_all
auc0 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
self.e <- subset(auc0,method=='self')
self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune","Liver"))
}

save(self.e,file='01_data/plot_data/F1E.RData')

stat.test <- compare_means(
  auc~level,data = self.e, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=1.2, to=1.65,length.out=3))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p2 <- ggboxplot(self.e, x = "level", y = "auc", fill = "level",
               width = 0.2,palette = c('#774ec7','#bd93cc','#a2c4b1'))+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  geom_hline(yintercept =0.9,color='#e60020')+
  ylim(0.05,1.68)+
  # stat_compare_means()+
  theme(legend.position="none")+    
  ylab("Internal AUC")+xlab('')+
  ggtitle('Data type')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p2 <- p2+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p2)

p=ggarrange(p1,p2,
            ncol = 2, nrow = 1,
            widths = c(4,3)
)
p
pdf("./02_figure/Fig1DE.pdf", height = 4.5, width = 5)
p
dev.off()



## F1FG -----
{
  a=get(load('./01_data/AUC_all_lasso.RData'))
  
  a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
  a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
  a$method[a$method == 'self'] <- 'internal'
  
  
  auc_self <- subset(a,method == 'internal')
  auc_external <- subset(a,method=='external')
  a1 <- a
  a1$group1 <- 'All'
  a_all <- rbind(a,a1)
  
  mean(auc_self$auc)
  max(auc_self$auc)
  min(auc_self$auc)
  sd(auc_self$auc)
  
  mean(auc_external$auc)
  max(auc_external$auc)
  min(auc_external$auc)
  sd(auc_external$auc)
  
  a_all$group1 <- factor(a_all$group1,levels = c("All","Intestinal","Metabolic","Mental","Autoimmune",'Liver'))
  
}
save(a_all,file='01_data/plot_data/F1FG.RData')
p <- ggboxplot(a_all, x = "method", y = "auc", fill = "method",
               palette = c('#1fb8b4','#ff7f0e'),width = 0.15)+ 
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  geom_hline(yintercept =0.9,color='#e60020')+
  facet_wrap(~group1,nrow = 1)+
  annotate('text',x=1:2,y=0.15,label=c('0.765','0.638'))+ #AUC
  geom_signif(comparisons =list(c('internal','external')),y_position = c(1.12, 1.32),test = 'wilcox.test',map_signif_level = function(x){ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')})+
  ylim(0.05,1.32)+
  theme(legend.position="top")+    
  xlab("") + ylab("AUC")+
  labs(fill = "AUC type")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),axis.text.x = element_blank(),axis.ticks=element_blank())
print(p)

pdf("./02_figure/Fig1FG.pdf", height = 4.2, width = 8.5)
p
dev.off()

