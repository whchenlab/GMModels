### Figure 3

## F3AB -----
auc0 <- get(load('01_data/AUC_all_lasso.RData'))
auc1 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
auc1 <- subset(auc1,method=='external')

auc1$data_type <- paste(auc1$data_type,auc1$taxon,sep = '_')
auc1$data_type <- factor(auc1$data_type)
mydata <- auc1[,c('group1','data_type','auc')]
colnames(mydata) <- c('disease_type','data_type','auc')

save(mydata,file = '01_data/plot_data/F3AB.RData')

#A
fit <- aov(auc ~ data_type*disease_type,data = mydata)
a<-summary(fit)
a

R2data <- data.frame('type'=c('Data type','Disease category','Interaction'),'R2'=c(a[[1]][["Sum Sq"]][1]/sum(a[[1]][["Sum Sq"]]),a[[1]][["Sum Sq"]][2]/sum(a[[1]][["Sum Sq"]]),a[[1]][["Sum Sq"]][3]/sum(a[[1]][["Sum Sq"]])),
                     p=c(paste0('p < ','2e-16'),paste0('p < ','2e-16'),paste0('p = ','0.109')))
R2data$type <- factor(R2data$type,levels=c('Data type','Disease category','Interaction'))

p <- ggplot(R2data, aes(x=type, y=R2))+
  geom_bar(data=R2data,mapping = aes(x=type, y=R2,fill=type),stat="identity",
           position = 'dodge',width = 0.4,alpha=0.3)+
  scale_fill_manual(values = c("#8c510a","#01665e",'black'))+
  geom_text(data=R2data,aes(label=round(R2,3)), vjust=-0.004)+
  geom_text(data=R2data,aes(label=p), vjust=-1.5)+
  theme(axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  labs( title = "ANOVA: ext−AUC~Data_type*Disease_category")+
  theme(text = element_text(size=13,face = 'plain',family =''))
print(p)

pdf("02_figure/Fig3A.pdf", height = 5, width = 5)
p
dev.off()


#B

fit <- aov(auc ~ disease_type*data_type,data = mydata)
a<-summary(fit)
a

R2data <- data.frame('type'=c('Disease category','Data type','Interaction'),'R2'=c(a[[1]][["Sum Sq"]][1]/sum(a[[1]][["Sum Sq"]]),a[[1]][["Sum Sq"]][2]/sum(a[[1]][["Sum Sq"]]),a[[1]][["Sum Sq"]][3]/sum(a[[1]][["Sum Sq"]])),
                     p=c(paste0('p < ','2e-16'),paste0('p = ','0.163'),paste0('p = ','0.109')))
R2data$type <- factor(R2data$type,levels=c('Disease category','Data type','Interaction'))
p <- ggplot(R2data, aes(x=type, y=R2))+
  geom_bar(data=R2data,mapping = aes(x=type, y=R2,fill=type),stat="identity",
           position = 'dodge',width = 0.4,alpha=0.3)+
  scale_fill_manual(values = c("#8c510a","#01665e",'black'))+
  geom_text(data=R2data,aes(label=round(R2,3)), vjust=-0.004)+
  geom_text(data=R2data,aes(label=p), vjust=-1.5)+
  theme(axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  
  labs( title = "ANOVA: ext−AUC~Disease_category*Data_type")+
  theme(text = element_text(size=13,face = 'plain',family =''))
print(p)


pdf("02_figure/Fig3B.pdf", height = 5, width = 5)
p
dev.off()

## F3C ----
{
  auc_all <- get(load('01_data/AUC_all_lasso.RData'))
  auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
  auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
  auc_all$level[auc_all$level=='Amplicon_genus'] <- '16S_genus'
  auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))
  self.e <- subset(auc_all,method=='external')
  self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune","Liver"))
  self.e <- subset(self.e,!disease%in%c('IBS','Adenoma'))
  self.e <- subset(self.e,level=='16S_genus')
}

save(self.e,file = '01_data/plot_data/F3C.RData')

stat.test <- compare_means(
  auc~group1,data = self.e, group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = c(seq(from=1, to=1.5,length.out=10)))
x <- stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')
# adj
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

pdf("./02_figure/Fig3C.pdf", height = 4.5, width = 3)
p
dev.off()

## F3D ------

auc_all <- get(load('data/auc_all_lasso.RData'))
auc_all$level[auc_all$level=='Metagenomics_species'] <- 'mNGS_species'
auc_all$level[auc_all$level=='Metagenomics_genus'] <- 'mNGS_genus'
auc_all$level <- factor(auc_all$level,levels=c("mNGS_species","mNGS_genus","Amplicon_genus"))

auc0 <- auc_all
auc0=auc0[auc0$group1!='Liver disease',]
auc0 <- subset(auc0,disease %in% names(which(table(unique(auc0[,c('disease','data_type')])$disease)>1)))
self.e <- subset(auc0,method=='external')
self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune"))
colnames(self.e)[8]='subtype'

data<-subset(self.e,(group2=='Intestinal')&(level%in%c('mNGS_species','Amplicon_genus')))

save(data,file = '01_data/plot_data/F3D.RData')

p <- ggboxplot(data, x = "level", y = "auc",color = 'level', palette = c('#774ec7','#a2c4b1'),
               facet.by = "disease",outlier.shape = NA,nrow = 1)+
  geom_jitter(data,mapping=aes(color=level),width =0.2,size=2)+
  scale_shape_manual(values = c(19, 15))+
  stat_compare_means(label.y = 1.05)+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank(),legend.position="top")+ 
  xlab("Data type") + ylab("ext-AUC")+
  labs(shape = "Data type")+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))  
print(p)

pdf("./02_figure/Fig3D.pdf", height = 4.5, width = 5.5)
p
dev.off()


