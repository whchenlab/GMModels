### Figure S7 ----


## FS7A ----

self.e <- get(load('./01_data/MSI_all.RData'))
self.e$data_type[self.e$data_type=='Metagenomics'] <- 'mNGS'
self.e$data_type[self.e$data_type=='Amplicon'] <- '16S'
self.e$level <- paste0(self.e$data_type,'_',self.e$taxon)
self.e$level <- factor(self.e$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))

save(self.e,file = '01_data/plot_data/FS7A.RData')

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



pdf("./02_figure/FigS7A.pdf", height = 5, width = 8)
p
dev.off()

## FS7B ----
# data

self.e <- get(load('./01_data/MSI_all.RData'))
self.e$data_type[self.e$data_type=='Metagenomics'] <- 'mNGS'
self.e$data_type[self.e$data_type=='Amplicon'] <- '16S'
self.e$level <- paste0(self.e$data_type,'_',self.e$taxon)
self.e$level <- factor(self.e$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))

self.e <- subset(self.e,disease %in% names(which(table(unique(self.e[,c('disease','data_type')])$disease)>1)))
self.e$group1 <- factor(self.e$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune"))
colnames(self.e)[colnames(self.e)=='group1']='subtype'

save(self.e,file = '01_data/plot_data/FS7B.RData')

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

pdf("./02_figure/FigS7B.pdf", height = 5, width = 8)
p
dev.off()

## FS7C ----

self.e <- get(load('./01_data/MSI_all.RData'))

self.e$data_type[self.e$data_type=='Amplicon'] <- '16S'
self.e$level <- paste0(self.e$data_type,'_',self.e$taxon)

self.e <- subset(self.e,level=='16S_genus'&disease!='IBS')

save(self.e,file = '01_data/plot_data/FS7C.RData')

stat.test <- compare_means(
  dis_value~group1,data = self.e, group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = c(seq(from=4, to=7.5,length.out=10)))
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



pdf("./02_figure/FigS7C.pdf", height = 5, width = 3)
p
dev.off()
