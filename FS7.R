### Figure S7 ----

# Fig S7A ----
result_MSI <- get(load('01_data/result_MSI_DA.RData'))
result_MSI <- subset(result_MSI,Group=='cutoff')
auc_all <- get(load("01_data/AUC_all_lasso.RData"))
auc_all$TR <- paste(auc_all$TR,auc_all$TE,sep='_')
result_MSI_mat <- unique(result_MSI[,c('disease','data_type','taxon')])
colnames(result_MSI)
result_MSI_auc_all <- data.frame()
result_MSI_auc_median <- data.frame()
for (i in seq_len(nrow(result_MSI_mat))) {
  result_MSI_temp <- subset(result_MSI,(disease==result_MSI_mat$disease[i]) & (data_type==result_MSI_mat$data_type[i]) & (taxon==result_MSI_mat$taxon[i]))
  result_MSI_temp$auc <- 0
  for (j in seq_len(nrow(result_MSI_temp))) {
    auc_all_temp <- subset(auc_all,(TR==result_MSI_temp[j,'TR'])&(disease==result_MSI_mat$disease[i]) & (data_type==result_MSI_mat$data_type[i]) & (taxon==result_MSI_mat$taxon[i]))
    result_MSI_temp$auc[j] <- auc_all_temp$auc
  }
  result_MSI_auc_all <- rbind(result_MSI_auc_all,result_MSI_temp)
  result_MSI_auc_median <- rbind(result_MSI_auc_median,(cbind(unique(result_MSI_temp[,c("disease","data_type","group1","taxon")]),
                                                              ALDEx2=median(result_MSI_temp$ALDEx2),
                                                              MaAslin2=median(result_MSI_temp$MaAslin2),
                                                              LEfSe=median(result_MSI_temp$LEfSe),
                                                              auc=median(result_MSI_temp$auc))))
}
save(result_MSI_auc_all,file = '01_data/result_MSI_auc_all.RData')
save(result_MSI_auc_median,file = '01_data/result_MSI_auc_median.RData')

#plot
auc_marker0 <- get(load('01_data/result_MSI_auc_median.RData'))
auc_marker0$disease_type2 <- ifelse(auc_marker0$group1=='Intestinal','Intestinal','non-Intestinal')
auc_marker0$level <- paste0(auc_marker0$data_type,'_',auc_marker0$taxon)
auc_marker0$level[auc_marker0$level=='Metagenomics_genus']	<- 'mNGS_genus'
auc_marker0$level[auc_marker0$level=='Metagenomics_species']	<- 'mNGS_species'
auc_marker0$level[auc_marker0$level=='Amplicon_genus']	<- '16S_genus'
auc_marker0$level <- factor(auc_marker0$level,levels = c('mNGS_species','mNGS_genus','16S_genus'))

#
index <- 'ALDEx2'
auc_marker <- auc_marker0[c(index,"disease","data_type","group1","taxon","auc","level" ,"disease_type2")]
colnames(auc_marker)[1] <- 'dis'
save(auc_marker,file = '01_data/plot_data/FS7A.RData')

{
  library(psych)
  cor <- corr.test(auc_marker[,c('dis','auc')],method = 'spearman')
  cor$p[1,2]
  cor$r[1,2]
}
{
  require(psych)
  require(ggcorrplot)
  library(showtext)
  library(ggpubr)
  library(ggrepel)
}
p1_1=ggscatter(data=auc_marker,x='dis',y='auc',shape='level',color='disease',palette = "igv",alpha=0.8,size=2)+
  geom_smooth(mapping=aes(x=dis,y=auc),data=auc_marker,formula =y~x, method=lm, inherit.aes = F,size=0.5)+
  geom_text_repel(data=auc_marker,aes(x=dis,y=auc,label = disease),inherit.aes = F,size=3)+
  xlab("MSI (median)") + ylab("External AUC (median)")+
  labs(color = "disease")+
  theme(legend.position = 'none')+
  annotate('text',x=15,y=1,label=paste0('r=',round(cor$r[1,2],2),'  ','p','=',format(cor$p[1,2],scientific=T,digits=3)),size=4)+
  border()
p1_1 <- p1_1+ggtitle('ALDEx2')

#MaAslin2
index <- 'MaAslin2'
auc_marker <- auc_marker0[c(index,"disease","data_type","group1","taxon","auc","level" )]
colnames(auc_marker)[1] <- 'dis'
{
  library(psych)
  cor <- corr.test(auc_marker[,c('dis','auc')],method = 'spearman')
  cor$p[1,2]
  cor$r[1,2]
}
{
  require(psych)
  require(ggcorrplot)
  library(showtext)
  library(ggpubr)
  library(ggrepel)
}
p2_1=ggscatter(data=auc_marker,x='dis',y='auc',shape='level',color='disease',palette = "igv",alpha=0.8,size=2)+
  geom_smooth(mapping=aes(x=dis,y=auc),data=auc_marker,formula =y~x, method=lm, inherit.aes = F,size=0.5)+
  geom_text_repel(data=auc_marker,aes(x=dis,y=auc,label = disease),inherit.aes = F,size=3)+
  xlab("MSI (median)") + ylab("External AUC (median)")+
  labs(color = "disease")+
  theme(legend.position = 'none')+
  annotate('text',x=40,y=1,label=paste0('r=',round(cor$r[1,2],2),'  ','p','=',format(cor$p[1,2],scientific=T,digits=3)),size=4)+
  border()
p2_1 <- p2_1+ggtitle('MaAslin2')

library(patchwork)
pA=p1_1+p2_1

pdf('02_figure/FigS7A.pdf',width = 11,height = 5)
pA
dev.off()

# Fig S7B ----
# ALDEx2
# disease category
result_MSI <- get(load('01_data/result_MSI_DA.RData'))
result_MSI <- subset(result_MSI,Group=='cutoff')
save(result_MSI,file = '01_data/plot_data/FS7B.RData')


result_MSI$level <- paste0(result_MSI$data_type,'_',result_MSI$t)

result_MSI$level[result_MSI$level=='Metagenomics_species'] <- 'mNGS_species'
result_MSI$level[result_MSI$level=='Metagenomics_genus'] <- 'mNGS_genus'
result_MSI$level[result_MSI$level=='Amplicon_genus'] <- '16S_genus'
result_MSI$level <- factor(result_MSI$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))

colnames(result_MSI)
# [1] "TR"        "ALDEx2"    "LEfSe"     "MaAslin2"  "Group"     "disease"   "data_type" "group1"    "taxon"     "level" 

{
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
}
{
  result_MSI$level[result_MSI$level=='Metagenomics_genus'] <- 'mNGS_genus'
  result_MSI$level[result_MSI$level=='Metagenomics_species'] <- 'mNGS_species'
}


stat.test <- compare_means(
  ALDEx2~group1,data = result_MSI, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=60, to=100,length.out=10))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p3_1 <- ggboxplot(result_MSI, x = "group1", y = "ALDEx2", fill = "group1",
                  palette = "jco",width = 0.2,outlier.shape = NA)+ 
  theme_bw() +
  stat_compare_means()+
  ylim(0,110)+
  theme(legend.position="none")+    
  ylab("MSI (ALDEx2)")+xlab('')+
  ggtitle('Disease category')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p3_1 <- p3_1+stat_pvalue_manual(stat.test,label = "p.adj.signif")

print(p3_1)




# data type
result_MSI <- get(load('01_data/result_MSI_DA.RData'))
result_MSI <- subset(result_MSI,Group=='cutoff')
result_MSI$level <- paste0(result_MSI$data_type,'_',result_MSI$t)
colnames(result_MSI)

result_MSI <- subset(result_MSI,disease %in% names(which(table(unique(result_MSI[,c('disease','data_type')])$disease)>1)))
result_MSI$group1_new <- factor(result_MSI$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune"))

result_MSI$level[result_MSI$level=='Metagenomics_species'] <- 'mNGS_species'
result_MSI$level[result_MSI$level=='Metagenomics_genus'] <- 'mNGS_genus'
result_MSI$level[result_MSI$level=='Amplicon_genus'] <- '16S_genus'
result_MSI$level <- factor(result_MSI$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))



stat.test <- compare_means(
  ALDEx2~level,data = result_MSI, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=70, to=90,length.out=3))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p4_1 <- ggboxplot(result_MSI, x = "level", y = "ALDEx2", fill = "level",
                  width = 0.2,palette = c('#774ec7','#bd93cc','#a2c4b1'),outlier.shape = NA)+
  ylim(0,110)+
  theme_bw() +
  stat_compare_means()+
  theme(legend.position="none")+    
  ylab("MSI (ALDEx2)")+xlab('')+
  ggtitle('Data type')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p4_1 <- p4_1+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p4_1)


library(patchwork)
pB <- p3+p4
pB

pB_1=ggarrange(p3_1,p4_1,
               ncol = 2, nrow = 1,
               widths = c(4,3)
)
pB_1

pdf("./02_figure/FigS7B.pdf", height = 4.5, width = 7)
pB_1
dev.off()



# Fig S7C ----

# MaAslin2

# disease category all
result_MSI <- get(load('01_data/result_MSI_DA.RData'))
result_MSI$ALDEx2 <- result_MSI$MaAslin2
result_MSI <- subset(result_MSI,Group=='cutoff')
save(result_MSI,file = '01_data/plot_data/FS7C.RData')


result_MSI$level <- paste0(result_MSI$data_type,'_',result_MSI$t)

colnames(result_MSI)

result_MSI$level[result_MSI$level=='Metagenomics_species'] <- 'mNGS_species'
result_MSI$level[result_MSI$level=='Metagenomics_genus'] <- 'mNGS_genus'
result_MSI$level[result_MSI$level=='Amplicon_genus'] <- '16S_genus'
result_MSI$level <- factor(result_MSI$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))

# [1] "TR"        "ALDEx2"    "LEfSe"     "MaAslin2"  "Group"     "disease"   "data_type" "group1"    "taxon"     "level" 

{
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(ggpubr)
}
{
  result_MSI$level[result_MSI$level=='Metagenomics_genus'] <- 'mNGS_genus'
  result_MSI$level[result_MSI$level=='Metagenomics_species'] <- 'mNGS_species'
}

stat.test <- compare_means(
  ALDEx2~group1,data = result_MSI, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=350, to=550,length.out=10))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p5_1 <- ggboxplot(result_MSI, x = "group1", y = "ALDEx2", fill = "group1",
                  palette = "jco",width = 0.2,outlier.shape = NA)+ 
  theme_bw() +
  stat_compare_means()+
  ylim(0,560)+
  theme(legend.position="none")+    
  ylab("MSI (MaAslin2)")+xlab('')+
  ggtitle('Disease category')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p5_1 <- p5_1+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p5_1)



# data type all
result_MSI <- get(load('01_data/result_MSI_DA.RData'))
result_MSI$ALDEx2 <- result_MSI$MaAslin2
result_MSI <- subset(result_MSI,Group=='cutoff')
result_MSI$level <- paste0(result_MSI$data_type,'_',result_MSI$t)
colnames(result_MSI)

result_MSI <- subset(result_MSI,disease %in% names(which(table(unique(result_MSI[,c('disease','data_type')])$disease)>1)))
result_MSI$group1_new <- factor(result_MSI$group1 , levels = c("Intestinal","Metabolic","Mental","Autoimmune"))


result_MSI$level[result_MSI$level=='Metagenomics_species'] <- 'mNGS_species'
result_MSI$level[result_MSI$level=='Metagenomics_genus'] <- 'mNGS_genus'
result_MSI$level[result_MSI$level=='Amplicon_genus'] <- '16S_genus'
result_MSI$level <- factor(result_MSI$level,levels=c("mNGS_species","mNGS_genus","16S_genus"))



stat.test <- compare_means(
  ALDEx2~level,data = result_MSI, 
  # group.by = "level",
  method = "wilcox.test"
)%>%mutate(y.position = seq(from=350, to=500,length.out=3))
x=stat.test$p.adj
stat.test$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

p6_1 <- ggboxplot(result_MSI, x = "level", y = "ALDEx2", fill = "level",
                  width = 0.2,palette = c('#774ec7','#bd93cc','#a2c4b1'),outlier.shape = NA)+
  ylim(0,560)+
  theme_bw() +
  stat_compare_means()+
  theme(legend.position="none")+    
  ylab("MSI (MaAslin2)")+xlab('')+
  ggtitle('Data type')+
  theme(axis.text.x=element_text(angle=20, hjust=0.8,face = 'plain',size=13),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p6_1 <- p6_1+stat_pvalue_manual(stat.test,label = "p.adj.signif")
print(p6_1)



pC <- p5+p6
pC

pC_1=ggarrange(p5_1,p6_1,
               ncol = 2, nrow = 1,
               widths = c(4,3)
)
pC_1

pdf("./02_figure/FigS7C.pdf", height = 4.5, width = 7)
pC_1
dev.off()



# Fig S7D ----
# MSI plot 
result_MSI <- get(load('01_data/result_MSI_DA.RData'))

result_MSI <- subset(result_MSI,Group=='cutoff')

colnames(result_MSI)

save(result_MSI,file = '01_data/plot_data/FS7D.RData')
# [1] "TR"       "ALDEx2"   "LEfSe"    "MaAslin2" "Group" 

library(ggplot2)
library(ggrepel)

r = paste("cor: ",round(cor(result_MSI$LEfSe,result_MSI$ALDEx2,method = 'spearman'),3), sep = "")
p_v = paste("p: ",format(cor.test(result_MSI$LEfSe,result_MSI$ALDEx2,method = 'spearman')$p.value,scientific=T,digits=3), sep = "")
title1=paste(r," ",p_v, sep = "")
title1

r = paste("cor: ",round(cor(result_MSI$LEfSe,result_MSI$MaAslin2,method = 'spearman'),3), sep = "")
p_v = paste("p: ",format(cor.test(result_MSI$LEfSe,result_MSI$MaAslin2,method = 'spearman')$p.value,scientific=T,digits=3), sep = "")
title2=paste(r," ",p_v, sep = "")
title2


library(ggpubr)
p1=ggplot(result_MSI,mapping = aes(x=LEfSe, y=ALDEx2),size=1)+
  geom_point(size=2)+
  labs(x='MSI (LEfSe)',y='MSI (ALDEx2)')+
  color_palette('igv')+
  stat_smooth(mapping =aes(x=LEfSe, y=ALDEx2) ,method=lm,inherit.aes = F,size=0.5)+
  # scale_shape_manual(values = c(4,18,20))+
  theme_classic()+
  annotate('text',x=0.5,y=70,label=title1,size=4)+
  border()
p1 <- p1+ggtitle('ALDEx2')

p2=ggplot(result_MSI,mapping = aes(x=LEfSe, y=MaAslin2),size=1)+
  geom_point(size=2)+
  labs(x='MSI (LEfSe)',y='MSI (MaAslin2)')+
  color_palette('igv')+
  stat_smooth(mapping =aes(x=LEfSe, y=MaAslin2) ,method=lm,inherit.aes = F,size=0.5)+
  # scale_shape_manual(values = c(4,18,20))+
  theme_classic()+
  annotate('text',x=1,y=400,label=title2,size=4)+
  border()
p2 <- p2+ggtitle('MaAslin2')

library(patchwork)
pB <- p1+p2


pdf('02_figure/FigS7D.pdf',width = 9,height = 5)
pB
dev.off()




# Fig S7E ----
library(UpSetR)        
library(VennDiagram) 

result_upset <- get(load('01_data/DA.RData'))
result_upset <- subset(result_upset,top=='N')
result_upset$BH <- as.numeric(result_upset$BH)
# result_upset_lefse <- subset(result_upset,(method=='LEfSe'))
result_upset <- subset(result_upset,(method=='LEfSe')|(BH<0.1))
result_upset$DAfeature <- paste(result_upset$disease,
                                result_upset$data_type,result_upset$taxon,
                                result_upset$project_id,result_upset$DA,
                                sep = '_')


result_upset <- result_upset[,c('method','DAfeature')]
unique(result_upset$method)
#[1] "ALDEX2"   "Maaslin2" "LEfSe"   
result_upset_1 <- subset(result_upset,method=='ALDEX2')
result_upset_2 <- subset(result_upset,method=='Maaslin2')
result_upset_3 <- subset(result_upset,method=='LEfSe')
# max_nrow <- max(nrow(result_upset_1),nrow(result_upset_2),nrow(result_upset_3))
# result_upset_1 <- c(result_upset_1$DAfeature,rep('',max_nrow - nrow(result_upset_1)))
# result_upset_2 <- c(result_upset_2$DAfeature,rep('',max_nrow - nrow(result_upset_2)))
# result_upset_3 <- c(result_upset_3$DAfeature)
# 
# upset_dat <- data.frame('ALDEX2'=result_upset_1,'Maaslin2'=result_upset_2,'LEfSe'=result_upset_3)

# upset_dat <- read.delim('https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/Venn/flower.txt')                      # 这里读取了网络上的demo数据，将此处换成你自己电脑里的文件
upset_list <- list(result_upset_1$DAfeature, result_upset_2$DAfeature, result_upset_3$DAfeature)   # 制作Upset图搜所需要的列表文件
names(upset_list) <- colnames(upset_dat[1:3])    # 把列名赋值给列表的key值
a <- fromList(upset_list)
save(upset_list,file = '01_data/plot_data/FS7E.RData')

pD <- upset(fromList(upset_list),  
            matrix.color = "gray23", main.bar.color = "#3B4992", sets.bar.color = "#3B4992",
            nsets = 10000,     
            nintersects = 10000, 
            order.by = "freq", 
            keep.order = F, 
            mb.ratio = c(0.6,0.4),   
            text.scale = 1.5, 
            mainbar.y.label = "Marker intersections", sets.x.label = "No. of markers per method" #坐标轴名称
)
pD
# library(ggplotify)
# pD <- as.ggplot(pD)+ggtitle('D')+theme(plot.title = element_text(size=18,hjust=-0.1))
# pD


pdf("./02_figure/FigS7E.pdf", height = 5, width = 5.3)
pD
dev.off()

