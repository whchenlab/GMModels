### Figure S9 ----


## FS9A ----
auc_marker_all <- get(load('01_data/MSI_all.cutoff.RData'))
auc_marker_all$dis[auc_marker_all$dis==Inf] <- 0

for (i in unique(auc_marker_all$cutoff)){
  dis_median <- median(subset(auc_marker_all,cutoff==i)$dis)
  if (i==unique(auc_marker_all$cutoff)[1]){
    auc_marker_all_median <- data.frame(dis=dis_median,a='t',cutoff=i)
  }else{
    auc_marker_all_median <- rbind(auc_marker_all_median,data.frame(dis=dis_median,a='t',cutoff=i))
  } 
}
auc_marker_all_median$cutoff <- seq_len(length(unique(auc_marker_all$cutoff)))

auc_marker_all$cutoff <- sprintf("%0.1f", auc_marker_all$cutoff)
class(auc_marker_all$cutoff)

save(auc_marker_all,file = '01_data/plot_data/FS9A.RData')

p<-ggboxplot( auc_marker_all,x='cutoff',y='dis',color=pal_npg("nrc")(10)[4],add = "jitter",width = 0.3)+
  xlab('Cutoff of LDA Score')+
  ylab('MSI (median)')+
  geom_line(data = auc_marker_all_median,aes(x=cutoff,y=dis),color=pal_npg("nrc")(10)[1])+
  geom_point(data = auc_marker_all_median,aes(x=cutoff,y=dis),color=pal_npg("nrc")(10)[1])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
p

pdf("./02_figure/FigS9A.pdf", height = 3.5, width = 5)
p
dev.off()

## FS9B ----
marker.lda.all <- get(load("01_data/marker.lda.all.0.RData"))
marker.lda.all <- subset(marker.lda.all,class=='adjust')
marker.lda.all$LDA <- as.numeric(marker.lda.all$LDA)
marker.lda.all$Class <- ifelse(marker.lda.all$LDA>=0,'Control','Case')
marker.lda.all$LDA <- abs(marker.lda.all$LDA)
marker.lda.all$group1 <- factor(marker.lda.all$group1 ,levels = c("Intestinal","Metabolic","Mental","Autoimmune","Liver"))
library(ggpubr)
KW <- (compare_means(LDA~group1,marker.lda.all,method = 'kruskal.test'))
KW$p.adj

library(ggsci)
p<-ggdensity(marker.lda.all, 'LDA', color="group1",palette = "jco",alpha = 0.1,add='mean',size=1,fill ="group1",rug = TRUE)+
  # stat_compare_means(data=marker.lda.all,aes(color="group1",y="LDA"),inherit.aes = F)+
  labs(x = 'LDA Score',y='Density')+
  theme(legend.position = 'right' )
  # annotate("text", label = paste0('Kruskal−Wallis, p = ',KW$p.adj),x=1.5,y=0.5,size=4)

p
pdf("./02_figure/FigS9B.pdf", height = 3.5, width = 5)
p
dev.off()



## FS9C ----
auc_marker_all_more2 <- get(load("01_data/auc_marker_all_more2.RData"))
auc_marker_all_less2 <- get(load("01_data/auc_marker_all_less2.RData"))
auc_marker_all <- rbind(auc_marker_all_less2,auc_marker_all_more2)
# auc_marker_all <- subset(auc_marker_all,cutoff==index_cutoff)
result_MSI_auc_median <- get(load("01_data/result_MSI_auc_median.RData"))
result_MSI_auc_median$level <- paste0(result_MSI_auc_median$data_type,'_',result_MSI_auc_median$taxon)
auc_marker_all$AUC <- 0
# result_MSI_auc_median
for (i in 1:nrow(auc_marker_all)) {
  auc_marker_all_t <- subset(result_MSI_auc_median,
                             (disease==auc_marker_all[i,"disease"])&(level==auc_marker_all[i,"level"]))
  auc_marker_all$AUC[i] <- auc_marker_all_t$auc
}

auc_marker_all$dis <- as.numeric(auc_marker_all$dis)
auc_marker_all$AUC <- as.numeric(auc_marker_all$AUC)
auc_marker_all$cutoff <- as.numeric(auc_marker_all$cutoff)
save(auc_marker_all,file = '01_data/plot_data/FS9C.RData')

# Cutoff 0.6 ----
auc_marker_all$cutoff[140]
auc_marker_all_temp <- subset(auc_marker_all,cutoff==auc_marker_all$cutoff[140])
nrow(auc_marker_all_temp)
library(ggplot2)
library(ggrepel)
library(ggpubr)

r = paste("cor: ",round(cor(auc_marker_all_temp$dis,auc_marker_all_temp$AUC,method = 'spearman'),3), sep = "")
p_v = paste("p: ",format(cor.test(auc_marker_all_temp$dis,auc_marker_all_temp$AUC,method = 'spearman')$p.value,scientific=T,digits=3), sep = "")
title1=paste(r," ",p_v, sep = "")
title1

p1=ggscatter(data=auc_marker_all_temp,x='dis',y='AUC',shape='level',color='disease',palette = "igv",alpha=0.8,size=2)+
  geom_smooth(mapping=aes(x=dis,y=AUC),data=auc_marker_all_temp,formula =y~x, method=lm, inherit.aes = F,size=0.5)+
  geom_text_repel(data=auc_marker_all_temp,aes(x=dis,y=AUC,label = disease),inherit.aes = F,size=3)+
  xlab("MSI (median)") + ylab("External AUC (median)")+
  labs(color = "disease")+
  theme(legend.position = 'none')+
  annotate('text',x=0.7,y=0.9,label=title1,size=4)+
  border()
p1 <- p1+ggtitle('0.6')
p1


# Cutoff 1.2 ----
auc_marker_all_temp <- subset(auc_marker_all,cutoff==1.2)
colnames(auc_marker_all_temp)
library(ggplot2)
library(ggrepel)
library(ggpubr)

r = paste("cor: ",round(cor(auc_marker_all_temp$dis,auc_marker_all_temp$AUC,method = 'spearman'),3), sep = "")
p_v = paste("p: ",format(cor.test(auc_marker_all_temp$dis,auc_marker_all_temp$AUC,method = 'spearman')$p.value,scientific=T,digits=3), sep = "")
title1=paste(r," ",p_v, sep = "")
title1

p2=ggscatter(data=auc_marker_all_temp,x='dis',y='AUC',shape='level',color='disease',palette = "igv",alpha=0.8,size=2)+
  geom_smooth(mapping=aes(x=dis,y=AUC),data=auc_marker_all_temp,formula =y~x, method=lm, inherit.aes = F,size=0.5)+
  geom_text_repel(data=auc_marker_all_temp,aes(x=dis,y=AUC,label = disease),inherit.aes = F,size=3)+
  xlab("MSI (median)") + ylab("External AUC (median)")+
  labs(color = "disease")+
  theme(legend.position = 'none')+
  annotate('text',x=0.7,y=0.9,label=title1,size=4)+
  border()
p2 <- p2+ggtitle('1.2')
p2

# Cutoff 1.8 ----
auc_marker_all_temp <- subset(auc_marker_all,cutoff==1.8)
colnames(auc_marker_all_temp)
library(ggplot2)
library(ggrepel)
library(ggpubr)

r = paste("cor: ",round(cor(auc_marker_all_temp$dis,auc_marker_all_temp$AUC,method = 'spearman'),3), sep = "")
p_v = paste("p: ",format(cor.test(auc_marker_all_temp$dis,auc_marker_all_temp$AUC,method = 'spearman')$p.value,scientific=T,digits=3), sep = "")
title1=paste(r," ",p_v, sep = "")
title1

p3=ggscatter(data=auc_marker_all_temp,x='dis',y='AUC',shape='level',color='disease',palette = "igv",alpha=0.8,size=2)+
  geom_smooth(mapping=aes(x=dis,y=AUC),data=auc_marker_all_temp,formula =y~x, method=lm, inherit.aes = F,size=0.5)+
  geom_text_repel(data=auc_marker_all_temp,aes(x=dis,y=AUC,label = disease),inherit.aes = F,size=3)+
  xlab("MSI (median)") + ylab("External AUC (median)")+
  labs(color = "disease")+
  theme(legend.position = 'none')+
  annotate('text',x=0.7,y=0.9,label=title1,size=4)+
  border()
p3 <- p3+ggtitle('1.8')
p3

library(patchwork)
p=p1+p2+p3
p

pdf("./02_figure/FigS9C.pdf",width = 11,height = 3.8)
p
dev.off()


