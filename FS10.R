### Figure S10  ----

library(ggplot2)
library(ggrepel)
library(ggsci)
library(dplyr)
# FS10A ----

load('01_data/curatedMetagenomicData/All_disease_info.RData')
All_disease_info$disease[All_disease_info$disease=='adenoma'] <- 'Adenoma'
All_disease_info_stat <- data.frame(table(All_disease_info$disease))
All_disease_info_stat[All_disease_info_stat$Var1=='adenoma','Var1'] <- 'Adenoma'
save(All_disease_info_stat,file='01_data/plot_data/FS10A.RData')


pA <- ggplot(All_disease_info_stat,aes(x=Var1,y=Freq)) +
  geom_col(aes(x=Var1,y=Freq,fill=Var1),width = 0.5) +
  geom_text(aes(label=Freq),position=position_stack(vjust = 0.5),size=7)+
  # facet_grid(~type, scales="free", space="free") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1,size=12),
        # axis.text.y=element_text(face = 'bold',size=12),
        plot.title=element_text(hjust=0.5)) +
  ylab("No. of project") +
  xlab('Disease') +
  labs(fill='Disease')+
  # coord_cartesian(ylim=c(0,11),expand=FALSE) +
  # scale_y_continuous(breaks=seq(0, 12, 2))+
  theme(panel.border = element_blank(), axis.line = element_line())+
  scale_fill_d3(alpha = 0.5)+
  theme(text = element_text(size=16,face = 'plain',family ='',colour = 'black'))
print(pA)


pdf("./02_figure/FigS10A.pdf", height = 3.8, width = 4.5)
pA
dev.off()



# FS10B ----
result <- data.frame()
# load('01_data/result_taxon_vs_pathway.RData')
load('01_data/result_taxon_vs_pathway.excluded.RData')
for (i in names(result_taxon_vs_pathway)){
  temp_result <- result_taxon_vs_pathway[[i]]
  # result <- rbind(result,temp_result$taxon$AUC, temp_result$pathway$AUC, temp_result$coverage$AUC)
  result <- rbind(result,temp_result$taxon$AUC, temp_result$pathway$AUC, temp_result$both$AUC)
}
result$method[result$method=='self'] <- 'Internal'
result$method[result$method=='external'] <- 'External'

table(result[,c('method','Label')])
result$Label <- factor(result$Label,levels = c('Taxon','Pathway','Both'))
result$method <- factor(result$method,levels = c('Internal','External'))
save(result,file='01_data/plot_data/FS10BC.RData')


stat.test1 <- compare_means(
  auc~Label,data = subset(result,method=='Internal'), 
  method = "wilcox.test",
  paired = T
) %>%mutate(y.position = rep(seq(from=1.1, to=1.5,length.out=3), times=1))
stat.test2 <- compare_means(
  auc~Label,data = subset(result,method=='External'), 
  method = "wilcox.test",
  paired = T
) %>%mutate(y.position = rep(seq(from=1.1, to=1.5,length.out=3), times=1))



x <- stat.test1$p.adj
stat.test1$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

x <- stat.test2$p.adj
stat.test2$p.adj.signif<-ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')

pB1 <- ggboxplot(subset(result,method=='Internal'), x = "Label", y = "auc",color = 'black',fill = "Label", palette = "aaas",
                 width = 0.3,add.params = list(size=1))+
  ylim(0.25,1.55)+
  # facet_grid(~method)+
  stat_compare_means(aes(label = paste0('p = ',..p.format..)),label.y = 1.5,hjust=-0.01)+
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1,size=12),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+ 
  xlab("") + ylab("Internal AUC")
# labs(color = "disease type")+
pB1 = pB1 + stat_pvalue_manual(stat.test1,label = "p.adj.signif")
print(pB1)

pB2 <- ggboxplot(subset(result,method=='External'), x = "Label", y = "auc",color = 'black',fill = "Label", palette = "aaas", 
                 width = 0.3,add.params = list(size=1))+
  ylim(0.25,1.55)+
  # facet_grid(~method)+
  stat_compare_means(aes(label = paste0('p = ',..p.format..)),label.y = 1.5,hjust=-0.01)+
  # axis.text.x = element_blank(),axis.ticks=element_blank()
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1,size=12),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+ 
  xlab("") + ylab("External AUC")
# labs(color = "disease type")+
pB2 = pB2 + stat_pvalue_manual(stat.test2,label = "p.adj.signif")
print(pB2)


library(patchwork)
pB <- pB1+pB2
# 
library(ggplotify)
pB <- as.ggplot(pB)

print(pB)
pdf("./02_figure/FigS10BC.pdf", height = 3.8, width = 7)
pB
dev.off()


