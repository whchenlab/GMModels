### Figure S8 ----


## FS8A ----
cross.siamcat.all <- get(load('01_data/Greengene_vs_Silva_all.RData'))
cross.siamcat.all$method[cross.siamcat.all$method=='self'] <- 'Internal'
cross.siamcat.all$method[cross.siamcat.all$method=='external'] <- 'External'
cross.siamcat.all$method <- factor(cross.siamcat.all$method,levels = c('Internal','External'))
cross.siamcat.all$Label <- factor(cross.siamcat.all$Label,levels = c('Silva','Greengene'))
cross.siamcat.all$auc <- as.numeric(cross.siamcat.all$auc)
save(cross.siamcat.all,file='01_data/plot_data/FS8A.RData')

library(ggpubr)
p <- ggboxplot(subset(cross.siamcat.all), x = 'Label', y = 'auc', fill = 'Label',palette = c("#F39B7FFF", "#91D1C2FF"), 
               line.color = "gray", line.size = 0.4,
               short.panel.labs = FALSE,width = 0.4)+
  facet_grid(~method)+
  stat_compare_means(paired = T,label.y=1.05)+
  # theme_bw()+
  theme(legend.position = 'top')+
  labs(y='AUC',x='')+
  theme(text = element_text(family = '',size = 12,face = 'plain'))

print(p)

pdf("./02_figure/FigS8A.pdf", height = 3.5, width = 2.1)
p
dev.off()




## FS8B ----
project_smote_last <- get(load('01_data/project_smote.RData'))
cross.siamcat.all <- get(load('01_data/SMOTE_AUC.RData'))
cross.siamcat.all$method[cross.siamcat.all$method=='self'] <- 'Internal'
cross.siamcat.all$method[cross.siamcat.all$method=='external'] <- 'External'
cross.siamcat.all$method <- factor(cross.siamcat.all$method,levels = c('Internal','External'))
cross.siamcat.all1 <- subset(cross.siamcat.all,(TE %in% project_smote_last$project))
save(cross.siamcat.all1,file='01_data/plot_data/FS8B.RData')

p <- ggboxplot(cross.siamcat.all1, x = 'Label', y = 'auc', fill = 'Label',palette = c("#FFB9007F","#5773CC7F"), 
               line.color = "gray", line.size = 0.4,
               short.panel.labs = FALSE,width = 0.4)+
  facet_grid(~method)+
  stat_compare_means(paired = T,label.y=1.05)+
  # theme_bw()+
  theme(legend.position = 'top')+
  labs(y='AUC',x='')+
  theme(text = element_text(family = '',size = 12,face = 'plain'))

print(p)

pdf("./02_figure/FigS8B.pdf", height = 3.5, width = 2.1)
p
dev.off()









