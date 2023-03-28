### Figure S1

## FS3 A ----
a=get(load('./01_data/AUC_all_lasso.RData'))
a$level[a$level=='Metagenomics_genus'] <- 'mNGS_genus'
a$level[a$level=='Metagenomics_species'] <- 'mNGS_species'
a$method[a$method == 'self'] <- 'internal'
a_d <- subset(a,group1=='Intestinal')
a_d$disease <- factor(a_d$disease,levels = c('CD','IBD','CRC','CDI','UC','Adenoma','IBS'))
save(a_d,file='01_data/plot_data/FS3A.RData')
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

pdf("./02_figure/FigS3A.pdf", height = 3.4, width = 6)
p
dev.off()

## FS3 B ----
file <- '04_model_raw_RData/single_LODO/02lasso/'
model.adj2 <- get(load(paste0(file,'D015179_Metagenomics.species_7.RData')))
FS1C <- model.adj2[["result"]]
save(FS1C,file='01_data/plot_data/FS3B.RData')
p2 <- my_auc.heatmap(datas=model.adj2$result,models='lasso_adj_nc')






