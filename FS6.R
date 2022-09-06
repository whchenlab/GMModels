### Figure S6 -----
auc_marker_lodo <- get(load('01_data/AUC_MSI_lodo_lasso.RData'))
auc_marker_lodo$level[auc_marker_lodo$level %in% 'Metagenomics_genus'] <- 'mNGS_genus'
auc_marker_lodo$level[auc_marker_lodo$level %in% 'Metagenomics_species'] <- 'mNGS_species'
auc_marker <- auc_marker_lodo

save(auc_marker,file = '01_data/plot_data/FS6.RData')


{
  cor <- corr.test(auc_marker[,1:2],method = 'spearman')
  cor$p[1,2]
  cor$r[1,2]
}
p1=ggscatter(data=auc_marker,x='dis',y='auc',shape='level',color='disease',palette = "igv",alpha=0.8,size=2)+
  geom_smooth(mapping=aes(x=dis,y=auc),data=auc_marker,formula =y~x, method=lm,inherit.aes = F,size=0.5)+
  # stat_poly_eq(data=auc_marker,mapping=aes(x=dis,y=auc,label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,inherit.aes = F,size=3)+
  geom_text_repel(data=auc_marker,aes(x=dis,y=auc,label = disease),inherit.aes = F,size=3)+
  xlab("MSI (median)") + ylab("External AUC (median)")+
  labs(color = "disease")+
  annotate('text',x=1.8,y=0.84,label=paste0('r=',round(cor$r[1,2],2),'  ','p','=',format(cor$p[1,2],scientific=T,digits=3)),size=3)+
  theme(legend.position = 'none')+
  border()
p2 <- ggdensity(auc_marker, 'dis', fill="disease_type2",palette = "npg")+
  clean_theme()+
  theme(legend.position = 'none')
wilcox.test(x=subset(auc_marker,disease_type2=='Intestinal')$dis,y=subset(auc_marker,disease_type2!='Intestinal')$dis)

p3 <- ggdensity(auc_marker, 'auc', fill="disease_type2", palette = "npg")+
  clean_theme()+
  ggpubr::rotate()+
  theme(legend.position = 'none')
wilcox.test(x=subset(auc_marker,disease_type2=='Intestinal')$auc,y=subset(auc_marker,disease_type2!='Intestinal')$auc)

p=ggarrange(p2, NULL, p1, p3, ncol = 2, nrow = 2, align = "hv", widths = c(3, 1), 
            heights = c(1, 3))

print(p)


pdf("./02_figure/FigS6.pdf", height = 3.8, width = 5)
p
dev.off()

