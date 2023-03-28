### Figure S6 -----

## FS6 A ----
auc_marker_lodo <- get(load('01_data/AUC_MSI_lodo_lasso.RData'))
auc_marker_lodo$level[auc_marker_lodo$level %in% 'Metagenomics_genus'] <- 'mNGS_genus'
auc_marker_lodo$level[auc_marker_lodo$level %in% 'Metagenomics_species'] <- 'mNGS_species'
auc_marker <- auc_marker_lodo

save(auc_marker,file = '01_data/plot_data/FS6A.RData')


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


pdf("./02_figure/FigS6A.pdf", height = 3.8, width = 5)
p
dev.off()

## FS6 B ----

self.e <- get(load('./01_data/MSI_all.RData'))

self.e$data_type[self.e$data_type=='Amplicon'] <- '16S'
self.e$level <- paste0(self.e$data_type,'_',self.e$taxon)

self.e <- subset(self.e,level=='16S_genus'&disease!='IBS')

save(self.e,file = '01_data/plot_data/FS6B.RData')

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



pdf("./02_figure/FigS6B.pdf", height = 5, width = 3)
p
dev.off()

## FS6 C ----
load('01_data/p_r_CRC_marker_new_species.RData')
p_r_test <- subset(p_r,method=='test')

save(p_r_test,file = '01_data/plot_data/FS6C_left.RData')
for (i in unique(p_r_test$num)){
  dis_median <- median(subset(p_r_test,num==i)$dis)
  if (i==1){
    p_r_test_median <- data.frame(dis=dis_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(dis=dis_median,method='test',num=i))
  } 
}
p3<-ggboxplot( p_r_test,x='num',y='dis',color=pal_npg("nrc")(10)[1],add = "jitter",width = 0.3)+
  xlab('No. of training cohorts')+
  ylab('MSI (median)')+
  labs(title = 'CRC')+
  geom_line(data = p_r_test_median,aes(x=num,y=dis),color=pal_npg("nrc")(10)[2])+
  geom_point(data = p_r_test_median,aes(x=num,y=dis),color=pal_npg("nrc")(10)[2])+
  # geom_hline(yintercept =0.5,color=pal_d3('category10')(10)[4])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))

p3


load('data/p_r_CD_marker_new_species.RData')
p_r_test <- subset(p_r,method=='test')

save(p_r_test,file = '01_data/plot_data/FS6C_right.RData')

for (i in unique(p_r_test$num)){
  dis_median <- median(subset(p_r_test,num==i)$dis)
  if (i==1){
    p_r_test_median <- data.frame(dis=dis_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(dis=dis_median,method='test',num=i))
  } 
}
p4<-ggboxplot( p_r_test,x='num',y='dis',color=pal_npg("nrc")(10)[1],add = "jitter",width = 0.3)+
  xlab('No. of training cohorts')+
  ylab('MSI (median)')+
  labs(title = 'CD')+
  geom_line(data = p_r_test_median,aes(x=num,y=dis),color=pal_npg("nrc")(10)[2])+
  geom_point(data = p_r_test_median,aes(x=num,y=dis),color=pal_npg("nrc")(10)[2])+
  # geom_hline(yintercept =0.5,color=pal_d3('category10')(10)[4])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))

p=p3+p4
p
pdf("./02_figure/FigS6C.pdf", height = 3.8, width = 8)
p
dev.off()




