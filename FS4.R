## FS4A ----
auc_external_mean_e <- get(load('01_data/AUC_ext_median.RData'))
auc_external_mean_e_i <- subset(auc_external_mean_e,group1=='Intestinal')
save(auc_external_mean_e_i,file='01_data/plot_data/FS4A.RData')

auc_external_mean_e1 <- subset(auc_external_mean_e_i,level=='Amplicon_genus')
auc_external_mean_e2 <- subset(auc_external_mean_e_i,level=='Metagenomics_genus')
auc_external_mean_e3 <- subset(auc_external_mean_e_i,level=='Metagenomics_species')

p1 <- ggboxplot(auc_external_mean_e1, x = "disease", color="method",y = "auc", alpha=0.3,
                palette = "npg",add = "jitter",width = 0.5)+ 
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  stat_compare_means(aes(group=method),label = "p.signif",label.y = 0.95,paired = T)+
  xlab("Amplicon_genus") + ylab("External AUC")+
  labs(fill ="disease type")+
  ylim(0.4,1)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")  
p2 <- ggboxplot(auc_external_mean_e2, x = "disease", color="method",y = "auc", alpha=0.3,
                palette = "npg",add = "jitter",width = 0.5)+ 
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  stat_compare_means(aes(group=method),label = "p.signif",label.y = 0.95,paired = T)+
  xlab("mNGS_genus") + ylab("")+
  labs(fill ="disease type")+
  ylim(0.4,1)+
  geom_hline(aes(yintercept=0.7),color='#8DD3C7')+
  geom_hline(aes(yintercept=0.6),color='orange')+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")  
p3 <- ggboxplot(auc_external_mean_e3, x = "disease", color="method",y = "auc", alpha=0.3,
                palette = "npg",add = "jitter",width = 0.5)+ 
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  # geom_hline(yintercept =0.9,color='#e60020')+
  stat_compare_means(aes(group=method),label = "p.signif",label.y = 0.95,paired = T)+
  xlab("mNGS_species") + ylab("")+
  ylim(0.4,1)+
  geom_hline(aes(yintercept=0.7),color='#8DD3C7')+
  geom_hline(aes(yintercept=0.6),color='orange')+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")  

p=ggarrange(p1,p2,p3,
            ncol = 3, nrow = 1, 
            widths = c(4,5,5)) 
print(p)

pdf("./02_figure/FigS4A.pdf", height = 2.8, width = 6.72)
p
dev.off()


## FS4B ----
load('./01_data/p_r_CRC_species.RData')
p_r_test <- subset(p_r,method=='test')

save(p_r_test,file='01_data/plot_data/FS4B_left.RData')

for (i in unique(p_r_test$num)){
  AUC_median <- median(subset(p_r_test,num==i)$AUC)
  if (i==1){
    p_r_test_median <- data.frame(AUC=AUC_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(AUC=AUC_median,method='test',num=i))
  } 
}
p_r_test_median
p3<-ggboxplot(p_r_test,x='num',y='AUC',color = pal_npg('nrc')(10)[5],add = "jitter",width = 0.2,add.params = list(size = 1))+
  xlab('No. of training cohorts')+
  labs(title = 'CRC')+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  geom_hline(yintercept =0.9,color='#e60020')+
  geom_line(data = p_r_test_median,aes(x=num,y=AUC),color='#8DD3C7')+
  geom_point(data = p_r_test_median,aes(x=num,y=AUC),color='#8DD3C7')+
  ylab("External AUC")+
  geom_hline(yintercept =0.7,color=pal_d3('category10')(10)[1])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylim(0.5,1)
p3

load('./data/p_r_CD_species.RData')
p_r_test <- subset(p_r,method=='test')

save(p_r_test,file='01_data/plot_data/FS4B_right.RData')
for (i in unique(p_r_test$num)){
  AUC_median <- median(subset(p_r_test,num==i)$AUC)
  if (i==1){
    p_r_test_median <- data.frame(AUC=AUC_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(AUC=AUC_median,method='test',num=i))
  } 
}
p_r_test_median

p4<-ggboxplot(p_r_test,x='num',y='AUC',color = pal_npg('nrc')(10)[5],add = "jitter",width = 0.2,add.params = list(size = 1))+
  xlab('No. of training cohorts')+
  labs(title = 'CD')+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  geom_hline(yintercept =0.9,color='#e60020')+
  geom_line(data = p_r_test_median,aes(x=num,y=AUC),color='#8DD3C7')+
  geom_point(data = p_r_test_median,aes(x=num,y=AUC),color='#8DD3C7')+
  ylab("External AUC")+
  geom_hline(yintercept =0.7,color=pal_d3('category10')(10)[1])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylim(0.5,1)
p4


p=ggarrange(p3,p4,
            ncol = 2, nrow = 1, 
            widths = c(6,5))

p
pdf("./02_figure/FigS4B.pdf", height = 2.8, width = 6.72)
p
dev.off()
## FS4C ----
#CRC
p_r0 <- get(load('./data/p_r_CRC_lasso_species_number_1_1.RData'))
p_r0 <- subset(p_r0,num!=0)
p_r0 <- p_r0[sort(p_r0$num,index.return =T)$ix,]


p_r <- p_r0
table(p_r$num)
s <- unique(p_r$num)
s<-s[!table(p_r$num)<max(table(p_r$num))]
p_r <- subset(p_r,num%in%s)

save(p_r,file = '01_data/plot_data/FS4C_left.RData')

l_max <- max(p_r$num)
for (i in sort(unique(p_r$num))){
  AUC_median <- median(subset(p_r,num==i)$AUC.test)
  if (i==sort(unique(p_r$num))[1]){
    p_r_test_median <- data.frame(AUC.test=AUC_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(AUC.test=AUC_median,method='test',num=i))
  }
}
class(p_r_test_median$num)


require(psych)
{
  cor <- corr.test(p_r_test_median[,c('num','AUC.test')],method = 'spearman')
  cor$p[1,2]
  cor$r[1,2]
}

label_x <- seq(100,l_max,100)
label_x <- c(16,40,label_x)

p1 <- ggplot(p_r,aes(x=num,y=AUC.test,group=num))+
  geom_boxplot(width = 2.5,color =pal_npg('nrc')(10)[6],fill =pal_npg('nrc')(10)[6])+
  xlab('No. of training samples')+
  ylab('External AUC')+
  theme_classic()+
  labs(title = 'CRC')+
  geom_line(data = p_r_test_median,aes(x=num,y=AUC.test),color=pal_npg('nrc')(10)[7],inherit.aes = F)+
  geom_point(data = p_r_test_median,aes(x=num,y=AUC.test),color=pal_npg('nrc')(10)[7],inherit.aes = F)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  geom_hline(yintercept =0.9,color='#e60020')+
  ylim(0.35,1)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+
  geom_smooth(data=p_r_test_median,aes(x=num, y=AUC.test), formula =y~x,se =T,method = 'lm',size=1,colour=pal_npg('nrc')(10)[8],inherit.aes = F)+
  annotate('text',x=200,y=0.95,label=paste0('r=',round(cor$r[1,2],2),'  ','p','=',format(cor$p[1,2],scientific=T,digits=3)),size=3)+
  
  # stat_poly_eq(data=p_r_test_median,mapping=aes(x=num,y=AUC.test,label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),method = 'lm',formula = y~x,parse=T,inherit.aes = F,size=3)+
  scale_x_discrete(limits = label_x)
p1



#CD
p_r0 <- get(load('./data/p_r_CD_lasso_species_number_1_1.RData'))
p_r0 <- get(load( 'data/p_r_CD_lasso_species_number_1_1_all.RData'))
p_r0 <- subset(p_r0,num!=0)
p_r0 <- p_r0[sort(p_r0$num,index.return =T)$ix,]


p_r <- p_r0
s <- unique(p_r$num)
table(p_r$num)
s<-s[-(19:22)]

p_r <- subset(p_r,num%in%s)

save(p_r,file = '01_data/plot_data/FS4C_right.RData')

l_max <- max(p_r$num)


for (i in sort(unique(p_r$num))){
  AUC_median <- median(subset(p_r,num==i)$AUC.test)
  if (i==sort(unique(p_r$num))[1]){
    p_r_test_median <- data.frame(AUC.test=AUC_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(AUC.test=AUC_median,method='test',num=i))
  }
}
# p_r_test_median$num <- 1:length(p_r_test_median$num)
class(p_r_test_median$num)
require(psych)
{
  cor <- corr.test(p_r_test_median[,c('num','AUC.test')],method = 'spearman')
  cor$p[1,2]
  cor$r[1,2]
}

label_x <- seq(100,l_max,100)
label_x <- c(16,40,label_x)
p2 <- ggplot(p_r,aes(x=num,y=AUC.test,group=num))+
  geom_boxplot(width = 2.5,color =pal_npg('nrc')(10)[6],fill =pal_npg('nrc')(10)[6])+
  xlab('No. of training samples')+
  ylab('External AUC')+
  theme_classic()+
  labs(title = 'CD')+
  geom_line(data = p_r_test_median,aes(x=num,y=AUC.test),color=pal_npg('nrc')(10)[7],inherit.aes = F)+
  geom_point(data = p_r_test_median,aes(x=num,y=AUC.test),color=pal_npg('nrc')(10)[7],inherit.aes = F)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  geom_hline(yintercept =0.9,color='#e60020')+
  ylim(0.35,1)+
  
  geom_smooth(data=p_r_test_median,aes(x=num, y=AUC.test), formula =y~x,se =T,method = 'lm',size=1,colour=pal_npg('nrc')(10)[8],inherit.aes = F)+
  annotate('text',x=150,y=0.95,label=paste0('r=',round(cor$r[1,2],2),'  ','p','=',format(cor$p[1,2],scientific=T,digits=3)),size=3)+
  
  # stat_poly_eq(data=p_r_test_median,mapping=aes(x=num,y=AUC.test,label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),method = 'lm',formula = y~x,parse=T,inherit.aes = F,size=3)+
  scale_x_discrete(limits = label_x)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))
p2

p <- ggarrange(p1,p2,nrow = 1,ncol = 2,widths = c(31,21))
print(p)
pdf("./02_figure/FigS4C.pdf", height = 3.4, width = 8.16)
p
dev.off()

