### Figure 4


## F34 ----
auc_external_mean_e <- get(load('01_data/AUC_ext_median.RData'))
auc_external_mean_e$group2 <- factor(auc_external_mean_e$group2,levels = c('non-Intestinal','Intestinal'))
save(auc_external_mean_e,file='01_data/plot_data/F4A.RData')

sig <- function(x){ifelse(x<0.05, ifelse(x<0.01, ifelse(x<0.001, ifelse(x<=0.0001, '****','***'),'**'),'*'),'ns')}
p1 <- ggviolin(subset(auc_external_mean_e,group2=='non-Intestinal'), x = "method", fill="method",y = "auc", alpha=0.3,
               palette = "npg",add = "boxplot",width = 0.4)+ 
  # facet_wrap(~group2)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  stat_compare_means(aes(label = ..p.signif..),comparisons=list(c('Single cohort','LODO')),paired = T,label.y = 1.05)+
  xlab("non-Intestinal") + ylab("External AUC (median)")+
  # labs(fill ="Method")+
  ylim(0.1,1.12)+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")  
p1

p2 <- ggviolin(subset(auc_external_mean_e,group2=='Intestinal'), x = "method", fill="method",y = "auc", alpha=0.3,
              palette = "npg",add = "boxplot",width = 0.4)+ 
  # facet_wrap(~group2)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  geom_hline(yintercept =0.8,color='#7b77ff')+
  stat_compare_means(aes(label = ..p.signif..),comparisons=list(c('Single cohort','LODO')),paired = T,label.y = 1.05)+
  xlab("Intestinal") + ylab("External AUC (median)")+
  # labs(fill ="Method")+
  ylim(0.1,1.12)+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")  
p2

p=ggarrange(p1,p2,
            ncol = 2, nrow = 1)
print(p)


pdf("./02_figure/Fig4A.pdf", height = 2.8, width = 5.5)
p
dev.off()


## F3B ----
auc_external_mean_e <- get(load('01_data/AUC_ext_median.RData'))
auc_external_mean_e_i <- subset(auc_external_mean_e,group1!='Intestinal')

save(auc_external_mean_e_i,file='01_data/plot_data/F4B.RData')

auc_external_mean_e1 <- subset(auc_external_mean_e_i,level=='Amplicon_genus')
auc_external_mean_e2 <- subset(auc_external_mean_e_i,level=='Metagenomics_genus')
auc_external_mean_e3 <- subset(auc_external_mean_e_i,level=='Metagenomics_species')

p1 <- ggboxplot(auc_external_mean_e1, x = "disease", color="method",y = "auc", alpha=0.3,
                palette = "npg",add = "jitter",width = 0.5)+ 
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  stat_compare_means(aes(group=method),label = "p.signif",label.y = 0.85,paired = T)+
  xlab("16S_genus") + ylab("External AUC")+
  labs(fill ="Disease type")+
  ylim(0.30,0.87)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")  #family = 字体 #p$layers[[1]]$aes_params$textsize <-13
p2 <- ggboxplot(auc_external_mean_e2, x = "disease", color="method",y = "auc", alpha=0.3,
                palette = "npg",add = "jitter",width = 0.5)+ 
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  stat_compare_means(aes(group=method),label = "p.signif",label.y = 0.85,paired = T)+
  xlab("mNGS_genus") + ylab("")+
  labs(fill ="disease type")+
  ylim(0.30,0.87)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")  #family = 字体 #p$layers[[1]]$aes_params$textsize <-13
p3 <- ggboxplot(auc_external_mean_e3, x = "disease", color="method",y = "auc", alpha=0.3,
                palette = "npg",add = "jitter",width = 0.5)+ 
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  stat_compare_means(aes(group=method),label = "p.signif",label.y = 0.85,paired = T)+
  xlab("mNGS_species") + ylab("")+
  ylim(0.30,0.87)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'),legend.position="none")  #family = 字体 #p$layers[[1]]$aes_params$textsize <-13


p=ggarrange(p1,p2,p3,
            ncol = 3, nrow = 1, 
            widths = c(2.2,1.15,1.15))
print(p)

pdf("./02_figure/Fig4B.pdf", height = 2.8, width = 7)
p
dev.off()

## F3C -----
#CCM ASD PD
load('./01_data/p_r_ASD.RData')
p_r_test <- subset(p_r,method=='test')

save(p_r_test,file = '01_data/plot_data/F4C_left.RData')


for (i in unique(p_r_test$num)){
  AUC_median <- median(subset(p_r_test,num==i)$AUC)
  if (i==1){
    p_r_test_median <- data.frame(AUC=AUC_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(AUC=AUC_median,method='test',num=i))
  }
}

p1<-ggboxplot(p_r_test,x='num',y='AUC',color = pal_npg('nrc')(10)[5],add = "jitter",width = 0.3)+
  xlab('No. of training dataset')+
  labs(title = 'ASD')+
  geom_line(data = p_r_test_median,aes(x=num,y=AUC),color='#8DD3C7')+
  geom_point(data = p_r_test_median,aes(x=num,y=AUC),color='#8DD3C7')+
  geom_hline(yintercept =0.5,color=pal_d3('category10')(10)[1])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylim(0.30,0.95)

load('./01_data/p_r_PD.RData')
p_r_test <- subset(p_r,method=='test')

save(p_r_test,file = '01_data/plot_data/F4C_right.RData')

for (i in unique(p_r_test$num)){
  AUC_median <- median(subset(p_r_test,num==i)$AUC)
  if (i==1){
    p_r_test_median <- data.frame(AUC=AUC_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(AUC=AUC_median,method='test',num=i))
  } 
}


p2<-ggboxplot(p_r_test,x='num',y='AUC',color = pal_npg('nrc')(10)[5],add = "jitter",width = 0.3)+
  xlab('No. of training dataset')+
  labs(title = 'PD')+
  geom_line(data = p_r_test_median,aes(x=num,y=AUC),color='#8DD3C7')+
  geom_point(data = p_r_test_median,aes(x=num,y=AUC),color='#8DD3C7')+
  geom_hline(yintercept =0.5,color=pal_d3('category10')(10)[1])+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+  
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylim(0.30,0.95)

p=p1+p2
p
pdf("./02_figure/Fig4C.pdf", height = 3.3, width = 7)
p
dev.off()


## F4D -----
#ASD
p_r0 <- get(load('./01_data/p_r_ASD_lasso_genus_number_1_1.RData'))
p_r0 <- subset(p_r0,num!=0)
p_r0 <- p_r0[sort(p_r0$num,index.return =T)$ix,]

p_r <- p_r0
s <- unique(p_r$num)
s<-s[!table(p_r$num)<max(table(p_r$num))]


p_r <- subset(p_r,num%in%s)
save(p_r,file = '01_data/plot_data/F4D_left.RData')

l_max <- max(p_r$num)
for (i in sort(unique(p_r$num))){
  AUC_median <- median(subset(p_r,num==i)$AUC.test)
  if (i==sort(unique(p_r$num))[1]){
    p_r_test_median <- data.frame(AUC.test=AUC_median,method='test',num=i)
  }else{
    p_r_test_median <- rbind(p_r_test_median,data.frame(AUC.test=AUC_median,method='test',num=i))
  }
}

require(psych)
{
  cor <- corr.test(p_r_test_median[,c('num','AUC.test')],method = 'spearman')
  cor$p[1,2]
  cor$r[1,2]
}
label_x <- seq(100,l_max,100)
label_x <- c(16,40,label_x)
library(ggpmisc)
p1 <- ggplot(p_r,aes(x=num,y=AUC.test,group=num))+
  geom_boxplot(width = 2.5,color =pal_npg('nrc')(10)[6],fill =pal_npg('nrc')(10)[6])+
  xlab('No. of training samples')+
  ylab('External AUC')+
  theme_classic()+
  labs(title = 'ASD')+
  geom_line(data = p_r_test_median,aes(x=num,y=AUC.test),color=pal_npg('nrc')(10)[7],inherit.aes = F)+
  geom_point(data = p_r_test_median,aes(x=num,y=AUC.test),color=pal_npg('nrc')(10)[7],inherit.aes = F)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  ylim(0.1,1)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+
  geom_smooth(data=p_r_test_median,aes(x=num, y=AUC.test), formula =y~x,se =T,method = 'lm',size=1,colour=pal_npg('nrc')(10)[8],inherit.aes = F)+
  annotate('text',x=150,y=0.85,label=paste0('r=',round(cor$r[1,2],2),'  ','p','=',format(cor$p[1,2],scientific=T,digits=3)),size=3)+
  
  stat_poly_eq(data=p_r_test_median,mapping=aes(x=num,y=AUC.test,label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),method = 'lm',formula = y~x,parse=T,inherit.aes = F,size=3)+
  scale_x_discrete(limits=label_x)
p1

#calculate 95CI
fit <- lm(AUC.test ~ num, p_r_test_median)
summary(fit)
confint(fit,level=0.95)
a <- data.frame(num = 1600)
result <-  predict(fit,a,interval = "prediction",level = 0.95)
print(result)
print(round(result,2))


#PD
p_r0 <- get(load('./01_data/p_r_PD_lasso_genus_number_1_1.RData'))
p_r0 <- subset(p_r0,num!=0)
p_r0 <- p_r0[sort(p_r0$num,index.return =T)$ix,]

p_r <- p_r0
s <- unique(p_r$num)
s<-s[!table(p_r$num)<max(table(p_r$num))]

p_r <- subset(p_r,num%in%s)

save(p_r,file = '01_data/plot_data/F4D_right.RData')

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
  labs(title = 'PD')+
  # theme_void()+
  # theme_linedraw()+
  geom_line(data = p_r_test_median,aes(x=num,y=AUC.test),color=pal_npg('nrc')(10)[7],inherit.aes = F)+
  geom_point(data = p_r_test_median,aes(x=num,y=AUC.test),color=pal_npg('nrc')(10)[7],inherit.aes = F)+
  geom_hline(yintercept =0.5,color='#dbdcdc')+
  geom_hline(yintercept =0.6,color='#ffd09a')+
  geom_hline(yintercept =0.7,color='#ffcbd8')+
  ylim(0.1,1)+
  theme(text = element_text(size=13,face = 'plain',family ='',colour = 'black'))+
  geom_smooth(data=p_r_test_median,aes(x=num, y=AUC.test), formula =y~x,se =T,method = 'lm',size=1,colour=pal_npg('nrc')(10)[8],inherit.aes = F)+
  annotate('text',x=150,y=0.85,label=paste0('r=',round(cor$r[1,2],2),'  ','p','=',format(cor$p[1,2],scientific=T,digits=3)),size=3)+
  
  stat_poly_eq(data=p_r_test_median,mapping=aes(x=num,y=AUC.test,label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),method = 'lm',formula = y~x,parse=T,inherit.aes = F,size=3)+
  scale_x_discrete(limits=label_x)
p2

#calculate 95CI
fit <- lm(AUC.test ~ num, p_r_test_median)
summary(fit)
confint(fit,level=0.95)
a <- data.frame(num = 2400)
result <-  predict(fit,a,interval = "prediction",level = 0.95)
print(result)
print(round(result,2))

p <- ggarrange(p1,p2,nrow = 1,ncol = 2,widths = c(21,45))
print(p)
pdf("./02_figure/Fig4D.pdf", height = 3.3, width = 7)
p
dev.off()

