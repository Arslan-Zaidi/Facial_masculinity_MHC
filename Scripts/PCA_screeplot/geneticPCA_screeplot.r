library(ggplot2)
library(data.table)

evec<-fread("../../Results/Genetic_pca/ADAPT_1000G_NoPalin_LDPrune_Euro_NoFins_ADAPTRef.pca.evec",header=T)
eval<-fread("../../Results/Genetic_pca/ADAPT_1000G_NoPalin_LDPrune_Euro_NoFins_ADAPTRef.pca.eval")
colnames(eval)<-c("Eigenvalue")
eval$gPC<-seq(1,nrow(eval),1)

scree<-ggplot(data=eval[which(eval$gPC<11),],aes(gPC,Eigenvalue))+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_x_continuous(breaks=seq(1,10))+
  ylim(c(1,4))

ggsave("../../Results/Genetic_pca/gPC_screeplot_1_10_09102018.pdf",scree,height=5,width=7)


