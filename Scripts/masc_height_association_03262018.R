#Author: AAZaidi
#Tests for association between facial masculinity and height
library(data.table)
library(ggplot2)

#read masculinity calculated per QL - unscaled
qlmasc.u<-as.matrix(fread('../Dataset/qlmasc_1233_noprop_03292018.txt',header=F,sep="\t"))

#read file containing data for sex, height, weight, etc.
eurofam.het<-read.table('../Dataset/euro_1233_masc_het_03292018.dat',header=T,sep="\t",stringsAsFactors = F)

#plot relationship between overall masculinity and height
masc.height.scatter<-ggplot(data=eurofam.het,aes(Height,avg.masc.unit,color=Sex))+
  geom_point(alpha=0.7)+
  stat_smooth(method="lm",se=T)+
  theme_bw()+
  scale_color_manual(values=c("#ff7f00","#377eb8"))+
  labs(x="Height in cm",y="Average facial masculinity")+
  theme(legend.position="none")

ggsave('../Results/Masc_v_height/masc_height_scatter_03292018.pdf',masc.height.scatter,height=6,width=4)

masc.density<-ggplot(data=eurofam.het,aes(avg.masc.unit,fill=Sex))+
  geom_density(alpha=0.6)+
  theme_bw()+
  scale_fill_manual(values=c("#ff7f00","#377eb8"))+
  coord_flip()+
  theme(axis.title.y=element_blank())

ggsave('../Results/masc_density_03292018.pdf',masc.density,height=7,width=4)


#center and scale data to generate standardized regression coefficients
eurofam.het[,c('avg.masc','Age','Height','Weight','phet_genome','phet_hla','gPC1','gPC2','gPC3','gPC4')]<-apply(eurofam.het[,c('avg.masc','Age','Height','Weight','phet_genome','phet_hla','gPC1','gPC2','gPC3','gPC4')],2,scale)

#testing relationship between overall masculinity and predictors
lm.masc.height.all<-lm(data=eurofam.het,avg.masc~Sex+Height+Age+Weight+gPC1+gPC2+gPC3)
capture.output(summary(lm.masc.height.all)$coefficients,file="../Results/Masc_v_height/lm_overallmasc_height_results_03292018.txt")

#separately within each sex
lm.masc.height.male<-lm(data=eurofam.het[which(eurofam.het$Sex=="Male"),], avg.masc~Height+Age+Weight+gPC1+gPC2+gPC3)
capture.output(summary(lm.masc.height.male)$coefficients,file="../Results/Masc_v_height/lm_overallmasc_height_results_male_03292018.txt")

lm.masc.height.female<-lm(data=eurofam.het[which(eurofam.het$Sex=="Female"),], avg.masc~Height+Age+Weight+gPC1+gPC2+gPC3)
capture.output(summary(lm.masc.height.female)$coefficients,file="../Results/Masc_v_height/lm_overallmasc_height_results_female_03292018.txt",fill=T)

#test whether the two slopes are different from each other
beta.m_f<-summary(lm.masc.height.male)$coefficients[2,1]-summary(lm.masc.height.female)$coefficients[2,1]
se.beta.m_f<-sqrt((summary(lm.masc.height.male)$coefficients[2,2]^2)+(summary(lm.masc.height.female)$coefficients[2,2]^2))
z.m_f<-beta.m_f/se.beta.m_f
p.m_f<-pnorm(abs(z.m_f),lower.tail = F)

qlmasc.u<-apply(qlmasc.u,2,scale)

#testing masculinity by height relationship
masc.height.b<-matrix(NA,7150,8)
masc.height.p<-matrix(NA,7150,8)
pb<-txtProgressBar(min=0,max=7150,style=3)
for(i in 1:7150){
  l1m<-lm(data=eurofam.het,qlmasc.u[,i]~Height+Sex+Age+Weight+gPC1+gPC2+gPC3)
  s1m<-summary(l1m)
  masc.height.b[i,]<-s1m$coefficients[c(1:8),1]
  masc.height.p[i,]<-s1m$coefficients[c(1:8),4]
  setTxtProgressBar(pb,i)
}

masc.height.q<-masc.height.p
masc.height.q[which(masc.height.p<0.05/7150)]<-1
masc.height.q[which(masc.height.p>=0.05/7150)]<-0

colnames(masc.height.b)<-colnames(masc.height.p)<-colnames(masc.height.q)<-c("Intercept","Height","Sex","Age","Weight","gPC1","gPC2","gPC3")
#write to file
write.table(masc.height.b,"../Results/Masc_v_height/qlmasc_height_b_noprop_sd_04142018.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(masc.height.p,"../Results/Masc_v_height/qlmasc_height_p_noprop_sd_0414018.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(masc.height.q,"../Results/Masc_v_height/qlmasc_height_q_noprop_sd_04142018.txt",sep="\t",col.names=T,row.names=F,quote=F)

#testing different slopes b/w masc and height in men and women
#separate ql masculinity table by sex for use later
qlmasc.m<-qlmasc.u[eurofam.het$Sex=="Male",]
qlmasc.f<-qlmasc.u[eurofam.het$Sex=="Female",]

qlmasc.m<-apply(qlmasc.m,2,scale)
qlmasc.f<-apply(qlmasc.f,2,scale)


masc.h.z<-matrix(NA,nrow=7150,ncol=5)
colnames(masc.h.z)<-c("slope_female","slope_male","diff","z","p.value")
pb<-txtProgressBar(min=0,max=7150,style=3)
for(i in 1:7150){
  datm<-cbind(qlmasc.m[,i],eurofam.het[eurofam.het$Sex=="Male",])
  colnames(datm)[1]<-"ql"
  datf<-cbind(qlmasc.f[,i],eurofam.het[eurofam.het$Sex=="Female",])
  colnames(datf)[1]<-"ql"
  #residualize
  l1m<-lm(data=datm,ql~Height+Age+Weight+gPC1+gPC2+gPC3)
  l1f<-lm(data=datf,ql~Height+Age+Weight+gPC1+gPC2+gPC3)
  
  b12<-summary(l1m)$coefficients[2,1]-summary(l1f)$coefficients[2,1]
  sb12<-sqrt((summary(l1m)$coefficients[2,2])+(summary(l1f)$coefficients[2,2]^2))
  z<-b12/sb12
  
  masc.h.z[i,1]<-summary(l1f)$coefficients[2,1]
  masc.h.z[i,2]<-summary(l1m)$coefficients[2,1]
  masc.h.z[i,3]<-b12
  masc.h.z[i,4]<-z
  masc.h.z[i,5]<-pnorm(z)
  setTxtProgressBar(pb,i)
}

masc.h.z<-as.data.frame(masc.h.z)
masc.h.z$q.value<-NA
masc.h.z$q.value[which(masc.h.z$p.value<(0.05/7150))]<-1
masc.h.z$q.value[which(masc.h.z$p.value>=(0.05/7150))]<-0

write.table(masc.h.z,'../Results/Masc_v_height/qlmasc_h_bysex_noprop_sd_04142018.txt',sep="\t",col.names=T,row.names=F,quote=F)
