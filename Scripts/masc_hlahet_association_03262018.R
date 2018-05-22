#Author: AAZaidi
#tests for association b/w facial masculinity and MHC heterozygosity
library(data.table)
library(ggplot2)

#read masculinity calculated per QL
qlmasc.u<-as.matrix(fread('../Dataset/qlmasc_1233_noprop_03292018.txt',header=F,sep="\t"))

#read file containing data for sex, height, weight, etc. 
eurofam.het<-read.table('../Dataset/euro_1233_masc_het_03292018.dat',header=T,sep="\t",stringsAsFactors = F)

#separate ql masculinity table by sex for use later
qlmasc.m<-qlmasc.u[eurofam.het$Sex=="Male",]
qlmasc.f<-qlmasc.u[eurofam.het$Sex=="Female",]

#center and scale data to generate standardized regression coefficients
eurofam.het[,c('avg.masc','Age','Height','Weight','phet_genome','phet_hla','gPC1','gPC2','gPC3','gPC4')]<-apply(eurofam.het[,c('avg.masc','Age','Height','Weight','phet_genome','phet_hla','gPC1','gPC2','gPC3','gPC4')],2,scale)

#testing relationship between overall masculinity and MHC heterozygosity without correcting for height
lm.masc.mhc.noh.all<-lm(data=eurofam.het,avg.masc~phet_hla+phet_genome+Sex+Age+Weight+gPC1+gPC2+gPC3)
capture.output(summary(lm.masc.mhc.noh.all)$coefficients,file="../Results/Masc_v_mhc/lm_overallmasc_mhc_noh_03292018.txt")

#testing relationship between overall masculinity and MHC heterozygosity with height correction
lm.masc.mhc.h.all<-lm(data=eurofam.het,avg.masc~phet_hla+phet_genome+Height+Sex+Age+Weight+gPC1+gPC2+gPC3)
capture.output(summary(lm.masc.mhc.h.all)$coefficients,file="../Results/Masc_v_mhc/lm_overallmasc_mhc_h_03292018.txt")

#separately within each sex
lm.masc.mhc.h.male<-lm(data=eurofam.het[which(eurofam.het$Sex=="Male"),], avg.masc~phet_hla+phet_genome+Height+Age+Weight+gPC1+gPC2+gPC3)
lm.masc.mhc.h.female<-lm(data=eurofam.het[which(eurofam.het$Sex=="Female"),], avg.masc~phet_hla+phet_genome+Height+Age+Weight+gPC1+gPC2+gPC3)
capture.output(summary(lm.masc.mhc.h.male)$coefficients,file="../Results/Masc_v_mhc/lm_overallmasc_mhc_h_male_03292018.txt")
capture.output(summary(lm.masc.mhc.h.female)$coefficients,file="../Results/Masc_v_mhc/lm_overallmasc_mhc_h_female_03292018.txt")



#test whether the two slopes are different from each other

beta.m_f<-summary(lm.masc.mhc.h.male)$coefficients[2,1]-summary(lm.masc.mhc.h.female)$coefficients[2,1]
se.beta.m_f<-sqrt((summary(lm.masc.mhc.h.male)$coefficients[2,2])+(summary(lm.masc.mhc.h.female)$coefficients[2,2]^2))
z.m_f<-beta.m_f/se.beta.m_f
p.m_f<-pnorm(abs(z.m_f),lower.tail=F)

#they are not significantly different from one another

#scale qlmasc
qlmasc.u<-apply(qlmasc.u,2,scale)


#testing masculinity by height relationship without interaction term with sex
#only going to save coefficients from HLA and genetic heterozygosity
masc.mhc.b<-matrix(NA,7150,2)
masc.mhc.p<-matrix(NA,7150,2)
pb<-txtProgressBar(min=0,max=7150,style=3)
for(i in 1:7150){
  l1m<-lm(data=eurofam.het,qlmasc.u[,i]~phet_hla+phet_genome+Height+Sex+Age+Weight+gPC1+gPC2+gPC3)
  s1m<-summary(l1m)
  masc.mhc.b[i,]<-s1m$coefficients[c(2:3),1]
  masc.mhc.p[i,]<-s1m$coefficients[c(2:3),4]
  setTxtProgressBar(pb,i)
}

masc.mhc.q<-masc.mhc.p
masc.mhc.q[which(masc.mhc.p<0.05/7150)]<-1
masc.mhc.q[which(masc.mhc.p>=0.05/7150)]<-0

colnames(masc.mhc.b)<-colnames(masc.mhc.p)<-colnames(masc.mhc.q)<-c("hla.het","gen.het")
#write to file
write.table(masc.mhc.b,"qlmasc_hlahet_b_noprop_04142018.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(masc.mhc.p,"qlmasc_hlahet_p_noprop_04142018.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(masc.mhc.q,"qlmasc_hlahet_q_noprop_04142018.txt",sep="\t",col.names=T,row.names=F,quote=F)

#testing different slopes b/w masc and height in men and women

qlmasc.m<-apply(qlmasc.m,2,scale)
qlmasc.f<-apply(qlmasc.f,2,scale)

masc.hla.z<-matrix(NA,nrow=7150,ncol=5)
colnames(masc.hla.z)<-c("slope_female","slope_male","diff","z","p.value")
pb<-txtProgressBar(min=0,max=7150,style=3)
for(i in 1:7150){
  datm<-cbind(qlmasc.m[,i],eurofam.het[eurofam.het$Sex=="Male",])
  colnames(datm)[1]<-"ql"
  datf<-cbind(qlmasc.f[,i],eurofam.het[eurofam.het$Sex=="Female",])
  colnames(datf)[1]<-"ql"
  #residualize
  l1m<-lm(data=datm,ql~phet_hla+phet_genome+Height+Age+Weight+gPC1+gPC2+gPC3)
  l1f<-lm(data=datf,ql~phet_hla+phet_genome+Height+Age+Weight+gPC1+gPC2+gPC3)
  
  b12<-summary(l1m)$coefficients[2,1]-summary(l1f)$coefficients[2,1]
  sb12<-sqrt((summary(l1m)$coefficients[2,2])+(summary(l1f)$coefficients[2,2]^2))
  z<-b12/sb12
  
  masc.hla.z[i,1]<-summary(l1f)$coefficients[2,1]
  masc.hla.z[i,2]<-summary(l1m)$coefficients[2,1]
  masc.hla.z[i,3]<-b12
  masc.hla.z[i,4]<-z
  masc.hla.z[i,5]<-pnorm(z,lower.tail=F)
  setTxtProgressBar(pb,i)
}

masc.hla.z<-as.data.frame(masc.hla.z)
masc.hla.z$q.value<-NA
masc.hla.z$q.value[which(masc.hla.z$p.value<(0.05/7150))]<-1
masc.hla.z$q.value[which(masc.hla.z$p.value>=(0.05/7150))]<-0

write.table(masc.hla.z,'qlmasc_hlahet_bysex_noprop_04142018.txt',sep="\t",col.names=T,row.names=F,quote=F)
