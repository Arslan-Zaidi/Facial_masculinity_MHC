#Author: AAZaidi
#This script calculates loads facial masculinity and calculates summary statistics 
#e.g. Cohen's D and degree of sexual dimorphism
library(data.table)
library(plyr)
library(car)

#read masculinity calculated per QL - unscaled - no size adjustment
qlmasc.u1<-as.matrix(fread('../Dataset/qlmasc_1233_sizenotcorrected_09302018.txt',header=F,sep="\t"))
qlmasc.u2<-as.matrix(fread('../Dataset/qlmasc_1233_sizecorrected_09302018.txt',header=F,sep="\t"))

#read file containing data for sex, height, weight, etc. 
eurofam.het<-read.table('../Dataset/euro_1233_masc_het_09302018.dat',header=T,sep="\t",stringsAsFactors = F)

#separate ql masculinity table by sex for use later
qlmasc.m<-qlmasc.u[eurofam.het$Sex=="Male",]
qlmasc.f<-qlmasc.u[eurofam.het$Sex=="Female",]


#calculating mean and standard deviation of overall facial masculinity for males and females

ddply(eurofam.het,.(Sex),summarize,mean=mean(avg.masc1),sd=sd(avg.masc1))

ddply(eurofam.het,.(Sex),summarize,mean=mean(Height),sd=sd(Height))


#check for equality of variance between males and females in height
leveneTest(Height~as.factor(Sex),data=eurofam.het)

#check for equality of variance between males and females in facial masculinity
#size not adjusted
leveneTest(avg.masc1~as.factor(Sex),data=eurofam.het)

#check for equality of variance between males and females in facial masculinity
#size adjusted
leveneTest(avg.masc2~as.factor(Sex),data=eurofam.het)


#calculate cohen's D for height and overall facial masculinity
male.index<-which(eurofam.het$Sex=="Male")
female.index<-which(eurofam.het$Sex=="Female")

cohen.d<-function(x){
  #x is a vector containing phenotype data in the SAME order as the male and female indices above
  #calculate mean difference
  male.df<-x[male.index]
  female.df<-x[female.index]
  
  mean.m<-mean(male.df)
  mean.f<-mean(female.df)
  mean.m_f<-mean.m-mean.f
  
  #calculate pooled sd
  var.m<-var(male.df)
  var.f<-var(female.df)
  pooled.sd<-sqrt((var.m+var.f)/2)

  #calculate D
  D<-mean.m_f/pooled.sd
  return(D)}

#cohenD - sex difference for height
cohen.d(eurofam.het$Height)

#cohenD - sex difference for masculinity
#size not adjusted
cohen.d(eurofam.het$avg.masc1)
#size adjusted
cohen.d(eurofam.het$avg.masc2)



#calculate cohen's D per QL
ql.cohen1<-apply(qlmasc.u1,2,cohen.d)
ql.cohen2<-apply(qlmasc.u2,2,cohen.d)

write.table(ql.cohen1,'../Results/Summary_dat/ql_sex_cohenD_sizenotcorrected_09302018.txt',col.names=T,row.names=F,quote=F,sep="\t")
write.table(ql.cohen2,'../Results/Summary_dat/ql_sex_cohenD_sizecorrected_09302018.txt',col.names=T,row.names=F,quote=F,sep="\t")

