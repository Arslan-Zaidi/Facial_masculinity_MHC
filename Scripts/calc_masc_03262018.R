#Author: AAZaidi
#This scripts calculates facial masculinity from 3D Procrustes coordinates
#Can be modified to work with 2D coordinates as well
library(data.table)
library(plyr)
library(ggplot2)

#read in symspace - size corrected
symql<-as.matrix(fread('~/Box Sync/MHC_paper/Data_from_git/SymShape_wosize_02102018.txt',header=F))
symid<-read.table('~/Box Sync/MHC_paper/Data_from_git/SymShape_IDs_wosize_02102018.txt',header=F)

#read in data file containing individual ddata, sex, height, weight, age etc. 
eurofam<-read.table('../Dataset/Euro_demographic_03292018.dat',header=T,sep="\t",stringsAsFactors=F)


#filter symql file for europeans
sym.euro.index<-which(symid$V1%in%eurofam$IID)
euroql<-symql[,sym.euro.index]

#transpose so individuals are rows and qls are columns
teuroql<-t(euroql)

#create new dataframe for europeans
euroid<-data.frame(IID=as.character(symid$V1[sym.euro.index]))
euroid<-join(euroid,eurofam,by="IID")

#list indices corresponding to order in symql file for males and females so they can be separated later
euro.m.index<-which(euroid$Sex=="Male")
euro.f.index<-which(euroid$Sex=="Female")

#create male consensus face
cons.males1<-apply(teuroql[euro.m.index,],2,mean)
cons.males<-matrix(cons.males1,7150,3,byrow=T)

#create female consensus face
cons.females<-matrix(apply(teuroql[euro.f.index,],2,mean),7150,3,byrow=T)

#subtract female face from male face
cons.diff<-cons.males-cons.females

#define function to calculate euclidean distance of vector
euclid<-function(x){sqrt(sum(x^2))} 

#calculate degree of sexual dimorphism per QL (geometric sexual dimorphism - gsd)
gsd<-apply(cons.diff,1,euclid)

write.table(gsd,'../Results/Summary_dat/ql_gsd_03292018.txt',sep="\t",col.names=T,row.names=F,quote=F)

#function to calculate centroid size
cal.centroid<-function(x){
  #use teurosym - where rows are individuals and columns are QLs
  indices<-seq(1,7150,1)
  x1<-3*indices - 2
  y1<-x1+1
  z1<-x1+2
  x2<-x[,x1]
  y2<-x[,y1]
  z2<-x[,z1]
  #calculate mean of x coordinates across landmarks
  x2.mean<-apply(x2,1,mean)
  y2.mean<-apply(y2,1,mean)
  z2.mean<-apply(z2,1,mean)
  x3<-apply(x2,2,function(y){(y-x2.mean)^2})
  y3<-apply(y2,2,function(y){(y-y2.mean)^2})
  z3<-apply(z2,2,function(y){(y-z2.mean)^2})
  overall<-sqrt(x3+y3+z3)
  overall<-apply(overall,1,mean)
  return(overall)
}

#plot histogram of centroids - make sure they are all around 1 - indicates size scaling
#hist(cal.centroid(t(euroql)))

#define function to project one vector x onto another y
proj.vec<-function(x,y){
  #projection of x onto y 
  compx.y<-(y%*%x)/euclid(y)
  #projx.y<-compx.y * (y)/(euclid(y))
  return(compx.y)
}

#calculate masculinity for each person per QL
qlmasc.mat<-matrix(NA,nrow(euroid),7150)
pb <- txtProgressBar(min = 0, max = 7150, style = 3)

for(i in 1:7150){
  qlindex<-3*i - 2
  target.df<-teuroql[,c(qlindex,qlindex+1,qlindex+2)]
  target.df<-t(apply(target.df,1,function(x){x-cons.females[i,]}))
  cons.ql<-cons.diff[i,]
  qlmasc.mat[,i]<-apply(target.df,1,function(x){
    result<-proj.vec(x,cons.ql)
    return(result)
    #result<-result/euclid(cons.ql)
  })
  setTxtProgressBar(pb,i)
}

#average masculinity over QL for each person and add to dataframe 
avg.masc<-apply(qlmasc.mat,1,mean)

euroid$avg.masc<-avg.masc

#scale average masculinity by the distance between male and female consensus face
#this to make masculinity more interpretable
euroid$avg.masc.unit<-euroid$avg.masc/mean(gsd)

#plot relationship between height and average masculinity
ggplot(euroid,aes(Height,avg.masc.unit,color=Sex))+geom_point()+stat_smooth(method="lm")

#write average masculinity data to file
write.table(euroid,'../Dataset/euro_1233_masc_het_03292018.dat',sep="\t",col.names=T,row.names=F,quote=F)

#write matrix of facial masculinity individual x QL to file
fwrite(as.data.table(qlmasc.mat),"../Dataset/qlmasc_1233_noprop_03292018.txt",sep="\t",col.names=F,row.names=F,quote=F)

