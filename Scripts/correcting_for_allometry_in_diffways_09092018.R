#Author: AAZaidi
#This scripts calculates facial masculinity from 3D Procrustes coordinates
#Can be modified to work with 2D coordinates as well
library(data.table)
library(plyr)
library(ggplot2)
library(GGally)

#read in symspace - size corrected - not available on github
symql<-as.matrix(fread('~/Box Sync/MHC_paper/Archived/Data_from_git/Symshape_wo_Size/SymShape_wosize_02102018.txt',header=F))
symid<-read.table('~/Box Sync/MHC_paper/Archived/Data_from_git/Symshape_wo_Size/SymShape_IDs_wosize_02102018.txt',header=F)

#read in data file containing individual ddata, sex, height, weight, age etc. 
eurofam<-read.table('~/Box Sync/MHC_paper/Euro_demographic_03262018.dat',header=T,sep="\t",stringsAsFactors=F)


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

##### 1ST METHOD: CALCULATE CONSENSUS FACES/MASCULINITY WITHOUT CORRECTING FOR HEIGHT #######

#create male consensus face
cons.males1<-matrix(apply(teuroql[euro.m.index,],2,mean),7150,3,byrow=T)

#create female consensus face
cons.females1<-matrix(apply(teuroql[euro.f.index,],2,mean),7150,3,byrow=T)

#subtract female face from male face
cons.diff1<-cons.males1-cons.females1

#define function to calculate euclidean distance of vector
euclid<-function(x){sqrt(sum(x^2))} 

#calculate degree of sexual dimorphism per QL (geometric sexual dimorphism - gsd)
gsd1<-apply(cons.diff1,1,euclid)

#define function to project one vector x onto another y
proj.vec<-function(x,y){
  #projection of x onto y 
  compx.y<-(y%*%x)/euclid(y)
  #projx.y<-compx.y * (y)/(euclid(y))
  return(compx.y)
}

#calculate masculinity for each person per QL
qlmasc.mat1<-matrix(NA,nrow(euroid),7150)
pb <- txtProgressBar(min = 0, max = 7150, style = 3)

for(i in 1:7150){
  qlindex<-3*i - 2
  target.df<-teuroql[,c(qlindex,qlindex+1,qlindex+2)]
  target.df<-t(apply(target.df,1,function(x){x-cons.females1[i,]}))
  cons.ql<-cons.diff1[i,]
  qlmasc.mat1[,i]<-apply(target.df,1,function(x){
    result<-proj.vec(x,cons.ql)
    return(result)
    #result<-result/euclid(cons.ql)
  })
  setTxtProgressBar(pb,i)
}

#average masculinity over QL for each person and add to dataframe 
avg.masc1<-apply(qlmasc.mat1,1,mean)

euroid$avg.masc1<-avg.masc1


################### 2ND METHOD: RESIDUALIZE QL COORDINATES ON HEIGHT BEFORE CONSTRUCTING CONSENSUS FACES ####################

#now residualize teuroql on height and repeat
#function to do so per QL - then we can iterate over QLs
residualize.height<-function(x){
  l1<-lm(data=euroid,x~Height)
  r1<-resid(l1)
  return(as.matrix(r1))
}

#residualize each of the 21,450 QL coordinates on height
reuroql<-apply(teuroql,2,residualize.height)

#create male consensus face
cons.males2<-matrix(apply(reuroql[euro.m.index,],2,mean),7150,3,byrow=T)

#create female consensus face
cons.females2<-matrix(apply(reuroql[euro.f.index,],2,mean),7150,3,byrow=T)

#subtract female face from male face
cons.diff2<-cons.males2-cons.females2

#calculate degree of sexual dimorphism per QL (geometric sexual dimorphism - gsd)
gsd2<-apply(cons.diff2,1,euclid)

#calculate masculinity for each person per QL
qlmasc.mat2<-matrix(NA,nrow(euroid),7150)
pb <- txtProgressBar(min = 0, max = 7150, style = 3)

for(i in 1:7150){
  qlindex<-3*i - 2
  target.df<-reuroql[,c(qlindex,qlindex+1,qlindex+2)]
  target.df<-t(apply(target.df,1,function(x){x-cons.females2[i,]}))
  cons.ql<-cons.diff2[i,]
  qlmasc.mat2[,i]<-apply(target.df,1,function(x){
    result<-proj.vec(x,cons.ql)
    return(result)
    #result<-result/euclid(cons.ql)
  })
  setTxtProgressBar(pb,i)
}

#average masculinity over QL for each person and add to dataframe 
avg.masc2<-apply(qlmasc.mat2,1,mean)

euroid$avg.masc2<-avg.masc2


###### 3RD METHOD: RESIDUALIZE MASCULINITY CALCULATED PER QL ON HEIGHT ##########


#a third way is to residualize on height after masculinity is calculated per QL 
qlmasc.mat3<-apply(qlmasc.mat1,2,residualize.height)

#average masculinity over QL for each person and add to dataframe 
avg.masc3<-apply(qlmasc.mat3,1,mean)

euroid$avg.masc3<-avg.masc3

####### METHOD 4: apply height correction to average masculinity from Method 1 ########

avg.masc4<-resid(lm(data=euroid,avg.masc1~Height))
euroid$avg.masc4<-avg.masc4

colnames(euroid)[c(12:15)]<-c("Method 1","Method 2","Method 3","Method 4")



#### PLOT ALL COMPARISONS #####

pm<-ggpairs(euroid,columns=c("Method 1","Method 2","Method 3","Method 4"),aes(color=Sex),
            lower=list(continuous=wrap("smooth",alpha=0.5)),
            upper=list(continuous=wrap("cor",size=4,alignPercent=0.8)),
            diag=list(continuous=wrap("densityDiag",alpha=0.5)))+
  theme_bw()+theme(axis.text.x = element_text(angle=90))

for(j in 1:pm$ncol){
  for(i in 1:pm$nrow){
    pm[i,j]<-pm[i,j]+scale_color_manual(values=c("#ff7f00","#377eb8"))+
      scale_fill_manual(values=c("#ff7f00","#377eb8"))
  }
}

ggsave("../Results/Summary_dat/masc_calc_method_comparison_09122018.pdf",pm,height=5,width=7)

