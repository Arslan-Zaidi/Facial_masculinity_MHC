library(ggplot2)
library(dplyr)
library(data.table)
library(WebPower)
library(cowplot)

#read file containing data for sex, height, weight, etc.
eurofam.het<-fread('../Dataset/euro_1233_masc_het_03292018.dat',header=T,sep="\t",stringsAsFactors = F)

#center and scale data to generate standardized regression coefficients
# eurofam.het[,c('avg.masc','Age','Height','Weight','phet_genome','phet_hla','gPC1','gPC2','gPC3','gPC4')]<-apply(eurofam.het[,c('avg.masc','Age','Height','Weight','phet_genome','phet_hla','gPC1','gPC2','gPC3','gPC4')],2,scale)

#full model with all covariates
lm.full<-lm(data=eurofam.het,avg.masc~phet_hla+phet_genome+Height+Sex+Age+Weight+gPC1+gPC2+gPC3)

#reduced model - without sex
lm.sex<-lm(data=eurofam.het,avg.masc~phet_hla+phet_genome+Height+Age+Weight+gPC1+gPC2+gPC3)

#reduced model - withut MHC heterozygosity
lm.mhc<-lm(data=eurofam.het,avg.masc~phet_genome+Height+Sex+Age+Weight+gPC1+gPC2+gPC3)

#get effect sizes for sex and MHC (Cohen's f2)
s.full=summary(lm.full)
s.sex=summary(lm.sex)
s.mhc=summary(lm.mhc)

f2sex=(s.full$r.squared - s.sex$r.squared)/(1-s.full$r.squared)
f2mhc=(s.full$r.squared - s.mhc$r.squared)/(1-s.full$r.squared)

#now let's do power analysis for varying effect sizes - scaled by effect size of sex
sim.f2=seq(0.01,1,0.01)
sim.f2=sim.f2*f2sex

#sample sizes
sim.n=seq(20,2000,10)

#create matrix and store power in it for alpha=0.05
mat1<-matrix(NA,nrow=length(sim.n),ncol=length(sim.f2))
for(i in 1:length(sim.n)){
  for(j in 1:length(sim.f2)){
    result=wp.regression(n=sim.n[i],p1=9,p2=8,f2=sim.f2[j])
    mat1[i,j]=result$power
  }
}

colnames(mat1)=sim.f2
dat1<-as.data.frame(mat1)
dat1.2<-dat1%>%
  mutate(N=sim.n)%>%
  melt(.,id.var="N",value.name="power")%>%
  mutate(f2=round(as.numeric(as.character(variable))/f2sex,3))

#plot power curves for different effect sizes
plt1<-ggplot(dat1.2,aes(N,power,color=f2,group=f2))+
  geom_line()+
  theme_bw()+
  scale_color_viridis_c()+
  labs(x="Sample size",
       y="Power",
       color="Effect size",
       title=expression(paste(alpha,"=0.05",sep="")))+
  geom_vline(xintercept=1233,color="red",linetype="dashed")+
  theme(panel.grid = element_blank())

#calculate power with sample size of 1230
pwr1<-dat1.2%>%
  filter(N==1230)%>%
  group_by(N,f2)%>%
  summarize(power=power*100)

#recompute power for alpha=0.05/9=0.006
mat2<-matrix(NA,nrow=length(sim.n),ncol=length(sim.f2))
for(i in 1:length(sim.n)){
  for(j in 1:length(sim.f2)){
    result=wp.regression(n=sim.n[i],p1=9,p2=8,f2=sim.f2[j],alpha=0.006)
    mat2[i,j]=result$power
  }
}

colnames(mat2)=sim.f2
dat2<-as.data.frame(mat2)
dat2.2<-dat2%>%
  mutate(N=sim.n)%>%
  melt(.,id.var="N",value.name="power")%>%
  mutate(f2=round(as.numeric(as.character(variable))/f2sex,3))

#plot power curves
plt2<-ggplot(dat2.2,aes(N,power,color=f2,group=f2))+
  geom_line()+
  theme_bw()+
  scale_color_viridis_c()+
  labs(x="Sample size",
       y="Power",
       color="Effect size",
       title=expression(paste(alpha,"=0.05/9 ~ 0.006",sep="")))+
  geom_vline(xintercept=1233,color="red",linetype="dashed")+
  theme(panel.grid = element_blank())

#calculate power with sample size of 1230
pwr2<-dat2.2%>%
  filter(N==1230)%>%
  group_by(N,f2)%>%
  summarize(power=power*100)

#combine both power curves
plt_all<-plot_grid(plt1+theme(legend.position = "none"),
                   plt2+theme(legend.position = "none"),
                   ncol=2,
                   align="h",
                   labels=c("A","B"))

legend<-get_legend(plt1)
plt_all2<-plot_grid(plt_all,
                    legend,
                    rel_widths = c(1,0.2))

#save plot to file
ggsave("power_calculations.pdf",plt_all2,height=3,width=7)