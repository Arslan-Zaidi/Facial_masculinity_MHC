#Author: AAZaidi
#tests for association between height and MHC heterozygosity
library(data.table)
library(ggplot2)

#read file containing data for sex, height, weight, etc. 
eurofam.het<-read.table('../Dataset/euro_1233_masc_het_03292018.dat',header=T,sep="\t",stringsAsFactors = F)

#center and scale data to generate standardized regression coefficients
eurofam.het[,c('phet_hla','phet_genome','gPC1','gPC2','gPC3','gPC4','Age','Height','Weight','avg.masc')]<-apply(eurofam.het[,c('phet_hla','phet_genome','gPC1','gPC2','gPC3','gPC4','Age','Height','Weight','avg.masc')],2,scale)

#test for an association betweeh height and hla heterozygosity
lm.height.hla.all<-lm(data=eurofam.het,Height~phet_hla+phet_genome+Sex+Age+gPC1+gPC2+gPC3)
capture.output(summary(lm.height.hla.all)$coefficients,file="../Results/Height_v_mhc/lm_height_hla_results_04142018.txt")


#separately within each sex
lm.height.hla.male<-lm(data=eurofam.het[which(eurofam.het$Sex=="Male"),], Height~phet_hla+phet_genome+Age+gPC1+gPC2+gPC3)
capture.output(summary(lm.height.hla.male)$coefficients,file="../Results/Height_v_mhc/lm_height_hla_results_male_04142018.txt")

lm.height.hla.female<-lm(data=eurofam.het[which(eurofam.het$Sex=="Female"),], Height~phet_hla+phet_genome+Age+gPC1+gPC2+gPC3)
capture.output(summary(lm.height.hla.female)$coefficients,file="../Results/lm_height_hla_results_female_04142018.txt")

#test whether the two slopes are different from each other
beta.m_f<-summary(lm.height.hla.male)$coefficients[2,1]-summary(lm.height.hla.female)$coefficients[2,1]
se.beta.m_f<-sqrt((summary(lm.height.hla.male)$coefficients[2,2]^2)+(summary(lm.height.hla.female)$coefficients[2,2]^2))
z.m_f<-beta.m_f/se.beta.m_f
p.m_f<-pnorm(abs(z.m_f))
