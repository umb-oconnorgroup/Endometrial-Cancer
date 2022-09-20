
######################################################################
############################# Regression #############################
######################################################################
#This script includes R code lines from a pipeline created by Dr. Douglas Loesch https://github.com/dloesch/leverage-ancestry

library(dplyr)
library(stringr)

Endodir<-"/Users/vborda/Documents/UMB/Endoseq_project/"
LA.file<-paste0(Endodir,"Local_Ancestry/MSP_all_chrs_Endometrial_4ancestries_edited.tsv") 
LA<- read.table(LA.file, head=TRUE)                 ## open local ancestry file

pos.inds<-seq(7, by = 2, len = (ncol(LA)/2)-3) 
subjects<-c()

##### Getting individual ids #####
for (i in pos.inds){
    ind.id<-gsub("N.0$","",colnames(LA[i]))
    subjects<-c(subjects,ind.id)
}

LA <- LA[order(LA$chm),]
final <- as.data.frame(subjects)
colnames(final) <- "ID"

segments <- paste0(LA[,1],":",LA[,2],"-", LA[,3])

########## Getting Local ancestry counts
#Subpopulation order/codes: AFR=0	EAS=1	EUR=2	NAT=3

#African ancestry
AFR <- LA[,7:ncol(LA)]
for(i in 1: ncol(AFR)){
  AFR[,i][AFR[,i] != 0] <- NA
  AFR[,i][AFR[,i] == 0] <- 1
  AFR[,i][is.na(AFR[,i])] <- 0
}

#AFR
afrtable <- as.data.frame(subjects)
colnames(afrtable) <- "ID"

for ( i in 1:nrow(AFR)){
  foo <- c()
  j <- 1
  while(j <= ncol(AFR)){
    foo2 <- AFR[i,j]+AFR[i,j+1]
    foo <- c(foo, foo2)
    j <- j+2
  }
  foo <- as.data.frame(foo)
  colnames(foo) <- segments[i]
  afrtable <- cbind(afrtable, foo)
}


#European ancestry
EUR <- LA[,7:ncol(LA)]
for(i in 1: ncol(EUR)){
  EUR[,i][EUR[,i] != 2] <- NA
  EUR[,i][EUR[,i] == 2] <- 1
  EUR[,i][is.na(EUR[,i])] <- 0
}

#EUR
eurtable <- as.data.frame(subjects)
colnames(eurtable) <- "ID"

for ( i in 1:nrow(EUR)){
  foo <- c()
  j <- 1
  while(j <= ncol(EUR)){
    foo2 <- EUR[i,j]+EUR[i,j+1]
    foo <- c(foo, foo2)
    j <- j+2
  }
  foo <- as.data.frame(foo)
  colnames(foo) <- segments[i]
  eurtable <- cbind(eurtable, foo)
}


write.table(eurtable,paste0(Endodir,"Local_ancestry_counts_EUR.txt"), sep='\t', col.names = T, row.names = F, quote = F)


#pca<- read.table("/Users/vborda/Documents/UMB/Endoseq_project/PCA//Principal_Components_Endometrial.txt", head=TRUE, fill = TRUE)
pca<- read.table("/Users/vborda/Documents/UMB/Endoseq_project/PCA/Principal_Components_Endometrial_SD.txt", head=TRUE, fill = TRUE)
clinic<- read.table("/Users/vborda/Documents/UMB/Endoseq_project/endoseq_survival_113021.txt", head=TRUE, fill = TRUE)
clinic_subset<-subset(clinic,select = c(1:3,8:9,19))
pca_subset<-subset(pca,select = c(1:6))

AFRsegments1<-merge(clinic_subset,afrtable,by.x="StudyID", by.y ="ID",all.y = TRUE)
AFRsegments_clinic<-merge(pca_subset,AFRsegments1,by.x="ID", by.y ="StudyID",all.y = TRUE)
AFRsegments_clinic_twohist<-AFRsegments_clinic[(AFRsegments_clinic$FINALHIST2=="endometrioid" | AFRsegments_clinic$FINALHIST2=="serous"),]
AFRsegments_clinic_twohist$FINALHIST[which(AFRsegments_clinic_twohist$FINALHIST2=="endometrioid")]<-0
AFRsegments_clinic_twohist$FINALHIST[which(AFRsegments_clinic_twohist$FINALHIST2=="serous")]<-1
AFRsegments_clinic_twohist<-as.data.frame(AFRsegments_clinic_twohist[(AFRsegments_clinic_twohist$DXAGE > 10),])
#AFRsegments_clinic_twohist<-as.data.frame(AFRsegments_clinic_twohist[(AFRsegments_clinic_twohist$PATBMI < 998),])
#AFRsegments_clinic_twohist<-as.data.frame(AFRsegments_clinic_twohist[(AFRsegments_clinic_twohist$DXAGE > 10),])

#AFRsegments_clinic_twohist$FINALHIST<-as.numeric(AFRsegments_clinic_twohist$FINALHIST)
AFRsegments_clinic_twohist$FINALHIST<-as.factor(AFRsegments_clinic_twohist$FINALHIST)
AFRsegments_clinic_twohist$FINALHIST2<-as.factor(AFRsegments_clinic_twohist$FINALHIST2)
AFRsegments_clinic_twohist$PTRACE<-as.factor(AFRsegments_clinic_twohist$PTRACE)
AFRsegments_clinic_twohist$DXAGE<-as.numeric(AFRsegments_clinic_twohist$DXAGE)
AFRsegments_clinic_twohist$PATBMI<-as.numeric(AFRsegments_clinic_twohist$PATBMI)
admixture<-"/home/victor/UMB/Women_Cancer/Admixture/Ancestry_proportions_Endometrial_samples.txt"
ADMIX<- subset(read.table(admixture, head=TRUE),select =c(1:5)) 
ADMIX$X.sample<-gsub("N$","",ADMIX$X.sample)
AFRsegments_clinic_twohist<-merge(ADMIX,AFRsegments_clinic_twohist,by.y="ID", by.x ="X.sample",all.y = TRUE)

EURsegments1<-merge(clinic_subset,eurtable,by.x="StudyID", by.y ="ID",all.y = TRUE)
EURsegments_clinic<-merge(pca_subset,EURsegments1,by.x="ID", by.y ="StudyID",all.y = TRUE)
EURsegments_clinic_twohist<-EURsegments_clinic[(EURsegments_clinic$FINALHIST2=="endometrioid" | EURsegments_clinic$FINALHIST2=="serous"),]
EURsegments_clinic_twohist$FINALHIST[which(EURsegments_clinic_twohist$FINALHIST2=="endometrioid")]<-0
EURsegments_clinic_twohist$FINALHIST[which(EURsegments_clinic_twohist$FINALHIST2=="serous")]<-1
EURsegments_clinic_twohist<-as.data.frame(EURsegments_clinic_twohist[(EURsegments_clinic_twohist$DXAGE > 10),])
#EURsegments_clinic_twohist<-as.data.frame(EURsegments_clinic_twohist[(EURsegments_clinic_twohist$PATBMI < 998),])
EURsegments_clinic_twohist$FINALHIST<-as.factor(EURsegments_clinic_twohist$FINALHIST)
EURsegments_clinic_twohist$FINALHIST2<-as.factor(EURsegments_clinic_twohist$FINALHIST2)
EURsegments_clinic_twohist$PTRACE<-as.factor(EURsegments_clinic_twohist$PTRACE)
EURsegments_clinic_twohist$DXAGE<-as.numeric(EURsegments_clinic_twohist$DXAGE)
EURsegments_clinic_twohist$PATBMI<-as.numeric(EURsegments_clinic_twohist$PATBMI)


boxplot(AFRsegments_clinic_twohist$AFR~AFRsegments_clinic_twohist$FINALHIST2)
compare_means(FINALHIST~AFR,  data = AFRsegments_clinic_twohist)

library("SKAT")
library("SNPRelate")
library("plyr")
library("logistf")
library("tidyr")

AFRsegments_clinic_twohist<-AFRsegments_clinic_twohist[(AFRsegments_clinic_twohist$PTRACE==2),]
EURsegments_clinic_twohist<-EURsegments_clinic_twohist[(EURsegments_clinic_twohist$PTRACE==2),]


p_snps<-c()
alpha=0.05
col.names.table.adjust<-c("TestedAllele","TestedAllele2","p-value","Beta","OR",paste("lower", 100 - 100 *alpha, "ci", sep = ""), paste("upper", 100 - 100 *alpha, "ci", sep = ""))
p_snps_adjust <- read.table(text = "", col.names = col.names.table.adjust)

for (i in c(12:ncol(EURsegments_clinic_twohist))){
  counting<-i-11
  adj<-logistf(data=EURsegments_clinic_twohist,FINALHIST ~ EURsegments_clinic_twohist[,i] + DXAGE + PC1 + PC2 + PC3 + PC4 +PC5,firth = TRUE, pl=TRUE) #  +PC3 + PC4 + PC5 + DXAGE
  p_snps_adjust[counting,1]<-colnames(EURsegments_clinic_twohist[i])
  p_snps_adjust[counting,2]<-colnames(EURsegments_clinic_twohist[i])
  p_snps_adjust[counting,3]<-as.numeric(adj$prob[2])
  p_snps_adjust[counting,4]<-as.numeric(adj$coefficients[2])                #Beta
  p_snps_adjust[counting,5]<-exp(as.numeric(adj$coefficients[2]))                #calculating OR
  p_snps_adjust[counting,6]<-exp(as.numeric(adj$ci.lower[2])) #calculating lower 
  p_snps_adjust[counting,7]<-exp(as.numeric(adj$ci.upper[2]))
}

p_snps_adjust<-p_snps_adjust %>% separate("TestedAllele2", c("start1", "end"), "-")
p_snps_adjust<-p_snps_adjust %>% separate("start1", c("chr", "start"), ":")
p_snps_adjust$start<-as.numeric(p_snps_adjust$start)
p_snps_adjust$chr<-as.numeric(p_snps_adjust$chr)

outputdir<-"/Users/vborda/Documents/UMB/Endoseq_project/plots/"
write.table(p_snps_adjust,paste0(outputdir,"Firth_Regression_counts_ALL_EUR_association_DXAGE+_SD_PCs12345_OR.txt"), sep='\t', col.names = T, row.names = F, quote = F)

extract1<-AFRsegments_clinic_twohist
observe<-extract1[,c(1,7:8,2439)]#2402 1878 2400 2402 1503 2438 #899 No signal
observe<-na.omit(observe)
table(observe$FINALHIST2,observe[,4])

library(lattice)
source("/home/victor/UMB/Women_Cancer/Association/Manhattan_plot.function.R")
ann<-rep(1, length(data.pvalues$pvalue))
ann[with(data.pvalues,chr==5 & start>=82317962 & start<88823722)]<-2
ann[with(data.pvalues, chr==9 & start>=112659325 & start<113969108)]<-3
ann[with(data.pvalues, chr==9 & start>=131432751 & start<132482243)]<-4

ann<-factor(ann, levels=1:4, labels=c("","chr5(q14.2-q14.3)","chr9(q32)","chr9(q34.13)"))
#draw plot with annotation

png("/home/victor/UMB/Women_Cancer/Association/Manhattan_SD_EUR_association_DXAGE+_SD_PCs12345.png", width=950, height=500)
manhattan.plot(data.pvalues$chr,data.pvalues$start,data.pvalues$pvalue,
               annotate=list(ann,"chr5(q14.2-q14.3)"=list(col="darkslateblue",cex=2,label=list(x=list(5,4.5),fontsize=15,offset=1)),
                                 "chr9(q32)"=list(col="indianred3",cex=2,label=list(x=list(9,4.5),fontsize=15)),
                                 "chr9(q34.13)"=list(col="greenyellow",cex=2,label=list(x=list(11,5),fontsize=15,offset=0.5))),sig.level=2.1e-05)
dev.off()
