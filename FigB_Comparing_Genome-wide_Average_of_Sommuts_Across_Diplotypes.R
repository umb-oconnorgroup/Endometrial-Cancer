### R script to estimate the genome-wide rate of somatic mutations per individual per diplotype
### packages required: dplyr, stringr, ggpubr, ggplot2

### Counting somatic mutations

Endodir<-"/Users/vborda/Documents/UMB/Endoseq_project/"
LA.file<-paste0(Endodir,"Local_Ancestry/MSP_all_chrs_Endometrial_4ancestries_edited.tsv")  ## remember the space between n snps

## Ancestry codes: AFR=0	EAS=1	EUR=2	NAT=3

LA<- read.table(LA.file, head=TRUE)                 ## open local ancestry file
pos.inds<-seq(7, by = 2, len = (ncol(LA)/2)-3)      ## getting the position of the individuals in the Local ancestry file

################################################################################################
######################### Counting per region per individual ###################################
######################### Difference between EUR and AFR dips###################################
################################################################################################

## 1) Counting the number of diplotypes and generating table of total sizes per individual

library(dplyr)
library(stringr)

ethnic<- read.table("/Users/vborda/Documents/UMB/Endoseq_project/EndoSeq.pop1", head=TRUE)
blackSD<-ethnic[(ethnic$Ethnicity=="Black"),]
pos.inds<-seq(7, by = 2, len = (ncol(LA)/2)-3) 
subjects<-c()
##### Getting individual ids #####
for (i in pos.inds){
  ind.id<-gsub("N.0$","",colnames(LA[i]))
  subjects<-c(subjects,ind.id)
}

SD_withindo<-intersect(subjects,blackSD$ID)
diplotypes<-c("EUR_EUR","EUR_AFR","EUR_NAT","EUR_EAS","AFR_AFR","AFR_NAT","AFR_EAS","NAT_NAT","NAT_EAS","EAS_EAS")
total.dips<-c()

data.total.size <- data.frame(matrix(ncol = length(diplotypes)+1, nrow = 0))
colnames(data.total.size) <- c("ID",diplotypes)

for (i in 1:length(pos.inds)){
  ind.dips<-c()

  ind.id<-gsub("N.0$","",colnames(LA[pos.inds[i]]))           ## Getting the ind ID / the "$" is to set that 0 is at the end
  ind.la<-subset(LA,select = c(1:3,pos.inds[i],pos.inds[i]+1))          ## Extracting both haplotypes for an individual
  
  message(paste0("counting dips in ", ind.id," (",i," of ",length(pos.inds),")"))
  for (column in 4:5) {
    ind.la[,column]<-gsub(0,"AFR",ind.la[,column])
    ind.la[,column]<-gsub(1,"EAS",ind.la[,column])
    ind.la[,column]<-gsub(2,"EUR",ind.la[,column])
    ind.la[,column]<-gsub(3,"NAT",ind.la[,column])
  }
  ind.la$diplotype<-str_c(ind.la[,4],"_",ind.la[,5])  ### Merging to column of Local ancestry to generate the diplotype
  ind.la <- ind.la %>% mutate(diplotype = gsub("AFR_EUR", "EUR_AFR", diplotype))    ### Standardizing diplotypes
  ind.la <- ind.la %>% mutate(diplotype = gsub("NAT_EUR", "EUR_NAT", diplotype))
  ind.la <- ind.la %>% mutate(diplotype = gsub("EAS_EUR", "EUR_EAS", diplotype))
  ind.la <- ind.la %>% mutate(diplotype = gsub("NAT_AFR", "AFR_NAT", diplotype))
  ind.la <- ind.la %>% mutate(diplotype = gsub("EAS_AFR", "AFR_EAS", diplotype))
  ind.la <- ind.la %>% mutate(diplotype = gsub("EAS_NAT", "NAT_EAS", diplotype))
  ind.dips<-ind.la$diplotype
  total.dips<-c(total.dips,ind.dips)
  data.total.size[i,"ID"]<-ind.id
  
  for (dip in diplotypes){
    ind.la.dip<-ind.la[(ind.la$diplotype==dip),]
    ind.la.dip$dif<-ind.la.dip$epos-ind.la.dip$spos
    data.total.size[i,dip]<-sum(ind.la.dip$dif)
  }
  
}

totalcounts<-as.data.frame(table(total.dips))
colnames(totalcounts) <- c("Diplotype","Frequency")

###########################################################
## 2) Extracting somatic mutations per diplotype
###########################################################

somaticdir<-paste0(Endodir,"/Somatic_mutation_list/list_somatic_samples_PASS/")
pos.inds<-seq(7, by = 2, len = (ncol(LA)/2)-3) 

cols<-c("EUR_EUR","EUR_AFR","EUR_NAT","EUR_EAS","AFR_AFR",
        "AFR_NAT","AFR_EAS","NAT_NAT","NAT_EAS","EAS_EAS",
        "Size_EUR_EUR","Size_EUR_AFR","Size_EUR_NAT","Size_EUR_EAS","Size_AFR_AFR",
        "Size_AFR_NAT","Size_AFR_EAS","Size_NAT_NAT","Size_NAT_EAS","Size_EAS_EAS")

data.density <- data.frame(matrix(ncol = length(cols)+1, nrow = 0))
data.density.perchr2 <- data.frame(matrix(ncol = length(cols)+2, nrow = 0))
colnames(data.density) <- c("ID",cols)
#colnames(data.density.perchr2) <- c("ID","chr",cols)

for (i in 1:length(pos.inds)){
  ind.la<-subset(LA,select = c(1:3,pos.inds[i],pos.inds[i]+1))          ## extracting both haplotypes for an individual
  ind.id<-gsub("N.0$","",colnames(LA[pos.inds[i]]))           ## getting the ind ID / the "$" is to set that 0 is at the end
  
  message(paste0("counting in ", ind.id," (",i," of ",length(pos.inds),")"))
  somaticmuts<- read.table(paste0(somaticdir,"Counting_PASS_", ind.id), head=FALSE)      ## Opening the somatic mutation file for the specific individual
  somaticmuts$V1<-gsub("chr","",somaticmuts$V1)
  colnames(somaticmuts)<-c("chromosome","position","REF","ALT")
  
  all.chrs.counting<-data.frame(matrix(ncol = 5, nrow = 0))
  colnames(all.chrs.counting) <- c("chr","Var1","Freq","diplotype","size")
  #data.density.perchr1 <- data.frame(matrix(ncol = length(cols)+2, nrow = 0))
  #colnames(data.density.perchr1) <- c("ID","chr",cols)
  
  for(chr in 1:22){
    
    counts<-c()
    ### save the chr and then look for the interval between tables
    ind.la.chr<-ind.la[ind.la$chm %in% chr, ]                                              #### Extracting the Chromosome from Local ancestry file
    somaticmuts.chr<-somaticmuts[somaticmuts$chromosome %in% chr, ]                        #### Extracting the Chromosome from somatic file
    counts<-as.data.frame(table(findInterval(somaticmuts.chr$position,ind.la.chr[,2])))    #### LINE to get the density
    counts$Var1<-as.numeric(levels(counts$Var1))
    counts<-counts[(counts$Var1>0),]
    if(isTRUE(nrow(counts)>0)){
      for (column in 4:5) {
        ind.la.chr[,column]<-gsub(0,"AFR",ind.la.chr[,column])
        ind.la.chr[,column]<-gsub(1,"EAS",ind.la.chr[,column])
        ind.la.chr[,column]<-gsub(2,"EUR",ind.la.chr[,column])
        ind.la.chr[,column]<-gsub(3,"NAT",ind.la.chr[,column])
      }                                                         ### transforming number to ancestry (i.e. 0 to AFR)
      
      counts$diplotype<-str_c(ind.la.chr[counts$Var1,4],"_",ind.la.chr[counts$Var1,5])  ### Merging to column of Local ancestry to generate the diplotype
      counts <- counts %>% mutate(diplotype = gsub("AFR_EUR", "EUR_AFR", diplotype))    ### Standardizing Dyplotypes
      counts <- counts %>% mutate(diplotype = gsub("NAT_EUR", "EUR_NAT", diplotype))
      counts <- counts %>% mutate(diplotype = gsub("EAS_EUR", "EUR_EAS", diplotype))
      counts <- counts %>% mutate(diplotype = gsub("NAT_AFR", "AFR_NAT", diplotype))
      counts <- counts %>% mutate(diplotype = gsub("EAS_AFR", "AFR_EAS", diplotype))
      counts <- counts %>% mutate(diplotype = gsub("EAS_NAT", "NAT_EAS", diplotype))
      
      counts$Var1 <- paste0(ind.la.chr[counts$Var1,1],":",ind.la.chr[counts$Var1,2],"-", ind.la.chr[counts$Var1,3])
      counts$chr<-chr
      all.chrs.counting<-rbind(all.chrs.counting,counts)
    }
  }
  
  data.density[i,1]<-ind.id
  diplotypes<-c("EUR_EUR","EUR_AFR","EUR_NAT","EUR_EAS","AFR_AFR",
                "AFR_NAT","AFR_EAS","NAT_NAT","NAT_EAS","EAS_EAS")
  for (dip in diplotypes){
    colsize<-paste0("Size_",dip)
    ind.dip<-all.chrs.counting[all.chrs.counting$diplotype %in% dip, ] 
    data.density[i,dip]<-sum(ind.dip$Freq)/(data.total.size[data.total.size$ID==ind.id,dip])
    data.density[i,colsize]<-(data.total.size[data.total.size$ID==ind.id,dip])/1000000
    
  }
  
}

data.density[,][is.na(data.density[,])] <- 0

write.table(data.density,paste0(Endodir,"Somatic_all_mutation_Density_by_individual.txt"), sep='\t', col.names = T, row.names = F, quote = F)


#############################################################
##  3) Keeping self-described black individuals
#############################################################

ethnic<- read.table("/Users/vborda/Documents/UMB/Endoseq_project/EndoSeq.pop1", head=TRUE)
#admix.path<-"/home/victor/UMB/Women_Cancer/Admixture/"
#ad.file<-paste0(admix.path,"Supervised_ancestry_proportions.txt") 
#admixture<- read.table(ad.file, head=TRUE)  
#admix.density<-merge(admixture,data.density,by.x="ID", by.y = "ID")
density.ethnicity<-merge(data.density,ethnic,by.x="ID", by.y = "ID",all.x = TRUE)

clinic<- read.table("/Users/vborda/Documents/UMB/Endoseq_project/endoseq_survival_113021.txt", head=TRUE, fill = TRUE)
admix.density.ethnicity2<-merge(density.ethnicity,clinic,by.x="ID", by.y = "StudyID",all.x = TRUE)

admix.density.selfdescribed<-admix.density.ethnicity2[(admix.density.ethnicity2$Ethnicity=="Black"),]
admix.density.ethnicity<-admix.density.selfdescribed
admix.density.ethnicity<-admix.density.ethnicity[!is.na(admix.density.ethnicity$ID),]

library(ggpubr)
library(ggplot2)

admix.density.ethnicity %>% count(FINALHIST2)

#write.table(admix.density.ethnicity,"/home/victor/UMB/Women_Cancer/Somatic_mutations/Agilent/SD_inds", sep='\t', col.names = T, row.names = F, quote = F)

###keeping serous

admix.density.1<-admix.density.ethnicity[(admix.density.ethnicity$FINALHIST2=="endometrioid"),]
admix.density.2<-admix.density.ethnicity[(admix.density.ethnicity$FINALHIST2=="clear"),]
admix.density.3<-admix.density.ethnicity[(admix.density.ethnicity$FINALHIST2=="mixed"),]
admix.density.4<-admix.density.ethnicity[(admix.density.ethnicity$FINALHIST2=="serous"),]

#############################################################
#############################################################

df<-c()
cols<-c("Difference","Diplotype")
sdf1 <- data.frame(ncol = 2, nrow = 0)
sdf2 <- data.frame(ncol = 2, nrow = 0)
sdf3 <- data.frame(ncol = 2, nrow = 0)
colnames(sdf1) <- cols
colnames(sdf2) <- cols
colnames(sdf3) <- cols

lower=30
#toplot1<-as.data.frame(admix.density.1[(admix.density.1$Size_EUR_EUR>lower  & admix.density.1$Size_EUR_AFR>lower),])
toplot1<-as.data.frame(admix.density.ethnicity[(admix.density.ethnicity$Size_EUR_EUR>lower  & admix.density.ethnicity$Size_EUR_AFR>lower),])
for (i in 1:nrow(toplot1)){
  sdf1[i,1]<-toplot1[i,"EUR_EUR"] - toplot1[i,"EUR_AFR"]  
  sdf1[i,2]<-"EUR_EUR - EUR_AFR"
}

#toplot2<-as.data.frame(admix.density.1[(admix.density.1$Size_EUR_EUR>lower  & admix.density.1$Size_AFR_AFR>lower),])
toplot2<-as.data.frame(admix.density.ethnicity[(admix.density.ethnicity$Size_EUR_EUR>lower  & admix.density.ethnicity$Size_AFR_AFR>lower),])
for (i in 1:nrow(toplot2)){
  sdf2[i,1]<-toplot2[i,"EUR_EUR"] - toplot2[i,"AFR_AFR"]  
  sdf2[i,2]<-"EUR_EUR - AFR_AFR"
}

#toplot3<-as.data.frame(admix.density.1[(admix.density.1$Size_EUR_AFR>lower  & admix.density.1$Size_AFR_AFR>lower),])
toplot3<-as.data.frame(admix.density.ethnicity[(admix.density.ethnicity$Size_EUR_AFR>lower  & admix.density.ethnicity$Size_AFR_AFR>lower),])
for (i in 1:nrow(toplot3)){
  sdf3[i,1]<-toplot3[i,"AFR_AFR"] - toplot3[i,"EUR_AFR"]
  sdf3[i,2]<-"AFR_AFR - EUR_AFR"
}

df<-rbind(sdf1,sdf2,sdf3)

wilcox.test(sdf1$Difference) #"EUR_EUR - EUR_AFR" ## 0.003
wilcox.test(sdf2$Difference) #"EUR_EUR - AFR_AFR" 0.16
wilcox.test(sdf3$Difference) #"AFR_AFR - EUR_AFR" 0.001775
mean(sdf1$Difference)
 ### Violin plot
df$Diplotype<-as.factor(df$Diplotype)
compare_means(Difference ~ 0, df, group.by="Diplotype",mu=0)
 cw_summary <- df %>% group_by(Diplotype) %>% tally()

 
library(ggforce)
 
jpeg(filename = "/Users/vborda/Documents/UMB/Endoseq_project/plots/09212022_Violinplot_dips_allsomatic_above30M_SD_PASS.jpg",
      width = 38, height = 15, units = "cm", pointsize =12,  res = 600)
 
ggplot(df, aes(x=Diplotype, y=Difference,color=Diplotype)) + 
 geom_violin(trim=FALSE) + geom_boxplot(width=0.1,alpha = 0.4) +
 scale_color_manual(values=c("blue", "gold4","firebrick"))+
 #ylim(c(-0.000006, 0.000007))+
 geom_text(data = cw_summary,
            aes(Diplotype, Inf, label = n), vjust = 1) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_zoom(ylim = c(-0.0000001,0.0000001),zoom.size = 1) +
 annotate(geom = "text", x = 1, y = 0.000007, label = "0.00177", parse = TRUE, color = "black", size = 4)  +    ### SD ALL 30M
 annotate(geom = "text", x = 2, y = 0.000007, label = "0.1692", parse = TRUE, color = "black", size = 4)   +      ### SD ALL 30M
 annotate(geom = "text", x = 3, y = 0.000007, label = "0.0032", parse = TRUE, color = "black", size = 4)  ### SD ALL 20M 
   #annotate(geom = "text", x = 1, y = 0.0000024, label = "0.3902", parse = TRUE, color = "black", size = 4)  +    ### SD ALL 30M
 
 
dev.off() 
