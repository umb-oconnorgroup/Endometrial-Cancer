#######################################################################################3
#######################################################################################3
### R script to determine the distribution of somatic variants per ancestry diplotype
## The result corresponds to Figure 3B
### Counting somatic mutations

Endodir<-"/home/victor/UMB/Women_Cancer/Somatic_mutations/"
LA.file<-paste0(Endodir,"/MSP_all_chrs_Endometrial_4ancestries_edited.tsv")  ## remember the space between n snps

## AFR=0	EAS=1	EUR=2	NAT=3

LA<- read.table(LA.file, head=TRUE)                 ## open local ancestry file
pos.inds<-seq(7, by = 2, len = (ncol(LA)/2)-3)      ## getting the position of the individuals in the Local ancestry file
diplotypes<-c("EUR_EUR","EUR_AFR","EUR_NAT","EUR_EAS","AFR_AFR","AFR_NAT","AFR_EAS","NAT_NAT","NAT_EAS","EAS_EAS")
matrix.cols<-c(1,2,3)
mat1.data<-c(rep(0,30))
final_matrix <- matrix(mat1.data,nrow=10,byrow=TRUE,dimnames=list(diplotypes,matrix.cols))

##################################################
################ to sum matrixes  ################

add_matrices_1 <- function(...) {
  a <- list(...)
  cols <- sort(unique(unlist(lapply(a, colnames))))
  rows <- sort(unique(unlist(lapply(a, rownames))))
  out <- array(0, dim=c(length(rows), length(cols)), dimnames=list(rows,cols))
  for(M in a) { out[rownames(M), colnames(M)] <- out[rownames(M), colnames(M)] + M }
  out
}

###############################################################################
######################  Counting the number of diplotypes #####################
#####################  and generating table of total sizes ####################
###############################################################################

library(dplyr)
library(stringr)

ethnic<- read.table("/home/victor/UMB/Women_Cancer/Somatic_mutations/EndoSeq.pop1", head=TRUE)
blackSD<-ethnic[(ethnic$Ethnicity=="Black"),]
pos.inds<-seq(7, by = 2, len = (ncol(LA)/2)-3) 
subjects<-c()
##### Getting individual ids #####
for (i in pos.inds){
  ind.id<-gsub("N.0$","",colnames(LA[i]))
  subjects<-c(subjects,ind.id)
}

SD_withindo<-intersect(subjects,blackSD$ID)                   ### getting selfdescribe individuals
diplotypes<-c("EUR_EUR","EUR_AFR","EUR_NAT","EUR_EAS","AFR_AFR","AFR_NAT","AFR_EAS","NAT_NAT","NAT_EAS","EAS_EAS")
total.dips<-c()

data.total.size <- data.frame(matrix(ncol = length(diplotypes)+1, nrow = 0))
colnames(data.total.size) <- c("ID",diplotypes)

for (i in 1:length(pos.inds)){
  ind.dips<-c()
  ind.id<-gsub("N.0$","",colnames(LA[pos.inds[i]]))           ## getting the ind ID / the "$" is to set that 0 is at the end
  if(isTRUE(ind.id %in% SD_withindo)){                        #### loop version for selfdescribed
    ind.la<-subset(LA,select = c(1:3,pos.inds[i],pos.inds[i]+1))          ## extracting both haplotypes for an individual
    
    message(paste0("counting dips in ", ind.id," (",i," of ",length(pos.inds),")"))
    for (column in 4:5) {
      ind.la[,column]<-gsub(0,"AFR",ind.la[,column])
      ind.la[,column]<-gsub(1,"EAS",ind.la[,column])
      ind.la[,column]<-gsub(2,"EUR",ind.la[,column])
      ind.la[,column]<-gsub(3,"NAT",ind.la[,column])
    }
    ind.la$diplotype<-str_c(ind.la[,4],"_",ind.la[,5])  ### Merging to column of Local ancestry to generate the diplotype
    ind.la <- ind.la %>% mutate(diplotype = gsub("AFR_EUR", "EUR_AFR", diplotype))    ### Standaring the dyplotypes
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
  
}

totalcounts<-as.data.frame(table(total.dips))
colnames(totalcounts) <- c("Diplotype","Frequency")

###############################################################################
################ Counting the somatic mutations per diplotype  ################
###############################################################################
library(dplyr)
library(stringr)
final_matrix<-c()
somaticdir<-paste0(Endodir,"/list_somatic_samples_PASS/")
#somaticdir<-paste0(Endodir,"/list_somatic_samples_PASS_INDELS/")
for (i in 1:length(pos.inds)){
  ind.id<-gsub("N.0$","",colnames(LA[pos.inds[i]]))           ## getting the ind ID / the "$" is to set that 0 is at the end
  if(isTRUE(ind.id %in% SD_withindo)){                        ## Getting Self-Described individuals
    ind.la<-subset(LA,select = c(1:3,pos.inds[i],pos.inds[i]+1))          ## extracting both haplotypes for an individual
    
    message(paste0("counting somatic muts in ", ind.id," (",i," of ",length(pos.inds),")"))
    #somaticmuts<- read.table(paste0(somaticdir,"Counting_PASS_INDELS_", ind.id), head=FALSE)     ## Opening the somatic mutation file for the specific individual
    somaticmuts<- read.table(paste0(somaticdir,"Counting_PASS_", ind.id), head=FALSE)     ## Opening the somatic mutation file for the specific individual
    somaticmuts$V1<-gsub("chr","",somaticmuts$V1)
    colnames(somaticmuts)<-c("chromosome","position","REF","ALT")
    for(chr in 1:22){
      
      counts<-c()
      ### save the chr and then look for the interval between tables
      ind.la.chr<-ind.la[ind.la$chm %in% chr, ]                                              #### Extracting the Chromosome from Local ancestry file
      somaticmuts.chr<-somaticmuts[somaticmuts$chromosome %in% chr, ]                        #### Extracting the Chromosome from somatic file
      counts<-as.data.frame(table(findInterval(somaticmuts.chr$position,ind.la.chr[,2])))    #### LINE to get the density
      counts$Var1<-as.numeric(levels(counts$Var1))
      counts<-counts[(counts$Var1>0),]                ## To remove any somatic variant that is no in any Local Ancestry interval
      if(isTRUE(nrow(counts)>0)){
      
        for (column in 4:5) {
          ind.la.chr[,column]<-gsub(0,"AFR",ind.la.chr[,column])
          ind.la.chr[,column]<-gsub(1,"EAS",ind.la.chr[,column])
          ind.la.chr[,column]<-gsub(2,"EUR",ind.la.chr[,column])
          ind.la.chr[,column]<-gsub(3,"NAT",ind.la.chr[,column])
        }
        
        counts$diplotype<-str_c(ind.la.chr[counts$Var1,4],"_",ind.la.chr[counts$Var1,5])  ### Merging to column of Local ancestry to generate the diplotype
        counts <- counts %>% mutate(diplotype = gsub("AFR_EUR", "EUR_AFR", diplotype))    ### Standaring the dyplotypes
        counts <- counts %>% mutate(diplotype = gsub("NAT_EUR", "EUR_NAT", diplotype))
        counts <- counts %>% mutate(diplotype = gsub("EAS_EUR", "EUR_EAS", diplotype))
        counts <- counts %>% mutate(diplotype = gsub("NAT_AFR", "AFR_NAT", diplotype))
        counts <- counts %>% mutate(diplotype = gsub("EAS_AFR", "AFR_EAS", diplotype))
        counts <- counts %>% mutate(diplotype = gsub("EAS_NAT", "NAT_EAS", diplotype))
        
        #### Realizar el conteo
        counting.table<-table(counts$diplotype,counts$Freq)
        final_matrix<-add_matrices_1(final_matrix,counting.table)
        }
    }
  }
}

m <- final_matrix[, order(as.integer(colnames(final_matrix)))]           #### sorting the column names
#write.table(m,"/home/victor/UMB/Women_Cancer/Somatic_mutations/Somatic_counting_exclusive.txt",sep="\t", row.names = TRUE, quote = F)

df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) <- c('Frequency', 'Counts', 'DIPLOTYPE')

### Generating a traverse table

for (dip in row.names(m)){
  toplot<-as.data.frame(m[dip,])
  colnames(toplot)<-"Frequency"
  toplot$Counts<-colnames(m)
  toplot$DIPLOTYPE<-dip
  df<-rbind(df,toplot)
}

df<-df[(df$Frequency>0),]    
df$Counts<-as.numeric(df$Counts)
df$LOG_Frequency<-log(df$Frequency)

for (j in 1:nrow(df)){
  for (rowdip in 1:nrow(totalcounts)) {
    if(isTRUE(df[j,3]==totalcounts[rowdip,1]))
      df$Relative_Frequency[j]<-(df$Frequency[j]/totalcounts[rowdip,2])#*100
  }
}

threedips<-c("AFR_AFR","EUR_AFR","EUR_EUR")
df_3dips<-df[df$DIPLOTYPE%in%threedips,]

twodips<-c("AFR_AFR","EUR_AFR")
df_2dips<-df[df$DIPLOTYPE%in%twodips,]

df_afr_afr<-subset(df[df$DIPLOTYPE%in%"AFR_AFR",],select =c(2,5))
df_eur_afr<-subset(df[df$DIPLOTYPE%in%"EUR_AFR",],select =c(2,5))
df_eur_eur<-subset(df[df$DIPLOTYPE%in%"EUR_EUR",],select =c(2,5))
colnames(df_afr_afr)<-c("COUNTS","AFR_AFR")
colnames(df_eur_afr)<-c("COUNTS","EUR_AFR")
colnames(df_eur_eur)<-c("COUNTS","EUR_EUR")

merged_afr_afr_eur_afr<-merge(df_afr_afr,df_eur_afr,all.x = TRUE,all.y = TRUE, by="COUNTS")
merged_afr_afr_eur_afr[is.na(merged_afr_afr_eur_afr[,])] <- 0

ks.test(merged_afr_afr_eur_afr$FREQ_AFR_AFR,merged_afr_afr_eur_afr$FREQ_EUR_AFR, exact = FALSE)

library(ggplot2)
library(tidyverse)
library(gapminder)

df_totest<-merged_afr_afr_eur_afr %>% pivot_longer(cols = AFR_AFR:EUR_AFR)

df_2dips$DIPLOTYPE <- factor(df_2dips$DIPLOTYPE,levels = c("AFR_AFR","EUR_AFR"))

df_totest$name <- factor(df_totest$name, levels = c("AFR_AFR","EUR_AFR"))

wilcox.test(df_totest$value~df_totest$name,paired=TRUE)
df_3dips_subset<-df_3dips[(df_3dips$Counts>24),] 

library(ggplot2)
library(tidyverse)
library(gapminder)

png(paste0(Endodir,"SD_Relative_Abundance_vs_SomaticMutations_PASS_3dips_ver4.png"),height=20,width=25,res = 90,units = "cm")
df_3dips %>%
  ggplot() + 
  geom_histogram(mapping = aes(x = Counts, y = Relative_Frequency, 
                               color = DIPLOTYPE, fill = DIPLOTYPE), 
                 stat = "identity",
                 alpha = 0.4,
                 position = "identity") + 
  scale_color_manual(values = alpha(c( "blue", "gold4","firebrick"), 1),
                     labels = c("AFR/AFR","EUR/AFR","EUR/EUR")) + 
  scale_fill_manual(values = alpha(c( "blue", "gold4","firebrick"), 1),
                    labels = c("AFR/AFR","EUR/AFR","EUR/EUR")) +
  labs(x = "Number of Somatic Mutatins per region", y = "-1/log(Frequency)", 
       fill = "DIPLOTYPE",
       title = "Comparing Diplotypes") 
       #subtitle = "Overall distribution shown in gray")
dev.off()

library(ggplot2)

png(paste0(Endodir,"SD_Abundance_vs_TotalSomaticMutations_rescaled_PASS_3dips_ver2.png"),height=20,width=20,res = 90,units = "cm")
ggplot(df_3dips, aes(x=Counts, y=Relative_Frequency,fill=DIPLOTYPE)) +
  geom_bar(stat='identity')+
  facet_wrap(~DIPLOTYPE,ncol=1) +
  xlab("Total number of somatic mutations in the Local ancestry region") + ylab("LOG of Abundance of Ancestry segments  ")+ xlim(c(0, 75)) + ylim(c(0, 12.5))
dev.off()
