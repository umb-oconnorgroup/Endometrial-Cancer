library("SKAT")
library(dplyr)

### creating SetID file

genotypes.path <-"/Users/vborda/Documents/UMB/Endoseq_project/SKAT/inputs/"
#######################################################
library("SKAT")
library(dplyr)

covariates_path<-"/Users/vborda/Documents/UMB/Endoseq_project/SKAT/inputs/"
pca_path<-"/Users/vborda/Documents/UMB/Endoseq_project/SKAT/inputs/"
cov_name<-paste0(covariates_path,"All_pheno_geno_SDBlack")
pca_name<-paste0(pca_path,"Principal_Components_Endometrial_SD.txt")


##################### RUNNING SKAT #####################

genotypes.path <-"/Users/vborda/Documents/UMB/Endoseq_project/SKAT/inputs/"
#clean.file<- "Nonsense_annotated_nonsynonymous_biallelic_SDBlack"
clean.file<- "Missense_Nonsense_annotated_nonsynonymous_biallelic_SDBlack"

File.Bed <- paste0(genotypes.path,clean.file,".bed")
File.Bim <- paste0(genotypes.path,clean.file,".bim")
File.Fam <- paste0(genotypes.path,clean.file,".fam")
#File.SetID <- paste0(genotypes.path,"nonsense_endo.SetID")
File.SetID <- paste0(genotypes.path,"missense_nonsense_endo.SetID")
File.SSD <- paste0(genotypes.path,clean.file,".SSD")
File.Info <- paste0(genotypes.path,clean.file,".Info")


##Close_SSD()

###### CZS ######
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)

## for the problem Error in if (SetID1 != SetID2) { : missing value where TRUE/FALSE needed
## you have to remove NA and change for any other values, in my case I changed for VOLXXX
###########################################################
covariates<-read.table(cov_name, header = T, sep = "", fill= TRUE)
pcs<-read.table(pca_name, header = T, sep = "", fill= TRUE)
covariates$FINALHIST[which(covariates$FINALHIST2=="endometrioid")]<-0
covariates$FINALHIST[which(covariates$FINALHIST2=="serous")]<-1
cov.pcs<-merge(covariates, pcs, by.x = "StudyID", by.y = "ID",sort = FALSE)

y<-cov.pcs$FINALHIST
Age<-cov.pcs$DXAGE
pc1<-cov.pcs$PC1
pc2<-cov.pcs$PC2
pc3<-cov.pcs$PC3
pc4<-cov.pcs$PC4
pc5<-cov.pcs$PC5

##  Open SNP set data file (SSD)

SSD_INFO<-Open_SSD(File.SSD, File.Info)


list.genes<-c()
gene.sets<-read.table(File.SetID, header = F, sep = "", fill= TRUE)
gene.sets.freqs<-as.data.frame(table(gene.sets[1]))
for(i in c(1:nrow(gene.sets.freqs))){
 if(gene.sets.freqs$Freq[i]>=6){
    list.genes<-c(list.genes,gene.sets.freqs$V1[i])
  }
}
length(list.genes)
bonferroni<-0.05/length(list.genes)
bonferroni<-0.05/11316
bonferroni
##40 Samples, 27169 Sets, 370061 Total SNPs
# When you have no covariates to adjust, we type 1 but if we have some covariates we change the 1 by X

obj_adjust.6cov<-SKAT_Null_Model(y ~ Age+pc1+pc2+pc3+pc4+pc5, out_type="D",n.Resampling = 1000)

####### BINARY

output.path<-"/Users/vborda/Documents/UMB/Endoseq_project/SKAT/outputs/"
out.SKAT0.05_6<-SKAT.SSD.All(SSD_INFO ,obj_adjust.6cov)
results<-out.SKAT0.05_6$results


results[(results$SetID=="PRRC2B" |results$SetID=="POMT1"|results$SetID=="UCK1"|results$SetID=="MED27"|results$SetID=="PRRT1B"|
           results$SetID=="RAPGEF1"|results$SetID=="XRCC4"|results$SetID=="VCAN"|results$SetID=="MIOS"|results$SetID=="RPA3" |
         results$SetID=="NDUFA4"|results$SetID=="ECPAS"|results$SetID=="OR2K2"|results$SetID=="ZNF483"|results$SetID=="PTGR1"|
           results$SetID=="SHOC1" | results$SetID=="UGCG"  ),]
