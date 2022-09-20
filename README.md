# Endoseq Project
# Genetic ancestries and their relationship with endometrial cancer outcomes

Local ancestry was performed using RFMIX with four super population groups as references 

#Three main analyses: 
1.  **Density distribution of somatic variants per diplotype:**<br/> 
    The script *Fig3B_DistributionOfDensityOfSomaticVariantWithinAncestryDiplotypes.R* generates the frequency distribution of somatic mutations for each         diplotype. This analysis requires for Local ancestry inferences from germline data and somatic mutation counts file with the format: **chr\tpos\tREF\tALT** (\t : TAB-delimited)
   
2. **Comparing Genome-wide average rate of somatic mutations across diplotypes in each individual:**<br/>
    The script *Fig3C_Comparing_Genome-wide_Average_of_Sommuts_Across_Diplotypes.R* determine (i) the total lenght for each diplotype (i.e. Eur/Eur, Afr/Afr, and Eur/Afr), (ii) the number of somatic mutations per diplotype per individual, and using both results will determine the genome-wide rate of somatic mutations per diplotype. After determine the average, the script will estimate the pairwise differences among diplotypes per individual and will plot the distribution of these differences. We determine if the difference is significantly away from zero as an indication that one average is greater than other averages.
3.  **Firth Regression to explore the relationship of local ancestry with total mutation burden (TMB) and histopathological outcomes:**<br/>
The script *AdmixtureMapping_FirthRegression_HistopathologicalOutcomes.R* generates a table of haplotype counts per locus for each individual from the local ancestry data (0,1, and 2 African counts in an specific loci). With this information coupled with clinical data (e.g. Age of Diagnosis,  histopathological type, etc) and principal components, this script will run a Firth regression and a manhattan plot with significant values annotated.



