# Endometrial Cancer

#Details: 
*  all commands in italic
*  paths in bold

**Raw dataset PATH: /local/chib/oconnor_genomes/EndometrialCancer/data/PAIRED**

Determining the read length

*zcat r_2006_FSFP192242716-1a_HWV27DSXX_L3_1.fq.gz | head -100 | awk '{if(NR%4==2) print length($1)}'*

- print the first 100th lines then every four lines select the second one (NR%4==2) and print the length of it.




