#!/bin/bash
##This is map_mutant_Portable.sh
##Created by F. BESNARD and Sana Dieudonné, 3 July 2015 - Last modif 2/12/2015
# History: original script designed on MacOSX (Speedy Felix Lab), adpated for PSMN (June 2016), adapted for portable version (Feb. 2017)
# last edit: 2017-02-27

### What is does:
#From a .vcf file containing 1 sample only, it computes for each reported variant position the frequency of reads for the ALT allele over the sum of (REF + ALT) read counts.
#Results are stored in a file named $argument2.freq.txt
#Then plots of frequency and coverage are geneated with a Rcript and stored in a new folder named "plots".

### Requirements:
#This scripts uses the R script Plot-fqcy_Portable.R. The info is given through the config file
#a file containing the scaf names and their sizes

### Usage:
##Run automatically in andalousian-map_Portable.sh (step 9)
#If run independantly:
#			     ($1st arg) ($2nd arg)    ($3rd arg)
#bash map_mutant_Portable.sh myvcf.vcf outfile-prefix config.file
#Take as 1st argument .vcf file restricted to the SNPs of the mapping strain (e.g. JU170)
#2nd argument is the outfile prefix

### Start of script
source $3

printf "${col} STARTING 'map_mutant_Portable.sh' ${col}\n"
printf "\e[0;32m computing allele frequencies from file $1 \e[0m \n"
printf "\e[0;32m Output files named $2.freq.txt and $2_CandidateScaff will be generated \e[0m \n"
printf "\e[0;32m Frequency and coverage plots will be generated in a new folder named '$2_plots' \e[0m \n"

mkdir $2_plots

printf "temporary file generated \n"
grep -v "#" $1 | grep -vP "GT\t" > filter.vcf
printf "Computing frequencies in progress. Progressing through: \n"
cut -f1 $ANDALOUSIAN_FOLDER/data/$Scaff_size | while read line; do echo $line ; grep "$line" filter.vcf | \
cut -f 1,2,4,5,10 | awk -F ":" '{print $1,$3,$2}' | awk -F "," '{print $1,$2}' | \
awk '{ if ($7+$8 > 0) print $1,$2,$3,$4,$5,$6,$7+$8,$8/($7+$8); else print $1,$2,$3,$4,$5,$6,$7+$8,0}' >> ./$2_plots/$2.freq.txt; done
printf "All frequencies computed. Cleaning temp. files \n"
rm filter.vcf

printf "Generating plots \n"
cd $2_plots
Rscript $ANDALOUSIAN_FOLDER/scripts/Plot-fqcy_Portable.R $2.freq.txt $2

printf "${col} END OF SCRIPT 'map_mutant_Portable.sh' ${col}\n"
### End of script
