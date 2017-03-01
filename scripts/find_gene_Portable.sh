#!/bin/bash
##This is find_gene_Portable.sh

##History:
##Original version created by F. BESNARD and Sana Dieudonné, 3 July 2015 
##Major edit 28/06/2016 to improve filtering
##Adaptation for PSMN then Portable version
##last modif 27/02/2017

#What is does:
#Calls the variants of the mutant only in the scaffolds where the mutation has been mapped
#Filters out all JU170 variations + quality filters
#Annotate the variations with snpeff

#### Requirements: ####
###Softwares and programs
# GATK (v3.6+)
# Snpeff (v4.1g+). 
###Datas 
#a config file (same as for andalousian-map_Portable.sh/GATK_fq-to-genotypeJU70.sh) for:
# 	-the reference genome
# 	-the file containing dbSNP_JU170 (the mapping strain) 
#a file containing the "Invariant" Scaffolds (ie containing no snps) between JU170 and CEW1
#	-'JU170_InvariantScaffold.interval_list'.
#The last .bam file of sequencing data (should end with .BQSR.bam)
#The file "mysample_CandidateScaff.tsv" containing the scaffolds where the mutation maps (located in a subfolder called "plots")
# gVCF files of mutagenized background strain (CEW1) and of the mapping strain
#	-CEW1.gVCF.vcf
#	-JU170.gVCF.vcf

### Usage:
#Run automatically in andalousian-map_Portable.sh (step 10)

#If it is trun independantly, It should be run as a 3rd step after:
#1. GATK_fq-to-genotypeJU170.sh
#2. map_mutant_Portable.sh
#Run from the folder corresponding to the mutant's analysis (named as "mysample" folder). This folder must contain mysample.RG.dedup.realign.BQSR.bam (generated after step1) and a "plots" folder containing "mysample_CandidateScaff.tsv"
#Run the following command: 
#bash find_gene_Portable.sh CONFIGFILE

CONFIGFILE=$1
source $CONFIGFILE

date
printf "${col} STARTING 'find_gene_Portable.sh' ${NC} \n"
printf "${col} find_gene_Portable.sh will find and annotate specific variations for the mutant $OUTPUT_NAME ${NC} \n"
printf "${col} Variation Call will be limited to intervalls defined in ${OUTPUT_NAME}_CandidateScaff.tsv ${NC} \n"

mkdir ${OUTPUT_NAME}_find_gene

#1. Format correct list of intervals for GATK
cut -f1 ./${OUTPUT_NAME}_plots/${OUTPUT_NAME}_CandidateScaff.tsv > ./${OUTPUT_NAME}_find_gene/${OUTPUT_NAME}.interval_list

#Add the scaffolds were JU170 is invariant (ie blind regions to the mapping approach)
cat $ANDALOUSIAN_FOLDER/data/$Invariant_Scaff >> ./${OUTPUT_NAME}_find_gene/${OUTPUT_NAME}.interval_list
#Sort and overwrite the list:
sort -o ./${OUTPUT_NAME}_find_gene/${OUTPUT_NAME}.interval_list ./${OUTPUT_NAME}_find_gene/${OUTPUT_NAME}.interval_list

cd ${OUTPUT_NAME}_find_gene

#2. Call variants
lastbam=$(echo `ls ../*.BQSR.bam`)
printf "${col} A variant call will be performed on $lastbam ${NC} \n"

printf "${col} ##### Call variants on selected intervals in a gVCF mode ##### ${NC} \n"
java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T HaplotypeCaller \
     -R $REFERENCE \
     --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
     -L ${OUTPUT_NAME}.interval_list\
     -I $lastbam \
     --dbsnp $ANDALOUSIAN_FOLDER/data/$dbSNP \
     -o $OUTPUT_NAME.gVCF.vcf

printf "${col} ##### Genotype the positions for $OUTPUT_NAME together with $Mapping_strain and $Background_strain ##### ${NC} \n"
java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T GenotypeGVCFs \
   -R $REFERENCE \
   --variant $ANDALOUSIAN_FOLDER/data/$Background_gVCF \
   --variant $ANDALOUSIAN_FOLDER/data/$Mapping_gVCF \
   --variant $OUTPUT_NAME.gVCF.vcf \
   -o ${OUTPUT_NAME}_raw.vcf

#3.Filter out bacground variants (those already present in the background or mapping strains), loci of coverage < 3 and hets
printf "${col} ##### Filter call set ##### ${NC} \n"
java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFERENCE \
   --variant ${OUTPUT_NAME}_raw.vcf \
   -L ${OUTPUT_NAME}.interval_list \
   -select 'DP >= 3 && !vc.hasAttribute("DB")' \
   -select 'vc.getGenotype("'$OUTPUT_NAME'").isHom() && vc.getGenotype("'$Background_strain'").isHom() && vc.getGenotype("'$Mapping_strain'").isHom()' \
   -select 'vc.getGenotype("'$OUTPUT_NAME'").alleles != vc.getGenotype("'$Background_strain'").alleles && vc.getGenotype("'$OUTPUT_NAME'").alleles != vc.getGenotype("'$Mapping_strain'").alleles' \
   -o ${OUTPUT_NAME}_filter.vcf

#4.Annotate Variants
printf "${col} ##### Annotate ##### ${NC} \n"
java -jar $snpEff/snpEff.jar $snpEff_database \
	-c $snpEff_config \
	-v ${OUTPUT_NAME}_filter.vcf > ${OUTPUT_NAME}_filter.snpeff.vcf
	
#5.Ouputs a small summary
printf "${col} ##### generate a small summary of the best candidates found (based on snpEff annotation) ##### ${NC} \n"
printf "##### RESULTS FOR gene_find_Portable.sh ##### \n" > ${OUTPUT_NAME}_summary_GeneFinding.txt
echo "#`date`"  >> ${OUTPUT_NAME}_summary_GeneFinding.txt

printf "# sample analyzed = $OUTPUT_NAME \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt

printf "# Variant analysis is restricted to intervals defined in `python -c 'import os; print(os.path.abspath("../${OUTPUT_NAME}_plots/${OUTPUT_NAME}_CandidateScaff.tsv"))'` \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
NumberOfIntervals=$(wc -l ../${OUTPUT_NAME}_plots/${OUTPUT_NAME}_CandidateScaff.tsv | cut -d " " -f 1)
printf "# Those Intervals defined by mapping contains $NumberOfIntervals scaffolds \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
printf "# 47 scaffolds with no JU170-SNPs were added to the analysis. Those scaffolds are stored in $ANDALOUSIAN_FOLDER/data/$Invariant_Scaff \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt

#Calculate the cumulative length of these intervals:
cat ${OUTPUT_NAME}.interval_list | while read line ; do grep "$line" $ANDALOUSIAN_FOLDER/data/$Scaff_size | cut -f 2 >>temp_size.tsv; done
CumulativeLength=$(awk 'END { print s } { s += $OUTPUT_NAME }' temp_size.tsv)
TotalLength=$(awk 'END { print s } { s += $2 }' $ANDALOUSIAN_FOLDER/data/$Scaff_size)
Fraction=`echo "$CumulativeLength $TotalLength" | awk '{printf "%.2f \n", $1/$2}'`
printf "# Altogether, these intervals are $CumulativeLength bp long ( $Fraction of Oscheius genome) \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
rm temp_size.tsv

#Warning if scaffolds from Warning_Scaff are also present (e.g. "False-Positive" scaffolds):
TEST=()
cat $ANDALOUSIAN_FOLDER/data/$Warning_Scaff | while read line ; do b=$(grep "$line" ${OUTPUT_NAME}.interval_list) ;	TEST+=($b) ; done
if [ ! -z "$TEST" ]; then 
	printf "# Warning: the selection of mapping intervals countains ${#TEST[@]} scaffolds that belong to the warning list $Warning_Scaff. \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
	printf "# Here are those scaffolds:\n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
	echo "#${TEST[@]}" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
else 
	printf "Your selection of mapping intervals does not contain any scaffolds from $Warning_Scaff. \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
fi 

NumberOfRawVariants=$(grep -vc "#" ${OUTPUT_NAME}_raw.vcf)
printf "#Variant analysis on the selected intervals calls $NumberOfRawVariants variants \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
NumberOfSpecificVariants=$(grep -vc "#" ${OUTPUT_NAME}_filter.vcf)
printf "# After filtering, Number of ${OUTPUT_NAME}-specific variants is = $NumberOfSpecificVariants \n \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt

#Sort candidate genes based on the functional annotation
#HIGH IMPACT mutations
printf "\n ## HIGH IMPACT mutations \n" >> ${OUTPUT_NAME}_summary_GeneFinding.txt
cat ${OUTPUT_NAME}_filter.snpeff.vcf | \
$snpEff/scripts/vcfEffOnePerLine.pl | \
java -jar $snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" | \
java -jar $snpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "LOF[*].GENE">> ${OUTPUT_NAME}_summary_GeneFinding.txt

NbHIGH=$(cat ${OUTPUT_NAME}_filter.snpeff.vcf | \
	$snpEff/scripts/vcfEffOnePerLine.pl | \
	java -jar $snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" | grep -vc "#")

#MODERATE IMPACT mutations
printf "\n ## MODERATE IMPACT mutations \n">> ${OUTPUT_NAME}_summary_GeneFinding.txt
cat ${OUTPUT_NAME}_filter.snpeff.vcf | \
$snpEff/scripts/vcfEffOnePerLine.pl | \
java -jar $snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'MODERATE'" | \
java -jar $snpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "LOF[*].GENE">> ${OUTPUT_NAME}_summary_GeneFinding.txt

NbMODERATE=$(cat ${OUTPUT_NAME}_filter.snpeff.vcf | \
$snpEff/scripts/vcfEffOnePerLine.pl | \
java -jar $snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'MODERATE'" | grep -vc "#")

#LOW IMPACT mutations
printf "\n ## LOW IMPACT mutations \n">> ${OUTPUT_NAME}_summary_GeneFinding.txt
cat ${OUTPUT_NAME}_filter.snpeff.vcf | \
$snpEff/scripts/vcfEffOnePerLine.pl | \
java -jar $snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'LOW'" | \
java -jar $snpEff/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "LOF[*].GENE">> ${OUTPUT_NAME}_summary_GeneFinding.txt

NbLOW=$(cat ${OUTPUT_NAME}_filter.snpeff.vcf | \
$snpEff/scripts/vcfEffOnePerLine.pl | \
java -jar $snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'LOW'" | grep -vc "#" )


##Counts of variants per functional class (in the file + STDOUT)
printf "${col} ##Counts of variants per functional class ${NC} \n"
echo "##Counts of variants per functional class">> ${OUTPUT_NAME}_summary_GeneFinding.txt

echo "Number of High-impact mutations: $NbHIGH"
echo "Number of High-impact mutations: $NbHIGH">> ${OUTPUT_NAME}_summary_GeneFinding.txt
echo "Number of Moderate-impact mutations: $NbMODERATE"
echo "Number of Moderate-impact mutations: $NbMODERATE">> ${OUTPUT_NAME}_summary_GeneFinding.txt
echo "Number of Low-impact mutations: $NbLOW"
echo "Number of Low-impact mutations: $NbLOW">> ${OUTPUT_NAME}_summary_GeneFinding.txt

echo "##end of file" >> ${OUTPUT_NAME}_summary_GeneFinding.txt

printf "${col} END OF SCRIPT 'find_gene_Portable.sh' ${NC} \n"
