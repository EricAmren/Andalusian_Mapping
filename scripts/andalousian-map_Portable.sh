#!/bin/bash

## script: 'andalousian-map_Portable.sh'
## Fabrice Besnard, 2017-02-27, v1 
# Requires a config file
# Derived from original script: andalousian-map.v2_PSMN.sh (2016-06-27)
# last edit: 2017-02-27

#What it does: 
# 1)map with bwa 
#		mem faster and more accurate algorithm (ideal for small genome and "long" 100-pb Illumina reads) ;
#		-a Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments ; 
#		-M=shorter split reads flagged as secondary (mem allows split alignment), necessary for PICARD compatibility ;
#		-can take several sequencing runs for the same sample.
#		Note: it does not filter the resulting BAM using FLAG: GATK-tools' filters will do it
#			  it is not necessary to use Picard/FixMateInformation or samtools/fixmate beacuse GATK/indelrealigner does it
# 2) do the preliminary procedures of GATK best practices before accurate calling 
#		->File formatting/compatibility: Add Read Groups, mark duplicated Reads, Index
#		->Sequencing data & Alignment improvement: Realignment around Indel + BQSR (given bootstrap and dbSNP_JU170)
#3) do a variant call with GATK Haplotype caller in a normal mode, RESTRICTED TO JU170 SNPs only.
#4) Use map_mutant_PSMN.sh
#		->From previous .vcf file containing 1 sample only, it computes for each reported variant position the frequency of reads for the ALT allele over the sum of (REF + ALT) read counts.
#		->Results are stored in a file named $argument2.freq.txt
#		->Then plots of frequency and coverage are geneated with a Rcript and stored in a new folder named "plots".
#5) Use find_gene_PSMN.sh (Here: tested with find_gene.v3_PSMN.sh)
#		->Calls the variants of the mutant only in the scaffolds where the mutation has been mapped
#		->Filters out all JU170 variations + quality filters
#		->Annotate the variations with snpeff

#### Required: 
##Softwares:
# bwa version: 0.7.5a-r405+
# samtools v0.1.18+
# Picard Version: 1.110+.
# GATK 3.7+ (java1.8 compatible)
# R (v3.3.0+)
# Snpeff (v4.1g+). 
##Note related to the use of snpEff: managing databases and snpEff.config
#If you are working with an organism (or assembly) which is not in a pre-built database released with snpEff:
#		#1/ Build the database according to http://snpeff.sourceforge.net/SnpEff_manual.html#databases +
#		#2/ Store this database in andalousian-mapping/data/snpEff-Databases in a dedicated folder with the same name + 
#		#3/ Place the corresponding genome.fa file (with the same name) in andalousian-mapping/data/snpEff-Databases/genomes the "genomes" folder +
#		#4/ Add the database entry to snpEff_andalousian.config
#See this example with the O.tipulae assembly n0t.2.0.
#The database version and snpEff.jar may need to be of the same version. Current O.tipulae.2.0 has been built with snpEff.v4.3. 

##Scripts:
# map_mutant_Portable.sh
# find_gene_Portable.sh
# Plot-fqcy_Portable.R

##Data
# The reference genome of the organism in FASTA format. 
# Corresponding bwas index + fasta file index + sequence dictionary of this reference genome, in the same folder. To generate those files, see http://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
# a vcf file containing the list of SNPs between the background and mapping strains (e.g. dbSNP_JU170_v2.0.vcf)
# a tab file containing the scaf names and their sizes in two columns, respectively (e.g. nOt.2.0_Scaff-size.tsv)
# a file containing the "Invariant" Scaffolds (ie containing no snps) between mapping/background strain
#	-e.g.: 'JU170_InvariantScaffold.interval_list'. 
# gVCF files of background strain (e.g. CEW1) and of the mapping strain (e.g. JU170)
#	-CEW1.gVCF.vcf
#	-JU170.gVCF.vcf

##Config file. a simple txt file containing variables.
#See in the folder andalousian-mapping the guides to format the config file:
#A detailed explanation: Config-file_Portable_TEMPLATE.txt
#An example: Config-file_Portable_MyLinux.txt

### Usage:
#bash andalousian-map_Portable.sh Config-file_mysample.txt

#From where it is started, it creates a subfolder named "mysample" to store output files if it does not exist already. 
#Suggestion: create an 'analysis' folder, go to it and start the analysis for all your samples from there.


###Start of script
CONFIGFILE=$1
source $CONFIGFILE

#Inform users of the parameters:
printf "${col} **Running andalousian-map_Portable.sh** ${NC} \n"
printf "${col}Reference genome selected: ${NC} $REFERENCE \n"
printf "${col}Numbers of sequencing runs used for mapping: ${NC} $NbRUNS \n"
for ((a=1; a <=$NbRUNS; a++))
	do
		echo "Sequencing Run $a, Fwd pair is: ${READS[$((2*$a-1))]}"
		echo "Sequencing Run $a, Reverse pair is: ${READS[$((2*$a))]}"
		echo "Sequencing Run $a, RGID is: ${RGID[$a]}"
		echo "Sequencing Run $a, RGPU is: ${RGPU[$a]}"
	done 
printf "${col}dbSNP selected:${NC} $dbSNP  \n"
printf "${col}Final file generated will be:${NC} $OUTPUT_NAME.vsMapStrain.vcf, where 'sample' is $RGSM \n"

#Create a subfolder if necessary
if [ ! -e $OUTPUT_NAME ]
then 
	mkdir $OUTPUT_NAME
else
	echo "a folder called $OUTPUT_NAME already exists. Output data will be stored in it."
fi

cd $OUTPUT_NAME

#STEP1. MAPPING
printf "${col} **STEP1** MAPPING pe-reads with bwa (mem, -aM, default) + sort bam output ${NC} \n"
#Note: bwa allowed for 4 threads. samtools sort syntax (with -T) fixed in new samtools version.
for ((a=1; a <=$NbRUNS; a++))
do
	echo "           ------1.$a Map Sequencing Run n°$a ------"
	bwa mem -t 4 -aM $REFERENCE ${READS[$((2*$a-1))]} ${READS[$((2*$a))]} | samtools view -buS -|samtools sort - ${OUTPUT_NAME}.SR$a
done

#Step2. READ GROUPS
printf "${col} **STEP2** Adding Read-Group Informations ${NC} \n"

for ((a=1; a <=$NbRUNS; a++))
do
	echo "           ------2.$a Read Groups for Sequencing Run n°$a ------"
	java -Xmx${RAM}g -jar $PICARD AddOrReplaceReadGroups \
	I=$OUTPUT_NAME.SR$a.bam \
	RGID=${RGID[$a]} RGSM=$RGSM RGLB=$RGLB RGPU=${RGPU[$a]} RGPL=$RGPL \
	O=$OUTPUT_NAME.RG.SR$a.bam
	printf "${col} Check & Clean files before next step ${NC} \n"
	if [ -e $OUTPUT_NAME.RG.SR$a.bam ]; then
	rm $OUTPUT_NAME.SR$a.bam
	else     
	printf "\e[0;31m ##ERROR## Expected bam file was not generated. Script will abort ${NC} \n"    
	exit 1
	fi
done

#Step3. Merge and Mark Duplicated Reads at once
printf "${col}**STEP3** Merge all Sequencing Runs and Mark duplicated reads ${NC} \n"
mkdir metrics

#This function allows to feed the MarkDuplicate tool of Picard with all bam previously generated:
function multipleDupe () { #$1 is the number of sequencing Runs=Nber of Bam inputs
	declare -A INPUTS
	for ((a=1; a <=$1; a++))
	do
	INPUTS[$a]="INPUT=$OUTPUT_NAME.RG.SR$a.bam"  #Note: this function only works if $OUTPUT_NAME is defined in the shell env
	done
	java -jar $PICARD MarkDuplicates \
	`echo "${INPUTS[@]}"`\
	OUTPUT=$OUTPUT_NAME.RG.dedup.SRmerged.bam \
	METRICS_FILE=./metrics/dedup.$OUTPUT_NAME.SRmerged.txt
}

multipleDupe $NbRUNS

printf "${col} Check & Clean files before next step ${NC} \n"
if [ -e $OUTPUT_NAME.RG.dedup.SRmerged.bam ]
then
	for ((a=1; a <=$NbRUNS; a++))
	do
		rm $OUTPUT_NAME.RG.SR$a.bam #Delete intermediate .bam files for separate SeqRuns.
	done	    
else
	printf "\e[0;31m ##ERROR##Expected merged+dedup bam file was not generated. Script will abort ${NC} \n"    
	exit 1
fi

#Step4. Index bam
printf "${col}**STEP4** Indexing BAM${NC} \n"
java -Xmx${RAM}g -jar $PICARD BuildBamIndex \
INPUT=$OUTPUT_NAME.RG.dedup.SRmerged.bam

#Step5. Get basics statistics
printf "${col} **STEP5** Generating basic stats on bam file ${NC} \n"
samtools flagstat $OUTPUT_NAME.RG.dedup.SRmerged.bam > ./metrics/$OUTPUT_NAME.RG.dedup.SRmerged.flagstat.txt

#Step6. Realign around Indels
printf "${col} **STEP6** Realign reads around Indel$ {NC} \n"
printf "${col}            ------6.1. Creating interval table------- ${NC} \n"
java -Xmx${RAM}g -jar $GATKHOME/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R $REFERENCE \
-I $OUTPUT_NAME.RG.dedup.SRmerged.bam \
-o ./metrics/indelrealigner.$OUTPUT_NAME.intervals
printf "${col}            ------6.2. Realigning Reads in bam------- ${NC} \n"
java -Xmx${RAM}g -jar $GATKHOME/GenomeAnalysisTK.jar -T IndelRealigner \
-R $REFERENCE \
-I $OUTPUT_NAME.RG.dedup.SRmerged.bam -o $OUTPUT_NAME.RG.dedup.SRmerged.realign.bam \
-targetIntervals ./metrics/indelrealigner.$OUTPUT_NAME.intervals --filter_bases_not_stored
printf "${col} Check & Clean files before next step ${NC} \n"
if [ -e $OUTPUT_NAME.RG.dedup.SRmerged.realign.bam ]; then    
rm $OUTPUT_NAME.RG.dedup.SRmerged.ba*
else printf "\e[0;31m ##ERROR## Expected realigned bam file was not generated. Script will abort ${NC} \n"    
exit 1
fi

#Step7. BQSR (With HC as caller, only one step is necessary when bootstrapping a first call)
printf "${col} **STEP7** BQSR ${NC} \n"
#Create a folder to contain files related to the BSQR
mkdir BQSR_files
printf "${col}            ------7.1. Perform a first variant call (HC)------- ${NC} \n"
java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $REFERENCE \
-stand_call_conf 10 \
-I $OUTPUT_NAME.RG.dedup.SRmerged.realign.bam  -o ./BQSR_files/$OUTPUT_NAME.HC0.vcf
printf "${col}            ------7.2. Analyze BaseQuality Covariates by bootstrapping the previous vcf + giving $dbSNP------- ${NC} \n"
java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T BaseRecalibrator \
-R $REFERENCE \
-I $OUTPUT_NAME.RG.dedup.SRmerged.realign.bam \
-knownSites ./BQSR_files/$OUTPUT_NAME.HC0.vcf \
-knownSites $ANDALOUSIAN_FOLDER/data/$dbSNP \
-o ./BQSR_files/$OUTPUT_NAME.BQSR.table
printf "${col}            ------7.3. Analyze Remaining covariates after BQSR------- ${NC} \n"
java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T BaseRecalibrator \
-R $REFERENCE \
-I $OUTPUT_NAME.RG.dedup.SRmerged.realign.bam -knownSites ./BQSR_files/$OUTPUT_NAME.HC0.vcf \
-BQSR ./BQSR_files/$OUTPUT_NAME.BQSR.table \
-o ./BQSR_files/$OUTPUT_NAME.post-BQSR.table
printf "${col}            ------7.4. Generate plots and stats of the BQSR ------- ${NC} \n"
java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T AnalyzeCovariates \
-R $REFERENCE \
-before ./BQSR_files/$OUTPUT_NAME.BQSR.table -after ./BQSR_files/$OUTPUT_NAME.post-BQSR.table \
-plots ./BQSR_files/$OUTPUT_NAME.BQSRplots.pdf
printf "${col}            ------7.5. Apply BQSR to the bam ------- ${NC} \n"
java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T PrintReads \
-R $REFERENCE \
-I $OUTPUT_NAME.RG.dedup.SRmerged.realign.bam -BQSR ./BQSR_files/$OUTPUT_NAME.BQSR.table \
-o $OUTPUT_NAME.RG.dedup.SRmerged.realign.BQSR.bam
printf "${col} Check & Clean files before next step ${NC} \n"
if [ -e $OUTPUT_NAME.RG.dedup.SRmerged.realign.BQSR.bam ]; then    
rm $OUTPUT_NAME.RG.dedup.SRmerged.realign.ba*
else     
printf "\e[0;31m ##ERROR## Expected BQSR-bam file was not generated. Script will abort ${NC} \n"
exit 1
fi

#Step8. Call (With HC as caller)
printf "${col} **STEP8** Variant Calling --Genotype $OUTPUT_NAME for $Mapping_strain positions-- ${NC} \n"
java -jar -Xmx${RAM}g  $GATKHOME/GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $REFERENCE \
-I $OUTPUT_NAME.RG.dedup.SRmerged.realign.BQSR.bam \
--genotyping_mode GENOTYPE_GIVEN_ALLELES \
--alleles $ANDALOUSIAN_FOLDER/data/$dbSNP \
-L $ANDALOUSIAN_FOLDER/data/$dbSNP \
--dbsnp $ANDALOUSIAN_FOLDER/data/$dbSNP \
-o $OUTPUT_NAME.vsMapStrain.vcf

if [ -e $OUTPUT_NAME.vsMapStrain.vcf ]; then    
printf "${col} **END OF Variant Calling** Variants of $OUTPUT_NAME called and stored into $OUTPUT_NAME.vsMapStrain.vcf ${NC} \n"
else     
printf "\e[0;31m ##ERROR## Expected VCF file was not generated. Last step may have failed ${NC} \n"
exit 1
fi

#Step 9. Mapping of the mutation
printf "${col} **STEP9** Mapping the mutation --Compute $Mapping_strain allele frequency for each scaffold and make plots in $OUTPUT_NAME dataset-- ${NC}\n"
bash $ANDALOUSIAN_FOLDER/scripts/map_mutant_Portable.sh $OUTPUT_NAME.vsMapStrain.vcf $OUTPUT_NAME $1

#Step 10 Finding the best condidate genes
printf "${col} **STEP10** find the best condidate genes thanks to annotations ${NC}\n"
bash $ANDALOUSIAN_FOLDER/scripts/find_gene_Portable.sh $1

