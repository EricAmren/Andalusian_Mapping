import subprocess
import os

nbRun = 1
readID = 'mf76'
#SAMPLE = ['mf76_1_sub', 'mf76_2_sub']
PATH = '/home/amren/Master1/Andalusian_Mapping/'
SAMPLES = ['/home/amren/Master1/Andalusian_Mapping/data/sample/mf76_1_sub.fastq.gz', '/home/amren/Master1/Andalusian_Mapping/data/sample/mf76_2_sub.fastq.gz']
#SAMPLES = ['/home/amren/Master1/Andalusian_Mapping/data/Full_data/mf76_1.fq.gz', '/home/amren/Master1/Andalusian_Mapping/data/Full_data/mf76_2.fq.gz']
genome = PATH + 'data/Reference_genomes/Otipulae/nOt.2.0/nOt.2.0.fna'

RAM = "4"
RGPL = "ILLUMINA"
RGSM = "mf76"
RGLB = "mf76_LB1"
RGID = "FCH3LFLBBXX_mf76"
RGPU = "FCH3LFLBBXX_L5_mf76"
dbSNP=PATH + "data/dbSNP_JU170_v2.0.vcf"

# Software Path
picard = "/home/amren/soft/picard.jar"
#picard = "/home/amren/soft/picard/build/libs/picard.jar"
#GATK = "~/soft/GenomeAnalysisTK.jar"
GATK = "~/soft/GenomeAnalysisTK-3.7-0-gcfedb67/GenomeAnalysisTK.jar"

# Main
os.chdir(PATH)
if not os.path.isdir(readID):
    os.mkdir(readID)
os.chdir(PATH + readID)

rule all:
    input:
        readID + ".SR1.bam",
        readID + ".SR1.fixed.bam",
        readID + ".RG.SR1.bam",
        readID + ".RG.dedup.SRmerged.bam",
        "metrics/dedup." + readID + ".SRmerged.txt",
        readID + ".RG.dedup.SRmerged.bai",
        "metrics/" + readID + ".RG.dedup.SRmerged.flagstat.txt",
        "metrics/indelrealigner." + readID + ".intervals",
        readID + ".RG.dedup.SRmerged.realign.bam",
        "BQSR_files/" + readID +".HC0.vcf",
        "BQSR_files/" + readID + ".BQSR.table",
        "BQSR_files/" + readID + ".post-BQSR.table",
        "BQSR_files/" + readID + ".BQSRplots.pdf",
        readID + ".RG.dedup.SRmerged.realign.BQSR.bam",
        readID + ".vsMapStrain.vcf"

rule mapping:
    input:
        genome,
        expand('{sample}', sample=SAMPLES)
    output:
        readID + ".SR1.bam"
    shell:
        "bwa mem -t 4 -aM {input} |\
        samtools view -buS -|\
        samtools sort - {readID}.SR1"

rule fix_mate_info:
    input:
        readID + ".SR1.bam"
    output:
        readID + ".SR1.fixed.bam"
    shell:
        "java -Xmx{RAM}g -jar {picard} FixMateInformation \
        I={input} \
        O={output} \
        SO=coordinate \
        VALIDATION_STRINGENCY=LENIENT"

# STEP 2
rule add_read_group_info:
    input:
        readID + ".SR1.fixed.bam",
    output:
        readID + ".RG.SR1.bam"
    shell:
        "java -Xmx{RAM}g -jar {picard} AddOrReplaceReadGroups \
        I={input} \
        RGID={RGID} \
        RGSM={RGSM} \
        RGLB={RGLB} \
        RGPU={RGPU} \
        RGPL={RGPL} \
        O={output}"

# STEP 3

if not os.path.isdir("metrics"):
    os.mkdir("metrics")

rule mark_duplicates:
    input:
        readID + ".RG.SR1.bam"
    output:
        merged_bam = readID + ".RG.dedup.SRmerged.bam",
        dedup = "metrics/dedup." + readID + ".SRmerged.txt"
    shell:
        "java -jar {picard} MarkDuplicates \
        INPUT={input} \
        OUTPUT={output.merged_bam} \
        METRICS_FILE={output.dedup}"

# STEP 4

rule build_bam_index:
    input:
        readID + ".RG.dedup.SRmerged.bam"
    output:
        readID + ".RG.dedup.SRmerged.bai"
    shell:
        "java -Xmx{RAM}g -jar {picard} BuildBamIndex I={input}"

# STEP 5

rule get_basic_stats:
    input:
        readID + ".RG.dedup.SRmerged.bai"
    output:
        "metrics/" + readID + ".RG.dedup.SRmerged.flagstat.txt"
    shell:
        "samtools flagstat {input} > {output}"

# STEP 6

rule creating_interval_table:
    input:
        bam_file = readID + ".RG.dedup.SRmerged.bam"
    output:
        "metrics/indelrealigner." + readID + ".intervals"
    shell:
        "java -Xmx{RAM}g -jar {GATK} \
        -T RealignerTargetCreator \
        -R {genome} \
        -I {input.bam_file} \
        -o {output}"

rule realign_reads:
    input:
        merged_bam = readID + ".RG.dedup.SRmerged.bam",
        intervals = "metrics/indelrealigner." + readID + ".intervals"
    output:
        readID + ".RG.dedup.SRmerged.realign.bam"
    shell:
        "java -Xmx{RAM}g -jar {GATK} \
        -T IndelRealigner -R {genome} \
        -I {input.merged_bam} -o {output} \
        -targetIntervals {input.intervals} \
        --filter_bases_not_stored"

# STEP 7

if not os.path.isdir("BQSR_files"):
    os.mkdir("BQSR_files")

rule variantCalling:
    input:
        realigned_bam = readID + ".RG.dedup.SRmerged.realign.bam"
    output:
        "BQSR_files/{readID}.HC0.vcf"
    shell:
        "java -jar -Xmx{RAM}g {GATK} \
        -T HaplotypeCaller \
        -R {genome} \
        -stand_call_conf 10 \
        -I {input} \
        -o {output}"

rule analyse_BaseQuality:
    input:
        realigned_bam = readID + ".RG.dedup.SRmerged.realign.bam",
        vcf_file = "BQSR_files/"+ readID + ".HC0.vcf"
    output:
        "BQSR_files/" + readID + ".BQSR.table"
    shell:
        "java -jar -Xmx{RAM}g {GATK} \
        -T BaseRecalibrator \
        -R {genome} \
        -I {input.realigned_bam} \
        -knownSites {input.vcf_file} \
        -knownSites {dbSNP} \
        -o {output}"

rule analyse_remaining_covariates:
    input:
        realigned_bam = readID + ".RG.dedup.SRmerged.realign.bam",
        vcf_file = "BQSR_files/"+ readID + ".HC0.vcf"
    output:
        "BQSR_files/" + readID + ".post-BQSR.table"
    shell:
        "java -jar -Xmx{RAM}g {GATK} \
        -T BaseRecalibrator \
        -R {genome} \
        -I {input.realigned_bam} \
        -knownSites {input.vcf_file} \
        -o {output}"

rule generate_covariates_plots:
    input:
        before = "BQSR_files/" + readID + ".BQSR.table",
        after = "BQSR_files/" + readID + ".post-BQSR.table"
    output:
        "BQSR_files/" + readID + ".BQSRplots.pdf"
    shell:
        "java -jar -Xmx{RAM}g {GATK} \
        -T AnalyzeCovariates \
        -R {genome} \
        -before {input.before} \
        -after {input.after} \
        -plots {output}"

rule apply_BQSR_to_BAM:
    input:
        realigned_bam = readID + ".RG.dedup.SRmerged.realign.bam",
        BQSR_table = "BQSR_files/" + readID + ".BQSR.table"
    output:
        readID + ".RG.dedup.SRmerged.realign.BQSR.bam"
    shell:
        "java -jar -Xmx{RAM}g {GATK} -T PrintReads \
        -R {genome} \
        -I {input.realigned_bam} \
        -BQSR {input.BQSR_table} \
        -o {output}"

# TODO : REMOVE FILES

# STEP 8:

rule variant_calling_for_strain_positions:
    input:
        readID + ".RG.dedup.SRmerged.realign.BQSR.bam"
    output:
        readID + ".vsMapStrain.vcf"
    shell:
        "java -jar -Xmx{RAM}g {GATK} -T HaplotypeCaller \
        -R {genome} \
        -I {input} \
        --genotyping_mode GENOTYPE_GIVEN_ALLELES \
        -L {dbSNP} \
        --dbsnp {dbSNP} \
        -o {output}"

# STEP 9 : map_mutant_Portable.sh on vsMapStrain.vcf readID configFile


# STEP 10 : find_gene_Portable configFile


