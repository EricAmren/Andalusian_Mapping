import subprocess
import os

nbRun = 1
readID = 'mf76'
SAMPLE = ['mf76_1_sub', 'mf76_2_sub']
PATH = '/home/amren/Master1/Andalusian_Mapping/'
SAMPLES = ['/home/amren/Master1/Andalusian_Mapping/data/sample/mf76_1_sub.fastq.gz', '/home/amren/Master1/Andalusian_Mapping/data/sample/mf76_2_sub.fastq.gz']
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
GATK = "~/soft/GenomeAnalysisTK.jar"

# Main
os.chdir(PATH)
#os.chdir("/home/amren/Master1/Andalusian_Mapping/test")
if not os.path.isdir(readID):
    os.mkdir(readID)
os.chdir(PATH + readID)

rule all:
    input:
        readID + ".SR1.bam",
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
        "BQSR_files/" + readID + ".BQSRplots.pdf"

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

# STEP 2
rule add_group_info:
    input:
        readID + ".SR1.bam",
    output:
        readID + ".RG.SR1.bam"
    shell:
        "java -Xmx{RAM}g -jar {picard} AddOrReplaceReadGroups \
        I={input} RGID={RGID} RGSM={RGSM} RGLB={RGLB} RGPU={RGPU} RGPL={RGPL} \
        O={readID}.RG.SR1.bam"

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
        "java -jar {picard} MarkDuplicates INPUT={input} \
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
        genome,
        bam_file = readID + ".RG.dedup.SRmerged.bam"
    output:
        "metrics/indelrealigner." + readID + ".intervals"
    shell:
        "java -Xmx{RAM}g -jar {GATK} -T RealignerTargetCreator \
        -R {genome} -I {input.bam_file} -o {output}"

rule realign_reads:
    input:
        genome,
        merged_bam = readID + ".RG.dedup.SRmerged.bam",
        intervals = "metrics/indelrealigner." + readID + ".intervals"
    output:
        readID + ".RG.dedup.SRmerged.realign.bam"
    shell:
        "java -Xmx{RAM}g -jar {GATK} -T IndelRealigner -R {genome} \
        -I {input.merged_bam} -o {output} \
        -targetIntervals {input.intervals} --filter_bases_not_stored"

# STEP 7

if not os.path.isdir("BQSR_files"):
    os.mkdir("BQSR_files")

rule variantCalling:
    input:
        genome,
        realigned_bam = readID + ".RG.dedup.SRmerged.realign.bam"
    output:
        "BQSR_files/{readID}.HC0.vcf"
    shell:
        "java -jar -Xmx{RAM}g {GATK} -T HaplotypeCaller \
        -R {genome} -stand_call_conf 10 -I {input.realigned_bam} \
        -o {output}"

rule analyse_BaseQuality:
    input:
        genome,
        realigned_bam = readID + ".RG.dedup.SRmerged.realign.bam",
        vcf_file = "BQSR_files/"+ readID + ".HC0.vcf"
    output:
        "BQSR_files/" + readID + ".BQSR.table"
    shell:
        "java -jar -Xmx{RAM}g {GATK} -T BaseRecalibrator \
        -R {genome} -I {input.realigned_bam} -knownSites {input.vcf_file} \
        -knownSites {dbSNP} -o {output}"

rule analyse_remaining_covariates:
    input:
        realigned_bam = readID + ".RG.dedup.SRmerged.realign.bam",
        vcf_file = "BQSR_files/"+ readID + ".HC0.vcf"
    output:
        "BQSR_files/" + readID + ".post-BQSR.table"
    shell:
        "java -jar -Xmx{RAM}g {GATK} -T BaseRecalibrator \
        -R {genome} -I {input.realigned_bam} -knownSites {input.vcf_file} \
        -o {output}"

rule generate_plots:
    input:
        before = "BQSR_files/" + readID + ".BQSR.table",
        after = "BQSR_files/" + readID + ".post-BQSR.table"
    output:
        "BQSR_files/" + readID + ".BQSRplots.pdf"
    shell:
        "java -jar -Xmx{RAM}g {GATK} \
        -T AnalyseCovariates \
        -R {genome} \
        -before {input.before} \
        -after {input.after} \
        -plots {output}"
