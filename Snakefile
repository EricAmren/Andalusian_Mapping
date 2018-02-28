import vcf
import subprocess
import os
import re

nbRun = 1
readID = 'mf76'
PATH = '/home/amren/Master1/Andalusian_Mapping/'
SAMPLES = ['/home/amren/Master1/Andalusian_Mapping/data/sample/mf76_1_sub.fastq.gz', '/home/amren/Master1/Andalusian_Mapping/data/sample/mf76_2_sub.fastq.gz']
genome = PATH + 'data/Reference_genomes/Otipulae/nOt.2.0/nOt.2.0.fna'

RAM = "8"
RGPL = "ILLUMINA"
RGSM = "mf76"
RGLB = "mf76_LB1"
RGID = "FCH3LFLBBXX_mf76"
RGPU = "FCH3LFLBBXX_L5_mf76"
dbSNP= PATH + "data/dbSNP_JU170_v2.0.vcf"
SCAFF_SIZE= PATH + "data/nOt.2.0_Scaff-size.tsv"
INVARIANT_SCAFF= PATH + "data/JU170_InvariantScaffold.interval_list"
BACKGROUND_GVCF = PATH + "data/CEW1.gVCF.vcf"
MAPPING_GVCF = PATH + "data/JU170.gVCF.vcf"
snpEff_database = "O.tipulae.2.0"
snpEff_contig= PATH + "data/snpEff-Databases/snpEff_andalousian.config"

# Software Path
picard = "/home/amren/soft/picard.jar"
#GATK = "~/soft/GenomeAnalysisTK-3.7-0-gcfedb67/GenomeAnalysisTK.jar"
GATK = "/home/amren/soft/GenomeAnalysisTK.jar"
newGATK = "~/Master1/ProjetS2/Log/gatk-4.0.1.1/gatk"
snpEff = "~/Master1/ProjetS2/Log/snpEff/snpEff.jar"
#snpEff = "~/soft/snpEff/snpEff.jar"


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
        readID + ".RG.dedup.SRmerged.bam",
        readID + ".RG.dedup.SRmerged.bai",
        "metrics/" + readID + ".RG.dedup.SRmerged.flagstat.txt",
        "metrics/indelrealigner." + readID + ".intervals",
        readID + ".RG.dedup.SRmerged.realign.bam",
        "BQSR_files/" + readID +".HC0.vcf",
        "BQSR_files/" + readID + ".BQSR.table",
        "BQSR_files/" + readID + ".post-BQSR.table",
        "BQSR_files/" + readID + ".BQSRplots.pdf",
        readID + ".RG.dedup.SRmerged.realign.BQSR.bam",
        readID + ".vsMapStrain.vcf",
        readID + "_plots/" + readID + ".freq.txt",
        readID + "_plots/" + readID + "_CandidateScaff.tsv",
        readID + "_plots/" + readID + "_summary.tsv",
        PATH + readID + "/" + readID + "_find_gene/" + readID + ".interval_list",
        readID + "_find_gene/" + readID + ".gVCF.vcf",
        #readID + "_find_gene/" + readID + "_raw.vcf",
        PATH + readID + "/" +readID + "_find_gene/" + readID + "_raw.vcf",
        #readID + "_find_gene/" + readID + "_filter.vcf"
        #readID + "_filter.vcf"
        PATH + readID + "/" + readID + "_find_gene/" + readID + "_filter.vcf",
        PATH + readID + "/" + readID + "_find_gene/" + readID + "_filter.snpeff.vcf"


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
        "java -Xmx{RAM}g -jar {picard} BuildBamIndex INPUT={input}"

# STEP 5

rule get_basic_stats:
    input:
        readID + ".RG.dedup.SRmerged.bam"
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
        --alleles {dbSNP} \
        -L {dbSNP} \
        --dbsnp {dbSNP} \
        -o {output}"

# STEP 9 : map_mutant_Portable.sh on vsMapStrain.vcf readID configFile

def build_freq_file(vcf_file, output_file):
    with open(output_file, 'w') as freq_file:
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        for record in vcf_reader:
            for sample in record.samples:
                line = ""
                GT = sample['GT']
                line += record.CHROM + ' ' + str(record.POS) + ' '
                line += record.REF[0] + ' ' + str(record.ALT[0]) + ' '
                line += GT
                try:
                    DP = sample['DP']
                    line += ' ' + str(DP)
                except AttributeError:
                    line += ' 0' 
                try: 
                    AD = sample['AD']
                    sept = str(int(AD[0]) + int(AD[1]))
                    huit = str(int(AD[1])/(int(AD[0]) + int(AD[1])))
                    line += ' ' + sept + ' ' + huit
                except AttributeError:
                    line += " 0 0.0"
                line += '\n'
                freq_file.write(line)

rule get_frequencies:
    input:
        vcf = readID + ".vsMapStrain.vcf"
    output:
        txt = readID + "_plots/" + readID + ".freq.txt"
    run:
        build_freq_file(input.vcf, output.txt)

rule generating_plot:
    input:
        txt = readID + "_plots/" + readID + ".freq.txt"
    output:
        readID + "_plots/" + readID + "_CandidateScaff.tsv",
        readID + "_plots/" + readID + "_summary.tsv"

    shell:
        "Rscript ../scripts/Plot-fqcy_Portable.R {input.txt} {readID}_plots/{readID}"

# STEP 10 : find_gene_Portable configFile

rule build_candidate_scaffs_list:
    input:
        candidate_scaff = readID + "_plots/" + readID + "_CandidateScaff.tsv"
    output:
        #interval_list = readID + "_find_gene/" + readID + ".interval_list"
        interval_list = PATH + readID + "/" + readID + "_find_gene/" + readID + ".interval_list"
    run:
        shell("cut -f1 {input} > {output}")
        shell("cat {INVARIANT_SCAFF} >> {output}")
        shell("sort -o {output} {output}")

rule variant_call:
    input:
        bam = readID + ".RG.dedup.SRmerged.realign.BQSR.bam",
        interval_list = PATH + readID + "/" + readID + "_find_gene/" + readID + ".interval_list"
        #interval_list = readID + "_find_gene/" + readID + ".interval_list"
    output:
        readID + "_find_gene/" + readID + ".gVCF.vcf"
    shell:
        "java -jar -Xmx{RAM}g {GATK} \
        -T HaplotypeCaller \
        -R {genome} \
        --emitRefConfidence GVCF \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        -L {input.interval_list} \
        -I {input.bam} \
        --dbsnp {dbSNP} \
        -o {output}"

rule genotype_GVCF:
    input:
        readID + "_find_gene/" + readID + ".gVCF.vcf"
    output:
        #readID + "_find_gene/" + readID + "_raw.vcf"
        vcf = PATH + readID + "/" +readID + "_find_gene/" + readID + "_raw.vcf",
    shell:
        "java -jar -Xmx{RAM}g {GATK} \
        -T GenotypeGVCFs \
        -R {genome} \
        --variant {BACKGROUND_GVCF} \
        --variant {MAPPING_GVCF} \
        --variant {input} \
        -o {output}"

rule filter_variants:
    input:
        vcf = PATH + readID + "/" +readID + "_find_gene/" + readID + "_raw.vcf",
        interval_list = PATH + readID + "/" + readID + "_find_gene/" + readID + ".interval_list"
    output:
        PATH + readID + "/" + readID + "_find_gene/" + readID + "_filter.vcf"
    shell:
        """
        cd {readID}_find_gene
        java -jar -Xmx{RAM}g {GATK} \
        -T SelectVariants \
        -R {genome} \
        --variant {input.vcf} -L {input.interval_list} \
        -select 'DP >= 3 && !vc.hasAttribute("DB")' \
        -select 'vc.getGenotype({readID}).isHom() && \
        vc.getGenotype("{BACKGROUND_GVCF}").isHom() && \
        vc.getGenotype("{MAPPING_GVCF}").isHom()' \
        -select 'vc.getGenotype({readID}).alleles != \
        vc.getGenotype("{BACKGROUND_GVCF}").alleles && \
        vc.getGenotype({readID}).alleles != vc.getGenotype("{MAPPING_GVCF}").alleles' \
        -o mf76_filter.vcf
        cd ..
        """

rule annotate_variants:
    input:
        PATH + readID + "/" + readID + "_find_gene/" + readID + "_filter.vcf"
    output:
        PATH + readID + "/" + readID + "_find_gene/" + readID + "_filter.snpeff.vcf"
    shell:
        """
        cd mf76_find_gene
        java -jar {snpEff} {snpEff_database} \
        -c {snpEff_contig} \
        -v {input} > {output}
        cd ..
        """
