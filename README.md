# Andalousian_Mapping
Pipeline to perform in a single step mapping-by-sequencing and mutation identification.

##General description
This pipeline has been originally developped for the nematode species *Oscheius tipulae* by the [Félix lab](http://www.ibens.ens.fr/?rubrique29).
This repository contains all data to perform mapping-by-sequencing with Oscheius tipulae using the reference strain CEW1 and the mapping strain JU170 (see Félix Lab [Strain database](http://www.justbio.com/worms/index.php)).
However, this is a versatile pipeline that allows working with any species and any mapping strains (See "How to adapt it to my purpose ?"). 

**Andalousian_Mapping** is mainly encoded in bash and it can be run on Linux (tested on Ubuntu 16.04 LTS) and MacOS (not Windows). It uses popular third-party NGS softwares (six exactly, see "Requirements") that must be installed in the computing machine.

Many bio-informatic steps are needed from the raw sequencing data to the final identification of candidate mutations. **Andalousian_Mapping** allows pipelining them all via one simple command line. This pipeline can then be used routinely with very basic knowledge in shell coding.

##Why using **Andalousian_Mapping** ?
Recently, [Cloudmap](https://usegalaxy.org/u/gm2123/p/cloudmap), a nice automatic pipeline of mapping-by-sequencing/mutant-identification has been developped and integrated to Galaxy, with an intuitive and user-friendly GUI. It works very well with model species like *Caenorhabditis elegans* or *Arabidopsis thaliana*.
However, all genomes are not available in Galaxy. We had to develop Anadalousian-Mapping to work with the non-model species Oscheius tipulae, which had no published genome. But we tried to keep simple procedures so that a large audiance of biologists (and not only bio-informaticians) could find it useful. 

We propose that **Andalousian-Mapping** will be an attractive alternative to other pipelines for people who need to work with non-reference genome assembly but with minimal efforts in computing/bio-informatics. Once the set-up is done, the procedure is really simple: one command line only.

Here are some advantages of **Andalousian-Mapping** that may match your needs:
-Most third-party software easy to install and so popular that they might already be on your system (e.g. R, bwa, samtools) ;
-The most complicated software to install are GATK and snpEff but they both have very clear tutorials and debug forum online. If it remains a problem for you, any people with intermediate computing skills should easily and quickly debug this step for you ;
-Once the set-up is done, the procedure is really simple: one command line only ;
-You can run it locally and don't need to queue your job in a distant server ;
-You can adapt it to any reference genomes. Hence, you can study non-model organisms or use an alternative assembly of a model-organism (different strain depedning on the design of your mapping cross).
-You manage the versions of your software: you can either benefit from the latest improvements of third-party software or freeze you set-up in the configuration you prefer.

##What it does:
**Andalousian-Mapping** takes as input paired-end reads from a bulk sequencing of a F2 population generated through a cross of a (ideally recessive) mutant strain with a mapping strain.
-> It outputs plots of allele frequency (mutant strain allele / mapping-strain allele) for each scaffolds/Chromosome present in the reference genome. This allows to visualy map the region linked to the mutation selected in the F2 population.
-> It outputs an annotated vcf file containing the mutations linked to the phenotype-causing mutations with their functional impacts
-> It outputs a text file that select the best candidate mutations based on the predicted functional impact.

##Requirements
###Software
* [bwa](http://bio-bwa.sourceforge.net/) version 0.7.5a-r405 or later
* [samtools](http://samtools.sourceforge.net/) version 0.1.18  or later
* [Picard](https://broadinstitute.github.io/picard/) Version 1.110  or later
* [GATK](https://software.broadinstitute.org/gatk/) Version 3.7  or later
* [R](https://www.r-project.org/). Version 3.3  or later
* [snpEff](http://snpeff.sourceforge.net/). Version 4.1g or later

##Usage
Run a simple stereotyped command line: 

`bash andalousian-map_Portable.sh Config-file_mysample.txt`

For each new analysis, you only change the configuration file (Config-file_mysample.txt), which is a simple texte file.

##How to adapt it to my purpose ?

###Notes:
We named our pipeline "Andalousian_Mapping" in reference to the "Hawaïan Mapping" procedure developped in C.elegans, involving the Hawaïan mapping strain CB4856. Our Oscheius tipuale mapping strain is a wild isolate sampled in Sevilla, Andalousia.
