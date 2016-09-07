#!/bin/bash 

display_usage() {
  echo -e ""
  echo -e "This script runs the Alignmet, sort and addRG steps of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir split_prefix RGID RGLB RGSM ncore mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  split_prefix: input file prefix (samplename + split_tag) (eg. TEST.MG01HX08_123_HK3GMCCXX_1). Input file names are assumed to be split_prefix.R[12].fastq.gz."
  echo -e "  RGID: lane name (eg. L1)."
  echo -e "  RGLB: lane ID (eg. MG01HX08_123_HK3GMCCXX_1)."
  echo -e "  RGSM: sample ID (eg. S1)."
  echo -e "  ncore: number of CPUs (recommended: 4)."
  echo -e "  mem: memory (recommended: 2G)."
  echo -e ""
}

# if no arguments supplied, display usage 
if [  $# -le 0 ] 
then
  display_usage
  exit 0
fi

# check whether user had supplied -h or --help . If yes display usage 
if [[ $# == 1 && ( $1 == "--help" ||  $1 == "-h" )]] 
then
  display_usage
  exit 0
fi

####################

project_indir=$1
project_outdir=$2
split_prefix=$3  # sample+split_tag
RGID=$4
RGLB=$5
RGSM=$6
ncore=$7
mem=$8

if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi




RGPU=unit1
fastq_extension=fq  

### step2 BWA alignment
/usr/local/bin/bwa-0.7.8/bin/bwa mem -M -t $ncore /resources/human_g1k_v37_decoy.fasta $project_indir/$split_prefix.R1.$fastq_extension $project_indir/$split_prefix.R2.$fastq_extension | /usr/local/bin/samtools-1.3/bin/samtools view -bSho $project_outdir/$split_prefix.bwa.bam - 

### step3 sort SAM
/usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -Djava.io.tmpdir=/tmp -jar /usr/local/bin/picard-1.130/picard.jar SortSam SORT_ORDER=coordinate INPUT=$project_outdir/$split_prefix.bwa.bam OUTPUT=$project_outdir/$split_prefix.sort.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

### step4 ADD read group
/usr/local/bin/jdk1.7.0_71/bin/java -jar -Xmx$mem -Xms$mem /usr/local/bin/picard-1.130/picard.jar AddOrReplaceReadGroups RGID=$RGID RGLB=$RGLB RGPL=illumina RGPU=$RGPU RGSM=$RGSM I=$project_outdir/$split_prefix.sort.bam O=$project_outdir/$split_prefix.addrg.bam 
/usr/local/bin/samtools-1.3/bin/samtools index $project_outdir/$split_prefix.addrg.bam 

