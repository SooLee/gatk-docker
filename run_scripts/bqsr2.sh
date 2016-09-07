#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the BQSR (2) step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_bqsr_indir project_outdir prefix chr ncore mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_bqsr_indir: directory in which input bqsr files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  prefix: prefix of the input/output bam/bqsr files. The input file names are assumed to be prefix.indel.chr.bam and prefix.bqsr and the output file name will be prefix.final.chr.bam."
  echo -e "  chr: chromosome (eg. 21, or 'decoy', or 'unmapped')."
  echo -e "  ncore: number of CPUs (recommended: 2)"
  echo -e "  mem: memory (recommended: 8G)."
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
project_bqsr_indir=$2
project_outdir=$3
prefix=$4
chr=$5
ncore=$6
mem=$7


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_bqsr_indir ]; then
  echo -e "project input bqsr directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi



### step7 BQSR (2)

if [ $chr == 'decoy' ]; then
  /usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T PrintReads -R /resources/human_g1k_v37_decoy.fasta -BQSR $project_bqsr_indir/$prefix.bqsr -I $project_indir/$prefix.indel.$chr.bam -o $project_outdir/$prefix.final.$chr.bam -nct $ncore -XL /resources/b37_nodecoy.bed
elif [ $chr == 'unmapped' ]; then
  cp $project_indir/$prefix.indel.$chr.bam $project_outdir/$prefix.final.$chr.bam
else  # 1,2,3, ..., X,MT
  /usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T PrintReads -R /resources/human_g1k_v37_decoy.fasta -BQSR $project_bqsr_indir/$prefix.bqsr -I $project_indir/$prefix.indel.$chr.bam -o $project_outdir/$prefix.final.$chr.bam -nct $ncore -L $chr
fi

# indexing
/usr/local/bin/samtools-1.3/bin/samtools index $project_outdir/$prefix.final.$chr.bam

