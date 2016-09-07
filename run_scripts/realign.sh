#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the Indel Realignment step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir prefix chr ncore mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  prefix: prefix of the input/output bam files. The input file name is assumed to be prefix.mkdup.bam and the output file name will be prefix.indel.chr.bam."
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
project_outdir=$2
prefix=$3
chr=$4
ncore=$5
mem=$6


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi



### step6 InDel realignment
if [ $chr == 'decoy' ]; then
  /usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /resources/human_g1k_v37_decoy.fasta -I $project_indir/$prefix.mkdup.bam -o $project_outdir/$prefix.indel.$chr.intervals -XL /resources/b37_nodecoy.bed
  /usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T IndelRealigner -R /resources/human_g1k_v37_decoy.fasta -I $project_indir/$prefix.mkdup.bam -targetIntervals $project_outdir/$prefix.indel.$chr.intervals -o $project_outdir/$prefix.indel.$chr.bam
elif [ $chr == 'unmapped' ]; then
  /usr/local/bin/samtools-1.3/bin/samtools view -bh -o $project_outdir/$prefix.indel.$chr.bam $project_indir/$prefix.mkdup.bam \* 
else # 1,2,...,X,MT
  /usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /resources/human_g1k_v37_decoy.fasta -I $project_indir/$prefix.mkdup.bam -o $project_outdir/$prefix.indel.$chr.intervals -L $chr
  /usr/local/bin/jdk1.7.0_71/bin/java -Xmx$mem -Xms$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T IndelRealigner -R /resources/human_g1k_v37_decoy.fasta -I $project_indir/$prefix.mkdup.bam -targetIntervals $project_outdir/$prefix.indel.$chr.intervals -o $project_outdir/$prefix.indel.$chr.bam
fi

# indexing
/usr/local/bin/samtools-1.3/bin/samtools index $project_outdir/$prefix.indel.$chr.bam

