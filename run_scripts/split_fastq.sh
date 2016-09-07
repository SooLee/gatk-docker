#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the Split_fastq step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 fastq_gz project_outdir mate"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  fastq_gz: path to the gzipped fastq file (eithe R1 or R2)"
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  mate: 'R1' or 'R2'. Must match the mate of the fastq file."
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


fastq=$1
project_outdir=$2
mate=$3



if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi

zcat $fastq | python /usr/local/bin/run_scripts/ali_split_fq.py -readRno $mate -wd $project_outdir

