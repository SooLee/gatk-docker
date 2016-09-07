#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script runs the Merge final bam step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir prefix chr_list_str"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  prefix: prefix of the output bam file. The output file name will be prefix.mkdup.bam."
  echo -e "  chr_list_str: comma-separated list of chromosomes to be merged. (1,2,3,...,MT,decoy,unmapped)"
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
chr_list_str=$4


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi



chr_list=(${chr_list_str//,/ }) # convert to array
bam_list_str=''
for tmpchr in "${chr_list[@]}"
do
  bam_list_str=$project_indir/$prefix.final.$tmpchr.bam' '$bam_list_str
done


### merging final bam (for storage purpose)
/usr/local/bin/samtools-1.3/bin/samtools merge $project_outdir/$prefix.final.bam $bam_list_str 
/usr/local/bin/samtools-1.3/bin/samtools index $project_outdir/$prefix.final.bam


