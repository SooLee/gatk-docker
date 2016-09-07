#!/bin/bash 
display_usage() {
  echo -e ""
  echo -e "This script run the Generate GVCF step of a variant calling pipeline based on bwa-gatk."  
  echo -e ""
  echo -e "Usage:"
  echo -e "  $0 project_indir project_outdir prefix chr region ncore mem"
  echo -e ""
  echo -e "Arguments:"
  echo -e "  project_indir: directory in which input files are placed."
  echo -e "  project_outdir: directory in which output files are placed."
  echo -e "  prefix: prefix of the input/output bam/vcf files. The input file names are assumed to be prefix.final.chr.bam and the output file name will be prefix.region.g.vcf."
  echo -e "  chr: chromosome (eg. 21)."
  echo -e "  region: region (eg. 21:1-50000). The region must match the chromosome."
  echo -e "  ncore: number of CPUs (recommended: 2)"
  echo -e "  mem: memory (recommended: 4G)."
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
region=$5
ncore=$6
mem=$7


if [ ! -d $project_indir ]; then
  echo -e "project input directory doesn't exist."
  exit 1
fi

if [ ! -d $project_outdir ]; then
  mkdir -p $project_outdir
fi


### HaplotypeCaller (Simple example)
### step1. generate g.vcf
java -Xms$mem -Xmx$mem -jar /usr/local/bin/GATK3.5n/GenomeAnalysisTK.jar -T HaplotypeCaller -R /resources/human_g1k_v37_decoy.fasta -I $project_indir/$prefix.final.$chr.bam --dbsnp /resources/dbsnp_138.b37.vcf -o $project_outdir/$prefix.$region.g.vcf --emitRefConfidence GVCF -l INFO -rf BadCigar -nct $ncore -dt NONE -L $region


